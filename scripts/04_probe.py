"""Ask each model for the next residues at every cut point -> results/probe/<model>.jsonl.

Each cut point is asked twice: on the real sequence, and on the same sequence
shuffled (same composition, no order). Real minus shuffled accuracy is what the
model knows about residue order beyond composition.

Dry by default: prints what a run would send. `--run` makes paid calls; resumable,
so finished (accession, cut, condition) probes are skipped.

    uv run python scripts/04_probe.py --models gpt-4.1-2025-04-14 --limit 20 --run
"""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from bioseq_mem.probe import K, parse, prompt, score, shuffled
from bioseq_mem.providers import MODELS, call

ap = argparse.ArgumentParser()
ap.add_argument("--models", default=",".join(MODELS))
ap.add_argument("--limit", type=int, default=None, help="first N probes only")
ap.add_argument("--workers", type=int, default=8)
ap.add_argument("--run", action="store_true", help="make paid API calls")
args = ap.parse_args()

seqs = pd.read_csv("data/sequences_pilot.tsv", sep="\t")
probes = [
    {"accession": r.accession, "pfam": r.pfam, "condition": cond, "cut": int(c), "prefix": s[: int(c)], "target": s[int(c) : int(c) + K]}
    for r in seqs.itertuples()
    for c in r.cuts.split(",")
    for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession)))
]
probes = probes[: args.limit]
chars = sum(len(prompt(p["prefix"])) for p in probes)
print(f"{len(probes)} probes per model, ~{chars // 4:,} input tokens per model (chars/4)")

if not args.run:
    print("dry run; pass --run to call:", args.models)
    raise SystemExit

out_dir = Path("results/probe")
out_dir.mkdir(parents=True, exist_ok=True)
for model in args.models.split(","):
    path = out_dir / f"{model}.jsonl"
    done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
    todo = [p for p in probes if (p["accession"], p["cut"], p["condition"]) not in done]

    def one(p, model=model):
        try:
            r = call(model, prompt(p["prefix"]))
            rec = {"text": r.text, "stop": r.stop, "top_logprobs": r.top_logprobs[:K], "usage": r.usage}
        except Exception as e:  # recorded, not scored as a miss
            rec = {"error": f"{type(e).__name__}: {e}"[:500]}
        pred = parse(rec.get("text", ""))
        return {"model": model, **{k: p[k] for k in ("accession", "pfam", "condition", "cut", "target")}, "pred": pred, **score(pred, p["target"]), **rec}

    with path.open("a") as f, ThreadPoolExecutor(args.workers) as pool:
        for i, rec in enumerate(pool.map(one, todo), 1):
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print(f"{model}: {i}/{len(todo)}", end="\r")
    print()
