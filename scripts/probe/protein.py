"""Ask each model for the next residues at every cut point -> results/probe/<model>.jsonl.

Each cut point is asked twice: on the real sequence, and on the same sequence
shuffled (same composition, no order). Real minus shuffled accuracy is what the
model knows about residue order beyond composition.

Dry by default: prints what a run would send. `--run` makes paid calls; resumable,
so finished (accession, cut, condition) probes are skipped.

    uv run python scripts/probe/protein.py --models gpt-4.1-2025-04-14 --limit 20 --run
"""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from bioseq_probe.probe import K, cut_off, parse, parse_echo, prompt, prompt_echo, score, shuffled
from bioseq_probe.providers import MODELS, call

ap = argparse.ArgumentParser()
ap.add_argument("--models", default=",".join(MODELS))
ap.add_argument("--limit", type=int, default=None, help="first N probes only")
ap.add_argument("--workers", type=int, default=8)
ap.add_argument("--run", action="store_true", help="make paid API calls")
ap.add_argument("--sequences", default="data/sequences_pilot.tsv")
ap.add_argument("--out", default="results/probe")
ap.add_argument("--guess", action="store_true", help="add the best-guess instruction (bioseq_probe.probe.GUESS)")
ap.add_argument("--max-tokens", type=int, default=4000)
ap.add_argument("--echo", action="store_true", help="copy-then-continue prompt (bioseq_probe.probe.prompt_echo)")
args = ap.parse_args()

seqs = pd.read_csv(args.sequences, sep="\t")
probes = [
    {"accession": r.accession, "pfam": r.pfam, "condition": cond, "cut": int(c), "prefix": s[: int(c)], "target": s[int(c) : int(c) + K]}
    for r in seqs.itertuples()
    for c in r.cuts.split(",")
    for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession)))
]
probes = probes[: args.limit]
chars = sum(len(prompt(p["prefix"], args.guess)) for p in probes)
print(f"{len(probes)} probes per model, ~{chars // 4:,} input tokens per model (chars/4)")

if not args.run:
    print("dry run; pass --run to call:", args.models)
    raise SystemExit

out_dir = Path(args.out)
out_dir.mkdir(parents=True, exist_ok=True)
for model in args.models.split(","):
    path = out_dir / f"{model}.jsonl"
    done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
    todo = [p for p in probes if (p["accession"], p["cut"], p["condition"]) not in done]

    def one(p, model=model):
        try:
            text = prompt_echo(p["prefix"]) if args.echo else prompt(p["prefix"], args.guess)
            r = call(model, text, args.max_tokens)
            rec = {"text": r.text, "stop": r.stop, "top_logprobs": r.top_logprobs[:K], "usage": r.usage}
        except Exception as e:  # recorded, not scored as a miss
            rec = {"error": f"{type(e).__name__}: {e}"[:500]}
        pred = parse_echo(rec.get("text", ""), p["prefix"]) if args.echo else parse(rec.get("text", ""), p["prefix"], cut=cut_off(rec))
        return {"model": model, **{k: p[k] for k in ("accession", "pfam", "condition", "cut", "target")}, "pred": pred, **score(pred, p["target"]), **rec}

    with path.open("a") as f, ThreadPoolExecutor(args.workers) as pool:
        for i, rec in enumerate(pool.map(one, todo), 1):
            f.write(json.dumps(rec) + "\n")
            f.flush()
            print(f"{model}: {i}/{len(todo)}", end="\r")
    print()
