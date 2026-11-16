"""Nucleotide probes: next 20 nt after a 200-nt prefix, real and dinucleotide-shuffled
-> results/dna/<model>.jsonl. Dry by default; --run makes paid calls; resumable.

    uv run python scripts/probe/dna.py --models claude-opus-5-5 --run
"""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from bioseq_probe.dna import K, PREFIX, dinucleotide_shuffle, parse, prompt
from bioseq_probe.providers import MODELS, call

ap = argparse.ArgumentParser()
ap.add_argument("--models", default=",".join(MODELS))
ap.add_argument("--workers", type=int, default=16)
ap.add_argument("--run", action="store_true")
ap.add_argument("--guess", action="store_true", help="add the best-guess instruction")
ap.add_argument("--out", default="results/dna")
ap.add_argument("--windows", default="data/dna_windows_*.tsv", help="glob of window files")
args = ap.parse_args()

wins = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(Path().glob(args.windows))])
probes = []
for w in wins.itertuples():
    for cond, s in (("real", w.window), ("shuffled", dinucleotide_shuffle(w.window, w.locus))):
        probes.append({"species": w.species, "region": w.region, "locus": w.locus, "condition": cond, "prefix": s[:PREFIX], "target": s[PREFIX : PREFIX + K]})
print(f"{len(probes)} probes per model")
if not args.run:
    raise SystemExit("dry run; pass --run")

out = Path(args.out)
out.mkdir(parents=True, exist_ok=True)
for model in args.models.split(","):
    path = out / f"{model}.jsonl"
    done = {(d["locus"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
    todo = [p for p in probes if (p["locus"], p["condition"]) not in done]

    def one(p, model=model):
        try:
            r = call(model, prompt(p["prefix"], args.guess))
            rec = {"text": r.text, "stop": r.stop, "usage": r.usage}
        except Exception as e:
            rec = {"error": f"{type(e).__name__}: {e}"[:500]}
        return {"model": model, **{k: p[k] for k in ("species", "region", "locus", "condition", "target")}, "pred": parse(rec.get("text", "")), **rec}

    with path.open("a") as f, ThreadPoolExecutor(args.workers) as pool:
        for rec in pool.map(one, todo):
            f.write(json.dumps(rec) + "\n")
            f.flush()
    print(model, "done", flush=True)
