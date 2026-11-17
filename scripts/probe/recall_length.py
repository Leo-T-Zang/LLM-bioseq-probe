"""How far can a model recite a protein letter for letter? -> results/recall_length/<model>.jsonl

Common case set: every pilot cut (real sequence) that any model continued with all 10 residues
right, the 10 not present in the prefix. Every model gets the same cases and is asked for the
next 200 residues; the score is the exact run against the true sequence.
One output file per model: results/recall_length/<model>.jsonl.

    uv run python scripts/probe/recall_length.py claude-opus-5-5,gemini-4-carbon --run
"""

import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from bioseq_probe import probe
from bioseq_probe.providers import call

N_AA = 200
OUT_DIR = Path("results/recall_length")


def prompt(prefix: str) -> str:
    return probe.prompt(prefix).replace(f"exactly the next {probe.K} amino acids", f"the next {N_AA} amino acids (or fewer if the protein ends)")


seqs = pd.read_csv("data/sequences_pilot.tsv", sep="\t").set_index("accession")
cases = {}
for path in sorted(Path("results/probe").glob("*.jsonl")):
    for d in map(json.loads, path.open()):
        s = seqs.loc[d["accession"]].sequence
        if d["condition"] != "real" or probe.parse(d.get("text", ""), s[: d["cut"]], cut=probe.cut_off(d))[:10] != d["target"]:
            continue
        if d["target"] not in s[: d["cut"]]:
            cases[(d["accession"], d["cut"])] = {"id": f"{d['accession']}:{d['cut']}", "prefix": s[: d["cut"]], "truth": s[d["cut"] :]}
cases = list(cases.values())
print(len(cases), "common cases")
if "--run" not in sys.argv:
    raise SystemExit("dry run; pass --run")

OUT_DIR.mkdir(parents=True, exist_ok=True)
for model in sys.argv[1].split(","):
    out = OUT_DIR / f"{model}.jsonl"
    done = {r["id"] for r in map(json.loads, out.open())} if out.exists() else set()
    todo = [c for c in cases if c["id"] not in done]

    def one(c, model=model):
        try:
            r = call(model, prompt(c["prefix"]), 8000)
            text = r.text
        except Exception as e:
            return {"model": model, "id": c["id"], "error": f"{type(e).__name__}: {e}"[:300]}
        pred = probe.parse(text, c["prefix"], prefer_len=20, cut=probe.cut_off({"stop": r.stop}))
        run = 0
        while run < min(len(pred), len(c["truth"])) and pred[run] == c["truth"][run]:
            run += 1
        return {"model": model, "id": c["id"], "truth_available": len(c["truth"]), "written": len(pred), "exact_run": run, "stop": r.stop, "text": text}

    with ThreadPoolExecutor(24) as pool, out.open("a") as f:
        for r in pool.map(one, todo):
            f.write(json.dumps(r) + "\n")
    print(model, "done", flush=True)
