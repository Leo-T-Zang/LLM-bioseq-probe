"""How far can a model recite a protein letter for letter? -> results/recall_length.jsonl

Common case set: every pilot cut (real sequence) that any model continued with all 10 residues
right, the 10 not present in the prefix. Every model gets the same cases and is asked for the
next 200 residues; the score is the exact run against the true sequence.

    uv run python scripts/23_recall_length.py claude-opus-5-5,gemini-4-carbon --run
"""

import json
import re
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from bioseq_mem import probe
from bioseq_mem.providers import call

N_AA = 200
STANDARD = "ACDEFGHIKLMNPQRSTVWY"
OUT = Path("results/recall_length.jsonl")


def prompt(prefix: str) -> str:
    return probe.prompt(prefix).replace(f"exactly the next {probe.K} amino acids", f"the next {N_AA} amino acids (or fewer if the protein ends)")


seqs = pd.read_csv("data/sequences_pilot.tsv", sep="\t").set_index("accession")
cases = {}
for path in sorted(Path("results/probe").glob("*.jsonl")):
    for d in map(json.loads, path.open()):
        if d["condition"] != "real" or probe.parse(d.get("text", ""))[:10] != d["target"]:
            continue
        s = seqs.loc[d["accession"]].sequence
        if d["target"] not in s[: d["cut"]]:
            cases[(d["accession"], d["cut"])] = {"id": f"{d['accession']}:{d['cut']}", "prefix": s[: d["cut"]], "truth": s[d["cut"] :]}
cases = list(cases.values())
print(len(cases), "common cases")
if "--run" not in sys.argv:
    raise SystemExit("dry run; pass --run")

done = {(r["model"], r["id"]) for r in map(json.loads, OUT.open())} if OUT.exists() else set()
for model in sys.argv[1].split(","):
    todo = [c for c in cases if (model, c["id"]) not in done]

    def one(c, model=model):
        try:
            text = call(model, prompt(c["prefix"]), 8000).text
        except Exception as e:
            return {"model": model, "id": c["id"], "error": f"{type(e).__name__}: {e}"[:300]}
        pred = re.sub(f"[^{STANDARD}]", "", re.sub(r"<[^>]*>", "", text).upper().split(">")[0])
        run = 0
        while run < min(len(pred), len(c["truth"])) and pred[run] == c["truth"][run]:
            run += 1
        return {"model": model, "id": c["id"], "truth_available": len(c["truth"]), "written": len(pred), "exact_run": run, "text": text[:2000]}

    with ThreadPoolExecutor(24) as pool, OUT.open("a") as f:
        for r in pool.map(one, todo):
            f.write(json.dumps(r) + "\n")
    print(model, "done", flush=True)
