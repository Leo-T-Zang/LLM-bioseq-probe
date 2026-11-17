"""The protein probe with no instruction at all: the chat message is only ">sequence" and the
FASTA prefix -> results/no_instruction/<model>.jsonl. Dry by default; --run makes paid calls.

A subset of the main set: 100 proteins drawn with a fixed seed, their first 2 cuts, real and
shuffled (400 questions per model). Scoring: if the reply repeats the end of the prefix (its last
15 residues), the residues after it; otherwise the first run of letters at the start of the reply.
Replies that start with prose count as non-answers.

    uv run python scripts/probe/no_instruction.py --models gpt-5.5,claude-opus-5-5 --run
"""

import argparse
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

from bioseq_probe.probe import STANDARD, fasta, shuffled
from bioseq_probe.providers import call

OUT = Path("results/no_instruction")


def parse(text: str, prefix: str) -> str:
    letters = re.sub(r"[^A-Z]", "", re.sub(r"<[^>]*>", "", text.upper().replace(">SEQUENCE", "")))
    i = letters.rfind(prefix[-15:])
    if i != -1:
        cont = letters[i + 15 :]
    else:
        first = re.sub(r"[\s>]", "", text.strip().split("\n\n")[0])
        cont = first if re.fullmatch(r"[A-Za-z]+", first or "") else ""
    cont = cont.upper()[:10]
    return cont if cont and set(cont) <= set(STANDARD) else ""


ap = argparse.ArgumentParser()
ap.add_argument("--models", required=True)
ap.add_argument("--run", action="store_true")
args = ap.parse_args()
seqs = pd.read_csv("data/confirm_sequences.tsv", sep="\t").sample(100, random_state=2026)
todo = [{"accession": r.accession, "pfam": r.pfam, "cut": int(c), "condition": cond, "prefix": s[: int(c)], "target": s[int(c) : int(c) + 10]}
        for r in seqs.itertuples() for c in r.cuts.split(",")[:2]
        for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession)))]
OUT.mkdir(parents=True, exist_ok=True)
for model in args.models.split(","):
    print(model, len(todo), "questions")
    if not args.run:
        continue

    def one(p, model=model):
        try:
            r = call(model, f">sequence\n{fasta(p['prefix'])}")
            rec = {"text": r.text, "stop": r.stop, "usage": r.usage}
        except Exception as e:
            rec = {"error": f"{type(e).__name__}: {e}"[:300]}
        pred = parse(rec.get("text", ""), p["prefix"])
        return {"model": model, **{k: p[k] for k in ("accession", "pfam", "cut", "condition", "target")}, "pred": pred,
                "valid": bool(pred), "first_correct": bool(pred) and pred[0] == p["target"][0], **rec}

    with (OUT / f"{model}.jsonl").open("w") as f, ThreadPoolExecutor(16) as pool:
        for rec in pool.map(one, todo):
            f.write(json.dumps(rec) + "\n")
