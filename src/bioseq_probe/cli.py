"""Score a new model on the protein test: `uv run probe-model --model <id> [--provider ...] --run`.

Asks the model for the next amino acids at every cut of the confirmation set (668 obscure
bacterial and archaeal proteins, 167 DUF families, 5 cuts each, real and shuffled: 6,680
questions), then reports the real-minus-shuffled gap in the next amino acid with a 95%
bootstrap CI over families, next to the models already measured. Resumable; dry by default.
"""

import argparse
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from bioseq_probe.probe import K, parse, prompt, score, shuffled
from bioseq_probe.providers import MODELS, call

PROVIDERS = ("anthropic", "openai_chat", "openai_responses", "gemini", "fireworks")


def probes(seqs: pd.DataFrame) -> list[dict]:
    return [
        {"accession": r.accession, "pfam": r.pfam, "condition": cond, "cut": int(c), "prefix": s[: int(c)], "target": s[int(c) : int(c) + K]}
        for r in seqs.itertuples()
        for c in r.cuts.split(",")
        for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession)))
    ]


def family_gap(records: list[dict], rng=np.random.default_rng(0)) -> tuple[float, float, float, int]:
    """Real-minus-shuffled next-residue accuracy (points), 95% bootstrap CI over families, valid replies."""
    rows = [{"pfam": d["pfam"], "condition": d["condition"], "hit": p[0] == d["target"][0]}
            for d in records if "error" not in d and d.get("stop") != "refusal" and (p := parse(d.get("text", "")))]
    fam = pd.DataFrame(rows).groupby(["pfam", "condition"]).hit.mean().unstack().dropna()
    gap = (fam["real"] - fam["shuffled"]).to_numpy() * 100
    boot = gap[rng.integers(0, len(gap), (10_000, len(gap)))].mean(1)
    return gap.mean(), *np.percentile(boot, [2.5, 97.5]), len(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, help="model id as the provider's API names it")
    ap.add_argument("--provider", choices=PROVIDERS, help="needed only for a model not in bioseq_probe.providers.MODELS")
    ap.add_argument("--limit", type=int, help="first N questions only, for a quick look")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--run", action="store_true", help="make the paid API calls")
    args = ap.parse_args()
    if args.model not in MODELS:
        if not args.provider:
            ap.error(f"{args.model} is not in MODELS; pass --provider")
        MODELS[args.model] = (args.provider, {})

    todo = probes(pd.read_csv("data/confirm_sequences.tsv", sep="\t"))[: args.limit]
    out = Path("results/confirm") / f"{args.model}.jsonl"
    done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, out.open())} if out.exists() else set()
    todo = [p for p in todo if (p["accession"], p["cut"], p["condition"]) not in done]
    print(f"{args.model}: {len(todo)} questions to ask ({len(done)} already done)")
    if todo and not args.run:
        raise SystemExit("dry run; pass --run to make the paid API calls")

    def one(p):
        try:
            r = call(args.model, prompt(p["prefix"]))
            rec = {"text": r.text, "stop": r.stop, "usage": r.usage}
        except Exception as e:
            rec = {"error": f"{type(e).__name__}: {e}"[:500]}
        pred = parse(rec.get("text", ""))
        return {"model": args.model, **{k: p[k] for k in ("accession", "pfam", "condition", "cut", "target")}, "pred": pred, **score(pred, p["target"]), **rec}

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("a") as f, ThreadPoolExecutor(args.workers) as pool:
        for rec in pool.map(one, todo):
            f.write(json.dumps(rec) + "\n")

    gap, lo, hi, n = family_gap([json.loads(line) for line in out.open()])
    print(f"\n{args.model}: next amino acid, real minus shuffled = {gap:+.1f} points (95% CI {lo:+.1f} to {hi:+.1f}), {n} valid replies")
    print("for reference (results/confirm_tests.tsv): Claude Opus 5.5 +27.4, Opus 5 +9.7, Opus 4.7 +8.9, Opus 4.8 +7.7, GPT-6 Astra +2.9, GPT-5.5 +2.1, Gemini 4 Carbon +0.9")
