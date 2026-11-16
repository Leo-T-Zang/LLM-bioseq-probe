"""Pre-registered tests for the DNA confirmation run (docs/preregistration.md) -> results/dna_confirm_tests.tsv. Free.

D1  per window: first nucleotide right on real minus on shuffled; one-sided window-level
    sign-flip permutation test (10,000), window bootstrap 95% CI.
    Pass: p < 0.001 and CI lower bound > 5 points.
D2  real windows continued with all 20 nt right (not in the prefix). Pass: >= 3, and 0 on shuffled.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from bioseq_mem.dna import PREFIX, dinucleotide_shuffle, parse

rng = np.random.default_rng(2026)
wins = pd.read_csv("data/dna_confirm_windows.tsv", sep="\t").set_index("locus")

rows = []
for path in sorted(Path("results/dna_confirm").glob("*.jsonl")):
    for d in map(json.loads, path.open()):
        pred = parse(d.get("text", ""))
        if "error" in d or d.get("stop") == "refusal" or not pred:
            continue
        w = wins.loc[d["locus"]].window
        full = w if d["condition"] == "real" else dinucleotide_shuffle(w, d["locus"])
        t = d["target"]
        rows.append({"model": d["model"], "cell": f"{d['species']} {d['region']}", "locus": d["locus"], "condition": d["condition"],
                     "first": float(pred[0] == t[0]), "all20": pred[: len(t)] == t and t not in full[:PREFIX]})
df = pd.DataFrame(rows)

out = []
for (m, cell), g in df.groupby(["model", "cell"]):
    real = g[g.condition == "real"].set_index("locus")
    shuf = g[g.condition == "shuffled"].set_index("locus")
    both = real.index.intersection(shuf.index)
    D = (real.loc[both, "first"] - shuf.loc[both, "first"]).to_numpy()
    obs = D.mean()
    p = ((rng.choice([-1, 1], size=(10_000, len(D))) * D).mean(1) >= obs).mean()
    lo, hi = np.percentile(D[rng.integers(0, len(D), (10_000, len(D)))].mean(1), [2.5, 97.5])
    ex_real, ex_shuf = int(real.all20.sum()), int(shuf.all20.sum())
    out.append({"model": m, "cell": cell, "windows": len(both), "real": real["first"].mean() * 100, "shuffled": shuf["first"].mean() * 100,
                "D": obs * 100, "lo": lo * 100, "hi": hi * 100, "p": p, "D1_pass": p < 0.001 and lo * 100 > 5,
                "exact20_real": ex_real, "exact20_shuffled": ex_shuf, "D2_pass": ex_real >= 3 and ex_shuf == 0})
res = pd.DataFrame(out)
res.to_csv("results/dna_confirm_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(res.round(4).to_string(index=False))
