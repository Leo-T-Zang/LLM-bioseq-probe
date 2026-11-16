"""Nucleotide probe results -> results/dna_summary.tsv. Free.

Per model, species and region: first-nucleotide accuracy on real vs dinucleotide-
shuffled windows (bootstrap CI over loci), and exact 20-nt continuations (absent
from the prefix). In-context baseline: an order-5 Markov guess from the prefix itself
(most frequent base after the last 5 nt, backing off to shorter contexts).
"""

import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from bioseq_mem.dna import PREFIX, dinucleotide_shuffle, parse

rng = np.random.default_rng(0)
wins = pd.concat([pd.read_csv(p, sep="\t") for p in sorted(Path("data").glob("dna_windows_*.tsv"))]).set_index("locus")


def markov_guess(prefix: str) -> str:
    for k in range(5, -1, -1):
        ctx = prefix[len(prefix) - k :] if k else ""
        nxt = Counter(prefix[i + k] for i in range(len(prefix) - k) if prefix[i : i + k] == ctx)
        if nxt:
            return nxt.most_common(1)[0][0]
    return "A"


rows = []
for path in sorted(p for d in ("results/dna", "results/dna_guess") for p in Path(d).glob("*.jsonl")):
    if path.parent.name == "dna" and path.stem == "deepseek-v4p1-flash":
        continue  # declined almost every plain prompt; its best-guess run is in dna_guess
    for d in map(json.loads, path.open()):
        pred = parse(d.get("text", ""))
        if "error" in d or d.get("stop") == "refusal" or not pred:
            continue
        w = wins.loc[d["locus"]].window
        full = w if d["condition"] == "real" else dinucleotide_shuffle(w, d["locus"])
        prefix, t = full[:PREFIX], d["target"]
        rows.append({"model": d["model"], "species": d["species"], "region": d["region"], "locus": d["locus"], "condition": d["condition"],
                     "first": pred[0] == t[0], "all20": pred[: len(t)] == t and t not in prefix, "markov": markov_guess(prefix) == t[0]})
df = pd.DataFrame(rows)
out = []
for (m, sp, rg), g in df.groupby(["model", "species", "region"]):
    real = g[g.condition == "real"].set_index("locus")["first"]
    shuf = g[g.condition == "shuffled"].set_index("locus")["first"]
    both = real.index.intersection(shuf.index)
    diff = (real[both].astype(float) - shuf[both].astype(float)).to_numpy()
    boot = diff[rng.integers(0, len(diff), (4000, len(diff)))].mean(1)
    out.append({"model": m, "species": sp, "region": rg, "loci": len(both), "real": real.mean() * 100, "shuffled": shuf.mean() * 100,
                "diff": diff.mean() * 100, "lo": np.percentile(boot, 2.5) * 100, "hi": np.percentile(boot, 97.5) * 100,
                "markov_real": g[g.condition == "real"].markov.mean() * 100,
                "exact20_real": int(g[g.condition == "real"].all20.sum()), "exact20_shuffled": int(g[g.condition == "shuffled"].all20.sum())})
res = pd.DataFrame(out)
res.to_csv("results/dna_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(res.round(1).to_string(index=False))
