"""Pre-registered tests for the confirmation run (docs/preregistration.md) -> results/confirm_tests.tsv. Free.

H1  per-family D = residue-1 accuracy real - shuffled; one-sided family sign-flip
    permutation test (10,000), family-bootstrap 95% CI. Pass: p < 0.001, CI low > 10 points.
H2  Spearman(D, log10 family size), one-sided p < 0.01.
H3  real cuts with all 10 residues right (not in the prefix) vs the count expected if
    positions were independent at the model's own per-position accuracy, per family.
    Pass: observed >= 10, observed / expected > 100, and 0 such runs on shuffled cuts.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from bioseq_mem.probe import K, parse, shuffled

rng = np.random.default_rng(2026)
seqs = pd.read_csv("data/confirm_sequences.tsv", sep="\t").set_index("accession")

rows = []
for path in sorted(Path("results/confirm").glob("*.jsonl")):
    for d in map(json.loads, path.open()):
        if "error" in d or d.get("stop") == "refusal":
            continue
        pred = parse(d.get("text", ""))
        if not pred:
            continue
        s = seqs.loc[d["accession"]]
        full = s.sequence if d["condition"] == "real" else shuffled(s.sequence, d["accession"])
        prefix, t = full[: d["cut"]], d["target"]
        hits = [k < len(pred) and pred[k] == t[k] for k in range(K)]
        rows.append({"model": d["model"], "pfam": d["pfam"], "proteins": s.family_proteins, "condition": d["condition"],
                     "first": hits[0], **{f"p{k}": hits[k] for k in range(K)}, "all10": all(hits) and t not in prefix})
df = pd.DataFrame(rows)

out = []
for model, g in df.groupby("model"):
    fam = g.groupby(["pfam", "condition"])["first"].mean().unstack().dropna()
    D = (fam["real"] - fam["shuffled"]).to_numpy()
    obs = D.mean()
    flips = rng.choice([-1, 1], size=(10_000, len(D)))
    p_perm = ((flips * D).mean(1) >= obs).mean()
    boot = D[rng.integers(0, len(D), (10_000, len(D)))].mean(1)
    lo, hi = np.percentile(boot, [2.5, 97.5])
    size = g.groupby("pfam").proteins.first().reindex(fam.index)
    rho, p_two = spearmanr(D, np.log10(size))
    p_rho = p_two / 2 if rho > 0 else 1 - p_two / 2
    real = g[g.condition == "real"]
    per_fam = real.groupby("pfam")[[f"p{k}" for k in range(K)]].mean()
    n_cuts = real.groupby("pfam").size()
    expected = float((per_fam.prod(axis=1) * n_cuts).sum())
    observed = int(real.all10.sum())
    shuffled_runs = int(g[g.condition == "shuffled"].all10.sum())
    out.append({
        "model": model, "families": len(D), "real_acc": real["first"].mean(), "shuffled_acc": g[g.condition == "shuffled"]["first"].mean(),
        "D_points": obs * 100, "D_ci_low": lo * 100, "D_ci_high": hi * 100, "H1_p_perm": p_perm,
        "H1_pass": p_perm < 0.001 and lo * 100 > 10,
        "H2_rho": rho, "H2_p": p_rho, "H2_pass": rho > 0 and p_rho < 0.01,
        "H3_observed": observed, "H3_expected": expected, "H3_ratio": observed / expected if expected else float("inf"),
        "H3_shuffled_runs": shuffled_runs, "H3_pass": observed >= 10 and (expected == 0 or observed / expected > 100) and shuffled_runs == 0,
    })
res = pd.DataFrame(out)
res.to_csv("results/confirm_tests.tsv", sep="\t", index=False)
pd.set_option("display.width", 250)
print(res.round(4).to_string(index=False))
