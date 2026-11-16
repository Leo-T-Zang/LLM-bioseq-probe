"""Pilot results for every model and prompt (plain and best-guess) -> results/pilot_all.tsv. Free.

Real-minus-shuffled residue-1 accuracy, all families and non-vertebrate DUF
families, with 95% family-bootstrap CIs; valid-reply counts; exact 10-residue runs.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

from bioseq_probe.probe import K, parse

rng = np.random.default_rng(0)
seqs = pd.read_csv("data/sequences_pilot.tsv", sep="\t").set_index("accession")
group = pd.read_csv("data/sequence_groups.tsv", sep="\t").set_index("accession").group


def boot(g):
    t = g.groupby(["pfam", "condition"])["first"].agg(["sum", "count"]).unstack("condition").fillna(0)
    if ("count", "shuffled") not in t or len(t) < 3:
        return np.nan, np.nan, np.nan
    s, c = t["sum"][["real", "shuffled"]].to_numpy(), t["count"][["real", "shuffled"]].to_numpy()
    idx = rng.integers(0, len(t), (4000, len(t)))
    with np.errstate(invalid="ignore", divide="ignore"):
        pt = s[:, 0].sum() / c[:, 0].sum() - s[:, 1].sum() / c[:, 1].sum()
        b = s[idx, 0].sum(1) / c[idx, 0].sum(1) - s[idx, 1].sum(1) / c[idx, 1].sum(1)
    return pt * 100, np.nanpercentile(b, 2.5) * 100, np.nanpercentile(b, 97.5) * 100


out = []
for prompt_kind, d in (("plain", "results/probe"), ("guess", "results/probe_guess")):
    for path in sorted(Path(d).glob("*.jsonl")):
        recs = [json.loads(line) for line in path.open()]
        rows = []
        for r in recs:
            pred = parse(r.get("text", ""))
            if "error" in r or r.get("stop") == "refusal" or not pred:
                continue
            t = r["target"]
            rows.append({"pfam": r["pfam"], "condition": r["condition"], "first": pred[0] == t[0],
                         "all10": pred[:K] == t, "nvduf": group[r["accession"]] != "Vertebrata" and seqs.loc[r["accession"], "lit_cell"] == "poor"})
        g = pd.DataFrame(rows)
        if g.empty:
            continue
        a, nv = boot(g), boot(g[g.nvduf])
        out.append({"model": path.stem, "prompt": prompt_kind, "probes": len(recs), "valid": len(g),
                    "diff_all": a[0], "lo_all": a[1], "hi_all": a[2], "diff_nvduf": nv[0], "lo_nvduf": nv[1], "hi_nvduf": nv[2],
                    "exact10": int(((g.condition == "real") & g.all10).sum())})
res = pd.DataFrame(out)
res.to_csv("results/pilot_all.tsv", sep="\t", index=False)
pd.set_option("display.width", 220)
print(res.round(1).to_string(index=False))
