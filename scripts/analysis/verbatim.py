"""Exact 10-residue continuations (pilot) per model and family group, against the
count expected if each position were guessed independently at the model's own
per-position accuracy (per family). -> results/verbatim.tsv. Free."""

import json
from pathlib import Path

import pandas as pd

from bioseq_probe.probe import K, cut_off, parse, shuffled

seqs = pd.read_csv("data/sequences_pilot.tsv", sep="\t").set_index("accession")
group = pd.read_csv("data/sequence_groups.tsv", sep="\t").set_index("accession").group
rows = []
for path in sorted(Path("results/probe").glob("*.jsonl")):
    if path.stem == "gpt-6-astra":
        continue
    for d in map(json.loads, path.open()):
        if "error" in d or d.get("stop") == "refusal":
            continue
        s = seqs.loc[d["accession"]]
        full = s.sequence if d["condition"] == "real" else shuffled(s.sequence, d["accession"])
        pred = parse(d.get("text", ""), full[: d["cut"]], cut=cut_off(d))
        if not pred:
            continue
        t = d["target"]
        hits = [k < len(pred) and pred[k] == t[k] for k in range(K)]
        vert = "vertebrate" if group[d["accession"]] == "Vertebrata" else "non-vertebrate"
        kind = "DUF" if s.lit_cell == "poor" else "well-studied"
        rows.append({"model": d["model"], "pfam": d["pfam"], "condition": d["condition"], "group": f"{vert} {kind}",
                     **{f"p{k}": hits[k] for k in range(K)}, "all10": all(hits) and t not in full[: d["cut"]]})
df = pd.DataFrame(rows)
out = []
for (m, g), x in df[df.condition == "real"].groupby(["model", "group"]):
    per_fam = x.groupby("pfam")[[f"p{k}" for k in range(K)]].mean()
    expected = float((per_fam.prod(axis=1) * x.groupby("pfam").size()).sum())
    sh = df[(df.model == m) & (df.group == g) & (df.condition == "shuffled")].all10.sum()
    out.append({"model": m, "group": g, "cuts": len(x), "families": x.pfam.nunique(), "observed": int(x.all10.sum()),
                "expected": expected, "shuffled_observed": int(sh)})
res = pd.DataFrame(out)
res.to_csv("results/verbatim.tsv", sep="\t", index=False)
pd.set_option("display.width", 200)
print(res.pivot(index="model", columns="group", values="observed").to_string())
print(res.pivot(index="model", columns="group", values="expected").round(3).to_string())
print("shuffled exact runs, all models:", res.shuffled_observed.sum())
