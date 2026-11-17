"""Next-residue accuracy per model and family cell, against free baselines. Free.

Baselines, each predicting residue 1 from the prefix alone:
  background  always the most common residue overall (L)
  unigram     the prefix's most common residue
  copy        the residue that followed the longest earlier repeat of the prefix's
              suffix (>= 3 residues), else unigram: what in-context copying alone gets
The main control is the model itself on the shuffled sequence (same composition,
no order): real minus shuffled is what it knows about order.
A model's score counts only valid replies; refusals and errors are reported apart.
"""

import json
import math
from collections import Counter
from pathlib import Path

import pandas as pd

from bioseq_probe.probe import cut_off, parse, score, shuffled

seqs = pd.read_csv("data/sequences_pilot.tsv", sep="\t").set_index("accession")


def copy_guess(prefix: str) -> str:
    for n in range(min(12, len(prefix) - 1), 2, -1):
        suf = prefix[-n:]
        i = prefix.rfind(suf, 0, len(prefix) - 1)
        if i != -1 and i + n < len(prefix):
            return prefix[i + n]
    return Counter(prefix).most_common(1)[0][0]


def first_residue_logprob(tops: list[dict[str, float]], residue: str) -> float | None:
    """P(next residue) from the first output token's top-20, summing tokens that start with it."""
    if not tops:
        return None
    p = sum(math.exp(lp) for tok, lp in tops[0].items() if tok.strip().upper().startswith(residue))
    return math.log(max(p, 1e-6))


rows = []
for path in sorted(Path("results/probe").glob("*.jsonl")):
    for d in map(json.loads, path.open()):
        s = seqs.loc[d["accession"]]
        full = s.sequence if d["condition"] == "real" else shuffled(s.sequence, d["accession"])
        prefix, t = full[: d["cut"]], d["target"][0]
        d.update(score(parse(d.get("text", ""), prefix, cut=cut_off(d)), d["target"]))  # re-scored, so parser fixes apply to all
        rows.append(
            {
                **{k: d[k] for k in ("model", "condition", "accession", "pfam", "cut", "valid", "first_correct", "exact_run")},
                "cell": f"seq-{s.seq_cell}/lit-{s.lit_cell}",
                "era": s.era,
                "first_seen": s.first_seen,
                "outcome": "error" if "error" in d else ("refusal" if d.get("stop") == "refusal" else ("valid" if d["valid"] else "invalid")),
                "logp": first_residue_logprob(d.get("top_logprobs") or [], t),
                "background": t == "L",
                "unigram": t == Counter(prefix).most_common(1)[0][0],
                "copy": t == copy_guess(prefix),
            }
        )
df = pd.DataFrame(rows)
if df.empty:
    raise SystemExit("no results yet")

print(pd.crosstab(df.model, df.outcome), "\n")
v = df[df.outcome == "valid"]
acc = v.groupby(["model", "cell", "condition"]).first_correct.mean().unstack()
acc["real - shuffled"] = acc["real"] - acc["shuffled"]
print("residue-1 accuracy (valid replies)\n", acc.unstack("cell").round(3))
base = df.drop_duplicates(["accession", "cut", "condition"]).groupby(["cell", "condition"])[["background", "unigram", "copy"]].mean()
print("\nbaselines\n", base.round(3))
era = v.groupby(["model", "era", "condition"]).first_correct.mean().unstack()
era["real - shuffled"] = era["real"] - era["shuffled"]
print("\nby era: old = entry before 2023, recent = entry from 2025-06\n", era.unstack("era").round(3))
print("\nmean log P(next residue), logprob models\n", v.dropna(subset=["logp"]).groupby(["model", "cell", "condition"]).logp.mean().unstack().round(2))
print("\nprobes with >= 5 residues exactly right\n", v.assign(hit=v.exact_run >= 5).groupby(["model", "condition"]).hit.sum().unstack())
df.to_csv("results/probe_scored.tsv", sep="\t", index=False)
