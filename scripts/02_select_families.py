"""Pick the pilot families: a 2x2 of sequence abundance x literature attention.

Sequence-rich: >= 5,000 proteins; sequence-poor: 300-1,500. Only 9 DUFs reach
20,000 proteins with no structure, so the rich threshold is 5,000.

Literature attention: a DUF ("unknown function") family with no PDB structure is
literature-poor; a named family with >= 100 structures is literature-rich. Europe PMC
hits for the Pfam accession are recorded as a check. Every family must have enough
UniProtKB entries created after CUTOFF to draw test sequences from. Free.
"""

from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_mem.uniprot import europepmc_hits

PER_CELL = 16  # spare families; 03 keeps the first 10 per cell that yield enough sequences
SEED = 0

df = pd.read_csv("data/pfam_families.tsv", sep="\t")
df["duf"] = df["name"].str.contains("unknown function", case=False, na=False)
df["seq"] = None
df.loc[df["proteins"].between(300, 1_500), "seq"] = "poor"
df.loc[df["proteins"] >= 5_000, "seq"] = "rich"
lit_poor = df["duf"] & (df["structures"] == 0)
lit_rich = ~df["duf"] & (df["structures"] >= 100)
df["lit"] = None
df.loc[lit_poor, "lit"] = "poor"
df.loc[lit_rich, "lit"] = "rich"
cands = df.dropna(subset=["seq", "lit"])
print(cands.groupby(["seq", "lit"]).size())


def check(row):
    hits = europepmc_hits(f'"{row.pfam}"')
    if row.lit == "poor" and hits > 3:
        return None
    return {**row.to_dict(), "epmc_hits": hits}


picked = []
with ThreadPoolExecutor(4) as pool:
    for (seq, lit), group in cands.groupby(["seq", "lit"]):
        rows = [r for _, r in group.sample(frac=1, random_state=SEED).iterrows()]
        got = []
        for i in range(0, len(rows), 24):
            got += [r for r in pool.map(check, rows[i : i + 24]) if r]
            print(seq, lit, f"checked {min(i + 24, len(rows))}, kept {len(got)}", flush=True)
            if len(got) >= PER_CELL:
                break
        picked += got[:PER_CELL]

pd.DataFrame(picked).to_csv("data/families_pilot.tsv", sep="\t", index=False)
print(pd.DataFrame(picked)[["seq", "lit", "pfam", "name", "proteins", "structures", "epmc_hits"]].to_string())
