"""Confirmation run sample, as pre-registered in docs/preregistration.md -> data/confirm_*.tsv. Free.

200 bacterial/archaeal-member DUF families (50 per size bin), 4 UniRef90-singleton
TrEMBL sequences each from Bacteria or Archaea, 5 cuts inside the Pfam domain.
"""

import math
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_mem.uniprot import europepmc_hits, pfam_spans, search, uniref90_singletons

BINS = [(300, 1_000), (1_000, 3_000), (3_000, 10_000), (10_000, 10**9)]
PER_BIN = 50
PER_FAMILY = 4
CUTS_PER_SEQ = 5
MIN_PREFIX = 30
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")
PROKARYOTES = "(taxonomy_id:2 OR taxonomy_id:2157)"


def entropy(s: str) -> float:
    c = Counter(s)
    return -sum(v / len(s) * math.log2(v / len(s)) for v in c.values())


def draw(fam) -> list[dict]:
    if europepmc_hits(f'"{fam.pfam}"') > 3:
        return []
    rng = random.Random(f"confirm-{fam.pfam}")
    entries = search(
        f"xref:pfam-{fam.pfam} AND reviewed:false AND {PROKARYOTES} AND length:[120 TO 800]",
        "accession,length,sequence,organism_name,lineage",
        size=120,
    )
    entries = [e for e in entries if set(e["sequence"]["value"]) <= STANDARD and entropy(e["sequence"]["value"]) > 3.9]
    rng.shuffle(entries)
    entries = entries[:60]
    single = uniref90_singletons([e["primaryAccession"] for e in entries])
    rows, md5s = [], set()
    for e in entries:
        s = e["sequence"]["value"]
        if e["sequence"]["md5"] in md5s or e["primaryAccession"] not in single:
            continue
        spans = [(a, b) for a, b in pfam_spans(e["primaryAccession"], fam.pfam) if b - a >= 40]
        if not spans:
            continue
        a, b = max(spans, key=lambda x: x[1] - x[0])
        lo = max(a - 1 + 10, MIN_PREFIX)
        if b - 10 - lo < CUTS_PER_SEQ:
            continue
        md5s.add(e["sequence"]["md5"])
        rows.append(
            {
                "pfam": fam.pfam, "family_proteins": fam.proteins, "size_bin": fam.size_bin,
                "accession": e["primaryAccession"], "upi": e["extraAttributes"]["uniParcId"],
                "organism": e["organism"]["scientificName"], "kingdom": e["organism"]["lineage"][0],
                "length": len(s), "domain": f"{a}-{b}",
                "cuts": ",".join(map(str, sorted(rng.sample(range(lo, b - 10), k=CUTS_PER_SEQ)))),
                "sequence": s,
            }
        )
        if len(rows) == PER_FAMILY:
            break
    print(fam.pfam, fam.size_bin, len(rows), flush=True)
    return rows if len(rows) == PER_FAMILY else []


fams = pd.read_csv("data/pfam_families.tsv", sep="\t")
pilot = set(pd.read_csv("data/families_pilot.tsv", sep="\t").pfam)
fams = fams[fams["name"].str.contains("unknown function", case=False, na=False) & (fams.structures == 0) & ~fams.pfam.isin(pilot)]
rows = []
with ThreadPoolExecutor(16) as pool:
    for i, (lo, hi) in enumerate(BINS):
        cand = fams[(fams.proteins >= lo) & (fams.proteins < hi)].sample(frac=1, random_state=i).assign(size_bin=f"{lo}-{hi if hi < 10**9 else 'max'}")
        kept = 0
        for start in range(0, len(cand), 32):
            batch = [f for _, f in cand.iloc[start : start + 32].iterrows()]
            for fam_rows in pool.map(draw, batch):
                if fam_rows and kept < PER_BIN:
                    rows += fam_rows
                    kept += 1
            print(f"bin {lo}-{hi}: {kept}/{PER_BIN}", flush=True)
            if kept >= PER_BIN:
                break

df = pd.DataFrame(rows)
df.to_csv("data/confirm_sequences.tsv", sep="\t", index=False)
print(df.groupby("size_bin").pfam.nunique(), len(df), "sequences")
