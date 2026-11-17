"""Small families, below the main set's 300-protein floor -> data/small_family_sequences.bin{0,1,2}.tsv. Free.

Same selection as scripts/data/04_confirm_sequences.py (DUF families, no structure, at most 3
Europe PMC hits, 4 UniRef90-singleton TrEMBL proteins from Bacteria or Archaea, 5 cuts inside the
Pfam domain), but families of 20 to 300 proteins: 30 per bin of 20-60, 60-150 and 150-300.
"""

import math
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_probe.uniprot import europepmc_hits, pfam_spans, search, uniref90_singletons

BINS = [(20, 60), (60, 150), (150, 300)]
PER_BIN = 30
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
    rng = random.Random(f"small-{fam.pfam}")
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


import sys
from pathlib import Path

# One process per size bin (argument: 0, 1 or 2), each appending a family's rows as soon as it has
# 4 proteins, so a run can stop and resume. Families already tried are skipped; data/small_families.log
# from the first, sequential run counts as tried.
b = int(sys.argv[1])
lo, hi = BINS[b]
out = Path(f"data/small_family_sequences.bin{b}.tsv")
tried = Path(f"data/small_family_tried.bin{b}.txt")
done_fams = set(pd.read_csv(out, sep="\t").pfam) if out.exists() else set()
seen = set(tried.read_text().split()) if tried.exists() else set()
log = Path("data/small_families.log")
if log.exists():
    for line in log.read_text().splitlines():
        parts = line.split()
        if len(parts) == 3 and parts[0].startswith("PF") and parts[1] == f"{lo}-{hi}" and parts[2] == "0":
            seen.add(parts[0])
fams = pd.read_csv("data/pfam_families.tsv", sep="\t")
exclude = set(pd.read_csv("data/families_pilot.tsv", sep="\t").pfam) | set(pd.read_csv("data/confirm_sequences.tsv", sep="\t").pfam)
fams = fams[fams["name"].str.contains("unknown function", case=False, na=False) & (fams.structures == 0) & ~fams.pfam.isin(exclude)]
cand = fams[(fams.proteins >= lo) & (fams.proteins < hi)].sample(frac=1, random_state=b).assign(size_bin=f"{lo}-{hi}")
cand = cand[~cand.pfam.isin(seen | done_fams)]
kept = len(done_fams)
with ThreadPoolExecutor(24) as pool:
    for start in range(0, len(cand), 48):
        if kept >= PER_BIN:
            break
        batch = [f for _, f in cand.iloc[start : start + 48].iterrows()]
        for fam, fam_rows in zip(batch, pool.map(draw, batch)):
            with tried.open("a") as t:
                t.write(fam.pfam + "\n")
            if fam_rows and kept < PER_BIN:
                pd.DataFrame(fam_rows).to_csv(out, sep="\t", index=False, mode="a", header=not out.exists())
                kept += 1
        print(f"bin {lo}-{hi}: {kept}/{PER_BIN}", flush=True)
