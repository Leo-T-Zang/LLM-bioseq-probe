"""Draw test sequences for each pilot family -> data/sequences_pilot.tsv. Free.

Per family, unreviewed (TrEMBL) entries from two eras: OLD entries created
before 2023-01-01, which every tested model could have seen, and RECENT entries
created on or after 2025-06-01. Each sequence's UniParc first-seen date is kept,
so each model can be scored separately on sequences before and after its own
cutoff: higher accuracy only before the cutoff = recall; the same on both = learned
regularities. A sequence qualifies when it uses only the 20 standard amino acids
and is not low-complexity; cut points fall inside the Pfam domain. Its UniRef90
cluster size (relatives at >= 90% identity) is recorded, not filtered on.
"""

import math
import random
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_probe.uniprot import pfam_spans, search, uniparc_first_seen, uniref90_size

ERAS = {"old": "date_created:[* TO 2022-12-31]", "recent": "date_created:[2025-06-01 TO *]"}
PER_ERA = 3
FAMILIES_PER_CELL = 10
CUTS_PER_SEQ = 5
MIN_PREFIX = 30
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def entropy(s: str) -> float:
    c = Counter(s)
    return -sum(v / len(s) * math.log2(v / len(s)) for v in c.values())


def draw(fam) -> list[dict]:
    rng = random.Random(fam.pfam)
    rows = []
    for era, window in ERAS.items():
        entries = search(
            f"xref:pfam-{fam.pfam} AND reviewed:false AND {window} AND length:[120 TO 800]",
            "accession,length,sequence,organism_name",
            size=200,
        )
        entries = [e for e in entries if set(e["sequence"]["value"]) <= STANDARD and entropy(e["sequence"]["value"]) > 3.9]
        rng.shuffle(entries)
        kept, md5s = 0, set()
        for e in entries:
            s = e["sequence"]["value"]
            if e["sequence"]["md5"] in md5s:
                continue
            spans = [(a, b) for a, b in pfam_spans(e["primaryAccession"], fam.pfam) if b - a >= 40]
            if not spans:
                continue
            a, b = max(spans, key=lambda x: x[1] - x[0])
            lo = max(a - 1 + 10, MIN_PREFIX)
            if b - 10 - lo < CUTS_PER_SEQ:
                continue
            cuts = sorted(rng.sample(range(lo, b - 10), k=CUTS_PER_SEQ))
            md5s.add(e["sequence"]["md5"])
            rows.append(
                {
                    "pfam": fam.pfam, "seq_cell": fam.seq, "lit_cell": fam.lit, "era": era,
                    "accession": e["primaryAccession"], "upi": e["extraAttributes"]["uniParcId"],
                    "organism": e["organism"]["scientificName"], "length": len(s),
                    "uniref90_size": uniref90_size(e["primaryAccession"]),
                    "domain": f"{a}-{b}", "cuts": ",".join(map(str, cuts)), "sequence": s,
                }
            )
            kept += 1
            if kept == PER_ERA:
                break
    seen = uniparc_first_seen([r["upi"] for r in rows])
    for r in rows:
        r["first_seen"] = seen.get(r["upi"])
    print(fam.pfam, fam.seq, fam.lit, len(rows), flush=True)
    return rows


fams = pd.read_csv("data/families_pilot.tsv", sep="\t")
with ThreadPoolExecutor(4) as pool:
    drawn = list(pool.map(draw, [f for _, f in fams.iterrows()]))
per_cell = Counter()
rows = []
for fam_rows in drawn:  # families keep 02's random order
    if len(fam_rows) < 2 * PER_ERA:
        continue
    cell = (fam_rows[0]["seq_cell"], fam_rows[0]["lit_cell"])
    if per_cell[cell] < FAMILIES_PER_CELL:
        per_cell[cell] += 1
        rows += fam_rows
print(dict(per_cell))
pd.DataFrame(rows).to_csv("data/sequences_pilot.tsv", sep="\t", index=False)
print(len(rows), "sequences")
