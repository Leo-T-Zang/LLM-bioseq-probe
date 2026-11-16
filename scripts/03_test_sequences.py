"""Draw post-cutoff test sequences for each pilot family -> data/sequences_pilot.tsv. Free.

A sequence qualifies when UniParc first saw it on or after CUTOFF (a UniProt entry
created after the cutoff can carry a sequence that is decades old), it has no
>= 90% identical relative (a UniRef90 singleton), it uses only the 20 standard
amino acids and is not low-complexity. Cut points fall inside the Pfam domain.
"""

import math
import random
from collections import Counter

from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_mem.uniprot import pfam_spans, search, uniparc_first_seen, uniref90_size

CUTOFF = "2026-07-01"
PER_FAMILY = 5
CUTS_PER_SEQ = 5
MIN_PREFIX = 30
STANDARD = set("ACDEFGHIKLMNPQRSTVWY")


def entropy(s: str) -> float:
    c = Counter(s)
    return -sum(v / len(s) * math.log2(v / len(s)) for v in c.values())


def draw(fam) -> list[dict]:
    rng = random.Random(fam.pfam)
    rows = []
    entries = search(
        f"xref:pfam-{fam.pfam} AND date_created:[{CUTOFF} TO *] AND reviewed:false AND length:[120 TO 800]",
        "accession,length,sequence,organism_name",
        size=200,
    )
    entries = [e for e in entries if set(e["sequence"]["value"]) <= STANDARD and entropy(e["sequence"]["value"]) > 3.9]
    seen = uniparc_first_seen([e["extraAttributes"]["uniParcId"] for e in entries])
    rng.shuffle(entries)
    kept, md5s = 0, set()
    for e in entries:
        upi, s = e["extraAttributes"]["uniParcId"], e["sequence"]["value"]
        if seen.get(upi, "0") < CUTOFF or e["sequence"]["md5"] in md5s:
            continue
        if uniref90_size(e["primaryAccession"]) != 1:
            continue
        spans = [(a, b) for a, b in pfam_spans(e["primaryAccession"], fam.pfam) if b - a >= 40]
        if not spans:
            continue
        a, b = max(spans, key=lambda x: x[1] - x[0])
        lo = max(a - 1 + 10, MIN_PREFIX)
        cuts = sorted(rng.sample(range(lo, b - 10), k=min(CUTS_PER_SEQ, b - 10 - lo)))
        md5s.add(e["sequence"]["md5"])
        rows.append(
            {
                "pfam": fam.pfam, "seq_cell": fam.seq, "lit_cell": fam.lit,
                "accession": e["primaryAccession"], "upi": upi, "first_seen": seen[upi],
                "organism": e["organism"]["scientificName"], "length": len(s),
                "domain": f"{a}-{b}", "cuts": ",".join(map(str, cuts)), "sequence": s,
            }
        )
        kept += 1
        if kept == PER_FAMILY:
            break
    print(fam.pfam, fam.seq, fam.lit, f"{kept}/{PER_FAMILY} from {len(entries)}", flush=True)
    return rows


fams = pd.read_csv("data/families_pilot.tsv", sep="\t")
with ThreadPoolExecutor(4) as pool:
    rows = [r for fam_rows in pool.map(draw, [f for _, f in fams.iterrows()]) for r in fam_rows]
pd.DataFrame(rows).to_csv("data/sequences_pilot.tsv", sep="\t", index=False)
print(len(rows), "sequences")
