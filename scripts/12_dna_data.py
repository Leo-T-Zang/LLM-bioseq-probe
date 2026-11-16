"""Nucleotide probe windows -> data/dna_windows_{genomic,bacterial}.tsv. Free.

    uv run python scripts/12_dna_data.py genomic|bacterial

Each window is 220 nt: a 200-nt prefix and the 20 nt to predict. Regions:
  human and mouse (Ensembl): CDS of a protein-coding gene (canonical transcript),
    intron (inside a protein-coding gene, no exon, no repeat), intergenic (no gene
    within 5 kb, no repeat); repeats from Ensembl's soft mask and repeat features.
  bacterial CDS: the coding sequence (ENA) of TrEMBL proteins from the confirmation
    set: obscure, UniRef90-singleton bacterial and archaeal proteins.
"""

import random
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_mem.http import get

ENSEMBL = "https://rest.ensembl.org"
WIN = 220
PER_REGION = 60
BACTERIAL = 100
rng = random.Random(7)


def ens(path: str, **params):
    return get(f"{ENSEMBL}{path}", {"content-type": "application/json", **params}).json()


def region_seq(species, chrom, start):
    r = get(f"{ENSEMBL}/sequence/region/{species}/{chrom}:{start}..{start + WIN - 1}:1", {"mask": "soft", "content-type": "text/plain"})
    return r.text.strip()


def sample_genomic(species: str, n_chrom: int) -> list[dict]:
    chroms = {c["name"]: c["length"] for c in ens(f"/info/assembly/{species}")["top_level_region"]
              if c["coord_system"] == "chromosome" and c["name"] in {str(i) for i in range(1, n_chrom + 1)}}
    names, weights = list(chroms), list(chroms.values())
    got = {"cds": [], "intron": [], "intergenic": []}
    while min(len(v) for v in got.values()) < PER_REGION:
        chrom = rng.choices(names, weights)[0]
        pos = rng.randrange(10_000, chroms[chrom] - 10_000)
        near = ens(f"/overlap/region/{species}/{chrom}:{pos - 5000}-{pos + WIN + 5000}", feature="gene")
        window = ens(f"/overlap/region/{species}/{chrom}:{pos}-{pos + WIN - 1}", feature=["exon", "repeat", "gene"])
        coding = [g for g in window if g["feature_type"] == "gene" and g.get("biotype") == "protein_coding"]
        has = lambda t: any(f["feature_type"] == t for f in window)
        if not near and not has("repeat") and len(got["intergenic"]) < PER_REGION:
            region, locus = "intergenic", f"{chrom}:{pos}"
        elif coding and not has("exon") and not has("repeat") and len(got["intron"]) < PER_REGION:
            region, locus = "intron", f"{chrom}:{pos} in {coding[0]['id']}"
        elif coding and len(got["cds"]) < PER_REGION:
            gene = ens(f"/lookup/id/{coding[0]['id']}", expand=1)
            tx = [t for t in gene.get("Transcript", []) if t.get("is_canonical")]
            if not tx:
                continue
            cds = get(f"{ENSEMBL}/sequence/id/{tx[0]['id']}", {"type": "cds", "content-type": "text/plain"}).text.strip()
            if len(cds) < WIN + 30 or set(cds) - set("ACGT"):
                continue
            start = rng.randrange(0, len(cds) - WIN)
            got["cds"].append({"species": species, "region": "cds", "locus": f"{tx[0]['id']}:{start}", "window": cds[start : start + WIN]})
            continue
        else:
            continue
        s = region_seq(species, chrom, pos)
        if len(s) == WIN and set(s) <= set("ACGT"):  # soft-masked repeats are lowercase
            got[region].append({"species": species, "region": region, "locus": locus, "window": s})
    print(species, {k: len(v) for k, v in got.items()}, flush=True)
    return [r for v in got.values() for r in v]


def bacterial_cds(acc: str) -> dict | None:
    entry = get(f"https://rest.uniprot.org/uniprotkb/{acc}", {"fields": "xref_embl"}).json()
    pids = [p["value"] for x in entry.get("uniProtKBCrossReferences", []) for p in x.get("properties", [])
            if x["database"] == "EMBL" and p["key"] == "ProteinId" and p["value"] != "-"]
    if not pids:
        return None
    fasta = get(f"https://www.ebi.ac.uk/ena/browser/api/fasta/{pids[0]}").text
    seq = "".join(fasta.splitlines()[1:]).upper()
    if len(seq) < WIN + 30 or set(seq) - set("ACGT"):
        return None
    start = random.Random(acc).randrange(0, len(seq) - WIN)
    return {"species": "bacteria/archaea", "region": "cds", "locus": f"{acc} {pids[0]}:{start}", "window": seq[start : start + WIN]}


part = sys.argv[1]  # "genomic" or "bacterial"
if part == "genomic":
    rows = sample_genomic("human", 22) + sample_genomic("mouse", 19)
else:
    accs = pd.read_csv("data/confirm_sequences.tsv", sep="\t").accession.drop_duplicates().sample(frac=1, random_state=7).tolist()
    with ThreadPoolExecutor(8) as pool:
        rows = [r for r in pool.map(bacterial_cds, accs[: BACTERIAL * 2]) if r][:BACTERIAL]
df = pd.DataFrame(rows)
df.to_csv(f"data/dna_windows_{part}.tsv", sep="\t", index=False)
print(df.groupby(["species", "region"]).size())
