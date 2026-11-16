"""Nucleotide probe windows -> data/dna_windows_{genomic,bacterial}.tsv. Free.

    uv run python scripts/12_dna_data.py genomic|bacterial|confirm

Each window is 220 nt (5'->3' on the feature's strand): a 200-nt prefix and the 20 nt
to predict.

genomic: offline, from UCSC hg38/mm39 .2bit (RepeatMasker soft-masking) and GENCODE
basic annotation (human v50, mouse vM39; Ensembl canonical transcripts), in
data/genomes/. Per species, PER_REGION windows of each region:
  cds          spliced CDS of protein-coding transcripts
  utr5, utr3   spliced UTRs of protein-coding transcripts
  promoter     1 kb upstream of protein-coding TSSs
  intron       introns of protein-coding transcripts
  lncrna       spliced exons of lncRNA transcripts
  intergenic   no annotated gene within 5 kb
Non-coding windows containing any repeat (lowercase) or N are rejected; CDS windows
only reject N.

bacterial: the coding sequence (ENA) of TrEMBL proteins from the confirmation set:
obscure, UniRef90-singleton bacterial and archaeal proteins.
"""

import bisect
import gzip
import random
import re
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

from bioseq_mem.http import get
from bioseq_mem.twobit import TwoBit

WIN = 220
PER_REGION = 100
BACTERIAL = 100
GENOMES = {
    "human": ("data/genomes/hg38.2bit", "data/genomes/gencode.v50.basic.annotation.gtf.gz", 22),
    "mouse": ("data/genomes/mm39.2bit", "data/genomes/gencode.vM39.basic.annotation.gtf.gz", 19),
}
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_gtf(path: str, chroms: set[str]):
    """Canonical transcripts: id -> (chrom, strand, gene_type, {feature: [(start0, end)]}); and all gene spans."""
    tx, genes = {}, defaultdict(list)
    attr = re.compile(r'(\S+) "([^"]*)"')
    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[0] not in chroms:
                continue
            start, end = int(c[3]) - 1, int(c[4])
            if c[2] == "gene":
                genes[c[0]].append((start, end))
                continue
            if 'tag "Ensembl_canonical"' not in c[8] or c[2] not in ("exon", "CDS", "UTR", "transcript"):
                continue
            a = dict(attr.findall(c[8]))
            t = tx.setdefault(a["transcript_id"], (c[0], c[6], a.get("gene_type"), defaultdict(list)))
            if c[2] == "transcript":
                t[3]["transcript"].append((start, end))
            else:
                t[3][c[2]].append((start, end))
    return tx, genes


def spliced(tb: TwoBit, chrom: str, strand: str, parts: list[tuple[int, int]]) -> str:
    s = "".join(tb.get(chrom, a, b) for a, b in sorted(parts))
    return revcomp(s) if strand == "-" else s


def window(s: str, rng: random.Random, coding: bool) -> str | None:
    if len(s) < WIN:
        return None
    for _ in range(20):
        i = rng.randrange(0, len(s) - WIN + 1)
        w = s[i : i + WIN]
        if "N" not in w.upper() and (coding or w.isupper()):
            return w.upper()
    return None


ALL_REGIONS = ("cds", "utr5", "utr3", "promoter", "intron", "lncrna", "intergenic")


def sample_species(species: str, regions=ALL_REGIONS, per=PER_REGION, seed="dna", exclude=frozenset()) -> list[dict]:
    """exclude: transcript ids (and intergenic loci) already used, e.g. by the pilot."""
    path2bit, gtf, n = GENOMES[species]
    tb = TwoBit(path2bit)
    chroms = {f"chr{i}" for i in range(1, n + 1)}
    tx, genes = read_gtf(gtf, chroms)
    rng = random.Random(f"{seed}-{species}")
    pc = [(t, v) for t, v in tx.items() if v[2] == "protein_coding" and v[3]["CDS"] and t not in exclude]
    lnc = [(t, v) for t, v in tx.items() if v[2] == "lncRNA" and t not in exclude]
    rows = []

    def take(region, items, make):
        if region not in regions:
            return
        got, tries = 0, 0
        order = items[:]
        rng.shuffle(order)
        for t, v in order:
            tries += 1
            w = make(t, v)
            if w:
                rows.append({"species": species, "region": region, "locus": f"{species}:{region}:{t}", "window": w})
                got += 1
                if got == per:
                    break
        print(species, region, got, f"from {tries} transcripts", flush=True)

    def utr(v, five):
        chrom, strand, _, f = v
        cds_lo, cds_hi = min(a for a, _ in f["CDS"]), max(b for _, b in f["CDS"])
        left = [(a, b) for a, b in f["UTR"] if b <= cds_lo]
        right = [(a, b) for a, b in f["UTR"] if a >= cds_hi]
        parts = (left if strand == "+" else right) if five else (right if strand == "+" else left)
        return spliced(tb, chrom, strand, parts) if parts else ""

    def promoter(v):
        chrom, strand, _, f = v
        a, b = f["transcript"][0]
        lo, hi = (a - 1000, a) if strand == "+" else (b, b + 1000)
        s = tb.get(chrom, max(lo, 0), hi)
        return revcomp(s) if strand == "-" else s

    def intron(v):
        chrom, strand, _, f = v
        ex = sorted(f["exon"])
        gaps = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1) if ex[i + 1][0] - ex[i][1] >= WIN + 100]
        if not gaps:
            return ""
        a, b = rng.choice(gaps)
        i = rng.randrange(a + 50, b - 50 - WIN + 1)
        s = tb.get(chrom, i, i + WIN)
        return revcomp(s) if strand == "-" else s

    take("cds", pc, lambda t, v: window(spliced(tb, v[0], v[1], v[3]["CDS"]), rng, coding=True))
    take("utr5", pc, lambda t, v: window(utr(v, True), rng, coding=False))
    take("utr3", pc, lambda t, v: window(utr(v, False), rng, coding=False))
    take("promoter", pc, lambda t, v: window(promoter(v), rng, coding=False))
    take("intron", pc, lambda t, v: window(intron(v), rng, coding=False))
    take("lncrna", lnc, lambda t, v: window(spliced(tb, v[0], v[1], v[3]["exon"]), rng, coding=False))

    if "intergenic" not in regions:
        return rows
    # intergenic: no gene within 5 kb
    padded = {}
    for c, spans in genes.items():
        merged = []
        for a, b in sorted((max(a - 5000, 0), b + 5000) for a, b in spans):
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        padded[c] = ([a for a, _ in merged], [b for _, b in merged])
    names = sorted(chroms)
    weights = [tb.size(c) for c in names]
    got = tries = 0
    while got < per:
        tries += 1
        c = rng.choices(names, weights)[0]
        i = rng.randrange(0, tb.size(c) - WIN)
        starts, ends = padded.get(c, ([], []))
        k = bisect.bisect_right(starts, i + WIN) - 1
        if k >= 0 and ends[k] > i:
            continue
        w = tb.get(c, i, i + WIN)
        if w.isupper() and "N" not in w and f"{species}:intergenic:{c}:{i}" not in exclude:
            rows.append({"species": species, "region": "intergenic", "locus": f"{species}:intergenic:{c}:{i}", "window": w})
            got += 1
    print(species, "intergenic", got, f"from {tries} positions", flush=True)
    return rows


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


part = sys.argv[1]
if part == "genomic":
    rows = sample_species("human") + sample_species("mouse")
    out = "data/dna_windows_genomic.tsv"
elif part == "bacterial":
    accs = pd.read_csv("data/confirm_sequences.tsv", sep="\t").accession.drop_duplicates().sample(frac=1, random_state=7).tolist()
    with ThreadPoolExecutor(8) as pool:
        rows = [r for r in pool.map(bacterial_cds, accs[: BACTERIAL * 2]) if r][:BACTERIAL]
    out = "data/dna_windows_bacterial.tsv"
else:  # "confirm": pre-registered DNA confirmation run, docs/preregistration.md
    pilot = pd.concat([pd.read_csv(p, sep="\t") for p in ("data/dna_windows_genomic.tsv", "data/dna_windows_bacterial.tsv")])
    used = {loc.split(":")[2] for loc in pilot.locus if loc.count(":") == 2} | set(pilot.locus)
    rows = []
    for sp in ("human", "mouse"):
        rows += sample_species(sp, ("cds", "intron", "promoter", "intergenic"), per=500, seed="dna-confirm", exclude=used)
    pilot_accs = {loc.split()[0] for loc in pilot.locus if " " in loc}
    accs = [a for a in pd.read_csv("data/confirm_sequences.tsv", sep="\t").accession.drop_duplicates() if a not in pilot_accs]
    with ThreadPoolExecutor(8) as pool:
        rows += [r for r in pool.map(bacterial_cds, accs) if r][:500]
    out = "data/dna_confirm_windows.tsv"
df = pd.DataFrame(rows)
df.to_csv(out, sep="\t", index=False)
print(df.groupby(["species", "region"]).size())
