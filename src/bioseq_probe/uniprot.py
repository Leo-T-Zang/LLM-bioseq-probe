"""Small wrappers over the UniProt, UniParc, UniRef, InterPro and Europe PMC REST APIs."""

from bioseq_probe.http import get

UNIPROT = "https://rest.uniprot.org"


def count(query: str, db: str = "uniprotkb") -> int:
    r = get(f"{UNIPROT}/{db}/search", {"query": query, "size": 0})
    return int(r.headers["x-total-results"])


def search(query: str, fields: str, size: int = 500, db: str = "uniprotkb") -> list[dict]:
    r = get(f"{UNIPROT}/{db}/search", {"query": query, "fields": fields, "size": size})
    return r.json()["results"]


def uniparc_first_seen(upis: list[str]) -> dict[str, str]:
    """UniParc id -> date the sequence first appeared in any source database."""
    out = {}
    for i in range(0, len(upis), 100):
        chunk = upis[i : i + 100]
        q = " OR ".join(f"upi:{u}" for u in chunk)
        for r in search(q, "upi,oldestCrossRefCreated", size=100, db="uniparc"):
            out[r["uniParcId"]] = r["oldestCrossRefCreated"]
    return out


def uniref90_size(accession: str) -> int | None:
    rs = search(f"uniprot_id:{accession} AND identity:0.9", "id,count", size=1, db="uniref")
    return rs[0]["memberCount"] if rs else None


def pfam_spans(accession: str, pfam: str) -> list[tuple[int, int]]:
    url = f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{pfam}/protein/uniprot/{accession}/"
    r = get(url, {"format": "json"})
    if r.status_code == 204 or not r.content:  # no match of this family on the protein
        return []
    d = r.json()
    locs = d["proteins"][0]["entry_protein_locations"]
    return [(f["start"], f["end"]) for loc in locs for f in loc["fragments"]]


def europepmc_hits(query: str) -> int:
    r = get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        {"query": query, "format": "json", "pageSize": 1},
    )
    return int(r.json()["hitCount"])


def uniref90_singletons(accessions: list[str]) -> set[str]:
    """The accessions that are alone in their UniRef90 cluster (no other sequence at >= 90% identity)."""
    out = set()
    for i in range(0, len(accessions), 50):
        q = " OR ".join(f"uniprot_id:{a}" for a in accessions[i : i + 50])
        for r in search(f"({q}) AND identity:0.9", "id,count", size=100, db="uniref"):
            if r["memberCount"] == 1:
                out.add(r["id"].removeprefix("UniRef90_"))
    return out
