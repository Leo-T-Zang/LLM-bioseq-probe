"""Small wrappers over the UniProt, UniParc, UniRef, InterPro and Europe PMC REST APIs."""

from bioseq_mem.http import get

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
    d = get(url, {"format": "json"}).json()
    locs = d["proteins"][0]["entry_protein_locations"]
    return [(f["start"], f["end"]) for loc in locs for f in loc["fragments"]]


def europepmc_hits(query: str) -> int:
    r = get(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
        {"query": query, "format": "json", "pageSize": 1},
    )
    return int(r.json()["hitCount"])
