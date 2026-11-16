"""Every Pfam family with its InterPro counters -> data/pfam_families.tsv. Free."""

import pandas as pd

from bioseq_probe.http import get

url = "https://www.ebi.ac.uk/interpro/api/entry/pfam/"
params = {"page_size": 200, "extra_fields": "counters", "format": "json"}
rows = []
while url:
    page = get(url, params).json()
    params = None  # the "next" link carries them
    for r in page["results"]:
        m, c = r["metadata"], r["extra_fields"]["counters"]
        rows.append(
            {
                "pfam": m["accession"],
                "name": m["name"],
                "type": m.get("type"),
                "proteins": c.get("proteins", 0),
                "structures": c.get("structures", 0),
                "taxa": c.get("taxa", 0),
            }
        )
    url = page.get("next")
    print(len(rows), end="\r", flush=True)

df = pd.DataFrame(rows)
df.to_csv("data/pfam_families.tsv", sep="\t", index=False)
print(f"\n{len(df)} families")
