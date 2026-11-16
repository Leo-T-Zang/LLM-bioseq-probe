#!/usr/bin/env bash
# UniParc id and Pfam annotations of every unreviewed UniProtKB entry created on or after the
# cutoff, length 120-800, in 8 parallel length bins -> data/new_entries/*.tsv. Free.
# Used to show that new UniProtKB entries almost never carry a new sequence:
# of ~3.8M entries created 2026-07-01..2026-11-04, ~10 have a sequence first seen after 2026-07-01.
set -euo pipefail
CUTOFF=2026-07-01
bins=("120 TO 180" "181 TO 230" "231 TO 280" "281 TO 340" "341 TO 410" "411 TO 490" "491 TO 600" "601 TO 800")
for i in "${!bins[@]}"; do
  q="date_created:[$CUTOFF TO *] AND reviewed:false AND length:[${bins[$i]}]"
  curl -sS --compressed --retry 5 -G "https://rest.uniprot.org/uniprotkb/stream" \
    --data-urlencode "query=$q" --data-urlencode "fields=accession,uniparc_id,xref_pfam" --data-urlencode "format=tsv" \
    -o "data/new_entries/bin$i.tsv" &
done
wait
