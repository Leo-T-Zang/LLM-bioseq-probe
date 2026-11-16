#!/usr/bin/env bash
# Launch the confirmation and DNA runs as their inputs appear; stop all probe runs past the spend cap.
set -u
cd "$(dirname "$0")/.."
CAP=290

spent() { uv run python -W ignore scripts/cost.py | tail -1 | sed -E 's/spent \$([0-9.]+).*/\1/'; }

until [ -s data/confirm_sequences.tsv ] && ! pgrep -f "scripts/[0]8_confirm_data" >/dev/null; do sleep 30; done
for m in claude-opus-5-5 claude-opus-5 claude-opus-4-6 gpt-5.5 gemini-4-carbon gpt-6-astra; do
  nohup uv run python -W ignore scripts/04_probe.py --sequences data/confirm_sequences.tsv --out results/confirm \
    --models "$m" --workers 24 --run > "results/logs/confirm_$m.log" 2>&1 &
done
echo "confirmation runs started $(date -u)"

uv run python -W ignore -u scripts/12_dna_data.py bacterial > data/12_bacterial.log 2>&1
until [ -s data/dna_windows_genomic.tsv ] && ! pgrep -f "scripts/[1]2_dna_data" >/dev/null; do sleep 30; done
for m in claude-opus-5-5 claude-opus-5 gemini-4-carbon gpt-5.5 deepseek-v4p1-flash; do
  nohup uv run python -W ignore scripts/13_dna_probe.py --models "$m" --workers 16 --run > "results/logs/dna_$m.log" 2>&1 &
done
echo "DNA runs started $(date -u)"

sleep 60
while pgrep -f "scripts/[01][0-9]_[a-z_]*probe.py" >/dev/null; do
  s=$(spent)
  if python3 -c "import sys; sys.exit(0 if float('$s') > $CAP else 1)"; then
    pkill -f "scripts/[01][0-9]_[a-z_]*probe.py"
    echo "spend cap hit at \$$s"
    exit 1
  fi
  sleep 120
done
echo "all runs finished $(date -u); spent \$$(spent)"
