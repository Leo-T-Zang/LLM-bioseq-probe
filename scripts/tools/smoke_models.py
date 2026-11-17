"""One tiny probe per model to check its settings are accepted. Costs cents."""

import sys
from concurrent.futures import ThreadPoolExecutor

from bioseq_probe.probe import cut_off, parse, prompt
from bioseq_probe.providers import MODELS, call

PREFIX = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF"  # GFP, famous


def one(model):
    try:
        r = call(model, prompt(PREFIX))
        return model, f"stop={r.stop} pred={parse(r.text, PREFIX, cut=cut_off({'stop': r.stop}))!r} raw={r.text[:60]!r} usage_out={r.usage.get('output_tokens') or r.usage.get('completion_tokens') or r.usage.get('candidates_token_count')}"
    except Exception as e:
        return model, f"ERROR {type(e).__name__}: {str(e)[:300]}"


models = sys.argv[1].split(",") if len(sys.argv) > 1 else list(MODELS)
with ThreadPoolExecutor(12) as pool:
    for m, out in pool.map(one, models):
        print(f"{m:30s} {out}")
