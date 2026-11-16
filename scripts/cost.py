"""Spend so far, and projected for the full run, from the recorded token usage. Free.

Prices are USD per 1M tokens (input, output). Anthropic's are from its pricing
table, OpenAI's from developers.openai.com/api/docs/pricing (2026-11-16); those
marked ASSUMED are guesses set high, so the projection errs on overspending.
"""

import json
from pathlib import Path

PRICE = {
    "claude-opus-5-5": (4, 20), "claude-opus-5": (5, 25), "claude-opus-4-6": (5, 25), "claude-sonnet-4-5-20250929": (3, 15),
    "gpt-4o-2024-11-20": (2.5, 10), "gpt-4.1-2025-04-14": (2, 8),
    "gpt-5.5": (5, 30), "gpt-6.1-sol": (2, 10), "gpt-6-astra": (10, 50),
    "gemini-2.5-pro": (1.25, 10), "gemini-3.5-flash": (1, 6), "gemini-3.8-flash": (1, 6), "gemini-4-carbon": (5, 30),  # flash/4 ASSUMED
}
FIREWORKS = (3, 12)  # ASSUMED for every open-weight model on Fireworks, set high
FULL = 2220  # probes per model in the pilot run


def tokens(u: dict) -> tuple[int, int]:
    t_in = u.get("input_tokens") or u.get("prompt_tokens") or u.get("prompt_token_count") or 0
    t_out = u.get("output_tokens") or u.get("completion_tokens") or ((u.get("candidates_token_count") or 0) + (u.get("thoughts_token_count") or 0))
    return t_in, t_out


total_now = total_full = 0.0
for p in sorted(p for d in ("results/probe", "results/confirm", "results/dna") for p in Path(d).glob("*.jsonl")):
    recs = [json.loads(line) for line in p.open()]
    t = [tokens(r.get("usage") or {}) for r in recs]
    pi, po = PRICE.get(p.stem, FIREWORKS)
    spent = sum(a * pi + b * po for a, b in t) / 1e6
    per = spent / len(recs)
    total_now += spent
    total_full += per * FULL
    print(f"{p.parent.name + '/' + p.stem:40s} n={len(recs):5d} spent ${spent:7.3f}  per probe ${per:.5f}  full ${per * FULL:7.2f}  mean out tokens {sum(b for _, b in t) / len(t):.0f}")
print(f"spent ${total_now:.2f}; projected full run ${total_full:.2f}")
