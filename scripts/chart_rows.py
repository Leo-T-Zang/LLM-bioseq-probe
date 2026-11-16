"""Rows for the report's native charts (pilot): JS object literals on stdout. Free."""

import numpy as np
import pandas as pd

from bioseq_mem.providers import MODELS

LABEL = {"gpt-4o-2024-11-20": "GPT-4o", "gpt-4.1-2025-04-14": "GPT-4.1", "gpt-5.5": "GPT-5.5", "gpt-6.1-sol": "GPT-6.1 Sol",
         "claude-sonnet-4-5-20250929": "Claude Sonnet 4.5", "claude-opus-4-6": "Claude Opus 4.6", "claude-opus-5": "Claude Opus 5",
         "claude-opus-5-5": "Claude Opus 5.5", "gemini-2.5-pro": "Gemini 2.5 Pro", "gemini-3.5-flash": "Gemini 3.5 Flash",
         "gemini-3.8-flash": "Gemini 3.8 Flash", "gemini-4-carbon": "Gemini 4 Carbon"}
VENDOR = {"openai_chat": "OpenAI", "openai_responses": "OpenAI", "anthropic": "Anthropic", "gemini": "Google"}
rng = np.random.default_rng(0)
v = pd.read_csv("results/probe_scored.tsv", sep="\t")
v = v[v.outcome == "valid"]


def boot(g, cond=None):
    t = g.groupby(["pfam", "condition"]).first_correct.agg(["sum", "count"]).unstack("condition").fillna(0)
    idx = rng.integers(0, len(t), (4000, len(t)))
    st = lambda c, ix=slice(None): t[("sum", c)].to_numpy()[ix].sum(-1) / t[("count", c)].to_numpy()[ix].sum(-1)
    if cond:
        pt, b = st(cond), st(cond, idx)
    else:
        pt, b = st("real") - st("shuffled"), st("real", idx) - st("shuffled", idx)
    lo, hi = np.percentile(b, [2.5, 97.5])
    return round(pt * 100, 1), round(lo * 100, 1), round(hi * 100, 1)


def js(rows):
    return "[" + ",".join("{" + ", ".join(f"{k}:{v!r}" if isinstance(v, str) else f"{k}:{v}" for k, v in r.items()) + "}" for r in rows) + "]"


a, b, c = [], [], []
for m in MODELS:
    g = v[v.model == m]
    base = {"model": LABEL[m], "vendor": VENDOR[MODELS[m][0]]}
    for cond in ("real", "shuffled"):
        pt, lo, hi = boot(g, cond)
        a.append({**base, "seq": cond, "acc": pt, "lo": lo, "hi": hi})
    for cell in ("seq-rich/lit-rich", "seq-rich/lit-poor", "seq-poor/lit-rich", "seq-poor/lit-poor"):
        pt, lo, hi = boot(g[g.cell == cell])
        b.append({**base, "cell": cell, "diff": pt, "lo": lo, "hi": hi})
    for era in ("old", "recent"):
        pt, lo, hi = boot(g[g.era == era])
        c.append({**base, "era": era, "diff": pt, "lo": lo, "hi": hi})
for name, rows in (("A", a), ("B", b), ("C", c)):
    print(name, js(rows))
