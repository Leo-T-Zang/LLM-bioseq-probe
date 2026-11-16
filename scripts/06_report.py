"""Self-contained HTML report from results/probe_scored.tsv -> results/report.html. Free.

Main figure: per model, residue-1 accuracy on the real sequence minus on the
shuffled one, in each family cell, with a 95% bootstrap CI over families
(cut points within a family are not independent).
"""

import base64
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bioseq_mem.providers import MODELS

VENDOR_COLOR = {"OpenAI": "#2a78d6", "Anthropic": "#eb6834", "Google": "#1baf7a"}  # categorical slots 1-3
CELLS = ["seq-rich/lit-poor", "seq-rich/lit-rich", "seq-poor/lit-poor", "seq-poor/lit-rich"]
CELL_TITLE = {
    "seq-rich/lit-poor": "Many sequences, few papers (DUFs)",
    "seq-rich/lit-rich": "Many sequences, many papers",
    "seq-poor/lit-rich": "Few sequences, many papers",
    "seq-poor/lit-poor": "Few sequences, few papers (DUFs)",
}
rng = np.random.default_rng(0)


def vendor(model: str) -> str:
    return {"openai_chat": "OpenAI", "openai_responses": "OpenAI", "anthropic": "Anthropic", "gemini": "Google"}[MODELS[model][0]]


def family_bootstrap(g: pd.DataFrame, n: int = 2000) -> tuple[float, float, float]:
    """Real-minus-shuffled accuracy, resampling families."""
    per_fam = g.groupby(["pfam", "condition"]).first_correct.agg(["sum", "count"]).unstack("condition").fillna(0)
    s, c = per_fam["sum"].to_numpy(), per_fam["count"].to_numpy()  # columns: real, shuffled
    point = s[:, 0].sum() / c[:, 0].sum() - s[:, 1].sum() / c[:, 1].sum()
    idx = rng.integers(0, len(per_fam), (n, len(per_fam)))
    boot = s[idx, 0].sum(1) / c[idx, 0].sum(1) - s[idx, 1].sum(1) / c[idx, 1].sum(1)
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return point, lo, hi


df = pd.read_csv("results/probe_scored.tsv", sep="\t")
models = [m for m in MODELS if m in set(df.model)]
valid = df[df.outcome == "valid"]

rows = []
for m in models:
    for cell in CELLS:
        g = valid[(valid.model == m) & (valid.cell == cell)]
        if g.condition.nunique() < 2:
            continue
        d, lo, hi = family_bootstrap(g)
        acc = g.groupby("condition").first_correct.mean()
        rows.append({"model": m, "vendor": vendor(m), "cell": cell, "real": acc["real"], "shuffled": acc["shuffled"], "diff": d, "lo": lo, "hi": hi, "n": len(g)})
res = pd.DataFrame(rows)
copy_base = df[df.condition == "real"].drop_duplicates(["accession", "cut"]).groupby("cell")["copy"].mean()

fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
ypos = {m: i for i, m in enumerate(reversed(models))}
for ax, cell in zip(axes.flat, CELLS):
    r = res[res.cell == cell]
    ax.axvline(0, color="#52514e", lw=1)
    for _, x in r.iterrows():
        y = ypos[x.model]
        ax.plot([x.lo * 100, x.hi * 100], [y, y], color=VENDOR_COLOR[x.vendor], lw=2, solid_capstyle="round")
        ax.plot(x["diff"] * 100, y, "o", ms=8, color=VENDOR_COLOR[x.vendor], mec="#fcfcfb", mew=2)
    ax.set_title(CELL_TITLE[cell], fontsize=11, loc="left")
    ax.grid(axis="x", color="#e5e4e0", lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
ax0 = axes[0, 0]
ax0.set_yticks(list(ypos.values()), list(ypos.keys()), fontsize=9)
axes[1, 0].set_yticks(list(ypos.values()), list(ypos.keys()), fontsize=9)
for ax in axes[1]:
    ax.set_xlabel("Residue-1 accuracy, real minus shuffled (percentage points)")
handles = [plt.Line2D([], [], marker="o", ls="", color=c, ms=8, label=v) for v, c in VENDOR_COLOR.items()]
fig.legend(handles=handles, loc="upper right", ncol=3, frameon=False)
fig.suptitle("How much better each model predicts the next residue when residue order is intact", x=0.02, ha="left", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.95))


def png(fig, name: str) -> str:
    fig.savefig(f"results/{name}.png", dpi=160, facecolor="#fcfcfb")
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, facecolor="#fcfcfb")
    return base64.b64encode(buf.getvalue()).decode()


img = png(fig, "fig_main")

# Era split: old entries (before 2023, seen by every model) vs recent (from 2025-06).
era_rows = []
for m in models:
    for era in ("old", "recent"):
        g = valid[(valid.model == m) & (valid.era == era)]
        if g.condition.nunique() == 2:
            d, lo, hi = family_bootstrap(g)
            era_rows.append({"model": m, "vendor": vendor(m), "era": era, "diff": d, "lo": lo, "hi": hi})
era_res = pd.DataFrame(era_rows)
fig2, ax = plt.subplots(figsize=(8, 5.5))
ax.axvline(0, color="#52514e", lw=1)
for _, x in era_res.iterrows():
    y = ypos[x.model] + (0.15 if x.era == "old" else -0.15)
    style = dict(color=VENDOR_COLOR[x.vendor], alpha=1.0 if x.era == "old" else 0.45)
    ax.plot([x.lo * 100, x.hi * 100], [y, y], lw=2, solid_capstyle="round", **style)
    ax.plot(x["diff"] * 100, y, "o" if x.era == "old" else "s", ms=7, mec="#fcfcfb", mew=1.5, **style)
ax.set_yticks(list(ypos.values()), list(ypos.keys()), fontsize=9)
ax.set_xlabel("Residue-1 accuracy, real minus shuffled (percentage points)")
ax.grid(axis="x", color="#e5e4e0", lw=0.8)
ax.set_axisbelow(True)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color="#52514e", label="old entries (before 2023)"),
                   plt.Line2D([], [], marker="s", ls="", color="#52514e", alpha=0.45, label="recent entries (from 2025-06)")],
          loc="lower right", frameon=False, fontsize=9)
ax.set_title("Old vs recent sequences, all families", loc="left", fontsize=11)
fig2.tight_layout()
img2 = png(fig2, "fig_era")

outcomes = pd.crosstab(df.model, df.outcome).reindex(models)
table = res.assign(**{k: (res[k] * 100).round(1) for k in ("real", "shuffled", "diff", "lo", "hi")})
table = table.pivot(index="model", columns="cell", values=["real", "shuffled", "diff"]).reindex(models)
tokens = {}
for p in Path("results/probe").glob("*.jsonl"):
    t_in = t_out = 0
    for d in map(json.loads, p.open()):
        u = d.get("usage") or {}
        t_in += u.get("input_tokens") or u.get("prompt_tokens") or u.get("prompt_token_count") or 0
        t_out += u.get("output_tokens") or u.get("completion_tokens") or ((u.get("candidates_token_count") or 0) + (u.get("thoughts_token_count") or 0))
    tokens[p.stem] = {"input tokens": t_in, "output tokens": t_out}

html = f"""<!doctype html><meta charset="utf-8"><title>bioseq-memorization: next-residue probe</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:1100px;margin:2em auto;padding:0 1em;color:#0b0b0b;background:#fcfcfb}}
table{{border-collapse:collapse;font-size:13px;margin:1em 0}}td,th{{border-bottom:1px solid #e5e4e0;padding:3px 8px;text-align:right}}
th:first-child,td:first-child{{text-align:left}}.muted{{color:#52514e}}img{{max-width:100%}}</style>
<h1>Do closed LLMs know protein sequence order?</h1>
<p class="muted">Each probe gives a model the start of an unreviewed (TrEMBL) protein from a Pfam family
and asks for the next residues. The same cut is asked on the sequence shuffled (same composition, no order).
Real minus shuffled is what the model knows about order beyond composition. Intervals: 95% bootstrap over families.</p>
<img src="data:image/png;base64,{img}" alt="real-minus-shuffled accuracy per model and family cell">
<h2>Old vs recent sequences</h2>
<p class="muted">Higher accuracy on old entries than on recent ones, in a model whose cutoff precedes the recent ones, points to recall.</p>
<img src="data:image/png;base64,{img2}" alt="real-minus-shuffled accuracy by era">
<h2>In-context copy baseline (real sequences)</h2>
<p class="muted">Accuracy from copying what followed an earlier repeat of the prefix's last residues; no training needed.</p>
{(copy_base * 100).round(1).to_frame("copy baseline, %").to_html()}
<h2>Residue-1 accuracy, % (valid replies)</h2>
{table.to_html()}
<h2>Reply outcomes</h2>
<p class="muted">Refusals and errors are not counted as misses.</p>
{outcomes.to_html()}
<h2>Tokens used</h2>
{pd.DataFrame(tokens).T.reindex(models).to_html()}
"""
Path("results/report.html").write_text(html)
print("results/report.html")
