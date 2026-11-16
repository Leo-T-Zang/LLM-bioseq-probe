"""Publication figures from results/probe_scored.tsv -> results/fig_*.{pdf,png}. Free.

Fig. 1  residue-1 accuracy on real vs shuffled sequences, per model, all families.
Fig. 2  real minus shuffled, per model, in each family cell (sequences x papers).
Fig. 3  real minus shuffled on old (before 2023) vs recent (from 2025-06) entries.
Fig. 4  real minus shuffled on non-vertebrate DUF families (post hoc: vertebrate DUF members
        have well-known human orthologs), sequence-rich vs sequence-poor.
Intervals: 95% bootstrap over families (cut points within a family are not independent).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bioseq_probe.providers import MODELS

FULL = 2220
LABEL = {
    "gpt-4o-2024-11-20": "GPT-4o", "gpt-4.1-2025-04-14": "GPT-4.1", "gpt-5.5": "GPT-5.5", "gpt-6.1-sol": "GPT-6.1 Sol", "gpt-6-astra": "GPT-6 Astra",
    "claude-sonnet-4-5-20250929": "Claude Sonnet 4.5", "claude-opus-4-6": "Claude Opus 4.6", "claude-opus-4-8": "Claude Opus 4.8", "claude-opus-5": "Claude Opus 5", "claude-opus-5-5": "Claude Opus 5.5",
    "gemini-2.5-pro": "Gemini 2.5 Pro", "gemini-3.5-flash": "Gemini 3.5 Flash", "gemini-3.8-flash": "Gemini 3.8 Flash", "gemini-4-carbon": "Gemini 4 Carbon",
}
VENDOR = {"openai_chat": "OpenAI", "openai_responses": "OpenAI", "anthropic": "Anthropic", "gemini": "Google"}
COLOR = {"OpenAI": "#2a78d6", "Anthropic": "#eb6834", "Google": "#1baf7a"}
CELLS = [  # 2x2 grid: rows = sequences in the family, columns = papers
    ("seq-rich/lit-rich", "Many sequences, many papers"),
    ("seq-rich/lit-poor", "Many sequences, few papers (DUFs)"),
    ("seq-poor/lit-rich", "Few sequences, many papers"),
    ("seq-poor/lit-poor", "Few sequences, few papers (DUFs)"),
]
# Titles state the finding; checked against results/analysis.txt after each run.
TITLE1 = "Newer Claude and Gemini models predict the next residue well above composition; GPT models do not"
TITLE2 = "Every signal is largest in well-studied families; in DUF families only Claude Opus 5 and 5.5 keep a large one"
TITLE3 = "Claude Opus 5.5 scores the same on old and newly added entries; Gemini models tend to drop on new ones"
TITLE4 = "Bacterial, plant and fungal DUFs: only Claude Opus 4.8 and later beat composition clearly"
INK, MUTED, GRID = "#1a1a1a", "#6b6a66", "#e6e5e1"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.7, "xtick.color": MUTED, "ytick.color": INK,
    "xtick.major.width": 0.7, "ytick.major.size": 0, "axes.spines.top": False, "axes.spines.right": False,
    "savefig.dpi": 300, "pdf.fonttype": 42,
})
rng = np.random.default_rng(0)

df = pd.read_csv("results/probe_scored.tsv", sep="\t")
valid = df[df.outcome == "valid"]
models = [m for m in MODELS if m in set(df.model) and m in LABEL]  # closed models; open-weight ones are in common_vs_obscure.py
n_done = df.groupby("model").size()
name = {m: LABEL[m] + ("  (partial)" if n_done[m] < FULL else "") for m in models}
vendor = {m: VENDOR[MODELS[m][0]] for m in models}
ypos = {m: len(models) - 1 - i + (2 - ["OpenAI", "Anthropic", "Google"].index(vendor[m])) * 0.6 for i, m in enumerate(models)}


def boot(g: pd.DataFrame, diff: bool, n: int = 4000):
    """Accuracy (or real - shuffled) and its 95% CI, resampling families."""
    t = g.groupby(["pfam", "condition"]).first_correct.agg(["sum", "count"]).unstack("condition").fillna(0)
    idx = rng.integers(0, len(t), (n, len(t)))

    def stat(cond, ix=slice(None)):
        s, c = t[("sum", cond)].to_numpy(), t[("count", cond)].to_numpy()
        return s[ix].sum(-1) / c[ix].sum(-1)

    if diff:
        point, draws = stat("real") - stat("shuffled"), stat("real", idx) - stat("shuffled", idx)
    else:
        point, draws = stat("real"), stat("real", idx)
    return point, *np.nanpercentile(draws, [2.5, 97.5])


def style_axis(ax, xlabel):
    ax.grid(axis="x", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(xlabel)


def yticks(ax, show=True):
    ax.set_yticks([ypos[m] for m in models], [name[m] for m in models] if show else [])
    for t, m in zip(ax.get_yticklabels(), models):
        t.set_color(COLOR[vendor[m]])


def save(fig, stem):
    Path("results/figures").mkdir(parents=True, exist_ok=True)
    fig.savefig(f"results/figures/{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"results/figures/{stem}.png", bbox_inches="tight", facecolor="white")


# ---- Fig. 1: absolute accuracy, real vs shuffled -------------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 4.4))
for m in models:
    g = valid[valid.model == m]
    sh = boot(g.assign(condition=g.condition.map({"shuffled": "real", "real": "x"})), diff=False)  # shuffled as "real"
    re = boot(g, diff=False)
    y, c = ypos[m], COLOR[vendor[m]]
    ax.plot([sh[0] * 100, re[0] * 100], [y, y], color=c, lw=1, alpha=0.35, zorder=1)
    ax.errorbar(re[0] * 100, y, xerr=[[100 * (re[0] - re[1])], [100 * (re[2] - re[0])]], fmt="o", ms=5.5, color=c, mec="white", mew=0.8, elinewidth=1.1, capsize=0, zorder=3)
    ax.errorbar(sh[0] * 100, y, xerr=[[100 * (sh[0] - sh[1])], [100 * (sh[2] - sh[0])]], fmt="o", ms=5.5, mfc="white", color=c, mew=1.1, elinewidth=1.1, capsize=0, zorder=2)
copy = df[df.condition == "real"].drop_duplicates(["accession", "cut"])["copy"].mean() * 100
ax.axvline(5, color=MUTED, lw=0.8, ls=":")
ax.axvline(copy, color=MUTED, lw=0.8, ls="--")
yticks(ax)
style_axis(ax, "Next-residue accuracy (%)")
ax.set_xlim(0, None)
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=MUTED, mec="white", ms=6, label="real sequence"),
                   plt.Line2D([], [], marker="o", ls="", mfc="white", color=MUTED, ms=6, label="same residues, shuffled"),
                   plt.Line2D([], [], ls=":", color=MUTED, lw=0.8, label="chance (1 in 20)"),
                   plt.Line2D([], [], ls="--", color=MUTED, lw=0.8, label=f"copy-from-prefix baseline ({copy:.0f}%)")],
          loc="upper right", frameon=False, fontsize=7.5)
ax.set_title(TITLE1, loc="left", fontweight="bold")
save(fig, "fig1_accuracy")

# ---- Fig. 2: real - shuffled by family cell ------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(7.2, 7.4), sharey=True, sharex=True)
for k, (ax, (cell, title)) in enumerate(zip(axes.flat, CELLS)):
    nfam = valid[valid.cell == cell].pfam.nunique()
    ax.axvline(0, color=MUTED, lw=0.8)
    for m in models:
        g = valid[(valid.model == m) & (valid.cell == cell)]
        if g.condition.nunique() < 2:
            continue
        d, lo, hi = boot(g, diff=True)
        ax.errorbar(d * 100, ypos[m], xerr=[[100 * (d - lo)], [100 * (hi - d)]], fmt="o", ms=4.5, color=COLOR[vendor[m]], mec="white", mew=0.6, elinewidth=1, capsize=0)
    ax.set_title(f"{'abcd'[k]}   {title}  ({nfam} families)", loc="left", fontsize=8.5)
    style_axis(ax, "Real minus shuffled (percentage points)" if k >= 2 else "")
    yticks(ax, show=k % 2 == 0)
fig.suptitle(TITLE2, x=0.01, ha="left", fontweight="bold", fontsize=9)
fig.tight_layout(rect=(0, 0, 1, 0.97), h_pad=2.0)
save(fig, "fig2_cells")

# ---- Fig. 3: old vs recent entries ---------------------------------------------------------
fig, ax = plt.subplots(figsize=(5.4, 4.2))
ax.axvline(0, color=MUTED, lw=0.8)
for m in models:
    for era, dy, filled in (("old", 0.17, True), ("recent", -0.17, False)):
        g = valid[(valid.model == m) & (valid.era == era)]
        if g.condition.nunique() < 2:
            continue
        d, lo, hi = boot(g, diff=True)
        c = COLOR[vendor[m]]
        ax.errorbar(d * 100, ypos[m] + dy, xerr=[[100 * (d - lo)], [100 * (hi - d)]], fmt="o", ms=4.5, color=c,
                    mfc=c if filled else "white", mec=c if not filled else "white", mew=0.9, elinewidth=1, capsize=0)
yticks(ax)
style_axis(ax, "Real minus shuffled accuracy (percentage points)")
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=MUTED, mec="white", ms=5.5, label="entries on UniProt before 2023"),
                   plt.Line2D([], [], marker="o", ls="", mfc="white", color=MUTED, ms=5.5, label="entries added from June 2025")],
          loc="lower right", frameon=False, fontsize=7.5)
ax.set_title(TITLE3, loc="left", fontweight="bold")
save(fig, "fig3_era")

# ---- Fig. 4: non-vertebrate DUF families, post hoc -----------------------------------------
group = pd.read_csv("data/sequence_groups.tsv", sep="\t").set_index("accession").group
nv = valid[valid.accession.map(group).ne("Vertebrata") & valid.cell.str.endswith("lit-poor")]
fig, ax = plt.subplots(figsize=(5.4, 4.2))
ax.axvline(0, color=MUTED, lw=0.8)
for cell, dy, filled in (("seq-rich/lit-poor", 0.17, True), ("seq-poor/lit-poor", -0.17, False)):
    for m in models:
        g = nv[(nv.model == m) & (nv.cell == cell)]
        d, lo, hi = boot(g, diff=True)
        c = COLOR[vendor[m]]
        ax.errorbar(d * 100, ypos[m] + dy, xerr=[[100 * (d - lo)], [100 * (hi - d)]], fmt="o", ms=4.5, color=c,
                    mfc=c if filled else "white", mec=c if not filled else "white", mew=0.9, elinewidth=1, capsize=0)
nfam = nv.groupby("cell").pfam.nunique()
yticks(ax)
style_axis(ax, "Real minus shuffled accuracy (percentage points)")
ax.legend(handles=[plt.Line2D([], [], marker="o", ls="", color=MUTED, mec="white", ms=5.5, label=f"many sequences ({nfam['seq-rich/lit-poor']} families)"),
                   plt.Line2D([], [], marker="o", ls="", mfc="white", color=MUTED, ms=5.5, label=f"few sequences ({nfam['seq-poor/lit-poor']} families)")],
          loc="lower right", frameon=False, fontsize=7.5)
ax.set_title(TITLE4, loc="left", fontweight="bold")
save(fig, "fig4_nonvertebrate_duf")
print("results/figures/fig1_accuracy, fig2_cells, fig3_era, fig4_nonvertebrate_duf (.pdf, .png)")
