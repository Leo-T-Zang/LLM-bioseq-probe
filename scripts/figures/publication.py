"""The five figures used in the README and the write-up -> figures/fig1..fig5 (.png, .pdf). Free.

Reads the summary tables written by scripts/analysis and scripts/figures:
results/confirm_tests.tsv, results/raw_summary.tsv, results/dna_confirm_tests.tsv,
results/common_vs_obscure.tsv, results/recall_length_summary.tsv, and for Fig 5 the replies
in results/confirm/ with results/confirm_family_opus55.tsv (Europe PMC hits per family).
"""

import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from bioseq_probe.figstyle import ACCENT, GREY, GRID, INK, LIGHT_GREY, MUTED, color, header, save

LABEL = {"claude-opus-4-6": "Claude Opus 4.6", "claude-opus-4-8": "Claude Opus 4.8", "claude-opus-5": "Claude Opus 5",
         "claude-opus-5-5": "Claude Opus 5.5", "gpt-5.5": "GPT-5.5", "gpt-6-astra": "GPT-6 Astra", "gemini-4-carbon": "Gemini 4 Carbon",
         "claude-sonnet-4-5-20250929": "Claude Sonnet 4.5", "gemini-3.5-flash": "Gemini 3.5 Flash", "gemini-3.8-flash": "Gemini 3.8 Flash",
         "gpt-6.1-sol": "GPT-6.1 Sol"}


def forest(ax, rows, y_of, xlabel):
    """rows: dicts with name, d, lo, hi, raw (bool)."""
    for r in rows:
        y, c = y_of[r["name"]], color(r["name"])
        ax.plot([r["lo"], r["hi"]], [y, y], color=c, lw=1.4, solid_capstyle="round", zorder=2)
        ax.plot(r["d"], y, "o", ms=4.2, color=c, mfc="white" if r["raw"] else c, mec=c, mew=1.0, zorder=3)
        ax.text(max(r["hi"], 0) + 0.8, y, f"{r['d']:+.1f}", va="center", fontsize=6.3,
                color=INK if r["name"] == "Claude Opus 5.5" else MUTED, fontweight="bold" if r["name"] == "Claude Opus 5.5" else "normal")
    ax.axvline(0, color=MUTED, lw=0.6, zorder=1)
    ax.set_xlabel(xlabel)


# Fig 1: every model on the obscure proteins ---------------------------------------------------
c = pd.read_csv("results/confirm_tests.tsv", sep="\t").set_index("model")
raw = pd.read_csv("results/raw_summary.tsv", sep="\t")
groups = [
    ("Anthropic", [("claude-opus-5-5", "chat"), ("claude-opus-5", "chat"), ("claude-opus-4-8", "chat"), ("claude-opus-4-6", "chat")]),
    ("OpenAI", [("gpt-6-astra", "chat"), ("gpt-5.5", "chat")]),
    ("Google", [("gemini-4-carbon", "chat"), ("Gemini 2.5 Pro (raw)", "raw")]),
    ("Open-weight", [(n, "raw") for n in raw.name if n.endswith("(raw)") and not n.startswith(("Gemini", "Claude"))]),
]
rows, y_of, ticks, y = [], {}, [], 0.0
for g, members in groups:
    ticks.append((g, y))
    y -= 0.85
    for key, how in members:
        if how == "chat":
            r = c.loc[key]
            name = LABEL[key]
            rows.append({"name": name, "d": r.D_points, "lo": r.D_ci_low, "hi": r.D_ci_high, "raw": False})
        else:
            r = raw[raw.name == key].iloc[0]
            name = key.replace(" (raw)", "")
            rows.append({"name": name, "d": r.d, "lo": r.lo, "hi": r.hi, "raw": True})
        y_of[name] = y
        y -= 1
    y -= 0.35
fig, ax = plt.subplots(figsize=(4.6, 4.3))
forest(ax, rows, y_of, "Next amino acid: real minus shuffled accuracy (percentage points)")
ax.set_yticks(list(y_of.values()), list(y_of.keys()))
for t, m in zip(ax.get_yticklabels(), y_of):
    if m == "Claude Opus 5.5":
        t.set_fontweight("bold")
for g, gy in ticks:
    ax.text(-0.02, gy - 0.25, g.upper(), transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6, color=MUTED, fontweight="bold")
ax.set_xlim(-3, 32)
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=GREY, ms=4.2, label="chat prompt"),
                   Line2D([], [], marker="o", ls="", color=GREY, mfc="white", ms=4.2, mew=1.0, label="raw continuation (completion endpoint or prefill)")],
          loc="lower right", handletextpad=0.3)
header(fig, ax, "Only Claude Opus 5.5 predicts obscure bacterial proteins well",
       "668 bacterial and archaeal proteins from 167 unknown-function (DUF) families, no relative at ≥90% identity;\n"
       "points: mean over families, bars: 95% bootstrap CI. Claude Sonnet 4.5 refused 97% of prefill requests (not shown).")
save(fig, "fig1_protein")

# Fig 2: DNA by region ------------------------------------------------------------------------
t = pd.read_csv("results/dna_confirm_tests.tsv", sep="\t")
cells = [("human cds", "Human coding (CDS)"), ("mouse cds", "Mouse coding (CDS)"), ("bacteria/archaea cds", "Bacterial coding (CDS)"),
         ("human promoter", "Human promoter"), ("mouse promoter", "Mouse promoter"), ("human intron", "Human intron"),
         ("mouse intron", "Mouse intron"), ("human intergenic", "Human intergenic"), ("mouse intergenic", "Mouse intergenic")]
models = [("claude-opus-5-5", 0.22), ("claude-opus-5", 0.0), ("gpt-5.5", -0.22)]
fig, ax = plt.subplots(figsize=(4.6, 3.6))
for i, (cell, name) in enumerate(cells):
    y = -i - (0.5 if i >= 3 else 0)
    for m, dy in models:
        r = t[(t.model == m) & (t.cell == cell)].iloc[0]
        col = color(LABEL[m])
        ax.plot([r.lo, r.hi], [y + dy, y + dy], color=col, lw=1.3, solid_capstyle="round")
        ax.plot(r.D, y + dy, "o", ms=3.6, color=col, mec="white", mew=0.4)
        if m == "claude-opus-5-5" and r.exact20_real:
            ax.text(r.hi + 1, y + dy, f"{int(r.exact20_real)} exact 20-nt runs", va="center", fontsize=6, color=ACCENT)
ax.set_yticks([-i - (0.5 if i >= 3 else 0) for i in range(len(cells))], [n for _, n in cells])
ax.axvline(0, color=MUTED, lw=0.6)
ax.set_xlim(-12, 46)
ax.set_xlabel("Next nucleotide: real minus dinucleotide-shuffled accuracy (percentage points)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color=color(LABEL[m]), ms=4, label=LABEL[m]) for m, _ in models],
          loc="lower right", handletextpad=0.3)
header(fig, ax, "On DNA, Claude Opus 5.5 predicts human and mouse coding sequence only",
       "500 repeat-free windows per region (hg38, mm39, and the genes of the proteins in Fig. 1); 200 nt shown, next nucleotide scored;\n"
       "95% bootstrap CI over windows. An exact 20-nt continuation has a chance probability of about 1e-12.")
save(fig, "fig2_dna")

# Fig 3: well-known vs obscure proteins -------------------------------------------------------
k = pd.read_csv("results/common_vs_obscure.tsv", sep="\t")
wide = k.pivot(index="model", columns="sub", values="d")
order = wide.sort_values("well-known").index.tolist()
fig, ax = plt.subplots(figsize=(4.6, 4.4))
for y, m in enumerate(order):
    a, b = wide.loc[m, "well-known"], wide.loc[m, "obscure"]
    col = color(m)
    ax.plot([b, a], [y, y], color=LIGHT_GREY, lw=1.2, zorder=1)
    ax.plot(a, y, "o", ms=4.2, color="white", mec=col, mew=1.1, zorder=2)
    ax.plot(b, y, "o", ms=4.2, color=col, mec=col, zorder=3)
ax.set_yticks(range(len(order)), [m + (" *" if k[k.model == m].guess.iloc[0] else "") for m in order])
for tl in ax.get_yticklabels():
    if tl.get_text().startswith("Claude Opus 5.5"):
        tl.set_fontweight("bold")
ax.axvline(0, color=MUTED, lw=0.6)
ax.set_xlabel("Next amino acid: real minus shuffled accuracy (percentage points)")
ax.legend(handles=[Line2D([], [], marker="o", ls="", color="white", mec=MUTED, mew=1.1, ms=4.2, label="well-known: vertebrate proteins, well-studied families"),
                   Line2D([], [], marker="o", ls="", color=MUTED, ms=4.2, label="obscure: non-vertebrate DUF proteins")],
          loc="lower right", handletextpad=0.3)
header(fig, ax, "Many models know famous proteins; only recent Claude Opus models know obscure ones",
       "Pilot set, 37 Pfam families (16 well-known, 17 obscure); means over families. * best-guess prompt (the model declined the plain one).")
save(fig, "fig3_known_vs_obscure")

# Fig 4: letter-for-letter recitation of well-known proteins ----------------------------------
s = pd.read_csv("results/recall_length_summary.tsv", sep="\t").set_index("model")
s.index = [LABEL.get(m, m) for m in s.index]
s = s.sort_values(">= 100")
fig, ax = plt.subplots(figsize=(4.6, 2.9))
shades = {">= 50": "#E6E9EF", ">= 100": "#C5CBD6", ">= 200": "#8B93A3"}
ys = range(len(s))
for col_name, shade in shades.items():
    ax.barh(list(ys), s[col_name], height=0.62, color=[color(m) if col_name == ">= 200" and m.startswith("Claude") else shade for m in s.index],
            edgecolor="white", linewidth=0.4, label=f"exact run {col_name.replace('>=', '≥')} residues", zorder=2)
for y, m in zip(ys, s.index):
    ax.text(s.loc[m, ">= 50"] + 0.8, y, f"longest {s.loc[m, 'longest']}", va="center", fontsize=6, color=MUTED)
ax.set_yticks(list(ys), list(s.index))
ax.set_xlabel("Cut points continued letter for letter (of 158)")
ax.legend(loc="lower right", handletextpad=0.4, reverse=True)
header(fig, ax, "Several models recite hundreds of residues of famous proteins",
       "Cut points that any model continued 10/10 in the pilot; asked for the next 200 residues. Bars overlap: each darker bar is a subset; Claude's \u2265200 segment in colour.")
save(fig, "fig4_recall")

# Fig 5: the gap by family size -----------------------------------------------------------------
seqs = pd.read_csv("data/confirm_sequences.tsv", sep="\t")
bins = ["300-1000", "1000-3000", "3000-10000", "10000-max"]
bin_of = seqs.groupby("pfam").size_bin.first()
papers = pd.read_csv("results/confirm_family_opus55.tsv", sep="\t").set_index("pfam").epmc


def family_gaps(model):
    d = pd.DataFrame(json.loads(line) for line in open(f"results/confirm/{model}.jsonl"))
    d = d[d.valid.astype(bool)]
    fam = d.groupby(["pfam", "condition"]).first_correct.mean().unstack().dropna()
    return (fam.real - fam.shuffled) * 100


def by_bin(g, rng=np.random.default_rng(0)):
    out = []
    for b in bins:
        x = g[bin_of.reindex(g.index) == b].to_numpy()
        boot = x[rng.integers(0, len(x), (10_000, len(x)))].mean(1)
        out.append((x.mean(), *np.percentile(boot, [2.5, 97.5])))
    return np.array(out)


fig, ax = plt.subplots(figsize=(4.6, 3.0))
xs = np.arange(len(bins))
ends = []
order = ["gpt-5.5", "gpt-6-astra", "gemini-4-carbon", "claude-opus-4-6", "claude-opus-4-8", "claude-opus-5", "claude-opus-5-5"]
for m in order:
    g = family_gaps(m)
    r = by_bin(g)
    col, top = color(LABEL[m]), m == "claude-opus-5-5"
    if top:
        ax.fill_between(xs, r[:, 1], r[:, 2], color=col, alpha=0.15, lw=0)
    ax.plot(xs, r[:, 0], "-o", color=col, lw=1.6 if top else 1.0, ms=3.6 if top else 2.8, zorder=3 if top else 2)
    ends.append([r[-1, 0], m])
    if top:
        z = by_bin(g[papers.reindex(g.index) == 0])
        ax.plot(xs, z[:, 0], "--", color=col, lw=1.0, zorder=3)
        ax.text(xs[-1] - 0.06, z[-1, 0] + 0.6, "zero-paper families", fontsize=6, color=ACCENT, ha="right", va="bottom")
ends.sort()
orig = [y for y, _ in ends]
for i in range(1, len(ends)):  # keep end labels at least 2.6 points apart ...
    ends[i][0] = max(ends[i][0], ends[i - 1][0] + 2.6)
start = 0
for i in range(1, len(ends) + 1):  # ... and centre each pushed-apart run on its lines
    if i == len(ends) or ends[i][0] - ends[i - 1][0] > 2.61:
        shift = np.mean(orig[start:i]) - np.mean([y for y, _ in ends[start:i]])
        for e in ends[start:i]:
            e[0] += shift
        start = i
for y, m in ends:
    top = m == "claude-opus-5-5"
    ax.text(xs[-1] + 0.12, y, LABEL[m], va="center", fontsize=6.2, color=INK if top else MUTED, fontweight="bold" if top else "normal")
ax.axhline(0, color=MUTED, lw=0.6)
ax.set_ylim(-6.5, 47)
ax.set_xticks(xs, ["300–1k", "1k–3k", "3k–10k", ">10k"])
ax.set_xlim(-0.25, len(bins) - 1 + 1.25)
ax.set_xlabel("Proteins in the family (Pfam)")
ax.set_ylabel("Real minus shuffled (points)")
ax.grid(axis="y", color=GRID, lw=0.5)
ax.grid(axis="x", visible=False)
header(fig, ax, "Claude Opus 5.5's gap grows with the number of sequences in a family",
       "Next amino acid on the 167 DUF families of Fig. 1, by family size (50, 50, 50, 17 families); mean over families,\n"
       "band: 95% bootstrap CI for Opus 5.5. Dashed: the 128 families with no Europe PMC hit.")
save(fig, "fig5_family_size")
print("figures/fig1_protein, fig2_dna, fig3_known_vs_obscure, fig4_recall, fig5_family_size")
