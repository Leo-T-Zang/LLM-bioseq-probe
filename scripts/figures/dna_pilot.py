"""DNA pilot figure: real minus dinucleotide-shuffled accuracy of the next nucleotide,
per model and region -> results/fig_dna.{png,pdf}. Free. Reads results/dna_summary.tsv."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

LABEL = {"claude-opus-5-5": "Claude Opus 5.5", "claude-opus-5": "Claude Opus 5", "gpt-5.5": "GPT-5.5",
         "deepseek-v4p1-flash": "DeepSeek V4.1 Flash *", "gemini-4-carbon": "Gemini 4 Carbon"}
REGIONS = [("human", "cds"), ("mouse", "cds"), ("bacteria/archaea", "cds"), ("human", "utr5"), ("mouse", "utr5"),
           ("human", "utr3"), ("mouse", "utr3"), ("human", "promoter"), ("mouse", "promoter"), ("human", "lncrna"),
           ("mouse", "lncrna"), ("human", "intron"), ("mouse", "intron"), ("human", "intergenic"), ("mouse", "intergenic")]
NAME = {"cds": "coding (CDS)", "utr5": "5' UTR", "utr3": "3' UTR", "promoter": "promoter (1 kb)", "lncrna": "lncRNA exon",
        "intron": "intron", "intergenic": "intergenic"}
BLUE, GREY, MUTED = "#2a78d6", "#a9a8a2", "#6b6a66"
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED})

s = pd.read_csv("results/dna_summary.tsv", sep="\t")
models = [m for m in LABEL if m in set(s.model)]
fig, axes = plt.subplots(1, len(models), figsize=(2.3 * len(models) + 1.6, 5.6), sharey=True, sharex=True)
ys = list(range(len(REGIONS)))[::-1]
for ax, m in zip(axes, models):
    for y, (sp, rg) in zip(ys, REGIONS):
        r = s[(s.model == m) & (s.species == sp) & (s.region == rg)]
        if r.empty or r.loci.iloc[0] < 30:
            ax.text(0, y, f"n={0 if r.empty else int(r.loci.iloc[0])}", ha="center", va="center", fontsize=7, color=MUTED)
            continue
        r = r.iloc[0]
        sig = r.lo > 0
        c = BLUE if sig else GREY
        ax.plot([r.lo, r.hi], [y, y], color=c, lw=1.8, solid_capstyle="round")
        ax.plot(r["diff"], y, "o", color=c, ms=5, mec="white", mew=0.6)
        if r.exact20_real:
            ax.text(r.hi + 2, y, f"{int(r.exact20_real)}x20", fontsize=7, va="center", color="#b23b3b")
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.grid(axis="x", color="#e6e5e1", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(LABEL[m], fontsize=9, fontweight="bold")
    ax.set_xlim(-35, 55)
axes[0].set_yticks(ys, [f"{sp.split('/')[0]}: {NAME[rg]}" for sp, rg in REGIONS])
for y in (11.5, 3.5):
    for ax in axes:
        ax.axhline(y, color="#d9d8d3", lw=0.6)
fig.supxlabel("Next nucleotide: real minus dinucleotide-shuffled accuracy (percentage points), 95% CI over ~100 windows", fontsize=8.5)
fig.suptitle("DNA (pilot): next-nucleotide prediction after a 200-nt prefix, repeat-free windows from hg38, mm39 and bacterial genes\n"
             "blue = CI above 0; red Nx20 = exact 20-nt continuations (chance ~1e-12); Gemini 4 Carbon: most prompts blocked; * best-guess prompt",
             x=0.01, ha="left", fontsize=9.5, fontweight="bold")
fig.tight_layout(rect=(0, 0.02, 1, 0.93))
Path("results/figures").mkdir(parents=True, exist_ok=True)
fig.savefig("results/figures/fig_dna.png", dpi=150, bbox_inches="tight")
fig.savefig("results/figures/fig_dna.pdf", bbox_inches="tight")
print("results/figures/fig_dna.png")
