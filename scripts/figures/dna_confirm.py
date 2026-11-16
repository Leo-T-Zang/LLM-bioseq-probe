"""DNA confirmation figure from results/dna_confirm_tests.tsv -> results/fig_dna_confirm.{png,pdf}. Free."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

LABEL = {"claude-opus-5-5": "Claude Opus 5.5", "claude-opus-5": "Claude Opus 5", "gpt-5.5": "GPT-5.5"}
CELLS = ["human cds", "mouse cds", "bacteria/archaea cds", "human promoter", "mouse promoter",
         "human intron", "mouse intron", "human intergenic", "mouse intergenic"]
NAME = {"cds": "coding (CDS)", "promoter": "promoter (1 kb)", "intron": "intron", "intergenic": "intergenic"}
BLUE, GREY, MUTED = "#2a78d6", "#a9a8a2", "#6b6a66"
plt.rcParams.update({"font.size": 8.5, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED})

t = pd.read_csv("results/dna_confirm_tests.tsv", sep="\t")
models = [m for m in LABEL if m in set(t.model)]
done = t.groupby("model").windows.sum()
final = all(done[m] >= 4400 for m in models)
fig, axes = plt.subplots(1, len(models), figsize=(2.9 * len(models) + 1.8, 4.2), sharey=True, sharex=True)
ys = list(range(len(CELLS)))[::-1]
for ax, m in zip(axes, models):
    for y, cell in zip(ys, CELLS):
        r = t[(t.model == m) & (t.cell == cell)]
        if r.empty:
            continue
        r = r.iloc[0]
        c = BLUE if r.D1_pass else GREY
        ax.plot([r.lo, r.hi], [y, y], color=c, lw=2, solid_capstyle="round")
        ax.plot(r.D, y, "o", color=c, ms=5.5, mec="white", mew=0.6)
        if r.exact20_real:
            ax.text(r.hi + 1.5, y, f"{int(r.exact20_real)}x20", fontsize=7.5, va="center", color="#b23b3b")
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.axvline(5, color=MUTED, lw=0.8, ls="--")
    ax.grid(axis="x", color="#e6e5e1", lw=0.6)
    ax.set_axisbelow(True)
    ax.set_title(f"{LABEL[m]}  ({int(done[m])} windows)", fontsize=9, fontweight="bold")
    ax.axhline(5.5, color="#d9d8d3", lw=0.6)
axes[0].set_yticks(ys, [f"{c.split()[0].split('/')[0]}: {NAME[c.split()[1]]}" for c in CELLS])
fig.supxlabel("Next nucleotide: real minus dinucleotide-shuffled accuracy (points), 95% CI over windows; dashed = threshold (5 points)", fontsize=8.5)
fig.suptitle(f"DNA (confirmation, {'final' if final else 'interim'}): 500 new repeat-free 220-nt windows per region, hg38 / mm39 / bacterial genes\n"
             "blue = passes D1 (p < 0.001, CI low > 5); red Nx20 = exact 20-nt continuations (chance ~1e-12)",
             x=0.01, ha="left", fontsize=9.5, fontweight="bold")
fig.tight_layout(rect=(0, 0.02, 1, 0.9))
Path("results/figures").mkdir(parents=True, exist_ok=True)
fig.savefig("results/figures/fig_dna_confirm.png", dpi=150, bbox_inches="tight")
fig.savefig("results/figures/fig_dna_confirm.pdf", bbox_inches="tight")
print("results/figures/fig_dna_confirm.png")
