"""Recall-length figure: how many common cases each model recites exactly for >= N residues
-> results/fig_recall_length.png. Free."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

LABEL = {"claude-opus-5-5": "Claude Opus 5.5", "claude-opus-5": "Claude Opus 5", "claude-opus-4-8": "Claude Opus 4.8",
         "gemini-4-carbon": "Gemini 4 Carbon", "gemini-3.8-flash": "Gemini 3.8 Flash", "gemini-3.5-flash": "Gemini 3.5 Flash",
         "gpt-6.1-sol": "GPT-6.1 Sol", "gpt-5.5": "GPT-5.5"}
THRESH = [20, 50, 100, 200]
rows = [json.loads(line) for p in Path("results/recall_length").glob("*.jsonl") for line in p.open()]
df = pd.DataFrame(rows)
df = df[df.exact_run.notna()]
n_cases = df.groupby("model").id.nunique()
models = [m for m in LABEL if m in set(df.model) and n_cases[m] >= 150]
tab = pd.DataFrame({f">= {t}": df.groupby("model").exact_run.apply(lambda s, t=t: int((s >= t).sum())) for t in THRESH}).loc[models]
tab["longest"] = df.groupby("model").exact_run.max().loc[models].astype(int)
print(tab.to_string())

fig, ax = plt.subplots(figsize=(8.6, 4.4))
w = 0.2
shades = ["#cfe0f5", "#86b6ef", "#3987e5", "#184f95"]
for i, (t, c) in enumerate(zip(THRESH, shades)):
    ys = [k + (i - 1.5) * w for k in range(len(models))][::-1]
    ax.barh(ys, tab[f">= {t}"], height=w, color=c, label=f"exact run >= {t} amino acids")
for k, m in enumerate(models[::-1]):
    ax.text(tab.loc[m, ">= 20"] + 1.5, k, f"longest {tab.loc[m, 'longest']}", va="center", fontsize=8, color="#52514e")
ax.set_yticks(range(len(models)), [LABEL[m] for m in models[::-1]])
ax.set_xlabel(f"Cases recited letter for letter (of {int(n_cases.loc[models].max())} common cut points)")
ax.spines[["top", "right"]].set_visible(False)
ax.legend(frameon=False, fontsize=8, loc="lower right")
ax.set_title("PROTEIN (pilot cut points any model got 10/10 right): asked for the next 200 amino acids,\n"
             "how far does each model continue exactly?", loc="left", fontsize=10, fontweight="bold")
Path("results/figures").mkdir(parents=True, exist_ok=True)
fig.savefig("results/figures/fig_recall_length.png", dpi=150, bbox_inches="tight")
tab.to_csv("results/recall_length_summary.tsv", sep="\t")
