"""Interim figures for a quick look (not the report): results/interim_*.png. Free."""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

LABEL = {"gpt-4o-2024-11-20": "GPT-4o", "gpt-4.1-2025-04-14": "GPT-4.1", "gpt-5.5": "GPT-5.5", "gpt-6.1-sol": "GPT-6.1 Sol", "gpt-6-astra": "GPT-6 Astra",
         "claude-sonnet-4-5-20250929": "Claude Sonnet 4.5", "claude-opus-4-6": "Claude Opus 4.6", "claude-opus-5": "Claude Opus 5", "claude-opus-5-5": "Claude Opus 5.5",
         "gemini-2.5-pro": "Gemini 2.5 Pro", "gemini-3.5-flash": "Gemini 3.5 Flash", "gemini-3.8-flash": "Gemini 3.8 Flash", "gemini-4-carbon": "Gemini 4 Carbon",
         "deepseek-v4p1-flash": "DeepSeek V4.1 Flash", "qwen3p8-2p4t-a95b": "Qwen 3.8", "kimi-k3p1-pro": "Kimi K3.1 Pro", "glm-5p3": "GLM-5.3",
         "minimax-m3": "MiniMax M3", "gpt-oss-120b": "gpt-oss-120b", "nemotron-3-ultra-nvfp4": "Nemotron 3 Ultra"}
OPEN = {"deepseek-v4p1-flash", "qwen3p8-2p4t-a95b", "kimi-k3p1-pro", "glm-5p3", "minimax-m3", "gpt-oss-120b", "nemotron-3-ultra-nvfp4"}
BLUE, GREY, INK, MUTED = "#2a78d6", "#a9a8a2", "#1a1a1a", "#6b6a66"
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED})


def forest(ax, labels, d, lo, hi, hot):
    y = range(len(labels))[::-1]
    for yi, l, a, b, c, h in zip(y, labels, d, lo, hi, hot):
        col = BLUE if h else GREY
        ax.plot([b, c], [yi, yi], color=col, lw=2, solid_capstyle="round")
        ax.plot(a, yi, "o", color=col, ms=6, mec="white", mew=0.8)
        ax.text(c + 0.6, yi, f"{a:+.1f}", va="center", fontsize=8, color=INK if h else MUTED)
    ax.set_yticks(list(y), labels)
    ax.axvline(0, color=MUTED, lw=0.8)
    ax.grid(axis="x", color="#e6e5e1", lw=0.6)
    ax.set_axisbelow(True)


# A: confirmation, interim
c = pd.read_csv("results/confirm_tests.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(8.5, 3.8))
labels = [f"{LABEL[m]}  ({f}/167 families)" for m, f in zip(c.model, c.families)]
forest(ax, labels, c.D_points, c.D_ci_low, c.D_ci_high, c.model.eq("claude-opus-5-5"))
ax.axvline(10, color=MUTED, lw=1, ls="--")
ax.text(10.3, -0.45, "pre-registered bar:\nCI lower bound > 10", fontsize=8, color=MUTED, va="bottom")
ax.set_xlabel("Next amino acid: real minus shuffled accuracy (percentage points), 95% CI over families")
ax.set_title("PROTEIN, confirmation run (interim): bacterial/archaeal DUF proteins, no 90% relative,\n"
             "cuts inside the Pfam domain", loc="left", fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig("results/interim_confirmation.png", dpi=150, bbox_inches="tight")

# B: pilot, all models, non-vertebrate DUFs
p = pd.read_csv("results/pilot_all.tsv", sep="\t")
p = p.sort_values("prompt").drop_duplicates("model", keep="first")  # "guess" sorts before "plain": prefer it where run
order = [m for m in LABEL if m in set(p.model)]
p = p.set_index("model").loc[order]
labels = [f"{LABEL[m]}{' *' if r.prompt == 'guess' else ''}{' (open)' if m in OPEN else ''}{' (partial)' if r.probes < 2220 else ''}" for m, r in p.iterrows()]
fig, ax = plt.subplots(figsize=(8.5, 6.4))
forest(ax, labels, p.diff_nvduf, p.lo_nvduf, p.hi_nvduf, [m.startswith("claude-opus-5") for m in p.index])
ax.set_xlabel("Next amino acid: real minus shuffled accuracy (percentage points), 95% CI over families")
ax.set_title("PROTEIN, pilot: 17 non-vertebrate DUF families (bacteria, plants, fungi),\n"
             "cuts inside the Pfam domain; * best-guess prompt", loc="left", fontsize=10, fontweight="bold")
fig.tight_layout()
fig.savefig("results/interim_pilot_models.png", dpi=150, bbox_inches="tight")
print("results/interim_confirmation.png results/interim_pilot_models.png")
