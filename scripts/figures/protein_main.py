"""Main protein figure: every model on the confirmation set, chat-prompt models and
raw-continuation models together -> results/figures/fig_protein_main.png. Free."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from bioseq_probe.probe import parse_raw

LABEL = {"claude-sonnet-4-5-20250929": "Claude Sonnet 4.5", "gemini-2.5-pro": "Gemini 2.5 Pro", "deepseek-v4p1-flash": "DeepSeek V4.1 Flash",
         "qwen3p8-2p4t-a95b": "Qwen 3.8", "kimi-k3p1-pro": "Kimi K3.1 Pro", "glm-5p3": "GLM-5.3", "minimax-m3": "MiniMax M3",
         "gpt-oss-120b": "gpt-oss-120b", "nemotron-3-ultra-nvfp4": "Nemotron 3 Ultra"}
CHAT = {"claude-opus-5-5": "Claude Opus 5.5", "claude-opus-5": "Claude Opus 5", "claude-opus-4-8": "Claude Opus 4.8", "claude-opus-4-6": "Claude Opus 4.6",
        "gpt-6-astra": "GPT-6 Astra", "gpt-5.5": "GPT-5.5", "gemini-4-carbon": "Gemini 4 Carbon"}
BLUE, ORANGE, GREY, MUTED = "#2a78d6", "#eb6834", "#a9a8a2", "#6b6a66"
rng = np.random.default_rng(0)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": MUTED})


def family_diff(df):
    fam = df.groupby(["pfam", "condition"])["first"].mean().unstack().dropna()
    D = (fam["real"] - fam["shuffled"]).to_numpy()
    lo, hi = np.percentile(D[rng.integers(0, len(D), (4000, len(D)))].mean(1), [2.5, 97.5])
    return D.mean() * 100, lo * 100, hi * 100, len(D)


rows = []
for m, name in LABEL.items():
    p = Path(f"results/raw/{m}.jsonl")
    if not p.exists():
        continue
    recs = [json.loads(line) for line in p.open()]
    for r in recs:  # re-scored from the text, so parser fixes apply
        r["pred"] = parse_raw(r.get("text", "") or "")
    ok = [r for r in recs if "error" not in r and r.get("stop") != "refusal" and r["pred"]]
    refused = sum(r.get("stop") == "refusal" for r in recs)
    if len(ok) < 200:
        rows.append({"name": f"{name} (raw)", "kind": "raw", "d": np.nan, "lo": np.nan, "hi": np.nan, "note": f"{refused}/{len(recs)} refused"})
        continue
    df = pd.DataFrame([{"pfam": r["pfam"], "condition": r["condition"], "first": r["pred"][0] == r["target"][0]} for r in ok])
    d, lo, hi, nfam = family_diff(df)
    rows.append({"name": f"{name} (raw)", "kind": "raw", "d": d, "lo": lo, "hi": hi, "note": f"{len(recs)}/6680 probes"})
c = pd.read_csv("results/confirm_tests.tsv", sep="\t").set_index("model")
for m, name in CHAT.items():
    r = c.loc[m]
    rows.append({"name": f"{name} (chat)", "kind": "chat", "d": r.D_points, "lo": r.D_ci_low, "hi": r.D_ci_high, "note": "6680/6680 probes"})
res = pd.DataFrame(rows)
res = pd.concat([res[res.kind == "chat"], res[res.kind == "raw"]], ignore_index=True)

fig, ax = plt.subplots(figsize=(8.8, 5.6))
ys = list(range(len(res)))[::-1]
for y, r in zip(ys, res.itertuples()):
    if np.isnan(r.d):
        ax.text(0.5, y, r.note, va="center", fontsize=8, color="#b23b3b")
        continue
    col = BLUE if r.name.startswith("Claude Opus 5.5") else (ORANGE if r.kind == "raw" else GREY)
    ax.plot([r.lo, r.hi], [y, y], color=col, lw=2, solid_capstyle="round")
    ax.plot(r.d, y, "o", color=col, ms=6, mec="white", mew=0.8)
    ax.text(r.hi + 0.6, y, f"{r.d:+.1f}" + ("" if r.note.startswith("6680/6680") else f"   ({r.note})"), va="center", fontsize=7.5, color=MUTED)
ax.set_yticks(ys, res.name)
ax.axvline(0, color=MUTED, lw=0.8)
ax.grid(axis="x", color="#e6e5e1", lw=0.6)
ax.set_axisbelow(True)
ax.set_xlabel("Next amino acid: real minus shuffled accuracy (percentage points), 95% CI over families")
ax.set_title("PROTEIN: next amino acid on 668 obscure bacterial/archaeal proteins (167 DUF families, no 90% relative)\n"
             "chat prompt (grey, blue) or raw continuation via completion endpoint or prefill (orange)", loc="left", fontsize=10, fontweight="bold")
Path("results/figures").mkdir(parents=True, exist_ok=True)
fig.savefig("results/figures/fig_protein_main.png", dpi=150, bbox_inches="tight")
fig.savefig("results/figures/fig_protein_main.pdf", bbox_inches="tight")
res.to_csv("results/raw_summary.tsv", sep="\t", index=False)
print(res.round(1).to_string(index=False))
