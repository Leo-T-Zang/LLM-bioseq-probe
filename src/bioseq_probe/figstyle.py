"""One look for every figure: journal-sized type, thin lines, and each lab's own colour
(Claude models darker for newer generations, open-weight models grey). Checked for
colour-blind separation; every model is also named on the chart."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt

INK = "#1D2433"
MUTED = "#6B7385"
GRID = "#ECEEF2"
GREY = "#A3AAB8"
LIGHT_GREY = "#D3D7DF"
ACCENT = "#D97757"  # Anthropic
ACCENT_TEXT = "#B4532F"  # the same hue, dark enough for text
OPENAI = "#4B4B55"
GPT = {"GPT-6.1 Sol": "#2C2C33", "GPT-6 Astra": OPENAI, "GPT-5.5": "#7C7C87"}  # newer darker
GOOGLE = "#4285F4"
PROTEIN_LM = "#9B59B6"  # protein language models, drawn as diamonds
CLAUDE = {  # older generations lighter
    "Claude Sonnet 4.5": "#F6DFD4", "Claude Opus 4.5": "#F3D3C4", "Claude Opus 4.6": "#EFC1AB",
    "Claude Opus 4.7": "#EAAD92", "Claude Opus 4.8": "#E59C7D", "Claude Opus 5": "#DF8A68", "Claude Opus 5.5": ACCENT, "Claude Fable 5.1": "#C0694C", "Claude Fable 5.2": "#9C3A1C",
}

for path in Path("/usr/share/fonts").rglob("LiberationSans-*.ttf"):
    fm.fontManager.addfont(str(path))

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 7, "legend.fontsize": 6.5,
    "axes.linewidth": 0.6, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": INK, "xtick.major.width": 0.6, "xtick.major.size": 2.5,
    "ytick.major.size": 0, "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
    "axes.grid": True, "axes.grid.axis": "x", "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "legend.frameon": False, "figure.dpi": 150, "savefig.dpi": 300, "pdf.fonttype": 42, "svg.fonttype": "none",
})


def color(model: str) -> str:
    """By display name: Claude by generation, GPT charcoal, Gemini blue, open-weight grey."""
    if model.startswith("GPT-"):
        return GPT.get(model, OPENAI)
    if model.startswith("Gemini"):
        return GOOGLE
    if model.startswith(("ProtGPT", "ESM")):
        return PROTEIN_LM
    return CLAUDE.get(model, GREY)


def header(fig, ax, title: str, subtitle: str):
    """Bold one-line finding above the axes; the grey note saying exactly what is plotted goes
    below the axes, right-aligned, like a caption."""
    fig.canvas.draw()
    h = fig.get_figheight() * 72  # points
    fig.text(0.0, ax.get_position().y1 + 10 / h, title, ha="left", va="bottom", fontsize=8.5, fontweight="bold", color=INK)
    box = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.transFigure.inverted())
    fig.text(box.x1, box.y0 - 6 / h, subtitle, ha="right", va="top", fontsize=6.2, color=MUTED, linespacing=1.35,
             multialignment="right")


def save(fig, stem: str):
    Path("figures").mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"figures/{stem}.{ext}", bbox_inches="tight", pad_inches=0.04, facecolor="white")
