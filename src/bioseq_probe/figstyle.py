"""One look for every figure: journal-sized type, thin lines, grey by default, indigo for the
Claude models (darker for newer generations)."""

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
ACCENT = "#3D5A98"  # indigo, shared with the method diagram
CLAUDE = {  # older generations lighter
    "Claude Sonnet 4.5": "#D6DEEE", "Claude Opus 4.6": "#BCC8E2", "Claude Opus 4.8": "#93A6CF",
    "Claude Opus 5": "#6A82B8", "Claude Opus 5.5": ACCENT,
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
    return CLAUDE.get(model, GREY)


def header(fig, ax, title: str, subtitle: str):
    """Bold one-line finding, then grey lines saying exactly what is plotted, set just above the axes."""
    fig.canvas.draw()
    top = ax.get_position().y1
    h = fig.get_figheight() * 72  # points
    n = subtitle.count("\n") + 1
    sub_y = top + 10 / h
    fig.text(0.0, sub_y, subtitle, ha="left", va="bottom", fontsize=6.5, color=MUTED, linespacing=1.35)
    fig.text(0.0, sub_y + (n * 6.5 * 1.35 + 4) / h, title, ha="left", va="bottom", fontsize=8.5, fontweight="bold", color=INK)


def save(fig, stem: str):
    Path("figures").mkdir(exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(f"figures/{stem}.{ext}", bbox_inches="tight", pad_inches=0.04, facecolor="white")
