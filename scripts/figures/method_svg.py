"""The test as one picture: figures/fig0_method.svg and fig0_method_mobile.svg. Free.

Two rows of residues: the end of a real protein, and the same residues shuffled (thin lines
show where each one went). The model fills in the next residue; the accuracies are the main
set's (results/confirm_tests.tsv). The letters are the residues before cut 109 of A0ABX2F111,
where Claude Opus 5.5 answered E (right); the shuffled row is schematic. Colours are CSS
variables with fallbacks, so a page can restyle them.
"""

from pathlib import Path

import pandas as pd

ACC, CUT = "A0ABX2F111", 109
seq = pd.read_csv("data/confirm_sequences.tsv", sep="\t").set_index("accession").loc[ACC, "sequence"]
assert seq[CUT] == "E"
t = pd.read_csv("results/confirm_tests.tsv", sep="\t").set_index("model").loc["claude-opus-5-5"]
REAL, SHUF, GAP = t.real_acc * 100, t.shuffled_acc * 100, t.D_points

SANS = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
MONO = "'SFMono-Regular', Menlo, Consolas, 'Liberation Mono', monospace"
C = {
    "ink": "#1d2433", "muted": "#6b7385", "faint": "#b4bac6", "line": "#d9dde4", "bead": "#ffffff",
    "accent": "#d97757", "accent-text": "#b4532f", "accent-soft": "#fbeee8",
    "good": "#2f8a5f", "bad": "#c2453d", "chip": "#1d2433", "chip-text": "#ffffff",
}
PERM = {9: [4, 7, 0, 8, 2, 5, 1, 6, 3], 4: [2, 0, 3, 1]}  # shuffled position j holds real residue PERM[j]


def v(name):
    return f"var(--bp-{name}, {C[name]})"


def text(x, y, s, size=13, color="ink", weight=400, anchor="start", font=SANS):
    return (f'<text x="{x:.1f}" y="{y:.1f}" text-anchor="{anchor}" style="font-family:{font};font-size:{size}px;'
            f'font-weight:{weight};fill:{v(color)}">{s}</text>')


def bead(x, y, letter, r, stroke="faint", fill="bead", dash=False, color="ink", weight=500):
    d = ";stroke-dasharray:3 2.5" if dash else ""
    return (f'<circle cx="{x:.1f}" cy="{y}" r="{r}" style="fill:{v(fill)};stroke:{v(stroke)};stroke-width:1.3{d}"/>'
            + text(x, y + r * 0.37, letter, round(r * 1.02, 1), color, weight, "middle", MONO))


def arrow(x0, x1, y):
    return f'<path d="M{x0:.1f},{y} H{x1:.1f}" style="fill:none;stroke:{v("faint")};stroke-width:1.3" marker-end="url(#bp-arr)"/>'


def figure(n, compact=False, step=30, r=12.5):
    """compact: the phone layout, with the model named over a plain arrow."""
    real = seq[CUT - n : CUT]
    shuf = "".join(real[p] for p in PERM[n])
    y1, y2, x0 = 30, 134, 108 if compact else 118
    xs = [x0 + i * step for i in range(n)]
    out = []
    for j, p in enumerate(PERM[n]):  # where each residue went
        a, b = xs[p], xs[j]
        out.append(f'<path d="M{a:.1f},{y1 + r + 4} C{a:.1f},{(y1 + y2) / 2} {b:.1f},{(y1 + y2) / 2} {b:.1f},{y2 - r - 4}" '
                   f'style="fill:none;stroke:{v("line")};stroke-width:1.1"/>')
    q = xs[-1] + step
    for y, letters in ((y1, real), (y2, shuf)):
        out += [bead(x, y, a, r) for x, a in zip(xs, letters)]
        out.append(bead(q, y, "?", r, "accent", "accent-soft", dash=True, color="accent-text", weight=700))
    out.append(text(24, y1 + 5, "Real", 15, "ink", 650))
    out.append(text(24, y2 + 5, "Shuffled", 15, "ink", 650))
    cx, chip = q + r + 14, 82
    for y, ans, ok, acc in ((y1, "E", True, REAL), (y2, "R", False, SHUF)):
        if compact:
            out.append(arrow(cx, cx + 62, y))
            out.append(text(cx + 30, y - 8, "Opus 5.5", 11.5, "muted", 600, "middle"))
            bx = cx + 68 + r
        else:
            out.append(arrow(cx, cx + 20, y))
            out.append(f'<rect x="{cx + 26:.1f}" y="{y - 15}" width="{chip}" height="30" rx="15" style="fill:{v("chip")}"/>')
            out.append(text(cx + 26 + chip / 2, y + 4.5, "Opus 5.5", 12.5, "chip-text", 600, "middle"))
            a = cx + 26 + chip + 6
            out.append(arrow(a, a + 20, y))
            bx = a + 26 + r
        mark = "good" if ok else "bad"
        out.append(bead(bx, y, ans, r, mark, color=mark, weight=700))
        out.append(text(bx + r + 1, y - r + 3, "✓" if ok else "✗", 12, mark, 700))
        out.append(text(bx + r + 18, y + 7.5, f"{acc:.1f}%", 21, "ink" if ok else "muted", 700))
    gx = bx + r + 96
    out.append(f'<path d="M{gx - 7:.1f},{y1} H{gx:.1f} V{y2} H{gx - 7:.1f}" style="fill:none;stroke:{v("accent")};stroke-width:1.4"/>')
    out.append(text(gx + 10, (y1 + y2) / 2 + 3, f"+{GAP:.1f}", 21, "accent-text", 700))
    out.append(text(gx + 10, (y1 + y2) / 2 + 20, "points", 12, "accent-text", 500))
    return gx + 80, y2 + 30, "".join(out)


DEFS = (f'<defs><marker id="bp-arr" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto"><path d="M1,1.5 L8.5,5 L1,8.5" style="fill:none;stroke:{v("faint")};stroke-width:1.6;'
        f'stroke-linecap:round;stroke-linejoin:round"/></marker></defs>')
LABEL = (f"Claude Opus 5.5 gets the next amino acid of a real protein right {REAL:.1f}% of the time, and of the "
         f"same protein shuffled {SHUF:.1f}% of the time: a gap of {GAP:.1f} points.")


def svg(w, h, body):
    return (f'<svg xmlns="http://www.w3.org/2000/svg" class="bp-figure" viewBox="0 0 {w:.0f} {h}" width="{w:.0f}" '
            f'height="{h}" role="img" aria-label="{LABEL}" focusable="false">{DEFS}{body}</svg>\n')


Path("figures").mkdir(exist_ok=True)
Path("figures/fig0_method.svg").write_text(svg(*figure(9)))
Path("figures/fig0_method_mobile.svg").write_text(svg(*figure(4, compact=True)))
print("figures/fig0_method.svg, figures/fig0_method_mobile.svg")
