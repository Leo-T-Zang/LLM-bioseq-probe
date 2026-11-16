"""The test, drawn by hand as SVG: figures/fig0_method.svg (wide) and fig0_method_mobile.svg. Free.

Colours are CSS variables with fallbacks, so a page can restyle them (dark mode).
Panels 2 and 3 show a real question from the main set and Claude Opus 5.5's real replies
(A0ABX2F111, Kibdelosporangium persicum, Pfam PF21962, cut 109; results/confirm/claude-opus-5-5.jsonl).
"""

from pathlib import Path

import pandas as pd

from bioseq_probe.probe import shuffled

ACC, CUT = "A0ABX2F111", 109
seq = pd.read_csv("data/confirm_sequences.tsv", sep="\t").set_index("accession").loc[ACC, "sequence"]
shuf = shuffled(seq, ACC)
REAL_TAIL, REAL_NEXT = seq[CUT - 10 : CUT], seq[CUT : CUT + 10]
SHUF_TAIL, SHUF_NEXT = shuf[CUT - 10 : CUT], shuf[CUT : CUT + 10]
REPLY_REAL, REPLY_SHUF = "ENNLSIANMD", "RVWVGSLDGF"
assert REAL_NEXT[0] == REPLY_REAL[0] and SHUF_NEXT[0] != REPLY_SHUF[0]

SANS = "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Helvetica, Arial, sans-serif"
MONO = "'SFMono-Regular', Menlo, Consolas, 'Liberation Mono', monospace"
C = {  # name: fallback (a page can set --bp-<name>, e.g. for dark mode)
    "panel": "#f6f7f9", "card": "#ffffff", "ink": "#1d2433", "body": "#3c4454", "muted": "#6b7385",
    "faint": "#a3aab8", "line": "#e1e4ea", "accent": "#3d5a98", "accent-text": "#2c4a8a", "accent-soft": "#eaf0fb",
    "good": "#2f8a5f", "good-soft": "#e3f3ea", "bad": "#c2453d", "bad-soft": "#fbe7e5",
    "hero": "#1d2433", "hero-text": "#ffffff", "tile": "#eceef2",
}


def v(name):
    return f"var(--bp-{name}, {C[name]})"


def text(x, y, s, size=11.5, color="body", weight=400, anchor="start", font=SANS):
    return (f'<text x="{x}" y="{y}" style="font-family:{font};font-size:{size}px;font-weight:{weight};'
            f'fill:{v(color)}" text-anchor="{anchor}">{s}</text>')


def rect(x, y, w, h, fill="card", stroke="line", r=8, sw=1):
    st = f"stroke:{v(stroke)};stroke-width:{sw}" if stroke else "stroke:none"
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" style="fill:{v(fill)};{st}"/>'


def arrow(x, y0, y1):
    return f'<path d="M{x},{y0} V{y1}" style="fill:none;stroke:{v("faint")};stroke-width:1.25" marker-end="url(#bp-arr)"/>'


def tiles(x, y, letters, first=None, w=15.8, gap=1.8):
    """One residue per tile; first='right'/'wrong' marks the scored residue."""
    out = []
    for i, a in enumerate(letters):
        fill, stroke = "tile", None
        if i == 0 and first == "right":
            fill, stroke = "good-soft", "good"
        if i == 0 and first == "wrong":
            fill, stroke = "bad-soft", "bad"
        out.append(rect(x + i * (w + gap), y, w, 20, fill, stroke, r=3, sw=1.25))
        out.append(text(x + i * (w + gap) + w / 2, y + 14.5, a, 11.5, "ink", 500, "middle", MONO))
    return "".join(out)


def badge(x, y, n, title):
    return (f'<circle cx="{x + 10}" cy="{y - 4}" r="10" style="fill:{v("accent")}"/>'
            + text(x + 10, y, n, 11.5, "hero-text", 700, "middle") + text(x + 28, y, title, 13, "ink", 650))


def panel1():
    p = [badge(0, 16, "1", "Proteins nobody studies")]
    rows = [(224, ["Pfam families named", "“domain of unknown function”"]),
            (200, ["no solved structure,", "≤3 papers, ≥300 proteins"]),
            (176, ["bacteria and archaea; no other", "UniProt protein ≥90% identical"])]
    y = 36
    for w, lines in rows:
        x = (224 - w) / 2
        p.append(rect(x, y, w, 44))
        p.append(text(112, y + 18, lines[0], 11.5, "body", 400, "middle"))
        p.append(text(112, y + 34, lines[1], 11.5, "body", 400, "middle"))
        p.append(arrow(112, y + 45, y + 57))
        y += 60
    p.append(rect(32, y, 160, 34, "accent-soft", "accent", sw=1.25))
    p.append(text(112, y + 22, "668 proteins, 167 families", 12, "ink", 650, "middle"))
    p.append(text(0, y + 62, "No paper describes them, so", 11.5, "muted"))
    p.append(text(0, y + 78, "reasoning and literature cannot", 11.5, "muted"))
    p.append(text(0, y + 94, "predict them. Training data can.", 11.5, "muted"))
    return "".join(p)


def panel2():
    p = [badge(0, 16, "2", "Ask for the next amino acids")]
    p.append(rect(0, 32, 224, 98))
    p.append(text(12, 52, "Predict how this protein continues.", 11, "muted"))
    p.append(text(12, 67, "Reply with the next 10 amino acids.", 11, "muted"))
    p.append(text(12, 93, "&gt;sequence", 11.5, "ink", 400, "start", MONO))
    p.append(text(12, 111, f"…{seq[CUT - 23 : CUT]}", 11.5, "ink", 400, "start", MONO))
    p.append(rect(184, 98, 16, 19, "accent-soft", "accent", r=3, sw=1.25))
    p.append(text(192, 112, "?", 11.5, "accent-text", 700, "middle", MONO))
    p.append(arrow(112, 131, 147))
    p.append(rect(52, 150, 120, 28, "hero", None, r=14))
    p.append(text(112, 168, "Claude Opus 5.5", 11.5, "hero-text", 650, "middle"))
    p.append(arrow(112, 179, 195))
    p.append(text(0, 213, "reply", 11, "muted"))
    p.append(text(36, 214, "\u2713", 13, "good", 700, "middle"))
    p.append(tiles(48, 199, REPLY_REAL, "right"))
    p.append(text(0, 239, "truth", 11, "muted"))
    p.append(tiles(48, 225, REAL_NEXT))
    p.append(text(0, 270, "Only the first amino acid is scored;", 11.5, "muted"))
    p.append(text(0, 286, "5 cuts in each protein’s domain.", 11.5, "muted"))
    p.append(text(0, 308, "Real reply; <tspan style=\"font-style:italic\">Kibdelosporangium persicum</tspan>,", 10.5, "faint"))
    p.append(text(0, 322, "family PF21962 (975 proteins).", 10.5, "faint"))
    return "".join(p)


def panel3():
    p = [badge(0, 16, "3", "Again, on a shuffled copy")]
    p.append(text(0, 46, "Same letters, random order:", 11.5, "muted"))
    p.append(text(0, 69, "real", 11, "muted"))
    p.append(tiles(48, 55, REAL_TAIL))
    p.append(text(0, 95, "shuffled", 11, "muted"))
    p.append(tiles(48, 81, SHUF_TAIL))
    p.append(arrow(136, 107, 121))
    p.append(text(0, 141, "reply", 11, "muted"))
    p.append(text(36, 142, "\u2717", 13, "bad", 700, "middle"))
    p.append(tiles(48, 127, REPLY_SHUF, "wrong"))
    p.append(text(0, 167, "truth", 11, "muted"))
    p.append(tiles(48, 153, SHUF_NEXT))
    p.append(rect(0, 190, 224, 96, "accent-soft", "accent", sw=1.25))
    p.append(text(112, 214, "gap = real − shuffled accuracy", 12.5, "ink", 650, "middle"))
    p.append(text(112, 236, "Opus 5.5 over 3,340 cuts:", 11.5, "body", 400, "middle"))
    p.append(text(112, 254, "34.9% − 7.5% = 27.4 points", 12.5, "accent-text", 700, "middle"))
    p.append(text(112, 274, "A guess at random scores 5%.", 11, "muted", 400, "middle"))
    p.append(text(0, 308, "From training data, the gap should", 11.5, "muted"))
    p.append(text(0, 322, "grow with the family’s size.", 11.5, "muted"))
    return "".join(p)


DEFS = (f'<defs><marker id="bp-arr" viewBox="0 0 10 10" refX="8.5" refY="5" markerWidth="7" markerHeight="7" '
        f'orient="auto-start-reverse"><path d="M1,1.5 L8.5,5 L1,8.5" style="fill:none;stroke:{v("faint")};'
        f'stroke-width:1.6;stroke-linecap:round;stroke-linejoin:round"/></marker></defs>')
LABEL = ("How the test works: ask for the next amino acid of proteins nobody studies, then of the same "
         "proteins shuffled; the gap in accuracy is what the model learned from sequence data.")


def svg(w, h, body, title_lines):
    head = "".join(text(24, 38 + 24 * i, t, 18, "ink", 700, "start", SANS) for i, t in enumerate(title_lines))
    return (f'<svg xmlns="http://www.w3.org/2000/svg" class="bp-figure" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
            f'role="img" aria-label="{LABEL}" focusable="false">{DEFS}'
            f'<rect x="0.5" y="0.5" width="{w - 1}" height="{h - 1}" rx="16" style="fill:{v("panel")};stroke:{v("line")};stroke-width:1"/>'
            f"{head}{body}</svg>\n")


panels = [panel1(), panel2(), panel3()]
rule = f'style="stroke:{v("line")};stroke-width:1"'
wide = "".join(f'<g transform="translate({24 + i * 248},66)">{p}</g>' for i, p in enumerate(panels))
wide += "".join(f'<line x1="{260 + i * 248}" y1="66" x2="{260 + i * 248}" y2="388" {rule}/>' for i in range(2))
narrow = "".join(f'<g transform="translate(24,{118 + i * 356})">{p}</g>' for i, p in enumerate(panels))
narrow += "".join(f'<line x1="24" y1="{118 + i * 356 - 22}" x2="248" y2="{118 + i * 356 - 22}" {rule}/>' for i in (1, 2))
Path("figures").mkdir(exist_ok=True)
Path("figures/fig0_method.svg").write_text(svg(760, 412, wide, ["Real protein or shuffled copy: the gap is what training taught"]))
Path("figures/fig0_method_mobile.svg").write_text(svg(272, 1168, narrow, ["Real protein or shuffled", "copy: the gap is what", "training taught"]))
print("figures/fig0_method.svg, figures/fig0_method_mobile.svg")
