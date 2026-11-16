"""The nucleotide probe: prompt, parsing, scoring and a dinucleotide-preserving shuffle."""

import random
import re
import zlib
from collections import defaultdict

from bioseq_mem.probe import GUESS

K = 20  # nucleotides asked for; 4^-20 ~ 1e-12, so an exact run is not luck
PREFIX = 200


def prompt(prefix: str, guess: bool = False) -> str:
    """guess=True adds bioseq_mem.probe.GUESS, for models that decline the plain prompt."""
    wrapped = "\n".join(prefix[i : i + 60] for i in range(0, len(prefix), 60))
    return (
        "Below is the beginning of a DNA sequence (5' to 3'). Predict how it continues. "
        f"{GUESS if guess else ''}Reply with exactly the next {K} nucleotides (A, C, G, T) on one line, and nothing else.\n\n"
        f">sequence\n{wrapped}"
    )


def parse(text: str) -> str:
    text = re.sub(r"<[^>]*>", "", text).strip()
    letters = re.sub(r"[^A-Za-z]", "", text.splitlines()[-1] if text else "").upper()
    return letters if letters and set(letters) <= set("ACGT") else ""


def dinucleotide_shuffle(seq: str, key: str) -> str:
    """Random sequence with exactly the same dinucleotide counts and first/last base
    (Altschul-Erickson, via a random Eulerian path), seeded by key."""
    rng = random.Random(zlib.crc32(key.encode()))
    edges = defaultdict(list)
    for a, b in zip(seq, seq[1:]):
        edges[a].append(b)
    last = seq[-1]
    while True:  # a random last exit per vertex, forming a tree into `last`
        exit_ = {v: rng.choice(out) for v, out in edges.items() if v != last}
        ok = True
        for v in exit_:
            seen, u = set(), v
            while u != last and ok:
                if u in seen:
                    ok = False
                seen.add(u)
                u = exit_.get(u, last)
            if not ok:
                break
        if ok:
            break
    order = {}
    for v, out in edges.items():
        out = out.copy()
        if v in exit_:
            out.remove(exit_[v])
        rng.shuffle(out)
        order[v] = out + ([exit_[v]] if v in exit_ else [])
    pos = dict.fromkeys(order, 0)
    walk, u = [seq[0]], seq[0]
    for _ in range(len(seq) - 1):
        u = order[u][pos[u]]
        pos[walk[-1]] += 1
        walk.append(u)
    return "".join(walk)
