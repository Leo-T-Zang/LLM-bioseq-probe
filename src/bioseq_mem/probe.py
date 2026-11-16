"""The next-residue probe: prompt construction and scoring."""

import random
import re
import zlib

K = 10  # residues asked for; residue 1 is the clean next-token measure
STANDARD = "ACDEFGHIKLMNPQRSTVWY"


def fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def prompt(prefix: str) -> str:
    return (
        "Below is the beginning of a protein sequence in one-letter amino acid codes. "
        f"Predict how it continues. Reply with exactly the next {K} amino acids as one-letter "
        "codes on one line, and nothing else.\n\n"
        f">sequence\n{fasta(prefix)}"
    )


def parse(text: str) -> str:
    letters = re.sub(r"[^A-Za-z]", "", text.strip().splitlines()[-1] if text.strip() else "").upper()
    return letters if set(letters) <= set(STANDARD) else ""


def score(pred: str, target: str) -> dict:
    run = 0
    while run < min(len(pred), len(target)) and pred[run] == target[run]:
        run += 1
    return {
        "valid": bool(pred),
        "first_correct": bool(pred) and pred[0] == target[0],
        "exact_run": run,
        "identity": sum(a == b for a, b in zip(pred, target)) / len(target),
    }


def shuffled(seq: str, accession: str) -> str:
    """The same residues in a random order, seeded by accession: composition kept, order lost."""
    residues = list(seq)
    random.Random(zlib.crc32(accession.encode())).shuffle(residues)
    return "".join(residues)
