"""The next-residue probe: prompt construction and scoring."""

import random
import re
import zlib

K = 10  # residues asked for; residue 1 is the clean next-token measure
STANDARD = "ACDEFGHIKLMNPQRSTVWY"


def fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


GUESS = "Give your best guess even if you are uncertain. Do not explain or refuse. "


def prompt(prefix: str, guess: bool = False) -> str:
    """guess=True adds GUESS, for models that decline the plain prompt."""
    return (
        "Below is the beginning of a protein sequence in one-letter amino acid codes. "
        f"Predict how it continues. {GUESS if guess else ''}Reply with exactly the next {K} amino acids as one-letter "
        "codes on one line, and nothing else.\n\n"
        f">sequence\n{fasta(prefix)}"
    )


ANCHOR = 15  # a reply that starts with at least this much of the prefix's end is restating it
LABEL = re.compile(r"\s*(?::|=|->|→)\s*")  # "Answer: ...", "I'll predict: ...", "... → ..."
REASONING = re.compile(r"<(think|thinking|reasoning)>.*?(?:</\1>|$)", re.S | re.I)
WORD_LEN = 8  # capitals before prose on the same line count as a sequence only from this length (not FASTA, WARNING)
# stop reasons for a reply that hit the output-token limit, across providers
CUT_OFF = {"length", "max_tokens", "MAX_TOKENS", "FinishReason.MAX_TOKENS", "max_output_tokens"}


def cut_off(record: dict) -> bool:
    """True if the reply hit the output-token limit before it finished."""
    return str(record.get("stop")) in CUT_OFF


def _read_line(line: str, alphabet: str, labelled: bool = False) -> tuple[str, bool]:
    """(the sequence a line gives, whether the whole line is sequence). Sequence letters are
    capitals of the alphabet; spaces and numbers between them are ignored. labelled: the text
    follows a label such as "Answer:", so it is meant to be a sequence."""
    toks = [t.strip("_") for t in line.split() if not re.fullmatch(r"[\d.:()#•-]+", t)]  # numbering or bullets: "1.", "61", "-"
    if not toks:
        return "", False
    bare = toks[:-1] + [toks[-1].rstrip(".;!")]  # a comma elsewhere marks a list ("E, D, K, H"), not a sequence
    if all(re.fullmatch(f"[{alphabet}]+", t) for t in bare):
        return "".join(bare), True
    if len(toks) == 1 and len(bare[0]) >= WORD_LEN and bare[0].isalpha() and set(bare[0].upper()) <= set(alphabet) \
            and 2 * sum(c.isupper() for c in bare[0]) >= len(bare[0]):
        return bare[0].upper(), True  # at least half capitals: LYPIGaalga, stratRCEIGA; not Sacrificing
    parts = LABEL.split(line)
    if len(parts) > 1:
        seq, whole = _read_line(parts[-1], alphabet, labelled=True)
        if seq:
            return seq, whole
    prose = any(c.islower() for c in line)  # in prose, only a run of at least WORD_LEN capitals is a sequence
    j = 0
    while re.fullmatch(f"[{alphabet}]+", bare[j]):
        j += 1
    seq, nxt = "".join(bare[:j]), toks[j]
    run = re.match(f"[{alphabet}]*", nxt).group()
    stop = nxt[len(run) : len(run) + 1]
    if len(nxt) <= 40 and re.fullmatch(r"[A-Za-z0-9]+(?:_[A-Za-z0-9]+)+", nxt.strip(".,;:!?()[]{}<>")):  # an identifier: VALIDATION_FAILED
        pass
    elif stop.islower():  # a word follows, perhaps glued on ("LSKLTLoading"): its capital is not a residue
        seq += run[:-1]
    elif stop in ":=,;" or nxt[len(run) :].startswith("->"):  # a label or a list follows ("ANSWER:", "E, D, K")
        seq += run
    elif stop.isascii() and stop.isupper():  # a capital that is not a residue: O, U, X, Z, B, J
        if not prose and (labelled or (len(run) >= 3 and sum(c.isalpha() for c in nxt) >= 10)):
            return seq + run, False  # the sequence ends there: "Answer: KLMNOPQRST", "...KNLPKFINALIZED"
        seq += run  # otherwise it counts only if long: not NOTE, PROTEIN, SEQUENCE
    elif nxt[len(run) + 1 : len(run) + 2].islower() and len(run) < WORD_LEN:
        pass  # a word with an apostrophe or a hyphen: I’m, C-terminal
    elif not prose:  # any other symbol ends the sequence: KLMN?PQRST, a letter of another script
        return seq + run, False
    else:
        seq += run  # in prose: ESM-2 and WNT5A are terms, a long run before a symbol is a sequence
    return (seq, False) if len(seq) >= WORD_LEN else ("", False)  # capitals before prose: not FASTA, WARNING


def _closed(line: str) -> str:
    """The line without an unfinished quotation at its end."""
    if line.count('"') % 2:
        line = line[: line.rfind('"')]
    if line.count("“") > line.count("”"):
        line = line[: line.rfind("“")]
    opens = [m.start() for m in re.finditer(r"(?<![A-Za-z0-9])'(?=[A-Za-z0-9])", line)]  # 'ACDEF, not I'll or proteins'
    if opens and not re.search(r"(?<=\S)'(?![A-Za-z0-9])", line[opens[-1] + 1 :]):  # closed: 'ACDEF.' too
        line = line[: opens[-1]]
    return line


def _items(text: str, alphabet: str) -> list[tuple[str, bool]]:
    """(sequence, opens the reply) for every sequence a reply states, in order; see sequences_in."""
    text = re.sub(r"<[^>]*>", "\n", REASONING.sub("\n", text or ""))  # drop reasoning, then markup
    items, block, opens, first = [], "", False, True

    def close():
        nonlocal block, opens
        if block:
            items.append((block, opens))
        block, opens = "", False

    for raw in text.splitlines():
        line = re.sub(r"[*`~'\"“”]", "", _closed(raw.strip())).strip()
        if raw.strip().startswith(("```", ">")) or not re.search("[A-Za-z]", line):
            close()  # a code fence, a FASTA header, or a line without letters ends a sequence
            continue
        seq, whole = _read_line(line, alphabet)
        parts = LABEL.split(line)
        label = " ".join(parts[:-1]).split() if len(parts) > 1 else []
        if len(parts) > 1:
            close()  # "Correction: ..." or "Contains: ..." starts afresh; it does not continue the lines above
        opening = first and bool(seq) and len(label) <= 3  # "Answer: ..." opens; a sentence that quotes does not
        if seq and whole:
            if not block:
                opens = opening
            block += seq
        else:
            close()
            if seq:
                items.append((seq, opening))
        first = False
    close()
    return items


def sequences_in(text: str, alphabet: str = STANDARD) -> list[str]:
    """Every sequence a reply states, in order. A sequence is:
    - lines holding only capital sequence letters (spaces, numbers, markdown ignored), joined
      when they wrap, also after a label ("Answer: ACDEF" then "GHIKL");
    - the capitals after a label, as in "I'll predict: DLKVLKEGVD";
    - capitals that open a line and end at a symbol that is not a residue (KLMN?PQRST, a letter
      of another script), or before prose on the same line if at least WORD_LEN long;
    - a single word that is mostly capitals (LYPIGaalga).
    Words are not residues: "I'll predict" is not ILLPREDICT, "Deep" is not DEEP. Text inside
    <think>/<reasoning> tags and an unfinished quotation at the end of a line are skipped."""
    return [s for s, _ in _items(text, alphabet)]


def restated(seq: str, prefix: str) -> str | None:
    """If seq starts by restating the end of the prefix (at least ANCHOR letters of it, or all of
    a shorter prefix), what follows the restatement ("" if nothing); else None."""
    for n in range(min(len(seq), len(prefix)), min(ANCHOR, len(prefix)) - 1, -1):
        if n and seq[:n] == prefix[-n:]:
            return seq[n:]
    return None


def parse(text: str, prefix: str = "", alphabet: str = STANDARD, prefer_len: int = 1, cut: bool = False) -> str:
    """The model's answer; "" if the reply gives none.
    A finished reply: the last sequence it states, since a final answer or a correction comes
    last (with prefer_len, the last one at least that long if there is one, else the last one).
    A reply cut off at the token limit (cut=True; see cut_off) has no final answer: only a
    sequence it opens with counts, because what follows is often a degenerate repeat or
    unfinished reasoning. With a prefix, a restatement of its end is removed (see restated);
    a sequence that is only a restatement gives nothing."""
    items = _items(text, alphabet)
    seqs = [s for s, opens in items[:1] if opens] if cut else [s for s, _ in items]
    if prefix:
        seqs = [c for s in seqs if (c := s if (r := restated(s, prefix)) is None else r)]
    long = [s for s in seqs if len(s) >= prefer_len]
    return (long or seqs or [""])[-1]


def parse_raw(text: str) -> str:
    """The first K residues of a raw continuation: sequence lines up to a '>' header, markup, a
    blank line or the first line of prose. A word's capital ("Loading", "This") is not a residue."""
    out = ""
    for line in re.split(r"[<>]|\n\s*\n", text or "", maxsplit=1)[0].split("\n"):  # markup, a header or a blank line ends it
        if not line.strip():
            continue
        seq, whole = _read_line(line.strip(), STANDARD)
        out += seq
        if not whole or len(out) >= K:
            break
    return out[:K]


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


def shown_prefix(seq: str, accession: str, condition: str, cut: int) -> str:
    """The prefix a model was shown: the real or the shuffled sequence up to the cut."""
    return (seq if condition == "real" else shuffled(seq, accession))[:cut]


def prompt_echo(prefix: str) -> str:
    """Copy-then-continue: the model writes the given part, then keeps going, so the
    residues it predicts follow its own copy of the sequence, as in plain next-token prediction."""
    return (
        "Below is the beginning of a protein sequence in one-letter amino acid codes. "
        f"Write out the sequence: copy the given part exactly, then continue it with the next {K} amino acids. "
        "Reply with the sequence only, on one line, with no spaces or other text.\n\n"
        f">sequence\n{fasta(prefix)}"
    )


def parse_echo(text: str, prefix: str) -> str:
    """The residues the model wrote after restating the prefix's end (see restated), in the
    last stated sequence that restates it; "" if none does."""
    found = echo(text, prefix)
    return found or ""


def echo(text: str, prefix: str) -> str | None:
    """What follows the reply's restatement of the prefix's end ("" if nothing), or None if the
    reply does not restate it."""
    for s in reversed(sequences_in(text)):
        if (r := restated(s, prefix)) is not None:
            return r
    return None
