"""How replies are read. Run with: uv run --with pytest pytest tests

Replies marked with a run, model and record are verbatim from results/outputs/ (some shortened);
the others are made up to pin down one rule each."""

from types import SimpleNamespace

from bioseq_probe import dna
from bioseq_probe.probe import cut_off, echo, parse, parse_echo, parse_raw
from bioseq_probe.providers import responses_stop

PREFIX = "MSTNKLLEAVRQGDLAAVKRLLDEGADPNARDNLGRTPLHLAAQAGHLEIVKLLLEHGAD"  # made up


def test_words_are_not_residues():
    # probe, Qwen 3.8, B4SB17:289 (end of a 2,374-character reply)
    assert parse('A reasonable prediction might be "DLKVLKEGVD" or "DLEALKAGID".\n\nI\'ll predict: DLKVLKEGVD') == "DLKVLKEGVD"
    # confirm, Claude Opus 5, A0AAJ5NE02:196
    assert parse("ASRDFEGNPI\n\nSacrificing\n\nDeep") == "ASRDFEGNPI"
    assert parse("The answer is: GQPNGRLPRS") == "GQPNGRLPRS"
    assert parse("ФКНЫЛЕТЕНН → FKHYLETENN") == "FKHYLETENN"
    assert parse("I cannot predict this sequence.") == ""
    # probe, gpt-oss-120b: a refusal with a curly apostrophe
    assert parse("I’m sorry, but I can’t provide that.") == "" and parse_raw("I’m sorry") == ""
    assert parse("C-terminal region") == ""
    # terms in an explanation do not replace the answer (confirm, Claude Opus 5; probe, Qwen 3.8)
    assert parse("VLAGAAIAFS\n\n- ESM-2 style predictions favour hydrophobic residues") == "VLAGAAIAFS"
    assert parse("NCSNPQPLQN\n\nLet me verify: WNT5A sequence is long") == "NCSNPQPLQN"
    # lists and identifiers (probe, Claude Opus 5, A0ABQ0EWJ1:114: an answer, then API documentation)
    assert parse("E, D, K, H.") == "" and parse("ACDEFGHIKL\nContains charged residues: E, D, K, H.") == "ACDEFGHIKL"
    assert parse("E , D , K , H.") == "" and parse("ACDEFGHIKL\nContains charged residues: E , D , K , H.") == "ACDEFGHIKL"
    assert parse("E; D; K; H.") == "" and parse("RSLGLTHGRA\nContains charged residues: E ; D ; K ; H.") == "RSLGLTHGRA"
    assert parse('RSLGLTHGRA\n\n{"code": "VALIDATION_FAILED"}') == "RSLGLTHGRA"
    assert parse("- ACDEFGHIKL") == "ACDEFGHIKL"  # a bullet is not part of the answer
    assert parse('RSLGLTHGRA\n\n```json\n{\n  "code": "VALIDATION_FAILED"\n}\n```') == "RSLGLTHGRA"
    assert parse("??????????") == ""


def test_the_last_answer_counts():
    # confirm, Claude Opus 5.5, A0A4R6QJY3:353
    assert parse("AQRLSQLARSL\n\nWait, need exactly 10: AQRLSQLARS") == "AQRLSQLARS"
    # confirm, Claude Opus 5.5, M0HYM5:274 (middle shortened)
    assert parse("MDRMLDELGVEAEIPKLTGLRLEN\n\nThe sequence is cut off ...\n\nThe first 10 are: MDRMLDELGV") == "MDRMLDELGV"
    # a reply that opens with prose quoting a family member, then answers (as Claude Fable 5.2 does)
    assert parse('This resembles the DUF2997 family proteins: "MKTIEITVSPDGQTRVETKG..."\n\nLEEALGKRTS') == "LEEALGKRTS"


def test_sequence_and_text_on_one_line():
    assert parse("ACDEFGHIKL Note: uncertain") == "ACDEFGHIKL"
    # recall_length, Gemini 3.5 Flash, A0A9B1I0J9:88 (shortened)
    assert parse("VESVKAKVPKIKVEALFNEQMVLQEAKALEELVKGKRERLK Standard", prefer_len=20) == "VESVKAKVPKIKVEALFNEQMVLQEAKALEELVKGKRERLK"
    assert parse("FASTA format:") == ""  # a capitalised word before prose is not a sequence
    for heading in ("ANSWER:", "NOTE: none", "PROTEIN SEQUENCE:", "SEQUENCE", "DNA sequence"):
        assert parse(heading) == "", heading


def test_formatting():
    # confirm, Claude Opus 4.8, G8RU09:264: a wrapped answer is one sequence
    assert parse("DRAKWDRLLI\nAQ") == "DRAKWDRLLIAQ"
    assert parse("Answer: ACDEF\nGHIKL") == "ACDEFGHIKL"
    assert parse("Answer: ACDEF\nCorrection: GHIKL") == "GHIKL"  # a new label starts afresh
    assert parse("**TKCVIVLSHC**") == "TKCVIVLSHC"
    assert parse("```\nLRLSSSLHGG\n```") == "LRLSSSLHGG"
    assert parse("<br>HLGLLQGSLL") == "HLGLLQGSLL"
    assert parse("QVAVGLLHRLAEAG\n\nA: QVAVGLLHRL") == "QVAVGLLHRL"
    # confirm, Claude Opus 4.5, A0A1H1A8Q4:225: mostly capitals
    assert parse("LYPIGaalga") == "LYPIGAALGA"
    # confirm, Claude Opus 5.5, A0A0B3RLQ6:143
    assert parse("stratRCEIGA") == "STRATRCEIGA"


def test_reasoning_and_unfinished_quotes_are_skipped():
    assert parse("<think>Just guess a plausible continuation.</think>LYLSLFEGLE") == "LYLSLFEGLE"
    assert parse("<think>\nACDEFGHIKL\n</think>") == ""
    assert parse('Unclosed "ACDEFGHIKL') == ""
    assert parse("Unclosed 'ACDEFGHIKL") == "" and parse("'ACDEFGHIKL'") == "ACDEFGHIKL"
    assert parse("The proteins' sequence: ACDEFGHIKL") == "ACDEFGHIKL"  # a possessive is not a quotation
    assert parse("'ACDEFGHIKL.'") == "ACDEFGHIKL" and parse("'ACDEFGHIKL?'") == "ACDEFGHIKL"


def test_an_invalid_symbol_ends_the_sequence():
    assert parse("KLMN?PQRST") == "KLMN"
    assert parse("Answer: KLMNOPQRST") == "KLMN"  # O is not one of the 20
    # dna, GPT-5.5, a human intergenic window: a letter of another script
    assert dna.parse("TCTGTTCCTTCAATCCTCאַט") == "TCTGTTCCTTCAATCCTC"


def test_restating_the_prefix():
    assert parse(PREFIX + "KVLAG", PREFIX) == "KVLAG"  # the whole prefix, then the continuation
    assert parse(PREFIX[-20:] + "\nKVLAG", PREFIX) == "KVLAG"  # its end, wrapped, then more
    assert parse(PREFIX, PREFIX) == ""  # a restatement only: no prediction
    assert parse("ACDEFGHIKL", "ACD") == "EFGHIKL"  # a prefix shorter than ANCHOR, restated whole
    # a guess that is a piece of the prefix, not its end, stays a guess
    # (probe, Qwen 3.8, A0ABX0IXB0:302: "I'll predict: KTKVFGVVISM", found inside the prefix)
    prefix = "MKV" + "KTKVFGVVISM" + "LLEAGRTPLHLAAQAGHLEIVKLLLEHGAD"
    assert parse("I'll predict: KTKVFGVVISM", prefix) == "KTKVFGVVISM"
    window = "ACGTTGCAACGTAGGCTAGCTAACGGATCCATGCAAGT" * 5  # made up
    assert dna.parse(window[-18:] + "GC", window) == "GC"  # GPT-5.5 restates the end of the window
    assert dna.parse(window[40:60], window) == window[40:60]  # a stretch from the middle: a guess


def test_replies_cut_off_at_the_token_limit():
    # confirm, Claude Opus 5: an answer, then unrelated text until the limit
    assert parse("LLSGPVTAQQ\n\nВыполнено\nПонравился ответ?\n\nAnswer: AAGGTTCCAA", cut=True) == "LLSGPVTAQQ"
    # probe, Kimi K3.1 Pro: reasoning that never reached an answer
    assert parse("I need predict next 10 aa. Could be NuMA1.\nWVDLQAGPEE... perhaps next", cut=True) == ""
    assert parse("```text\nACDEFGHIKL\n```\nand then", cut=True) == "ACDEFGHIKL"
    assert parse('Close to the DUF2997 proteins: "MKTIEITVSPDG" so', cut=True) == ""
    assert cut_off({"stop": "max_output_tokens"}) and not cut_off({"stop": "end_turn"})


def test_openai_responses_stop_reason():
    incomplete = SimpleNamespace(status="incomplete", incomplete_details=SimpleNamespace(reason="max_output_tokens"))
    assert responses_stop(incomplete) == "max_output_tokens" and cut_off({"stop": responses_stop(incomplete)})
    assert responses_stop(SimpleNamespace(status="completed", incomplete_details=None)) == "completed"


def test_dna():
    assert dna.parse("The next 20 nucleotides are: ACGTACGTAC") == "ACGTACGTAC"
    # dna_confirm, Claude Opus 5: the answer, then a lowercase fragment
    assert dna.parse("GAAATGGAAAATGGCTGA\n\nGcagacaacag") == "GAAATGGAAAATGGCTGA"
    assert dna.parse("<answer>GTGCTGCTGGCCACCGCCT</answer>\n\nI need exactly 20.\n\n<answer>GTGCTGCTGGCCACCGCCTG</answer>") == "GTGCTGCTGGCCACCGCCTG"


def test_raw_continuation():
    # raw, GLM-5.3, A0A495EE07:355: the "L" of "Loading" is not a residue
    assert parse_raw("LSKLT\nKQKQ\nLoading peptide info... please wait.\nPeptide:  323-") == "LSKLTKQKQ"
    assert parse_raw("LHSYF\n\nOf course. Here are a few ways") == "LHSYF"
    assert parse_raw("QGQGQGQGQGQGQG\n```\n\nThis is a") == "QGQGQGQGQG"
    assert parse_raw("This is a short protein") == ""
    assert parse_raw("\nMKVLAKVLAGRT") == "MKVLAKVLAG"
    # raw, MiniMax M3: residues, then the markup of a table
    assert parse_raw("GPGIPR</td></tr><tr><td>SEQ ID NO: 2</td><td>") == "GPGIPR"
    # raw, Gemini 2.5 Pro: a repeat cut off mid-word
    assert parse_raw("LSVLSVLSVLSVLSVLSV-t") == "LSVLSVLSVL"
    assert parse_raw("IGQYQGTRLADGKLQAGLGIQXDLAAVG and so on") == "IGQYQGTRLA"  # an X, then prose, on one line
    assert parse_raw("QGRLQGRLQGRLQGRLQGRL_QGRL*") == "QGRLQGRLQG"  # stray symbols far along a long repeat
    assert parse_raw("QGRL" * 12 + "_" + "QGRL" * 12) == "QGRLQGRLQG"  # (raw, Gemini 2.5 Pro: 2,666 letters, one _)
    # no_instruction, GPT-5.5: metadata, not residues
    assert parse_raw("FASTA format:\n\n```fasta\n>sequence\nAAVELGGRAGYALVSD") == ""


def test_echo():
    # echo, GPT-5.5, A0A8J7G199:68 (the prefix here is the reply's copy of it)
    reply = "MNEQFIKLGHGYGDIYELFELMKTNQTRIHNAFLLITEHQGKAVCSFAVALKPAGSSQFMPIYICLEGKNLPKFINALIZED\n"
    assert parse_echo(reply, reply[: reply.index("KNLPK")]).startswith("KNLPK")
    assert parse_echo("KVLAG", PREFIX) == ""
    # no_instruction, GPT-5.5, A0A286RE48:83: a copy of the prefix and nothing more
    assert echo(PREFIX, PREFIX) == "" and echo("KVLAG", PREFIX) is None
    # no_instruction, GPT-5.5, A0A078L199:31: a copy, then a list of residues that is not a continuation
    assert echo(f"Sequence: {PREFIX}\n- Contains charged residues: E, D, K, H.", PREFIX) == ""


def test_recall_takes_the_recited_sequence():
    # recall_length, Claude Opus 5.5, A0A4W2CF04:53 (sequence shortened)
    reply = ("This looks like the N-terminus of an arrestin (beta-arrestin), so the continuation below is "
             "reconstructed from memory of beta-arrestin-1 rather than derived from the fragment. Individual "
             "residues are likely to be wrong.\n\nYVTLTCAFRYGREDLDVLGLTFRKDLFVANVQSFPPAPEDKKPLTRLQERLIKK")
    assert parse(reply, PREFIX, prefer_len=20).startswith("YVTLTCAFRYGREDLDVLG")
    assert parse("ACD", prefer_len=20) == "ACD"  # a preference, not a threshold
