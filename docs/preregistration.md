# Confirmation run: pre-registration

Written 2026-11-16, before any confirmation data was drawn or any model was called.
The git commit that adds this file is the timestamp. The pilot
(`data/sequences_pilot.tsv`, `results/probe/`) motivated it; none of its families
are reused.

## Question

Was Claude Opus 5.5 (and Opus 5) trained on protein sequence data at a scale that
reaches obscure bacterial and archaeal families no paper describes?

## Sample

- **Families:** Pfam families named "unknown function" (DUFs) with no PDB
  structure, at most 3 Europe PMC hits for the accession, at least 300 proteins,
  not in the pilot. 200 families, 50 from each size bin of proteins in the family:
  300-1k, 1k-3k, 3k-10k, >= 10k (fewer if a bin runs out).
- **Sequences:** per family, 4 unreviewed (TrEMBL) entries from Bacteria or
  Archaea, length 120-800, only the 20 standard residues, not low-complexity, and
  a UniRef90 singleton (no other sequence at >= 90% identity in UniProt).
- **Cuts:** 5 per sequence, uniformly at random inside the longest Pfam domain
  span (>= 40 residues), at least 30 residues into the sequence.
- Seeds fixed in code. Families yielding fewer than 4 sequences are dropped.

## Probes and models

The pilot's prompt, parser and settings, unchanged (`src/bioseq_mem/probe.py`,
`src/bioseq_mem/providers.py` at this commit). Each cut asked on the real sequence
and on the same sequence shuffled. Models: claude-opus-5-5, claude-opus-5
(primary); claude-opus-4-6, gpt-5.5, gemini-4-carbon (comparison).

## Outcomes and tests (per model; refusals, errors and unparseable replies excluded)

1. **H1, order signal.** D = residue-1 accuracy on real minus on shuffled cuts.
   Test: one-sided family-level sign-flip permutation test on per-family D
   (10,000 permutations); 95% CI by family bootstrap.
   Pass: p < 0.001 and CI lower bound > 10 percentage points.
2. **H2, dose-response.** Spearman correlation between per-family D and
   log10(proteins in family), one-sided, p < 0.01.
3. **H3, verbatim recall.** Count of real cuts with all 10 residues exactly right
   (none present in the prefix), against the count expected if residues were
   independent at the model's own per-position accuracy: sum over cuts of the
   product of position-wise accuracies, estimated per family. Pass: observed
   >= 10 and observed / expected > 100; also 0 such runs on shuffled cuts.

**Verdict for a model:** H1 and H3 pass = trained on these sequences (verbatim)
and on their regularities; H1 and H2 pass = systematic sequence training. All
three are reported whatever they show, for every model listed.

## Exploratory (not part of the verdict)

Accuracy by UniParc first-seen year; exact-run lengths; refusal rates.

## Amendment 1 (2026-11-16, before any confirmation model call)

- Add gpt-6-astra (reasoning effort low, its lowest; $10/$50 per 1M tokens) as a
  comparison model. It first runs on the pilot set; it joins the confirmation run
  only if its pilot real-minus-shuffled gain on DUF families exceeds 5 points.
- Batched UniRef90 singleton lookup (same criterion, fewer API calls).
- Protein language models (ProGen2, ESM-2), run locally on the same probes, are a
  reference level only; they do not change any test or verdict.

## Amendment 2 (2026-11-16, before any confirmation model call)

- gpt-6-astra joins the confirmation run at the user's request, without the pilot gate.
- Open-weight models on Fireworks (DeepSeek V4.1 Flash, Qwen3.8 2.4T, Kimi K3.1 Pro,
  GLM-5.3, MiniMax M3, gpt-oss-120b, Nemotron 3 Ultra) run on the pilot set as
  exploratory comparisons; thinking off where the model allows, else lowest effort.
- Protein language model reference dropped.

---

# DNA confirmation run: pre-registration

Written 2026-11-16 after the DNA pilot (`results/dna_summary.tsv`) and before any
confirmation window was drawn or any model was called. The pilot suggested a signal
for Claude Opus 5.5 in human and mouse coding sequence only.

## Sample

New windows, none overlapping a pilot transcript or locus, 220 nt each (200-nt
prefix, 20 nt to predict), from hg38/mm39 (.2bit, RepeatMasker soft-masking) and
GENCODE basic v50/vM39 canonical transcripts, sampled as in `scripts/12_dna_data.py`
with a new seed. 500 windows per cell:

- human and mouse: CDS, intron, promoter (1 kb upstream of the TSS), intergenic
  (no gene within 5 kb); non-coding windows contain no repeat and no N;
- bacterial/archaeal CDS: ENA coding sequences of proteins in
  `data/confirm_sequences.tsv` not used in the pilot (as many as available, up to 500).

## Probes and models

The pilot DNA prompt and parser, unchanged; each window asked on the real sequence
and on a dinucleotide-preserving shuffle of the same 220 nt. Models: claude-opus-5-5
(primary), claude-opus-5, gpt-5.5. (Gemini 4 Carbon blocks most DNA prompts.)

## Tests (per model and cell; errors, refusals, unparseable replies excluded)

1. **D1, order signal.** D = first-nucleotide accuracy on real minus shuffled windows.
   One-sided window-level sign-flip permutation test (10,000), window bootstrap 95% CI.
   Pass: p < 0.001 and CI lower bound > 5 percentage points.
2. **D2, verbatim recall.** Real windows continued with all 20 nt exactly right, the 20 nt
   not present in the prefix. Chance is about 1e-12 per window. Pass: >= 3 such
   windows and 0 on shuffled windows.

Predictions: Opus 5.5 passes D1 and D2 in human and mouse CDS, and fails both in
introns, promoters, intergenic DNA and bacterial CDS. Every cell is reported.
