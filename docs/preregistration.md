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
