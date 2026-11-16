# Findings

Run through each provider's API in November 2026. Figures: `figures/fig0_method` to `fig5_family_size`.

## Claude Opus 5.5 was trained on protein sequence data

On 668 bacterial and archaeal proteins from 167 unknown-function families (no papers to speak
of, no relative at 90% identity), Claude Opus 5.5 gets the next amino acid right 34.9% of the
time, against 7.5% on the same residues shuffled and 5% for a blind guess: a gap of 27.4
points (95% CI 24.7 to 30.1, sign-flip permutation p < 0.0001).

- The gap appears in 95% of families (median 25 points); the 10 largest families give 14% of it.
- It grows with family size (Spearman 0.27, p = 0.0002): 22 points for families of 300-1k
  proteins, 37 for families over 10k. In the 128 families with zero Europe PMC hits it is
  26.8 points and still grows with size (rho 0.28, p = 0.001); it does not track paper counts
  (rho 0.05).
- It is not reasoning at test time: on the 55% of questions where Opus 5.5 wrote 15 tokens or
  fewer, it scores 35.6% vs 7.3%. GPT-6 Astra thinks for 200+ tokens a question and reaches
  about 10%.
- It is not the prompt: asking Opus 5.5 to copy the sequence and continue gives the same result.
- It is not recall of UniProt pages: entries on uniprot.org before 2023 and entries added from
  June 2025 score alike (37.2 vs 37.6 points, pilot).
- Exact recall of these obscure proteins is weak: 12 cuts continued 10/10 where 0.56 are
  expected from per-position accuracy; re-asked for 200 residues, Opus 5.5 rarely goes past 10.
  It learned the families' statistics more than each protein's letters.

## Within Claude, the signal is recent

Claude Opus 4.6: +0.9 points. Opus 4.8: +7.7 (5.6 to 9.9). Opus 5: +9.7 (7.4 to 12.2). Opus
5.5: +27.4. Opus 4.8 and Opus 5 also grow with family size (rho 0.31 and 0.29), so the data
was there from Opus 4.8 on and much stronger in Opus 5.5.

## No other model shows it

GPT-5.5 (+2.1), GPT-6 Astra (+2.9), Gemini 4 Carbon (+0.9) through the chat API; Gemini 2.5
Pro (+1.1) and seven open-weight models (DeepSeek V4.1 Flash, Qwen 3.8, Kimi K3.1 Pro, GLM-5.3,
MiniMax M3, gpt-oss-120b, Nemotron 3 Ultra: -0.0 to +1.9) as raw next-token continuation of
the FASTA text. Absence of a signal does not prove absence of training data.

## Famous proteins prove little

On well-studied vertebrate proteins most models do well (Gemini 4 Carbon +47.6 points, Opus
4.8 to 5.5 +44 to +52), and several recite hundreds of residues letter for letter (Opus 5
longest 325, Opus 5.5 305, Gemini 4 Carbon 201). These sequences fill papers and web pages.
GPT models are the exception: GPT-5.5's longest exact run on these proteins is 2 residues.

## DNA: coding sequence only

On 500 repeat-free windows per region, Opus 5.5 predicts human coding sequence (+30.5, 25
exact 20-nucleotide continuations) and mouse coding sequence (+28.5, 5 exact), and nothing
clear in promoters, introns or intergenic DNA (-0 to +5). The DNA of the bacterial genes
behind the proteins in the main set gives +6.0 (p = 0.026), against +27.4 as protein: the
training data held protein sequences, not the genomes they come from.

## Safety behaviour seen along the way

- Gemini 4 Carbon blocks about 75% of DNA prompts before generation (`block_reason OTHER`).
- Claude Sonnet 4.5 refused 97% of protein prefill requests.
- OpenAI rejected 23 GPT-6 requests as "flagged for possible biological" misuse.
- GPT-4o, DeepSeek V4.1, gpt-oss and others often declined to guess a sequence ("I can't
  reliably predict…") until asked for a best guess.

## What this does not show

Whether the sequences came from a deliberate sequence corpus or crawled database pages; in
which training phase they entered; whether any lab trained on protein-model outputs; that the
other models never saw sequences.

## Related work

[PG-LLM](https://www.proteingymllm.com/) (Arora, Chen, Du, Marks and Church, 2026) ranks protein
variants with general-purpose LLMs; Claude Opus 5 approached ESM-2, with accuracy that barely
changes with alignment depth. A [post by Atharv Jayprakash](https://x.com/atharvjay872/status/2121477118450340342)
found GPT-6 Astra ranking SOD1 variants about as well as ESM-2. Variant ranking can draw on
biology and the literature; the test here removes that route, and GPT-6 shows almost no
sequence knowledge on it.
