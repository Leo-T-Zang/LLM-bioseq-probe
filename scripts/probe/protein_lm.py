"""Protein language models on the main set, as a reference -> results/plm_confirm/<model>__<mode>.jsonl.

Free on CPU (slow for ESM C and ESM3); scripts/probe/protein_lm_modal.py runs it on a GPU.

Same questions as the LLMs (data/confirm_sequences.tsv: every cut, real and shuffled), scored on the
residue right after the cut: the most likely of the 20 standard residues.

  prefix  only the residues before the cut: the same information the LLMs get. ESM-2 and ESM C are
          masked (bidirectional) models, so they see the prefix, a mask and the end-of-sequence token;
          ProtGPT2 is autoregressive and predicts the next token (a residue's probability = the summed
          probability of the BPE tokens that start with it).
  full    the whole sequence with that one residue masked (ESM-2 and ESM C only): more information
          than the LLMs get, an upper bound for these models.

Masking is done on token ids: the ESM SDK's "_" maps to a non-mask token in ESM C.
Sanity checks (full context, every residue masked in turn): ESM-2 650M ubiquitin 80%, hemoglobin
beta 62%, lysozyme 82%; ESM C 600M ubiquitin 100%, hemoglobin beta 80%; ESM3-open 1.4B ubiquitin 100%, hemoglobin beta 82%. ProtGPT2, next
residue from the left only: ubiquitin 65%, hemoglobin beta 37%.

    .venv-plm/bin/python scripts/probe/protein_lm.py esm2-650M:prefix,esm2-650M:full,esmc-600m:prefix,esmc-600m:full,esm3-open:prefix,esm3-open:full,protgpt2:prefix
"""

import json
import random
import sys
import zlib
from pathlib import Path

import pandas as pd
import torch

AA = "ACDEFGHIKLMNPQRSTVWY"
OUT = Path("results/plm_confirm")
torch.set_grad_enabled(False)
torch.set_num_threads(16)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def shuffled(seq: str, accession: str) -> str:  # same as bioseq_probe.probe.shuffled
    residues = list(seq)
    random.Random(zlib.crc32(accession.encode())).shuffle(residues)
    return "".join(residues)


def esm2():
    from transformers import AutoTokenizer, EsmForMaskedLM

    name = "facebook/esm2_t33_650M_UR50D"
    model, tok = EsmForMaskedLM.from_pretrained(name).to(DEVICE).eval(), AutoTokenizer.from_pretrained(name)
    ids = tok.convert_tokens_to_ids(list(AA))

    def logp(seq: str, cut: int, mode: str):
        x = tok(seq[:cut] + "A" if mode == "prefix" else seq, return_tensors="pt").input_ids
        x[0, cut + 1] = tok.mask_token_id  # +1 for <cls>
        return torch.log_softmax(model(input_ids=x.to(DEVICE)).logits[0, cut + 1, ids].float().cpu(), -1)

    return logp


def esmc():
    from esm.models.esmc import ESMC
    from esm.sdk.api import ESMProtein
    from esm.tokenization import EsmSequenceTokenizer

    model, tok = ESMC.from_pretrained("esmc_600m").to(DEVICE).eval(), EsmSequenceTokenizer()
    ids = tok.convert_tokens_to_ids(list(AA))

    def logp(seq: str, cut: int, mode: str):
        x = model.encode(ESMProtein(sequence=seq[:cut] + "A" if mode == "prefix" else seq)).sequence.clone()
        x[cut + 1] = tok.mask_token_id  # +1 for BOS
        return torch.log_softmax(model(sequence_tokens=x.unsqueeze(0).to(DEVICE)).sequence_logits[0, cut + 1, ids].float().cpu(), -1)

    return logp


def esm3():
    from esm.models.esm3 import ESM3
    from esm.sdk.api import ESMProtein
    from esm.tokenization import EsmSequenceTokenizer

    model, tok = ESM3.from_pretrained("esm3-open").to(DEVICE).float().eval(), EsmSequenceTokenizer()
    ids = tok.convert_tokens_to_ids(list(AA))

    def logp(seq: str, cut: int, mode: str):
        x = model.encode(ESMProtein(sequence=seq[:cut] + "A" if mode == "prefix" else seq)).sequence.clone()
        x[cut + 1] = tok.mask_token_id  # +1 for BOS
        return torch.log_softmax(model(sequence_tokens=x.unsqueeze(0).to(DEVICE)).sequence_logits[0, cut + 1, ids].float().cpu(), -1)

    return logp


def protgpt2():
    from transformers import AutoTokenizer, GPT2LMHeadModel

    model, tok = GPT2LMHeadModel.from_pretrained("nferruz/ProtGPT2").to(DEVICE).eval(), AutoTokenizer.from_pretrained("nferruz/ProtGPT2")
    vocab = tok.convert_ids_to_tokens(list(range(len(tok))))
    first = torch.tensor([AA.find(t[0]) if t and t[0] in AA else -1 for t in vocab])

    def logp(seq: str, cut: int, mode: str):
        prefix = seq[:cut]  # ProtGPT2 was trained on FASTA-like text, 60 residues a line
        text = "<|endoftext|>" + "\n".join(prefix[i : i + 60] for i in range(0, cut, 60)) + ("\n" if cut % 60 == 0 else "")
        p = torch.softmax(model(tok(text, return_tensors="pt").input_ids[:, -1024:].to(DEVICE)).logits[0, -1].float().cpu(), -1)
        mass = torch.zeros(len(AA)).index_add_(0, first.clamp(min=0), p * (first >= 0))
        return torch.log(mass / mass.sum() + 1e-12)

    return logp


LOADERS = {"esm2-650M": esm2, "esmc-600m": esmc, "esm3-open": esm3, "protgpt2": protgpt2}



def main(spec: str):
    seqs = pd.read_csv("data/confirm_sequences.tsv", sep="\t")
    OUT.mkdir(parents=True, exist_ok=True)
    runs = [r.split(":") for r in spec.split(",")]
    for mname in dict.fromkeys(m for m, _ in runs):
        logp_fn = LOADERS[mname]()
        for mode in [md for m, md in runs if m == mname]:
            path = OUT / f"{mname}__{mode}.jsonl"
            done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
            with path.open("a") as f:
                for r in seqs.itertuples():
                    for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession))):
                        for c in map(int, r.cuts.split(",")):
                            if (r.accession, c, cond) in done:
                                continue
                            lp = logp_fn(s, c, mode)
                            pred = AA[int(lp.argmax())]
                            f.write(json.dumps({"model": mname, "mode": mode, "accession": r.accession, "pfam": r.pfam,
                                                "condition": cond, "cut": c, "target": s[c : c + 10], "pred": pred,
                                                "first_correct": pred == s[c], "logp": float(lp[AA.index(s[c])])}) + "\n")
                    f.flush()
            print(mname, mode, "done", flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
