"""Protein language models on the same probes, as a reference level -> <out>/<model>.jsonl. Free, CPU.

Only architectures built into transformers (no remote code):
  ESM-2 650M (masked): prefix + <mask>, nothing after it; logits at the mask.
  ProtGPT2 (autoregressive, BPE): next-token distribution after the prefix; a residue's
    probability = the summed probability of tokens that start with it.
Prediction = the most likely of the 20 standard residues; scored on residue 1 only.

    .venv-plm/bin/python scripts/10_plm.py data/sequences_pilot.tsv results/plm
"""

import json
import random
import sys
import zlib
from pathlib import Path

import pandas as pd
import torch
from transformers import AutoTokenizer, EsmForMaskedLM, GPT2LMHeadModel

AA = "ACDEFGHIKLMNPQRSTVWY"
torch.set_grad_enabled(False)


def shuffled(seq: str, accession: str) -> str:  # same as bioseq_probe.probe.shuffled
    residues = list(seq)
    random.Random(zlib.crc32(accession.encode())).shuffle(residues)
    return "".join(residues)


def esm2(name):
    model = EsmForMaskedLM.from_pretrained(name).eval()
    tok = AutoTokenizer.from_pretrained(name)
    ids = tok.convert_tokens_to_ids(list(AA))

    def predict(prefix):
        enc = tok(prefix + tok.mask_token, return_tensors="pt")
        pos = int((enc.input_ids[0] == tok.mask_token_id).nonzero()[0])
        return torch.log_softmax(model(**enc).logits[0, pos, ids], -1)

    return predict


def protgpt2(name):
    model = GPT2LMHeadModel.from_pretrained(name).eval()
    tok = AutoTokenizer.from_pretrained(name)
    vocab = tok.convert_ids_to_tokens(list(range(len(tok))))
    first = torch.tensor([AA.find(t[0]) if t and t[0] in AA else -1 for t in vocab])

    def predict(prefix):
        # ProtGPT2 was trained on FASTA-like text with newlines every 60 residues.
        text = "<|endoftext|>" + "\n".join(prefix[i : i + 60] for i in range(0, len(prefix), 60))
        if len(prefix) % 60 == 0:
            text += "\n"
        x = tok(text, return_tensors="pt").input_ids[:, -1024:]
        p = torch.softmax(model(x).logits[0, -1], -1)
        mass = torch.zeros(len(AA)).index_add_(0, first.clamp(min=0), p * (first >= 0))
        return torch.log(mass / mass.sum() + 1e-12)

    return predict


MODELS = {"esm2-650M": lambda: esm2("facebook/esm2_t33_650M_UR50D"), "protgpt2": lambda: protgpt2("nferruz/ProtGPT2")}

seqs = pd.read_csv(sys.argv[1], sep="\t")
out_dir = Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
torch.set_num_threads(16)
for mname, load in MODELS.items():
    path = out_dir / f"{mname}.jsonl"
    done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
    predict = load()
    with path.open("a") as f:
        for r in seqs.itertuples():
            for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession))):
                for c in map(int, r.cuts.split(",")):
                    if (r.accession, c, cond) in done:
                        continue
                    logp = predict(s[:c])
                    pred, t = AA[int(logp.argmax())], s[c]
                    f.write(json.dumps({"model": mname, "accession": r.accession, "pfam": r.pfam, "condition": cond, "cut": c,
                                        "target": s[c : c + 10], "pred": pred, "first_correct": pred == t, "logp": float(logp[AA.index(t)])}) + "\n")
            f.flush()
    print(mname, "done", flush=True)
