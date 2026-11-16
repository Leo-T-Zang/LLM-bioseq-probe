"""ESM-family protein language models on the same probes -> <out>/<model>__<mode>.jsonl. Free, CPU.

Models: ESM C 600M, ESM3 open (1.4B) via EvolutionaryScale's `esm` package; ESM-2 650M via transformers.
Modes, both scored on the residue right after the cut (argmax over the 20 standard residues):
  prefix  only the residues before the cut, then a mask: the same information the LLMs get
  full    the whole sequence with that one residue masked: more information than the LLMs
          get, an upper bound for these bidirectional models

    .venv-plm/bin/python scripts/22_esm_models.py data/confirm_sequences.tsv results/plm_confirm esmc-600m
"""

import json
import random
import sys
import zlib
from pathlib import Path

import pandas as pd
import torch

AA = "ACDEFGHIKLMNPQRSTVWY"
torch.set_grad_enabled(False)
torch.set_num_threads(16)


def shuffled(seq: str, accession: str) -> str:  # same as bioseq_probe.probe.shuffled
    residues = list(seq)
    random.Random(zlib.crc32(accession.encode())).shuffle(residues)
    return "".join(residues)


def esm_sdk(kind: str):
    from esm.sdk.api import ESMProtein, LogitsConfig
    from esm.tokenization import EsmSequenceTokenizer

    if kind == "esmc-600m":
        from esm.models.esmc import ESMC

        model = ESMC.from_pretrained("esmc_600m").to("cpu").eval()
    else:
        from esm.models.esm3 import ESM3

        model = ESM3.from_pretrained("esm3-open").to("cpu").float().eval()
    ids = EsmSequenceTokenizer().convert_tokens_to_ids(list(AA))

    def logp(seq_with_mask: str, pos: int):  # "_" marks the masked residue; +1 for BOS
        t = model.encode(ESMProtein(sequence=seq_with_mask))
        lg = model.logits(t, LogitsConfig(sequence=True)).logits.sequence[0, pos + 1, ids]
        return torch.log_softmax(lg.float(), -1)

    return logp


def esm2():
    from transformers import AutoTokenizer, EsmForMaskedLM

    name = "facebook/esm2_t33_650M_UR50D"
    model, tok = EsmForMaskedLM.from_pretrained(name).eval(), AutoTokenizer.from_pretrained(name)
    ids = tok.convert_tokens_to_ids(list(AA))

    def logp(seq_with_mask: str, pos: int):
        enc = tok(seq_with_mask.replace("_", tok.mask_token), return_tensors="pt")
        return torch.log_softmax(model(**enc).logits[0, pos + 1, ids], -1)

    return logp


LOADERS = {"esmc-600m": lambda: esm_sdk("esmc-600m"), "esm3-open": lambda: esm_sdk("esm3-open"), "esm2-650M": esm2}
MODES = {"prefix": lambda s, c: s[:c] + "_", "full": lambda s, c: s[:c] + "_" + s[c + 1 :]}

seqs = pd.read_csv(sys.argv[1], sep="\t")
out_dir = Path(sys.argv[2])
out_dir.mkdir(parents=True, exist_ok=True)
for mname in sys.argv[3].split(","):
    logp_fn = LOADERS[mname]()
    for mode, make in MODES.items():
        if mname == "esm2-650M" and mode == "prefix":
            continue  # scripts/10_plm.py already scores ESM-2 on the prefix
        path = out_dir / f"{mname}__{mode}.jsonl"
        done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
        with path.open("a") as f:
            for r in seqs.itertuples():
                for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession))):
                    for c in map(int, r.cuts.split(",")):
                        if (r.accession, c, cond) in done:
                            continue
                        lp = logp_fn(make(s, c), c)
                        pred = AA[int(lp.argmax())]
                        f.write(json.dumps({"model": mname, "mode": mode, "accession": r.accession, "pfam": r.pfam, "condition": cond,
                                            "cut": c, "target": s[c : c + 10], "pred": pred, "first_correct": pred == s[c], "logp": float(lp[AA.index(s[c])])}) + "\n")
                f.flush()
        print(mname, mode, "done", flush=True)
