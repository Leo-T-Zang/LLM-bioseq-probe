"""scripts/probe/protein_lm.py on Modal GPUs, one model per GPU -> results/plm_confirm/. About $1.

The weights are public on Hugging Face (biohub/ESMC-600M, biohub/esm3-sm-open-v1, AI4PD/ProtGPT3-*),
so no token is needed.

    uv run --with modal modal run scripts/probe/protein_lm_modal.py --runs esmc-600m:prefix,esm3-open:prefix,protgpt3-1.3B:prefix,protgpt3-10B:prefix
"""

from pathlib import Path

import modal

app = modal.App("llm-bioseq-probe-plm")
image = (
    modal.Image.debian_slim(python_version="3.12")
    .pip_install("esm==3.4.3", "transformers", "pandas")
    .add_local_file("scripts/probe/protein_lm.py", "/work/protein_lm.py")
    .add_local_file("data/confirm_sequences.tsv", "/work/data/confirm_sequences.tsv")
)


@app.function(image=image, gpu=["L40S", "A100", "A100-80GB", "H100"], timeout=3 * 3600)
def run(spec: str) -> dict[str, str]:
    import os
    import sys

    os.chdir("/work")
    sys.path.insert(0, "/work")
    import protein_lm

    protein_lm.main(spec)
    return {p.name: p.read_text() for p in Path("results/plm_confirm").glob("*.jsonl")}


@app.local_entrypoint()
def main(runs: str):
    out = Path("results/plm_confirm")
    out.mkdir(parents=True, exist_ok=True)
    for files in run.map(runs.split(",")):
        for name, text in files.items():
            (out / name).write_text(text)
            print(name, text.count("\n"), "rows")
