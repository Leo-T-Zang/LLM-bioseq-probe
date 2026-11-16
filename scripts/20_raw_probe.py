"""Protein probes as raw continuation, the closest these APIs get to plain next-token
prediction -> results/raw/<model>.jsonl. Dry by default; --run makes paid calls; resumable.

The model receives ">sequence\\n" + the FASTA-wrapped prefix as text to continue:
  fireworks  /completions endpoint, no chat template (open-weight models)
  anthropic  assistant-turn prefill (Claude Sonnet 4.5; Opus 5 and 5.5 reject prefill)
  gemini     request ending with a model turn (Gemini 2.5 Pro; 3.5 and 4 reject it)
The prediction is the first 10 residues it writes, newlines ignored, stopping at '>'.
GPT models have no completion endpoint or prefill, so they are not here.

    uv run python scripts/20_raw_probe.py --models claude-sonnet-4-5-20250929 --run
"""

import argparse
import json
import os
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import anthropic
import httpx
import pandas as pd
from dotenv import load_dotenv
from google import genai
from google.genai import types

from bioseq_probe.probe import K, STANDARD, fasta, score, shuffled

load_dotenv()
ASK = "Continue this protein sequence."
MODELS = {
    "claude-sonnet-4-5-20250929": "anthropic",
    "gemini-2.5-pro": "gemini",
    **{m: "fireworks" for m in ("deepseek-v4p1-flash", "qwen3p8-2p4t-a95b", "kimi-k3p1-pro", "glm-5p3",
                                "minimax-m3", "gpt-oss-120b", "nemotron-3-ultra-nvfp4")},
}


def text_of(prefix: str) -> str:
    t = ">sequence\n" + fasta(prefix)
    return t + "\n" if len(prefix) % 60 == 0 else t


def continue_raw(model: str, prefix: str) -> tuple[str, dict, str | None]:
    """(continuation text, usage, stop reason); a refusal is kept apart from a miss."""
    kind, t = MODELS[model], text_of(prefix)
    if kind == "fireworks":
        for attempt in range(8):  # back off on rate limits (429)
            r = httpx.post("https://api.fireworks.ai/inference/v1/completions", timeout=120,
                           headers={"Authorization": f"Bearer {os.environ['FIREWORKS_API_KEY']}"},
                           json={"model": f"accounts/fireworks/models/{model}", "prompt": t, "max_tokens": 24, "temperature": 0})
            if r.status_code != 429:
                break
            time.sleep(min(2**attempt, 60))
        r.raise_for_status()
        j = r.json()
        return j["choices"][0]["text"], j.get("usage", {}), j["choices"][0].get("finish_reason")
    if kind == "anthropic":
        r = anthropic.Anthropic(timeout=120).messages.create(model=model, max_tokens=32, messages=[
            {"role": "user", "content": ASK}, {"role": "assistant", "content": t.rstrip("\n")}])
        return "".join(b.text for b in r.content if b.type == "text"), r.usage.model_dump(), r.stop_reason
    c = genai.Client(api_key=os.environ["GEMINI_API_KEY"], http_options=types.HttpOptions(api_version="v1beta", timeout=120_000))
    r = c.models.generate_content(model=model, contents=[
        types.Content(role="user", parts=[types.Part(text=ASK)]), types.Content(role="model", parts=[types.Part(text=t)])],
        config=types.GenerateContentConfig(temperature=0, max_output_tokens=2000, thinking_config=types.ThinkingConfig(thinking_budget=128)))
    stop = str(r.candidates[0].finish_reason) if r.candidates else f"blocked: {r.prompt_feedback}"
    return r.text or "", r.usage_metadata.model_dump() if r.usage_metadata else {}, stop


def parse_raw(text: str) -> str:
    """Residues written before any '>' or blank line, newlines and spaces dropped."""
    body = re.split(r">|\n\s*\n", text, maxsplit=1)[0]
    letters = re.sub(r"[\s]", "", body).upper()
    letters = re.match(r"[A-Z]*", letters).group(0)
    return letters[:K] if letters and set(letters[:K]) <= set(STANDARD) else ""


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default=",".join(MODELS))
    ap.add_argument("--sequences", default="data/confirm_sequences.tsv")
    ap.add_argument("--out", default="results/raw")
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()

    seqs = pd.read_csv(args.sequences, sep="\t")
    probes = [
        {"accession": r.accession, "pfam": r.pfam, "condition": cond, "cut": int(c), "prefix": s[: int(c)], "target": s[int(c) : int(c) + K]}
        for r in seqs.itertuples()
        for c in r.cuts.split(",")
        for cond, s in (("real", r.sequence), ("shuffled", shuffled(r.sequence, r.accession)))
    ]
    print(f"{len(probes)} probes per model")
    if not args.run:
        raise SystemExit("dry run; pass --run")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for model in args.models.split(","):
        path = out / f"{model}.jsonl"
        done = {(d["accession"], d["cut"], d["condition"]) for d in map(json.loads, path.open())} if path.exists() else set()
        todo = [p for p in probes if (p["accession"], p["cut"], p["condition"]) not in done]

        def one(p, model=model):
            try:
                text, usage, stop = continue_raw(model, p["prefix"])
                rec = {"text": text, "usage": usage, "stop": stop}
            except Exception as e:
                rec = {"error": f"{type(e).__name__}: {e}"[:500]}
            pred = parse_raw(rec.get("text", ""))
            return {"model": model, **{k: p[k] for k in ("accession", "pfam", "condition", "cut", "target")}, "pred": pred, **score(pred, p["target"]), **rec}

        with path.open("a") as f, ThreadPoolExecutor(args.workers) as pool:
            for rec in pool.map(one, todo):
                f.write(json.dumps(rec) + "\n")
                f.flush()
        print(model, "done", flush=True)
