"""Which APIs can continue a FASTA as raw text (prefill or completion endpoints)? Costs cents.

Each test gives the GFP sequence up to ...LVTTF and checks whether the model continues it
(truth: SYGVQCFSRY). Prints one line per (provider, model, method).
"""

import os
from concurrent.futures import ThreadPoolExecutor

import anthropic
import httpx
import openai
from dotenv import load_dotenv
from google import genai
from google.genai import types

from bioseq_mem.probe import fasta

load_dotenv()
GFP = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTF"
FA = ">sp|P42212|GFP_AEQVI Green fluorescent protein\n" + fasta(GFP)
ASK = "Continue this protein sequence."
T = 90.0


def anth(model):
    try:
        r = anthropic.Anthropic(timeout=T).messages.create(model=model, max_tokens=64, messages=[
            {"role": "user", "content": ASK}, {"role": "assistant", "content": FA}])
        return "anthropic prefill", model, "".join(b.text for b in r.content if b.type == "text")
    except Exception as e:
        return "anthropic prefill", model, f"ERR {e}"


def oai_completions(model):
    try:
        r = openai.OpenAI(timeout=T).completions.create(model=model, prompt=FA, max_tokens=24, temperature=0)
        return "openai /v1/completions", model, r.choices[0].text
    except Exception as e:
        return "openai /v1/completions", model, f"ERR {e}"


def oai_chat_prefill(model):
    try:
        r = openai.OpenAI(timeout=T).chat.completions.create(model=model, max_completion_tokens=2000, messages=[
            {"role": "user", "content": ASK}, {"role": "assistant", "content": FA}])
        return "openai chat, trailing assistant msg", model, r.choices[0].message.content
    except Exception as e:
        return "openai chat, trailing assistant msg", model, f"ERR {e}"


def gem(model):
    try:
        c = genai.Client(api_key=os.environ["GEMINI_API_KEY"], http_options=types.HttpOptions(api_version="v1beta", timeout=int(T * 1000)))
        r = c.models.generate_content(model=model, contents=[
            types.Content(role="user", parts=[types.Part(text=ASK)]), types.Content(role="model", parts=[types.Part(text=FA)])],
            config=types.GenerateContentConfig(temperature=0, max_output_tokens=2000))
        return "gemini, trailing model turn", model, r.text
    except Exception as e:
        return "gemini, trailing model turn", model, f"ERR {e}"


def fw_completions(model):
    try:
        r = httpx.post("https://api.fireworks.ai/inference/v1/completions", timeout=T,
                       headers={"Authorization": f"Bearer {os.environ['FIREWORKS_API_KEY']}"},
                       json={"model": f"accounts/fireworks/models/{model}", "prompt": FA, "max_tokens": 24, "temperature": 0})
        j = r.json()
        return "fireworks /completions (raw)", model, j["choices"][0]["text"] if "choices" in j else f"ERR {j}"
    except Exception as e:
        return "fireworks /completions (raw)", model, f"ERR {e}"


def fw_chat_prefill(model):
    try:
        r = httpx.post("https://api.fireworks.ai/inference/v1/chat/completions", timeout=T,
                       headers={"Authorization": f"Bearer {os.environ['FIREWORKS_API_KEY']}"},
                       json={"model": f"accounts/fireworks/models/{model}", "max_tokens": 24, "temperature": 0,
                             "messages": [{"role": "user", "content": ASK}, {"role": "assistant", "content": FA}]})
        j = r.json()
        return "fireworks chat, trailing assistant msg", model, j["choices"][0]["message"]["content"] if "choices" in j else f"ERR {j}"
    except Exception as e:
        return "fireworks chat, trailing assistant msg", model, f"ERR {e}"


jobs = [(anth, m) for m in ("claude-opus-5-5", "claude-opus-5", "claude-sonnet-4-5-20250929", "claude-haiku-4-5-20251001")]
jobs += [(oai_completions, m) for m in ("gpt-5.5", "gpt-4.1", "gpt-4o", "gpt-3.5-turbo-instruct", "davinci-002", "babbage-002")]
jobs += [(oai_chat_prefill, m) for m in ("gpt-4.1", "gpt-5.5")]
jobs += [(gem, m) for m in ("gemini-4-carbon", "gemini-3.5-flash", "gemini-2.5-pro")]
open_models = ("deepseek-v4p1-flash", "qwen3p8-2p4t-a95b", "kimi-k3p1-pro", "glm-5p3", "minimax-m3", "gpt-oss-120b", "nemotron-3-ultra-nvfp4")
jobs += [(fw_completions, m) for m in open_models] + [(fw_chat_prefill, m) for m in open_models]
with ThreadPoolExecutor(16) as pool:
    for how, m, out in pool.map(lambda j: j[0](j[1]), jobs):
        print(f"{how:38s} {m:28s} {out.strip()[:110]!r}")
print("truth continues: SYGVQCFSRYPDHMKQ")
