"""List the models each provider key in .env can reach. Listing is free."""

import os

import anthropic
import httpx
import openai
from dotenv import load_dotenv
from google import genai

load_dotenv()


def show(name, ids):
    print(f"== {name} ({len(ids)})")
    for i in sorted(ids):
        print("  ", i)


try:
    show("anthropic", [m.id for m in anthropic.Anthropic().models.list()])
except Exception as e:
    print("anthropic:", type(e).__name__, e)

try:
    ids = [m.id for m in openai.OpenAI().models.list()]
    show("openai", [i for i in ids if i.startswith(("gpt", "o", "davinci", "babbage"))])
except Exception as e:
    print("openai:", type(e).__name__, e)

try:
    client = genai.Client(api_key=os.environ["GEMINI_API_KEY"])
    show("gemini", [m.name for m in client.models.list() if "gemini" in m.name])
except Exception as e:
    print("gemini:", type(e).__name__, e)

try:
    r = httpx.get("https://openrouter.ai/api/v1/models", timeout=30)
    ids = [m["id"] for m in r.json()["data"]]
    show("openrouter (closed-source labs only)", [i for i in ids if i.split("/")[0] in {"openai", "anthropic", "google", "x-ai"}])
except Exception as e:
    print("openrouter:", type(e).__name__, e)
