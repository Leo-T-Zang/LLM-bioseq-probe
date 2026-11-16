"""One call shape across Anthropic, OpenAI and Gemini: prompt in, text (+ logprobs) out.

Each model's settings come from MODELS. Sampling is as close to greedy as each
API allows: temperature 0 where accepted; reasoning at its lowest setting where
it cannot be turned off.
"""

import os
from dataclasses import dataclass, field

import anthropic
import openai
from dotenv import load_dotenv
from google import genai
from google.genai import types

load_dotenv()

# name -> (provider, settings). Proposed pilot list; confirm before a paid run.
MODELS: dict[str, tuple[str, dict]] = {
    "gpt-4o-2024-11-20": ("openai_chat", {}),
    "gpt-4.1-2025-04-14": ("openai_chat", {}),
    "gpt-5.5": ("openai_responses", {"reasoning": {"effort": "none"}}),
    "gpt-6.1-sol": ("openai_responses", {"reasoning": {"effort": "low"}}),  # lowest it accepts
    "claude-sonnet-4-5-20250929": ("anthropic", {}),
    "claude-opus-4-6": ("anthropic", {"output_config": {"effort": "low"}}),
    "claude-opus-5": ("anthropic", {"thinking": {"type": "disabled"}, "output_config": {"effort": "low"}}),
    "claude-opus-5-5": ("anthropic", {"output_config": {"effort": "low"}}),
    "gemini-2.5-pro": ("gemini", {"thinking_config": {"thinking_budget": 128}}),
    "gemini-3.5-flash": ("gemini", {"thinking_config": {"thinking_level": "low"}}),
    "gemini-3.8-flash": ("gemini", {"thinking_config": {"thinking_level": "low"}}),
    "gemini-4-carbon": ("gemini", {"thinking_config": {"thinking_level": "low"}}),
}


@dataclass
class Reply:
    text: str
    stop: str | None = None  # provider stop reason; "refusal" is kept apart from a miss
    top_logprobs: list[dict[str, float]] = field(default_factory=list)  # per output token
    usage: dict = field(default_factory=dict)


def call(model: str, prompt: str, max_tokens: int = 4000) -> Reply:
    provider, cfg = MODELS[model]
    return _CALLS[provider](model, prompt, max_tokens, cfg)


def _anthropic(model, prompt, max_tokens, cfg):
    r = anthropic.Anthropic().messages.create(
        model=model, max_tokens=max_tokens, messages=[{"role": "user", "content": prompt}], **cfg
    )
    text = "".join(b.text for b in r.content if b.type == "text")
    return Reply(text, r.stop_reason, usage=r.usage.model_dump())


def _openai_chat(model, prompt, max_tokens, cfg):
    r = openai.OpenAI().chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
        temperature=0,
        max_tokens=64,
        logprobs=True,
        top_logprobs=20,
        **cfg,
    )
    c = r.choices[0]
    tops = [{t.token: t.logprob for t in tok.top_logprobs} for tok in (c.logprobs.content or [])]
    return Reply(c.message.content or "", c.finish_reason, tops, r.usage.model_dump())


def _openai_responses(model, prompt, max_tokens, cfg):
    r = openai.OpenAI().responses.create(model=model, input=prompt, max_output_tokens=max_tokens, **cfg)
    return Reply(r.output_text or "", r.status, usage=r.usage.model_dump() if r.usage else {})


def _gemini(model, prompt, max_tokens, cfg):
    client = genai.Client(api_key=os.environ["GEMINI_API_KEY"], http_options=types.HttpOptions(api_version="v1beta"))
    conf = types.GenerateContentConfig(
        temperature=0, max_output_tokens=max_tokens, thinking_config=types.ThinkingConfig(**cfg["thinking_config"])
    )
    r = client.models.generate_content(model=model, contents=prompt, config=conf)
    cand = r.candidates[0] if r.candidates else None
    stop = str(cand.finish_reason) if cand else "no_candidate"
    return Reply(r.text or "", stop, usage=r.usage_metadata.model_dump() if r.usage_metadata else {})


_CALLS = {
    "anthropic": _anthropic,
    "openai_chat": _openai_chat,
    "openai_responses": _openai_responses,
    "gemini": _gemini,
}
