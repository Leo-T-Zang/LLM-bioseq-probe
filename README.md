# bioseq-memorization

Do closed-source LLMs (Anthropic, OpenAI, Google) carry biological sequence
databases such as UniProt in their training data? This project probes them by
asking for verbatim recall of sequences that appear only in the databases,
with sequences released after each model's cutoff as the negative control.

The study design is still under discussion. Status: scaffold only.

## Setup

```bash
uv sync
uv run python scripts/list_models.py   # which models the keys reach; free
```

API keys come from `.env`, a symlink to `../phylo-evals/.env`.
