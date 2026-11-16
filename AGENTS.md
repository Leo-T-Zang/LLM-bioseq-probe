# LLM-bioseq-probe for coding agents

- API keys live only in `.env` (see `.env.example`). Never read, print, copy or commit them.
- Every model call costs money. Scripts are dry by default; get a human's go-ahead, with
  the estimated cost (`scripts/tools/cost.py`), before passing `--run`.
- Model replies go to `results/<run>/` and are never edited by hand; every number in the
  docs comes from a script over them. Publish with `scripts/tools/export_outputs.py`.
- Run scripts from the repo root. Plots use matplotlib.
