# bioseq-memorization for coding agents

- Credentials only in `.env` (a symlink to `../phylo-evals/.env`). Never read,
  print, copy or commit their values.
- A model call costs money. Listing models is free; get a human's go-ahead
  before any run that sends prompts, with its estimated cost.
- Raw model outputs go to `results/`, never edited by hand. Every number in a
  report comes from a script over them.
- Plots use matplotlib.
