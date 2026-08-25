# evades-search

Command-line tool to run the same two searches as the EVADES website's
"Analyse" tab — an HMM profile search ([HMMER](http://hmmer.org/)
`hmmsearch`) and a structural search ([Foldseek](https://github.com/steineggerlab/foldseek))
— against the EVADES anti-defence protein database, locally and offline
after a one-time database download. No account, no upload, no size
limits beyond your own machine.

This exists in response to a reviewer question on the EVADES manuscript
about whether users could run their own searches against the database
outside the website (in the spirit of tools like
[DefenseFinder](https://github.com/mdmparis/defense-finder)). It ports
the exact search logic from the website's backend
(`backend/app/tasks.py` in the `evades-webapp` repo), so a local run
against the same database produces the same hits as uploading the same
file on the website.

## Prerequisites

- Python 3.9+
- [HMMER](http://hmmer.org/) (`hmmsearch`, `hmmpress`):
  `brew install hmmer` or `conda install -c bioconda hmmer`
- [Foldseek](https://github.com/steineggerlab/foldseek):
  `brew install brewsci/bio/foldseek` or `conda install -c bioconda foldseek`

## Install

From a clone of this repo:

```bash
pipx install .        # isolated install, puts `evades-search` on your PATH
# or:
pip install .
```

For development (editable install + tests):

```bash
pip install -e '.[dev]'
pytest
```

## Quickstart

```bash
# One-time: download the EVADES HMM profile library and structure
# database (the same files served at /downloads/ on the EVADES website)
# and build the local HMMER/Foldseek indexes from them (~50 MB, cached
# under ~/.cache/evades-search).
evades-search fetch-db

# Search protein sequences against the HMM profile library.
evades-search hmm my_proteins.fasta

# Search a predicted structure against the Foldseek structure database.
evades-search structure my_protein.pdb
```

Check what's installed and where the local database lives:

```bash
evades-search info
```

## Output

`hmm` and `structure` both accept `-f/--format table|tsv|json` (default
`table`) and `-o/--output FILE` (default: stdout).

- `hmm` columns: `query_name, adp, moa, defence, evalue, score, hmm_from, hmm_to, ali_from, ali_to`
- `structure` columns: `query, adp, moa, defence, seq_identity, aln_len, prob, tm_score`

`adp` is the matched EVADES protein ID; `moa`/`defence` are its curated
mode-of-action and inhibited-defence-system annotations, joined in from
the same `metadata.tsv` the website uses.

## Matching the website's results

Thresholds are the same as the website by default: E-value ≤ 1e-5 for
`hmm`, TM-score ≥ 0.5 for `structure` (the standard "same fold"
convention, Zhang & Skolnick 2004). Both are overridable
(`--evalue`, `--tm-score-min`) if you want a looser or stricter cutoff
than the website uses. For structural search, chains/models belonging
to the same multimeric or NMR-ensemble target are collapsed to one hit
per protein, keeping the highest-confidence one — same as the website.

## Updating the database

The website's database occasionally changes as new proteins are added.
Re-sync your local copy with:

```bash
evades-search fetch-db --force
```

By default this pulls from the live EVADES server; point it elsewhere
with `--base-url` (or the `EVADES_SEARCH_BASE_URL` environment variable)
if the site moves.

## License

MIT — see `LICENSE`.
