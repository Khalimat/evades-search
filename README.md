# evades-search

Command-line search of the [EVADES](https://github.com/Khalimat/evades-webapp)
anti-defence protein database: an HMM profile search
([HMMER](http://hmmer.org/) `hmmsearch`) and a structural search
([Foldseek](https://github.com/steineggerlab/foldseek)), run locally and
offline against a downloaded copy of the database. No account, no
upload, no size limits beyond your own machine — useful for batch
searches, pipeline integration, or working with sequences/structures you
don't want to upload anywhere.

Results match the EVADES website's "Analyse" tab: same database, same
default thresholds, same hit-collapsing rules for multimeric/NMR
structures.

## Prerequisites

- Python 3.9+
- [HMMER](http://hmmer.org/) (`hmmsearch`, `hmmpress`):
  `brew install hmmer` or `conda install -c bioconda hmmer`
- [Foldseek](https://github.com/steineggerlab/foldseek):
  `brew install brewsci/bio/foldseek` or `conda install -c bioconda foldseek`

## Install

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
# database, and build the local HMMER/Foldseek indexes from them
# (~50 MB, cached under ~/.cache/evades-search).
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
`metadata.tsv`.

## Search thresholds

| | default | flag |
|---|---|---|
| `hmm` E-value cutoff | `1e-3` | `-E/--evalue` |
| `structure` minimum TM-score | `0.5` | `--tm-score-min` |

`1e-3` and `0.5` are the defaults used by the EVADES website itself.
TM-score `0.5` is the standard structural-biology convention for
"generally the same fold" (Zhang & Skolnick, 2004); below ~0.17 is
essentially random similarity. Raise or lower either flag for a
stricter or looser cutoff than the website's default.

For structural search, chains/models belonging to the same multimeric
or NMR-ensemble target are collapsed to one hit per protein, keeping
the highest-confidence one.

## Updating the database

The EVADES database occasionally changes as new proteins are added.
Re-sync your local copy with:

```bash
evades-search fetch-db --force
```

By default this pulls from the live EVADES server; point it elsewhere
with `--base-url` (or the `EVADES_SEARCH_BASE_URL` environment variable)
if the site moves.

## Suggesting new entries

Know of an anti-defence protein that isn't in the database yet? Open an
[issue](https://github.com/Khalimat/evades-search/issues) with the
protein and a link to the paper describing it.

## Cite us

EVADES: Encyclopaedia of bacterial virus anti-defence systems

Khalimat Murtazalieva, Evangelos Karatzas, Jiawei Wang, Robert D. Finn

## License

MIT — see `LICENSE`.
