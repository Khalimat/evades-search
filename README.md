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
- [HMMER](http://hmmer.org/) (`hmmsearch`, `hmmpress`)
- [Foldseek](https://github.com/steineggerlab/foldseek)

Neither requires conda — pick whichever of these fits your setup:

**macOS (Homebrew)**
```bash
brew install hmmer
brew install brewsci/bio/foldseek   # foldseek isn't in homebrew-core, needs this tap
```

**Linux (Debian/Ubuntu)**
```bash
sudo apt install hmmer              # packaged directly; foldseek is not
```
Foldseek doesn't have a Linux package; grab its static binary instead (no
root needed — this just unpacks a tarball and adds it to your `PATH`):
```bash
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz   # or foldseek-linux-arm64.tar.gz on ARM
tar xvzf foldseek-linux-avx2.tar.gz
export PATH="$(pwd)/foldseek/bin:$PATH"   # add to your shell profile to keep it
```

**conda/mamba (any OS)**
```bash
conda install -c conda-forge -c bioconda hmmer foldseek
```

**HPC / shared cluster**: both tools are also published as
[biocontainer](https://biocontainers.pro/) images
(`quay.io/biocontainers/hmmer`, `quay.io/biocontainers/foldseek`) if
Singularity/Apptainer is what's available to you instead of a package
manager.

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
# database from Zenodo, and build the local HMMER/Foldseek indexes
# from them (~50 MB, cached under ~/.cache/evades-search).
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

By default, `fetch-db` pulls from a permanent Zenodo deposit of the
bulk-download files (DOI: [10.5281/zenodo.22096345](https://doi.org/10.5281/zenodo.22096345)),
independent of the EVADES website's own server. Point it elsewhere with
`--base-url` (or the `EVADES_SEARCH_BASE_URL` environment variable) —
e.g. to pull from a different EVADES deployment's own `/downloads/`
directory instead.

When the underlying dataset changes (new proteins added, corrected
structures, etc.), a new version gets deposited to Zenodo and this
default is updated to point at it; re-sync your local copy with:

```bash
evades-search fetch-db --force
```

## Suggesting new entries

Know of an anti-defence protein that isn't in the database yet? Open an
[issue](https://github.com/Khalimat/evades-search/issues) with the
protein and a link to the paper describing it.

## Cite us

EVADES: Encyclopaedia of bacterial virus anti-defence systems

Khalimat Murtazalieva, Evangelos Karatzas, Jiawei Wang, Robert D. Finn

The underlying database is separately archived on Zenodo:
[10.5281/zenodo.22096345](https://doi.org/10.5281/zenodo.22096345).

## License

MIT — see `LICENSE`.
