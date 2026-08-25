# evades-search

A command-line tool to search EVADES (Encyclopaedia of bacterial virus
anti-defence systems) locally, using
[HMMER](http://hmmer.org/) `hmmsearch` for protein sequences and
[Foldseek](https://github.com/steineggerlab/foldseek) for structures.

- `hmm` takes a protein FASTA (single- or multi-sequence).
- `structure` takes one or more `.pdb`/`.cif` files, or a directory of
  them, including multi-chain/multimer structures.

## Prerequisites

- Python 3.9+
- [HMMER](http://hmmer.org/) (`hmmsearch`, `hmmpress`)
- [Foldseek](https://github.com/steineggerlab/foldseek)

**macOS**
```bash
brew install hmmer
brew install brewsci/bio/foldseek
```

**Linux**
```bash
sudo apt install hmmer
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz   # foldseek-linux-arm64.tar.gz on ARM
tar xvzf foldseek-linux-avx2.tar.gz
export PATH="$(pwd)/foldseek/bin:$PATH"
```

No root access (e.g. a shared cluster)? Build HMMER from source into your
home directory instead:
```bash
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar zxf hmmer.tar.gz && cd hmmer-*
./configure --prefix="$HOME/.local"
make && make install
export PATH="$HOME/.local/bin:$PATH"
```
Foldseek's static binary above already needs no root.

**conda/mamba** (creates a ready-to-use environment: Python, HMMER, Foldseek)
```bash
conda create -n evades-search -c conda-forge -c bioconda python=3.11 hmmer foldseek
conda activate evades-search
pip install .   # from a clone of this repo — see Install below
```

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

The database downloads to `~/.cache/evades-search` by default. To use a
different location instead, set `EVADES_SEARCH_CACHE_DIR` *before*
running `fetch-db` — a plain `export` only lasts for the current
terminal session, so add it to your shell's startup file to make it
stick:

```bash
echo 'export EVADES_SEARCH_CACHE_DIR=/path/to/somewhere' >> ~/.zshrc && source ~/.zshrc
```

(zsh shown above — it's the default shell on current macOS; run
`echo $SHELL` if you're not sure which you're on. bash: same line, but
in `~/.bashrc`. fish: no file to edit, just run
`set -Ux EVADES_SEARCH_CACHE_DIR /path/to/somewhere` once.)

```bash
# One-time: download the EVADES HMM profile library and structure
# database, and build the local HMMER/Foldseek indexes from them.
evades-search fetch-db

# Search protein sequences against the HMM profile library.
evades-search hmm my_proteins.fasta

# Search a single predicted structure.
evades-search structure my_protein.pdb

# Or search every structure in a directory at once, in one batched run.
evades-search structure my_structures_dir/
```

Check what's installed and where the local database lives:

```bash
evades-search info
```

## Output

`hmm` and `structure` both accept `-f/--format tsv|json|table` (default
`tsv`) and `-o/--output FILE`.

By default, results are written to a file next to the input — same
name, extension swapped for the format (e.g. `my_protein.pdb` →
`my_protein.tsv`; a directory input becomes `<dirname>.tsv` alongside
it). Pass `-o` to name the file yourself, or `-o -` for stdout.

- `hmm` columns: `query_name, adp, moa, defence, evalue, score, hmm_from, hmm_to, ali_from, ali_to`
- `structure` columns: `query, query_file, adp, moa, defence, seq_identity, aln_len, prob, tm_score`

`adp` is the matched EVADES protein ID; `moa`/`defence` are its curated
mode-of-action and inhibited-defence-system annotations. `query_file`
is which input file a hit came from; `query` adds the chain.

## Search thresholds

| | default | flag |
|---|---|---|
| `hmm` E-value cutoff | `1e-3` | `-E/--evalue` |
| `structure` minimum TM-score | `0.5` | `--tm-score-min` |

These match the EVADES website's own defaults; raise or lower either
flag for a stricter or looser cutoff.

Chains/models belonging to the same multimeric or NMR-ensemble query
structure are collapsed to one hit per protein, keeping the
highest-confidence one. When searching many structures at once, this
collapsing is per query file — two different files matching the same
EVADES protein each keep their own hit.

## Updating the database

By default, `fetch-db` pulls three files from a permanent Zenodo deposit
(DOI: [10.5281/zenodo.22096345](https://doi.org/10.5281/zenodo.22096345)):
`hmm_profiles.tar.gz`, `metadata.tsv`, and
`foldseek_monomer_structures.tar.gz`. Point it elsewhere with
`--base-url` (or the `EVADES_SEARCH_BASE_URL` environment variable).

Re-sync your local copy after the dataset updates:

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

The database this tool searches against is separately archived on
Zenodo: [10.5281/zenodo.22096345](https://doi.org/10.5281/zenodo.22096345).

## License

MIT — see `LICENSE`.
