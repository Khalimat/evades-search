"""Structural search (Foldseek) against the EVADES structure database —
the same query the web app's Analyse tab runs server-side (see
`backend/app/tasks.py::run_foldseek` in the evades-webapp repo). Ported
here so a local run against the same database produces the same hits,
including the TM-score filter and the multi-model/multi-chain
collapsing logic.
"""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from .metadata import Metadata

# TM-score >= 0.5 is the standard structural-biology convention for
# "generally the same fold" (Zhang & Skolnick, 2004); below ~0.17 is
# essentially random similarity. Matches the web app's default.
DEFAULT_TM_SCORE_MIN = 0.5

# Some source PDBs (e.g. NMR structures) contain many models of the same
# chain; `foldseek createdb` indexes each model as its own entry, e.g.
# "klca_MODEL_16_A". Strip that suffix back down to the base name.
_MODEL_SUFFIX_RE = re.compile(r"_MODEL_\d+_[A-Za-z0-9]+$")


class FoldseekError(RuntimeError):
    pass


def foldseek_available() -> bool:
    return shutil.which("foldseek") is not None


def _base_protein_name(target: str, known_ids: frozenset[str]) -> str:
    """Foldseek also appends "_<chain id>" for each chain of a
    multi-chain (multimer) structure, e.g. "acrib4_A", "acrib4_B" — on
    top of the "_MODEL_n_<chain>" suffix for multi-model NMR ensembles.
    Some protein IDs themselves contain underscores (e.g. "gp5_9",
    "acrvia2_plus", even a literal trailing "_" as in "acrif18_"), so a
    blind regex strip of the last "_..." segment would mangle those —
    match against the known ID set instead."""
    target = _MODEL_SUFFIX_RE.sub("", target)
    if target in known_ids:
        return target
    if "_" in target:
        prefix = target.rsplit("_", 1)[0]
        if prefix in known_ids:
            return prefix
    return target


def run_foldseek(
    structure_path: Path,
    foldseek_db: Path,
    metadata: Metadata,
    *,
    tm_score_min: float = DEFAULT_TM_SCORE_MIN,
    threads: int = 2,
    timeout: int = 250,
) -> list[dict]:
    if not foldseek_available():
        raise FoldseekError(
            "foldseek not found on PATH. Install Foldseek first "
            "(`brew install brewsci/bio/foldseek` or `conda install -c bioconda foldseek`)."
        )
    if not foldseek_db.with_suffix("").exists() and not Path(str(foldseek_db)).exists():
        raise FoldseekError(
            f"Foldseek database not found at {foldseek_db}. Run `evades-search fetch-db` first."
        )
    if not structure_path.exists():
        raise FoldseekError(f"Structure file not found: {structure_path}")

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        results_tsv = tmp_path / "results.tsv"
        cmd = [
            "foldseek", "easy-search",
            "--threads", str(threads),
            str(structure_path),
            str(foldseek_db),
            str(results_tsv),
            str(tmp_path / "tmp"),
            "--alignment-type", "1",   # exact TM-align (Kabsch superposition)
            "--tmscore-threshold", str(tm_score_min),
            # NOTE: no -e filter — under --alignment-type 1, Foldseek
            # redefines "evalue" as (qTMscore+tTMscore)/2 (higher is
            # better), not a real e-value. --tmscore-threshold is the
            # correct filter here.
            "--format-output",
            "query,target,fident,alnlen,prob,alntmscore",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            raise FoldseekError(f"foldseek failed: {proc.stderr[-2000:]}")

        return _parse_results(results_tsv, metadata, tm_score_min)


def _parse_results(results_tsv: Path, metadata: Metadata, tm_score_min: float) -> list[dict]:
    hits = []
    if results_tsv.exists():
        with results_tsv.open() as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 6:
                    continue
                hits.append({
                    "query": row[0],
                    "adp": row[1],
                    "seq_identity": float(row[2]),
                    "aln_len": int(row[3]),
                    "prob": float(row[4]),
                    # TM-score is mathematically bounded to (0, 1];
                    # floating-point rounding on near-perfect matches
                    # can occasionally push it a hair above 1.0.
                    "tm_score": min(float(row[5]), 1.0),
                })

    # Belt-and-braces re-check against the alntmscore field (already
    # enforced by --tmscore-threshold at the tool level).
    hits = [h for h in hits if h["tm_score"] >= tm_score_min]

    known_ids = metadata.known_ids()
    best_by_protein: dict[str, dict] = {}
    for hit in hits:
        base_name = _base_protein_name(hit["adp"], known_ids)
        hit["adp"] = base_name
        existing = best_by_protein.get(base_name)
        if existing is None or hit["prob"] > existing["prob"]:
            best_by_protein[base_name] = hit

    hits = sorted(best_by_protein.values(), key=lambda h: -h["prob"])
    result = []
    for h in hits:
        result.append({
            "query": h["query"],
            "adp": h["adp"],
            "moa": metadata.moa(h["adp"]),
            "defence": metadata.defence(h["adp"]),
            "seq_identity": h["seq_identity"],
            "aln_len": h["aln_len"],
            "prob": h["prob"],
            "tm_score": h["tm_score"],
        })
    return result
