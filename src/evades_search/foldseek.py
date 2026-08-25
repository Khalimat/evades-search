"""Structural search (Foldseek) against the EVADES structure database —
the same query the EVADES website's Analyse tab runs server-side.
Ported here so a local run against the same database produces the same
hits, including the TM-score filter and the multi-model/multi-chain
collapsing logic.
"""
from __future__ import annotations

import csv
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

from .metadata import Metadata

# TM-score >= 0.5 is the standard structural-biology convention for
# "generally the same fold" (Zhang & Skolnick, 2004); below ~0.17 is
# essentially random similarity. Matches the web app's default.
DEFAULT_TM_SCORE_MIN = 0.5

_STRUCTURE_SUFFIXES = (".pdb", ".cif")

# Some source PDBs (e.g. NMR structures) contain many models of the same
# chain; `foldseek createdb` indexes each model as its own entry, e.g.
# "klca_MODEL_16_A". Strip that suffix back down to the base name.
_MODEL_SUFFIX_RE = re.compile(r"_MODEL_\d+_[A-Za-z0-9]+$")


class FoldseekError(RuntimeError):
    pass


def foldseek_available() -> bool:
    return shutil.which("foldseek") is not None


def _strip_known_suffix(name: str, known_ids: frozenset[str]) -> str:
    """Strip Foldseek's appended "_MODEL_n_<chain>" and/or "_<chain>"
    suffixes, matching against a known set of base identifiers — used
    both for resolving a target hit back to its EVADES protein ID, and
    for resolving a query hit back to its source structure file. Some
    identifiers themselves contain underscores (e.g. "gp5_9",
    "acrvia2_plus", even a literal trailing "_" as in "acrif18_"), so a
    blind regex strip of the last "_..." segment would mangle those —
    match against the known ID set instead."""
    name = _MODEL_SUFFIX_RE.sub("", name)
    if name in known_ids:
        return name
    if "_" in name:
        prefix = name.rsplit("_", 1)[0]
        if prefix in known_ids:
            return prefix
    return name


def _base_protein_name(target: str, known_ids: frozenset[str]) -> str:
    return _strip_known_suffix(target, known_ids)


def _resolve_structure_files(paths: Sequence[Path]) -> list[Path]:
    """Expand a mix of structure files and directories into a flat,
    de-duplicated, sorted list of .pdb/.cif files — directories are
    searched recursively."""
    files: list[Path] = []
    for p in paths:
        if not p.exists():
            raise FoldseekError(f"Path not found: {p}")
        if p.is_dir():
            found = sorted(q for q in p.rglob("*") if q.suffix.lower() in _STRUCTURE_SUFFIXES)
            if not found:
                raise FoldseekError(f"No .pdb/.cif files found under directory: {p}")
            files.extend(found)
        elif p.suffix.lower() in _STRUCTURE_SUFFIXES:
            files.append(p)
        else:
            raise FoldseekError(f"Not a .pdb/.cif file or directory: {p}")
    # De-duplicate while preserving order (the same file could be
    # reachable twice if paths overlap, e.g. a file and its parent dir
    # both passed).
    seen = set()
    deduped = []
    for f in files:
        resolved = f.resolve()
        if resolved not in seen:
            seen.add(resolved)
            deduped.append(f)
    return deduped


def run_foldseek(
    structure_paths: Path | Sequence[Path],
    foldseek_db: Path,
    metadata: Metadata,
    *,
    tm_score_min: float = DEFAULT_TM_SCORE_MIN,
    threads: int = 2,
    timeout: int = 250,
) -> list[dict]:
    """Search one or more query structures (files and/or directories of
    files) against the EVADES Foldseek database in a single batched
    Foldseek run."""
    if not foldseek_available():
        raise FoldseekError(
            "foldseek not found on PATH. Install Foldseek first "
            "(`brew install brewsci/bio/foldseek` or `conda install -c bioconda foldseek`)."
        )
    if not foldseek_db.with_suffix("").exists() and not Path(str(foldseek_db)).exists():
        raise FoldseekError(
            f"Foldseek database not found at {foldseek_db}. Run `evades-search fetch-db` first."
        )

    paths = [structure_paths] if isinstance(structure_paths, Path) else list(structure_paths)
    files = _resolve_structure_files(paths)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)

        # Stage every query file under a query/ subdir with a
        # collision-proof name — Foldseek's "query" output column is
        # derived from filename stem + chain, and it drops the
        # directory a file came from, so two different files with the
        # same stem (e.g. two "candidate.pdb"s from different folders)
        # would otherwise become indistinguishable in the results.
        query_dir = tmp_path / "query"
        query_dir.mkdir()
        stem_to_source: dict[str, Path] = {}
        for i, f in enumerate(files):
            staged_stem = f"{i:04d}_{f.stem}"
            (query_dir / f"{staged_stem}{f.suffix}").symlink_to(f.resolve())
            stem_to_source[staged_stem] = f

        results_tsv = tmp_path / "results.tsv"
        cmd = [
            "foldseek", "easy-search",
            "--threads", str(threads),
            str(query_dir),
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

        return _parse_results(results_tsv, metadata, tm_score_min, stem_to_source)


def _parse_results(
    results_tsv: Path,
    metadata: Metadata,
    tm_score_min: float,
    stem_to_source: dict[str, Path],
) -> list[dict]:
    hits = []
    if results_tsv.exists():
        with results_tsv.open() as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                if len(row) < 6:
                    continue
                hits.append({
                    "raw_query": row[0],
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
    staged_stems = frozenset(stem_to_source)
    query_order = {stem: i for i, stem in enumerate(stem_to_source)}

    for hit in hits:
        staged_stem = _strip_known_suffix(hit["raw_query"], staged_stems)
        chain_suffix = hit["raw_query"][len(staged_stem):]
        source = stem_to_source[staged_stem]
        hit["query_file"] = str(source)
        hit["query"] = f"{source.stem}{chain_suffix}"
        hit["_query_order"] = query_order[staged_stem]

    # Collapse multi-model/multi-chain target duplicates down to one
    # hit per protein — scoped *per source query file*, so two
    # different query structures that both happen to match the same
    # EVADES protein each keep their own hit rather than one clobbering
    # the other.
    best_by_protein: dict[tuple[str, str], dict] = {}
    for hit in hits:
        base_name = _base_protein_name(hit["adp"], known_ids)
        hit["adp"] = base_name
        key = (hit["query_file"], base_name)
        existing = best_by_protein.get(key)
        if existing is None or hit["prob"] > existing["prob"]:
            best_by_protein[key] = hit

    ordered = sorted(best_by_protein.values(), key=lambda h: (h["_query_order"], -h["prob"]))
    result = []
    for h in ordered:
        result.append({
            "query": h["query"],
            "query_file": h["query_file"],
            "adp": h["adp"],
            "moa": metadata.moa(h["adp"]),
            "defence": metadata.defence(h["adp"]),
            "seq_identity": h["seq_identity"],
            "aln_len": h["aln_len"],
            "prob": h["prob"],
            "tm_score": h["tm_score"],
        })
    return result
