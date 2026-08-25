"""HMM profile search (HMMER `hmmsearch`) against the EVADES profile
library — the same query the EVADES website's Analyse tab runs
server-side. Ported here so a local run against the same database
produces the same hits, with the same E-value cutoff and column set.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path

from .metadata import Metadata

DEFAULT_EVALUE = 1e-3  # matches the web app; not every profile has a
                        # curated gathering (GA) threshold, so -E is
                        # used instead of --cut_ga
_PROFILE_SUFFIX_RE = re.compile(r"\.aln$")


class HmmSearchError(RuntimeError):
    pass


def hmmsearch_available() -> bool:
    return shutil.which("hmmsearch") is not None


def run_hmmsearch(
    fasta_path: Path,
    hmm_db: Path,
    metadata: Metadata,
    *,
    evalue: float = DEFAULT_EVALUE,
    cpu: int = 2,
    timeout: int = 600,
) -> list[dict]:
    if not hmmsearch_available():
        raise HmmSearchError(
            "hmmsearch not found on PATH. Install HMMER first "
            "(`brew install hmmer` or `conda install -c bioconda hmmer`)."
        )
    if not hmm_db.exists():
        raise HmmSearchError(
            f"HMM database not found at {hmm_db}. Run `evades-search fetch-db` first."
        )
    if not fasta_path.exists():
        raise HmmSearchError(f"FASTA file not found: {fasta_path}")

    with tempfile.TemporaryDirectory() as tmp:
        domtblout = Path(tmp) / "hits.domtblout"
        cmd = [
            "hmmsearch",
            "--domtblout", str(domtblout),
            "-E", str(evalue),
            "--cpu", str(cpu),
            str(hmm_db),
            str(fasta_path),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        if proc.returncode != 0:
            raise HmmSearchError(f"hmmsearch failed: {proc.stderr[-2000:]}")

        return parse_domtblout(domtblout, metadata)


def parse_domtblout(path: Path, metadata: Metadata) -> list[dict]:
    hits = []
    with path.open() as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 22:
                continue
            # HMM profiles are named "<protein id>.aln" after the
            # alignment file they were built from.
            adp = _PROFILE_SUFFIX_RE.sub("", fields[3])
            hits.append({
                "query_name": fields[0],
                "adp": adp,
                "moa": metadata.moa(adp),
                "defence": metadata.defence(adp),
                "evalue": float(fields[6]),
                "score": float(fields[7]),
                "hmm_from": int(fields[15]),
                "hmm_to": int(fields[16]),
                "ali_from": int(fields[17]),
                "ali_to": int(fields[18]),
            })
    hits.sort(key=lambda h: h["evalue"])
    return hits
