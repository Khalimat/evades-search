"""Fetches and builds the local copy of the EVADES search databases from
the same bulk-download files the EVADES website itself serves
(`data/downloads/{hmm_profiles.tar.gz,predicted_structures.tar.gz,metadata.tsv}`
in the evades-webapp repo, published at `<base_url>/`). Caches the built
HMMER and Foldseek indexes under a local cache directory so repeat runs
don't re-download or re-build anything.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.error import URLError

DEFAULT_BASE_URL = "https://167-233-198-65.sslip.io/downloads"
# ^ the live EVADES server's bulk-download directory (see OPERATIONS.md
# in evades-webapp). Override with --base-url / EVADES_SEARCH_BASE_URL
# if the site moves to a different host.

_HMM_ARCHIVE = "hmm_profiles.tar.gz"
_STRUCTURES_ARCHIVE = "predicted_structures.tar.gz"
_METADATA_FILE = "metadata.tsv"


@dataclass(frozen=True)
class Paths:
    cache_dir: Path
    hmm_db: Path            # hmmpress'd HMM library
    foldseek_db: Path       # foldseek createdb output (base name, no suffix)
    metadata_tsv: Path
    raw_dir: Path           # downloaded archives, kept for re-builds
    structures_dir: Path    # extracted .pdb/.cif files


def default_cache_dir() -> Path:
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "evades-search"


def paths_for(cache_dir: Path) -> Paths:
    return Paths(
        cache_dir=cache_dir,
        hmm_db=cache_dir / "hmm" / "evades_profiles.hmm",
        foldseek_db=cache_dir / "foldseek" / "evades_structures_db",
        metadata_tsv=cache_dir / "metadata.tsv",
        raw_dir=cache_dir / "raw",
        structures_dir=cache_dir / "structures",
    )


def is_built(paths: Paths) -> bool:
    return (
        paths.hmm_db.exists()
        and paths.hmm_db.with_suffix(".hmm.h3p").exists()
        and Path(str(paths.foldseek_db)).exists()
        and paths.metadata_tsv.exists()
    )


class FetchError(RuntimeError):
    pass


def _download(url: str, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        with urllib.request.urlopen(url, timeout=60) as resp, dest.open("wb") as f:
            shutil.copyfileobj(resp, f)
    except URLError as exc:
        raise FetchError(f"failed to download {url}: {exc}") from exc


def _prune_junk_files(directory: Path) -> None:
    """Remove macOS extraction artifacts (AppleDouble `._name` resource-fork
    files, `.DS_Store`) that can end up as real tar members when a tar.gz
    was built on macOS with extended attributes/Finder metadata attached —
    `tar tzf` hides them but Python's tarfile extracts them as regular
    files, which `foldseek createdb` would otherwise index as bogus
    structures."""
    for path in directory.rglob("._*"):
        if path.is_file():
            path.unlink()
    for path in directory.rglob(".DS_Store"):
        if path.is_file():
            path.unlink()


def _safe_extract(archive: Path, dest: Path) -> None:
    """Extract a tar.gz, refusing any member that would land outside
    `dest` (defends against a malicious/corrupt archive using `../`
    paths — Python's tarfile doesn't guard against this by default)."""
    dest.mkdir(parents=True, exist_ok=True)
    dest_resolved = dest.resolve()
    with tarfile.open(archive) as tf:
        for member in tf.getmembers():
            member_path = (dest / member.name).resolve()
            inside_dest = member_path == dest_resolved or str(member_path).startswith(
                str(dest_resolved) + os.sep
            )
            if not inside_dest:
                raise FetchError(f"refusing to extract unsafe path from archive: {member.name}")
        tf.extractall(dest)


def fetch(
    *,
    base_url: str = DEFAULT_BASE_URL,
    cache_dir: Path | None = None,
    force: bool = False,
    progress=print,
) -> Paths:
    """Download the bulk EVADES data files and build the local HMMER and
    Foldseek indexes from them. Idempotent — skips the build if it's
    already present unless `force=True`."""
    cache_dir = cache_dir or default_cache_dir()
    paths = paths_for(cache_dir)

    if is_built(paths) and not force:
        progress(f"Database already built at {cache_dir} (use --force to rebuild).")
        return paths

    if not shutil.which("hmmpress"):
        raise FetchError(
            "hmmpress not found on PATH. Install HMMER first "
            "(`brew install hmmer` or `conda install -c bioconda hmmer`)."
        )
    if not shutil.which("foldseek"):
        raise FetchError(
            "foldseek not found on PATH. Install Foldseek first "
            "(`brew install brewsci/bio/foldseek` or `conda install -c bioconda foldseek`)."
        )

    base_url = base_url.rstrip("/")
    paths.raw_dir.mkdir(parents=True, exist_ok=True)

    hmm_archive = paths.raw_dir / _HMM_ARCHIVE
    structures_archive = paths.raw_dir / _STRUCTURES_ARCHIVE

    progress(f"Downloading {_HMM_ARCHIVE} from {base_url} ...")
    _download(f"{base_url}/{_HMM_ARCHIVE}", hmm_archive)

    progress(f"Downloading {_STRUCTURES_ARCHIVE} from {base_url} ...")
    _download(f"{base_url}/{_STRUCTURES_ARCHIVE}", structures_archive)

    progress(f"Downloading {_METADATA_FILE} from {base_url} ...")
    _download(f"{base_url}/{_METADATA_FILE}", paths.metadata_tsv)

    progress("Building HMMER profile database (hmmpress) ...")
    hmm_extract_dir = paths.cache_dir / "_hmm_extract"
    if hmm_extract_dir.exists():
        shutil.rmtree(hmm_extract_dir)
    _safe_extract(hmm_archive, hmm_extract_dir)
    _prune_junk_files(hmm_extract_dir)
    hmm_files = list(hmm_extract_dir.glob("*.hmm"))
    if not hmm_files:
        raise FetchError(f"no .hmm file found inside {_HMM_ARCHIVE}")
    paths.hmm_db.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(hmm_files[0], paths.hmm_db)
    shutil.rmtree(hmm_extract_dir)
    _run_hmmpress(paths.hmm_db)

    progress("Building Foldseek structure database (foldseek createdb) ...")
    if paths.structures_dir.exists():
        shutil.rmtree(paths.structures_dir)
    _safe_extract(structures_archive, paths.structures_dir)
    _prune_junk_files(paths.structures_dir)
    structure_files = [
        p for p in paths.structures_dir.rglob("*") if p.suffix.lower() in (".pdb", ".cif")
    ]
    if not structure_files:
        raise FetchError(f"no .pdb/.cif files found inside {_STRUCTURES_ARCHIVE}")
    paths.foldseek_db.parent.mkdir(parents=True, exist_ok=True)
    _run_foldseek_createdb(paths.structures_dir, paths.foldseek_db)

    progress(f"Done. {len(structure_files)} structures indexed under {cache_dir}")
    return paths


def _run_hmmpress(hmm_db: Path) -> None:
    # -f overwrites any stale index left over from a previous fetch --force
    proc = subprocess.run(
        ["hmmpress", "-f", str(hmm_db)], capture_output=True, text=True
    )
    if proc.returncode != 0:
        raise FetchError(f"hmmpress failed: {proc.stderr[-2000:]}")


def _run_foldseek_createdb(structures_dir: Path, foldseek_db: Path) -> None:
    proc = subprocess.run(
        ["foldseek", "createdb", str(structures_dir), str(foldseek_db)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise FetchError(f"foldseek createdb failed: {proc.stderr[-2000:]}")
