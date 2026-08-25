from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from . import __version__, db, foldseek, hmm, output
from .metadata import Metadata


def _add_db_override_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--cache-dir", type=Path, default=None,
        help="Where the local database lives (default: ~/.cache/evades-search, "
             "or $XDG_CACHE_HOME/evades-search).",
    )
    parser.add_argument(
        "--hmm-db", type=Path, default=None,
        help="Override path to the hmmpress'd HMM library.",
    )
    parser.add_argument(
        "--foldseek-db", type=Path, default=None,
        help="Override path to the Foldseek database.",
    )
    parser.add_argument(
        "--metadata", type=Path, default=None,
        help="Override path to metadata.tsv.",
    )


def _resolve_paths(args: argparse.Namespace) -> db.Paths:
    cache_dir = args.cache_dir or db.default_cache_dir()
    paths = db.paths_for(cache_dir)
    if args.hmm_db:
        paths = db.Paths(**{**paths.__dict__, "hmm_db": args.hmm_db})
    if args.foldseek_db:
        paths = db.Paths(**{**paths.__dict__, "foldseek_db": args.foldseek_db})
    if args.metadata:
        paths = db.Paths(**{**paths.__dict__, "metadata_tsv": args.metadata})
    return paths


def cmd_fetch_db(args: argparse.Namespace) -> int:
    try:
        db.fetch(
            base_url=args.base_url,
            cache_dir=args.cache_dir or db.default_cache_dir(),
            force=args.force,
        )
    except db.FetchError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


def cmd_hmm(args: argparse.Namespace) -> int:
    paths = _resolve_paths(args)
    try:
        metadata = Metadata.load(paths.metadata_tsv)
        hits = hmm.run_hmmsearch(
            args.fasta, paths.hmm_db, metadata,
            evalue=args.evalue, cpu=args.cpu,
        )
    except (hmm.HmmSearchError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    with _open_out(args.output) as out:
        output.write(hits, output.HMM_COLUMNS, args.format, out)
    return 0


def cmd_structure(args: argparse.Namespace) -> int:
    paths = _resolve_paths(args)
    try:
        metadata = Metadata.load(paths.metadata_tsv)
        hits = foldseek.run_foldseek(
            args.structure, paths.foldseek_db, metadata,
            tm_score_min=args.tm_score_min, threads=args.threads,
        )
    except (foldseek.FoldseekError, FileNotFoundError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    with _open_out(args.output) as out:
        output.write(hits, output.FOLDSEEK_COLUMNS, args.format, out)
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    paths = _resolve_paths(args)
    def status(found: bool) -> str:
        return "[found]" if found else "[missing]"

    print(f"cache dir:    {paths.cache_dir}")
    print(f"HMM db:       {paths.hmm_db}  {status(paths.hmm_db.exists())}")
    print(f"foldseek db:  {paths.foldseek_db}  {status(Path(str(paths.foldseek_db)).exists())}")
    print(f"metadata.tsv: {paths.metadata_tsv}  {status(paths.metadata_tsv.exists())}")
    print()
    print(f"hmmsearch on PATH: {'yes' if hmm.hmmsearch_available() else 'no'}")
    print(f"foldseek on PATH:  {'yes' if foldseek.foldseek_available() else 'no'}")
    if paths.metadata_tsv.exists():
        try:
            n = len(Metadata.load(paths.metadata_tsv))
            print(f"proteins in metadata.tsv: {n}")
        except Exception:
            pass
    return 0


class _StdoutCtx:
    def __enter__(self):
        return sys.stdout

    def __exit__(self, *exc):
        return False


def _open_out(path: Path | None):
    if path is None:
        return _StdoutCtx()
    return path.open("w")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="evades-search",
        description="Run EVADES anti-defence protein searches (HMM + structural) locally — "
                     "the same searches as the EVADES web app's Analyse tab.",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    default_base_url = os.environ.get("EVADES_SEARCH_BASE_URL", db.DEFAULT_BASE_URL)
    p_fetch = sub.add_parser(
        "fetch-db", help="Download and build the local EVADES search databases.",
    )
    p_fetch.add_argument(
        "--base-url", default=default_base_url,
        help=(
            f"Base URL to fetch bulk-download files from (default: {default_base_url}; "
            "also settable via $EVADES_SEARCH_BASE_URL)."
        ),
    )
    p_fetch.add_argument(
        "--cache-dir", type=Path, default=None,
        help="Where to store the built database (default: ~/.cache/evades-search).",
    )
    p_fetch.add_argument(
        "--force", action="store_true",
        help="Re-download and rebuild even if already present.",
    )
    p_fetch.set_defaults(func=cmd_fetch_db)

    p_hmm = sub.add_parser(
        "hmm", help="Search a protein FASTA against the EVADES HMM profile library.",
    )
    p_hmm.add_argument("fasta", type=Path, help="Multi-FASTA of protein sequences to search.")
    p_hmm.add_argument(
        "-E", "--evalue", type=float, default=hmm.DEFAULT_EVALUE,
        help=f"E-value cutoff (default: {hmm.DEFAULT_EVALUE}).",
    )
    p_hmm.add_argument("--cpu", type=int, default=2, help="Threads for hmmsearch (default: 2).")
    p_hmm.add_argument("-f", "--format", choices=["table", "tsv", "json"], default="table")
    p_hmm.add_argument(
        "-o", "--output", type=Path, default=None, help="Write to a file instead of stdout.",
    )
    _add_db_override_args(p_hmm)
    p_hmm.set_defaults(func=cmd_hmm)

    p_struct = sub.add_parser(
        "structure", help="Search a predicted structure against the EVADES Foldseek database.",
    )
    p_struct.add_argument("structure", type=Path, help="Query structure file (.pdb or .cif).")
    p_struct.add_argument(
        "--tm-score-min", type=float, default=foldseek.DEFAULT_TM_SCORE_MIN,
        help=f"Minimum TM-score to keep a hit (default: {foldseek.DEFAULT_TM_SCORE_MIN}).",
    )
    p_struct.add_argument(
        "--threads", type=int, default=2, help="Threads for foldseek (default: 2).",
    )
    p_struct.add_argument("-f", "--format", choices=["table", "tsv", "json"], default="table")
    p_struct.add_argument(
        "-o", "--output", type=Path, default=None, help="Write to a file instead of stdout.",
    )
    _add_db_override_args(p_struct)
    p_struct.set_defaults(func=cmd_structure)

    p_info = sub.add_parser("info", help="Show database location, presence, and tool availability.")
    _add_db_override_args(p_info)
    p_info.set_defaults(func=cmd_info)

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
