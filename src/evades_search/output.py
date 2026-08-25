"""Renders hit lists (from hmm.run_hmmsearch / foldseek.run_foldseek) as
a human-readable table, TSV, or JSON.
"""
from __future__ import annotations

import csv
import json
import sys
from typing import IO, Sequence

HMM_COLUMNS = [
    "query_name", "adp", "moa", "defence", "evalue", "score",
    "hmm_from", "hmm_to", "ali_from", "ali_to",
]
FOLDSEEK_COLUMNS = [
    "query", "query_file", "adp", "moa", "defence", "seq_identity", "aln_len", "prob", "tm_score",
]


def write(hits: list[dict], columns: Sequence[str], fmt: str, out: IO[str] = sys.stdout) -> None:
    if fmt == "json":
        json.dump(hits, out, indent=2)
        out.write("\n")
    elif fmt == "tsv":
        _write_tsv(hits, columns, out)
    elif fmt == "table":
        _write_table(hits, columns, out)
    else:
        raise ValueError(f"unknown output format: {fmt}")


def _write_tsv(hits: list[dict], columns: Sequence[str], out: IO[str]) -> None:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(columns)
    for hit in hits:
        writer.writerow([hit.get(c, "") for c in columns])


def _write_table(hits: list[dict], columns: Sequence[str], out: IO[str]) -> None:
    if not hits:
        out.write("No hits.\n")
        return

    rows = [[str(hit.get(c, "")) for c in columns] for hit in hits]
    widths = [max(len(columns[i]), *(len(row[i]) for row in rows)) for i in range(len(columns))]

    def fmt_row(cells: Sequence[str]) -> str:
        return "  ".join(cell.ljust(widths[i]) for i, cell in enumerate(cells))

    out.write(fmt_row(columns) + "\n")
    out.write("  ".join("-" * w for w in widths) + "\n")
    for row in rows:
        out.write(fmt_row(row) + "\n")
