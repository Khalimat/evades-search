"""Loads the curated EVADES protein metadata table (metadata.tsv) used to
annotate search hits with mode-of-action / inhibited-defence columns, and
to resolve Foldseek's chain/model-suffixed target names back to plain
protein IDs. Same source file (and same "_" == empty convention) as the
EVADES web app's `backend/app/tasks.py::_load_metadata`.
"""
from __future__ import annotations

import csv
from pathlib import Path


def _clean(value: str | None) -> str:
    return value if value not in (None, "", "_") else ""


class Metadata:
    """protein ID -> {"defence": ..., "moa": ...} lookup, plus the set of
    known IDs (needed to resolve Foldseek's chain-suffixed target names)."""

    def __init__(self, rows: dict[str, dict[str, str]]):
        self._rows = rows

    @classmethod
    def load(cls, path: Path) -> "Metadata":
        if not path.exists():
            raise FileNotFoundError(
                f"metadata.tsv not found at {path}. Run `evades-search fetch-db` first."
            )
        with path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            rows = {
                row["ID"]: {
                    "defence": _clean(row.get("Defences")),
                    "moa": _clean(row.get("MoA")),
                }
                for row in reader
                if row.get("ID")
            }
        return cls(rows)

    @classmethod
    def empty(cls) -> "Metadata":
        return cls({})

    def defence(self, protein_id: str) -> str:
        return self._rows.get(protein_id, {}).get("defence", "")

    def moa(self, protein_id: str) -> str:
        return self._rows.get(protein_id, {}).get("moa", "")

    def known_ids(self) -> frozenset[str]:
        return frozenset(self._rows)

    def __len__(self) -> int:
        return len(self._rows)
