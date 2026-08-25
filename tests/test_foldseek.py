from pathlib import Path

from evades_search.foldseek import _base_protein_name, _parse_results
from evades_search.metadata import Metadata


def _metadata(tmp_path, ids):
    path = tmp_path / "metadata.tsv"
    body = "ID\tDefences\tMoA\n" + "".join(f"{i}\t_\t_\n" for i in ids)
    path.write_text(body)
    return Metadata.load(path)


def _stems(*names):
    return {name: Path(f"{name}.pdb") for name in names}


def test_base_protein_name_strips_nmr_model_suffix():
    known = frozenset({"klca"})
    assert _base_protein_name("klca_MODEL_16_A", known) == "klca"
    assert _base_protein_name("klca_MODEL_1_B", known) == "klca"


def test_base_protein_name_leaves_plain_names_unchanged():
    assert _base_protein_name("some_protein", frozenset()) == "some_protein"


def test_base_protein_name_collapses_multimer_chains():
    # Protein IDs that themselves contain underscores (or even a
    # trailing one) — a blind "strip the last _segment" would mangle
    # these, so the known-ID set has to be consulted.
    known = frozenset({"acrib4", "gp5_9", "acrvia2_plus", "acrif18_"})

    assert _base_protein_name("acrib4_A", known) == "acrib4"
    assert _base_protein_name("acrib4_B", known) == "acrib4"
    assert _base_protein_name("gp5_9_B", known) == "gp5_9"
    assert _base_protein_name("acrvia2_plus_A", known) == "acrvia2_plus"
    assert _base_protein_name("acrif18__A", known) == "acrif18_"
    assert _base_protein_name("acrib4_MODEL_2_A", known) == "acrib4"
    assert _base_protein_name("acrib4", known) == "acrib4"


def test_parse_results_filters_by_tm_score(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "query\tprotA\t0.5\t100\t0.9\t0.7\n"
        "query\tprotB\t0.4\t90\t0.6\t0.3\n"
    )
    metadata = _metadata(tmp_path, ["protA", "protB"])

    hits = _parse_results(results, metadata, tm_score_min=0.5, stem_to_source=_stems("query"))

    assert [h["adp"] for h in hits] == ["protA"]


def test_parse_results_collapses_multi_chain_keeping_best_prob(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "query\tacrib4_A\t0.5\t100\t0.7\t0.8\n"
        "query\tacrib4_B\t0.5\t100\t0.9\t0.8\n"
    )
    metadata = _metadata(tmp_path, ["acrib4"])

    hits = _parse_results(results, metadata, tm_score_min=0.5, stem_to_source=_stems("query"))

    assert len(hits) == 1
    assert hits[0]["adp"] == "acrib4"
    assert hits[0]["prob"] == 0.9  # kept the higher-confidence chain


def test_parse_results_sorted_by_prob_descending(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text(
        "query\tprotA\t0.5\t100\t0.6\t0.9\n"
        "query\tprotB\t0.5\t100\t0.9\t0.9\n"
    )
    metadata = _metadata(tmp_path, ["protA", "protB"])

    hits = _parse_results(results, metadata, tm_score_min=0.5, stem_to_source=_stems("query"))

    assert [h["adp"] for h in hits] == ["protB", "protA"]


def test_parse_results_clamps_tm_score_to_one(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text("query\tprotA\t0.5\t100\t0.9\t1.0000003\n")
    metadata = _metadata(tmp_path, ["protA"])

    hits = _parse_results(results, metadata, tm_score_min=0.5, stem_to_source=_stems("query"))

    assert hits[0]["tm_score"] == 1.0


def test_parse_results_keeps_same_protein_hit_separate_per_query_file(tmp_path):
    # Two different query structures (batch search) both matching the
    # same EVADES protein must NOT clobber each other — each query's
    # best hit against that protein should survive independently.
    results = tmp_path / "results.tsv"
    results.write_text(
        "0000_fileA\tprotA\t0.5\t100\t0.9\t0.8\n"
        "0001_fileB\tprotA\t0.5\t100\t0.6\t0.8\n"
    )
    metadata = _metadata(tmp_path, ["protA"])
    stems = {"0000_fileA": Path("fileA.pdb"), "0001_fileB": Path("fileB.pdb")}

    hits = _parse_results(results, metadata, tm_score_min=0.5, stem_to_source=stems)

    assert len(hits) == 2
    assert {h["query_file"] for h in hits} == {"fileA.pdb", "fileB.pdb"}
    assert all(h["adp"] == "protA" for h in hits)


def test_parse_results_resolves_query_file_and_chain_suffix(tmp_path):
    results = tmp_path / "results.tsv"
    results.write_text("0000_gp5_9\tprotA\t0.5\t100\t0.9\t0.8\n")
    metadata = _metadata(tmp_path, ["protA"])
    stems = {"0000_gp5_9": Path("/some/dir/gp5_9.pdb")}

    hits = _parse_results(results, metadata, tm_score_min=0.5, stem_to_source=stems)

    assert hits[0]["query_file"] == "/some/dir/gp5_9.pdb"
    assert hits[0]["query"] == "gp5_9"  # original stem, no numeric staging prefix leaked
