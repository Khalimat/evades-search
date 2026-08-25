from pathlib import Path

from evades_search.cli import _default_output_path, _resolve_output_path


def test_default_output_path_swaps_extension_for_single_file():
    assert _default_output_path([Path("my_proteins.fasta")], "tsv") == Path("my_proteins.tsv")
    assert _default_output_path([Path("my_protein.pdb")], "json") == Path("my_protein.json")
    assert _default_output_path([Path("my_protein.pdb")], "table") == Path("my_protein.txt")


def test_default_output_path_names_after_directory(tmp_path):
    d = tmp_path / "my_structures"
    d.mkdir()

    result = _default_output_path([d], "tsv")

    assert result == tmp_path / "my_structures.tsv"


def test_default_output_path_falls_back_to_generic_name_for_multiple_inputs():
    result = _default_output_path([Path("a.pdb"), Path("b.cif")], "tsv")

    assert result == Path("evades_search_results.tsv")


def test_resolve_output_path_dash_means_stdout():
    assert _resolve_output_path(Path("-"), [Path("a.pdb")], "tsv") is None


def test_resolve_output_path_explicit_path_wins():
    explicit = Path("custom_name.tsv")

    result = _resolve_output_path(explicit, [Path("a.pdb")], "tsv")

    assert result == explicit


def test_resolve_output_path_auto_derives_when_omitted():
    result = _resolve_output_path(None, [Path("a.pdb")], "json")

    assert result == Path("a.json")
