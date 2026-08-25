from evades_search.metadata import Metadata


def test_load_reads_defence_and_moa(tmp_path):
    path = tmp_path / "metadata.tsv"
    path.write_text("ID\tDefences\tMoA\npnk\tPrrC\tphosphorylates tRNA\nno_data\t_\t_\n")

    metadata = Metadata.load(path)

    assert metadata.defence("pnk") == "PrrC"
    assert metadata.moa("pnk") == "phosphorylates tRNA"
    assert metadata.defence("no_data") == ""
    assert metadata.defence("unknown_id") == ""
    assert metadata.known_ids() == frozenset({"pnk", "no_data"})


def test_load_missing_file_raises(tmp_path):
    try:
        Metadata.load(tmp_path / "does_not_exist.tsv")
        assert False, "expected FileNotFoundError"
    except FileNotFoundError:
        pass
