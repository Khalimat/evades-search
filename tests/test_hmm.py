from evades_search.hmm import parse_domtblout
from evades_search.metadata import Metadata

# domtblout columns used by parse_domtblout (0-indexed):
#   0 target(seq) name, 3 query(profile) name, 6 e-value, 7 score,
#   15 hmm-from, 16 hmm-to, 17 ali-from, 18 ali-to
DOMTBLOUT_HEADER_COMMENT = "# this is a comment line and should be ignored\n"


def _line(name, profile, evalue, score, hmm_from, hmm_to, ali_from, ali_to):
    fields = [name, "-", "100", profile, "-", "50", str(evalue), str(score), "0.1"]
    fields += ["1", "1", str(evalue), str(evalue), str(score), "0.1"]
    fields += [str(hmm_from), str(hmm_to), str(ali_from), str(ali_to)]
    fields += ["1", "60", "0.9", "some description"]
    return " ".join(fields) + "\n"


def test_parse_domtblout_sorts_by_evalue_ascending(tmp_path):
    domtblout = tmp_path / "hits.domtblout"
    domtblout.write_text(
        DOMTBLOUT_HEADER_COMMENT
        + _line("query1", "profileA", "1e-10", "55.2", 10, 60, 5, 55)
        + _line("query2", "profileB", "1e-20", "80.1", 20, 70, 15, 65)
    )

    hits = parse_domtblout(domtblout, Metadata.empty())

    assert len(hits) == 2
    assert hits[0]["query_name"] == "query2"  # more significant e-value comes first
    assert hits[0]["evalue"] == 1e-20
    assert hits[1]["query_name"] == "query1"
    assert hits[1]["adp"] == "profileA"
    assert hits[1]["hmm_from"] == 10
    assert hits[1]["ali_to"] == 55


def test_parse_domtblout_skips_short_and_comment_lines(tmp_path):
    domtblout = tmp_path / "hits.domtblout"
    domtblout.write_text(
        "# comment\n"
        "too short a line\n"
        + _line("query1", "profileA", "1e-10", "55.2", 10, 60, 5, 55)
    )

    hits = parse_domtblout(domtblout, Metadata.empty())

    assert len(hits) == 1
    assert hits[0]["query_name"] == "query1"


def test_parse_domtblout_strips_aln_suffix_from_profile_name(tmp_path):
    domtblout = tmp_path / "hits.domtblout"
    domtblout.write_text(
        DOMTBLOUT_HEADER_COMMENT
        + _line("query1", "pnk.aln", "1e-10", "55.2", 10, 60, 5, 55)
    )

    hits = parse_domtblout(domtblout, Metadata.empty())

    assert hits[0]["adp"] == "pnk"


def test_parse_domtblout_joins_metadata(tmp_path):
    metadata = Metadata.load(_write_metadata(tmp_path, "ID\tDefences\tMoA\npnk\tPrrC\tsome moa\n"))
    domtblout = tmp_path / "hits.domtblout"
    domtblout.write_text(_line("query1", "pnk.aln", "1e-10", "55.2", 10, 60, 5, 55))

    hits = parse_domtblout(domtblout, metadata)

    assert hits[0]["defence"] == "PrrC"
    assert hits[0]["moa"] == "some moa"


def _write_metadata(tmp_path, content):
    path = tmp_path / "metadata.tsv"
    path.write_text(content)
    return path
