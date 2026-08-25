import tarfile

from evades_search import db


def _make_archive_with_junk(tmp_path):
    archive = tmp_path / "structures.tar.gz"
    real = tmp_path / "protA.pdb"
    real.write_text("ATOM ...")
    junk = tmp_path / "._protA.pdb"
    junk.write_bytes(b"\x00\x05\x16\x07AppleDouble junk")
    ds_store = tmp_path / ".DS_Store"
    ds_store.write_bytes(b"junk")

    with tarfile.open(archive, "w:gz") as tf:
        tf.add(real, arcname="EVADES_v1/protA.pdb")
        tf.add(junk, arcname="EVADES_v1/._protA.pdb")
        tf.add(ds_store, arcname="EVADES_v1/.DS_Store")
    return archive


def test_prune_junk_files_removes_appledouble_and_ds_store(tmp_path):
    archive = _make_archive_with_junk(tmp_path)
    dest = tmp_path / "extracted"
    db._safe_extract(archive, dest)

    # Sanity check: tarfile really does extract the AppleDouble member —
    # this is the bug being guarded against (macOS's `tar tzf` hides these
    # entries from a listing, but Python's tarfile extracts them anyway).
    assert (dest / "EVADES_v1" / "._protA.pdb").exists()

    db._prune_junk_files(dest)

    remaining = sorted(p.name for p in dest.rglob("*") if p.is_file())
    assert remaining == ["protA.pdb"]
