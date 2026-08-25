import io
import tarfile
from unittest.mock import patch
from urllib.error import HTTPError, URLError

from evades_search import db


def test_default_cache_dir_prefers_evades_search_cache_dir_env(monkeypatch, tmp_path):
    monkeypatch.setenv("EVADES_SEARCH_CACHE_DIR", str(tmp_path / "custom"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))  # lower priority, ignored

    assert db.default_cache_dir() == tmp_path / "custom"


def test_default_cache_dir_falls_back_to_xdg_cache_home(monkeypatch, tmp_path):
    monkeypatch.delenv("EVADES_SEARCH_CACHE_DIR", raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))

    assert db.default_cache_dir() == tmp_path / "xdg" / "evades-search"


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


def _http_error(code):
    return HTTPError("http://example.test/f", code, "error", None, None)


def test_download_retries_transient_server_error_then_succeeds(tmp_path, monkeypatch):
    # Reproduces a Zenodo 504 on the first attempt, succeeding on retry —
    # what actually happened on a user's first real `fetch-db` run.
    monkeypatch.setattr(db, "_DOWNLOAD_BACKOFF_SECONDS", 0)
    dest = tmp_path / "f.tsv"
    responses = [_http_error(504), io.BytesIO(b"content")]

    with patch("evades_search.db.urllib.request.urlopen", side_effect=responses):
        db._download("http://example.test/f", dest)

    assert dest.read_bytes() == b"content"


def test_download_does_not_retry_client_error(tmp_path, monkeypatch):
    monkeypatch.setattr(db, "_DOWNLOAD_BACKOFF_SECONDS", 0)
    dest = tmp_path / "f.tsv"

    with patch(
        "evades_search.db.urllib.request.urlopen", side_effect=_http_error(404)
    ) as mock_urlopen:
        try:
            db._download("http://example.test/f", dest)
            assert False, "expected FetchError"
        except db.FetchError:
            pass

    assert mock_urlopen.call_count == 1  # no retry on a 404 — retrying can't fix it


def test_download_gives_up_after_max_retries(tmp_path, monkeypatch):
    monkeypatch.setattr(db, "_DOWNLOAD_BACKOFF_SECONDS", 0)
    dest = tmp_path / "f.tsv"

    with patch(
        "evades_search.db.urllib.request.urlopen", side_effect=URLError("connection refused")
    ) as mock_urlopen:
        try:
            db._download("http://example.test/f", dest)
            assert False, "expected FetchError"
        except db.FetchError:
            pass

    assert mock_urlopen.call_count == db._DOWNLOAD_RETRIES
