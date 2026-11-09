"""Tests for what the flux page's footer says about its data: which FEMM, and which day (flux_explorer_build.provenance).

The footer used to say "FEMM 4.2 solves ... 2026-10-28" whatever the data was. These tests pin down that it now comes from
facts: the date of the sweeps' own tables, and the version the installed FEMM states about itself, named only when the data
cannot be older than that FEMM. No FEMM is started. Run from P:   uv run pytest tests/test_flux_explorer_build.py -q
"""
import datetime
import os
import re
import time
from pathlib import Path

import pytest

import femm_session
import flux_explorer_build as build

DAY = 24 * 3600


def make_tables(folder: Path, written: float):
    """The two tables provenance() looks at, with a chosen modification time (seconds since 1970)."""
    for name in ("states.csv", "states_p120.csv"):
        (folder / name).write_text("gid\n0\n")
        os.utime(folder / name, (written, written))


def fake_femm(monkeypatch, tmp_path, version, built):
    """Make femm_session say: this FEMM has this version, and its program file was made at time `built`."""
    program = tmp_path / "femm_headlessw.exe"
    program.write_bytes(b"not a real program")
    os.utime(program, (built, built))
    monkeypatch.setattr(femm_session, "femm_version", lambda: version)
    monkeypatch.setattr(femm_session, "femm_program", lambda: program)


def test_the_footer_names_the_installed_femm_and_the_day_the_sweeps_finished(tmp_path, monkeypatch):
    now = time.time()
    make_tables(tmp_path, now)
    fake_femm(monkeypatch, tmp_path, "5.0.20261101", now - 3 * DAY)
    engine, day = build.provenance(str(tmp_path))
    assert engine == "FEMM 5.0.20261101"
    assert day == datetime.date.fromtimestamp(now).isoformat()


def test_data_older_than_the_installed_femm_does_not_get_that_femms_name(tmp_path, monkeypatch, capsys):
    now = time.time()
    make_tables(tmp_path, now - 10 * DAY)                        # the tables predate this FEMM's program file
    fake_femm(monkeypatch, tmp_path, "5.0.20261101", now - 3 * DAY)
    engine, day = build.provenance(str(tmp_path))
    assert engine == "FEMM" and "5.0" not in engine              # an earlier FEMM made them: we do not know which
    assert day == datetime.date.fromtimestamp(now - 10 * DAY).isoformat()
    assert "older than the installed FEMM" in capsys.readouterr().out    # and the builder says so on the console


def test_one_old_table_is_enough_to_withhold_the_version(tmp_path, monkeypatch):
    now = time.time()
    make_tables(tmp_path, now)
    os.utime(tmp_path / "states.csv", (now - 10 * DAY, now - 10 * DAY))     # the 120 deg table is new, the main one old
    fake_femm(monkeypatch, tmp_path, "5.0.20261101", now - 3 * DAY)
    assert build.provenance(str(tmp_path))[0] == "FEMM"


def test_no_femm_installed_means_no_version_but_still_a_day(tmp_path, monkeypatch):
    make_tables(tmp_path, time.time())

    def not_installed():
        raise FileNotFoundError("FEMM_DIR is not set")
    monkeypatch.setattr(femm_session, "femm_version", not_installed)
    assert build.provenance(str(tmp_path))[0] == "FEMM"


def test_a_program_that_states_no_version_gives_no_version(tmp_path, monkeypatch):
    now = time.time()
    make_tables(tmp_path, now)
    fake_femm(monkeypatch, tmp_path, None, now - 3 * DAY)
    assert build.provenance(str(tmp_path))[0] == "FEMM"


def test_the_footer_text_has_nothing_typed_in_about_the_engine_or_the_day():
    """The old footer said FEMM 4.2 and 2026-10-28 in the source. Neither may be written in the builder again."""
    source = Path(build.__file__).read_text(encoding="utf-8")
    assert "FEMM 4.2" not in source
    footer = source[source.index('"footer": ['):]                # (a date elsewhere, in a comment about a test, is fine)
    assert not re.search(r"20\d\d-\d\d-\d\d", footer)
    assert "{engine}" in footer and "{day}" in footer and "{n_solves" in footer


def test_the_real_data_on_this_pc_gets_the_real_femm(tmp_path):
    """On this PC the sweeps' tables are newer than FEMM 5.0 (installed 2026-11-01): the footer names it."""
    out = Path(build.OUT)
    if not (out / "states.csv").is_file():
        pytest.skip("no sweep data in out/flux_explorer")
    try:
        version = femm_session.femm_version()
    except FileNotFoundError:
        pytest.skip("FEMM is not installed")
    engine, day = build.provenance(str(out))
    assert re.fullmatch(r"\d{4}-\d\d-\d\d", day)
    assert engine in ("FEMM", f"FEMM {version}")
