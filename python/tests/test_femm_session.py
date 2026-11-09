"""Tests for femm_session: one private FEMM 5.0 per Python process.

The first tests start real FEMMs (about 0.3 s each), so they are marked `femm`. The last ones start none: they read a
program file's version, check how FEMM_DIR is found, or probe a port with a listener made for the test. The files that
FEMM writes go to P\\out\\tests\\ (an ordinary folder of P, which git ignores), never to a temp folder.
Run from P:   uv run pytest tests/test_femm_session.py -q

The shell this was written in started before FEMM 5.0 was installed, so it has no FEMM_DIR: the tests are
meant to pass there too, because femm_session takes FEMM's folder from config.toml first (the config names it), and
reads FEMM_DIR from the registry when the config names none and the shell lacks it.
"""
import multiprocessing
import os
import socket
import subprocess
import time
from pathlib import Path

import pytest

import femm_session as fs
import femm_test_tools as tools
import paths

OUT = Path(__file__).resolve().parent.parent / "out" / "tests"
DEFAULT_PORT = 48015                       # where openfemm() always goes, and where an interactive FEMM window listens
BELOW_NORMAL = 0x4000                      # BELOW_NORMAL_PRIORITY_CLASS
NORMAL = 0x20                              # NORMAL_PRIORITY_CLASS


def port_in_use(number):
    """True if something holds this port. Tries to bind it; it never connects.

    Not connecting matters: a FEMM that has no client yet takes the first connection as its client, and a
    headless FEMM quits when its client hangs up. A connect-and-leave probe would quit a FEMM that is not ours.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        try:
            s.bind(("127.0.0.1", number))
        except OSError:
            return True
    return False


def femm_is_running(pid):
    info = tools.process_info(pid)
    return info is not None and info[0].lower().startswith("femm_headless")


@pytest.mark.femm
def test_femm_dir_found_without_env(monkeypatch):
    monkeypatch.delenv("FEMM_DIR", raising=False)
    assert (fs.femm_dir() / "femm_headlessw.exe").is_file()


@pytest.mark.femm
def test_free_port_is_free():
    number = fs.free_port()
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", number))                   # raises OSError if somebody else has the port
    assert number != DEFAULT_PORT


@pytest.mark.femm
def test_open_femm_starts_its_own_femm():
    session = fs.open_femm()
    try:
        child = fs.pid(session)
        assert fs.port(session) != DEFAULT_PORT
        assert child is not None and tools.process_info(child)[0].lower() == "femm_headlessw.exe"
    finally:
        fs.close_femm()
    deadline = time.monotonic() + 5
    while femm_is_running(child) and time.monotonic() < deadline:
        time.sleep(0.05)
    assert not femm_is_running(child), "the FEMM was still running 5 s after close_femm()"


@pytest.mark.femm
def test_open_femm_starts_the_femm_of_the_config_even_if_femm_dir_names_another_folder(monkeypatch, tmp_path):
    """R4: femm.FemmSession looks in the environment's FEMM_DIR first, so open_femm() has to put femm_dir()'s answer there.

    Without that, a FEMM_DIR left over from an earlier install would decide which FEMM starts, and config.toml would only
    decide which one femm_dir() reports. tmp_path is a folder with no FEMM in it; afterwards FEMM_DIR names the config's.
    """
    monkeypatch.setenv("FEMM_DIR", str(tmp_path))
    session = fs.open_femm()
    try:
        assert os.environ["FEMM_DIR"] == str(paths.tool("femm"))
        child = fs.pid(session)
        assert child is not None and tools.process_info(child)[0].lower() == "femm_headlessw.exe"
    finally:
        fs.close_femm()


@pytest.mark.femm
def test_never_joins_an_instance_on_48015():
    standin = None
    if not port_in_use(DEFAULT_PORT):      # nothing there, so start a stand-in for an interactive FEMM window
        exe = fs.femm_dir() / "femm_headlessw.exe"
        standin = subprocess.Popen([str(exe), "--port", str(DEFAULT_PORT)],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # else: something already holds 48015, perhaps an interactive FEMM window. It is the stand-in. We only ever try
    # to bind the port to see that it is held; we never connect, send a command, close or kill it.
    try:
        deadline = time.monotonic() + 10
        while not port_in_use(DEFAULT_PORT) and time.monotonic() < deadline:
            time.sleep(0.1)
        assert port_in_use(DEFAULT_PORT), "no stand-in on 48015"
        session = fs.open_femm()
        try:
            child = fs.pid(session)
            assert fs.port(session) != DEFAULT_PORT
            assert child is not None                                  # None would mean "joined an existing FEMM"
            if standin is not None:
                assert child != standin.pid
            exe_name, parent = tools.process_info(child)
            assert exe_name.lower() == "femm_headlessw.exe"
            assert parent == os.getpid(), "pid(session) is a FEMM this test did not start"
        finally:
            fs.close_femm()
        assert port_in_use(DEFAULT_PORT), "close_femm() quit the FEMM on 48015, which was not ours"
        if standin is not None:
            assert standin.poll() is None, "the stand-in on 48015 was quit"
    finally:
        if standin is not None:
            standin.terminate()
            standin.wait(timeout=10)


@pytest.mark.femm
def test_two_processes_get_two_femms():
    folder = OUT / "two"
    ctx = multiprocessing.get_context("spawn")
    both_alive = ctx.Barrier(2)
    results = ctx.Queue()
    jobs = [(folder / "with_magnet.fem", 1_030_272), (folder / "no_magnet.fem", 0)]   # private file per worker
    workers = [ctx.Process(target=tools.magnet_worker, args=(path, hc, both_alive, results)) for path, hc in jobs]
    for w in workers:
        w.start()
    got = [results.get(timeout=120) for _ in workers]
    for w in workers:
        w.join(timeout=60)
    assert [r[0] for r in got] == ["ok", "ok"], got
    by_hc = {r[1]: r for r in got}
    with_magnet, no_magnet = by_hc[1_030_272], by_hc[0]
    assert with_magnet[2] != no_magnet[2]          # the ports differ
    assert with_magnet[3] != no_magnet[3]          # the FEMM processes differ
    assert with_magnet[5] != no_magnet[5]          # and so do the two Python workers that own them
    assert with_magnet[4] > 0.1                    # By [T] with the magnet
    assert abs(no_magnet[4]) < 1e-6                # and without it: the worker did not read the other's answer
    assert all(w.exitcode == 0 for w in workers)


def test_a_block_that_lowers_the_pytest_process_leaves_its_priority_class_as_it_found_it():
    """B6 of the final review. femm_session.lower_priority() (open_femm calls it) lowers the WHOLE process, so a test that opens a FEMM inside the pytest process left every test
    after it at below normal, timing tests included. tools.priority_class_restored() is the block the shared fixture of conftest.py wraps every test in: whatever happens inside,
    the class the process had when the block began is put back (and a block that lowers nothing leaves it alone)."""
    pid, original = os.getpid(), tools.priority_class(os.getpid())
    try:
        tools.set_own_priority_class(NORMAL)
        with tools.priority_class_restored():
            fs.lower_priority()
            assert tools.priority_class(pid) == BELOW_NORMAL
        assert tools.priority_class(pid) == NORMAL
        with tools.priority_class_restored():
            pass
        assert tools.priority_class(pid) == NORMAL
    finally:
        tools.set_own_priority_class(original)


_lowered_by_the_test_before = {}                    # what the pair of tests below tell each other: the class the fixture saw, then the class the next test finds


def test_a_test_that_lowers_the_process_is_followed_by_one_that_finds_it_as_it_was():
    """Part 1 of the pair that shows the shared fixture (conftest.py, autouse) at work: this test lowers the whole pytest process, as every test that opens a FEMM does."""
    _lowered_by_the_test_before["class"] = tools.priority_class(os.getpid())            # what the fixture saw before this test (and will put back)
    tools.set_own_priority_class(NORMAL)
    fs.lower_priority()
    assert tools.priority_class(os.getpid()) == BELOW_NORMAL


def test_the_next_test_finds_the_process_at_the_class_it_had_before_the_one_that_lowered_it():
    """Part 2: it runs right after part 1 (pytest keeps a file's order), and finds the class the process had before part 1 lowered it. It is skipped when run alone."""
    if "class" not in _lowered_by_the_test_before:
        pytest.skip("part 1 did not run before it (the two are meant to run one after the other)")
    assert tools.priority_class(os.getpid()) == _lowered_by_the_test_before["class"]


@pytest.mark.femm
def test_femm_runs_below_normal():
    # Start from normal, and show first that without the request FEMM really is normal. (Earlier tests that lowered this process no longer
    # make the check below pass by itself: the shared fixture of conftest.py puts the class back after every test.)
    tools.set_own_priority_class(NORMAL)
    session = fs.open_femm(below_normal=False)
    try:
        assert tools.priority_class(fs.pid(session)) == NORMAL
    finally:
        fs.close_femm()
    session = fs.open_femm(below_normal=True)
    try:
        assert tools.priority_class(fs.pid(session)) == BELOW_NORMAL
    finally:
        fs.close_femm()


@pytest.mark.femm
def test_paths_with_spaces():
    fem = OUT / "with space" / "magnet.fem"
    fs.open_femm()
    try:
        by = tools.solve_magnet(fem, 1_030_272)
    finally:
        fs.close_femm()
    assert by > 0.1
    assert fem.with_suffix(".ans").is_file()


# ------------------------------------------------------------------------- the version a program file states
def test_a_file_without_a_version_resource_states_no_version(tmp_path):
    """file_version() reads the text block inside an .exe. A text file, or a file that is not there, has none."""
    text = tmp_path / "notes.exe"
    text.write_text("this is not a program")
    assert fs.file_version(text) is None
    assert fs.file_version(tmp_path / "missing.exe") is None


def test_the_version_calls_are_declared_once_and_as_pointers():
    """argtypes / restype are set when version.dll is first used and kept (femm_session._version_calls): a pointer is a c_void_p."""
    import ctypes
    dll = fs._version_calls()
    assert dll is fs._version_calls()                                       # the same object: nothing is declared again per call
    assert dll.VerQueryValueW.argtypes[0] is ctypes.c_void_p and dll.GetFileVersionInfoW.argtypes[3] is ctypes.c_void_p
    assert fs.file_version(fs.femm_program()) == fs.file_version(fs.femm_program())          # and it still reads the same string


def test_the_installed_femm_states_its_version():
    """FEMM 5.0's programs carry "5.0.<date>" as their file version. The flux page's footer names it.

    Not marked `femm`: it starts no FEMM, it only reads the version text inside the installed program files.
    """
    version = fs.femm_version()
    assert version is not None and version.startswith("5.0."), version
    assert fs.file_version(fs.femm_dir() / "fkn.exe") == version            # the solver says the same
    assert fs.femm_program().name in ("femm_headlessw.exe", "femm_headless.exe", "femm_headless")


# ------------------------------------------------------------------------- where FEMM is, and which ports count as taken
def test_a_bogus_femm_dir_raises_file_not_found_and_names_femm_dir(monkeypatch, tmp_path):
    """femm_dir() refuses a FEMM_DIR that does not lead to FEMM, and the message names FEMM_DIR (no FEMM is started).

    open_femm() calls femm_dir() first, so a wrong setting is reported before anything is started or lowered. Three
    wrong settings: a folder that does not exist, a folder that exists but holds none of FEMM's programs (tmp_path has
    only the config in it), and no setting at all (the lookup of the environment and the registry is replaced by one
    that finds none).

    Since the tools plan's R4 the config comes first and FEMM_DIR is only the fall-back, so this test hands femm_dir() a
    config that names no FEMM. (tests/test_paths.py has the other half: a config that does name one wins over FEMM_DIR.)
    """
    (tmp_path / "config.toml").write_text("[tools]\nfemm = ''\n", encoding="utf-8")
    monkeypatch.setattr(paths, "CONFIG", tmp_path / "config.toml")
    for bogus in (tmp_path / "no such folder", tmp_path):
        monkeypatch.setenv("FEMM_DIR", str(bogus))
        with pytest.raises(FileNotFoundError, match="FEMM_DIR"):
            fs.femm_dir()
    monkeypatch.setattr(fs, "_femm_dir_setting", lambda: (None, None))
    with pytest.raises(FileNotFoundError, match="FEMM_DIR is not set"):
        fs.femm_dir()


@pytest.mark.skipif(not hasattr(socket, "SO_EXCLUSIVEADDRUSE"), reason="SO_EXCLUSIVEADDRUSE is a Windows socket option")
@pytest.mark.parametrize("exclusive", [True, False])
def test_port_taken_sees_a_listener_on_127_0_0_1_exclusive_or_not(exclusive):
    """_port_taken() sees a listener on 127.0.0.1 whether or not it holds the number exclusively: FEMM 5.0 listens there.

    This pins the first sentence of "which owners this probe sees" in _port_taken()'s docstring. The rest of that
    paragraph (an exclusive listener on 0.0.0.0 is seen; a plain one on 0.0.0.0 and an IPv6 one are not) was checked by
    hand on 2026-11-05 and is not tested here: a listener on every address can make Windows' firewall ask a question,
    and a listener on the loopback address never does.

    "Exclusive" is SO_EXCLUSIVEADDRUSE: the owner claims the number for itself alone. The listener is made with port 0
    (any free number), so no number is guessed, and it never accepts a connection: the probe only tries to bind.
    """
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as owner:
        if exclusive:
            owner.setsockopt(socket.SOL_SOCKET, socket.SO_EXCLUSIVEADDRUSE, 1)
        owner.bind(("127.0.0.1", 0))
        owner.listen(1)
        number = owner.getsockname()[1]
        assert fs._port_taken(number), f"a listener on 127.0.0.1:{number} (exclusive: {exclusive}) was not seen"
    assert not fs._port_taken(number), "the number is free again once its owner has closed the socket"
