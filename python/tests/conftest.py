"""pytest set-up for the tests.

Three jobs:

1. Put P (the folder above this one), this folder and P\\validation on sys.path, so the tests can `import femm_session`,
   `import femm_test_tools` and `import same_model` however pytest was started. The worker processes the tests spawn
   get a copy of the parent's sys.path, so this is enough for them too.
2. After every test, quit the FEMM that the test opened, if it forgot. A test that fails half way must not leave a
   femm_headlessw.exe running on the machine. And put the pytest process's priority class back: opening a FEMM lowers the whole
   process, and the tests after it would run at below normal.
3. Keep the one copy of what two test files need to see the PATH that Gmsh cuts: `process_path()` reads the PATH of the
   running PROCESS from Windows (Python's os.environ is a separate copy and is never cut), `cut_process_path()` cuts it
   the way a library does, and the `path_is_put_back` fixture restores os.environ["PATH"] whatever a test did.
   tests/test_check_tools.py and the bench's tests (tests/test_bench*.py) import them from here
   (`from conftest import process_path`).
"""
import os
import sys
from pathlib import Path

import pytest

HERE = Path(__file__).resolve().parent
P = HERE.parent
for folder in (str(P), str(HERE), str(P / "validation")):
    if folder not in sys.path:
        sys.path.insert(0, folder)

import paths  # noqa: E402  (config.toml: where out\ is)


def temp_root_under_out() -> Path | None:
    """Make pytest take its temporary folders from out\\.pytest (on D:, git-ignored) and not from the user's Temp on C:, where every run left ~300 files (1.2 MB) and pytest
    keeps the last three runs (149 of the final review). pytest reads PYTEST_DEBUG_TEMPROOT, when it first needs a temporary folder, as the folder in which it makes
    pytest-of-<user>\\pytest-<n>; that folder has to exist. setdefault: a user's own setting wins, and then None is returned; else the root that was set. Not --basetemp: it
    empties its folder at every start and would clobber a second pytest run going at the same time, where the numbered folders of the temp root do not.

    The name starts with a dot so that the stale scan (rulings.results, which stale.py runs on out\\) skips it at any depth, as it skips out\\.femm: the tests' own csv and
    json files are not studies' results, and out\\pytest was scanned and listed about 300 of them as unstamped."""
    if "PYTEST_DEBUG_TEMPROOT" in os.environ:
        return None
    root = paths.path("out") / ".pytest"
    root.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("PYTEST_DEBUG_TEMPROOT", str(root))
    return root


temp_root_under_out()                       # (conftest.py is imported before any test asks for a temporary folder)


@pytest.fixture(autouse=True)
def _quit_leftover_femm():
    yield
    import femm_session
    femm_session.close_femm()


@pytest.fixture(autouse=True)
def _the_priority_class_is_restored():
    """After every test, the pytest process has the priority class it had before it. A test that opens a FEMM here (open_femm lowers the process, and femm_session.lower_priority
    is process-wide), or calls anything else that lowers it, would otherwise leave every test after it at below normal, the timing-based ones of test_femm_queue.py among them
    (B6 of the final review of the tools plan). One shared fixture for all of them: femm_test_tools.priority_class_restored()."""
    import femm_test_tools
    with femm_test_tools.priority_class_restored():
        yield


# ----------------------------------------------------------------------------------------------- the PATH of the process
# Found by the first full run (2026-11-06): after gmsh.initialize() a tool row could not start powershell.exe by name.
# The Python module cuts the PATH of the running PROCESS short (os.environ, Python's own copy, still shows all of it).
def process_path() -> str:
    """The PATH of this process as Windows has it (what a program started by name is looked up on), not Python's copy."""
    import ctypes
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.GetEnvironmentVariableW.argtypes = [ctypes.c_wchar_p, ctypes.c_wchar_p, ctypes.c_uint]
    kernel.GetEnvironmentVariableW.restype = ctypes.c_uint
    size = kernel.GetEnvironmentVariableW("PATH", None, 0)
    buffer = ctypes.create_unicode_buffer(size + 1)
    kernel.GetEnvironmentVariableW("PATH", buffer, size + 1)
    return buffer.value


def cut_process_path(value: str) -> None:
    """What a library does behind Python's back: set the process's PATH without touching os.environ."""
    import ctypes
    kernel = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel.SetEnvironmentVariableW.argtypes = [ctypes.c_wchar_p, ctypes.c_wchar_p]
    assert kernel.SetEnvironmentVariableW("PATH", value)


@pytest.fixture
def path_is_put_back():
    """Yields the PATH as Python has it; afterwards puts it back into the process whatever the test did to it (os.environ's
    assignment calls SetEnvironmentVariable, so this restores the process's PATH too)."""
    saved = os.environ["PATH"]
    yield saved
    os.environ["PATH"] = saved
