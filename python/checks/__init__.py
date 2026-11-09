"""checks: the cases of check_tools.py, one module per tool. What is in this file is what every case shares.

Why it exists
    Before the test bench uses a tool on the motor's geometry, each tool has to give one known answer, driven from Python the way
    the pipeline will drive it (a ring that spins, a shaft that vibrates, a wire that carries a current). check_tools.py runs
    all of them. This package holds one module per tool, and this file holds the parts they have in common, so that a case is
    only its own physics.

What a case module looks like
    def run(folder: Path) -> Iterator[Row]:
        ...
        yield Row(tool, case, expected, measured, verdict, seconds)
    `folder` is out\\check_tools\\<name>\\, the case's own scratch folder. A case gives its rows one by one (a generator), so that
    the rows it made are kept if a later part cannot run. A tool that is not there is `raise Missing(name, where)`: the
    framework turns it into the row "missing: <name> (<where>)", with no traceback. Any other error becomes a failed row that
    says what it was, with the traceback in <folder>\\traceback.txt.

The modules (the name is what `check_tools.py --only` takes)
    femm        femm_check.py      FEMM 5.0: check_femm.py's 10 A solve, the torque by two methods
    gmsh        gmsh_check.py      Gmsh: motorgeo.to_gmsh on the default rotor (the Python module), and gmsh -version (the app)
    calculix    ccx_check.py       CalculiX: the spinning ring against Lame's formula
    prepomax    prepomax_check.py  PrePoMax: its version, and the version of the solver its settings name (nrbf.py reads settings.bin)
    ross        ross_check.py      ROSS, in its own environment: a shaft's first natural frequency (ross_shaft_case.py runs there)
    elmer       elmer_check.py     Elmer: a round conductor in MagnetoDynamics2D against Ampere's law
    salib       salib_check.py     SALib: Sobol indices of the Ishigami function
    pymoo       pymoo_check.py     pymoo: NSGA-II on ZDT1
    matlab      matlab_check.py    MATLAB: matlab -batch loads a copy of the plant, then simulates a small model; what the job held

The rules every case follows (R10 of the tools plan)
    * A program is started only through femm_queue.run_solver: one heavy run at a time on this PC, at below-normal priority, under
      a timeout, in a job that ends with its caller. tests/test_check_tools.py scans this folder, and check_tools.py, for anything
      else, and allows no exception.
    * Python that runs in-process (Gmsh's API, SALib, pymoo) is seconds of CPU and takes no slot; check_tools.py lowers its own process
      to below normal at the start, so that runs below normal too.
    * A whole run holds the lock of out\\check_tools (hold_the_run) from its first row to its last: a case empties its folder before it asks
      for the solver slot, so two runs at once would delete each other's files. A second run waits, and says who it waits for.
    * Nothing here ends a process by its number: no taskkill, no os.kill. The job that run_solver puts a program in ends
      everything below it, at the program's timeout and when its caller dies, however it dies.
    * Every path comes from config.toml through paths.py; a case folder is emptied only inside out\\check_tools.
    * A case compares with a number that is worked out from a formula in the module (not typed in from a previous run), and the
      formula is tested in tests/test_check_tools.py.

Self-check (the answers follow).
  1. Why does a case yield its rows one at a time instead of returning a list?
  2. Why is fresh() so careful about which folders it will empty?
  3. Why is there no module called ross.py, gmsh.py or femm.py here?
Answers.
  1. If the second part of a case cannot run (the Gmsh app is missing, say), the row of the first part (the Python module meshed
     the rotor) must still be in the table. With a list, the exception would throw it away.
  2. out\\check_tools also holds the scratch of other tasks (matlab_mcp\\, verify\\). A function that empties a folder must never
     be handed the wrong one by a typo or a stray "..": it resolves the path and refuses anything that is not inside
     out\\check_tools (or is that folder itself), and it refuses those two folders by name, matlab_mcp\\ (but for its run\\) and
     verify\\, so that no call, however it was built, can empty them.
  3. A script in this folder that is started directly puts this folder first on Python's search path, so `import ross` inside
     it would find our file instead of the library. The modules are named <tool>_check.py for that reason.
"""
from __future__ import annotations

import hashlib
import importlib
import os
import shutil
import time
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import Callable, Iterator, NamedTuple

import femm_queue
import paths

COLUMNS = ("tool", "case", "expected", "measured", "verdict", "seconds")
PASS, FAIL = "PASS", "FAIL"
NOT_STARTED = "(not started)"


@contextmanager
def intact_path():
    """Keep the PATH of this process as it was when the block was entered, whatever a library does to it inside.

    Why: on Windows, the Python gmsh module's initialize() cuts the PATH of the running process down to its first few entries,
    and finalize() does not put it back (measured 2026-11-06 with gmsh 4.15.2: 1,864 characters became 313; the PATH that Python
    shows in os.environ stays whole, only the process's own is cut). A program started BY NAME afterwards, such as
    powershell.exe, which lives only on the PATH, is then not found (FileNotFoundError). Programs started by their full path,
    as every case here does, and children, which get a copy of os.environ, are not affected. Assigning os.environ["PATH"] sets the
    process's variable again. check_tools.py wraps every check in this, and the cases that call gmsh wrap the call.
    """
    saved = os.environ.get("PATH")
    try:
        yield
    finally:
        if saved is not None:
            os.environ["PATH"] = saved


class Row(NamedTuple):
    """One line of the table and of out\\check_tools.csv. `seconds` is the time this row's part took."""
    tool: str
    case: str
    expected: str
    measured: str
    verdict: str
    seconds: float


class Missing(Exception):
    """A tool (or a file) the case needs is not there; check_tools shows it as the row `missing: <name> (<where>)`."""

    def __init__(self, name: str, where, case: str = NOT_STARTED):
        super().__init__(f"missing: {name} ({where})")
        self.name, self.where, self.case = name, str(where), case
        self.also: list[Missing] = []               # other tools of the same case that are missing too (see need_tools): a row each


class Check(NamedTuple):
    """What check_tools runs: the label shown in the tool column, and run(folder), which yields Rows."""
    label: str
    run: Callable[[Path], Iterator[Row]]


def verdict(ok: bool) -> str:
    return PASS if ok else FAIL


def need_tool(name: str, case: str = NOT_STARTED) -> Path:
    """The [tools] entry `name` of config.toml as a Path, or Missing that says where the config looks (or that it names nothing)."""
    try:
        return paths.tool(name)
    except FileNotFoundError:
        pass
    try:
        value = paths.setting("tools", name)
    except (KeyError, ValueError, FileNotFoundError):
        raise Missing(name, "not named in config.toml", case) from None
    if not value:
        raise Missing(name, "empty in config.toml", case) from None
    raise Missing(name, paths.P / value, case) from None


def need_tools(*names: str, case: str = NOT_STARTED) -> list[Path]:
    """need_tool for several tools at once. If some are missing, ALL of them are named: the first is raised as Missing, with the
    others in its .also, so that one run tells the whole list and not one tool per run."""
    found, missing = [], []
    for name in names:
        try:
            found.append(need_tool(name, case))
        except Missing as error:
            missing.append(error)
    if missing:
        missing[0].also = missing[1:]
        raise missing[0]
    return found


def need_file(name: str, path, case: str = NOT_STARTED) -> Path:
    """`path` if it is a file, else Missing. For what is not a [tools] entry: ROSS's python, a reference model, settings.bin."""
    if not Path(path).is_file():
        raise Missing(name, path, case)
    return Path(path)


def save_log(folder: Path, name: str, result) -> None:
    """Keep what a program printed in <folder>\\<name>.out.txt and .err.txt (a SolverRun's stdout and stderr), for the look afterwards."""
    (Path(folder) / f"{name}.out.txt").write_text(result.stdout, encoding="utf-8", errors="replace")
    (Path(folder) / f"{name}.err.txt").write_text(result.stderr, encoding="utf-8", errors="replace")


def sha256(path) -> str:
    """The SHA-256 of a file, as hex. The reference files are hashed before and after every run that copies them."""
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_root() -> Path:
    """out\\check_tools: the one folder the cases write in (it also holds the scratch of other tasks, which nobody here touches)."""
    return paths.path("out") / "check_tools"


def case_folder(name: str) -> Path:
    """out\\check_tools\\<name>, made if it is not there (and left as it is if it is: a case empties its own scratch with fresh())."""
    folder = check_root() / name
    folder.mkdir(parents=True, exist_ok=True)
    return folder


RUN_LOCK = ".lock"                          # inside out\\check_tools: the folder of the lock one whole run holds (hold_the_run); not a case, and never emptied


@contextmanager
def hold_the_run(label: str, *, say=None, poll_s: float = 30.0, look_s: float = 0.5):
    """Hold the lock of out\\check_tools for the `with` block, a whole check_tools run; yields femm_queue's holder note {"label", "pid", "since"}.

    Why (the final review, N2): every case empties its folder under out\\check_tools at the start (fresh), before it asks for the solver slot, so two runs at once would
    delete the files the other is using. With this lock a second run waits first. It says "waiting for another check_tools run: held by '<label>' (pid N) since HH:MM" at the
    first look and then every `poll_s` seconds, looks every `look_s` seconds (so that a Ctrl-C is noticed within that), and goes on when the first lets go.

    The lock is femm_queue's own, an operating-system lock on a file, taken with solver_slot(lock_dir=out\\check_tools\\.lock) and not the machine's solver slot (that is
    run_solver's, taken row by row): bench.hold_root does the same for a bench root. The system gives it up when the run ends however it ends, a killed run included."""
    say = say if say is not None else (lambda text: print(text, flush=True))
    lock_dir = check_root() / RUN_LOCK
    stack = ExitStack()
    last_said = None
    while True:
        try:
            note = stack.enter_context(femm_queue.solver_slot(label, wait=False, lock_dir=lock_dir))
            break
        except femm_queue.QueueBusy as busy:
            now = time.monotonic()
            if last_said is None or now - last_said >= poll_s:
                say(f"waiting for another check_tools run: held by {femm_queue._who(busy.holder)}")
                last_said = now
            time.sleep(look_s)
    try:
        yield note
    finally:
        stack.close()


RESERVED = ("matlab_mcp", "verify")        # folders of out\check_tools that hold files of earlier work: they are never emptied
ALLOWED_INSIDE_RESERVED = ("matlab_mcp", "run")     # ... except this sub-folder (the rest of matlab_mcp\ is never touched)


def fresh(folder) -> Path:
    """Empty `folder` (delete it and make it again) and return it. Only a folder inside out\\check_tools may be emptied, never that
    folder itself and never anything outside it: the path is resolved first, so a ".." or a junction cannot lead out. The folders
    of earlier work, out\\check_tools\\verify and out\\check_tools\\matlab_mcp, are refused too; only matlab_mcp\\run (and below) may
    be emptied."""
    root = check_root().resolve()
    target = Path(folder).resolve()
    if target == root or root not in target.parents:
        raise ValueError(f"{target} is not a case folder inside {root}: refusing to empty it")
    inside = tuple(part.lower() for part in target.relative_to(root).parts)
    if inside[0] in RESERVED and inside[:2] != ALLOWED_INSIDE_RESERVED:
        raise ValueError(f"{target} holds the files of an earlier task (out\\check_tools\\{inside[0]}): refusing to empty it "
                         f"(only out\\check_tools\\{ALLOWED_INSIDE_RESERVED[0]}\\{ALLOWED_INSIDE_RESERVED[1]} may be emptied)")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    return Path(folder)


# What `check_tools.py --only` takes, the label of the tool column, and the module that holds the case (all under checks/).
_CASES = (
    ("femm", "FEMM", "femm_check"),
    ("gmsh", "Gmsh", "gmsh_check"),
    ("calculix", "CalculiX", "ccx_check"),
    ("prepomax", "PrePoMax", "prepomax_check"),
    ("ross", "ROSS", "ross_check"),
    ("elmer", "Elmer", "elmer_check"),
    ("salib", "SALib", "salib_check"),
    ("pymoo", "pymoo", "pymoo_check"),
    ("matlab", "MATLAB", "matlab_check"),
)


def _lazy(module: str) -> Callable[[Path], Iterator[Row]]:
    """run(folder) that imports its module only when it is called: `--only gmsh` must not import ROSS's or MATLAB's code."""
    def run(folder: Path) -> Iterator[Row]:
        yield from importlib.import_module(f".{module}", __package__).run(folder)
    return run


def registry() -> dict[str, Check]:
    """The nine rows of check_tools.py in the order of the brief's table: {name: Check(label, run)}."""
    return {name: Check(label, _lazy(module)) for name, label, module in _CASES}
