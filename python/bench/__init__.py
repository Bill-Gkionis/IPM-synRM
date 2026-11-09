"""bench: each tool the pipeline will use, tried on the motor's own geometry, driven from Python the way the pipeline will drive it.

Why it exists
    Task 4 (check_tools.py) proved every tool on a KNOWN answer: a ring, a shaft, a wire. That is necessary and not enough. The
    pipeline will not draw rings; it will hand the motor to FEMM, to Gmsh, to Elmer, and read numbers back. The bench does exactly
    that, one chain at a time (Tasks 11 to 13 of the tools plan), and holds every number it gets against an expectation whose source
    is named: an earlier result, a hand estimate, or another tool. A result is "sound" when it makes physical sense, not when it is
    merely a number.

The two rules behind every case (R10 and R11 of the tools plan)
    R11  The geometry is motorgeo.Params(barrier="none"): the air band behind each magnet is left out, in Python (a design decision of
         2026-11-02). The CAD and FEM files (reference\\ and the Desktop folders) are read-only inputs; a case that needs a
         band-free version of one builds it in memory.
    R10  One heavy run at a time, machine-wide: every solver goes through femm_queue (femm_pool for FEMM, run_solver for every
         other program). Each run is at below-normal priority, under a timeout sized to the case, with at most 8 workers or threads
         (config.toml, [femm] workers), in a Windows job that ends with its caller. Nothing in bench/ starts a program itself or ends a
         process (tests/test_bench.py scans for both, and for the time limit of every call).

The cases (this folder: one module per case; bench/run_bench.py runs them)
    B1  b1_femm.py   Python -> FEMM 5.0 (femm_pool): no-load flux linkage, MTPA torque, cogging, on the no-band 90-degree sector
    B2  b2_gmsh.py   Python -> Gmsh (the API) -> .msh: motorgeo.to_gmsh on the sector, and the CAD drawings through ezdxf, shapely, Gmsh
    B3  b3_elmer.py  Gmsh -> Elmer (ElmerGrid, ElmerSolver, MagnetoDynamics2D): the no-load field, against FEMM's
    B4  b4_ccx.py    Gmsh -> CalculiX 2.23 (ccx, plane stress, 6-node triangles): the band-free rotor under centrifugal load at 10,000 and 12,000 rpm; the peak in the bridges
    B5  b5_prepomax.py  PrePoMax 3.0.0 from the command line with no window: B4's deck, CalculiX 2.23, B4's peak (b5_prepomax_script.py runs inside PrePoMax)
    B6  b6_ross.py   Python -> ROSS (its own environment envs\\ross; b6_ross_case.py runs there): the rotor on its shaft and two bearings, critical speeds, Campbell diagram
    B7  b7_salib.py  Python -> SALib -> FEMM (femm_pool): Morris screening of four rotor dimensions on the coarse-mesh torque at 10 A, run twice (femm_points.py is the FEMM side)
    B8  b8_matlab.py Python -> MATLAB (matlab -batch, b8_mtpa_map.m): a flux map from FEMM saved as .mat, interpolated in MATLAB, the MTPA angle read back
    B9  b9_rules.py  The rules themselves: two bench runs at once, in two processes, each in a root of its own (and two on one root, which the root lock serializes)

Where a case writes (Task 13, B9's rulings)
    The root is out\\bench, or the ROOT of `run_bench.py --out ROOT` (inside out\\; use_root() sets it, root() says it). Every case writes under ROOT\\<case>\\, and fresh() refuses
    anything outside the active root. A run holds a lock on its root for the whole run (hold_root(): femm_queue's own operating-system lock, taken with solver_slot(lock_dir=
    ROOT\\.lock)), so a second run on the SAME root waits, printing who holds it, and cannot empty the first run's folders. Two runs on two roots do not meet there; they meet at
    the solver slot (femm_queue), which is machine-wide.

What every case leaves, in ROOT\\<case>\\   (ROOT is out\\bench unless told otherwise; out\\ is git-ignored; no script of the bench writes anywhere else, and none into the project notes)
    results.json   the inputs, the model and mesh sizes, every number, and `checks`: each number against its expectation, with the
                   source of the expectation, a verdict (PASS or FAIL) and `gating`. A gating check (the default) decides the exit
                   status of run_bench.py. One that is not gating is made, shown as REPORTED with its own pass or fail, and kept
                   here, but does not decide it; its expected text says why (B3's plain-curve run: the verdict is the laminated
                   run's). One that cannot be made on this machine is SKIPPED, with the reason. So the list of checks is the same
                   on every PC, and nothing is hidden
    figure.png     the band-free motor inside that tool: the tool's own render where it has one (Gmsh), else the tool's result plotted
                   from Python in the style of learn/tools_figures.py (the whole cross-section from the sector, turbo for |B|, flux
                   lines as contours of A). It is copied into the project notes by hand; no script here does
    run.log        what the run said as it happened (the waiting line of the solver slot, the pool's start-up line, every program)
    runs.json      one record per run: its kind (femm_pool or run_solver), label, the priority class it was SEEN at (read back from the
                   operating system, not assumed), and its worker or thread count. B9 reads these
    plus the case's own files (meshes, model files, solver input and output), which the next run of the case empties first

Run (from the Python to FEMM folder)
    uv run python bench/run_bench.py --only B1,B2,B3        (or: uv run python -m bench.run_bench --only B1; B4 ... B9 likewise)
    uv run python bench/run_bench.py --out out\\bench\\mine  --only B6      the same, written under another root (inside out\\)
    uv run pytest tests/test_bench*.py -q                   the tests: no solver is run (test_bench.py: this folder's shared part, run_bench.py
                                                            and the scans; one file for each case)

Self-check (the answers follow).
  1. A case's folder is emptied before the case runs. Why does fresh() refuse out\\bench\\B1\\fem?
     (And, since Task 13: why must a root be inside out\\, and why does a second run on the same root wait?)
  2. The log names the priority class a program was SEEN at, read back by looking at it while it ran. Why not just write
     "below normal" because that is what was asked for?
  3. A number is 0.4 % from its expectation and the tolerance is 0.5 %. Another is 0.6 %. What does the bench do with the second?
Answers.
  1. A function that empties folders must never be handed a wrong one by a typo or a stray "..". It resolves the path first and
     accepts only a folder that is a direct child of the active root and is named like a case (B1, B2, ...); everything else is a
     ValueError, so a mistake costs an error message and not somebody's files. A root is only accepted inside out\\ for the same reason:
     the bench empties ROOT\\B1, and a root typed wrong would aim that at somebody's folder. A second run on the same root waits because
     its first act is to empty ROOT\\B<n>, which is where the first run is working: the lock on the root makes the two take turns.
  2. A request is not a fact. The class is set when the program is created; if something went wrong (a launcher that does not pass it
     on, a different Windows rule) the program would run at normal priority while the log said otherwise. Reading it back
     is the only way the log can be evidence for B9's "every run was below normal".
  3. It reports the second as FAILED, with the difference and its best explanation, and does not touch the tolerance or the expected
     value. A bench that adjusts its yardstick until everything passes measures nothing.
"""
from __future__ import annotations

import dataclasses
import json
import math
import re
import shutil
import sys
import textwrap
import threading
import time
from contextlib import ExitStack, contextmanager
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

import femm_queue
import paths

__all__ = ["root", "default_root", "root_allowed", "use_root", "root_lock_dir", "hold_root", "case_folder", "fresh", "compare", "check", "skipped", "is_gating", "describe", "params", "write_json", "fig_field", "legend_and_model_figure", "wrap_title", "priority_class",
           "class_name", "process_record", "BELOW_NORMAL", "OWN_LIMIT", "ChildWatcher", "watch_children", "HolderSampler", "Log", "run_program"]

OPTIONS: dict = {}                                                  # what run_bench.py was asked beyond the cases: {"mesh_study": True} for B1's finer-mesh point
BELOW_NORMAL = 0x4000                                               # Windows BELOW_NORMAL_PRIORITY_CLASS
CLASS_NAMES = {0x20: "normal", 0x40: "idle", 0x80: "high", 0x100: "realtime", 0x4000: "below_normal", 0x8000: "above_normal"}
CASE_NAME = re.compile(r"B[1-9][0-9]*")                             # B1 ... B9: the names of the bench's cases
NOT_GATING = "reported, not gating"                                 # what a check that does not decide the verdict says about itself
OWN_LIMIT = "the bench's own limit; the brief gives none"           # said, in brackets, next to every limit that the bench chose and the brief did not give
TITLE_CHARS = 92                                                    # a picture's title is wrapped to this many characters a line: it fits the 9-inch figure
POLL_S = 0.025                                                      # how often the watcher looks at the processes under this one


# ============================================================================================== where a case writes
_ACTIVE_ROOT: Path | None = None                                    # set by use_root() (run_bench.py --out ROOT); None means the default, out\bench
ROOT_LOCK_DIR = ".lock"                                             # inside a root: the folder of its lock (femm_queue's solver.lock and holder.json live in it)


def default_root() -> Path:
    """out\\bench: where the bench writes unless a run is told another root (config.toml [paths] out, so a test can point it somewhere else)."""
    return paths.path("out") / "bench"


def root() -> Path:
    """The root of this run: the folder every case writes under (ROOT\\<case>\\) and fresh() obeys. out\\bench, or the ROOT of `run_bench.py --out ROOT` (use_root)."""
    return _ACTIVE_ROOT if _ACTIVE_ROOT is not None else default_root()


def root_allowed(candidate) -> Path:
    """`candidate` as an absolute, resolved path, if it may be a root: strictly INSIDE out\\ (config.toml [paths] out), and not out\\ itself. ValueError otherwise.

    Why: the bench empties ROOT\\B<n> (fresh) and writes its logs in ROOT. A root typed wrong, or aimed with a stray "..", must not reach the project notes, reference\\ or somebody's folder
    that happens to hold a B1. A relative candidate starts at the project folder, as everywhere in paths.py. The path is resolved first, so ".." and junctions cannot lead out."""
    out = paths.path("out").resolve()
    target = Path(candidate)
    if not target.is_absolute():
        target = paths.P / target
    target = target.resolve()
    if target == out or out not in target.parents:
        raise ValueError(f"{target} is not a root the bench may write in: it has to be a folder inside {out} (and not that folder itself)")
    return target


@contextmanager
def use_root(candidate):
    """`with use_root(ROOT) as root:` the bench writes under ROOT (root_allowed: inside out\\) until the block ends, and then under what it wrote under before."""
    global _ACTIVE_ROOT
    target = root_allowed(candidate)
    previous, _ACTIVE_ROOT = _ACTIVE_ROOT, target
    try:
        yield target
    finally:
        _ACTIVE_ROOT = previous


def case_folder(case: str) -> Path:
    """ROOT\\<case> (out\\bench\\<case> by default), made if it is not there."""
    folder = root() / case
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def fresh(folder) -> Path:
    """Empty `folder` (delete it and make it again) and return it. Only ROOT\\B<n> of the ACTIVE root may be emptied: the path is resolved first, so a ".." or a junction cannot
    lead out, and anything that is not a direct child of the active root named like a case is a ValueError (another root's case folder included)."""
    base = root().resolve()
    target = Path(folder).resolve()
    if target.parent != base or not CASE_NAME.fullmatch(target.name):
        raise ValueError(f"{target} is not a case folder ({base}\\B<n>) of the active root: refusing to empty it")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    return Path(folder)


def root_lock_dir(root_folder=None) -> Path:
    """The folder of a root's lock: ROOT\\.lock (the active root's if none is given). fresh() empties case folders only, so the lock file is never deleted (femm_queue's rule)."""
    return (Path(root_folder) if root_folder is not None else root()) / ROOT_LOCK_DIR


def _holder_words(info: dict | None) -> str:
    """"'label' (pid N) since HH:MM" for a holder note, as femm_queue words its own waiting line; a holder that cannot be named says so."""
    if not info or info.get("pid") is None:
        return "a run whose note (holder.json) cannot be read"
    try:
        since = datetime.fromisoformat(info["since"]).strftime("%H:%M")
    except (TypeError, ValueError):
        since = "?"
    return f"'{info['label']}' (pid {info['pid']}) since {since}"


@contextmanager
def hold_root(label: str, *, folder=None, say=None, poll_s: float = 30.0, look_s: float = 0.5):
    """Hold the lock of a root for the `with` block (a whole bench run): yields femm_queue's holder note {"label", "pid", "since"}.

    The lock is femm_queue's own: an operating-system lock on a file, taken with solver_slot(lock_dir=ROOT\\.lock). Reused, not copied: the kernel gives it up when the holding
    process ends however it ends, so a killed run frees its root at once. While another run holds it this waits, looking every `look_s` seconds, and says (through `say`, which
    prints by default) 'waiting for the bench root <ROOT>: held by '<label>' (pid N) since HH:MM' at the first look and then every `poll_s`. Ctrl-C while waiting takes nothing and
    leaves nothing. `folder` is the ROOT (the active root if none is given). RuntimeError, as for the slot, if this thread already holds it."""
    say = say if say is not None else (lambda text: print(text, flush=True))
    target = Path(folder) if folder is not None else root()
    stack = ExitStack()
    last_said = None
    while True:
        try:
            note = stack.enter_context(femm_queue.solver_slot(label, wait=False, lock_dir=root_lock_dir(target)))
            break
        except femm_queue.QueueBusy as busy:
            now = time.monotonic()
            if last_said is None or now - last_said >= poll_s:
                say(f"waiting for the bench root {target}: held by {_holder_words(busy.holder)}")
                last_said = now
            time.sleep(look_s)
    try:
        yield note
    finally:
        stack.close()


# ============================================================================================== a number against its expectation
def _not_gating(name: str, reason: str) -> str:
    """The words a check that does not decide the verdict carries: "reported, not gating (<reason>)". A check that does not gate has to say why."""
    if not reason:
        raise ValueError(f"{name}: a check that does not gate the verdict has to say why (reason=...)")
    return f"{NOT_GATING} ({reason})"


def compare(name: str, measured: float, expected: float, tolerance_percent: float, source: str, unit: str = "", *, gating: bool = True,
            reason: str = "") -> dict:
    """One number held against its expectation: {name, unit, expected, source, measured, difference_percent, tolerance_percent, verdict, gating}.

    The difference is (measured - expected) / |expected| in percent, and the verdict is PASS when its size is within the tolerance
    (a number exactly on the tolerance passes: "within 0.5 %"). `source` says where the expectation comes from: an earlier result,
    a hand estimate, another tool. A zero expectation is a ValueError: a relative difference from zero means nothing, so the case has
    to say what it means (an absolute limit, with check()).

    `gating` is True unless the case says otherwise. With gating=False (and the `reason`) the check keeps its own verdict and is shown as REPORTED, but
    it does not decide run_bench.py's exit status; the item carries `note`, "reported, not gating (<reason>)", which describe() puts before the expectation."""
    if expected == 0:
        raise ValueError(f"{name}: the expected value is 0, a percentage difference from it is meaningless")
    difference = 100.0 * (measured - expected) / abs(expected)
    item = {"name": name, "unit": unit, "expected": expected, "source": source, "measured": measured, "difference_percent": difference,
            "tolerance_percent": tolerance_percent, "verdict": "PASS" if abs(difference) <= tolerance_percent + 1e-9 else "FAIL", "gating": gating}
    if not gating:
        item["note"] = _not_gating(name, reason)
    return item


def check(name: str, ok: bool, expected: str, measured: str, source: str, *, gating: bool = True, reason: str = "", **extra) -> dict:
    """A check that is not one number: {name, expected, measured, source, verdict, gating, ...extra} with the words of both sides.

    With gating=False (and the `reason`) the check keeps its own verdict but does not decide the exit status, and its expected text begins
    "reported, not gating (<reason>): "."""
    if not gating:
        expected = f"{_not_gating(name, reason)}: {expected}"
    return {"name": name, "expected": expected, "measured": measured, "source": source, "verdict": "PASS" if ok else "FAIL", "gating": gating, **extra}


def skipped(name: str, expected: str, reason: str, source: str) -> dict:
    """A check that cannot be made on this machine (a file it reads is not there): verdict SKIPPED, never gating, the reason where the measurement would be.
    It is in the list all the same, so that the list of checks is the same on every PC and the reader sees what was not looked at."""
    return {"name": name, "expected": expected, "measured": reason, "source": source, "verdict": "SKIPPED", "gating": False}


def is_gating(item: dict) -> bool:
    """Does this check decide the exit status? Yes unless it says gating False (an item made without the field, by an older case, gates)."""
    return bool(item.get("gating", True))


def describe(item: dict) -> tuple[str, str]:
    """(expected, measured) of a check or a compare() result as words, for the table the run prints. A compare() result that does not gate says so first."""
    if "difference_percent" in item:
        unit = f" {item['unit']}" if item.get("unit") else ""
        said = f"{item['note']}: " if item.get("note") else ""
        return (f"{said}{item['expected']:.6g}{unit} within {item['tolerance_percent']:g} %  [{item['source']}]",
                f"{item['measured']:.6g}{unit}  ({item['difference_percent']:+.3f} %)")
    return f"{item['expected']}  [{item['source']}]", str(item["measured"])


# ============================================================================================== the results file
def _plain(value: Any) -> Any:
    """`value` made of what JSON has: numpy numbers and arrays, tuples, sets and paths turned into plain ones, NaN and infinity into None."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [_plain(item) for item in (sorted(value) if isinstance(value, (set, frozenset)) else value)]
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if math.isfinite(number) else None
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path, data) -> Path:
    """Write `data` as UTF-8 JSON (indented, non-ASCII text kept as it is). NaN is not JSON: it is written as null."""
    target = Path(path)
    target.write_text(json.dumps(_plain(data), indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    return target


def params(p) -> dict:
    """The parameters of the motor (a motorgeo.Params) as a plain dict, for the `params` entry of a case's results.json. The one copy: B1, B2 and B3 use it."""
    return dataclasses.asdict(p)


# ============================================================================================== the field picture
def wrap_title(title: str) -> str:
    """A picture's title with every line wrapped to TITLE_CHARS characters (the line breaks the case wrote are kept): a picture saved with a tight box is as wide as
    its widest line, so a long one would leave most of the figure blank. fig_field and B2's bands.png both use it."""
    return "\n".join(textwrap.fill(line, TITLE_CHARS) for line in title.split("\n"))


def fig_field(xy, triangles, potential, flux_density, title: str, target) -> None:
    """The whole cross-section from the sector, as B1 and B3 show it: |B| on every triangle (turbo, 0 to 2 T) and the flux lines as contours of the vector
    potential, the 90-degree sector repeated four times (two poles per sector, periodic). The style of learn/tools_figures.py.

    xy (n, 2) the nodes in mm; triangles (m, 3) their indices; potential (n,) A at the nodes; flux_density (m,) |B| of each triangle [T]. The title may have
    several lines; each is wrapped to TITLE_CHARS characters, so that the picture, saved with a tight box, is not widened by a long line (B1's was 45 % blank)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    xy, triangles, potential = np.asarray(xy, float), np.asarray(triangles), np.asarray(potential, float)
    x, y = xy[:, 0], xy[:, 1]
    fig, ax = plt.subplots(figsize=(9, 9))
    pc = None
    for k in range(4):
        c, s = math.cos(math.radians(90 * k)), math.sin(math.radians(90 * k))
        t = mtri.Triangulation(c * x - s * y, s * x + c * y, triangles)
        pc = ax.tripcolor(t, facecolors=np.asarray(flux_density, float), cmap="turbo", vmin=0, vmax=2.0, shading="flat", rasterized=True)
        ax.tricontour(t, potential, levels=np.linspace(potential.min(), potential.max(), 26)[1:-1], colors="#111111", linewidths=0.4)
    fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02).set_label("|B| [T]")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(wrap_title(title), fontsize=10.5, color="#1d1d1f", loc="left")
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def legend_and_model_figure(legend: np.ndarray, model: np.ndarray, caption: str, target) -> None:
    """A tool's own picture, put together: its colour legend (an RGB image) at the left, the model (an RGB image) at the right, both at their own pixel scale and starting right under
    the caption (every line of it wrapped to TITLE_CHARS characters). cgx's frames (B4) and PrePoMax's viewport (B5) both come this way."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    title_in = 0.25 * (caption.count("\n") + 1) + 0.1                          # room for the caption, so that the two images start right under it
    width_in = 11.0
    height_in = width_in * max(legend.shape[0], model.shape[0]) / (legend.shape[1] + model.shape[1]) + title_in
    fig, axes = plt.subplots(1, 2, figsize=(width_in, height_in), gridspec_kw={"width_ratios": [legend.shape[1], model.shape[1]], "wspace": 0.0})
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1 - title_in / height_in)
    axes[0].imshow(legend)
    axes[1].imshow(model)
    for ax in axes:
        ax.axis("off")
    fig.text(0.0, 1 - 0.05 / height_in, wrap_title(caption), ha="left", va="top", fontsize=10.5, color="#1d1d1f")
    fig.savefig(target, dpi=130, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


# ============================================================================================== priority classes
def priority_class(pid: int) -> int | None:
    """The Windows priority class of a running process, asked of the operating system; None if the process cannot be looked at
    (it has ended, or there is no such process). The look is read-only: it can neither end nor change anything."""
    from check_queue import priority_class as ask         # check_queue.py's own reader (the live check of the queue): one definition
    try:
        return int(ask(pid))
    except OSError:
        return None


def class_name(value: int) -> str:
    """0x4000 -> "below_normal". A value that is not one of Windows' six classes is shown as the hex number, never hidden."""
    return CLASS_NAMES.get(value, f"{value:#x}")


GMSH_THREAD_OPTIONS = ("General.NumThreads", "Mesh.MaxNumThreads1D", "Mesh.MaxNumThreads2D", "Mesh.MaxNumThreads3D")


def gmsh_options() -> dict | None:
    """The thread options of the Gmsh library in this process, asked of Gmsh itself with gmsh.option.getNumber: {option name: value}. None if the gmsh module
    cannot be imported. A Gmsh that is already running is asked and left running; otherwise a session is opened for the question and closed (inside
    checks.intact_path(): gmsh.initialize() cuts the PATH of this process; tests/test_bench.py scans for the guard)."""
    try:
        import gmsh
    except ImportError:
        return None
    from checks import intact_path
    running = bool(gmsh.isInitialized())
    with intact_path():
        if not running:
            gmsh.initialize()
        try:
            if not running:
                gmsh.option.setNumber("General.Terminal", 0)
            return {name: int(gmsh.option.getNumber(name)) for name in GMSH_THREAD_OPTIONS}
        finally:
            if not running:
                gmsh.finalize()


def process_record() -> dict:
    """The record of the bench's own Python for runs.json: it is a run too (Gmsh's library meshes in it, numpy and matplotlib work in it). The class is the one
    the operating system reports for this process now (run_bench lowers it when it starts). The thread count is read back from Gmsh, not written here:
    `gmsh_options` holds what Gmsh answers for General.NumThreads and Mesh.MaxNumThreads1D/2D/3D, and `threads` is what its 2-D mesher uses, Mesh.MaxNumThreads2D
    or, when that is 0, General.NumThreads (1 and 0 on 2026-11-06, Gmsh 4.15.2). `threads_limit` is config.toml's [femm] workers, the most any run may use, so a
    thread count above it shows as above it. If gmsh cannot be imported both are None. What is NOT read back: numpy's own BLAS threads (threadpoolctl, which
    would say, is not installed). The bench's one matrix product is a rotation of at most a few hundred points by a 2 x 2 matrix, so none is expected; the
    record claims only what Gmsh says."""
    import os
    value = priority_class(os.getpid())
    options = gmsh_options()
    threads = (options["Mesh.MaxNumThreads2D"] or options["General.NumThreads"]) if options else None
    return {"kind": "in_process", "label": "the bench's own Python (the Gmsh library, numpy, matplotlib)", "pid": os.getpid(), "threads": threads,
            "threads_limit": paths.setting("femm", "workers"), "gmsh_options": options,
            "priority_class_seen": class_name(value) if value is not None else "unreadable", "time": time.strftime("%Y-%m-%d %H:%M:%S")}


class ChildWatcher(threading.Thread):
    """While a program runs, look every POLL_S seconds at the processes that descend from this one and read each one's priority class.

    run_solver starts its program as a child of the caller; what that program starts is below it. So "what did the operating system say the
    class of ElmerSolver was while it ran" is answered by reading the classes of this process's descendants. Read-only: it never opens a
    process for more than querying it, and it never ends one. A program that lives for less than POLL_S may be missed (ElmerGrid is
    about 0.2 s: caught); `seen()` then lacks it, and the log says so instead of guessing.

    It also counts, at every look, how many processes of each program name are alive under this one, and keeps the largest count (`peaks()`): B9's 'at most 8 workers' read
    from outside, not from what the pool says about itself."""

    def __init__(self, poll_s: float = POLL_S):
        super().__init__(daemon=True, name="child-watcher")
        import os
        self._me = os.getpid()
        self._poll_s = poll_s
        self._seen: dict[str, set[int]] = {}
        self._peak: dict[str, int] = {}
        self._stop_now = threading.Event()
        self.looks = 0

    def run(self) -> None:
        from checks.matlab_check import process_table             # a read-only Toolhelp snapshot: {pid: (parent pid, program name)}
        while not self._stop_now.is_set():
            self.looks += 1
            try:
                table = process_table()
            except OSError:
                table = {}
            below = {self._me}
            changed = True
            while changed:                                         # every process whose parent chain reaches this one
                changed = False
                for pid, (ppid, _) in table.items():
                    if ppid in below and pid not in below:
                        below.add(pid)
                        changed = True
            alive: dict[str, int] = {}
            for pid in below - {self._me}:
                value = priority_class(pid)
                if value is not None:
                    self._seen.setdefault(table[pid][1], set()).add(value)
                    alive[table[pid][1].lower()] = alive.get(table[pid][1].lower(), 0) + 1
            for name, count in alive.items():
                self._peak[name] = max(self._peak.get(name, 0), count)
            self._stop_now.wait(self._poll_s)

    def stop(self) -> None:
        self._stop_now.set()
        self.join(5)

    def seen(self) -> dict[str, set[int]]:
        """{program name: the priority classes it was seen at}."""
        return {name: set(values) for name, values in self._seen.items()}

    def peaks(self) -> dict[str, int]:
        """{program name in lower case: the most processes of that name alive under this one at any single look}."""
        return dict(sorted(self._peak.items()))

    def summary(self) -> dict[str, list[str]]:
        """{program name in lower case: its class names, sorted}: what goes into runs.json."""
        merged: dict[str, set[str]] = {}
        for name, values in self._seen.items():
            merged.setdefault(name.lower(), set()).update(class_name(value) for value in values)
        return {name: sorted(values) for name, values in sorted(merged.items())}


class HolderSampler(threading.Thread):
    """Looks at who holds a lock every `poll_s` seconds and keeps (seconds since `zero`, None for a free lock or (label, pid)): femm_queue.holder() for the solver slot (no
    `lock_dir`) or for another lock (a root's, hold_root). Read-only (asking takes the lock for an instant when it is free: femm_queue's own way of asking). B8 looks at the slot
    while MATLAB runs, and B9 looks at the slot and at a root's lock while its two runs take turns."""

    def __init__(self, zero: float, poll_s: float = 0.2, lock_dir=None):
        super().__init__(daemon=True, name="holder-sampler")
        self.zero, self.poll_s, self.lock_dir = zero, poll_s, lock_dir
        self.samples: list[tuple[float, tuple | None]] = []
        self.errors = 0
        self._stop_now = threading.Event()

    def run(self) -> None:
        while not self._stop_now.is_set():
            try:
                info = femm_queue.holder() if self.lock_dir is None else femm_queue.holder(self.lock_dir)
                seen = None if info is None else (info["label"], info["pid"])
            except OSError:
                self.errors += 1
                seen = ("unreadable", None)
            self.samples.append((time.perf_counter() - self.zero, seen))
            self._stop_now.wait(self.poll_s)

    def stop(self) -> None:
        self._stop_now.set()
        self.join(5)


@contextmanager
def watch_children():
    """`with watch_children() as watcher:` around a run_solver call; afterwards watcher.seen() / .summary() say what was seen."""
    watcher = ChildWatcher()
    watcher.start()
    try:
        yield watcher
    finally:
        watcher.stop()


# ============================================================================================== the log of a run
class _Tee:
    """A stand-in for sys.stdout that passes everything on and keeps the complete lines in the log."""

    def __init__(self, log: "Log", stream):
        self._log, self._stream, self._partial = log, stream, ""

    def write(self, text: str) -> int:
        written = self._stream.write(text)
        self._partial += text
        while "\n" in self._partial:
            line, self._partial = self._partial.split("\n", 1)
            self._log.lines.append(line)
        return len(text) if written is None else written

    def flush(self) -> None:
        self._stream.flush()

    def __getattr__(self, name):                                    # encoding, isatty, ...: whatever the real stream has
        return getattr(self._stream, name)


class Log:
    """What a case says while it runs (run.log), and one record for each run (runs.json), kept in the case's folder.

    `say(text)` prints a line and keeps it. `tee()` is a `with` block inside which everything printed (a pool's start-up line, the waiting
    line of the solver slot) is kept too. `add_run(record)` keeps the record of a run: its kind, label, the priority class it was seen at
    and its worker or thread count (B9 reads them). `save()` writes both files."""

    def __init__(self, folder):
        self.folder = Path(folder)
        self.lines: list[str] = []
        self.runs: list[dict] = []
        self._teeing = False

    def say(self, text: str = "") -> None:
        print(text, flush=True)
        if not self._teeing:                                        # inside tee() the tee has already kept it
            self.lines.extend(text.split("\n"))

    @contextmanager
    def tee(self):
        stream = sys.stdout
        tee = _Tee(self, stream)
        self._teeing, sys.stdout = True, tee
        try:
            yield self
        finally:
            sys.stdout, self._teeing = stream, False
            if tee._partial:
                self.lines.append(tee._partial)

    def add_run(self, record: dict) -> dict:
        self.runs.append(record)
        return record

    def save(self) -> None:
        self.folder.mkdir(parents=True, exist_ok=True)
        (self.folder / "run.log").write_text("\n".join(self.lines) + "\n", encoding="utf-8")
        write_json(self.folder / "runs.json", self.runs)


def run_program(log: Log, cmd, *, label: str, cwd, timeout_s: float, threads: int, name: str):
    """femm_queue.run_solver, watched and logged: the program runs under the solver slot at below-normal priority, and the log gets a record of it.

    The record: its kind ("run_solver"), label, the program, the thread count it was given (OMP_NUM_THREADS), the timeout, the exit code, the
    seconds, and the priority classes the watcher SAW the program (and anything it started) at. What the program printed is kept in
    <cwd>\\<name>.out.txt and .err.txt. Returns the SolverRun."""
    import checks
    with log.tee():
        with watch_children() as watcher:
            result = femm_queue.run_solver(cmd, label=label, cwd=cwd, timeout_s=timeout_s, threads=threads)
    checks.save_log(Path(cwd), name, result)
    configured = paths.setting("femm", "workers")
    log.add_run({"kind": "run_solver", "label": label, "program": Path(str(cmd[0])).name, "threads": threads, "threads_limit": configured,
                 "timeout_s": timeout_s, "returncode": result.returncode, "timed_out": result.timed_out, "seconds": round(result.seconds, 3),
                 "priority_classes_seen": watcher.summary(), "looks_at_the_processes": watcher.looks,
                 "time": time.strftime("%Y-%m-%d %H:%M:%S")})
    return result
