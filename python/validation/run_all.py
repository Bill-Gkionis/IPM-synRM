r"""run_all.py: re-run every script of this folder on FEMM 5.0, one at a time, and keep a table of what happened.

Why it exists
    The motor scripts were written and run on FEMM 4.2. FEMM 5.0 is another program with another Python binding
    (femm_session.py), so a script could break, or a number could move, without anybody noticing. This runner starts
    every script again, in a fixed order, from the folder it expects and with the arguments it needs. It keeps one log
    per script and a table (run_all.csv) with how each one ended. It does not judge the numbers: every script prints
    its own, and the README's "Expected numbers" and the FEMM 4.2 outputs in femm42-out (the reference: a checked copy of
    the `out` folder that FEMM 4.2 left in D:\dev\scripts\automation\proc\out) are what they are compared with afterwards.

Run (from the Python to FEMM folder; one runner at a time, and nothing else heavy on FEMM meanwhile)
    uv run python validation/run_all.py                 every step, in order: about 9 minutes on FEMM 5.0
    uv run python validation/run_all.py --list          the steps, what each one needs from earlier ones; runs nothing
    uv run python validation/run_all.py --only check_geometry,check_vs_fem,check_femm        a few (the smoke test)
    uv run python validation/run_all.py --from validate_full                 that step and every step after it
    uv run python validation/run_all.py --only barrier_study --timeout-factor 10     more patience for a slow step
  --only and --from always run in table order, whatever order is typed. A step that needs an earlier step's files
  (--list shows which) finds them only if that step ran before, in this run or an older one.
  Exit status: 0 if every step that ran has status ok. 1 if one did not. 2 if the runner could not start at all
  (a name it does not know, FEMM_DIR not found, or another runner already running).

How one step is run
    One step is one child process, started in the step's own folder with the environment's own python.exe
    (sys.executable). That is what `uv run python` runs, without uv's start-up check, so uv cannot change the
    environment in the middle of a run. Steps run strictly one after the other. Each child gets this environment:
      FEMM_DIR           from femm_session.femm_dir(), so a terminal that predates FEMM 5.0 works too;
      MOTOR_FIGURES_DIR  out\materials\figures, so the materials scripts keep their figures out of the repository's docs figures folder;
      PYTHONUNBUFFERED=1 so a step's log grows while it runs (watch it with Get-Content -Wait);
      PYTHONIOENCODING=utf-8 so the logs are UTF-8 whatever the console's code page is.
    The runner puts itself at below-normal priority first (femm_session.lower_priority()). Every child, and every
    FEMM a child starts, inherits that, so a long run does not take the desktop away. It also switches off Windows'
    "has stopped working" window (SetErrorMode), which the children inherit too: a crashing solver would otherwise wait
    for a click nobody is there to give. The crash is still written to the Application event log. The step's standard
    output and standard error go straight into out\validation\logs\NN-name.log (NN = the step's number); a step run
    again keeps the log of the run before as NN-name.log.prev.
    While a step runs the runner looks, twice a second, at which processes are below it, and remembers each as
    (id, program, creation time). That memory is the only thing it ever acts on: see "Whom the runner ends" below.
    Stopped by Ctrl-C, a step is ended with everything seen below it. If the runner is killed outright, Windows kills
    the steps too: the runner puts itself into a "job object" that ends every process in it when the runner ends, and
    every step, worker and FEMM is born inside it (join_kill_on_exit_job).

Time limits
    A step that runs longer than max(5 min, 5 x the time it took in the first FEMM 5.0 sweep, 2026-11-04) is ended. Its
    row says `timeout` and the run goes on with the next step. --timeout-factor changes the 5. (FEMM 4.2 took 2 h 13
    min for all the steps and FEMM 5.0 took 9, so the limits are in minutes, not hours; --list shows both times.)

Whom the runner ends
    Only processes it saw below the step while it ran, each only after its id, its program and its creation time have
    been checked again on the very handle it ends (an id is given out again after its owner ends; a creation time is
    not), newest first. That covers a step over its time limit, Ctrl-C, and a FEMM the step left running after it ended.
    Never touched, whatever they are: a process that existed before the step; a FEMM whose parent the runner never saw
    (a FEMM from somebody's debugging session that has since exited looks exactly like that); a new FEMM that hangs below
    a process outside the step, such as the fkn.exe that an interactive FEMM window starts when its user presses Analyze.

How to read run_all.csv (out\validation\run_all.csv; one row per step, rewritten after every step)
    step               the step's position in the table (and the NN of its log file)
    name               what --only and --from take
    command            what ran, as it would be typed after `cd <cwd>` (the runner uses python.exe, not uv)
    exit_code          the script's exit status. 0 = it finished normally. Anything else = it stopped with an error.
                       (A step the runner ended shows whatever that left, usually 1: read `status`.) Blank: not run.
    seconds            wall-clock time of the step. Does not include the up to 5 s the runner waits for FEMM to go.
    femm_left_running  how many FEMM programs (femm_headless.exe, femm_headlessw.exe, fkn.exe) that were NOT running
                       before the step are still running 5 s after it ends. Must be 0. Counted: the ones the runner saw
                       below the step (it then ends them, so one leak cannot disturb the next step), and new ones that
                       look like the step's but it never saw (orphans, or children of a process of the step), which it
                       only reports. Never counted: see "Whom the runner ends".
    cwd                the folder the step ran in, relative to Python to FEMM
    status             ok = exit 0 and nothing left running | failed = exit code not 0 | timeout = ended for taking too
                       long | leak = exit 0 but FEMMs still running | interrupted = Ctrl-C | runner-error = a bug in
                       this runner (the traceback is on the console) | not-run = chosen for this run, not run (yet)
    finished           when the row was written
  A step is good when status is ok. The file keeps the latest row of every step. A run starts by marking every step it
  will run `not-run`, so a run that stops early (Ctrl-C, a crash of the machine) leaves an honest table, never the
  "ok" of an older run. A partial run (--only, --from) replaces only its own rows, so after fixing a step and running
  it again with --only, the file can read "all ok" again. Rows of names that are no longer in the table are dropped.

The order, and what each step needs
    This is the plan's list with two changes. (1) barrier_study comes before reference_open: reference_open saves its
    models into out\barrier\fem, a folder that only barrier_study makes, and it fails at once without it. (2) A new
    step, hand_reference_check, follows validate_full (see below). Inputs that come from earlier steps:
      reference_open        out\barrier\fem (the folder)                      <- barrier_study
      hand_reference_check  out\validate\hand_000.fem                         <- validate_full
      validate_fix_zeros    out\validate\validate.csv                         <- validate_full
      barrier_plots         out\barrier\*.csv and fields\*.ans                <- barrier_study
                            out\validate\validate.csv (one optional figure)   <- validate_fix_zeros
      optimize_plot         out\opt\evaluations.csv, pareto_X.npy             <- optimize_demo
      flux_explorer_sweep120 out\flux_explorer\fem (the folder)               <- flux_explorer_sweep
      flux_explorer_build   out\flux_explorer\states.csv, *.npz               <- flux_explorer_sweep
                            out\flux_explorer\states_p120.csv, *_p120.npz     <- flux_explorer_sweep120
      pf_harmonics          out\barrier\*.csv                                 <- barrier_study
                            out\barrier\reference_open.csv                    <- reference_open
                            out\flux_explorer\states.csv                      <- flux_explorer_sweep
      materials_fields      out\barrier\fields\*.ans                          <- barrier_study
      materials_magnetic    out\materials\fields.json, current.json           <- materials_fields, materials_current
    If a step fails, the steps that read its files are still run (they will probably fail too); the summary says
    "its input comes from X, which also failed" so one cause is not mistaken for several.
    Outside out\ the steps read: the hand-drawn .FEM (check_vs_fem, validate_full, hand_reference_check and materials_magnetic:
    only read, never saved; validate_full works on copies), FEMM's matlib.dat and materials\research\m19_dc_maker.csv.
    Outside out\ they write: learn\ (tiny.dxf, magnet.fem, magnet.ans, out\: git ignores them), runs\ (evaluate's private
    files: an ordinary folder of Python to FEMM, which git ignores) and pages\Rotor Flux Explorer.html
    (flux_explorer_build), which git keeps as a binary file.
    Nothing is written into the repository's docs figures folder. After a run, check that it shows no change under git:
      git -C "<repository folder>" status --short -- docs/figures

What each step proves
    The numbers are the README's, FEMM 5.0 ones where it has them. "x42" is the FEMM 4.2 time.
     1 check_geometry   x42 ~15 s. The geometry kernel works without FEMM: the baseline sector has 40 regions and no rib
                        under 1 mm, kw1 = 0.9452, the DXF round trip is within +0.030 %, Gmsh meshes it, and two
                        impossible candidates (v_half_deg 55, mag_x 38.4) are refused before any solver runs.
     2 check_vs_fem     x42 ~2 s. motorgeo draws the hand-drawn model: rotor steel, bands and magnets 0.000 %, stator
                        steel -0.009 %, slot 54.127 against 54.152 mm2, winding the mirror of the .FEM's.
     3 check_femm       x42 2 s. FEMM 5.0 solves the baseline model, coarse mesh (9,580 nodes): no-load psi_d / psi_q
                        11.991 / 0.031 mWb (d on the magnet), T_age = T_wst = 0.8087 / 0.8085 N.m at iq = 10 A, and
                        0.9104 (id -6, iq 8 A) above 0.3944 N.m (+6, 8 A): Ld < Lq.
     4 check_femm_fine  x42 3 s. The same on the fine mesh (17,934 nodes): 12.032 / 0.031 mWb, 0.8082 N.m, 0.9203, 0.3826.
     5 example          x42 ~11 s. The one-screen tour (sector plot, DXF, SolidWorks equations, Gmsh mesh) and 30 rotor
                        positions: mean 0.729 N.m, ripple 0.158 N.m peak to peak.
     6 evaluate         x42 ~6 s. One optimizer candidate = one FEMM in a pool worker: 5 candidates on 4 workers; the
                        baseline T_mean 0.7289 N.m, T_ripple 0.1579; mag_y 5.0 gives 1.0013; web 1.2 gives 0.6635; the
                        last two candidates are refused without FEMM.
     7 mesh_study       x42 ~2 min. The torque over one 60 degree ripple period on 6 workers, coarse against fine
                        mesh, with the 6th/12th/18th/24th harmonics.
     8 worker_scaling   x42 ~4 min. The same 48 solves on 1, 2, 4, 6, 8 and 10 workers: bit-identical results
                        (max |T - T(1 worker)| = 0) and a speed-up of about 5.9x on 8 workers (x42), 6.4x on 10.
     9 barrier_study    x42 ~60 min (10 workers; the script runs 8). The air-band study on the fine mesh, 1,070 solves:
                        no load, cogging, torque against current angle at 4 currents, ripple at MTPA, field
                        snapshots. At 10 A and its own MTPA angle: with the band psi_PM 12.09 mWb, 0.857 N.m;
                        without it 16.57 mWb, 1.125 N.m (x42).
    10 reference_open   x42 ~7 min. The magnet inside the band (the 10-22 prototype): 19.18 mWb, 1.266 N.m (x42).
    11 validate_full    x42 ~30 min. The 90 degree slice is the same machine as the hand-drawn full 360 degree model: the flux
                        linkage differs by at most 1.13 % (the worst rotor position, as a share of the full model's
                        peak) and the peak torque at 10 A by 1.43 % (the slice's peak against the full model's)
                        (FEMM 5.0; x42: 1.15 % and 1.48 %). 96 solves, 8 workers, 27 s. Every job has files of its
                        own: until 2026-11-04 the no-load and the 10 A slice jobs shared four file names, which
                        crashed, hung or corrupted runs of this step.
    12 hand_reference_check  (new, seconds). FEMM 5.0's mi_moverotate turns each magnet's direction with the rotor, as 4.2's
                        did. validate_full's copy hand_000.fem was made by FEMM 5.0 from the hand-drawn .FEM with three changes
                        (air-gap label to group 0, rotor turned, currents set). This applies the same three changes to
                        the .FEM as plain arithmetic and compares. Exit 0 = the same model. Anything else means the
                        rotated copies are not the model validate_full thinks they are, so its numbers are wrong:
                        the step fails, and with it the run.
    13 validate_fix_zeros  x42 ~2 min. Re-solves the slice rows whose flux linkage came back as exactly 0.0 (a 4.2
                        glitch that was the shared file names above). On FEMM 5.0 it finds none ("0 rows to re-solve").
    14 barrier_plots    x42 ~1 min. The figures and summary.csv of the air-band study, from the study's own files.
    15 optimize_demo    x42 ~12 min. NSGA-II end to end (pymoo, motorgeo, shapely, FEMM, CSV): 8 candidates, 7
                        generations. It checks itself ("pymoo's Pareto scores match the evaluation log") and sends
                        the knee design back to a DXF (round trip within +0.030 % at 4.2).
    16 optimize_plot    seconds. The Pareto figure and the knee design's outline, from the demo's log.
    17 flux_explorer_sweep  x42 ~12 min. The fields for the interactive page: 34 models, 2,040 solves, 8 workers; the
                        script reports "flux-linkage retries: 0".
    18 flux_explorer_sweep120  x42 ~2 min. The 120 degree end position for the same fields: 306 solves.
    19 flux_explorer_build  seconds. Packs the sweep into pages\Rotor Flux Explorer.html (about 17 MB); it asserts that
                        every file of a rotor has the same mesh. The page's footer sentence, which says which FEMM
                        made the solves and when, is not typed into the script: since 708a98f the builder reads the
                        version from the installed FEMM and the day from the sweep's tables, so a re-run on another
                        day rewrites the tracked page with a new date in its footer.
    20 pf_harmonics     x42 ~10 s. Power factor and field harmonics from the study data on disk (no new solves); its
                        check T = 1.5 p |psi| |i| PF against FEMM's torque must agree to ~0.
    21 materials_fields, 22 materials_current, 23 materials_magnetic, 24 materials_thermal, 25 materials_structural,
    26 materials_crosssection  (seconds each; BUILD.md's order).
                        Every number and figure of the materials notes can be recomputed. Only fields (and the F3/F4
                        parts of magnetic) read FEMM results, from barrier_study's saved .ans files; current, thermal
                        and structural read none and should equal FEMM 4.2 era values. magnetic.py reads FEMM 5.0's
                        library for the ASM "M-19" curve and checks that the model's "M-19 Steel" is FEMM 4.2's.
    27 learn_01_shapely_basics  Prints "faces found: 2 areas: [8.0, 8.0]" and the other shapely answers. No FEMM.
    28 learn_02_ezdxf_basics    Writes learn\tiny.dxf and reads it back. No FEMM.
    29 learn_03_pyfemm_magnet   One magnet in a box of air, solved by its own FEMM on its own port: B above and inside
                        the magnet, the node count.
    30 learn_04_workers         Four process workers get four FEMMs (it counts femm_headlessw.exe on the whole PC, so run
                        it when no other headless FEMM runs); threads are refused a second FEMM; a shared file name
                        gives a wrong answer, private names give the right ones.
    31 learn_05_pymoo_toy       NSGA-II on the two-cooks toy problem; no FEMM.
    32 learn_geometry_figure, 33 learn_toy_front_figure   The two figures of the project notes, written to out\learn.

Self-check
    1. Why is the table's `femm_left_running` measured with a list of process ids taken before the step, and not by
       asking "is any femm_headlessw.exe running?" after it?
    2. A step times out. Which processes does the runner end, and which does it leave alone?
    3. barrier_plots was fixed and its row should read ok without a full run. What is typed? (Its inputs are
       in out\ from the last run.)
    4. A run is stopped half way with Ctrl-C. What do the rows of the steps that had not run yet say?
    Answers: 1. An interactive FEMM window may be open; its processes were there before the step, so they are not new and
    are never counted. 2. The step's python.exe and everything the runner saw below it while it ran (pool workers, their
    FEMMs), each after checking its id, program and creation time again. Nothing that existed before the step, and no
    FEMM it did not see below the step. 3. `uv run python validation/run_all.py --only barrier_plots`. 4. `not-run`,
    whatever an older run left there.
"""
from __future__ import annotations

import argparse
import csv
import ctypes
import datetime
import difflib
import os
import subprocess
import sys
import time
import traceback
from ctypes import wintypes
from dataclasses import dataclass
from pathlib import Path

P = Path(__file__).resolve().parents[1]                 # Python to FEMM
OUT = P / "out"                                         # a plain folder in P (git ignores it); FEMM 4.2's run is P/femm42-out
LOGS = OUT / "validation" / "logs"
LEDGER = OUT / "validation" / "run_all.csv"
LOCK = OUT / "validation" / "run_all.lock"
FIGURES = OUT / "materials" / "figures"                 # MOTOR_FIGURES_DIR for every child

FEMM_IMAGES = ("femm_headless.exe", "femm_headlessw.exe", "fkn.exe")   # the programs femm_left_running counts
MIN_TIMEOUT = 300.0                                     # seconds: no step gets less than 5 minutes
TIMEOUT_FACTOR = 5.0                                    # times the seconds the step took in the 2026-11-04 FEMM 5.0 sweep
GRACE = 5.0                                             # seconds a finished step gets for its FEMMs to quit
COLUMNS = ("step", "name", "command", "exit_code", "seconds", "femm_left_running", "cwd", "status", "finished")


# ================================================================================================ the step table
@dataclass(frozen=True)
class Step:
    """One script run. Everything the runner knows about it is here, and nowhere else."""
    name: str                                       # what --only and --from take; part of the log's file name
    script: str                                     # relative to `folder`
    args: tuple[str, ...] = ()
    folder: str = "."                               # working folder, relative to P: the scripts' relative paths need it
    needs: tuple[tuple[str, str], ...] = ()         # (earlier step that makes it, what it makes) for each input
    femm42: str = ""                                # the README's FEMM 4.2 time, for people ("~12 min")
    measured: int | None = None                     # seconds the step took in the 2026-11-04 FEMM 5.0 sweep (rounded up)
    mkdirs: tuple[str, ...] = ()                    # folders (relative to P) the script expects to exist already


STEPS: tuple[Step, ...] = (
    Step("check_geometry", "check_geometry.py", femm42="~15 s", measured=2),
    Step("check_vs_fem", "check_vs_fem.py", femm42="~2 s", measured=2),
    Step("check_femm", "check_femm.py", femm42="2 s", measured=2),
    Step("check_femm_fine", "check_femm.py", ("--fine",), femm42="3 s", measured=3),
    Step("example", "example.py", femm42="~11 s", measured=10),
    Step("evaluate", "evaluate.py", femm42="~6 s", measured=6),
    Step("mesh_study", "mesh_study.py", femm42="~2 min", measured=7),
    Step("worker_scaling", "worker_scaling.py", (r"out\scaling\worker_scaling.png",), femm42="~4 min", measured=34),
    Step("barrier_study", "barrier_study.py", femm42="~60 min", measured=131),
    Step("reference_open", "reference_open.py", femm42="~7 min", measured=13,
         needs=(("barrier_study", r"out\barrier\fem (the folder)"),)),
    Step("validate_full", "validate_full.py", femm42="~30 min", measured=27),
    Step("hand_reference_check", r"validation\hand_reference_check.py", (r"out\validate\hand_000.fem", "--job", "0"),
         femm42="new, seconds", measured=1, needs=(("validate_full", r"out\validate\hand_000.fem"),)),
    Step("validate_fix_zeros", "validate_fix_zeros.py", femm42="~2 min", measured=1,
         needs=(("validate_full", r"out\validate\validate.csv"),)),
    Step("barrier_plots", "barrier_plots.py", (r"out\barrier\figs",), femm42="~1 min", measured=5,
         needs=(("barrier_study", r"out\barrier\phase1.csv, phase2.csv, gap_*.csv, fields\*.ans"),
                ("validate_fix_zeros", r"out\validate\validate.csv (optional: one figure)"))),
    Step("optimize_demo", "optimize_demo.py", femm42="~12 min", measured=99),
    Step("optimize_plot", "optimize_plot.py", (r"out\opt\pareto.png",), femm42="seconds", measured=1,
         needs=(("optimize_demo", r"out\opt\evaluations.csv, pareto_X.npy"),)),
    Step("flux_explorer_sweep", "flux_explorer_sweep.py", femm42="~12 min", measured=119),
    Step("flux_explorer_sweep120", "flux_explorer_sweep120.py", femm42="~2 min", measured=20,
         needs=(("flux_explorer_sweep", r"out\flux_explorer\fem (the folder)"),)),
    Step("flux_explorer_build", "flux_explorer_build.py", femm42="seconds", measured=15,
         needs=(("flux_explorer_sweep", r"out\flux_explorer\states.csv, g*_air.npz, g*_steel.npz"),
                ("flux_explorer_sweep120", r"out\flux_explorer\states_p120.csv, g*_p120.npz"))),
    Step("pf_harmonics", "pf_harmonics.py", (r"out\pf_harmonics",), femm42="~10 s", measured=2,
         needs=(("barrier_study", r"out\barrier\phase1.csv, phase2.csv, gap_*.csv"),
                ("reference_open", r"out\barrier\reference_open.csv"),
                ("flux_explorer_sweep", r"out\flux_explorer\states.csv"))),
    Step("materials_fields", r"materials\fields.py", femm42="seconds", measured=2,
         needs=(("barrier_study", r"out\barrier\fields\*.ans"),)),
    Step("materials_current", r"materials\current.py", femm42="seconds", measured=2),
    Step("materials_magnetic", r"materials\magnetic.py", femm42="seconds", measured=2,
         needs=(("materials_fields", r"out\materials\fields.json"),
                ("materials_current", r"out\materials\current.json"))),
    Step("materials_thermal", r"materials\thermal.py", femm42="seconds", measured=1),
    Step("materials_structural", r"materials\structural.py", femm42="seconds", measured=1),
    Step("materials_crosssection", r"materials\crosssection.py", femm42="seconds", measured=2),
    Step("learn_01_shapely_basics", "01_shapely_basics.py", folder="learn", femm42="seconds", measured=1),
    Step("learn_02_ezdxf_basics", "02_ezdxf_basics.py", folder="learn", femm42="seconds", measured=1),
    Step("learn_03_pyfemm_magnet", "03_pyfemm_magnet.py", folder="learn", femm42="seconds", measured=1),
    Step("learn_04_workers", "04_workers.py", folder="learn", femm42="seconds", measured=17),
    Step("learn_05_pymoo_toy", "05_pymoo_toy.py", folder="learn", femm42="seconds", measured=1),
    Step("learn_geometry_figure", "geometry_figure.py", (r"..\out\learn\geometry_figure.png",), folder="learn",
         femm42="seconds", measured=2, mkdirs=(r"out\learn",)),
    Step("learn_toy_front_figure", "toy_front_figure.py", (r"..\out\learn\toy_front.png",), folder="learn",
         femm42="seconds", measured=2, mkdirs=(r"out\learn",)),
)


def command_text(step: Step) -> str:
    """The step as it would be typed after `cd <folder>`."""
    return "python " + subprocess.list2cmdline([step.script, *step.args])


def timeout_for(step: Step, factor: float = TIMEOUT_FACTOR) -> float:
    """Seconds a step may run: `factor` times the seconds it took in the first FEMM 5.0 sweep, never less than 5 min."""
    return max(MIN_TIMEOUT, factor * (step.measured or 0))


def select_steps(only: list[str] | None, from_step: str | None) -> list[Step]:
    """The steps to run, in table order. Raises ValueError, naming the nearest valid names, for an unknown name."""
    names = [s.name for s in STEPS]
    wanted = [n.strip() for chunk in (only or []) for n in chunk.split(",") if n.strip()]
    if only is not None and not wanted:
        raise ValueError("--only needs at least one step name")
    for name in wanted + ([from_step] if from_step else []):
        if name not in names:
            near = difflib.get_close_matches(name, names, n=3)
            raise ValueError(f"no step called {name!r}" + (f" (closest names: {', '.join(near)})" if near else "")
                             + ". --list shows the names.")
    if from_step:
        return list(STEPS[names.index(from_step):])
    if wanted:
        return [s for s in STEPS if s.name in wanted]
    return list(STEPS)


# ================================================================================================ processes (Windows)
# Windows keeps a list of running processes. Taking a "snapshot" of it is one call; each row has the process id,
# the id of its parent and the name of its program. A process id is NOT a name for ever: Windows gives an id to a new
# process as soon as the old one is gone, and the list keeps naming the old parent in its children. So whenever the
# runner must be sure which process it means, it uses an identity: (program, creation time). Two processes never
# share one, whatever happens to the ids.
_TH32CS_SNAPPROCESS = 0x2
_KILL_ON_JOB_CLOSE = 0x2000
_JOB_OBJECT_EXTENDED_LIMIT_INFORMATION = 9
_PROCESS_TERMINATE, _PROCESS_SET_QUOTA, _PROCESS_QUERY_LIMITED_INFORMATION = 0x1, 0x100, 0x1000
_SEM_FAILCRITICALERRORS, _SEM_NOGPFAULTERRORBOX = 0x1, 0x2

Identity = tuple[str, int]              # (program name in lower case, creation time as a Windows FILETIME integer)


class _ProcessEntry(ctypes.Structure):
    """PROCESSENTRY32W from the Windows API: one row of the list of running processes."""
    _fields_ = [("dwSize", wintypes.DWORD), ("cntUsage", wintypes.DWORD), ("th32ProcessID", wintypes.DWORD),
                ("th32DefaultHeapID", ctypes.c_size_t), ("th32ModuleID", wintypes.DWORD),
                ("cntThreads", wintypes.DWORD), ("th32ParentProcessID", wintypes.DWORD),
                ("pcPriClassBase", wintypes.LONG), ("dwFlags", wintypes.DWORD),
                ("szExeFile", wintypes.WCHAR * 260)]


class _FileTime(ctypes.Structure):
    """FILETIME: a moment, as a count of 100 ns steps split over two 32-bit halves."""
    _fields_ = [("dwLowDateTime", wintypes.DWORD), ("dwHighDateTime", wintypes.DWORD)]


class _IoCounters(ctypes.Structure):
    _fields_ = [(n, ctypes.c_ulonglong) for n in ("ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
                                                  "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class _BasicLimits(ctypes.Structure):
    _fields_ = [("PerProcessUserTimeLimit", ctypes.c_longlong), ("PerJobUserTimeLimit", ctypes.c_longlong),
                ("LimitFlags", wintypes.DWORD), ("MinimumWorkingSetSize", ctypes.c_size_t),
                ("MaximumWorkingSetSize", ctypes.c_size_t), ("ActiveProcessLimit", wintypes.DWORD),
                ("Affinity", ctypes.c_size_t), ("PriorityClass", wintypes.DWORD), ("SchedulingClass", wintypes.DWORD)]


class _ExtendedLimits(ctypes.Structure):
    """JOBOBJECT_EXTENDED_LIMIT_INFORMATION: a job's limits. The only one used here is "kill all on close"."""
    _fields_ = [("BasicLimitInformation", _BasicLimits), ("IoInfo", _IoCounters), ("ProcessMemoryLimit", ctypes.c_size_t),
                ("JobMemoryLimit", ctypes.c_size_t), ("PeakProcessMemoryUsed", ctypes.c_size_t),
                ("PeakJobMemoryUsed", ctypes.c_size_t)]


_k32 = None


def _kernel32():
    """kernel32.dll with the argument and result types of the few calls used here (a 64-bit handle is not an int)."""
    global _k32
    if _k32 is None:
        k = ctypes.WinDLL("kernel32", use_last_error=True)
        k.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
        k.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
        k.Process32FirstW.argtypes = [wintypes.HANDLE, ctypes.POINTER(_ProcessEntry)]
        k.Process32FirstW.restype = wintypes.BOOL
        k.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(_ProcessEntry)]
        k.Process32NextW.restype = wintypes.BOOL
        k.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
        k.OpenProcess.restype = wintypes.HANDLE
        k.CloseHandle.argtypes = [wintypes.HANDLE]
        k.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(_FileTime)] * 4
        k.GetProcessTimes.restype = wintypes.BOOL
        k.QueryFullProcessImageNameW.argtypes = [wintypes.HANDLE, wintypes.DWORD, wintypes.LPWSTR,
                                                 ctypes.POINTER(wintypes.DWORD)]
        k.QueryFullProcessImageNameW.restype = wintypes.BOOL
        k.TerminateProcess.argtypes = [wintypes.HANDLE, wintypes.UINT]
        k.TerminateProcess.restype = wintypes.BOOL
        k.SetErrorMode.argtypes = [wintypes.UINT]
        k.SetErrorMode.restype = wintypes.UINT
        k.CreateJobObjectW.argtypes = [ctypes.c_void_p, wintypes.LPCWSTR]
        k.CreateJobObjectW.restype = wintypes.HANDLE
        k.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, ctypes.c_void_p, wintypes.DWORD]
        k.SetInformationJobObject.restype = wintypes.BOOL
        k.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]
        k.AssignProcessToJobObject.restype = wintypes.BOOL
        _k32 = k
    return _k32


def running_processes() -> dict[int, tuple[int, str]]:
    """{process id: (parent's id, program name in lower case)} for every process running now."""
    k = _kernel32()
    snapshot = k.CreateToolhelp32Snapshot(_TH32CS_SNAPPROCESS, 0)
    if snapshot is None or snapshot == ctypes.c_void_p(-1).value:           # INVALID_HANDLE_VALUE
        raise ctypes.WinError(ctypes.get_last_error())
    table = {}
    try:
        entry = _ProcessEntry()
        entry.dwSize = ctypes.sizeof(entry)
        more = k.Process32FirstW(snapshot, ctypes.byref(entry))
        while more:
            table[entry.th32ProcessID] = (entry.th32ParentProcessID, entry.szExeFile.lower())
            more = k.Process32NextW(snapshot, ctypes.byref(entry))
    finally:
        k.CloseHandle(snapshot)
    return table


def _identity_of(handle) -> Identity | None:
    """(program, creation time) of the process behind an open handle, or None if Windows will not say."""
    k = _kernel32()
    made, ended, kernel, user = (_FileTime() for _ in range(4))
    if not k.GetProcessTimes(handle, ctypes.byref(made), ctypes.byref(ended), ctypes.byref(kernel), ctypes.byref(user)):
        return None
    buffer = ctypes.create_unicode_buffer(1024)
    size = wintypes.DWORD(len(buffer))
    if not k.QueryFullProcessImageNameW(handle, 0, buffer, ctypes.byref(size)):
        return None
    return os.path.basename(buffer.value).lower(), (made.dwHighDateTime << 32) | made.dwLowDateTime


def process_identity(pid: int) -> Identity | None:
    """The identity of the process that has this id right now, or None if there is none (or Windows will not say)."""
    k = _kernel32()
    handle = k.OpenProcess(_PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
    if not handle:
        return None
    try:
        return _identity_of(handle)
    finally:
        k.CloseHandle(handle)


def terminate_if_same(pid: int, identity: Identity) -> bool:
    """End the process with this id, but only if it is still the one with this identity. Returns whether it was ended.

    The identity is read from the same handle that is then used to end the process, so there is no gap in which the
    id could be handed to another process between the check and the kill.
    """
    k = _kernel32()
    handle = k.OpenProcess(_PROCESS_TERMINATE | _PROCESS_QUERY_LIMITED_INFORMATION, False, pid)
    if not handle:
        return False
    try:
        if _identity_of(handle) != identity:
            return False
        return bool(k.TerminateProcess(handle, 1))
    finally:
        k.CloseHandle(handle)


def femm_pids(table: dict[int, tuple[int, str]]) -> set[int]:
    """The ids of the FEMM programs (femm_headless.exe, femm_headlessw.exe, fkn.exe) in a process table."""
    return {pid for pid, (_, image) in table.items() if image in FEMM_IMAGES}


def descendants(table: dict[int, tuple[int, str]], root: int, created=None) -> set[int]:
    """Every process below `root` in the parent-child tree: its children, theirs, and so on.

    The table names each process's parent by id, and an id can be handed to a new process after its owner ended. If
    `created` (a function: process id -> creation time or None) is given, a process counts as a child only if it was
    created after the parent it names; a link that fails this is a stale one and is ignored.
    """
    children: dict[int, list[int]] = {}
    for pid, (parent, _) in table.items():
        children.setdefault(parent, []).append(pid)
    found, todo = set(), [root]
    while todo:
        parent = todo.pop()
        for child in children.get(parent, []):
            if child in found or child == root:
                continue
            if created is not None:
                born_parent, born_child = created(parent), created(child)
                if born_parent is None or born_child is None or born_child < born_parent:
                    continue
            found.add(child)
            todo.append(child)
    return found


def forget_gone(own: dict[int, Identity], table: dict[int, tuple[int, str]]) -> None:
    """Drop from `own` every process that has ended, or whose id now belongs to another process."""
    for pid, identity in list(own.items()):
        if pid not in table or process_identity(pid) != identity:
            del own[pid]


def look_below(own: dict[int, Identity], table: dict[int, tuple[int, str]], root: int) -> None:
    """Remember in `own` every process now below `root` (the step's python.exe); forget the ones that are gone."""
    seen: dict[int, Identity | None] = {}

    def identity(pid):
        if pid not in seen:
            seen[pid] = process_identity(pid)
        return seen[pid]

    def created(pid):
        found = identity(pid)
        return found[1] if found else None
    for pid in descendants(table, root, created):
        if pid not in own and identity(pid) is not None:
            own[pid] = identity(pid)
    forget_gone(own, table)


def belongs_to_step(pid: int, table: dict[int, tuple[int, str]], own) -> bool:
    """Does this process look like the step's? `own` holds the ids of the step's python.exe and what was seen below it.

    Yes if its parent has gone (an orphan: whatever started it ended and left it behind), or if any of its
    ancestors is in `own`. No if it hangs below a process that is still running outside the step.
    """
    parent = table[pid][0]
    if parent not in table:
        return True
    guard = {pid}
    while parent in table and parent not in guard:
        if parent in own:
            return True
        guard.add(parent)
        parent = table[parent][0]
    return parent in own                                # the chain ended at a process that has gone: was it ours?


def leaked_femms(before: set[int], table: dict[int, tuple[int, str]], own: dict[int, Identity],
                 identity=process_identity) -> tuple[dict[int, Identity], set[int]]:
    """The FEMM programs a step left running, as (ours, others). Both are new since `before`.

    ours: seen below the step while it ran and still the very same process. These are the only FEMMs the runner ends.
    others: never seen below the step, but they look like its: the parent has gone (orphans), or an ancestor is a
      process of the step. They are counted and reported, never touched: a FEMM from somebody's debugging session
      that has since exited looks exactly like this too.
    A new FEMM that hangs below a process outside the step (the fkn.exe of an interactive FEMM window) is neither.
    """
    new = femm_pids(table) - before
    ours = {pid: own[pid] for pid in new if pid in own and identity(pid) == own[pid]}
    others = {pid for pid in new - set(ours) if belongs_to_step(pid, table, own)}
    return ours, others


def kill_step(proc: subprocess.Popen, own: dict[int, Identity]) -> None:
    """End a step that must stop: its python.exe and everything remembered below it, and nothing else.

    Each process is ended only after its identity has been checked again (terminate_if_same), newest first: a child is
    always newer than its parent, so the children go before the process that started them. That order does not stop a
    parent that is still alive from starting one more child while the others are being ended, so it is not the net.
    The net is the repetition: the round is repeated, with a fresh look at what is below the step, until nothing of
    it is left.
    """
    deadline = time.monotonic() + 20
    while own and time.monotonic() < deadline:
        try:
            look_below(own, running_processes(), proc.pid)
        except OSError:
            pass                                        # cannot read the list: end what is remembered
        for pid, identity in sorted(own.items(), key=lambda item: item[1][1], reverse=True):
            terminate_if_same(pid, identity)
        time.sleep(0.2)
        try:
            forget_gone(own, running_processes())
        except OSError:
            pass
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        proc.kill()                                     # our own child, through the handle Popen has held since it began


_job = None                     # the job object's handle; it stays open until this process ends


def join_kill_on_exit_job() -> None:
    """Put THIS process into a Windows job object that ends every process in it when the process ends.

    A job object is a group of processes that Windows can treat as one. This one has a single rule, "kill them all when
    the last handle to the job closes", and the only handle is ours, so the group dies exactly when this runner does,
    however it ends (finished, Ctrl-C, or killed from outside). A process belongs to its parent's job from the moment
    it is born, so every step, pool worker and FEMM started from here on is in the group, with no gap in which one
    could slip out. Without it, killing the runner would leave a 60-minute study running on its own.
    Raises OSError if Windows refuses (the runner is then on its own).
    """
    global _job
    k = _kernel32()
    job = k.CreateJobObjectW(None, None)
    if not job:
        raise ctypes.WinError(ctypes.get_last_error())
    limits = _ExtendedLimits()
    limits.BasicLimitInformation.LimitFlags = _KILL_ON_JOB_CLOSE
    if not k.SetInformationJobObject(job, _JOB_OBJECT_EXTENDED_LIMIT_INFORMATION, ctypes.byref(limits),
                                     ctypes.sizeof(limits)):
        raise ctypes.WinError(ctypes.get_last_error())
    me = k.OpenProcess(_PROCESS_TERMINATE | _PROCESS_SET_QUOTA, False, os.getpid())
    if not me:
        raise ctypes.WinError(ctypes.get_last_error())
    try:
        if not k.AssignProcessToJobObject(job, me):
            raise ctypes.WinError(ctypes.get_last_error())
    finally:
        k.CloseHandle(me)
    _job = job


def suppress_crash_dialogs() -> None:
    """Switch off Windows' "has stopped working" window for this process, and so for every step (children inherit it).

    A solver that crashes, as FEMM 5.0's fkn.exe did on 2026-11-04, would otherwise wait for a click that nobody is
    there to give, and the run would wait with it until its time limit. The crash is still logged by Windows
    (Event Viewer, Application, "Application Error"). The flags are added to whatever the process had.
    """
    k = _kernel32()
    k.SetErrorMode(k.SetErrorMode(0) | _SEM_FAILCRITICALERRORS | _SEM_NOGPFAULTERRORBOX)


# ================================================================================================ one step
@dataclass
class Result:
    number: int
    step: Step
    exit_code: int
    seconds: float
    femm_left: int
    status: str                 # ok | failed | timeout | leak | interrupted | runner-error
    log: Path
    finished: str = ""


def now_text() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def fmt_duration(seconds: float | None) -> str:
    if seconds is None:
        return "-"
    if seconds < 60:
        return f"{seconds:.1f} s" if seconds != int(seconds) else f"{int(seconds)} s"
    seconds = int(round(seconds))
    hours, rest = divmod(seconds, 3600)
    minutes, secs = divmod(rest, 60)
    return f"{hours} h {minutes:02d} min" if hours else f"{minutes} min {secs:02d} s"


def keep_previous_log(path: Path) -> None:
    """A step run again keeps the log of the run before it as NN-name.log.prev (one generation, the older one goes)."""
    if path.exists():
        os.replace(path, path.with_name(path.name + ".prev"))


def run_step(step: Step, number: int, env: dict[str, str], logs: Path, timeout: float,
             grace: float = GRACE) -> Result:
    """Run one step to its end (or to its time limit) and return how it ended. Never raises for a failing script.

    The child's output goes straight into its log file (the child holds the file handle, so nothing passes through
    this process and the log grows as the child prints). The runner only waits, twice a second, notes which
    processes are below the step, and checks the clock. When the step has ended it gives its FEMMs `grace` seconds to
    quit, counts the ones that did not, and ends the ones it saw below the step.
    """
    cwd = P / step.folder
    for folder in step.mkdirs:
        (P / folder).mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    log_path = logs / f"{number:02d}-{step.name}.log"
    keep_previous_log(log_path)
    command = [sys.executable, step.script, *step.args]
    before = femm_pids(running_processes())                 # the FEMMs that exist now are none of this step's business
    status = None
    own: dict[int, Identity] = {}                           # the step's python.exe and what was seen below it
    with open(log_path, "wb") as log:
        log.write((f"# run_all step {number}: {step.name}\n# command: {command_text(step)}\n# folder:  {cwd}\n"
                   f"# started: {now_text()}   time limit {fmt_duration(timeout)}\n{'#' * 78}\n").encode("utf-8"))
        log.flush()
        started = time.monotonic()
        proc = subprocess.Popen(command, cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=log,
                                stderr=subprocess.STDOUT)
        root = process_identity(proc.pid)
        if root is not None:
            own[proc.pid] = root
        try:
            while True:
                try:
                    code = proc.wait(timeout=0.5)
                    break
                except subprocess.TimeoutExpired:
                    try:
                        look_below(own, running_processes(), proc.pid)
                    except OSError:
                        pass                                    # one missed look is not worth a failed step
                    if time.monotonic() - started > timeout:
                        status = "timeout"
                        kill_step(proc, own)
                        code = proc.wait()
                        break
        except KeyboardInterrupt:
            status = "interrupted"
            kill_step(proc, own)
            code = proc.wait()
        except BaseException:
            kill_step(proc, own)                        # a fault in this runner must not leave a study running unattended
            raise
        seconds = time.monotonic() - started
    # The step has ended. Give its FEMMs a moment to quit (FEMM 5.0 quits when its python hangs up), then count the
    # ones that did not. A Ctrl-C in here still gets the step its row.
    ours: dict[int, Identity] = {}
    others: set[int] = set()
    table: dict[int, tuple[int, str]] = {}
    try:
        deadline = time.monotonic() + grace
        while True:
            table = running_processes()
            forget_gone(own, table)
            ours, others = leaked_femms(before, table, own)
            if not (ours or others) or time.monotonic() >= deadline:
                break
            time.sleep(0.25)
    except KeyboardInterrupt:
        status = status or "interrupted"
    lines = []
    for pid, identity in sorted(ours.items()):
        ended = terminate_if_same(pid, identity)        # one leak must not disturb the steps after it
        lines.append(f"\n#   still running: {identity[0]} (pid {pid}), "
                     + ("ended by the runner" if ended else "could not be ended"))
    for pid in sorted(others):
        lines.append(f"\n#   still running: {table.get(pid, (0, '?'))[1]} (pid {pid}), not seen below the step: left alone")
    left = len(ours) + len(others)
    if status is None:
        status = "ok" if code == 0 and not left else "failed" if code != 0 else "leak"
    with open(log_path, "ab") as log:
        log.write((f"{'#' * 78}\n# finished: {now_text()}   exit code {code}   {seconds:.1f} s   status {status}   "
                   f"FEMM programs left running: {left}" + "".join(lines) + "\n").encode("utf-8"))
    return Result(number, step, code, seconds, left, status, log_path, now_text())


# ================================================================================================ the table of results
def ledger_row(result: Result) -> dict[str, object]:
    return {"step": result.number, "name": result.step.name, "command": command_text(result.step),
            "exit_code": result.exit_code, "seconds": f"{result.seconds:.1f}", "femm_left_running": result.femm_left,
            "cwd": result.step.folder, "status": result.status, "finished": result.finished}


def pending_row(step: Step, number: int) -> dict[str, object]:
    """The row of a step that this run will run but has not run yet."""
    return {"step": number, "name": step.name, "command": command_text(step), "exit_code": "", "seconds": "",
            "femm_left_running": "", "cwd": step.folder, "status": "not-run", "finished": ""}


def update_ledger(path: Path, new_rows: list[dict], steps: tuple[Step, ...] = STEPS) -> list[dict]:
    """Put these rows into run_all.csv: replace the rows of the same steps, keep the others, drop names that are no
    longer steps, renumber, and write the whole file in table order. Returns the rows written."""
    kept: dict[str, dict] = {}
    try:
        with open(path, newline="", encoding="utf-8") as f:
            kept = {row["name"]: row for row in csv.DictReader(f) if row.get("name")}
    except (OSError, ValueError, csv.Error, KeyError):  # ValueError includes UnicodeDecodeError
        pass                                            # no table yet, or one we cannot read: start a new one
    kept.update({row["name"]: row for row in new_rows})
    rows = []
    for number, step in enumerate(steps, start=1):
        if step.name in kept:
            rows.append({**kept[step.name], "step": number})
    for attempt in range(5):                            # an editor that has the file open can refuse a replace
        try:
            path.parent.mkdir(parents=True, exist_ok=True)
            temporary = path.with_name(path.name + ".tmp")
            with open(temporary, "w", newline="", encoding="utf-8") as f:
                writer = csv.DictWriter(f, fieldnames=COLUMNS, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)
            os.replace(temporary, path)
            break
        except OSError as exc:
            if attempt == 4:
                print(f"warning: could not write {path}: {exc}", flush=True)
            time.sleep(0.5)
    return rows


# ================================================================================================ one runner at a time
def take_lock(path: Path) -> int | None:
    """Create the lock file holding this runner's process id. Returns None if it is ours now, or the process id of a
    runner that still holds it. A lock left by a runner that is gone is replaced."""
    path.parent.mkdir(parents=True, exist_ok=True)
    for _ in range(3):
        try:
            fd = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            try:
                holder = int(path.read_text().strip())
            except (OSError, ValueError):
                holder = 0
            entry = running_processes().get(holder)
            if entry is not None and entry[1].startswith("python") and holder != os.getpid():
                return holder
            path.unlink(missing_ok=True)                # stale: its runner is gone
            continue
        with os.fdopen(fd, "w") as f:
            f.write(str(os.getpid()))
        return None
    raise OSError(f"could not take the lock {path}")


def drop_lock(path: Path) -> None:
    try:
        if path.read_text().strip() == str(os.getpid()):
            path.unlink()
    except OSError:
        pass


# ================================================================================================ the command line
def print_list(steps: list[Step], factor: float) -> None:
    """The steps in table order: number, name, command, folder, both times, time limit, and what each needs."""
    number = {s.name: i for i, s in enumerate(STEPS, start=1)}
    name_w = max(len(s.name) for s in steps)
    cmd_w = max(len(command_text(s)) for s in steps)
    print(f"{'#':>2}  {'name':<{name_w}}  {'command':<{cmd_w}}  {'folder':<6}  {'FEMM 4.2':<13}  {'FEMM 5.0':<10}  time limit")
    for s in steps:
        print(f"{number[s.name]:>2}  {s.name:<{name_w}}  {command_text(s):<{cmd_w}}  {s.folder:<6}  {s.femm42:<13}  "
              f"{fmt_duration(s.measured):<10}  {fmt_duration(timeout_for(s, factor))}")
        for provider, what in s.needs:
            print(f"{'':>2}  {'':<{name_w}}  needs {what}  <- {provider}")


def child_environment() -> dict[str, str]:
    """The environment each step gets. Imports femm_session here, not at the top, so that --list needs no FEMM."""
    if str(P) not in sys.path:
        sys.path.insert(0, str(P))
    import femm_session
    env = dict(os.environ)
    env["FEMM_DIR"] = str(femm_session.femm_dir())
    env["MOTOR_FIGURES_DIR"] = str(FIGURES)
    env["PYTHONUNBUFFERED"] = "1"
    env["PYTHONIOENCODING"] = "utf-8"
    return env


def log_tail(path: Path, lines: int = 8) -> str:
    try:
        with open(path, "rb") as f:
            f.seek(0, os.SEEK_END)
            f.seek(max(0, f.tell() - 4000))
            text = f.read().decode("utf-8", errors="replace")
    except OSError:
        return ""
    kept = [ln for ln in text.splitlines() if ln.strip() and not ln.startswith("#")]
    return "\n".join("        | " + ln for ln in kept[-lines:])


def print_summary(results: list[Result], rows: list[dict], planned: int | None = None) -> None:
    print("\n" + "=" * 100)
    print(f"{'#':>2}  {'name':<26} {'status':<12} {'exit':>5} {'seconds':>9} {'FEMM left':>9}")
    for r in results:
        print(f"{r.number:>2}  {r.step.name:<26} {r.status:<12} {r.exit_code:>5} {r.seconds:>9.1f} {r.femm_left:>9}")
    bad = [r for r in results if r.status != "ok"]
    failed_names = {r.step.name for r in bad}
    for r in bad:
        print(f"\nNOT OK: step {r.number} {r.step.name}: {r.status}, exit code {r.exit_code}, "
              f"{r.femm_left} FEMM program(s) left running")
        print(f"        log: {r.log}")
        blame = [provider for provider, _ in r.step.needs if provider in failed_names]
        if blame:
            print(f"        its input comes from {', '.join(blame)}, which also failed: fix that first")
        tail = log_tail(r.log)
        if tail:
            print(tail)
    skipped = (planned - len(results)) if planned is not None else 0
    print(f"\n{len(results)} step(s) run, {len(results) - len(bad)} ok, {len(bad)} not ok"
          + (f", {skipped} not run (the run was stopped)." if skipped else "."))
    ok_rows = sum(1 for row in rows if row["status"] == "ok")
    not_run = sum(1 for row in rows if row["status"] == "not-run")
    print(f"{LEDGER}: {len(rows)} of {len(STEPS)} steps recorded: {ok_rows} ok, "
          f"{len(rows) - ok_rows - not_run} not ok, {not_run} not run.")


def main(argv: list[str] | None = None) -> int:
    for stream in (sys.stdout, sys.stderr):                 # a log line with a character the console's code page lacks
        try:                                                # (a redirected stdout is often cp1252) must not end the run
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
    ap = argparse.ArgumentParser(prog="run_all.py", description="Re-run every script of this folder on FEMM 5.0, one at "
                                 "a time. See the module's docstring for the order, the table and what each step proves.")
    ap.add_argument("--only", action="append", metavar="NAME[,NAME...]", help="run just these steps (in table order)")
    ap.add_argument("--from", dest="from_step", metavar="NAME", help="run this step and every step after it")
    ap.add_argument("--list", action="store_true", help="print the steps and what each needs; run nothing")
    ap.add_argument("--timeout-factor", type=float, default=TIMEOUT_FACTOR, metavar="N",
                    help=f"a step may run N x the time it took in the first FEMM 5.0 sweep, at least 5 min "
                         f"(default {TIMEOUT_FACTOR:g})")
    args = ap.parse_args(argv)
    if args.only is not None and args.from_step:
        ap.error("use --only or --from, not both")
    try:
        steps = select_steps(args.only, args.from_step)
    except ValueError as exc:
        ap.error(str(exc))
    if args.list:
        print_list(steps, args.timeout_factor)
        return 0

    holder = take_lock(LOCK)
    if holder is not None:
        print(f"another run_all.py is running (process {holder}); one FEMM-heavy run at a time.\n"
              f"If it is not, delete {LOCK}.")
        return 2
    try:
        try:
            env = child_environment()
        except (FileNotFoundError, ImportError, ValueError) as exc:     # no FEMM_DIR, the FEMM binding cannot be imported, or config.toml is not valid TOML (ValueError)
            print(f"cannot start: {exc}")
            return 2
        import femm_session
        try:
            femm_session.lower_priority()                   # children and their FEMMs inherit below-normal priority
        except OSError as exc:
            print(f"warning: could not lower the priority ({exc}); the run will compete with the desktop", flush=True)
        try:
            suppress_crash_dialogs()                        # children inherit this too
        except OSError as exc:
            print(f"warning: could not switch off the crash window ({exc})", flush=True)
        try:
            join_kill_on_exit_job()
        except OSError as exc:
            print(f"warning: no kill-on-exit job ({exc}); stop the runner with Ctrl-C, not by killing it", flush=True)
        number = {s.name: i for i, s in enumerate(STEPS, start=1)}
        print(f"run_all: {len(steps)} of {len(STEPS)} steps, in {P}\n         python {sys.executable}\n"
              f"         FEMM_DIR {env['FEMM_DIR']}   MOTOR_FIGURES_DIR {env['MOTOR_FIGURES_DIR']}\n"
              f"         logs {LOGS}\n         table {LEDGER}", flush=True)
        results: list[Result] = []
        # every step of this run is "not-run" until it has run, so a run that stops early cannot leave the old "ok" behind
        rows = update_ledger(LEDGER, [pending_row(s, number[s.name]) for s in steps])
        try:
            for step in steps:
                limit = timeout_for(step, args.timeout_factor)
                print(f"\n[{number[step.name]:>2}/{len(STEPS)}] {step.name}: {command_text(step)}"
                      f"{'' if step.folder == '.' else '   (in ' + step.folder + ')'}\n"
                      f"        FEMM 4.2 {step.femm42 or '-'}, FEMM 5.0 {fmt_duration(step.measured)}; "
                      f"time limit {fmt_duration(limit)}; started {now_text()}", flush=True)
                try:
                    result = run_step(step, number[step.name], env, LOGS, limit)
                except Exception:                           # a bug in this runner must not stop the whole run
                    traceback.print_exc()
                    result = Result(number[step.name], step, -1, 0.0, 0, "runner-error",
                                    LOGS / f"{number[step.name]:02d}-{step.name}.log", now_text())
                results.append(result)
                rows = update_ledger(LEDGER, [ledger_row(result)])
                print(f"        -> {result.status}: exit {result.exit_code}, {fmt_duration(result.seconds)}, "
                      f"{result.femm_left} FEMM program(s) left running", flush=True)
                if result.status == "interrupted":
                    print("interrupted: the remaining steps were not run", flush=True)
                    break
        finally:
            print_summary(results, rows, planned=len(steps))
        return 0 if results and all(r.status == "ok" for r in results) else 1
    finally:
        drop_lock(LOCK)


if __name__ == "__main__":
    sys.exit(main())
