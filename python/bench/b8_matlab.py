"""B8: Python -> MATLAB (matlab -batch): a flux map of the band-free motor from FEMM, interpolated in MATLAB, the MTPA angle found there and read back by Python.

What it shows, and why
    MATLAB and Simulink are where the controller and the inverter will be modelled, and the flux map of the motor, psi_d and psi_q over the d and q currents, is what hands the
    motor to them (a controller looks the flux up, it does not run FEMM). Task 4 proved that `matlab -batch` starts, runs a script, gets an answer back and lets go. B8 sends the
    MOTOR through it, the way the pipeline will: FEMM computes the map, Python saves it as a .mat file, a MATLAB script loads it, interpolates it and walks the circle of 10 A in the
    current angle to find the angle of maximum torque per ampere (MTPA), and Python reads the answer back and holds it against its own (scipy, the same map).
    The case is sound when the two angles agree within 0.5 degrees (the brief), when no project file and no MATLAB preference changes (the brief's "no file outside the scratch
    folder changes", narrowed by a design decision: MathWorks' own files change whoever starts MATLAB, and are listed), and when exactly one MATLAB ran, under the solver slot,
    with every process of its job below normal.

The map (FEMM, through femm_pool, coarse mesh)
    GRID_POINTS x GRID_POINTS = 6 x 6 points over id in [-10, 0] A and iq in [0, 10] A, so that the 10 A circle (id = -10 sin(gamma), iq = 10 cos(gamma), gamma 0 to 90 degrees)
    lies inside it. Each point is solved at the rotor positions theta_e = 0 and 30 electrical degrees (the pair that cancels the 6th, 18th and 30th harmonics of the ripple, B7's
    choice) and psi_d and psi_q are the means: 36 x 2 = 72 FEMM solves (bench/femm_points.py, the worker B7 shares), in 8 tasks. FEMM is driven through the dq currents
    (motorgeo.FemmModel.solve: the inverse Park transform to the three phases), and the flux linkages come back through the Park transform of the phase flux linkages.

How Python passes the map to MATLAB and reads the answer (the files and the protocol)
    * flux_map.mat (scipy.io.savemat): id_A and iq_A (the axes), psi_d_Wb and psi_q_Wb ((id x iq): one row per id, one column per iq, the layout interp2(iq, id, V, ...) wants),
      a CONTROL map ctrl_psi_d_Wb and ctrl_psi_q_Wb (psi_d = psi_m + Ld id, psi_q = Lq iq, whose MTPA angle has a closed form, so that a map read the wrong way shows at once),
      pole_pairs, current_A and gamma_deg (the 9,001 angles from 0 to 90 degrees in steps of 0.01 that both sides search).
    * b8_mtpa_map.m (bench/, tracked; copied into the scratch folder, so MATLAB runs the copy), started by `matlab -sd <folder> -logfile matlab.log -singleCompThread
      -noFigureWindows -batch b8_mtpa_map` through bench.run_program (the solver slot, a timeout, below-normal priority). `-batch` is MathWorks' non-interactive mode: "Exits
      automatically with exit code 0 if statement executes successfully. Otherwise, MATLAB terminates with a non-zero exit code"; `-sd` sets the start folder; `-logfile` copies
      the command window; `-singleCompThread` limits MATLAB to one computational thread; `-noFigureWindows` keeps -batch, which otherwise "displays figure windows", from opening
      one (https://www.mathworks.com/help/matlab/ref/matlabwindows.html and .../commonly-used-startup-options.html, read 2026-11-06).
    * answer.mat (gamma_deg, torque_Nm, mtpa_deg, mtpa_Nm and the control's) and answer.json (the scalars and the facts about this MATLAB: version, prefdir, threads), which MATLAB
      writes in the same folder and Python reads (scipy.io.loadmat, json). figure.png is MATLAB's own: exportgraphics of a figure it draws with no window, in MATLAB's LIGHT theme
      (theme(fig, "light"): "sets the theme of the specified figure", introduced in R2025a, https://www.mathworks.com/help/matlab/ref/theme.html, read 2026-11-06; that figure only, no
      setting is changed). The picture goes to a white page, and exportgraphics follows the theme MATLAB runs with: a probe of 2026-11-06 with the graphics theme set to dark for a
      session (in scratch preferences) gave a black picture with white text without the call and the white one with it; with MATLAB's usual theme the two are the same.
    * The interpolation is LINEAR on both sides (MATLAB interp2 'linear' = scipy RegularGridInterpolator 'linear'), and the search is the same: the angle of the largest torque
      among the 9,001 angles. The case proves the path Python -> .mat -> MATLAB -> answer -> Python, not the interpolation; so it also reports, not gating, what a cubic spline
      gives in Python, and how the linear map's angle and torque curve compare with B1's 22.9 degrees and its curve (a coarse mesh and a 6 x 6 map will not match them exactly).

"No project file and no MATLAB preference changes" (the brief's "no file outside the scratch folder changes"): what was looked at, and what MATLAB does write
    Before and after the run, SHA-256 of: the inputs in the scratch folder (flux_map.mat and the script's copy), every file of reference\\ (the drawings, the plant: the CAD and FEM
    files are read-only), and the repository's copy of the script; and of every file of MATLAB's PREFERENCES FOLDER (prefdir), which is outside the scratch folder too. The real
    prefdir is <APPDATA>\\MathWorks\\MATLAB\\<release>: the form of the path that MathWorks' prefdir page shows as the example of what prefdir returns (C:\\Users\\<name>\\AppData\\Roaming\\
    MathWorks\\MATLAB\\R2025a: an example output, not a statement of a default), and what a `matlab -batch "disp(prefdir)"` with no redirect printed on this PC on 11-06
    (C:\\Users\\<name>\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b). That trial, run twice, changed exactly ONE of the folder's 7,640 files (197.4 MB): epfwk_cache-<build>-<hash>.json,
    rewritten each time with the same size (and, in the run that was compared byte for byte, other content), so MATLAB does write there. So B8 runs MATLAB with the environment variable
    MATLAB_PREFDIR on a folder inside the scratch folder ("To override the default location, create an environment variable named MATLAB_PREFDIR. MATLAB reads this value at startup":
    MathWorks Support Team, https://www.mathworks.com/matlabcentral/answers/93696, updated 2021-07-02; the `prefdir` reference page does not mention it). MATLAB then reports that folder
    as its prefdir (answer.json), creates 8 files (0.5 MB) in it, and the real prefdir, listed and hashed before and after (key and certificate files by size and time only: their
    contents are never read), is unchanged. A file that cannot be read (locked by the user's own MATLAB, say) is named, and the check FAILS: it cannot be shown unchanged.
    What MATLAB writes that it does not let us move is MathWorks' own: licensing and usage files under AppData\\Roaming\\MathWorks, the ServiceHost's logs and its credentials service's
    config under AppData\\Local\\MathWorks. B8 lists every such path with its size in the check "MathWorks' own files" (and in outside_writes.txt), reported and not gating: some of it
    may be the logon service host's own.

Exactly one MATLAB, under the slot (what was asked to be seen)
    checks.matlab_check's JobTracker follows every process that descends from the program run_solver started, from the process list read twenty times a second; the case counts
    the matlab.exe launcher (a child of this process), the MATLAB.exe engine below it, and any other MATLAB-named process in the job (none may be). The slot's holder is looked at
    every 0.2 s while MATLAB runs (femm_queue.holder()): every look between the first and the last sighting of matlab.exe must name this run's label and pid. After run_solver
    returns, nothing of the job is alive, and how soon after the last matlab.exe the job was empty is compared with Task 4's measurement (task-4-report.md section 5: within 0.05 s, the
    MathWorksServiceHost clients of the job aside). MATLAB processes this run did not start (the service host started at logon, its monitor and its crash reporter: they run all
    day) are listed before and after and never touched; only this job's processes count.

Run (from the Python to FEMM folder; MATLAB starts once, the pool once; about 1 minute)
    uv run python bench/run_bench.py --only B8
Writes out\\bench\\B8\\: results.json, figure.png (MATLAB's own), flux_map.mat, b8_mtpa_map.m, answer.mat, answer.json, matlab.log, matlab.out.txt, job_members.txt,
outside_writes.txt, prefdir\\ (the scratch preferences folder), fem\\ (the FEMM models), run.log, runs.json.

Numbers to expect (this PC, 2026-11-06, FEMM 5.0, MATLAB R2025b Update 5, warm; four real runs gave these numbers, the last after the solver slot's job began to hold its processes at
their class; the case takes 35 to 42 s warm, 110 s in the whole bench run, right after PrePoMax and CalculiX: the 7,640 files of the real prefdir were hashed in 57 s instead of 4 s and
MATLAB ran 44 s instead of 29)
    the map     6 x 6 points, 72 solves in 8 tasks (3.3 to 3.4 s), psi_d(0, 0) = 16.558 mWb (B1: 16.629 on the fine mesh, -0.43 %), Ld 1.17 mH and Lq 2.11 mH across the range (Lq above Ld), psi_q
                at iq = 0 is 0.21 % of its largest value, the torque at gamma = 0 is 1.00 N.m
    the angle   MATLAB 22.69 degrees = Python 22.69 (linear, the same 9,001 angles: difference 0.000), the largest torque 1.1223 N.m; the torque curves agree to 6e-16 of the peak;
                a cubic spline in Python gives 22.78; B1's MTPA angle on the fine mesh is 22.94 and its torque curve is matched to 0.65 % of the peak
    the control psi_m 16.6 mWb, Ld 0.1 mH, Lq 1.0 mH: the closed form gives 22.523 degrees, MATLAB 22.52, Python 22.52 (a transposed read would give 15.2)
    MATLAB      exit 0 in 22.8 to 29.5 s warm (44 s in the whole bench run): the script itself takes 5 s, and MATLAB then runs on for 18 to 24 s, until the MATLABWindow.exe processes of the figure end (without a figure the whole run
                is 4 s; that the figure's helpers are what it waits for is not isolated); one computational thread; the scratch prefdir gets 8 files (0.5 MB); the real prefdir (7,640 files, 197.4 MB) is unchanged
    the job     matlab.exe (the launcher, a child of this run) and MATLAB.exe below it; 3 MathWorksServiceHost.exe clients that MATLAB starts (one for the whole run, two for 0.1 to 0.2 s); and,
                because of the figure, 8 MATLABWindow.exe (MATLAB's web-graphics window host: one with seven below it) from about 3.5 s to the end. Empty 0.11 to 0.16 s after the last matlab.exe
                (Task 4, without a figure: 0.00 to 0.05 s; 0.19 s in the whole bench run). The slot's holder was this run at every one of 110 to 215 looks. Already running and never touched:
                MathWorksServiceHost.exe,
                MathWorksCrashReporter.exe and MathWorksServiceHost-Monitor.exe, started at logon
    priority    every process of the job, MATLABWindow.exe included: below normal at every one of 770 looks (matlab.exe, MATLAB.exe, the service hosts and the 8 MATLABWindow.exe), the rule of
                the gating check "every process of MATLAB's job was seen below normal". The solver slot's job holds them there (femm_queue._Job, JOB_OBJECT_LIMIT_PRIORITY_CLASS): before the
                job did, in the first real runs of 11-06, about half of the 8 MATLABWindow.exe were seen at normal and above normal (MATLAB's own helper: whether it set its class itself or
                was created with another one was not isolated, and the job that only killed its members on close did not hold it)
    elsewhere   MATLAB wrote 6 or 7 files (4 to 9 MB in all) outside the scratch folder: the service host's logs (3) and its credentials service's config under AppData\\Local\\MathWorks, the
                license file and the usage log under AppData\\Roaming\\MathWorks\\licensing; every path is listed in the check "MathWorks' own files", and no file of the project changed

Self-check (the answers follow).
  1. MATLAB's angle and Python's agree to 0.01 degrees. Does that show that the map is right?
  2. Why is there a control map, and why does it have a closed form?
  3. MATLAB is started with MATLAB_PREFDIR set. Why is the variable set in this process and removed afterwards, and not given to MATLAB some other way?
Answers.
  1. No: it shows that the path works. Both sides read the same file, with the same method; a wrong map would give the same wrong angle on both. The physics of the map is held
     against B1 (its psi_PM, its torque curve) and against what an IPM must do (Lq above Ld, no q flux without q current).
  2. If MATLAB took the file's rows for columns, its angle would be wrong and so would a scipy that made the same mistake, and they would still agree. A map whose answer is known in
     closed form (psi_d = psi_m + Ld id, psi_q = Lq iq: T = 1.5 p I (psi_m cos g + (Lq - Ld) I sin g cos g), whose maximum is at sin g = (-psi_m + sqrt(psi_m^2 + 8 D^2)) / (4 D),
     D = (Lq - Ld) I) tests the layout against something that is neither of them.
  3. run_solver starts the program with a copy of this process's environment and has no argument for another. The variable is set around that one call, in this process, and put back
     whatever happens, so nothing else this process starts afterwards inherits it.
"""
from __future__ import annotations

import contextlib
import json
import math
import os
import re
import shutil
import threading
import time
import zlib
from pathlib import Path

import numpy as np

import femm_queue
import motorgeo as mg
import paths
from checks import need_tool, sha256
from checks.matlab_check import JobTracker, Watcher, process_created, process_table

from . import OWN_LIMIT, HolderSampler, Log, check, compare, params, run_program, skipped
from . import femm_points
from .b1_femm import pool_record
from .femm_points import POSITIONS_EL, TASK_TIMEOUT_S

CASE = "B8"
TITLE = "Python -> MATLAB (matlab -batch): a FEMM flux map, interpolated and searched for the MTPA angle in MATLAB, read back"

# ---- the brief's numbers and the rulings
ANGLE_TOLERANCE_DEG = 0.5                                  # the brief: MATLAB's MTPA angle equals Python's own (scipy, the same map) within 0.5 degrees
GRID_POINTS = 6                                            # per axis: at least 5 x 5 (the brief)
ID_RANGE_A, IQ_RANGE_A = (-10.0, 0.0), (0.0, 10.0)         # the ruling: so that the 10 A circle lies inside the map
CURRENT_A = 10.0                                           # the brief: the MTPA angle at 10 A
POLE_PAIRS = 4
GAMMA_STEP_DEG = 0.01                                      # both sides search the 9,001 angles 0, 0.01, ..., 90
INTERPOLATION = "linear"                                   # MATLAB interp2 'linear' = scipy RegularGridInterpolator 'linear' (the ruling)
SPLINE = "cubic"                                           # what Python gets with a spline: reported, not gating
CONTROL_PSI_M, CONTROL_LD, CONTROL_LQ = 0.0166, 0.1e-3, 1.0e-3      # the control map: Wb, H, H. A salient machine of about this motor's size; its closed-form MTPA angle, 22.52 degrees,
                                                                   # is 7.3 degrees away from what a map read with rows for columns gives, so the control can tell
SCRIPT_NAME = "b8_mtpa_map"
MAT_VARIABLES = ("id_A", "iq_A", "psi_d_Wb", "psi_q_Wb", "ctrl_psi_d_Wb", "ctrl_psi_q_Wb", "pole_pairs", "current_A", "gamma_deg")
LABEL = "bench B8 MATLAB"                                  # what the solver slot says about the holder while MATLAB runs
FEMM_LABEL = "bench B8 FEMM"
TIMEOUT_S = 300                                            # MATLAB's own run is about 10 to 30 s; the time limit of Task 4's row
THREADS = 2                                                # OMP_NUM_THREADS and the like in MATLAB's environment; MATLAB itself is held to one thread (-singleCompThread)
THREAD_LIMIT = 8                                           # config.toml [femm] workers, the machine's limit for threads
CHUNK = 9                                                  # points per FEMM task (one model per task): 72 points in 8 tasks
POLL_S = 0.2                                               # how often the slot's holder is looked at while MATLAB runs
# ---- the bench's own limits (the brief gives none)
CURVE_TOLERANCE = 1e-9                                     # MATLAB's torque curve against scipy's, relative to the peak: the same maths on the same numbers
CONTROL_TOLERANCE_DEG = 0.05                               # either side's angle on the control map against its closed form: 5 steps of the search
B1_PSI_TOLERANCE_PERCENT = 3.0                             # psi_d at no current (coarse mesh, two positions) against B1's psi_PM (fine mesh, 120 positions)
B1_ANGLE_TOLERANCE_DEG = 2.0                               # the map's MTPA angle against B1's: the torque is flat near its peak, so the angle is not sharp
B1_CURVE_TOLERANCE_PERCENT = 3.0                           # the map's torque curve against B1's, at B1's angles, in per cent of the peak
CURRENT_ERROR_TOLERANCE_A = 1e-9                           # the currents FEMM was driven at against the ones asked for
Q_AXIS_TOLERANCE_PERCENT = 5.0                             # psi_q with no q current, in per cent of its largest value: zero by the symmetry of the machine about its d axis
EMPTY_TOLERANCE_S = 0.25                                   # the job empty within this long of the last matlab.exe: five looks of the watcher (Task 4: 0.00 to 0.05 s)
KEY_SUFFIXES = (".p12", ".pem", ".pfx", ".key")            # key and certificate files: listed by size and time, never opened
UNREADABLE = "unreadable"                                  # what the hash of a file that could not be read says (then: "unreadable: PermissionError")
# What the MATLAB script (its code, comments cut off) may not contain. Substrings that name a change outside the scratch folder or another program, and whole words that would
# change MATLAB's path, preferences or settings for good ("no permanent change"): \b so that fullpath( or a name that merely ends in 'settings' is not one. tests/test_bench_b8.py
# reads the script against both, and the words against examples.
SCRIPT_FORBIDDEN_TEXT = ("delete(", "rmdir", "movefile", "copyfile", "system(", "dos(", "unix(", "!", "cd(", "addpath", "mkdir", "websave", "urlwrite", "diary", "rethrow(", "exit")
SCRIPT_FORBIDDEN_WORDS = (r"\bsavepath\b", r"\bsetpref\b", r"\brmpref\b", r"\bsettings\b", r"\bpath\(")


# ============================================================================================== the map and Python's own search (no solver, no MATLAB)
def grid_axes(points: int = GRID_POINTS) -> tuple[np.ndarray, np.ndarray]:
    """The axes of the map: id from -10 to 0 A and iq from 0 to 10 A, `points` each, increasing."""
    return np.linspace(ID_RANGE_A[0], ID_RANGE_A[1], points), np.linspace(IQ_RANGE_A[0], IQ_RANGE_A[1], points)


def gamma_grid() -> np.ndarray:
    """The angles both sides search: 0 to 90 degrees in steps of GAMMA_STEP_DEG (9,001 of them; the ends exact)."""
    return np.linspace(0.0, 90.0, int(round(90.0 / GAMMA_STEP_DEG)) + 1)


def circle(gamma_deg, current: float = CURRENT_A) -> tuple[np.ndarray, np.ndarray]:
    """The current (id, iq) at each angle gamma of the circle of `current` amperes: id = -I sin(gamma), iq = I cos(gamma)."""
    g = np.radians(np.asarray(gamma_deg, float))
    return -current * np.sin(g), current * np.cos(g)


def torque_curve(id_axis, iq_axis, psi_d, psi_q, gamma_deg, current: float = CURRENT_A, pole_pairs: float = POLE_PAIRS, method: str = INTERPOLATION) -> np.ndarray:
    """The torque T = 1.5 p (psi_d iq - psi_q id) [N.m] on the circle of `current` amperes, from a flux map interpolated with scipy's RegularGridInterpolator (`method`: 'linear' is
    MATLAB's interp2 'linear'; 'cubic' is the spline). psi_d and psi_q are (id x iq): rows follow id_axis, columns iq_axis."""
    from scipy.interpolate import RegularGridInterpolator
    id_c, iq_c = circle(gamma_deg, current)
    points = np.column_stack([id_c, iq_c])
    flux_d = RegularGridInterpolator((np.asarray(id_axis, float), np.asarray(iq_axis, float)), np.asarray(psi_d, float), method=method)
    flux_q = RegularGridInterpolator((np.asarray(id_axis, float), np.asarray(iq_axis, float)), np.asarray(psi_q, float), method=method)
    return 1.5 * pole_pairs * (flux_d(points) * iq_c - flux_q(points) * id_c)


def mtpa(gamma_deg, torque) -> tuple[float, float]:
    """(the angle, the torque) of the largest torque: the first one if two tie, as MATLAB's max does."""
    k = int(np.argmax(np.asarray(torque, float)))
    return float(np.asarray(gamma_deg)[k]), float(np.asarray(torque)[k])


def analytic_map(id_axis, iq_axis, psi_m: float, ld: float, lq: float) -> tuple[np.ndarray, np.ndarray]:
    """The control map: psi_d = psi_m + Ld id and psi_q = Lq iq, as (id x iq) arrays."""
    id_axis, iq_axis = np.asarray(id_axis, float), np.asarray(iq_axis, float)
    return psi_m + ld * id_axis[:, None] + 0.0 * iq_axis[None, :], 0.0 * id_axis[:, None] + lq * iq_axis[None, :]


def analytic_mtpa_deg(psi_m: float, ld: float, lq: float, current: float = CURRENT_A) -> float:
    """The MTPA angle of the control map in closed form. T = 1.5 p I (psi_m cos g + (Lq - Ld) I sin g cos g); dT/dg = 0 gives 2 D s^2 + psi_m s - D = 0 with s = sin g and
    D = (Lq - Ld) I, so s = (-psi_m + sqrt(psi_m^2 + 8 D^2)) / (4 D)."""
    d = (lq - ld) * current
    return math.degrees(math.asin((-psi_m + math.sqrt(psi_m ** 2 + 8 * d ** 2)) / (4 * d)))


def agree(matlab_deg: float, python_deg: float, tolerance: float = ANGLE_TOLERANCE_DEG) -> bool:
    """The brief's check: |MATLAB - Python| <= 0.5 degrees (on the tolerance passes: 'within'). NaN and infinity never agree."""
    return bool(math.isfinite(matlab_deg) and math.isfinite(python_deg) and abs(matlab_deg - python_deg) <= tolerance + 1e-12)


def forbidden_in_script(code: str) -> list[str]:
    """What the code of the MATLAB script (comments already cut off) contains that it may not: the substrings of SCRIPT_FORBIDDEN_TEXT, and the patterns of SCRIPT_FORBIDDEN_WORDS
    that match as whole words. [] for a clean script."""
    return [text for text in SCRIPT_FORBIDDEN_TEXT if text in code] + [word for word in SCRIPT_FORBIDDEN_WORDS if re.search(word, code)]


_TRANSPOSE_AFTER = frozenset("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_)]}.'")      # a ' right after one of these is MATLAB's transpose operator


def strip_matlab_comments(text: str) -> str:
    """The code of a MATLAB script with its comments cut off, string-aware, for the static scan (forbidden_in_script) to read.

    A % starts a comment to the end of the line, unless it is inside a '...' or "..." literal (so `fprintf('%s\\n', x); savepath` keeps its savepath). A ' right after a letter,
    a digit, _, ), ], }, . or another ' is the transpose operator, not the start of a string (so `x = y'; % savepath` ends in a comment). A doubled quote inside a literal is an
    escaped quote. What follows `...` (a line continuation) is ignored by MATLAB, and is cut. A line that is only %{ opens a block comment and one that is only %} closes it (they
    nest); every line in between is cut. Strings stay: a forbidden word inside one is still seen, which is the careful side. The line breaks stay, so the lines keep their numbers."""
    lines: list[str] = []
    depth = 0                                                  # block comments opened and not yet closed
    for line in text.splitlines():
        marker = line.strip()
        if depth:
            depth += 1 if marker == "%{" else -1 if marker == "%}" else 0
            lines.append("")
            continue
        if marker == "%{":
            depth = 1
            lines.append("")
            continue
        kept: list[str] = []
        quote = ""                                             # the quote of the literal we are inside, or ""
        i = 0
        while i < len(line):
            c = line[i]
            if quote:
                kept.append(c)
                if c == quote:
                    if line[i + 1:i + 2] == quote:             # a doubled quote: an escaped one, still inside
                        kept.append(quote)
                        i += 1
                    else:
                        quote = ""
            elif c == "%" or line.startswith("...", i):
                break
            elif c == '"' or (c == "'" and not (i > 0 and line[i - 1] in _TRANSPOSE_AFTER)):
                quote = c
                kept.append(c)
            else:
                kept.append(c)
            i += 1
        lines.append("".join(kept))
    return "\n".join(lines)


def write_flux_map(path, id_axis, iq_axis, psi_d, psi_q, ctrl_psi_d, ctrl_psi_q, pole_pairs, current, gamma_deg) -> Path:
    """flux_map.mat for MATLAB (scipy.io.savemat): MAT_VARIABLES. The maps are (id x iq): one row per id, one column per iq."""
    import scipy.io
    scipy.io.savemat(str(path), {"id_A": np.asarray(id_axis, float), "iq_A": np.asarray(iq_axis, float), "psi_d_Wb": np.asarray(psi_d, float), "psi_q_Wb": np.asarray(psi_q, float),
                                 "ctrl_psi_d_Wb": np.asarray(ctrl_psi_d, float), "ctrl_psi_q_Wb": np.asarray(ctrl_psi_q, float), "pole_pairs": float(pole_pairs),
                                 "current_A": float(current), "gamma_deg": np.asarray(gamma_deg, float)})
    return Path(path)


def read_answer(folder) -> dict:
    """What MATLAB wrote in `folder`: answer.mat (the curve and the angles) and answer.json (the scalars and the facts). {"error": None, ...} or {"error": why} if a file is missing
    or lacks something: a failed MATLAB step is a failed check that says which, never a KeyError."""
    import scipy.io
    from scipy.io.matlab import MatReadError
    folder = Path(folder)
    if not (folder / "answer.mat").is_file():
        return {"error": "answer.mat is not there"}
    try:
        mat = scipy.io.loadmat(str(folder / "answer.mat"), squeeze_me=True)
    except (OSError, ValueError, TypeError, EOFError, NotImplementedError, MatReadError, zlib.error) as error:     # empty, cut short, damaged, locked, a v7.3 file: whatever it is, it is named
        return {"error": f"answer.mat cannot be read: {type(error).__name__}: {error}"}
    wanted = ("gamma_deg", "torque_Nm", "mtpa_deg", "mtpa_Nm", "ctrl_torque_Nm", "ctrl_mtpa_deg", "ctrl_mtpa_Nm")
    if [k for k in wanted if k not in mat]:
        return {"error": f"answer.mat lacks {[k for k in wanted if k not in mat]}"}
    if not (folder / "answer.json").is_file():
        return {"error": "answer.json is not there"}
    try:
        info = json.loads((folder / "answer.json").read_text(encoding="utf-8"))
    except ValueError as error:
        return {"error": f"answer.json is not JSON: {error}"}
    values = {}
    for k in wanted:
        try:
            values[k] = np.asarray(mat[k], float) if k.endswith(("deg", "Nm")) and np.ndim(mat[k]) else float(mat[k])
        except (TypeError, ValueError) as error:
            return {"error": f"answer.mat's {k} is not a number: {error}"}
    return {"error": None, "info": info, **values}


# ============================================================================================== the command, the environment, the processes
def matlab_command(matlab, folder) -> list:
    """The command line (a list, never one string): MATLAB starts in `folder`, copies its command window to matlab.log, is held to one computational thread, opens no figure
    window, and runs the script b8_mtpa_map as a statement; -batch makes it exit by itself, with a nonzero code if the script fails."""
    folder = Path(folder)
    return [matlab, "-sd", folder, "-logfile", folder / "matlab.log", "-singleCompThread", "-noFigureWindows", "-batch", SCRIPT_NAME]


@contextlib.contextmanager
def temporary_environment(**values):
    """Set environment variables of this process for the block and put them back (an old value restored, a new one removed) however the block ends."""
    old = {name: os.environ.get(name) for name in values}
    os.environ.update({name: str(value) for name, value in values.items()})
    try:
        yield
    finally:
        for name, value in old.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


def scratch_prefdir(folder) -> Path:
    """The preferences folder MATLAB is given for this run: prefdir\\ inside the case's folder."""
    return Path(folder) / "prefdir"


def matlab_release(matlab_root) -> str | None:
    """The release of the MATLAB installed in `matlab_root` ('R2025b'), from its VersionInfo.xml; None if the file is not there or says none."""
    try:
        text = (Path(matlab_root) / "VersionInfo.xml").read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    found = re.search(r"<release>(R\d{4}[ab])</release>", text)
    return found.group(1) if found else None


def release_name(release) -> str:
    """MathWorks' name of a release: 'R2025b'. MATLAB's version('-release') returns '2025b' (without the R), the name of the folder in VersionInfo.xml has it; neither is doubled."""
    if not release:
        return "(release not reported)"
    return str(release) if str(release).startswith("R") else f"R{release}"


def real_prefdir(matlab_root, appdata=None) -> Path | None:
    """MATLAB's own preferences folder: <APPDATA>\\MathWorks\\MATLAB\\<release>. That is the form of the path that MathWorks' prefdir page shows as an example of what prefdir
    returns (C:\\Users\\<name>\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025a: an example output, not a statement of a default), and what a `matlab -batch "disp(prefdir)"` with no
    redirect printed on this PC on 2026-11-06 (C:\\Users\\<name>\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b). None if the release or APPDATA is not known: never a guess."""
    release = matlab_release(matlab_root)
    base = appdata if appdata is not None else os.environ.get("APPDATA")
    return Path(base) / "MathWorks" / "MATLAB" / release if release and base else None


SlotSampler = HolderSampler                                # looks at femm_queue.holder() while MATLAB runs (the one copy is in bench/__init__.py)
GRAPHICS_HOST = "matlabwindow.exe"                         # MATLAB's web-graphics window host: it starts when a figure is drawn and exported (the first real run saw 8 of them)


def classify_job(members, my_pid: int) -> dict:
    """What the processes of the job are: {"launchers" (matlab.exe started by this process), "engines" (a MATLAB.exe below a launcher), "hosts" (MathWorks' service hosts, started by
    MATLAB), "graphics" (MATLABWindow.exe: the window host MATLAB starts to draw and export a figure), "other_matlab" (names of any other MATLAB-named process in the job: a second
    MATLAB would be one), "others" (names of whatever else)}. Windows names are not case-sensitive: launcher and engine are told apart by their parents."""
    launchers = [m for m in members if m.name.lower() == "matlab.exe" and m.ppid == my_pid]
    launcher_pids = {m.pid for m in launchers}
    engines = [m for m in members if m.name.lower() == "matlab.exe" and m.ppid in launcher_pids and m not in launchers]
    hosts = [m for m in members if m.name.lower().startswith("mathworks")]
    graphics = [m for m in members if m.name.lower() == GRAPHICS_HOST]
    named = [m for m in members if m not in launchers and m not in engines and m not in hosts and m not in graphics]
    return {"launchers": len(launchers), "engines": len(engines), "hosts": len(hosts), "graphics": len(graphics), "other_matlab": [m.name for m in named if "matlab" in m.name.lower()],
            "others": [m.name for m in named if "matlab" not in m.name.lower()]}


def exactly_one_matlab(job: dict) -> bool:
    """One launcher, one engine below it, and no other MATLAB-named process in the job."""
    return job["launchers"] == 1 and job["engines"] == 1 and not job["other_matlab"]


class ClassSampler(threading.Thread):
    """Looks at the priority class of every process the tracker knows, every `poll_s` seconds, and keeps for each (pid, program) the classes in the order they were seen, with the
    time each was first seen: {(pid, name): [(seconds, class name)]}. It is how a process that raises its own priority shows (inherited below normal, then normal). Read-only.

    The tracker keeps its members after they end, and Windows can give the pid of an ended process to a stranger. So a member is looked at only while the process on its pid has
    the creation time the member carries (checks.matlab_check.process_created); a member that carries none cannot be checked and is not looked at (review M3)."""

    def __init__(self, tracker, zero: float, poll_s: float = 0.1):
        super().__init__(daemon=True, name="class-sampler")
        self.tracker, self.zero, self.poll_s = tracker, zero, poll_s
        self.sequences: dict[tuple[int, str], list[tuple[float, str]]] = {}
        self._stop_now = threading.Event()

    def run(self) -> None:
        from . import class_name, priority_class
        while not self._stop_now.is_set():
            now = time.perf_counter() - self.zero
            for member in self.tracker.members():
                if member.created is None or process_created(member.pid) != member.created:
                    continue                                       # another process now (or one we cannot tell): its class is not this member's
                value = priority_class(member.pid)
                if value is None:
                    continue
                seen = self.sequences.setdefault((member.pid, member.name), [])
                if not seen or seen[-1][1] != class_name(value):
                    seen.append((round(now, 2), class_name(value)))
            self._stop_now.wait(self.poll_s)

    def stop(self) -> None:
        self._stop_now.set()
        self.join(5)

    def left_below_normal(self) -> dict[str, list[str]]:
        """{"<program> (pid N)": ["<class> at <seconds> s", ...]} for every process that was ever seen at a class other than below normal."""
        return {f"{name} (pid {pid})": [f"{c} at {t:g} s" for t, c in seq] for (pid, name), seq in sorted(self.sequences.items()) if any(c != "below_normal" for _, c in seq)}


def holders_in_window(samples, first: float, last: float) -> dict:
    """{"looks": how many looks were taken between `first` and `last` seconds, "holders": the set of what they saw (None for a free slot)}."""
    inside = [seen for t, seen in samples if first <= t <= last]
    return {"looks": len(inside), "holders": set(inside)}


# ============================================================================================== no file outside the scratch folder changes
def tree_snapshot(folder) -> dict[str, tuple[int, int, str | None]]:
    """{path relative to `folder`: (size, modification time in ns, SHA-256)} of every file under it. A key or certificate file (KEY_SUFFIXES) has None for its hash: it is listed
    by size and time and never opened.

    The folder can be the user's own MATLAB's, so a file can be locked, or can vanish between the listing and the look (review M4): a file that vanished is not in the snapshot
    (a diff against a snapshot that had it shows it removed); one that cannot be read has "unreadable: <error name>" for its hash (and size and time -1 if even those could not
    be read). unreadable_files() lists them, and the check that needs the folder unchanged FAILS naming them: they cannot be shown unchanged."""
    folder = Path(folder)
    out: dict[str, tuple[int, int, str | None]] = {}
    for path in sorted(folder.rglob("*")):
        name = str(path.relative_to(folder))
        try:
            if not path.is_file():
                continue
            stat = path.stat()
        except FileNotFoundError:
            continue
        except OSError as error:
            out[name] = (-1, -1, f"{UNREADABLE}: {type(error).__name__}")
            continue
        digest = None
        if path.suffix.lower() not in KEY_SUFFIXES:
            try:
                digest = sha256(path)
            except FileNotFoundError:
                continue
            except OSError as error:
                digest = f"{UNREADABLE}: {type(error).__name__}"
        out[name] = (stat.st_size, stat.st_mtime_ns, digest)
    return out


def unreadable_files(snapshot: dict) -> list[str]:
    """The names in a tree_snapshot whose file could not be read (their hash says so), sorted."""
    return sorted(name for name, (_, _, digest) in snapshot.items() if isinstance(digest, str) and digest.startswith(f"{UNREADABLE}:"))


def tree_diff(before: dict, after: dict) -> dict[str, list[str]]:
    """What differs between two snapshots: {"added", "removed", "changed"} (sorted names). Changed means any of size, time or hash: a file rewritten with the same content is
    still touched, and the real preferences folder must not even be touched."""
    return {"added": sorted(set(after) - set(before)), "removed": sorted(set(before) - set(after)), "changed": sorted(k for k in set(before) & set(after) if before[k] != after[k])}


def changed_since(roots, since_ns: int, exclude=(), skip_names=()) -> list[dict]:
    """The files under `roots` whose modification time is `since_ns` or later: [{"path", "bytes", "mtime_ns"}] sorted by path. Folders in `exclude`, and folders named in
    `skip_names`, are not entered; a root that is not there is skipped. Read-only (a walk and a stat)."""
    skip = {os.path.normcase(os.path.abspath(p)) for p in exclude}
    found: list[dict] = []

    def walk(folder: str) -> None:
        try:
            entries = list(os.scandir(folder))
        except OSError:
            return
        for entry in entries:
            try:
                if entry.is_dir(follow_symlinks=False):
                    if entry.name not in skip_names and os.path.normcase(os.path.abspath(entry.path)) not in skip:
                        walk(entry.path)
                elif entry.is_file(follow_symlinks=False):
                    stat = entry.stat(follow_symlinks=False)
                    if stat.st_mtime_ns >= since_ns:
                        found.append({"path": entry.path, "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns})
            except OSError:
                continue
    for r in roots:
        if Path(r).is_dir():
            walk(str(r))
    return sorted(found, key=lambda f: f["path"])


def processes_like(table: dict, pattern: str = r"matlab|mathworks") -> list[tuple[int, str]]:
    """[(pid, name)] of the processes in a process table whose program name matches the pattern: what else is running that has to do with MATLAB."""
    return sorted((pid, name) for pid, (_, name) in table.items() if re.search(pattern, name, re.I))


# ============================================================================================== the FEMM side
def map_tasks(id_axis, iq_axis, fem_dir) -> tuple[list[dict], list[list[tuple[int, int, int]]]]:
    """The tasks that solve the map: every (id, iq) at every rotor position, in chunks of CHUNK points (one model per task), and for each task which (i, j, k) each of its points is."""
    points: list[tuple[tuple[int, int, int], dict]] = []
    for i, a in enumerate(id_axis):
        for j, b in enumerate(iq_axis):
            for k, theta in enumerate(POSITIONS_EL):
                points.append(((i, j, k), {"theta_e": float(theta), "id": float(a), "iq": float(b)}))
    tasks, index = [], []
    for n, start in enumerate(range(0, len(points), CHUNK)):
        chunk = points[start:start + CHUNK]
        tasks.append({"id": n, "tag": f"map_c{n:02d}", "overrides": {}, "fem_dir": str(fem_dir), "points": [p for _, p in chunk]})
        index.append([ijk for ijk, _ in chunk])
    return tasks, index


def assemble_map(results: list[dict], tasks: list[dict], index, shape: tuple[int, int]) -> dict:
    """psi_d and psi_q as (id x iq) arrays, the mean over the rotor positions, from the pool's results; also the largest difference between the currents FEMM was driven at
    (the ones motorgeo.FemmModel.solve reports back) and the ones the tasks asked for, and how many rotor positions each point was averaged over."""
    sums = {key: np.zeros(shape) for key in ("psi_d", "psi_q", "T_age", "T_dq")}
    count = np.zeros(shape)
    worst = 0.0
    for result, task, ijk in zip(results, tasks, index):
        for row, asked, (i, j, _) in zip(result["rows"], task["points"], ijk):
            for key in sums:
                sums[key][i, j] += row[key]
            count[i, j] += 1
            worst = max(worst, abs(row["id"] - asked["id"]), abs(row["iq"] - asked["iq"]))
    return {"psi_d": sums["psi_d"] / count, "psi_q": sums["psi_q"] / count, "T_age": sums["T_age"] / count, "T_dq": sums["T_dq"] / count, "positions": int(count.min()),
            "current_error": worst}


# ============================================================================================== the case
def run(folder: Path, log: Log) -> dict:
    """Run B8 in `folder` (already empty). Returns the results dict, whose "checks" list every verdict."""
    folder = Path(folder)
    started = time.perf_counter()
    matlab = need_tool("matlab")
    matlab_root = matlab.parent.parent
    reference_dir = paths.path("reference_fem").parent
    script_source = paths.P / "bench" / f"{SCRIPT_NAME}.m"
    # ---- what must not change: listed and hashed before anything runs
    t = time.perf_counter()
    reference_before = tree_snapshot(reference_dir)
    script_before = sha256(script_source)
    real = real_prefdir(matlab_root)
    prefdir_before = tree_snapshot(real) if real is not None and real.is_dir() else None
    log.say(f"B8: before: reference\\ {len(reference_before)} files, the repository's script, MATLAB's real prefdir {real} "
            f"({'not found' if prefdir_before is None else f'{len(prefdir_before):,} files, {sum(v[0] for v in prefdir_before.values()) / 1e6:.1f} MB'}): hashed in {time.perf_counter() - t:.1f} s")
    # ---- the flux map from FEMM, through femm_pool, on the coarse mesh
    id_axis, iq_axis = grid_axes()
    tasks, index = map_tasks(id_axis, iq_axis, folder / "fem")
    with log.tee():
        with femm_queue.femm_pool(FEMM_LABEL, task_timeout_s=TASK_TIMEOUT_S) as pool:
            map_started = time.perf_counter()
            results = femm_points.run_tasks(pool, tasks, "B8: the flux map", log)
            record = pool_record(FEMM_LABEL, pool, results, time.perf_counter() - map_started, task_timeout_s=TASK_TIMEOUT_S)
    log.add_run(record)
    flux = assemble_map(results, tasks, index, (len(id_axis), len(iq_axis)))
    psi_d, psi_q = flux["psi_d"], flux["psi_q"]
    # ---- Python's own answer, on the same map (scipy): linear (the check) and a spline (reported)
    gamma = gamma_grid()
    torque_py = torque_curve(id_axis, iq_axis, psi_d, psi_q, gamma, CURRENT_A, POLE_PAIRS, INTERPOLATION)
    python_deg, python_nm = mtpa(gamma, torque_py)
    spline_deg, spline_nm = mtpa(gamma, torque_curve(id_axis, iq_axis, psi_d, psi_q, gamma, CURRENT_A, POLE_PAIRS, SPLINE))
    ctrl_d, ctrl_q = analytic_map(id_axis, iq_axis, CONTROL_PSI_M, CONTROL_LD, CONTROL_LQ)
    ctrl_deg, _ = mtpa(gamma, torque_curve(id_axis, iq_axis, ctrl_d, ctrl_q, gamma, CURRENT_A, POLE_PAIRS, INTERPOLATION))
    # ---- the inputs, in the scratch folder
    write_flux_map(folder / "flux_map.mat", id_axis, iq_axis, psi_d, psi_q, ctrl_d, ctrl_q, POLE_PAIRS, CURRENT_A, gamma)
    shutil.copy2(script_source, folder / f"{SCRIPT_NAME}.m")
    inputs = [folder / "flux_map.mat", folder / f"{SCRIPT_NAME}.m"]
    inputs_before = {p.name: sha256(p) for p in inputs}
    prefdir = scratch_prefdir(folder)
    prefdir.mkdir()
    # ---- MATLAB: one run, under the slot, watched from outside
    table_before = process_table()
    others_before = processes_like(table_before)
    tracker = JobTracker(os.getpid(), table_before)
    zero = time.perf_counter()
    watcher, sampler, classes = Watcher(tracker, zero), SlotSampler(zero), ClassSampler(tracker, zero)
    run_started_ns = time.time_ns()
    watcher.start()
    sampler.start()
    classes.start()
    try:
        with temporary_environment(MATLAB_PREFDIR=prefdir):
            matlab_run = run_program(log, matlab_command(matlab, folder), label=LABEL, cwd=folder, timeout_s=TIMEOUT_S, threads=THREADS, name="matlab")
        returned_at = time.perf_counter() - zero
    finally:
        watcher.stop()
        sampler.stop()
        classes.stop()
    table_after = process_table()
    job_text = tracker.summary(returned_at, looks=watcher.looks, failed=watcher.failed, last_error=watcher.last_error)
    members = tracker.members()
    (folder / "job_members.txt").write_text(job_text + "\n\n" + "\n".join(f"{m.pid:>7}  {m.name:<28} parent {m.ppid:>7}  first {m.first_seen:7.2f} s  last {m.last_seen:7.2f} s"
                                                                       for m in members) + "\n", encoding="utf-8")
    job = {**classify_job(members, os.getpid()), "alive_after": [m.name for m in tracker.alive(table_after)],
           "members": [(m.name, round(m.first_seen, 2), round(m.last_seen, 2)) for m in members], "empty_after_s": None, "looks": watcher.looks, "failed_looks": watcher.failed}
    last_matlab = tracker.last_seen("matlab.exe")
    if last_matlab is not None:
        job["empty_after_s"] = round(returned_at - last_matlab, 3)
    first_matlab = min((m.first_seen for m in members if m.name.lower() == "matlab.exe"), default=None)
    window = holders_in_window(sampler.samples, first_matlab, last_matlab) if first_matlab is not None and last_matlab is not None else {"looks": 0, "holders": set()}
    answer = read_answer(folder)
    # ---- after: the same lists and hashes
    inputs_after = {p.name: sha256(p) if p.is_file() else None for p in inputs}
    reference_after = tree_snapshot(reference_dir)
    reference_diff = tree_diff(reference_before, reference_after)
    prefdir_after = tree_snapshot(real) if prefdir_before is not None and real is not None else None
    prefdir_diff = tree_diff(prefdir_before, prefdir_after) if prefdir_after is not None else {"added": [], "removed": [], "changed": []}
    scratch_files = [p for p in prefdir.rglob("*") if p.is_file()]
    outside = changed_since([Path(os.environ[v]) / "MathWorks" for v in ("LOCALAPPDATA", "APPDATA", "ProgramData") if os.environ.get(v)],
                            run_started_ns - 2_000_000_000, exclude=[real] if real is not None else [])
    in_project = changed_since([paths.P], run_started_ns - 2_000_000_000, exclude=[paths.path("out"), paths.P / "runs", paths.P / "femm42-out", paths.P / ".venv"],
                               skip_names={"__pycache__", ".pytest_cache", ".git"})
    (folder / "outside_writes.txt").write_text(
        f"Files under MathWorks' folders (AppData\\Local, AppData\\Roaming without the real prefdir, ProgramData) changed since the run began: {len(outside)}\n"
        + "\n".join(f"{f['bytes']:>10,}  {f['path']}" for f in outside) + f"\n\nFiles changed under the project folder (without out\\, runs\\, femm42-out\\, .venv): {len(in_project)}\n"
        + "\n".join(f"{f['bytes']:>10,}  {f['path']}" for f in in_project) + "\n", encoding="utf-8")
    # ---- the numbers, against their expectations
    b1, b1_reason = femm_points.b1_reference()
    figure = folder / "figure.png"
    context = make_context(flux=flux, id_axis=id_axis, iq_axis=iq_axis, gamma=gamma, torque_py=torque_py, python_deg=python_deg, spline_deg=spline_deg, ctrl_deg=ctrl_deg,
                           answer=answer, matlab_run=matlab_run, figure=figure, inputs_before=inputs_before, inputs_after=inputs_after, reference_before=reference_before,
                           reference_after=reference_after, reference_diff=reference_diff, script_before=script_before, script_after=sha256(script_source), real=real,
                           prefdir_before=prefdir_before, prefdir_after=prefdir_after, prefdir_diff=prefdir_diff, prefdir=prefdir, scratch_files=scratch_files, job=job, window=window, sampler=sampler, others_before=others_before,
                           others_after=processes_like(table_after), run_record=log.runs[-1], raised=classes.left_below_normal(), b1=b1, b1_reason=b1_reason, outside=outside, in_project=in_project)
    results_dict = {"case": CASE, "title": TITLE, "params": params(mg.Params(barrier="none")), "matlab": {"command": [str(c) for c in matlab_command(matlab, folder)],
                    "prefdir_env": str(prefdir), "returncode": matlab_run.returncode, "timed_out": matlab_run.timed_out, "seconds": round(matlab_run.seconds, 2),
                    "release": (answer.get("info") or {}).get("matlab_release"), "version": (answer.get("info") or {}).get("matlab_version")},
                    "map": {"id_A": id_axis, "iq_A": iq_axis, "psi_d_Wb": psi_d, "psi_q_Wb": psi_q, "positions_el_deg": list(POSITIONS_EL), "solves": sum(len(t["points"]) for t in tasks),
                            "tasks": len(tasks), "nodes": sorted({r["nodes"] for r in results if r["nodes"]}), "interpolation": INTERPOLATION},
                    "answer": {"matlab_mtpa_deg": answer.get("mtpa_deg"), "matlab_mtpa_Nm": answer.get("mtpa_Nm"), "python_mtpa_deg": python_deg, "python_mtpa_Nm": python_nm,
                               "python_spline_mtpa_deg": spline_deg, "python_spline_mtpa_Nm": spline_nm, "control_closed_form_deg": analytic_mtpa_deg(CONTROL_PSI_M, CONTROL_LD, CONTROL_LQ),
                               "control_matlab_deg": answer.get("ctrl_mtpa_deg"), "control_python_deg": ctrl_deg, "b1": b1, "info": answer.get("info")},
                    "files": {"inputs_before": inputs_before, "inputs_after": inputs_after, "reference_files": len(reference_before), "reference_diff": reference_diff,
                              "script_before": script_before},
                    "prefdir": {"real": str(real), "files": None if prefdir_before is None else len(prefdir_before), "diff": prefdir_diff, "scratch_files": len(scratch_files),
                                "scratch_bytes": sum(p.stat().st_size for p in scratch_files)},
                    "priority_classes_over_time": {f"{name} (pid {pid})": seq for (pid, name), seq in sorted(classes.sequences.items())},
                    "job": {**job, "summary": job_text, "slot_looks_inside": window["looks"], "slot_holders_inside": [list(h) if h else None for h in window["holders"]],
                            "slot_looks_total": len(sampler.samples)},
                    "processes_like_matlab": {"before": others_before, "after": processes_like(table_after)},
                    "outside_writes": {"count": len(outside), "bytes": sum(f["bytes"] for f in outside), "files": [{"path": f["path"], "bytes": f["bytes"]} for f in outside], "in_project": [f["path"] for f in in_project]},
                    "pool": record, "seconds": {"matlab": round(matlab_run.seconds, 2)}}
    log.say(f"B8: MATLAB {release_name(results_dict['matlab']['release'])}: MTPA {answer.get('mtpa_deg')} deg (Python {python_deg:.2f}, spline {spline_deg:.2f}); job {job['members']}; "
            f"exit {matlab_run.returncode} in {matlab_run.seconds:.1f} s")
    results_dict["checks"] = make_checks(context)
    results_dict["seconds_total"] = round(time.perf_counter() - started, 1)
    return results_dict


def make_context(*, flux, id_axis, iq_axis, gamma, torque_py, python_deg, spline_deg, ctrl_deg, answer, matlab_run, figure, inputs_before, inputs_after, reference_before,
                 reference_after, reference_diff, script_before, script_after, real, prefdir_before, prefdir_after, prefdir_diff, prefdir, scratch_files, job, window, sampler,
                 others_before, others_after, run_record, raised, b1, b1_reason, outside, in_project) -> dict:
    """Everything make_checks needs, from what run() measured (see good_context in tests/test_bench_b8.py for its shape)."""
    info = answer.get("info") or {}
    psi_d, psi_q = flux["psi_d"], flux["psi_q"]
    matlab_deg = answer.get("mtpa_deg", float("nan")) if answer.get("error") is None else float("nan")
    torque_matlab = answer.get("torque_Nm") if answer.get("error") is None else None
    peak = float(np.max(np.abs(torque_py)))
    curve_difference = float(np.max(np.abs(torque_matlab - torque_py)) / peak) if torque_matlab is not None and np.shape(torque_matlab) == np.shape(torque_py) else float("inf")
    mid = len(iq_axis) // 2
    ld = (psi_d[-1, mid] - psi_d[0, mid]) / (id_axis[-1] - id_axis[0])                  # d(psi_d)/d(id) across the whole id range at mid iq [H]
    lq = (psi_q[len(id_axis) // 2, -1] - psi_q[len(id_axis) // 2, 0]) / (iq_axis[-1] - iq_axis[0])          # d(psi_q)/d(iq) across the whole iq range at mid id [H]
    id_c, iq_c = circle(gamma)
    context = {
        "angle": {"matlab_deg": matlab_deg, "python_deg": python_deg, "curve_difference": curve_difference, "angles": len(gamma),
                  "torque_matlab_Nm": float(np.max(torque_matlab)) if torque_matlab is not None else float("nan"), "torque_python_Nm": float(np.max(torque_py))},
        "control": {"closed_form_deg": analytic_mtpa_deg(CONTROL_PSI_M, CONTROL_LD, CONTROL_LQ), "matlab_deg": answer.get("ctrl_mtpa_deg", float("nan")) if answer.get("error") is None
                    else float("nan"), "python_deg": ctrl_deg},
        "run": {"returncode": matlab_run.returncode, "timed_out": matlab_run.timed_out, "seconds": matlab_run.seconds, "answer_error": answer.get("error"),
                "figure_bytes": figure.stat().st_size if figure.is_file() else 0, "files": {n: (figure.parent / n).is_file() for n in ("answer.mat", "answer.json", "figure.png")}},
        "files": {"inputs": len(inputs_before), "inputs_changed": [n for n in inputs_before if inputs_after.get(n) != inputs_before[n]], "reference_files": len(reference_before),
                  "reference_changed": reference_diff["added"] + reference_diff["removed"] + reference_diff["changed"], "script_changed": script_after != script_before,
                  "unreadable": sorted(set(unreadable_files(reference_before)) | set(unreadable_files(reference_after)))},
        "prefdir": {"real": str(real), "files": 0 if prefdir_before is None else len(prefdir_before), "bytes": 0 if prefdir_before is None else sum(v[0] for v in prefdir_before.values()),
                    "added": prefdir_diff["added"], "removed": prefdir_diff["removed"], "changed": prefdir_diff["changed"], "used_by_matlab": info.get("prefdir", ""), "expected": str(prefdir),
                    "scratch_files": len(scratch_files), "found": prefdir_before is not None,
                    "unreadable": sorted(set(unreadable_files(prefdir_before or {})) | set(unreadable_files(prefdir_after or {})))},
        "job": job,
        "slot": {"expected": (LABEL, os.getpid()), "looks": window["looks"], "holders": sorted(window["holders"], key=str), "looks_total": len(sampler.samples),
                 "unreadable_looks": sampler.errors},
        "threads": {"matlab": info.get("max_comp_threads", -1), "limit": THREAD_LIMIT},
        "priority": {"classes_seen": run_record.get("priority_classes_seen", {}), "looks": run_record.get("looks_at_the_processes", 0), "raised": raised},
        "map": {"points": (len(id_axis), len(iq_axis)), "solves": None, "max_current_error": flux["current_error"], "circle_inside": bool(id_c.min() >= id_axis[0] and id_c.max() <= id_axis[-1]
                and iq_c.min() >= iq_axis[0] and iq_c.max() <= iq_axis[-1]), "psi_m_mWb": 1e3 * float(psi_d[-1, 0]), "ld_mH": 1e3 * float(ld), "lq_mH": 1e3 * float(lq),
                "psi_q_at_iq0_mWb": 1e3 * float(np.max(np.abs(psi_q[:, 0]))), "psi_q_max_mWb": 1e3 * float(np.max(np.abs(psi_q))), "torque_at_gamma0_Nm": float(torque_py[0]),
                "torque_at_gamma90_Nm": float(torque_py[-1]), "positions": flux["positions"]},
        "b1": b1, "b1_reason": b1_reason, "spline": {"angle_deg": spline_deg},
        "outside": {"count": len(outside), "bytes": sum(f["bytes"] for f in outside), "by_root": _by_root(outside), "in_repo": [f["path"] for f in in_project],
                    "files": [{"path": f["path"], "bytes": f["bytes"]} for f in outside]},
        "others": {"before": others_before, "after": others_after},
    }
    if b1 is not None:
        angles = sorted(g for g in b1["curve"] if 0.0 <= g <= 90.0)
        mapped = torque_curve(id_axis, iq_axis, psi_d, psi_q, np.array(angles), CURRENT_A, POLE_PAIRS, INTERPOLATION)
        context["map_vs_b1"] = {"max_difference_percent_of_peak": 100 * float(np.max(np.abs(mapped - np.array([b1["curve"][g] for g in angles])))) / b1["peak_Nm"], "angles": angles}
    return context


def _by_root(files: list[dict]) -> dict:
    """How many of the files are under each of the MathWorks folders (the first three path parts below the drive and the user: AppData\\Local\\MathWorks and so on)."""
    counts: dict[str, int] = {}
    for f in files:
        parts = Path(f["path"]).parts
        key = str(Path(*parts[:parts.index("MathWorks") + 1])) if "MathWorks" in parts else str(Path(f["path"]).parent)
        counts[key] = counts.get(key, 0) + 1
    return counts


# ============================================================================================== the checks
def make_checks(context: dict) -> list[dict]:
    """Every number of B8 against its expectation, with the source of the expectation (the shape of `context` is good_context's in tests/test_bench_b8.py)."""
    out: list[dict] = []
    angle, control, run_, files, prefdir, job = (context[k] for k in ("angle", "control", "run", "files", "prefdir", "job"))
    # ---- the brief: the two angles, and what the case adds to make them mean something
    ok = agree(angle["matlab_deg"], angle["python_deg"])
    out.append(check(f"MTPA angle: MATLAB's equals Python's (the same map, linear interpolation, the same search) within {ANGLE_TOLERANCE_DEG:g} degrees", ok,
                     f"|MATLAB - Python| <= {ANGLE_TOLERANCE_DEG:g} degrees (the brief); both interpolate the 6 x 6 map linearly (interp2 'linear' = scipy 'linear') and search the {angle['angles']:,} "
                     f"angles of the 10 A circle (a design decision)",
                     f"MATLAB {angle['matlab_deg']:.2f} degrees, Python {angle['python_deg']:.2f} degrees (difference {abs(angle['matlab_deg'] - angle['python_deg']):.3f}); the largest torque "
                     f"{angle['torque_matlab_Nm']:.5f} (MATLAB) and {angle['torque_python_Nm']:.5f} N.m (Python)", "answer.mat read back, and scipy on the same map"))
    out.append(check("torque curve: MATLAB's equals Python's at every angle", angle["curve_difference"] <= CURVE_TOLERANCE,
                     f"the largest difference of the {angle['angles']:,} torques is below {CURVE_TOLERANCE:g} of the peak: the same maths on the same numbers, so the layout of the file and the "
                     f"interpolation are the same on both sides ({OWN_LIMIT})", f"{angle['curve_difference']:.2e} of the peak", "answer.mat's torque_Nm against scipy's"))
    out.append(check("control map: MATLAB and Python both find the closed-form MTPA angle of psi_d = psi_m + Ld id, psi_q = Lq iq",
                     agree(control["matlab_deg"], control["closed_form_deg"], CONTROL_TOLERANCE_DEG) and agree(control["python_deg"], control["closed_form_deg"], CONTROL_TOLERANCE_DEG),
                     f"sin g = (-psi_m + sqrt(psi_m^2 + 8 D^2)) / (4 D), D = (Lq - Ld) I, within {CONTROL_TOLERANCE_DEG:g} degrees (the search steps {GAMMA_STEP_DEG:g}); a map read with rows for columns "
                     f"would miss it by degrees ({OWN_LIMIT})",
                     f"closed form {control['closed_form_deg']:.3f}, MATLAB {control['matlab_deg']:.2f}, Python {control['python_deg']:.2f} degrees", "the closed form, derived in the docstring"))
    ran = (run_["returncode"] == 0 and not run_["timed_out"] and all(run_["files"].values()) and run_["answer_error"] is None and run_["figure_bytes"] > 20_000)
    out.append(check("MATLAB ran: exit code 0, no timeout, answer.mat, answer.json and MATLAB's own figure.png written", bool(ran),
                     f"matlab -batch exits 0 when the script succeeds (MathWorks); under a {TIMEOUT_S} s timeout; the three files MATLAB writes are there and figure.png is a real picture",
                     f"exit code {run_['returncode']}, timed out {run_['timed_out']}, {run_['seconds']:.1f} s, files {run_['files']}, figure.png {run_['figure_bytes']:,} bytes"
                     + (f", {run_['answer_error']}" if run_["answer_error"] else ""), "MathWorks: matlab (Windows), -batch"))
    # ---- no project file and no MATLAB preference changes
    changed = files["inputs_changed"] + files["reference_changed"] + (["bench\\b8_mtpa_map.m"] if files["script_changed"] else [])
    prefdir_changed = prefdir["added"] + prefdir["removed"] + prefdir["changed"]
    unreadable = files.get("unreadable", []) + prefdir.get("unreadable", [])
    found = bool(prefdir.get("found", True))
    out.append(check("no project file and no MATLAB preference changes: the inputs, reference\\, the repository's script and MATLAB's real preferences folder are what they were",
                     not changed and not prefdir_changed and not unreadable and found and files["inputs"] >= 2 and files["reference_files"] >= 1,
                     "SHA-256 before and after: flux_map.mat and the script's copy (the inputs in the scratch folder), every file of reference\\, and bench\\b8_mtpa_map.m; and every file of MATLAB's "
                     "real preferences folder, listed with its size, time and SHA-256 (key and certificate files by size and time only), with no file added, removed or touched (the brief: 'No "
                     "file outside the scratch folder changes', narrowed by a design decision to these: MathWorks' own files change whoever starts MATLAB and are listed below); a file that "
                     "cannot be read cannot be shown unchanged",
                     f"{files['inputs']} inputs, {files['reference_files']} reference files, the script: " + ("unchanged" if not changed else f"CHANGED: {changed}")
                     + f"; MATLAB's real preferences folder {prefdir['real']} "
                     + ("was not found" if not found else f"({prefdir['files']:,} files, {prefdir['bytes'] / 1e6:.1f} MB): "
                        + ("unchanged" if not prefdir_changed else f"CHANGED: added {prefdir['added'][:5]}, removed {prefdir['removed'][:5]}, changed {prefdir['changed'][:5]}"))
                     + (f"; COULD NOT BE READ (so not shown unchanged): {unreadable[:6]}" if unreadable else ""),
                     "checks.sha256 before and after, and the listing of the real folder"))
    redirected = os.path.normcase(os.path.abspath(prefdir["used_by_matlab"] or "?")) == os.path.normcase(os.path.abspath(prefdir["expected"]))
    out.append(check("prefdir: MATLAB used the scratch preferences folder it was given (MATLAB_PREFDIR), not its own", bool(redirected),
                     "MATLAB_PREFDIR on prefdir\\ in the scratch folder (MathWorks Support Team: 'To override the default location, create an environment variable named MATLAB_PREFDIR'): "
                     "MATLAB reports that folder as its prefdir; without it the real folder was written to (a trial of 11-06 changed one of its 7,640 files, epfwk_cache-*.json)",
                     f"MATLAB's prefdir was {prefdir['used_by_matlab'] or 'not reported'} ({prefdir['scratch_files']} files there); expected {prefdir['expected']}",
                     "answer.json's prefdir"))
    # ---- one MATLAB, under the slot, gone afterwards
    already = ", ".join(f"{name} (pid {pid})" for pid, name in context["others"]["before"]) or "no MATLAB or MathWorks process"
    out.append(check("exactly one MATLAB ran: one matlab.exe launcher with its MATLAB.exe, started by this run, and no other MATLAB in the job", exactly_one_matlab(job),
                     "the job's processes, followed from the process list: one matlab.exe (a child of this process), one MATLAB.exe below it, no other MATLAB-named process (the brief: 'Exactly one MATLAB "
                     "runs'); the MATLAB processes already running (the logon service host, its monitor, its crash reporter) are not the job's and are not touched",
                     f"{job['launchers']} launcher, {job['engines']} engine, {job['hosts']} MathWorks service-host clients and {job.get('graphics', 0)} MATLABWindow.exe graphics hosts "
                     f"(MATLAB's own helpers, for the figure) started by MATLAB, other MATLAB-named: {job['other_matlab'] or 'none'}, other: {job['others'] or 'none'}; "
                     f"already running before the run: {already}",
                     "checks.matlab_check.JobTracker"))
    slot, holders = context["slot"], set(context["slot"]["holders"])
    out.append(check("under the solver slot: at every look while MATLAB ran, the holder was this run's label and pid", slot["looks"] >= 3 and holders == {tuple(slot["expected"])},
                     f"femm_queue.holder() every {POLL_S:g} s between the first and the last sighting of matlab.exe: at least 3 looks, every one naming '{slot['expected'][0]}' (pid {slot['expected'][1]})",
                     f"{slot['looks']} looks inside the run ({slot['looks_total']} in all); the holders seen: " + (", ".join(str(h) for h in sorted(holders, key=str)) or "none")
                     + (f"; {slot['unreadable_looks']} looks could not be made" if slot.get("unreadable_looks") else ""), "femm_queue.holder()"))
    out.append(check("the job is empty afterwards: nothing MATLAB started is still running when run_solver returns", job["alive_after"] == [] and job.get("failed_looks", 0) == 0,
                     "run_solver returns when the job is empty (femm_queue); no member of the job is alive in a process list taken right after, and every look at the list worked",
                     "nothing left" if not job["alive_after"] else f"STILL RUNNING: {job['alive_after']}", "checks.matlab_check.JobTracker.alive"))
    members = ", ".join(f"{n} ({a:.1f}-{b:.1f} s)" for n, a, b in job["members"])
    empties = job["empty_after_s"]
    names = {n.lower() for n, _, _ in job["members"]}
    out.append(check("the job empties at once after MATLAB (Task 4: within 0.05 s of the last matlab.exe)", empties is not None and empties <= EMPTY_TOLERANCE_S and "matlab.exe" in names
                     and names <= {"matlab.exe", "mathworksservicehost.exe", GRAPHICS_HOST},
                     f"Task 4 (task-4-report.md section 5, 8 runs, warm, no figure): matlab.exe, MATLAB.exe and MathWorksServiceHost.exe clients of MATLAB's own; empty within 0.05 s of the "
                     f"last matlab.exe ({EMPTY_TOLERANCE_S:g} s here, five looks of the watcher: {OWN_LIMIT}). This run also draws a figure, which adds MATLABWindow.exe",
                     f"empty {empties if empties is not None else 'n/a'} s after the last matlab.exe; members: {members}",
                     "Task 4's measurement, repeated", gating=False, reason="a comparison with Task 4's measurement, which was made on another run and without a figure"))
    threads = context["threads"]
    out.append(check("at most 8 threads: MATLAB was held to one computational thread", 1 <= threads["matlab"] <= threads["limit"],
                     f"-singleCompThread (MathWorks: 'Limit MATLAB to a single computational thread'); maxNumCompThreads, which the script reports, between 1 and the machine's limit of {threads['limit']} "
                     f"(config.toml [femm] workers); run_program also passes {THREADS} in OMP_NUM_THREADS", f"MATLAB reports {threads['matlab']} computational thread(s)", "answer.json: maxNumCompThreads"))
    classes = context["priority"]["classes_seen"]
    core = {name: values for name, values in classes.items() if name != GRAPHICS_HOST}
    core_seen = {c for values in core.values() for c in values}
    out.append(check("priority: the matlab.exe launcher, the MATLAB.exe engine and the service hosts were seen below normal", bool(core) and core_seen == {"below_normal"}
                     and context["priority"]["looks"] > 0,
                     "run_solver starts the program below normal and what the program starts inherits that; the class of each process was read back from the operating system while it ran",
                     "; ".join(f"{n}: {', '.join(v)}" for n, v in sorted(core.items())) + f" ({context['priority']['looks']} looks)", "bench.run_program's watcher (runs.json)"))
    raised = context["priority"].get("raised", {})
    every_seen = {c for values in classes.values() for c in values}
    shown = "; ".join(f"{n}: {', '.join(v)}" for n, v in sorted(classes.items())) + ("" if GRAPHICS_HOST in classes else f"; {GRAPHICS_HOST}: not seen (no figure window was started)")
    rule_ok = bool(classes) and every_seen == {"below_normal"} and context["priority"]["looks"] > 0
    out.append(check("priority: every process of MATLAB's job was seen below normal, the figure's MATLABWindow.exe included", rule_ok,
                     "R10: every process below normal. run_solver's job holds every process in it at the class the call asked for (JOB_OBJECT_LIMIT_PRIORITY_CLASS, femm_queue._Job): a "
                     "process that asks for a higher class, or starts a child with one, stays where it is. So no process of MATLAB's job may be seen at another class; the evidence is each "
                     "image name's classes, read from the operating system while it ran (MATLABWindow.exe is MATLAB's web-graphics window host, started when the script draws the figure)",
                     shown + f" ({context['priority']['looks']} looks)" + ("" if rule_ok or not raised else " | in the order seen: " + "; ".join(f"{k}: {' -> '.join(v)}" for k, v in list(raised.items())[:8])),
                     "bench.run_program's watcher (runs.json), and a look at each process every 0.1 s (ClassSampler)"))
    # ---- the map itself: the brief's size, FEMM's currents, the physics, B1
    m = context["map"]
    out.append(check("the map's grid: at least 5 x 5 points over id -10 to 0 A and iq 0 to 10 A, and the 10 A circle lies inside it", m["points"][0] >= 5 and m["points"][1] >= 5 and m["circle_inside"],
                     "the brief: 'at least 5 x 5'; the ruling: id in [-10, 0] A and iq in [0, 10] A, so that the 10 A circle is inside", f"{m['points'][0]} x {m['points'][1]} points; circle inside: {m['circle_inside']}",
                     "the grid, and the circle at 9,001 angles"))
    out.append(check("FEMM was driven at the requested currents: id and iq came back as asked", m["max_current_error"] <= CURRENT_ERROR_TOLERANCE_A,
                     f"the inverse Park transform and its inverse return the currents to {CURRENT_ERROR_TOLERANCE_A:g} A ({OWN_LIMIT})", f"largest difference {m['max_current_error']:.2e} A",
                     "motorgeo.FemmModel.solve's id and iq"))
    b1 = context["b1"]
    if b1 is None:
        reason = context["b1_reason"]
        out.append(skipped("magnet flux: psi_d at no current against B1's psi_PM", f"within {B1_PSI_TOLERANCE_PERCENT:g} % ({OWN_LIMIT})", reason, "B1"))
        out.append(skipped("the map's torque curve against B1's torque at each current angle", f"within {B1_CURVE_TOLERANCE_PERCENT:g} % of the peak ({OWN_LIMIT})", reason, "B1"))
        out.append(skipped("the map's MTPA angle against B1's", f"within {B1_ANGLE_TOLERANCE_DEG:g} degrees ({OWN_LIMIT})", reason, "B1"))
    else:
        out.append(compare("magnet flux: psi_d at no current against B1's psi_PM (fine mesh, 120 positions)", m["psi_m_mWb"], b1["psi_pm_mWb"], B1_PSI_TOLERANCE_PERCENT,
                           f"{b1['path']}: summary.none.psi_pm_mWb; {B1_PSI_TOLERANCE_PERCENT:g} % is {OWN_LIMIT} (a coarse mesh and two positions)", "mWb"))
        diff = context["map_vs_b1"]["max_difference_percent_of_peak"]
        out.append(check("the map's torque curve against B1's: the same machine at every current angle B1 solved", diff <= B1_CURVE_TOLERANCE_PERCENT,
                         f"the map's linear torque at B1's {len(context['map_vs_b1']['angles'])} angles differs from B1's mean torque (fine mesh, six positions) by at most {B1_CURVE_TOLERANCE_PERCENT:g} % "
                         f"of the peak ({OWN_LIMIT}); the torque of the flux map is the dq torque, B1's the air-gap element's, and their means over the positions agree",
                         f"the largest difference is {diff:.2f} % of B1's {b1['peak_Nm']:.4f} N.m peak", f"{b1['path']}: summary.none.mtpa.curve"))
        out.append(check("the map's MTPA angle (linear) against B1's (22.9 degrees on the fine mesh)", abs(angle["python_deg"] - b1["gamma_deg"]) <= B1_ANGLE_TOLERANCE_DEG,
                         f"within {B1_ANGLE_TOLERANCE_DEG:g} degrees of B1's {b1['gamma_deg']:.2f} ({OWN_LIMIT}): the torque is flat near its peak, so the angle is not sharp, and a coarse mesh and a "
                         f"6 x 6 map will not match it exactly", f"{angle['python_deg']:.2f} degrees against {b1['gamma_deg']:.2f} (difference {angle['python_deg'] - b1['gamma_deg']:+.2f})",
                         f"{b1['path']}: summary.none.mtpa.bench_method.gamma_deg", gating=False, reason="the brief compares MATLAB with Python; B1 is the physical cross-check"))
    out.append(check("saliency: Lq is above Ld, so there is a reluctance torque to find an angle for", m["lq_mH"] > m["ld_mH"] > 0,
                     "an IPM machine has Lq > Ld > 0 (the d axis crosses the magnets, the q axis stays in steel); the slopes of the map across its whole range", f"Ld {m['ld_mH']:.3f} mH, Lq {m['lq_mH']:.3f} mH",
                     "the map's psi_d along id and psi_q along iq"))
    out.append(check("the q axis: psi_q is zero without q current, and so is the torque at gamma = 90 degrees", 100 * m["psi_q_at_iq0_mWb"] / m["psi_q_max_mWb"] <= Q_AXIS_TOLERANCE_PERCENT,
                     f"the machine is symmetric about its d axis, so with iq = 0 there is no q flux: within {Q_AXIS_TOLERANCE_PERCENT:g} % of the largest psi_q ({OWN_LIMIT})",
                     f"the largest |psi_q| at iq = 0 is {m['psi_q_at_iq0_mWb']:.3f} mWb ({100 * m['psi_q_at_iq0_mWb'] / m['psi_q_max_mWb']:.2f} % of {m['psi_q_max_mWb']:.2f} mWb); "
                     f"the torque at 90 degrees is {m['torque_at_gamma90_Nm']:.4f} N.m, at 0 degrees {m['torque_at_gamma0_Nm']:.4f}", "the map"))
    # ---- reported, not gating
    out.append(check("spline: Python's MTPA angle with a cubic spline on the same map", True,
                     "what a smoother interpolation gives (the ruling: reported, with how it compares with the linear angle and with B1's)",
                     f"{context['spline']['angle_deg']:.2f} degrees with the spline, {angle['python_deg']:.2f} with the linear interpolation (difference {context['spline']['angle_deg'] - angle['python_deg']:+.2f}), "
                     + (f"B1's {b1['gamma_deg']:.2f}" if b1 else "B1 not available"), "scipy RegularGridInterpolator('cubic')", gating=False,
                     reason="the brief compares MATLAB with Python on one method; this is the spline's number"))
    outside = context["outside"]
    listing = "; ".join(f"{f['path']} ({f['bytes']:,} bytes)" for f in outside["files"]) or "none"
    out.append(check("MathWorks' own files: what changed under AppData\\Local\\MathWorks, AppData\\Roaming\\MathWorks and ProgramData\\MathWorks during the run, and nothing in the project", not outside["in_repo"],
                     "MATLAB writes its licensing and usage files and the service host's logs outside the scratch folder, where MATLAB_PREFDIR does not reach; every path is listed, and nothing of "
                     "the project (outside out\\) may change", f"{outside['count']} files, {outside['bytes'] / 1e6:.2f} MB, in " + (", ".join(f"{k}: {v}" for k, v in outside["by_root"].items()) or "none")
                     + ": " + listing + ("; project files changed: " + ", ".join(outside["in_repo"]) if outside["in_repo"] else "; no project file changed") + " (the list: outside_writes.txt)",
                     "a walk of those folders for files newer than the run's start", gating=False,
                     reason="MathWorks' own files cannot be moved, and part of the list may be the logon service host's (the brief's 'no file outside the scratch folder changes' is held on the "
                            "inputs, reference\\, the repository and the real prefdir)"))
    return out
