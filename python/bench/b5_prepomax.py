"""B5: PrePoMax 3.0.0, driven from the command line with no window: B4's deck goes in, CalculiX 2.23 runs it, B4's peak stress comes out.

What it shows, and why
    PrePoMax is the pre- and post-processor of the structural chain: the window in which a CalculiX model is drawn, loaded and looked at. Task 4 could only check its version and the
    solver its settings name, because a window cannot be driven from a script the way the solvers are. B5 asks what can be done WITHOUT the window, and does it: PrePoMax 3.0.0 has
    a command line, and a Python API behind a `--script` switch that runs in a session with no window (`-g No`). B5 hands it B4's case (the same deck B4's builder writes), lets it
    import, run CalculiX 2.23 and read the result, and holds its answer against B4's. It is sound when PrePoMax starts, takes B4's case, runs CalculiX 2.23 and reproduces B4's peak
    within 1 %, and closes cleanly; or, where a step does need the window, when the report says exactly which. None did.

PrePoMax's command-line options, with the documentation each comes from (PrePoMax's own documentation is the only source: prepomax.fs.um.si, and what ships in its folder)
    The program is started as PrePoMax.com (a small console launcher beside PrePoMax.exe: it starts PrePoMax.exe, waits, and passes its exit code on).
    From the v3.0.0 manual (https://prepomax.fs.um.si/wp-content/uploads/2026/11/PrePoMax-v3.0.0-manual.pdf, section 3 "Workflow, general options", pages 7-8):
      -f, --file             the file to open or import (a .pmx, or a mesh or geometry file: .inp, .stp, .unv ...)
      -g, --showGui          Yes (default) or No: show the window; No can only be used for regeneration (the manual), for regeneration, a script or the automation server (--help)
      -o, --overwrite        Yes or No (default): overwrite the .pmx after a regeneration
      -p, --parameters       overwrite the model's parameters (p1=1.2, "a=1.2; b=31.4", or a .pmp file); regeneration only
          --parametersExport the .pmp file that the parameters of the .pmx given with -r are exported to
      -r, --regenerate       regenerate a .pmx: rerun all the commands of its history (a form of automation); the work directory holds the geometry and mesh files; with r1, r2 or r3
                             (--r1, --r2, --r3 in the program's own help; --r0 for all) only the pre-processing, only the analysis, only the post-processing
      -u, --unitSystem       the unit system for an import (M_KG_S_C, MM_TON_S_C, M_TON_S_C, IN_LB_S_F, UNIT_LESS; also for a .frd)
      -w, --workDirectory    the work directory
      -x, --exitAfterRegeneration   Yes (default) or No: exit after regeneration
          help, version      the help screen, the version
    From the program's own help screen (`PrePoMax.com --help`; the manual's table is older than the program) and the shipped Help\\Python\\Python_API_Overview.html ("Unattended script
    execution"; the web site's news of v3.0.0 names the feature: "Complete Python support: script execution and converting the command history to a Python script"):
          --script           a Python file to run after the application starts (after the -f file is opened), recorded as one command that a regeneration reruns
          --scriptArgs       the arguments that reach the script as sys.argv[1:] (split on spaces; sys.argv[0] is the script)
          --save             the .pmx to save after the script (so that the saved file regenerates the scripted session)
          --onPrompt         Fail or Accept: what to do when a dialog asks a question and nobody is there (default Fail for a script, Accept for a regeneration)
          --report           a machine-readable JSON record of the run (every switch, the commands and their durations, messages, prompts, the analyses and their status, the files made)
          --journal          the session's command history written out as a Python script when the run ends
          --timeout          the longest the whole run may take, in seconds; on expiry the running analyses are killed and the program exits (0 = no limit)
          --rpc, --rpcPort, --rpcTimeout   the automation (MCP) server that an MCP client can drive a session through (Off, On, ReadOnly); off in the shipped settings
          --settings, --noSettingsSave     read the settings from another file; do not write the settings file (a script run and a -g No session with --rpc never write it)
          --r0 ... --r3      regenerate everything, or only the pre-processing, the analysis, the post-processing
    Exit codes (the API overview): 0 completed; 1 command line error (or the automation server could not start); 2 model or regeneration error; 3 script error; 4 timed out; 5 a prompt
    was refused (--onPrompt Fail); 6 the PrePoMax.com launcher could not start the program; 7 cancelled with Ctrl+C.

Method (the files and the protocol)
    1. `PrePoMax.com --version` and `--help` (through bench.run_program, like every program): the version, and the options the program lists. B5 holds the list against the manual's.
    2. The case. B5 builds B4's pressure model on B4's main mesh (b4_ccx.build_model and write_deck: Gmsh in this process, the deck written in B5's own folder; B5 never writes into
       B4's) and solves the same deck with ccx directly, B4's way, as the control. The deck is PrePoMax-friendly by construction (B4's writer puts *BOUNDARY inside the step and writes
       the axis of *DLOAD CENTRIF as whole numbers, which PrePoMax's importer needs; ccx does not care).
    3. PrePoMax: `PrePoMax.com --script b5_prepomax_script.py --scriptArgs case.json -g No -x Yes -w <folder> --report report.json --journal journal.py --save model.pmx --timeout 700
       --onPrompt Fail`, through run_program (the solver slot, below normal, 4 CPUs, a timeout). The script (bench/b5_prepomax_script.py, which runs inside PrePoMax) imports the
       deck, reads back what the import made, replaces the section by one with the stack's 20 mm thickness, adds the centrifugal load, runs CalculiX through PrePoMax's own analysis
       object, opens the result, reads the largest von Mises stress, and plots it in PrePoMax's own viewport. case.json goes in; answer.json, PrePoMax's report.json, its journal
       and its CalculiX files (job\\Job-1.inp, .frd, .dat) come out; the case reads them.
    4. What the import does and does not take. It takes the mesh (6-node triangles), the node sets, the element set, the material, the two boundary conditions, the nodal forces (one load
       for each component of each loaded node, 1,600 here), the outputs. It does not take the *DLOAD CENTRIF line ("Failed to import distributed load": PrePoMax's own export writes the
       same line) and the section's thickness: the script gives both by API calls (2 calls). Without the thickness the solver uses 1 mm and the stress is wrong (145 MPa on the
       0.2 mm mesh, where the right answer is 32).
    5. settings.bin. It is hashed before and after (Task 1 patched its solver fields: ccx_static.exe and Spooles; PrePoMax rewrites it on exit in an interactive session, but a script run
       never writes it) and its solver fields are read back with checks/nrbf.py. B5 never changes it.
    6. The processes: before PrePoMax starts, any PrePoMax that is already running (somebody's own window) is recorded and left alone; afterwards no process of this run is left. Only
       run_solver's job ends anything, at its timeout. A process table that cannot be read is not an empty one: nothing can then be said about what was left, and the check fails
       with that reason (a check that cannot look does not pass). If a step of the script fails, answer.json (written in a `finally`) says which, and one failed check names it.
    7. The picture: PrePoMax's own render (the viewport renders with the window hidden): the von Mises stress, true-scale deformation (the displacement is a few micrometres: the sector
       looks as it is), saved with a transparent background; Python crops the legend and the model and sets them under a caption. PrePoMax's scale bar is left out: with true-scale
       deformation it reads 0.9 mm for the 45 mm sector (a fault of the picture, not of the numbers).

Run (from the Python to FEMM folder; the solver slot is taken for PrePoMax and for ccx, a couple of minutes)
    uv run python bench/run_bench.py --only B4,B5          (B5 compares with out\\bench\\B4\\results.json when it is there; alone it says so and skips that one comparison)
Writes out\\bench\\B5\\: results.json, figure.png, run.log, runs.json, prepomax\\ (case.json, answer.json, report.json, journal.py, model.pmx, the deck job.inp, job\\ with PrePoMax's own
CalculiX run), direct\\ (the control: B4's own way of running the same deck), and the output of --version and --help.

Numbers to expect (this PC, 2026-11-06, PrePoMax 3.0.0, CalculiX 2.23; the run's own are in results.json)
    --version, --help   "PrePoMax 3.0.0.0", exit code 0; 27 options listed: the manual's 11 and 16 more (--script ... --r3)
    the case            B4's pressure model at 10,000 rpm on the 0.1 mm mesh: the steel's 45,271 nodes and 22,130 6-node triangles (the mesh's other 25,263 nodes are the slots', which this
                        model has no elements for); the control, ccx run directly: 33.4540 MPa, 6 s, the same as B4's number to the digit
    the import          45,271 nodes, 22,130 elements, the element set STEEL, the material M19, 10 named node sets and 1,600 made for the forces, 2 boundary conditions, 1,600 concentrated
                        forces; not taken: the *DLOAD CENTRIF line and the section's thickness
    PrePoMax's run      exit code 0 ("completed"), status Ok, about 234 s in all, of which: the import 74 s, adding the centrifugal load 73 s, the analysis 80 s (CalculiX itself 6 s: the rest
                        is PrePoMax writing its input file and reading the result back; the cost goes with the 1,600 separate loads), the other five steps 3 s; no prompt; 4 CPUs
    the peak            33.4541 MPa at node 527 (+0.0003 % from the control and from B4; the .frd read by the bench: 33.4541 MPa); the smallest 0.854 MPa
    settings.bin        10,654 bytes, SHA-256 df4cff967880fec1..., the same after the run; solver ccx_static.exe of CalculiX 2.23, type 3 (Spooles)
    priority            ccx_static.exe, prepomax.exe, PrePoMax.com and conhost.exe all below normal; no process of the run is left
    the picture         PrePoMax's own render of the stress with its window hidden, a PNG of about 300 kB with a transparent background

Self-check (the answers follow).
  1. Why is the model's thickness set by the script, when the deck had it?
  2. PrePoMax's peak stress agrees with ccx's own run to 0.0003 %. Does that show that PrePoMax is right?
  3. settings.bin has the same hash after the run. Is that a property of PrePoMax or of how it was started?
Answers.
  1. The importer reads the section but not its thickness (found on 2026-11-06); a plane-stress element without one is 1 mm thick, so the body force, which scales with the thickness, falls
     by 20 and the nodal forces, which do not, are then 20 times too strong: 145 MPa instead of 32. The thickness is the one number of the deck that cannot be left to the import.
  2. It shows that PrePoMax hands CalculiX the same problem: the same mesh, material, supports and loads, so the same answer. It is the translation that is tested, not CalculiX (B4 and
     Task 4 test that). A wrong import (a missing thickness) gave 145 MPa: the check would have caught it.
  3. Of how it was started: an interactive session rewrites the file on exit (the window's position), a `--script` run never does (the API overview says so and --noSettingsSave is implied).
     The case would also pass a changed file whose solver fields are unchanged, and says which it is.
"""
from __future__ import annotations

import json
import re
import shutil
import time
from pathlib import Path

import numpy as np

import checks
import motorgeo as mg
import paths
from checks import ccx_check, nrbf, prepomax_check
from checks.matlab_check import process_table

from . import OWN_LIMIT, Log, check, compare, legend_and_model_figure, params, root, run_program, skipped
from . import b4_ccx as b4

CASE = "B5"
TITLE = "PrePoMax 3.0.0 from the command line (no window): B4's deck, CalculiX 2.23, B4's peak"

TOLERANCE_PERCENT = 1.0                 # the brief: PrePoMax reproduces B4's peak within 1 %
VARIANT = b4.VERDICT_VARIANT            # B4's pressure model ...
RPM = b4.SPEEDS_RPM[0]                  # ... at 10,000 rpm ...
H_MIN = b4.H_MIN_MAIN                   # ... on its main mesh
UNIT_SYSTEM = "MM_TON_S_C"              # PrePoMax's names (the unit system of the deck: mm, N, tonne, s)
MODEL_SPACE = "PlaneStress"
STEP_NAME = "Step-1"                    # what the importer calls the deck's step
MATERIAL, ELEMENT_SET = "M19", "STEEL"  # B4's deck
SCRIPT_NAME = "b5_prepomax_script.py"
CPUS = b4.CCX_THREADS
SMALL_TIMEOUT_S = 60                    # --version and --help
PMX_TIMEOUT_S = 900                     # run_solver's limit for the PrePoMax run
PMX_OWN_TIMEOUT_S = 700                 # PrePoMax's own --timeout: it ends its analyses and exits first
SOLVER_TIMEOUT_S = 500                  # the analysis object's own timeout
MAGNIFICATION = 2
EXPECTED_SOLVER_TYPE = 3                # Task 1: settings.bin's _calculixDefaultSolverType, 3 = Spooles
SAME_PEAK_PERCENT = 0.01                # PrePoMax's own number against the .frd it wrote (the bench's own limit)
GAP_PX = 30                             # the runs of a PrePoMax picture (legend, model, scale bar, info box) are separated by at least this many empty pixels

MANUAL = "the v3.0.0 manual (prepomax.fs.um.si/wp-content/uploads/2026/11/PrePoMax-v3.0.0-manual.pdf), section 3, pp. 7-8"
HELP = "PrePoMax.com --help, and Help\\Python\\Python_API_Overview.html ('Unattended script execution') in PrePoMax's folder"
OPTIONS = (
    {"option": "-f, --file", "does": "the file to open or import (a .pmx, or a mesh or geometry file)", "source": MANUAL},
    {"option": "-g, --showGui", "does": "Yes (default) or No: show the window; No runs without one, for regeneration, a script or the automation server", "source": MANUAL + "; " + HELP},
    {"option": "-o, --overwrite", "does": "Yes or No (default): overwrite the .pmx after a regeneration", "source": MANUAL},
    {"option": "-p, --parameters", "does": "overwrite the model's parameters (p1=1.2, \"a=1.2; b=31.4\", or a .pmp file); regeneration only", "source": MANUAL},
    {"option": "--parametersExport", "does": "the .pmp file the parameters of the .pmx given with -r are exported to", "source": MANUAL},
    {"option": "-r, --regenerate", "does": "regenerate a .pmx: rerun its command history (r1 pre-processing only, r2 analysis only, r3 post-processing only; --r0 all)", "source": MANUAL + "; " + HELP},
    {"option": "-u, --unitSystem", "does": "the unit system of an import: M_KG_S_C, MM_TON_S_C, M_TON_S_C, IN_LB_S_F, UNIT_LESS (also for a .frd)", "source": MANUAL},
    {"option": "-w, --workDirectory", "does": "the work directory (a bare file name of -f, --script, --save, --report, --journal, --settings is taken from it)", "source": MANUAL + "; " + HELP},
    {"option": "-x, --exitAfterRegeneration", "does": "Yes (default) or No: exit after the regeneration or the script", "source": MANUAL},
    {"option": "--help, --version", "does": "the help screen; the version", "source": MANUAL},
    {"option": "--script", "does": "a Python file to run after start-up, recorded as one command of the history", "source": HELP},
    {"option": "--scriptArgs", "does": "the arguments that reach the script as sys.argv[1:]", "source": HELP},
    {"option": "--save", "does": "the .pmx saved after the script, so that it regenerates the scripted session", "source": HELP},
    {"option": "--onPrompt", "does": "Fail or Accept: what to do with a dialog's question when nobody is there", "source": HELP},
    {"option": "--report", "does": "a machine-readable JSON record of the run", "source": HELP},
    {"option": "--journal", "does": "the session's history written out as a Python script when the run ends", "source": HELP},
    {"option": "--timeout", "does": "the longest the whole run may take (s); on expiry the analyses are killed and it exits", "source": HELP},
    {"option": "--rpc, --rpcPort, --rpcTimeout", "does": "the automation (MCP) server for an MCP client: Off (default), On, ReadOnly; its port; its request timeout", "source": HELP},
    {"option": "--settings, --noSettingsSave", "does": "read the settings from another file; never write the settings file (implied by a script run)", "source": HELP},
    {"option": "--r0, --r1, --r2, --r3", "does": "regenerate everything, the pre-processing, the analysis, the post-processing", "source": HELP},
)
MANUAL_OPTIONS = ("file", "showGui", "overwrite", "parameters", "parametersExport", "regenerate", "unitSystem", "workDirectory", "exitAfterRegeneration", "help", "version")      # the manual's table
HELP_ONLY_OPTIONS = ("script", "scriptArgs", "save", "onPrompt", "report", "journal", "timeout", "rpc", "rpcPort", "rpcTimeout", "settings", "noSettingsSave", "r0", "r1", "r2", "r3")
EXIT_CODES = {0: "completed", 1: "command line error", 2: "model or regeneration error", 3: "script error", 4: "timed out", 5: "a prompt was refused", 6: "the launcher could not start the program",
              7: "cancelled with Ctrl+C"}


# ============================================================================================== what the program says about itself
def parse_help(text: str) -> list[str]:
    """The long option names a help screen lists (`  -f, --file ...` and `  --script ...` lines), in the order they appear, without the dashes."""
    return re.findall(r"^\s{2}(?:-\w,\s+)?--(\w+)", text, re.M)


def parse_version(text: str) -> str | None:
    """The version of the first line of `--version` ("PrePoMax 3.0.0.0": the version resource's four numbers), or None."""
    found = re.search(r"PrePoMax\s+(\d+\.\d+\.\d+(?:\.\d+)?)", text)
    return found.group(1) if found else None


def options_missing(listed: list[str]) -> dict:
    """What the help screen lacks: {"manual": the manual's options it does not list, "help_only": the help-only options it does not list}."""
    have = set(listed)
    return {"manual": [o for o in MANUAL_OPTIONS if o not in have], "help_only": [o for o in HELP_ONLY_OPTIONS if o not in have]}


# ============================================================================================== settings.bin, the processes
def enum_value(member):
    """The number inside a .NET enum as nrbf.py reads it (an object with one member, value__, for an enum such as SolverTypeEnum); a plain number is returned as it is."""
    if hasattr(member, "names") and hasattr(member, "values"):
        return dict(zip(member.names, member.values)).get("value__")
    return member


def settings_state(program: Path) -> dict:
    """settings.bin next to PrePoMax.exe: its SHA-256 and size, and its solver fields read with checks/nrbf.py (nothing is started):
    {"path", "sha256", "bytes", "calculix_exe" (as stored), "solver" (resolved), "default_solver_type" (the number of the enum: 3 is Spooles)}."""
    settings = Path(program).parent / "settings.bin"
    members = nrbf.read_class_members(settings.read_bytes(), prepomax_check.CLASS)
    return {"path": str(settings), "sha256": checks.sha256(settings), "bytes": settings.stat().st_size, "calculix_exe": members.get(prepomax_check.KEY),
            "solver": str(prepomax_check.solver_path(settings, program)), "default_solver_type": enum_value(members.get("_calculixDefaultSolverType"))}


def prepomax_processes() -> dict[int, str] | None:
    """{pid: name} of the PrePoMax programs running now (PrePoMax.exe, PrePoMax.com), from a read-only snapshot of the process table; None if the table could not be read (an empty
    table is not an unreadable one: with None nothing can be said about what a run left, and the check of it fails with that reason)."""
    try:
        table = process_table()
    except OSError:
        return None
    return {pid: name for pid, (_, name) in table.items() if "prepomax" in name.lower()}


def processes_left(before: dict | None, after: dict | None) -> dict:
    """What a run left behind: {"readable", "reason", "started_and_left" (the processes in `after` that were not in `before`: this run's own, still running), "bills" (those of `before`
    still running)}. If a snapshot is None (the process table could not be read) nothing can be said: readable is False, `reason` says which snapshot failed, and both dicts are empty."""
    if before is None or after is None:
        when = "before and after the run" if before is None and after is None else "before the run" if before is None else "after the run"
        return {"readable": False, "reason": f"the process table (a read-only snapshot, checks.matlab_check.process_table) could not be read {when}", "started_and_left": {}, "bills": {}}
    return {"readable": True, "reason": "", "started_and_left": {pid: name for pid, name in after.items() if pid not in before},
            "bills": {pid: name for pid, name in before.items() if pid in after}}


def snapshot_words(snapshot: dict | None) -> str:
    """One snapshot of the PrePoMax processes in words: the pids and names, "none", or that the process table could not be read."""
    return "unknown: the process table could not be read" if snapshot is None else str(snapshot) if snapshot else "none"


def process_words(processes: dict, key: str) -> str:
    """The processes under `key` ("started_and_left" or "bills") of a processes_left result in words: the pids and names, "none", or that the table could not be read, and why."""
    if not processes["readable"]:
        return f"unknown: {processes['reason']}"
    return str(processes[key]) if processes[key] else "none"


# ============================================================================================== B4's results, the case file, the answers
def b4_reference(mesh_nodes: int, results_path: Path | None = None) -> dict:
    """B4's peak for this case, from out\\bench\\B4\\results.json: {"ok": True, "peak", "file", "nodes"} or {"ok": False, "reason"} when the file is not there, is not B4's pressure
    model at this speed, or was made on a mesh with another number of nodes (then it is not the same problem)."""
    path = Path(results_path) if results_path is not None else root() / "B4" / "results.json"
    if not path.is_file():
        return {"ok": False, "reason": f"{path} is not there: B4 has not been run on this PC (uv run python bench/run_bench.py --only B4)"}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        run = data["runs"][f"{VARIANT}_{RPM}_h{H_MIN:g}"]
        nodes, peak = int(data["mesh"]["nodes"]), float(run["peak"]["value_mpa"])
    except (KeyError, TypeError, ValueError, json.JSONDecodeError) as error:
        return {"ok": False, "reason": f"{path} has no {VARIANT} run at {RPM} rpm on the {H_MIN} mm mesh with a peak ({type(error).__name__}: {error})"}
    if nodes != mesh_nodes:
        return {"ok": False, "reason": f"{path} is for a mesh of {nodes:,} nodes, this one has {mesh_nodes:,}: not the same problem"}
    return {"ok": True, "peak": peak, "file": str(path), "nodes": nodes}


def case_inputs(folder: Path, deck: Path, omega: float, thickness: float) -> dict:
    """The content of case.json: what PrePoMax's side needs (absolute paths as text, the units, the step, the load's speed in rad/s, the work folder of the analysis, where the answer and
    the picture go)."""
    folder = Path(folder)
    return {"deck": str(deck), "unit_system": UNIT_SYSTEM, "model_space": MODEL_SPACE, "step": STEP_NAME, "material": MATERIAL, "element_set": ELEMENT_SET, "thickness": thickness,
            "omega": omega, "cpus": CPUS, "job_dir": str(folder / "job"), "solver_timeout_s": SOLVER_TIMEOUT_S, "answer": str(folder / "answer.json"), "image": str(folder / "view.png"),
            "magnification": MAGNIFICATION}


def read_json(path, what: str) -> dict:
    """A JSON file as a dict; ValueError, saying which file and why, if it is not there or is not JSON."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"{what} {path} is not there")
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"{what} {path} is not JSON: {error}") from error


def import_findings(imported: dict, expected: dict) -> dict:
    """What PrePoMax's importer made of B4's deck against what the deck holds. `imported` is answer.json's "imported", `expected` {"nodes", "elements", "node_sets" (names),
    "loaded_nodes", "boundary_conditions"}. Returns {"mesh_ok", "sets_ok", "forces_ok", "bcs_ok", "gaps"}: `gaps` names what the importer did not take (the centrifugal load, the thickness)."""
    part = imported["parts"][0] if imported.get("parts") else {"nodes": None, "elements": None}
    forces = imported.get("loads", {}).get("concentrated_force", 0)
    gaps = []
    if "centrifugal" not in imported.get("loads", {}):
        gaps.append("the *DLOAD CENTRIF line (the importer reports 'Failed to import distributed load')")
    if any(s.get("thickness") in (None, 0) for s in imported.get("sections", [])) or not imported.get("sections"):
        gaps.append("the section's thickness (the imported section has none: 1 mm would be used)")
    return {"mesh_ok": part["nodes"] == expected["nodes"] and part["elements"] == expected["elements"],
            "sets_ok": set(expected["node_sets"]) <= set(imported.get("node_sets", {})) and ELEMENT_SET in imported.get("element_sets", {}) and MATERIAL in imported.get("materials", []),
            "forces_ok": forces == 2 * expected["loaded_nodes"], "bcs_ok": len(imported.get("boundary_conditions", [])) == expected["boundary_conditions"], "gaps": gaps}


def findings_of(answer: dict, expected: dict) -> dict | None:
    """import_findings of the read-back in answer.json; None if the script never got that far (answer.json is written in a `finally`, so a failed import step leaves an answer
    without "imported", and the failed step is in its "steps")."""
    return import_findings(answer["imported"], expected) if answer.get("imported") else None


def peak_check(prepomax: float, reference: float, name: str, source: str, **kw) -> dict:
    """PrePoMax's peak von Mises stress against a reference peak, within the brief's 1 %."""
    return compare(name, prepomax, reference, TOLERANCE_PERCENT, source, "MPa", **kw)


# ============================================================================================== the picture
def _runs(flags, gap: int = GAP_PX) -> list[tuple[int, int]]:
    """The runs of True in `flags` as (first, last) index pairs, runs closer than `gap` positions merged into one."""
    out: list[list[int]] = []
    for i in np.flatnonzero(np.asarray(flags)):
        if out and i - out[-1][1] <= gap:
            out[-1][1] = int(i)
        else:
            out.append([int(i), int(i)])
    return [(a, b) for a, b in out]


def prepomax_picture(png, caption: str, target) -> dict:
    """PrePoMax's own picture (a PNG with a transparent background, as pmx.view.save_image writes it) cut into what is wanted: the colour legend (the run of columns at the left), the model
    (the widest run of columns, and in it the tallest run of rows: the scale bar below it, which reads 0.9 mm for the 45 mm sector, is left out), set on white under `caption`.
    Returns {"made", "size", "legend", "model"} (the pixel boxes); made is False if the picture does not have a legend and a model."""
    from PIL import Image
    image = Image.open(png).convert("RGBA")
    alpha = np.asarray(image)[..., 3] > 0
    columns = _runs(alpha.any(axis=0))
    if len(columns) < 2:
        return {"made": False, "size": image.size}
    legend_cols, model_cols = columns[0], max(columns, key=lambda r: r[1] - r[0])
    if legend_cols == model_cols:
        return {"made": False, "size": image.size}
    legend_rows = _runs(alpha[:, legend_cols[0]:legend_cols[1] + 1].any(axis=1))
    model_rows = max(_runs(alpha[:, model_cols[0]:model_cols[1] + 1].any(axis=1)), key=lambda r: r[1] - r[0])
    white = Image.new("RGBA", image.size, "white")
    white.alpha_composite(image)
    flat = np.asarray(white.convert("RGB"))
    pad = lambda a, b, limit: (max(a - 8, 0), min(b + 9, limit))
    r0, r1 = pad(legend_rows[0][0], legend_rows[-1][1], flat.shape[0])
    c0, c1 = pad(*legend_cols, flat.shape[1])
    legend = flat[r0:r1, c0:c1]
    r0, r1 = pad(*model_rows, flat.shape[0])
    c0, c1 = pad(*model_cols, flat.shape[1])
    model = flat[r0:r1, c0:c1]
    legend_and_model_figure(legend, model, caption, target)
    return {"made": True, "size": image.size, "legend": list(legend.shape[:2]), "model": list(model.shape[:2])}


def picture_caption(nodes: int, peak: float, peak_node: int) -> str:
    """The caption of figure.png: four lines, none longer than bench.TITLE_CHARS (a longer one is wrapped, and leaves a word alone on a line). `nodes` is the number of nodes of the model
    PrePoMax solved (the steel's: the pressure model has no slots' elements), `peak` the largest von Mises stress in MPa and `peak_node` where PrePoMax found it."""
    return (f"PrePoMax 3.0.0, run from the command line without its window: the band-free rotor,\n"
            f"its 90-degree sector at {RPM:,} rpm. von Mises stress [MPa] of B4's {VARIANT} model,\n"
            f"{nodes:,} nodes, solved by CalculiX 2.23 through PrePoMax. Peak {peak:.2f} MPa at node {peak_node}.\n"
            f"PrePoMax's own render; its scale bar (it reads 0.9 mm here) is left out.")


# ============================================================================================== the case
def run(folder: Path, log: Log) -> dict:
    """Run B5 in `folder` (already empty). Returns the results dict whose "checks" list every verdict."""
    folder = Path(folder)
    started = time.perf_counter()
    com = paths.tool("prepomax").with_name("PrePoMax.com")
    program = paths.tool("prepomax")
    if not com.is_file():
        raise checks.Missing("PrePoMax.com", com)
    ccx = paths.tool("ccx")
    before_hash = settings_state(program)
    before_processes = prepomax_processes()
    log.say(f"B5: settings.bin {before_hash['sha256'][:16]}..., solver {before_hash['solver']}, type {before_hash['default_solver_type']}; PrePoMax processes already running (left alone): "
            f"{snapshot_words(before_processes)}")
    # ---- the program's own words: its version and its options
    version = run_program(log, [com, "--version"], label="bench B5 PrePoMax --version", cwd=folder, timeout_s=SMALL_TIMEOUT_S, threads=1, name="prepomax_version")
    helped = run_program(log, [com, "--help"], label="bench B5 PrePoMax --help", cwd=folder, timeout_s=SMALL_TIMEOUT_S, threads=1, name="prepomax_help")
    listed = parse_help(helped.stdout + helped.stderr)
    found_version = parse_version(version.stdout + version.stderr)
    # ---- B4's case in B5's own folder, and the control: ccx run on it directly, B4's way
    p = mg.Params(barrier="none")
    model = b4.build_model(p, H_MIN)
    mesh = model["mesh"]
    pmx_dir = folder / "prepomax"
    deck = b4.write_deck(pmx_dir, model, VARIANT, RPM)
    at = b4.loads_at(model, RPM)
    steel_nodes = {n for t in mesh.steel for n in t}
    direct = b4.solve(model, VARIANT, RPM, folder / "direct", log, f"bench B5 ccx direct {VARIANT} {RPM} rpm (the control)")
    direct_peak = b4.peak_of(direct["stress"], steel_nodes)[0] if direct["finished"] else None
    same_deck = checks.sha256(deck) == checks.sha256(folder / "direct" / "job.inp")
    reference = b4_reference(len(mesh.nodes))
    log.say(f"B5: B4's {VARIANT} model: the steel's {len(steel_nodes):,} nodes (of the mesh's {len(mesh.nodes):,}: the slots' are not in this model); the control (ccx directly) "
            f"{f'{direct_peak:.4f} MPa' if direct_peak is not None else 'did not finish'}, "
            f"{direct['seconds']:.1f} s; B4's own number: {reference['peak'] if reference['ok'] else reference['reason']}")
    # ---- PrePoMax: the script run, no window
    shutil.copyfile(paths.P / "bench" / SCRIPT_NAME, pmx_dir / SCRIPT_NAME)
    omega = b4.omega(RPM)
    (pmx_dir / "case.json").write_text(json.dumps(case_inputs(pmx_dir, deck, omega, p.stack), indent=1), encoding="utf-8")
    cmd = [com, "--script", SCRIPT_NAME, "--scriptArgs", "case.json", "-g", "No", "-x", "Yes", "-w", str(pmx_dir), "--report", "report.json", "--journal", "journal.py",
           "--save", "model.pmx", "--timeout", str(PMX_OWN_TIMEOUT_S), "--onPrompt", "Fail"]
    pmx = run_program(log, cmd, label="bench B5 PrePoMax --script", cwd=pmx_dir, timeout_s=PMX_TIMEOUT_S, threads=CPUS, name="prepomax")
    after_processes = prepomax_processes()
    after_hash = settings_state(program)
    results = {"case": CASE, "title": TITLE, "params": params(p), "variant": VARIANT, "rpm": RPM, "mesh_nodes": len(mesh.nodes), "options": list(OPTIONS), "exit_codes": EXIT_CODES,
               "version": {"text": (version.stdout + version.stderr).strip(), "parsed": found_version, "returncode": version.returncode},
               "help": {"options_listed": listed, "missing": options_missing(listed), "returncode": helped.returncode},
               "settings": {"before": before_hash, "after": after_hash, "expected_solver": str(ccx)},
               "processes": {"before": before_processes, "after": after_processes, **processes_left(before_processes, after_processes)},
               "direct": {"finished": direct["finished"], "peak_mpa": direct_peak, "seconds": direct["seconds"], "same_deck_as_prepomax": same_deck}, "reference": reference,
               "prepomax": {"returncode": pmx.returncode, "timed_out": pmx.timed_out, "seconds": round(pmx.seconds, 2), "command": [str(c) for c in cmd]}}
    try:
        answer = read_json(pmx_dir / "answer.json", "PrePoMax's answer")
        report = read_json(pmx_dir / "report.json", "PrePoMax's report")
    except ValueError as error:
        results["error"] = str(error)
        results["checks"] = make_checks(results, None, None, None, before_hash, after_hash, pmx, found_version, listed, False)
        results["seconds_total"] = round(time.perf_counter() - started, 1)
        return results
    frd = pmx_dir / "job" / "Job-1.frd"
    frd_peak = None
    if frd.is_file():
        stress = ccx_check.read_frd_stress(frd)
        frd_peak = b4.peak_of(stress, set(stress))[0]
    expected = {"nodes": len(steel_nodes), "elements": len(mesh.steel), "node_sets": sorted(b4.wall_sets(at["loads"])) + ["EDGE0", "EDGE90"],
                "loaded_nodes": len(b4.nodal_forces(at["loads"])), "boundary_conditions": 2}
    findings = findings_of(answer, expected)
    results.update({"answer": answer, "report": {k: report.get(k) for k in ("version", "exitCode", "exitStatus", "error", "durationSeconds", "messages", "prompts", "jobs", "artifacts", "commands")},
                    "invocation": report.get("invocation"), "import": {"expected": expected, **(findings or {})}, "frd_peak_mpa": frd_peak,
                    "peak_mpa": answer["results"]["peak"]["value"] if "results" in answer else None})
    picture = {"made": False}
    if "image" in answer and "results" in answer:
        peak = answer["results"]["peak"]
        picture = prepomax_picture(answer["image"]["path"], picture_caption(len(steel_nodes), peak["value"], peak["node"]), folder / "figure.png")
    results["picture"] = picture
    results["checks"] = make_checks(results, answer, report, findings, before_hash, after_hash, pmx, found_version, listed, picture["made"], log=log)
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    log.say(f"B5 PrePoMax: exit {pmx.returncode}, {pmx.seconds:.1f} s; job {answer.get('job', {}).get('status')}; peak {results['peak_mpa']} MPa (the .frd: {frd_peak}); the control {direct_peak}; "
            f"settings.bin unchanged: {before_hash['sha256'] == after_hash['sha256']}; processes left: {process_words(results['processes'], 'started_and_left')}")
    return results


def make_checks(results: dict, answer: dict | None, report: dict | None, findings: dict | None, before: dict, after: dict, pmx, version: str | None, listed: list[str],
                picture_made: bool, log: Log | None = None) -> list[dict]:
    """Every number of B5 against its expectation, with the source of the expectation. When PrePoMax gave no answer (answer is None) only the checks that need none are made, and the
    check of the run itself says what happened."""
    out: list[dict] = []
    out.append(check("PrePoMax starts from the command line and says its version", version is not None and version.startswith("3.0.0") and results["version"]["returncode"] == 0,
                     "PrePoMax.com --version prints the version, 3.0.0 (Task 1 installed 3.0.0; the same as checks/prepomax_check's version resource), and exits 0",
                     f"{results['version']['text']!r}, exit code {results['version']['returncode']}", "PrePoMax.com --version (documented: the manual's 'version')"))
    missing = results["help"]["missing"]
    out.append(check("the options: the manual's table and the program's own help agree, and the script options are there",
                     not missing["manual"] and not missing["help_only"] and results["help"]["returncode"] == 0,
                     f"--help lists every option of the manual's table ({len(MANUAL_OPTIONS)}) and the {len(HELP_ONLY_OPTIONS)} that only the program's help and the API overview describe "
                     f"(--script, --scriptArgs, --save, --report, --journal, --timeout, ...)",
                     f"{len(listed)} options listed; the manual's missing from --help: {missing['manual'] or 'none'}; the script options missing: {missing['help_only'] or 'none'}; "
                     f"the manual's table lacks (the program lists, the web manual does not): {', '.join(HELP_ONLY_OPTIONS)}",
                     MANUAL + "; " + HELP))
    same_hash = before["sha256"] == after["sha256"]
    fields_same = (before["solver"], before["default_solver_type"]) == (after["solver"], after["default_solver_type"])
    solver_ok = Path(after["solver"]) == Path(results["settings"]["expected_solver"]) and after["default_solver_type"] == EXPECTED_SOLVER_TYPE
    out.append(check("settings.bin is not changed, and its solver fields are the ones Task 1 set", (same_hash or fields_same) and solver_ok,
                     f"the hash of settings.bin is the same before and after the run (a script run never writes it); if it were not, its solver fields (read with checks/nrbf.py) must be "
                     f"unchanged: ccx_static.exe of CalculiX 2.23 ({results['settings']['expected_solver']}) and the default solver type {EXPECTED_SOLVER_TYPE} (Spooles)",
                     f"SHA-256 {before['sha256'][:16]}... before, {after['sha256'][:16]}... after ({'the same' if same_hash else 'CHANGED'}); {after['bytes']:,} bytes; solver {after['solver']}, "
                     f"type {after['default_solver_type']}; solver fields unchanged: {fields_same}", "checks/nrbf.py, checks/prepomax_check.py; Task 1's INSTALLED.txt"))
    if answer is None:
        out.append(check("PrePoMax runs the script and answers", False, "PrePoMax.com --script b5_prepomax_script.py ... writes answer.json and report.json",
                         f"exit code {pmx.returncode}, timed out {pmx.timed_out}; {results.get('error', '')}", "PrePoMax's --report and the script's answer.json"))
        return out
    exit_code = report.get("exitCode")
    steps = answer["steps"]
    invocation = report.get("invocation") or {}
    processes = results["processes"]
    out.append(check("PrePoMax runs without its window, and closes cleanly", invocation.get("showGui") == "No" and exit_code == 0 and report.get("exitStatus") == "Ok" and pmx.returncode == 0
                     and not pmx.timed_out and all(s["ok"] for s in steps) and processes["readable"] and not processes["started_and_left"] and not report.get("prompts"),
                     "-g No (no window), exit code 0 and status Ok in PrePoMax's own report, no prompt refused, every step of the script worked, and no process of this run is left (a process "
                     "table that cannot be read fails this: a check that cannot look does not pass); any PrePoMax that was already running (somebody's own window) is untouched",
                     f"showGui {invocation.get('showGui')}, exit code {exit_code} ({EXIT_CODES.get(exit_code)}), status {report.get('exitStatus')}, {report.get('durationSeconds', 0):.1f} s; "
                     f"{sum(s['ok'] for s in steps)} of {len(steps)} steps: " + ", ".join(f"{s['step']} ({s['seconds']:.1f} s)" if "seconds" in s else s["step"] for s in steps) + f"; prompts {len(report.get('prompts') or [])}; settings saved by the run: "
                     f"{invocation.get('settingsSaved')}; processes left by this run: {process_words(processes, 'started_and_left')}; PrePoMax processes that were already running "
                     f"and still are: {process_words(processes, 'bills')}", "PrePoMax's --report; the process table before and after (checks.matlab_check.process_table, read only)"))
    needed = [key for key in ("imported", "job", "results", "image") if not answer.get(key)]
    if needed:                                          # a step failed: answer.json (written in a `finally`) says which, and has none of what the later steps would have made
        failed = [s for s in steps if not s["ok"]]
        out.append(check("PrePoMax's script made what the checks below read", False,
                         "answer.json holds the read-back of the import, the job, the results and the picture's file: every step of the script worked",
                         "missing from answer.json: " + ", ".join(needed) + "; " + (f"failed step: {failed[0]['step']}: {failed[0].get('error')}" if failed else
                                                                                       f"no failed step is recorded (the script's error: {answer.get('error')})"),
                         "answer.json (b5_prepomax_script.py writes it in a `finally`)"))
        return out
    imported = answer["imported"]
    expected = results["import"]["expected"]
    out.append(check("PrePoMax takes B4's case: the mesh, the sets, the material, the supports and the nodal forces",
                     findings["mesh_ok"] and findings["sets_ok"] and findings["forces_ok"] and findings["bcs_ok"],
                     f"the imported part has B4's {expected['nodes']:,} nodes and {expected['elements']:,} 6-node triangles; the node sets {', '.join(expected['node_sets'][:2])}... and the "
                     f"wall sets, the element set {ELEMENT_SET}, the material {MATERIAL}, {expected['boundary_conditions']} boundary conditions and {2 * expected['loaded_nodes']} nodal forces "
                     f"(two components for each of {expected['loaded_nodes']} nodes) are there",
                     f"part {imported['parts'][0]['name']}: {imported['parts'][0]['nodes']:,} nodes, {imported['parts'][0]['elements']:,} elements ({imported['parts'][0]['type']}); node sets "
                     f"{len(imported['node_sets'])} named + {imported.get('auto_created_node_sets', 0)} made for the forces; element sets {imported['element_sets']}; materials "
                     f"{imported['materials']}; boundary conditions {len(imported['boundary_conditions'])}; loads {imported['loads']}; outputs {imported['field_outputs']} and "
                     f"{imported['history_outputs']}", "the script's read-back of the model after pmx.import_file"))
    out.append(check("the importer reads the whole deck", not findings["gaps"],
                     "the import makes a complete model: nothing the deck says has to be said again through the API",
                     "what the importer did not take, found on 2026-11-06: " + ("; ".join(findings["gaps"]) if findings["gaps"] else "nothing") + ". The script gives both through the API "
                     "(sections.add_solid with the thickness, step.loads.add_centrifugal): without the thickness the result is 145 MPa, not 32 (the 0.2 mm mesh)",
                     "pmx.import_file's read-back", gating=False, reason="PrePoMax's importer has two gaps; the script bridges them and the result shows it"))
    job = answer["job"]
    settings_solver = Path(after["solver"])
    banner = re.search(r"CalculiX Version\s+(\d+\.\d+)", pmx.stdout)
    out.append(check("PrePoMax runs CalculiX 2.23, the solver of its settings", job["status"] == "OK" and Path(job["executable"]) == settings_solver and banner is not None
                     and banner.group(1) == "2.23" and job["number_of_cpus"] <= paths.setting("femm", "workers"),
                     f"the analysis ends OK; the program it ran is the one settings.bin names (ccx_static.exe in CalculiX 2.23's folder) and its banner says 2.23; at most "
                     f"{paths.setting('femm', 'workers')} CPUs",
                     f"status {job['status']}; executable {job['executable']}; CalculiX {banner.group(1) if banner else 'no banner'}; {job['number_of_cpus']} CPUs; argument {job['argument']}",
                     "pmx.jobs, PrePoMax's printout of the solver"))
    peak = answer["results"]["peak"]["value"]
    out.append(peak_check(peak, results["direct"]["peak_mpa"], "PrePoMax's peak von Mises stress reproduces the same deck solved by ccx directly",
                          "the control: B4's own way of running this deck, in B5's folder (the two decks are identical: " + str(results["direct"]["same_deck_as_prepomax"]) + ")")
               if results["direct"]["peak_mpa"] is not None else check("PrePoMax's peak against the control", False, "the control run finished", "the control did not finish", "b4_ccx.solve"))
    if results["reference"]["ok"]:
        out.append(peak_check(peak, results["reference"]["peak"], "PrePoMax reproduces B4's peak von Mises stress",
                              f"the brief: B4's pressure model at 10,000 rpm on its {H_MIN} mm mesh, {results['reference']['file']}"))
    else:
        out.append(skipped("PrePoMax reproduces B4's peak von Mises stress", f"within {TOLERANCE_PERCENT} % of B4's number", results["reference"]["reason"], "out\\bench\\B4\\results.json"))
    frd = results["frd_peak_mpa"]
    out.append(compare("PrePoMax's own number agrees with the .frd it wrote", peak, frd, SAME_PEAK_PERCENT, f"the .frd of PrePoMax's CalculiX run, read by the bench ({OWN_LIMIT}: the same numbers "
                                                                                                            "through two readers)", "MPa") if frd else
               check("PrePoMax's own number agrees with the .frd it wrote", False, "job\\Job-1.frd is there", "no .frd", "PrePoMax's job folder"))
    seen = log.runs[-1]["priority_classes_seen"] if log is not None and log.runs else {}
    classes = sorted({c for cs in seen.values() for c in cs})
    out.append(check("everything PrePoMax started ran below normal priority", classes == ["below_normal"],
                     "R10: every program the bench starts, and every program that program starts, runs at below-normal priority (the job inherits it)",
                     "classes seen while the run was on: " + ", ".join(f"{name} {cs}" for name, cs in seen.items()), "bench.run_program's watcher (the classes read back from the operating system)"))
    out.append(check("the picture: figure.png is PrePoMax's own render of the stress", picture_made and Path(answer["image"]["path"]).is_file(),
                     "PrePoMax's viewport (with its window hidden) draws the von Mises stress and saves a PNG; figure.png is its legend and its model under a caption",
                     f"saved {answer['image']['bytes']:,} bytes; legend and model cut out: {picture_made}", "pmx.view.save_image"))
    out.append(skipped("an automation (MCP) session of PrePoMax: --rpc On", "an MCP client drives a session through PrePoMax's automation server",
                       "not made: it needs an MCP client on the other end and the session's access token, and it is off in the shipped settings; the script run needs neither, and does the same work",
                       HELP))
    return out
