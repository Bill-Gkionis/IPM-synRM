"""gmsh_check: Gmsh meshes the default rotor from Python, and the Gmsh app prints the version we installed.

What it does, and why
    Gmsh is the mesh maker for everything that is not FEMM: Elmer and CalculiX read its meshes. It comes twice here, and a row
    for each, because the two are separate installs that can drift apart:
      * the Python module `gmsh` (pinned to 4.15.2 in pyproject.toml), which motorgeo.to_gmsh drives to mesh the motor's cross
        section: every region of the sector is a surface with a physical group of its own, rotor and stator are meshed apart (the
        sliding interface of the air gap), and the file is written in the MSH 2.2 format that ElmerGrid reads;
      * the app gmsh.exe (D:\\Apps\\Gmsh), which the pipeline will use for pictures and for meshing from the command line.
    The first row passes when the default rotor (motorgeo.Params()) meshes with more than 1,000 nodes (it is 9,176 at
    the mesh sizes of the regions); the second when `gmsh -version` prints 4.15.2.

Method
    1. mg.build(Params(), gap="sliding") draws the sector with the air gap split into two sliding arcs, mg.to_gmsh meshes it and
       writes out\\check_tools\\gmsh\\rotor.msh. No program is started: the Python module is a library, a second of CPU.
    2. The app is started through femm_queue.run_solver with the one argument -version. Its exit status is not looked at (many
       programs use it for something else); the text it prints is.

A side effect of the Python module that this row found (2026-11-06)
    gmsh.initialize() cuts the PATH of the running process down to its first few entries and finalize() does not restore it
    (1,864 characters became 313). os.environ still shows the whole PATH, so nothing looks wrong, but a program started by
    its bare name afterwards (powershell.exe, which lives only on the PATH) is no longer found. Programs started by their full
    path, as run_solver does, and children (they get a copy of os.environ) are not affected. checks.intact_path() puts the
    PATH back; this row and the Elmer row use it round their Gmsh calls, and check_tools.py round every check. A script of the
    pipeline that meshes with Python and then starts a program by name needs the same.

Numbers to expect (2026-11-06): the module 4.15.2 meshes the sector into 9,176 nodes (a 947 kB MSH 2.2 file, 5 boundary groups)
    in 0.3 s; the app prints 4.15.2 in 0.1 s.

Self-check (the answers follow).
  1. Why does the first row not start gmsh.exe?
  2. What would a Gmsh app of 4.16 and a Python module of 4.15.2 break, and where would this row show it?
Answers.
  1. The Python module is the Gmsh library inside the Python process; there is no separate program to start (the pip wheel
     carries the library with it, and does not need gmsh.exe).
  2. Nothing at once: both read and write the same files. But a mesh made by one could differ from a mesh made by the other, and
     a pipeline that mixed them would be hard to reproduce. The second row fails (4.16 is not 4.15.2) and says so.
"""
from __future__ import annotations

import re
import time
from pathlib import Path

import femm_queue

from . import FAIL, PASS, Row, fresh, intact_path, need_tool, save_log

EXPECT_VERSION = "4.15.2"
MIN_NODES = 1000


def run(folder: Path):
    fresh(folder)
    import gmsh
    import motorgeo as mg

    started = time.perf_counter()
    params = mg.Params()
    mesh_file = folder / "rotor.msh"
    with intact_path():                             # gmsh.initialize() cuts the process's PATH short (see checks.intact_path)
        nodes, boundaries = mg.to_gmsh(params, mg.build(params, gap="sliding"), str(mesh_file))
    seconds = round(time.perf_counter() - started, 2)
    header = mesh_file.read_text(encoding="ascii", errors="replace")[:40].split()
    msh_version = header[1] if len(header) > 1 else "?"
    ok = nodes > MIN_NODES and mesh_file.stat().st_size > 0 and msh_version == "2.2"
    yield Row("Gmsh", "motorgeo.to_gmsh on the default rotor (the Python module)", f"meshes with more than {MIN_NODES:,} nodes",
              f"gmsh {gmsh.__version__}: {nodes:,} nodes, MSH {msh_version} file of {mesh_file.stat().st_size / 1e3:,.0f} kB, "
              f"{len(boundaries)} boundary groups", PASS if ok else FAIL, seconds)

    exe = need_tool("gmsh", case="gmsh -version (the app)")
    started = time.perf_counter()
    result = femm_queue.run_solver([exe, "-version"], label="check_tools gmsh -version", cwd=folder, timeout_s=60, threads=1)
    seconds = round(time.perf_counter() - started, 2)
    save_log(folder, "gmsh_version", result)
    text = (result.stdout + result.stderr).strip()
    found = re.search(r"\b\d+\.\d+\.\d+\b", text)
    version = found.group(0) if found else (f"no version in: {text[:80]!r}" if not result.timed_out else "timed out")
    yield Row("Gmsh", "gmsh -version (the app)", f"the app prints {EXPECT_VERSION}", f"the app prints {version}",
              PASS if version == EXPECT_VERSION else FAIL, seconds)
