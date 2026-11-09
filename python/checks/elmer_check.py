"""elmer_check: Elmer solves the magnetic field of a round wire (MagnetoDynamics2D), and Ampere's law says what it must be.

What it does, and why
    Elmer is the open solver for the fields FEMM is not used for. Its 2D magnetic solver, MagnetoDynamics2D, finds the vector
    potential Az of a field that has only an x-y plane (the field of the motor's cross-section) from the currents in it. The
    known-answer case is the simplest conductor there is: a long round copper wire of radius R carrying a current I, alone in
    air. Ampere's law gives its field exactly. B points around the wire, and
        inside  (r <= R):  B(r) = mu0 I r / (2 pi R^2)    (the current enclosed grows as r^2, so B grows as r)
        outside (r >= R):  B(r) = mu0 I / (2 pi r)        (all of I is enclosed, B falls as 1/r)
    With R = 10 mm and I = 100 A (mu0 = 4 pi 1e-7): B(R) = 2.0 mT, B(R/2) = 1.0 mT.

The chain the pipeline will use (every program through femm_queue.run_solver, in the case folder)
    1. Gmsh (the Python module, in-process) draws a disk of radius R (the conductor, physical group 1) inside a disk of radius
       10 R (air, group 2; its edge is the boundary, group 1) and meshes them in the MSH 2.2 format that ElmerGrid reads, finest
       (0.2 mm) at the wire's surface and growing to 8 mm at the edge: 21,153 nodes, 42,225 triangles.
    2. ElmerGrid 14 2 case.msh turns the mesh into Elmer's own files (the folder case\\: mesh.header, .nodes, .elements, .boundary).
    3. ElmerSolver case.sif solves it. The input file (case_sif()) has three solvers:
         MagnetoDynamics2D            solves  -div( (1/mu) grad Az ) = Jz   for Az, with Az = 0 on the outer edge;
         MagnetoDynamicsCalcFields    computes B = curl(Az z-hat) = (dAz/dy, -dAz/dx) at the nodes from Az;
         SaveLine                     writes B along the x axis from 0 to 2R at 40 equal steps (points every R/20, so r = R/2 and
                                      r = R are points 10 and 20) to line.dat, with the names of its columns in line.dat.names.
       The current density is the same in the whole wire, Jz = I / (pi R^2) = 318,309.886 A/m^2 (the current is distributed
       evenly because the wire is a solid DC conductor; a copper wire at 50 Hz would crowd it to the skin, which a magnetostatic
       solver does not model). Az = 0 on a circle round a single wire is the exact free-space solution up to a constant, so the
       finite domain costs no accuracy.
    4. This module reads line.dat and compares |B| = sqrt(Bx^2 + By^2) at r = R/2 and at r = R with Ampere's law. Both must be
       within 2 %.

Why 0.2 mm and not coarser
    Elmer's elements are linear, so B = curl(A) is constant on each triangle and the nodal value is an average of the triangles
    around the node. B(r) has a corner at r = R (it grows inside, falls outside), and an average rounds a corner off. The error at
    r = R falls in proportion to the element size (measured 2026-11-06): -2.21 % at 0.8 mm, -1.09 % at 0.4 mm, -0.48 % at 0.2 mm,
    -0.22 % at 0.1 mm; at r = R/2, where B is smooth, it is 0.25 % or less at every size. 0.2 mm leaves three quarters of the
    2 % unused and still solves in under a second (0.8 mm would fail the row).

Numbers to expect (2026-11-06, Elmer 26.2-devel rev 11c118114, MUMPS and Lua linked in, no MPI): |B(R)| 1.9904 mT (-0.48 %),
    |B(R/2)| 1.0019 mT (+0.19 %); the whole row 1.5 s (ElmerSolver's own total 0.6 s). The printouts of ElmerGrid and ElmerSolver are
    kept in the case folder (elmergrid.out.txt, elmersolver.out.txt) beside line.dat.

Self-check (the answers follow).
  1. Where would the field of the same wire be 0.5 mT? (Both answers.)
  2. The conductor's relative permeability is set to 1. If it were steel (500), what would change inside and what outside?
  3. Why is Az = 0 on the outer circle not an approximation of "zero field at infinity"?
Answers.
  1. At r = R/4 = 2.5 mm inside (B = 2.0 mT * r/R) and at r = 4 R = 40 mm outside (B = 2.0 mT * R/r).
  2. Inside, B would be 500 times larger (B = mu I r / (2 pi R^2)), and so would the energy in the wire; outside nothing changes:
     the field outside depends on the enclosed current only, and a 2D wire of permeability mu with the same current has the same
     H outside and the same B, because the air around it is the same.
  3. For one conductor the potential is Az = -(mu0 I / 2 pi) ln(r / r0) + const, and r0 is whatever sets the constant: choosing
     Az = 0 at r0 = 10 R fixes the constant and changes B by nothing. A second wire would break this, which is why the motor's
     model has its boundary far from the stator and does not rely on it.
"""
from __future__ import annotations

import math
import re
import time
from pathlib import Path

import femm_queue

from . import FAIL, PASS, Row, fresh, intact_path, need_tools, save_log

RADIUS = 0.010                  # m: the wire
OUTER = 0.100                   # m: the edge of the air round it (10 R)
CURRENT = 100.0                 # A
MU0 = 4e-7 * math.pi            # H/m
H_MIN, H_MAX = 0.0002, 0.008    # m: element size at the wire's surface and at the outer edge
TOLERANCE_PERCENT = 2.0
TIMEOUT_S = 180


def field_of_conductor(current: float, radius: float, r: float) -> float:
    """B [T] at distance r [m] from the axis of a long round wire of this radius [m] with this current [A] spread evenly over it."""
    if r <= radius:
        return MU0 * current * r / (2 * math.pi * radius ** 2)
    return MU0 * current / (2 * math.pi * r)


def current_density(current: float, radius: float) -> float:
    """Jz [A/m^2] of a current spread evenly over the wire's cross-section."""
    return current / (math.pi * radius ** 2)


def case_sif(current: float = CURRENT, radius: float = RADIUS, mesh: str = "case") -> str:
    """The ElmerSolver input file: the wire (body 1, carries the current) in air (body 2), Az = 0 on the outer circle (boundary 1)."""
    return f"""Header
  CHECK KEYWORDS Warn
  Mesh DB "." "{mesh}"
End

Simulation
  Max Output Level = 4
  Coordinate System = Cartesian 2D
  Simulation Type = Steady State
  Steady State Max Iterations = 1
End

Constants
  Permeability of Vacuum = {MU0!r}
End

Body 1
  Name = "Wire"
  Target Bodies(1) = 1
  Equation = 1
  Material = 1
  Body Force = 1
End

Body 2
  Name = "Air"
  Target Bodies(1) = 2
  Equation = 1
  Material = 2
End

Equation 1
  Active Solvers(3) = 1 2 3
End

Material 1
  Name = "Copper"
  Relative Permeability = 1.0
  Relative Permittivity = 1.0
End

Material 2
  Name = "Air"
  Relative Permeability = 1.0
  Relative Permittivity = 1.0
End

Body Force 1
  Name = "Current"
  Current Density = Real {current_density(current, radius)!r}
End

Boundary Condition 1
  Name = "OuterEdge"
  Target Boundaries(1) = 1
  Az = Real 0.0
End

Solver 1
  Equation = "MgDyn2D"
  Procedure = "MagnetoDynamics2D" "MagnetoDynamics2D"
  Variable = Az
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
End

Solver 2
  Equation = "MgDynPost"
  Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"
  Potential Variable = String "Az"
  Calculate Nodal Fields = Logical True
  Calculate Elemental Fields = Logical False
  Linear System Solver = Direct
  Linear System Direct Method = UMFPack
End

Solver 3
  Exec Solver = After All
  Equation = "SaveLine"
  Procedure = "SaveData" "SaveLine"
  Filename = "line.dat"
  Polyline Coordinates(2,3) = 0.0 0.0 0.0  {2 * radius!r} 0.0 0.0
  Polyline Divisions(1) = 40
End
"""


def write_mesh(path: Path) -> tuple[int, int]:
    """Mesh the wire in its air with Gmsh (the Python module) and write path in the MSH 2.2 format. Returns (nodes, triangles)."""
    import gmsh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("wire")
        occ = gmsh.model.occ
        wire = occ.addDisk(0, 0, 0, RADIUS, RADIUS)
        air = occ.addDisk(0, 0, 0, OUTER, OUTER)
        occ.fragment([(2, air)], [(2, wire)])            # the air disk is cut by the wire: two surfaces that share the wire's circle
        occ.synchronize()
        area = {tag: occ.getMass(2, tag) for _, tag in gmsh.model.getEntities(2)}
        wire_surface, air_surface = min(area, key=area.get), max(area, key=area.get)

        def reach(curve: int) -> float:              # how far from the axis the curve goes (its bounding box is a hair larger than it)
            x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(1, curve)
            return max(abs(x0), abs(x1), abs(y0), abs(y1))

        curves = [tag for _, tag in gmsh.model.getEntities(1)]
        edge = [c for c in curves if abs(reach(c) - OUTER) < 1e-6]
        skin = [c for c in curves if abs(reach(c) - RADIUS) < 1e-6]
        if len(edge) != 1 or len(skin) != 1 or len(curves) != 2:
            raise RuntimeError(f"expected the wire's circle and the outer circle, found curves {curves} (edge {edge}, wire {skin})")
        gmsh.model.addPhysicalGroup(2, [wire_surface], 1, name="wire")
        gmsh.model.addPhysicalGroup(2, [air_surface], 2, name="air")
        gmsh.model.addPhysicalGroup(1, edge, 1, name="outer edge")
        distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(distance, "CurvesList", skin)
        gmsh.model.mesh.field.setNumber(distance, "Sampling", 400)
        size = gmsh.model.mesh.field.add("Threshold")                 # H_MIN within 4 mm of the wire's surface, growing to H_MAX at 60 mm
        for key, value in (("InField", distance), ("SizeMin", H_MIN), ("SizeMax", H_MAX), ("DistMin", 0.004), ("DistMax", 0.06)):
            gmsh.model.mesh.field.setNumber(size, key, value)
        gmsh.model.mesh.field.setAsBackgroundMesh(size)
        for key, value in (("Mesh.MeshSizeExtendFromBoundary", 0), ("Mesh.MeshSizeFromPoints", 0), ("Mesh.MeshSizeFromCurvature", 0),
                           ("Mesh.MshFileVersion", 2.2)):
            gmsh.option.setNumber(key, value)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(path))
        nodes = len(gmsh.model.mesh.getNodes()[0])
        triangles = sum(len(tags) for tags in gmsh.model.mesh.getElements(2)[1])
    finally:
        gmsh.finalize()
    return nodes, triangles


def read_line(folder: Path):
    """(x, |B|) along the saved line: the columns are found from the names file, not assumed."""
    import numpy as np

    columns = {}
    for text in (folder / "line.dat.names").read_text(encoding="utf-8", errors="replace").splitlines():
        found = re.match(r"\s*(\d+):\s*(.+?)\s*$", text)
        if found:
            columns[found.group(2).lower()] = int(found.group(1)) - 1
    needed = ("coordinate 1", "magnetic flux density 1", "magnetic flux density 2")
    missing = [name for name in needed if name not in columns]
    if missing:
        raise RuntimeError(f"line.dat.names does not list {', '.join(missing)} (it lists: {', '.join(columns)})")
    data = np.loadtxt(folder / "line.dat")
    return data[:, columns["coordinate 1"]], np.hypot(data[:, columns["magnetic flux density 1"]], data[:, columns["magnetic flux density 2"]])


def run(folder: Path):
    solver, grid = need_tools("elmer", "elmergrid")
    fresh(folder)
    case = f"MagnetoDynamics2D: round copper wire, R = {RADIUS * 1e3:g} mm, {CURRENT:g} A, in air (meshed with Gmsh)"
    expected = (f"|B| within {TOLERANCE_PERCENT:g} % of {field_of_conductor(CURRENT, RADIUS, RADIUS) * 1e3:.1f} mT at r = R "
                f"and of {field_of_conductor(CURRENT, RADIUS, RADIUS / 2) * 1e3:.1f} mT at r = R/2")
    started = time.perf_counter()
    with intact_path():                             # gmsh.initialize() cuts the process's PATH short (see checks.intact_path)
        nodes, triangles = write_mesh(folder / "case.msh")
    (folder / "case.sif").write_text(case_sif(), encoding="ascii")
    converted = femm_queue.run_solver([grid, "14", "2", "case.msh"], label="check_tools Elmer ElmerGrid", cwd=folder, timeout_s=TIMEOUT_S, threads=1)
    save_log(folder, "elmergrid", converted)
    if converted.timed_out or converted.returncode != 0 or not (folder / "case" / "mesh.header").is_file():
        yield Row("Elmer", case, expected, f"ElmerGrid failed (exit code {converted.returncode}); see elmergrid.out.txt", FAIL,
                  round(time.perf_counter() - started, 2))
        return
    solved = femm_queue.run_solver([solver, "case.sif"], label="check_tools Elmer ElmerSolver", cwd=folder, timeout_s=TIMEOUT_S, threads=2)
    save_log(folder, "elmersolver", solved)
    seconds = round(time.perf_counter() - started, 2)
    if solved.timed_out or solved.returncode != 0 or "ALL DONE" not in solved.stdout or not (folder / "line.dat").is_file():
        why = f"timed out after {TIMEOUT_S} s" if solved.timed_out else f"exit code {solved.returncode}"
        tail = " | ".join(line for line in solved.stdout.splitlines() if "ERROR" in line or "FATAL" in line.upper())[:200]
        yield Row("Elmer", case, expected, f"ElmerSolver failed: {why}. {tail} (see elmersolver.out.txt)", FAIL, seconds)
        return
    x, b = read_line(folder)
    errors, shown = [], []
    for label, r in (("R", RADIUS), ("R/2", RADIUS / 2)):
        i = int(abs(x - r).argmin())
        if abs(x[i] - r) > 1e-6:
            raise RuntimeError(f"line.dat has no point at r = {label} (nearest x = {x[i]:.6f} m)")
        exact = field_of_conductor(CURRENT, RADIUS, r)
        errors.append(100 * (b[i] / exact - 1))
        shown.append(f"|B({label})| = {b[i] * 1e3:.4f} mT ({errors[-1]:+.2f} %)")
    banner = re.search(r"Version:\s*(\S+)\s*\(Rev:\s*(\w+)", solved.stdout)
    version = f"Elmer {banner.group(1)} rev {banner.group(2)}" if banner else "Elmer ?"
    measured = f"{version}, {nodes:,} nodes: " + ", ".join(shown)
    yield Row("Elmer", case, expected, measured, PASS if all(abs(e) < TOLERANCE_PERCENT for e in errors) else FAIL, seconds)
