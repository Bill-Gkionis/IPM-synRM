"""B3: Gmsh -> Elmer (ElmerGrid, ElmerSolver, MagnetoDynamics2D) on the motor's own mesh: the no-load field, held against FEMM's.

What it shows, and why
    Elmer is the open solver the pipeline can use where FEMM cannot go (a field that is not 2D, a coupling to another physics). Task 4 proved the
    chain Gmsh -> ElmerGrid -> ElmerSolver on a wire, where Ampere's law gives the answer. B3 sends the MOTOR through it: B2's mesh of the band-free
    sector, N45 magnets, M-19 steel on its nonlinear B-H curve, the sector's periodic edges, no current. At no load the field is the magnets'
    alone, and FEMM (B1's model, the same circle in the air gap) gives the number to hold Elmer against. Elmer is sound here if it reproduces FEMM's
    air-gap flux density (its fundamental within 3 %, its harmonic pattern alike) and the flux the magnets put through the winding (psi_PM within 3 %),
    and converges.

Method (the files and the protocol)
    1. The mesh. B2's path 1: motorgeo.to_gmsh on motorgeo.build(Params(barrier="none"), gap="sliding"), size_max 1.5 mm, written as case.msh (MSH 2.2).
       The sector's two periodic edges (0 and 90 degrees) carry the same nodes at the same radii, and the two sides of the air gap (rotor side
       and stator side, meshed apart so that the rotor can turn later) carry the same nodes at the same angles; interface_facts() checks both before
       Elmer is asked, because ElmerSolver refuses a periodic pair that does not match one to one (a finer mesh of 0.7 mm does not: Elmer says "This is
       not conforming!"; that size is not used).
    2. ElmerGrid 14 2 case.msh -scale 0.001 0.001 0.001 -merge 1e-9 (through run_solver). 14 -> 2 is "Gmsh to ElmerSolver format"; -scale turns the mesh
       from millimetres into metres (Elmer's equations are in SI); -merge fuses the nodes that coincide, which closes the sliding air-gap interface at the
       fixed rotor position (269 nodes) and leaves ElmerGrid saying "a single piece in the mesh, so the mesh is conforming". The physical group
       names and tags of the Gmsh file are the bodies and boundaries of the case; the input file is written from them (sif_text).
    3. ElmerSolver case.sif (through run_solver, 4 threads). Solver 1 MagnetoDynamics2D solves for the vector potential Az; solver 2
       MagnetoDynamicsCalcFields gives B = curl(Az) at the nodes; solver 3 writes the result as ASCII VTU (case\\case_t0001.vtu), which this module reads.
    4. The numbers. B_r on the mid-gap circle (r = 45.4 mm, the arc where the two sides of the gap meet; FEMM's air-gap element is centred on the
       same circle) is sampled every 0.25 degrees from the nodes on that circle. The flux linkage of each phase is computed from Az the way FEMM computes
       a circuit's: for every coil side, turns x (the mean of Az over the side) x the stack length, summed over the phase's sides and the 4 sectors.
       psi_PM is the length of the (d, q) flux vector at the rotor position of the model (an amplitude-invariant Park transform of the three phases).
    5. FEMM's side: B1's model (barrier="none", fine mesh) solved at the same drawn rotor position (B1's solve_chunk, on a femm_pool of 8 workers),
       its 360 samples of B_r on the same circle, and psi_PM as B1 defines it (the fundamental of phase A's flux linkage over 120 rotor positions).
       B3 builds no FEMM model of its own: it calls b1_femm.reference(), which runs B1's own model builder and solver. That is the point: B3 holds Elmer
       against exactly the model that B1 validated against the README's row, so what B3 says of Elmer is said of that model.

Two things the physics needed that the brief did not say (found on 2026-11-06)
    * The magnet. FEMM's N45 (mu_r 1.05, coercivity Hc = 1,030,272 A/m) is B = mu0 mu_r (H + Hc). Elmer's 2D equation is curl(nu curl A) = curl M with
      nu = 1/(mu0 mu_r) (Models Manual, 22.2 Theory, equation 22.1), that is B = mu0 mu_r (H + M): the magnetisation M of Elmer's `Magnetization 1/2` is FEMM's
      Hc, not Hc times mu_r. Check: the 10 x 2 mm N45 magnet in a 60 mm box of air (learn/03, check_queue.py) gives By = 0.1582 T at its centre in FEMM,
      0.1588 T in Elmer.
    * The lamination. FEMM's "M-19 Steel" is the 47 points of its B-H curve AND a lamination fill factor of 0.98 (motorgeo.M19_STEEL_FEMM42["LamFill"]),
      and FEMM applies it in magnetostatics as B_eff(H) = f B(H) + (1 - f) mu0 H: FEMM with fill 1.0 and the curve changed by that formula gives, to
      every printed digit, FEMM with fill 0.98 (psi_A = -15.6797 mWb), while FEMM with fill 1.0 and the curve as it is gives psi_A = -15.308 mWb: 2.4 % less
      (the likely reason, not tested: a stack that is 2 % weaker saturates earlier, its bridges leak less, and more of the magnets' flux reaches the
      winding). The probe was made once, outside the bench: a femm_pool task that sets motorgeo.M19_STEEL_FEMM42["LamFill"] and ["BH"] in its own process
      before to_femm (the bench changes no file of motorgeo). An Elmer model that takes the 47 points alone is therefore another steel. B3 runs both:
      `laminated` (the points and the fill, the same steel as FEMM's: the verdict of the case; results.json says "verdict_run": "laminated") and `plain`
      (the 47 points exactly as the brief says, no fill: reported with its own pass or fail and the difference, but not gating: every check about it, and
      the one that compares the two runs, has "gating": false and says so in its expected text, so a plain miss cannot change the exit status).

The keywords, and where each comes from (Elmer Models Manual 26.1 of 20 January 2026, ElmerSolver Manual of the same date, and elmerfem's regression tests;
the sections are cited by number and by the heading the PDF prints: Models Manual "22.2 Theory", "22.2.1 Boundary Conditions", "22.3 Keywords", "63.2 Keywords" (of
chapter 63, "Result output in different formats"), ElmerSolver Manual "9.3 Conforming conditions" (of chapter 9, "Enforcing shared conditions between boundaries"))
    Material:  H-B Curve = Variable "dummy" / Real Cubic / <B H pairs> / End      Models Manual 22.3 Keywords ("H-B Curve ... must be given as a cubic spline");
               syntax: elmerfem fem/tests/mgdyn_bh/case.sif (h-b curve = Variable "dummy" Real ... End, a table of B and H) and ModelDescription.F90
               (the table reader: "Real Cubic", "Real Cubic Monotone", "Real"). The solver evaluates H at |B|, so the table is B first, H second.
    Material:  Relative Permeability, Magnetization 1, Magnetization 2           Models Manual 22.3; elmerfem fem/tests/mgdyn2D_pm/ptest.sif (a magnet with
               Relative Permeability = 1.05 and Magnetization 1/2 = Real ...). One material per magnet: each has its own direction.
    Solver 1:  Nonlinear System Max Iterations, Nonlinear System Convergence Tolerance      Models Manual 22.3: "If a nonlinear H-B curve is given, then
               Newton's linearization is applied after the first iteration" (MagnetoDynamics2D.F90 switches to Newton-Raphson at iteration 2).
    Solver 1:  Apply Conforming BCs = Logical True                                ElmerSolver Manual 9.3 Conforming conditions
    BC:        Conforming BC = Integer <the paired BC>, Radial Projector = Logical True    ElmerSolver Manual 9.3 Conforming conditions, quoted whole:
               "a direct node-to-node mapping, or edge-to-edge mapping on the two boundaries. In practice we use a simple search and identify node pairs by their
               coordinates" (the radial projector matches by radius); elmerfem
               fem/tests/mgdyn2D_em_conforming/konesec.sif (a PM machine sector, its edges eliminated at assembly); MortarUtils.F90 PeriodicPermutation:
               Radial Projector is periodic, Anti Radial Projector is antiperiodic. Two poles per sector: periodic.
    BC:        Az = Real 0.0 on the shaft and the stator's outer edge                Models Manual 22.2.1 Boundary Conditions (Dirichlet condition for A3), and Task 4's elmer_check
    Solver 3:  ResultOutputSolve with Vtu Format, Ascii Output, Save Geometry Ids   Models Manual 63.2 Keywords
    Constants: Permeability of Vacuum = 4 pi e-7                                    Models Manual 22.3 (SI units; hence ElmerGrid -scale 0.001)
    (CHECK KEYWORDS Warn in the header makes ElmerSolver name any keyword it does not know; the run prints none.)

Run (from the Python to FEMM folder; it takes the solver slot for about a minute: FEMM's reference first, then ElmerGrid and ElmerSolver)
    uv run python bench/run_bench.py --only B3
Writes out\\bench\\B3\\: results.json, figure.png (Elmer's field), profile.png (B_r on the mid-gap circle and its harmonics, Elmer against FEMM), run.log,
runs.json, and one folder per Elmer run (elmer_laminated\\, elmer_plain\\, elmer_s1\\, elmer_s0.8\\, elmer_s0.5\\) with case.sif, the converted mesh, the solver's
own printout and the result file.

Numbers to expect (2026-11-06, Elmer 26.2-devel rev 11c118114, FEMM 5.0.20261101; the whole case is about 30 s, of which FEMM's reference 12 s)
    mesh     B2's: 8,819 nodes in the saved file (that is what ElmerGrid reads; motorgeo.to_gmsh returns 8,906, which also counts the 87 arc centres that no element
             uses and the file does not hold); ElmerGrid merges the 269 nodes of the gap's two sides and leaves Elmer 8,550 nodes, 16,914 triangles, one piece; 0.06 s
    solver   ElmerSolver 1.4 s on 4 threads: 11 nonlinear iterations (relative change 2.1e-7), no warning; Az on the two periodic edges equal to the last digit
    FEMM     B_r fundamental 0.32248 T on the mid-gap circle; psi_PM 16.629 mWb (over 120 positions), |psi_dq| 16.571 mWb at the drawn position
    laminated   (the default mesh) B_r fundamental 0.32141 T (-0.33 %); psi_PM 16.521 mWb (-0.65 %; -0.30 % against FEMM at the same position); the largest harmonics
             (order: FEMM / Elmer, % of the fundamental) 5: 24.4 / 22.8, 3: 20.4 / 20.2, 9: 13.7 / 12.1; shapes 7.5 % rms
    plain    B_r fundamental 0.31386 T (-2.67 %); psi_PM 16.133 mWb (-2.98 %; -2.64 % at the same position): the 47 points alone, without FEMM's lamination
    meshes   laminated, B_r and psi_PM against FEMM: size_max 1.0 -0.84 % and -1.11 % (4.5 % rms), 0.8 -1.04 % and -1.35 % (1.7 %), 0.5 -1.08 % and -1.48 % (1.2 %)

Self-check (the answers follow).
  1. Why is the sector's periodic condition "periodic" here and "anti-periodic" for a machine with one pole in its sector?
  2. A first-order Elmer triangle has a constant B. How can the field along the mid-gap circle be a smooth curve?
  3. Why is psi_PM taken as the length of the (d, q) vector and not as the flux linkage of phase A?
  4. The plain run is 2.7 % below FEMM and the laminated one 0.3 % (on the default mesh). Which is the better solver, and why does the question not make sense?
Answers.
  1. A sector of one pole pitch ends where the field has the opposite sign (north pole on one edge, south on the other): A(edge 2) = - A(edge 1). The
     90-degree sector holds two poles, a whole pair: the pattern repeats exactly, A(edge 2) = A(edge 1).
  2. It is not: B_r comes from the nodal fields, which MagnetoDynamicsCalcFields makes by averaging the triangles around each node, and it is
     interpolated linearly between the nodes of the circle (0.33 degrees apart). The waveform therefore has small steps and spikes where FEMM's
     air-gap element, which carries the gap field as a smooth series, has none: that is the 7.5 % rms between the two shapes on B2's mesh, and it
     falls to 1.2 % when the mesh is finer (0.5 mm).
  3. At one rotor position phase A's flux linkage is psi_PM cos(angle): it depends on where the rotor happens to be. The length of the (d, q) vector does
     not (apart from the small 6th-harmonic ripple), so it is the number that can be held against FEMM's psi_PM whatever the position.
  4. Neither: the plain run solves a different steel. FEMM's M-19 has a 2 % lamination factor and Elmer's plain run has none; once Elmer is given
     the same steel, the two tools agree to a third of a percent on the default mesh (size_max 1.5 mm: B_r fundamental -0.33 %, psi_PM -0.30 % at the same
     position). That is not the whole story, and it is not a figure that improves with the mesh: on finer meshes the difference grows, to -1.08 % in the B_r
     fundamental and -1.48 % in psi_PM at size_max 0.5 mm (still inside the 3 % of the brief). One contributor is known: the N45 area is 0.172 % short at
     every mesh size (the chords of the magnet fillets, B2). The rest is not isolated. The two programs share the 47 points of the B-H curve but not the
     interpolation between them, and the saturated bridges set the leakage, so a small difference in how each treats the knee of the curve would show there:
     candidates, not shown. (The chords on the rotor's surface that B2 reports shrink with size_max while the difference grows, so
     at most they offset it on the default mesh.) Comparing two solvers means giving them the same problem first.
"""
from __future__ import annotations

import math
import re
import shutil
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

import motorgeo as mg
import paths
from checks import intact_path

from . import OWN_LIMIT, Log, check, compare, fig_field, params, run_program
from . import b1_femm, b2_gmsh

CASE = "B3"
TITLE = "Gmsh -> Elmer (ElmerGrid, ElmerSolver, MagnetoDynamics2D): the no-load field of the band-free sector, against FEMM's"
VERDICT_RUN = "laminated"               # the Elmer run that carries the case's verdict (the ledger's ruling); the plain-curve run is reported beside it

HC = 1_030_272.0                        # A/m: N45's coercivity as FEMM has it (motorgeo.to_femm: mi_addmaterial("N45", 1.05, 1.05, 1030272))
MU_R_MAGNET = 1.05
MU0 = 4e-7 * math.pi
FUNDAMENTAL_TOLERANCE_PERCENT = 3.0     # the brief: the air-gap flux density's fundamental within 3 % of FEMM's
PSI_TOLERANCE_PERCENT = 3.0             # the brief: psi_PM within 3 % of FEMM's
# "The harmonic pattern alike" has no number in the brief. The three limits below are the bench's own, chosen after a first exploratory run on the default mesh
# (differences of up to 2.2 points per harmonic, 7.7 % rms) and before the case's own run; the finer-mesh rows show the shapes converging (4.5 %, 1.7 % rms).
PATTERN_POINTS = 3.0                    # every harmonic 2 ... 20, as a percentage of its own fundamental, within this many points
PATTERN_RMS_PERCENT = 10.0              # ... the two waveforms, each divided by its fundamental, within this rms (percent of the fundamental)
PATTERN_TOP = 3                         # ... and the three largest harmonics are the same orders
PATTERN_ORDERS = range(2, 21)
MAX_ITERATIONS = 60
CONVERGENCE_TOLERANCE = 1e-6
MESH_SCALE = "0.001"                    # ElmerGrid -scale: the mesh is in millimetres, Elmer works in metres
MERGE_DISTANCE = "1e-9"                 # ElmerGrid -merge: far below the smallest edge of the mesh (0.07 mm), far above the coincidence of the two sides (1e-13)
GRID_TIMEOUT_S = 120
SOLVER_TIMEOUT_S = 600
SOLVER_THREADS = 4
GRID_THREADS = 1
SENSITIVITY_SIZES = (1.0, 0.8, 0.5)     # finer meshes of path 1 for the laminated run (a size whose periodic edges do not match, 0.7 for one, is skipped and said)
NODE_TOLERANCE = 1e-9                   # nodes of two sides that "coincide": radii in mm, angles in degrees


# ============================================================================================== the input file
def physical_names(path) -> dict[str, tuple[int, int]]:
    """{name: (dimension, tag)} of the physical groups of an MSH 2.2 file (its $PhysicalNames section). The tags are Elmer's body and boundary numbers."""
    text = Path(path).read_text(encoding="ascii", errors="replace")
    block = text.split("$PhysicalNames")[1].split("$EndPhysicalNames")[0].strip().splitlines()
    out: dict[str, tuple[int, int]] = {}
    for line in block[1:]:
        dim, tag, name = line.split(None, 2)
        out[name.strip().strip('"')] = (int(dim), int(tag))
    return out


def steel_table(fill: float) -> list[tuple[float, float]]:
    """The (B, H) pairs of FEMM 4.2's "M-19 Steel" as a laminated stack of fill factor `fill`: B_eff(H) = fill B(H) + (1 - fill) mu0 H, the way FEMM
    applies its lamination fill in magnetostatics (the module docstring: FEMM with fill 1.0 and the curve changed this way reproduces FEMM with fill 0.98).
    With fill = 1 these are the 47 points of motorgeo.M19_STEEL_FEMM42 exactly."""
    return [(fill * float(b) + (1.0 - fill) * MU0 * float(h), float(h)) for b, h in mg.M19_STEEL_FEMM42["BH"]]


def sif_text(names: dict, regions, *, fill: float, mesh: str = "case", max_iterations: int = MAX_ITERATIONS, tolerance: float = CONVERGENCE_TOLERANCE) -> str:
    """The ElmerSolver input file of the no-load motor: one body for every region (named, targeting the physical group of its name), the steel on its H-B
    curve (laminated with `fill`), the magnets with their coercivity as magnetisation along their own direction, air and coils at mu_r = 1 (no current), A = 0 on
    the shaft and the outer edge, the two periodic edges as a conforming pair, and the solvers. ASCII, one keyword to a line, nothing but what the
    module docstring's table names."""
    magnets = [r for r in regions if r.name.startswith("magnet_")]
    material_of: dict[str, int] = {r.name: 3 + i for i, r in enumerate(magnets)}
    bodies = sorted(((tag, name) for name, (dim, tag) in names.items() if dim == 2))
    region = {r.name: r for r in regions}
    out = ["Header", "  CHECK KEYWORDS Warn", f'  Mesh DB "." "{mesh}"', "End", "",
           "Simulation", "  Max Output Level = 5", "  Coordinate System = Cartesian 2D", "  Simulation Type = Steady State",
           "  Steady State Max Iterations = 1", "End", "",
           "Constants", f"  Permeability of Vacuum = {MU0!r}", "End", ""]
    for k, (tag, name) in enumerate(bodies, start=1):
        r = region[name]
        material = 1 if r.material == "M-19 Steel" else material_of[name] if r.material == "N45" else 2
        out += [f"Body {k}", f'  Name = "{name}"', f"  Target Bodies(1) = {tag}", "  Equation = 1", f"  Material = {material}", "End", ""]
    out += ["Equation 1", "  Active Solvers(3) = 1 2 3", "End", ""]
    rows = [f"      {b!r} {h!r}" for b, h in steel_table(fill)]
    permittivity = "  Relative Permittivity = 1.0"            # not used by a magnetostatic solve; MagnetoDynamicsCalcFields warns if a material has none
    out += ["Material 1", '  Name = "M-19 Steel"', '  H-B Curve = Variable "dummy"', "    Real Cubic", *rows, "    End", permittivity, "End", ""]
    out += ["Material 2", '  Name = "Air and coil"', "  Relative Permeability = 1.0", permittivity, "End", ""]
    for r in magnets:
        angle = math.radians(r.magdir)
        out += [f"Material {material_of[r.name]}", f'  Name = "{r.name}"', f"  Relative Permeability = {MU_R_MAGNET!r}",
                f"  Magnetization 1 = Real {HC * math.cos(angle)!r}", f"  Magnetization 2 = Real {HC * math.sin(angle)!r}", permittivity, "End", ""]
    out += ["Boundary Condition 1", '  Name = "A0"', f"  Target Boundaries(1) = {names['A0'][1]}", "  Az = Real 0.0", "End", "",
            "Boundary Condition 2", '  Name = "edge_0"', f"  Target Boundaries(1) = {names['edge_0'][1]}", "  Conforming BC = Integer 3",
            "  Radial Projector = Logical True", "End", "",
            "Boundary Condition 3", '  Name = "edge_sector"', f"  Target Boundaries(1) = {names['edge_sector'][1]}", "End", ""]
    out += ["Solver 1", '  Equation = "MgDyn2D"', '  Procedure = "MagnetoDynamics2D" "MagnetoDynamics2D"', "  Variable = Az",
            f"  Nonlinear System Max Iterations = {max_iterations}", f"  Nonlinear System Convergence Tolerance = {tolerance!r}",
            "  Apply Conforming BCs = Logical True", "  Linear System Solver = Direct", "  Linear System Direct Method = UMFPack", "End", "",
            "Solver 2", '  Equation = "MgDynPost"', '  Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"', '  Potential Variable = String "Az"',
            "  Calculate Nodal Fields = Logical True", "  Calculate Elemental Fields = Logical False",
            "  Linear System Solver = Direct", "  Linear System Direct Method = UMFPack", "End", "",
            "Solver 3", "  Exec Solver = After All", '  Equation = "ResultOutput"', '  Procedure = "ResultOutputSolve" "ResultOutputSolver"',
            '  Output File Name = "case"', "  Vtu Format = Logical True", "  Ascii Output = Logical True", "  Save Geometry Ids = Logical True", "End", ""]
    return "\n".join(out) + "\n"


# ============================================================================================== what the programs said
def parse_grid_log(text: str) -> dict:
    """What ElmerGrid says about the merge and the result: nodes read, nodes after the merge and how many were merged, whether the mesh is one piece
    (conforming), the pieces, the coordinates and elements saved. A number it does not say is None."""
    merged = re.search(r"(\d+) original nodes merged to (\d+) new nodes", text)
    single = re.search(r"There is a single piece in the mesh, so the mesh is conforming", text)
    several = re.search(r"There are (\d+) pieces", text)
    coordinates, elements = re.search(r"Saving (\d+) coordinates", text), re.search(r"Saving (\d+) element topologies", text)
    return {"nodes_read": int(merged.group(1)) if merged else None, "nodes_merged_to": int(merged.group(2)) if merged else None,
            "nodes_merged": int(merged.group(1)) - int(merged.group(2)) if merged else None,
            "conforming": True if single else False if several else None, "pieces": 1 if single else int(several.group(1)) if several else None,
            "coordinates_saved": int(coordinates.group(1)) if coordinates else None, "elements_saved": int(elements.group(1)) if elements else None}


def parse_solver_log(text: str) -> dict:
    """What ElmerSolver says: the version, the threads, the nonlinear iterations of MagnetoDynamics2D (one `ComputeChange: NS (ITER=n)` line each, with the
    relative change of the solution), whether it finished (ALL DONE), the field energy, the time, the warnings and errors."""
    changes = [float(m.group(1)) for m in re.finditer(r"ComputeChange: NS \(ITER=\d+\) \(NRM,RELC\):\s*\(\s*\S+\s+(\S+)\s*\)\s*::\s*mgdyn2d", text)]
    version = re.search(r"Version:\s*(\S+)\s*\(Rev:\s*(\w+)", text)
    threads = re.search(r"Running in parallel with (\d+) threads", text)
    energy = re.search(r"Field Energy:\s*(\S+)", text)
    seconds = re.search(r"SOLVER TOTAL TIME\(CPU,REAL\):\s*(\S+)\s+(\S+)", text)
    return {"iterations": len(changes), "changes": changes, "last_change": changes[-1] if changes else None, "finished": "ALL DONE" in text,
            "threads": int(threads.group(1)) if threads else None, "version": version.group(1) if version else None,
            "revision": version.group(2) if version else None, "energy_J": float(energy.group(1)) if energy else None,
            "solver_seconds_real": float(seconds.group(2)) if seconds else None,
            "warnings": [m.group(1).strip() for m in re.finditer(r"^WARNING::\s*(.+)$", text, re.M)],
            "errors": [m.group(1).strip() for m in re.finditer(r"^ERROR::\s*(.+)$", text, re.M)]}


# ============================================================================================== the result file
def read_vtu(path) -> dict:
    """Elmer's ASCII VTU as arrays: {"points_mm" (n, 2), "az" (n,), "b" (n, 2) (Bx, By), "triangles" (m, 3), "body_ids" (m,), "lines" (k, 2), "line_ids" (k,)}.
    The file is in metres (the mesh was scaled), so the points come back in millimetres, like everything else here."""
    piece = ET.parse(str(path)).getroot().find(".//Piece")

    def array(element):
        kind = element.get("type", "")
        data = np.array(element.text.split(), dtype=float if kind.startswith("Float") else np.int64)
        width = int(element.get("NumberOfComponents", "1"))
        return data.reshape(-1, width) if width > 1 else data

    point = {da.get("Name"): array(da) for da in piece.find("PointData")}
    cell = {da.get("Name"): array(da) for da in piece.find("CellData")}
    points = array(piece.find("Points")[0])
    parts = {da.get("Name"): array(da) for da in piece.find("Cells")}
    types, offsets, connectivity = parts["types"], parts["offsets"], parts["connectivity"]
    starts = np.concatenate([[0], offsets[:-1]])
    triangle_cells, line_cells = np.where(types == 5)[0], np.where(types == 3)[0]
    cut = lambda cells: np.array([connectivity[starts[i]:offsets[i]] for i in cells], dtype=np.int64).reshape(len(cells), -1)
    return {"points_mm": points[:, :2] * 1e3, "az": point["az"], "b": point["magnetic flux density"][:, :2], "triangles": cut(triangle_cells),
            "body_ids": cell["GeometryIds"][triangle_cells], "lines": cut(line_cells), "line_ids": cell["GeometryIds"][line_cells]}


def _triangle_areas(points: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    a, b, c = points[triangles[:, 0]], points[triangles[:, 1]], points[triangles[:, 2]]
    return 0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))


def flux_linkage(vtu: dict, names: dict, regions, stack_mm: float, sectors: int) -> dict[str, float]:
    """Flux linkage [Wb] of each phase of the whole machine from the vector potential, as FEMM computes a circuit's: for every coil side of the phase,
    turns x (the integral of Az over the side / its area) x the stack length, summed over the sides and multiplied by the number of sectors."""
    area = _triangle_areas(vtu["points_mm"], vtu["triangles"])
    mean_a = vtu["az"][vtu["triangles"]].mean(axis=1)
    psi = {"A": 0.0, "B": 0.0, "C": 0.0}
    for r in regions:
        if r.material == "Coil" and r.circuit in psi:
            side = vtu["body_ids"] == names[r.name][1]
            psi[r.circuit] += r.turns * float((mean_a[side] * area[side]).sum() / area[side].sum())
    return {phase: value * stack_mm * 1e-3 * sectors for phase, value in psi.items()}


def park(psi_abc, theta_e: float) -> tuple[float, float]:
    """(psi_d, psi_q) of three phase flux linkages at the electrical angle theta_e [degrees]: motorgeo.FemmModel's amplitude-invariant transform."""
    angles = [math.radians(theta_e + s) for s in (0, -120, 120)]
    return (2 / 3 * sum(v * math.cos(t) for v, t in zip(psi_abc, angles)), -2 / 3 * sum(v * math.sin(t) for v, t in zip(psi_abc, angles)))


def gap_profile(vtu: dict, radius_mm: float, angles_deg) -> tuple[np.ndarray, np.ndarray]:
    """(B_r, B_t) [T] on the circle of radius `radius_mm` at `angles_deg`: the nodal field of the nodes ON that circle, interpolated linearly in the angle.
    ValueError if the mesh has no node on it."""
    r = np.hypot(vtu["points_mm"][:, 0], vtu["points_mm"][:, 1])
    on = np.where(np.abs(r - radius_mm) < 1e-6)[0]
    if len(on) < 2:
        raise ValueError(f"the mesh has {len(on)} node(s) on the circle of radius {radius_mm} mm: no profile can be sampled")
    angle = np.degrees(np.arctan2(vtu["points_mm"][on, 1], vtu["points_mm"][on, 0]))
    order = np.argsort(angle)
    angle, bx, by = angle[order], vtu["b"][on, 0][order], vtu["b"][on, 1][order]
    br, bt = bx * np.cos(np.radians(angle)) + by * np.sin(np.radians(angle)), -bx * np.sin(np.radians(angle)) + by * np.cos(np.radians(angle))
    return np.interp(angles_deg, angle, br), np.interp(angles_deg, angle, bt)


# ============================================================================================== comparing the two tools
def compare_waveforms(femm_br, other_br) -> dict:
    """Two samples of B_r over one electrical period (the 90-degree sector), compared the way B3 judges them.

    fundamental_percent   the other's fundamental against FEMM's, in percent (the number of the 3 % test)
    harmonics             for each order 2 ... 20: each waveform's amplitude as a percentage of ITS OWN fundamental, and the difference in points
    rms_percent_of_fundamental   the two waveforms, each divided by its own fundamental, their rms difference as a percentage of that fundamental
    alike                 "the harmonic pattern alike", set before the run: every harmonic within PATTERN_POINTS points, the rms within PATTERN_RMS_PERCENT, and
                          the PATTERN_TOP largest harmonics the same orders. It judges the shape, so a difference in amplitude alone does not spoil it."""
    femm, other = np.asarray(femm_br, float), np.asarray(other_br, float)
    hf, ho = b1_femm.harmonic_amplitudes(femm), b1_femm.harmonic_amplitudes(other)
    rows = [{"order": n, "femm_percent": 100 * hf[n] / hf[1], "other_percent": 100 * ho[n] / ho[1],
             "difference_points": 100 * ho[n] / ho[1] - 100 * hf[n] / hf[1]} for n in PATTERN_ORDERS]
    top = lambda h: [n for n in sorted(PATTERN_ORDERS, key=lambda n: -h[n])[:PATTERN_TOP]]
    rms = 100 * float(np.sqrt(np.mean((other / ho[1] - femm / hf[1]) ** 2)))
    alike = (all(abs(row["difference_points"]) <= PATTERN_POINTS for row in rows) and rms <= PATTERN_RMS_PERCENT and set(top(hf)) == set(top(ho)))
    return {"fundamental_femm": float(hf[1]), "fundamental_other": float(ho[1]), "fundamental_ratio": float(ho[1] / hf[1]),
            "fundamental_percent": 100 * float(ho[1] / hf[1] - 1), "harmonics": rows, "rms_percent_of_fundamental": rms,
            "top_orders_femm": top(hf), "top_orders_other": top(ho), "alike": bool(alike)}


def nodes_match(first, second, tolerance: float) -> dict:
    """Do two sets of positions (radii of the two periodic edges, angles of the two sides of the gap) coincide one to one? {"match", "count": (n1, n2),
    "max_difference"}: the positions are sorted, and each is compared with its partner (max_difference is None when the counts differ)."""
    a, b = np.sort(np.asarray(first, float)), np.sort(np.asarray(second, float))
    if len(a) != len(b):
        return {"match": False, "count": (len(a), len(b)), "max_difference": None}
    worst = float(np.abs(a - b).max()) if len(a) else 0.0
    return {"match": worst <= tolerance, "count": (len(a), len(b)), "max_difference": worst}


def interface_facts(msh_path) -> dict:
    """The mesh facts Elmer needs, read from the MSH file in a Gmsh session of its own: {"periodic_edges": the radii of the nodes of edge_0 against edge_sector,
    "gap_interface": the angles of the nodes of sliding_rotor against sliding_stator}, each as nodes_match() says. The nodes of a group are those of its line
    elements (the corner nodes included)."""
    import gmsh
    with intact_path():
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.open(str(msh_path))
            tags, coords, _ = gmsh.model.mesh.getNodes()
            xy = {int(t): c[:2] for t, c in zip(tags, np.asarray(coords, float).reshape(-1, 3))}
            group = {gmsh.model.getPhysicalName(d, t): (d, t) for d, t in gmsh.model.getPhysicalGroups(1)}

            def nodes(name):
                found = set()
                for entity in gmsh.model.getEntitiesForPhysicalGroup(*group[name]):
                    for _, _, node_tags in zip(*gmsh.model.mesh.getElements(1, entity)):
                        found.update(int(n) for n in node_tags)
                return np.array([xy[n] for n in sorted(found)])

            radius = lambda name: np.hypot(*nodes(name).T)
            angle = lambda name: np.degrees(np.arctan2(nodes(name)[:, 1], nodes(name)[:, 0]))
            return {"periodic_edges": nodes_match(radius("edge_0"), radius("edge_sector"), NODE_TOLERANCE),
                    "gap_interface": nodes_match(angle("sliding_rotor"), angle("sliding_stator"), NODE_TOLERANCE)}
        finally:
            gmsh.finalize()


# ============================================================================================== one Elmer run
def run_elmer(workdir: Path, msh: Path, names: dict, regions, fill: float, log: Log, label: str) -> dict:
    """ElmerGrid and ElmerSolver on one mesh in a folder of their own (case.msh, case.sif, case\\, the printouts), both through run_solver. Returns
    {"grid", "solver", "vtu", "seconds", "returncodes"}; "vtu" is None if the solver did not leave a result."""
    workdir = Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(msh, workdir / "case.msh")
    started = time.perf_counter()
    grid_run = run_program(log, [paths.tool("elmergrid"), "14", "2", "case.msh", "-scale", MESH_SCALE, MESH_SCALE, MESH_SCALE, "-merge", MERGE_DISTANCE],
                           label=f"{label} ElmerGrid", cwd=workdir, timeout_s=GRID_TIMEOUT_S, threads=GRID_THREADS, name="elmergrid")
    grid = parse_grid_log(grid_run.stdout)
    out = {"grid": grid, "solver": None, "vtu": None, "returncodes": {"elmergrid": grid_run.returncode}, "seconds": None}
    if grid_run.timed_out or grid_run.returncode != 0 or not (workdir / "case" / "mesh.header").is_file():
        out["seconds"] = time.perf_counter() - started
        return out
    (workdir / "case.sif").write_text(sif_text(names, regions, fill=fill), encoding="ascii", newline="\n")
    solver_run = run_program(log, [paths.tool("elmer"), "case.sif"], label=f"{label} ElmerSolver", cwd=workdir, timeout_s=SOLVER_TIMEOUT_S,
                             threads=SOLVER_THREADS, name="elmersolver")
    out["solver"] = parse_solver_log(solver_run.stdout)
    out["solver"]["timed_out"] = solver_run.timed_out
    out["returncodes"]["elmersolver"] = solver_run.returncode
    vtu_path = workdir / "case" / "case_t0001.vtu"
    if vtu_path.is_file() and out["solver"]["finished"]:
        out["vtu"] = read_vtu(vtu_path)
    out["seconds"] = time.perf_counter() - started
    return out


def analyse_run(run: dict, names: dict, regions, p: mg.Params, theta_drawn: float, femm_samples: np.ndarray, angles: np.ndarray) -> dict:
    """The numbers of one Elmer run: the profile on the mid-gap circle against FEMM's, the flux linkages, psi_PM, and how well the periodic condition holds."""
    vtu = run["vtu"]
    br, bt = gap_profile(vtu, p.r_ro + p.airgap / 2, angles)
    comparison = compare_waveforms(femm_samples, br)
    psi = flux_linkage(vtu, names, regions, p.stack, p.n_sectors)
    d, q = park([psi["A"], psi["B"], psi["C"]], theta_drawn)
    x, y = vtu["points_mm"][:, 0], vtu["points_mm"][:, 1]
    edge_0, edge_90 = np.where(np.abs(y) < 1e-9)[0], np.where(np.abs(x) < 1e-6)[0]
    r = np.hypot(x, y)
    a0, a90 = vtu["az"][edge_0[np.argsort(r[edge_0])]], vtu["az"][edge_90[np.argsort(r[edge_90])]]
    periodic = float(np.abs(a0 - a90).max() / np.abs(vtu["az"]).max()) if len(a0) == len(a90) else None
    return {"br_T": br, "bt_T": bt, "comparison": comparison, "psi_abc_Wb": psi, "psi_d_Wb": d, "psi_q_Wb": q, "psi_pm_mWb": 1e3 * math.hypot(d, q),
            "periodic_relative_difference": periodic, "az_max": float(np.abs(vtu["az"]).max())}


# ============================================================================================== pictures
def read_field(vtu: dict) -> tuple:
    """(nodes in mm, triangles, Az at the nodes, |B| of each triangle) of an Elmer result: the nodal field averaged over each triangle. What bench.fig_field draws."""
    b = np.hypot(vtu["b"][:, 0], vtu["b"][:, 1])[vtu["triangles"]].mean(axis=1)
    return vtu["points_mm"], vtu["triangles"], vtu["az"], b


def fig_profile(angles, femm, lam, plain, target: Path) -> None:
    """B_r on the mid-gap circle over the 90-degree sector, FEMM against Elmer's two runs, and the harmonics of the three."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.5, 1]})
    for values, label, colour, style in ((femm, "FEMM 5.0 (air-gap element)", "#1f6feb", "-"), (lam["br_T"], "Elmer, laminated steel (FEMM's M-19)", "#e0782f", "-"),
                                         (plain["br_T"], "Elmer, the 47 points alone", "#2a9d5c", "--")):
        axes[0].plot(angles, values, style, color=colour, lw=1.4, label=label)
    axes[0].set(xlabel="angle in the stator frame [degrees] (one electrical period = 90)", ylabel="radial flux density B_r [T]", title="No load: B_r on the mid-gap circle (r = 45.4 mm)")
    orders = np.arange(1, 21)
    hf = b1_femm.harmonic_amplitudes(femm)
    for k, (values, label, colour) in enumerate(((femm, "FEMM", "#1f6feb"), (lam["br_T"], "Elmer laminated", "#e0782f"), (plain["br_T"], "Elmer plain", "#2a9d5c"))):
        h = b1_femm.harmonic_amplitudes(values)
        axes[1].bar(orders + (k - 1) * 0.27, [100 * h[n] / hf[1] for n in orders], width=0.27, color=colour, label=label)
    axes[1].set(xlabel="harmonic order (1 = the fundamental, 4 pole pairs)", ylabel="amplitude [% of FEMM's fundamental]", title="Harmonics of B_r", yscale="log", ylim=(0.05, 130),
                xticks=list(range(1, 21)))
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(target, dpi=130)
    plt.close(fig)


# ============================================================================================== the case
def run(folder: Path, log: Log) -> dict:
    """Run B3 in `folder` (already empty). Returns the results dict whose "checks" list every verdict."""
    folder = Path(folder)
    p = mg.Params(barrier="none")
    regions = mg.build(p, gap="sliding")
    started = time.perf_counter()
    # ---- the mesh: B2's path 1, the same function and size
    made = b2_gmsh.make_path1_mesh(folder / "path1.msh")
    names = physical_names(folder / "path1.msh")
    facts = interface_facts(folder / "path1.msh")
    log.say(f"B3: B2's mesh, {made['nodes_in_file']:,} nodes in the file ({made['to_gmsh_nodes']:,} counted by to_gmsh, {made['nodes_without_element']} of them arc centres "
            f"that no element uses); periodic edges {facts['periodic_edges']['count']} nodes, match {facts['periodic_edges']['match']}; "
            f"the gap's two sides {facts['gap_interface']['count']} nodes, match {facts['gap_interface']['match']}")
    # ---- FEMM: B1's model at the drawn position, and psi_PM as B1 defines it
    reference = b1_femm.reference(folder, log, "bench B3 FEMM reference")
    angles = np.asarray(reference["drawn"]["gap_deg"], float)
    femm_br = np.asarray(reference["drawn"]["gap_br"], float)
    theta_drawn = reference["drawn"]["theta_e"]
    femm_psi_here = 1e3 * math.hypot(reference["drawn"]["psi_d"], reference["drawn"]["psi_q"])
    log.say(f"B3: FEMM at the drawn position (theta_e {theta_drawn:.2f}): B_r fundamental {b1_femm.fundamental(femm_br):.5f} T, |psi_dq| {femm_psi_here:.4f} mWb; "
            f"psi_PM over the 120 positions {reference['psi_pm_mWb']:.4f} mWb ({reference['drawn']['nodes']:,} nodes)")
    # ---- Elmer: the laminated steel (the verdict) and the plain 47 points (the brief's words, with its own verdicts)
    runs: dict[str, dict] = {}
    for key, fill in (("laminated", mg.M19_STEEL_FEMM42["LamFill"]), ("plain", 1.0)):
        runs[key] = run_elmer(folder / f"elmer_{key}", folder / "path1.msh", names, regions, fill, log, f"bench B3 Elmer {key}")
        if runs[key]["vtu"] is None:
            raise RuntimeError(f"Elmer's {key} run left no result: grid {runs[key]['grid']}, solver {runs[key]['solver']}; see {folder / ('elmer_' + key)}")
        runs[key]["numbers"] = analyse_run(runs[key], names, regions, p, theta_drawn, femm_br, angles)
        n = runs[key]["numbers"]
        g = runs[key]["grid"]
        log.say(f"B3 Elmer {key}: ElmerGrid read {_fmt(g['nodes_read'], ',')} nodes, merged {g['nodes_merged']} across the gap, {_fmt(g['nodes_merged_to'], ',')} for Elmer; "
                f"{runs[key]['solver']['iterations']} nonlinear iterations, B_r fundamental {n['comparison']['fundamental_other']:.5f} T "
                f"({n['comparison']['fundamental_percent']:+.2f} % against FEMM), psi_PM {n['psi_pm_mWb']:.4f} mWb, {runs[key]['seconds']:.1f} s")
    # ---- the same, on finer meshes (laminated)
    sensitivity: list[dict] = []
    for size in SENSITIVITY_SIZES:
        scratch = folder / f"path1_s{size:g}.msh"
        made_s = b2_gmsh.make_path1_mesh(scratch, size)
        edge = interface_facts(scratch)
        finer_run = finer_numbers = None
        if edge["periodic_edges"]["match"] and edge["gap_interface"]["match"]:
            r = run_elmer(folder / f"elmer_s{size:g}", scratch, physical_names(scratch), regions, mg.M19_STEEL_FEMM42["LamFill"], log, f"bench B3 Elmer size {size:g}")
            if r["vtu"] is not None:
                finer_run, finer_numbers = r, analyse_run(r, physical_names(scratch), regions, p, theta_drawn, femm_br, angles)
        sensitivity.append(sensitivity_row(size, made_s["nodes_in_file"], edge, finer_run, finer_numbers, reference["psi_pm_mWb"]))
        scratch.unlink()
    results = {"case": CASE, "title": TITLE, "verdict_run": VERDICT_RUN, "params": params(p),
               "mesh": {"made": made, "interfaces": facts, "size_max_mm": made["size_max"], "node_counts": node_counts(made, runs[VERDICT_RUN]["grid"])},
               "elmer": {k: {"grid": v["grid"], "solver": v["solver"], "seconds": round(v["seconds"], 2), "returncodes": v["returncodes"]} for k, v in runs.items()},
               "femm": {"theta_e": theta_drawn, "psi_pm_mWb": reference["psi_pm_mWb"], "psi_dq_here_mWb": femm_psi_here, "nodes": reference["drawn"]["nodes"],
                        "psi_abc_mWb": [1e3 * reference["drawn"][k] for k in ("psi_a", "psi_b", "psi_c")]},
               "profile": {"angles_deg": angles, "femm_br_T": femm_br, "elmer_laminated_br_T": runs["laminated"]["numbers"]["br_T"],
                           "elmer_plain_br_T": runs["plain"]["numbers"]["br_T"]},
               "numbers": {k: {n: v for n, v in runs[k]["numbers"].items() if n not in ("br_T", "bt_T")} for k in runs}, "sensitivity": sensitivity}
    results["checks"] = make_checks(runs, facts, reference, femm_psi_here, sensitivity, made)
    lam = runs[VERDICT_RUN]                                    # a title of short lines (bench.fig_field wraps any line over its width, but a wrap can leave an orphan word)
    fig_field(*read_field(lam["vtu"]), f"Elmer 26.2 (ElmerGrid, ElmerSolver, MagnetoDynamics2D): the band-free rotor at no load\n"
              f"the drawn position; flux density |B| and flux lines (contours of Az)\n"
              f"B2's mesh as Elmer has it: {lam['grid']['nodes_merged_to']:,} nodes ({made['nodes_in_file']:,} in the file, {lam['grid']['nodes_merged']} merged across the gap)\n"
              f"{lam['solver']['iterations']} nonlinear iterations; psi_PM {lam['numbers']['psi_pm_mWb']:.2f} mWb (FEMM {reference['psi_pm_mWb']:.2f}), "
              f"B_r fundamental {lam['numbers']['comparison']['fundamental_other']:.4f} T ({lam['numbers']['comparison']['fundamental_percent']:+.2f} %)", folder / "figure.png")
    fig_profile(angles, femm_br, runs["laminated"]["numbers"], runs["plain"]["numbers"], folder / "profile.png")
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    return results


def _pair_text(row: dict) -> str:
    return f"{row['order']}: {row['femm_percent']:.1f} / {row['other_percent']:.1f}"


def _fmt(value, spec: str) -> str:
    """`value` formatted with `spec`, or "n/a" if the program did not give it: a failed check has to say so, not raise while it is being worded."""
    return "n/a" if value is None else format(value, spec)


def node_counts(made: dict, grid: dict) -> dict:
    """The node counts of B2's mesh as B3 meets them, each with what it counts (they are three different numbers, and the 8,906 that to_gmsh returns is not 'the nodes in
    the file'). `made` is b2_gmsh.make_path1_mesh's answer, `grid` parse_grid_log's reading of ElmerGrid's printout of the run that is the verdict's.
      nodes_in_file          the nodes of the saved path1.msh: B2's mesh (8,819)
      to_gmsh_nodes          what motorgeo.to_gmsh returned: the nodes of its Gmsh session (8,906)
      nodes_without_element  the difference: the centres of the arcs, geometry that no element uses and the file does not hold (87)
      elmergrid_read         the nodes ElmerGrid says it read ('N original nodes'); the same file, so the same number as nodes_in_file
      merged_in_the_gap      the nodes ElmerGrid merged: the two sides of the sliding air gap, closed at the fixed rotor position (269)
      elmer_nodes            what is left for ElmerSolver (8,550)
    elmergrid_read_b2s_file says whether B3 and B2 mean the same mesh by 'B2's mesh'; meaning has one line for each count."""
    return {"nodes_in_file": made["nodes_in_file"], "to_gmsh_nodes": made["to_gmsh_nodes"], "nodes_without_element": made["nodes_without_element"],
            "elmergrid_read": grid["nodes_read"], "merged_in_the_gap": grid["nodes_merged"], "elmer_nodes": grid["nodes_merged_to"],
            "elmergrid_read_b2s_file": grid["nodes_read"] == made["nodes_in_file"],
            "meaning": {"nodes_in_file": "the nodes of the saved path1.msh, B2's mesh",
                        "to_gmsh_nodes": "what motorgeo.to_gmsh returned: the nodes of its Gmsh session",
                        "nodes_without_element": "the centres of the arcs: geometry that no element uses, not in the file",
                        "elmergrid_read": "the nodes ElmerGrid read from the file ('original nodes')",
                        "merged_in_the_gap": "the nodes ElmerGrid merged: the gap's rotor side onto its stator side",
                        "elmer_nodes": "the nodes ElmerSolver works with, after the merge"}}


def sensitivity_row(size: float, nodes_in_file: int, edge: dict, elmer_run: dict | None = None, numbers: dict | None = None, femm_psi_mwb: float | None = None) -> dict:
    """One row of the finer-mesh table, as results.json keeps it. The mesh of path 1 at `size` mm has `nodes_in_file` nodes in its file (what B2's refinement table
    counts, 10,711 at size_max 1.0); `edge` is interface_facts() of it. If Elmer could use the mesh (`elmer_run`, run_elmer's answer, and its analyse_run `numbers`)
    the row also has Elmer's side, each count named for what it counts: elmergrid_read (the file's nodes as ElmerGrid read them), merged_in_the_gap, and elmer_nodes,
    what is left after the merge (10,442 at size_max 1.0). There is no key called `nodes`: it was Elmer's count and looked like the file's (P3 of the re-review)."""
    row = {"size_max_mm": size, "nodes_in_file": nodes_in_file, "periodic_edges_match": edge["periodic_edges"]["match"], "gap_interface_match": edge["gap_interface"]["match"]}
    if elmer_run is not None and numbers is not None:
        grid, comparison = elmer_run["grid"], numbers["comparison"]
        row.update({"elmergrid_read": grid["nodes_read"], "merged_in_the_gap": grid["nodes_merged"], "elmer_nodes": grid["nodes_merged_to"],
                    "iterations": elmer_run["solver"]["iterations"], "fundamental_percent": comparison["fundamental_percent"], "psi_pm_mWb": numbers["psi_pm_mWb"],
                    "psi_percent": 100 * (numbers["psi_pm_mWb"] / femm_psi_mwb - 1), "rms_percent_of_fundamental": comparison["rms_percent_of_fundamental"],
                    "pattern_alike": comparison["alike"]})
    return row


def _verdict_only_for(key: str) -> dict:
    """The keywords of bench.check and bench.compare for a check about the Elmer run `key`: it gates when it is the run that carries the verdict (VERDICT_RUN, the
    laminated steel: the same steel as FEMM's), and is reported, not gating, when it is the other (the 47 points alone, a different steel)."""
    return {} if key == VERDICT_RUN else {"gating": False, "reason": f"the verdict is the {VERDICT_RUN} run's"}


def make_checks(runs: dict, facts: dict, reference: dict, femm_psi_here: float, sensitivity: list[dict], made: dict) -> list[dict]:
    """Every number of B3 against its expectation, with the source of the expectation.

    The checks about the `laminated` run, the mesh and the finer meshes gate the verdict. Every check about the `plain` run (and the one that compares the two
    runs) is reported beside them with its own pass or fail and does not gate: the plain run solves another steel, so its miss is not Elmer's.
    `made` is b2_gmsh.make_path1_mesh's answer for the mesh B3 used: the mesh check names its node counts and holds ElmerGrid to having read that file."""
    out: list[dict] = []
    lam = runs["laminated"]
    grid = lam["grid"]
    out.append(check("mesh: the periodic edges and the two sides of the gap carry matching nodes, and ElmerGrid makes one conforming mesh",
                     facts["periodic_edges"]["match"] and facts["gap_interface"]["match"] and grid["conforming"] is True
                     and grid["nodes_merged"] == facts["gap_interface"]["count"][0] and grid["nodes_read"] == made["nodes_in_file"],
                     f"edge_0 and edge_sector: the same radii; sliding_rotor and sliding_stator: the same angles, to within {NODE_TOLERANCE:g} ({OWN_LIMIT}); "
                     f"ElmerGrid: reads the nodes of B2's file, one piece, the gap's nodes merged",
                     f"edges {facts['periodic_edges']['count']} nodes, largest difference {_fmt(facts['periodic_edges']['max_difference'], '.1e')} mm; gap sides "
                     f"{facts['gap_interface']['count']} nodes, {_fmt(facts['gap_interface']['max_difference'], '.1e')} deg. B2's mesh: {made['nodes_in_file']:,} nodes in the file "
                     f"({made['to_gmsh_nodes']:,} counted by to_gmsh, {made['nodes_without_element']} of them arc centres that no element uses); ElmerGrid read "
                     f"{_fmt(grid['nodes_read'], ',')}, merged {grid['nodes_merged']} (the gap) and left {_fmt(grid.get('nodes_merged_to'), ',')} for Elmer, pieces: {grid['pieces']}",
                     "interface_facts(), ElmerGrid's own printout, make_path1_mesh()"))
    for key, label in (("laminated", "laminated steel"), ("plain", "plain 47 points")):
        s, n = runs[key]["solver"], runs[key]["numbers"]
        out.append(check(f"{key}: ElmerSolver converges", s["finished"] and not s["errors"] and s["last_change"] is not None and s["last_change"] < CONVERGENCE_TOLERANCE
                         and s["iterations"] < MAX_ITERATIONS and not s["timed_out"],
                         f"ALL DONE with the relative change below {CONVERGENCE_TOLERANCE:g} in fewer than {MAX_ITERATIONS} nonlinear iterations ({OWN_LIMIT}); "
                         f"the brief says only that the solver converges and the iteration count is reported",
                         f"{s['iterations']} nonlinear iterations (Picard, then Newton from the second), last relative change {_fmt(s['last_change'], '.2e')}; "
                         f"{s['solver_seconds_real']} s on {s['threads']} thread(s); Elmer {s['version']} rev {s['revision']}; warnings: {len(s['warnings'])}",
                         "ElmerSolver's printout (ComputeChange lines)", **_verdict_only_for(key)))
        out.append(check(f"{key}: the periodic condition holds in the solution", n["periodic_relative_difference"] is not None and n["periodic_relative_difference"] < 1e-9,
                         f"Az on edge_0 equals Az on edge_sector at every radius (the conforming condition is applied), to a relative difference below 1e-9 ({OWN_LIMIT})",
                         f"largest difference {_fmt(n['periodic_relative_difference'], '.1e')} of the largest |Az| ({n['az_max']:.3e} Wb/m)", "case\\case_t0001.vtu", **_verdict_only_for(key)))
    for key, label in (("laminated", "Elmer with FEMM's steel: the 47 points and the lamination fill 0.98"), ("plain", "Elmer with the 47 points alone, no lamination fill")):
        c, n = runs[key]["numbers"]["comparison"], runs[key]["numbers"]
        out.append(compare(f"{key}: B_r on the mid-gap circle, fundamental, Elmer against FEMM", c["fundamental_other"], c["fundamental_femm"], FUNDAMENTAL_TOLERANCE_PERCENT,
                           f"FEMM 5.0, B1's model at the drawn position, mo_getgapb on the same circle ({label})", "T", **_verdict_only_for(key)))
        top = c["top_orders_femm"]
        out.append(check(f"{key}: the harmonic pattern of B_r is alike", c["alike"],
                         f"every harmonic 2 to 20 within {PATTERN_POINTS:g} points of its own fundamental, the shapes within {PATTERN_RMS_PERCENT:g} % rms, and the "
                         f"{PATTERN_TOP} largest harmonics the same orders (the bench's own limits: the brief gives none)",
                         f"largest differences {max(c['harmonics'], key=lambda r: abs(r['difference_points']))['difference_points']:+.2f} points at order "
                         f"{max(c['harmonics'], key=lambda r: abs(r['difference_points']))['order']}; shapes {c['rms_percent_of_fundamental']:.2f} % rms; the largest harmonics "
                         f"(order: FEMM / Elmer, % of the fundamental): " + ", ".join(_pair_text(next(r for r in c['harmonics'] if r['order'] == o)) for o in top)
                         + f"; same orders: {c['top_orders_femm']} / {c['top_orders_other']}", "FEMM's gap field, 360 samples", **_verdict_only_for(key)))
        out.append(compare(f"{key}: psi_PM from Az in the coil sides, against FEMM's psi_PM (the fundamental over 120 positions)", n["psi_pm_mWb"], reference["psi_pm_mWb"],
                           PSI_TOLERANCE_PERCENT, "FEMM 5.0, B1's model: the fundamental of phase A's flux linkage over 120 rotor positions", "mWb", **_verdict_only_for(key)))
        out.append(compare(f"{key}: psi_PM against FEMM's flux vector at the same rotor position (the tighter comparison)", n["psi_pm_mWb"], femm_psi_here, PSI_TOLERANCE_PERCENT,
                           "FEMM 5.0, the same model and the same position: the length of (psi_d, psi_q)", "mWb", **_verdict_only_for(key)))
    shift = runs["laminated"]["numbers"]["comparison"]["fundamental_percent"] - runs["plain"]["numbers"]["comparison"]["fundamental_percent"]
    out.append(check("the lamination fill is what separates the plain run from FEMM", shift > 1.5,
                     "the same Elmer, given FEMM's lamination fill (0.98), moves toward FEMM by about what the fill does in FEMM itself (2.4 % in psi, found by solving FEMM with "
                     f"fill 1.0): by more than 1.5 points ({OWN_LIMIT})",
                     f"the B_r fundamental moves by {shift:+.2f} points (from {runs['plain']['numbers']['comparison']['fundamental_percent']:+.2f} % to "
                     f"{runs['laminated']['numbers']['comparison']['fundamental_percent']:+.2f} % against FEMM)", "the two Elmer runs", **_verdict_only_for("plain")))
    done = [row for row in sensitivity if "fundamental_percent" in row]
    out.append(check("laminated: the agreement holds on finer meshes", bool(done) and all(abs(r["fundamental_percent"]) < FUNDAMENTAL_TOLERANCE_PERCENT
                                                                                           and abs(r["psi_percent"]) < PSI_TOLERANCE_PERCENT for r in done),
                     f"B_r fundamental and psi_PM within {FUNDAMENTAL_TOLERANCE_PERCENT:g} % of FEMM's at every mesh size whose periodic edges match",
                     "; ".join(f"size_max {r['size_max_mm']:g}: " + (f"{r['elmer_nodes']:,} nodes in Elmer after ElmerGrid's merge ({r['nodes_in_file']:,} in the file, "
                                                                   f"{r['merged_in_the_gap']} merged across the gap), B_r {r['fundamental_percent']:+.2f} %, psi {r['psi_percent']:+.2f} %, shapes "
                                                                   f"{r['rms_percent_of_fundamental']:.1f} % rms" if "fundamental_percent" in r else
                                                                   "skipped, the periodic edges or the gap's sides do not match") for r in sensitivity), "the same chain on finer meshes"))
    return out
