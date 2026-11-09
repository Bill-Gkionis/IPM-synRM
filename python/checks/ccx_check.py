"""ccx_check: CalculiX solves a ring that spins, and the stress at its bore is compared with Lame's exact answer.

What it does, and why
    A structural solver is trusted when it reproduces a case whose answer is known on paper. A steel ring that spins is such a
    case: every bit of steel is pulled outward by the centrifugal force, the ring stretches, and the stretching puts a hoop
    (tangential) stress in it that is largest at the inner edge, the bore. Lame worked that stress out exactly in 1852 for a flat
    ring with free edges. The motor's rotor is a ring of this size (shaft radius 7.5 mm, rotor radius 45 mm) that spins at up to
    10,000 rpm, so the same case is the first check of the bridge-stress work that comes later. It is the ring of the appendix
    of the project note on finite-element tools (FEA), with CalculiX 2.23 in place of 2.22.

Method (what the code does, in the order it does it)
    1. Lame's formula, for a thin flat ring, inner radius a, outer radius b, density rho, angular speed w, Poisson ratio nu:
           k      = (3 + nu) / 8 * rho * w^2
           hoop   = k * (a^2 + b^2 + a^2 b^2 / r^2 - (1 + 3 nu) / (3 + nu) * r^2)
           radial = k * (a^2 + b^2 - a^2 b^2 / r^2 - r^2)
       At the bore (r = a) with steel at 10,000 rpm the hoop stress is 14.098 MPa; at the rim 3.362 MPa. The radial stress is
       zero at both edges, because nothing pushes on a free edge.
    2. A quarter of the ring is meshed by this code (structured, 8-node quadrilaterals, plane stress, 96 x 96 elements, 28,033
       nodes). The other three quarters are not needed: the cut edges are held by symmetry (a node on the x axis may not move
       in y, a node on the y axis may not move in x). The load is the centrifugal force, *DLOAD ... CENTRIF.
    3. The input deck is written in units of mm, N, tonne, s (so stress comes out in MPa and density in tonne/mm^3).
    4. ccx_static.exe solves it through femm_queue.run_solver (below-normal priority, a timeout, one at a time on the PC). It
       writes ring.frd; this code reads the stress at every node of the bore, turns each into its hoop component
       (sigma_t = sxx sin^2 + syy cos^2 - 2 sxy sin cos, with the angle of the node) and averages them.
    5. The row passes when the mean bore hoop stress is within 0.5 % of Lame's.
    CalculiX 2.23's zip has no plain ccx.exe: ccx_static.exe (SPOOLES) is the build that runs on this PC (config.toml says so).

Numbers to expect (this PC, 2026-11-06, CalculiX 2.23): Lame 14.0977 MPa at the bore; CalculiX 14.1060 MPa (+0.06 %), the worst node
    0.06 % off; the rim 3.362 MPa; about 5 s. The coarser meshes of that note give 1.5 % (133 nodes), 0.40 % (1,825), 0.16 %
    (7,105): the error falls with the mesh, which is what a correct solver does.

Self-check (the answers follow).
  1. Why is it enough to mesh a quarter of the ring?
  2. The hoop stress at the bore is 3.5 times the one at the rim, although the steel at the rim moves faster. Why?
  3. What would the stress be at twice the speed, and why?
Answers.
  1. The ring and its load are the same after a quarter turn and after a reflection, so the stress field repeats. A node on a cut
     edge can slide along the edge but not across it; that is all the missing three quarters would do.
  2. The rim is free to stretch while the bore is held back by the ring's inner steel, which wants to stretch less: the bore
     carries the hoop load of the whole ring. The hoop stress is a difference of terms that is largest at the smallest radius.
  3. Four times: every term of Lame's formula is proportional to w^2 (the centrifugal load is rho w^2 r), and (2w)^2 = 4 w^2.
"""
from __future__ import annotations

import math
import re
import time
from pathlib import Path

import femm_queue

from . import PASS, FAIL, Row, fresh, need_tool, save_log

BORE, OUTER = 7.5, 45.0                     # mm: the shaft's radius and the rotor's
E, NU, RHO = 200e3, 0.3, 7.65e-9            # MPa, -, tonne/mm^3 (steel)
RPM = 10_000
OMEGA = 2 * math.pi * RPM / 60              # rad/s
MESH = (96, 96)                             # elements across the ring, and along the quarter arc
TOLERANCE_PERCENT = 0.5
JOB = "ring"
TIMEOUT_S = 180


def lame_hoop(r: float) -> float:
    """Hoop stress [MPa] at radius r [mm] in the spinning ring (plane stress, free edges)."""
    k = (3 + NU) / 8 * RHO * OMEGA ** 2
    return k * (BORE ** 2 + OUTER ** 2 + BORE ** 2 * OUTER ** 2 / r ** 2 - (1 + 3 * NU) / (3 + NU) * r ** 2)


def write_deck(path: Path, nr: int, nt: int) -> dict[int, tuple[float, float]]:
    """Write the quarter ring as a CalculiX input deck (nr x nt CPS8 elements). Returns {node number: (x, y)}."""
    def nid(i: int, j: int) -> int:
        return i * (2 * nt + 1) + j + 1

    lines = ["*HEADING", f"spinning ring, plane stress, {nr}x{nt} CPS8", "*NODE"]
    coords: dict[int, tuple[float, float]] = {}
    for i in range(2 * nr + 1):
        r = BORE + (OUTER - BORE) * i / (2 * nr)
        for j in range(2 * nt + 1):
            if i % 2 and j % 2:                                     # the centre of an 8-node element is not a node
                continue
            angle = math.pi / 2 * j / (2 * nt)
            x, y = r * math.cos(angle), r * math.sin(angle)
            coords[nid(i, j)] = (x, y)
            lines.append(f"{nid(i, j)}, {x:.9f}, {y:.9f}, 0.0")
    lines.append("*ELEMENT, TYPE=CPS8, ELSET=EALL")
    count = 0
    for ie in range(nr):
        for je in range(nt):
            i0, j0 = 2 * ie, 2 * je
            count += 1
            corner_and_mid = (nid(i0, j0), nid(i0 + 2, j0), nid(i0 + 2, j0 + 2), nid(i0, j0 + 2),
                              nid(i0 + 1, j0), nid(i0 + 2, j0 + 1), nid(i0 + 1, j0 + 2), nid(i0, j0 + 1))
            lines.append(f"{count}, " + ", ".join(map(str, corner_and_mid)))

    def node_set(name: str, ids: list[int]) -> list[str]:
        return [f"*NSET, NSET={name}"] + [", ".join(map(str, ids[k:k + 8])) for k in range(0, len(ids), 8)]

    lines += node_set("EDGE0", [nid(i, 0) for i in range(2 * nr + 1)])            # on the x axis: no motion in y
    lines += node_set("EDGE90", [nid(i, 2 * nt) for i in range(2 * nr + 1)])      # on the y axis: no motion in x
    lines += ["*MATERIAL, NAME=STEEL", "*ELASTIC", f"{E}, {NU}", "*DENSITY", f"{RHO}",
              "*SOLID SECTION, ELSET=EALL, MATERIAL=STEEL", "1.0",
              "*BOUNDARY", "EDGE0, 2, 2", "EDGE90, 1, 1",
              "*STEP", "*STATIC",
              "*DLOAD", f"EALL, CENTRIF, {OMEGA * OMEGA:.6e}, 0., 0., 0., 0., 0., 1.",
              "*NODE FILE", "U", "*EL FILE", "S", "*END STEP"]
    path.write_text("\n".join(lines) + "\n", encoding="ascii")
    return coords


def read_frd_stress(path: Path) -> dict[int, list[float]]:
    """The nodal stresses of an ASCII .frd file: {node: [sxx, syy, szz, sxy, syz, szx]} (fixed-width columns of 12 characters)."""
    stress: dict[int, list[float]] = {}
    block = None
    with open(path, encoding="ascii", errors="replace") as handle:
        for line in handle:
            if line.startswith(" -4"):
                block = line.split()[1]
            elif line.startswith(" -3"):
                block = None
            elif block == "STRESS" and line.startswith(" -1"):
                stress[int(line[3:13])] = [float(line[13 + 12 * k:25 + 12 * k]) for k in range(6)]
    return stress


def bore_hoop(coords: dict[int, tuple[float, float]], stress: dict[int, list[float]]) -> list[float]:
    """The hoop stress [MPa] at every node of the bore."""
    values = []
    for node, (x, y) in coords.items():
        r = math.hypot(x, y)
        if abs(r - BORE) < 1e-6 and node in stress:
            sxx, syy, _, sxy, _, _ = stress[node]
            c, s = x / r, y / r
            values.append(sxx * s * s + syy * c * c - 2 * sxy * s * c)
    return values


def run(folder: Path):
    ccx = need_tool("ccx")
    fresh(folder)
    nr, nt = MESH
    coords = write_deck(folder / f"{JOB}.inp", nr, nt)
    case = f"spinning ring, {RPM:,} rpm, {nr}x{nt} CPS8 elements (the FEA note's ring)"
    expected = f"bore hoop stress within {TOLERANCE_PERCENT} % of Lame's {lame_hoop(BORE):.3f} MPa"
    started = time.perf_counter()
    run_ = femm_queue.run_solver([ccx, "-i", JOB], label="check_tools CalculiX ring", cwd=folder, timeout_s=TIMEOUT_S, threads=2)
    seconds = round(time.perf_counter() - started, 2)
    save_log(folder, "ccx", run_)                       # ccx.out.txt: what the solver said (its version banner, the step it ran)
    frd = folder / f"{JOB}.frd"
    if run_.timed_out or run_.returncode != 0 or not frd.is_file():
        why = "timed out after %d s" % TIMEOUT_S if run_.timed_out else f"exit code {run_.returncode}, no {frd.name}" if not frd.is_file() else f"exit code {run_.returncode}"
        tail = " | ".join((run_.stdout + run_.stderr).strip().splitlines()[-3:])
        yield Row("CalculiX", case, expected, f"ccx failed: {why}. {tail}", FAIL, seconds)
        return
    hoop = bore_hoop(coords, read_frd_stress(frd))
    if not hoop:
        yield Row("CalculiX", case, expected, "no stress was found at the bore nodes of the result file", FAIL, seconds)
        return
    exact = lame_hoop(BORE)
    mean = sum(hoop) / len(hoop)
    error = 100 * (mean - exact) / exact
    worst = 100 * max(abs(h - exact) for h in hoop) / exact
    banner = re.search(r"Version\s+(\d+\.\d+)", run_.stdout)
    version = banner.group(1) if banner else "?"
    measured = f"ccx {version}: {mean:.4f} MPa at the bore ({error:+.2f} %; worst of {len(hoop)} nodes {worst:.2f} %), {len(coords):,} nodes"
    yield Row("CalculiX", case, expected, measured, PASS if abs(error) < TOLERANCE_PERCENT else FAIL, seconds)
