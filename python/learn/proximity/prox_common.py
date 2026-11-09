"""learn/proximity/prox_common.py: the small FEMM models and the helpers the proximity-effect demos share.

What it is, and why
    The three FEMM demos of this folder (01_one_wire.py, 02_two_wires.py, 03_slot_strands.py) draw round copper wires in air and read FEMM's time-harmonic solution back as a picture of the
    current density. The drawing, the solving and the reading are the same for all, so they are here, once, as plain functions of a few lines each. Nothing here is specific to a motor: a wire is a
    circle of copper (conductivity 1 / rho), a circuit is a way of saying how much current it carries, and the model's outer boundary is far enough away to leave the wires alone.

How FEMM's time-harmonic solver is used (the vocabulary of the three demos)
    frequency    mi_probdef(f, ...): f in Hz; the solver then finds a COMPLEX vector potential A (the z component of the magnetic potential, B = curl A) at every node. Every phasor here is a PEAK value.
    conductor    a block with a conductivity: FEMM adds the eddy currents J = -j w sigma A (plus a constant that makes the total current what the circuit says).
    circuit      mi_addcircprop(name, current, type): 'series' (type 1) makes the block carry exactly that current; with current 0 the wire carries no net current and only eddy currents flow in it
                 (the proximity effect alone); 'parallel' (type 0) joins several wires to a common voltage and lets the current divide.
    field B      the outer boundary can be given A = -B x, a uniform field B along y far from the wires (units: FEMM takes x in mm, so A1 = -B * 1e-3).
    reading      the current density in a conductor, from the nodal A of the solution file: J = I / S - j w sigma (A - mean A over the wire) with S the wire's area: the part I / S is the uniform
                 DC distribution and the other part is the eddy current, which has no net current. (FEMM's own block integrals agree to 1e-4; optimization/ac_slot.py checks it.)

The functions
    wire_task(task)        the pool function (runs in a worker on its own FEMM): draws the wires, solves, returns each wire's loss and current and the current density of every triangle in them.
    solve(tasks)           runs tasks through the project's solver slot (femm_queue.femm_pool: one heavy run at a time, below-normal priority, a timeout): a few seconds each.
    plot_wires(ax, ...)    draws the triangles of the wires coloured by a value (matplotlib's PolyCollection), the surface and ink colours of the project.
    The task is a dict: {'tag', 'f', 'wires': [{'x', 'y', 'd', 'i' (peak, complex), 'circuit' (a name)}], 'circuits': {name: (current, 'series' | 'parallel')}, 'B' (tesla or None),
    'depth_mm', 'radius_mm' (of the outer boundary), 'mesh_over_d', 'arc_deg'}.

Self-check (the answers follow).
  1. A wire in a series circuit of 0 A has current density everywhere in it. Where does it come from, and what is its total?
  2. Why is the outer boundary of the model drawn at 40 or 100 radii of the wire and not at 5?
Answers.
  1. From the field: a changing B across the wire drives an axial electric field -y dB/dt (y across the wire), and the conductor answers with a current density that is plus on one side and minus on
     the other. The circuit of 0 A says the total is zero, and the current density is the loop that closes through the wire's own two halves.
  2. A wire alone carries a field that falls as 1 / r outside it; with the boundary at A = 0 five radii away the model would hold the field back and change the current in the wire. At 40 radii the change
     is under a part in a thousand (optimization/ac_calibration.py: its skin-effect check at 40 radii agrees with the exact answer to 8e-4 or better). For the uniform-field case the wire's own response,
     a dipole falling as 1 / r^2, has to be held at the boundary too, so the box is 100 radii wide; its convergence stage finds that going from 40 to 100 radii moves the loss by less than 1e-4.
"""
from __future__ import annotations

import math
import os
import sys
import time
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))          # the code folder, wherever this is run from

import numpy as np

import femm_queue
import paths
from optimization import ac_theory

OUT = paths.path("out") / "learn" / "proximity"
T_COPPER = 100.0
RHO = float(ac_theory.resistivity(T_COPPER))
SURFACE, INK, INK2, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#898781"


def frequency_for(ratio: float, d_mm: float, t_c: float = T_COPPER) -> float:
    """The frequency [Hz] at which a wire of d_mm has d / delta = ratio in copper at t_c: f = rho / (pi mu0 delta^2) with delta = d / ratio."""
    return float(ac_theory.resistivity(t_c)) / (math.pi * ac_theory.MU0 * (d_mm * 1e-3 / ratio) ** 2)


def wire_task(task: dict) -> dict:
    """One pool task: wires in air in a round box, solved at one frequency (module docstring). Returns {'status', 'loss_W' (per wire), 'current_A' (complex, per wire), 'area_m2',
    'tri' (n, 3, 2) the corners [mm] of every triangle inside a wire, 'wire' (n,) which wire it belongs to, 'J' (n, 3) the complex current density [A/m^2] at its corners, 'nodes'}."""
    import femm
    started = time.perf_counter()
    fem = os.path.join(task["fem_dir"], f"{task['tag']}.fem")
    os.makedirs(task["fem_dir"], exist_ok=True)
    ans = fem[:-4] + ".ans"
    try:
        rho = float(ac_theory.resistivity(T_COPPER))
        arc = task.get("arc_deg", 2.0)
        femm.newdocument(0)
        femm.mi_probdef(task["f"], "millimeters", "planar", 1e-10, task.get("depth_mm", 20.0), 30)
        femm.mi_getmaterial("Air")
        femm.mi_addmaterial("Cu", 1, 1, 0, 0, 1.0 / rho / 1e6, 0, 0, 1, 0, 0, 0, 0, 0)
        radius = task["radius_mm"]
        B = task.get("B")
        femm.mi_addboundprop("OUT", 0, 0.0 if B is None else -B * 1e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        femm.mi_addnode(radius, 0)
        femm.mi_addnode(-radius, 0)
        femm.mi_addarc(radius, 0, -radius, 0, 180, 5.0)
        femm.mi_addarc(-radius, 0, radius, 0, 180, 5.0)
        femm.mi_selectarcsegment(0, radius)
        femm.mi_selectarcsegment(0, -radius)
        femm.mi_setarcsegmentprop(5.0, "OUT", 0, 0)
        femm.mi_clearselected()
        for name, (current, kind) in task["circuits"].items():
            femm.mi_addcircprop(name, current, 1 if kind == "series" else 0)
        for w in task["wires"]:
            r = w["d"] / 2.0
            n_seg = 2 * math.ceil(180.0 / arc)
            r *= math.sqrt(2.0 * math.pi / (n_seg * math.sin(2.0 * math.pi / n_seg)))      # a polygon of the circle's area
            x, y = w["x"], w["y"]
            femm.mi_addnode(x + r, y)
            femm.mi_addnode(x - r, y)
            femm.mi_addarc(x + r, y, x - r, y, 180, arc)
            femm.mi_addarc(x - r, y, x + r, y, 180, arc)
            femm.mi_addblocklabel(x, y)
            femm.mi_selectlabel(x, y)
            femm.mi_setblockprop("Cu", 0, w["d"] * task.get("mesh_over_d", 1.0 / 12.0), w["circuit"], 0, 0, 1)
            femm.mi_clearselected()
        femm.mi_addblocklabel(0.0, radius / 2.0)
        femm.mi_selectlabel(0.0, radius / 2.0)
        femm.mi_setblockprop("Air", 1, 0, "", 0, 0, 0)
        femm.mi_clearselected()
        femm.mi_saveas(fem)
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        loss, cur, area = [], [], []
        for w in task["wires"]:
            femm.mo_selectblock(w["x"], w["y"])
            loss.append(float(np.real(femm.mo_blockintegral(4))))
            cur.append(complex(femm.mo_blockintegral(7)))
            area.append(float(np.real(femm.mo_blockintegral(5))))
            femm.mo_clearblock()
        nodes = int(femm.mo_numnodes())
        femm.mo_close()
        tri, wire, jj = _read_wires(ans, task["wires"], cur, area, task["f"], rho)
    except Exception as error:
        return {"status": "error", "tag": task["tag"], "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4)}
    finally:
        try:
            femm.mi_close()
        except Exception:
            pass
        for path in (ans, fem):
            try:
                os.remove(path)
            except OSError:
                pass
    return {"status": "ok", "tag": task["tag"], "loss_W": loss, "current_A": cur, "area_m2": area, "tri": tri, "wire": wire, "J": jj, "nodes": nodes, "seconds": time.perf_counter() - started}


def _read_wires(ans: str, wires: list[dict], cur: list, area: list, f: float, rho: float):
    """The triangles inside each wire and the current density at their corners, from the solution file (J = I / S - j w sigma (A - mean A over the wire))."""
    from optimization import ac_slot
    sol = ac_slot.read_ans(ans)
    centres = np.array([[w["x"], w["y"]] for w in wires])
    from scipy.spatial import cKDTree
    dist, nearest = cKDTree(centres).query(np.array(sol["labels"], dtype=float))
    lab_to_wire = np.where(dist < 1e-4, nearest, -1)
    tri_wire = lab_to_wire[sol["label"]]
    ok = tri_wire >= 0
    g = ac_slot.triangle_geometry(sol)
    a_tri = g["A"][ok]
    k = tri_wire[ok]
    wgt = g["area_mm2"][ok]
    out_j = np.empty(a_tri.shape, dtype=complex)
    omega = 2.0 * math.pi * f
    for i in range(len(wires)):
        m = k == i
        a_mean = np.sum(wgt[m] * a_tri[m].mean(axis=1)) / wgt[m].sum()
        out_j[m] = cur[i] / area[i] - 1j * omega / rho * (a_tri[m] - a_mean)
    return sol["xy"][sol["tri"][ok]], k, out_j


def solve(tasks: list[dict], label: str, workers: int = 4, timeout_s: int = 300) -> list[dict]:
    """Run the tasks on FEMM workers under the project's solver slot (one heavy run at a time, below-normal priority); a task that failed raises."""
    with femm_queue.femm_pool(label, workers=min(workers, len(tasks)), task_timeout_s=timeout_s) as pool:
        res = pool.map(wire_task, tasks)
    bad = [r for r in res if not (isinstance(r, dict) and r.get("status") == "ok")]
    if bad:
        raise RuntimeError(f"{label}: " + "; ".join(str(b.get('error', b)) for b in bad[:3]))
    return res


def plot_wires(ax, res: dict, value, cmap, norm, d_mm: float | None = None, lim: float | None = None, edge: bool = False):
    """Draw every triangle of the wires of a result coloured by value (one number a triangle) on the axes ax; returns the collection (for a colour bar)."""
    from matplotlib.collections import PolyCollection
    pc = PolyCollection(res["tri"], array=value, cmap=cmap, norm=norm, edgecolors="face" if not edge else "none", linewidths=0.3)
    ax.add_collection(pc)
    ax.set_aspect("equal")
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK2, labelsize=8)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    if lim is not None:
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    return pc


def style_axes(ax, xlabel: str, ylabel: str) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
    ax.tick_params(colors=INK2, labelsize=8.5)
    ax.grid(True, color=GRID, linewidth=0.75)
    ax.set_axisbelow(True)
    ax.set_xlabel(xlabel, color=INK2, fontsize=9)
    ax.set_ylabel(ylabel, color=INK2, fontsize=9)
