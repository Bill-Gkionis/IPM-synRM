"""optimization.maps_femm: the FEMM side of the flux maps and the iron-loss map: the tasks a femm_pool worker runs, one per (set of) operating point(s).

What it is, and why
    The optimized motor (design.params(final.json's overrides)) has its geometry fixed. What the torque-speed map needs of FEMM is its behaviour over the
    CURRENT plane: the flux linkages psi_d(i_d, i_q), psi_q(i_d, i_q) (with the saturation of the steel and its cross-coupling), the torque, and what the
    flux density does in the steel (for the iron loss). This file holds the three kinds of question, each a module-level function that
    femm_queue.femm_pool.map() can hand to a worker (a worker is a new Python with a FEMM of its own, which imports this module by name):
      flux_task    the flux linkages and torques of a list of (i_d, i_q) points at a list of rotor positions: the raw rows of the flux maps and of their checks
      frozen_task  one point with FEMM's FROZEN PERMEABILITY: the magnets' flux and the current's flux in the saturated steel of that point, told apart
      loss_task    one point over ONE FULL ELECTRICAL PERIOD (360 degrees, N positions): the flux density waveform of every triangle of the steel, reduced
                   to harmonic amplitudes and one time-domain mean (iron_loss.waveform_terms), for the iron-loss sums
    Nothing here starts FEMM or ends a process: the pool does (R10, one heavy run at a time). Every task writes model and solution files of its own (two FEMMs
    never share a file name) and removes them. A task that fails returns {"status": "error", ...} instead of raising, so one bad point does not end a batch.
    FEMM's M-19 block lists a lamination thickness of 0.635 mm; it enters only FEMM's AC eddy-current solver and changes nothing in these magnetostatic
    numbers (the loss model's thickness is iron_loss.py's: the datasheet's 0.356 mm for the fit, the CAD's 0.3 mm for the headline).
    Driven by maps_run.py; the numbers are analysed by maps_report.py. FEMM 5.0 (D:\\Apps\\femm-qt, build 1144 at the time of writing), the FINE mesh
    (smart mesh on, about 17,000 nodes) for every number, a 2D model 20 mm deep, the 10-turn placeholder winding, currents in A PEAK (the length of the
    amplitude-invariant (i_d, i_q) vector is the phase peak current).

The method, step by step
    The flux linkage and the torque at one rotor position are motorgeo.FemmModel.solve's (the same code as design.score, so the same numbers):
      rotor position theta_e [electrical degrees]; d-q currents -> phase currents by the inverse Park transform with the rotor angle (d is the magnet axis);
      FEMM solves the magnetostatic problem (nonlinear M-19, N45 as a straight line, periodic air-gap element) on the 90 degree sector; the flux linkage
      of each phase is the sector's circuit flux times the 4 sectors; psi_d, psi_q by the Park transform; the torque T_age from the air-gap element's
      integral (the study's torque), T_wst from the weighted stress tensor on the rotor, T_dq = 1.5 p (psi_d i_q - psi_q i_d).
    The frozen-permeability split (frozen_task), because the flux of a saturated machine is NOT the magnets' flux plus L times the current:
      a. solve the loaded point (nonlinear): the steel of every triangle has a permeability that depends on the flux in it;
      b. freeze them: FEMM's "previous solution, frozen permeability" (mi_setprevious(file, 2)) solves a LINEAR problem whose every triangle has the
         permeability (B / H) that the loaded solution had;
      c. in that linear problem superposition holds exactly, so four solves give everything: F0 (no current: the magnets' flux psi_PM in the loaded steel),
         Fd and Fq (a test current on one axis: the columns of the frozen inductance matrix L_fp, with the magnets on and F0 taken away), and Fref (the
         loaded point's own currents: it must give back the nonlinear solution, which proves that FEMM freezes the SECANT permeability);
      d. psi(loaded) = psi_PM,fp + L_fp i must then hold: the check of the whole procedure (maps_report.py prints it).
    A FEMM quirk found on the way, and what is done about it (2026-11-09): FEMM's frozen permeability IGNORES the lamination fill of the steel. The M-19 block has
    a fill of 0.98, and with it the frozen solve with the loaded point's own sources misses the loaded flux linkage by 0.9 % (psi_d at the rated point; 3.8 % if the
    fill is set to 0.9: it goes as 1 - fill); with a fill of 1.0 it returns it to 1e-7. FEMM applies a fill f < 1 as steel and insulation in parallel under the same H:
    B_eff(H) = f B_steel(H) + (1 - f) mu0 H (a fill-1 material with that curve reproduced the fill-0.98 flux linkages and torque to 1e-13). So frozen_task
    solves the EQUIVALENT model (use_effective_steel: fill 1, the effective curve in every steel block): the same machine, to 1e-13 (checked again at each point),
    on which the frozen method closes. The same formula says what the steel's own flux density is, for the loss model: B_steel = (B - (1 - f) mu0 H) / f, which is
    B / f to within (1 - f) mu0 H / f, 0.3 per cent at the most saturated triangle and 0.01 per cent in the middle of the steel: iron_loss.py uses B / f.
    The loss waveform (loss_task): the model is solved at N positions equally spaced over 360 electrical degrees (the 90 degree sector is one electrical period
    of a motor with 4 pole pairs: the field of the sector repeats after it exactly), the currents turning with the rotor. After each solve the .ans file is read
    (femfile.read_ans: first-order triangles, one B = curl A per triangle), and the steel triangles' (Bx, By) are kept. The mesh must be the same at every
    position (only the air-gap element turns); it is checked. A rotor triangle is a fixed piece of the rotor, a stator triangle a fixed piece of the stator:
    each sees, over the period, exactly the waveform its steel sees.
    Regions: rotor steel = triangles whose centroid lies inside the air gap's mean radius; stator steel is split at the radius of the deepest point of the
    slots: inside it the TEETH (tooth tips, tooth bodies and the root fans), outside it the YOKE.

Run
    Not on its own: maps_run.py opens the pool, builds the task lists and stores the results (it waits for the solver slot while another run holds it).

Numbers to expect (FEMM 5.0, fine mesh, the optimized design; the same to the last digit on repeated runs, FEMM 5.0 is deterministic here)
    15,693 nodes and 30,290 triangles, 17,208 of them steel (teeth 7,773, yoke 1,564, rotor 7,871; one 90 degree sector); a flux solve takes 0.5 ... 1 s (15 positions: 9 ... 24 s on a
    worker of 8, more at the heaviest corner of the grid); a frozen point of 6 positions (6 solves each: native, equivalent, four frozen) 16 ... 18 s; a loss point of 72 positions 63 s,
    of 144 positions 126 s (a solve 0.5 s, read_ans 0.15 s, the flux density of 30,290 triangles and the harmonics of 17,208). The regression of flux_task on compare.csv agrees to 1e-14
    relative (the same code path as design.score). maps_run.py's log has the timings of the run.

Self-check (the answers follow).
  1. Why must frozen_task solve the loaded point and the frozen ones in two different model files (and the native model in a third)?
  2. The frozen solve with the loaded point's own currents (Fref) returns the loaded flux linkage. What does it show, and what would it look like if FEMM
     froze the incremental permeability instead?
  3. A rotor triangle's B is read in the model's coordinates at every rotor position. Why is that the waveform the rotor's steel sees?
Answers.
  1. The solver reads the previous solution (the .ans of the loaded point) while it writes its own .ans under the model's name. One file for both would
     overwrite the solution it is reading from. The native model (fill 0.98) is kept apart so that each position can show the equivalent model is the same problem.
  2. That the frozen permeability is the secant B / H: the same sources in a linear problem with those permeabilities reproduce the solution that defined them.
     With the incremental (tangent) permeability the loaded point's own currents would NOT give the loaded flux linkage back.
  3. FEMM turns the rotor by rotating the rotor side of the air-gap element, not the rotor's mesh: the rotor's triangles stay where they were drawn, in the
     rotor's own frame, and the field they carry is the field in that frame.
"""
from __future__ import annotations

import math
import os
import shutil
import time
import traceback
from pathlib import Path

import numpy as np

import femfile
import motorgeo as mg
from optimization import design, iron_loss

STEEL = "M-19 Steel"


# ============================================================================================== helpers (they run inside a worker)
def _files(task: dict) -> tuple[str, str]:
    """The model file and the solution file of this task: names of its own, so that two FEMMs never share one."""
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
    return fem, fem[:-4] + ".ans"


def _cleanup(femm, *files: str) -> None:
    try:
        femm.mi_close()
    except Exception:
        pass
    for path in files:
        try:
            os.remove(path)
        except OSError:
            pass


def _open_model(task: dict, fem: str):
    """Draw the design on the fine mesh, save it as `fem` and open it: (Params, regions, FemmModel). A design motorgeo rejects is an error."""
    import femm
    p = design.params(task["overrides"])
    regions = mg.build(p)
    _, violations = mg.check(p, regions)
    if violations:
        raise RuntimeError(f"{task['tag']}: motorgeo.check rejects this design {task['overrides']}: {'; '.join(violations)}")
    mg.to_femm(p, regions, fem, fine=True)
    femm.mi_close()
    return p, regions, mg.FemmModel(p, fem)


def _error(error: Exception, started: float) -> dict:
    return {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4), "seconds": time.perf_counter() - started}


# ============================================================================================== flux_task
def flux_task(task: dict) -> dict:
    """A list of (i_d, i_q) points at a list of rotor positions.

    task = {"tag", "overrides", "fem_dir", "points": [[i_d, i_q], ...] [A peak], "positions": [theta_e, ...] [electrical degrees]}
    Returns {"status": "ok", "rows": array (n, 9): i_d, i_q, theta_e, psi_d [Wb], psi_q [Wb], T_age, T_wst, T_dq [N.m], psi_retries, "nodes", "seconds"}.
    The model is drawn once and reused for every point and position (FemmModel.solve turns the rotor and sets the currents each time)."""
    import femm
    started = time.perf_counter()
    fem, ans = _files(task)
    try:
        p, regions, model = _open_model(task, fem)
        rows = []
        for i_d, i_q in task["points"]:
            for theta in task["positions"]:
                r = model.solve(float(theta), float(i_d), float(i_q))
                rows.append([i_d, i_q, theta, r["psi_d"], r["psi_q"], r["T_age"], r["T_wst"], r["T_dq"], r["psi_retries"]])
        femm.mi_loadsolution()
        nodes = int(femm.mo_numnodes())
        femm.mo_close()
    except Exception as error:
        return _error(error, started)
    finally:
        _cleanup(femm, fem, ans)
    return {"status": "ok", "kind": "flux", "rows": np.array(rows, float), "nodes": nodes, "seconds": time.perf_counter() - started}


# ============================================================================================== lamfill_task
def lamfill_task(task: dict) -> dict:
    """Does FEMM apply the steel's lamination fill factor (0.98) in this model? One question, a few solves.

    task = {"tag", "overrides", "fem_dir", "fills": [0.98, 1.0, 0.5], "points": [[i_d, i_q], ...], "theta_e"}.
    The loss model needs the steel's own flux density: if FEMM's block average B is fill x B_steel (laminated steel in parallel with its insulation), the steel
    carries B / fill; if FEMM ignores the fill for the lamination type of this model ("not laminated, or laminated in the plane": LamType 0), B is the steel's.
    The test sets the fill of "M-19 Steel" (mi_modifymaterial, property 8) to each value and solves the same points: a fill that is applied moves the flux
    linkage (0.5 moves it a lot), one that is ignored does not. Returns {"rows": array (n_fills x n_points, 7): fill, i_d, i_q, psi_d, psi_q, T_age, T_dq}."""
    import femm
    started = time.perf_counter()
    fem, ans = _files(task)
    try:
        p, regions, model = _open_model(task, fem)
        rows = []
        for fill in task["fills"]:
            femm.mi_modifymaterial(STEEL, 8, float(fill))
            for i_d, i_q in task["points"]:
                r = model.solve(float(task["theta_e"]), float(i_d), float(i_q))
                rows.append([fill, i_d, i_q, r["psi_d"], r["psi_q"], r["T_age"], r["T_dq"]])
    except Exception as error:
        return _error(error, started)
    finally:
        _cleanup(femm, fem, ans)
    return {"status": "ok", "kind": "lamfill", "rows": np.array(rows, float), "seconds": time.perf_counter() - started}


# ============================================================================================== frozen_task
MU0 = 4e-7 * math.pi
FROZEN_COLUMNS = ["theta_e", "ref_psi_d", "ref_psi_q", "F0_psi_d", "F0_psi_q", "Fd_psi_d", "Fd_psi_q", "Fq_psi_d", "Fq_psi_q", "Fref_psi_d", "Fref_psi_q",
                  "T_age_ref", "T_age_Fref", "T_dq_ref", "native_psi_d", "native_psi_q"]


def use_effective_steel(femm, regions, name: str = "M-19 eff") -> None:
    """In the OPEN model: replace the lamination-filled M-19 by a fill-1 steel with FEMM's own effective B-H curve, in every steel block.

    FEMM treats a lamination fill f < 1 ("not laminated, or laminated in the plane") as steel and insulation side by side under the same H, so the block's
    curve is B_eff(H) = f B_steel(H) + (1 - f) mu0 H. A fill-1 material with that curve is the same problem for FEMM: the two models gave the same flux linkages
    and torque to 1e-13 relative at the rated point (checked again at every frozen point, in the rows' native_psi columns). The reason it is needed:
    FEMM's frozen permeability ignores the fill. With f = 0.98 the frozen solve with the loaded point's own sources misses the loaded flux linkage by 0.9 %
    (3.8 % for f = 0.9, so it goes as 1 - f); with the fill-1 equivalent it returns it to 1e-7, as the secant-permeability method must."""
    m = mg.M19_STEEL_FEMM42
    f = m["LamFill"]
    femm.mi_addmaterial(name, m["Mu_x"], m["Mu_y"], m["H_c"], m["J_re"], m["Sigma"], m["d_lam"], m["Phi_h"], 1.0, m["LamType"], m["Phi_hx"], m["Phi_hy"],
                        m["NStrands"], m["WireD"])
    for b, h in m["BH"]:
        femm.mi_addbhpoint(name, f * b + (1.0 - f) * MU0 * h, h)
    for r in regions:
        if r.material == STEEL:
            group = 3 if r.name.startswith("gap_") else 1 if r.side == "rotor" else 2          # motorgeo.to_femm's own groups
            femm.mi_selectlabel(*r.label)
            femm.mi_setblockprop(name, 0, r.mesh, r.circuit, r.magdir, group, r.turns)
            femm.mi_clearselected()


def frozen_task(task: dict) -> dict:
    """One loaded point with FEMM's frozen permeability, at a list of rotor positions.

    task = {"tag", "overrides", "fem_dir", "i_d", "i_q", "positions", "test_current"}: the loaded point [A peak], the positions [electrical degrees] and the test
    current [A] of the two axis solves (any value: the frozen problem is linear; one that makes a flux linkage well above round-off, by default 10 A).
    The model is the design with the lamination fill turned into an equivalent fill-1 steel (use_effective_steel: FEMM's frozen permeability does not handle a
    fill below 1). Per position (see the module docstring, "the frozen-permeability split"): the nonlinear solve of the loaded point on the NATIVE model (the
    check that the equivalent one is the same problem), the same on the equivalent model, then four frozen solves on a copy of it whose "previous solution" is
    that one: F0 (no current), Fd (test current on d), Fq (on q), Fref (the loaded point's currents again).
    Returns {"status": "ok", "rows": array (n_positions, 16), columns FROZEN_COLUMNS, "test_current", "i_d", "i_q", "seconds"}."""
    import femm
    started = time.perf_counter()
    fem_n, ans_n = _files(task)                                # the native model (fill 0.98)
    fem_a = fem_n[:-4] + "_eq.fem"                             # the equivalent model: the loaded point, nonlinear
    ans_a = fem_a[:-4] + ".ans"
    fem_b = fem_n[:-4] + "_frozen.fem"                         # the same with a previous solution: frozen permeability
    ans_b = fem_b[:-4] + ".ans"
    ans_ref = fem_n[:-4] + "_loaded.ans"                       # a copy of the loaded point's solution: the previous solution of the frozen model
    delta = float(task.get("test_current", 10.0))
    try:
        p, regions, model_n = _open_model(task, fem_n)         # fem_n is saved by motorgeo.to_femm and stays as it is
        use_effective_steel(femm, regions)
        femm.mi_saveas(fem_a)                                  # the equivalent model, saved under its own name
        femm.mi_setprevious(ans_ref, 2)                        # the frozen model: the same drawing with a previous solution (type 2: frozen permeability)
        femm.mi_saveas(fem_b)
        femm.mi_close()                                         # closes fem_b; fem_a, saved before the change, has no previous solution
        rows = []
        for theta in task["positions"]:
            theta = float(theta)
            i_d, i_q = float(task["i_d"]), float(task["i_q"])
            native = mg.FemmModel(p, fem_n).solve(theta, i_d, i_q)         # the real model
            femm.mi_close()
            ref = mg.FemmModel(p, fem_a).solve(theta, i_d, i_q, keep_ans=ans_ref)    # the same problem with the fill carried by the curve: the loaded point
            femm.mi_close()
            model_b = mg.FemmModel(p, fem_b)                   # same position, the permeabilities of ans_ref frozen
            f0 = model_b.solve(theta, 0.0, 0.0)
            fd = model_b.solve(theta, delta, 0.0)
            fq = model_b.solve(theta, 0.0, delta)
            fr = model_b.solve(theta, i_d, i_q)
            femm.mi_close()
            rows.append([theta, ref["psi_d"], ref["psi_q"], f0["psi_d"], f0["psi_q"], fd["psi_d"], fd["psi_q"], fq["psi_d"], fq["psi_q"],
                         fr["psi_d"], fr["psi_q"], ref["T_age"], fr["T_age"], ref["T_dq"], native["psi_d"], native["psi_q"]])
    except Exception as error:
        return _error(error, started)
    finally:
        _cleanup(femm, fem_n, ans_n, fem_a, ans_a, fem_b, ans_b, ans_ref)
    return {"status": "ok", "kind": "frozen", "rows": np.array(rows, float), "test_current": delta, "i_d": task["i_d"], "i_q": task["i_q"],
            "seconds": time.perf_counter() - started}


# ============================================================================================== loss_task
def slot_radius(regions) -> float:
    """The radius of the deepest point of the slots [mm]: the stator's teeth end there, the yoke begins (motorgeo.check's "yoke" is r_so minus this)."""
    return max(math.hypot(*xy) for r in regions if r.name.startswith(("opening_", "coil_")) for xy in r.polygon().exterior.coords)


def classify_steel(sol: dict, p, regions) -> dict:
    """Which triangles of a solved .ans are steel, which region each one belongs to, and their areas and centroids. Returns
    {"steel": bool per triangle, "region": int8 per STEEL triangle (0 teeth, 1 yoke, 2 rotor), "area_mm2", "centroid" (per steel triangle), "r_split"}.
    Steel is the material "M-19 Steel" (by the .ans's own block-label table, which the triangles' label numbers index); rotor/stator by the centroid's radius
    against the air gap's mean radius; teeth/yoke by the slots' deepest radius."""
    labels = sol["labels"]
    material = np.array([labels[k]["material"] for k in sol["label"]])
    corners = sol["xy"][sol["tri"]]                                         # (n_tri, 3, 2)
    centroid = corners.mean(axis=1)
    radius = np.hypot(centroid[:, 0], centroid[:, 1])
    d1, d2 = corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]
    area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])
    steel = material == STEEL
    r_gap = 0.5 * (p.r_ro + p.r_si)
    r_split = slot_radius(regions)
    region = np.where(radius[steel] < r_gap, 2, np.where(radius[steel] <= r_split, 0, 1)).astype(np.int8)
    return {"steel": steel, "region": region, "area_mm2": area[steel], "centroid": centroid[steel], "r_split": r_split}


def loss_task(task: dict) -> dict:
    """One operating point over ONE ELECTRICAL PERIOD: the flux density waveform of every steel triangle, reduced for the iron-loss sums.

    task = {"tag", "overrides", "fem_dir", "i_d", "i_q", "n_positions"}: the point [A peak] and N (even): positions theta_e = 360 k / N, k = 0 ... N - 1.
    Returns {"status": "ok", "Ak" (n_steel, N/2) float32 [T], "xdot15" (n_steel,) float32 [T^1.5 per period^1.5] (see iron_loss.waveform_terms), "psi_d", "psi_q",
    "T_age", "T_wst" (N each), "mesh": {"region", "area_mm2", "centroid", "r_split", "nodes", "triangles", "steel_area_mm2" per region, "polygon_area_mm2" per region},
    "seconds"}. The mesh must be the same at every position: it is compared triangle by triangle with the first one's."""
    import femm
    started = time.perf_counter()
    fem, ans = _files(task)
    n = int(task["n_positions"])
    try:
        p, regions, model = _open_model(task, fem)
        positions = [360.0 * k / n for k in range(n)]
        psi_d, psi_q, t_age, t_wst = [], [], [], []
        static = None
        bx = by = None
        for k, theta in enumerate(positions):
            row = model.solve(theta, float(task["i_d"]), float(task["i_q"]))
            sol = femfile.read_ans(ans)
            b_x, b_y = femfile.element_b(sol["xy"], sol["A"], sol["tri"])
            if static is None:
                static = classify_steel(sol, p, regions)
                static["tri"], static["xy"] = sol["tri"], sol["xy"]
                bx = np.empty((n, int(static["steel"].sum())))
                by = np.empty_like(bx)
            elif not (np.array_equal(sol["tri"], static["tri"]) and np.array_equal(sol["xy"], static["xy"])):
                raise RuntimeError(f"the mesh at theta_e = {theta} differs from the one at {positions[0]}: the triangles are not comparable")
            bx[k], by[k] = b_x[static["steel"]], b_y[static["steel"]]
            psi_d.append(row["psi_d"]); psi_q.append(row["psi_q"]); t_age.append(row["T_age"]); t_wst.append(row["T_wst"])
        Ak, xdot15 = iron_loss.waveform_terms(bx, by)
        region = static["region"]
        steel_area = [float(static["area_mm2"][region == r].sum()) for r in range(3)]
        visible = mg.visible(regions)                           # what each region shows once the regions drawn after it (slots, magnets) are painted over it
        poly_area = [float(sum(visible[r.name].area for r in regions if r.material == STEEL and r.side != "rotor")),
                     float(sum(visible[r.name].area for r in regions if r.material == STEEL and r.side == "rotor"))]
        mesh = {"region": region, "area_mm2": static["area_mm2"], "centroid": static["centroid"], "r_split": static["r_split"], "nodes": int(len(static["xy"])),
                "triangles": int(len(static["tri"])), "steel_area_mm2": steel_area, "polygon_area_mm2_stator_rotor": poly_area}
    except Exception as error:
        return _error(error, started)
    finally:
        _cleanup(femm, fem, ans)
    return {"status": "ok", "kind": "loss", "i_d": task["i_d"], "i_q": task["i_q"], "positions": positions, "psi_d": psi_d, "psi_q": psi_q, "T_age": t_age,
            "T_wst": t_wst, "Ak": Ak.astype(np.float32), "xdot15": xdot15.astype(np.float32), "mesh": mesh, "seconds": time.perf_counter() - started}
