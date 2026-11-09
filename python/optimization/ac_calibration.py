"""optimization.ac_calibration: FEMM's time-harmonic solver against the cases whose answer is known, before it is trusted on the slot (item 1 of the copper-loss work).

What it is, and why
    The slot study (optimization/ac_slot.py) draws hundreds of round strands and lets FEMM find the eddy currents in each. Whether those numbers mean anything depends on three things this
    module tests one at a time, each on a case with an exact answer (optimization/ac_theory.py), at the SAME discretisation the slot model uses (strand outlines of 10 degree segments,
    mesh size d / 6, the area of the polygon made equal to the circle's) and at a much finer one, so that the discretisation's own error is known:
      (a) SKIN EFFECT: one round wire alone, carrying a current (a series circuit of 1 A), at d / delta from 0.1 to 5. FEMM's loss over the DC loss of the same wire (made from FEMM's own
          area of the wire) against the exact Bessel-function result R_ac / R_dc = Re[(z / 2) I_0(z) / I_1(z)], z = (1 + j) d / (2 delta).
      (b) PROXIMITY EFFECT: one round wire with no net current (a series circuit of 0 A) in a uniform transverse field (the outer boundary of the model prescribes A = -B x, so far from the
          wire the field is B along y). FEMM's loss against the exact result P' = (2 pi a B^2 / (mu0^2 sigma)) Im[gamma conj(I_1(gamma a) / I_0(gamma a))], whose low-frequency limit is the
          winding study's pi sigma w^2 B^2 d^4 / 128 (ac_theory.py derives both). The ratio to the low-frequency formula is what 'the formula holds' means: 0.993 at d / delta = 1,
          0.898 at 2, 0.638 at 3.
      (c) THE PARALLEL CIRCUIT: two identical wires side by side across a uniform field, joined in ONE parallel circuit of 0 A total (FEMM gives them a common voltage and lets the currents
          divide), against the same two wires in separate series circuits of 0 A each. The circulating current between them is jw sigma S B s / 2, so the extra loss is w^2 sigma S B^2 s^2 / 4
          (two_wire_circulating) and the pair's total is that plus twice the single wire's proximity loss. This is the case the strands of one turn in parallel are made of.
    The wire is d = 0.63 mm (the thick strand of the study: its d / delta runs from 0.1 to 5 between 145 Hz and 361 kHz) in copper at 100 degC, 20 mm deep, FEMM 5.0.20261101; every result
    carries its mesh (the number of nodes).

Pass criteria (set BEFORE the run, from the discretisation: the polygon of a 10 degree outline is made area-exact, so what is left is the mesh)
    d / delta from 0.1 to 3 (the range of the winding): skin |FEMM / exact - 1| <= 1e-3 and proximity <= 5e-3 at the slot's settings; <= 3e-4 and <= 1e-3 at the fine ones; the parallel pair's
    total loss within 5e-3 of the analytic. d / delta = 5 is reported and not a criterion (beyond the range of the winding). A failure is reported, not repaired by widening a tolerance.

Run (from the Python to FEMM folder; one pool session of a few seconds, waiting for the solver slot if another run holds it)
    uv run python -m optimization.ac_calibration
Writes out/optimization/copper_ac/calibration/calibration.json (every number), calibration.csv and fig_calibration.png. Exit status 0 if every criterion passes (1 on this PC: the proximity
criteria fail above d / delta 1.5 at the slot's settings and above 1 at the fine ones, see below).
    uv run python -m optimization.ac_calibration --redraw        only redraws the figure from calibration.json (no FEMM)

Numbers to expect (FEMM 5.0.20261101; the run takes 7.9 s of pool time, 4,109 / 8,213 nodes for the pair, 3,894 / 5,868 for the skin case at the slot's / the fine settings)
    skin       |FEMM / exact - 1| at most 8.4e-4 at the slot's settings and 8.1e-5 at the fine ones for d / delta 0.1 ... 3: PASS (1.2e-2 at d / delta 5, not a criterion).
    proximity  FEMM's loss over the exact one, minus 1 (FEMM above), slot's settings: 5e-6 at d / delta 0.1, 1.4e-4 at 0.7, 5.6e-4 at 1, 2.7e-3 at 1.5, 8.0e-3 at 2, 1.7e-2 at 2.5, 2.9e-2 at 3;
               fine settings: 2.6e-4 at 1, 1.3e-3 at 1.5, 3.7e-3 at 2, 7.9e-3 at 2.5, 1.3e-2 at 3. The criteria (5e-3 and 1e-3) FAIL above d / delta 1.5 (slot's settings) and 1 (fine): reported,
               the tolerances not widened. The winding's own range is d / delta <= 1.2 at the 20 kHz carrier for its thickest strand (0.63 mm; 0.59 for the 0.315 mm strand), where the slot-settings error is about 1e-3 or less (5.6e-4 at 1).
    pair       the parallel pair's total loss within 4.9e-5 of the analytic at 200 Hz and 1.2e-3 at 1,000 Hz at both settings (the circulating current 17.285 A against 17.286 A, and 86.38 A
               against 86.43 A): PASS. The 1.2e-3 at 1 kHz is the analytic formula's own low-frequency approximation (it is the same at both settings).
    why the proximity criteria fail: the convergence stage (22 solves) shows the error is the mesh of the AIR round the wire. At d / delta 2 (box 40 radii) it is 2.7e-3 with the automatic mesh, 2.1e-3
               with 0.2 mm elements in the air, 1.3e-3 at 0.1 mm and 7.0e-4 at 0.05 mm; at d / delta 3 it is 9.6e-3, 7.4e-3, 4.7e-3 and 2.5e-3. The box radius matters less: 40 against 100 radii at a
               fixed 0.2 mm air mesh differs by 2e-5 (d / delta 2) and 8e-5 (d / delta 3). The slot model's strands are surrounded by other strands, which is not this case: its own mesh convergence
               (ac_slot_study.py, sensitivity stage) is 0.05 % at 20 kHz.

Self-check (the answers follow).
  1. At d / delta = 0.1 FEMM's R_ac / R_dc is 1.0000001 and the exact one 1.0000001. What does the check test at that point, and what does it not?
  2. The skin-effect result does not depend on the radius of the model's outer boundary. Why, and why does the proximity result?
  3. Why is the loss of two parallel wires in the field so much more than that of the same two wires in separate circuits?
Answers.
  1. It tests that FEMM's solve, its polygon of the wire and the DC reference are consistent (the ratio is 1 to 1e-7 where nothing happens). It says nothing about the eddy currents,
     which are 1e-7 of the loss: the checks at d / delta of 1 and above are the ones that do.
  2. A wire alone carries a current that is the same all round, so the field inside depends only on the current enclosed at each radius (Ampere's law), whatever lies beyond. A uniform
     external field must be brought in from the boundary, and the boundary is a finite distance: the wire's own response (a dipole, falling as 1 / r^2) has to be held at the boundary,
     which is why the model is 100 radii wide for the proximity case and 40 for the skin case (the convergence stage: from 40 to 100 radii the result moves by less than 1e-4).
  3. In separate circuits each wire's eddy currents close within the wire. In a parallel circuit a current can also run along one wire and come back along the other, driven by the
     difference of A at their two positions (jw (A_1 - A_2) over the resistance of the loop), and its loop is the whole separation wide: its loss goes as s^2, the eddy loss as d^2.
"""
from __future__ import annotations

import csv
import json
import math
import time
import traceback
from pathlib import Path

import numpy as np

import femm_queue
import paths
from optimization import ac_slot, ac_theory

OUT = paths.path("out") / "optimization" / "copper_ac" / "calibration"
D_MM = 0.63
DEPTH_MM = 20.0
T_COPPER = 100.0
RATIOS = (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0)           # d / delta
RANGE_MAX = 3.0                                                              # the criteria apply up to here
SETTINGS = {"slot": {"arc_deg": ac_slot.DEFAULT_OPTIONS["arc_deg"], "mesh_cu_over_d": ac_slot.DEFAULT_OPTIONS["mesh_cu_over_d"], "area_correct": ac_slot.DEFAULT_OPTIONS["area_correct"]},
            "fine": {"arc_deg": 1.0, "mesh_cu_over_d": 1.0 / 20.0, "area_correct": True}}
CRITERIA = {"slot": {"skin": 1e-3, "prox": 5e-3}, "fine": {"skin": 3e-4, "prox": 1e-3}, "pair": 5e-3}
PAIR_S_MM = 2.0                                                              # the two wires' centre distance
PAIR_FREQS = (200.0, 1000.0)
B_PROX = 1.0                                                                 # T peak: the models are linear, so the value is only a scale


def frequency_for(ratio: float, d_mm: float = D_MM, t_c: float = T_COPPER) -> float:
    """The frequency at which a wire of d_mm has d / delta = ratio in copper at t_c: f = rho / (pi mu0 delta^2), delta = d / ratio."""
    rho = float(ac_theory.resistivity(t_c))
    return rho / (math.pi * ac_theory.MU0 * (d_mm * 1e-3 / ratio) ** 2)


# ---------------------------------------------------------------------------------------------- the pool's tasks
def _draw_wire(femm, x: float, y: float, d_mm: float, settings: dict, circuit: str) -> None:
    r = d_mm / 2.0
    if settings["area_correct"]:
        n_seg = 2 * math.ceil(180.0 / settings["arc_deg"])
        r *= math.sqrt(2.0 * math.pi / (n_seg * math.sin(2.0 * math.pi / n_seg)))
    femm.mi_addnode(x + r, y)
    femm.mi_addnode(x - r, y)
    femm.mi_addarc(x + r, y, x - r, y, 180, settings["arc_deg"])
    femm.mi_addarc(x - r, y, x + r, y, 180, settings["arc_deg"])
    femm.mi_addblocklabel(x, y)
    femm.mi_selectlabel(x, y)
    femm.mi_setblockprop("Cu", 0, d_mm * settings["mesh_cu_over_d"], circuit, 0, 0, 1)
    femm.mi_clearselected()


def run_task(task: dict) -> dict:
    """One pool task: a model of kind 'skin', 'prox' or 'pair' (parallel or series), built, solved and read in this worker's FEMM; the files are removed. Errors come back as
    {'status': 'error'} (one bad task must not end a run)."""
    import femm
    import os
    started = time.perf_counter()
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
    ans = fem[:-4] + ".ans"
    try:
        d, s, settings = task["d_mm"], task.get("s_mm", 0.0), task["settings"]
        rho = float(ac_theory.resistivity(task["t_c"]))
        a = d / 2.0
        kind = task["kind"]
        radius = task.get("radius_over_a", 40.0 if kind == "skin" else 100.0) * a if kind != "pair" else 60.0
        femm.newdocument(0)
        femm.mi_probdef(task["f"], "millimeters", "planar", 1e-10, DEPTH_MM, 30)
        femm.mi_getmaterial("Air")
        femm.mi_addmaterial("Cu", 1, 1, 0, 0, 1.0 / rho / 1e6, 0, 0, 1, 0, 0, 0, 0, 0)
        if kind == "skin":
            femm.mi_addboundprop("OUT", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        else:
            femm.mi_addboundprop("OUT", 0, -task["B"] * 1e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0)       # A = A1 x, x in mm: B_y = -A1 * 1e3 tesla
        femm.mi_addnode(radius, 0)
        femm.mi_addnode(-radius, 0)
        femm.mi_addarc(radius, 0, -radius, 0, 180, 5.0)
        femm.mi_addarc(-radius, 0, radius, 0, 180, 5.0)
        femm.mi_selectarcsegment(0, radius)
        femm.mi_selectarcsegment(0, -radius)
        femm.mi_setarcsegmentprop(5.0, "OUT", 0, 0)
        femm.mi_clearselected()
        if kind == "pair":
            centres = [(-s / 2.0, 0.0), (s / 2.0, 0.0)]
            if task["circuit"] == "parallel":
                femm.mi_addcircprop("w", 0.0, 0)
                circ = ["w", "w"]
            else:
                femm.mi_addcircprop("w1", 0.0, 1)
                femm.mi_addcircprop("w2", 0.0, 1)
                circ = ["w1", "w2"]
        else:
            centres = [(0.0, 0.0)]
            femm.mi_addcircprop("w", 1.0 if kind == "skin" else 0.0, 1)
            circ = ["w"]
        for (x, y), c in zip(centres, circ):
            _draw_wire(femm, x, y, d, settings, c)
        femm.mi_addblocklabel(0.0, radius / 2.0)
        femm.mi_selectlabel(0.0, radius / 2.0)
        if task.get("air_mesh_mm"):                                      # a fixed mesh size in the air (the convergence test of the exterior); by default FEMM chooses
            femm.mi_setblockprop("Air", 0, task["air_mesh_mm"], "", 0, 0, 0)
        else:
            femm.mi_setblockprop("Air", 1, 0, "", 0, 0, 0)
        femm.mi_clearselected()
        femm.mi_saveas(fem)
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        loss, cur, area = [], [], []
        for (x, y) in centres:
            femm.mo_selectblock(x, y)
            loss.append(float(np.real(femm.mo_blockintegral(4))))
            cur.append(complex(femm.mo_blockintegral(7)))
            area.append(float(np.real(femm.mo_blockintegral(5))))
            femm.mo_clearblock()
        nodes = int(femm.mo_numnodes())
        femm.mo_close()
    except Exception as error:
        return {"status": "error", "tag": task["tag"], "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4), "seconds": time.perf_counter() - started}
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
    return {"status": "ok", "tag": task["tag"], "loss_W": loss, "current_A": cur, "area_m2": area, "nodes": nodes, "seconds": time.perf_counter() - started}


def make_tasks(fem_dir: str) -> list[dict]:
    tasks = []
    for name, st in SETTINGS.items():
        for ratio in RATIOS:
            f = frequency_for(ratio)
            tasks.append({"kind": "skin", "tag": f"cal_skin_{name}_{ratio:g}".replace(".", "p"), "fem_dir": fem_dir, "d_mm": D_MM, "f": f, "t_c": T_COPPER, "settings": st, "ratio": ratio, "set": name})
            tasks.append({"kind": "prox", "tag": f"cal_prox_{name}_{ratio:g}".replace(".", "p"), "fem_dir": fem_dir, "d_mm": D_MM, "f": f, "t_c": T_COPPER, "settings": st, "ratio": ratio,
                          "set": name, "B": B_PROX})
        for f in PAIR_FREQS:
            for circuit in ("series", "parallel"):
                tasks.append({"kind": "pair", "tag": f"cal_pair_{name}_{circuit}_{f:g}", "fem_dir": fem_dir, "d_mm": D_MM, "s_mm": PAIR_S_MM, "f": f, "t_c": T_COPPER, "settings": st,
                              "set": name, "B": B_PROX, "circuit": circuit})
    return tasks


CONVERGENCE = [(ratio, rad, air) for ratio in (2.0, 3.0) for rad in (20.0, 40.0, 100.0) for air in (None, 0.2, 0.1, 0.05, 0.025)
               if not (rad == 100.0 and air is not None and air < 0.2) and not (rad == 40.0 and air is not None and air < 0.05)]       # (d / delta, model radius over the wire's, air mesh [mm]: None = FEMM's own)


def make_convergence_tasks(fem_dir: str) -> list[dict]:
    """The proximity case at d / delta of 2 and 3 with the wire meshed finely (d / 40, 2 degree outline) and the exterior's mesh and size varied: it shows where the remaining proximity error of
    the main run comes from (the exterior: the field of the wire's own eddy currents, a dipole falling as 1 / r^2, is badly resolved by the mesh FEMM chooses)."""
    tasks = []
    for ratio, rad, air in CONVERGENCE:
        tasks.append({"kind": "prox", "tag": f"cal_conv_{ratio:g}_{rad:g}_{air}".replace(".", "p"), "fem_dir": fem_dir, "d_mm": D_MM, "f": frequency_for(ratio), "t_c": T_COPPER,
                      "settings": {"arc_deg": 2.0, "mesh_cu_over_d": 1.0 / 40.0, "area_correct": True}, "ratio": ratio, "set": "conv", "B": B_PROX, "radius_over_a": rad, "air_mesh_mm": air})
    return tasks


def analyse_convergence(tasks: list[dict], results: list[dict]) -> list[dict]:
    rho = float(ac_theory.resistivity(T_COPPER))
    L, d = DEPTH_MM * 1e-3, D_MM * 1e-3
    out = []
    for t, r in zip(tasks, results):
        femm_val = r["loss_W"][0] / (float(ac_theory.proximity_lowfreq(d, t["f"], t["B"], rho)) * L)
        out.append({"d_over_delta": t["ratio"], "radius_over_a": t["radius_over_a"], "air_mesh_mm": t["air_mesh_mm"], "nodes": r["nodes"], "seconds": r["seconds"],
                    "error": femm_val / float(ac_theory.proximity_ratio(d, t["f"], rho)) - 1.0})
    return out


# ---------------------------------------------------------------------------------------------- the comparison (no FEMM)
def analyse(tasks: list[dict], results: list[dict]) -> dict:
    rho = float(ac_theory.resistivity(T_COPPER))
    L = DEPTH_MM * 1e-3
    d = D_MM * 1e-3
    rows, pairs = [], []
    by = {t["tag"]: r for t, r in zip(tasks, results)}
    for t, r in zip(tasks, results):
        if r.get("status") != "ok":
            continue
        if t["kind"] == "skin":
            r_dc = rho * L / r["area_m2"][0]
            femm_val = r["loss_W"][0] / (0.5 * 1.0 ** 2 * r_dc)
            exact = float(ac_theory.skin_ratio(d, t["f"], rho))
            rows.append({"case": "skin", "set": t["set"], "d_over_delta": t["ratio"], "f_Hz": t["f"], "femm": femm_val, "exact": exact, "error": femm_val / exact - 1.0, "nodes": r["nodes"]})
        elif t["kind"] == "prox":
            femm_val = r["loss_W"][0] / (float(ac_theory.proximity_lowfreq(d, t["f"], t["B"], rho)) * L)           # over the low-frequency formula
            exact = float(ac_theory.proximity_ratio(d, t["f"], rho))
            rows.append({"case": "prox", "set": t["set"], "d_over_delta": t["ratio"], "f_Hz": t["f"], "femm": femm_val, "exact": exact, "error": femm_val / exact - 1.0, "nodes": r["nodes"],
                         "current_leak_A": abs(r["current_A"][0])})
    for name in SETTINGS:
        for f in PAIR_FREQS:
            ser, par = by.get(f"cal_pair_{name}_series_{f:g}"), by.get(f"cal_pair_{name}_parallel_{f:g}")
            if not ser or not par or ser.get("status") != "ok" or par.get("status") != "ok":
                continue
            p_single = float(ac_theory.proximity_exact(d, f, B_PROX, rho)) * L
            p_circ = float(ac_theory.two_wire_circulating(d, PAIR_S_MM * 1e-3, f, B_PROX, rho)) * L
            i_circ_exact = 2.0 * math.pi * f * (1.0 / rho) * (math.pi * d ** 2 / 4.0) * B_PROX * PAIR_S_MM * 1e-3 / 2.0
            pairs.append({"set": name, "f_Hz": f, "series_total_W": sum(ser["loss_W"]), "series_exact_W": 2.0 * p_single, "series_error": sum(ser["loss_W"]) / (2.0 * p_single) - 1.0,
                          "parallel_total_W": sum(par["loss_W"]), "parallel_exact_W": 2.0 * p_single + p_circ, "parallel_error": sum(par["loss_W"]) / (2.0 * p_single + p_circ) - 1.0,
                          "circulating_femm_W": sum(par["loss_W"]) - sum(ser["loss_W"]), "circulating_exact_W": p_circ,
                          "i_circ_femm_A": abs(par["current_A"][0]), "i_circ_exact_A": i_circ_exact, "nodes": par["nodes"]})
    return {"rows": rows, "pairs": pairs}


def criteria_table(an: dict) -> tuple[list[dict], bool]:
    out, ok_all = [], True
    for name in SETTINGS:
        for case in ("skin", "prox"):
            sel = [r for r in an["rows"] if r["set"] == name and r["case"] == case and r["d_over_delta"] <= RANGE_MAX]
            worst = max(abs(r["error"]) for r in sel)
            ok = worst <= CRITERIA[name][case]
            ok_all &= ok
            out.append({"set": name, "case": case, "worst_error": worst, "criterion": CRITERIA[name][case], "ok": ok, "n": len(sel)})
    for pr in an["pairs"]:
        worst = max(abs(pr["parallel_error"]), abs(pr["series_error"]))
        ok = worst <= CRITERIA["pair"]
        ok_all &= ok
        out.append({"set": pr["set"], "case": f"pair at {pr['f_Hz']:g} Hz", "worst_error": worst, "criterion": CRITERIA["pair"], "ok": ok, "n": 2})
    return out, ok_all


def figure(an: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    surface, ink, ink2, grid, axis = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    rho = float(ac_theory.resistivity(T_COPPER))
    d = D_MM * 1e-3
    x = np.logspace(-1.05, 0.75, 200)
    fq = np.array([frequency_for(v) for v in x])
    fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.2), facecolor=surface)
    series = {"slot": (orange, "o", "slot-model settings (10 deg outline, mesh d/6)"), "fine": (blue, "s", "fine settings (1 deg outline, mesh d/20)")}
    ax = axes[0, 0]
    ax.plot(x, ac_theory.skin_ratio(d, fq, rho), color=gray, lw=2.0, label="exact (Bessel functions)")
    for name, (col, mk, lab) in series.items():
        sel = [r for r in an["rows"] if r["case"] == "skin" and r["set"] == name]
        ax.plot([r["d_over_delta"] for r in sel], [r["femm"] for r in sel], mk, color=col, ms=6, mec=surface, mew=1.5, label=f"FEMM, {lab}")
    ax.set_xscale("log")
    ax.set_xlabel("d / delta  (wire diameter over skin depth)", fontsize=8.5, color=ink2)
    ax.set_ylabel("R_ac / R_dc of one wire alone", fontsize=8.5, color=ink2)
    ax.set_title("Skin effect: one wire alone", fontsize=10, color=ink, loc="left", fontweight="bold")
    ax = axes[0, 1]
    ax.plot(x, ac_theory.proximity_ratio(d, fq, rho), color=gray, lw=2.0, label="exact (Bessel functions)")
    for name, (col, mk, lab) in series.items():
        sel = [r for r in an["rows"] if r["case"] == "prox" and r["set"] == name]
        ax.plot([r["d_over_delta"] for r in sel], [r["femm"] for r in sel], mk, color=col, ms=6, mec=surface, mew=1.5, label=f"FEMM, {lab}")
    ax.set_xscale("log")
    ax.set_xlabel("d / delta", fontsize=8.5, color=ink2)
    ax.set_ylabel("proximity loss over the low-frequency formula", fontsize=8.5, color=ink2)
    ax.set_title("Proximity effect: no net current, uniform transverse field", fontsize=10, color=ink, loc="left", fontweight="bold")
    for k, (case, title) in enumerate((("skin", "Skin effect: FEMM over exact, minus 1"), ("prox", "Proximity effect: FEMM over exact, minus 1"))):
        ax = axes[1, k]
        for name, (col, mk, lab) in series.items():
            sel = [r for r in an["rows"] if r["case"] == case and r["set"] == name]
            ax.plot([r["d_over_delta"] for r in sel], [100 * r["error"] for r in sel], "-" + mk, color=col, lw=1.5, ms=6, mec=surface, mew=1.5, label=lab)
        ax.axhline(0, color=axis, lw=1.0)
        ax.axvline(RANGE_MAX, color=axis, lw=1.0, ls=":")
        ax.text(RANGE_MAX / 1.04, ax.get_ylim()[0] + 0.62 * (ax.get_ylim()[1] - ax.get_ylim()[0]), "criteria apply\nup to here", fontsize=7.5, color=ink2, va="bottom", ha="right")
        ax.set_xscale("log")
        ax.set_xlabel("d / delta", fontsize=8.5, color=ink2)
        ax.set_ylabel("error [%]", fontsize=8.5, color=ink2)
        ax.set_title(title, fontsize=10, color=ink, loc="left", fontweight="bold")
    for ax in axes.ravel():
        ax.set_facecolor(surface)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(axis)
        ax.tick_params(colors=ink2, labelsize=8.5)
        ax.grid(True, color=grid, lw=0.75)
        ax.set_axisbelow(True)
        leg = ax.legend(frameon=False, fontsize=7.8, loc="best")
        for t in leg.get_texts():
            t.set_color(ink2)
    fig.text(0.012, 0.008, f"A copper wire of {D_MM} mm at {T_COPPER:g} degC, {DEPTH_MM:g} mm deep, FEMM 5.0 time-harmonic; skin: series circuit of 1 A, proximity: 0 A in a uniform field "
             "(prescribed A on a boundary 100 radii out). Source: optimization/ac_calibration.py.", fontsize=7.5, color=ink2)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(path, dpi=130, facecolor=surface)
    plt.close(fig)


def main(argv=None) -> int:
    import argparse
    ap = argparse.ArgumentParser(description="FEMM's time-harmonic solver against the exact cases (skin, proximity, parallel pair); --redraw only redraws fig_calibration.png from calibration.json")
    ap.add_argument("--redraw", action="store_true", help="redraw the figure from the saved calibration.json (no FEMM)")
    args = ap.parse_args(argv)
    if args.redraw:
        saved = json.loads((OUT / "calibration.json").read_text(encoding="utf-8"))
        figure({"rows": saved["rows"], "pairs": saved["pairs"]}, OUT / "fig_calibration.png")
        print(f"redrawn {OUT / 'fig_calibration.png'}")
        return 0
    import femm_session
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = make_tasks(str(OUT / "fem"))
    conv_tasks = make_convergence_tasks(str(OUT / "fem"))
    t0 = time.perf_counter()
    with femm_queue.femm_pool("copper ac: calibration", workers=8, task_timeout_s=400) as pool:
        waited = time.perf_counter() - t0
        t1 = time.perf_counter()
        results_all = pool.map(run_task, tasks + conv_tasks)
        run_s = time.perf_counter() - t1
    results, conv_results = results_all[:len(tasks)], results_all[len(tasks):]
    bad = [(t["tag"], r) for t, r in zip(tasks + conv_tasks, results_all) if not (isinstance(r, dict) and r.get("status") == "ok")]
    if bad:
        raise RuntimeError("calibration tasks failed: " + "; ".join(f"{tag}: {r.get('error', r) if isinstance(r, dict) else r}" for tag, r in bad))
    conv = analyse_convergence(conv_tasks, conv_results)
    an = analyse(tasks, results)
    crit, ok_all = criteria_table(an)
    print(f"{len(tasks)} FEMM solves in {run_s:.1f} s ({waited:.1f} s waiting for the slot), FEMM {femm_session.femm_version()}\n")
    print("skin effect: R_ac / R_dc of one wire alone, FEMM against the Bessel solution")
    print("   d/delta   f [Hz]      exact     FEMM(slot)  err(slot)    FEMM(fine)  err(fine)    nodes slot / fine")
    for ratio in RATIOS:
        a = next(r for r in an["rows"] if r["case"] == "skin" and r["set"] == "slot" and r["d_over_delta"] == ratio)
        b = next(r for r in an["rows"] if r["case"] == "skin" and r["set"] == "fine" and r["d_over_delta"] == ratio)
        print(f"  {ratio:6.1f} {a['f_Hz']:10.1f}  {a['exact']:10.7f}  {a['femm']:10.7f} {a['error']:+9.2e}   {b['femm']:10.7f} {b['error']:+9.2e}   {a['nodes']:6d} / {b['nodes']:6d}")
    print("\nproximity effect: loss of one wire in a uniform field over the low-frequency formula pi sigma w^2 B^2 d^4 / 128")
    print("   d/delta   exact       FEMM(slot)  err(slot)    FEMM(fine)  err(fine)")
    for ratio in RATIOS:
        a = next(r for r in an["rows"] if r["case"] == "prox" and r["set"] == "slot" and r["d_over_delta"] == ratio)
        b = next(r for r in an["rows"] if r["case"] == "prox" and r["set"] == "fine" and r["d_over_delta"] == ratio)
        print(f"  {ratio:6.1f}   {a['exact']:10.7f}  {a['femm']:10.7f} {a['error']:+9.2e}   {b['femm']:10.7f} {b['error']:+9.2e}")
    print("\nthe parallel circuit: two wires 2 mm apart across a uniform field of 1 T (peak)")
    for pr in an["pairs"]:
        print(f"  {pr['set']:5s} {pr['f_Hz']:6g} Hz  separate circuits {pr['series_total_W']:.6g} W (exact {pr['series_exact_W']:.6g}, {pr['series_error']:+.1e});  one parallel circuit "
              f"{pr['parallel_total_W']:.6g} W (exact {pr['parallel_exact_W']:.6g}, {pr['parallel_error']:+.1e});  circulating current {pr['i_circ_femm_A']:.4f} A (exact {pr['i_circ_exact_A']:.4f})")
    print("\nwhere the remaining proximity error comes from: the exterior's mesh (the wire meshed at d/40 with a 2 degree outline; model radius over the wire's, air mesh, FEMM / exact - 1)")
    print("   d/delta   R/a   air mesh [mm]   nodes   FEMM / exact - 1")
    for c in conv:
        print(f"   {c['d_over_delta']:6.1f}  {c['radius_over_a']:5.0f}   {str(c['air_mesh_mm']) if c['air_mesh_mm'] else 'FEMM auto':>12s}   {c['nodes']:7d}   {c['error']:+.3e}")
    print("\ncriteria (set before the run):")
    for c in crit:
        print(f"  {'ok  ' if c['ok'] else 'FAIL'} {c['set']:5s} {c['case']:18s} worst error {c['worst_error']:.2e} against {c['criterion']:.0e} ({c['n']} points)")
    (OUT / "calibration.json").write_text(json.dumps({"femm_version": femm_session.femm_version(), "d_mm": D_MM, "depth_mm": DEPTH_MM, "t_copper_C": T_COPPER, "settings": SETTINGS,
                                                      "criteria": crit, "rows": an["rows"], "pairs": an["pairs"], "convergence": conv, "seconds": run_s}, indent=1), encoding="utf-8")
    with open(OUT / "calibration.csv", "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=["case", "set", "d_over_delta", "f_Hz", "femm", "exact", "error", "nodes"], extrasaction="ignore")
        w.writeheader()
        w.writerows(an["rows"])
    figure(an, OUT / "fig_calibration.png")
    print("\nall criteria pass" if ok_all else "\nSOME CRITERIA FAILED")
    return 0 if ok_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
