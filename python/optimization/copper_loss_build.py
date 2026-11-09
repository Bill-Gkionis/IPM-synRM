"""optimization.copper_loss_build: packages the studies' outputs into the tables of copper_loss.py (optimization/copper_loss_tables.json).

What it is, and why
    copper_loss.py is fast because its physics is closed-form and its data is a few small tables. The data comes from FEMM runs and Python studies that are slow and live in
    out/optimization/copper_ac/ (git-ignored, on this PC): the machine's field (ac_machine.py), the slot's AC resistance (ac_slot_study.py), the PWM ripple (ac_pwm.py). This module reads
    their outputs once and writes the one JSON file that sits next to copper_loss.py, so that the optimizer needs neither FEMM nor the studies' folders. The file records where each
    table came from.

What goes in (each table, its source and what the fast model does with it)
    geometry   the optimized stator's coil layer (36.218 mm^2), the radii of the two layers' centroids (50.08 and 56.30 mm), the chord between the two coil sides of a coil (36.85 mm for the pitch of
               4 slots) and the first study's span at the mid radius: motorgeo, via ac_slot.layer_polygons. For the end turns.
    machine    the fits of G (the machine's field at the strands) and of H (the spread of the vector potential) over the seven harmonic groups (ac_machine.fit_split: the no-load part plus
               AT^2 times a bilinear total times linear shares), the quadratic fits of the flux linkages and the torque, with the range of the points; the inductances of the 2-D model
               at the 120 series turns of its 10 turns a coil side (L_d from the fit's own flux linkages, held against the brief's 1.390 mH). ac_machine.py.
    pwm        h(m), the split of the ripple's power over the carrier's harmonics. ac_pwm.py.
    slot       the slot constant K of the static field for layers in phase (15.95 mm^2, the mean of 12 layouts, spread 3.4 %) and 60 degrees apart (13.01 mm^2), and c_iron(f): FEMM's extra
               loss over the static law as a function of frequency (the steel's screening); the circulating-current model's constants. ac_slot_study.py (sweep, val and circ stages).

Run
    uv run python -m optimization.copper_loss --build          (no FEMM, 0.1 s; needs the studies' outputs)
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

import motorgeo as mg
import paths
from optimization import ac_slot, copper_loss, design

AC = paths.path("out") / "optimization" / "copper_ac"
L_D_BRIEF, L_Q_BRIEF = 1.390e-3, 2.112e-3          # H at 10 turns a coil side, the brief's fine-mesh values at the rated point (checked against the fit below)


def geometry() -> dict:
    ov = json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]
    p = design.params(ov)
    polys = ac_slot.layer_polygons(p)
    r0, r1 = polys[0].centroid.y, polys[1].centroid.y
    span = math.radians(p.coil_pitch * 360.0 / p.slots)
    from optimization import winding
    first = winding.end_turn_geometry(p)                 # the first study's stand-in: the radius of the centroid of the two layers together, and the span at it
    return {"A_layer_mm2": float(polys[0].area), "r_layer_mm": [float(r0), float(r1)], "coil_pitch_slots": p.coil_pitch, "span_deg": math.degrees(span),
            "chord_mm": float(math.sqrt(r0 ** 2 + r1 ** 2 - 2 * r0 * r1 * math.cos(span))), "span_mid_mm": float(first["span_mm"]), "r_mid_mm": float(first["r_mean_mm"]),
            "stack_ref_mm": p.stack, "slots": p.slots, "pole_pairs": p.pole_pairs}


def machine_block() -> dict:
    tab = json.loads((AC / "machine" / "tables.json").read_text(encoding="utf-8"))
    summ = json.loads((AC / "machine" / "summary.json").read_text(encoding="utf-8"))
    m = {k: tab[k] for k in ("i_ref_turns", "at_rated", "groups", "ranges", "G", "H", "psi_T_coef", "points")}
    pts = summ["points"]
    stats = {st["name"]: st for st in summ["stats"]}
    psi_pm = stats["no load"]["psi_d"]
    rated = stats["MTPA 1 x rated"]
    l_d = (rated["psi_d"] - psi_pm) / rated["i_d"]
    l_q = rated["psi_q"] / rated["i_q"]
    m["psi_pm_ref_Wb"], m["L_d_fit_rated_H"], m["L_q_fit_rated_H"] = psi_pm, l_d, l_q
    m["L_d_ref_H"], m["L_q_ref_H"] = L_D_BRIEF, L_Q_BRIEF
    m["note"] = ("fits are per 120 series turns (10 turns a coil side, one path) and a 20 mm stack; the field's G_k depends on the ampere-turns only; L_d_ref_H is the brief's value "
                 f"(the fit's secant at the rated point is {1e3 * l_d:.4f} mH, L_q {1e3 * l_q:.4f} mH)")
    return m


def pwm_block() -> dict:
    t = json.loads((AC / "pwm" / "ripple_table.json").read_text(encoding="utf-8"))
    return {k: t[k] for k in ("m", "h_rms", "pp_max", "frac", "tail", "n_harm")}


def slot_block() -> dict:
    from optimization import ac_theory
    sweep = json.loads((AC / "slot" / "sweep.json").read_text(encoding="utf-8"))
    rows = sweep["rows"]
    rho = float(ac_theory.resistivity(100.0))
    low = [r for r in rows if r["f_Hz"] <= 1090.0]
    k_vals = [r["K_static"] for r in rows if r["case"] == "same" and r["f_Hz"] == 50.0]
    val_path = AC / "slot" / "val.json"
    if val_path.exists():
        k_vals += [r["K_static"] for r in json.loads(val_path.read_text(encoding="utf-8")) if r["f_Hz"] == 5000.0]
    k_same = float(np.mean(k_vals))
    k_diff = float(np.mean([r["K_static"] for r in low if r["case"] == "diff"]))
    fs = sorted({r["f_Hz"] for r in rows})

    def c_of(winding, case):
        out = []
        for f in fs:
            sel = [r for r in rows if r["winding"] == winding and r["case"] == case and r["f_Hz"] == f]
            vals = []
            for r in sel:
                d_m = sweep["windings"][winding]["d"] * 1e-3
                skin_part = r["p_dc"] * (float(ac_theory.skin_ratio(d_m, f, rho)) - 1.0)
                prox_part = r["formula_exact"] - skin_part
                vals.append((r["excess"] - skin_part) / prox_part)
            out.append(float(np.mean(vals)))
        return out
    c_w1 = c_of("W1", "same")
    c_w2 = c_of("W2", "same")
    return {"K_same_m2": k_same, "K_diff_m2": k_diff, "w_same": 6.0 / 9.0, "K_spread_same": float(np.std(k_vals) / k_same), "n_layouts": len(k_vals),
            "c_f": [1.0] + fs, "c": [1.0] + c_w1, "c_w2": [1.0] + c_w2,
            "note": "K from the static field of every layout (sweep: W1, W2 x 3 layouts, val: 6 more windings); c(f) = FEMM's proximity loss over the static-field law with the lone wire's exact factors, "
                    "W1 (0.315 mm strands), layers in phase; c_w2 is the same for W2 (0.63 mm) as a check"}


def circ_block() -> dict | None:
    path = AC / "slot" / "circ_model.json"
    if not path.exists():
        return None
    c = json.loads(path.read_text(encoding="utf-8"))
    return {"v_A": c["v_A"], "ns": c["ns"], "kappa": c["kappa"], "rolloff_fit": c["rolloff_fit"]}


def build(out_path: Path | None = None, log=print) -> Path:
    out_path = out_path or copper_loss.TABLES_PATH
    t0 = time.perf_counter()
    tables = {"provenance": {"built": time.strftime("%Y-%m-%d %H:%M:%S"), "by": "optimization/copper_loss_build.py", "sources": str(AC)},
              "geometry": geometry(), "machine": machine_block(), "pwm": pwm_block(), "slot": slot_block()}
    circ = circ_block()
    if circ is not None:
        tables["circ"] = circ
    out_path.write_text(json.dumps(tables, indent=1), encoding="utf-8")
    log(f"written {out_path} ({out_path.stat().st_size / 1e3:.1f} kB, {time.perf_counter() - t0:.1f} s)")
    return out_path
