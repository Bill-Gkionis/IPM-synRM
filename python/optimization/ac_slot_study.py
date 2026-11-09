"""optimization.ac_slot_study: the slot-model studies of the copper-loss work: R_ac / R_dc of two real windings from the fundamental to the PWM frequency, the strand-by-strand
check of the winding study's formula, the current-density pictures, the circulating currents of strands in parallel, and how much each modelling choice matters.

What it is, and why
    optimization/ac_slot.py draws one slot of the optimized stator with every strand in it and solves it in FEMM's time-harmonic solver. This module decides which solves to make, runs
    them (all the stages asked for share ONE pool session, resumable: a stage's results are kept in out/optimization/copper_ac/slot/raw/<stage>.pkl and only the missing tasks are run) and turns the
    per-strand numbers into the answers the brief asks for. The two windings are the winding study's example (section 9.5 of the Optimization note) and its thick-strand twin, both at the
    rated 118.6 ampere-turns a coil side (the same N / a = 6.5, so the same 18.25 A peak phase current and the same field):
        W1   26 turns a coil side x 8 strands of 0.315 mm (0.367 mm overall), 4 parallel paths: 208 strands a layer, fill 0.448, 0.570 A peak a strand
        W2   13 turns a coil side x 4 strands of 0.630 mm (0.704 mm overall, interpolated), 2 paths: 52 strands a layer, fill 0.448, 2.281 A peak a strand
    The layers of a slot carry the same current in 6 of the 9 slots of the sector ('same': both layers of one phase) and currents 60 degrees apart in the other 3 ('diff'); both are made.

The stages (python -m optimization.ac_slot_study --stage <name>; 'all' runs them in turn)
    sweep    R_ac / R_dc against frequency, 50 Hz to 200 kHz, series circuits (every strand carries its share: equal sharing), three random layouts of each winding for 'same' and
             one for 'diff', with the magnetostatic solves of the same layouts that give the field the formula needs (the strands' own, and the homogenised slot's). It also gives the
             formula-against-FEMM comparison, strand by strand and by depth. Writes sweep.json, sweep.csv and the figures fig_slot_rac, fig_slot_collapse and fig_slot_formula.
    maps     the current density inside the strands (|J| over the DC density of the same current) and the loss of each strand over its DC loss: fig_slot_maps.png.
    circ     the strands of one turn in PARALLEL circuits, placed as 'line' (the worst: every turn spans the slot's depth), 'random', 'cluster' (compact blobs) and 'row' (side by side
             at one depth), against the same layouts with equal sharing: the extra loss at the fundamental, at the slot harmonic and at the switching frequency. fig_slot_circ.png.
    sens     what each modelling choice does to R_ac / R_dc at 5 and 20 kHz: the cell's width and boundaries, the steel's lamination and permeability, the rotor's surface, the mesh.
    rotor    a uniform field applied from the rotor side (A = -B x on the cell's outer boundary), strands with no current: the strands near the opening in the rotor's field.
    val      other strand sizes and fills (the examples of the winding study's table): the scaling law R_ac / R_dc - 1 = K (d k / delta^2)^2 of the next section, tested.
    Each stage's FEMM time is under 3 minutes (the pool session ends by itself); another run holding the solver slot may make it wait.

What the studies find, in the order they are made (the numbers are printed by the run and quoted in the report; the reasoning is here)
    1. The formula as the winding study used it: P_k = pi sigma w^2 d^4 L |B_k|^2 / 128 for each strand, B_k the local peak field. With B_k from the magnetostatic solve of the same strands
       it reproduces FEMM's time-harmonic excess loss (loss above DC) to 0.1 % wherever the steel is lossless (sensitivity stage, sigma_steel = 0: FEMM / formula 1.0006 at 5 kHz and
       0.9992 at 20 kHz, d / delta 0.59 for W1; the law of learn/proximity/04_rac_vs_frequency.py holds to 0.2 % up to 40 kHz and 1.3 % at 100 kHz), strand by strand (at 50 Hz the ratio is
       1.0009 for the strands next to the opening and 0.9993 in layer 1): the strands do NOT need to be drawn for the formula to hold, only the field has to be right.
    2. What the formula misses is the steel: the laminations' eddy currents (FEMM's complex effective permeability, 0.356 mm sheets) reduce the field the strands see, and the loss with it, by 4.1 %
       at 5 kHz and 11.4 % at 20 kHz; 0.635 mm sheets by 9.8 % and 21.2 %, 0.2 mm sheets by 0.8 % and 5.3 %. The static field is the field of infinitely thin sheets. The correction is a
       function of frequency alone (c_iron), made in the sweep.
    3. A slot of ANY winding then follows one law: with the copper fill k = N n_s (pi d^2 / 4) / A_layer and the skin depth delta, the extra loss over DC, divided by the DC loss, is
           (R_ac - R_dc) / R_dc = K (d k / delta^2)^2 c_iron(f) g(d / delta)
       where K (a slot-geometry constant [m^2]) comes from ONE magnetostatic solve (its mean-square field per unit strand current), g is the exact lone-wire proximity factor of ac_theory, and the
       corrections c_iron and a high-frequency correction are read from the sweep: the fast model uses this (optimization/copper_loss.py); the val stage tests it on other windings.
    4. In parallel circuits the strands of a turn share a voltage, and their currents divide as their positions in the field make them: the loss against equal sharing is the circulating-current
       loss. It goes as the square of the spread of the vector potential A over the strands of a turn: large when the strands of a turn are spread in depth (line, random), small when they lie
       side by side at one depth (row: the field of a slot points across it, so A depends on depth only).

Run (from the Python to FEMM folder)
    uv run python -m optimization.ac_slot_study --stage sweep       (about 1 minute of FEMM)
    uv run python -m optimization.ac_slot_study                     (every stage)
Everything goes to out/optimization/copper_ac/slot/ . The results of a stage are tables (json, csv) and figures.

Numbers to expect: see the printed summaries; the headline (FEMM 5.0.20261101, mesh 60 thousand nodes for 416 strands in a three-slot cell, same-phase layers, mean of three random layouts,
the range in brackets): W1's R_ac / R_dc is 1.00037 at 200 Hz, 1.0109 at 1,090 Hz (the fundamental at 16,330 rpm), 1.855 at 10 kHz and 4.27 (4.17 ... 4.40) at 20 kHz; W2's is 1.0016 at 200 Hz,
1.0472 at 1,090 Hz, 4.67 at 10 kHz and 14.8 (14.5 ... 15.0) at 20 kHz. Layers 60 degrees apart: 1.0003 and 3.80 for W1, 1.0013 and 12.5 for W2 (the sweep table). The whole study (245 tasks,
8 workers) takes 106 s of FEMM. The slot constant K is 15.28 mm^2 for W1 (W2 16.48) with the layers in phase, 12.72 (W2 13.31) with the layers 60 degrees apart.

Self-check (the answers follow).
  1. W1 and W2 have the same copper fill, the same ampere-turns and the same slot. Why is W2's R_ac / R_dc larger at the same frequency?
  2. FEMM's loss falls 4 % below the formula's at 5 kHz and 11 % below it at 20 kHz, while the formula is exact to 0.1 % in the lossless-steel run. Which part of the model is responsible?
  3. Why does the circulating loss of a 'row' layout stay small when that of a 'line' layout is large, at the same frequency and the same strands?
Answers.
  1. The proximity loss of a winding goes as d^2 at a given copper area and field, and W2's strands have twice the diameter: four times the extra loss, until d / delta is large enough for the
     strands' self-shielding (the exact factor of ac_theory.proximity_ratio) to take over.
  2. The steel: the laminations' eddy currents (their effective complex permeability, whose imaginary part grows with frequency) reduce the field the strands see, more the higher the
     frequency and the thicker the sheet (0.2 mm: 0.8 % and 5.3 % at 5 and 20 kHz; 0.356 mm: 4.1 % and 11.4 %; 0.635 mm: 9.8 % and 21.2 %). The formula's inputs are right (copper),
     the field's source is not (iron), and it is the one number of the slot model that rests on an assumed sheet thickness.
  3. In a row the strands of a turn are at the same depth, where A is nearly the same (a slot's field points across it, A varies with the depth); in a line they sample every depth, so the EMF between
     two of them is the whole potential difference across the layer, and the current that EMF drives is the loss.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import time
from pathlib import Path

import numpy as np

import femm_queue
import femm_session
import paths
from optimization import ac_slot, ac_theory, design, winding

OUT = paths.path("out") / "optimization" / "copper_ac" / "slot"
RAW = OUT / "raw"
FEM = OUT / "fem"
AT_RATED = 118.6                                      # ampere-turns a coil side, the brief's rated point
T_COPPER = 100.0
RHO = float(ac_theory.resistivity(T_COPPER))
SIGMA = 1.0 / RHO
PHASE_DIFF = complex(math.cos(math.pi / 3.0), -math.sin(math.pi / 3.0))        # the second layer's current lags the first by 60 deg in the slots of two phases
_CAT = winding.catalogue()
D_OVER = {round(float(d), 3): float(o) for d, o in zip(_CAT["d_mm"], _CAT["d_over_mm"])}
WINDINGS = {"W1": {"N": 26, "a": 4, "n_s": 8, "d": 0.315}, "W2": {"N": 13, "a": 2, "n_s": 4, "d": 0.630}}
for _w in WINDINGS.values():
    _w["d_over"] = D_OVER[round(_w["d"], 3)]
    _w["n_layer"] = _w["N"] * _w["n_s"]
    _w["i_strand"] = AT_RATED / _w["n_layer"]                 # A peak in each strand: the layer's total is the 118.6 ampere-turns
FREQS = (50.0, 100.0, 200.0, 500.0, 1090.0, 2000.0, 5000.0, 10000.0, 20000.0, 40000.0, 60000.0, 100000.0, 200000.0)
SEEDS = (1, 2, 3)
KEY_FREQS = {"fundamental, 3,000 rpm": 200.0, "fundamental, 16,330 rpm": 1090.0, "PWM 10 kHz": 10000.0, "PWM 20 kHz": 20000.0}
SENS_FREQS = (5000.0, 20000.0)
SENS_VARIANTS = {
    "baseline: 3 pitches, A=0 cuts, 0.356 mm laminations, mu_r 4416": {},
    "1 pitch wide": {"n_pitch": 1},
    "5 pitches wide": {"n_pitch": 5},
    "rotor surface A = 0 (shielded rotor)": {"rotor": "dirichlet"},
    "steel without eddy currents (sigma = 0)": {"steel": {"mu_r": 4416.0, "sigma_MS": 0.0, "d_lam_mm": 0.0, "lam_fill": 1.0}},
    "steel mu_r = 1000, no eddy currents": {"steel": {"mu_r": 1000.0, "sigma_MS": 0.0, "d_lam_mm": 0.0, "lam_fill": 1.0}},
    "steel mu_r = 100000, no eddy currents": {"steel": {"mu_r": 1e5, "sigma_MS": 0.0, "d_lam_mm": 0.0, "lam_fill": 1.0}},
    "laminations 0.200 mm": {"steel": {"mu_r": 4416.0, "sigma_MS": 1.9, "d_lam_mm": 0.20, "lam_fill": 0.98}},
    "laminations 0.635 mm (the project's FEMM steel)": {"steel": {"mu_r": 4416.0, "sigma_MS": 1.9, "d_lam_mm": 0.635, "lam_fill": 0.98}},
    "steel mu_r 10000, 0.356 mm laminations": {"steel": {"mu_r": 1e4, "sigma_MS": 1.9, "d_lam_mm": 0.356, "lam_fill": 0.98}},
    "finer mesh (arcs 5 deg, strands d/12, air 0.08)": {"arc_deg": 5.0, "mesh_cu_over_d": 1.0 / 12.0, "mesh_air_mm": 0.08},
    "coarser mesh (arcs 15 deg, strands d/4, air 0.25)": {"arc_deg": 15.0, "mesh_cu_over_d": 0.25, "mesh_air_mm": 0.25},
}
CIRC_FREQS = (200.0, 1090.0, 5000.0, 20000.0, 60000.0)
CIRC_ARRANGEMENTS = ("line", "random", "cluster", "row")
VAL_WINDINGS = {                                       # other strand sizes and fills, from the winding study's table of section 9.5 and two of its own
    "V1: 20 x 13 x 0.280": {"N": 20, "n_s": 13, "d": 0.280},
    "V2: 27 x 19 x 0.200": {"N": 27, "n_s": 19, "d": 0.200},
    "V3: 14 x 29 x 0.224": {"N": 14, "n_s": 29, "d": 0.224},
    "V4: 20 x 7 x 0.315 (fill 0.30)": {"N": 20, "n_s": 7, "d": 0.315},
    "V5: 9 x 9 x 0.500": {"N": 9, "n_s": 9, "d": 0.500},
    "V6: 5 x 4 x 1.000": {"N": 5, "n_s": 4, "d": 1.000},
}
VAL_FREQS = (5000.0, 20000.0, 100000.0)


def overrides() -> dict:
    return json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]


# ---------------------------------------------------------------------------------------------- running stages (one pool session each, resumable)
def run_stages(names: tuple, workers: int = 8, timeout_s: int = 600, log=print) -> dict:
    """Run the tasks of the stages in `names` that are not in raw/<stage>.pkl yet, ALL IN ONE pool session (one wait for the solver slot, one pool start-up), keep each stage's results
    in its own file, and return {stage: {tag: result}}. A task that fails stops the run after the others' results are saved (so a re-run only repeats what is missing)."""
    RAW.mkdir(parents=True, exist_ok=True)
    stage_tasks = {name: STAGE_TASKS[name]() for name in names}
    stores: dict = {}
    todo: list[tuple[str, dict]] = []
    for name in names:
        path = RAW / f"{name}.pkl"
        stores[name] = pickle.load(open(path, "rb")) if path.exists() else {}
        todo += [(name, t) for t in stage_tasks[name] if t["tag"] not in stores[name]]
    if todo:
        t0 = time.perf_counter()
        with femm_queue.femm_pool("copper ac: slot " + "+".join(names), workers=workers, task_timeout_s=timeout_s) as pool:
            waited = time.perf_counter() - t0
            t1 = time.perf_counter()
            res = pool.map(ac_slot.run_task, [t for _, t in todo])
        bad = []
        for (name, t), r in zip(todo, res):
            if isinstance(r, dict) and r.get("status") == "ok":
                stores[name][t["tag"]] = r
            else:
                bad.append((t["tag"], r))
        for name in names:
            with open(RAW / f"{name}.pkl", "wb") as handle:
                pickle.dump(stores[name], handle)
        log(f"[{'+'.join(names)}] {len(todo) - len(bad)} of {len(todo)} tasks done in {time.perf_counter() - t1:.1f} s of FEMM ({waited:.1f} s waiting for the slot)")
        if bad:
            raise RuntimeError("; ".join(f"{tag}: {r.get('error', r) if isinstance(r, dict) else r}" for tag, r in bad[:3]))
    else:
        log(f"[{'+'.join(names)}] all {sum(len(v) for v in stage_tasks.values())} tasks are saved already")
    return {name: {t["tag"]: stores[name][t["tag"]] for t in stage_tasks[name]} for name in names}


def layers_for(w: dict, case: str) -> tuple:
    """The two layers' strand currents [A peak] of the winding for the 'same' or 'diff' phase relation."""
    i = w["i_strand"]
    return (i, i) if case == "same" else (i, i * PHASE_DIFF)


# ---------------------------------------------------------------------------------------------- per-strand numbers from the three solves of a layout
def strand_table(ac: dict, ms: dict, hom: dict | None, i_layers: tuple, f: float, d_mm: float) -> dict:
    """Every strand's numbers at frequency f: the FEMM loss, the DC loss, and the formula's predictions from the magnetostatic solves (module docstring). All in W for the stack of the
    model; arrays over strands. `ms` and `hom` are the layout's own magnetostatic results (unit current a strand in layer 0 alone, and in layer 1 alone): the complex layer currents combine them.
    P_formula = pi sigma w^2 d_eff^4 L |B|^2 / 128 with B the mean field over the strand (and with the homogenised slot's field at the strand's centre for 'homog'); formula_exact adds the lone
    wire's proximity factor and its skin loss (it predicts the EXCESS over DC: p_excess_exact); low_f is the exact low-frequency eddy loss from the potential A (sigma w^2 L / 2 times the integral of |A - mean|^2 over the strand)."""
    w = 2.0 * math.pi * f
    L = ac["stack_mm"] * 1e-3
    area = ac["area_m2"]
    d_eff = np.sqrt(4.0 * area / math.pi)
    i0, i1 = complex(i_layers[0]), complex(i_layers[1])
    bx = i0 * ms["Bx_mean"][0] + i1 * ms["Bx_mean"][1]
    by = i0 * ms["By_mean"][0] + i1 * ms["By_mean"][1]
    b2 = np.abs(bx) ** 2 + np.abs(by) ** 2
    p_dc = np.array([abs(i_layers[li]) for li in ac["layer"]]) ** 2 * RHO * L / (2.0 * area)
    i_sq = abs(i0) ** 2 * ms["A_dev_integral"][0] + abs(i1) ** 2 * ms["A_dev_integral"][1] + 2.0 * np.real(i0 * np.conj(i1)) * ms["A_dev_cross"]
    out = {"p_femm": ac["loss_W"], "p_dc": p_dc, "d_eff": d_eff, "xy": ac["xy"], "layer": ac["layer"], "b2": b2,
           "p_formula": math.pi * SIGMA * w ** 2 * d_eff ** 4 * L * b2 / 128.0, "p_low_f": SIGMA * w ** 2 * L / 2.0 * i_sq}
    if hom is not None:
        bxh = i0 * hom["Bx"][0] + i1 * hom["Bx"][1]
        byh = i0 * hom["By"][0] + i1 * hom["By"][1]
        out["p_homog"] = math.pi * SIGMA * w ** 2 * d_eff ** 4 * L * (np.abs(bxh) ** 2 + np.abs(byh) ** 2) / 128.0
    ratio = ac_theory.proximity_ratio(d_eff, f, RHO)
    skin = ac_theory.skin_ratio(d_eff, f, RHO)
    out["p_excess_exact"] = out["p_formula"] * ratio + p_dc * (skin - 1.0)
    out["prox_ratio"], out["skin_ratio"] = ratio, skin
    return out


DEPTH_GROUPS = (("0 - 1 mm below the coil's top edge", 0.0, 1.0), ("1 - 2 mm", 1.0, 2.0), ("2 mm to the layers' interface", 2.0, 6.55), ("layer 1 (the bottom half)", 6.55, 99.0))


def depth_below_top(p, xy: np.ndarray) -> np.ndarray:
    top = ac_slot.layer_polygons(p)[0].bounds[1]
    return xy[:, 1] - top


def summarise(tab: dict, p) -> dict:
    """Totals and ratios of one layout at one frequency: R_ac / R_dc, the excess loss (above DC) against the formulas, and the same by depth group. 'exact' is the whole lone-wire theory:
    the formula times the exact proximity factor, plus the skin loss of each strand."""
    exc = tab["p_femm"] - tab["p_dc"]
    out = {"R": float(tab["p_femm"].sum() / tab["p_dc"].sum()), "excess": float(exc.sum()), "p_dc": float(tab["p_dc"].sum()), "formula": float(tab["p_formula"].sum()),
           "formula_exact": float(tab["p_excess_exact"].sum()), "low_f": float(tab["p_low_f"].sum()), "ratio_formula": float(exc.sum() / tab["p_formula"].sum()),
           "ratio_exact": float(exc.sum() / tab["p_excess_exact"].sum())}
    if "p_homog" in tab:
        out["homog"] = float(tab["p_homog"].sum())
        out["ratio_homog"] = float(exc.sum() / tab["p_homog"].sum())
    depth = depth_below_top(p, tab["xy"])
    groups = []
    for name, lo, hi in DEPTH_GROUPS:
        m = (depth >= lo) & (depth < hi)
        groups.append({"group": name, "n": int(m.sum()), "excess": float(exc[m].sum()), "formula": float(tab["p_formula"][m].sum()),
                       "ratio": float(exc[m].sum() / tab["p_formula"][m].sum()) if m.any() else float("nan"), "share_of_excess": float(exc[m].sum() / exc.sum())})
    out["groups"] = groups
    return out


# ---------------------------------------------------------------------------------------------- stage: sweep
def sweep_tasks() -> list[dict]:
    ov = overrides()
    tasks = []
    for wn, w in WINDINGS.items():
        for seed in SEEDS:
            tasks.append(ac_slot.task_ms(f"ms_{wn}_s{seed}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], seed=seed, t_c=T_COPPER))
            tasks.append(ac_slot.task_hom(f"hom_{wn}_s{seed}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], seed=seed))
            for f in FREQS:
                tasks.append(ac_slot.task_ac(f"ac_{wn}_same_s{seed}_{f:g}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], f=f, i_layers=layers_for(w, "same"),
                                             seed=seed, t_c=T_COPPER))
        for f in FREQS:
            tasks.append(ac_slot.task_ac(f"ac_{wn}_diff_s1_{f:g}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], f=f, i_layers=layers_for(w, "diff"),
                                         seed=1, t_c=T_COPPER))
    return tasks


def analyse_sweep(res: dict, log=print) -> dict:
    p = design.params(overrides())
    rows, tables = [], {}
    for wn, w in WINDINGS.items():
        for case, seeds in (("same", SEEDS), ("diff", (1,))):
            for seed in seeds:
                ms, hom = res[f"ms_{wn}_s{seed}"], res[f"hom_{wn}_s{seed}"]
                for f in FREQS:
                    ac = res[f"ac_{wn}_{case}_s{seed}_{f:g}"]
                    tab = strand_table(ac, ms, hom, layers_for(w, case), f, w["d"])
                    s = summarise(tab, p)
                    delta = float(ac_theory.skin_depth(f, RHO))
                    k_fill = w["n_layer"] * math.pi * w["d"] ** 2 / 4.0 / 36.218
                    law = (w["d"] * 1e-3 * k_fill / delta ** 2) ** 2                          # (d k / delta^2)^2 [1 / m^2]
                    s.update({"winding": wn, "case": case, "seed": seed, "f_Hz": f, "d_over_delta": w["d"] * 1e-3 / delta, "fill": k_fill,
                              "K": (s["R"] - 1.0) / law, "K_static": s["formula"] / s["p_dc"] / law})
                    rows.append(s)
                    tables[(wn, case, seed, f)] = tab
    return {"rows": rows, "tables": tables}


def print_sweep(an: dict, log=print) -> None:
    rows = an["rows"]
    log("\nR_ac / R_dc of the slot (series circuits: equal sharing), mean over the 3 layouts for 'same', one layout for 'diff'")
    log("   f [Hz]   d/delta W1 | W1 same    (3 layouts min ... max)   W1 diff | d/delta W2 | W2 same    (min ... max)   W2 diff")
    for f in FREQS:
        cells = []
        for wn in WINDINGS:
            same = [r for r in rows if r["winding"] == wn and r["case"] == "same" and r["f_Hz"] == f]
            diff = [r for r in rows if r["winding"] == wn and r["case"] == "diff" and r["f_Hz"] == f][0]
            vals = [r["R"] for r in same]
            cells.append(f"{same[0]['d_over_delta']:6.3f} | {np.mean(vals):9.5f} ({min(vals):8.5f} ... {max(vals):8.5f})  {diff['R']:9.5f}")
        log(f"  {f:8g}   {cells[0]} | {cells[1]}")
    log("\nthe slot constant K = (R_ac / R_dc - 1) / (d k / delta^2)^2 [mm^2] from the STATIC field (the formula) and from FEMM; same-phase layers (mean of 3 layouts) and layers 60 deg apart")
    log("   f [Hz]   W1: K static   K FEMM   | 60 deg: K static   K FEMM   ||  W2: K static   K FEMM   | 60 deg: K static   K FEMM")
    for f in FREQS:
        cells = []
        for wn in WINDINGS:
            same = [r for r in rows if r["winding"] == wn and r["case"] == "same" and r["f_Hz"] == f]
            diff = [r for r in rows if r["winding"] == wn and r["case"] == "diff" and r["f_Hz"] == f][0]
            cells.append(f"{1e6 * np.mean([r['K_static'] for r in same]):10.3f} {1e6 * np.mean([r['K'] for r in same]):9.3f}  | {1e6 * diff['K_static']:12.3f} {1e6 * diff['K']:9.3f}")
        log(f"  {f:8g}   {cells[0]}  ||  {cells[1]}")
    log("\nthe winding study's formula against FEMM (excess loss over DC, FEMM / formula), same-phase layers, mean of 3 layouts")
    log("   f [Hz]   W1: FEMM/formula  FEMM/homog  FEMM/low-f   near opening (0-1 mm)  1-2 mm   rest of layer 0   layer 1 |  W2: FEMM/formula  FEMM/homog  FEMM/low-f   near opening")
    for f in FREQS:
        cells = []
        for wn in WINDINGS:
            same = [r for r in rows if r["winding"] == wn and r["case"] == "same" and r["f_Hz"] == f]
            fm = lambda key: float(np.mean([r[key] for r in same]))
            gr = lambda j: float(np.mean([r["groups"][j]["ratio"] for r in same]))
            if wn == "W1":
                cells.append(f"{fm('ratio_formula'):8.4f}  {fm('ratio_homog'):10.4f}  {fm('excess') / fm('low_f'):10.4f}   {gr(0):14.4f}  {gr(1):9.4f}  {gr(2):9.4f}  {gr(3):9.4f}")
            else:
                cells.append(f"{fm('ratio_formula'):8.4f}  {fm('ratio_homog'):10.4f}  {fm('excess') / fm('low_f'):10.4f}   {gr(0):14.4f}")
        log(f"  {f:8g}   {cells[0]} |  {cells[1]}")


def write_sweep(an: dict) -> None:
    rows = an["rows"]
    keep = ["winding", "case", "seed", "f_Hz", "d_over_delta", "fill", "R", "excess", "p_dc", "formula", "formula_exact", "low_f", "ratio_formula", "ratio_homog", "K"]
    with open(OUT / "sweep.csv", "w", newline="", encoding="utf-8") as handle:
        w = csv.DictWriter(handle, fieldnames=keep, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    (OUT / "sweep.json").write_text(json.dumps({"windings": WINDINGS, "freqs": FREQS, "t_copper_C": T_COPPER, "rows": rows, "key_freqs": KEY_FREQS}, indent=1, default=float), encoding="utf-8")


# ---------------------------------------------------------------------------------------------- stage: maps
MAP_CASES = (("W1", 20000.0), ("W2", 20000.0), ("W2", 100000.0))


def maps_tasks() -> list[dict]:
    ov = overrides()
    return [ac_slot.task_ac(f"map_{wn}_{f:g}", str(FEM), ov, n_turns=WINDINGS[wn]["N"], n_s=WINDINGS[wn]["n_s"], d_mm=WINDINGS[wn]["d"], d_over_mm=WINDINGS[wn]["d_over"], f=f,
                            i_layers=layers_for(WINDINGS[wn], "same"), seed=1, t_c=T_COPPER, want_field=True) for wn, f in MAP_CASES]


# ---------------------------------------------------------------------------------------------- stage: circ
def circ_tasks() -> list[dict]:
    ov = overrides()
    tasks = []
    for wn, w in WINDINGS.items():
        tasks.append(ac_slot.task_ms(f"circ_ms_{wn}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], seed=1, t_c=T_COPPER))
        for f in CIRC_FREQS:
            tasks.append(ac_slot.task_ac(f"circ_{wn}_series_{f:g}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], f=f, i_layers=layers_for(w, "same"),
                                         seed=1, t_c=T_COPPER, circuit="series"))
            for arr in CIRC_ARRANGEMENTS:
                tasks.append(ac_slot.task_ac(f"circ_{wn}_{arr}_{f:g}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], f=f, i_layers=layers_for(w, "same"),
                                             seed=1, t_c=T_COPPER, circuit="parallel", arrangement=arr))
    return tasks


def circulating_lowf(par: dict, ms: dict, i_layers: tuple, f: float) -> float:
    """The low-frequency circulating-current loss [W] predicted from the static field: strand i of a turn carries, besides its share, the current j w sigma S (mean A of the turn - A_i), so the loss is
    sum over strands of w^2 sigma S |A_i - A_turn|^2 L / 2 (peak phasors; the strand length L is the same for every strand in the 2-D model). A_i is the mean over the strand of the vector
    potential of the magnetostatic solve (ms), combined for the complex layer currents. Valid while the circulating currents' own fields are small (low frequency)."""
    w = 2.0 * math.pi * f
    L = par["stack_mm"] * 1e-3
    i0, i1 = complex(i_layers[0]), complex(i_layers[1])
    a = i0 * ms["A_mean"][0] + i1 * ms["A_mean"][1]
    area = par["area_m2"]
    total = 0.0
    for li in range(2):
        for t in np.unique(par["turn"][par["layer"] == li]):
            sel = (par["layer"] == li) & (par["turn"] == t)
            total += float(np.sum(w ** 2 * SIGMA * area[sel] * np.abs(a[sel] - a[sel].mean()) ** 2 * L / 2.0))
    return total


def analyse_circ(res: dict, log=print) -> list[dict]:
    rows = []
    for wn, w in WINDINGS.items():
        ms = res[f"circ_ms_{wn}"]
        for f in CIRC_FREQS:
            base = res[f"circ_{wn}_series_{f:g}"]
            p_base = float(base["loss_W"].sum())
            p_dc = float(ac_slot.dc_loss_W(base, layers_for(w, "same")).sum())
            for arr in CIRC_ARRANGEMENTS:
                par = res[f"circ_{wn}_{arr}_{f:g}"]
                p_par = float(par["loss_W"].sum())
                cur = par["current_A"]
                i_share = w["i_strand"]
                dev = np.abs(cur - np.array([complex(layers_for(w, "same")[li]) for li in par["layer"]]))
                lowf = circulating_lowf(par, ms, layers_for(w, "same"), f)
                rows.append({"winding": wn, "arrangement": arr, "f_Hz": f, "d_over_delta": w["d"] * 1e-3 / float(ac_theory.skin_depth(f, RHO)), "p_series_W": p_base, "p_parallel_W": p_par,
                             "extra_W": p_par - p_base, "extra_over_dc": (p_par - p_base) / p_dc, "R_series": p_base / p_dc, "R_parallel": p_par / p_dc,
                             "lowf_prediction_W": lowf, "extra_over_lowf": (p_par - p_base) / lowf if lowf > 0 else float("nan"),
                             "i_dev_rms_over_share": float(np.sqrt(np.mean(dev ** 2)) / i_share)})
    log("\ncirculating currents: the loss of strands in PARALLEL circuits (one circuit a turn) minus the same layouts with equal sharing, over the DC loss; same-phase layers, seed 1")
    for wn in WINDINGS:
        log(f"  {wn}:  f [Hz]    " + "   ".join(f"{arr:>22s}" for arr in CIRC_ARRANGEMENTS) + "      (R_ac/R_dc with equal sharing)")
        for f in CIRC_FREQS:
            sel = {r["arrangement"]: r for r in rows if r["winding"] == wn and r["f_Hz"] == f}
            log(f"      {f:8g}   " + "   ".join(f"{sel[a]['extra_over_dc']:12.4g} (i dev {sel[a]['i_dev_rms_over_share']:5.2f})" for a in CIRC_ARRANGEMENTS) + f"      {sel['line']['R_series']:.5f}")
        log(f"  {wn}: FEMM's extra over the low-frequency prediction from the static potential (1 = the prediction holds):")
        for f in CIRC_FREQS:
            sel = {r["arrangement"]: r for r in rows if r["winding"] == wn and r["f_Hz"] == f}
            log(f"      {f:8g}   " + "   ".join(f"{sel[a]['extra_over_lowf']:22.4f}" for a in CIRC_ARRANGEMENTS))
    return rows


# ---------------------------------------------------------------------------------------------- the circulating-current model, from the circ stage and the sweep's static solves
CIRC_NS = (2, 4, 8, 13, 16, 26, 52)                      # strands in hand for which the within-turn spread of A is measured (divisors of the 208 strands of W1's layer)


def build_circ_model(sweep_res: dict, circ_rows: list[dict], log=print) -> dict:
    """The constants of the fast model of the circulating currents (optimization/copper_loss.py), from what the studies measured:
      v_A      the spread of the vector potential over the strands of one layer, per (the layer's total current)^2 [Wb^2 / m^2 / A^2]: with all strands in phase,
               V_A = mean over a layer's strands of (A_i - mean A)^2 from the magnetostatic solves of the sweep (3 layouts of W1 and 3 of W2: they must agree, the slot is the same);
      kappa    the share of that spread the strands of one TURN see, for each placement: (the mean over turns of the sum over the turn's strands of (A_i - A_turn)^2) over (the sum over the
               layer's strands of (A_i - A_layer)^2), for turns of n_s strands cut from W1's layouts (CIRC_NS); 'random' follows 1 - 1 / n_s, 'line' 1 - 1 / n_s^2, 'row' and 'cluster' are
               small: the field of a slot points across it, so A depends on the depth, and strands side by side at one depth see nearly the same A;
      rolloff  FEMM's circulating loss over the low-frequency prediction from the static potential (circulating_lowf), against x = (d / delta)^2 for each placement: 1 at low frequency; the currents
               the strands exchange make a field of their own that opposes the EMF between them, so the loss grows slower than f^2 and then saturates.
    The loss of the strands of a layer at frequency f, per DC loss of the same strands, is then  Q = (w sigma k A_layer)^2 kappa v_A rolloff(x)  (the low-frequency law times the roll-off)."""
    out: dict = {"ns": list(CIRC_NS)}
    # ---- v_A
    vals = {}
    for wn, w in WINDINGS.items():
        v = []
        for seed in SEEDS:
            ms = sweep_res[f"ms_{wn}_s{seed}"]
            a = ms["A_mean"][0] + ms["A_mean"][1]                       # 1 A in every strand of both layers: the layer total is n_layer amperes
            for li in range(2):
                sel = ms["layer"] == li
                v.append(float(np.mean((a[sel] - a[sel].mean()) ** 2)) / w["n_layer"] ** 2)
        vals[wn] = v
    v_a = float(np.mean(vals["W1"] + vals["W2"]))
    out["v_A"] = v_a
    out["v_A_by_winding"] = {wn: [float(np.mean(v)), float(np.std(v))] for wn, v in vals.items()}
    log(f"\nthe spread of A over the strands of a layer per (layer current)^2: W1 {np.mean(vals['W1']):.4e} (std {np.std(vals['W1']):.1e}), W2 {np.mean(vals['W2']):.4e} (std {np.std(vals['W2']):.1e}) Wb^2/m^2/A^2")
    # ---- kappa
    kappa = {}
    w1 = WINDINGS["W1"]
    for arr in ("random", "line", "cluster", "row"):
        ks = []
        for n_s in CIRC_NS:
            per_seed = []
            for seed in SEEDS:
                ms = sweep_res[f"ms_W1_s{seed}"]
                a = ms["A_mean"][0] + ms["A_mean"][1]
                n_turns = w1["n_layer"] // n_s
                turn = ac_slot.assign_turns(ms["xy"], ms["layer"], n_turns, n_s, arr, seed + 1000)
                within, whole = 0.0, 0.0
                for li in range(2):
                    sel = ms["layer"] == li
                    whole += float(np.sum((a[sel] - a[sel].mean()) ** 2))
                    for t in range(n_turns):
                        m = sel & (turn == t)
                        within += float(np.sum((a[m] - a[m].mean()) ** 2))
                per_seed.append(within / whole)
            ks.append(float(np.mean(per_seed)))
        kappa[arr] = ks
    out["kappa"] = kappa
    log("kappa (the within-turn share of A's spread) against the strands in hand, from W1's layouts (mean of 3):")
    log("   n_s   " + "  ".join(f"{n:6d}" for n in CIRC_NS))
    for arr, ks in kappa.items():
        log(f"   {arr:8s}" + "  ".join(f"{k:6.3f}" for k in ks))
    # ---- the roll-off: r(u) = 1 / (1 + (u / u_c)^p), u the frequency for line, random and row and the bundle size n_s (d / delta)^2 for cluster; fitted to FEMM's points where r > 0.02
    from scipy.optimize import least_squares
    n_s_of = {wn: w["n_s"] for wn, w in WINDINGS.items()}
    roll, fits = {}, {}
    for arr in ("line", "random", "cluster", "row"):
        pts = [((r["d_over_delta"] ** 2) * n_s_of[r["winding"]] if arr == "cluster" else r["f_Hz"], r["extra_over_lowf"], r["winding"], r["f_Hz"]) for r in circ_rows if r["arrangement"] == arr]
        pts.sort()
        roll[arr] = {"u": [p[0] for p in pts], "r": [p[1] for p in pts]}
        use = [(u, rr) for u, rr, _, _ in pts if rr > 0.02]
        uu, rr = np.array([p[0] for p in use]), np.array([p[1] for p in use])
        fit = least_squares(lambda q: np.log(1.0 / (1.0 + (uu / q[0]) ** q[1])) - np.log(rr), x0=[np.median(uu), 1.5], bounds=([1e-6, 0.3], [1e9, 6.0]))
        fits[arr] = {"u_c": float(fit.x[0]), "p": float(fit.x[1]), "rms_log": float(np.sqrt(np.mean(fit.fun ** 2))), "n": int(len(uu))}
    out["rolloff_points"] = roll
    out["rolloff_fit"] = fits
    log("roll-off r = FEMM's circulating loss over the low-frequency prediction, r(u) = 1 / (1 + (u / u_c)^p), u = f [Hz] (line, random, row) or n_s (d / delta)^2 (cluster):")
    for arr in ("line", "random", "cluster", "row"):
        f = fits[arr]
        log(f"   {arr:8s} u_c = {f['u_c']:10.4g}, p = {f['p']:.2f}, rms error of ln r {f['rms_log']:.2f} over {f['n']} points; points: " + "  ".join(f"{u:.4g}:{rr:.3f}" for u, rr in zip(roll[arr]["u"], roll[arr]["r"])))
    return out


# ---------------------------------------------------------------------------------------------- stage: sens
def sens_tasks() -> list[dict]:
    ov = overrides()
    w = WINDINGS["W1"]
    tasks = []
    for idx, (name, opt) in enumerate(SENS_VARIANTS.items()):
        key = f"v{idx:02d}"
        tasks.append(ac_slot.task_ms(f"sens_ms_{key}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], seed=1, options=opt))
        for f in SENS_FREQS:
            tasks.append(ac_slot.task_ac(f"sens_ac_{key}_{f:g}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], f=f, i_layers=layers_for(w, "same"), seed=1,
                                         options=opt))
    return tasks


def analyse_sens(res: dict, log=print) -> list[dict]:
    w = WINDINGS["W1"]
    rows = []
    for idx, name in enumerate(SENS_VARIANTS):
        key = f"v{idx:02d}"
        ms = res[f"sens_ms_{key}"]
        row = {"variant": name}
        for f in SENS_FREQS:
            ac = res[f"sens_ac_{key}_{f:g}"]
            tab = strand_table(ac, ms, None, layers_for(w, "same"), f, w["d"])
            exc = (tab["p_femm"] - tab["p_dc"]).sum()
            row[f"R_{f:g}"] = float(tab["p_femm"].sum() / tab["p_dc"].sum())
            row[f"formula_{f:g}"] = float(exc / tab["p_formula"].sum())
            row[f"nodes_{f:g}"] = int(ac["nodes"])
        rows.append(row)
    base = rows[0]
    log("\nsensitivity of W1's R_ac / R_dc (same-phase layers, seed 1) to the modelling choices, and the formula's ratio in each")
    log("  " + f"{'variant':58s}" + "".join(f" R at {f / 1e3:g} kHz   excess vs baseline   FEMM/formula  " for f in SENS_FREQS))
    for r in rows:
        cells = "".join(f" {r[f'R_{f:g}']:9.4f}   {100 * ((r[f'R_{f:g}'] - 1) / (base[f'R_{f:g}'] - 1) - 1):+9.1f} %        {r[f'formula_{f:g}']:8.3f}     " for f in SENS_FREQS)
        log(f"  {r['variant']:58s}{cells}  ({r['nodes_20000']} nodes)")
    return rows


# ---------------------------------------------------------------------------------------------- stage: rotor (a uniform field from the rotor side)
def rotor_tasks() -> list[dict]:
    ov = overrides()
    tasks = []
    for wn, w in WINDINGS.items():
        tasks.append(ac_slot.task_ms(f"rot_ms_{wn}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], seed=1, drive_field_T=1.0))
        for f in (2000.0, 10000.0, 20000.0, 50000.0, 100000.0, 200000.0):
            tasks.append(ac_slot.task_ac(f"rot_ac_{wn}_{f:g}", str(FEM), ov, n_turns=w["N"], n_s=w["n_s"], d_mm=w["d"], d_over_mm=w["d_over"], f=f, i_layers=(0.0, 0.0), seed=1,
                                         drive_field_T=1.0))
    return tasks


def analyse_rotor(res: dict, log=print) -> list[dict]:
    p = design.params(overrides())
    rows = []
    for wn, w in WINDINGS.items():
        ms = res[f"rot_ms_{wn}"]
        for f in (2000.0, 10000.0, 20000.0, 50000.0, 100000.0, 200000.0):
            ac = res[f"rot_ac_{wn}_{f:g}"]
            om = 2.0 * math.pi * f
            L = ac["stack_mm"] * 1e-3
            d_eff = np.sqrt(4.0 * ac["area_m2"] / math.pi)
            b2 = ms["Bx_mean"][0] ** 2 + ms["By_mean"][0] ** 2
            p_form = math.pi * SIGMA * om ** 2 * d_eff ** 4 * L * b2 / 128.0
            p_exact = p_form * ac_theory.proximity_ratio(d_eff, f, RHO)
            p_lf = SIGMA * om ** 2 * L / 2.0 * ms["A_dev_integral"][0]
            depth = depth_below_top(p, ac["xy"])
            groups = []
            for name, lo, hi in DEPTH_GROUPS:
                m = (depth >= lo) & (depth < hi)
                groups.append(float(ac["loss_W"][m].sum() / p_form[m].sum()))
            rows.append({"winding": wn, "f_Hz": f, "d_over_delta": w["d"] * 1e-3 / float(ac_theory.skin_depth(f, RHO)), "p_femm_W": float(ac["loss_W"].sum()), "femm_over_formula": float(ac["loss_W"].sum() / p_form.sum()),
                         "femm_over_formula_exact": float(ac["loss_W"].sum() / p_exact.sum()), "femm_over_low_f": float(ac["loss_W"].sum() / p_lf.sum()), "groups": groups,
                         "net_current_max_A": float(np.abs(ac["current_A"]).max())})
    log("\nthe rotor's field: a uniform field applied from the rotor side, strands with no current. FEMM's loss over the formula (static local field), by depth below the coil's top edge")
    for wn in WINDINGS:
        log(f"  {wn}:  f [Hz]  d/delta   FEMM/formula   FEMM/formula-exact   FEMM/low-f   by depth group: " + " | ".join(g[0] for g in DEPTH_GROUPS))
        for r in [r for r in rows if r["winding"] == wn]:
            log(f"      {r['f_Hz']:8g}  {r['d_over_delta']:6.3f}   {r['femm_over_formula']:9.4f}   {r['femm_over_formula_exact']:12.4f}   {r['femm_over_low_f']:9.4f}     " + "  ".join(f"{g:8.4f}" for g in r["groups"]))
    return rows


# ---------------------------------------------------------------------------------------------- stage: val (the scaling law on other strands)
def val_tasks() -> list[dict]:
    ov = overrides()
    tasks = []
    for name, v in VAL_WINDINGS.items():
        key = name.split(":")[0]
        d_over = D_OVER[round(v["d"], 3)]
        i_s = AT_RATED / (v["N"] * v["n_s"])
        tasks.append(ac_slot.task_ms(f"val_ms_{key}", str(FEM), ov, n_turns=v["N"], n_s=v["n_s"], d_mm=v["d"], d_over_mm=d_over, seed=1))
        for f in VAL_FREQS:
            tasks.append(ac_slot.task_ac(f"val_ac_{key}_{f:g}", str(FEM), ov, n_turns=v["N"], n_s=v["n_s"], d_mm=v["d"], d_over_mm=d_over, f=f, i_layers=(i_s, i_s), seed=1))
    return tasks


def analyse_val(res: dict, log=print) -> list[dict]:
    rows = []
    for name, v in VAL_WINDINGS.items():
        key = name.split(":")[0]
        d_over = D_OVER[round(v["d"], 3)]
        i_s = AT_RATED / (v["N"] * v["n_s"])
        ms = res[f"val_ms_{key}"]
        k_fill = v["N"] * v["n_s"] * math.pi * v["d"] ** 2 / 4.0 / 36.218
        for f in VAL_FREQS:
            ac = res[f"val_ac_{key}_{f:g}"]
            tab = strand_table(ac, ms, None, (i_s, i_s), f, v["d"])
            delta = float(ac_theory.skin_depth(f, RHO))
            R = float(tab["p_femm"].sum() / tab["p_dc"].sum())
            static = float(tab["p_formula"].sum() / tab["p_dc"].sum())
            rows.append({"winding": name, "fill": k_fill, "f_Hz": f, "d_over_delta": v["d"] * 1e-3 / delta, "R": R, "static_formula_over_dc": static,
                         "K_static": static / ((v["d"] * 1e-3 * k_fill / delta ** 2) ** 2), "K_femm": (R - 1.0) / ((v["d"] * 1e-3 * k_fill / delta ** 2) ** 2),
                         "ratio_formula": float((tab["p_femm"].sum() - tab["p_dc"].sum()) / tab["p_formula"].sum())})
    log("\nthe scaling law (R_ac - R_dc) / R_dc = K (d k / delta^2)^2: K from the STATIC field of each layout (one magnetostatic solve) and from FEMM's time-harmonic loss, same-phase layers")
    log("  winding                          fill    f [Hz]   d/delta   R_ac/R_dc   K static [mm^2]   K FEMM [mm^2]   FEMM/formula")
    for r in rows:
        log(f"  {r['winding']:30s} {r['fill']:6.3f}  {r['f_Hz']:8g}  {r['d_over_delta']:7.3f}   {r['R']:9.4f}   {1e6 * r['K_static']:14.3f}   {1e6 * r['K_femm']:13.3f}   {r['ratio_formula']:9.4f}")
    return rows


# ---------------------------------------------------------------------------------------------- the driver
STAGE_TASKS = {"sweep": sweep_tasks, "maps": maps_tasks, "circ": circ_tasks, "sens": sens_tasks, "rotor": rotor_tasks, "val": val_tasks}


def main() -> int:
    ap = argparse.ArgumentParser(description="the slot-model studies (module docstring)")
    ap.add_argument("--stage", choices=("all", *STAGE_TASKS), default="all")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    summary: dict = {"femm_version": femm_session.femm_version(), "t_copper_C": T_COPPER, "windings": WINDINGS, "options": ac_slot.DEFAULT_OPTIONS}
    stages = tuple(STAGE_TASKS) if args.stage == "all" else (args.stage,)
    results = run_stages(stages, args.workers)
    from optimization import ac_slot_figures
    for stage in stages:
        res = results[stage]
        if stage == "sweep":
            an = analyse_sweep(res)
            print_sweep(an)
            write_sweep(an)
            summary["sweep_rows"] = len(an["rows"])
            ac_slot_figures.fig_rac(an, OUT / "fig_slot_rac.png")
            ac_slot_figures.fig_collapse(an, OUT / "fig_slot_collapse.png")
            ac_slot_figures.fig_formula(an, res, OUT / "fig_slot_formula.png")
        elif stage == "maps":
            ac_slot_figures.fig_maps(res, OUT / "fig_slot_maps.png")
        elif stage == "circ":
            rows = analyse_circ(res)
            (OUT / "circ.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
            if "sweep" in results:
                (OUT / "circ_model.json").write_text(json.dumps(build_circ_model(results["sweep"], rows), indent=1), encoding="utf-8")
            ac_slot_figures.fig_circ(rows, OUT / "fig_slot_circ.png")
        elif stage == "sens":
            rows = analyse_sens(res)
            (OUT / "sens.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
        elif stage == "rotor":
            rows = analyse_rotor(res)
            (OUT / "rotor.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
        elif stage == "val":
            rows = analyse_val(res)
            (OUT / "val.json").write_text(json.dumps(rows, indent=1), encoding="utf-8")
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=str), encoding="utf-8")
    print(f"\ndone in {time.perf_counter() - started:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
