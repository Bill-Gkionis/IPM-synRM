"""optimization.thermal_study: what the thermal network says about the optimized motor: how much loss it can shed, the binding limit, the time constants, the overloads.

What it is, and why
    thermal.py is the network; this file asks it the design questions and writes the answers down (no solver is started: the network takes
    microseconds, and its finite-element check is thermal_fem.py). For each cooling (natural convection, forced air over a finned housing, a liquid
    jacket), each stack length (20, 40, 60, 80 mm: one to four rows of magnets) and each insulation class (F 155 C, H 180 C at the winding's hot spot),
    at the rated 3,000 rpm and a 40 C ambient (40 C coolant inlet):
      1. THE SUSTAINABLE LOSS. The continuous current density at which the first limit is reached: the winding's hot spot against the class, or the
         magnets against N45's 80 C (a hard limit). The copper loss then is the DC loss rho(T) J^2 V_cu at the temperatures the network
         finds (copper 0.45 filled, the winding study's mean turn 2 L + pi span + 4 e), with the PROVISIONAL iron-loss placeholder, the bearing and
         windage estimates, and no magnet eddy loss (the input defaults to 0). The ampere-turns per coil side are J times the copper area of a coil
         side (peak = sqrt 2 x rms; the project's 100 A-turns at 10 A and 10 turns are peak).
      2. THE BINDING LIMIT in each case: winding or magnets.
      3. THE TIME CONSTANTS and the overloads: how long 1.5x and 2x the continuous current can be carried, from cold and from hot.
      4. SENSITIVITY: the sustainable current when each assumed coefficient is moved to the ends of its range, and when the iron loss or the magnets'
         eddy loss is not what the placeholder says. The binding limit is only as good as these.
      5. The speed: what the maximum operating speed does to the paths that depend on it (gap, end space, windage, bearings).
    It also draws the four figures (thermal_study draws them with the palette of winding_figures.py): the network, the node temperatures at the rated
    point, the sustainable current against stack length, and the transient heating curves.

Run (from the Python to FEMM folder; a few seconds)
    uv run python -m optimization.thermal_study              everything: tables, sensitivity, transients, figures
    uv run python -m optimization.thermal_study --no-figures
Writes out\\optimization\\thermal\\: sustainable.csv, reference_current.csv, time_constants.csv, overload.csv, sensitivity.csv, loss_inputs.csv, speed.csv, study.json,
study_output.txt, and the figures.

Numbers to expect (this PC, 2026-11-09; a run takes about 10 to 20 s, figures included). The magnets bind in all 24 rows, class F and class H alike. Continuous current density
[A/mm^2 rms], copper loss [W] and ampere-turns per coil side [peak] for the stacks of 20 / 40 / 60 / 80 mm:
    natural convection  3.02 / 2.39 / 1.85 / 1.32 A/mm^2;  18.7 / 14.6 / 10.3 / 6.1 W;  70 / 55 / 43 / 30 A-turns
    forced air          6.81 / 6.89 / 6.88 / 6.86 A/mm^2;   99.9 / 124.9 / 147.9 / 170.4 W;  157 / 159 / 159 / 158
    liquid jacket       11.54 / 13.50 / 14.44 / 14.98 A/mm^2;  315 / 528 / 714 / 887 W;  266 / 311 / 333 / 345
At 60 mm the slowest mode is 4,940 s (natural), 1,500 s (forced air) and 980 s (liquid), and the winding's fast mode 551 / 417 / 131 s; the hot spot reaches 63 % of its rise in
4,170 / 720 / 130 s and the magnets in 5,040 / 1,950 / 1,140 s. 2x the continuous current from cold takes the hot spot to class F in 406 s (forced air) and 28 s (liquid); never
(natural). The winding study's reference current (118.6 ampere-turns a coil side, peak: 5.15 A/mm^2) passes the magnets' limit by 48 to 75 K in natural convection (magnets 128 to 155 C;
the winding's hot spot 145 / 150 / 157 / 163 C, so class F too from 60 mm), and leaves 15 to 17 K of magnet margin in forced air (magnets 63 to 65 C) and 32 to 34 K in the jacket (46 to
48 C). The sensitivity run (134 rows) moves 22 coefficients to the ends of their ranges, one at a time, plus the copper-loss module's shorter diamond end turn (48.55 mm against
61.15 mm: J -6 % natural, +1 % forced air, +6 % liquid, magnets binding in all three), and the binding limit is the magnets' in every row (three natural-convection
rows cannot even carry the non-copper heat, J = 0: a bare-metal housing at emissivity 0.1, a rotor iron loss of half the stator's, and every rotor path unfavourable together);
with every rotor-path coefficient at its favourable end together the magnets still bind (J +18 % natural, +5 % forced air, +3.5 % liquid). With a magnet limit of 150 C (SH
grade) the winding binds in natural and forced air (J 5.09 and 9.75), and with 100 C or more in the jacket (J 15.67).

Self-check (the answers follow).
  1. Class H allows 25 K more at the hot spot than class F, yet the sustainable current is the same in most rows of the table. Why?
  2. The sustainable copper loss is much smaller than the total heat the cooling can carry. Which losses use the rest?
  3. Doubling the stack length does not double the sustainable copper loss. Which parts of the network do not grow with the stack?
Answers.
  1. The magnets bind first: they have a hard 80 C limit and the winding's hot spot is still far below 155 C when they reach it, so the class never
     comes into play. Only where the winding binds (the liquid jacket on a short stack) do F and H differ.
  2. The iron loss (a provisional placeholder at 200 Hz), the windage and the bearings: they are heat the cooling must also carry. The sustainable
     TOTAL is a property of the cooling and the limit; the copper gets what is left.
  3. The housing's end bells and end space, the shaft and the bearings: their areas and conduction paths belong to the ends, not the stack.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time

import numpy as np

import paths
from optimization import thermal
from optimization.thermal import (CLASS_F_C, CLASS_H_C, MAGNET_LIMIT_C, Assumptions, Cooling, ThermalModel)

STACKS = (20.0, 40.0, 60.0, 80.0)
COOLINGS = (("natural", Cooling.natural()), ("forced air", Cooling.forced_air()), ("liquid", Cooling.liquid()))
CLASSES = (("F", CLASS_F_C), ("H", CLASS_H_C))
OUT = lambda: paths.path("out") / "optimization" / "thermal"


def sustainable_table() -> list[dict]:
    """The headline table: every cooling x stack x class at 3,000 rpm."""
    rows = []
    for cname, cooling in COOLINGS:
        for stack in STACKS:
            model = ThermalModel(stack_mm=stack, cooling=cooling)
            for klass, limit in CLASSES:
                s = model.sustainable(wind_limit_c=limit)
                r = s["result"]
                non_cu = sum(v for k, v in r.sources.items() if not k.startswith("copper"))
                rows.append({"cooling": cname, "stack_mm": stack, "class": klass, "J_rms_A_mm2": s["J_rms"], "P_cu_W": s["P_cu"], "P_cu_slot_W": r.copper_slot,
                             "P_cu_end_W": r.copper_end, "P_other_W": non_cu, "P_total_W": s["P_total"], "AT_pk_per_coil_side": s["AT_pk"],
                             "AT_rms_per_coil_side": s["AT_rms"], "binding": s["binding"], "J_winding_limit": s["J_winding_limit"],
                             "J_magnet_limit": s["J_magnet_limit"], "T_hot_spot_C": s["hot_spot"], "T_hot_slot_C": s["hot_spot_slot"],
                             "T_hot_end_C": s["hot_spot_end"], "T_magnet_C": s["magnet"], **{f"T_{n}_C": v for n, v in s["T"].items()}})
    return rows


REFERENCE_AT_PEAK = 118.6      # ampere-turns per coil side [peak]: the winding study's reference operating point (the example of winding.py and copper_loss.py)


def reference_point_table() -> list[dict]:
    """What the winding study's reference current does in every cooling and stack: the temperatures at 118.6 ampere-turns per coil side (peak), 3,000 rpm,
    the provisional iron loss and no magnet eddy loss. The same loss model as the sustainable table, with the current given instead of found: the
    margins are the class F limit (155 C) and the magnets' 80 C minus the temperature, so a negative margin is a limit passed."""
    rows = []
    for cname, cooling in COOLINGS:
        for stack in STACKS:
            model = ThermalModel(stack_mm=stack, cooling=cooling)
            j = REFERENCE_AT_PEAK / math.sqrt(2.0) / model.ampere_turns(1.0)["a_cu_side_mm2"]
            r = model.steady_state(model.default_losses(copper=model.copper(j)))
            rows.append({"cooling": cname, "stack_mm": stack, "AT_pk_per_coil_side": REFERENCE_AT_PEAK, "J_rms_A_mm2": j, "P_cu_W": r.copper_loss,
                         "T_hot_spot_C": r.hot_spot, "T_magnet_C": r.magnet, "margin_winding_F_K": r.margin_winding, "margin_magnet_K": r.margin_magnet,
                         "inside_both_limits": bool(r.margin_winding >= 0 and r.margin_magnet >= 0), "runaway": bool(r.runaway)})
    return rows


def dynamics_table(rows: list[dict]) -> tuple[list[dict], list[dict]]:
    """Time constants and overloads at each case's continuous current (class F)."""
    tcs, ovs = [], []
    for row in rows:
        if row["class"] != "F":
            continue
        cooling = dict(COOLINGS)[row["cooling"]]
        model = ThermalModel(stack_mm=row["stack_mm"], cooling=cooling)
        j = row["J_rms_A_mm2"]
        tc = model.time_constants(j)
        tcs.append({"cooling": row["cooling"], "stack_mm": row["stack_mm"], "J_rms": j, "tau_slowest_s": tc["tau_s"][0] if tc["tau_s"] else float("nan"),
                    "tau_list_s": [round(x, 1) for x in tc["tau_s"][:6]], "weights": [round(x, 3) for x in tc["winding_weight"][:6]],
                    "t63_hot_spot_s": tc["t63_hot_spot_s"], "t63_magnet_s": tc["t63_magnet_s"], "rise_hot_spot_K": tc["rise_hot_spot_K"],
                    "rise_magnet_K": tc["rise_magnet_K"]})
        for factor in (1.5, 2.0):
            for start in ("cold", "hot"):
                o = model.overload(factor, j, start=start)
                ovs.append({"cooling": row["cooling"], "stack_mm": row["stack_mm"], "factor": factor, "start": start, "J_cont": j, "J_overload": factor * j,
                            "t_classF_s": o["winding_155C_s"], "t_classH_s": o["winding_180C_s"], "t_magnet_s": o["magnet_80C_s"],
                            "hot_spot_start_C": o["hot_spot_start"], "magnet_start_C": o["magnet_start"], "magnet_rise_60s_K": o["magnet_rise_60s"]})
    return tcs, ovs


# the assumed coefficients, each moved to the ends of its range: (name, low, high)
SENSITIVITY = [
    ("contact_factor", 1.0, 8.0), ("k_matrix", 0.15, 1.0), ("liner_mm", 0.15, 0.40), ("housing_wall_mm", 3.0, 10.0), ("flange_mm2kw", 100.0, 3000.0),
    ("rotor_shaft_mm2kw", 100.0, 2000.0), ("bearing_r_kw", 0.5, 5.0), ("shaft_stub_mm", 4.0, 20.0), ("rotor_k_fraction", 0.5, 1.0),
    ("magnet_gap_mm", 0.05, 0.3), ("k_adhesive", 0.2, 1.5), ("end_area_factor", 0.6, 1.4), ("end_clearance_mm", 3.0, 15.0),
    ("b_tooth_t", 1.3, 1.8), ("b_yoke_t", 1.0, 1.6), ("rotor_iron_fraction", 0.0, 0.5), ("wind_gap_factor", 1.0, 2.0), ("t_gap_film_c", 50.0, 100.0),
    ("emissivity", 0.10, 0.95), ("natural_fin_gain", 1.0, 1.8), ("air_speed_ms", 2.0, 12.0), ("flow_lpm", 1.0, 10.0),
]
H_END_SCALE = (0.5, 2.0)       # the end-space coefficient moved as a whole: both k1 and k2 (and the cap) scaled
# every assumption that decides the rotor's cooling moved to the same end together (the rotor's paths: the shaft, bearings, end space, magnet bond, rotor iron loss)
ROTOR_PATHS = {
    "all rotor paths favourable": dict(shaft_stub_mm=4.0, bearing_r_kw=0.5, rotor_shaft_mm2kw=100.0, rotor_k_fraction=1.0, end_area_factor=1.4, h_end_k1=26.58,
                                       h_end_k2=3.386, h_end_max=200.0, magnet_gap_mm=0.05, k_adhesive=1.5, rotor_iron_fraction=0.0),
    "all rotor paths unfavourable": dict(shaft_stub_mm=20.0, bearing_r_kw=5.0, rotor_shaft_mm2kw=2000.0, rotor_k_fraction=0.5, end_area_factor=0.6, h_end_k1=6.645,
                                         h_end_k2=0.8465, h_end_max=50.0, magnet_gap_mm=0.3, k_adhesive=0.2, rotor_iron_fraction=0.5),
}


def sensitivity_table(stack_mm: float = 60.0) -> list[dict]:
    """One-at-a-time: the sustainable current density (class F) with each assumed coefficient at the ends of its range. Only the coefficients that
    the cooling uses are reported for it (air speed for forced air, flow for the liquid, emissivity and fin gain for natural convection)."""
    rows = []
    relevant = {"natural": {"emissivity", "natural_fin_gain"}, "forced air": {"air_speed_ms"}, "liquid": {"flow_lpm", "housing_wall_mm"}}
    for cname, cooling in COOLINGS:
        base = ThermalModel(stack_mm=stack_mm, cooling=cooling).sustainable()
        row0 = {"cooling": cname, "stack_mm": stack_mm, "parameter": "(base)", "setting": "", "J_rms": base["J_rms"], "binding": base["binding"],
                "J_winding_limit": base["J_winding_limit"], "J_magnet_limit": base["J_magnet_limit"], "change_%": 0.0}
        rows.append(row0)
        trials = [(name, lo, hi) for name, lo, hi in SENSITIVITY if name not in {"emissivity", "natural_fin_gain", "air_speed_ms", "flow_lpm"} or name in relevant[cname]]
        for name, lo, hi in trials:
            for setting in (lo, hi):
                a = Assumptions(**{name: type(getattr(Assumptions(), name))(setting)})
                s = ThermalModel(stack_mm=stack_mm, cooling=cooling, assume=a).sustainable()
                rows.append({"cooling": cname, "stack_mm": stack_mm, "parameter": name, "setting": setting, "J_rms": s["J_rms"], "binding": s["binding"],
                             "J_winding_limit": s["J_winding_limit"], "J_magnet_limit": s["J_magnet_limit"],
                             "change_%": 100.0 * (s["J_rms"] / base["J_rms"] - 1.0)})
        for name, kw in ROTOR_PATHS.items():
            s = ThermalModel(stack_mm=stack_mm, cooling=cooling, assume=Assumptions(**kw)).sustainable()
            rows.append({"cooling": cname, "stack_mm": stack_mm, "parameter": name, "setting": "", "J_rms": s["J_rms"], "binding": s["binding"],
                         "J_winding_limit": s["J_winding_limit"], "J_magnet_limit": s["J_magnet_limit"], "change_%": 100.0 * (s["J_rms"] / base["J_rms"] - 1.0)})
        # the copper-loss module (copper_loss.py) draws a diamond end turn of 48.55 mm per coil side (both ends) against the winding study's half circle (61.15 mm here);
        # its overhang is about 13.6 mm against 20.2 mm
        s = ThermalModel(stack_mm=stack_mm, cooling=cooling, assume=Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6)).sustainable()
        rows.append({"cooling": cname, "stack_mm": stack_mm, "parameter": "end turn: diamond 48.55 mm", "setting": "", "J_rms": s["J_rms"],
                     "binding": s["binding"], "J_winding_limit": s["J_winding_limit"], "J_magnet_limit": s["J_magnet_limit"], "change_%": 100.0 * (s["J_rms"] / base["J_rms"] - 1.0),
                     "P_cu_W": s["P_cu"]})
        for scale in H_END_SCALE:
            a = Assumptions(h_end_k1=13.29 * scale, h_end_k2=1.693 * scale, h_end_max=100.0 * scale)
            s = ThermalModel(stack_mm=stack_mm, cooling=cooling, assume=a).sustainable()
            rows.append({"cooling": cname, "stack_mm": stack_mm, "parameter": "end-space h (all of it)", "setting": scale, "J_rms": s["J_rms"],
                         "binding": s["binding"], "J_winding_limit": s["J_winding_limit"], "J_magnet_limit": s["J_magnet_limit"],
                         "change_%": 100.0 * (s["J_rms"] / base["J_rms"] - 1.0)})
    return rows


def loss_inputs_table(stack_mm: float = 60.0) -> list[dict]:
    """The sustainable current density against the unknown or provisional losses: the magnets' eddy loss (0 ... 20 W), the iron loss (0, 0.5, 1, 2 x the
    placeholder) and the magnet temperature limit (80 N, 100 M, 120 H, 150 SH grade)."""
    rows = []
    for cname, cooling in COOLINGS:
        model = ThermalModel(stack_mm=stack_mm, cooling=cooling)
        ir = model.iron_placeholder()
        for pm in (0.0, 2.0, 5.0, 10.0, 20.0):
            s = model.sustainable(magnet_loss=pm)
            rows.append({"cooling": cname, "stack_mm": stack_mm, "input": "magnet eddy loss [W]", "value": pm, "J_rms": s["J_rms"], "P_cu_W": s["P_cu"],
                         "binding": s["binding"], "T_magnet": s["magnet"], "T_hot_spot": s["hot_spot"]})
        for scale in (0.0, 0.5, 1.0, 2.0):
            iron = {k: (v * scale if k != "f_hz" else v) for k, v in ir.items()}
            s = model.sustainable(iron=iron)
            rows.append({"cooling": cname, "stack_mm": stack_mm, "input": "iron loss x placeholder", "value": scale, "J_rms": s["J_rms"], "P_cu_W": s["P_cu"],
                         "binding": s["binding"], "T_magnet": s["magnet"], "T_hot_spot": s["hot_spot"]})
        for limit, grade in ((80.0, "N"), (100.0, "M"), (120.0, "H"), (150.0, "SH")):
            s = model.sustainable(magnet_limit_c=limit)
            rows.append({"cooling": cname, "stack_mm": stack_mm, "input": f"magnet limit [C] ({grade} grade)", "value": limit, "J_rms": s["J_rms"],
                         "P_cu_W": s["P_cu"], "binding": s["binding"], "T_magnet": s["magnet"], "T_hot_spot": s["hot_spot"]})
    return rows


def speed_table(stack_mm: float = 60.0) -> list[dict]:
    """What the speed does to the paths that depend on it, the iron loss HELD at the 3,000 rpm placeholder (isolating the thermal paths)."""
    rows = []
    for cname, cooling in COOLINGS:
        base_iron = ThermalModel(stack_mm=stack_mm, cooling=cooling).iron_placeholder()
        for rpm in (500.0, 1000.0, 3000.0, 8000.0, thermal.MAX_SPEED_RPM):
            m = ThermalModel(stack_mm=stack_mm, cooling=cooling, speed_rpm=rpm)
            s = m.sustainable(iron=base_iron)
            rows.append({"cooling": cname, "stack_mm": stack_mm, "rpm": rpm, "h_gap": m.gap["h"], "Nu_gap": m.gap["Nu"], "regime": m.gap["regime"],
                         "h_end": m.h_end, "P_windage_W": m.p_windage, "P_bearing_W": m.p_bearing, "J_rms": s["J_rms"], "binding": s["binding"],
                         "T_magnet": s["magnet"]})
    return rows


def _fmt_table(rows: list[dict], cols: list[tuple[str, str, str]]) -> str:
    """A plain-text table: cols = (key, heading, format)."""
    widths = [max(len(h), 8) for _, h, _ in cols]
    lines = ["  ".join(h.rjust(w) for (_, h, _), w in zip(cols, widths))]
    for r in rows:
        cells = []
        for (k, _, f), w in zip(cols, widths):
            v = r[k]
            cells.append((f.format(v) if not isinstance(v, str) else v).rjust(w))
        lines.append("  ".join(cells))
    return "\n".join(lines)


def write_csv(name: str, rows: list[dict]) -> None:
    OUT().mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT() / f"{name}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()} for r in rows])


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.perf_counter()
    OUT().mkdir(parents=True, exist_ok=True)
    lines: list[str] = []

    def say(text: str = "") -> None:
        print(text)
        lines.append(text)
    rows = sustainable_table()
    write_csv("sustainable", rows)
    say("SUSTAINABLE LOSS at 3,000 rpm, 40 C ambient / coolant inlet, copper fill 0.45, iron loss = the PROVISIONAL placeholder, no magnet eddy loss")
    say(_fmt_table(rows, [("cooling", "cooling", ""), ("stack_mm", "stack mm", "{:.0f}"), ("class", "class", ""), ("P_cu_W", "P_cu W", "{:.1f}"),
                          ("J_rms_A_mm2", "J A/mm2", "{:.2f}"), ("AT_pk_per_coil_side", "AT pk", "{:.0f}"), ("AT_rms_per_coil_side", "AT rms", "{:.0f}"),
                          ("binding", "binding", ""), ("T_hot_spot_C", "hot spot", "{:.1f}"), ("T_magnet_C", "magnets", "{:.1f}"), ("P_other_W", "P other W", "{:.1f}"),
                          ("P_total_W", "P total W", "{:.1f}")]))
    ref = reference_point_table()
    write_csv("reference_current", ref)
    say(f"\nTHE WINDING STUDY'S REFERENCE CURRENT ({REFERENCE_AT_PEAK} ampere-turns a coil side, peak) in each cooling: temperatures at that current (a negative margin is a limit passed)")
    for r in ref:
        say(f"  {r['cooling']:10s} {r['stack_mm']:3.0f} mm  J {r['J_rms_A_mm2']:5.2f}  P_cu {r['P_cu_W']:6.1f} W  hot spot {r['T_hot_spot_C']:6.1f} C (margin to F {r['margin_winding_F_K']:+6.1f} K)"
            f"  magnets {r['T_magnet_C']:6.1f} C (margin {r['margin_magnet_K']:+6.1f} K)  {'inside both limits' if r['inside_both_limits'] else 'LIMIT PASSED'}")
    tcs, ovs = dynamics_table(rows)
    write_csv("time_constants", tcs)
    write_csv("overload", ovs)
    say("\nTIME CONSTANTS at each case's continuous current (class F): the modes of C^-1 G above 2 s (slowest first), their weight in the winding's rise, and 63 % rise times from cold")
    for r in tcs:
        say(f"  {r['cooling']:10s} {r['stack_mm']:3.0f} mm  J {r['J_rms']:5.2f}  tau {r['tau_list_s']}  weights {r['weights']}  t63: hot spot {r['t63_hot_spot_s']:.0f} s, magnets {r['t63_magnet_s']:.0f} s"
            f"  (rises {r['rise_hot_spot_K']:.0f} K / {r['rise_magnet_K']:.0f} K)")
    say("\nOVERLOADS: seconds until the winding's hot spot passes class F (155 C), class H (180 C), and until the magnets pass 80 C (None: not within 4 h)")
    for r in ovs:
        def f(x):
            return "never" if x is None else f"{x:.0f}"
        say(f"  {r['cooling']:10s} {r['stack_mm']:3.0f} mm  {r['factor']:.1f}x ({r['J_overload']:5.2f} A/mm2) from {r['start']:4s}:  F {f(r['t_classF_s']):>6s}  H {f(r['t_classH_s']):>6s}  magnets {f(r['t_magnet_s']):>6s}"
            f"   start: hot spot {r['hot_spot_start_C']:.0f} C, magnets {r['magnet_start_C']:.1f} C; magnets +{r['magnet_rise_60s_K']:.2f} K after 60 s")
    sens = sensitivity_table()
    write_csv("sensitivity", sens)
    say("\nSENSITIVITY (60 mm, class F): sustainable J with each assumed coefficient at the ends of its range (change against the base)")
    for cname, _ in COOLINGS:
        say(f"  {cname}:")
        for r in [x for x in sens if x["cooling"] == cname]:
            flag = "" if r["parameter"] == "(base)" else ("   <- binding limit changes" if r["binding"] != next(x["binding"] for x in sens if x["cooling"] == cname and x["parameter"] == "(base)") else "")
            say(f"    {r['parameter']:26s} {str(r['setting']):>8s}  J {r['J_rms']:6.2f}  ({r['change_%']:+6.1f} %)  binding {r['binding']}{flag}")
    li = loss_inputs_table()
    write_csv("loss_inputs", li)
    say("\nTHE LOSS INPUTS (60 mm, class F): sustainable J against the magnets' eddy loss, the iron loss and the magnet limit")
    for r in li:
        say(f"  {r['cooling']:10s} {r['input']:28s} {r['value']:7.1f}  J {r['J_rms']:6.2f}  P_cu {r['P_cu_W']:7.1f} W  binding {r['binding']:7s}  magnets {r['T_magnet']:.1f} C  hot spot {r['T_hot_spot']:.1f} C")
    sp = speed_table()
    write_csv("speed", sp)
    say("\nSPEED (60 mm, class F, iron loss held at its 3,000 rpm placeholder): the speed-dependent paths")
    for r in sp:
        say(f"  {r['cooling']:10s} {r['rpm']:7.0f} rpm  gap h {r['h_gap']:6.1f} (Nu {r['Nu_gap']:.2f}, {r['regime']})  end-space h {r['h_end']:5.1f}  windage {r['P_windage_W']:.2f} W  bearings {r['P_bearing_W']:.2f} W"
            f"  J {r['J_rms']:.2f} ({r['binding']})")
    (OUT() / "study.json").write_text(json.dumps({"sustainable": rows, "time_constants": tcs, "overload": ovs, "sensitivity": sens, "loss_inputs": li, "speed": sp, "reference_current": ref},
                                                 indent=1, default=float), encoding="utf-8")
    if not args.no_figures:
        from optimization import thermal_figures
        paths_out = thermal_figures.make_all(rows)
        say("\nfigures: " + ", ".join(str(p) for p in paths_out))
    say(f"\nstudy done in {time.perf_counter() - t0:.1f} s")
    (OUT() / "study_output.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
