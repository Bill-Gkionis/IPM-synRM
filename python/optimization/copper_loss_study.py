"""optimization.copper_loss_study: the answers of the copper-loss work, from the fast model (optimization/copper_loss.py): where the proximity loss matters in this motor and where it does not.

What it is, and why
    The slot model (ac_slot.py), the machine's field (ac_machine.py) and the PWM ripple (ac_pwm.py) are the evidence; copper_loss.py is their one fast summary. This module uses it the way the winding
    optimizer will, and asks the design question: for the optimized motor, 3,000 rpm to the 16,330 rpm stress limit, does the AC loss of the copper matter, which part of it, for which winding?
    Five things, all printed and written to out/optimization/copper_ac/final/:
      1. THE MODEL AGAINST FEMM. The slot's R_ac / R_dc from copper_loss.slot_r_ratio against every FEMM slot solve of the sweep stage (two windings, 13 frequencies, three layouts) and of the val stage
         (six more strands and fills): the largest and the mean error. The surrogate is checked on what made it (the sweep) and on what did not (the val stage).
      2. THE LOSS OVER THE SPEED RANGE: for the winding study's example (W1: 26 turns x 8 x 0.315 mm, 4 paths) and its thick-strand twin (W2: 13 x 4 x 0.630 mm, 2 paths), at 3,000, 6,000, 10,000 and
         16,330 rpm at the rated 118.6 ampere-turns, MTPA angle: every part of the loss and its share of the DC loss.
      3. THE PROXIMITY LOSS over the strand size and the speed (at a constant copper fill), at the rated current and at twice it, as a map.
      4. THE PWM RIPPLE: for the six windings of the winding study's inverter table (10 ... 100 A limits), at 10 and 20 kHz on 48 and 60 V: the ripple, its loss, its share.
      5. THE CIRCULATING CURRENTS: what each placement of the strands in hand costs, over the speed range, for W1 and W2.
      6. THE VOLTAGE LIMIT: the speed at which the 48 V bus can no longer make the winding's voltage at the rated current and the MTPA angle.
      7. THE FIGURES (optimization/copper_loss_figures.py): the loss parts over the speed, the proximity map, the PWM bars, the circulating currents.
    The operating points above the corner speed are in field weakening in reality; the table takes the rated current at the MTPA angle at every speed (an upper bound on the armature's field) and, in a second
    table, the field-weakening angles 45 and 70 degrees; the voltage need is reported ('modulation_index'; above 1 the inverter cannot make it, and the ripple is at 1).

Run
    uv run python -m optimization.copper_loss_study            (no FEMM: 2 s, after the studies and copper_loss_tables.json; also draws the four figures of copper_loss_figures.py)

Numbers to expect (the run prints them all; the winding study's example W1 and its thick-strand twin W2, 118.6 ampere-turns a coil side, 100 degC, 20 mm stack, 20 kHz, 48 V)
    1. the model's slot R_ac / R_dc against FEMM, the extra loss over DC at f >= 1 kHz: W1 with the layers in phase +4.5 % (+0.3 ... +7.7 % over the three layouts: the layouts' own spread of the slot
       constant K is 3.4 %), W2 -3.4 % at low frequency and -2.4 % at d / delta 1.2, then +3.3 %, +10.7 % and +21.6 % at d / delta 2.0, 2.6 and 3.7 (60 ... 200 kHz); six other windings (not used to
       make the model): within 3.2 % up to d / delta 2.1 at 5, 20 and 100 kHz, and +28.6 % at d / delta 4.2 (outside the range of any winding here).
    2. AC loss over DC (skin + proximity + PWM), MTPA angle 23.0 deg: W1 0.031 % at 3,000 rpm, 0.073 % at 6,000, 0.175 % at 10,000 and 0.44 % at 16,330 rpm (the proximity loss 0.0069 ... 0.206 W
       of 48.5 W of DC loss, 14.16 W in the slots and 34.37 W in the end turns); W2 0.10, 0.27, 0.68 and 1.73 % (proximity 0.028 ... 0.82 W).
    3. the map at 16,330 rpm, proximity + skin over DC at the rated ampere-turns: 0.17 % for 0.2 mm strands, 1.7 % for 0.63 mm, 8.2 % for 1.4 mm (0.20, 1.97 and 9.5 % at twice the rated).
    4. the six windings of the inverter table: ripple 0.043 ... 1.26 A rms at 20 kHz, its loss 0.006 ... 0.61 % of the DC loss (0.014 ... 0.87 % at 10 kHz), the largest for the 4-turn 100 A winding.
    5. strands in parallel at random: circulating loss 0.67 % of DC at 3,000 rpm and 13.8 % at 16,330 rpm for W1 (W2 0.57 and 11.7 %); side by side in a row 0.012 and 0.27 % (W2 0.009 and 0.21 %).
    6. the 48 V bus reaches modulation index 1 at 3,320 rpm for both windings at this current and angle: above it the tables of 2 are hypothetical (field weakening).

Self-check (the answers follow).
  1. W2 has the same DC loss as W1 (same copper, same ampere-turns). Why is its proximity loss four times larger?
  2. A winding's p_pwm is 1 % of its DC loss at 20 kHz. What does doubling the switching frequency do?
Answers.
  1. The proximity loss of a strand goes as d^4 and a given copper area holds a number of strands that goes as 1 / d^2: the winding's loss goes as d^2, and W2's strands are twice as thick.
  2. The ripple halves (a quarter of the power), the strands' AC resistance rises (as f^2 while the strands are thin): the DC-like part falls to a quarter, the proximity part of the ripple loss does not change
     (the ripple squared falls as f^-2, the extra resistance rises as f^2). So p_pwm falls, but by less than the factor of 4, as far as proximity dominates it. The table of 4 shows it: halving f_sw
     from 20 to 10 kHz raises the loss share from 0.0058 to 0.0141 % for the 10 A winding (2.4 times) and from 0.61 to 0.87 % for the 100 A winding (1.4 times).
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

import paths
from optimization import ac_theory, copper_loss as cl, copper_loss_figures as figs

OUT = figs.OUT
SLOT = paths.path("out") / "optimization" / "copper_ac" / "slot"
AT_RATED = figs.AT_RATED
V_DC, F_SW = 48.0, 20000.0
SPEEDS = (3000.0, 6000.0, 10000.0, 16330.0)
GAMMA_MTPA = figs.mtpa_angle(AT_RATED)                 # deg: 23.02 at the rated 118.6 ampere-turns, from the angle search of ac_machine.py (the brief's 22.47 is 0.55 deg away: the torque is flat)
W1, W2, INVERTER_TABLE, i_for = figs.W1, figs.W2, figs.INVERTER_TABLE, figs.i_for       # the two example windings, the winding study's inverter table, N i / a = ampere-turns


def row(w: dict, rpm: float, gamma: float, at: float, f_sw: float = F_SW, v_dc: float = V_DC, arrangement: str = "transposed", t_c: float = 100.0, stack: float = 20.0) -> dict:
    return cl.losses(dict(w, arrangement=arrangement), {"i_peak": i_for(w, at), "gamma_deg": gamma, "rpm": rpm, "f_sw": f_sw, "v_dc": v_dc}, t_c, stack)


# ---------------------------------------------------------------------------------------------- 1. the model against FEMM
def check_slot_surrogate(log=print) -> dict:
    sweep = json.loads((SLOT / "sweep.json").read_text(encoding="utf-8"))
    rho = float(ac_theory.resistivity(100.0))
    errs = []
    for r in sweep["rows"]:
        w = sweep["windings"][r["winding"]]
        pred = float(cl.slot_r_ratio(w["d"] * 1e-3, r["fill"], r["f_Hz"], rho, "same" if r["case"] == "same" else "diff"))
        errs.append((r["winding"], r["case"], r["f_Hz"], r["R"], pred, (pred - 1) / (r["R"] - 1) - 1 if r["R"] > 1 + 1e-9 else float("nan")))
    out = {"sweep": errs}
    log("\n1. THE SURROGATE AGAINST FEMM (the slot's R_ac / R_dc: the model's extra loss over FEMM's)")
    for wn in ("W1", "W2"):
        for case in ("same", "diff"):
            sel = [e for e in errs if e[0] == wn and e[1] == case and e[2] >= 1000.0 and not math.isnan(e[5])]
            if sel:
                vals = np.array([e[5] for e in sel])
                log(f"   {wn} {case:4s}: f >= 1 kHz: mean error of the extra loss {100 * vals.mean():+6.2f} %, range {100 * vals.min():+6.2f} ... {100 * vals.max():+6.2f} % over {len(sel)} solves")
    vp = SLOT / "val.json"
    if vp.exists():
        val = json.loads(vp.read_text(encoding="utf-8"))
        ve = []
        for r in val:
            d = {"V1: 20 x 13 x 0.280": 0.280, "V2: 27 x 19 x 0.200": 0.200, "V3: 14 x 29 x 0.224": 0.224, "V4: 20 x 7 x 0.315 (fill 0.30)": 0.315, "V5: 9 x 9 x 0.500": 0.500, "V6: 5 x 4 x 1.000": 1.000}[r["winding"]]
            pred = float(cl.slot_r_ratio(d * 1e-3, r["fill"], r["f_Hz"], rho, "same"))
            ve.append((r["winding"], r["f_Hz"], r["d_over_delta"], r["R"], pred, (pred - 1) / (r["R"] - 1) - 1))
        out["val"] = ve
        log("   the other strands and fills (not used to make the model):")
        for e in ve:
            log(f"     {e[0]:32s} {e[1]:8g} Hz  d/delta {e[2]:6.3f}  FEMM R_ac/R_dc {e[3]:8.4f}  model {e[4]:8.4f}  error of the extra loss {100 * e[5]:+6.1f} %")
    return out


# ---------------------------------------------------------------------------------------------- 2. the loss over the speed range
def speed_tables(log=print) -> dict:
    res = {}
    log("\n2. THE LOSS OVER THE SPEED RANGE (118.6 ampere-turns a coil side, 20 mm stack, copper at 100 degC; 20 kHz, 48 V bus; W)")
    for name, w in (("W1: 26 turns x 8 x 0.315 mm, 4 paths", W1), ("W2: 13 turns x 4 x 0.630 mm, 2 paths", W2)):
        for gamma, label in ((GAMMA_MTPA, "MTPA angle"), (45.0, "field weakening 45 deg"), (70.0, "field weakening 70 deg")):
            log(f"  {name}, {label}")
            log("     rpm    f_e [Hz]  d/delta   DC slots   DC ends     skin      prox       pwm   (ripple A)   (m)    AC / DC [%]   total [W]")
            rows = []
            for rpm in SPEEDS:
                r = row(w, rpm, gamma, AT_RATED)
                dc = r["p_dc_slot"] + r["p_dc_end"]
                ac = r["p_skin"] + r["p_prox"] + r["p_pwm"] + r["p_circ"]
                rows.append({"rpm": rpm, **r})
                log(f"  {rpm:7.0f} {r['f_e_hz']:9.1f}  {r['d_over_delta_fund']:7.4f} {r['p_dc_slot']:9.3f} {r['p_dc_end']:9.3f} {r['p_skin']:9.5f} {r['p_prox']:9.4f} {r['p_pwm']:9.5f}   "
                    f"({r['ripple_rms_a']:6.3f})  {r['modulation_index']:5.2f}{'*' if r['overmodulated'] else ' '}  {100 * ac / dc:10.4f}   {r['p_total']:9.3f}")
            res[(name, label)] = rows
    log("  (* the inverter cannot make this voltage at that speed: the ripple is taken at modulation index 1; above the corner speed the real operating point is in field weakening)")
    return res


# ---------------------------------------------------------------------------------------------- 3. the proximity loss map
def proximity_map(log=print) -> dict:
    """Proximity + skin loss over DC at constant fill 0.448 for strand diameters 0.2 ... 1.4 mm (strands in hand continuous), 13 turns a coil side, 2 paths, over the speed, at 1 x and 2 x rated."""
    ds = np.array([0.2, 0.28, 0.4, 0.5, 0.63, 0.8, 1.0, 1.25, 1.4])
    speeds = np.linspace(1000.0, 16330.0, 12)
    out = {"d_mm": ds, "rpm": speeds}
    log("\n3. THE PROXIMITY + SKIN LOSS over the DC loss [%] at the MTPA angle of each current, constant copper fill 0.448 (13 turns a coil side, 2 paths)")
    for fac in (1.0, 2.0):
        n_s = 0.448 * 36.218 / (13 * math.pi * ds ** 2 / 4.0)
        r = cl.losses_arrays(13, 2, n_s[:, None], ds[:, None], ds[:, None] * 1.12, fac * AT_RATED * 2 / 13, figs.mtpa_angle(fac * AT_RATED), speeds[None, :], 0.0, 48.0, 100.0, 20.0)
        share = 100.0 * (r["p_prox"] + r["p_skin"]) / (r["p_dc_slot"] + r["p_dc_end"])
        out[f"share_{fac:g}x"] = share
        log(f"  {fac:g} x rated ampere-turns:  rows d [mm] = " + " ".join(f"{d:6.3f}" for d in ds))
        for j in (0, 3, 7, 11):
            log(f"     {speeds[j]:7.0f} rpm: " + " ".join(f"{share[i, j]:6.3f}" for i in range(len(ds))) + "   %")
    return out


# ---------------------------------------------------------------------------------------------- 4. PWM
def pwm_table(log=print) -> list[dict]:
    log("\n4. THE PWM RIPPLE for the winding study's inverter table (rated 118.6 ampere-turns, MTPA angle, 3,000 rpm; loss of the ripple over the DC loss)")
    log("   limit [A]   N  a  n_s   d [mm]   I peak [A]   L [uH]   f_sw 10 kHz: ripple rms [A]   share [%]   | 20 kHz: ripple [A]   share [%]   (48 V bus; 60 V: 20 kHz ripple, share)")
    rows = []
    for lim, N, a, n_s, d in INVERTER_TABLE:
        w = {"N": N, "a": a, "n_s": n_s, "d_mm": d, "d_over_mm": figs.d_over_for(d)}
        cells = []
        for f_sw, v in ((10000.0, 48.0), (20000.0, 48.0), (20000.0, 60.0)):
            r = row(w, 3000.0, GAMMA_MTPA, AT_RATED, f_sw, v)
            cells.append((r["ripple_rms_a"], 100.0 * r["p_pwm"] / (r["p_dc_slot"] + r["p_dc_end"]), r))
        rows.append({"limit": lim, "winding": w, "ripple_10k": cells[0][0], "share_10k": cells[0][1], "ripple_20k": cells[1][0], "share_20k": cells[1][1], "ripple_20k_60V": cells[2][0],
                     "share_20k_60V": cells[2][1], "i_peak": i_for(w), "l_uH": 1e6 * cells[1][2]["l_phase_h"]})
        log(f"   {lim:6d}  {N:3d} {a:2d} {n_s:3d}  {d:7.3f}   {i_for(w):9.2f}   {1e6 * cells[1][2]['l_phase_h']:7.1f}   {cells[0][0]:12.3f}  {cells[0][1]:12.4f}   | {cells[1][0]:10.3f}  {cells[1][1]:12.4f}   ({cells[2][0]:.3f}, {cells[2][1]:.4f})")
    return rows


# ---------------------------------------------------------------------------------------------- 5. circulating currents
def circulating_table(log=print) -> dict:
    out = {}
    log("\n5. CIRCULATING CURRENTS between the strands of a turn (118.6 ampere-turns, MTPA angle; W over the DC loss [%])")
    for name, w in (("W1", W1), ("W2", W2)):
        log(f"  {name}:   rpm    " + "   ".join(f"{arr:>10s}" for arr in cl.ARRANGEMENTS[1:]))
        for rpm in SPEEDS:
            cells = []
            for arr in cl.ARRANGEMENTS[1:]:
                r = row(w, rpm, GAMMA_MTPA, AT_RATED, arrangement=arr)
                cells.append(100.0 * r["p_circ"] / (r["p_dc_slot"] + r["p_dc_end"]))
            out[(name, rpm)] = cells
            log(f"        {rpm:7.0f}   " + "   ".join(f"{c:10.4f}" for c in cells))
    return out


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    cl.tables()
    check_slot_surrogate()
    speed_tables()
    proximity_map()
    pwm_table()
    circulating_table()
    print("\n6. THE VOLTAGE LIMIT of the 48 V bus at the rated current and the MTPA angle (above it the machine is in field weakening and the tables above are hypothetical):")
    for name, w in (("W1", W1), ("W2", W2)):
        c = figs.corner_rpm(w, GAMMA_MTPA)
        print(f"   {name}: the inverter reaches modulation index 1 at {c:.0f} rpm" if c else f"   {name}: modulation index stays below 1 up to 16,330 rpm")
    print("\n7. THE FIGURES")
    for name, fn in (("fig_final_speed.png", figs.fig_speed), ("fig_final_map.png", figs.fig_map), ("fig_final_pwm.png", figs.fig_pwm), ("fig_final_circ.png", figs.fig_circ)):
        fn(OUT / name)
        print(f"   written {OUT / name}")
    print(f"\ndone in {time.perf_counter() - t0:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
