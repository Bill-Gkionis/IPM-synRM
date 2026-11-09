"""optimization.thermal_fan_study: the shaft fan: its model and sizing, the continuous limit against speed, where the heat goes, and the capability of each cooling.

What it is, and why
    The cooling was decided: forced air from a fan on the shaft, like a simple totally enclosed motor (liquid cooling is parked). A shaft fan is not the fixed-speed
    fan of the first study: its air speed is proportional to the shaft speed, so it cools least where the motor needs it most (low speed, high torque, the largest copper
    loss) and its power grows with the cube of the speed. This file asks the network (thermal.py, Cooling.shaft_fan) and the motor's own loss and torque models
    (thermal_machine.py: flux maps, iron-loss map, copper-loss model) five questions, and writes the answers down. No solver is started.
      1. THE FAN. The air speed v(n) = v_ref n / n_ref in the fin channels of the forced-air preset (36 fins, 15 x 3 mm), the pressure the channels and the cowl ask for,
         the air power and the shaft power over a fan efficiency (P ~ n^3). Two sizings: (a) 6 m/s at 3,000 rpm (the preset's speed), and (b) the fan whose power at the
         top speed, 16,330 rpm, is held to a budget (50 W), which makes it a much weaker fan at 3,000 rpm. The budget is chosen below with its justification. A table of
         budgets shows what each buys.
      2. THE CONTINUOUS LIMIT AGAINST SPEED, from standstill to 16,330 rpm, for the shaft fan (both sizings), natural convection, the fixed-speed forced air and, for
         reference, the liquid jacket: the ampere-turns per coil side (peak) at which the winding's hot spot (class F) or the magnets (N45 80 C; also the 150 C grade)
         reach their limit, at MTPA currents, with the real iron loss at that point, the copper-loss model as the network's copper function, windage and bearings.
         LABEL: MTPA at every speed. Above the corner speed the real motor weakens the field and has less iron loss; the motor map handles that. The low-speed part is exact.
      3. WHERE THE HEAT GOES at each cooling's continuous limit (60 mm): into which parts, and out through which surfaces.
      4. THE CAPABILITY TABLE per cooling at 3,000 rpm, 40 C ambient: ampere-turns, current density, MTPA torque from the flux maps, the mechanical power T x omega (an
         ESTIMATE), the losses and efficiency, the temperatures, the binding limit, and the short-term torque at 2x the current and how long it lasts from cold. For
         every stack (20, 40, 60, 80 mm) and for both magnet grades.
      5. The sensitivity of the shaft fan's numbers to its assumed coefficients.

The estimates and their assumptions (every power in this file is an estimate)
    T x omega is the MTPA torque of the flux maps (the d-q torque, scaled linearly with the stack: the 2D maps have no end effects) times the shaft speed. It assumes
    sinusoidal currents at the MTPA angle, ignores the inverter's voltage limit (the number of turns is free: at a given ampere-turn the voltage scales with the turns),
    and is the air-gap (electromagnetic) power. Efficiency = (T omega - windage - bearings - shaft fan) / (T omega + copper + iron [+ the external fan of the forced_air
    preset]); it leaves out the inverter. The copper loss is the copper-loss model's (DC plus AC, example winding, 20 kHz, 48 V) at the copper temperature the network
    finds; the iron loss the iron-loss map's at the MTPA point; windage and bearings the network's estimates; no magnet eddy loss (an input, 0 here).
    The 150 C grade is an open decision: the maps are those of the N45 magnets; a 150 C grade has a lower remanence (a few per cent; check before citing), so its torque
    would be a little lower than shown.

Run (from the Python to FEMM folder; 3 minutes for everything, 40 s with --quick --no-sensitivity)
    uv run python -m optimization.thermal_fan_study               everything and the figures
    uv run python -m optimization.thermal_fan_study --no-figures
    uv run python -m optimization.thermal_fan_study --quick       the 60 mm stack only (a third of the work)
    uv run python -m optimization.thermal_fan_study --figures-only   redraw the figures from the saved fan_study.json
Writes out\\optimization\\thermal\\: fan_model.csv, fan_sizing.csv, fan_speed_limit.csv, fan_heat_paths.csv, fan_capability.csv, fan_sensitivity.csv, fan_study.json,
fan_study_output.txt, fan_preheat.csv, and the four figures (thermal_fan_figures.py).

Numbers to expect (this PC, 2026-11-09; 182 s for everything; 60 mm stack, N45 80 C, class F, 40 C ambient, MTPA, the example winding)
    The fan: (a) 6 m/s at 3,000 rpm takes 12.4 W of shaft power at 3,000 rpm, 449 W at 10,000 and 1,941 W at 16,330 (it warms its air by 0.34 / 3.7 / 9.8 K); (b), held to 50 W at 16,330 rpm,
    blows 9.58 m/s at the top speed and 1.76 m/s at 3,000 rpm, where it takes 0.32 W. The budget is 2.2 % of the 2,297 W (an estimate) that (a) allows at 3,000 rpm. The fan's own heat in
    its air costs (a) 2.1 % of its ampere-turns at 6,000 rpm, 8.3 % at 10,000 (142.4 against 155.2) and 63 % at 16,330 (33.1 against 88.3); (b) 0.3 and 1.7 % at 6,000 and 10,000 rpm.
    The limit against speed, ampere-turns a coil side (peak), at 0 / 3,000 / 10,000 / 16,330 rpm: natural convection 80.8 / 48.6 / none / none; shaft fan (a) 95.1 / 158.1 / 142.4 / 33.1;
    shaft fan (b) 95.1 / 96.6 / 59.1 / none; fixed-speed forced air 176.1 / 157.8 / 62.3 / none; liquid jacket 355.9 (the flux maps' edge) / 350.5 / 286.4 / 183.8. 'None': the iron,
    windage and bearing heat alone passes the magnets' 80 C (natural convection from 5,000 rpm, fixed forced air from 12,000). At standstill the shaft fan has stopped and carries 54 % of
    what the fixed-speed fan carries; with the 150 C grade the limits are 132 / 117 / none / none, 152 / 228 / 260 / 222 (a), 152 / 164 / 177 / 144 (b), 234 / 226 / 177 / 102 and 356 / 356 / 356 / 331.
    Where the heat goes (3,000 rpm): natural convection leaves 15 % by convection and 22 % by radiation from the housing and 63 % from the end bells; the fixed forced air and (a) leave 62 %
    through the housing and 38 % through the bells, (b) 53 and 47 %; the slot winding gives 76 to 78 % of its heat to the teeth and 22 to 24 % to the yoke (69 and 31 % in natural convection,
    where the iron loss warms the teeth); the end turns send 24 to 31 % of the heat made through the slot (11 % in natural convection, 39 % in the jacket) and 4 to 5 % into the end-space air;
    the stator heats the rotor across the air gap by 1 to 2 % of the total, and the rotor gives it to the end-space air and through the shaft.
    The capability at 3,000 rpm, 60 mm: continuous ampere-turns, MTPA torque, T x omega (estimate), efficiency: natural 48.6, 2.15 N.m, 0.67 kW, 95.2 %; (a) 158.1, 7.31 N.m, 2.30 kW, 92.8 %;
    (b) 96.6, 4.37 N.m, 1.37 kW, 94.9 %; fixed forced air 157.8, 7.30 N.m, 2.29 kW, 92.9 %; liquid 350.5, 15.10 N.m, 4.74 kW, 86.3 %. At 2x the current from cold the first limit comes after
    3,563 s (natural, magnets), 462 s ((a), winding), 1,527 s ((b), magnets), 466 s (fixed, winding), 25 s (liquid, winding; 2x is beyond the flux maps). The 150 C grade raises the
    continuous torque of natural convection to 5.36 N.m and of the fans to 10.54 N.m ((a)) and 7.57 N.m ((b)).

Self-check (the answers follow).
  1. The shaft fan blows only 0.5 m/s at 250 rpm and already carries 19 % more current than at standstill. Why does so little air help so much?
  2. Sizing (a) carries less current at 16,330 rpm (33 A-turns) than at 6,000 rpm (169), although its air is five times faster. Name the two reasons.
  3. Why is the fan's budget a limit on peak power and not on energy, and what does a larger budget buy?
Answers.
  1. The fins are a large area (3.5 times the plain housing) and even laminar flow between them (Nu = 4, h = 9 W/m2K) beats the natural convection of the same fins; the stream's capacity
     m cp (T_wall - T_air), 3 W/K at 0.5 m/s, is what limits it, and that is already more than the 0.6 W/K of still-air convection.
  2. The iron loss at MTPA currents grows with the speed (284 W at the limit against 88 W at 6,000 rpm), and the fan's own 1.9 kW warms its air by 9.8 K, which keeps about
     17 W/K x 9.8 K = 169 W of cooling from the housing; at 6,000 rpm the fan costs 98 W and warms the air by 1.35 K.
  3. P goes as n^3, so the power is large only near the top speed: a 200 W fan at 16,330 rpm costs 10 W at 6,000 rpm and 1.2 W at 3,000. A larger budget buys current between 6,000 and
     10,000 rpm (118 A-turns at 6,000 rpm for 200 W against 56 for 50 W) and nothing at 3,000 rpm (97.4 against 96.6).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import time

import numpy as np

import paths
from optimization import thermal, thermal_machine as tm
from optimization.thermal import CLASS_F_C, MAX_SPEED_RPM, RATED_RPM, Cooling

STACKS = (20.0, 40.0, 60.0, 80.0)
MAIN_STACK = 60.0
SPEEDS = (0, 250, 500, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 8000, 10000, 12000, 14000, MAX_SPEED_RPM)
FAN_SPEEDS = (0, 500, 1000, 2000, 3000, 5000, 8000, 10000, 12000, MAX_SPEED_RPM)
BUDGET_W = 50.0                                          # the shaft fan's power at 16,330 rpm in sizing (b)
BUDGETS_W = (10.0, 25.0, 50.0, 100.0, 200.0, 400.0, 800.0, 1600.0)
TRADE_SPEEDS = (1000, 3000, 6000, 10000, MAX_SPEED_RPM)
MAGNETS = (("N45", 80.0), ("150C", 150.0))
KEYS = ("natural", "shaft_a", "shaft_b", "forced", "liquid")
NAMES = {"natural": "natural convection (plain housing)", "shaft_a": "shaft fan (a): 6 m/s at 3,000 rpm", "shaft_b": f"shaft fan (b): {BUDGET_W:.0f} W at 16,330 rpm",
         "forced": "fixed-speed forced air (6 m/s)", "liquid": "liquid jacket (reference)"}
HEAT_CASES = (("natural", 3000), ("shaft_a", 1000), ("shaft_a", 3000), ("shaft_a", 10000), ("shaft_b", 1000), ("shaft_b", 3000), ("shaft_b", 10000),
              ("forced", 3000), ("liquid", 3000))
BIND_CHAR = {"winding": "w", "magnet": "m"}
OUT = lambda: paths.path("out") / "optimization" / "thermal"


def coolings(stack_mm: float, assume) -> dict:
    """The five coolings for one stack. Sizing (b) is made for the stack (the finned length sets the channel's friction a little)."""
    b = thermal.size_shaft_fan(BUDGET_W, stack_mm=stack_mm, assume=assume)
    return {"natural": Cooling.natural(), "shaft_a": Cooling.shaft_fan(6.0, float(RATED_RPM)), "shaft_b": b["cooling"], "forced": Cooling.forced_air(), "liquid": Cooling.liquid()}


# ---------------------------------------------------------------------------------------------------- 1. the fan
def fan_table(mach, assume, stack_mm: float = MAIN_STACK) -> list[dict]:
    """For the two shaft fans at each speed: the air speed, the channel's Re and friction, the pressure, the flow, the air and shaft power, the warming of the air by the
    fan's loss, and the outside's conductance (housing + bells, to the ambient at 80 C surfaces) next to the same for natural convection."""
    cs = coolings(stack_mm, assume)
    rows = []
    for key in ("shaft_a", "shaft_b"):
        for rpm in FAN_SPEEDS:
            m = tm.model_for(cs[key], stack_mm, rpm, mach, assume)
            fs = m.fan_state()
            p = m.outside_parts(80.0, 80.0)
            nat = tm.model_for(cs["natural"], stack_mm, rpm, mach, assume).outside_parts(80.0, 80.0)
            rows.append({"sizing": key, "v_ref_ms": cs[key].v_ref_ms, "rpm": rpm, "v_ms": fs["v"], "re": fs["re"], "f_darcy": fs["f_darcy"], "dp_pa": fs["dp_pa"],
                         "q_m3s": fs["q_m3s"], "p_air_w": fs["p_air_w"], "p_shaft_w": fs["p_shaft_w"], "dT_air_k": fs["dT_air_k"],
                         "g_housing_w_per_k": p["h_conv"] + p["h_rad"], "g_bells_w_per_k": p["b_conv"] + p["b_rad"],
                         "g_outside_w_per_k": p["h_conv"] + p["h_rad"] + p["b_conv"] + p["b_rad"],
                         "g_outside_natural_plain_w_per_k": nat["h_conv"] + nat["h_rad"] + nat["b_conv"] + nat["b_rad"]})
    return rows


def sizing_table(mach, assume, stack_mm: float = MAIN_STACK) -> list[dict]:
    """The budgets of the fan's power at 16,330 rpm: the air speed they give at the top speed and at 3,000 rpm, the power at 3,000 rpm and at the top, and the continuous
    ampere-turns (N45, class F, real losses) at five speeds. Budget 0 is the finned housing with the fan removed (still-air convection at every speed); the last row
    is sizing (a) (6 m/s at 3,000 rpm), whose power at the top speed is what it is."""
    rows = []
    a_top = thermal.shaft_fan_power(MAX_SPEED_RPM, Cooling.shaft_fan(6.0, float(RATED_RPM)), stack_mm, assume=assume)
    for label, budget in [("no fan", 0.0)] + [(f"{b:g} W", b) for b in BUDGETS_W] + [("(a) 6 m/s at 3,000 rpm", a_top)]:
        if budget == 0.0:
            c = Cooling.shaft_fan(0.0, float(RATED_RPM))
            v_top, v_ref = 0.0, 0.0
        elif label.startswith("(a)"):
            c = Cooling.shaft_fan(6.0, float(RATED_RPM))
            v_ref = 6.0
            v_top = 6.0 * MAX_SPEED_RPM / RATED_RPM
        else:
            s = thermal.size_shaft_fan(budget, stack_mm=stack_mm, assume=assume)
            c, v_ref, v_top = s["cooling"], s["v_ref"], s["v_top"]
        row = {"label": label, "budget_w": budget, "v_ref_ms": v_ref, "v_top_ms": v_top}
        for rpm in TRADE_SPEEDS:
            m = tm.model_for(c, stack_mm, rpm, mach, assume)
            cl = tm.continuous_limit(m, mach)
            row[f"AT_{rpm}"] = cl["AT_pk"]
            row[f"feasible_{rpm}"] = cl["feasible"]
            row[f"P_fan_{rpm}_w"] = m.fan_power()
        rows.append(row)
    return rows


def preheat_rows(mach, assume, stack_mm: float = MAIN_STACK) -> list[dict]:
    """What the fan's own heat does: the continuous ampere-turns of the two sizings at three speeds with the stream warmed by the fan's loss and without it."""
    cs = coolings(stack_mm, assume)
    off = tm.Assumptions(**{**assume.__dict__, "fan_air_preheat": False})
    rows = []
    for key in ("shaft_a", "shaft_b"):
        for rpm in (6000, 10000, MAX_SPEED_RPM):
            m_on = tm.model_for(cs[key], stack_mm, rpm, mach, assume)
            on = tm.continuous_limit(m_on, mach)
            offv = tm.continuous_limit(tm.model_for(cs[key], stack_mm, rpm, mach, off), mach)
            rows.append({"sizing": key[-1], "rpm": rpm, "p_fan_w": m_on.fan_power(), "dT_air_k": m_on.dT_fan, "AT_with": on["AT_pk"], "AT_without": offv["AT_pk"]})
    return rows


# ---------------------------------------------------------------------------------------------------- 2. the continuous limit against speed
def speed_rows(mach, assume, stacks=STACKS) -> list[dict]:
    rows = []
    for stack in stacks:
        cs = coolings(stack, assume)
        for key in KEYS:
            for mag, mlim in MAGNETS:
                for rpm in SPEEDS:
                    m = tm.model_for(cs[key], stack, rpm, mach, assume)
                    cl = tm.continuous_limit(m, mach, magnet_limit_c=mlim)
                    rows.append(flat(key, mag, cl))
    return rows


def flat(key: str, mag: str, cl: dict) -> dict:
    """One row for a CSV from a continuous_limit description."""
    binding = cl["binding"] if (cl["feasible"] or cl["limited_by_maps"]) else "none: limit passed at zero current"
    return {"cooling": key, "magnets": mag, "stack_mm": cl["stack_mm"], "rpm": cl["rpm"], "AT_pk": cl["AT_pk"], "AT_rms": cl["AT_rms"], "J_rms": cl["J_rms"],
            "feasible": cl["feasible"], "limited_by_maps": cl["limited_by_maps"], "binding": binding, "gamma_deg": cl["gamma_deg"], "torque_Nm": cl["torque"],
            "p_em_w": cl["p_em_w"], "p_shaft_net_w": cl["p_em_w"] - cl["windage"] - cl["bearing"] - (cl["fan"] if cl["fan_on_shaft"] else 0.0), "copper_w": cl["copper"], "copper_slot_w": cl["copper_slot"], "copper_end_w": cl["copper_end"], "iron_teeth_w": cl["iron_teeth"],
            "iron_yoke_w": cl["iron_yoke"], "iron_rotor_w": cl["iron_rotor"], "windage_w": cl["windage"], "bearing_w": cl["bearing"], "fan_w": cl["fan"],
            "fan_on_shaft": cl["fan_on_shaft"], "efficiency": cl["efficiency"], "hot_spot_c": cl["hot_spot"], "magnet_c": cl["magnet"], "v_air_ms": cl["v_air"],
            "J_winding_limit": cl["J_winding_limit"], "J_magnet_limit": cl["J_magnet_limit"]}


# ---------------------------------------------------------------------------------------------------- 3. where the heat goes
def heat_rows(mach, assume, stack_mm: float = MAIN_STACK) -> list[dict]:
    cs = coolings(stack_mm, assume)
    rows = []
    for key, rpm in HEAT_CASES:
        m = tm.model_for(cs[key], stack_mm, rpm, mach, assume)
        cl = tm.continuous_limit(m, mach)
        hf = cl["heat"]
        tot = hf["total_made"]
        r = {"cooling": key, "rpm": rpm, "stack_mm": stack_mm, "AT_pk": cl["AT_pk"], "feasible": cl["feasible"], "binding": cl["binding"], "total_made_w": tot,
             "balance_w": hf["balance_w"], "magnet_c": cl["magnet"], "hot_spot_c": cl["hot_spot"], "fan_w": cl["fan"]}
        for k, v in hf["made"].items():
            r[f"made_{k}_w"] = v
        for k, v in hf["exits"].items():
            r[f"out_{k}_w"] = v
        for k, v in hf["groups"].items():
            r[f"grp_{k}_w"] = v
        for k, v in hf["links"].items():
            r[f"link_{k}_w"] = v
        for k in [c for c in list(r) if c.startswith(("made_", "out_", "grp_", "link_"))]:
            r[k[:-2] + "_pct"] = 100.0 * r[k] / tot if tot > 0 else float("nan")
        rows.append(r)
    return rows


# ---------------------------------------------------------------------------------------------------- 4. the capability table
def capability_rows(mach, assume, stacks=STACKS) -> list[dict]:
    """Per cooling, stack and magnet grade at 3,000 rpm: the continuous point and the 2x short-term point."""
    rows = []
    for stack in stacks:
        cs = coolings(stack, assume)
        for key in KEYS:
            for mag, mlim in MAGNETS:
                m = tm.model_for(cs[key], stack, RATED_RPM, mach, assume)
                cl = tm.continuous_limit(m, mach, magnet_limit_c=mlim)
                sh = tm.short_term(m, mach, cl, 2.0)
                r = flat(key, mag, cl)
                r.update({"two_x_AT_pk": sh["AT_pk"], "two_x_torque_Nm": sh["torque"], "two_x_beyond_maps": sh["beyond_maps"], "two_x_t_winding_s": sh["t_winding_s"],
                          "two_x_t_magnet_s": sh["t_magnet_s"], "two_x_t_first_s": sh["t_first_s"], "two_x_first": sh["first"], "two_x_iron_w": sh.get("iron_total_w")})
                rows.append(r)
    return rows


# ---------------------------------------------------------------------------------------------------- 5. the sensitivity of the shaft fan to its coefficients
SENS = [("fan_efficiency", 0.10, 0.35), ("fan_k_system", 2.0, 6.0), ("duct_f_re", 54.0, 64.0), ("fan_fin_gain_conv", 1.5, 3.5), ("fan_fin_gain_rad", 1.0, 1.3),
        ("mixed_exponent", 2.0, 4.0), ("air_bulk_warming", False, True), ("fan_air_preheat", False, True), ("emissivity", 0.3, 0.95), ("contact_factor", 1.0, 8.0), ("housing_wall_mm", 3.0, 10.0),
        ("t_air_film_c", 40.0, 80.0)]
SENS_SPEEDS = (0, 1000, 3000, 10000)


def sensitivity_rows(mach, stack_mm: float = MAIN_STACK) -> list[dict]:
    """Each assumed coefficient of the shaft fan (and the housing's) at the ends of its range, one at a time: the continuous AT (N45, real losses) of both sizings at three
    speeds, and the fan's power. Sizing (b) is made again for each setting (its budget is the design constraint), so a coefficient that changes the fan's power also
    changes the fan's speed."""
    rows = []
    base = tm.assumptions()
    settings = [("(base)", None, base)]
    for name, lo, hi in SENS:
        for val in (lo, hi):
            cur = getattr(base, name)
            settings.append((name, val, tm.Assumptions(**{**base.__dict__, name: type(cur)(val)})))
    for name, val, asm in settings:
        cs = coolings(stack_mm, asm)
        row = {"coefficient": name, "setting": val}
        for key in ("shaft_a", "shaft_b"):
            for rpm in SENS_SPEEDS:
                m = tm.model_for(cs[key], stack_mm, rpm, mach, asm)
                cl = tm.continuous_limit(m, mach)
                row[f"{key}_AT_{rpm}"] = cl["AT_pk"]
            row[f"{key}_v_ref_ms"] = cs[key].v_ref_ms
            row[f"{key}_P_fan_3000_w"] = thermal.shaft_fan_power(RATED_RPM, cs[key], stack_mm, assume=asm)
            row[f"{key}_P_fan_top_w"] = thermal.shaft_fan_power(MAX_SPEED_RPM, cs[key], stack_mm, assume=asm)
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------------------------------- checks
def checks(mach, assume, heat: list[dict], stack_mm: float = MAIN_STACK) -> list[tuple[str, bool, str]]:
    """Consistency checks of the pieces against each other (printed and kept in the output)."""
    out = []
    # 1. energy balance of every heat-path case
    worst = max(abs(r["balance_w"]) for r in heat)
    out.append(("every heat-path case: heat out = heat made (to 1e-6 W: the state is solved to 1e-10 K)", worst < 1e-6, f"largest imbalance {worst:.2e} W"))
    # 2. the network's DC copper loss against the copper-loss model's DC loss at one point
    m = tm.model_for(Cooling.forced_air(), stack_mm, RATED_RPM, mach, assume)
    at = 118.6
    pt = mach.point(at, RATED_RPM, stack_mm)
    f, share, r100 = mach.copper(at, pt["gamma_deg"], RATED_RPM, stack_mm)
    j = at / (math.sqrt(2.0) * m.ampere_turns(1.0)["a_cu_side_mm2"])
    p_net = m.copper(j).total(100.0)
    p_mod = r100["p_dc_slot"] + r100["p_dc_end"]
    out.append(("DC copper loss: network (end turn 48.55 mm, the winding's fill) = copper-loss model to 1e-5 (the slot area is carried to 5 digits in the two)", abs(p_net / p_mod - 1.0) < 1e-5,
                f"{p_net:.6f} W against {p_mod:.6f} W at 118.6 A-turns, 100 C, 60 mm; the AC parts add {100 * (r100['p_total'] / p_mod - 1):.3f} %"))
    # 3. the fan laws: P(n) / n^3 from 5,000 to 16,330 rpm
    ca = Cooling.shaft_fan(6.0, float(RATED_RPM))
    p3 = thermal.shaft_fan_power(RATED_RPM, ca, stack_mm, assume=assume)
    p16 = thermal.shaft_fan_power(MAX_SPEED_RPM, ca, stack_mm, assume=assume)
    ratio = p16 / (p3 * (MAX_SPEED_RPM / RATED_RPM) ** 3)
    out.append(("fan laws: P(16,330) / (P(3,000) (16,330/3,000)^3) close to 1 (the friction factor falls slowly with Re)", 0.9 < ratio < 1.05, f"{ratio:.3f}"))
    # 4. standstill: the shaft fan has no forced convection
    ms = tm.model_for(ca, stack_mm, 0, mach, assume)
    out.append(("standstill: the shaft fan's forced conductance is zero (natural convection and radiation only)", ms.g_forced_h == 0.0 and ms.g_forced_b == 0.0,
                f"v {ms.v_air} m/s"))
    # 5. shaft fan (a) at 3,000 rpm against the fixed forced air preset (the same fins and speed)
    a_lim = tm.continuous_limit(tm.model_for(ca, stack_mm, RATED_RPM, mach, assume), mach)
    f_lim = tm.continuous_limit(tm.model_for(Cooling.forced_air(), stack_mm, RATED_RPM, mach, assume), mach)
    d = 100.0 * (a_lim["AT_pk"] / f_lim["AT_pk"] - 1.0)
    out.append(("shaft fan (a) at 3,000 rpm against the fixed 6 m/s forced air: within 2 % (radiation and mixing add, the air's warming takes away)", abs(d) < 2.0,
                f"{a_lim['AT_pk']:.2f} against {f_lim['AT_pk']:.2f} A-turns ({d:+.2f} %)"))
    return out


# ---------------------------------------------------------------------------------------------------- output
def write_csv(name: str, rows: list[dict]) -> None:
    OUT().mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(OUT() / f"{name}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows([{k: (json.dumps(v) if isinstance(v, (list, dict)) else v) for k, v in r.items()} for r in rows])


def f1(x, n=1):
    return "-" if x is None or (isinstance(x, float) and math.isnan(x)) else f"{x:.{n}f}"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-figures", action="store_true")
    ap.add_argument("--quick", action="store_true", help="the 60 mm stack only")
    ap.add_argument("--no-sensitivity", action="store_true")
    ap.add_argument("--figures-only", action="store_true", help="redraw the figures from the saved fan_study.json (no calculation)")
    args = ap.parse_args(argv)
    if args.figures_only:
        from optimization import thermal_fan_figures
        d = json.loads((OUT() / "fan_study.json").read_text(encoding="utf-8"))
        for p in thermal_fan_figures.make_all(d["fan"], d["sizing"], d["speed"], d["heat"], d["capability"]):
            print(p)
        return 0
    t0 = time.perf_counter()
    lines: list[str] = []

    def say(text: str = "") -> None:
        print(text)
        lines.append(text)
    mach = tm.Machine()
    assume = tm.assumptions()
    stacks = (MAIN_STACK,) if args.quick else STACKS
    say(f"Machine: the flux maps cover up to {mach.at_max:.1f} ampere-turns a coil side (peak); the example winding {mach.winding} has fill {mach.fill:.4f}; "
        f"network end turn {assume.end_turn_mm} mm, overhang {assume.end_overhang_mm} mm")

    # ---- 1. the fan
    ftab = fan_table(mach, assume)
    write_csv("fan_model", ftab)
    say("\n1. THE SHAFT FAN, 60 mm stack: air speed in the fin channels, the channel's Reynolds number, the pressure the fan builds, the flow, the air power and the shaft power, the warming of the "
        "air by the fan's own loss, and the outside's conductance to the ambient (housing + bells at 80 C surfaces; natural convection on the plain housing for comparison)")
    for r in ftab:
        say(f"  ({r['sizing'][-1]}) v_ref {r['v_ref_ms']:5.2f} m/s  {r['rpm']:6d} rpm: v {r['v_ms']:6.2f} m/s  Re {r['re']:7.0f}  f {r['f_darcy']:.4f}  dp {r['dp_pa']:7.1f} Pa  Q {r['q_m3s']:.4f} m3/s  "
            f"P_air {r['p_air_w']:8.2f} W  P_shaft {r['p_shaft_w']:9.2f} W  dT_air {r['dT_air_k']:5.2f} K  G_out {r['g_outside_w_per_k']:6.3f} W/K (plain natural {r['g_outside_natural_plain_w_per_k']:.3f})")
    sz = sizing_table(mach, assume)
    write_csv("fan_sizing", sz)
    pa3 = thermal.shaft_fan_power(RATED_RPM, Cooling.shaft_fan(6.0, float(RATED_RPM)), MAIN_STACK, assume=assume)
    pa16 = thermal.shaft_fan_power(MAX_SPEED_RPM, Cooling.shaft_fan(6.0, float(RATED_RPM)), MAIN_STACK, assume=assume)
    say(f"\n   THE SIZING TRADE (60 mm, N45, class F, MTPA, real losses): the continuous ampere-turns a coil side (peak) at five speeds, and the fan's power there, for a fan "
        f"whose power at 16,330 rpm is held to a budget. Sizing (a), 6 m/s at 3,000 rpm, takes {pa3:.1f} W at 3,000 rpm and {pa16:.0f} W at 16,330 rpm.")
    say("   budget / sizing        v_ref   v_top |" + "".join(f" {s:>6d} rpm: AT (P_fan W) |" for s in TRADE_SPEEDS))
    for r in sz:
        say(f"   {r['label']:22s} {r['v_ref_ms']:6.2f} {r['v_top_ms']:7.2f} |" + "".join(
            f" {r[f'AT_{s}']:6.1f}{'*' if not r[f'feasible_{s}'] else ' '} ({r[f'P_fan_{s}_w']:7.1f}) |" for s in TRADE_SPEEDS))
    say("   (* a limit is passed with no current at all: the iron, windage and bearing heat alone)")
    say("   sizing (b) for each stack (the finned length sets the channel's friction a little): " + ", ".join(
        f"{st:.0f} mm {coolings(st, assume)['shaft_b'].v_ref_ms:.3f} m/s at 3,000 rpm" for st in STACKS))
    pre = preheat_rows(mach, assume)
    write_csv("fan_preheat", pre)
    say("   THE FAN'S OWN HEAT warms its air (all of its shaft power ends as heat in the stream before the fins): continuous ampere-turns with / without that warming, 60 mm, N45")
    for r in pre:
        say(f"      ({r['sizing']}) at {r['rpm']:6d} rpm: fan {r['p_fan_w']:7.1f} W, the air warmed by {r['dT_air_k']:5.2f} K: {r['AT_with']:6.1f} A-turns with the warming, {r['AT_without']:6.1f} without "
            + (f"({100 * (r['AT_with'] / r['AT_without'] - 1):+.1f} %)" if r["AT_without"] > 0 else "(no current either way: the iron, windage and bearing heat alone passes the limit)"))
    p_em_a = tm.continuous_limit(tm.model_for(Cooling.shaft_fan(6.0, float(RATED_RPM)), MAIN_STACK, RATED_RPM, mach, assume), mach)["p_em_w"]
    mtop = tm.model_for(Cooling.natural(), MAIN_STACK, MAX_SPEED_RPM, mach, assume)
    say(f"   THE BUDGET CHOSEN: {BUDGET_W:.0f} W at 16,330 rpm. Why: it is {100 * BUDGET_W / p_em_a:.1f} % of the continuous mechanical power (an estimate, T x omega = {p_em_a:.0f} W) that "
        f"the 6 m/s shaft fan allows at 3,000 rpm, about what a small motor's fan costs as a share of its output (check before citing); it is {BUDGET_W / (mtop.p_windage + mtop.p_bearing):.1f} "
        f"times the motor's other friction losses at the top speed (windage {mtop.p_windage:.1f} W + bearings {mtop.p_bearing:.1f} W); sizing (a) would take {pa16:.0f} W there. "
        f"A larger budget buys more between 6,000 and 10,000 rpm and nothing at 3,000 rpm (the table above): the budget is a limit on the peak power, and the energy it costs below "
        f"the top speed is small (P goes as n^3).")

    # ---- 2. the continuous limit against speed
    sp = speed_rows(mach, assume, stacks)
    write_csv("fan_speed_limit", sp)
    say("\n2. THE CONTINUOUS LIMIT AGAINST SPEED (60 mm, N45 80 C, class F, MTPA, real iron loss, the copper-loss model, windage, bearings, the fan): ampere-turns a coil side, peak "
        "(* = a limit is passed with no current at all; ^ = the flux maps' edge, 355.9, is reached first; w / m = the winding / the magnets bind); torque N.m at that point")
    say("   rpm     " + "".join(f"{NAMES[k][:26]:>28s}" for k in KEYS))
    for rpm in SPEEDS:
        cells = []
        for k in KEYS:
            r = next(x for x in sp if x["cooling"] == k and x["magnets"] == "N45" and x["stack_mm"] == MAIN_STACK and x["rpm"] == rpm)
            flag = "*" if not r["feasible"] else ("^" if r["limited_by_maps"] else " ")
            cells.append(f"{r['AT_pk']:6.1f}{flag} T {f1(r['torque_Nm'], 2):>6s} {BIND_CHAR.get(r['binding'], '-')}")
        say(f"   {rpm:6d}  " + "".join(f"{c:>28s}" for c in cells))
    say("   the same with the 150 C magnet grade (class F, 60 mm):")
    for rpm in SPEEDS:
        cells = []
        for k in KEYS:
            r = next(x for x in sp if x["cooling"] == k and x["magnets"] == "150C" and x["stack_mm"] == MAIN_STACK and x["rpm"] == rpm)
            flag = "*" if not r["feasible"] else ("^" if r["limited_by_maps"] else " ")
            cells.append(f"{r['AT_pk']:6.1f}{flag} T {f1(r['torque_Nm'], 2):>6s} {BIND_CHAR.get(r['binding'], '-')}")
        say(f"   {rpm:6d}  " + "".join(f"{c:>28s}" for c in cells))

    # ---- 3. where the heat goes
    ht = heat_rows(mach, assume)
    write_csv("fan_heat_paths", ht)
    say("\n3. WHERE THE HEAT GOES at each cooling's continuous limit (60 mm, N45, class F): watts and per cent of the heat made in the motor (the fan's power is a shaft load that "
        "heats the air and is not in the sources)")
    for r in ht:
        say(f"  {NAMES[r['cooling']][:34]:34s} {r['rpm']:6d} rpm: {r['AT_pk']:6.1f} A-turns, heat made {r['total_made_w']:7.1f} W (fan {r['fan_w']:7.1f} W), hot spot {r['hot_spot_c']:5.1f} C, "
            f"magnets {r['magnet_c']:5.1f} C" + ("" if r["feasible"] else "   [a limit is passed with no current]"))
        say("      made:  " + ", ".join(f"{k[5:-2]} {r[k]:.1f}" for k in r if k.startswith("made_") and k.endswith("_w") and abs(r[k]) > 0.005))
        say("      out:   " + ", ".join(f"{k[4:-2]} {r[k]:.1f} W ({r[k[:-2] + '_pct']:.1f} %)" for k in r if k.startswith("out_") and k.endswith("_w") and abs(r[k]) > 0.005))
        say("      inside:" + ", ".join(f" {k[4:-2]} {r[k]:.1f} W ({r[k[:-2] + '_pct']:.1f} %)" for k in r if k.startswith("grp_") and k.endswith("_w")))

    # ---- 4. the capability table
    cap = capability_rows(mach, assume, stacks)
    write_csv("fan_capability", cap)
    for mag, mlim in MAGNETS:
        say(f"\n4. CAPABILITY at 3,000 rpm, 40 C ambient, class F, magnets {mag} ({mlim:.0f} C): continuous point (MTPA, real losses; T x omega is an ESTIMATE of the mechanical power), "
            f"then 2x the current (torque from the maps, time from cold to the first limit)")
        say("   cooling                          stack |    AT pk (rms)  J rms | T N.m   P_em W | Pcu   Fe(t/y/r)        Wind Brg  Fan | eff % | hot spot  magnets  limit | 2x: AT    T N.m   first limit (s)")
        for r in cap:
            if r["magnets"] != mag:
                continue
            first = "-" if r["two_x_t_first_s"] is None else f"{r['two_x_first']} {r['two_x_t_first_s']:.0f}"
            if not r["feasible"]:
                first = "no current"
            tq2 = "beyond maps" if r["two_x_torque_Nm"] is None and r["feasible"] else f1(r["two_x_torque_Nm"], 2)
            say(f"   {NAMES[r['cooling']][:30]:30s} {r['stack_mm']:4.0f} | {r['AT_pk']:6.1f} ({r['AT_rms']:5.1f}) {r['J_rms']:5.2f} | {f1(r['torque_Nm'], 2):>6s} {r['p_em_w']:7.0f} | "
                f"{r['copper_w']:5.1f} {r['iron_teeth_w']:4.1f}/{r['iron_yoke_w']:4.1f}/{r['iron_rotor_w']:3.1f} {r['windage_w']:4.2f} {r['bearing_w']:4.2f} {r['fan_w']:5.1f} | "
                f"{100 * r['efficiency']:5.1f} | {r['hot_spot_c']:6.1f}  {r['magnet_c']:6.1f}  {r['binding'][:7]:7s}| {r['two_x_AT_pk']:6.1f} {tq2:>11s}   {first}")

    # ---- 5. sensitivity
    if not args.no_sensitivity:
        sens = sensitivity_rows(mach)
        write_csv("fan_sensitivity", sens)
        say("\n5. SENSITIVITY of the shaft fan to its coefficients (60 mm, N45, class F, MTPA): continuous ampere-turns of sizings (a) and (b) at 0 / 1,000 / 3,000 / 10,000 rpm and the fans' power; "
            "sizing (b) is made again for each setting (its 50 W budget is the constraint)")
        base = sens[0]
        for r in sens:
            tag = f"{r['coefficient']} {r['setting']}" if r["setting"] is not None else "(base)"
            cells = " ".join(f"{r[f'{k}_AT_{n}']:6.1f}" for k in ("shaft_a", "shaft_b") for n in SENS_SPEEDS)
            say(f"   {tag:30s} (a) / (b) AT at 0 / 1,000 / 3,000 / 10,000 rpm: {cells}   v_ref (b) {r['shaft_b_v_ref_ms']:.2f} m/s, P_fan (a) {r['shaft_a_P_fan_3000_w']:.1f} W at 3,000 and {r['shaft_a_P_fan_top_w']:.0f} W at the top")
    else:
        sens = []

    # ---- checks
    say("\nCHECKS")
    chk = checks(mach, assume, ht)
    for name, ok, detail in chk:
        say(f"  {'PASS' if ok else 'FAIL'}  {name}  [{detail}]")
    (OUT() / "fan_study.json").write_text(json.dumps({"fan": ftab, "sizing": sz, "speed": sp, "heat": ht, "capability": cap, "sensitivity": sens,
                                                      "checks": [{"name": n, "ok": ok, "detail": d} for n, ok, d in chk]}, indent=1, default=float), encoding="utf-8")
    if not args.no_figures:
        from optimization import thermal_fan_figures
        figs = thermal_fan_figures.make_all(ftab, sz, sp, ht, cap)
        say("\nfigures: " + ", ".join(str(p) for p in figs))
    say(f"\nstudy done in {time.perf_counter() - t0:.1f} s")
    (OUT() / "fan_study_output.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0 if all(ok for _, ok, _ in chk) else 1


if __name__ == "__main__":
    raise SystemExit(main())
