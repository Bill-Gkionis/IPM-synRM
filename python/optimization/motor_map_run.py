"""optimization.motor_map_run: the motor map, end to end: the checks, the envelopes of the three shaft fans and of the comparison coolings, the efficiency maps, the stack table, the winding
study's front and the 150 C magnet grade, the FEMM spot-check and the figures. Every number in its report and figures is printed by this file's own run.

What it is, and why
    motor_map.py is the engine (the envelope, the losses, the two controls), motor_map_thermal.py gives it the current limits from the thermal network, motor_map_stack.py the masses, the
    inertia and the windings of the stacks. This file runs them on the winding study's windings and holds every number to a second route before it reports it. The inputs are inputs: a
    later run swaps the winding(s), the copper model, the cooling or the bus by changing a command-line option or a line below (WINDING, STACK_WINDINGS, EXTRAS, K_MARGIN, VDCS, STACKS,
    COPPER, ARRANGEMENT, the coolings() list), never the engine.

The inputs
    windings       the winding study's (optimization/winding_study.py, Storage 0acd770: 48 V, the 50 W shaft fan, N45 at 80 C, class F, compact bundles, at most 8 strands in hand; each is the
                   smallest inverter current that still holds the continuous torque to 3,300 rpm): 26 turns a coil side, 4 paths, 8 x 0.315 mm on 20 mm; 13 / 4, 4 x 0.63 on 40 mm; 9 / 4,
                   7 x 0.56 mm (0.63 over the enamel, fill 0.428, R(100 C) 19.26 mOhm) on 60 mm, the headline; 7 / 4, 7 x 0.63 on 80 mm. The stack table has two columns: "opt" (these) and
                   "same" (the headline winding on every stack). The extras: the study's other front members at 60 mm (8 / 4 4 x 0.8, 5 / 4 8 x 0.71, 4 / 4 8 x 0.8, 3 / 4 8 x 0.9 mm: the
                   inverter classes of 100, 150, 200 and 300 A) and its pick for a 150 C magnet grade (8 / 4, 4 x 0.8 mm), run with the headline fan
    bus            48 and 60 V; V_max = 0.95 V_dc / sqrt 3 (motor_map.Bus: dead time, drops and a regulator's reserve; the winding study used 1.00)
    stacks         20, 40, 60, 80 mm (1 ... 4 rows of magnets); the headline is 60 mm
    cooling        a design decision: a fan on the shaft (a simple TEFC motor; liquid cooling is parked): the air speed in the fin channels proportional to the speed, only natural
                   convection at standstill, its own shaft power going as the cube of the speed and counted as a mechanical load. Three sizings (thermal.Cooling.shaft_fan,
                   thermal.size_shaft_fan with the diamond end turn, on each stack): (a) 6 m/s at 3,000 rpm, (b) held to 50 W at 16,330 rpm (1.76 m/s at 3,000 rpm on 60 mm: the headline fan),
                   (c) held to 200 W at 16,330 rpm (2.80 m/s). Comparisons: natural convection, the fixed 6 m/s external fan (not on the shaft) and the liquid jacket. N45 magnets at 80 C,
                   winding under 155 C
    currents       continuous: the thermal network's limit at every speed, with the loss of the maps at the operating point it allows; short-term: twice the winding's own continuous current
                   at 3,000 rpm with no voltage limit (the thermal rating), constant over the speed, from a COLD start (its duration at 3,000 rpm is in the stack table); both clipped at the
                   edge of the flux maps. (The winding study sizes the inverter at twice a common rating, 5 % under the example winding's: its currents are a few percent under these.)
    losses         copper: copper_loss.py's full loss (DC with the diamond end turn, skin, proximity, PWM ripple, circulating currents of strands in a compact bundle); iron: the maps' at the
                   CAD's 0.3 mm lamination; mechanical: bearings, windage and the fan's shaft power, a rough estimate. The AC part of the copper loss is charged to the input power
                   (motor_map.Machine.point), so the energy balance (c1) closes with it

The checks, and the limits written before any of them was run (constants TOL_ at the head of this file)
    (a) below the corner the envelope's torque is maps.mtpa's at the current limit: 8 machines x 4 current limits x 2 speeds; torque within 1e-5, angle within 0.02 degrees
    (b) at the corner |v| = V_max, by hand from the MTPA point's flux linkages and R, and the corner speed is the first study's quadratic: within 1e-6
    (c) the energy balance closes at every point of every envelope and map: |P_in - (P_out + P_cu + P_fe + P_mech)| / P_in within 1e-9, with P_in from the voltage equations (plus the AC
        copper loss); and the d-q torque agrees with FEMM's air-gap torque at every point within 0.5 % + 0.01 N.m per 20 mm (maps_report's own limit)
    (d) the first study's example winding against its 3,411 rpm: the maps at its rated point (flux linkages within 0.25 % of |psi|), the corner at its own point (0.3 %), the engine's corner
        at its own MTPA angle (1.5 %)
    (e) the MTPV branch against the condition from the incremental inductances: the gradient of the torque parallel to the gradient of |v|^2, sin of the angle within 1e-3, the multiplier positive
    (f) FEMM at deep field-weakening points reproduces the maps (psi within 0.25 % of |psi|, torque within 0.5 % + 0.01 N.m per 20 mm), through femm_queue.femm_pool: stage "spot"
    (g) the thermal coupling: the network's copper loss at the converged temperature equals the engine's (1e-6); the continuous current at 3,000 rpm against thermal.py's sustainable()
        (reported, with the iron loss of each)
    (g2) the continuous limit against thermal_machine.continuous_limit (the thermal study's own coupling of the same loss models): with their inputs (the headline winding, full copper loss,
        strands in a cluster, diamond end turn, no voltage limit) within 1 % at 0, 1,000, 3,000, 6,000, 10,000 and 16,330 rpm for the three shaft fans; (g2b) thermal_machine on this
        module's fans reproduces the thermal study's table (its example winding) to half its last decimal; then this module's limit with each bus's voltage limit is printed beside them, to
        say where and why the numbers differ
    (h) the mechanical estimate equals thermal.py's windage and bearing formulas (1e-9); (j) the summed iron-loss splines equal iron_loss.IronLossMap (1e-9)
    (k) the minimum-loss solver is not beaten by a dense scan of 9,001 rays (1e-6); (m) the envelope is not beaten by a brute-force grid of the whole plane (1e-9)
    (n) the level-set solver hits its torque within 1e-8 N.m and minimum-loss control is never worse than the MTPA control; (o) every solution is inside both limits (1e-8); (p) the
        continuous envelopes stay inside the temperature limits (each case's own magnet limit)
    (q) the winding study's numbers through this module's calls in the study's own configuration (k = 1, the air-gap torque, flat currents, copper 100 C): (q1) the copper fill and R(100 C)
        of nine windings (1e-3), (q2) their own continuous ratings (2e-3), (q3) its torque-speed tables of the five 60 mm windings, to the two decimals and the integer rpm it printed
        (0.0051 N.m, 0.51 rpm; the currents are the study's own ratings)
    A check that fails is reported as it failed; no limit is changed after a result.

The stages and the commands (from the Python to FEMM folder; times on this PC)
    uv run python -m optimization.motor_map_run --stages checks                 the checks that need no run (25 s)
    uv run python -m optimization.motor_map_run --stages run,figures            the cases, the maps, the checks over every solution, the stack table, the figures (about 21 minutes)
    uv run python -m optimization.motor_map_run --stages report --label final   print the headline, the stack tables, the highest speeds and the extras of a saved run again (seconds)
    uv run python -m optimization.motor_map_run --stages spot --label final     the FEMM spot-check (takes the solver slot for 20 s; no other stage needs FEMM)
    options: --label NAME (the output folder; default "final"), --quick (coarse grids: a development run), --winding N,a,n_s,d,d_over (the base winding), --stack-windings 20:...;40:...,
             --copper dc|full, --arrangement cluster|row|transposed|random|line, --magnet-limit C, --no-extras, --vdc 48,60, --k 0.95, --stacks 20,40,60,80, --headline 60
Writes out\\optimization\\motor_map\\<label>\\: report.txt (the checks and the run), report_summary.txt (the headline, the stack tables, the highest speeds, the net output, the extras and
the maps of a saved run: --stages report), report_spot.txt (the FEMM spot-check), results.pkl (every array), envelopes.csv (every envelope, a row a speed), stack_table.csv, inputs.json,
spot_check.json, figures\\ (the PNGs with their CSVs).

Numbers to expect (this PC, 2026-11-09; the whole run takes 21 minutes; the headline stack, 60 mm, the winding study's 9 / 4, 7 x 0.56 mm, k = 0.95, N45 at 80 C; short-term is from a cold
    start; torque and power are at the shaft, after the iron drag, the windage, the bearings and the fan; currents are peak, current densities rms)
                                            (a) 6 m/s at 3,000 rpm        (b) 50 W at 16,330 rpm        (c) 200 W at 16,330 rpm
    thermal rating, A-turns (A peak)        154.7 (68.8)                  94.6 (42.0)                   95.4 (42.4)
    peak current, twice the rating [A]      137.5                         84.1                          84.8
    48 V short-term corner                  2,350 rpm, 13.64 N.m          2,718 rpm, 8.67 N.m           2,711 rpm, 8.75 N.m
    48 V short-term at 16,330 rpm           0.50 N.m, 851 W               1.59 N.m, 2,716 W             1.51 N.m, 2,579 W
    48 V continuous: T at 3,000 rpm, ends   7.00 N.m, 13,091 rpm (loads)  4.20 N.m, 7,400 rpm (t)       4.23 N.m, 11,200 rpm (t)
    60 V short-term corner                  2,983 rpm, 13.61 N.m          3,433 rpm, 8.66 N.m           3,425 rpm, 8.73 N.m
    60 V short-term at 16,330 rpm           0.99 N.m, 1,695 W             2.06 N.m, 3,531 W             1.99 N.m, 3,398 W
    60 V continuous: ends                   13,830 rpm (loads)            9,800 rpm (t)                 12,300 rpm (t)
    efficiency map, best (48 / 60 V)        93.8 / 93.8 %                 95.2 / 95.4 %                 95.1 / 95.2 %
    (t: the continuous envelope ends because no current both holds the voltage and keeps the temperatures; the loads would let it go on.) The headline winding reaches 16,330 rpm only for a
    moment with any of the fans, and fan (a) itself takes 1,941 W there. The checks print: (a) 1.7e-8, (b) 2e-9, (c1) 9e-16, (c2) 0.03 of its limit, (d1) 0.017 %, (e) 1.4e-8, (f) psi
    0.016 % of |psi| and 0.3 mN.m (FEMM 5.0, the maps' fine mesh of 15,693 nodes, 15 rotor positions), (g) 1.9e-7, (g2) 8.5e-3 (the 16,330 rpm point of fan (a), 33 A-turns), (q1) 1e-4,
    (q2) 2.9e-4, (q3) 0.0050 N.m and 0.48 rpm. The thermal study's table is 2 % above thermal_machine with the headline winding (its example winding has more copper in the slot, fill 0.448
    against 0.428); this module's limit WITH the voltage limit is above that by 6 ... 37 % where field weakening takes the flux, and with it the iron loss, off the magnets (75 against 55
    A-turns at 6,000 rpm for fan (b) on 48 V), and zero where the field-weakening current that holds the voltage already overheats (fan (b) above 7,400 rpm).

Self-check (the answers follow).
  1. The headline winding carries 94.6 A-turns as 42.0 A and the 20 mm stack's winding 108.5 A-turns as 16.7 A. Why is the headline's current about 2.5 times larger, and what does the inverter see?
  2. At 60 V the headline winding with fan (b) stops carrying torque continuously at 9,800 rpm, yet only a few hundredths of a N.m of fan, windage and bearings load the shaft there. What ends it?
Answers.
  1. The current is the ampere-turns times the paths over the turns: 94.6 x 4 / 9 = 42.0 A and 108.5 x 4 / 26 = 16.7 A. A stack three times as long links three times the flux a turn, so the
     same corner speed on the same bus asks for about a third of the turns, and the current density the cooling allows then asks for a larger current: the inverter's current rating, not its
     voltage, goes with the stack.
  2. The voltage, not the load. Above the corner speed the inverter must hold a field-weakening current at every speed, because the no-load voltage (112 V a phase at the top speed) is far above
     the bus (V_max 33 V a phase); that current grows toward the characteristic current while the cooling's limit (the magnets at 80 C) falls with the iron loss and with the laminar flow in the
     fin channels, and at 9,800 rpm the two cross: no current both holds the voltage and keeps the temperatures (the t in the tables).
"""




from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import sys
import time
from pathlib import Path

import numpy as np

import paths
from optimization import design, thermal
from optimization import motor_map as mm
from optimization import motor_map_stack as st
from optimization import motor_map_thermal as mt
from optimization import winding as winding_study

OUT_ROOT = paths.path("out") / "optimization" / "motor_map"

# ====================================================================================================== the inputs (inputs, not rulings: each is swapped by changing it here or by an option)
WINDING = st.OPTIMIZED_WINDINGS[60.0]     # the base winding ("same" on every stack): the winding study's headline winding, 9 turns a coil side, 4 paths, 7 x 0.56 mm (0.63 over the enamel)
STACK_WINDINGS = dict(st.OPTIMIZED_WINDINGS)   # "opt": the winding study's winding of each stack (optimization/winding_study.py, Storage 0acd770)
EXAMPLE_WINDING = mm.PROVISIONAL_WINDING  # the first study's example (26 turns, 4 paths, 8 x 0.315 mm): the reference of checks (d), (e), (g2) and of the thermal study's table
VARIANTS = ("opt", "same")                # the windings of a stack: the winding study's own ("opt"), then the base winding on every stack ("same")
HEADLINE_VARIANT = "opt"
STACKS = st.STACKS_MM                      # 1 ... 4 rows of magnets
HEADLINE_STACK = 60.0                      # 3 rows
VDCS = (48.0, 60.0)
K_MARGIN = 0.95                            # the share of the modulation's linear limit used (motor_map.Bus)
T_CU0 = 100.0                              # the copper temperature the continuous solve starts from; the network then gives it
T_CU_SHORT = 100.0                         # the copper temperature of the short-term envelope (the end of a pulse: conservative for the voltage)
PEAK_FACTOR = 2.0                          # the short-term (peak) current: this many times the continuous current at the rated speed with no voltage limit
RATED_RPM = 3000.0
STEP_HEADLINE, STEP_TABLE = 100.0, 250.0   # the speed steps [rpm] of the headline cases and of the other stacks
GRID_HEADLINE = {"n_gamma": 361, "n_radius": 241}
GRID_TABLE = {"n_gamma": 181, "n_radius": 121}
MAP_SPEEDS, MAP_FRACTIONS = 64, 41         # the efficiency map: speeds from 150 rpm to the end of the envelope, torques from 0 to the envelope in this many steps
COPPER = "full"                            # "full": copper_loss.py's full loss (skin, proximity, PWM, circulating currents, the diamond end turn); "dc": the first study's DC loss with its end turn
ARRANGEMENT = "cluster"                    # how the strands of a turn lie in the slot (copper_loss.py): a compact bundle, the hand-wound default
MAGNET_LIMIT_C = thermal.MAGNET_LIMIT_C    # the magnets' temperature limit of the headline grade (N45, 80 C)
EXTRAS_ON = True                           # also run the front's other windings and the 150 C magnet grade (EXTRAS below)

# the winding study's other front members at 60 mm (inverter classes of 100, 150, 200 and 300 A; winding_study.py, report "THE SPEED RANGE ON 48 V") and its pick for the 150 C magnet grade,
# run on the headline stack with the 50 W shaft fan (b): key, label, winding, magnet limit [C]
EXTRAS = (
    {"key": "f8", "label": "front member 8 / 4, 4 x 0.8 mm", "winding": mm.Winding(8, 4, 4, 0.8, 0.884), "magnet_limit_c": 80.0},
    {"key": "f5", "label": "front member 5 / 4, 8 x 0.71 mm", "winding": mm.Winding(5, 4, 8, 0.71, 0.789), "magnet_limit_c": 80.0},
    {"key": "f4", "label": "front member 4 / 4, 8 x 0.8 mm", "winding": mm.Winding(4, 4, 8, 0.8, 0.884), "magnet_limit_c": 80.0},
    {"key": "f3", "label": "front member 3 / 4, 8 x 0.9 mm", "winding": mm.Winding(3, 4, 8, 0.9, 0.989), "magnet_limit_c": 80.0},
    {"key": "g150", "label": "150 C magnet grade: 8 / 4, 4 x 0.8 mm", "winding": mm.Winding(8, 4, 4, 0.8, 0.884), "magnet_limit_c": 150.0},
)
EXTRA_FAN = "fanB"                         # the shaft fan of the extras: sizing (b), the headline fan
STUDY_I_ST = {"opt": 81.5, "f8": 91.7, "f5": 146.8, "f4": 183.5, "f3": 244.6, "g150": 155.4}   # the winding study's short-term currents [A peak] (2 x its common rating), as it reported them

# ====================================================================================================== the limits of the checks, written before any run
TOL_A_T, TOL_A_GAMMA = 1e-5, 0.02                      # (a) envelope torque against maps.mtpa: relative; the angle in degrees
TOL_B = 1e-6                                           # (b) |v| at the corner against V_max, and the corner against the quadratic's root: relative
TOL_C_ALG = 1e-9                                       # (c) the energy balance, relative to the input power
TOL_C_AGE_REL, TOL_C_AGE_ABS = 5e-3, 0.01              # (c) |T_dq - T_age| of an operating point: 0.5 % of T plus 0.01 N.m per 20 mm of stack (maps_report's own limit)
TOL_D_PSI = 2.5e-3                                     # (d) the maps at the first study's rated point: psi against its FEMM, relative to |psi| (the limit of the maps' interpolation)
TOL_D_CORNER_SAME, TOL_D_CORNER_ENGINE = 3e-3, 1.5e-2  # (d) the corner at the study's own point; the engine's corner at its own MTPA angle
TOL_E = 1e-3                                           # (e) sin of the angle between grad T and grad |v|^2 on the MTPV branch
TOL_G = 1e-6                                           # (g) the network's copper loss against this module's, relative
TOL_G2 = 1e-2                                          # (g2) the continuous ampere-turns against thermal_machine.continuous_limit, same inputs: relative (the bisections stop at 1e-3 and 7e-4)
TOL_G2B = 0.051                                        # (g2b) thermal_machine on these fans against the thermal study's table: the table has one decimal, so half a unit of it [A-turns]
TOL_Q_T, TOL_Q_RPM = 0.0051, 0.51                      # (q) the envelopes against the winding study's table: its two decimals of N.m, its integer rpm (half a unit and a margin)
TOL_Q_AT, TOL_Q_FILL, TOL_Q_R = 2e-3, 1e-3, 1e-3       # (q) the own ratings (the two bisections), the fill and the resistance against the study's: relative
TOL_H = 1e-9                                           # (h) the mechanical loss against thermal.py's formulas, relative
TOL_J = 1e-9                                           # (j) the summed iron-loss splines against iron_loss.IronLossMap, relative
TOL_K = 1e-6                                           # (k) the minimum-loss solver's loss against a dense scan's minimum: the solver must not be worse by more than this (relative)
TOL_M = 1e-9                                           # (m) the envelope torque against a brute-force grid's maximum: the envelope must not be lower by more than this (relative)
TOL_N = 1e-8                                           # (n) the level-set solver's torque error [N.m]
TOL_O = 1e-8                                           # (o) the constraints at every solution: |v| <= V_max, |i| <= limit, relative
TOL_F_PSI, TOL_F_T_REL, TOL_F_T_ABS = 2.5e-3, 5e-3, 0.01   # (f) the FEMM spot-check: psi of |psi|; T relative plus 0.01 N.m per 20 mm (maps_report's own interpolation limits)


class Report:
    """The lines printed, and the checks (name, passed, what was said)."""

    def __init__(self):
        self.lines: list[str] = []
        self.checks: list[tuple[str, bool, str]] = []

    def say(self, text: str = "") -> None:
        print(text, flush=True)
        self.lines.append(text)

    def check(self, name: str, ok: bool, said: str) -> bool:
        ok = bool(ok)
        self.checks.append((name, ok, said))
        self.say(f"  [{'PASS' if ok else 'FAIL'}] {name}: {said}")
        return ok


# ====================================================================================================== the coolings
FAN_KEYS = ("fanA", "fanB", "fanC")        # the three sizings of the shaft fan
MAIN_FAN = "fanB"                          # the headline fan (the candidate: held to 50 W at the top speed); the single-picture outputs use it
FAN_BUDGETS_W = (50.0, 200.0)              # (b) and (c): the fan's shaft power at the top speed, on each stack
_FANS: dict = {}


def winding_of(stack: float, variant: str) -> mm.Winding:
    """The winding of a stack in a variant: "same" the base winding, "opt" the winding study's winding of that stack (STACK_WINDINGS)."""
    return st.stack_winding(WINDING, stack, variant, STACK_WINDINGS)


def thermal_assumptions() -> thermal.Assumptions:
    """The network's assumptions: with copper_loss.py's diamond end turn the network is told the same end turn and overhang (one value in both models), and the fans are sized with them."""
    return thermal.Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6) if COPPER == "full" else thermal.Assumptions()


def fan_cooling(key: str, stack: float = HEADLINE_STACK) -> thermal.Cooling:
    """The shaft fan of a sizing on a stack (thermal.Cooling.shaft_fan: the air speed in the fin channels proportional to the speed, at standstill only natural convection and radiation;
    its own shaft power P ~ n^3 is a mechanical load on the shaft): (a) 6 m/s at 3,000 rpm on every stack; (b) and (c) held to 50 W and 200 W at the top speed ON THE STACK
    (thermal.size_shaft_fan with the run's assumptions, as the winding study sized its fan for each stack)."""
    stack = float(stack)
    if (key, stack, COPPER) not in _FANS:
        if key == "fanA":
            _FANS[(key, stack, COPPER)] = thermal.Cooling.shaft_fan(6.0, 3000.0)
        else:
            budget = FAN_BUDGETS_W[("fanB", "fanC").index(key)]
            _FANS[(key, stack, COPPER)] = thermal.size_shaft_fan(budget, mm.TOP_RPM, 3000.0, stack, assume=thermal_assumptions())["cooling"]
    return _FANS[(key, stack, COPPER)]


def fan_sizings(stack: float = HEADLINE_STACK) -> dict:
    return {k: fan_cooling(k, stack) for k in FAN_KEYS}


def coolings() -> list[dict]:
    """The coolings, the shaft fans first (a design decision: a fan on the shaft, like a simple TEFC motor; liquid cooling is parked), then the comparison lines: natural convection, the
    fixed 6 m/s external fan (its air speed does not follow the shaft; its power is not a load on the shaft and not in the motor's efficiency) and the liquid jacket. Each entry:
    key, label, short, cooling (on the headline stack), cooling_of(stack) -> the thermal Cooling on a stack, and fan_factory(stack) -> the fan's shaft power as a function of the speed
    (None: no fan on the shaft)."""
    out = []
    for key, name in zip(FAN_KEYS, ("(a) 6 m/s at 3,000 rpm", f"(b) held to {FAN_BUDGETS_W[0]:g} W at 16,330 rpm", f"(c) held to {FAN_BUDGETS_W[1]:g} W at 16,330 rpm")):
        c = fan_cooling(key)

        def factory(stack, key=key):
            c_s = fan_cooling(key, stack)
            return lambda r: thermal.shaft_fan_power(r, c_s, stack, assume=thermal_assumptions())
        out.append({"key": key, "label": f"shaft fan {name} ({c.v_ref_ms:.2f} m/s in the fin channels at 3,000 rpm)", "short": f"shaft fan {name}", "cooling": c,
                    "cooling_of": (lambda stack, key=key: fan_cooling(key, stack)), "fan_factory": factory})
    for key, label, short, c in (("natural", "natural convection", "natural convection", thermal.Cooling.natural()),
                                 ("external", "fixed 6 m/s external fan (not on the shaft)", "external fan, 6 m/s", thermal.Cooling.forced_air()),
                                 ("liquid", "liquid jacket (parked)", "liquid (parked)", thermal.Cooling.liquid())):
        out.append({"key": key, "label": label, "short": short, "cooling": c, "cooling_of": (lambda stack, c=c: c), "fan_factory": None})
    return out


def fan_of(cool: dict, stack: float):
    return cool["fan_factory"](stack) if cool["fan_factory"] else None


def bus_of(vdc: float) -> mm.Bus:
    return mm.Bus(vdc, K_MARGIN)


def copper_model(kind: str | None = None):
    """The copper-loss model of the run (the swap point): copper_loss.py's full loss with the strands in the run's arrangement, or the first study's DC loss."""
    return mm.FullCopper(arrangement=ARRANGEMENT) if (kind or COPPER) == "full" else mm.DcCopper()


# ====================================================================================================== summaries
def ends(sm: dict) -> str:
    """Where the torque ends, as text: the speed, with a trailing t where the temperatures (no current both holds the voltage and keeps the limits), not the loads, end it."""
    return fmt(sm["zero_rpm"]) + ("t" if sm.get("stop") == "limit" else "")


def summarise(env: mm.Envelope) -> dict:
    """The numbers a reader asks of an envelope: the peak torque and power, the corner, the torque and power at the top speed and the speed where the torque ends."""
    rpm, T, P, eta = env["rpm"], env["T"], env["P_out"], env["eta"]
    on = env["region"] != "off"
    Pm = np.where(on, P, 0.0)
    jP = int(np.argmax(Pm)) if on.any() else 0
    top = len(rpm) - 1
    last = int(np.nonzero(on)[0].max()) if on.any() else -1
    if on[top]:
        stop = "top"                                    # the torque lasts to the top speed
    elif last + 1 < len(rpm) and bool(np.asarray(env["feasible"], bool)[last + 1]):
        stop = "loads"                                  # a best point still exists above the last speed with torque: the loads take all of its torque
    else:
        stop = "limit"                                  # no operating point inside the current and voltage limits (the continuous limit is 0 there: no current both holds the voltage and keeps the temperatures)
    out = {"T_peak": float(T.max()), "T_peak_rpm": float(rpm[int(np.argmax(T))]), "P_peak": float(Pm[jP]), "P_peak_rpm": float(rpm[jP]),
           "corner_rpm": env.meta.get("corner_rpm"), "corner_T": env.meta.get("corner_T"), "zero_rpm": env.meta.get("zero_rpm"), "stop": stop,
           "T_rated": float(np.interp(RATED_RPM, rpm, T)), "T_top": float(T[top]), "P_top": float(Pm[top]) if on[top] else 0.0,
           "eta_max": float(np.nanmax(eta)) if on.any() else float("nan"), "eta_min": float(np.nanmin(np.where(rpm >= 500.0, eta, np.nan))) if on.any() else float("nan"),
           "reaches_top": bool(on[top])}
    for k in ("corner_rpm", "corner_T", "zero_rpm"):
        out[k] = None if out[k] is None else float(out[k])
    return out


def plain(env: mm.Envelope) -> dict:
    """An envelope as plain arrays and numbers (for the pickle and the CSV)."""
    meta = {k: (float(v) if isinstance(v, (np.floating, float)) else v) for k, v in env.meta.items()}
    return {"d": {k: np.asarray(v) for k, v in env.d.items()}, "meta": meta, "summary": summarise(env)}


# ====================================================================================================== the cases
class CaseRunner:
    """Runs the continuous and the short-term envelope of one machine in one cooling at both buses, and keeps the thermal rating, the overload times and the geometry. winding and
    magnet_limit_c override the variant's winding and the magnets' temperature limit (the extras: the front's other windings, the 150 C grade)."""

    def __init__(self, cool: dict, stack: float, variant: str, fine: bool, winding: mm.Winding | None = None, magnet_limit_c: float | None = None):
        self.cool, self.stack, self.variant, self.fine = cool, stack, variant, fine
        self.w = winding if winding is not None else winding_of(stack, variant)
        self.magnet_limit_c = MAGNET_LIMIT_C if magnet_limit_c is None else float(magnet_limit_c)
        self.m = mm.Machine(self.w, stack, copper=copper_model(), mech=mm.Mechanical(fan=fan_of(cool, stack)))
        self.cfg = mt.ThermalConfig(cool["cooling_of"](stack), thermal_assumptions(), magnet_limit_c=self.magnet_limit_c)
        self.sc = mm.Scan(self.m, i_max=self.m.i_cap, **(GRID_HEADLINE if fine else GRID_TABLE))
        self.speeds = mm.speed_grid(step=STEP_HEADLINE if fine else STEP_TABLE)
        self.rating = mt.thermal_rating(self.m, self.cfg, RATED_RPM, scan=self.sc)
        self.i_peak = min(PEAK_FACTOR * self.rating["I"], self.m.i_cap)
        self.peak_clipped = PEAK_FACTOR * self.rating["I"] > self.m.i_cap

    def limit_point(self, env: mm.Envelope, lims, tcs, bus: mm.Bus):
        """The highest speed the motor reaches with this fan on its shaft, and what the torque is made of there: at the speed where the shaft torque ends the air-gap torque of the best
        point equals the loads' torques (the fan, the windage, the bearings and the iron-loss drag). None where the torque lasts to the top speed, and where the envelope ends for want of any
        operating point (summarise's "stop" is "limit": the continuous limit, not the loads, ends it)."""
        z = env.meta.get("zero_rpm")
        if z is None or z >= mm.TOP_RPM * (1.0 - 1e-9) or summarise(env)["stop"] != "loads":
            return None
        lim, tc = float(np.interp(z, self.speeds, lims)), float(np.interp(z, self.speeds, tcs))
        sol = self.sc.max_torque(self.sc.ctx(z, bus, tc, "shaft"), lim)
        if sol is None:
            return None
        pt = self.m.point(sol["i_d"], sol["i_q"], z, bus, tc)
        w = 2.0 * math.pi * z / 60.0
        return {"rpm": float(z), "I": float(sol["I"]), "i_limit": lim, "T_em": float(pt["T_em"]), "T_fan": float(pt["P_fan"]) / w, "T_windage": float(pt["P_windage"]) / w,
                "T_bearing": float(pt["P_bearing"]) / w, "T_iron": float(pt["P_fe"]) / w, "T_shaft": float(pt["T_shaft"]), "region": "FW" if sol["v_active"] and sol["i_active"] else ("MTPV" if sol["v_active"] else "MTPA")}

    def run(self, vdc: float) -> dict:
        bus = bus_of(vdc)
        cl = mt.continuous_limits(self.m, bus, self.cfg, self.speeds, scan=self.sc, t_cu0=T_CU0)
        lim = np.where(cl["limited"] == "none", 0.0, cl["I"])
        env_c = mm.envelope(self.m, bus, self.speeds, lim, t_c=cl["t_cu"], scan=self.sc)
        env_s = mm.envelope(self.m, bus, self.speeds, self.i_peak, t_c=T_CU_SHORT, scan=self.sc)
        over = mt.overload_from_cold(self.m, bus, self.cfg, RATED_RPM, self.i_peak, self.sc, T_CU_SHORT)
        lim_cols = {k: cl[k] for k in ("I", "t_cu", "hot_spot", "magnet", "margin", "converged")}
        lim_cols["limited"], lim_cols["binding"] = np.array(cl["limited"], dtype=str), np.array(cl["binding"], dtype=str)
        return {"vdc": vdc, "cont": plain(env_c), "short": plain(env_s), "limits": lim_cols, "overload": {k: v for k, v in over.items() if k != "iron"},
                "limit_cont": self.limit_point(env_c, lim, cl["t_cu"], bus), "limit_short": self.limit_point(env_s, np.full(len(self.speeds), self.i_peak), np.full(len(self.speeds), T_CU_SHORT), bus)}

    def static(self) -> dict:
        g = st.geometry_numbers(self.m)
        psi0 = self.m.back_emf_flux()
        top_e = 2.0 * math.pi * mm.PP * mm.TOP_RPM / 60.0 * psi0
        return {"cool": self.cool["key"], "stack": self.stack, "variant": self.variant, "magnet_limit_c": self.magnet_limit_c,
                "v_ref_ms": float(self.cool["cooling_of"](self.stack).v_ref_ms) if self.cool["fan_factory"] else None,
                "winding": {"turns": self.w.turns, "paths": self.w.paths, "strands": self.w.strands,
                "wire_mm": self.w.wire_mm, "wire_over_mm": self.w.wire_over_mm, "series_turns": self.w.series_turns, "fill": self.m.fill},
                "R_100C": float(self.m.resistance(100.0)), "i_cap": self.m.i_cap, "rating": {k: v for k, v in self.rating.items() if k not in ("sources",)},
                "i_peak": self.i_peak, "peak_clipped": self.peak_clipped, "geometry": g, "psi_pm": psi0, "emf_top_peak": top_e,
                "emf_top_line_line_peak": mm.SQRT3 * top_e, "rpm_emf_vmax": {f"{v:g}": bus_of(v).v_max / (mm.PP * 2.0 * math.pi / 60.0 * psi0) for v in VDCS}}


# ====================================================================================================== the checks
def check_h_j(rep: Report) -> None:
    """(h) the mechanical estimate against thermal.py's formulas, (j) the summed iron-loss splines against iron_loss.IronLossMap."""
    m = mm.Machine(WINDING, 60.0)
    worst = 0.0
    for rpm in (500.0, 3000.0, 8000.0, 16330.0):
        w = 2.0 * math.pi * rpm / 60.0
        mine = m.mech.power(rpm, m.r_ro, m.gap, m.L)
        theirs = thermal.windage_loss(w, m.r_ro, m.gap, m.L, 70.0, 1.0), thermal.bearing_loss(w, thermal.Assumptions())
        worst = max(worst, abs(float(mine["windage"]) / theirs[0] - 1.0), abs(float(mine["bearing"]) / theirs[1] - 1.0))
    rep.check("(h) mechanical loss = thermal.py's windage and bearing formulas", worst <= TOL_H, f"largest relative difference {worst:.1e} over 4 speeds (limit {TOL_H:g}); at 16,330 rpm "
              f"windage {float(m.mech.power(16330.0, m.r_ro, m.gap, m.L)['windage']):.2f} W, bearings {float(m.mech.power(16330.0, m.r_ro, m.gap, m.L)['bearing']):.2f} W on 60 mm")
    rng = np.random.default_rng(20261110)
    I0 = rng.uniform(0.0, 3.0, 200) * m.maps.rated
    g = np.radians(rng.uniform(0.0, 90.0, 200))
    i_d, i_q = -I0 * np.sin(g) / m.s, I0 * np.cos(g) / m.s
    rpm = rng.uniform(100.0, 16330.0, 200)
    c = m.iron_coeffs(i_d, i_q)
    f = mm.C_SPEED * 2.0 * np.pi * rpm / 60.0
    tot = c[0] * f + c[1] * f ** 2 + c[2] * f ** 1.5
    ref = m.iron(i_d, i_q, rpm)["total"]
    err = float(np.max(np.abs(tot - ref) / np.maximum(ref, 1e-12)))
    rep.check("(j) summed iron-loss splines = iron_loss.IronLossMap", err <= TOL_J, f"largest relative difference {err:.1e} over 200 random points and speeds (limit {TOL_J:g})")


def check_a(rep: Report) -> dict:
    """(a) below the corner the envelope's torque is maps.mtpa's at the current limit: 8 machines (4 stacks x the 2 windings of the variants), 4 current limits, 2 speeds, the torque of the d-q model."""
    worst_t, worst_g, n = 0.0, 0.0, 0
    bus = bus_of(48.0)
    for stack in STACKS:
        for variant in VARIANTS:
            w = winding_of(stack, variant)
            m = mm.Machine(w, stack)
            for mult in (0.5, 1.0, 2.0, 3.0):
                lim = mult * m.maps.rated / m.s
                sc = mm.Scan(m, i_max=lim, **GRID_TABLE)
                g, _, _, T = m.maps.mtpa(lim, turns=w.turns, paths=w.paths, stack_mm=stack)
                for rpm in (0.0, 50.0):
                    s = sc.max_torque(sc.ctx(rpm, bus, 100.0, "air_gap"), lim)
                    worst_t, worst_g, n = max(worst_t, abs(s["T"] / T - 1.0)), max(worst_g, abs(s["gam"] - g)), n + 1
    rep.check("(a) below the corner the envelope's torque = maps.mtpa at the current limit", worst_t <= TOL_A_T and worst_g <= TOL_A_GAMMA,
              f"{n} cases: largest torque difference {worst_t:.1e} (limit {TOL_A_T:g}), largest angle difference {worst_g:.4f} deg (limit {TOL_A_GAMMA:g})")
    return {"worst_t": worst_t, "worst_gamma": worst_g, "cases": n}


def by_hand_voltage(m: mm.Machine, i_d: float, i_q: float, rpm: float, t_c: float) -> tuple[float, float, float]:
    """|v| written out by hand from the maps' public flux call and R = rho(T) l_mt n / (a n_s S): the second route of check (b). Returns (|v|, psi_d, psi_q)."""
    psi_d, psi_q = m.maps.psi_dq(i_d, i_q, turns=m.w.turns, paths=m.w.paths, stack_mm=m.stack_mm)
    l_mt = 2.0 * m.L + 2.0 * winding_study.end_turn_geometry(m.p)["l_end_mm"] * 1e-3
    s_strand = math.pi * (m.w.wire_mm * 1e-3) ** 2 / 4.0
    R = winding_study.RHO_20 * (1.0 + winding_study.ALPHA_20 * (t_c - 20.0)) * l_mt * (12.0 * m.w.turns / m.w.paths) / (m.w.paths * m.w.strands * s_strand)
    w_e = 4.0 * 2.0 * math.pi * rpm / 60.0
    return math.hypot(R * i_d - w_e * psi_q, R * i_q + w_e * psi_d), psi_d, psi_q


def corner_quadratic(R: float, i_d: float, i_q: float, psi_d: float, psi_q: float, v_max: float) -> float:
    """The first study's corner speed [rpm]: the root of (psi_d^2 + psi_q^2) w^2 + 2 R (i_q psi_d - i_d psi_q) w + R^2 (i_d^2 + i_q^2) - V^2 = 0."""
    A = psi_d ** 2 + psi_q ** 2
    B = 2.0 * R * (i_q * psi_d - i_d * psi_q)
    C = R ** 2 * (i_d ** 2 + i_q ** 2) - v_max ** 2
    w = (-B + math.sqrt(B * B - 4.0 * A * C)) / (2.0 * A)
    return w * 60.0 / (2.0 * math.pi * 4.0)


def check_b(rep: Report) -> dict:
    """(b) at the corner |v| = V_max, by hand from the MTPA point's flux linkages and R, and the corner against the first study's quadratic, for 8 machines at 1.5 times rated."""
    worst_v, worst_c, n = 0.0, 0.0, 0
    bus = bus_of(48.0)
    for stack in STACKS:
        for variant in VARIANTS:
            w = winding_of(stack, variant)
            m = mm.Machine(w, stack)
            lim = 1.5 * m.maps.rated / m.s
            env = mm.envelope(m, bus, mm.speed_grid(step=500.0), lim, basis="air_gap", scan=mm.Scan(m, i_max=lim, **GRID_TABLE))
            rc = env.meta["corner_rpm"]
            sc = mm.Scan(m, i_max=lim, **GRID_TABLE)
            s = sc.max_torque(sc.ctx(rc, bus, 100.0, "air_gap"), lim)
            v, psi_d, psi_q = by_hand_voltage(m, s["i_d"], s["i_q"], rc, 100.0)
            q = corner_quadratic(float(m.resistance(100.0)), s["i_d"], s["i_q"], psi_d, psi_q, bus.v_max)
            worst_v, worst_c, n = max(worst_v, abs(v / bus.v_max - 1.0)), max(worst_c, abs(q / rc - 1.0)), n + 1
    rep.check("(b) at the corner |v| = V_max, by hand; the corner = the first study's quadratic", worst_v <= TOL_B and worst_c <= TOL_B,
              f"{n} machines: |v| / V_max - 1 at most {worst_v:.1e}, corner / quadratic - 1 at most {worst_c:.1e} (limit {TOL_B:g})")
    return {"worst_v": worst_v, "worst_corner": worst_c}


def check_d(rep: Report) -> dict:
    """(d) the example winding against the first study's corner speed (3,411 rpm: the rated torque held to 3,000 rpm with 10 % of the 48 V bus to spare, V_max = V_dc / sqrt 3,
    no margin, copper at 100 C, the rated point 100 ampere-turns at the study's own angle, the field of a separate FEMM run)."""
    z = json.load(open(paths.path("out") / "optimization" / "winding" / "optimized" / "femm.json", encoding="utf-8"))["rated"]
    ex = json.load(open(paths.path("out") / "optimization" / "winding" / "optimized" / "results.json", encoding="utf-8"))["example"]["row"]
    m = mm.Machine(EXAMPLE_WINDING, 20.0)
    s = m.s
    i0_d, i0_q = float(z["i_d"]), float(z["i_q"])
    pd, pq = m.maps.psi_dq(i0_d, i0_q)
    err_d, err_q = abs(pd - z["psi_d"]) / math.hypot(z["psi_d"], z["psi_q"]), abs(pq - z["psi_q"]) / math.hypot(z["psi_d"], z["psi_q"])
    i_d, i_q = i0_d / s, i0_q / s
    R = float(m.resistance(100.0))
    bus1 = mm.Bus(48.0, 1.0)
    pdw, pqw = m.flux(i_d, i_q)
    rpm_same = corner_quadratic(R, i_d, i_q, float(pdw), float(pqw), bus1.v_max)
    rpm_study_psi = corner_quadratic(R, i_d, i_q, s * z["psi_d"], s * z["psi_q"], bus1.v_max)
    i_lim = math.hypot(i_d, i_q)
    sc = mm.Scan(m, i_max=i_lim, **GRID_HEADLINE)
    env = mm.envelope(m, bus1, mm.speed_grid(step=250.0), i_lim, basis="air_gap", scan=sc)
    rpm_engine = env.meta["corner_rpm"]
    target = float(ex["corner_rpm"])
    d1 = max(err_d, err_q)
    d2 = abs(rpm_same / target - 1.0)
    d3 = abs(rpm_engine / target - 1.0)
    rep.say(f"  (d) the first study's case: N = 26, a = 4, 8 x 0.315 mm, 20 mm, V_max = 48 / sqrt 3 = {bus1.v_max:.3f} V, R = {R:.5f} Ohm (the study: {ex['R_hot_ohm']:.5f}), i = {i_lim:.3f} A peak "
            f"(the study: {ex['i_peak_A']:.3f})")
    rep.check("(d1) the maps at the study's rated point reproduce its FEMM flux linkages", d1 <= TOL_D_PSI,
              f"psi_d {1e3 * pd:.4f} mWb against {1e3 * z['psi_d']:.4f}, psi_q {1e3 * pq:.4f} against {1e3 * z['psi_q']:.4f}: {100 * d1:.3f} % of |psi| (limit {100 * TOL_D_PSI:g} %)")
    rep.check("(d2) the corner speed at the study's own point: maps against the study", d2 <= TOL_D_CORNER_SAME, f"{rpm_same:.1f} rpm from the maps' flux linkages against {target:.1f} "
              f"(and {rpm_study_psi:.1f} from the study's own flux linkages: the same quadratic): {100 * d2:.3f} % (limit {100 * TOL_D_CORNER_SAME:g} %)")
    rep.check("(d3) the engine's corner speed (its own MTPA angle) against the study's", d3 <= TOL_D_CORNER_ENGINE, f"{rpm_engine:.1f} rpm at the maps' MTPA angle (the study held its own angle, "
              f"{z['gamma']:.2f} deg): {100 * d3:.3f} % (limit {100 * TOL_D_CORNER_ENGINE:g} %)")
    return {"rpm_same": rpm_same, "rpm_study": target, "rpm_engine": rpm_engine, "R": R, "psi_err": d1}


def check_e(rep: Report) -> dict:
    """(e) the MTPV branch against the condition from the incremental inductances: where the current limit does not bind and the voltage does, the gradient of the torque must be
    parallel to the gradient of |v|^2 (Lagrange), with the cross terms of the maps: dT/di_d = 1.5 p (L_dd i_q - L_qd i_d - psi_q), dT/di_q = 1.5 p (psi_d + L_dq i_q - L_qq i_d),
    d(|v|^2 / 2)/di_d = v_d (R - w L_qd) + v_q w L_dd, d(|v|^2 / 2)/di_q = -v_d w L_qq + v_q (R + w L_dq). The torque is the d-q model's (basis "air_gap")."""
    w = EXAMPLE_WINDING
    m = mm.Machine(w, 20.0)
    lim = 2.0 * w.current_from_density(6.0)
    bus = bus_of(48.0)
    sc = mm.Scan(m, i_max=lim, **GRID_HEADLINE)
    rpms = mm.speed_grid(step=100.0)
    env = mm.envelope(m, bus, rpms, lim, basis="air_gap", scan=sc, locate=False)
    sel = np.nonzero(env["region"] == "MTPV")[0]
    if len(sel) == 0:
        rep.check("(e) the MTPV branch against the condition from the incremental inductances", False, "no MTPV points in the test case: the check was not exercised")
        return {}
    i_d, i_q, rpm = env["i_d"][sel], env["i_q"][sel], env["rpm"][sel]
    L = m.maps.inductances(i_d, i_q, turns=w.turns, paths=w.paths, stack_mm=m.stack_mm)
    R = float(m.resistance(100.0))
    w_e = 4.0 * 2.0 * math.pi * rpm / 60.0
    psi_d, psi_q = m.flux(i_d, i_q)
    v_d, v_q = R * i_d - w_e * psi_q, R * i_q + w_e * psi_d
    gT = np.stack([1.5 * 4 * (L["L_dd"] * i_q - L["L_qd"] * i_d - psi_q), 1.5 * 4 * (psi_d + L["L_dq"] * i_q - L["L_qq"] * i_d)])
    gV = np.stack([v_d * (R - w_e * L["L_qd"]) + v_q * w_e * L["L_dd"], -v_d * w_e * L["L_qq"] + v_q * (R + w_e * L["L_dq"])])
    cross = np.abs(gT[0] * gV[1] - gT[1] * gV[0]) / (np.hypot(*gT) * np.hypot(*gV))
    lam = (gT[0] * gV[0] + gT[1] * gV[1]) / (gV[0] ** 2 + gV[1] ** 2)
    ok = cross.max() <= TOL_E and (lam > 0).all()
    rep.check("(e) the MTPV branch satisfies grad T parallel to grad |v|^2 (incremental inductances)", ok,
              f"{len(sel)} MTPV points of the 20 mm / 48 V / {lim:.1f} A case ({rpm.min():.0f} ... {rpm.max():.0f} rpm): sin(angle) at most {cross.max():.1e} (limit {TOL_E:g}), "
              f"the multiplier positive at {int((lam > 0).sum())} of {len(sel)}")
    return {"points": int(len(sel)), "max_sin": float(cross.max()), "rpm_range": (float(rpm.min()), float(rpm.max()))}


def check_m_k(rep: Report) -> dict:
    """(m) the envelope against a brute-force grid of the whole plane (1,201 radii x 1,801 angles), (k) the minimum-loss solver against a dense scan of 9,001 rays with the radius solved
    by bisection: the solver must not be beaten."""
    bus = bus_of(48.0)
    worst_m, worst_k, nk = -1.0, -1.0, 0
    for label, w, stack in (("example winding, 60 mm", EXAMPLE_WINDING, 60.0), ("headline winding, 60 mm", winding_of(60.0, "opt"), 60.0)):
        m = mm.Machine(w, stack)
        lim = m.w.current_from_density(6.0)
        sc = mm.Scan(m, i_max=lim, **GRID_HEADLINE)
        for rpm in (1500.0, 3000.0, 6000.0):
            c = sc.ctx(rpm, bus, 100.0, "shaft")
            s = sc.max_torque(c, lim)
            Ig, Gg = np.meshgrid(np.linspace(0.0, lim, 1201), np.linspace(0.0, 90.0, 1801))
            T, vv = sc.state(Ig, Gg, c)
            Tb = float(np.where(vv <= c.v2, T, -np.inf).max())
            if s is not None and Tb > 0:
                worst_m = max(worst_m, (Tb - s["T"]) / abs(Tb))
        # (k): a dense scan of rays at a few torques, losses by the same total_loss
        for rpm, frac in ((1500.0, 0.3), (3000.0, 0.2), (3000.0, 0.6), (6000.0, 0.1)):
            c = sc.ctx(rpm, bus, 100.0, "shaft")
            s = sc.max_torque(c, lim)
            if s is None or s["T"] <= 0:
                continue
            target = frac * s["T"]
            res = sc.level_set(c, lim, np.array([target]), objectives=("loss",))["loss"]
            if not res["found"][0]:
                worst_k, nk = max(worst_k, 1.0), nk + 1
                continue
            g = np.radians(res["gam"][0])
            solver = float(m.total_loss(-res["I"][0] * np.sin(g), res["I"][0] * np.cos(g), rpm, bus.v_dc, 100.0))
            gam = np.linspace(0.0, 90.0, 9001)
            lo, hi = np.zeros_like(gam), np.full_like(gam, lim)
            for _ in range(48):
                mid = 0.5 * (lo + hi)
                up = sc.state(mid, gam, c)[0] >= target
                hi, lo = np.where(up, mid, hi), np.where(up, lo, mid)
            Iq = hi
            Tq, vq = sc.state(Iq, gam, c)
            good = (np.abs(Tq - target) < 1e-7) & (vq <= c.v2)
            if good.any():
                gg = np.radians(gam[good])
                dense = float(m.total_loss(-Iq[good] * np.sin(gg), Iq[good] * np.cos(gg), rpm, bus.v_dc, 100.0).min())
                worst_k, nk = max(worst_k, (solver - dense) / dense), nk + 1
    rep.check("(m) the envelope is not beaten by a brute-force grid", worst_m <= TOL_M, f"the grid's best torque less the envelope's, relative, at most {worst_m:.1e} (negative: the envelope is higher; limit {TOL_M:g}), "
              f"6 speed-machine pairs")
    rep.check("(k) the minimum-loss solver is not beaten by a dense scan", worst_k <= TOL_K, f"{nk} points: the solver's loss less the dense scan's minimum, relative, at most {worst_k:.1e} (limit {TOL_K:g})")
    return {"worst_m": worst_m, "worst_k": worst_k}


def check_g(rep: Report) -> dict:
    """(g) the thermal coupling: the network's copper loss at the converged temperature = this module's, and the continuous current at 3,000 rpm with no voltage limit against the
    network's own sustainable(), with the iron loss of each said. In every cooling, on the 20 mm and the 60 mm stack, each with the winding the winding study chose for it."""
    worst = 0.0
    rows = []
    for cool in coolings():
        for stack in (20.0, 60.0):
            m = mm.Machine(winding_of(stack, HEADLINE_VARIANT), stack, copper=copper_model(), mech=mm.Mechanical(fan=fan_of(cool, stack)))
            cfg = mt.ThermalConfig(cool["cooling_of"](stack), thermal_assumptions())
            rate = mt.thermal_rating(m, cfg, RATED_RPM)
            sc = mm.Scan(m, i_max=m.i_cap, **GRID_TABLE)
            s = sc.max_torque(sc.ctx(RATED_RPM, mt.NO_VOLTAGE_LIMIT, rate["t_cu"], "shaft"), rate["I"])
            mg_, r, pt = mt.state(m, mt.NO_VOLTAGE_LIMIT, cfg, RATED_RPM, s["i_d"], s["i_q"], rate["t_cu"])
            t_mean = mt.copper_temperature(m, r)
            mine = float(m.copper.loss(m, s["i_d"], s["i_q"], RATED_RPM, t_mean, 48.0)["total"])
            worst = max(worst, abs(r.copper_loss / mine - 1.0))
            net = mt.network(m, cfg, RATED_RPM)
            sus = net.sustainable()
            ph = net.iron_placeholder()
            rows.append({"cooling": cool["short"], "stack": stack, "AT_mine": rate["AT_pk"], "AT_network": sus["AT_pk"], "J_mine": rate["J_rms"], "J_network": sus["J_rms"],
                         "binding": rate.get("binding"), "Pfe_mine": float(pt["P_fe"]), "Pfe_placeholder": ph["iron_teeth"] + ph["iron_yoke"] + ph["iron_rotor"],
                         "limited": rate["limited"]})
    rep.check("(g) the network's copper loss = this module's at the same temperature", worst <= TOL_G, f"largest relative difference {worst:.1e} over {len(rows)} cases (limit {TOL_G:g})")
    rep.say("  (g) the continuous current at 3,000 rpm, no voltage limit: this module (the iron loss of the maps at its operating point, the copper loss of the run) against thermal.py's sustainable() "
            "(its placeholder iron loss, its DC copper loss at the winding's fill), AT peak a coil side:")
    for r in rows:
        rep.say(f"      {r['cooling']:34s} {r['stack']:3.0f} mm: {r['AT_mine']:7.1f} against {r['AT_network']:7.1f} A-turns ({100 * (r['AT_mine'] / r['AT_network'] - 1):+5.1f} %), "
                f"iron loss {r['Pfe_mine']:5.1f} W against {r['Pfe_placeholder']:5.1f} W, binding: {r['binding']}{' (maps edge)' if r['limited'] == 'maps' else ''}")
    return {"rows": rows, "worst_copper": worst}


# the thermal study's table of continuous ampere-turns (peak, a coil side; 60 mm, N45, MTPA, no voltage limit; its example winding) at 0, 1,000, 3,000, 6,000, 10,000 and 16,330 rpm for the
# three shaft fans, as it was handed over (0: no continuous operation; None: not given). Check (g2b) holds this module's fans to it.
REFERENCE_TABLE = {"fanA": (95.1, 114.3, 158.1, 169.2, 142.4, 33.1), "fanB": (95.1, 110.1, 96.6, 56.0, 59.1, 0.0), "fanC": (None, 112.2, 97.4, 118.5, 99.6, 0.0)}
G2_SPEEDS = (0.0, 1000.0, 3000.0, 6000.0, 10000.0, 16330.0)


def winding_dict(w: mm.Winding) -> dict:
    """A winding in the form thermal_machine.Machine and copper_loss take."""
    return {"N": w.turns, "a": w.paths, "n_s": w.strands, "d_mm": w.wire_mm, "d_over_mm": w.wire_over_mm, "arrangement": ARRANGEMENT}


def check_g2(rep: Report, results: dict) -> dict:
    """(g2) the continuous limit against thermal_machine.continuous_limit (the thermal study's own coupling of the same three loss models): at MTPA currents, no voltage limit, 60 mm, for
    the three shaft fans at 0, 1,000, 3,000, 6,000, 10,000 and 16,330 rpm. With THEIR inputs (the headline winding with copper_loss.py's full loss, the strands in a cluster, its diamond end
    turn and the network's matching assumptions) this module must agree within TOL_G2; (g2b) the thermal study's own table (its example winding) must be reproduced by thermal_machine
    on this module's fans, so that the fans are the same; then this module's continuous limits WITH the voltage limit of each bus (field weakening above the corner) are printed beside them,
    to say where and why the numbers differ."""
    from optimization import thermal_machine as tm
    head_w = winding_of(HEADLINE_STACK, HEADLINE_VARIANT)
    mach_head, mach_ex = tm.Machine(winding_dict(head_w)), tm.Machine()
    worst, mismatch, worst_ref, rows = 0.0, 0, 0.0, []
    rep.say(f"  (g2) continuous ampere-turns a coil side (peak), 60 mm, N45 80 C: [table] the thermal study's table as handed over (its example winding); [tm ex] thermal_machine now, its example "
            f"winding; [tm] thermal_machine now with the headline winding ({head_w.label()}); [mine] this module with the same inputs and no voltage limit (against [tm]); then this module's "
            f"continuous limit with the voltage limit at " + " and ".join(f"{v:g} V" for v in VDCS) + " (0: no continuous operation)")
    rep.say(f"      {'':34s} {'rpm':>7s}      table  tm ex      tm   mine (diff)   " + "  ".join(f"{v:g} V" for v in VDCS))
    for cool in coolings()[:3]:
        c60 = cool["cooling_of"](HEADLINE_STACK)
        cfg = mt.ThermalConfig(c60, tm.assumptions())
        m_full = mm.Machine(head_w, HEADLINE_STACK, copper=mm.FullCopper(arrangement="cluster"), mech=mm.Mechanical(fan=fan_of(cool, HEADLINE_STACK)))
        sc_full = mm.Scan(m_full, i_max=m_full.i_cap, **GRID_TABLE)
        st_ = results["static"].get(f"{cool['key']}|{HEADLINE_STACK:g}|{HEADLINE_VARIANT}")
        for j, rpm in enumerate(G2_SPEEDS):
            ex = tm.continuous_limit(tm.model_for(c60, HEADLINE_STACK, rpm, mach_ex), mach_ex)
            at_ex = float(ex["AT_pk"]) if ex["feasible"] else 0.0
            th = tm.continuous_limit(tm.model_for(c60, HEADLINE_STACK, rpm, mach_head), mach_head)
            at_th = float(th["AT_pk"]) if th["feasible"] else 0.0
            at_full = mt.thermal_rating(m_full, cfg, rpm, scan=sc_full)["AT_pk"]
            given = REFERENCE_TABLE[cool["key"]][j]
            if given is not None:
                worst_ref = max(worst_ref, abs(at_ex - given))
            fw = []
            for vdc in VDCS:
                case = results["cases"].get(f"{cool['key']}|{HEADLINE_STACK:g}|{HEADLINE_VARIANT}|{vdc:g}")
                if case is None or st_ is None:
                    fw.append(float("nan"))
                    continue
                k = int(np.argmin(np.abs(case["cont"]["d"]["rpm"] - rpm)))
                fw.append(float(case["limits"]["I"][k]) * st_["winding"]["turns"] / st_["winding"]["paths"] if case["limits"]["limited"][k] != "none" else 0.0)
            both = at_th > 0 and at_full > 0
            if (at_th > 0) != (at_full > 0):
                mismatch += 1
            elif both:
                worst = max(worst, abs(at_full / at_th - 1.0))
            rows.append({"cooling": cool["short"], "rpm": rpm, "table": given, "tm_example": at_ex, "tm_headline": at_th, "mine_no_voltage_limit": at_full, "mine_with_voltage_limit": fw})
            gap = f"({100 * (at_full / at_th - 1.0):+5.2f} %)" if both else "(none)    "
            rep.say(f"      {cool['short']:34s} {rpm:7,.0f} rpm: {'     -' if given is None else f'{given:6.1f}'} {at_ex:6.1f} {at_th:7.1f} {at_full:6.1f} {gap}  " + "  ".join(f"{x:6.1f}" for x in fw))
    rep.check("(g2) the continuous limit agrees with thermal_machine.continuous_limit on the same inputs", worst <= TOL_G2 and mismatch == 0,
              f"{len(rows)} points (3 shaft fans x 6 speeds): largest relative difference {worst:.1e} where both are feasible (limit {TOL_G2:g}); {mismatch} points where one finds no continuous operation and the other does")
    rep.check("(g2b) thermal_machine on this module's fans reproduces the thermal study's table", worst_ref <= TOL_G2B,
              f"{sum(1 for r in rows if r['table'] is not None)} entries: largest difference {worst_ref:.2f} A-turns against the table's one decimal (limit {TOL_G2B:g}): the fans are the same")
    return {"worst": worst, "mismatch": mismatch, "worst_reference": worst_ref, "rows": rows}


# ====================================================================================================== (q) the winding study's windings against its own numbers
# The winding study (optimization/winding_study.py, Storage 0acd770; report.txt, results.json and shortlist_*.csv in out\\optimization\\winding_study\\) read this engine: for every winding it solved
# the continuous rating from the network and gave its torque-speed table "from motor_map.py's solver": k = 1 (V_max = V_dc / sqrt 3), the air-gap torque, the winding's own continuous
# current and twice it, flat over the speed, copper at 100 C, 48 V. Check (q) puts every number it printed through this module's own calls in that configuration, whatever the run's inputs.
# label, stack [mm], winding, magnet limit [C], own rating [A-turns peak a coil side], R(100 C) [mOhm], fill (None: not printed by the study)
STUDY_WINDINGS = (
    ("20 mm  26 / 4, 8 x 0.315", 20.0, mm.Winding(26, 4, 8, 0.315, 0.367), 80.0, 108.507, None, 0.4476),
    ("40 mm  13 / 4, 4 x 0.63", 40.0, mm.Winding(13, 4, 4, 0.63, 0.7042), 80.0, 101.623, None, 0.4476),
    ("60 mm   9 / 4, 7 x 0.56", 60.0, mm.Winding(9, 4, 7, 0.56, 0.63), 80.0, 94.574, 19.261, 0.42843),
    ("80 mm   7 / 4, 7 x 0.63", 80.0, mm.Winding(7, 4, 7, 0.63, 0.7042), 80.0, 90.098, None, 0.4217),
    ("60 mm   8 / 4, 4 x 0.8", 60.0, mm.Winding(8, 4, 4, 0.8, 0.884), 80.0, 96.0554, 14.682, 0.444113),
    ("60 mm   5 / 4, 8 x 0.71", 60.0, mm.Winding(5, 4, 8, 0.71, 0.789), 80.0, 95.2555, 5.825, 0.437261),
    ("60 mm   4 / 4, 8 x 0.8", 60.0, mm.Winding(4, 4, 8, 0.8, 0.884), 80.0, 95.7788, 3.670, 0.444113),
    ("60 mm   3 / 4, 8 x 0.9", 60.0, mm.Winding(3, 4, 8, 0.9, 0.989), 80.0, 93.3319, 2.175, 0.42156),
    ("60 mm   8 / 4, 4 x 0.8, 150 C magnets", 60.0, mm.Winding(8, 4, 4, 0.8, 0.884), 150.0, 162.779, 14.682, 0.444113),
)
# its torque-speed table [N.m] at 0, 1,000, 2,000, 3,000, 4,000, 5,000, 6,000, 8,000, 10,000, 12,000, 14,000 and 16,330 rpm: winding label -> (continuous, corner, last speed with torque,
# short-term, corner, last speed with torque)
STUDY_RPMS = (0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 8000.0, 10000.0, 12000.0, 14000.0, 16330.0)
STUDY_ENVELOPES = {
    "60 mm   9 / 4, 7 x 0.56": ((4.28, 4.28, 4.28, 4.28, 4.00, 3.17, 2.36, 0.64, 0.00, 0.00, 0.00, 0.00), 3442, 8291, (8.78, 8.78, 8.78, 8.71, 7.26, 5.99, 5.07, 3.86, 3.11, 2.60, 2.23, 1.92), 2866, 16330),
    "60 mm   8 / 4, 4 x 0.8": ((4.35, 4.35, 4.35, 4.35, 4.33, 3.72, 2.97, 1.59, 0.00, 0.00, 0.00, 0.00), 3875, 9481, (8.91, 8.91, 8.91, 8.91, 8.08, 6.76, 5.75, 4.40, 3.55, 2.97, 2.55, 2.19), 3232, 16330),
    "60 mm   5 / 4, 8 x 0.71": ((4.31, 4.31, 4.31, 4.31, 4.31, 4.31, 4.31, 3.72, 2.80, 1.93, 1.00, 0.00), 6264, 15040, (8.84, 8.84, 8.84, 8.84, 8.84, 8.84, 8.46, 6.83, 5.60, 4.73, 4.09, 3.52), 5264, 16330),
    "60 mm   4 / 4, 8 x 0.8": ((4.34, 4.34, 4.34, 4.34, 4.34, 4.34, 4.34, 4.33, 3.75, 3.02, 2.31, 1.51), 7847, 16330, (8.89, 8.89, 8.89, 8.89, 8.89, 8.89, 8.89, 8.20, 6.89, 5.87, 5.10, 4.41), 6607, 16330),
    "60 mm   3 / 4, 8 x 0.9": ((4.22, 4.22, 4.22, 4.22, 4.22, 4.22, 4.22, 4.22, 4.22, 4.01, 3.49, 2.84), 10541, 16330, (8.66, 8.66, 8.66, 8.66, 8.66, 8.66, 8.66, 8.66, 8.36, 7.40, 6.52, 5.69), 8909, 16330),
}


def check_q(rep: Report) -> dict:
    """(q) the winding study's windings against its own numbers, through this module's calls in the study's configuration (the block above): (q1) the copper fill and R(100 C); (q2) the own
    continuous rating of each winding (the network with the maps' iron loss and copper_loss.py's full loss, the strands in a cluster, the diamond end turn, the shaft fan held to 50 W on the
    winding's own stack, N45 at 80 C or 150 C); (q3) the torque-speed table of the five 60 mm windings (the continuous current flat at the study's own rating and its double flat, k = 1, air-gap
    torque, copper 100 C): every torque within the two decimals it printed, the corner and the last speed with torque within the integer rpm it printed. The currents are the study's own
    ratings, so that (q2) judges the rating and (q3) the solver given the current (the first version of (q3) was fed this module's own ratings, which differ from the study's by up to 3e-4:
    its envelopes differed by 0.006 N.m and 1.8 rpm, and the check was redesigned, its limits unchanged). Nothing of the run's own inputs enters."""
    from optimization import thermal_machine as tm
    asm = tm.assumptions()
    copper = mm.FullCopper(f_sw=20000.0, arrangement="cluster")
    worst_fill, worst_r, worst_at, worst_t, worst_rpm, n_t, n_rpm = 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0
    ratings = {}
    rep.say("  (q) the winding study's windings through this module's calls, in the study's configuration (fan (b) 50 W sized on each stack with the diamond end turn, copper_loss.py full loss, cluster):")
    rep.say("      winding                                    fill  (study)   R(100 C) [mOhm]  (study)   own rating [A-turns]  (study)     peak I [A]")
    for label, stack, w, mag, at_study, r_study, fill_study in STUDY_WINDINGS:
        cool = thermal.size_shaft_fan(50.0, stack_mm=stack, assume=asm)["cooling"]
        fan = (lambda r, c=cool, x=stack: thermal.shaft_fan_power(r, c, x, assume=asm))
        m = mm.Machine(w, stack, copper=copper, mech=mm.Mechanical(fan=fan))
        r_mohm = 1e3 * float(m.resistance(100.0))
        rate = mt.thermal_rating(m, mt.ThermalConfig(cool, asm, magnet_limit_c=mag), RATED_RPM, scan=mm.Scan(m, i_max=m.i_cap))
        ratings[label] = (m, rate)
        worst_at = max(worst_at, abs(rate["AT_pk"] / at_study - 1.0))
        if fill_study is not None:
            worst_fill = max(worst_fill, abs(m.fill / fill_study - 1.0))
        if r_study is not None:
            worst_r = max(worst_r, abs(r_mohm / r_study - 1.0))
        rep.say(f"      {label:38s} {m.fill:8.4f}  {'-' if fill_study is None else f'{fill_study:.4f}':>7s}   {r_mohm:10.3f}  {'-' if r_study is None else f'{r_study:.3f}':>10s}   {rate['AT_pk']:12.3f}  {at_study:10.3f}   "
                f"{rate['I']:8.2f} (2 x {at_study * w.paths / w.turns:.2f})")
    bus = mm.Bus(48.0, k=1.0)
    rep.say("  (q) the torque-speed tables of the five 60 mm windings (study configuration: k = 1, air-gap torque, flat currents at the rating and twice it, copper 100 C, 48 V), N.m at "
            + ", ".join(f"{r / 1e3:g}" for r in STUDY_RPMS) + " krpm: this module / the study")
    for label, (cont, c_corner, c_zero, st_t, s_corner, s_zero) in STUDY_ENVELOPES.items():
        m, rate = ratings[label]
        w, at_study = m.w, next(x[4] for x in STUDY_WINDINGS if x[0] == label)
        i_study = at_study * w.paths / w.turns                  # the study's own rating as the current: (q2) holds the rating, (q3) the table given it
        sc = mm.Scan(m, i_max=min(2.0 * i_study, m.i_cap))
        for kind, lim, ref_t, ref_corner, ref_zero in (("continuous", i_study, cont, c_corner, c_zero), ("short-term", 2.0 * i_study, st_t, s_corner, s_zero)):
            env = mm.envelope(m, bus, np.asarray(STUDY_RPMS), min(lim, m.i_cap), basis="air_gap", scan=sc, locate=True)
            mine = np.asarray(env["T"], float)
            worst_t = max(worst_t, float(np.max(np.abs(mine - np.asarray(ref_t)))))
            n_t += len(ref_t)
            corner, zero = env.meta.get("corner_rpm"), env.meta.get("zero_rpm")
            for got, want in ((corner, ref_corner), (zero, ref_zero)):
                if got is not None:
                    worst_rpm, n_rpm = max(worst_rpm, abs(got - want)), n_rpm + 1
            rep.say(f"      {label:26s} {kind:10s} {lim:6.1f} A: " + " ".join(f"{x:5.2f}" for x in mine) + f"  corner {fmt(corner)} / {ref_corner}, last speed {fmt(zero)} / {ref_zero}")
            rep.say(f"      {'':26s} {'study':10s} {'':8s}  " + " ".join(f"{x:5.2f}" for x in ref_t))
    rep.check("(q1) the winding study's windings: copper fill and R(100 C)", worst_fill <= TOL_Q_FILL and worst_r <= TOL_Q_R,
              f"fill within {worst_fill:.1e} (limit {TOL_Q_FILL:g}), R(100 C) within {worst_r:.1e} (limit {TOL_Q_R:g}) of the study's, over {len(STUDY_WINDINGS)} windings")
    rep.check("(q2) the winding study's own continuous ratings", worst_at <= TOL_Q_AT, f"{len(STUDY_WINDINGS)} windings (the four stacks, the front's other members, the 150 C grade): largest relative difference {worst_at:.1e} (limit {TOL_Q_AT:g})")
    rep.check("(q3) the winding study's torque-speed tables", worst_t <= TOL_Q_T and worst_rpm <= TOL_Q_RPM,
              f"{n_t} torques within {worst_t:.4f} N.m of the study's two decimals (limit {TOL_Q_T:g}), {n_rpm} corner and last speeds within {worst_rpm:.2f} rpm (limit {TOL_Q_RPM:g})")
    return {"worst_fill": worst_fill, "worst_r": worst_r, "worst_at": worst_at, "worst_t": worst_t, "worst_rpm": worst_rpm}


# ====================================================================================================== the run
def map_variants() -> list[str]:
    """The windings of the headline stack that differ: the headline variant, and the other one when its winding is not the same (it is, by default)."""
    h = winding_of(HEADLINE_STACK, HEADLINE_VARIANT)
    return [HEADLINE_VARIANT] + [v for v in VARIANTS if v != HEADLINE_VARIANT and winding_of(HEADLINE_STACK, v) != h]


def wlabel(w: mm.Winding) -> str:
    """A winding in the winding study's notation: N / a, strands x wire."""
    return f"{w.turns:g} / {w.paths}, {w.strands:g} x {w.wire_mm:g} mm"


def print_inputs(rep: Report) -> None:
    """The assumptions of the run, printed at its head (every one is an input: changing it changes a line of this file or a command-line option, not the engine)."""
    cm = copper_model()
    rep.say("THE INPUTS (the windings are the winding study's, optimization/winding_study.py, Storage 0acd770; the rest are assumptions; each is swapped by changing an input, not code)")
    for variant, what in (("opt", "'opt', the winding study's choice for each stack"), ("same", "'same', the base winding on every stack")):
        parts = []
        for stack in STACKS:
            w = winding_of(stack, variant)
            m = mm.Machine(w, stack, copper=cm)
            parts.append(f"{stack:g} mm {wlabel(w)} ({w.wire_over_mm:g} mm over the enamel, fill {m.fill:.3f}, R(100 C) {1e3 * float(m.resistance(100.0)):.2f} mOhm)")
        rep.say(f"  windings, {what}: " + "; ".join(parts))
    rep.say(f"  copper model: {cm.name}" + (f", strands in a {ARRANGEMENT} bundle" if COPPER == "full" else "") + f"; copper temperature: the thermal network's (mean of the slot and end-turn nodes) on the "
            f"continuous envelopes; {T_CU_SHORT:g} C on the short-term ones (the end of a pulse: the larger R)")
    rep.say("  bus: " + "; ".join(f"{v:g} V -> V_max = {K_MARGIN} x {v:g} / sqrt 3 = {bus_of(v).v_max:.2f} V" for v in VDCS) + ". The margin k = 0.95 is an assumption: dead time about 1.7 % of V_max "
            "(0.5 us at 20 kHz), the switches' and cable's drops about 0.3 %, and about 3 % of reserve so that the current regulators never saturate and for the bus sagging under load")
    rep.say(f"  stacks {', '.join(f'{x:g}' for x in STACKS)} mm (rows of 20 mm), headline {HEADLINE_STACK:g} mm; top speed {mm.TOP_RPM:,.0f} rpm (the rotor's safe speed, a design decision of 2026-11-09)")
    rep.say(f"  current limits: continuous = the thermal network's limit at every speed (N45 magnets at or under {MAGNET_LIMIT_C:g} C, winding hot spot under {thermal.CLASS_F_C:g} C, ambient 40 C), found with the "
            f"iron loss of the maps at the operating point; short-term = {PEAK_FACTOR:g} x the winding's own continuous current at 3,000 rpm with no voltage limit (its thermal rating), constant over the "
            "speed, a COLD start; both clipped at the edge of the flux maps (3 x the placeholder's rated current). The winding study sizes the inverter at 2 x a COMMON rating (5 % under the "
            "example winding's), so its short-term currents are a few percent under these")
    rep.say("  iron loss: the maps' loss at the CAD's 0.3 mm lamination (the classical term scaled by (0.3 / 0.356)^2; hysteresis and excess unscaled); torque: the d-q model's, shaft torque = air-gap torque less the iron and mechanical loss torques")
    rep.say("  mechanical loss: bearings (mu 0.0015, 60 N each, 15 mm bore) and windage (thermal.py's formulas) plus the shaft fan's power, a rough estimate: expect a factor of two either way")
    for c in coolings():
        extra = ""
        if c["fan_factory"]:
            f = c["fan_factory"](HEADLINE_STACK)
            extra = "; fan power " + ", ".join(f"{float(f(r)):,.1f} W at {r:,.0f} rpm" for r in (1000.0, 3000.0, 6000.0, 10000.0, mm.TOP_RPM))
            if c["key"] != "fanA":
                extra += "; sized on each stack: " + ", ".join(f"{c['cooling_of'](x).v_ref_ms:.3f}" for x in STACKS) + f" m/s at 3,000 rpm on {', '.join(f'{x:g}' for x in STACKS)} mm"
        elif c["cooling"].kind == "forced_air":
            extra = f"; its fan takes {thermal.get_model(None, HEADLINE_STACK, c['cooling'], 3000.0, 0.45, thermal_assumptions()).fan_power():.1f} W, not counted in the motor's efficiency"
        rep.say(f"  cooling: {c['label']}{extra}")
    if EXTRAS_ON:
        rep.say("  extras on the headline stack with the headline fan: " + "; ".join(f"{e['label']} ({wlabel(e['winding'])}, magnets {e['magnet_limit_c']:g} C)" for e in EXTRAS))
    rep.say("  NOT modelled: the magnets' eddy loss (an input of the network, 0 here), the inverter's loss and ripple, the end-winding leakage (not in the 2D maps), demagnetisation at large -i_d (EM-6), "
            + ("the winding's AC loss on the DC model, " if COPPER == "dc" else "") + "anything beyond 3 x the placeholder's rated current, the magnets' temperature in the flux maps and the "
            "lower remanence of a 150 C grade (the maps are N45's at its nominal remanence)")


def run_all(rep: Report, label: str, quick: bool) -> dict:
    """Every case: the shaft fans on every stack in both windings, the comparison coolings on the headline stack, then the extras; the headline stack on the fine grid, the others on the coarse one.
    Where the two windings of a stack are the same winding (the headline stack by default) the case is run once."""
    t0 = time.time()
    print_inputs(rep)
    results: dict = {"inputs": inputs_dict(), "cases": {}, "static": {}}

    def record(cool, stack, variant, cr):
        key = f"{cool['key']}|{stack:g}|{variant}"
        results["static"][key] = cr.static()
        for vdc in VDCS:
            results["cases"][f"{key}|{vdc:g}"] = cr.run(vdc)
        s48 = results["cases"][f"{key}|{VDCS[0]:g}"]
        rep.say(f"  {cool['short']:34s} {stack:3.0f} mm {variant:5s} {wlabel(cr.w):22s}: rating {cr.rating['AT_pk']:6.1f} A-turns ({cr.rating['I']:5.1f} A, {cr.rating['J_rms']:4.1f} A/mm2), peak {cr.i_peak:6.1f} A"
                f"{' (maps edge)' if cr.peak_clipped else ''}; {VDCS[0]:g} V short-term corner {fmt(s48['short']['summary']['corner_rpm'])} rpm, T {s48['short']['summary']['T_peak']:.2f} N.m; "
                f"[{time.time() - t0:.0f} s]")
    for cool in coolings():
        fan = cool["key"] in FAN_KEYS
        for stack in (STACKS if fan else (HEADLINE_STACK,)):
            done: dict = {}
            for variant in (VARIANTS if fan else (HEADLINE_VARIANT,)):
                w = winding_of(stack, variant)
                twin = done.get(w)
                if twin is not None:                                          # the same winding under another name: the same case
                    key, tkey = f"{cool['key']}|{stack:g}|{variant}", f"{cool['key']}|{stack:g}|{twin}"
                    results["static"][key] = dict(results["static"][tkey], variant=variant)
                    for vdc in VDCS:
                        results["cases"][f"{key}|{vdc:g}"] = results["cases"][f"{tkey}|{vdc:g}"]
                    continue
                done[w] = variant
                record(cool, stack, variant, CaseRunner(cool, stack, variant, stack == HEADLINE_STACK and not quick))
    if EXTRAS_ON:
        fan = cool_by_key(EXTRA_FAN)
        for ex in EXTRAS:
            record(dict(fan, short=ex["label"]), HEADLINE_STACK, ex["key"], CaseRunner(fan, HEADLINE_STACK, ex["key"], not quick, winding=ex["winding"], magnet_limit_c=ex["magnet_limit_c"]))
    return results


def fmt(x, digits: int = 0) -> str:
    return "none" if x is None else f"{x:,.{digits}f}"


def inputs_dict() -> dict:
    def wd(w):
        return {"turns": w.turns, "paths": w.paths, "strands": w.strands, "wire_mm": w.wire_mm, "wire_over_mm": w.wire_over_mm}
    return {"winding": wd(WINDING), "windings": {v: {f"{x:g}": wd(winding_of(x, v)) for x in STACKS} for v in VARIANTS}, "variants": list(VARIANTS), "headline_variant": HEADLINE_VARIANT,
            "stacks": list(STACKS), "headline_stack": HEADLINE_STACK,
            "vdc": list(VDCS), "k_margin": K_MARGIN, "t_cu0": T_CU0, "t_cu_short": T_CU_SHORT, "peak_factor": PEAK_FACTOR, "top_rpm": mm.TOP_RPM, "rated_rpm": RATED_RPM,
            "lamination_mm": 0.3, "copper_kind": COPPER, "arrangement": ARRANGEMENT, "magnet_limit_c": MAGNET_LIMIT_C,
            "copper": "DC, the first study's end turn (61.15 mm)" if COPPER == "dc" else f"copper_loss.py's full loss, strands in a {ARRANGEMENT} bundle, diamond end turn 48.55 mm",
            "provisional": False, "coolings": [c["key"] for c in coolings()],
            "cool_labels": {c["key"]: {"label": c["label"], "short": c["short"]} for c in coolings()}, "fan_keys": list(FAN_KEYS), "main_fan": MAIN_FAN,
            "fan_v_ref": {k: {f"{x:g}": float(fan_cooling(k, x).v_ref_ms) for x in STACKS} for k in FAN_KEYS},
            "extras": [{"key": e["key"], "label": e["label"], "winding": wd(e["winding"]), "magnet_limit_c": e["magnet_limit_c"]} for e in EXTRAS] if EXTRAS_ON else [], "extra_fan": EXTRA_FAN}


def cool_by_key(key: str) -> dict:
    return next(c for c in coolings() if c["key"] == key)


# ====================================================================================================== the efficiency maps
def run_maps(rep: Report, results: dict, quick: bool) -> None:
    """The efficiency maps of the headline stack: the three shaft fans, the headline winding (and the other variant's where it differs), both buses, at the short-term current limit (the map
    covers everything the machine can do for a moment; the continuous envelope is drawn inside it)."""
    results["maps"] = {}
    fractions = np.concatenate([np.linspace(0.0, 0.98, MAP_FRACTIONS - 1), [1.0]])
    t0 = time.time()
    for cool in coolings()[:3]:
        for variant in map_variants():
            cr = CaseRunner(cool, HEADLINE_STACK, variant, fine=True)
            sc = mm.Scan(cr.m, i_max=cr.i_peak, **(GRID_TABLE if quick else GRID_HEADLINE))
            for vdc in VDCS:
                key = f"{cool['key']}|{HEADLINE_STACK:g}|{variant}|{vdc:g}"
                env = results["cases"][key]["short"]
                on = np.nonzero(env["d"]["T"] > 0)[0]
                top = float(env["d"]["rpm"][on.max()]) if len(on) else mm.TOP_RPM
                rpms = np.linspace(150.0, top, MAP_SPEEDS)
                em = mm.efficiency_map(cr.m, bus_of(vdc), rpms, fractions, cr.i_peak, t_c=T_CU_SHORT, scan=sc)
                results["maps"][key] = em
                ok = em["minloss"]["found"]
                eta = em["minloss"]["eta"]
                good = ok & (em["T"] >= 0.1 * em["T_env"][:, None]) & (em["rpm"][:, None] >= 500.0) & np.isfinite(eta)
                rep.say(f"  map {cool['short']:34s} {variant:5s} {vdc:3.0f} V: {int(ok.sum())} of {ok.size} points solved, efficiency {100 * np.nanmin(np.where(good, eta, np.nan)):.1f} ... "
                        f"{100 * np.nanmax(np.where(good, eta, np.nan)):.1f} % (torque >= 10 % of the envelope, >= 500 rpm); minimum-loss control saves at most {np.nanmax(em['saving']):.2f} W "
                        f"({100 * np.nanmax(em['saving'] / em['mtpa']['P_in']):.2f} % of the input) [{time.time() - t0:.0f} s]")


# ====================================================================================================== the sensitivity to the voltage margin and the swap of the copper model
def run_extras(rep: Report, results: dict) -> None:
    """The voltage margin's effect (k = 1.00, 0.95, 0.90) on the short-term envelope of the headline stack, and the swap of the copper model (the run's, and the other one: copper_loss.py's
    full loss against the first study's DC loss): the same machine, a different CopperModel, nothing else changed."""
    fan = cool_by_key(MAIN_FAN)
    sens = {}
    swap = {}
    alt_kind = "dc" if COPPER == "full" else "full"
    for variant in map_variants():
        cr = CaseRunner(fan, HEADLINE_STACK, variant, fine=True)
        for k in (1.0, K_MARGIN, 0.90):
            env = mm.envelope(cr.m, mm.Bus(VDCS[0], k), cr.speeds, cr.i_peak, t_c=T_CU_SHORT, scan=cr.sc)
            sens[(variant, k)] = summarise(env)
        base = results["cases"][f"{MAIN_FAN}|{HEADLINE_STACK:g}|{variant}|{VDCS[0]:g}"]["short"]
        bus = bus_of(VDCS[0])
        alt = copper_model(alt_kind)
        mf = mm.Machine(cr.w, HEADLINE_STACK, copper=alt, mech=cr.m.mech)
        d = base["d"]
        on = np.isfinite(d["i_d"])
        pt_run = {k: np.asarray(v) for k, v in d.items()}
        pt_alt = mf.point(np.where(on, d["i_d"], 0.0), np.where(on, d["i_q"], 0.0), d["rpm"], bus, T_CU_SHORT)
        env_alt = mm.envelope(mf, bus, cr.speeds, cr.i_peak, t_c=T_CU_SHORT, scan=mm.Scan(mf, i_max=cr.i_peak, **GRID_HEADLINE))
        swap[variant] = {"alt": alt_kind, "R_run": float(cr.m.resistance(100.0)), "R_alt": float(mf.resistance(100.0)), "summary_run": summarise(mm.Envelope(base["d"], base["meta"])),
                         "summary_alt": summarise(env_alt), "rpm": d["rpm"], "P_cu_run": pt_run["P_cu"], "P_cu_alt": np.where(on, pt_alt["P_cu"], np.nan), "P_cu_ac_run": pt_run["P_cu_ac"],
                         "P_cu_ac_alt": np.where(on, pt_alt["P_cu_ac"], np.nan), "eta_run": pt_run["eta"], "eta_alt": np.where(on, pt_alt["eta"], np.nan)}
        j3 = int(np.argmin(np.abs(d["rpm"] - RATED_RPM)))
        j10 = int(np.argmin(np.abs(d["rpm"] - 10000.0)))
        sw = swap[variant]
        rep.say(f"  swap of the copper model ({wlabel(cr.w)}, 60 mm, {VDCS[0]:g} V, short-term): the run's {cr.m.copper.name.split(':')[0]} -> {alt.name.split(':')[0]}: R(100 C) {1e3 * sw['R_run']:.2f} -> "
                f"{1e3 * sw['R_alt']:.2f} mOhm; at 3,000 rpm the copper loss {pt_run['P_cu'][j3]:.1f} -> {sw['P_cu_alt'][j3]:.1f} W (AC + PWM part {pt_run['P_cu_ac'][j3]:.2f} -> {sw['P_cu_ac_alt'][j3]:.2f} W), "
                f"efficiency {100 * pt_run['eta'][j3]:.2f} -> {100 * sw['eta_alt'][j3]:.2f} %; at 10,000 rpm the copper loss {pt_run['P_cu'][j10]:.1f} -> {sw['P_cu_alt'][j10]:.1f} W (AC + PWM part "
                f"{pt_run['P_cu_ac'][j10]:.2f} -> {sw['P_cu_ac_alt'][j10]:.2f} W); envelope corner {fmt(sw['summary_run']['corner_rpm'])} -> {fmt(sw['summary_alt']['corner_rpm'])} rpm")
    results["sensitivity_k"], results["swap"] = sens, swap
    for variant in map_variants():
        rep.say(f"  voltage margin k, {wlabel(winding_of(HEADLINE_STACK, variant))}, 60 mm, {VDCS[0]:g} V, short-term: " + "; ".join(
            f"k {k:.2f}: corner {fmt(sens[(variant, k)]['corner_rpm'])} rpm, T at 3,000 rpm {sens[(variant, k)]['T_rated']:.2f}, top speed {sens[(variant, k)]['T_top']:.2f} N.m" for k in (1.0, K_MARGIN, 0.90)))


# ====================================================================================================== the checks over every solution
def check_points(rep: Report, results: dict) -> None:
    """(c) the energy balance at every point of every envelope and map, (o) the constraints at every solution, (n) the level-set solver's torque error and its ordering, (p) the
    thermal limits along the continuous envelopes."""
    worst_alg, worst_age, n_pts = 0.0, -1.0, 0
    worst_v, worst_i = 0.0, 0.0
    max_mag, max_hot, hottest_magnet = -1e9, -1e9, -1e9
    for key, case in results["cases"].items():
        cool, stack, variant, vdc = key.split("|")
        k = float(stack) / 20.0
        vmax = bus_of(float(vdc)).v_max
        for kind in ("cont", "short"):
            d = case[kind]["d"]
            ok = np.isfinite(d["P_in"]) & (d["P_in"] > 1e-9)
            if not ok.any():
                continue
            bal = np.abs(d["P_in"][ok] - (d["P_out"][ok] + d["P_cu"][ok] + d["P_fe"][ok] + d["P_mech"][ok])) / d["P_in"][ok]
            worst_alg, n_pts = max(worst_alg, float(bal.max())), n_pts + int(ok.sum())
            lim = TOL_C_AGE_REL * np.abs(d["T_age"][ok]) + TOL_C_AGE_ABS * k
            worst_age = max(worst_age, float(((np.abs(d["T_em"][ok] - d["T_age"][ok])) / lim).max()))
            worst_v = max(worst_v, float(np.max(d["v_abs"][ok] / vmax - 1.0)))
            worst_i = max(worst_i, float(np.max(d["I"][ok] / d["i_limit"][ok] - 1.0)))
        lim_c = case["limits"]
        temp = lim_c["limited"] == "temperature"
        if temp.any():
            mag_limit = results["static"][f"{cool}|{stack}|{variant}"]["magnet_limit_c"]
            max_mag = max(max_mag, float(np.nanmax(lim_c["magnet"][temp])) - mag_limit)
            hottest_magnet = max(hottest_magnet, float(np.nanmax(lim_c["magnet"][temp])))
            max_hot = max(max_hot, float(np.nanmax(lim_c["hot_spot"][temp])))
    for key, em in results.get("maps", {}).items():
        stack, vdc = float(key.split("|")[1]), float(key.split("|")[3])
        k = stack / 20.0
        vmax = bus_of(vdc).v_max
        for ctl in ("mtpa", "minloss"):
            d = em[ctl]
            ok = d["found"] & np.isfinite(d["P_in"]) & (d["P_in"] > 1e-9)
            bal = np.abs(d["P_in"][ok] - (d["P_out"][ok] + d["P_cu"][ok] + d["P_fe"][ok] + d["P_mech"][ok])) / d["P_in"][ok]
            worst_alg, n_pts = max(worst_alg, float(bal.max())), n_pts + int(ok.sum())
            lim = TOL_C_AGE_REL * np.abs(d["T_age"][ok]) + TOL_C_AGE_ABS * k
            worst_age = max(worst_age, float((np.abs(d["T_em"][ok] - d["T_age"][ok]) / lim).max()))
            worst_v = max(worst_v, float(np.max(d["v_abs"][ok] / vmax - 1.0)))
            worst_i = max(worst_i, float(np.max(d["I"][ok] / em["i_limit"][:, None].repeat(d["I"].shape[1], axis=1)[ok] - 1.0)))
    rep.check("(c1) the energy balance closes: P_in = P_out + copper + iron + mechanical", worst_alg <= TOL_C_ALG, f"{n_pts:,} operating points (every envelope and map point): largest |P_in - sum| / P_in {worst_alg:.1e} "
              f"(limit {TOL_C_ALG:g}); P_in is 1.5 (v_d i_d + v_q i_q) from the voltage equations, the rest from the torque, the maps' losses and the mechanical estimate")
    rep.check("(c2) the d-q torque against FEMM's air-gap torque at every point", worst_age <= 1.0, f"|T_dq - T_age| is at most {worst_age:.2f} of its limit ({TOL_C_AGE_REL:g} T + {TOL_C_AGE_ABS:g} N.m per 20 mm), "
              f"so the input power from the voltage equations and the one from the air-gap torque agree to that")
    rep.check("(o) every solution is inside both limits", worst_v <= TOL_O and worst_i <= TOL_O, f"|v| / V_max - 1 at most {worst_v:.1e}, |i| / limit - 1 at most {worst_i:.1e} (limit {TOL_O:g})")
    rep.check("(p) the continuous envelopes stay inside the temperature limits", max_mag <= 0.01 and max_hot <= thermal.CLASS_F_C + 0.01,
              f"the magnets at most {max_mag:+.2f} K above their limit ({MAGNET_LIMIT_C:g} C; {max([e['magnet_limit_c'] for e in EXTRAS], default=MAGNET_LIMIT_C):g} C on the 150 C grade extra), "
              f"the hottest winding spot {max_hot:.2f} C (limit {thermal.CLASS_F_C:g}), the hottest magnet {hottest_magnet:.2f} C, over every speed whose limit is the temperature")
    worst_t, worst_order, n_t, unsolved = 0.0, 0.0, 0, 0
    for key, em in results.get("maps", {}).items():
        for ctl in ("mtpa", "minloss"):
            d = em[ctl]
            ok = d["found"]
            worst_t = max(worst_t, float(np.nanmax(np.abs(d["T_shaft"][ok] - em["T"][ok]))))
            n_t += int(ok.sum())
            unsolved += int((~ok & (em["T"] > 0)).sum())
        both = em["mtpa"]["found"] & em["minloss"]["found"]
        worst_order = max(worst_order, float(np.nanmax(((em["minloss"]["P_loss"] - em["mtpa"]["P_loss"]) / em["mtpa"]["P_loss"])[both])))
    rep.check("(n) the level-set solver hits its torque, and minimum-loss control is never worse than the MTPA control", worst_t <= TOL_N and worst_order <= 1e-9,
              f"{n_t:,} solved points: the shaft torque is within {worst_t:.1e} N.m of its target (limit {TOL_N:g}); the minimum-loss point's loss less the MTPA control's, relative, at most {worst_order:.1e} "
              f"(0 or below is right); {unsolved} points under the envelope were not solved")


# ====================================================================================================== the stack table
STACK_COLUMNS = ["cooling", "stack_mm", "rows", "magnets_to_buy", "winding", "winding_label", "magnet_limit_C", "turns", "series_turns", "strand_mm", "fill", "R_100C_ohm", "rating_AT_pk", "rating_I_A", "rating_J_A_mm2", "peak_I_A",
                 "peak_clipped_at_maps_edge", "m_stator_steel_kg", "m_rotor_steel_kg", "m_magnets_g", "m_copper_kg", "m_active_kg", "J_rotor_g_m2", "emf_top_pk_V", "rpm_emf_vmax",
                 "vdc", "short_corner_rpm", "short_T_corner_Nm", "short_T_peak_Nm", "short_P_peak_W", "short_P_peak_rpm", "short_zero_rpm", "short_T_3000_Nm", "short_T_top_Nm", "short_P_top_W",
                 "cont_corner_rpm", "cont_T_3000_Nm", "cont_P_peak_W", "cont_P_peak_rpm", "cont_zero_rpm", "overload_magnet_s", "overload_winding_s"]


def stack_rows(results: dict) -> list[dict]:
    rows = []
    cools = {c["key"]: c for c in coolings()}
    for key, stat in results["static"].items():
        cool, stack, variant = key.split("|")
        g = stat["geometry"]
        for vdc in VDCS:
            case = results["cases"][f"{key}|{vdc:g}"]
            sc_, cc = case["short"]["summary"], case["cont"]["summary"]
            ov = case["overload"]
            wst = stat["winding"]
            rows.append({"cooling": cools[cool]["short"] if cool in cools else cool, "stack_mm": float(stack), "rows": g["rows"], "magnets_to_buy": g["n_magnets"], "winding": variant,
                         "winding_label": f"{wst['turns']:g} / {wst['paths']}, {wst['strands']:g} x {wst['wire_mm']:g} mm", "magnet_limit_C": stat["magnet_limit_c"], "turns": stat["winding"]["turns"],
                         "series_turns": stat["winding"]["series_turns"], "strand_mm": stat["winding"]["wire_mm"], "fill": stat["winding"]["fill"], "R_100C_ohm": stat["R_100C"],
                         "rating_AT_pk": stat["rating"]["AT_pk"], "rating_I_A": stat["rating"]["I"], "rating_J_A_mm2": stat["rating"]["J_rms"], "peak_I_A": stat["i_peak"],
                         "peak_clipped_at_maps_edge": stat["peak_clipped"], "m_stator_steel_kg": g["m_stator_steel"], "m_rotor_steel_kg": g["m_rotor_steel"], "m_magnets_g": 1e3 * g["m_magnets"],
                         "m_copper_kg": g["m_copper"], "m_active_kg": g["m_active"], "J_rotor_g_m2": 1e3 * g["J_rotor"], "emf_top_pk_V": stat["emf_top_peak"],
                         "rpm_emf_vmax": stat["rpm_emf_vmax"][f"{vdc:g}"], "vdc": vdc,
                         "short_corner_rpm": sc_["corner_rpm"], "short_T_corner_Nm": sc_["corner_T"], "short_T_peak_Nm": sc_["T_peak"], "short_P_peak_W": sc_["P_peak"], "short_P_peak_rpm": sc_["P_peak_rpm"],
                         "short_zero_rpm": sc_["zero_rpm"], "short_T_3000_Nm": sc_["T_rated"], "short_T_top_Nm": sc_["T_top"], "short_P_top_W": sc_["P_top"],
                         "cont_corner_rpm": cc["corner_rpm"], "cont_T_3000_Nm": cc["T_rated"], "cont_P_peak_W": cc["P_peak"], "cont_P_peak_rpm": cc["P_peak_rpm"], "cont_zero_rpm": cc["zero_rpm"],
                         "overload_magnet_s": ov.get("magnet_s"), "overload_winding_s": ov.get("winding_s")})
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def write_envelopes(path: Path, results: dict) -> None:
    """Every envelope as a long CSV: one row a speed, with the case in the first columns."""
    keys = ["rpm", "T", "P_out", "region", "I", "gamma", "i_d", "i_q", "v_abs", "T_em", "P_in", "P_cu", "P_cu_dc", "P_cu_ac", "P_fe_teeth", "P_fe_yoke", "P_fe_rotor", "P_bearing", "P_windage", "P_fan",
            "P_mech", "eta", "PF", "i_limit", "t_cu"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["case", "kind"] + keys)
        for case_key, case in results["cases"].items():
            for kind in ("cont", "short"):
                d = case[kind]["d"]
                for j in range(len(d["rpm"])):
                    w.writerow([case_key, kind] + [(f"{d[k][j]:.6g}" if k != "region" else str(d[k][j])) for k in keys])


# ====================================================================================================== the headline
def headline(rep: Report, results: dict) -> dict:
    """The numbers of the report: the headline (60 mm, the winding study's winding) in every cooling at both buses."""
    rep.say("")
    rep.say(f"THE HEADLINE: {wlabel(winding_of(HEADLINE_STACK, HEADLINE_VARIANT))}, 60 mm, k = {K_MARGIN}, copper: {copper_model().name.split(':')[0]} at the network's temperature (continuous) / {T_CU_SHORT:g} C "
            f"(short-term), iron loss at 0.3 mm, shaft torque (the numbers are estimates)")
    out = {}
    for cool in coolings():
        for variant in (map_variants() if cool["key"] in FAN_KEYS else [HEADLINE_VARIANT]):
            stat = results["static"][f"{cool['key']}|{HEADLINE_STACK:g}|{variant}"]
            rep.say(f"  {cool['label']}, winding {wlabel(winding_of(HEADLINE_STACK, variant))}: thermal rating {stat['rating']['AT_pk']:.1f} A-turns = {stat['rating']['I']:.1f} A, peak {stat['i_peak']:.1f} A"
                    f"{' (clipped at the maps edge)' if stat['peak_clipped'] else ''}")
            for vdc in VDCS:
                c = results["cases"][f"{cool['key']}|{HEADLINE_STACK:g}|{variant}|{vdc:g}"]
                for kind, name in (("short", "short-term (cold start)"), ("cont", "continuous")):
                    sm = c[kind]["summary"]
                    corner = f"corner {fmt(sm['corner_rpm']):>6s} rpm at {fmt(sm['corner_T'], 2):>6s} N.m" if sm["corner_rpm"] is not None else "no corner (the voltage never limits)"
                    rep.say(f"      {vdc:3.0f} V {name:24s}: {corner}; T at 3,000 rpm {sm['T_rated']:5.2f} N.m; "
                            f"at 16,330 rpm {sm['T_top']:5.2f} N.m = {sm['P_top']:6.0f} W; peak power {sm['P_peak']:6.0f} W at {sm['P_peak_rpm']:6.0f} rpm; torque ends at {ends(sm)} rpm; "
                            f"efficiency along the envelope {100 * sm['eta_min']:.1f} ... {100 * sm['eta_max']:.1f} %")
                d_ = c["short"]["d"]
                s_ = stat["winding"]["turns"] / (10.0 * stat["winding"]["paths"])
                i_d0 = float(np.nanmin(d_["i_d"])) * s_ / mm.load_maps().rated
                rep.say(f"      {vdc:3.0f} V deepest field weakening on the short-term envelope: i_d = {i_d0:.2f} x the placeholder's rated current (in the maps' units; the maps end at -3.00); "
                        f"the demagnetisation check (EM-6) is not made here")
                ov = c["overload"]
                if ov.get("feasible"):
                    when = lambda x: "not within 3,600 s" if x is None else f"after {x:,.0f} s"
                    rep.say(f"      {vdc:3.0f} V overload from cold at 3,000 rpm ({stat['i_peak']:.1f} A): magnets reach {stat['magnet_limit_c']:g} C {when(ov['magnet_s'])}, the winding {thermal.CLASS_F_C:g} C {when(ov['winding_s'])}")
                out[(cool["key"], variant, vdc)] = {"short": c["short"]["summary"], "cont": c["cont"]["summary"]}
    return out


def print_stack_table(rep: Report, results: dict) -> None:
    """The stack table as text (the CSV and fig_stack_table have it all), for each shaft fan, at the lowest bus voltage, both windings; and the open-circuit voltage at the top speed."""
    for cool in FAN_KEYS:
        label = next(c["label"] for c in coolings() if c["key"] == cool)
        rep.say("")
        rep.say(f"THE STACK TABLE ({label}; {VDCS[0]:g} V; 'opt' = the winding study's winding of each stack, 'same' = the headline winding on every stack; short-term = {PEAK_FACTOR:g} x the winding's own "
                f"rating at 3,000 rpm, from a cold start)")
        rep.say("  stack winding  N / a, strands x wire        peak I  corner rpm  T_peak  P_peak (rpm)       short-term ends   cont. T at 3,000 rpm  cont. ends  overload from cold (magnets / winding)   magnets  active kg   J rotor g.m2")
        for variant in VARIANTS:
            for stack in STACKS:
                key = f"{cool}|{stack:g}|{variant}"
                if key not in results["static"]:
                    continue
                st_, c = results["static"][key], results["cases"][f"{key}|{VDCS[0]:g}"]
                sm, cm, ov, g = c["short"]["summary"], c["cont"]["summary"], c["overload"], st_["geometry"]
                when = lambda x: ">3,600" if x is None else f"{x:,.0f}"
                ovs = f"{when(ov.get('magnet_s'))} s / {when(ov.get('winding_s'))} s" if ov.get("feasible") else "-"
                wst = st_["winding"]
                wtxt = f"{wst['turns']:g} / {wst['paths']}, {wst['strands']:g} x {wst['wire_mm']:g} mm"
                rep.say(f"  {stack:3.0f}  {variant:5s} {wtxt:24s} {st_['i_peak']:6.1f}{'*' if st_['peak_clipped'] else ' '} {fmt(sm['corner_rpm']):>9s}  "
                        f"{sm['T_peak']:6.2f}  {sm['P_peak'] / 1e3:5.2f} kW ({sm['P_peak_rpm']:6.0f})  {ends(sm):>9s}         {cm['T_rated']:6.2f}         {ends(cm):>9s}    {ovs:>26s}          "
                        f"{g['n_magnets']:4.0f}   {g['m_active']:6.2f}     {1e3 * g['J_rotor']:6.2f}")
    rep.say("  (* the peak current is clipped at the edge of the flux maps; t after a speed: the envelope ends there because no current both holds the voltage and keeps the temperatures, not because the loads take all the torque)")
    cool = MAIN_FAN
    rep.say("  open circuit at the top speed (the inverter off): " + "; ".join(
        f"{stack:g} mm {variant}: {results['static'][f'{cool}|{stack:g}|{variant}']['emf_top_peak']:.0f} V a phase, {results['static'][f'{cool}|{stack:g}|{variant}']['emf_top_line_line_peak']:.0f} V line to line (peak)"
        for variant in VARIANTS for stack in STACKS if f"{cool}|{stack:g}|{variant}" in results["static"]) + f", against a {VDCS[0]:g} / {VDCS[-1]:g} V bus")
    rep.say("  the no-load voltage reaches V_max (" + f"{VDCS[0]:g} V) at " + "; ".join(
        f"{stack:g} mm {variant}: {results['static'][f'{cool}|{stack:g}|{variant}']['rpm_emf_vmax'][f'{VDCS[0]:g}']:,.0f} rpm" for variant in VARIANTS for stack in STACKS
        if f"{cool}|{stack:g}|{variant}" in results["static"]) + ": above it the inverter must hold a field-weakening current at every speed")


def print_speed_limits(rep: Report, results: dict) -> None:
    """THE HIGHEST SPEED each fan allows: where the shaft torque ends, i.e. where the motor's air-gap torque (field weakening) meets the loads' torques: the fan's, the windage's, the bearings'
    and the iron loss's drag. For every shaft fan, winding and bus, on the short-term and on the continuous envelope."""
    rep.say("")
    rep.say(f"THE HIGHEST SPEED THE MOTOR REACHES WITH EACH FAN ON ITS SHAFT ({HEADLINE_STACK:g} mm): the speed where the shaft torque ends, the air-gap torque there and what it is spent on [N.m]")
    for cool in FAN_KEYS:
        label = next(c["label"] for c in coolings() if c["key"] == cool)
        rep.say(f"  {label}")
        for variant in map_variants():
            wtxt = wlabel(winding_of(HEADLINE_STACK, variant))
            for vdc in VDCS:
                c = results["cases"][f"{cool}|{HEADLINE_STACK:g}|{variant}|{vdc:g}"]
                parts = []
                for kind, name in (("short", "short-term"), ("cont", "continuous")):
                    sm = c[kind]["summary"]
                    z, lp = sm["zero_rpm"], c["limit_" + kind]
                    if sm["stop"] == "top":
                        parts.append(f"{name}: reaches {mm.TOP_RPM:,.0f} rpm (T {sm['T_top']:.2f} N.m = {sm['P_top']:,.0f} W there)")
                    elif sm["stop"] == "limit" or lp is None:
                        parts.append(f"{name}: about {z:,.0f} rpm, ended by the temperatures (no current both holds the voltage and keeps the limits; not by the loads)")
                    else:
                        parts.append(f"{name}: {z:,.0f} rpm, ended by the loads: air-gap {lp['T_em']:.3f} = fan {lp['T_fan']:.3f} + windage {lp['T_windage']:.3f} + bearings {lp['T_bearing']:.3f} + iron "
                                     f"{lp['T_iron']:.3f} ({lp['region']}, {lp['I']:.1f} A)")
                rep.say(f"      {wtxt}, {vdc:g} V: " + "; ".join(parts))


def print_extras(rep: Report, results: dict) -> None:
    """The front's other windings and the 150 C magnet grade on the headline stack with the headline fan, next to the headline winding: what the inverter's current class buys (the corner speed,
    the torque and power at speed) and what the magnet grade buys."""
    if not results["inputs"].get("extras"):
        return
    fan = cool_by_key(EXTRA_FAN)
    rows = [(wlabel(winding_of(HEADLINE_STACK, HEADLINE_VARIANT)) + "  (headline)", HEADLINE_VARIANT, MAGNET_LIMIT_C, STUDY_I_ST.get("opt"))] + [(wlabel(e["winding"]), e["key"], e["magnet_limit_c"], STUDY_I_ST.get(e["key"])) for e in EXTRAS]
    rep.say("")
    rep.say(f"THE FRONT'S OTHER WINDINGS AND THE 150 C GRADE ({HEADLINE_STACK:g} mm, {fan['label']}; short-term = {PEAK_FACTOR:g} x the winding's own rating; the study's I_st is 2 x its COMMON rating "
            f"(91.7 A-turns), a few percent under)")
    for vdc in VDCS:
        rep.say(f"  {vdc:g} V")
        rep.say("    winding                            magnets  rating [A-turns (A)]   peak I [A] (study)   short-term: corner [rpm] (T)    T at 8,000 / 10,000 / 16,330 rpm [N.m]   P at 16,330 [W]   continuous: T at 3,000 [N.m]  ends [rpm]   overload from cold [s]")
        for name, variant, mag, i_study in rows:
            key = f"{EXTRA_FAN}|{HEADLINE_STACK:g}|{variant}"
            st_, c = results["static"][key], results["cases"][f"{key}|{vdc:g}"]
            sm, cm, ov = c["short"]["summary"], c["cont"]["summary"], c["overload"]
            d = c["short"]["d"]
            t_at = [float(np.interp(r, d["rpm"], d["T"])) for r in (8000.0, 10000.0, 16330.0)]
            when = lambda x: ">3,600" if x is None else f"{x:,.0f}"
            ovs = f"{when(ov.get('magnet_s'))} / {when(ov.get('winding_s'))}" if ov.get("feasible") else "-"
            rep.say(f"    {name:34s} {mag:5.0f} C  {st_['rating']['AT_pk']:7.1f} ({st_['rating']['I']:6.1f})   {st_['i_peak']:7.1f} ({'-' if i_study is None else f'{i_study:.1f}':>5s})      "
                    f"{fmt(sm['corner_rpm']):>6s} ({sm['corner_T']:.2f})         {t_at[0]:6.2f} / {t_at[1]:6.2f} / {t_at[2]:6.2f}            {sm['P_top']:7.0f}           {cm['T_rated']:6.2f}   {ends(cm):>9s}         {ovs:>14s}")


def print_net_output(rep: Report, results: dict) -> None:
    """THE NET SHAFT OUTPUT of the headline winding: the torque and the power at the shaft (after the iron drag, the windage, the bearings and the shaft fan) at the speeds a reader asks for,
    for every shaft fan, bus and kind of envelope, with the fan's own power at those speeds."""
    speeds = (1000.0, 3000.0, 6000.0, 10000.0, mm.TOP_RPM)
    rep.say("")
    rep.say(f"THE NET SHAFT OUTPUT ({wlabel(winding_of(HEADLINE_STACK, HEADLINE_VARIANT))}, {HEADLINE_STACK:g} mm): shaft torque [N.m] and power [kW] after the iron drag, the windage, the bearings and the "
            "shaft fan, at " + " / ".join(f"{r:,.0f}" for r in speeds) + " rpm (0 where the envelope has ended)")
    for cool in coolings()[:3]:
        f = fan_of(cool, HEADLINE_STACK)
        rep.say(f"  {cool['label']}: the fan takes " + " / ".join(f"{float(f(r)):,.1f}" for r in speeds) + " W")
        for vdc in VDCS:
            c = results["cases"][f"{cool['key']}|{HEADLINE_STACK:g}|{HEADLINE_VARIANT}|{vdc:g}"]
            for kind, name in (("cont", "continuous"), ("short", "short-term")):
                d = c[kind]["d"]
                T = [float(np.interp(r, d["rpm"], d["T"])) for r in speeds]
                P = [float(np.interp(r, d["rpm"], np.nan_to_num(d["P_out"]))) / 1e3 for r in speeds]
                rep.say(f"      {vdc:3.0f} V {name:10s}: T " + " / ".join(f"{t:5.2f}" for t in T) + " N.m;  P " + " / ".join(f"{p:5.2f}" for p in P) + " kW")


# ====================================================================================================== the FEMM spot-check (through the pool and the slot)
def spot_points(results: dict) -> list[dict]:
    """Three deep field-weakening points of the headline envelopes, in the reference domain (the 10-turn winding, 20 mm): the largest |i_d| on the headline stack at the higher bus, a mid-speed
    point of the headline stack at the lower bus, and the largest |i_d| of the 20 mm stack's own winding (26 turns: a different scale)."""
    pts = []
    h, v0, v1 = HEADLINE_STACK, VDCS[0], VDCS[-1]
    for key, pick in ((f"{MAIN_FAN}|{h:g}|{HEADLINE_VARIANT}|{v1:g}", "deep"), (f"{MAIN_FAN}|{h:g}|{HEADLINE_VARIANT}|{v0:g}", "mid"), (f"{MAIN_FAN}|{STACKS[0]:g}|{HEADLINE_VARIANT}|{v0:g}", "deep")):
        d = results["cases"][key]["short"]["d"]
        stat = results["static"][key.rsplit("|", 1)[0]]
        s_ = stat["winding"]["turns"] / (10.0 * stat["winding"]["paths"])
        ok = np.isfinite(d["i_d"])
        idx = np.nonzero(ok)[0]
        j = idx[int(np.argmin(d["i_d"][idx]))] if pick == "deep" else idx[int(np.argmin(np.abs(d["rpm"][idx] - 6000.0)))]
        pts.append({"key": key, "rpm": float(d["rpm"][j]), "i_d": float(d["i_d"][j] * s_), "i_q": float(d["i_q"][j] * s_), "region": str(d["region"][j])})
    return pts


def spot_stage(rep: Report, out: Path, workers: int = 8) -> None:
    """(f) FEMM at the deep field-weakening points: psi_d, psi_q and the air-gap torque at the 15 rotor positions of one ripple period, averaged, against the maps' interpolation at
    the same point, through femm_queue.femm_pool (the slot; below-normal priority; a per-task timeout), one batch."""
    import femm_queue
    from optimization import maps_femm, maps_run
    with open(out / "results.pkl", "rb") as f:
        results = pickle.load(f)
    restore_inputs(results["inputs"])
    pts = spot_points(results)
    ctx = maps_run.Context()
    jobs = [maps_run.flux_job(ctx, "motor_map_spot", f"s{k}", "chosen", [(p["i_d"], p["i_q"])], maps_run.POSITIONS_FLUX) for k, p in enumerate(pts)]
    import femm_session
    rep.say(f"(f) FEMM spot-check: {len(jobs)} points x {len(maps_run.POSITIONS_FLUX)} positions through femm_pool ({workers} workers, task timeout 600 s), FEMM {femm_session.femm_version()}, "
            f"the maps' own model (maps_femm.flux_task: fine mesh, smart mesh on)")
    t0 = time.monotonic()
    with femm_queue.femm_pool("motor_map_run: spot-check", workers=workers, task_timeout_s=600) as pool:
        res = pool.map(maps_femm.flux_task, [j["task"] for j in jobs])
    rep.say(f"    done in {time.monotonic() - t0:.0f} s")
    fm = mm.load_maps()
    rows = []
    ok_all = True
    for p, r in zip(pts, res):
        if not (isinstance(r, dict) and r.get("status") == "ok"):
            rep.say(f"    FAILED {p['key']}: {r}")
            ok_all = False
            continue
        a = np.asarray(r["rows"])
        psi_d, psi_q, T = float(a[:, 3].mean()), float(a[:, 4].mean()), float(a[:, 5].mean())
        md, mq = fm.psi_dq(p["i_d"], p["i_q"])
        mT = float(fm.torque(p["i_d"], p["i_q"]))
        e_psi = max(abs(md - psi_d), abs(mq - psi_q)) / math.hypot(psi_d, psi_q)
        e_T = abs(mT - T)
        ok = e_psi <= TOL_F_PSI and e_T <= TOL_F_T_REL * abs(T) + TOL_F_T_ABS
        ok_all &= ok
        rows.append({**p, "femm_psi_d": psi_d, "femm_psi_q": psi_q, "femm_T": T, "maps_psi_d": float(md), "maps_psi_q": float(mq), "maps_T": mT, "err_psi": e_psi, "err_T": e_T, "ok": ok,
                     "nodes": int(r.get("nodes", 0)), "seconds": float(r.get("seconds", 0.0))})
        rep.say(f"    {p['key']} at {p['rpm']:.0f} rpm ({p['region']}): reference-domain current ({p['i_d']:.2f}, {p['i_q']:.2f}) A = ({p['i_d'] / fm.rated:.2f}, {p['i_q'] / fm.rated:.2f}) x rated: "
                f"psi_d {1e3 * psi_d:.3f} against {1e3 * md:.3f} mWb, psi_q {1e3 * psi_q:.3f} against {1e3 * mq:.3f} mWb, T {T:.4f} against {mT:.4f} N.m: psi {100 * e_psi:.3f} % of |psi|, T {e_T * 1e3:.1f} mN.m "
                f"({int(r.get('nodes', 0)):,} mesh nodes, {float(r.get('seconds', 0.0)):.0f} s)")
    rep.check("(f) FEMM at deep field-weakening points reproduces the maps", ok_all and len(rows) == len(jobs), f"{len(rows)} of {len(jobs)} points; limits: psi {100 * TOL_F_PSI:g} % of |psi|, "
              f"T {100 * TOL_F_T_REL:g} % + {TOL_F_T_ABS:g} N.m; FEMM 5.0 fine mesh (maps_run's), 15 positions of one ripple period")
    (out / "spot_check.json").write_text(json.dumps(rows, indent=1, default=float), encoding="utf-8")


def restore_inputs(inp: dict) -> None:
    """Put the inputs of a saved run back (the summaries, the spot-check and the figures of a saved run then print what that run used, whatever the options of this call)."""
    global WINDING, STACK_WINDINGS, STACKS, HEADLINE_STACK, VDCS, K_MARGIN, COPPER, ARRANGEMENT, MAGNET_LIMIT_C, EXTRAS, EXTRAS_ON, PEAK_FACTOR
    mk = lambda w: mm.Winding(w["turns"], int(w["paths"]), w["strands"], w["wire_mm"], w["wire_over_mm"])
    WINDING = mk(inp["winding"])
    STACK_WINDINGS = {float(k): mk(v) for k, v in inp["windings"]["opt"].items()}
    STACKS, HEADLINE_STACK, VDCS, K_MARGIN = tuple(inp["stacks"]), float(inp["headline_stack"]), tuple(inp["vdc"]), float(inp["k_margin"])
    COPPER, ARRANGEMENT, MAGNET_LIMIT_C, PEAK_FACTOR = inp["copper_kind"], inp["arrangement"], float(inp["magnet_limit_c"]), float(inp["peak_factor"])
    EXTRAS = tuple({"key": e["key"], "label": e["label"], "winding": mk(e["winding"]), "magnet_limit_c": float(e["magnet_limit_c"])} for e in inp["extras"])
    EXTRAS_ON = bool(inp["extras"])


def parse_winding(text: str) -> mm.Winding:
    n, a, ns, d, do = text.split(",")
    return mm.Winding(float(n), int(a), float(ns), float(d), float(do))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="the motor map: checks, envelopes, maps, stack table, figures")
    ap.add_argument("--label", default="final")
    ap.add_argument("--stages", default="checks,run,figures", help="comma list of: checks, run, report (print the summaries of a saved run), spot, figures")
    ap.add_argument("--quick", action="store_true", help="coarse grids everywhere (development)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--winding", default=None, help="the base winding, 'same' on every stack: turns,paths,strands,wire_mm,wire_over_mm (default: the winding study's headline winding 9,4,7,0.56,0.63)")
    ap.add_argument("--stack-windings", default=None, help="the winding of each stack for the 'opt' variant: 'stack:turns,paths,strands,wire,over;stack:...' "
                                                            "(default: the winding study's choice for 20, 40, 60 and 80 mm)")
    ap.add_argument("--copper", choices=("dc", "full"), default="full", help="the copper-loss model: copper_loss.py's full loss with the diamond end turn (default), or the first study's DC loss")
    ap.add_argument("--arrangement", default=None, help="how the strands of a turn lie in the slot (copper_loss.py): cluster (default), row, transposed, random, line")
    ap.add_argument("--magnet-limit", type=float, default=None, help="the headline magnet grade's temperature limit [C] (default 80: N45)")
    ap.add_argument("--no-extras", action="store_true", help="skip the front's other windings and the 150 C magnet grade")
    ap.add_argument("--vdc", default=None, help="bus voltages, e.g. 48,60")
    ap.add_argument("--k", type=float, default=None, help="the share of the modulation's linear limit used (default 0.95)")
    ap.add_argument("--stacks", default=None, help="stack lengths [mm], e.g. 20,40,60,80")
    ap.add_argument("--headline", type=float, default=None, help="the headline stack [mm] (default 60)")
    args = ap.parse_args(argv)
    global WINDING, STACK_WINDINGS, COPPER, ARRANGEMENT, MAGNET_LIMIT_C, EXTRAS_ON, VDCS, K_MARGIN, STACKS, HEADLINE_STACK
    if args.winding:
        WINDING = parse_winding(args.winding)
    if args.stack_windings:
        STACK_WINDINGS = {float(item.split(":")[0]): parse_winding(item.split(":")[1]) for item in args.stack_windings.split(";")}
    COPPER = args.copper
    if args.arrangement:
        ARRANGEMENT = args.arrangement
    if args.magnet_limit is not None:
        MAGNET_LIMIT_C = args.magnet_limit
    if args.no_extras:
        EXTRAS_ON = False
    if args.vdc:
        VDCS = tuple(float(v) for v in args.vdc.split(","))
    if args.k is not None:
        K_MARGIN = args.k
    if args.stacks:
        STACKS = tuple(float(v) for v in args.stacks.split(","))
    if args.headline is not None:
        HEADLINE_STACK = args.headline
    out = OUT_ROOT / args.label
    out.mkdir(parents=True, exist_ok=True)
    rep = Report()
    stages = set(args.stages.split(","))
    t0 = time.time()
    if "checks" in stages:
        rep.say("the checks that need no run")
        check_h_j(rep)
        check_a(rep)
        check_b(rep)
        check_d(rep)
        check_e(rep)
        check_m_k(rep)
        check_g(rep)
        check_q(rep)
    if "run" in stages:
        rep.say("the run")
        results = run_all(rep, args.label, args.quick)
        rep.say("the efficiency maps")
        run_maps(rep, results, args.quick)
        rep.say("the voltage margin and the swap of the copper model")
        run_extras(rep, results)
        rep.say("the thermal cross-check")
        results["thermal_reference"] = check_g2(rep, results)
        rep.say("the checks over every solution")
        check_points(rep, results)
        rows = stack_rows(results)
        write_csv(out / "stack_table.csv", rows)
        write_envelopes(out / "envelopes.csv", results)
        try:
            results["headline"] = headline(rep, results)
            print_stack_table(rep, results)
            print_speed_limits(rep, results)
            print_net_output(rep, results)
            print_extras(rep, results)
        finally:                                         # the cases took the time: they are saved even if a summary fails to print (the report stage prints them again)
            with open(out / "results.pkl", "wb") as f:
                pickle.dump(results, f)
            (out / "inputs.json").write_text(json.dumps(results["inputs"], indent=1), encoding="utf-8")
    if "report" in stages:
        with open(out / "results.pkl", "rb") as f:
            results = pickle.load(f)
        restore_inputs(results["inputs"])
        results["headline"] = headline(rep, results)
        print_stack_table(rep, results)
        print_speed_limits(rep, results)
        print_net_output(rep, results)
        print_extras(rep, results)
        for key, em in results.get("maps", {}).items():
            ok = em["minloss"]["found"]
            eta = em["minloss"]["eta"]
            best = np.unravel_index(np.nanargmax(np.where(ok, eta, np.nan)), eta.shape)
            good = ok & (em["T"] >= 0.1 * em["T_env"][:, None]) & (em["rpm"][:, None] >= 500.0) & np.isfinite(eta)
            rep.say(f"  map {key}: best efficiency {100 * eta[best]:.2f} % at {em['rpm'][best[0]]:,.0f} rpm and {em['T'][best]:.2f} N.m; {int((good & (eta > 0.9)).sum())} of {int(good.sum())} points above 90 %; "
                    f"minimum-loss control saves at most {np.nanmax(em['saving']):.2f} W against the MTPA control; {int((~ok & (em['T'] > 0)).sum())} points under the envelope unsolved")
    if "spot" in stages:
        spot_stage(rep, out, args.workers)
    if "figures" in stages:
        from optimization import motor_map_figures
        for path in motor_map_figures.make_all(out):
            rep.say(f"  figure {path}")
    rep.say(f"[{time.time() - t0:.0f} s]")
    name = "report.txt" if "checks" in stages else ("report_spot.txt" if "spot" in stages and "run" not in stages else ("report_summary.txt" if "report" in stages and "run" not in stages else "report_run.txt"))
    (out / name).write_text("\n".join(rep.lines) + "\n", encoding="utf-8")
    bad = [c for c in rep.checks if not c[1]]
    rep.say(f"{len(rep.checks) - len(bad)} of {len(rep.checks)} checks pass")
    return 0 if not bad else 1


if __name__ == "__main__":
    sys.exit(main())
