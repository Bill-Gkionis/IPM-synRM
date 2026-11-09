"""B6: Python -> ROSS (envs\\ross): the motor's rotor on its shaft and two ball bearings: its critical speeds, the Campbell diagram, the margins at 10,000 and 12,000 rpm.

What it shows, and why
    A rotor that spins at one of its natural frequencies shakes itself apart, so the speed range of the motor has to stay clear of them. ROSS is the library that finds them
    for a shaft with discs and bearings; Task 4 proved it on a bar whose answer is a formula. B6 sends the MOTOR's rotor through it: the 15 mm shaft, the 90 mm x 60 mm stack
    as a disc with its inertias (laminated M-19 and the magnets' mass), two ball bearings, and asks for the critical speeds and the Campbell diagram. It is sound when ROSS's
    first critical speed is within 10 % of a Jeffcott hand estimate made from the same numbers, and when the margins at 10,000 and 12,000 rpm are reported (there is no speed
    target yet; Steps holds the bridges at 12,000 rpm to start).

Method (the files and the protocol)
    1. The numbers (bench/b6_ross.py, P's side). Each comes with its source:
         shaft      C45E, 15 mm, E 212 GPa, nu 0.29, 7,850 kg/m3: 04 Structural 6.1 (the note's shaft steel). G = E / (2 (1 + nu)).
         stack      90 mm x 60 mm (the brief; the hub's "200 laminations of 0.3 mm make a 60 mm stack with 3 rows of magnets"). The steel is the ring between the 15 mm shaft and the
                    rim, less the 16 magnet slots (261.6 mm2, motorgeo's), at M-19's 7,650 kg/m3 times the stacking factor 0.98 (FEMM's lamination fill of the project's M-19). The factor is
                    the note's own rule, not a departure from its density: the table of 04 Structural 6.1 gives 7,650 kg/m3 "for a sheet" and, for the rotor and stator, "a stack: x its
                    stacking factor"; the magnets are N45, 16 of 9.08 x 1.57 mm in three rows of 19.14 mm (the hub), 7,600 kg/m3.
                    That is 2.664 kg of steel and 0.0995 kg of magnets. As a disc: Ip = m (ro^2 + ri^2) / 2 and Id = m (3 (ro^2 + ri^2) + L^2) / 12 of a hollow cylinder of that mass
                    (the slots take steel away at r = 33-43 mm and the magnets put almost the same mass back, so the error is below 1 %).
         span       105 mm between the bearings' centres: the 60 mm stack, and on each side 15 mm of end winding (the half chord of a 4-slot coil, 18 mm at the slots' mean radius,
                    times tan 30 degrees, plus 5 mm of straight: an estimate, nothing in the repository says), a 3 mm clearance and half of a 9 mm wide bearing (a 6002, 15 x 32 x 9 mm).
                    No overhang of the shaft beyond the bearings is modelled (no load, no coupling yet). The stack adds mass and inertia, not stiffness: a stack of laminations on a
                    shaft is far less stiff than a solid block (the note, 5), so the shaft's 15 mm carries the bending alone.
         bearings   k = 3.7e7 N/m each (37 N/um), no damping. From Hertz's sphere on a plane: each ball squeezes between two raceways, the approach of the two bodies is 2 delta,
                    delta = (9 Q^2 / (16 E*^2 R))^(1/3) (E* = E / (2 (1 - nu^2)), E = 200 GPa, nu 0.3 for the bearing steel, R the ball's radius); the ball's stiffness is
                    dQ/d(approach) = (3/2) Q / approach; Z balls all under the preload Q (no radial clearance) give a radial stiffness (Z/2) of that. The bearing is a 6002 (15 x 32 x 9 mm)
                    and its balls are the maker's: NSK's page for the 6002 (https://www.nsk.com/engineering/6002-apn.html, read on 2026-11-06) lists Z = 9 balls of Dw = 4.762 mm (3/16 inch),
                    and the load ratings Cr = 6,150 N and C0r = 2,840 N, which the calculation does not use. (The case first assumed 8 balls of 4.763 mm and 5.6-6.3 kN from catalogues it
                    did not name; the page replaced that on the review's request: the stiffness rose by 12 %.) Assumed: the preload, Q = 25 N each (an axial preload of about 45 N from a wave
                    spring, at a contact angle of about 12 degrees). The estimate is a LOWER bound: the grooves fit the ball better than a plane does (a factor of 1.3 to 1.6 in stiffness).
                    The preload is a guess; so the case asks what the first critical speed does when the stiffness, or the span, is off by a factor of 2 (the variants below).
    2. case.json (in) and result.json (out) in the case folder: the python of envs\\ross (checks.ross_check.ross_python, built from paths.P, never typed) runs
       bench/b6_ross_case.py under bench.run_program (the solver slot, below normal, a timeout), with -I (isolated mode: nothing of P's environment can leak into it).
       The script builds the ROSS rotor and writes its numbers; nothing is parsed from what it prints (ROSS prints a REFPROP notice when it is imported).
    3. What ROSS is asked: run_modal at rest (the natural frequencies, in pairs: x and y), run_critical_speed (every speed at which a natural frequency equals the speed, with its
       whirl; unbalance excites FORWARD whirl only), run_campbell over 0 to 36,000 rpm in 500 rpm steps (the damped natural frequencies of the first 8 modes). Six variants of the
       model (base, bearing stiffness x2 and x0.5, span x2 and x0.5, rigid bearings) give the first frequency and critical speed again.
    4. The Jeffcott hand estimate, from the same numbers: a disc at mid-span of a shaft of bending stiffness k_shaft = 48 E I / L^3 on two bearings that are springs in parallel under the
       load, 2 k_b, in series with the shaft: 1/k = 1/k_shaft + 1/(2 k_b); the moving mass is the disc and 17/35 of the shaft's own; f = sqrt(k/m) / 2 pi, and the critical speed is
       60 f rpm. The 17/35 = 0.486 is Rayleigh's method for a simply supported beam under a central load: take the static deflection as the shape, y ~ 3x - 4x^3 for 0 <= x <= 1/2
       (x in units of the span), and the beam's kinetic energy at the speed of its middle is that of a mass of 2 (the integral from 0 to 1/2 of (3x - 4x^3)^2 dx) = 17/35 of the beam's
       mass at the middle (a test does the integral). Blevins (Formulas for Natural Frequency and Mode Shape, Van Nostrand Reinhold, 1979, Table 8-8, case 4, p. 158: a pinned-pinned
       beam with a mass at mid-span) gives the same share rounded to 0.49. The table, the case and the page are as the Ansys LS-DYNA verification manual (VM-LSDYNA-SOLVE-046)
       quotes them; the book itself was not opened.
    5. The numbers B6 reports: ROSS's first forward critical speed against Jeffcott's (the brief's 10 %), the first critical speeds (forward and backward, in rpm), the Campbell
       diagram's crossing with the 1x line (an independent way to the same speed, from the arrays), the margins (critical / speed, and the distance in %) at the two speeds, the
       effect of each variant.
    6. The picture: ROSS's own Campbell figure is a plotly figure, and saving a plotly figure as a PNG needs the package kaleido, which envs\\ross does not have (and which must not
       be installed there). So ROSS's numbers (its Campbell arrays, passed back as JSON) are plotted with matplotlib in P's own environment, in the usual form: frequency against
       speed, forward branches solid and backward ones dashed, the 1x line, the two speeds, the critical speed marked.

Run (from the Python to FEMM folder; the solver slot is taken for about 15 s)
    uv run python bench/run_bench.py --only B6
Writes out\\bench\\B6\\: results.json, figure.png (the Campbell diagram), case.json and result.json (the protocol), ross.out.txt and ross.err.txt (what ROSS printed), run.log, runs.json.

Numbers to expect (this PC, 2026-11-06, ROSS 3.1.0; the whole case takes 20.6 s alone and 23.4 s inside the whole bench run, nearly all of it ROSS's start-up and its six variants)
    rotor       the stack 2.764 kg (steel 2.664, magnets 0.0995; Ip 2.88e-3 and Id 2.27e-3 kg m2) on 0.1457 kg of shaft: 2.910 kg, 138 degrees of freedom
    bearings    3.66e7 N/m each (Hertz: NSK's 9 balls of 4.762 mm at 25 N); the shaft alone, 48 E I / L^3, is 2.18e7 N/m, and the two bearings in series with it 1.68e7 N/m
    ROSS        the first natural frequency 380.57 Hz at rest (twice: x and y), the next pair 704.8 Hz (the tilting of the stack), the shaft's own bending 7,554 Hz; the first forward
                critical speed 22,838 rpm (the backward one 22,830, and a tilting mode at 28,140 that ROSS labels mixed, in all 16 solo runs of the repeatability trial and in the
                bench's own runs); the Campbell diagram's crossing with the 1x line 22,838 rpm, the same to 1e-8
    Jeffcott    23,265 rpm (387.8 Hz): ROSS is 1.84 % below it (not isolated: the likely reason is the shear deformation of ROSS's Timoshenko shaft, which Jeffcott's 48 E I / L^3 leaves
                out: it is about 4 % of the shaft's compliance here, so about 2 % of the frequency; the same offset, 2.1 %, with rigid bearings)
    margins     the first critical speed is 2.28 x 10,000 rpm (+128 %) and 1.90 x 12,000 rpm (+90 %)
    variants    bearing stiffness x2: 24,245 rpm (+6 %); x0.5: 20,627 (-10 %); rigid bearings 25,944 (+14 %); span x2: 9,038 rpm (-60 %: BELOW 10,000 rpm); span x0.5: 39,820 (+74 %; a span of
                52.5 mm is shorter than the 60 mm stack, so it cannot be built: the formula's number); against their Jeffcott estimates ROSS is 0.6 to 2.9 % below in every variant
    run to run  not exactly repeatable. In 16 solo runs of this step on 2026-11-06 (8 with THREADS = 2, this case's own setting, and 8 with threads=1, one after another through the
                slot) the result took two values in each group: ROSS's modes near 450,000 rpm (critical speeds of 446,165, 447,547 and 451,958 rpm) differ by up to 4.9e-8 (4.3e-8 with one
                thread) and the critical speeds that matter here (22,830 to 28,140 rpm) by 1.5e-11 at most. Of the Campbell diagram's whirl flags, 106 (128 with one thread) change from run
                to run: 46 (64) between 0.0 and 0.5, 57 (63) between 0.5 and 1.0 and 3 (1) between 0.0 and 1.0. The case's own `critical.whirl` label of the highest mode (451,958 rpm)
                changes between Forward and Mixed (Mixed in all 8 runs with THREADS = 2; Forward in 4 and Mixed in 4 with one thread); the labels of the other five modes (22,830 to
                447,547 rpm) never change. Setting OMP_NUM_THREADS=1 (run_solver's threads=1) did not remove the variation; the cause is not isolated. B9 compares B6's numbers between two
                runs to 1e-6 (B3's exactly) for that reason, and compares no labels
Self-check (the answers follow).
  1. The rotor is a disc in the middle of a shaft; why does ROSS give its first frequency twice?
  2. The bearings' stiffness is doubled. Does the critical speed double?
  3. The second pair of branches (the tilting of the stack, 705 Hz at rest) parts as the speed rises: the forward branch climbs, the backward one falls, and the 1x line crosses the
     backward one at 28,140 rpm. Is that crossing a critical speed to design around, like the one at 22,838 rpm?
Answers.
  1. The rotor is round and the bearings are the same in x and in y: it bends in the x direction exactly as in y, two modes with the same frequency. A keyway, or bearings of
     different stiffness in x and y, would split the pair (and the pair's forward and backward whirl split with the speed, by the gyroscopic effect of the disc).
  2. No: the shaft and the two bearings are springs in series, 1/k = 1/k_shaft + 1/(2 k_b). Doubling k_b makes the bearing part of the compliance half, not the shaft's. With the shaft
     about as stiff as the bearings the critical speed rises by only 6 % (22,838 to 24,245 rpm); rigid bearings are the limit, +14 %, where only the shaft's compliance is left.
  3. Not for unbalance: unbalance is a force that turns with the rotor, in the direction of the rotation, so it excites the whirl that goes round the same way (frequency = speed, the
     1x line). The 1x line crosses the backward branch too, but a force that goes forward cannot drive a mode that goes backward. The first pair does not part in the diagram: the
     stack moves sideways in it without tilting, so the disc's gyroscopic effect hardly touches it (22,830 backward and 22,838 forward).
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

import motorgeo as mg
import paths
from checks import Missing
from checks import ross_check

from . import OWN_LIMIT, Log, check, compare, params, run_program, wrap_title
from .b4_ccx import M19, N45

CASE = "B6"
TITLE = "Python -> ROSS: the motor's rotor on its shaft and two ball bearings, critical speeds and Campbell diagram"
VERSION = "3.1.0"                       # the pin of envs\ross\pyproject.toml
TOLERANCE_PERCENT = 10.0                # the brief: the first critical speed within 10 % of a Jeffcott hand estimate from the same numbers
SPEEDS_RPM = (10_000, 12_000)           # the brief: margins at these two (there is no speed target yet; Steps holds the bridges at 12,000 rpm to start)

# ---- the shaft: C45E (04 Structural 6.1)
SHAFT_DIAMETER_MM = 15.0
SHAFT_E = 212e9                         # Pa
SHAFT_NU = 0.29
SHAFT_RHO = 7850.0                      # kg/m3
# ---- the stack: 90 mm x 60 mm (the brief)
STACK_OD_MM = 90.0
STACK_LENGTH_MM = 60.0
STACKING_FACTOR = mg.M19_STEEL_FEMM42["LamFill"]      # the share of the stack that is steel: FEMM's lamination fill of the project's M-19 (0.98)
RHO_STEEL = M19.density_kg_m3           # 7,650 kg/m3 as a sheet (the note, 6.4)
RHO_MAGNET = N45.density_kg_m3          # 7,600 kg/m3
MAGNET_ROWS, MAGNET_ROW_LENGTH_MM = 3, 19.14           # the hub: a 60 mm stack with 3 rows of magnets, each 19.14 mm long
# ---- the layout between the bearings (an estimate: nothing in the repository gives it)
END_WINDING_MM, CLEARANCE_MM, BEARING_WIDTH_MM = 15.0, 3.0, 9.0       # the stator's end winding beyond the core, the clearance to the bearing, a 6002's width
# ---- the bearings: a 6002 as Hertz sees it (the balls are the maker's, the preload is assumed: see the module docstring)
BALLS, BALL_DIAMETER_MM = 9, 4.762      # NSK's page for the 6002 (https://www.nsk.com/engineering/6002-apn.html, read 2026-11-06): z = 9, Dw = 4.762 mm
PRELOAD_PER_BALL_N = 25.0               # assumed: an axial preload of about 45 N at a contact angle of about 12 degrees (nothing of the maker's)
BEARING_E_MPA, BEARING_NU = 200e3, 0.3
RIGID_STIFFNESS = 1e11                  # N/m: a wall (Task 4's check uses the same)
# The share of the shaft's mass that moves with the disc: Rayleigh with the static deflection of a central load, 2 int_0^(1/2) (3x - 4x^3)^2 dx = 17/35 (Blevins 1979, Table 8-8
# case 4, rounds it to 0.49, as the Ansys verification manual VM-LSDYNA-SOLVE-046 quotes it; the book was not opened)
JEFFCOTT_SHAFT_SHARE = 17 / 35
PLOT_CUT_HZ = 1500                      # the Campbell picture shows the branches under this frequency; the others lie above its top
ELEMENTS = 22                           # even, so that the disc sits on a node at mid-span
CAMPBELL_MAX_RPM, CAMPBELL_POINTS, CAMPBELL_MODES = 36_000, 73, 8
CRITICAL_MODES = 12
TIMEOUT_S = 300
THREADS = 2
CAMPBELL_AGREE_PERCENT = 1.0            # the critical speed of ROSS's solver against the 1x crossing of its own Campbell arrays: the bench's own limit
PAIR_TOLERANCE = 1e-3                   # the two modes of a pair have the same frequency within this share (isotropic rotor)
MASS_TOLERANCE_PERCENT = 0.1            # the rotor ROSS built has the mass that was asked for
MARGIN_RULE_PERCENT = 25.0              # a rule of thumb for a motor: its first critical speed at least 25 % above its top speed (not in the brief; reported, not gating)


# ============================================================================================== the numbers
def shear_modulus() -> float:
    """G = E / (2 (1 + nu)) of the shaft [Pa]."""
    return SHAFT_E / (2 * (1 + SHAFT_NU))


def bearing_span_mm(stack_mm: float, end_winding_mm: float, clearance_mm: float, bearing_width_mm: float) -> float:
    """The distance between the bearings' centres: the stack and, on each side, the end winding, the clearance and half a bearing."""
    return stack_mm + 2 * (end_winding_mm + clearance_mm + bearing_width_mm / 2)


SPAN_MM = bearing_span_mm(STACK_LENGTH_MM, END_WINDING_MM, CLEARANCE_MM, BEARING_WIDTH_MM)


def disk_inertias(mass_kg: float, inner_radius_m: float, outer_radius_m: float, length_m: float) -> tuple[float, float]:
    """(Ip, Id) [kg m2] of a hollow cylinder: Ip about the axis, m (ro^2 + ri^2) / 2; Id about a diameter through the centre, m (3 (ro^2 + ri^2) + L^2) / 12."""
    squares = outer_radius_m ** 2 + inner_radius_m ** 2
    return 0.5 * mass_kg * squares, mass_kg / 12 * (3 * squares + length_m ** 2)


def rotor_stack(p: mg.Params) -> dict:
    """The rotor's stack as a disc: the steel's mass (the ring between the shaft and the rim less the magnet slots of the geometry of `p`, at the laminated density), the magnets'
    mass, the total, the inertias, and the slot area. {"slots", "slots_area_mm2", "steel_mass_kg", "magnet_mass_kg", "mass_kg", "Ip_kgm2", "Id_kgm2"}."""
    regions = mg.build(p, gap="sliding")
    slot_regions = [r for r in regions if r.name.startswith("magnet_")]
    slots = len(slot_regions) * p.n_sectors
    slots_area = sum(r.polygon().area for r in slot_regions) * p.n_sectors                         # mm2
    ring_area = math.pi / 4 * (STACK_OD_MM ** 2 - (2 * p.r_shaft) ** 2)                            # mm2
    steel = (ring_area - slots_area) * 1e-6 * STACK_LENGTH_MM * 1e-3 * RHO_STEEL * STACKING_FACTOR
    magnets = slots * (p.mag_w * 1e-3) * (p.mag_t * 1e-3) * (MAGNET_ROWS * MAGNET_ROW_LENGTH_MM * 1e-3) * RHO_MAGNET
    mass = steel + magnets
    ip, id_ = disk_inertias(mass, p.r_shaft * 1e-3, STACK_OD_MM / 2 * 1e-3, STACK_LENGTH_MM * 1e-3)
    return {"slots": slots, "slots_area_mm2": slots_area, "steel_mass_kg": steel, "magnet_mass_kg": magnets, "mass_kg": mass, "Ip_kgm2": ip, "Id_kgm2": id_}


def ball_bearing_radial_stiffness(balls: int, ball_diameter_mm: float, preload_per_ball_n: float, youngs_mpa: float = BEARING_E_MPA, poisson: float = BEARING_NU) -> float:
    """The radial stiffness [N/m] of a ball bearing whose `balls` balls all carry the same preload: Hertz's sphere on a plane at both raceways, (Z/2) x (3/2) Q / (2 delta).
    A lower bound: the grooves fit the ball better than a plane. It grows as the cube root of the preload and of the ball's diameter, and in proportion to the number of balls."""
    e_star = youngs_mpa / (2 * (1 - poisson ** 2))                                                  # MPa, two identical steels
    delta_mm = (9 * preload_per_ball_n ** 2 / (16 * e_star ** 2 * (ball_diameter_mm / 2))) ** (1 / 3)
    return balls / 2 * 1.5 * preload_per_ball_n / (2 * delta_mm * 1e-3)


BEARING_STIFFNESS = ball_bearing_radial_stiffness(BALLS, BALL_DIAMETER_MM, PRELOAD_PER_BALL_N)


# ============================================================================================== the hand estimate (Jeffcott)
def shaft_bending_stiffness(youngs_pa: float, diameter_m: float, span_m: float) -> float:
    """k = 48 E I / L^3 [N/m]: the stiffness of a round shaft simply supported at both ends, loaded at the middle."""
    return 48 * youngs_pa * (math.pi * diameter_m ** 4 / 64) / span_m ** 3


def series_stiffness(shaft_stiffness: float, bearing_stiffness: float) -> float:
    """The shaft in series with the two bearings (parallel springs under the central load, 2 k_b): 1/k = 1/k_shaft + 1/(2 k_b)."""
    return 1 / (1 / shaft_stiffness + 1 / (2 * bearing_stiffness))


def jeffcott(disk_mass: float, shaft_mass: float, shaft_stiffness: float, bearing_stiffness: float) -> dict:
    """Jeffcott's rotor: a disc of `disk_mass` at mid-span of a massless-but-for-17/35 shaft; f = sqrt(k / m) / 2 pi with k the shaft in series with the bearings.
    {"stiffness_n_per_m", "mass_kg", "frequency_hz", "critical_rpm"}."""
    k = series_stiffness(shaft_stiffness, bearing_stiffness)
    m = disk_mass + JEFFCOTT_SHAFT_SHARE * shaft_mass
    f = math.sqrt(k / m) / (2 * math.pi)
    return {"stiffness_n_per_m": k, "mass_kg": m, "frequency_hz": f, "critical_rpm": 60 * f}


def hand_estimate(p: mg.Params, bearing_stiffness: float | None = None, span_mm: float | None = None) -> dict:
    """The Jeffcott estimate of the motor's rotor from the numbers of the case (the disc's mass, the shaft, the bearing and the span, any of the last two replaced for a variant)."""
    span = (SPAN_MM if span_mm is None else span_mm) * 1e-3
    k_b = BEARING_STIFFNESS if bearing_stiffness is None else bearing_stiffness
    stack = rotor_stack(p)
    diameter = SHAFT_DIAMETER_MM * 1e-3
    shaft_mass = SHAFT_RHO * math.pi / 4 * diameter ** 2 * span
    k_shaft = shaft_bending_stiffness(SHAFT_E, diameter, span)
    return {"span_mm": span * 1e3, "bearing_stiffness_n_per_m": k_b, "shaft_stiffness_n_per_m": k_shaft, "disk_mass_kg": stack["mass_kg"], "shaft_mass_kg": shaft_mass,
            **jeffcott(stack["mass_kg"], shaft_mass, k_shaft, k_b)}


def jeffcott_check(hand_rpm: float, ross_rpm: float, **kw) -> dict:
    """The brief's: ROSS's first critical speed within 10 % of the Jeffcott estimate of the same rotor."""
    return compare("the first critical speed against Jeffcott's estimate from the same numbers", ross_rpm, hand_rpm, TOLERANCE_PERCENT,
                   "the brief: a Jeffcott hand estimate (the shaft's bending stiffness in series with the bearings, the disc and 17/35 of the shaft's mass)", "rpm", **kw)


# ============================================================================================== the case file, and what comes back
def case_inputs(p: mg.Params) -> dict:
    """The content of case.json: everything ROSS's environment needs, in SI units, with the variants (the base, the bearing's stiffness x2 and x0.5, the span x2 and x0.5, rigid bearings)."""
    stack = rotor_stack(p)
    span = SPAN_MM * 1e-3
    k = BEARING_STIFFNESS
    variants = [{"name": "base", "k": k, "length": span}, {"name": "bearing stiffness x2", "k": 2 * k, "length": span}, {"name": "bearing stiffness x0.5", "k": k / 2, "length": span},
                {"name": "span x2", "k": k, "length": 2 * span}, {"name": "span x0.5", "k": k, "length": span / 2}, {"name": "rigid bearings", "k": RIGID_STIFFNESS, "length": span}]
    return {"shaft": {"diameter": SHAFT_DIAMETER_MM * 1e-3, "length": span, "elements": ELEMENTS, "E": SHAFT_E, "rho": SHAFT_RHO, "G_s": shear_modulus()},
            "disk": {"node": ELEMENTS // 2, "m": stack["mass_kg"], "Id": stack["Id_kgm2"], "Ip": stack["Ip_kgm2"]},
            "bearing": {"k": k, "c": 0.0},
            "campbell": {"max_rpm": CAMPBELL_MAX_RPM, "points": CAMPBELL_POINTS, "frequencies": CAMPBELL_MODES},
            "critical": {"num_modes": CRITICAL_MODES}, "speeds_rpm": list(SPEEDS_RPM), "variants": variants}


REQUIRED = ("ross_version", "rotor_mass_kg", "ndof", "wn_hz_at_rest", "critical", "campbell", "variants")


def read_result(path) -> dict:
    """result.json as a dict. ValueError, saying what is wrong, if the file is not there, is not JSON, or lacks a key ROSS's side is to give."""
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"{path} is not there: ROSS's side did not write its answer")
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"{path} is not JSON: {error}") from error
    missing = [key for key in REQUIRED if key not in data]
    if missing:
        raise ValueError(f"{path} lacks {', '.join(missing)}")
    return data


def rpm_of(rad_s: float) -> float:
    return rad_s * 30 / math.pi


def all_critical_speeds_rpm(result: dict) -> list[tuple[float, str]]:
    """Every critical speed ROSS found, in rpm, with its whirl ("Forward" or "Backward"), the lowest first."""
    found = result["critical"]
    return sorted((rpm_of(w), whirl) for w, whirl in zip(found["wn_rad_s"], found["whirl"]))


def forward_critical_speeds_rpm(result: dict) -> list[float]:
    """The critical speeds of FORWARD whirl in rpm, the lowest first: the ones a synchronous unbalance can reach."""
    return [rpm for rpm, whirl in all_critical_speeds_rpm(result) if whirl == "Forward"]


def first_forward_critical_rpm(result: dict) -> float:
    """The lowest critical speed of forward whirl [rpm]; ValueError if ROSS found none."""
    forward = forward_critical_speeds_rpm(result)
    if not forward:
        raise ValueError("ROSS's answer has no forward critical speed (it found only backward ones)")
    return forward[0]


def campbell_crossings(result: dict, order: float = 1.0) -> list[dict]:
    """Where a forward branch of the Campbell diagram meets the line frequency = order x speed (the 1x line for synchronous unbalance), found by linear interpolation between the
    diagram's points: [{"mode", "rpm", "frequency_hz", "order"}], by speed. A branch is a column of the diagram; it counts as forward where ROSS's whirl value is 0."""
    camp = result["campbell"]
    speed = np.asarray(camp["speed_rad_s"], float) * 30 / math.pi                                  # rpm
    frequency = np.asarray(camp["wd_rad_s"], float) / (2 * math.pi)                                # Hz
    whirl = np.asarray(camp["whirl_values"], float)
    gap = frequency - order * speed[:, None] / 60                                                  # branch minus the line
    out = []
    for j in range(frequency.shape[1]):
        for i in range(len(speed) - 1):
            if whirl[i, j] == 0 and whirl[i + 1, j] == 0 and gap[i, j] * gap[i + 1, j] < 0:
                t = gap[i, j] / (gap[i, j] - gap[i + 1, j])
                out.append({"mode": j, "rpm": float(speed[i] + t * (speed[i + 1] - speed[i])), "frequency_hz": float(frequency[i, j] + t * (frequency[i + 1, j] - frequency[i, j])),
                            "order": order})
    return sorted(out, key=lambda c: c["rpm"])


def margins(critical_rpm: float, speeds=SPEEDS_RPM) -> dict:
    """{speed: {"ratio": critical / speed, "margin_percent": (critical - speed) / speed in %}}: negative when the critical speed is below the speed, and never hidden."""
    return {s: {"ratio": critical_rpm / s, "margin_percent": 100 * (critical_rpm - s) / s} for s in speeds}


# ============================================================================================== the picture
def branch_style(whirl) -> str:
    """The line style of one Campbell branch from the whirl ROSS gives along it (0 forward, 1 backward, 0.5 'mixed': it cannot say): solid when the points where it says forward outnumber
    those where it says backward, dashed for the reverse, dotted when it never says or the two tie. The mixed points do not vote: ROSS gives 0.5 at a share of the tilting pair's points
    that changes from one run to the next (a median over all the points drew the pair's upper branch solid in one run of the case and dotted in the next)."""
    values = np.asarray(whirl, float)
    forward, backward = int(np.count_nonzero(values == 0.0)), int(np.count_nonzero(values == 1.0))
    return "-" if forward > backward else "--" if backward > forward else ":"


def branch_split(result: dict) -> dict:
    """The branches of ROSS's Campbell diagram (the columns of its array) that the picture shows, those whose frequency stays under PLOT_CUT_HZ at every speed, and those above its top:
    {"total", "shown", "above", "above_min_hz", "above_max_hz"} (the last two are the lowest and the highest frequency of the branches above over all the speeds; None if there are none)."""
    frequency = np.asarray(result["campbell"]["wd_rad_s"], float) / (2 * math.pi)
    shown = (frequency < PLOT_CUT_HZ).all(axis=0)
    above = frequency[:, ~shown]
    return {"total": int(frequency.shape[1]), "shown": int(shown.sum()), "above": int((~shown).sum()),
            "above_min_hz": float(above.min()) if above.size else None, "above_max_hz": float(above.max()) if above.size else None}


def campbell_title(result: dict, first_critical_rpm: float, margins_by_speed: dict, disk_mass_kg: float, stiffness: float, span_mm: float) -> str:
    """The title of figure.png: it says what the plot shows (the branches under PLOT_CUT_HZ, how many of the first CAMPBELL_MODES, and where the others lie), the rotor, and the first
    forward critical speed with its margins. Every line is at most bench.TITLE_CHARS characters."""
    split = branch_split(result)
    lines = [f"ROSS {result['ross_version']}: Campbell diagram of the motor's rotor, no damping: {split['shown']} of the first {split['total']} branches"]
    if split["above"]:
        lines.append(f"(the other {split['above']} lie between {split['above_min_hz']:,.0f} and {split['above_max_hz']:,.0f} Hz, above the plot)")
    lines += [f"{disk_mass_kg:.2f} kg stack on a 15 mm shaft; two bearings of {stiffness:.2g} N/m, {span_mm:.0f} mm apart",
              f"first forward critical speed {first_critical_rpm:,.0f} rpm: {margins_by_speed[10_000]['ratio']:.2f} x 10,000 rpm, {margins_by_speed[12_000]['ratio']:.2f} x 12,000 rpm",
              "matplotlib plot of ROSS's arrays (its plotly figure needs kaleido, which envs\\ross lacks)"]
    return "\n".join(lines)


def fig_campbell(result: dict, critical_rpm: float, crossing: dict | None, target, title: str) -> None:
    """The Campbell diagram from ROSS's arrays, in matplotlib: frequency [Hz] against speed [rpm], forward branches solid and backward ones dashed, the 1x line, the two speeds,
    the first forward critical speed marked. The branches above PLOT_CUT_HZ are off the plot (the title says how many, and where)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    camp = result["campbell"]
    speed = np.asarray(camp["speed_rad_s"], float) * 30 / math.pi
    frequency = np.asarray(camp["wd_rad_s"], float) / (2 * math.pi)
    whirl = np.asarray(camp["whirl_values"], float)
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    styles = [branch_style(whirl[:, j]) for j in range(frequency.shape[1])]       # a column is one branch; ROSS's whirl value of it is 0 (forward), 1 (backward) or 0.5 (it cannot say) at each speed
    for j, style in enumerate(styles):
        ax.plot(speed, frequency[:, j], style, color=f"C{j // 2}", lw=1.7)
    top = float(np.nanmax(np.where(frequency < PLOT_CUT_HZ, frequency, np.nan))) * 1.08
    ax.plot(speed, speed / 60, "-.", color="#222222", lw=1.3, label="1x (synchronous unbalance)")
    for s in SPEEDS_RPM:
        ax.axvline(s, color="#c0392b", lw=1.2, alpha=0.8)
        ax.text(s, top * 0.97, f" {s:,} rpm", color="#c0392b", rotation=90, va="top", ha="left", fontsize=9)
    ax.plot([], [], "-", color="#555555", label="forward whirl")
    ax.plot([], [], "--", color="#555555", label="backward whirl")
    if ":" in styles:
        ax.plot([], [], ":", color="#555555", label="whirl ROSS does not class")
    if crossing is not None:
        ax.plot(crossing["rpm"], crossing["frequency_hz"], "o", ms=9, mfc="none", mec="#1f1f1f", mew=1.8)
        ax.annotate(f"first critical speed\n{critical_rpm:,.0f} rpm", (crossing["rpm"], crossing["frequency_hz"]), xytext=(15, -45), textcoords="offset points", fontsize=9.5,
                    arrowprops={"arrowstyle": "-", "color": "#444444"})
    ax.set(xlim=(0, float(speed.max())), ylim=(0, top), xlabel="speed [rpm]", ylabel="natural frequency [Hz]")
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9)             # the corner the branches and the 1x line leave empty (at the upper left it covered the label of the 10,000 rpm line)
    ax.set_title(wrap_title(title), fontsize=10.5, loc="left", color="#1d1d1f")
    fig.tight_layout()
    fig.savefig(target, dpi=130)
    plt.close(fig)


# ============================================================================================== the case
def run(folder: Path, log: Log) -> dict:
    """Run B6 in `folder` (already empty). Returns the results dict whose "checks" list every verdict."""
    folder = Path(folder)
    started = time.perf_counter()
    p = mg.Params(barrier="none")
    python = ross_check.ross_python()
    if not python.is_file():
        raise Missing("ross", python)
    case = case_inputs(p)
    (folder / "case.json").write_text(json.dumps(case, indent=1), encoding="utf-8")
    hands = {v["name"]: hand_estimate(p, bearing_stiffness=v["k"], span_mm=v["length"] * 1e3) for v in case["variants"]}
    script = paths.P / "bench" / "b6_ross_case.py"
    # -I: isolated mode, so nothing of the caller's environment (PYTHONPATH, P's own packages) can leak into ROSS's python
    ross = run_program(log, [python, "-I", script, folder], label="bench B6 ROSS", cwd=folder, timeout_s=TIMEOUT_S, threads=THREADS, name="ross")
    seconds = round(time.perf_counter() - started, 2)
    answered = not ross.timed_out and ross.returncode == 0 and (folder / "result.json").is_file()
    results = {"case": CASE, "title": TITLE, "params": params(p), "inputs": case, "stack": rotor_stack(p), "hand_estimate": hands, "bearing": {
        "stiffness_n_per_m": BEARING_STIFFNESS, "balls": BALLS, "ball_diameter_mm": BALL_DIAMETER_MM, "preload_per_ball_n": PRELOAD_PER_BALL_N, "youngs_mpa": BEARING_E_MPA,
        "poisson": BEARING_NU, "model": "Hertz sphere on a plane at both raceways, Z/2 balls effective; a lower bound; Z and D_w are NSK's for the 6002, the preload is assumed"},
        "layout": {"span_mm": SPAN_MM, "stack_mm": STACK_LENGTH_MM, "end_winding_mm": END_WINDING_MM, "clearance_mm": CLEARANCE_MM, "bearing_width_mm": BEARING_WIDTH_MM},
        "ross": {"returncode": ross.returncode, "timed_out": ross.timed_out, "seconds": seconds}}
    if not answered:
        results["checks"] = [check("ROSS runs in its own environment and answers", False, "envs\\ross's python runs bench/b6_ross_case.py and writes result.json",
                                   f"exit code {ross.returncode}, timed out {ross.timed_out}, result.json there: {(folder / 'result.json').is_file()}; "
                                   f"{' | '.join(ross.stderr.strip().splitlines()[-3:])}", "bench/b6_ross_case.py through run_program")]
        results["seconds_total"] = round(time.perf_counter() - started, 1)
        return results
    result = read_result(folder / "result.json")
    analysis = analyse(result, hands, p)
    results["result"] = {k: v for k, v in result.items() if k != "campbell"}
    results["analysis"] = analysis
    fig_campbell(result, analysis["first_forward_critical_rpm"], analysis["campbell_first_crossing"], folder / "figure.png",
                 campbell_title(result, analysis["first_forward_critical_rpm"], analysis["margins"], analysis["disk_mass_kg"], BEARING_STIFFNESS, SPAN_MM))
    log.say(f"B6 ROSS {result['ross_version']}: rotor mass {result['rotor_mass_kg']:.3f} kg, {result['ndof']} dof; f1 at rest {result['wn_hz_at_rest'][0]:.2f} Hz; first forward critical speed "
            f"{analysis['first_forward_critical_rpm']:,.0f} rpm (Jeffcott {hands['base']['critical_rpm']:,.0f}, {analysis['jeffcott_difference_percent']:+.2f} %); Campbell 1x crossing "
            f"{(analysis['campbell_first_crossing'] or {}).get('rpm', float('nan')):,.0f} rpm; margins {', '.join(f'{s:,}: x{m['ratio']:.2f}' for s, m in analysis['margins'].items())}; {seconds:.1f} s")
    results["checks"] = make_checks(result, analysis, hands, ross, folder / "figure.png")
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    return results


def analyse(result: dict, hands: dict, p: mg.Params) -> dict:
    """The numbers of B6 from ROSS's answer: the first forward critical speed and the Campbell diagram's crossing with the 1x line, against the hand estimate; the margins; the
    variants against their own hand estimates."""
    first = first_forward_critical_rpm(result)
    crossings = campbell_crossings(result, 1.0)
    stack = rotor_stack(p)
    variants = []
    for v in result["variants"]:
        hand = hands[v["name"]]
        variants.append({**v, "hand_critical_rpm": hand["critical_rpm"], "hand_f_hz": hand["frequency_hz"],
                         "difference_percent": None if v["critical_forward_rpm"] is None else 100 * (v["critical_forward_rpm"] - hand["critical_rpm"]) / hand["critical_rpm"]})
    return {"first_forward_critical_rpm": first, "forward_critical_speeds_rpm": forward_critical_speeds_rpm(result), "all_critical_speeds_rpm": all_critical_speeds_rpm(result),
            "campbell_crossings_1x": crossings, "campbell_first_crossing": crossings[0] if crossings else None, "margins": margins(first),
            "jeffcott_difference_percent": 100 * (first - hands["base"]["critical_rpm"]) / hands["base"]["critical_rpm"], "disk_mass_kg": stack["mass_kg"],
            "mass_expected_kg": stack["mass_kg"] + hands["base"]["shaft_mass_kg"], "variants": variants,
            "f1_pair_ratio": result["wn_hz_at_rest"][1] / result["wn_hz_at_rest"][0]}


def make_checks(result: dict, analysis: dict, hands: dict, ross, figure) -> list[dict]:
    """Every number of B6 against its expectation, with the source of the expectation. `ross` is the SolverRun of ROSS's side, `figure` the path of the picture."""
    out: list[dict] = []
    out.append(check("ROSS runs in its own environment and answers", result["ross_version"] == VERSION and ross.returncode == 0 and not ross.timed_out,
                     f"envs\\ross's python (never P's) runs bench/b6_ross_case.py through run_program and writes result.json; ROSS {VERSION} (the pin of envs\\ross\\pyproject.toml)",
                     f"exit code {ross.returncode}, {ross.seconds:.1f} s, ROSS {result['ross_version']}, "
                     f"{result['ndof']} degrees of freedom, {result['elements']} shaft elements, every variant solved: {all(v['error'] is None for v in result['variants'])}",
                     "bench/b6_ross_case.py, result.json"))
    out.append(compare("the rotor ROSS built has the mass that was asked for", result["rotor_mass_kg"], analysis["mass_expected_kg"], MASS_TOLERANCE_PERCENT,
                       f"the disc ({analysis['disk_mass_kg']:.3f} kg: the stack's steel and magnets) and the shaft (7,850 kg/m3 x its volume), the case's own sum ({OWN_LIMIT})", "kg"))
    out.append(jeffcott_check(hands["base"]["critical_rpm"], analysis["first_forward_critical_rpm"]))
    out.append(check("the first natural frequency comes twice (x and y)", abs(analysis["f1_pair_ratio"] - 1) <= PAIR_TOLERANCE,
                     f"an isotropic rotor bends in x as in y: the first two natural frequencies at rest are equal within {PAIR_TOLERANCE:g} of themselves ({OWN_LIMIT})",
                     f"{result['wn_hz_at_rest'][0]:.3f} Hz and {result['wn_hz_at_rest'][1]:.3f} Hz (ratio {analysis['f1_pair_ratio']:.6f}); the next pair {result['wn_hz_at_rest'][2]:.1f} Hz, "
                     f"the bending of the shaft in the second mode {result['wn_hz_at_rest'][4]:.0f} Hz", "ROSS run_modal at rest; the same as Task 4's shaft row"))
    cross = analysis["campbell_first_crossing"]
    difference = None if cross is None else 100 * (cross["rpm"] - analysis["first_forward_critical_rpm"]) / analysis["first_forward_critical_rpm"]
    out.append(check("the Campbell diagram and the critical-speed solver agree", cross is not None and abs(difference) <= CAMPBELL_AGREE_PERCENT,
                     f"the first crossing of a forward branch with the 1x line, found by interpolation in ROSS's Campbell arrays, is the first forward critical speed of ROSS's own solver "
                     f"within {CAMPBELL_AGREE_PERCENT} % ({OWN_LIMIT}): two ways to the same speed, which also tests the whirl convention of the arrays (0 forward)",
                     "no forward branch crosses the 1x line" if cross is None else f"{cross['rpm']:,.0f} rpm (mode {cross['mode']}, {cross['frequency_hz']:.1f} Hz) against "
                     f"{analysis['first_forward_critical_rpm']:,.0f} rpm ({difference:+.3f} %); all the 1x crossings: " + ", ".join(f"{c['rpm']:,.0f}" for c in analysis["campbell_crossings_1x"]) + " rpm",
                     "ross.Rotor.run_campbell and run_critical_speed, the same rotor"))
    m = analysis["margins"]
    out.append(check("the first critical speed clears the speeds with a margin (a rule of thumb)", all(v["margin_percent"] >= MARGIN_RULE_PERCENT for v in m.values()),
                     f"reported at 10,000 and 12,000 rpm: the first forward critical speed at least {MARGIN_RULE_PERCENT:g} % above each ({OWN_LIMIT}: a rule of thumb for a motor; the brief only "
                     f"asks that the margins are reported, and no speed target exists yet)",
                     f"first forward critical speed {analysis['first_forward_critical_rpm']:,.0f} rpm: " + "; ".join(f"at {s:,} rpm the ratio is {v['ratio']:.2f} (a margin of {v['margin_percent']:+.1f} %)" for s, v in m.items())
                     + "; all critical speeds (rpm, whirl): " + ", ".join(f"{r:,.0f} {w[0]}" for r, w in analysis["all_critical_speeds_rpm"] if r < 100_000),
                     "ROSS run_critical_speed", gating=False, reason="the brief asks for the margins to be reported, not to meet a number"))
    base = next(v for v in analysis["variants"] if v["name"] == "base")
    by = {v["name"]: v for v in analysis["variants"]}
    ordered = all(by[a]["critical_forward_rpm"] is not None for a in by)
    sane = ordered and by["bearing stiffness x0.5"]["critical_forward_rpm"] < base["critical_forward_rpm"] < by["bearing stiffness x2"]["critical_forward_rpm"] < by["rigid bearings"]["critical_forward_rpm"] \
        and by["span x2"]["critical_forward_rpm"] < base["critical_forward_rpm"] < by["span x0.5"]["critical_forward_rpm"]
    effect = lambda name: f"{by[name]['critical_forward_rpm']:,.0f} rpm ({100 * (by[name]['critical_forward_rpm'] / base['critical_forward_rpm'] - 1):+.0f} %, Jeffcott {by[name]['hand_critical_rpm']:,.0f})" \
        if by[name]["critical_forward_rpm"] is not None else f"no answer ({by[name]['error']})"
    out.append(check("the critical speed moves the way the physics says when the bearing or the span is changed", sane,
                     "stiffer bearings raise it (x0.5 < base < x2 < rigid), a longer span lowers it (span x2 < base < span x0.5): the shaft and the bearings are springs in series, the shaft's "
                     "stiffness goes as 1/L^3",
                     f"base {base['critical_forward_rpm']:,.0f} rpm; bearing stiffness x2 {effect('bearing stiffness x2')}, x0.5 {effect('bearing stiffness x0.5')}; span x2 {effect('span x2')}, "
                     f"x0.5 {effect('span x0.5')} (52.5 mm is shorter than the 60 mm stack: it cannot be built, it is the formula's number); rigid bearings {effect('rigid bearings')}",
                     "the six variants of the same ROSS model"))
    saved = Path(figure).is_file() and Path(figure).stat().st_size > 20_000
    out.append(check("the picture: figure.png is saved", saved, "out\\bench\\B6\\figure.png is the Campbell diagram, plotted with matplotlib from ROSS's arrays",
                      f"saved: {saved}" + (f", {Path(figure).stat().st_size:,} bytes" if saved else ""), "bench/b6_ross.py"))
    return out
