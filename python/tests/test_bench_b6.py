"""Tests for B6 (bench/b6_ross.py and bench/b6_ross_case.py): Python -> ROSS. Everything that needs no ROSS: the rotor's numbers, the bearing, the Jeffcott hand estimate, the
JSON that ROSS's environment hands back, the critical speeds, the Campbell diagram's crossings and the margins.

Run from P:   uv run pytest tests/test_bench_b6.py -q      (a second; ROSS lives in its own environment and is not started, no solver slot is taken)

In groups:
  1. The inputs: the shaft and the stack with their sources, the rotor's mass and inertias, the bearing (a 6002: the maker's balls, an assumed preload) and the span, the brief's 10 % pin.
  2. The hand estimate: Jeffcott's disc on a shaft whose bending stiffness is in series with the bearings, from the same numbers; the shaft's 17/35 is Rayleigh's integral.
  3. What ROSS hands back: the case file in, the result file out, the critical speeds (forward ones only), the Campbell diagram's crossing with the 1x line, the margins, the line style
     of a branch and the title that says which branches the plot shows.
"""
import ast
import json
import math

import numpy as np
import pytest

import motorgeo as mg
from bench_test_tools import P


# ----------------------------------------------------------------------------------------------- 1. the inputs
def test_the_briefs_ten_percent_and_the_two_speeds_are_pinned():
    """The brief: the first critical speed within 10 % of a Jeffcott hand estimate; margins at 10,000 and 12,000 rpm. Change a constant and this goes red."""
    from bench import b6_ross as b6
    assert b6.TOLERANCE_PERCENT == 10.0 and b6.SPEEDS_RPM == (10_000, 12_000) and b6.VERSION == "3.1.0"
    item = b6.jeffcott_check(22_000.0, 22_000.0 * 1.0999)
    assert item["verdict"] == "PASS" and item["tolerance_percent"] == 10.0 and item["gating"] is True
    assert b6.jeffcott_check(22_000.0, 22_000.0 * 1.0999 / 1.0)["verdict"] == "PASS"
    assert b6.jeffcott_check(22_000.0, 22_000.0 * 1.101)["verdict"] == "FAIL" and b6.jeffcott_check(22_000.0, 22_000.0 * 0.899)["verdict"] == "FAIL"


def test_the_shaft_is_the_notes_c45e_and_the_stack_is_the_notes_m19_and_n45():
    """04 Structural 6.1: the shaft is C45E, E 212 GPa, nu 0.29, 7,850 kg/m3, 15 mm; the stack's steel is M-19, 7,650 kg/m3 as a sheet (6.4), times the stacking factor (the note, 4.1:
    a stack weighs a little less than its sheets); the magnets are N45, 7,600 kg/m3. The stack is 90 mm x 60 mm (the brief)."""
    from bench import b4_ccx as b4
    from bench import b6_ross as b6
    assert (b6.SHAFT_E, b6.SHAFT_NU, b6.SHAFT_RHO, b6.SHAFT_DIAMETER_MM) == (212e9, 0.29, 7850.0, 15.0)
    assert b6.shear_modulus() == pytest.approx(212e9 / (2 * 1.29))
    assert (b6.STACK_OD_MM, b6.STACK_LENGTH_MM) == (90.0, 60.0)
    assert b6.STACKING_FACTOR == mg.M19_STEEL_FEMM42["LamFill"] == 0.98                    # FEMM's lamination fill of the project's M-19: the stack's share of steel
    assert b6.RHO_STEEL == b4.M19.density_kg_m3 and b6.RHO_MAGNET == b4.N45.density_kg_m3


def test_the_rotor_stack_is_the_steel_ring_minus_the_slots_plus_the_magnets():
    """The disc of the model: the steel between the 15 mm shaft and the 90 mm rim over 60 mm, less the 16 slots (261.6 mm2 in all, motorgeo's), times the laminated density, plus the
    magnets (16 of 9.08 x 1.57 mm in three rows of 19.14 mm: the hub's 'magnet 1.57 x 9.08 x 19.14 mm' and 'a 60 mm stack with 3 rows of magnets')."""
    from bench import b6_ross as b6
    s = b6.rotor_stack(mg.Params(barrier="none"))
    ring = math.pi / 4 * (0.090 ** 2 - 0.015 ** 2)
    slots = 16 * 16.348342e-6
    assert s["slots_area_mm2"] == pytest.approx(261.57, abs=0.05)
    assert s["steel_mass_kg"] == pytest.approx((ring - slots) * 0.060 * 7650.0 * 0.98, rel=1e-4) and s["steel_mass_kg"] == pytest.approx(2.664, abs=0.005)
    assert s["magnet_mass_kg"] == pytest.approx(16 * 9.08e-3 * 1.57e-3 * 3 * 19.14e-3 * 7600.0, rel=1e-9) and s["magnet_mass_kg"] == pytest.approx(0.0995, abs=0.0005)
    assert s["mass_kg"] == pytest.approx(s["steel_mass_kg"] + s["magnet_mass_kg"]) and s["mass_kg"] == pytest.approx(2.764, abs=0.006)
    assert s["Ip_kgm2"] == pytest.approx(0.5 * s["mass_kg"] * (0.045 ** 2 + 0.0075 ** 2)) and s["Id_kgm2"] == pytest.approx(s["mass_kg"] / 12 * (3 * (0.045 ** 2 + 0.0075 ** 2) + 0.060 ** 2))


def test_the_inertias_of_a_hollow_cylinder_are_the_textbook_ones():
    """Ip (about the axis) = m (ro^2 + ri^2) / 2; Id (about a diameter through the centre) = m (3 (ro^2 + ri^2) + L^2) / 12. A solid disc (ri = 0) has m r^2 / 2 and m (3 r^2 + L^2) / 12;
    a thin ring has Ip = m r^2 and Id = m r^2 / 2."""
    from bench import b6_ross as b6
    assert b6.disk_inertias(2.0, 0.0, 0.1, 0.02) == pytest.approx((2.0 * 0.01 / 2, 2.0 * (3 * 0.01 + 0.0004) / 12))
    ip, id_ = b6.disk_inertias(1.0, 0.0999999, 0.1, 1e-9)
    assert ip == pytest.approx(0.01, rel=1e-6) and id_ == pytest.approx(0.005, rel=1e-6)
    ip, id_ = b6.disk_inertias(3.0, 0.0075, 0.045, 0.060)
    assert ip > id_ * 0.5 and id_ > 0 and ip == pytest.approx(3.0 * (0.045 ** 2 + 0.0075 ** 2) / 2)


def test_the_bearing_stiffness_is_hertz_s_sphere_on_a_plane_twice_per_ball_and_the_balls_in_parallel():
    """Each ball squeezes between two raceways, so the approach of the two bodies is 2 x delta, delta = (9 Q^2 / (16 E*^2 R))^(1/3) (Hertz, sphere on a plane, E* = E / (2 (1 - nu^2)),
    R the ball's radius). The ball's stiffness is dQ/d(approach) = (3/2) Q / approach; Z balls all under the preload Q (no radial clearance) give a radial stiffness (Z/2) x that.
    It grows as the cube root of the preload and of the ball's diameter. The estimate is a lower bound: the grooves fit the ball better than a plane does."""
    from bench import b6_ross as b6
    z, d_w, q, e, nu = 8, 4.763, 25.0, 200e3, 0.3
    e_star = e / (2 * (1 - nu ** 2))
    delta = (9 * q ** 2 / (16 * e_star ** 2 * (d_w / 2))) ** (1 / 3)                      # mm, one contact
    k = z / 2 * 1.5 * q / (2 * delta * 1e-3)
    got = b6.ball_bearing_radial_stiffness(z, d_w, q, e, nu)
    assert got == pytest.approx(k, rel=1e-12) and got == pytest.approx(3.25e7, rel=0.02)
    assert b6.ball_bearing_radial_stiffness(z, d_w, 8 * q, e, nu) == pytest.approx(2 * got, rel=1e-9)               # 8 times the preload: twice the stiffness (cube root)
    assert b6.ball_bearing_radial_stiffness(2 * z, d_w, q, e, nu) == pytest.approx(2 * got, rel=1e-12)              # twice the balls: twice the stiffness
    assert b6.ball_bearing_radial_stiffness(z, 8 * d_w, q, e, nu) == pytest.approx(2 * got, rel=1e-9)               # eight times the diameter: twice


def test_the_bearing_is_a_6002_with_the_makers_balls_and_an_assumed_preload():
    """Fix round 1, item 9: NSK's page for the 6002 (https://www.nsk.com/engineering/6002-apn.html, read 2026-11-06) lists z = 9 balls of Dw = 4.762 mm; the case had assumed 8 of 4.763
    mm. The preload (25 N a ball) is the one assumed number, and the case still asks what a factor of 2 in the stiffness does."""
    from bench import b6_ross as b6
    assert (b6.BALLS, b6.BALL_DIAMETER_MM, b6.PRELOAD_PER_BALL_N) == (9, 4.762, 25.0)
    assert b6.BEARING_STIFFNESS == pytest.approx(b6.ball_bearing_radial_stiffness(9, 4.762, 25.0), rel=1e-12) and b6.BEARING_STIFFNESS == pytest.approx(3.662e7, rel=2e-3)
    assert b6.BEARING_STIFFNESS / b6.ball_bearing_radial_stiffness(8, 4.763, 25.0) == pytest.approx(1.125, abs=0.001)                  # the first assumption was 12 % softer
    source = (P / "bench" / "b6_ross.py").read_text(encoding="utf-8")
    assert "https://www.nsk.com/engineering/6002-apn.html" in source and "Cr = 6,150 N" in source and "C0r = 2,840 N" in source and "Z = 9 balls of Dw = 4.762 mm" in source


def test_the_bearing_span_is_the_stack_plus_the_end_winding_the_clearance_and_a_bearing_width():
    """The layout: the stack, then the stator's end winding (15 mm beyond the core), a 3 mm clearance and half a bearing (a 6002 is 9 mm wide) on each side:
    60 + 2 (15 + 3 + 4.5) = 105 mm between the bearings' centres."""
    from bench import b6_ross as b6
    assert b6.bearing_span_mm(60.0, 15.0, 3.0, 9.0) == pytest.approx(105.0)
    assert (b6.END_WINDING_MM, b6.CLEARANCE_MM, b6.BEARING_WIDTH_MM) == (15.0, 3.0, 9.0) and b6.SPAN_MM == pytest.approx(105.0)
    assert b6.SPAN_MM > b6.STACK_LENGTH_MM                                                  # the stack fits between the bearings


# ----------------------------------------------------------------------------------------------- 2. the hand estimate
def test_the_jeffcott_estimate_puts_the_shafts_bending_stiffness_in_series_with_the_two_bearings():
    """The shaft carries the disc at mid-span: k_shaft = 48 E I / L^3. The bearings, each of stiffness k_b, are two springs in parallel under the load, 2 k_b; the shaft and them are in
    series: 1/k = 1/k_shaft + 1/(2 k_b). The moving mass is the disc and 17/35 of the shaft's (the shaft's own share, for the shape the central load bends it into). f = sqrt(k/m) / 2 pi."""
    from bench import b6_ross as b6
    e, d, span = 212e9, 0.015, 0.105
    inertia = math.pi * d ** 4 / 64
    assert b6.shaft_bending_stiffness(e, d, span) == pytest.approx(48 * e * inertia / span ** 3, rel=1e-12) and b6.shaft_bending_stiffness(e, d, span) == pytest.approx(2.18e7, rel=0.01)
    assert b6.series_stiffness(2.0e7, 3.0e7) == pytest.approx(1 / (1 / 2.0e7 + 1 / 6.0e7))
    assert b6.series_stiffness(2.0e7, 1e15) == pytest.approx(2.0e7, rel=1e-6)           # a rigid bearing: the shaft alone
    h = b6.jeffcott(disk_mass=2.764, shaft_mass=0.1457, shaft_stiffness=2.184e7, bearing_stiffness=3.25e7)
    k = 1 / (1 / 2.184e7 + 1 / 6.5e7)
    m = 2.764 + 17 / 35 * 0.1457
    assert h["stiffness_n_per_m"] == pytest.approx(k, rel=1e-12) and h["mass_kg"] == pytest.approx(m, rel=1e-12)
    assert h["frequency_hz"] == pytest.approx(math.sqrt(k / m) / (2 * math.pi), rel=1e-12) and h["critical_rpm"] == pytest.approx(60 * h["frequency_hz"], rel=1e-12)


def test_the_shaft_s_share_of_the_moving_mass_is_rayleighs_17_over_35():
    """Fix round 1, item 9c: the beam's share is 2 x the integral from 0 to 1/2 of (3x - 4x^3)^2 dx, the static deflection of a simply supported beam under a central load (x in units
    of the span, the middle's deflection scaled to 1): 17/35 = 0.486 (Blevins, Table 8-8 case 4, rounds it to 0.49, as the Ansys LS-DYNA verification manual VM-LSDYNA-SOLVE-046 quotes it; the book was not opened)."""
    from numpy.polynomial import Polynomial
    from bench import b6_ross as b6
    shape = Polynomial([0.0, 3.0, 0.0, -4.0])
    assert shape(0.5) == pytest.approx(1.0)                                                      # the middle's deflection is the scale
    squared = (shape * shape).integ()
    assert 2 * (squared(0.5) - squared(0.0)) == pytest.approx(17 / 35, rel=1e-12) == pytest.approx(b6.JEFFCOTT_SHAFT_SHARE, rel=1e-12)
    assert round(b6.JEFFCOTT_SHAFT_SHARE, 2) == 0.49


def test_the_jeffcott_critical_speed_falls_with_a_longer_span_and_rises_with_stiffer_bearings():
    """A factor of 2 in the span is a factor of 8 in the shaft's stiffness (L^3), a factor of 2 in a bearing's stiffness is 2 in that spring: the first moves the critical speed
    more. These are the hand numbers that the sensitivity table of the case is read against."""
    from bench import b6_ross as b6
    base = b6.hand_estimate(mg.Params(barrier="none"))
    assert base["span_mm"] == pytest.approx(105.0) and 340 < base["frequency_hz"] < 420
    stiffer = b6.hand_estimate(mg.Params(barrier="none"), bearing_stiffness=2 * base["bearing_stiffness_n_per_m"])
    softer = b6.hand_estimate(mg.Params(barrier="none"), bearing_stiffness=base["bearing_stiffness_n_per_m"] / 2)
    longer = b6.hand_estimate(mg.Params(barrier="none"), span_mm=210.0)
    assert softer["frequency_hz"] < base["frequency_hz"] < stiffer["frequency_hz"]
    assert longer["frequency_hz"] < 0.45 * base["frequency_hz"]


# ----------------------------------------------------------------------------------------------- 3. what ROSS hands back
def ross_result(speeds_rpm=(0, 10_000, 20_000, 30_000)):
    """A result.json in miniature: one isotropic rotor whose first pair splits with the speed (the forward branch rises, the backward one falls), ROSS's
    own critical speeds (rad/s, with the whirl each belongs to), and the Campbell arrays (speed_rad_s, wd_rad_s (speeds x modes), whirl_values, 0 forward, 1 backward)."""
    speed = [s * math.pi / 30 for s in speeds_rpm]
    wd = [[2 * math.pi * f for f in (370.0 - 0.001 * s / 100, 370.0 + 0.001 * s / 100, 700.0 - 0.01 * s / 100, 700.0 + 0.04 * s / 100)] for s in speeds_rpm]
    return {"ross_version": "3.1.0", "rotor_mass_kg": 2.91, "ndof": 92, "elements": 22,
            "wn_hz_at_rest": [370.0, 370.0, 700.0, 700.0, 7000.0, 7000.0],
            "critical": {"wn_rad_s": [2300.0, 2365.0, 2900.0], "whirl": ["Backward", "Forward", "Backward"]},
            "campbell": {"speed_rad_s": speed, "wd_rad_s": wd, "whirl_values": [[1.0, 0.0, 1.0, 0.0]] * len(speeds_rpm)},
            "variants": [{"name": "base", "f1_hz": 370.0, "critical_forward_rpm": 22_000.0}]}


def test_the_critical_speeds_are_the_forward_ones_in_rpm_and_the_first_is_the_lowest():
    """A synchronous unbalance excites forward whirl only: ROSS also lists the backward crossings (the isotropic rotor's backward mode crosses the 1x line too), which no
    unbalance reaches. ROSS gives rad/s; the case speaks in rpm."""
    from bench import b6_ross as b6
    r = ross_result()
    forward = b6.forward_critical_speeds_rpm(r)
    assert forward == pytest.approx([2365.0 * 30 / math.pi]) and b6.first_forward_critical_rpm(r) == pytest.approx(2365.0 * 30 / math.pi)
    assert b6.all_critical_speeds_rpm(r) == [(pytest.approx(2300.0 * 30 / math.pi), "Backward"), (pytest.approx(2365.0 * 30 / math.pi), "Forward"), (pytest.approx(2900.0 * 30 / math.pi), "Backward")]
    r["critical"] = {"wn_rad_s": [2900.0, 2300.0], "whirl": ["Backward", "Backward"]}
    with pytest.raises(ValueError, match="no forward"):
        b6.first_forward_critical_rpm(r)


def test_a_missing_or_wrong_result_is_a_value_error_that_says_what_is_missing(tmp_path):
    from bench import b6_ross as b6
    good = tmp_path / "result.json"
    good.write_text(json.dumps(ross_result()), encoding="utf-8")
    assert b6.read_result(good)["ross_version"] == "3.1.0"
    broken = ross_result()
    del broken["critical"]
    (tmp_path / "broken.json").write_text(json.dumps(broken), encoding="utf-8")
    with pytest.raises(ValueError, match="critical"):
        b6.read_result(tmp_path / "broken.json")
    with pytest.raises(ValueError, match="not there"):
        b6.read_result(tmp_path / "nope.json")


def test_the_campbell_diagram_crosses_the_one_times_line_where_a_forward_branch_meets_the_speed():
    """The forward branch of the first pair (the one that rises with the speed) meets the line frequency = speed at a speed that the interpolation between the diagram's points gives.
    Here: a branch f(N) = 370 Hz + 0.001 N / 100 meets N / 60: N = 370 / (1/60 - 1e-5) = 22,168 rpm (to the interpolation's own error)."""
    from bench import b6_ross as b6
    speeds = np.arange(0, 40_001, 500)
    wd = np.array([[2 * math.pi * (370.0 + 1e-5 * s), 2 * math.pi * (370.0 - 1e-5 * s), 2 * math.pi * (700.0 + 4e-4 * s)] for s in speeds])
    campbell = {"speed_rad_s": (speeds * math.pi / 30).tolist(), "wd_rad_s": wd.tolist(), "whirl_values": [[0.0, 1.0, 0.0]] * len(speeds)}
    crossings = b6.campbell_crossings({"campbell": campbell}, order=1)
    assert len(crossings) == 1 and crossings[0]["rpm"] == pytest.approx(370.0 / (1 / 60 - 1e-5), rel=2e-3) and crossings[0]["mode"] == 0
    twice = b6.campbell_crossings({"campbell": campbell}, order=2)                                   # the 2x line, N / 30: both rising branches meet it inside the range
    assert [c["mode"] for c in twice] == [0, 2]
    assert twice[0]["rpm"] == pytest.approx(370.0 / (2 / 60 - 1e-5), rel=2e-3) and twice[1]["rpm"] == pytest.approx(700.0 / (2 / 60 - 4e-4), rel=2e-3)


def test_a_branch_is_drawn_by_the_whirl_ross_gives_it_where_it_gives_one_and_the_mixed_points_do_not_vote():
    """ROSS gives 'mixed' (0.5) at a changing share of the points of the tilting pair (the first run of the case had 37 of 73 points forward on the upper branch and drew it solid, the
    next had 33 and, with a median over all the points, drew it dotted). Only the points where ROSS says forward (0) or backward (1) vote; none or a tie is dotted."""
    from bench import b6_ross as b6
    assert b6.branch_style([0.0] * 73) == "-" and b6.branch_style([1.0] * 73) == "--"
    assert b6.branch_style([0.0] * 33 + [0.5] * 40) == "-"                      # the second run's upper branch of the tilting pair: 33 forward, 40 mixed (the median said mixed)
    assert b6.branch_style([1.0] * 37 + [0.5] * 36) == "--"
    assert b6.branch_style([0.0] * 14 + [0.5] * 59) == "-" and b6.branch_style([1.0] * 72 + [0.0]) == "--"
    assert b6.branch_style([0.5] * 10) == ":" and b6.branch_style([0.0, 1.0]) == ":" and b6.branch_style([]) == ":"
    assert b6.branch_style(np.array([0.0, 0.5, 1.0, 1.0])) == "--"


def campbell_arrays():
    """A Campbell diagram in miniature: 8 branches over 5 speeds, the first four under 1,500 Hz (two pairs) and the last four near 7.3 kHz (the shaft's own bending)."""
    speeds = np.linspace(0.0, 3770.0, 5)
    low = [[2 * math.pi * f] * 5 for f in (375.0, 375.0, 695.0, 695.0)]
    high = [[2 * math.pi * (7270.0 + 6 * j + 2 * i) for i in range(5)] for j in range(4)]
    wd = np.array(low + high).T
    return {"ross_version": "3.1.0", "campbell": {"speed_rad_s": speeds.tolist(), "wd_rad_s": wd.tolist(), "whirl_values": [[0.0] * 8] * 5}}


def test_the_campbell_title_says_which_branches_the_plot_shows_and_where_the_others_are():
    """Fix round 1, item 9d: the title said 'the first 8 modes' over a plot of 4 branches. It counts the branches under the plot's top and says where the rest lie."""
    import bench
    from bench import b6_ross as b6
    result = campbell_arrays()
    split = b6.branch_split(result)
    assert (split["total"], split["shown"], split["above"]) == (8, 4, 4) and split["above_min_hz"] == pytest.approx(7270.0) and split["above_max_hz"] == pytest.approx(7296.0)
    title = b6.campbell_title(result, 23_000.0, b6.margins(23_000.0), 2.764, 3.66e7, 105.0)
    lines = title.split("\n")
    assert lines[0] == "ROSS 3.1.0: Campbell diagram of the motor's rotor, no damping: 4 of the first 8 branches"
    assert lines[1] == "(the other 4 lie between 7,270 and 7,296 Hz, above the plot)"
    assert "2.76 kg stack on a 15 mm shaft; two bearings of 3.7e+07 N/m, 105 mm apart" in lines and "2.30 x 10,000 rpm, 1.92 x 12,000 rpm" in lines[3] and "23,000 rpm" in lines[3]
    assert all(len(line) <= bench.TITLE_CHARS for line in lines) and bench.wrap_title(title) == title
    result["campbell"]["wd_rad_s"] = [row[:4] for row in result["campbell"]["wd_rad_s"]]
    only = b6.campbell_title(result, 23_000.0, b6.margins(23_000.0), 2.764, 3.66e7, 105.0).split("\n")
    assert only[0].endswith("4 of the first 4 branches") and not any("above the plot" in line for line in only)               # nothing above: no sentence about it


def test_margins_are_the_distance_of_the_critical_speed_from_each_speed():
    """margin = (critical - speed) / speed in percent, and the ratio critical / speed. At 10,000 and 12,000 rpm for a critical speed of 22,500 rpm: 125 % and 87.5 %."""
    from bench import b6_ross as b6
    got = b6.margins(22_500.0, (10_000, 12_000))
    assert got[10_000] == pytest.approx({"ratio": 2.25, "margin_percent": 125.0})
    assert got[12_000] == pytest.approx({"ratio": 1.875, "margin_percent": 87.5})
    below = b6.margins(9_000.0, (10_000,))
    assert below[10_000]["margin_percent"] == pytest.approx(-10.0)                                  # a critical speed below the running speed is a negative margin, never hidden


def test_the_case_file_carries_every_number_ross_needs_and_the_sensitivity_variants():
    """case.json is what ROSS's environment reads: the shaft, the disc (node, mass, Id, Ip), the bearing, the Campbell sweep and the variants (the bearing's stiffness x2 and x0.5, the
    span x2 and x0.5, rigid bearings). The disc sits at the middle node of an even number of elements."""
    from bench import b6_ross as b6
    c = b6.case_inputs(mg.Params(barrier="none"))
    assert c["shaft"]["elements"] % 2 == 0 and c["disk"]["node"] == c["shaft"]["elements"] // 2
    assert c["shaft"]["length"] == pytest.approx(0.105) and c["shaft"]["diameter"] == pytest.approx(0.015)
    assert c["disk"]["m"] == pytest.approx(2.764, abs=0.006) and c["disk"]["Id"] > 0 and c["disk"]["Ip"] > c["disk"]["Id"]
    assert c["bearing"]["k"] == pytest.approx(3.66e7, rel=0.01) and c["bearing"]["c"] == 0.0                # 9 balls of 4.762 mm at 25 N each (NSK's page for the 6002)
    names = [v["name"] for v in c["variants"]]
    assert names == ["base", "bearing stiffness x2", "bearing stiffness x0.5", "span x2", "span x0.5", "rigid bearings"]
    base, double, half, wide, narrow, rigid = c["variants"]
    assert (base["k"], base["length"]) == (c["bearing"]["k"], c["shaft"]["length"]) and double["k"] == 2 * base["k"] and half["k"] == base["k"] / 2
    assert wide["length"] == 2 * base["length"] and narrow["length"] == base["length"] / 2 and rigid["k"] >= 1e11
    json.dumps(c)                                                                                    # it can be written
    assert c["campbell"]["max_rpm"] >= 2 * max(b6.SPEEDS_RPM) and c["speeds_rpm"] == [10_000, 12_000]


def test_the_ross_side_script_runs_in_ross_s_environment_and_imports_nothing_of_the_projects():
    """bench/b6_ross_case.py is started by the python of envs\\ross, which has none of P's modules (R3 of the tools plan): it may import the standard library and ROSS only,
    and ROSS only inside main() (so that this environment can import the file for a look). It never starts a program."""
    source = (P / "bench" / "b6_ross_case.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    imported = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    imported |= {(n.module or "").split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert imported <= {"json", "math", "sys", "pathlib", "ross", "numpy", "__future__"}, imported
    top = {a.name for n in tree.body if isinstance(n, ast.Import) for a in n.names}
    assert "ross" not in top and "numpy" not in top, "ROSS is imported inside main(), after the case file is read"
    import importlib.util
    spec = importlib.util.spec_from_file_location("b6_ross_case_for_test", P / "bench" / "b6_ross_case.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)                                                                  # importing it needs no ROSS
    assert callable(module.main)


def test_the_case_is_registered_with_run_bench_after_b5():
    from bench import run_bench
    assert list(run_bench.CASES)[5] == "B6" and run_bench.CASES["B6"] == "b6_ross"          # after B5 (Task 13 added B7 to B9 after it)


def test_the_documentation_says_where_the_17_over_35_was_read_how_long_the_case_takes_and_what_ross_labels_the_tilting_mode():
    """T12-P1, T12-P2 and T13-P5 of the final review, wording that the data and the sources decide. (1) Blevins' table, case and page are second-hand: quoted by the Ansys
    LS-DYNA verification manual VM-LSDYNA-SOLVE-046; the book was not opened, and the SolidWorks manual was not checked, so it is not cited. (2) The case takes 20.6 s alone (the
    first fix round's 23 s was a busier run). (3) The tilting mode at 28,140 rpm is labelled Mixed by ROSS in all 16 solo runs of the repeatability trial and in the bench's own
    runs; only the highest mode (about 452,000 rpm) changes its label."""
    from bench import b6_ross
    doc = b6_ross.__doc__
    assert "VM-LSDYNA-SOLVE-046" in doc and "the book itself was not opened" in doc and "SolidWorks" not in doc and "SolidWorks" not in (P / "bench" / "b6_ross.py").read_text(encoding="utf-8")
    assert "the whole case takes 23 s" not in doc and "20.6 s" in doc
    assert "labels backward in some runs and mixed in others" not in doc and "that ROSS labels mixed" in doc
