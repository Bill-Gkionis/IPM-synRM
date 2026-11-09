"""Tests for optimization/thermal_fem.py: the exact solutions the solvers are held to, the CalculiX deck writer and reader, the slice's segments, and
the rule that turns the slice's numbers into a verdict. No solver and no FEMM are started.

Run from P:   uv run pytest tests/test_thermal_fem.py -q      (a few seconds)

The point of testing the exact solutions here is that a calibration is only as good as its reference: if composite_exact were wrong, a solver that agreed with
it would be "calibrated" onto a mistake. So the reference is checked against hand formulas first.
"""
import math

import numpy as np
import pytest
from shapely.geometry import MultiPolygon, Polygon

from optimization import thermal, thermal_fem as tf


def test_cylinder_exact_is_the_textbook_case():
    c = tf.cylinder_exact(radius=0.040, k=1.0, q=1.0e5, h=20.0, t_inf=0.0)
    assert c["T_surface"] == pytest.approx(100.0)                         # q R / (2 h)
    assert c["T_centre"] == pytest.approx(140.0)                          # + q R^2 / (4 k)
    assert c["q_total"] == pytest.approx(1.0e5 * math.pi * 0.040 ** 2)


def test_composite_exact_matches_hand_formulas():
    layers = tf.COMPOSITE
    ex = tf.composite_exact(layers, tf.COMPOSITE_H, tf.COMPOSITE_TINF)
    T = ex["T"]
    # the surface: all the heat made leaves through the outer radius
    q_total = sum(L.q * math.pi * (L.r_out ** 2 - L.r_in ** 2) for L in layers)
    assert ex["q_total"] == pytest.approx(q_total, rel=1e-12)
    assert T(layers[-1].r_out) == pytest.approx(tf.COMPOSITE_TINF + q_total / (2 * math.pi * layers[-1].r_out * tf.COMPOSITE_H), rel=1e-12)
    # a layer with no heat source passes a constant heat flow, so its drop is Q' ln(r2/r1) / (2 pi k) (the liner and the contact layer)
    for name in ("liner", "contact"):
        L = next(x for x in layers if x.name == name)
        q_at = sum(x.q * math.pi * (min(L.r_in, x.r_out) ** 2 - x.r_in ** 2) for x in layers if L.r_in > x.r_in)
        assert T(L.r_in) - T(L.r_out) == pytest.approx(q_at * math.log(L.r_out / L.r_in) / (2 * math.pi * L.k), rel=1e-8)
    # the hole's wall is the hottest point, and the temperature falls outward
    radii = np.linspace(layers[0].r_in, layers[-1].r_out, 40)
    values = [T(float(r)) for r in radii]
    assert all(a >= b for a, b in zip(values[:-1], values[1:]))


def test_ccx_deck_structure(tmp_path):
    meta = tf.write_ccx_annulus(tmp_path / "t.inp", tf.COMPOSITE, [3, 1, 4, 1, 2], 15.0, 0.0)
    text = (tmp_path / "t.inp").read_text()
    assert meta["n_elements"] == 11
    assert "*ELEMENT, TYPE=CPS8" in text and "*HEAT TRANSFER, STEADY STATE" in text
    elements = [line for line in text.split("*ELEMENT")[1].split("*ELSET")[0].splitlines()[1:] if line.strip()]
    assert len(elements) == 11 and all(len(line.split(",")) == 9 for line in elements)       # an id and eight nodes
    film = text.split("*FILM")[1].split("*NODE FILE")[0].strip().splitlines()
    assert len(film) == 1 and film[0].startswith("11, F2,")                                  # the outermost element's outer edge only
    assert "0.015" in film[0]                                                                  # 15 W/(m2 K) = 0.015 mW/(mm2 K)
    assert text.count("BF") == 2                                                              # winding and iron have heat sources
    assert "L0, BF, 0.2" in text                                                              # 2e5 W/m3 = 0.2 mW/mm3


def test_ccx_numbers_are_never_longer_than_twenty_characters():
    for value in (6.627128124470999e-09, 1.0 / 3.0, 2.2250738585072014e-308, 123456789.123456789):
        assert len(tf.ccx_number(value)) <= 19


def test_read_frd_temperature(tmp_path):
    frd = tmp_path / "a.frd"
    frd.write_text(" -4  NDTEMP      1    1\n -5  T           1    1    0    0\n -1         1 1.40000E+02\n -1         2 1.00000E+02\n -3\n -4  DISP        4    1\n"
                   " -1         1 9.99999E+09\n -3\n", encoding="ascii")
    assert tf.read_frd_temperature(frd) == {1: 140.0, 2: 100.0}


def test_slice_segments_close_every_region():
    p = thermal.design_params()
    task, _ = tf.build_sector_task(tf.SECTOR_VARIANTS["base"], p, 1.0, tf.OUT(), "test", split_teeth=True)
    regions = {}
    for name, d in task["regions"].items():
        parts = [Polygon(r) for r in d["parts"]]
        regions[name] = parts[0] if len(parts) == 1 else MultiPolygon(parts)
    nodes, segs = tf.noded_segments(regions)
    degree = {}
    for i, j in segs:
        degree[i] = degree.get(i, 0) + 1
        degree[j] = degree.get(j, 0) + 1
    assert min(degree.values()) >= 2                                      # no dangling edge: every region is a closed loop
    # the heat made in the slice, from the regions' own sources, is what the network is given
    areas = {k: v.area * 1e-6 for k, v in regions.items()}
    made = sum(task["regions"][k]["q"] * areas[k] for k in regions)
    assert made == pytest.approx(20.0 + 1.19e5 * (areas["teeth_shoe"] + areas["teeth_fin"] + areas["teeth_root"]) + 7.5e4 * areas["yoke"], rel=1e-9)


def _row(variant, hot_drop, mean_drop, hot_amb=0.1, mean_amb=0.1, h=60.0, copper=True, teeth_drop=0.0, yoke_drop=0.0, teeth_amb=0.0, yoke_amb=0.0):
    return {"variant": variant, "h_outer": h, "has_copper": copper, "T_hot_drop_err_%": hot_drop, "T_mean_drop_err_%": mean_drop, "T_hot_amb_err_%": hot_amb,
            "T_mean_amb_err_%": mean_amb, "T_teeth_drop_err_%": teeth_drop, "T_yoke_drop_err_%": yoke_drop, "T_teeth_amb_err_%": teeth_amb, "T_yoke_amb_err_%": yoke_amb}


def test_the_verdict_follows_the_stated_tolerances():
    rows = [_row("ok", 1.0, 2.0),
            _row("hot spot too far", 10.5, 2.0),                                     # 10.5 % > 10 % internal
            _row("ambient too far", 1.0, 1.0, hot_amb=5.5),                          # 5.5 % > 5 % above the ambient
            _row("liquid: only the internal criterion", 1.0, 1.0, hot_amb=6.0, mean_amb=6.0, h=2000.0),
            _row("iron only, tooth off by 18 %", 0.0, 0.0, copper=False, teeth_drop=18.0),
            _row("iron only, fine", 0.0, 0.0, copper=False, teeth_drop=5.0, yoke_drop=3.0)]
    res = tf.evaluate_sector_rows(rows)
    verdicts = {r["variant"]: r["pass"] for r in res["rows"]}
    assert verdicts == {"ok": True, "hot spot too far": False, "ambient too far": False, "liquid: only the internal criterion": True,
                        "iron only, tooth off by 18 %": False, "iron only, fine": True}
    assert res["ok"] is False


def test_tolerances_are_the_ones_written_down_in_advance():
    assert (tf.TOL_CYLINDER_PERCENT, tf.TOL_COMPOSITE_PERCENT, tf.TOL_THIN_LAYER_PERCENT) == (0.5, 1.0, 3.0)
    assert (tf.TOL_SECTOR_INTERNAL_PERCENT, tf.TOL_SECTOR_AMBIENT_PERCENT) == (10.0, 5.0)
