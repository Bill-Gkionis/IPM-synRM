"""Tests for validation/provenance.py: the small analyses behind the report's hand-made numbers, on made-up data with known answers.

No FEMM is started. Run from P:   uv run pytest tests/test_provenance.py -q
"""
import math

import numpy as np
import pandas as pd
import pytest
from PIL import Image

import provenance as pv

DEG = np.arange(0, 90, 0.25)                                         # the 360 samples of the air-gap files


def waveform(deg):
    """A smooth waveform with the sector's period (90 degrees): one cycle of 0.25 T and a little third harmonic."""
    return 0.25 * np.sin(np.radians(4 * deg)) + 0.03 * np.sin(np.radians(12 * deg))


# ----------------------------------------------------------------------------------------------- F4: the air-gap flux density
def test_the_slope_is_a_central_difference_in_tesla_per_degree():
    b = 0.25 * np.sin(np.radians(4 * DEG))
    slope = pv.slope_per_degree(DEG, b)
    assert slope.max() == pytest.approx(0.25 * math.radians(4), rel=1e-3)           # the steepest slope of A sin(k x) is A k
    assert slope[0] == pytest.approx(0.25 * math.radians(4), rel=1e-3)               # at 0 degrees, where the sine crosses zero
    assert slope[90] < 1e-3                                                         # at 22.5 degrees the sine is at its peak: flat


def test_the_best_shift_is_the_angle_a_waveform_was_moved_by():
    old = waveform(DEG)
    new = waveform(DEG + 0.12)                                                       # new(theta) = old(theta + 0.12)
    assert pv.best_shift(DEG, old, new) == pytest.approx(-0.12, abs=0.011)           # so new has to be moved back by 0.12 degrees
    assert pv.best_shift(DEG, old, old) == 0.0                                       # equal waveforms need no shift
    noisy = old + 0.002 * np.random.default_rng(1).standard_normal(len(DEG))
    assert pv.best_shift(DEG, old, noisy) == pytest.approx(0.0, abs=0.05)            # noise is not a shift


def test_the_worst_sample_its_slope_and_the_steep_overlap_are_reported():
    old = pd.DataFrame({"deg": DEG, "br": waveform(DEG)})
    new = old.copy()
    k = int(np.argmax(pv.slope_per_degree(DEG, old.br.values)))                      # the steepest sample
    new.loc[k, "br"] += 0.01                                                         # one difference, there
    r = pv.gap_case(old, new)
    peak = np.abs(old.br).max()
    assert r["worst_deg"] == DEG[k] and r["max_pct"] == pytest.approx(100 * 0.01 / peak)
    assert len(r["top"]) == 5 and r["top"][0][0] == DEG[k] and r["top"][0][1] == pytest.approx(100 * 0.01 / peak)      # the bump, signed, first
    assert r["slope_at_worst"] == pytest.approx(r["steepest"])                       # the worst sample IS the steepest
    assert r["steep_overlap"] >= 1 and r["quarter"] == 90
    assert r["shift_deg"] == pytest.approx(0.0, abs=0.011)                           # one bumped sample is not a shift
    assert r["pct_of_peak_per_0p02_deg_at_steepest"] == pytest.approx(100 * 0.02 * r["steepest"] / peak)


def test_the_ladder_comparison_needs_its_files_and_says_so(tmp_path):
    assert pv.gap_ladder(tmp_path, tmp_path / "mesh") is None                        # nothing has run: None, not a crash


# ----------------------------------------------------------------------------------------------- F1: the offset of psi_d
def write_phase1(folder, offset_wb):
    """A phase1.csv with a 3 x 3 torque-vs-angle grid, two positions each, whose psi_d is `offset_wb` higher than the base."""
    rows = []
    for variant in ("with band", "no band"):
        for i in (5, 10, 15):
            for g in (0, 30, 60):
                for th in (0, 10):
                    base = 0.012 - 0.0002 * i * math.sin(math.radians(g)) + 0.0001 * th / 10
                    rows.append({"study": "gamma", "variant": variant, "theta_e": th, "I": i, "gamma": g, "psi_d": base + offset_wb,
                                 "psi_q": 0.02})
    (folder / "barrier").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(folder / "barrier" / "phase1.csv", index=False)


def test_a_constant_offset_in_psi_d_is_found_with_its_size_and_its_signs(tmp_path):
    write_phase1(tmp_path / "old", 0.0)
    write_phase1(tmp_path / "new", 0.00004)                                          # +0.04 mWb in every state
    pd.DataFrame({"with band": [12.0], "no band": [16.0]}, index=["psi_pm_mWb"]).to_csv(tmp_path / "old" / "barrier" / "summary.csv")
    rows = pv.psi_d_offset(tmp_path / "old", tmp_path / "new")
    assert [r["variant"] for r in rows] == ["with band", "no band"]
    r = rows[0]
    assert r["states"] == 9 and r["positive"] == 9
    assert r["mean"] == pytest.approx(0.04, abs=1e-9) and r["std"] == pytest.approx(0.0, abs=1e-9)
    assert r["pct_of_magnet_flux"] == pytest.approx(100 * 0.04 / 12.0)
    assert r["pct_of_largest"] == pytest.approx(100 * 0.04 / 20.0)                   # psi_q = 20 mWb is the largest in the grid


# ----------------------------------------------------------------------------------------------- the pandas trap
def test_pandas_misreads_are_counted_against_python_floats(tmp_path):
    path = tmp_path / "evaluations.csv"
    path.write_text("gen,mag_y,mag_x,v_half_deg,web,feasible\n"
                    "0,5.0,39.227554466990604,60.0,1.2227404896806187,True\n"
                    "0,5.5,39.0,61.0,1.0,False\n", encoding="utf-8")                  # the infeasible row is not counted
    wrong, total = pv.pandas_misreads(path)
    assert total == 4 and 1 <= wrong <= 2                                            # the two long numbers of the optimizer's own CSV


# ----------------------------------------------------------------------------------------------- the picture
def test_the_pictures_are_compared_by_halves_and_by_channel(tmp_path):
    a = np.zeros((4, 8, 4), dtype=np.uint8)
    b = a.copy()
    b[1, 6] = (40, 10, 0, 0)                                                         # the right half is the second half of the columns
    b[2, 7] = (0, 0, 0, 0)                                                           # unchanged
    Image.fromarray(a, "RGBA").save(tmp_path / "a.png")
    Image.fromarray(b, "RGBA").save(tmp_path / "b.png")
    r = pv.picture_difference(tmp_path / "a.png", tmp_path / "b.png")
    assert (r["differ"], r["total"], r["left"], r["right"]) == (1, 32, 0, 1)
    assert r["max_channel"] == 40 and r["brighter"] == 1 and r["darker"] == 0
    Image.fromarray(np.zeros((4, 4, 4), dtype=np.uint8), "RGBA").save(tmp_path / "c.png")
    assert pv.picture_difference(tmp_path / "a.png", tmp_path / "c.png") == {"same_size": False}
