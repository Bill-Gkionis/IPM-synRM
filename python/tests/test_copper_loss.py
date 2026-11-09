"""Tests for the copper-loss work (optimization/ac_theory.py, ac_pwm.py, ac_slot.py's geometry, copper_loss.py): everything that needs no FEMM.

  test_ac_theory_*       the exact wire results agree with each other (Bessel against Kelvin functions, closed form against quadrature) and with the first study's own proximity_ratio
  test_pwm_*             the ripple's closed forms against Parseval, the one-leg textbook case and an independent time-domain simulation
  test_packing_*         the random strand packing: no two strands overlap, all inside their layer, the layer areas of the optimized design; the turn assignments hold n_s strands a turn
  test_copper_loss_*     the interface (keys, positivity, units), the DC loss against winding.py with the first study's end-turn stand-in, the scaling laws (p_prox with the speed squared,
                         zero at standstill), the temperature, the proximity loss at the first study's own point (its formula, within the fit's 0.5 % residual), the ripple's 1 / f_sw, no
                         circulating current for a single strand, strands in a line against in a row, the array function against the scalar one, and the speed per call. Skipped if
                         copper_loss_tables.json has not been built (python -m optimization.copper_loss --build).

Run from the Python to FEMM folder:   uv run pytest tests/test_copper_loss.py -q
"""
import json
import math
import time

import numpy as np
import pytest

import paths
from optimization import ac_pwm, ac_slot, ac_theory, copper_loss, design, winding


# ------------------------------------------------------------------------------------------------------ ac_theory
def test_ac_theory_closed_forms_agree():
    rho = float(ac_theory.resistivity(100.0))
    d = 0.63e-3
    for ratio in (0.1, 0.5, 1.0, 2.0, 3.0, 5.0):
        f = rho / (math.pi * ac_theory.MU0 * (d / ratio) ** 2)
        assert float(ac_theory.skin_ratio(d, f, rho)) == pytest.approx(float(ac_theory.skin_ratio_kelvin(d, f, rho)), rel=1e-10)
        assert float(ac_theory.proximity_ratio(d, f, rho)) == pytest.approx(ac_theory.proximity_ratio_quadrature(d, f, rho), rel=1e-9)
        assert float(ac_theory.proximity_ratio(d, f, rho)) == pytest.approx(winding.proximity_ratio(d, f, rho), rel=1e-9)


def test_ac_theory_known_values():
    rho = float(ac_theory.resistivity(100.0))
    d = 0.63e-3
    f1 = rho / (math.pi * ac_theory.MU0 * d ** 2)                  # d / delta = 1
    assert float(ac_theory.skin_ratio(d, f1, rho)) == pytest.approx(1.0013007, abs=2e-7)
    assert float(ac_theory.proximity_ratio(d, f1, rho)) == pytest.approx(0.9928916, abs=2e-7)
    assert float(ac_theory.two_wire_circulating(d, d, 1000.0, 1.0, rho)) == pytest.approx(4.0 * 2.0 * float(ac_theory.proximity_lowfreq(d, 1000.0, 1.0, rho)))


# ------------------------------------------------------------------------------------------------------ ac_pwm
def test_pwm_parseval_and_leg():
    theta = 2.0 * math.pi * (np.arange(360) + 0.5) / 360
    r = ac_pwm.ripple_period(ac_pwm.duty_cycles(theta, 0.8), 0, 200)
    assert float(r["harm"].sum(axis=1).mean()) == pytest.approx(float(r["ms"].mean()), rel=2e-6)
    assert ac_pwm.ripple_vs_m([0.0], 360)["h_rms"][0] == 0.0


def test_pwm_matches_time_domain_simulation():
    for m in (0.4, 1.0):
        analytic = ac_pwm.ripple_vs_m([m], 360)["h_rms"][0]
        assert ac_pwm.simulate_rms(m, ratio=120, steps=1500) == pytest.approx(analytic, rel=3e-3)


# ------------------------------------------------------------------------------------------------------ the strand packing
@pytest.fixture(scope="module")
def optimized():
    ov = json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]
    return design.params(ov)


def test_packing_no_overlap_inside(optimized):
    polys = ac_slot.layer_polygons(optimized)
    assert polys[0].area == pytest.approx(36.218, abs=1e-3) and polys[1].area == pytest.approx(36.218, abs=1e-3)
    lay = ac_slot.strand_layout(optimized, 208, 0.367, seed=3)
    from scipy.spatial import cKDTree
    for li in range(2):
        xy = lay["xy"][lay["layer"] == li]
        assert len(xy) == 208
        assert len(cKDTree(xy).query_pairs(0.367 - 1e-6)) == 0
        import shapely
        assert shapely.contains_xy(polys[li].buffer(-0.367 / 2 + 1e-7), xy[:, 0], xy[:, 1]).all()


@pytest.mark.parametrize("arrangement", ac_slot.ARRANGEMENTS)
def test_turn_assignment_holds_n_s_strands_a_turn(optimized, arrangement):
    lay = ac_slot.strand_layout(optimized, 13 * 4, 0.7042, seed=1)
    turn = ac_slot.assign_turns(lay["xy"], lay["layer"], 13, 4, arrangement, seed=5)
    for li in range(2):
        counts = np.bincount(turn[lay["layer"] == li], minlength=13)
        assert counts.tolist() == [4] * 13


def test_line_arrangement_spans_the_depth(optimized):
    lay = ac_slot.strand_layout(optimized, 26 * 8, 0.367, seed=2)
    spread = {}
    for arr in ("line", "row"):
        turn = ac_slot.assign_turns(lay["xy"], lay["layer"], 26, 8, arr, seed=1)
        sel = lay["layer"] == 0
        spread[arr] = np.mean([np.ptp(lay["xy"][sel & (turn == t), 1]) for t in range(26)])
    assert spread["line"] > 4.0 * spread["row"]


# ------------------------------------------------------------------------------------------------------ copper_loss
tables_missing = not copper_loss.TABLES_PATH.exists()
needs_tables = pytest.mark.skipif(tables_missing, reason="copper_loss_tables.json is not built")

W1 = {"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367}
W2 = {"N": 13, "a": 2, "n_s": 4, "d_mm": 0.630, "d_over_mm": 0.7042}
OP = {"i_peak": 18.25, "gamma_deg": 22.47, "rpm": 3000.0, "f_sw": 20000.0, "v_dc": 48.0}
KEYS = ("p_dc_slot", "p_dc_end", "p_prox", "p_skin", "p_pwm", "p_circ")


@needs_tables
def test_copper_loss_interface():
    r = copper_loss.losses(W1, OP, 100.0, 20.0)
    for k in KEYS:
        assert k in r and r[k] >= 0.0 and math.isfinite(r[k]), k
    assert r["p_total"] == pytest.approx(sum(r[k] for k in KEYS))
    assert r["p_circ"] == 0.0                                  # the default arrangement is transposed


@needs_tables
def test_copper_loss_dc_matches_the_first_study():
    """With the first study's end turn (the half circle, 2 mm legs) the DC loss is its R i^2 formula: winding.Model.dc_loss_at_fill at the example's fill, at 100 ampere-turns."""
    p = design.params(json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"])
    w = dict(W1, end_model="arc", end_ext_mm=2.0)
    op = dict(OP, i_peak=100.0 * 4 / 26, rpm=0.0)
    r = copper_loss.losses(w, op, 100.0, 20.0)
    g = winding.end_turn_geometry(p)
    rho = winding.resistivity(100.0)
    a_s = math.pi * (0.315e-3) ** 2 / 4
    R = rho * g["l_mt_mm"] * 1e-3 * 12 * 26 / (4 ** 2 * 8 * a_s)
    assert r["p_dc_slot"] + r["p_dc_end"] == pytest.approx(1.5 * R * op["i_peak"] ** 2, rel=1e-9)


@needs_tables
def test_copper_loss_scalings():
    base = copper_loss.losses(W1, dict(OP, rpm=3000.0), 100.0, 20.0)
    fast = copper_loss.losses(W1, dict(OP, rpm=6000.0), 100.0, 20.0)
    still = copper_loss.losses(W1, dict(OP, rpm=0.0), 100.0, 20.0)
    assert still["p_prox"] == 0.0 and still["p_skin"] == 0.0
    assert fast["p_prox"] / base["p_prox"] == pytest.approx(4.0, rel=0.02)             # the speed squared (d / delta is small)
    hot = copper_loss.losses(W1, OP, 150.0, 20.0)
    assert hot["p_dc_slot"] / base["p_dc_slot"] == pytest.approx((1 + 0.00393 * 130) / (1 + 0.00393 * 80), rel=1e-9)
    long_stack = copper_loss.losses(W1, OP, 100.0, 40.0)
    assert long_stack["p_dc_end"] == pytest.approx(base["p_dc_end"], rel=1e-9)         # the end turns do not depend on the stack
    assert long_stack["p_dc_slot"] == pytest.approx(2.0 * base["p_dc_slot"], rel=1e-9)


@needs_tables
def test_copper_loss_speed_per_call():
    t0 = time.perf_counter()
    n = 400
    for _ in range(n):
        copper_loss.losses(W1, OP, 100.0, 20.0)
    per_call = (time.perf_counter() - t0) / n
    assert per_call < 5e-3                                      # the optimizer needs microseconds to milliseconds


@needs_tables
@pytest.mark.parametrize("w", [W1, W2], ids=["W1", "W2"])
def test_copper_loss_proximity_matches_the_first_study_at_its_point(w):
    """At the first study's own point (100 ampere-turns a coil side, 19.657 deg, 3,000 rpm: G = 4.8222292e-4 T^2/rad^2 over the coil area, from out/optimization/winding/optimized/femm.json)
    the proximity loss is the first study's low-frequency formula 72 N n_s L (pi sigma d^4 w^2 / 64) G; the fast model's fit reproduces that point to the fit's residual (0.5 % here)."""
    g_first = 4.8222292381e-04
    rho = float(ac_theory.resistivity(100.0))
    d, omega, length = w["d_mm"] * 1e-3, 2.0 * math.pi * 200.0, 20e-3
    formula = 72 * w["N"] * w["n_s"] * length * (math.pi / rho * d ** 4 * omega ** 2 / 64.0) * g_first
    op = dict(OP, i_peak=100.0 * w["a"] / w["N"], gamma_deg=19.656668507770743, rpm=3000.0)
    r = copper_loss.losses(w, op, 100.0, 20.0)
    assert not r["extrapolated"]
    assert r["p_prox"] == pytest.approx(formula, rel=0.01)


@needs_tables
def test_copper_loss_ripple_goes_as_one_over_the_switching_frequency():
    a = copper_loss.losses(W1, dict(OP, f_sw=10000.0), 100.0, 20.0)
    b = copper_loss.losses(W1, dict(OP, f_sw=20000.0), 100.0, 20.0)
    assert a["ripple_rms_a"] / b["ripple_rms_a"] == pytest.approx(2.0, rel=1e-9)


@needs_tables
@pytest.mark.parametrize("arrangement", ["row", "cluster", "random", "line"])
def test_copper_loss_no_circulating_current_for_a_single_strand(arrangement):
    """A single strand has nobody to share the turn's current with: no circulating current, whatever the placement (the module's 'self-check' 3)."""
    r = copper_loss.losses(dict(W1, n_s=1, arrangement=arrangement), dict(OP, rpm=16330.0), 100.0, 20.0)
    assert r["p_circ"] == 0.0


@needs_tables
def test_copper_loss_parallel_strands_in_a_line_lose_more_than_in_a_row():
    row = copper_loss.losses(dict(W1, arrangement="row"), dict(OP, rpm=10000.0), 100.0, 20.0)["p_circ"]
    line = copper_loss.losses(dict(W1, arrangement="line"), dict(OP, rpm=10000.0), 100.0, 20.0)["p_circ"]
    assert line > 10.0 * row > 0.0


@needs_tables
def test_copper_loss_arrays_agree_with_the_scalar_function():
    speeds = np.array([1000.0, 5000.0, 16330.0])
    arr = copper_loss.losses_arrays(W2["N"], W2["a"], W2["n_s"], W2["d_mm"], W2["d_over_mm"], OP["i_peak"], OP["gamma_deg"], speeds, OP["f_sw"], OP["v_dc"], 100.0, 20.0)
    for i, rpm in enumerate(speeds):
        one = copper_loss.losses(W2, dict(OP, rpm=float(rpm)), 100.0, 20.0)
        for k in KEYS:
            assert float(arr[k][i]) == pytest.approx(one[k], rel=1e-12, abs=1e-15), (k, rpm)
