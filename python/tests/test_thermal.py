"""Tests for optimization/thermal.py (the thermal network) and optimization/thermal_geometry.py (the slice for its finite-element check).

Run from P:   uv run pytest tests/test_thermal.py -q      (about a second; no solver and no FEMM are started)

What is held to what:
  1. The building blocks against textbook cases (the rectangle's series, the air-gap correlation's worked example, the single-resistor limit).
  2. The network against itself: heat out equals heat made in every cooling, a copper loss that rises with temperature is solved the same way in the
     matrix and by iteration, the transient ends at the steady state, a limit is a root (the margin there is zero and negative beyond it).
  3. The network against the project's earlier numbers: the reference design's copper loss is the winding study's 56.6 W.
  4. Directions that must hold: a longer stack carries more copper loss, class H is never worse than F, an overload of 2x is shorter than 1.5x.
"""
import math

import numpy as np
import pytest

from optimization import thermal, thermal_geometry
from optimization.thermal import Cooling, Losses, ThermalModel

COOLINGS = (Cooling.natural(), Cooling.forced_air(), Cooling.liquid())


def test_selfcheck_passes():
    failed = [name for name, ok, _ in thermal.selfcheck(verbose=False) if not ok]
    assert not failed, failed


@pytest.mark.parametrize("cooling", COOLINGS, ids=lambda c: c.kind)
def test_energy_balance(cooling):
    m = ThermalModel(stack_mm=60.0, cooling=cooling)
    r = m.steady_state(m.default_losses(copper=m.copper(5.0)))
    made = sum(r.sources.values())
    out = sum(r.heat_out.values())
    assert abs(out - made) < 1e-4 * made, (made, out)
    assert r.converged and not r.runaway


def test_copper_feedback_closed_form_single_node():
    """With huge internal conductances the network is one node: T = T_amb + P(T)/G with P(T) = P20 (1 + alpha (T - 20)), so
    T = (T_amb G + P20 (1 - 20 alpha)) / (G - alpha P20). The matrix solution must give exactly that."""
    m = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air())
    m.G0 = m.G0 * 1e6
    s = m._slot_branch(0.45)
    for k in ("g_wt", "g_wy", "g_we"):
        s[k] *= 1e6
    cl = m.copper(6.0)
    p20 = sum(cl.p20())
    r = m.steady_state(Losses(copper=cl, bearing=0.0, windage=0.0))
    g = sum(m.external(r.T["H"], r.T["B"], 0.0)[:2])
    t_exact = (m.amb * g + p20 * (1 - 20 * thermal.ALPHA_CU)) / (g - thermal.ALPHA_CU * p20)
    assert abs(r.T["W"] - t_exact) < 5e-3, (r.T["W"], t_exact)


def test_reference_design_copper_loss_is_the_winding_studys():
    """The winding study: 1.5 R i^2 = 56.604 W at fill 0.45, 100 C, 100 A-turns peak on the reference design (20 mm stack)."""
    g = thermal.build_geom(thermal.design_params({}))
    m = ThermalModel(overrides={}, stack_mm=20.0)
    s = m._slot_branch(0.45)
    j_pk = 100.0 / (0.45 * g.a_layer * 1e6)
    p = thermal.copper_resistivity(100.0) * (j_pk / math.sqrt(2) * 1e6) ** 2 * (s["v_slot"] + s["v_end"])
    assert abs(p - 56.604) < 0.1, p


def test_end_turn_length_can_be_replaced():
    """Assumptions(end_turn_mm, end_overhang_mm) replace the winding study's half circle (61.15 mm a coil side, 20.2 mm of overhang): the copper in the
    end turns follows the length, the copper in the slots does not, and the overhang sets the end space's length."""
    base = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air())
    short = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air(), assume=thermal.Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6))
    assert base.l_end_side == pytest.approx(0.06115, rel=1e-3) and base.l_overhang == pytest.approx(0.0202, rel=1e-2)
    assert short.l_end_side == pytest.approx(0.04855, rel=1e-12) and short.l_overhang == pytest.approx(0.0136, rel=1e-12)
    assert short.l_es == pytest.approx(0.0136 + short.a.end_clearance_mm * 1e-3, rel=1e-12)
    cb, cs = base.copper(5.0), short.copper(5.0)
    assert cs.slot(100.0) == pytest.approx(cb.slot(100.0), rel=1e-12)
    assert cs.end(100.0) / cb.end(100.0) == pytest.approx(48.55 / base.l_end_side * 1e-3, rel=1e-9)


def test_geometry_of_the_optimized_design():
    g = thermal.build_geom(thermal.design_params())
    assert abs(g.a_slot * 1e6 - 72.436) < 0.01
    assert abs((g.r_so - g.r_yb) * 1e3 - 4.1608) < 1e-3
    assert abs(g.a_layer * 1e6 - 36.218) < 0.01


def test_sustainable_limit_is_a_root():
    for cooling in COOLINGS:
        m = ThermalModel(stack_mm=40.0, cooling=cooling)
        s = m.sustainable()
        ir = m.iron_placeholder()
        base = dict(iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"])
        at = m.steady_state(Losses(copper=m.copper(s["J_rms"]), **base))
        beyond = m.steady_state(Losses(copper=m.copper(s["J_rms"] * 1.01), **base))
        margin = min(at.margin_winding, at.margin_magnet)
        assert 0.0 <= margin < 0.05, (cooling.kind, margin)
        assert min(beyond.margin_winding, beyond.margin_magnet) < 0.0
        assert s["binding"] in ("winding", "magnet")


def test_class_h_is_never_worse_and_stack_helps_where_the_ends_are_cooled():
    for cooling in (Cooling.forced_air(), Cooling.liquid()):
        previous = 0.0
        for stack in (20.0, 40.0, 60.0, 80.0):
            m = ThermalModel(stack_mm=stack, cooling=cooling)
            f, h = m.sustainable(wind_limit_c=155.0), m.sustainable(wind_limit_c=180.0)
            assert h["J_rms"] >= f["J_rms"] - 1e-9
            assert f["P_cu"] > previous            # more copper in a longer stack, more heat shed
            previous = f["P_cu"]


def test_magnet_limit_is_the_owners_hard_limit():
    m = ThermalModel(stack_mm=60.0, cooling=Cooling.forced_air())
    s = m.sustainable()
    assert s["magnet"] <= thermal.MAGNET_LIMIT_C + 1e-3
    assert s["hot_spot"] <= thermal.CLASS_F_C + 1e-3
    hotter = m.sustainable(magnet_limit_c=150.0)
    assert hotter["J_rms"] > s["J_rms"]


def test_callable_copper_and_fixed_number_agree_when_the_function_is_constant():
    m = ThermalModel(stack_mm=40.0, cooling=Cooling.natural())
    a = m.steady_state(m.default_losses(copper=12.0))
    b = m.steady_state(m.default_losses(copper=lambda t: 12.0))
    assert abs(a.T["W"] - b.T["W"]) < 1e-3


def test_transient_ends_at_the_steady_state_and_overloads_order():
    m = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air())
    s = m.sustainable()
    ir = m.iron_placeholder()
    ss = m.steady_state(Losses(copper=m.copper(s["J_rms"]), iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"]))
    tr = m.transient(10 * 3600.0, s["J_rms"], start="cold", dt_max_s=120.0)
    assert abs(tr.T[-1, thermal.IDX["W"]] - ss.T["W"]) < 0.1
    assert np.all(np.diff(tr.T[:, thermal.IDX["W"]]) > -1e-6)            # a heating curve from cold only rises
    o15 = m.overload(1.5, s["J_rms"], start="cold")
    o20 = m.overload(2.0, s["J_rms"], start="cold")
    assert o20["winding_155C_s"] < o15["winding_155C_s"]
    hot = m.overload(2.0, s["J_rms"], start="hot")
    assert hot["winding_155C_s"] < o20["winding_155C_s"]


def test_radial_2d_energy_balance_and_ordering():
    m = ThermalModel(stack_mm=1000.0, cooling=Cooling.natural())
    r = m.radial_2d(p_slot=700.0, p_teeth=40.0, p_yoke=30.0, h_outer=15.0, t_amb=0.0)
    assert abs(r["G_ext"] * r["T_H"] - 770.0) < 1e-6 * 770.0
    assert r["T_hot"] > r["T_W"] > r["T_T"] > r["T_Y"] > r["T_H"] > 0.0


def test_slice_geometry_tiles_without_overlap():
    p = thermal.design_params()
    sl = thermal_geometry.build_slice(p)
    info = thermal_geometry.check_slice(sl, p)
    assert abs(info["slot_mm2"] - info["slot_motorgeo_mm2"]) < 0.02
    assert abs(info["steel_slice_mm2"] - info["steel_motorgeo_mm2_per_pitch"]) < 0.01
    assert abs(info["overlap_mm2"]) < 1e-6 and abs(info["uncovered_mm2"]) < 1e-3


def test_transient_is_exact_for_a_single_capacity():
    """With huge internal conductances the network is one capacity C (the sum of every node's) behind the outside's conductance G, heated by P(T) =
    P20 (1 + alpha (T - 20)) + P_other: dT/dt = (P(T) - G (T - T_amb)) / C, a linear equation with the exact solution T_ss - (T_ss - T0) exp(-t / tau),
    tau = C / (G - alpha P20). The stepping (an exponential per step) must reproduce it to round-off, at any step."""
    m = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air())
    m.G0 = m.G0 * 1e7
    s = m._slot_branch(0.45)
    for k in ("g_wt", "g_wy", "g_we"):
        s[k] *= 1e7
    g_out = 3.0
    m.external = lambda t_h, t_b, q: (2.0, 1.0, 0.0, m.amb)                   # a constant outside: 3 W/K in all
    j = 5.0
    p20 = sum(m.copper(j).p20())
    base = m._iron_losses(None, None)
    other = base.iron_teeth + base.iron_yoke + base.iron_rotor + 0.5 * m.p_windage * 2 + m.p_bearing
    c_tot = float(sum(m.C0) - m.C0[thermal.IDX["W"]] - m.C0[thermal.IDX["E"]] + s["c_w"] + s["c_e"])
    alpha = thermal.ALPHA_CU
    t_ss = (g_out * m.amb + p20 * (1 - 20 * alpha) + other) / (g_out - alpha * p20)
    tau = c_tot / (g_out - alpha * p20)
    r = m.transient(3.0 * tau, j, start="cold", dt_s=tau / 40, dt_max_s=tau / 10)
    expected = t_ss - (t_ss - m.amb) * np.exp(-r.t / tau)
    assert np.max(np.abs(r.T[:, thermal.IDX["R"]] - expected)) < 2e-3, np.max(np.abs(r.T[:, thermal.IDX["R"]] - expected))


def test_correlations_against_known_values():
    # air at 300 K: density 1.177 kg/m3, viscosity 1.846e-5 Pa s, conductivity 0.0263 W/(m K), Pr 0.707
    a = thermal.air_props(26.85)
    assert abs(a["rho"] - 1.177) < 0.01 and abs(a["mu"] - 1.846e-5) < 5e-8 and abs(a["k"] - 0.0263) < 3e-4 and abs(a["Pr"] - 0.707) < 0.01
    # a fin with m (H + t/2) = 1 has the efficiency tanh(1) = 0.7616
    k, t, height = 154.0, 3e-3, 15e-3
    m_h = 1.0 / (height + t / 2)
    h = m_h ** 2 * k * t / 2
    assert abs(thermal.fin_efficiency(h, k, t, height) - math.tanh(1.0)) < 1e-12
    # Gnielinski turbulent value at Re = 1e4, Pr = 0.71 is Nu about 33 (textbook)
    assert 28.0 < thermal.nusselt_duct(1.0e4, 0.71) < 38.0
    assert thermal.nusselt_duct(2000.0, 0.71, 4.0) == 4.0
    # the maker's typical loss for 29-gauge M-19: 3.19 W/kg at 1.5 T and 60 Hz, 1.36 W/kg at 1.0 T (materials note 02 Current)
    assert abs(thermal.iron_w_per_kg(1.5, 60.0) - 3.19) < 0.01
    assert abs(thermal.iron_w_per_kg(1.0, 60.0) - 1.36) < 0.01
    # copper's resistivity: 1.7241e-8 at 20 C and 39 % more at 120 C (+0.393 %/K)
    assert abs(thermal.copper_resistivity(120.0) / thermal.copper_resistivity(20.0) - 1.393) < 1e-3


def test_air_gap_regimes_and_speed_dependence():
    still = thermal.air_gap_h(0.0, 0.045, 0.0008)
    assert still["regime"].startswith("laminar") and abs(still["Nu"] - 2.018) < 0.01
    rated = thermal.air_gap_h(2 * math.pi * 3000 / 60, 0.045, 0.0008)
    fast = thermal.air_gap_h(2 * math.pi * 16330 / 60, 0.045, 0.0008)
    assert rated["regime"] == "Taylor vortices" and fast["regime"].startswith("turbulent")
    assert still["h"] < rated["h"] < fast["h"]
