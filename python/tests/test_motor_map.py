"""Tests of the motor map's engine (optimization/motor_map.py, motor_map_thermal.py, motor_map_stack.py). No FEMM and no solver slot: the saved flux and iron-loss maps
(out\\optimization\\maps\\chosen\\, which maps_report.py writes; the tests that need them are skipped without them).

Run from the code folder:   uv run pytest tests/test_motor_map.py -q      (about 4 s)

In groups:
  1. The inputs: the winding family, the mechanical estimate against thermal.py's formulas, the polygon formula of the inertia.
  2. The machine: the scaling rule, the summed iron-loss splines against the map, one operating point's energy balance, the swap of the copper model (DC against copper_loss.py).
  3. The solver: the envelope against the MTPA table, the corner against its closed form, the MTPV branch against the Lagrange condition, the level sets against their targets.
  4. The thermal coupling: the continuous current at the rated speed against the network's own sustainable().
  5. What ends an envelope: the top speed, the loads or the temperatures (motor_map_run.summarise).
  6. The winding study's windings: the fill, the resistance, the rating and the torque-speed table of the headline winding against the study's own numbers, the shaft fan sized on its stack,
     the inputs of a saved run, and the energy balance with the full copper loss.
"""
import math

import numpy as np
import pytest

from optimization import maps
from optimization import motor_map as mm
from optimization import motor_map_stack as st
from optimization import motor_map_thermal as mt
from optimization import thermal

HAVE = maps.FLUX_MAP.exists() and maps.iron_loss.__module__ and (maps.FLUX_MAP.parent / "iron_map.npz").exists()
pytestmark = pytest.mark.skipif(not HAVE, reason="the saved flux and iron-loss maps are not there")

W = mm.PROVISIONAL_WINDING
BUS = mm.Bus(48.0, 0.95)


# ---------------------------------------------------------------------------------------------- 1. the inputs
def test_scaled_winding_keeps_the_fill_and_the_copper_in_the_slot():
    for f in (0.5, 1.0 / 3.0, 3.0):
        w = W.scaled(f)
        m0, m1 = mm.Machine(W, 20.0), mm.Machine(w, 20.0)
        assert w.turns == pytest.approx(W.turns * f)
        assert m1.fill == pytest.approx(m0.fill, rel=1e-12)
        assert w.turns * w.copper_area_m2 == pytest.approx(W.turns * W.copper_area_m2, rel=1e-12)


def test_current_density_round_trip_and_the_example_winding():
    i = W.current_from_density(6.0)
    assert i == pytest.approx(math.sqrt(2) * 6e6 * 4 * 8 * math.pi * (0.315e-3) ** 2 / 4, rel=1e-12)
    assert W.density_from_current(i) == pytest.approx(6.0, rel=1e-12)
    assert W.series_turns == 78.0


def test_mechanical_loss_is_thermals_formula():
    mech = mm.Mechanical()
    m = mm.Machine(W, 60.0)
    for rpm in (500.0, 3000.0, 16330.0):
        w = 2 * math.pi * rpm / 60
        p = mech.power(rpm, m.r_ro, m.gap, m.L)
        assert float(p["windage"]) == pytest.approx(thermal.windage_loss(w, m.r_ro, m.gap, m.L, 70.0, 1.0), rel=1e-12)
        assert float(p["bearing"]) == pytest.approx(thermal.bearing_loss(w, thermal.Assumptions()), rel=1e-12)
    assert float(mech.power(0.0, m.r_ro, m.gap, m.L)["total"]) == 0.0


def test_shaft_fan_is_a_load_on_the_shaft():
    fan = lambda r: 1e-8 * np.asarray(r, float) ** 3
    mech = mm.Mechanical(fan=fan)
    p = mech.power(np.array([1000.0, 2000.0]), 0.045, 0.0008, 0.06)
    assert p["fan"] == pytest.approx([10.0, 80.0])
    assert p["total"] == pytest.approx(p["bearing"] + p["windage"] + p["fan"])


def test_polar_moment_of_a_ring_is_the_closed_form():
    from shapely.geometry import Point
    ring = Point(0, 0).buffer(45.0, 1000).difference(Point(0, 0).buffer(7.5, 1000))
    assert st.polar_moment_mm4(ring) == pytest.approx(math.pi / 2 * (45.0 ** 4 - 7.5 ** 4), rel=1e-5)


def test_bus_limit_overrides_the_modulation_limit_and_keeps_the_bus_voltage():
    assert mm.Bus(48.0, 0.95).v_max == pytest.approx(0.95 * 48.0 / math.sqrt(3))
    free = mm.Bus(48.0, 1.0, limit=1e9)
    assert free.v_max == 1e9 and free.v_dc == 48.0
    assert mt.NO_VOLTAGE_LIMIT.v_dc == 48.0 and mt.NO_VOLTAGE_LIMIT.v_max == 1e9


# ---------------------------------------------------------------------------------------------- 2. the machine
def test_flux_follows_the_scaling_rule():
    m = mm.Machine(W, 60.0)
    i_d, i_q = -9.0, 17.0
    s, k = m.s, 3.0
    pd0, pq0 = m.maps.psi_dq(s * i_d, s * i_q)                       # the 10-turn, 20 mm placeholder at the scaled current
    pd, pq = m.flux(i_d, i_q)
    assert pd == pytest.approx(s * k * pd0, rel=1e-12) and pq == pytest.approx(s * k * pq0, rel=1e-12)
    assert m.torque_dq(i_d, i_q) == pytest.approx(k * m.maps.torque_dq(s * i_d, s * i_q), rel=1e-12)


def test_summed_iron_splines_equal_the_regions_of_the_map():
    m = mm.Machine(W, 60.0)
    rng = np.random.default_rng(1)
    I0 = rng.uniform(0.0, 3.0, 50) * m.maps.rated
    g = np.radians(rng.uniform(0.0, 90.0, 50))
    i_d, i_q, rpm = -I0 * np.sin(g) / m.s, I0 * np.cos(g) / m.s, rng.uniform(50.0, 16330.0, 50)
    c = m.iron_coeffs(i_d, i_q)
    f = mm.C_SPEED * 2 * math.pi * rpm / 60
    assert c[0] * f + c[1] * f ** 2 + c[2] * f ** 1.5 == pytest.approx(m.iron(i_d, i_q, rpm)["total"], rel=1e-9)


def test_an_operating_point_closes_its_energy_balance():
    m = mm.Machine(W, 60.0)
    i_d, i_q = np.array([-6.0, -15.0, -2.0]), np.array([12.0, 8.0, 3.0])
    rpm = np.array([800.0, 6000.0, 12000.0])
    pt = m.point(i_d, i_q, rpm, BUS, 100.0)
    assert np.all(m.balance_error(pt) < 1e-12)
    assert pt["P_out"] == pytest.approx(pt["T_shaft"] * 2 * math.pi * rpm / 60)
    assert pt["P_in"] == pytest.approx(1.5 * (pt["v_d"] * i_d + pt["v_q"] * i_q))


def test_the_dc_copper_model_is_copper_loss_with_the_first_studys_end_turn():
    from optimization import copper_loss
    m = mm.Machine(W, 20.0)
    dc = m.copper.loss(m, -4.0, 14.0, 3000.0, 100.0, 48.0)["total"]
    r = copper_loss.losses({"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367, "end_model": "arc", "end_ext_mm": 2.0},
                           {"i_peak": math.hypot(4.0, 14.0), "gamma_deg": math.degrees(math.atan2(4.0, 14.0)), "rpm": 3000.0, "f_sw": 0.0, "v_dc": 48.0}, 100.0, 20.0)
    assert dc == pytest.approx(r["p_dc_slot"] + r["p_dc_end"], rel=1e-9)
    assert float(m.resistance(100.0)) == pytest.approx(r["r_dc_ohm"], rel=1e-9)


def test_full_copper_adds_to_the_dc_loss():
    m = mm.Machine(W, 20.0, copper=mm.FullCopper())
    cu = m.copper.loss(m, -4.0, 14.0, 6000.0, 100.0, 48.0)
    assert cu["total"] == pytest.approx(cu["dc"] + cu["ac"]) and cu["ac"] > 0.0 and cu["dc"] > 0.0


# ---------------------------------------------------------------------------------------------- 3. the solver
def test_envelope_equals_the_mtpa_table_below_the_corner():
    m = mm.Machine(W, 60.0)
    lim = W.current_from_density(6.0)
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    g, _, _, T = m.maps.mtpa(lim, turns=W.turns, paths=W.paths, stack_mm=60.0)
    s = sc.max_torque(sc.ctx(100.0, BUS, 100.0, "air_gap"), lim)
    assert s["T"] == pytest.approx(T, rel=1e-6) and s["gam"] == pytest.approx(g, abs=0.02)
    assert not s["v_active"] and s["i_active"]


def test_corner_speed_is_the_roots_of_the_first_studys_quadratic():
    m = mm.Machine(W, 20.0)
    lim = 20.0
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    env = mm.envelope(m, BUS, mm.speed_grid(step=500.0), lim, basis="air_gap", scan=sc)
    rc = env.meta["corner_rpm"]
    s = sc.max_torque(sc.ctx(rc, BUS, 100.0, "air_gap"), lim)
    pd, pq = m.flux(s["i_d"], s["i_q"])
    R = float(m.resistance(100.0))
    A, B, C = pd ** 2 + pq ** 2, 2 * R * (s["i_q"] * pd - s["i_d"] * pq), R ** 2 * s["I"] ** 2 - BUS.v_max ** 2
    root = (-B + math.sqrt(B * B - 4 * A * C)) / (2 * A) * 60 / (2 * math.pi * 4)
    assert rc == pytest.approx(root, rel=1e-6)


def test_regions_follow_each_other_in_order_and_torque_never_rises_above_the_corner():
    m = mm.Machine(W, 20.0)
    lim = 2.0 * W.current_from_density(6.0)
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    env = mm.envelope(m, BUS, mm.speed_grid(step=250.0), lim, scan=sc)
    seq = [r for r in env["region"] if r != "off"]
    rank = {"MTPA": 0, "FW": 1, "MTPV": 2}
    assert [rank[r] for r in seq] == sorted(rank[r] for r in seq)
    assert set(seq) == {"MTPA", "FW", "MTPV"}
    T = env["T"]
    above = env["rpm"] > env.meta["corner_rpm"]
    assert np.all(np.diff(T[above]) <= 1e-9)


def test_the_mtpv_branch_satisfies_the_lagrange_condition():
    w = W
    m = mm.Machine(w, 20.0)
    lim = 2.0 * w.current_from_density(6.0)
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    env = mm.envelope(m, BUS, mm.speed_grid(step=500.0), lim, basis="air_gap", scan=sc, locate=False)
    sel = np.nonzero(env["region"] == "MTPV")[0]
    assert len(sel) > 5
    i_d, i_q, rpm = env["i_d"][sel], env["i_q"][sel], env["rpm"][sel]
    L = m.maps.inductances(i_d, i_q, turns=w.turns, paths=w.paths, stack_mm=20.0)
    R = float(m.resistance(100.0))
    we = 4 * 2 * math.pi * rpm / 60
    pd, pq = m.flux(i_d, i_q)
    vd, vq = R * i_d - we * pq, R * i_q + we * pd
    gT = np.stack([6.0 * (L["L_dd"] * i_q - L["L_qd"] * i_d - pq), 6.0 * (pd + L["L_dq"] * i_q - L["L_qq"] * i_d)])
    gV = np.stack([vd * (R - we * L["L_qd"]) + vq * we * L["L_dd"], -vd * we * L["L_qq"] + vq * (R + we * L["L_dq"])])
    sin = np.abs(gT[0] * gV[1] - gT[1] * gV[0]) / (np.hypot(*gT) * np.hypot(*gV))
    assert sin.max() < 1e-4


def test_level_sets_hit_their_torque_and_minimum_loss_is_not_worse_than_mtpa():
    m = mm.Machine(st.stack_winding(W, 60.0, "hold"), 60.0)
    lim = m.w.current_from_density(6.0)
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    for rpm in (1500.0, 5000.0):
        c = sc.ctx(rpm, BUS, 100.0, "shaft")
        top = sc.max_torque(c, lim)["T"]
        targets = np.array([0.0, 0.1, 0.4, 0.8, 0.99]) * top
        res = sc.level_set(c, lim, targets)
        loss = {}
        for name in ("loss", "current"):
            r = res[name]
            assert r["found"].all()
            g = np.radians(r["gam"])
            i_d, i_q = -r["I"] * np.sin(g), r["I"] * np.cos(g)
            pt = m.point(i_d, i_q, rpm, BUS, 100.0)
            assert pt["T_shaft"] == pytest.approx(targets, abs=1e-8)
            assert np.all(pt["v_abs"] <= BUS.v_max * (1 + 1e-8)) and np.all(pt["I"] <= lim * (1 + 1e-8))
            loss[name] = pt["P_loss"]
        assert np.all(loss["loss"] <= loss["current"] * (1 + 1e-9))


def test_a_thin_feasible_lens_at_the_current_limit_is_found():
    m = mm.Machine(W, 60.0)
    lim = W.current_from_density(6.0)
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    c = sc.ctx(4500.0, BUS, 100.0, "shaft")                             # the feasible region is a lens a fraction of an ampere deep
    top = sc.max_torque(c, lim)["T"]
    res = sc.level_set(c, lim, np.array([0.2, 0.5]) * top)
    assert res["loss"]["found"].all() and res["current"]["found"].all()


def test_a_level_set_whose_best_point_is_the_exit_of_the_voltage_curve_is_found():
    """The target found missing in the first run (the 60 mm / 60 V short-term map at 1,948 rpm, half the envelope's torque, on the coarse grid): its best point is where the level
    curve LEAVES the voltage curve (the exit), which neither the rays' nor the circles' crossings reach alone."""
    fan = lambda r: thermal.shaft_fan_power(r, thermal.Cooling.shaft_fan(6.0 * 3000.0 / mm.TOP_RPM, 3000.0), 60.0)
    m = mm.Machine(W, 60.0, mech=mm.Mechanical(fan=fan))
    bus = mm.Bus(60.0, 0.95)
    lim = 29.393754585275182
    sc = mm.Scan(m, i_max=lim, n_gamma=181, n_radius=121)
    c = sc.ctx(1947.7777777777776, bus, 100.0, "shaft")
    top = sc.max_torque(c, lim)["T"]
    target = np.array([0.5025641025641026 * top])
    res = sc.level_set(c, lim, target)
    for name in ("loss", "current"):
        assert res[name]["found"].all()
        g = np.radians(res[name]["gam"])
        pt = m.point(-res[name]["I"] * np.sin(g), res[name]["I"] * np.cos(g), 1947.7777777777776, bus, 100.0)
        assert pt["T_shaft"] == pytest.approx(target, abs=1e-8) and pt["v_abs"][0] <= bus.v_max * (1 + 1e-8)


# ---------------------------------------------------------------------------------------------- 4. the thermal coupling
def test_thermal_rating_follows_the_network_within_the_iron_loss_difference():
    m = mm.Machine(W, 20.0)
    cool = thermal.Cooling.forced_air()
    mine = mt.thermal_rating(m, mt.ThermalConfig(cool), 3000.0)
    net = thermal.get_model(None, 20.0, cool, 3000.0, m.fill, thermal.Assumptions()).sustainable()
    assert mine["AT_pk"] == pytest.approx(net["AT_pk"], rel=0.03)
    assert mine["binding"] == net["binding"] == "magnet"


def test_the_continuous_limit_falls_where_the_iron_loss_grows_and_stays_inside_the_limits():
    m = mm.Machine(W, 20.0)
    cfg = mt.ThermalConfig(thermal.Cooling.forced_air())
    rpms = np.array([0.0, 1000.0, 2000.0, 3000.0])
    cl = mt.continuous_limits(m, BUS, cfg, rpms, n_gamma=181, n_radius=121)
    assert np.all(cl["limited"] == "temperature")
    assert np.all(np.diff(cl["I"]) < 0)
    assert np.all(cl["magnet"] <= thermal.MAGNET_LIMIT_C + 0.05) and np.all(cl["hot_spot"] <= thermal.CLASS_F_C + 0.05)


def test_stack_geometry_scales_with_the_stack():
    g20, g60 = st.geometry_numbers(mm.Machine(W, 20.0)), st.geometry_numbers(mm.Machine(W, 60.0))
    assert g60["n_magnets"] == 3 * g20["n_magnets"] == 48
    assert g60["m_rotor_steel"] == pytest.approx(3 * g20["m_rotor_steel"], rel=1e-12)
    assert g60["m_stator_steel"] == pytest.approx(3 * g20["m_stator_steel"], rel=1e-12)
    assert 0.0008 < g20["J_rotor"] < 0.0012

def test_overload_from_cold_is_the_networks_overload_with_the_same_inputs():
    m = mm.Machine(W, 20.0)
    cool = thermal.Cooling.shaft_fan()
    cfg = mt.ThermalConfig(cool)
    th = thermal.get_model(None, 20.0, cool, 3000.0, m.fill, thermal.Assumptions())
    j = th.sustainable()["J_rms"]
    ref = th.overload(2.0, j, start="cold", rpm=3000.0)
    mine = mt.overload_from_cold(m, BUS, cfg, 3000.0, 0.0, iron={k: v for k, v in th.iron_placeholder(3000.0).items() if k != "f_hz"}, j_rms=2.0 * j)
    assert mine["magnet_s"] == pytest.approx(ref["magnet_80C_s"], rel=1e-9)
    assert mine["winding_s"] == ref["winding_155C_s"] or mine["winding_s"] == pytest.approx(ref["winding_155C_s"], rel=1e-9)


# ---------------------------------------------------------------------------------------------- 5. what ends an envelope
def _fake_envelope(region, feasible, zero_rpm):
    n = len(region)
    on = np.array([r != "off" for r in region])
    d = {"rpm": 600.0 + 100.0 * np.arange(n), "T": np.where(on, 1.0, 0.0), "P_out": np.where(on, 100.0, 0.0), "eta": np.where(on, 0.9, 0.0), "region": np.array(region, dtype=object),
         "feasible": np.array(feasible, bool)}
    return mm.Envelope(d, {"corner_rpm": None, "corner_T": None, "zero_rpm": zero_rpm})


def test_an_envelope_ends_at_the_top_by_the_loads_or_by_the_limit():
    from optimization import motor_map_run as run
    top = run.summarise(_fake_envelope(["MTPA", "FW", "FW"], [1, 1, 1], 200.0))
    loads = run.summarise(_fake_envelope(["MTPA", "FW", "off"], [1, 1, 1], 150.0))        # a best point exists above the last speed with torque: the loads take all of its torque
    limit = run.summarise(_fake_envelope(["MTPA", "FW", "off"], [1, 1, 0], 150.0))        # no operating point above it: the limits end it
    assert (top["stop"], loads["stop"], limit["stop"]) == ("top", "loads", "limit")
    assert run.ends(limit).endswith("t") and not run.ends(loads).endswith("t") and not run.ends(top).endswith("t")


# ---------------------------------------------------------------------------------------------- 6. the winding study's windings
def test_the_optimized_windings_of_each_stack_and_the_headline_winding_against_the_study():
    table = st.OPTIMIZED_WINDINGS
    assert sorted(table) == [20.0, 40.0, 60.0, 80.0]
    assert st.stack_winding(table[60.0], 20.0, "opt") == table[20.0] and st.stack_winding(table[60.0], 80.0, "same") == table[60.0]
    with pytest.raises(ValueError):
        st.stack_winding(table[60.0], 30.0, "opt")
    w = table[60.0]
    assert (w.turns, w.paths, w.strands, w.wire_mm, w.wire_over_mm) == (9, 4, 7, 0.56, 0.63)
    m = mm.Machine(w, 60.0, copper=mm.FullCopper(arrangement="cluster"))
    assert m.fill == pytest.approx(0.42843, rel=1e-4)                       # the study's fill and R(100 C) for 9 turns, 4 paths, 7 x 0.56 mm
    assert 1e3 * float(m.resistance(100.0)) == pytest.approx(19.261, rel=1e-3)


def test_the_headline_winding_reproduces_the_studys_rating_and_envelope():
    asm = thermal.Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6)
    cool = thermal.size_shaft_fan(50.0, stack_mm=60.0, assume=asm)["cooling"]
    w = st.OPTIMIZED_WINDINGS[60.0]
    m = mm.Machine(w, 60.0, copper=mm.FullCopper(arrangement="cluster"))
    rate = mt.thermal_rating(m, mt.ThermalConfig(cool, asm), 3000.0)
    assert rate["AT_pk"] == pytest.approx(94.574, rel=2e-3) and rate["limited"] == "temperature" and rate["binding"] == "magnet"
    i = 94.574 * w.paths / w.turns                                          # the study's own rating as the current
    env = mm.envelope(m, mm.Bus(48.0, k=1.0), np.array([0.0, 3000.0, 6000.0, 10000.0]), 2.0 * i, basis="air_gap", scan=mm.Scan(m, i_max=2.0 * i), locate=True)
    assert np.asarray(env["T"]) == pytest.approx([8.78, 8.71, 5.07, 3.11], abs=0.0051)       # the study's short-term table (two decimals)
    assert env.meta["corner_rpm"] == pytest.approx(2866.0, abs=0.51)


def test_the_shaft_fan_is_sized_on_its_stack_with_the_run_assumptions():
    from optimization import motor_map_run as run
    for stack in (20.0, 60.0):
        c = run.fan_cooling("fanB", stack)
        p_top = thermal.shaft_fan_power(mm.TOP_RPM, c, stack, assume=run.thermal_assumptions())
        assert p_top == pytest.approx(50.0, rel=1e-6)                         # held to its budget at the top speed ON THAT STACK, with the diamond end turn
    assert run.fan_cooling("fanB", 20.0).v_ref_ms > run.fan_cooling("fanB", 60.0).v_ref_ms > 1.7


def test_the_energy_balance_closes_with_the_full_copper_loss_too():
    m = mm.Machine(st.OPTIMIZED_WINDINGS[60.0], 60.0, copper=mm.FullCopper(arrangement="cluster"))
    i_d, i_q, rpm = np.array([-5.0, -30.0, -20.0]), np.array([20.0, 40.0, 10.0]), np.array([1000.0, 10000.0, 16330.0])
    pt = m.point(i_d, i_q, rpm, BUS, 100.0)
    assert np.all(pt["P_cu_ac"] > 0.0)                                  # the AC, PWM and circulating part is not zero, and it is in the input power
    assert np.max(m.balance_error(pt)) < 1e-12
    assert np.all(pt["P_in"] == pytest.approx(1.5 * (pt["v_d"] * i_d + pt["v_q"] * i_q) + pt["P_cu_ac"], rel=1e-14))


def test_the_inputs_of_a_saved_run_come_back_as_they_were():
    from optimization import motor_map_run as run
    before = (run.WINDING, dict(run.STACK_WINDINGS), run.STACKS, run.HEADLINE_STACK, run.VDCS, run.K_MARGIN, run.COPPER, run.ARRANGEMENT, run.MAGNET_LIMIT_C, run.EXTRAS, run.EXTRAS_ON)
    rec = run.inputs_dict()
    run.restore_inputs(rec)
    after = (run.WINDING, dict(run.STACK_WINDINGS), run.STACKS, run.HEADLINE_STACK, run.VDCS, run.K_MARGIN, run.COPPER, run.ARRANGEMENT, run.MAGNET_LIMIT_C, run.EXTRAS, run.EXTRAS_ON)
    assert after == before
    assert rec["windings"]["opt"]["60"]["turns"] == 9 and rec["windings"]["opt"]["20"]["turns"] == 26 and [e["key"] for e in rec["extras"]] == ["f8", "f5", "f4", "f3", "g150"]
    assert rec["extras"][-1]["magnet_limit_c"] == 150.0 and rec["fan_v_ref"]["fanB"]["60"] == pytest.approx(1.76035, rel=1e-4)
