"""Tests for the shaft fan of optimization/thermal.py and for optimization/thermal_machine.py (the network driven by the flux maps, the iron-loss map and the
copper-loss model).

Run from P:   uv run pytest tests/test_thermal_fan.py -q      (about 15 seconds; no solver and no FEMM are started)

What is held to what:
  1. Nothing that was there moves: the three earlier coolings give the same continuous current density, copper loss and hot spot as before the shaft fan was added
     (twelve numbers kept from the first study's output), and their external conductances equal the sum of the convection and radiation parts.
  2. The fan: the air speed follows the shaft speed; the power goes as the cube of it (the fan laws); a budget gives the speed that costs it; standstill is natural convection only.
  3. The heat flows: heat out equals heat made, and every node's flows add up to its source, in every cooling.
  4. The motor's models behind the network: the MTPA point and its iron loss, the DC copper loss of the network against the copper-loss model's, the continuous limit's binding
     margin, a limit that cannot be met (no current at all), and the flux maps' edge.
"""
import math

import numpy as np
import pytest

from optimization import thermal
from optimization.thermal import Assumptions, Cooling, Losses, ThermalModel

# the first study's output (sustainable.csv, class F, 3,000 rpm, the placeholder iron loss): J_rms [A/mm2], copper loss [W], hot spot [C]
GOLDEN = {("natural", 20): (3.01544189453125, 18.658130842065574, 84.6187808400261), ("natural", 40): (2.394866943359375, 14.554720444639681, 81.5615773297975),
          ("natural", 60): (1.847381591796875, 10.336007403924484, 80.24255639626173), ("natural", 80): (1.31805419921875, 6.11697338924269, 79.48870825484107),
          ("forced_air", 20): (6.8084716796875, 99.87970839001656, 104.77772626616809), ("forced_air", 40): (6.8902587890625, 124.88560569975954, 95.89212435858136),
          ("forced_air", 60): (6.883544921875, 147.87180286156033, 91.96925993041268), ("forced_air", 80): (6.8646240234375, 170.39216245821444, 90.54400468585145),
          ("liquid", 20): (11.536865234375, 315.32472969580806, 150.0372104412033), ("liquid", 40): (13.50341796875, 527.5358347115559, 139.1653368278582),
          ("liquid", 60): (14.443359375, 713.6041820297071, 132.8550767091743), ("liquid", 80): (14.98046875, 886.8896314444946, 132.95855832198953)}
COOLING = {"natural": Cooling.natural, "forced_air": Cooling.forced_air, "liquid": Cooling.liquid}


@pytest.mark.parametrize("key", sorted(GOLDEN))
def test_earlier_coolings_are_unchanged(key):
    kind, stack = key
    s = ThermalModel(stack_mm=float(stack), cooling=COOLING[kind]()).sustainable()
    j, p, hot = GOLDEN[key]
    assert s["J_rms"] == pytest.approx(j, rel=1e-9)
    assert s["P_cu"] == pytest.approx(p, rel=1e-9)
    assert s["hot_spot"] == pytest.approx(hot, rel=1e-9)


@pytest.mark.parametrize("kind", ["natural", "forced_air", "liquid", "shaft_fan"])
def test_outside_parts_add_up_to_external(kind):
    c = Cooling.shaft_fan(3.0, 3000.0) if kind == "shaft_fan" else COOLING[kind]()
    m = ThermalModel(stack_mm=60.0, cooling=c)
    th, tb = 75.0, 65.0
    p = m.outside_parts(th, tb)
    g_ha, g_ba, g_hc, _ = m.external(th, tb, 100.0)
    r_h2 = math.log(m.r_ho / m.r_hmid) / (2 * math.pi * thermal.MATERIALS["aluminium"]["k"] * m.L)
    g_h = p["h_conv"] + p["h_rad"]
    assert g_ba == pytest.approx(p["b_conv"] + p["b_rad"], rel=1e-12)
    assert g_ha == pytest.approx(1.0 / (r_h2 + 1.0 / g_h) if g_h > 0 else 0.0, rel=1e-12)
    if kind == "liquid":
        assert p["h_cool"] == pytest.approx(g_hc, rel=1e-12)
    xs = m.external_full(th, tb, 100.0)
    if kind != "shaft_fan":
        assert xs[4] == 0.0 and xs[5] == 0.0
    else:
        assert xs[4] > 0.0 and xs[5] > 0.0            # the fan's loss warms the stream: a little less heat leaves


def test_air_speed_follows_the_shaft_and_power_goes_as_the_cube():
    c = Cooling.shaft_fan(6.0, 3000.0)
    for rpm, v in ((0.0, 0.0), (1500.0, 3.0), (3000.0, 6.0), (6000.0, 12.0)):
        assert ThermalModel(stack_mm=60.0, cooling=c, speed_rpm=rpm).v_air == pytest.approx(v)
    m = ThermalModel(stack_mm=60.0, cooling=c, speed_rpm=3000.0)
    p = [m.fan_state(v)["p_shaft_w"] for v in (6.0, 12.0, 24.0)]
    assert 7.4 < p[1] / p[0] < 8.2 and 7.4 < p[2] / p[1] < 8.2           # twice the speed, 8 x the power, the friction factor bending it a little
    st = m.fan_state(6.0)
    assert st["p_shaft_w"] == pytest.approx(st["dp_pa"] * st["q_m3s"] / m.a.fan_efficiency, rel=1e-12)
    assert m.fan_state(0.0)["p_shaft_w"] == 0.0
    # a Cooling.shaft_fan(...) is a value: hashable, comparable
    assert Cooling.shaft_fan(6.0, 3000.0) == Cooling.shaft_fan(6.0, 3000.0) and hash(Cooling.shaft_fan(6.0)) == hash(Cooling.shaft_fan(6.0, 3000.0))


def test_shaft_fan_power_function_and_the_budget_sizing():
    c = Cooling.shaft_fan(6.0, 3000.0)
    arr = thermal.shaft_fan_power(np.array([0.0, 1500.0, 3000.0]), c, 60.0)
    assert arr.shape == (3,) and arr[0] == 0.0 and arr[1] < arr[2]
    assert thermal.shaft_fan_power(3000.0, c, 60.0) == pytest.approx(arr[2], rel=1e-12)
    for budget in (10.0, 50.0, 400.0):
        s = thermal.size_shaft_fan(budget, stack_mm=60.0)
        assert thermal.shaft_fan_power(16330.0, s["cooling"], 60.0) == pytest.approx(budget, rel=1e-6)
        assert s["v_ref"] == pytest.approx(s["v_top"] * 3000.0 / 16330.0, rel=1e-12)
        assert s["cooling"].kind == "shaft_fan"
    with pytest.raises(ValueError):
        thermal.shaft_fan_power(3000.0, Cooling.forced_air(), 60.0)


def test_standstill_is_natural_convection_of_the_finned_housing():
    m = ThermalModel(stack_mm=60.0, cooling=Cooling.shaft_fan(6.0, 3000.0), speed_rpm=0.0)
    assert m.v_air == 0.0 and m.g_forced_h == 0.0 and m.g_forced_b == 0.0 and m.dT_fan == 0.0
    plain = ThermalModel(stack_mm=60.0, cooling=Cooling.natural(), speed_rpm=0.0)
    p_fan, p_plain = m.outside_parts(80.0, 80.0), plain.outside_parts(80.0, 80.0)
    a = m.a
    assert p_fan["h_conv"] == pytest.approx(p_plain["h_conv"] * a.fan_fin_gain_conv, rel=1e-12)
    assert p_fan["h_rad"] == pytest.approx(p_plain["h_rad"] * a.fan_fin_gain_rad, rel=1e-12)


def test_mixed_convection_lies_between_the_blend_and_its_limits():
    m = ThermalModel(stack_mm=60.0, cooling=Cooling.shaft_fan(6.0, 3000.0), speed_rpm=3000.0)
    p = m._fan_parts(80.0, 80.0)
    f, n = m.g_forced_h, p["h_natural"]
    assert max(f, n) <= p["h_conv"] <= f + n                           # a blend of exponent 3 is between the larger and the sum
    off = ThermalModel(stack_mm=60.0, cooling=Cooling.shaft_fan(6.0, 3000.0), speed_rpm=3000.0, assume=Assumptions(air_bulk_warming=False))
    assert off.g_forced_h > m.g_forced_h                                # the air's warming takes conductance away
    assert m.g_forced_h + m.g_forced_b < m.mdot_cp_air                  # the stream cannot carry more than m cp per kelvin


@pytest.mark.parametrize("kind", ["natural", "forced_air", "liquid", "shaft_fan"])
def test_heat_flows_balance_node_by_node(kind):
    c = Cooling.shaft_fan(3.0, 3000.0) if kind == "shaft_fan" else COOLING[kind]()
    m = ThermalModel(stack_mm=40.0, cooling=c, speed_rpm=3000.0)
    res = m.steady_state(Losses(copper=40.0, iron_teeth=6.0, iron_yoke=4.0, iron_rotor=1.5, magnet=0.8), tol=1e-10, max_iter=200)
    hf = m.heat_flows(res)
    assert abs(hf["balance_w"]) < 1e-6
    assert sum(hf["exits"].values()) == pytest.approx(hf["total_made"], abs=1e-6)
    L, X, g, src = hf["links"], hf["exits"], hf["groups"], hf["made"]
    wind_half = 0.5 * src["windage"]
    # node by node: what leaves = what is made (links are from>to; E>W is the end turns' heat into the winding)
    assert L["W>T"] + L["W>Y"] - L["E>W"] == pytest.approx(src["copper_slot"], abs=1e-6)                                     # W
    assert L["E>W"] + L["E>A"] + L["E>B"] == pytest.approx(src["copper_end"], abs=1e-6)                                      # E
    assert L["T>Y"] + L["T>A"] + L["T>R"] - L["W>T"] == pytest.approx(src["iron_teeth"] + wind_half, abs=1e-6)               # T
    assert L["Y>H"] - L["T>Y"] - L["W>Y"] == pytest.approx(src["iron_yoke"], abs=1e-6)                                       # Y
    assert L["R>M"] + L["R>S"] - L["T>R"] - L["A>R"] == pytest.approx(src["iron_rotor"] + wind_half, abs=1e-6)               # R
    assert -L["R>M"] == pytest.approx(src["magnet"], abs=1e-6)                                                               # M
    assert L["S>G"] - L["R>S"] - L["A>S"] == pytest.approx(0.0, abs=1e-6)                                                    # S
    assert L["G>B"] - L["S>G"] == pytest.approx(src["bearing"], abs=1e-6)                                                    # G
    assert g["rotor_made"] == pytest.approx(g["rotor_across_gap"] + g["rotor_to_end_space_air"] + g["rotor_through_shaft"], abs=1e-6)
    # the exits take what reaches the housing and the bells
    h_in = L["Y>H"] + L["B>H"]                                  # the yoke's heat, and the bells' (B>H negative: the housing heats the bells)
    b_out = L["A>B"] + L["E>B"] + L["G>B"] - L["B>H"]           # what reaches the bells through the air, the end turns and the bearings, less what goes to the housing
    assert X["housing_convection"] + X["housing_radiation"] + X["housing_to_coolant"] == pytest.approx(h_in, abs=1e-6)
    assert X["bells_convection"] + X["bells_radiation"] == pytest.approx(b_out, abs=1e-6)
    assert X["housing_radiation"] >= 0.0 and X["bells_radiation"] >= 0.0


def test_a_function_of_the_current_density_for_the_iron_loss_equals_the_dict_when_constant():
    m = ThermalModel(stack_mm=60.0, cooling=Cooling.forced_air())
    ir = m.iron_placeholder()
    a = m.sustainable(iron=ir)
    b = m.sustainable(iron=lambda j: ir)
    c = m.sustainable(losses_at=lambda j: Losses(copper=m.copper(j), iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"]))
    assert a["J_rms"] == pytest.approx(b["J_rms"], rel=1e-12) and a["J_rms"] == pytest.approx(c["J_rms"], rel=1e-12)
    # an iron loss that grows with the current lowers the limit
    d = m.sustainable(iron=lambda j: {k: v * (1.0 + 0.2 * j) for k, v in ir.items() if k != "f_hz"})
    assert d["J_rms"] < a["J_rms"]
    # j_hi: a limit beyond the largest current searched is reported at j_hi
    e = m.sustainable(j_hi=2.0)
    assert e["J_rms"] == pytest.approx(2.0) and e["feasible"]


def test_overload_takes_the_iron_loss_of_the_overload_and_of_the_state_before():
    m = ThermalModel(stack_mm=60.0, cooling=Cooling.forced_air())
    ir = m.iron_placeholder()
    hot_iron = {k: 3.0 * v for k, v in ir.items() if k != "f_hz"}
    j = 6.0
    base = m.overload(2.0, j, start="cold", iron=ir)
    more = m.overload(2.0, j, start="cold", iron=hot_iron)
    assert more["winding_155C_s"] < base["winding_155C_s"]
    start_a = m.overload(2.0, j, start="hot", iron=ir, iron_start=ir)
    start_b = m.overload(2.0, j, start="hot", iron=ir, iron_start=hot_iron)
    assert start_b["hot_spot_start"] > start_a["hot_spot_start"]


# ------------------------------------------------------------------------------------------------ the motor's models
@pytest.fixture(scope="module")
def mach():
    pytest.importorskip("optimization.maps")
    from optimization import thermal_machine as tm
    try:
        return tm.Machine()
    except FileNotFoundError:
        pytest.skip("the flux maps or the iron-loss map are not in out/optimization/maps/chosen")


def test_machine_point_and_the_units(mach):
    from optimization import maps
    assert mach.at_rated == pytest.approx(10.0 * mach.fm.rated) and mach.at_max == pytest.approx(3.0 * mach.at_rated)
    pt = mach.point(mach.at_rated, 3000.0, 20.0)
    g, i_d, i_q, t = mach.fm.mtpa(mach.fm.rated, stack_mm=20.0)
    assert pt["torque"] == pytest.approx(t) and pt["gamma_deg"] == pytest.approx(g) and i_d == pytest.approx(pt["i_d"])
    assert pt["iron_teeth"] > 0.0 and not pt["beyond_maps"]
    assert mach.point(0.0, 0.0, 60.0)["iron_teeth"] == 0.0                         # standstill: no iron loss
    with pytest.raises(ValueError):
        mach.point(mach.at_max * 1.01, 3000.0, 60.0)
    edge = mach.point(mach.at_max * 1.5, 3000.0, 60.0, clip=True)
    assert edge["beyond_maps"] and edge["torque"] is None
    # the torque scales with the stack
    assert mach.point(100.0, 3000.0, 60.0)["torque"] == pytest.approx(3.0 * mach.point(100.0, 3000.0, 20.0)["torque"], rel=1e-9)
    assert maps.winding_scale(26, 4) == pytest.approx(0.65)


def test_dc_copper_of_the_network_is_the_copper_loss_models(mach):
    from optimization import thermal_machine as tm
    m = tm.model_for(Cooling.forced_air(), 60.0, 3000.0, mach)
    at = 118.6
    pt = mach.point(at, 3000.0, 60.0)
    f, share, r100 = mach.copper(at, pt["gamma_deg"], 3000.0, 60.0)
    j = at / (math.sqrt(2.0) * m.ampere_turns(1.0)["a_cu_side_mm2"])
    assert m.copper(j).total(100.0) == pytest.approx(r100["p_dc_slot"] + r100["p_dc_end"], rel=1e-5)
    assert 0.5 < share < 0.6 and f(100.0) == pytest.approx(r100["p_total"], rel=1e-12)
    assert f(20.0) < f(100.0)                                                       # the copper's resistivity follows its temperature


def test_continuous_limit_margin_is_zero_at_the_binding_limit(mach):
    from optimization import thermal_machine as tm
    for cooling, rpm in ((Cooling.forced_air(), 3000.0), (Cooling.shaft_fan(6.0, 3000.0), 1000.0), (Cooling.shaft_fan(6.0, 3000.0), 10000.0)):
        m = tm.model_for(cooling, 60.0, rpm, mach)
        cl = tm.continuous_limit(m, mach)
        assert cl["feasible"] and not cl["limited_by_maps"]
        margin = cl["margin_winding"] if cl["binding"] == "winding" else cl["margin_magnet"]
        assert abs(margin) < 0.02, (cooling.label, rpm, margin)                      # the bisection's tolerance on J, a few hundredths of a kelvin
        assert min(cl["margin_winding"], cl["margin_magnet"]) > -0.02
        assert cl["torque"] > 0 and 0 < cl["efficiency"] < 1 and cl["p_em_w"] == pytest.approx(cl["torque"] * 2 * math.pi * rpm / 60.0)
        # more current than the limit passes a limit
        over = m.steady_state(tm.losses_at_factory(m, mach)(cl["J_rms"] * 1.02))
        assert min(over.margin_winding, over.margin_magnet) < 0.0


def test_a_limit_that_cannot_be_met_and_the_maps_edge(mach):
    from optimization import thermal_machine as tm
    hot = tm.continuous_limit(tm.model_for(Cooling.natural(), 60.0, 10000.0, mach), mach)
    assert not hot["feasible"] and hot["J_rms"] == 0.0 and hot["AT_pk"] == 0.0 and hot["torque"] == 0.0
    assert min(hot["margin_winding"], hot["margin_magnet"]) < 0.0                    # the iron, windage and bearing heat alone passes the magnets' limit
    liq = tm.continuous_limit(tm.model_for(Cooling.liquid(), 60.0, 0.0, mach), mach)
    assert liq["limited_by_maps"] and liq["AT_pk"] == pytest.approx(mach.at_max, rel=1e-6) and liq["binding"] == "maps"


def test_short_term_torque_and_time(mach):
    from optimization import thermal_machine as tm
    m = tm.model_for(Cooling.forced_air(), 60.0, 3000.0, mach)
    cl = tm.continuous_limit(m, mach)
    sh = tm.short_term(m, mach, cl, 2.0)
    assert sh["torque"] > cl["torque"] and sh["AT_pk"] == pytest.approx(2.0 * cl["AT_pk"])
    assert sh["t_first_s"] is not None and 0.0 < sh["t_first_s"] < 4 * 3600.0 and sh["first"] in ("winding", "magnet")
    liq = tm.model_for(Cooling.liquid(), 60.0, 3000.0, mach)
    sh2 = tm.short_term(liq, mach, tm.continuous_limit(liq, mach), 2.0)
    assert sh2["torque"] is None and sh2["beyond_maps"] and sh2["t_first_s"] is not None               # 2x lies beyond the maps: no torque, a time
    # a longer overload from cold takes longer at a lower factor
    assert tm.short_term(m, mach, cl, 1.5)["t_first_s"] > sh["t_first_s"]


def test_the_shaft_fan_cools_least_at_low_speed(mach):
    from optimization import thermal_machine as tm
    sa = Cooling.shaft_fan(6.0, 3000.0)
    at = {rpm: tm.continuous_limit(tm.model_for(sa, 60.0, float(rpm), mach), mach)["AT_pk"] for rpm in (0, 3000, 6000)}
    fixed = tm.continuous_limit(tm.model_for(Cooling.forced_air(), 60.0, 0.0, mach), mach)["AT_pk"]
    assert at[0] < 0.7 * fixed                          # at standstill the fixed-speed fan keeps blowing, the shaft fan has stopped
    assert at[3000] > 1.3 * at[0] and at[6000] > at[0]


def test_operating_point_at_any_angle(mach):
    from optimization import thermal_machine as tm
    m = tm.model_for(Cooling.shaft_fan(6.0, 3000.0), 60.0, 6000.0, mach)
    cl = tm.continuous_limit(m, mach)
    res = tm.operating_point(m, mach, cl["AT_pk"])                                   # the MTPA point at the continuous limit: the binding margin is zero
    assert abs(min(res.margin_winding, res.margin_magnet)) < 0.02
    # a field-weakened point (a larger angle at the same ampere-turns): less torque, the iron loss of that point, another temperature
    pt_mtpa, pt_fw = mach.point(cl["AT_pk"], 6000.0, 60.0), mach.point(cl["AT_pk"], 6000.0, 60.0, gamma_deg=60.0)
    assert pt_fw["torque"] < pt_mtpa["torque"] and pt_fw["gamma_deg"] == 60.0
    assert pt_fw["iron_teeth"] != pt_mtpa["iron_teeth"]
    assert mach.point(cl["AT_pk"], 6000.0, 60.0, gamma_deg=pt_mtpa["gamma_deg"])["torque"] == pytest.approx(pt_mtpa["torque"], rel=2e-3)   # the d-q torque at the MTPA angle
    res_fw = tm.operating_point(m, mach, cl["AT_pk"], gamma_deg=60.0)
    assert res_fw.magnet != pytest.approx(res.magnet, abs=0.05)
    assert m.heat_flows(res_fw)["balance_w"] == pytest.approx(0.0, abs=1e-3)


def test_the_fans_own_heat_warms_the_air_and_can_be_switched_off():
    c = Cooling.shaft_fan(6.0, 3000.0)
    on = ThermalModel(stack_mm=60.0, cooling=c, speed_rpm=10000.0)
    off = ThermalModel(stack_mm=60.0, cooling=c, speed_rpm=10000.0, assume=Assumptions(fan_air_preheat=False))
    assert on.dT_fan == pytest.approx(on.fan_state()["dT_air_k"]) and on.dT_fan > 3.0 and off.dT_fan == 0.0
    assert off.external_full(80.0, 70.0, 0.0)[4:] == (0.0, 0.0)
    r_on, r_off = (m.steady_state(Losses(copper=100.0, iron_teeth=20.0)) for m in (on, off))
    assert r_on.T["H"] > r_off.T["H"] + 1.0 and r_on.hot_spot > r_off.hot_spot            # a warmer stream: a hotter housing and winding
    # the fan's power itself is not a source: the heat made is the same
    assert sum(r_on.sources.values()) == pytest.approx(sum(r_off.sources.values()), rel=1e-12)
