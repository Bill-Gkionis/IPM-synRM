"""Tests for the winding study (optimization/winding_study.py): everything that needs no FEMM.

  test_candidates_*          the enumeration holds exactly the windings the constraints allow (fill, the slot opening, the strands in hand, N, a), counted a second way
  test_resistance_*          the closed-form resistance against the copper-loss module's DC loss (1.5 I^2 R at the three phases)
  test_corner_speed_*        the closed-form corner speed solves |v| = V_max, goes with the bus and is zero where the resistive drop alone is over the limit
  test_hs_curve_*            the high-speed point's solver: its answer holds the demanded torque and the voltage limit, no point of a finer search has a smaller current, and the
                             motor-map engine's own solver agrees (skipped if that module cannot be loaded)
  test_rating_*              the continuous rating is the root of the thermal margin (a lower current is cool, a higher one is hot), a circulating-current layout rates lower than a compact
                             one, a pruned winding (fill under the DC rating's floor) cannot carry the common rating
  test_short_term_*          the short-term current is twice the continuous one, clipped at the maps' range
  test_front_*               the front of the headline case is the exact Pareto front of its allowed windings (independent loops and winding.verify_front) and holds the recommended winding
  test_pick_*                the choice of a winding for an inverter class and the recommendation agree with the allowed set
  test_rule_table_*          the voltage rule's margin table: the 10 % row is the recommendation, a larger margin needs more current
  test_arrangement_*         at the high-speed point random and line cost more than the compact arrangements
  test_manufacturing_note    strands abreast, room in the opening, stiffness
  test_engine_*              the torque-speed envelope of the motor-map engine against the closed-form corner speed (skipped if the engine cannot be loaded)
  test_front_figure_*        the front's figure is drawn (the other figures need the whole run's results: they are looked at, not tested)

The headline case (60 mm, 48 V, about twenty thousand windings) runs once for the whole module: 5 s.

Run from the Python to FEMM folder:   uv run pytest tests/test_winding_study.py -q
"""
import math

import numpy as np
import pytest

from optimization import copper_loss as cl
from optimization import maps, thermal, winding
from optimization import winding_study as ws


# ------------------------------------------------------------------------------------------------------ the shared run
@pytest.fixture(scope="module")
def case():
    return ws.Case()


@pytest.fixture(scope="module")
def ctx(case):
    return ws.get_context(case)


@pytest.fixture(scope="module")
def run(case):
    return ws.run_case(case)


def _index(cand, n, a, n_s, d):
    sel = np.flatnonzero((cand["N"] == n) & (cand["a"] == a) & (cand["n_s"] == n_s) & np.isclose(cand["d_mm"], d))
    assert len(sel) == 1, (n, a, n_s, d, len(sel))
    return int(sel[0])


# ------------------------------------------------------------------------------------------------------ the candidates
def test_candidates_respect_the_constraints():
    c = ws.candidates(8)
    assert c["fill"].max() <= ws.FILL_MAX + 1e-9
    assert (c["d_over_mm"] + 2.0 * ws.CLEARANCE_MM <= ws.SLOT_OPEN_MM + 1e-9).all()
    assert c["n_s"].min() == 1 and c["n_s"].max() == 8
    assert c["N"].min() >= 1 and c["N"].max() <= 60
    assert set(np.unique(c["a"]).tolist()) <= set(ws.PATHS)
    # the fill is the bare copper area over the layer's
    j = _index(c, 9, 4, 7, 0.56)
    assert c["fill"][j] == pytest.approx(9 * 7 * math.pi * 0.56 ** 2 / 4.0 / ws.layer_area_mm2(), rel=1e-12)


def test_candidates_count_matches_a_second_count():
    cat = winding.catalogue()
    area = ws.layer_area_mm2()
    wires = [(d, o) for d, o in zip(cat["d_mm"], cat["d_over_mm"]) if o + 2.0 * ws.CLEARANCE_MM <= ws.SLOT_OPEN_MM + 1e-9]
    count = 0
    for n in range(1, 61):
        for d, _ in wires:
            for n_s in range(1, 7):
                if n * n_s * math.pi * d * d / 4.0 / area <= ws.FILL_MAX + 1e-12:
                    count += len(ws.PATHS)
    assert count == len(ws.candidates(6)["N"])
    assert len(ws.candidates(4)["N"]) < len(ws.candidates(8)["N"]) < len(ws.candidates(30)["N"])


def test_the_wire_that_does_not_pass_the_opening_is_out():
    w = ws.wires()
    assert w["d_mm"].max() <= 1.0 + 1e-9
    assert 1.12 not in np.round(w["d_mm"], 3)


# ------------------------------------------------------------------------------------------------------ the resistance
@pytest.mark.parametrize("n,a,n_s,d,d_over", [(9, 4, 7, 0.56, 0.63), (26, 4, 8, 0.315, 0.367), (4, 2, 8, 0.8, 0.884)])
def test_resistance_matches_the_copper_loss_module(ctx, n, a, n_s, d, d_over):
    r = float(ws.series_resistance(ctx, n, a, n_s, d))
    i = 30.0
    lo = cl.losses({"N": n, "a": a, "n_s": n_s, "d_mm": d, "d_over_mm": d_over, "arrangement": "cluster"}, {"i_peak": i, "gamma_deg": 20.0, "rpm": 3000.0, "f_sw": 20000.0, "v_dc": 48.0}, 100.0, 60.0)
    assert float(lo["p_dc_slot"] + lo["p_dc_end"]) == pytest.approx(1.5 * i * i * r, rel=1e-5)


def test_resistance_goes_with_the_turns_and_the_copper(ctx):
    base = float(ws.series_resistance(ctx, 8, 4, 4, 0.8))
    assert float(ws.series_resistance(ctx, 16, 4, 4, 0.8)) == pytest.approx(2.0 * base, rel=1e-12)       # twice the turns
    assert float(ws.series_resistance(ctx, 8, 2, 4, 0.8)) == pytest.approx(4.0 * base, rel=1e-12)       # half the paths: twice the length in series, half the copper in parallel
    assert float(ws.series_resistance(ctx, 8, 4, 8, 0.8)) == pytest.approx(0.5 * base, rel=1e-12)       # twice the strands
    assert float(ws.series_resistance(ctx, 8, 4, 4, 0.8, 100.0)) > float(ws.series_resistance(ctx, 8, 4, 4, 0.8, 20.0))


# ------------------------------------------------------------------------------------------------------ the corner speed
def test_corner_speed_solves_the_voltage_equation():
    rng = np.random.default_rng(3)
    v_max = 48.0 / ws.SQRT3
    for _ in range(40):
        psi_d, psi_q = rng.uniform(0.002, 0.03), rng.uniform(0.002, 0.03)
        i_d, i_q = -rng.uniform(0.0, 60.0), rng.uniform(5.0, 80.0)
        r = rng.uniform(0.001, 0.02)
        rpm = float(ws.corner_speed(psi_d, psi_q, i_d, i_q, r, v_max))
        if r * math.hypot(i_d, i_q) >= v_max:
            assert rpm == 0.0
        else:
            assert float(ws.voltage(psi_d, psi_q, i_d, i_q, r, rpm)) == pytest.approx(v_max, rel=1e-9)


def test_corner_speed_goes_with_the_bus_without_resistance():
    a = float(ws.corner_speed(0.01, 0.02, -20.0, 40.0, 0.0, 48.0 / ws.SQRT3))
    b = float(ws.corner_speed(0.01, 0.02, -20.0, 40.0, 0.0, 60.0 / ws.SQRT3))
    assert b / a == pytest.approx(60.0 / 48.0, rel=1e-12)
    assert a == pytest.approx(60.0 * (48.0 / ws.SQRT3) / (2.0 * math.pi * ws.PP * math.hypot(0.01, 0.02)), rel=1e-12)        # |v| = w |psi|


def test_corner_speed_is_zero_when_the_resistive_drop_alone_is_over_the_limit():
    assert float(ws.corner_speed(0.01, 0.01, -50.0, 60.0, 0.5, 48.0 / ws.SQRT3)) == 0.0


# ------------------------------------------------------------------------------------------------------ the high-speed point
def _hs_state(ctx, c, j, sol):
    """The torque [N.m] and the voltage [V] at the solver's answer for candidate j, computed from the maps directly (placeholder scaling)."""
    s = c["N"][j] / (10.0 * c["a"][j])
    g = math.radians(float(sol["gamma"][0]))
    i_d, i_q = -float(sol["I"][0]) * math.sin(g), float(sol["I"][0]) * math.cos(g)
    pd, pq = ctx.fm.psi_dq(i_d * s, i_q * s, stack_mm=ctx.stack)
    psi_d, psi_q = s * pd, s * pq
    r = float(ws.series_resistance(ctx, c["N"][j], c["a"][j], c["n_s"][j], c["d_mm"][j]))
    return float(1.5 * ws.PP * (psi_d * i_q - psi_q * i_d)), float(ws.voltage(psi_d, psi_q, i_d, i_q, r, ws.HS_RPM))


@pytest.mark.parametrize("n,a,n_s,d", [(9, 4, 7, 0.56), (5, 4, 8, 0.71), (4, 4, 8, 0.8)])
def test_hs_curve_holds_the_torque_and_the_voltage(ctx, n, a, n_s, d):
    c = ws.candidates(8)
    j = _index(c, n, a, n_s, d)
    req = ws.requirement(ctx)
    curve = ws.hs_curve(ctx)
    s = c["N"][j] / (10.0 * c["a"][j])
    r = float(ws.series_resistance(ctx, c["N"][j], c["a"][j], c["n_s"][j], c["d_mm"][j]))
    i_st = 2.0 * req["AT"] * c["a"][j] / c["N"][j]
    sol = curve.solve(np.array([s]), np.array([r]), ws.HS_RPM, 48.0 / ws.SQRT3, np.array([i_st]))
    assert sol["feasible"][0]
    t, v = _hs_state(ctx, c, j, sol)
    assert t == pytest.approx(curve.target, rel=2e-3)
    assert v <= 48.0 / ws.SQRT3 * (1.0 + 1e-3)
    assert sol["I"][0] <= i_st + 1e-9


def test_hs_curve_has_no_better_point_in_a_finer_search(ctx):
    """An independent search: on a grid twice as fine in angle and radius, the smallest current on the demanded-torque contour that holds the voltage."""
    c = ws.candidates(8)
    req = ws.requirement(ctx)
    curve = ws.hs_curve(ctx)
    fm = ctx.fm
    gam = np.linspace(0.0, 90.0, 1441)
    i0 = np.linspace(0.0, min(-fm.range_d[0], fm.range_q), 960)
    g = np.radians(gam)[:, None]
    d0, q0 = -i0[None, :] * np.sin(g), i0[None, :] * np.cos(g)
    pd, pq = fm.psi_dq(d0, q0, stack_mm=ctx.stack)
    t0 = 1.5 * ws.PP * (pd * q0 - pq * d0)
    want = curve.target
    reach = (t0 >= want).any(axis=1)
    first = np.argmax(t0 >= want, axis=1)
    rows = np.flatnonzero(reach & (first > 0))
    f = (want - t0[rows, first[rows] - 1]) / (t0[rows, first[rows]] - t0[rows, first[rows] - 1])
    i_c = i0[first[rows] - 1] + f * (i0[first[rows]] - i0[first[rows] - 1])
    th = np.radians(gam[rows])
    pd_c, pq_c = fm.psi_dq(-i_c * np.sin(th), i_c * np.cos(th), stack_mm=ctx.stack)
    for n, a, n_s, d in ((9, 4, 7, 0.56), (6, 4, 4, 0.9), (4, 4, 8, 0.8)):
        j = _index(c, n, a, n_s, d)
        s = c["N"][j] / (10.0 * c["a"][j])
        r = float(ws.series_resistance(ctx, c["N"][j], c["a"][j], c["n_s"][j], c["d_mm"][j]))
        i_st = 2.0 * req["AT"] * c["a"][j] / c["N"][j]
        v = ws.voltage(s * pd_c, s * pq_c, -i_c * np.sin(th) / s, i_c * np.cos(th) / s, r, ws.HS_RPM)
        cur = i_c / s
        ok = (v <= 48.0 / ws.SQRT3) & (cur <= i_st)
        assert ok.any()
        sol = curve.solve(np.array([s]), np.array([r]), ws.HS_RPM, 48.0 / ws.SQRT3, np.array([i_st]))
        assert sol["I"][0] == pytest.approx(cur[ok].min(), rel=6e-3)          # the polar grid and the refinement at the voltage boundary: a few tenths of a percent


def test_hs_curve_infeasible_windings_carry_no_answer(ctx):
    curve = ws.hs_curve(ctx)
    sol = curve.solve(np.array([0.5]), np.array([0.5]), ws.HS_RPM, 48.0 / ws.SQRT3, np.array([5.0]))       # a huge resistance and a 5 A limit: cannot give the torque
    assert not sol["feasible"][0] and np.isnan(sol["I"][0])


def test_hs_solver_agrees_with_the_motor_map_engine(ctx, run):
    mm = pytest.importorskip("optimization.motor_map")
    c, o = run.cand, run.obj
    j = _index(c, 9, 4, 7, 0.56)
    try:
        w = mm.Winding(float(c["N"][j]), int(c["a"][j]), float(c["n_s"][j]), float(c["d_mm"][j]), float(c["d_over_mm"][j]))
        m = mm.Machine(w, 60.0, copper=mm.FullCopper(f_sw=20000.0, arrangement="cluster"))
        assert float(m.resistance(100.0)) == pytest.approx(float(o["R100"][j]), rel=1e-4)
        bus = mm.Bus(48.0, k=1.0)
        sc = mm.Scan(m, i_max=float(o["I_st"][j]))
        cx = sc.ctx(ws.HS_RPM, bus, 100.0, "air_gap")
        lv = sc.level_set(cx, float(o["I_st"][j]), np.array([o["T_hs"]]), objectives=("current",))["current"]
    except (AttributeError, TypeError, KeyError) as error:                     # the engine is another module's and in progress: an interface change is not this study's failure
        pytest.skip(f"motor_map's interface changed: {error}")
    assert lv["found"][0]
    assert lv["I"][0] == pytest.approx(float(o["I_hs"][j]), rel=1e-2)
    assert lv["gam"][0] == pytest.approx(float(o["gamma_hs"][j]), abs=0.5)


# ------------------------------------------------------------------------------------------------------ the rating
def test_rating_is_the_root_of_the_thermal_margin(ctx, case):
    w = dict(ws.REFERENCE)
    fill = ws.REFERENCE["N"] * ws.REFERENCE["n_s"] * math.pi * ws.REFERENCE["d_mm"] ** 2 / 4.0 / ctx.a_layer
    r = ws.continuous_rating(ctx, w, case.arrangement)
    assert r["feasible"] and not r["limited_by_maps"]
    m0, _, _ = ws.margin_at(ctx, w, case.arrangement, fill, r["AT"])
    assert abs(m0) < 0.02                                                      # the root
    assert ws.margin_at(ctx, w, case.arrangement, fill, 0.95 * r["AT"])[0] > 0.0
    assert ws.margin_at(ctx, w, case.arrangement, fill, 1.05 * r["AT"])[0] < 0.0
    assert r["binding"] == "magnet"                                           # the N45 at 80 C binds before the class F winding, in this cooling


def test_rating_falls_with_the_circulating_currents(ctx):
    w = {"N": 9, "a": 4, "n_s": 7, "d_mm": 0.56, "d_over_mm": 0.63}
    at = {arr: ws.continuous_rating(ctx, w, arr)["AT"] for arr in ("transposed", "cluster", "random", "line")}
    assert at["transposed"] >= at["cluster"] > at["random"] > at["line"]


def test_rating_of_the_requirement_is_a_share_of_the_reference(ctx, case):
    req = ws.requirement(ctx)
    assert req["AT"] == pytest.approx((1.0 - case.rating_margin) * req["AT_ref"], rel=1e-12)
    assert req["T"] < req["T_ref"] and req["T"] > 0.9 * req["T_ref"]


def test_pruned_windings_cannot_carry_the_common_rating(ctx, case):
    """The fill floor is a safe prune: windings just below it fail the network's check at the common rating."""
    req = ws.requirement(ctx)
    floor = ctx.fill_floor(req["AT"])
    c = ws.candidates(8)
    below = np.flatnonzero((c["fill"] < floor) & (c["fill"] > floor - 0.02))
    assert len(below) > 20
    pick = below[:: len(below) // 8][:8]
    rows = [(int(c["N"][k]), int(c["a"][k]), int(c["n_s"][k]), float(c["d_mm"][k]), float(c["d_over_mm"][k])) for k in pick]
    for rr in ws.check_rows({"case": case, "arrangement": "cluster", "rows": rows}):
        assert rr["margin"] < 0.0


def test_rating_rises_with_the_fill(ctx):
    fills = [0.1, 0.2, 0.3, 0.4, 0.45]
    ats = [ctx.at_guess(f) for f in fills]
    assert all(b > a for a, b in zip(ats, ats[1:]))


def test_thermal_check_runs_the_whole_enumeration(run):
    k, o = run.check, run.obj
    assert k["solved"].sum() < len(k["solved"])                                # some windings are pruned
    assert (k["margin"][~k["solved"]] < 0).all() and not o["thermal_ok"][~k["solved"]].any()
    assert o["thermal_ok"].sum() == ((k["margin"] >= 0) & ~k["runaway"]).sum()
    assert o["thermal_ok"].sum() > 100


# ------------------------------------------------------------------------------------------------------ the objectives
def test_objectives_are_consistent_with_the_requirement(run, ctx):
    req = ws.requirement(ctx)
    o, c = run.obj, run.cand
    idx = np.flatnonzero(run.allowed(True))
    assert len(idx) > 5
    assert np.allclose(o["I_cont"][idx], req["AT"] * c["a"][idx] / c["N"][idx], rtol=1e-12)
    assert np.allclose(o["I_st"][idx], 2.0 * o["I_cont"][idx], rtol=1e-12)               # not clipped at 60 mm on 48 V
    assert np.allclose(o["corner_bus"][idx], o["corner48"][idx], rtol=1e-12)
    assert (o["corner60"][idx] > o["corner48"][idx]).all()
    assert (o["P_hs"][idx] >= o["P_hs_dc"][idx]).all() and np.isfinite(o["P_hs"][idx]).all()
    assert (o["corner_bus"][idx] >= ws.RULE_RPM).all() and (c["fill"][idx] <= ws.FILL_MAX + 1e-9).all() and (c["n_s"][idx] <= 8).all()


def test_windings_of_one_n_over_a_share_the_current(run):
    """psi and i scale with N / a (maps.py's rule): 8 turns on 4 paths and 4 turns on 2 paths have the same short-term current, and nearly the same corner speed."""
    c, o = run.cand, run.obj
    j1, j2 = _index(c, 8, 4, 4, 0.8), _index(c, 4, 2, 8, 0.8)
    assert o["I_st"][j1] == pytest.approx(o["I_st"][j2], rel=1e-12)
    assert o["corner48"][j1] == pytest.approx(o["corner48"][j2], rel=2e-2)               # only the resistance differs (n_s a, N / a)


def test_short_term_current_is_clipped_at_the_maps_range():
    """With a high rating (fan (a), 150 C magnets: about 217 A-turns) twice it is past the maps' 3 x rated: the short-term ampere-turns stop there."""
    c2 = ws.Case(cooling="fan_a", magnet_limit_c=150.0)
    ctx2 = ws.get_context(c2)
    cand = ws.candidates(2)
    n = 4
    sub = {k: v[:n] for k, v in cand.items()}
    ob = ws.objectives(c2, sub, {"margin": np.ones(n), "runaway": np.zeros(n, bool)})
    assert (ob["AT"] * 2.0 > ctx2.at_max).all()
    assert np.allclose(ob["AT_st"], ctx2.at_max, rtol=1e-12)
    assert np.allclose(ob["I_st"], ctx2.at_max * sub["a"] / sub["N"], rtol=1e-12)


# ------------------------------------------------------------------------------------------------------ the front
def _pareto_by_loops(f):
    keep = []
    for i in range(len(f)):
        beaten = any((f[j] <= f[i]).all() and (f[j] < f[i]).any() for j in range(len(f)) if j != i)
        if not beaten:
            keep.append(i)
    return keep


def test_front_is_the_exact_pareto_front(run):
    o = run.obj
    idx = np.flatnonzero(run.allowed(True))
    f = winding.round_sig(np.stack([o["I_st"][idx], -o["corner_bus"][idx], o["P_hs"][idx]], axis=1), 9)
    expect = set(idx[_pareto_by_loops(f)].tolist())
    got = set(run.front(True).tolist())
    assert got == expect and 3 <= len(got) <= len(idx)
    verdict = winding.verify_front(np.stack([o["I_st"][idx], -o["corner_bus"][idx], o["P_hs"][idx]], axis=1), np.flatnonzero(np.isin(idx, run.front(True))))
    assert verdict["ok"]


def test_front_members_trade_current_against_speed_and_loss(run):
    o = run.obj
    fr = run.front(True)
    fr = fr[np.lexsort((o["corner_bus"][fr], np.round(o["I_st"][fr], 6)))]
    assert (np.diff(o["corner_bus"][fr]) >= -1e-6).all()                      # more current, more speed range
    assert o["P_hs"][fr][0] > o["P_hs"][fr][-1]                                # and less copper loss at the high-speed point


# ------------------------------------------------------------------------------------------------------ the choices
def test_pick_and_recommend_agree_with_the_allowed_set(run):
    o, c = run.obj, run.cand
    ok = np.flatnonzero(run.allowed(True))
    rec = ws.recommend(run)
    assert rec in set(ok.tolist())
    assert o["I_st"][rec] == pytest.approx(o["I_st"][ok].min(), rel=1e-9)
    assert rec in set(run.front(True).tolist()) or any(np.isclose(o["I_st"][j], o["I_st"][rec]) and np.isclose(o["P_hs"][j], o["P_hs"][rec], rtol=1e-3) for j in run.front(True))
    assert ws.pick_for_class(run, 0.99 * o["I_st"][ok].min()) is None            # below the smallest current that satisfies the rule: no winding
    top = ws.pick_for_class(run, 1e9)
    assert o["corner_bus"][top] >= (1 - 2e-3) * o["corner_bus"][ok].max()
    for lim in (100.0, 150.0, 200.0, 300.0):
        j = ws.pick_for_class(run, lim)
        assert j is not None and o["I_st"][j] <= lim + 1e-9 and j in set(ok.tolist())


def test_cap_only_narrows_the_set(run):
    a4, a8, a30 = (set(np.flatnonzero(run.allowed(True, cap)).tolist()) for cap in (4, 8, 30))
    assert a4 <= a8 <= a30 and len(a4) < len(a8) < len(a30)
    assert all(run.cand["n_s"][j] <= 4 for j in a4)


def test_the_rule_and_the_demand_narrow_the_set(run):
    all_th = np.flatnonzero(run.allowed(False, None, False))
    no_rule = np.flatnonzero(run.allowed(False, None, True))
    with_rule = np.flatnonzero(run.allowed(True, None, True))
    assert len(with_rule) <= len(no_rule) <= len(all_th)
    assert (run.obj["corner_bus"][with_rule] >= ws.RULE_RPM).all()


# ------------------------------------------------------------------------------------------------------ the arrangements
def test_arrangement_costs_at_the_high_speed_point(run):
    c, o = run.cand, run.obj
    j = _index(c, 9, 4, 7, 0.56)
    w = {"N": 9, "a": 4, "n_s": 7, "d_mm": 0.56, "d_over_mm": 0.63}
    p = {}
    for arr in ("transposed", "row", "cluster", "random", "line"):
        lh = cl.losses(dict(w, arrangement=arr), {"i_peak": float(o["I_hs"][j]), "gamma_deg": float(o["gamma_hs"][j]), "rpm": ws.HS_RPM, "f_sw": 20000.0, "v_dc": 48.0}, 100.0, 60.0)
        p[arr] = float(lh["p_total"])
    assert p["transposed"] < p["row"] <= p["cluster"] < p["random"] < p["line"]
    assert p["random"] / p["cluster"] > 1.1                                    # the layouts to avoid cost tenths of the loss
    assert p["cluster"] == pytest.approx(float(o["P_hs"][j]), rel=1e-9)         # the run's number is this arrangement's


def test_manufacturing_note():
    m = ws.manufacturing({"N": 9, "a": 4, "n_s": 7, "d_mm": 0.56, "d_over_mm": 0.63})
    assert m["strands_abreast"] == 1 and m["opening_room_mm"] == pytest.approx(0.47, abs=1e-9)
    assert m["stiffness_vs_0315"] == pytest.approx((0.56 / 0.315) ** 4, rel=1e-12)
    assert m["strands_per_coil_side"] == 63 and m["coils_per_path"] == 3
    thin = ws.manufacturing({"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367})
    assert thin["strands_abreast"] == 2 and thin["bundle_stiffness_vs_example"] == pytest.approx(1.0, rel=1e-12)


# ------------------------------------------------------------------------------------------------------ the case, the cooling, the engine
def test_case_is_hashable_and_labelled(case):
    assert hash(case) == hash(ws.Case())
    assert case != ws.Case(stack_mm=40.0)
    assert "60 mm" in case.label() and "48 V" in case.label()


def test_the_shaft_fan_is_held_to_its_power_at_the_top_speed(ctx, case):
    p = thermal.shaft_fan_power(ws.TOP_RPM, ctx.cooling, ctx.stack)
    assert p == pytest.approx(case.fan_budget_w, rel=0.03)
    assert thermal.shaft_fan_power(ws.RATED_RPM, ctx.cooling, ctx.stack) < 1.0


def test_engine_corner_speed_agrees_with_the_closed_form(run):
    pytest.importorskip("optimization.motor_map")
    c, o = run.cand, run.obj
    j = _index(c, 9, 4, 7, 0.56)
    try:
        sr = ws.speed_range(run.case, ws.winding_dict(c, j), run.arrangement, rpms=(0.0, 3000.0, 6000.0))
    except (AttributeError, TypeError, KeyError) as error:
        pytest.skip(f"motor_map's interface changed: {error}")
    if sr is None:
        pytest.skip("motor_map could not be loaded")
    d = ws.detail(run.case, ws.winding_dict(c, j), run.arrangement)
    assert sr["corner_cont"] == pytest.approx(d["corner48_rpm"], rel=1e-2)
    assert sr["T_cont"][0] == pytest.approx(d["T_cont"], rel=2e-2)
    assert sr["T_st"][0] > 1.9 * sr["T_cont"][0] * 0.95                       # the short-term current gives about twice the torque at standstill


def test_rule_table_ten_percent_is_the_recommendation(run):
    rows = ws.rule_table(ws.Report(), run)
    row10 = next(r for r in rows if abs(r["margin_pct"] - 10.0) < 1e-9)
    rec = ws.recommend(run)
    c = run.cand
    assert (row10["N"], row10["a"], row10["n_s"], row10["d_mm"]) == (c["N"][rec], c["a"][rec], c["n_s"][rec], c["d_mm"][rec])
    assert all(b["I_st"] >= a["I_st"] - 1e-9 for a, b in zip(rows, rows[1:]))          # a larger margin needs more inverter current
    assert all(r["corner"] >= r["rule_rpm"] for r in rows)


def test_front_figure_is_drawn(run, case, tmp_path):
    from optimization import winding_study_figures as wf
    path = tmp_path / "front.png"
    wf.fig_front({"run": run, "case": case, "rec_index": ws.recommend(run)}, path)
    assert path.exists() and path.stat().st_size > 20000
