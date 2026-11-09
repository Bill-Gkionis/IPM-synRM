"""Tests for B7 (bench/b7_salib.py, and bench/femm_points.py, the FEMM side it shares with B8): Python -> SALib -> FEMM. Everything that needs no FEMM: the fixed choices of the
brief and of the rulings, the four parameters and their ranges against motorgeo's checks, the Morris sample and its seed, the SALib wiring on a known answer, the
exact comparison of two runs, the ranking, the figure, the checks.

Run from P:   uv run pytest tests/test_bench_b7.py -q      (about 15 s: the grid of 256 designs is built and checked by motorgeo; FEMM is not started, no slot is taken)

In groups:
  1. The fixed choices: gamma, the current, the levels, the seed, the trajectories, the solves per run, the four parameters.
  2. The ranges: around the design, and every level of the grid inside motorgeo's checks (the widest rung of the ladder that passes, and the next one up does not).
  3. The Morris sample: deterministic with its seed, on the grid, one parameter moved at a time; the columns are the named parameters.
  4. SALib on a known answer (the wiring), the exact comparison of two runs, the one-at-a-time swings, the ranking and its physics.
  5. The picture (two Morris views, no blank margin) and the checks.
  6. The FEMM side (femm_points): a design that motorgeo rejects never reaches FEMM; the worker can be sent to a pool by name.
  7. What B7 and B8 share (fix round 1 of Task 13): the rotor positions and the task limit written once, B1's numbers read once from the active root (and why they are missing, if
     they are), the pool record's limit given as a parameter; SALib's defaults with the lines of its installed source that say them, read back; the spread of each index.
"""
import ast
import copy
import json
import pickle

import numpy as np
import pytest

import motorgeo as mg
from bench_test_tools import P, use_config


# ----------------------------------------------------------------------------------------------- 1. the fixed choices
def test_the_current_angle_is_b1s_mtpa_angle_on_the_fine_mesh_and_says_where_it_comes_from():
    """A design decision: gamma = 22.9 degrees, B1's MTPA angle on the fine mesh, typed as a constant whose source names B1's results.json of 11-06. At gamma = 0 only the
    magnet torque counts, which would hide the saliency parameters; at the MTPA angle both kinds count."""
    from bench import b7_salib as b7
    assert b7.GAMMA_DEG == 22.9 and b7.AMPS == 10.0
    assert "out\\bench\\B1\\results.json" in b7.GAMMA_SOURCE and "11-06" in b7.GAMMA_SOURCE and "22.9" in b7.GAMMA_SOURCE and "fine mesh" in b7.GAMMA_SOURCE
    point = b7.current_point(30.0)
    assert point["theta_e"] == 30.0
    assert point["id"] == pytest.approx(-10.0 * np.sin(np.radians(22.9))) and point["iq"] == pytest.approx(10.0 * np.cos(np.radians(22.9)))
    assert np.hypot(point["id"], point["iq"]) == pytest.approx(10.0)


def test_the_size_of_the_run_is_stated_and_within_the_sixty_solves_of_the_ruling():
    """r x (k + 1) samples, each solved at the rotor positions of POSITIONS_EL: at most about 60 FEMM solves per run, on the coarse mesh."""
    from bench import b7_salib as b7
    assert b7.LEVELS == 4 and b7.SEED == 1 and b7.TRAJECTORIES == 6 and len(b7.PARAMETERS) == 4
    assert b7.samples() == b7.TRAJECTORIES * (len(b7.PARAMETERS) + 1) == 30
    assert b7.POSITIONS_EL == (0.0, 30.0) and b7.solves_per_run() == 30 * 2 == 60 <= b7.SOLVE_BUDGET == 60


def test_the_four_parameters_are_rotor_dimensions_and_each_changes_the_band_free_drawing():
    """Four rotor dimensions of motorgeo.Params. With barrier="none" some of its rotor fields draw nothing (rib, bridge, web, bar_t, bar_apex belong to the band or only
    guard the build): a Morris study on them would have an index of exactly 0 by construction, which says nothing about the motor. Each chosen one moves the drawing."""
    from bench import b7_salib as b7
    base = mg.build(mg.Params(barrier="none"))
    assert set(b7.PARAMETERS) <= set(mg.Params.__dataclass_fields__) and len(set(b7.PARAMETERS)) == 4
    assert not set(b7.PARAMETERS) & {"airgap", "stack", "slots", "pole_pairs", "r_so", "r_si", "turns"}                  # rotor dimensions, not the stator's or the winding's

    def drawing(**over):
        regions = mg.build(mg.Params(barrier="none", **over))
        return [(r.name, [tuple(np.round(pt, 9)) for e in r.edges for pt in e.points()]) for r in regions]
    for name in b7.PARAMETERS:
        moved = drawing(**{name: getattr(mg.Params(), name) * (1 + b7.RANGE_PERCENT / 100)})
        assert moved != drawing(), f"{name} does not change the drawing"
    for inert in ("rib", "bridge", "web", "bar_t", "bar_apex"):                                                            # the control: these draw nothing without the band
        assert drawing(**{inert: getattr(mg.Params(), inert) * 1.02}) == drawing(), inert
    assert [r.name for r in base]                                                                                          # (the baseline builds)


def test_every_parameter_has_a_physics_line_in_one_line_and_the_ranking_prints_it():
    from bench import b7_salib as b7
    assert set(b7.PHYSICS) == set(b7.PARAMETERS)
    for name, line in b7.PHYSICS.items():
        assert "\n" not in line and 40 < len(line) < 400, name


# ----------------------------------------------------------------------------------------------- 2. the ranges
def test_the_bounds_are_the_design_value_plus_and_minus_the_range_and_the_design_is_inside():
    from bench import b7_salib as b7
    problem = b7.problem()
    base = mg.Params(barrier="none")
    assert problem["num_vars"] == 4 and problem["names"] == list(b7.PARAMETERS)
    for name, (low, high) in zip(problem["names"], problem["bounds"]):
        value = getattr(base, name)
        assert low == pytest.approx(value * (1 - b7.RANGE_PERCENT / 100)) and high == pytest.approx(value * (1 + b7.RANGE_PERCENT / 100)) and low < value < high
    assert b7.RANGE_PERCENT == 1.5 and 1.5 in b7.RANGE_LADDER_PERCENT and b7.RANGE_LADDER_PERCENT == tuple(sorted(b7.RANGE_LADDER_PERCENT, reverse=True))


def test_every_level_of_the_morris_grid_passes_motorgeos_checks():
    """The ruling: 'check every corner of the box and every level before any FEMM run'. 4 levels in 4 parameters is 256 designs, the 16 corners among them: none is rejected,
    neither by build() (GeometryError) nor by check() (a violation)."""
    from bench import b7_salib as b7
    grid = b7.grid_overrides()
    assert len(grid) == 256 == b7.LEVELS ** 4
    screened = b7.screen(grid)
    assert screened["points"] == 256 and screened["count"] == 0 and screened["rejected"] == []
    corners = b7.grid_overrides(corners_only=True)
    assert len(corners) == 16 and all(c in grid for c in corners)


def test_the_chosen_range_is_the_widest_rung_of_the_ladder_that_passes():
    """The next rung up has designs that motorgeo rejects, so the chosen range cannot be widened without breaking a trajectory. (The design sits 0.158 mm from the build guard
    of the outer bridge: three of the four dimensions push the slot's outer end toward the rotor surface, and together at their corners they use it up. The corners are where
    it breaks first.) Checked on the 16 corners of the next rung: one at least is rejected."""
    from bench import b7_salib as b7
    ladder = list(b7.RANGE_LADDER_PERCENT)
    above = ladder[ladder.index(b7.RANGE_PERCENT) - 1]
    assert above > b7.RANGE_PERCENT
    screened = b7.screen(b7.grid_overrides(percent=above, corners_only=True))
    assert screened["count"] >= 1 and all(isinstance(reason, str) and reason for _, reason in screened["rejected"])


def test_a_design_that_motorgeo_rejects_is_named_with_its_reason():
    from bench import b7_salib as b7
    assert b7.rejection({}) is None
    assert "reaches the bridge" in b7.rejection({"slot_w": 9.28 * 1.10})                                                   # build() refuses (GeometryError)
    assert "centre_rib" in b7.rejection({"mag_y": 4.8})                                                                    # check() reports a violation
    screened = b7.screen([{}, {"slot_w": 12.0}])
    assert screened["points"] == 2 and screened["count"] == 1 and screened["rejected"][0][0] == {"slot_w": 12.0}


# ----------------------------------------------------------------------------------------------- 3. the Morris sample
def test_the_sample_is_the_same_every_time_for_one_seed_and_another_for_another():
    from bench import b7_salib as b7
    first, second = b7.sample(), b7.sample()
    assert first.shape == (30, 4) and np.array_equal(first, second)                                                        # exactly, not nearly
    assert not np.array_equal(first, b7.sample(seed=b7.SEED + 1))


def test_the_sample_is_on_the_four_level_grid_inside_the_bounds_and_moves_one_parameter_at_a_time():
    from bench import b7_salib as b7
    problem = b7.problem()
    low, high = np.array([b[0] for b in problem["bounds"]]), np.array([b[1] for b in problem["bounds"]])
    unit = (b7.sample() - low) / (high - low)
    assert unit.min() >= -1e-12 and unit.max() <= 1 + 1e-12
    assert np.allclose(unit * 3, np.round(unit * 3), atol=1e-9)                                                            # levels 0, 1/3, 2/3, 1
    for t in range(b7.TRAJECTORIES):                                                                                       # a trajectory: k + 1 points, each one step from the last in one parameter
        block = unit[t * 5:(t + 1) * 5]
        moved = [int(np.count_nonzero(~np.isclose(a, b))) for a, b in zip(block[:-1], block[1:])]
        assert moved == [1, 1, 1, 1], (t, moved)
        assert sorted(int(np.argmax(~np.isclose(a, b))) for a, b in zip(block[:-1], block[1:])) == [0, 1, 2, 3]            # every parameter once


def test_the_columns_of_the_sample_are_the_named_parameters_in_the_order_of_the_problem():
    from bench import b7_salib as b7
    row = b7.sample()[3]
    over = b7.overrides_of(row)
    assert list(over) == list(b7.PARAMETERS) and all(isinstance(v, float) for v in over.values())
    assert [over[n] for n in b7.PARAMETERS] == [float(v) for v in row]
    assert mg.Params(barrier="none", **over).barrier == "none"


def test_no_sample_of_this_run_is_rejected_by_motorgeo():
    from bench import b7_salib as b7
    screened = b7.screen([b7.overrides_of(row) for row in b7.sample()])
    assert screened["points"] == 30 and screened["count"] == 0


# ----------------------------------------------------------------------------------------------- 4. SALib on a known answer, two runs, the swings, the ranking
def test_salib_reproduces_a_known_linear_answer_on_this_very_sample():
    """The wiring of sample -> evaluate -> analyze: for a linear function every elementary effect is the coefficient times the range, so mu* = |c| x range, mu = c x range
    (signed) and sigma = 0, whichever points the sample holds. A mixed-up column, or a range taken the wrong way, shows here and not in the motor's numbers."""
    from bench import b7_salib as b7
    X = b7.sample()
    coef = np.array([4.0, -2.5, 0.0, 10.0])
    ranges = np.array([hi - lo for lo, hi in b7.problem()["bounds"]])
    result = b7.analyse(X, X @ coef)
    assert result["names"] == list(b7.PARAMETERS)
    assert np.allclose(result["mu_star"], np.abs(coef) * ranges, rtol=1e-9, atol=1e-12)
    assert np.allclose(result["mu"], coef * ranges, rtol=1e-9, atol=1e-12)
    assert np.allclose(result["sigma"], 0.0, atol=1e-9) and np.allclose(result["mu_star_conf"], 0.0, atol=1e-9)
    assert b7.ranking(result) == [b7.PARAMETERS[i] for i in np.argsort(-np.abs(coef) * ranges, kind="stable")]


def test_two_runs_must_agree_exactly_and_one_unit_in_the_last_place_is_a_difference():
    """'Two runs give identical indices': identical, not close. FEMM 5.0 was deterministic in B1 to 1e-13; if it is not exactly so here, the check FAILS with the difference
    and is never loosened into a tolerance. The comparison is on every array of the analysis and on the torques."""
    from bench import b7_salib as b7
    X = b7.sample()
    Y = np.linspace(1.0, 1.2, len(X))
    a, b = b7.analyse(X, Y), b7.analyse(X, Y.copy())
    assert b7.exactly_equal(a, b) and b7.exactly_equal(Y, Y.copy())
    Z = Y.copy()
    Z[7] = np.nextafter(Z[7], 2.0)                                                                                         # one unit in the last place
    assert not b7.exactly_equal(Y, Z)
    c = dict(a)
    c["sigma"] = [np.nextafter(v, 1.0) for v in a["sigma"]]
    assert not b7.exactly_equal(a, c) and not b7.exactly_equal(a, {k: v for k, v in a.items() if k != "sigma"})
    assert not b7.exactly_equal(np.array([1.0, 2.0]), np.array([1.0, 2.0, 3.0]))


def test_the_difference_between_two_runs_is_measured_so_that_a_fail_can_say_how_large():
    from bench import b7_salib as b7
    X = b7.sample()
    Y = np.linspace(1.0, 1.2, len(X))
    Z = Y.copy()
    Z[3] += 3e-13
    assert b7.largest_difference(Y, Z) == pytest.approx(3e-13, rel=1e-3) and b7.largest_difference(Y, Y.copy()) == 0.0
    assert b7.largest_difference(b7.analyse(X, Y), b7.analyse(X, Y.copy())) == 0.0


def test_the_one_at_a_time_designs_are_the_centre_and_each_end_of_each_range():
    from bench import b7_salib as b7
    designs = b7.oat_designs()
    assert len(designs) == 1 + 2 * 4 and designs[0] == {"name": "centre", "param": None, "side": None, "overrides": {}}
    base = mg.Params(barrier="none")
    for d in designs[1:]:
        name, side = d["param"], d["side"]
        factor = 1 - b7.RANGE_PERCENT / 100 if side == "low" else 1 + b7.RANGE_PERCENT / 100
        assert d["overrides"] == {name: pytest.approx(getattr(base, name) * factor)} and side in ("low", "high")
    assert [(d["param"], d["side"]) for d in designs[1:]] == [(n, s) for n in b7.PARAMETERS for s in ("low", "high")]
    assert b7.screen([d["overrides"] for d in designs])["count"] == 0                                                      # and motorgeo accepts all nine


def test_the_swing_of_a_parameter_is_its_torque_at_the_high_end_minus_at_the_low_end():
    from bench import b7_salib as b7
    torque = {"centre": 1.12, "slot_t low": 1.10, "slot_t high": 1.13, "slot_w low": 1.00, "slot_w high": 1.25, "mag_y low": 1.16, "mag_y high": 1.08,
              "v_half_deg low": 1.121, "v_half_deg high": 1.119}
    swings = b7.swings(torque)
    assert swings == {"slot_t": pytest.approx(0.03), "slot_w": pytest.approx(0.25), "mag_y": pytest.approx(-0.08), "v_half_deg": pytest.approx(-0.002)}
    assert b7.swing_ranking(swings) == ["slot_w", "mag_y", "slot_t", "v_half_deg"]                                          # by size, not by sign


def test_the_ranking_is_by_mu_star_with_the_largest_first():
    from bench import b7_salib as b7
    result = {"names": ["slot_t", "slot_w", "mag_y", "v_half_deg"], "mu_star": [0.02, 0.12, 0.04, 0.005], "mu": [0.02, 0.12, -0.04, -0.005],
              "sigma": [0.001, 0.01, 0.002, 0.0004], "mu_star_conf": [0.001, 0.01, 0.003, 0.0005]}
    assert b7.ranking(result) == ["slot_w", "mag_y", "slot_t", "v_half_deg"]
    lines = b7.ranking_lines(result, {"slot_w": 0.12, "mag_y": -0.04, "slot_t": 0.02, "v_half_deg": -0.005})
    assert len(lines) == 4 and lines[0].startswith("1. slot_w") and "mu* 0.1200" in lines[0] and b7.PHYSICS["slot_w"] in lines[0]
    assert lines[1].startswith("2. mag_y") and "mu -0.0400" in lines[1] and "sigma" in lines[1]


# ----------------------------------------------------------------------------------------------- 5. the picture and the checks
def synthetic_analysis():
    return {"names": ["slot_t", "slot_w", "mag_y", "v_half_deg"], "mu_star": [0.0177, 0.125, 0.040, 0.0075], "mu": [0.0177, 0.125, -0.040, -0.0075],
            "sigma": [0.0008, 0.006, 0.002, 0.0006], "mu_star_conf": [0.001, 0.008, 0.003, 0.0007]}


def test_the_figure_draws_the_two_morris_views_the_ruling_asks_for():
    """Left: mu* per parameter as horizontal bars, the largest on top, each with SALib's confidence interval as an error bar. Right: sigma against mu*, each parameter a labelled
    point, with the line sigma = mu* drawn. In N.m, the torque's units."""
    import matplotlib
    matplotlib.use("Agg")
    from bench import b7_salib as b7
    fig = b7.draw_morris(synthetic_analysis(), "a title")
    try:
        left, right = fig.axes
        bars = [p for p in left.patches]
        assert len(bars) == 4
        widths = {round(b.get_width(), 6) for b in bars}
        assert widths == {0.0177, 0.125, 0.04, 0.0075}                                                                     # the values of mu*
        labels_on_y = left.get_yticklabels()
        ys = [t.get_position()[1] for t in labels_on_y]
        top_label = labels_on_y[int(np.argmin(ys) if left.yaxis_inverted() else np.argmax(ys))].get_text()
        bottom_label = labels_on_y[int(np.argmax(ys) if left.yaxis_inverted() else np.argmin(ys))].get_text()
        assert top_label == "slot_w" and bottom_label == "v_half_deg"                                                      # the largest is on top, the smallest at the bottom
        assert any(getattr(c, "errorbar", None) is not None or type(c).__name__ == "ErrorbarContainer" for c in left.containers)      # the confidence interval as an error bar
        bar_errors = [c for c in left.containers if type(c).__name__ == "ErrorbarContainer"]
        assert bar_errors and len(bar_errors[0].lines[2]) >= 1                                                           # (the bars of the interval are drawn)
        assert "N" in left.get_xlabel() and "mu*" in left.get_xlabel()
        labels = {t.get_text() for t in right.texts}
        assert {"slot_t", "slot_w", "mag_y", "v_half_deg"} <= labels                                                       # each point is labelled
        diagonal = [ln for ln in right.lines if len(ln.get_xdata()) == 2 and np.allclose(ln.get_xdata(), ln.get_ydata())]
        assert diagonal, "the line sigma = mu* is not drawn"
        assert "sigma" in right.get_ylabel() and "mu*" in right.get_xlabel() and "N" in right.get_ylabel()
    finally:
        import matplotlib.pyplot as plt
        plt.close(fig)


def test_the_picture_is_saved_tight_with_no_blank_margin_and_a_title_that_fits(tmp_path):
    import matplotlib.image as image
    from bench import b7_salib as b7
    target = tmp_path / "figure.png"
    title = "Morris screening (SALib 1.6.0) of four rotor dimensions " * 4 + "\nsecond line"
    b7.fig_morris(synthetic_analysis(), title, target)
    picture = image.imread(str(target))
    height, width = picture.shape[:2]
    assert width <= 1.1 * 11.5 * 130 and height > 0.3 * width                                                              # the title did not widen it; the plots are most of it
    white = (picture[..., :3] > 0.97).all(axis=2)

    def blank_border(lines) -> int:
        """How many lines (rows or columns) at the start are completely white."""
        count = 0
        for line in lines:
            if not line.all():
                break
            count += 1
        return count
    borders = {"top": blank_border(white), "bottom": blank_border(white[::-1]), "left": blank_border(white.T), "right": blank_border(white.T[::-1])}
    assert all(b <= 0.04 * min(height, width) for b in borders.values()), f"a blank margin (white lines at each edge): {borders}"      # the 0.12 inch of padding is 16 pixels; a blank half would be hundreds


def good_context():
    """What make_checks needs, for a run in which everything is as the brief says."""
    from bench import b7_salib as b7
    X = b7.sample()
    coef = np.array([0.5, 4.0, -1.2, -0.2])
    Y = 1.12 + (X - X.mean(axis=0)) @ coef
    analysis = b7.analyse(X, Y)
    swings = {n: float(c * (hi - lo)) for n, c, (lo, hi) in zip(b7.PARAMETERS, coef, b7.problem()["bounds"])}
    return {"X1": X, "X2": X.copy(), "Y1": Y, "Y2": Y.copy(), "analysis1": analysis, "analysis2": b7.analyse(X, Y.copy()), "rejected_grid": b7.screen([]),
            "grid_points": 256, "rejected_samples": b7.screen([]), "solves": 60, "solves_run2": 60, "oat_solves": 18, "swings": swings, "known": b7.known_answer(X), "figure": None,
            "centre_torque": 1.1249, "b1": {"best_sample_T_Nm": 1.12914, "best_sample_gamma_deg": 23.0, "path": "out/bench/B1/results.json"}, "b1_reason": None, "mesh_nodes": (8500, 9000)}


def test_a_run_in_which_everything_is_as_the_brief_says_passes_every_gating_check(tmp_path):
    from bench import b7_salib as b7
    figure = tmp_path / "figure.png"
    figure.write_bytes(b"x" * 30000)
    context = good_context()
    context["figure"] = figure
    checks = b7.make_checks(context)
    names = [c["name"] for c in checks]
    assert len(names) == len(set(names))
    assert all(c["verdict"] == "PASS" for c in checks if c.get("gating", True)), [c["name"] for c in checks if c["verdict"] != "PASS"]
    for needle in ("rejected", "identical", "ranking", "known answer", "solves"):
        assert any(needle in n for n in names), needle


def test_each_failure_the_brief_names_turns_its_check_red():
    from bench import b7_salib as b7

    def verdicts(**changes):
        context = good_context()
        context["figure"] = None
        context.update(changes)
        return {c["name"]: c["verdict"] for c in b7.make_checks(context) if c.get("gating", True)}
    base = verdicts()
    assert "FAIL" not in {v for k, v in base.items() if "picture" not in k}
    rejected = {"points": 30, "count": 1, "rejected": [({"slot_w": 12.0}, "reaches the bridge")]}
    assert [k for k, v in verdicts(rejected_samples=rejected).items() if v == "FAIL" and "rejected" in k]
    Y2 = good_context()["Y1"].copy()
    Y2[5] = np.nextafter(Y2[5], 2.0)
    assert [k for k, v in verdicts(Y2=Y2).items() if v == "FAIL" and "identical" in k]                                       # one ulp is enough to fail "identical"
    assert [k for k, v in verdicts(solves=61).items() if v == "FAIL" and "solves" in k]
    assert [k for k, v in verdicts(solves_run2=59).items() if v == "FAIL" and "solves" in k]                                 # the second session counts too
    assert [k for k, v in verdicts(oat_solves=16).items() if v == "FAIL" and "solves" in k]                                  # and so do the one-at-a-time designs of the first
    swapped = {"slot_t": 9.0, "slot_w": 0.1, "mag_y": -0.2, "v_half_deg": -0.1}
    assert [k for k, v in verdicts(swings=swapped).items() if v == "FAIL" and "ranking" in k]


# ----------------------------------------------------------------------------------------------- 6. the FEMM side
def test_a_design_that_motorgeo_rejects_never_reaches_femm(tmp_path):
    """The worker builds the design and checks it before it touches FEMM: a rejected sample is an error that names the reason, never a number."""
    from bench import femm_points
    task = {"id": 0, "tag": "bad", "overrides": {"slot_w": 12.0}, "fem_dir": str(tmp_path), "points": [{"theta_e": 0.0, "id": -4.0, "iq": 9.0}]}
    with pytest.raises(Exception, match="reaches the bridge"):
        femm_points.solve_points(task)
    assert list(tmp_path.iterdir()) == []                                                                                  # no model file was written


def test_the_worker_is_a_module_level_function_that_a_pool_can_send_by_name():
    from bench import femm_points
    blob = pickle.dumps(femm_points.solve_points)
    assert pickle.loads(blob) is femm_points.solve_points and femm_points.solve_points.__module__ == "bench.femm_points"
    assert pickle.loads(pickle.dumps(femm_points.run_tasks)) is femm_points.run_tasks


def test_b7_opens_a_femm_pool_in_one_place_with_a_time_limit_and_calls_it_for_each_of_the_two_runs():
    """'Run the whole sample and evaluate step twice, in two femm_pool sessions': one function opens the pool (so that the two runs cannot differ in how they are made), and the
    case calls it twice. Every femm_pool call has its task_timeout_s (tests/test_bench.py scans all of bench/ for that)."""
    from test_bench import queue_calls
    from test_check_tools import scopes_calling
    source = (P / "bench" / "b7_salib.py").read_text(encoding="utf-8")
    found, missing = queue_calls(source)
    assert found == 1 and missing == []
    assert scopes_calling(source, None, "run_session") == ["run", "run"]
    assert scopes_calling(source, "femm_queue", "femm_pool") == ["run_session"]


def test_the_figure_title_says_what_the_ruling_asks_and_every_line_fits_the_figure():
    """The title says what was varied, the output (the mean torque at gamma = 22.9 degrees, 10 A), r x (k + 1) = the number of solves, the levels and the seed; each line is at
    most bench.TITLE_CHARS long, so that wrap_title leaves it as written (no orphan word, no widened figure)."""
    import bench
    from bench import b7_salib as b7
    title = b7.figure_title("1.6.0")
    for needle in ("slot_t, slot_w, mag_y, v_half_deg", "+-1.5 %", "mean torque", "22.9", "10 A", "6 x 5 = 30 samples", "60 solves", "4 levels", "seed 1", "SALib 1.6.0", "coarse mesh"):
        assert needle in title, needle
    assert all(len(line) <= bench.TITLE_CHARS for line in title.split("\n")) and bench.wrap_title(title) == title



# ----------------------------------------------------------------------------------------------- 7. what B7 and B8 share (fix round 1 of Task 13)
def module_level_assignments(source: str) -> set[str]:
    """The names a module assigns at its top level (so that a constant typed in two modules is found)."""
    return {t.id for node in ast.parse(source).body if isinstance(node, ast.Assign) for t in node.targets if isinstance(t, ast.Name)}


def test_the_rotor_positions_and_the_task_limit_of_the_femm_side_are_written_once_in_femm_points():
    """POSITIONS_EL (the pair of rotor positions 30 electrical degrees apart) and TASK_TIMEOUT_S were typed in B7 and again in B8 (review M6). They belong to the FEMM side that
    both share: femm_points holds them, and B7 and B8 import them."""
    from bench import b7_salib as b7
    from bench import b8_matlab as b8
    from bench import femm_points
    assert femm_points.POSITIONS_EL == (0.0, 30.0) and femm_points.TASK_TIMEOUT_S == 120
    for module in (b7, b8):
        assert module.POSITIONS_EL is femm_points.POSITIONS_EL and module.TASK_TIMEOUT_S == femm_points.TASK_TIMEOUT_S
        source = (P / "bench" / f"{module.__name__.split('.')[-1]}.py").read_text(encoding="utf-8")
        assert not {"POSITIONS_EL", "TASK_TIMEOUT_S"} & module_level_assignments(source), module.__name__


def test_the_pool_record_says_the_task_limit_it_is_given_and_b7_and_b8_do_not_patch_it_afterwards():
    """pool_record is B1's and said B1's limit (300 s); B7 and B8 overwrote the field after the call (review M6). The limit is a parameter now: B1's own by default."""
    from bench import b1_femm
    from bench import b7_salib, b8_matlab                                                                                   # noqa: F401  (the modules whose source is read)

    class Pool:
        workers = 8
    results = [{"worker": {"pid": 11, "class": 0x4000, "femm_pid": 12, "femm_class": 0x4000}}]
    assert b1_femm.pool_record("x", Pool(), results, 1.0)["task_timeout_s"] == b1_femm.TASK_TIMEOUT_S == 300                 # B1's own, as before
    assert b1_femm.pool_record("x", Pool(), results, 1.0, task_timeout_s=120)["task_timeout_s"] == 120
    for name in ("b7_salib.py", "b8_matlab.py"):
        source = (P / "bench" / name).read_text(encoding="utf-8")
        assert 'record["task_timeout_s"]' not in source, name
        calls = [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "pool_record"]
        assert calls and all("task_timeout_s" in {k.arg for k in c.keywords} for c in calls), name


B1_RESULTS = {"summary": {"none": {"psi_pm_mWb": 16.629, "mtpa": {"bench_method": {"gamma_deg": 22.936, "best_sample_gamma_deg": 23.0, "best_sample_T_Nm": 1.12914},
                                                                  "curve": {"0.0": {"T_age_Nm": 0.998}, "22.5": {"T_age_Nm": 1.129}, "90.0": {"T_age_Nm": 0.0}}}}}}


def put_b1(root, data, text=None):
    """Write a B1 results.json under `root` (the file the case reads); `text` is written as it is when given."""
    folder = root / "B1"
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "results.json").write_text(text if text is not None else json.dumps(data), encoding="utf-8")
    return folder / "results.json"


def test_b1s_numbers_are_read_by_one_function_from_the_active_root_and_the_other_roots_are_not_this_runs(monkeypatch, tmp_path):
    """B7 and B8 each wrote b1_reference() (review M5), with different fields and a path typed as out\\bench. One function in femm_points reads what both need, from the file in
    the ACTIVE root (a run on another root has its own B1 or none): (the numbers, None) or (None, why not)."""
    import bench
    from bench import femm_points
    out = use_config(monkeypatch, tmp_path)
    put_b1(out / "bench", B1_RESULTS)                                                                                         # the default root's B1
    with bench.use_root(out / "bench" / "B9" / "run1") as root:
        data, why = femm_points.b1_reference()
        assert data is None and why.startswith("no B1 results in this run's root") and str(root / "B1" / "results.json") in why        # the default root's file is not this run's
        path = put_b1(root, B1_RESULTS)
        data, why = femm_points.b1_reference()
    assert why is None and data["path"] == str(path)
    assert (data["psi_pm_mWb"], data["gamma_deg"], data["best_sample_T_Nm"], data["best_sample_gamma_deg"]) == (16.629, 22.936, 1.12914, 23.0)
    assert data["curve"] == {0.0: 0.998, 22.5: 1.129, 90.0: 0.0} and data["peak_Nm"] == 1.129
    for name in ("b7_salib.py", "b8_matlab.py"):
        source = (P / "bench" / name).read_text(encoding="utf-8")
        assert "def b1_reference" not in source and "femm_points.b1_reference()" in source and "is not in this run's root" not in source, name


@pytest.mark.parametrize("dotted", ["summary", "summary.none", "summary.none.psi_pm_mWb", "summary.none.mtpa", "summary.none.mtpa.curve", "summary.none.mtpa.bench_method",
                                    "summary.none.mtpa.bench_method.gamma_deg", "summary.none.mtpa.bench_method.best_sample_T_Nm", "summary.none.mtpa.bench_method.best_sample_gamma_deg"])
def test_b1_results_that_lack_a_number_say_which_one_and_are_not_taken_for_missing_results(monkeypatch, tmp_path, dotted):
    """'No B1 results in this root' and 'B1's results lack key K' are different things to tell a reader: the first says to run B1, the second that its file is not what was expected."""
    import bench
    from bench import femm_points
    out = use_config(monkeypatch, tmp_path)
    broken = copy.deepcopy(B1_RESULTS)
    *parents, last = dotted.split(".")
    node = broken
    for key in parents:
        node = node[key]
    del node[last]
    with bench.use_root(out / "bench" / "mine") as root:
        put_b1(root, broken)
        data, why = femm_points.b1_reference()
    assert data is None and "lack" in why and dotted in why and "no B1 results" not in why, why


def test_b1_results_that_are_not_json_or_hold_a_value_that_is_not_a_number_say_so(monkeypatch, tmp_path):
    import bench
    from bench import femm_points
    out = use_config(monkeypatch, tmp_path)
    with bench.use_root(out / "bench" / "mine") as root:
        put_b1(root, None, text="{not json")
        data, why = femm_points.b1_reference()
        assert data is None and "not JSON" in why and "no B1 results" not in why
        words = copy.deepcopy(B1_RESULTS)
        words["summary"]["none"]["psi_pm_mWb"] = "sixteen"
        put_b1(root, words)
        data, why = femm_points.b1_reference()
        assert data is None and "not a number" in why and "no B1 results" not in why


def test_salibs_defaults_are_cited_with_the_lines_of_its_installed_source_and_each_line_says_it():
    """B7 cites where its choices come from (review M10): SALib 1.6.0's own source in the environment, the file and the line, with the text on that line. Read back here, so that
    a SALib that moved a line, or changed a default, shows. The docstring of B7 names every file:line."""
    from pathlib import Path

    import SALib

    from bench import b7_salib as b7
    folder = Path(SALib.__file__).parent
    assert b7.SALIB_SOURCES
    wanted = ("4 levels", "100 resamples", "95 % confidence", "the step 2/3", "mu* is the mean of", "sigma")
    assert all(any(what.startswith(w) for what, *_ in b7.SALIB_SOURCES) for w in wanted)
    for what, file, line, text in b7.SALIB_SOURCES:
        lines = (folder / file).read_text(encoding="utf-8").splitlines()
        assert text in lines[line - 1], (what, file, line, lines[line - 1])
        assert f"{file}:{line}" in b7.__doc__, (what, file, line)


def test_each_index_is_given_with_its_spread_as_a_share_of_mu_star_and_the_docstring_states_the_measured_range():
    """The docstring said 'every sigma is a few per cent of its mu*' (review C3); the measured shares (sigma_over_mu_star of the real run's results.json, unrounded) are 1.9 % (slot_w),
    2.8 % (mag_y), 10 % (slot_t) and 26.5 % (v_half_deg), and the docstring states them."""
    from bench import b7_salib as b7
    analysis = {"names": ["slot_t", "slot_w", "mag_y", "v_half_deg"], "mu_star": [0.0186, 0.1261, 0.0414, 0.0074], "sigma": [0.0019, 0.0024, 0.0011, 0.0020]}
    ratios = b7.sigma_over_mu_star(analysis)
    assert list(ratios) == ["slot_t", "slot_w", "mag_y", "v_half_deg"]
    assert ratios == {"slot_t": pytest.approx(0.1022, abs=1e-4), "slot_w": pytest.approx(0.0190, abs=1e-4), "mag_y": pytest.approx(0.0266, abs=1e-4), "v_half_deg": pytest.approx(0.2703, abs=1e-4)}
    doc = b7.__doc__
    assert "a few per cent" not in doc and all(share in doc for share in ("1.9 %", "2.8 %", "10 %", "26.5 %"))


def test_the_solves_check_states_both_counts_the_morris_solves_of_a_session_and_the_one_at_a_time_ones():
    """60 Morris solves a session, two sessions, plus 18 one-at-a-time solves in the first: 138 in the case (a design decision: the 18 stay, as the physics cross-check)."""
    from bench import b7_salib as b7
    solves = next(c for c in b7.make_checks(good_context()) if "solves" in c["name"])
    text = solves["expected"] + " | " + solves["measured"]
    assert solves["verdict"] == "PASS" and all(needle in text for needle in ("60", "18", "138")), text
    assert b7.oat_solves() == 18 and b7.solves_in_the_case() == 138


def test_without_b1s_results_b7_skips_the_comparison_and_says_why_they_are_missing():
    from bench import b7_salib as b7
    context = good_context()
    context["b1"], context["b1_reason"] = None, "B1's results (D:\\x\\B1\\results.json) lack summary.none.mtpa.bench_method.best_sample_T_Nm"
    skipped = [c for c in b7.make_checks(context) if c["verdict"] == "SKIPPED"]
    assert len(skipped) == 1 and "lack summary.none.mtpa.bench_method.best_sample_T_Nm" in skipped[0]["measured"] and skipped[0]["gating"] is False
