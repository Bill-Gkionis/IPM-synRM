"""Tests for B1 (bench/b1_femm.py): Python -> FEMM 5.0. Everything that needs no solver.

Run from P:   uv run pytest tests/test_bench_b1.py -q      (a few seconds; no FEMM is started)

In groups:
  1. The arithmetic: the positions and angles of the sweeps (the same ones as barrier_study.py), the estimators (the fundamental of a
     flux-linkage waveform, the peak of a torque-vs-angle parabola), the drawn rotor position, and the geometry check that barrier="none"
     and barrier="steel" are one machine (the bench's geometry and the README row's).
  2. The analysis on synthetic rows (formulas with known answers): the numbers, the verdicts, a miss, a cogging wave that does not repeat, the
     three checks of the study's own numbers that are skipped (not dropped) where its file is not there.
  3. A stage on the pool (a stub pool): a chunk that did not come back is an error that names it, whatever the first id of the stage.
  4. The limits the bench chose, which say so next to their number, and the brief's, which do not.
"""
import json
import math

import numpy as np
import pytest

import motorgeo as mg
from bench_test_tools import use_config


# ----------------------------------------------------------------------------------------------- 1. the arithmetic
def test_b1_sweeps_use_the_positions_and_angles_of_barrier_study():
    """The README's band-free row comes from barrier_study.py, so B1 has to sample the same positions or the numbers are not comparable."""
    from bench import b1_femm
    import barrier_study
    assert b1_femm.noload_positions() == [float(t) for t in range(0, 360, 3)]
    assert b1_femm.coarse_gammas() == [float(g) for g in barrier_study.GAMMAS] and len(b1_femm.coarse_gammas()) == 13
    assert b1_femm.mean_positions() == [float(t) for t in barrier_study.MEAN_POS] == [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]
    assert b1_femm.I_RATED == barrier_study.I_RATED == 10
    assert list(b1_femm.GAP_DEG) == list(barrier_study.GAP_DEG) and len(b1_femm.GAP_DEG) == 360
    assert b1_femm.VARIANTS["steel"] == barrier_study.VARIANTS["no band"] == "steel"          # the row's band-free variant
    assert b1_femm.VARIANTS["none"] == "none" and b1_femm.VARIANTS["air"] == barrier_study.VARIANTS["with band"]


def test_the_cogging_period_is_5_mechanical_degrees_because_lcm_of_36_and_8_is_72():
    from bench import b1_femm
    p = mg.Params()
    assert math.lcm(p.slots, 2 * p.pole_pairs) == 72
    assert b1_femm.cogging_period_mech(p) == pytest.approx(5.0)
    assert b1_femm.cogging_period_el(p) == pytest.approx(20.0)           # 5 mechanical degrees x 4 pole pairs
    positions = b1_femm.cogging_positions()
    assert positions[0] == 0.0 and positions[-1] == 40.0 and len(positions) == 81      # two periods, so that "it repeats" is tested
    assert len([t for t in positions if t <= 20.0]) >= 11                 # the brief: at least 11 positions in the first period
    assert all(b - a == pytest.approx(0.5) for a, b in zip(positions, positions[1:]))
    assert 20.0 in positions                                              # the end of the first period is a solved point


def test_current_angle_is_measured_from_the_q_axis_towards_minus_d():
    from bench import b1_femm
    i_d, i_q = b1_femm.currents(0.0, 10.0)
    assert (abs(i_d), i_q) == (0.0, 10.0)
    i_d, i_q = b1_femm.currents(90.0, 10.0)
    assert (i_d, abs(i_q)) == (pytest.approx(-10.0), pytest.approx(0.0, abs=1e-9))
    i_d, i_q = b1_femm.currents(30.0, 10.0)
    assert i_d == pytest.approx(-5.0) and i_q == pytest.approx(10 * math.cos(math.radians(30)))
    assert math.hypot(*b1_femm.currents(22.8, 10.0)) == pytest.approx(10.0)


def test_the_fundamental_of_a_sampled_period_is_its_cosine_amplitude_whatever_the_offset_and_the_harmonics():
    from bench import b1_femm
    n = 120
    k = np.arange(n)
    wave = 3.0 + 2.5 * np.cos(2 * np.pi * k / n + 0.7) + 0.4 * np.cos(3 * 2 * np.pi * k / n) + 0.1 * np.cos(5 * 2 * np.pi * k / n + 1.0)
    assert b1_femm.fundamental(wave) == pytest.approx(2.5, abs=1e-12)
    assert b1_femm.fundamental(np.full(n, 4.0)) == pytest.approx(0.0, abs=1e-12)
    spectrum = b1_femm.harmonic_amplitudes(wave)
    assert spectrum[0] == pytest.approx(3.0) and spectrum[1] == pytest.approx(2.5) and spectrum[3] == pytest.approx(0.4)
    assert spectrum[5] == pytest.approx(0.1) and len(spectrum) == n // 2 + 1


def test_the_peak_of_a_torque_curve_comes_from_a_parabola_through_the_best_sample_and_its_neighbours():
    from bench import b1_femm
    gammas = np.arange(0, 90.1, 7.5)
    exact = 1.1 - 0.0004 * (gammas - 22.3) ** 2                            # a parabola: the estimate must be exact
    g, t = b1_femm.parabola_peak(gammas, exact)
    assert g == pytest.approx(22.3, abs=1e-9) and t == pytest.approx(1.1, abs=1e-12)
    sine = np.sin(np.radians(2 * gammas + 20))                              # not a parabola: close, and never above the true peak by much
    g, t = b1_femm.parabola_peak(gammas, sine)
    assert abs(g - 35.0) < 1.0 and abs(t - 1.0) < 0.01
    edge = 1.0 - 0.001 * (gammas - 1.0) ** 2                                # a peak at the first sample: the neighbours are taken from inside
    g, t = b1_femm.parabola_peak(gammas, edge)
    assert g == pytest.approx(1.0, abs=1e-9)
    with pytest.raises(ValueError):
        b1_femm.parabola_peak([0.0, 1.0], [0.0, 1.0])                       # three points at least


def test_the_one_degree_steps_surround_the_coarse_peak_and_stay_inside_0_to_90():
    from bench import b1_femm
    assert b1_femm.refine_gammas(22.4) == [18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0]
    assert b1_femm.refine_gammas(22.6)[4] == 23.0
    assert b1_femm.refine_gammas(1.0) == [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
    assert b1_femm.refine_gammas(89.0) == [85.0, 86.0, 87.0, 88.0, 89.0, 90.0]
    assert all(b - a == 1.0 for a, b in zip(b1_femm.refine_gammas(40.0), b1_femm.refine_gammas(40.0)[1:]))


def test_the_drawn_position_is_the_one_where_the_rotor_angle_of_the_model_is_zero():
    """No-load FEMM and Elmer are compared at the position the geometry is drawn in: the AGE inner angle is then 0 and the stator
    frame and the rotor frame are the same frame (B3 depends on it)."""
    from bench import b1_femm
    p = mg.Params(barrier="none")
    axes, _ = mg.phase_axes(p)
    theta_e = b1_femm.drawn_theta_e(p, axes["A"])
    theta_d0 = 90 / p.pole_pairs
    assert theta_e / p.pole_pairs + axes["A"] - theta_d0 == pytest.approx(0.0, abs=1e-9)       # FemmModel.solve_abc's rotor_deg
    assert -360.0 < theta_e < 360.0                                       # (it is negative here: the drawn rotor sits "before" electrical 0)


def test_every_job_is_a_plain_dict_with_its_own_name_so_that_no_two_solves_share_a_file_and_the_table_can_be_rebuilt():
    from bench import b1_femm
    jobs = b1_femm.stage_one_jobs()
    assert {j["variant"] for j in jobs} == {"none", "steel", "air"}
    kinds = {(j["variant"], j["study"]) for j in jobs}
    assert ("none", "noload") in kinds and ("none", "gamma") in kinds and ("none", "cogging") in kinds and ("none", "drawn") in kinds
    assert ("air", "cogging") in kinds and ("air", "noload") not in kinds and ("air", "gamma") not in kinds    # the README's 25.7 only
    per = {v: sum(1 for j in jobs if j["variant"] == v) for v in ("none", "steel", "air")}
    assert per["none"] == per["steel"] == 120 + 13 * 6 + 81 + 1 and per["air"] == 81
    json.dumps(jobs)                                                       # picklable and JSON-able: they cross a pipe to a worker
    assert len({(j["variant"], j["study"], j["theta_e"], j.get("gamma"), j["id"], j["iq"]) for j in jobs}) == len(jobs)
    drawn = [j for j in jobs if j["study"] == "drawn"]
    assert len(drawn) == 2 and all(j["id"] == 0.0 and j["iq"] == 0.0 and len(j["gap_deg"]) == 360 for j in drawn)


def test_chunks_keep_one_variant_and_cover_every_job_exactly_once():
    from bench import b1_femm
    jobs = b1_femm.stage_one_jobs()
    chunks = b1_femm.make_chunks(jobs, size=14, mesh_scale=1.0, fem_dir="x")
    assert sum(len(c["jobs"]) for c in chunks) == len(jobs)
    assert all(len({j["variant"] for j in c["jobs"]}) == 1 and c["variant"] == c["jobs"][0]["variant"] for c in chunks)
    assert all(len(c["jobs"]) <= 14 for c in chunks)
    assert len({c["id"] for c in chunks}) == len(chunks)
    assert all(c["mesh_scale"] == 1.0 for c in chunks)


def test_the_bench_geometry_and_the_readme_rows_geometry_are_the_same_machine_with_a_different_mesh():
    """The README's band-free row is barrier="steel" (the band's outline stays in the drawing, filled with steel), not barrier="none"
    (R11: the band is not drawn at all). Same steel and magnets, so the physics is the same; the mesh is not."""
    from bench import b1_femm
    facts = b1_femm.same_machine()
    assert facts["only_in_steel"] and all(name.startswith("barrier_") for name in facts["only_in_steel"])
    assert facts["only_in_none"] == []
    assert facts["n_band_regions"] == 4                                  # two poles x two magnets in the 90-degree sector
    assert facts["steel_area_none_mm2"] == pytest.approx(facts["steel_area_steel_mm2"] + facts["band_area_mm2"], abs=1e-6)
    assert facts["other_regions_identical"] is True                       # every other region has the same outline and label
    assert facts["magnet_area_none_mm2"] == pytest.approx(facts["magnet_area_steel_mm2"], abs=1e-9)
    assert facts["mesh_size_band_mm"] == 0.5 and facts["mesh_size_rotor_steel_mm"] == 1.5      # why the meshes differ: 0.5 mm blocks in the steel


# ----------------------------------------------------------------------------------------------- 2. the analysis, on synthetic rows
def synthetic_rows(psi_mwb=16.61, t_peak=1.127, g_peak=22.8, cog_pp_mnm=52.0, air_pp_mnm=25.7, scale=1.0):
    """Rows as solve_chunk makes them, from formulas whose answers are known: psi_a = psi cos(theta_e), a parabola of torque over gamma that
    peaks at (g_peak, t_peak), a sine of cogging with period 20 electrical degrees. No FEMM: this is how the analysis is tested."""
    from bench import b1_femm
    rows = []

    def row(variant, study, theta_e, gamma=None, torque=0.0, psi=0.0, i_rated=0.0, **extra):
        r = {"study": study, "variant": variant, "mesh_scale": scale, "theta_e": theta_e, "I": i_rated, "gamma": gamma, "id": 0.0, "iq": 0.0,
             "psi_a": psi * math.cos(math.radians(theta_e)), "psi_b": 0.0, "psi_c": 0.0, "psi_d": psi, "psi_q": 0.0, "T_age": torque, "T_wst": torque,
             "T_dq": torque, "rotor_deg": 0.0, "psi_retries": 0, "seconds": 0.5}
        r.update(extra)
        return r

    coarse = b1_femm.coarse_gammas()
    fine = [float(g) for g in range(19, 28)]
    for variant in ("none", "steel"):
        for th in b1_femm.noload_positions():
            rows.append(row(variant, "noload", th, psi=psi_mwb * 1e-3))
        for g in sorted(set(coarse) | set(fine)):
            for th in b1_femm.mean_positions():
                rows.append(row(variant, "gamma", th, gamma=g, i_rated=10.0, torque=t_peak - 0.0004 * (g - g_peak) ** 2))
        for th in b1_femm.cogging_positions():
            rows.append(row(variant, "cogging", th, torque=0.5 * cog_pp_mnm * 1e-3 * math.sin(2 * math.pi * th / 20.0)))
        rows.append(row(variant, "drawn", -170.1, psi=psi_mwb * 1e-3, nodes=17934, elements=35000, ans="x.ans", gap_deg=[0.0, 0.25],
                        gap_br=[0.1, 0.2], gap_bt=[0.0, 0.0]))
    for th in b1_femm.cogging_positions():
        rows.append(row("air", "cogging", th, torque=0.5 * air_pp_mnm * 1e-3 * math.sin(2 * math.pi * th / 20.0)))
    return rows


def test_the_analysis_turns_rows_into_the_numbers_and_the_verdicts_of_the_brief(monkeypatch):
    from bench import b1_femm
    monkeypatch.setattr(b1_femm, "stored_row", lambda: None)
    summary, checks = b1_femm.analyse(synthetic_rows(), mg.Params(), b1_femm.same_machine())
    for v in ("none", "steel"):
        assert summary[v]["psi_pm_mWb"] == pytest.approx(16.61, rel=1e-9)
        assert summary[v]["mtpa"]["bench_method"]["gamma_deg"] == pytest.approx(22.8, abs=1e-6)
        assert summary[v]["mtpa"]["bench_method"]["T_Nm"] == pytest.approx(1.127, abs=1e-9)
        assert summary[v]["mtpa"]["row_method"]["T_Nm"] == pytest.approx(1.127, abs=1e-9)       # the row's estimator, on the coarse samples alone
        assert summary[v]["mtpa"]["best_coarse_gamma_deg"] == 22.5
        assert summary[v]["cogging"]["pp_mNm"] == pytest.approx(52.0, rel=1e-3) and summary[v]["cogging"]["repeat_percent_of_pp"] < 1e-6
        assert summary[v]["cogging"]["period_mech_deg"] == 5.0 and summary[v]["cogging"]["positions_first_period"] == 41
        assert summary[v]["nodes"] == 17934
    assert summary["air"]["cogging"]["pp_mNm"] == pytest.approx(25.7, rel=1e-3)
    assert checks and all(c["verdict"] == "PASS" for c in checks if c["gating"]), [c for c in checks if c["gating"] and c["verdict"] != "PASS"]
    assert all(c["verdict"] == "SKIPPED" for c in checks if not c["gating"])                    # without the study's file the three checks of it are skipped
    names = [c["name"] for c in checks]
    assert any(n.startswith("none: psi_PM") for n in names) and any(n.startswith("steel: MTPA") for n in names) and any(n.startswith("air: cogging") for n in names)
    assert all(c["source"] for c in checks)                      # every number says where its expectation comes from


def test_a_number_outside_the_brief_s_tolerance_is_a_failed_check_named_for_its_variant_and_the_tolerance_is_never_moved(monkeypatch):
    from bench import b1_femm
    monkeypatch.setattr(b1_femm, "stored_row", lambda: None)
    assert b1_femm.TOLERANCE_PERCENT == 0.5 and b1_femm.ROW == {"psi_pm_mWb": 16.61, "T_mtpa_Nm": 1.127}      # the brief's numbers, as written
    _, checks = b1_femm.analyse(synthetic_rows(psi_mwb=16.61 * 1.01, t_peak=1.127 * 0.99), mg.Params(), b1_femm.same_machine())
    failed = [c for c in checks if c["verdict"] == "FAIL"]
    assert {c["name"].split(":")[0] for c in failed} == {"none", "steel"}
    psi = next(c for c in failed if c["name"].startswith("none: psi_PM"))
    assert psi["difference_percent"] == pytest.approx(1.0, abs=1e-6) and psi["tolerance_percent"] == 0.5 and psi["expected"] == 16.61


def test_a_cogging_wave_that_does_not_repeat_is_a_failed_check(monkeypatch):
    from bench import b1_femm
    monkeypatch.setattr(b1_femm, "stored_row", lambda: None)
    rows = synthetic_rows()
    for r in rows:
        if r["variant"] == "none" and r["study"] == "cogging" and r["theta_e"] > 20.0:
            r["T_age"] *= 0.5                                      # the second period is half of the first
    _, checks = b1_femm.analyse(rows, mg.Params(), b1_femm.same_machine())
    bad = next(c for c in checks if c["name"].startswith("none: cogging torque repeats"))
    assert bad["verdict"] == "FAIL" and "mN.m" in bad["measured"]


def test_the_stored_row_of_the_study_adds_the_same_model_checks_for_steel_only(monkeypatch):
    from bench import b1_femm
    stored = {"no band": {"psi_pm_mWb": 16.607303314702392, "T_mtpa_10A": 1.1274468002754214, "cogging_pp_mNm": 52.428870872152295},
              "with band": {"psi_pm_mWb": 12.126829312310509}}
    monkeypatch.setattr(b1_femm, "stored_row", lambda: stored)
    _, checks = b1_femm.analyse(synthetic_rows(psi_mwb=16.607303314702392, t_peak=1.1274468002754214, cog_pp_mnm=52.428870872152295), mg.Params(),
                                b1_femm.same_machine())
    same = [c for c in checks if "the study's own" in c["name"]]
    assert len(same) == 3 and all(c["name"].startswith("steel:") for c in same) and all(c["tolerance_percent"] == 0.01 for c in same)
    assert all(c["verdict"] == "PASS" for c in same)


def test_without_the_study_s_summary_the_three_same_model_checks_are_still_in_the_list_as_skipped_with_the_reason(monkeypatch):
    """M9 of the review of Task 11: out\\barrier\\summary.csv is git-ignored, so on a PC that never ran barrier_study.py the three checks vanished without a
    word (48 checks here, 45 elsewhere). They are now always in the list: made when the file is there, SKIPPED (never gating) with the reason when not."""
    import bench
    from bench import b1_femm
    stored = {"no band": {"psi_pm_mWb": 16.607303314702392, "T_mtpa_10A": 1.1274468002754214, "cogging_pp_mNm": 52.428870872152295}, "with band": {}}
    rows = synthetic_rows(psi_mwb=16.607303314702392, t_peak=1.1274468002754214, cog_pp_mnm=52.428870872152295)
    monkeypatch.setattr(b1_femm, "stored_row", lambda: stored)
    _, made = b1_femm.analyse(rows, mg.Params(), b1_femm.same_machine())
    monkeypatch.setattr(b1_femm, "stored_row", lambda: None)
    _, not_made = b1_femm.analyse(rows, mg.Params(), b1_femm.same_machine())
    assert len(made) == len(not_made) and [c["name"] for c in made] == [c["name"] for c in not_made]            # the same list on every PC
    skipped = [c for c in not_made if "the study's own" in c["name"]]
    assert len(skipped) == 3 and all(c["verdict"] == "SKIPPED" and c["gating"] is False for c in skipped)
    assert all(c["measured"] == "out\\barrier\\summary.csv is not there: barrier_study.py makes it" for c in skipped)
    assert all(bench.is_gating(c) and c["verdict"] == "PASS" for c in made if "the study's own" in c["name"])        # and when it is there they gate, as before
    assert all("0.01 %" in c["expected"] for c in skipped)                                                    # a skipped check still says what it would have asked


def test_the_stored_row_is_read_from_the_study_s_summary_when_it_is_there(monkeypatch, tmp_path):
    from bench import b1_femm
    out = use_config(monkeypatch, tmp_path)
    assert b1_femm.stored_row() is None
    (out / "barrier").mkdir(parents=True)
    (out / "barrier" / "summary.csv").write_text(",with band,no band,change_pct\npsi_pm_mWb,12.1,16.6,36.9\nT_mtpa_10A,0.859,1.127,31.2\n", encoding="utf-8")
    row = b1_femm.stored_row()
    assert row["no band"]["psi_pm_mWb"] == 16.6 and row["with band"]["T_mtpa_10A"] == 0.859


# ----------------------------------------------------------------------------------------------- 3. a stage on the pool: a chunk that did not come back
class StubPool:
    """A pool that returns what it was given for each chunk, in the order of the chunks (femm_pool.map's contract), and records what it was asked."""

    def __init__(self, answers):
        self.answers, self.asked = answers, []

    def map(self, fn, items):
        items = list(items)
        self.asked.append([c["id"] for c in items])
        return [self.answers(c) for c in items]


def one_job_each(n):
    """n x 14 + 2 jobs of one variant: n + 1 chunks of the pool (CHUNK = 14)."""
    return [{"variant": "none", "study": "noload", "theta_e": float(k), "id": 0.0, "iq": 0.0} for k in range(14 * n + 2)]


def test_a_chunk_that_did_not_come_back_is_an_error_that_names_it_and_what_the_pool_said_whatever_the_first_id_of_the_stage(tmp_path):
    """M3 of the review of Task 11: the message was built with the chunk's ID as an index into the results, which are indexed by position. In the second stage
    of a run (first_id > 0) that raised an IndexError and hid the real news."""
    from bench import Log, b1_femm
    log = Log(tmp_path)
    good = lambda c: {"chunk": c["id"], "rows": [], "meshes": [], "worker": {}}
    for first_id in (0, 7):
        lost_id = first_id + 1
        pool = StubPool(lambda c: {"status": "timeout"} if c["id"] == lost_id else good(c))
        with pytest.raises(RuntimeError) as error:
            b1_femm.run_stage(pool, one_job_each(2), tmp_path, "stage 2", log, first_id=first_id)
        assert pool.asked == [[first_id, first_id + 1, first_id + 2]]
        message = str(error.value)
        assert f"chunk(s) [{lost_id}]" in message and "'status': 'timeout'" in message and "stage 2" in message, message
    pool = StubPool(lambda c: {"status": "died"} if c["id"] in (8, 9) else good(c))                 # two lost: both are named, each with its own news
    with pytest.raises(RuntimeError) as error:
        b1_femm.run_stage(pool, one_job_each(2), tmp_path, "stage 2", log, first_id=7)
    assert "chunk(s) [8, 9]" in str(error.value) and str(error.value).count("'status': 'died'") == 2


def test_a_stage_whose_chunks_all_came_back_returns_them_in_order(tmp_path):
    from bench import Log, b1_femm
    pool = StubPool(lambda c: {"chunk": c["id"], "rows": [{"theta_e": c["id"]}], "meshes": [], "worker": {}})
    results = b1_femm.run_stage(pool, one_job_each(2), tmp_path, "stage 1", Log(tmp_path), first_id=3)
    assert [r["chunk"] for r in results] == [3, 4, 5] and [row["theta_e"] for row in b1_femm.table(results)] == [3, 4, 5]


# ----------------------------------------------------------------------------------------------- 4. the limits the bench chose are labelled as its own
def test_every_limit_that_the_bench_chose_says_so_next_to_its_number(monkeypatch):
    """M1 of the review of Task 11: only the brief's numbers (0.5 % for the row, the 11 positions and the 5 degrees) are the brief's. The cogging repeat of 5 %, the
    mesh-study criterion, the 0.5 % that the README's cogging number is held to (the brief gives that number none), the 0.01 % of the same model against the study's
    own numbers and the 0.1 % between the two torque methods are the bench's own, and the check that uses each one says so in the words of bench.OWN_LIMIT."""
    import bench
    from bench import b1_femm
    stored = {"no band": {"psi_pm_mWb": 16.61, "T_mtpa_10A": 1.127, "cogging_pp_mNm": 52.0}, "with band": {}}
    monkeypatch.setattr(b1_femm, "stored_row", lambda: stored)
    rows = synthetic_rows()
    summary, checks = b1_femm.analyse(rows, mg.Params(), b1_femm.same_machine())
    text = lambda c: f"{c['expected']} {c['source']}"                               # a number's limit is in its source, a check's in its expected text
    own = lambda name: next(c for c in checks if c["name"].startswith(name))
    assert "(" + bench.OWN_LIMIT + ")" in own("geometry:")["expected"]                      # P4 of the re-review: the 1e-6 mm2 the band's area must add up to
    assert "(" + bench.OWN_LIMIT + ")" in own("none: cogging torque repeats")["expected"]
    assert "(" + bench.OWN_LIMIT + ")" in own("steel: cogging torque repeats")["expected"]
    assert "(" + bench.OWN_LIMIT + ")" in own("air: cogging")["source"] and own("air: cogging")["tolerance_percent"] == 0.5          # the number is not moved, only said
    assert "(" + bench.OWN_LIMIT + ")" in own("none: torque from the air-gap element")["source"]
    for c in checks:
        if "the study's own" in c["name"]:
            assert "(" + bench.OWN_LIMIT + ")" in text(c), c["name"]
    for name in ("none: psi_PM", "none: MTPA", "steel: psi_PM (fundamental", "steel: MTPA torque at 10 A"):                     # the brief's own 0.5 % is not labelled as the bench's
        assert bench.OWN_LIMIT not in text(own(name)), name
    extra = synthetic_rows(scale=0.5)
    study = b1_femm.mesh_study_checks(summary, rows + extra)
    assert len(study) == 3 and all("(" + bench.OWN_LIMIT + ")" in c["expected"] for c in study)
