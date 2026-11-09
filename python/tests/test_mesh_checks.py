"""Tests for validation/mesh_checks.py: the model it solves is the model's own, only denser, and its measures are the study's.

The first tests need no FEMM. The tests marked `femm` start real FEMMs (a few seconds each): they pin down what the density
factor does to the node count, and that FEMM 5.0 reads a FEMM 4.2 solution the way FEMM 4.2 did (that one needs FEMM 4.2's old
out folder and is skipped where it is not). Run from P:   uv run pytest tests/test_mesh_checks.py -q
"""
import dataclasses
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import mesh_checks as mc
import motorgeo as mg
import paths

REAL_OLD = Path(os.environ.get("FEMM42_OUT") or paths.path("femm42_out"))      # FEMM 4.2's run: config.toml, unless the variable names another
OUT = Path(__file__).resolve().parent.parent / "out" / "tests" / "mesh_checks"


# ----------------------------------------------------------------------------------------------- the model at density f
def test_density_multiplies_every_mesh_size_and_changes_nothing_else():
    p = mg.Params()
    base, scaled = mg.build(p), mc.scaled_regions(p, 2.0)
    assert len(base) == len(scaled) > 30
    for a, b in zip(base, scaled):
        assert b.mesh == pytest.approx(2.0 * a.mesh)
        assert dataclasses.replace(b, mesh=a.mesh) == a                  # name, material, edges, label, circuit, turns: all the same
    assert [r.mesh for r in mc.scaled_regions(p, 1.0)] == [r.mesh for r in base]


def test_the_density_does_not_touch_the_model_object_it_started_from():
    p = mg.Params()
    before = [r.mesh for r in mg.build(p)]
    mc.scaled_regions(p, 0.5)
    assert [r.mesh for r in mg.build(p)] == before


# ----------------------------------------------------------------------------------------------- the points that are solved
def test_the_job_list_is_the_barrier_studys_points():
    jobs = mc.slice_jobs()
    sets = pd.Series([j["set"] for j in jobs]).value_counts().to_dict()
    assert sets == {"noload": 24, "cogging": 41, "gamma0": 6, "gamma30": 6, "gamma90": 6, "ripple": 61, "check0": 1, "check1": 1,
                    "check2": 1, "check3": 1}
    assert len({(j["set"], j["theta"], j["id"], j["iq"]) for j in jobs}) == len(jobs)          # nothing is solved twice
    assert sum(j["gap"] for j in jobs) == 1                                                    # the air-gap B is sampled once
    ripple = [j for j in jobs if j["set"] == "ripple"]
    assert {(j["id"], j["iq"]) for j in ripple} == {mc.IQ_ID_42}                               # at FEMM 4.2's own best angle
    assert math.degrees(math.atan2(-mc.IQ_ID_42[0], mc.IQ_ID_42[1])) == pytest.approx(mc.GAMMA_MTPA_42, abs=1e-9)


def frame(sets):
    rows = []
    for name, thetas, fn in sets:
        rows += [{"set": name, "theta": t, **fn(t)} for t in thetas]
    return pd.DataFrame(rows)


def test_the_measures_of_a_known_waveform():
    t18 = np.arange(0, 60)                                           # 60 positions: harmonic k has 6 k cycles per electrical turn
    ripple = lambda t: {"psi_a": 0.0, "psi_d": 0.012, "psi_q": 0.02, "T_age": 0.8 + 0.05 * math.sin(2 * math.pi * 3 * t / 60)}
    d = frame([("noload", range(0, 360, 15), lambda t: {"psi_a": 0.012 * math.cos(math.radians(t)), "psi_d": 0.012, "psi_q": 0.0, "T_age": 0.0}),
               ("cogging", np.arange(0, 20.01, 0.5), lambda t: {"psi_a": 0, "psi_d": 0, "psi_q": 0, "T_age": 0.01 * math.sin(math.radians(18 * t))}),
               ("gamma0", range(0, 60, 10), lambda t: {"psi_a": 0, "psi_d": 0.0121, "psi_q": 0.0196, "T_age": 0.73}),
               ("gamma30", range(0, 60, 10), lambda t: {"psi_a": 0, "psi_d": 0.0064, "psi_q": 0.017, "T_age": 0.86}),
               ("gamma90", range(0, 60, 10), lambda t: {"psi_a": 0, "psi_d": 0.0008, "psi_q": 0.0, "T_age": 0.002}),
               ("ripple", range(0, 61), ripple)])
    m = mc.measures(d)
    assert m["psi_pm_mWb"] == pytest.approx(12.0)                    # the fundamental of psi_a, in mWb
    assert m["psi_d_noload_mean_mWb"] == pytest.approx(12.0)
    assert m["cogging_pp_mNm"] == pytest.approx(20.0, rel=0.02)      # 0.01 N.m peak: 20 mN.m peak to peak
    assert m["T_mean_gamma0_Nm"] == pytest.approx(0.73) and m["T_mean_gamma30_Nm"] == pytest.approx(0.86)
    assert m["psi_d_gamma90_mWb"] == pytest.approx(0.8) and m["psi_q_gamma0_mWb"] == pytest.approx(19.6)
    assert m["T_mean_mtpa_Nm"] == pytest.approx(0.8)
    assert m["ripple_pp_Nm"] == pytest.approx(0.1, rel=0.02)
    assert m["harmonic18_mNm"] == pytest.approx(50.0) and m["harmonic6_mNm"] == pytest.approx(0, abs=1e-9)     # k = 3 is the 18th


def test_old_rows_pick_exactly_the_points_that_are_solved_again(tmp_path):
    """FEMM 4.2's phase1.csv and phase2.csv hold many more rows; old_rows keeps the ones the ladder solves."""
    (tmp_path / "barrier").mkdir()
    rows1 = []
    for variant in ("with band", "no band"):
        rows1 += [dict(study="noload", variant=variant, theta_e=float(t), I=np.nan, gamma=np.nan, psi_a=1, psi_d=1, psi_q=0, T_age=0) for t in range(0, 360, 3)]
        rows1 += [dict(study="cogging", variant=variant, theta_e=float(t), I=np.nan, gamma=np.nan, psi_a=1, psi_d=1, psi_q=0, T_age=0) for t in np.arange(0, 20.01, 0.5)]
        rows1 += [dict(study="gamma", variant=variant, theta_e=float(t), I=float(i), gamma=float(g), psi_a=1, psi_d=1, psi_q=0, T_age=0)
                  for i in (5, 10) for g in (0, 7.5, 30, 90) for t in range(0, 60, 10)]
    pd.DataFrame(rows1).to_csv(tmp_path / "barrier" / "phase1.csv", index=False)
    rows2 = [dict(study="ripple", variant=v, theta_e=float(t), psi_a=1, psi_d=1, psi_q=0, T_age=0) for v in ("with band", "no band") for t in range(0, 61)]
    pd.DataFrame(rows2).to_csv(tmp_path / "barrier" / "phase2.csv", index=False)
    d = mc.old_rows(tmp_path)
    counts = d["set"].value_counts().to_dict()
    assert counts == {"noload": 24, "cogging": 41, "gamma0": 6, "gamma30": 6, "gamma90": 6, "ripple": 61}
    assert set(d.columns) == {"set", "theta", "psi_a", "psi_d", "psi_q", "T_age"}


def test_the_air_gap_deviation_is_a_share_of_the_old_waveforms_peak():
    deg = np.arange(0, 90, 0.25)
    old = pd.DataFrame({"deg": deg, "br": 0.25 * np.sin(np.radians(4 * deg))})
    new = old.copy()
    new.loc[10, "br"] += 0.0025                                       # one sample, 1 % of the peak 0.25 T
    largest, rms = mc.gap_deviation(old, new)
    assert largest == pytest.approx(1.0) and rms == pytest.approx(1.0 / math.sqrt(len(deg)))


def test_the_check_femm_ladder_prints_femm_42s_numbers_from_compare_and_writes_its_csv(tmp_path, monkeypatch, capsys):
    """The four check_femm points at every density, against FEMM 4.2's coarse numbers, which live in compare.README_CHECK_FEMM.
    (A full run once stopped here with a KeyError, after 7 minutes, because that dictionary's keys had changed.)"""
    monkeypatch.setattr(mc, "OUT", tmp_path)
    rows = [{"tag": tag, "set": f"check{k}", "psi_d": 0.012 + 0.0001 * n, "psi_q": 1e-5, "T_age": t, "T_wst": t}
            for n, tag in enumerate(("c-natural", "f-natural")) for k, t in enumerate((0.0004, 0.8087, 0.9104, 0.3944))]
    cal = pd.DataFrame([{"tag": "c-natural", "nodes": 9580}, {"tag": "f-natural", "nodes": 17934}])
    out = mc.summarise_check_femm(pd.DataFrame(rows), cal)
    text = capsys.readouterr().out
    assert "11,644 nodes: psi_d 11.941, T 0.8061 / 0.911 / 0.3886" in text                # FEMM 4.2's, as compare.py has them
    assert "9,580" in text and "17,934" in text and "0.8087" in text
    assert list(out.tag) == ["c-natural", "f-natural"] and (tmp_path / "check_femm_ladder.csv").is_file()
    assert out.loc[0, "T_age_1"] == 0.8087 and out.loc[1, "psi_d_0_mWb"] == pytest.approx(12.1)


# ----------------------------------------------------------------------------------------------- with FEMM
@pytest.mark.femm
def test_density_one_is_the_models_own_mesh_and_a_smaller_f_gives_more_nodes():
    import femm_session
    OUT.mkdir(parents=True, exist_ok=True)
    femm_session.open_femm()
    try:
        p = mg.Params()
        nodes = {f: mc.count_nodes(mc.open_model(p, f, False, OUT / f"density_{f}.fem")) for f in (1.0, 0.8166)}
    finally:
        femm_session.close_femm()
    assert nodes[1.0] == 9580                                          # what check_femm's default model has on FEMM 5.0
    assert nodes[0.8166] == pytest.approx(11722, rel=0.01)             # the calibrated factor that gives FEMM 4.2's count, within 1 %
    assert nodes[0.8166] > nodes[1.0]


@pytest.mark.femm
@pytest.mark.skipif(not (REAL_OLD / "barrier" / "fields" / "noload_with_band.ans").is_file(), reason="FEMM 4.2's old out folder is not on this PC")
def test_femm_5_reads_a_femm_4_2_solution_the_way_femm_4_2_did():
    """The post-processor of FEMM 5.0 on FEMM 4.2's own saved solutions: flux linkages and the weighted-stress-tensor torque must be
    the old numbers to rounding, the air-gap B within a tenth of a percent of its peak. (The air-gap torque integral is the one
    quantity that is not: it reads about 0.1 % higher, on the very same field. The report, section 5.2 (project note "4.2 vs 5.0"),
    says what is known of it.)"""
    df = mc.run_post42(REAL_OLD)
    assert len(df) == 6
    assert (df.psi_max_rel_diff < 1e-12).all()
    assert (df.T_wst_diff.abs() < 1e-12).all()
    assert (df.gap_br_max_diff_pct_of_peak < 0.1).all()
    loaded = df[df.T_age_old.abs() > 0.1]
    assert ((loaded.T_age_diff / loaded.T_age_old).abs() < 0.002).all()


# ----------------------------------------------------------------------------------------------- one bit of a design
def test_a_design_read_back_from_a_csv_must_be_read_exactly(tmp_path):
    """pandas' default CSV reader is one unit in the last place off for some of the optimizer's numbers, and for this model one
    unit in the last place is enough to give FEMM a different mesh (see the 'noise' part). Python's float() is exact, and so is
    pandas with float_precision="round_trip". This pins the trap down with two numbers taken from the optimizer's own CSV."""
    import io
    values = ["39.227554466990604", "1.2227404896806187"]
    path = tmp_path / "designs.csv"
    path.write_text("mag_x,web\n" + ",".join(values) + "\n", encoding="utf-8")
    rows = mc.exact_csv_rows(path)
    assert [float(rows[0]["mag_x"]), float(rows[0]["web"])] == [float(v) for v in values]
    exact = pd.read_csv(io.StringIO(path.read_text()), float_precision="round_trip").iloc[0]
    assert [exact.mag_x, exact.web] == [float(v) for v in values]
    default = pd.read_csv(io.StringIO(path.read_text())).iloc[0]
    assert [default.mag_x, default.web] != [float(v) for v in values]          # the trap: the default reader moved at least one of them


def test_one_unit_in_the_last_place_is_a_different_number_to_the_model():
    """The noise part moves a parameter by math.nextafter. That is a different float, and a valid design."""
    base = mg.Params()
    moved = mg.Params(mag_y=math.nextafter(base.mag_y, math.inf))
    assert moved.mag_y != base.mag_y and abs(moved.mag_y - base.mag_y) < 1e-14
    assert mc.NOISE_VARIANTS["mag_y + 1 ulp"]["mag_y"] == moved.mag_y
    assert len(mg.build(moved)) == len(mg.build(base))                         # the same drawing, to any eye


# ----------------------------------------------------------------------------------------------- the current, the ripple, the full model
def test_the_current_vector_is_the_studies_convention():
    """gamma is measured from the q axis toward -d: all q at 0 degrees, all -d at 90 degrees."""
    assert mc.current_vector(10.0, 0) == (-0.0, 10.0)
    i_d, i_q = mc.current_vector(20.0, 90)
    assert i_d == -20.0 and i_q == pytest.approx(0.0, abs=1e-14)
    i_d, i_q = mc.current_vector(10.0, 30)
    assert (i_d, i_q) == pytest.approx((-5.0, 10 * math.cos(math.radians(30))))
    assert math.hypot(*mc.current_vector(7.0, 41.3)) == pytest.approx(7.0)


def test_the_ripple_measures_of_a_known_period():
    t = np.arange(60)
    torque = 0.8 + 0.05 * np.sin(2 * np.pi * 3 * t / 60) + 0.002 * np.sin(2 * np.pi * 4 * t / 60)    # harmonics 18 and 24 of the electrical angle
    m = mc.ripple_measures(torque)
    assert m["T_mean_mtpa_Nm"] == pytest.approx(0.8) and m["ripple_pp_Nm"] == pytest.approx(np.ptp(torque))
    assert m["harmonic18_mNm"] == pytest.approx(50.0) and m["harmonic24_mNm"] == pytest.approx(2.0)
    assert m["harmonic6_mNm"] == pytest.approx(0.0, abs=1e-9) and set(m) >= {f"harmonic{6 * k}_mNm" for k in range(1, 7)}


def test_the_full_models_jobs_are_validate_fulls_48_solves_in_its_order():
    jobs = mc.full_model_jobs()
    assert len(jobs) == 48 and [k for k, _, _ in jobs] == list(range(48))
    assert jobs[0] == (0, 0.0, 0.0) and jobs[23] == (23, 345.0, 0.0) and jobs[24] == (24, 0.0, 10.0) and jobs[47] == (47, 345.0, 10.0)


def test_the_full_model_is_read_from_the_copy_beside_the_code():
    """validate_full.py and mesh_checks.py read one file: config.toml's reference_fem, the copy in P/reference (the reference
    commit). validate_full.py read the IPM-synRM clone on the Desktop until the tools plan's Task 2."""
    import hashlib
    import validate_full as vf
    assert mc.FULL_FEM == paths.path("reference_fem") == Path(vf.FEM)
    assert mc.FULL_FEM.is_file() and mc.FULL_FEM.parent.name == "reference"
    assert hashlib.sha256(mc.FULL_FEM.read_bytes()).hexdigest().startswith("3ac2f628")


def synthetic_cycle(shift, own_mean=0.0, n=24):
    """A cogging-like waveform with a 18th harmonic (n = 24 positions per electrical cycle), shifted by `shift` mN.m."""
    k = np.arange(n)
    return (own_mean + shift + 20.0 * np.sin(2 * np.pi * 18 * k / 360 * (360 / n))) / 1e3      # N.m


def test_a_common_offset_is_told_from_noise():
    """FEMM 4.2: a waveform with mean +1.7 mN.m. FEMM 5.0: the same waveform 2.4 mN.m lower on three meshes whose noise is +-0.3 mN.m."""
    rng = np.random.default_rng(3)
    old = synthetic_cycle(1.7)
    meshes = np.column_stack([synthetic_cycle(-0.7) + rng.uniform(-0.0003, 0.0003, 24) for _ in range(3)])
    s = mc.offset_analysis(old, meshes[:, 0], meshes)
    assert s["mean_old"] == pytest.approx(1.7, abs=1e-6) and s["means_new"][0] == pytest.approx(-0.7, abs=0.3)
    assert s["shift"] == pytest.approx(-2.4, abs=0.3) and s["n_negative"] >= 22 and s["n_positions"] == 24
    assert s["median_abs_diff"] == pytest.approx(2.4, abs=0.4) and s["noise_median"] < 0.7
    assert s["within_noise"] <= 2                                               # the offset is 4 times the noise: no position is within it
    assert s["median_abs_diff_no_offset"] < 0.4 and s["within_noise_no_offset"] >= 18      # without the offset the two agree to the noise
    assert s["spread_of_means"] < 0.3 and abs(s["shift"]) / s["spread_of_means"] > 5       # and the three meshes agree with one another


def test_the_dense_grid_separates_the_mean_from_the_sampling():
    """A 24th harmonic of amplitude A on a mean m: a 24-position grid sees m + A cos(90 degrees x shift/3.75), the 96-position mean sees m."""
    k = np.arange(96)
    a, m = 0.8, 0.3
    torque = m + a * np.cos(2 * np.pi * 24 * k / 96) + 20.0 * np.sin(2 * np.pi * 18 * k / 96)
    r = mc.dense_analysis(torque)
    assert r["mean"] == pytest.approx(m)
    assert r["grid_means"] == pytest.approx([m + a, m, m - a, m])               # shifted by 0, 3.75, 7.5, 11.25 degrees
    assert r["harmonics"][24] == pytest.approx(a) and r["harmonics"][18] == pytest.approx(20.0) and r["harmonics"][6] == pytest.approx(0.0, abs=1e-9)


def test_the_flagged_parts_read_the_old_studys_rows(tmp_path):
    (tmp_path / "barrier").mkdir()
    rows = [{"study": "ripple", "variant": v, "theta_e": t, "id": -4.0 - (v == "no band"), "iq": 8.0 + 0.5 * (v == "no band"), "T_age": 1.0 + 0.001 * t + (v == "no band")}
            for v in ("with band", "no band") for t in range(0, 61)]
    pd.DataFrame(rows).to_csv(tmp_path / "barrier" / "phase2.csv", index=False)
    assert mc.ripple_current_no_band(tmp_path) == (-5.0, 8.5)                   # the no-band row at position 0, as floats
    m = mc.sweep_ripple_no_band(tmp_path)                                       # the 61st sample (60 degrees) repeats the first and is dropped
    assert m["T_mean_mtpa_Nm"] == pytest.approx(2.0295)                         # the mean of 2.000 ... 2.059 (60 samples)
    phase1 = pd.DataFrame([{"study": "gamma", "variant": "with band", "I": 20.0, "gamma": 90.0, "theta_e": t, "T_age": 0.001 * t} for t in (0.0, 40.0)])
    phase1.to_csv(tmp_path / "barrier" / "phase1.csv", index=False)
    assert mc.sweep_state_torque(tmp_path, 20.0, 90.0, 40.0) == pytest.approx(0.04)


@pytest.mark.skipif(not (REAL_OLD / "barrier" / "phase2.csv").is_file(), reason="FEMM 4.2's old out folder is not on this PC")
def test_the_ripple_current_typed_in_mesh_checks_is_the_one_femm_42_solved():
    """slice_jobs() types FEMM 4.2's (i_d, i_q) of the with-band ripple sweep; its phase2.csv is where it comes from."""
    p2 = pd.read_csv(REAL_OLD / "barrier" / "phase2.csv", float_precision="round_trip")
    r = p2[(p2.study == "ripple") & (p2.variant == "with band")].sort_values("theta_e").iloc[0]
    assert (float(r["id"]), float(r["iq"])) == mc.IQ_ID_42
