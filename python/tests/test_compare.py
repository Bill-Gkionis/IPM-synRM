"""Tests for validation/compare.py: the rules it judges FEMM 5.0's numbers by, and the things it must notice.

No FEMM is started. The last test runs the whole comparison on this PC's real folders and is skipped where they are not there
(FEMM 4.2's old out folder is deleted when it is no longer needed). Run from P:   uv run pytest tests/test_compare.py -q
"""
import csv
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import compare as cmp_
import paths

REAL_OLD = Path(os.environ.get("FEMM42_OUT") or paths.path("femm42_out"))      # FEMM 4.2's run: config.toml, unless the variable names another
REAL_NEW = Path(__file__).resolve().parent.parent / "out"


# ----------------------------------------------------------------------------------------------- the rules
def test_a_percentage_is_new_minus_old_over_old():
    assert cmp_.plain_percent(0.8061, 0.8087) == pytest.approx(0.3225, abs=1e-3)
    assert cmp_.plain_percent(-2.0, -1.0) == pytest.approx(50.0)           # relative to the size of the old value
    assert cmp_.plain_percent(0.0, 1.0) is None                            # a percentage of nothing is not a number


def test_a_small_old_value_is_judged_against_a_tenth_of_the_scale():
    """The 'Self-check' of compare.py: psi_q at no load, 0.023 mWb, becomes 0.031 mWb; the machine's flux is 12 mWb."""
    assert cmp_.plain_percent(0.023, 0.031) == pytest.approx(34.8, abs=0.1)
    judged = cmp_.judged_percent(0.023, 0.031, "flux", scale=12.0)
    assert judged == pytest.approx(100 * 0.008 / 1.2, abs=1e-6)            # 0.008 mWb of 10 % of 12 mWb
    assert abs(judged) < cmp_.limit_for("flux", "coarse")                  # fair: well inside 2 %


def test_a_value_above_the_floor_is_judged_by_itself():
    assert cmp_.judged_percent(5.0, 5.1, "flux", scale=12.0) == pytest.approx(2.0)
    assert cmp_.judged_percent(5.0, 5.1, "flux", scale=12.0) == cmp_.plain_percent(5.0, 5.1)


def test_a_ripple_or_a_harmonic_has_no_floor_however_small_it_is():
    """0.00066 -> 0.00104 N.m is +56 %: flagged, then explained in the report. The scale must not hide it."""
    assert cmp_.judged_percent(0.00066, 0.00104, "harmonic", scale=0.05) == pytest.approx(57.6, abs=0.1)
    assert cmp_.verdict_for("harmonic", "fine", 57.6, 0.00066, 0.00104) == "FLAG"


@pytest.mark.parametrize("kind, mesh, limit", [("mean_torque", "fine", 1.0), ("mean_torque", "coarse", 2.0), ("flux", "fine", 1.0),
                                               ("flux", "coarse", 2.0), ("ripple", "fine", 10.0), ("cogging", "coarse", 10.0),
                                               ("harmonic", "fine", 10.0)])
def test_the_tolerances_are_the_plans(kind, mesh, limit):
    assert cmp_.limit_for(kind, mesh) == limit


def test_verdicts():
    assert cmp_.verdict_for("mean_torque", "fine", 0.99, 1.0, 1.0099) == "ok"
    assert cmp_.verdict_for("mean_torque", "fine", 1.01, 1.0, 1.0101) == "FLAG"
    assert cmp_.verdict_for("mean_torque", "coarse", 1.5, 1.0, 1.015) == "ok"          # 2 % on the coarse mesh
    assert cmp_.verdict_for("mean_torque", "fine", -1.5, 1.0, 0.985) == "FLAG"         # a drop counts as much as a rise
    assert cmp_.verdict_for("flux", "fine", 0.0, 3.0, 3.0) == "identical"
    assert cmp_.verdict_for("info", "fine", 17.0, 1, 2) == "info"                      # no tolerance: shown, never judged
    assert cmp_.verdict_for("exact", "fine", None, 3, 3) == "identical"
    assert cmp_.verdict_for("exact", "fine", None, 3, 4) == "DIFFERS"


# ----------------------------------------------------------------------------------------------- groups of numbers
def test_a_group_reports_median_p95_max_and_counts_what_is_beyond_the_tolerance():
    rep = cmp_.Report()
    old = np.full(100, 10.0)
    new = old * (1 + np.concatenate([np.full(97, 0.001), [0.03, 0.03, 0.05]]))       # 97 at +0.1 %, three of 3, 3 and 5 %
    row = rep.group("s", "x", old, new, "mean_torque", "coarse", scale=10.0)
    assert row.stats[0] == pytest.approx(0.1) and row.stats[2] == pytest.approx(5.0)
    assert row.stats[1] == pytest.approx(0.1, abs=0.2)                                # p95 is still among the small ones
    assert (row.n, row.n_flagged, row.verdict) == (100, 3, "FLAG")                    # 3, 3 and 5 % are beyond 2 %
    assert len(rep.detail) == 100


def test_the_floor_does_not_forgive_a_difference_that_is_large_against_a_tenth_of_the_scale():
    rep = cmp_.Report()
    old = np.array([10.0, 0.001, -10.0])
    new = np.array([10.0, 0.05, -10.0])                      # the middle one: 0.049 off. A tenth of the scale is 1.0: 4.9 % of that
    row = rep.group("s", "x", old, new, "flux", "fine", scale=10.0)
    assert (row.verdict, row.n_flagged, row.n_plain_flagged) == ("FLAG", 1, 1)


def test_the_floor_forgives_a_difference_that_is_small_against_a_tenth_of_the_scale():
    rep = cmp_.Report()
    old = np.array([10.0, 0.001, -10.0])
    new = np.array([10.0, 0.0105, -10.0])                    # 0.0095 off a value of 0.001: +950 % of itself, 0.95 % of the floor
    row = rep.group("s", "x", old, new, "flux", "fine", scale=10.0)
    assert row.verdict == "ok" and row.n_flagged == 0
    assert row.n_plain_flagged == 1                          # the number is not hidden: it is counted by the plain rule


def test_a_waveform_is_judged_against_its_own_size():
    rep = cmp_.Report()
    t = np.linspace(0, 2 * np.pi, 120, endpoint=False)
    old = np.sin(t)
    new = np.sin(t) + 0.004 * np.cos(3 * t)                  # 0.4 % of the peak everywhere, huge near the zero crossings
    row = rep.group("s", "x", old, new, "flux", "fine", scale=1.0, wave=True)
    assert row.stats[2] <= 0.4 + 1e-9 and row.verdict == "ok"


def test_equal_groups_are_identical_and_a_group_with_no_tolerance_is_only_info():
    rep = cmp_.Report()
    a = np.array([1.0, 2.0, 3.0])
    assert rep.group("s", "x", a, a.copy(), "flux", "fine", scale=3.0).verdict == "identical"
    assert rep.group("s", "y", a, a * 2, "info", "fine").verdict == "info"


def test_csv_rows_have_three_lines_for_a_group_and_the_verdict_is_on_the_max_line():
    rep = cmp_.Report()
    rep.group("s", "x", np.full(10, 1.0), np.full(10, 1.5), "mean_torque", "coarse", scale=1.0)
    rows = list(rep.csv_rows())
    assert [r[1] for r in rows] == ["x: median |delta|", "x: p95 |delta|", "x: max |delta|"]
    assert [r[5] for r in rows] == ["info", "info", "FLAG"]
    assert all(len(r) == len(cmp_.CSV_HEADER) for r in rows)


# ----------------------------------------------------------------------------------------------- small helpers
def test_harmonic_amplitudes_of_a_known_waveform():
    t = np.arange(60)
    y = 0.5 + 2.0 * np.sin(2 * np.pi * 3 * t / 60) + 0.25 * np.cos(2 * np.pi * 6 * t / 60)
    h = cmp_.harmonic_amplitudes(y)
    assert h[0] == pytest.approx(0.5) and h[3] == pytest.approx(2.0) and h[6] == pytest.approx(0.25) and h[1] == pytest.approx(0, abs=1e-12)


@pytest.mark.parametrize("name, kind", [("psi_pm_mWb", "flux"), ("Ld_10A_mH", "flux"), ("saliency_5A", "flux"), ("B1_gap_T", "flux"),
                                        ("T_mtpa_10A", "mean_torque"), ("T_q_20A", "mean_torque"), ("T_pm_10A_Nm", "mean_torque"),
                                        ("mtpa_gamma_15A", "angle"), ("cogging_pp_mNm", "cogging"), ("ripple_pp_Nm", "ripple"),
                                        ("ripple_pct", "ripple"), ("emf_thd_pct", "harmonic"), ("B_gap_thd_pct", "harmonic")])
def test_every_summary_quantity_has_its_kind(name, kind):
    assert cmp_.summary_kind(name) == kind


def test_old_logs_are_read_whatever_their_encoding(tmp_path):
    (tmp_path / "a.log").write_bytes("96 solves in 1721 s".encode("utf-16"))               # a BOM, then two bytes per character
    (tmp_path / "b.log").write_bytes(b"\xef\xbb\xbf" + "phase 1: 946 solves in 3116 s".encode("utf-8"))
    (tmp_path / "c.log").write_text("plain 7 s", encoding="utf-8")
    assert cmp_.read_text(tmp_path / "a.log") == "96 solves in 1721 s"
    assert cmp_.read_text(tmp_path / "b.log").startswith("phase 1")
    assert cmp_.read_text(tmp_path / "c.log") == "plain 7 s"


def test_check_femm_lines_are_read():
    line = "th_e=   0 id= -6 iq=  8 | psi_d=  5.188 psi_q= 15.883 mWb | T_age= 0.9104 T_wst= 0.9100 T_dq= 0.8208 N.m | 0.3 s"
    m = cmp_.CHECK_FEMM_LINE.search(line)
    assert (int(m.group(1)), int(m.group(2))) == (-6, 8)
    assert [float(x) for x in m.groups()[2:]] == [5.188, 15.883, 0.9104, 0.9100, 0.8208]


def test_the_check_femm_log_is_read_by_one_function_and_every_quantity_has_a_name():
    lines = ["# run_all header",
             "th_e=   0 id=  0 iq=  0 | psi_d= 11.991 psi_q=  0.031 mWb | T_age= 0.0000 T_wst= 0.0000 T_dq= 0.0000 N.m | 0.3 s",
             "th_e=   0 id=  0 iq= 10 | psi_d= 12.000 psi_q= 19.500 mWb | T_age= 0.8087 T_wst= 0.8085 T_dq= 0.8081 N.m | 0.3 s",
             "something else"]
    got = cmp_.parse_check_femm(lines)
    assert set(got) == {(0, 0), (0, 10)}
    assert got[(0, 10)] == {"psi_d": 12.0, "psi_q": 19.5, "T_age": 0.8087, "T_wst": 0.8085, "T_dq": 0.8081}
    assert list(cmp_.CHECK_FEMM_UNITS) == list(got[(0, 10)])                  # the units are listed in the order of the log line
    assert cmp_.parse_check_femm([]) == {}
    assert {k for point in cmp_.README_CHECK_FEMM.values() for k in point} <= set(cmp_.CHECK_FEMM_UNITS)


def test_a_log_header_and_footer_are_not_part_of_the_output(tmp_path):
    folder = tmp_path / "validation" / "logs"
    folder.mkdir(parents=True)
    (folder / "03-check_femm.log").write_text("# run_all step 3\n####\nthe numbers\n# finished: 0\n", encoding="utf-8")
    assert cmp_.log_lines(cmp_.find_log(tmp_path, "check_femm")) == ["the numbers"]
    assert cmp_.find_log(tmp_path, "nothing") is None


# ----------------------------------------------------------------------------------------------- repeated runs
def write_logs(folder, name, prev, cur):
    folder.mkdir(parents=True, exist_ok=True)
    (folder / f"{name}.log.prev").write_text("# header\n" + prev + "\n# finished\n", encoding="utf-8")
    (folder / f"{name}.log").write_text("# header\n" + cur + "\n# finished\n", encoding="utf-8")


def test_two_runs_that_differ_only_in_times_and_process_ids_are_the_same_run(tmp_path):
    logs = tmp_path / "validation" / "logs"
    write_logs(logs, "03-a", "T 0.8087 N.m | 0.3 s\n2 workers: 7.1 s, speed-up 1.92x\nworker PIDs: [1, 2]\nthread finished, PID 17\ncheck 27.7 ms",
               "T 0.8087 N.m | 0.4 s\n2 workers: 7.3 s, speed-up 1.89x\nworker PIDs: [9, 8]\nthread finished, PID 99\ncheck 27.2 ms")
    assert cmp_.same_logs(tmp_path) == (1, 0)


def test_two_runs_that_print_different_numbers_are_counted(tmp_path):
    logs = tmp_path / "validation" / "logs"
    write_logs(logs, "03-a", "T 0.8087 N.m | 0.3 s", "T 0.8088 N.m | 0.3 s")
    write_logs(logs, "04-b", "same 1.0 s", "same 1.1 s")
    assert cmp_.same_logs(tmp_path) == (2, 1)


# ----------------------------------------------------------------------------------------------- the slice against the full model
def test_the_slice_against_the_full_model():
    rows = []
    for th in range(0, 360, 15):
        full = 10 * math.sin(math.radians(th))
        rows.append({"model": "hand-drawn 360", "theta_e": th, "I": 10.0, "psi_a": full, "T": full})
        rows.append({"model": "motorgeo 90", "theta_e": th, "I": 10.0, "psi_a": 0.99 * full, "T": 0.98 * full})
        rows.append({"model": "hand-drawn 360", "theta_e": th, "I": 0.0, "psi_a": full, "T": 0.0})
        rows.append({"model": "motorgeo 90", "theta_e": th, "I": 0.0, "psi_a": 0.995 * full, "T": 0.0})
    got = cmp_.slice_vs_full(pd.DataFrame(rows))
    assert got["flux, 10 A"] == pytest.approx(1.0, abs=1e-6)             # worst gap: 1 % of the full model's peak
    assert got["flux, no load"] == pytest.approx(0.5, abs=1e-6)
    assert got["peak torque, 10 A"] == pytest.approx(-2.0, abs=1e-6)     # the slice's peak is 2 % lower


def write_validate(folder, flux_gap, torque_gap):
    """validate/validate.csv with a slice that is `flux_gap` (a share) low in flux linkage and `torque_gap` low in torque."""
    rows = []
    for current in (0.0, 10.0):
        for th in range(0, 360, 15):
            full = 10 * math.sin(math.radians(th))
            rows.append({"model": "hand-drawn 360", "theta_e": th, "I": current, "psi_a": full, "T": full})
            rows.append({"model": "motorgeo 90", "theta_e": th, "I": current, "psi_a": (1 - flux_gap) * full, "T": (1 - torque_gap) * full})
    (folder / "validate").mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(folder / "validate" / "validate.csv", index=False)


def slice_checks(old, new):
    rep = cmp_.Report()
    cmp_.study_self(rep, old, new)
    return [r for r in rep.rows if r.quantity.startswith("the slice against the full model")]


def test_the_slice_checks_show_the_limit_that_their_verdict_used(tmp_path):
    """'1.1 % / 1.5 %, or better': with FEMM 4.2's own figures (1.15 % and 1.48 % here) the slice must be no worse than they are,
    and the limit_pct column says so; without them the plan's 1.1 and 1.5 are the limits."""
    write_validate(tmp_path / "old", 0.0115, 0.0148)
    write_validate(tmp_path / "new", 0.0113, 0.0143)
    flux, peak = slice_checks(tmp_path / "old", tmp_path / "new")
    assert (flux.verdict, peak.verdict) == ("ok", "ok")
    assert flux.limit_pct == pytest.approx(1.15) and peak.limit_pct == pytest.approx(1.48)
    assert (flux.femm42, peak.femm42) == (pytest.approx(1.15), pytest.approx(1.48))          # the limit is the FEMM 4.2 column
    write_validate(tmp_path / "new", 0.0117, 0.0149)                                         # worse than FEMM 4.2 (the torque is still under 1.5)
    assert [r.verdict for r in slice_checks(tmp_path / "old", tmp_path / "new")] == ["FAIL", "FAIL"]
    flux, peak = slice_checks(tmp_path / "nothing-old", tmp_path / "new")                    # no FEMM 4.2 numbers: the plan's
    assert (flux.limit_pct, peak.limit_pct) == (1.1, 1.5) and (flux.verdict, peak.verdict) == ("FAIL", "ok")


# ----------------------------------------------------------------------------------------------- the audit of the two extra rules
def test_the_audit_names_the_rows_that_pass_only_through_the_floor_or_the_waveform_rule(monkeypatch):
    t = np.linspace(0, 2 * np.pi, 120, endpoint=False)

    def fake_compare_all(old, new, page=None):
        rep = cmp_.Report()
        rep.value("s", "small flux", 0.02, 0.03, "flux", "fine", scale=12.0)               # +50 % of itself, 0.83 % of a tenth of 12
        rep.value("s", "small flux, close", 0.02, 0.03188, "flux", "fine", scale=12.0)     # 0.99 % of a limit of 1 %
        rep.value("s", "big flux", 5.0, 5.01, "flux", "fine", scale=12.0)                  # ok with or without the floor
        rep.value("s", "flag until 0.15", 0.02, 0.035, "flux", "fine", scale=12.0)         # FLAG (1.25 % of the floor), until the floor is 0.15
        rep.group("s", "wave", np.sin(t), np.sin(t) + 0.004 * np.cos(3 * t), "flux", "fine", scale=1.0, wave=True)
        rep.group("s", "both", np.array([10.0, 0.001, -10.0, 10.0]), np.array([10.5, 0.0105, -10.0, 10.0]), "flux", "fine", scale=10.0)
        return rep

    monkeypatch.setattr(cmp_, "compare_all", fake_compare_all)
    lines = []
    got = cmp_.audit(Path("old"), Path("new"), out=lines.append)
    assert [x["quantity"] for x in got["floor"]] == ["small flux", "small flux, close"]
    first, second = got["floor"]
    assert first["margin"] == pytest.approx(100 * (1 - 0.0100 / 1.2 * 100), abs=0.1) and not first["close"]
    assert second["close"] and second["margin"] < cmp_.NEAR
    assert first["plain"] == pytest.approx(50.0)
    assert [(x["quantity"], x["n_beyond"], x["n_beyond_plain"]) for x in got["more"]] == [("both", 1, 2)]      # the floor hides one member
    assert [x["quantity"] for x in got["wave"]] == ["wave"]
    assert "2 rows pass only through the floor" in "\n".join(lines) and "1 of them are within 4 %" in "\n".join(lines)
    assert cmp_.FLOOR == 0.10 and cmp_.WAVE_RULE is True                                      # the rules are put back
    # the floor's own value: "small flux, close" (0.0119 on a floor of 1.2 = 0.99 %) is FLAG as soon as the floor drops to 0.09 (1.1 %), and
    # "flag until 0.15" (0.015 on a floor of 1.2 = 1.25 %) stops being FLAG when the floor rises to 0.15 (0.83 %)
    assert set(got["floors"]) == set(cmp_.FLOORS_TRIED)
    assert got["floors"][0.09]["became"] == ["s | small flux, close"] and got["floors"][0.09]["stopped"] == []
    assert got["floors"][0.12]["became"] == [] and got["floors"][0.12]["stopped"] == []
    assert got["floors"][0.15]["became"] == [] and got["floors"][0.15]["stopped"] == ["s | flag until 0.15"]
    assert any("FLOOR = 0.09: 3 FLAG rows" in ln for ln in lines) and any("becomes FLAG: s | small flux, close" in ln for ln in lines)


def test_the_rules_are_put_back_even_when_the_block_fails():
    with pytest.raises(RuntimeError):
        with cmp_.rules(floor=0.0, wave=False):
            assert cmp_.FLOOR == 0.0 and cmp_.WAVE_RULE is False
            raise RuntimeError("boom")
    assert cmp_.FLOOR == 0.10 and cmp_.WAVE_RULE is True


# ----------------------------------------------------------------------------------------------- whole studies on made-up folders
def make_scaling(folder, rows):
    (folder / "scaling").mkdir(parents=True, exist_ok=True)
    (folder / "scaling" / "scaling.csv").write_text("workers,wall_s,processes,max_diff_Nm\n" + "".join(f"{a},{b},{c},{d}\n" for a, b, c, d in rows))


def test_worker_scaling_must_be_bit_identical_on_the_new_side(tmp_path):
    make_scaling(tmp_path / "old", [(1, 94.8, 1, 0), (2, 49.8, 2, 0)])
    make_scaling(tmp_path / "new", [(1, 13.5, 1, 0), (2, 7.1, 2, 1e-12)])           # one count differs in the 12th decimal
    rep = cmp_.Report()
    cmp_.study_scaling(rep, tmp_path / "old", tmp_path / "new")
    check = next(r for r in rep.rows if r.kind == "must hold")
    assert check.verdict == "FAIL"
    make_scaling(tmp_path / "new", [(1, 13.5, 1, 0), (2, 7.1, 2, 0)])
    rep = cmp_.Report()
    cmp_.study_scaling(rep, tmp_path / "old", tmp_path / "new")
    assert next(r for r in rep.rows if r.kind == "must hold").verdict == "ok"


def make_opt(folder, rows):
    (folder / "opt").mkdir(parents=True, exist_ok=True)
    cols = "gen,mag_y,mag_x,v_half_deg,web,feasible,reason,gamma,T_mean,ripple,centre_rib,q_rib,outer_bridge,web_mm,seconds\n"
    (folder / "opt" / "evaluations.csv").write_text(cols + "".join(",".join(str(v) for v in r) + "\n" for r in rows))


def opt_row(gen, y, T, ripple=0.15, gamma=25.0):
    return (gen, y, 39.4, 60.0, 1.5, True, "", gamma, T, ripple, 1.5, 1.8, 1.4, 1.46, 10.0)


def test_the_optimizer_compares_generations_0_and_1_row_by_row_and_not_the_rest(tmp_path):
    old = [opt_row(0, 6.0, 0.85), opt_row(0, 5.5, 1.08), opt_row(1, 5.6, 1.05), opt_row(2, 5.4, 1.03)]
    new = [opt_row(0, 6.0, 0.8461), opt_row(0, 5.5, 1.0821), opt_row(1, 5.6, 1.0561), opt_row(2, 5.9, 1.09)]    # generation 2: other designs
    make_opt(tmp_path / "old", old)
    make_opt(tmp_path / "new", new)
    rep = cmp_.Report()
    cmp_.study_optimize(rep, tmp_path / "old", tmp_path / "new")
    rows = [r for r in rep.rows if r.stats is not None]
    assert {r.study for r in rows} == {"optimize_demo: generation 0, row by row", "optimize_demo: generation 1, row by row"}
    later = [r for r in rep.rows if "generation 2: the same" in r.quantity]
    assert later and later[0].note.startswith("no:")                                  # the searches parted ways at generation 2


def test_a_design_that_differs_in_generation_0_is_reported_not_compared(tmp_path):
    make_opt(tmp_path / "old", [opt_row(0, 6.0, 0.85)])
    make_opt(tmp_path / "new", [opt_row(0, 6.1, 0.85)])
    rep = cmp_.Report()
    cmp_.study_optimize(rep, tmp_path / "old", tmp_path / "new")
    first = next(r for r in rep.rows if "generation 0" in r.study)
    assert first.verdict == "DIFFERS"
    assert not [r for r in rep.rows if r.stats is not None and "generation 0" in r.study]


# ----------------------------------------------------------------------------------------------- the whole thing
def test_empty_folders_give_missing_rows_a_csv_and_a_failing_exit_code_not_a_crash(tmp_path):
    (tmp_path / "old").mkdir()
    (tmp_path / "new").mkdir()
    out_csv, detail = tmp_path / "t.csv", tmp_path / "d.csv"
    code = cmp_.main([str(tmp_path / "old"), str(tmp_path / "new"), "--csv", str(out_csv), "--detail", str(detail),
                      "--page", str(tmp_path / "no-page.html"), "--quiet"])
    assert code == 1
    rows = list(csv.reader(out_csv.open(encoding="utf-8")))
    assert rows[0] == cmp_.CSV_HEADER and any(r[5] == "MISSING" for r in rows[1:])
    assert detail.is_file()


def test_a_folder_that_is_not_there_is_a_usage_error(tmp_path, capsys):
    assert cmp_.main([str(tmp_path / "nope"), str(tmp_path)]) == 2
    assert "OLD_OUT is not a folder" in capsys.readouterr().err


@pytest.mark.skipif(not (REAL_OLD / "barrier" / "phase1.csv").is_file() or not (REAL_NEW / "barrier" / "phase1.csv").is_file(),
                    reason="FEMM 4.2's old out folder or FEMM 5.0's out folder is not on this PC")
def test_the_real_comparison_on_this_pc(tmp_path):
    rep = cmp_.compare_all(REAL_OLD, REAL_NEW)
    by = {(r.study, r.quantity): r for r in rep.rows}
    t = by[("barrier: summary.csv", "with band: T_mtpa_10A")]
    assert (t.femm42, t.verdict) == (pytest.approx(0.857377, abs=1e-5), "ok") and t.delta_pct == pytest.approx(0.217, abs=0.01)
    assert by[("barrier: phase 2 (ripple at the best current angle, snapshots)", "with band: torque harmonic 24 per electrical cycle [N.m]")].verdict == "FLAG"
    assert by[("worker_scaling", "largest |T - T(1 worker)| over every worker count [N.m]")].verdict == "ok"
    nodes = by[("mesh: nodes and triangles", "nodes, coarse mesh, the default model")]
    assert (nodes.femm42, nodes.femm50) == (11644, 9580)
    nonfemm = [r for r in rep.rows if r.study == "not FEMM: identical" and r.verdict not in ("identical",)]
    assert [r.quantity for r in nonfemm] == ["file check_vs_fem.png (pixels that differ)"]           # the one picture that is not pixel-identical
    assert not [r for r in rep.rows if r.verdict in ("FAIL", "MISSING")], [(r.study, r.quantity) for r in rep.rows if r.verdict in ("FAIL", "MISSING")]
