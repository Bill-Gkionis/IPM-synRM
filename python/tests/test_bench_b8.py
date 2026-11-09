"""Tests for B8 (bench/b8_matlab.py and bench/b8_mtpa_map.m): Python -> MATLAB (matlab -batch): a flux map from FEMM, interpolated in MATLAB, the MTPA angle read back. Everything
that needs neither FEMM nor MATLAB: the brief's half degree, the map's grid and the circle in it, Python's own search against a closed form, the file the script reads, what the
script may do, the command line, the hashes and the prefdir listing, the MATLAB processes of the job, who holds the slot, the checks.

Run from P:   uv run pytest tests/test_bench_b8.py -q      (a few seconds; MATLAB is not started, no slot is taken)

In groups:
  1. The brief's numbers: 0.5 degrees, a map of at least 5 x 5 over id -10..0 and iq 0..10 A, the 10 A circle inside it.
  2. Python's side: the linear search on a map with a known answer (the closed form), the layout of the map, the spline, the torque curve.
  3. The file Python writes and the script reads, and what the script is allowed to do.
  4. The MATLAB command, the environment (MATLAB_PREFDIR), the processes of the job and who holds the slot.
  5. No project file and no MATLAB preference changes: hashes, the prefdir listing (a file that cannot be read is said, not raised), the scan of what MATLAB wrote elsewhere.
  6. The checks (fix round 1 of Task 13: the job holds every process at its class, so the graphics host gates like the rest).
"""
import json
import math
import os
import re
import time
from pathlib import Path

import numpy as np
import pytest
import scipy.io

from bench_test_tools import P


# ----------------------------------------------------------------------------------------------- 1. the brief's numbers
def test_the_briefs_half_degree_is_pinned_and_it_is_matlab_against_python():
    """MATLAB's MTPA angle equals Python's own (scipy, the same map) within 0.5 degrees. On the tolerance passes ('within'); a hundredth beyond it fails. Change the constant and
    this goes red."""
    from bench import b8_matlab as b8
    assert b8.ANGLE_TOLERANCE_DEG == 0.5
    assert b8.agree(23.0, 23.5) and b8.agree(23.5, 23.0) and b8.agree(23.0, 23.0)
    assert not b8.agree(23.0, 23.51) and not b8.agree(23.51, 23.0) and not b8.agree(23.0, 22.49)
    assert not b8.agree(float("nan"), 23.0) and not b8.agree(23.0, float("inf"))


def test_the_map_is_at_least_five_by_five_over_id_minus_ten_to_zero_and_iq_zero_to_ten_and_holds_the_ten_amp_circle():
    from bench import b8_matlab as b8
    id_axis, iq_axis = b8.grid_axes()
    assert len(id_axis) >= 5 and len(iq_axis) >= 5 and b8.GRID_POINTS >= 5
    assert id_axis[0] == -10.0 and id_axis[-1] == 0.0 and iq_axis[0] == 0.0 and iq_axis[-1] == 10.0
    assert np.all(np.diff(id_axis) > 0) and np.all(np.diff(iq_axis) > 0)
    gamma = b8.gamma_grid()
    id_c, iq_c = b8.circle(gamma)
    assert id_c.min() >= id_axis[0] and id_c.max() <= id_axis[-1] and iq_c.min() >= iq_axis[0] and iq_c.max() <= iq_axis[-1]          # the circle is inside the map
    assert np.allclose(np.hypot(id_c, iq_c), b8.CURRENT_A, rtol=1e-12) and b8.CURRENT_A == 10.0
    assert (id_c[0], iq_c[0]) == (-0.0, 10.0) and id_c[-1] == -10.0 and iq_c[-1] == pytest.approx(0.0, abs=1e-14)                  # gamma 0: all q current; 90: all -d


def test_the_angles_searched_are_the_quarter_circle_in_hundredths_of_a_degree():
    from bench import b8_matlab as b8
    gamma = b8.gamma_grid()
    assert len(gamma) == 9001 and gamma[0] == 0.0 and gamma[-1] == 90.0 and np.allclose(np.diff(gamma), 0.01) and b8.GAMMA_STEP_DEG == 0.01


# ----------------------------------------------------------------------------------------------- 2. Python's side
def test_pythons_linear_search_finds_the_closed_form_mtpa_angle_of_a_map_with_a_known_answer():
    """T = 1.5 p I (psi_m cos g + (Lq - Ld) I sin g cos g): the maximum is at sin g = (-psi_m + sqrt(psi_m^2 + 8 D^2)) / (4 D), D = (Lq - Ld) I. A map with psi_d = psi_m + Ld id
    and psi_q = Lq iq is bilinear, so the linear interpolation is exact and the search finds the closed form to half a step of the angles (0.005 degrees)."""
    from bench import b8_matlab as b8
    id_axis, iq_axis = b8.grid_axes()
    psi_d, psi_q = b8.analytic_map(id_axis, iq_axis, b8.CONTROL_PSI_M, b8.CONTROL_LD, b8.CONTROL_LQ)
    gamma = b8.gamma_grid()
    angle, torque = b8.mtpa(gamma, b8.torque_curve(id_axis, iq_axis, psi_d, psi_q, gamma, b8.CURRENT_A, 4.0))
    closed = b8.analytic_mtpa_deg(b8.CONTROL_PSI_M, b8.CONTROL_LD, b8.CONTROL_LQ)
    assert closed == pytest.approx(22.52, abs=0.01) and abs(angle - closed) <= 0.0051
    assert torque == pytest.approx(1.5 * 4 * 10 * (b8.CONTROL_PSI_M * math.cos(math.radians(closed)) + (b8.CONTROL_LQ - b8.CONTROL_LD) * 10 * math.sin(math.radians(closed))
                                                  * math.cos(math.radians(closed))), rel=1e-7)           # the angle sampled is up to 0.005 degrees from the exact maximum, a flat top: less than 1e-7


def test_the_closed_form_is_the_maximum_of_the_torque_and_the_search_would_notice_a_map_read_the_wrong_way():
    """The control exists to catch a map read with rows and columns swapped. Here: read the transposed map and the angle is far from the closed form."""
    from bench import b8_matlab as b8
    id_axis, iq_axis = b8.grid_axes()
    psi_d, psi_q = b8.analytic_map(id_axis, iq_axis, b8.CONTROL_PSI_M, b8.CONTROL_LD, b8.CONTROL_LQ)
    gamma = b8.gamma_grid()
    closed = b8.analytic_mtpa_deg(b8.CONTROL_PSI_M, b8.CONTROL_LD, b8.CONTROL_LQ)
    # the map read the wrong way: a file whose rows are taken for columns is the transposed map under the same axes
    angle_swapped, _ = b8.mtpa(gamma, b8.torque_curve(id_axis, iq_axis, psi_d.T, psi_q.T, gamma, b8.CURRENT_A, 4.0))
    assert abs(angle_swapped - closed) > 5.0, angle_swapped
    brute = np.linspace(0, 90, 90001)
    t = np.cos(np.radians(brute)) * b8.CONTROL_PSI_M + (b8.CONTROL_LQ - b8.CONTROL_LD) * 10 * np.sin(np.radians(brute)) * np.cos(np.radians(brute))
    assert brute[np.argmax(t)] == pytest.approx(closed, abs=2e-3)


def test_scipys_linear_method_reproduces_a_bilinear_function_and_the_spline_is_a_different_interpolation():
    from scipy.interpolate import RegularGridInterpolator
    from bench import b8_matlab as b8
    id_axis, iq_axis = b8.grid_axes()
    f = lambda a, b: 3.0 + 2.0 * a - 5.0 * b + 0.7 * a * b                                # bilinear: the linear interpolation is exact
    values = f(id_axis[:, None], iq_axis[None, :])
    point = np.array([[-3.7, 6.2]])
    assert RegularGridInterpolator((id_axis, iq_axis), values, method="linear")(point)[0] == pytest.approx(f(-3.7, 6.2), rel=1e-12)
    assert b8.INTERPOLATION == "linear" and b8.SPLINE == "cubic"
    curved = np.exp(id_axis[:, None] / 6.0) * (1 + iq_axis[None, :] ** 2 / 50.0)
    gamma = b8.gamma_grid()
    lin = b8.torque_curve(id_axis, iq_axis, curved, curved, gamma, 10.0, 4.0, method="linear")
    cub = b8.torque_curve(id_axis, iq_axis, curved, curved, gamma, 10.0, 4.0, method=b8.SPLINE)
    assert lin.shape == cub.shape == gamma.shape and not np.allclose(lin, cub, rtol=1e-6)


def test_the_torque_of_a_flux_map_is_one_and_a_half_p_psi_d_iq_minus_psi_q_id():
    from bench import b8_matlab as b8
    id_axis, iq_axis = np.array([-10.0, -5.0, 0.0]), np.array([0.0, 5.0, 10.0])
    psi_d = np.full((3, 3), 0.012)
    psi_q = np.full((3, 3), 0.02)
    gamma = np.array([0.0, 30.0, 90.0])
    expected = 1.5 * 4 * (0.012 * 10 * np.cos(np.radians(gamma)) - 0.02 * (-10 * np.sin(np.radians(gamma))))
    assert b8.torque_curve(id_axis, iq_axis, psi_d, psi_q, gamma, 10.0, 4.0) == pytest.approx(expected, rel=1e-12)
    assert b8.mtpa(np.array([0.0, 1.0, 2.0]), np.array([1.0, 3.0, 2.0])) == (1.0, 3.0)
    assert b8.mtpa(np.array([0.0, 1.0, 2.0]), np.array([3.0, 3.0, 2.0])) == (0.0, 3.0)                  # a tie: the first, as MATLAB's max does


# ----------------------------------------------------------------------------------------------- 3. the file and the script
def test_the_flux_map_file_has_the_names_and_the_layout_the_matlab_script_reads(tmp_path):
    """psi_d_Wb and psi_q_Wb are (id x iq): one row per id, one column per iq, the layout that interp2(iq, id, V, ...) wants. The axes are rows; the scalars are 1 x 1."""
    from bench import b8_matlab as b8
    id_axis, iq_axis = b8.grid_axes()
    psi_d = np.arange(36.0).reshape(6, 6)
    psi_q = -psi_d
    gamma = b8.gamma_grid()
    b8.write_flux_map(tmp_path / "flux_map.mat", id_axis, iq_axis, psi_d, psi_q, psi_d + 1, psi_q + 1, 4, 10.0, gamma)
    back = scipy.io.loadmat(tmp_path / "flux_map.mat")
    assert set(b8.MAT_VARIABLES) <= set(back)
    assert back["psi_d_Wb"].shape == (6, 6) and np.array_equal(back["psi_d_Wb"], psi_d) and back["id_A"].shape == (1, 6) and back["iq_A"].shape == (1, 6)
    assert back["gamma_deg"].shape == (1, 9001) and back["pole_pairs"].item() == 4.0 and back["current_A"].item() == 10.0
    assert np.array_equal(back["ctrl_psi_d_Wb"], psi_d + 1) and np.array_equal(back["ctrl_psi_q_Wb"], psi_q + 1)


def test_an_answer_file_that_cannot_be_read_is_an_error_that_names_it_and_not_a_traceback(tmp_path, monkeypatch):
    """MATLAB can leave a half-written or damaged answer.mat (or a locked one). read_answer says what is wrong with it, as it does for a missing file or a key (review M4,
    Task 12's lesson 3): a FAIL naming the step, never a scipy exception out of the case."""
    import zlib

    from bench import b8_matlab as b8
    answer = tmp_path / "answer.mat"
    for bad in (b"", os.urandom(300), b"not a MAT file at all\n" * 20):                                                # empty (truncated), random bytes, text
        answer.write_bytes(bad)
        said = b8.read_answer(tmp_path)
        assert said["error"] and "answer.mat" in said["error"] and "cannot be read" in said["error"], (bad[:10], said)
    good = tmp_path / "good.mat"
    scipy.io.savemat(good, {"a": np.arange(200.0)}, do_compression=True)
    data = good.read_bytes()
    answer.write_bytes(data[: len(data) // 2])                                                                        # cut in the middle: OSError from scipy
    assert "cannot be read" in b8.read_answer(tmp_path)["error"]
    for error in (PermissionError("locked by another program"), zlib.error("corrupt compressed data"), NotImplementedError("HDF reader needed")):
        def refuse(*args, error=error, **kwargs):
            raise error
        monkeypatch.setattr(scipy.io, "loadmat", refuse)
        said = b8.read_answer(tmp_path)
        assert "answer.mat cannot be read" in said["error"] and type(error).__name__ in said["error"]
        monkeypatch.undo()
    scipy.io.savemat(answer, {"gamma_deg": np.linspace(0, 90, 5), "torque_Nm": np.ones(5), "mtpa_deg": "twenty", "mtpa_Nm": 1.1, "ctrl_torque_Nm": np.ones(5), "ctrl_mtpa_deg": 22.5,
                              "ctrl_mtpa_Nm": 1.0})
    (tmp_path / "answer.json").write_text("{}", encoding="utf-8")
    assert "mtpa_deg" in b8.read_answer(tmp_path)["error"] and "not a number" in b8.read_answer(tmp_path)["error"]     # a variable that is text, not a number


def test_the_answer_is_read_back_from_matlabs_two_files_and_a_missing_one_is_said_not_raised(tmp_path):
    from bench import b8_matlab as b8
    assert "answer.mat" in b8.read_answer(tmp_path)["error"]                                                          # nothing there: an error text, not an exception
    gamma = np.linspace(0, 90, 5)
    scipy.io.savemat(tmp_path / "answer.mat", {"gamma_deg": gamma, "torque_Nm": gamma * 0 + 1.0, "mtpa_deg": 22.5, "mtpa_Nm": 1.1, "ctrl_torque_Nm": gamma * 0 + 2.0,
                                               "ctrl_mtpa_deg": 22.52, "ctrl_mtpa_Nm": 1.0})
    assert "answer.json" in b8.read_answer(tmp_path)["error"]
    (tmp_path / "answer.json").write_text(json.dumps({"matlab_release": "2025b", "prefdir": "C:\\x", "max_comp_threads": 1}), encoding="utf-8")
    answer = b8.read_answer(tmp_path)
    assert answer["error"] is None and answer["mtpa_deg"] == 22.5 and answer["ctrl_mtpa_deg"] == 22.52 and answer["info"]["max_comp_threads"] == 1
    assert np.array_equal(answer["gamma_deg"], gamma) and answer["torque_Nm"].shape == (5,)


def test_the_matlab_script_reads_what_python_writes_interpolates_linearly_and_writes_only_where_it_may():
    """A static look at bench/b8_mtpa_map.m (MATLAB is not started): it asks for every variable Python writes; it interpolates with interp2 'linear', the method the case names; it
    writes only answer.mat, answer.json and figure.png, in the current folder; it deletes nothing, runs no other program, and changes no path."""
    from bench import b8_matlab as b8
    text = (P / "bench" / "b8_mtpa_map.m").read_text(encoding="utf-8")
    code = b8.strip_matlab_comments(text)                                                                               # the code without its comments (string-aware: see the test below)
    for name in b8.MAT_VARIABLES:
        assert f"'{name}'" in code, name
    assert "interp2(" in code and "'linear'" in code and b8.INTERPOLATION == "linear"
    written = set(re.findall(r"fopen\('([^']+)', *'w'", code)) | set(re.findall(r"save\('([^']+)'", code)) | set(re.findall(r"exportgraphics\(\w+, '([^']+)'", code))
    assert written == {"answer.json", "answer.mat", "figure.png"}, written
    assert b8.forbidden_in_script(code) == [], b8.forbidden_in_script(code)                                            # nothing outside the folder, nothing permanent (see the next test)
    assert "flux_map.mat" in code and "load('flux_map.mat')" in code
    assert "sprintf('MATLAB R%s" in text and "version('-release')" in text                                              # the figure's title says MATLAB R2025b, MathWorks' own name (the release is '2025b')
    # The figure goes to the project notes' MATLAB chapter, a white page: it gets MATLAB's light theme, for this figure only (theme(fig, "light"): R2025a and later; MathWorks' page for
    # it is cited in the script). A global setting would be `settings`, which stays forbidden.
    assert 'theme(fig, "light")' in code and 'isMATLABReleaseOlderThan("R2025a")' in code
    assert code.index('theme(fig, "light")') < code.index("tiledlayout(") < code.index("exportgraphics(")                 # before anything is drawn into the figure (the layout comes first), and so before it is exported
    assert "https://www.mathworks.com/help/matlab/ref/theme.html" in text and "R2025a" in text and not re.search(r"\bsettings\b", code)
    assert b8.SCRIPT_NAME == "b8_mtpa_map" and re.fullmatch(r"[A-Za-z]\w*", b8.SCRIPT_NAME)                                # a name MATLAB can run as a script
    assert "-v7" in code                                                                                                # answer.mat is a format scipy.io.loadmat reads


def test_the_script_may_not_change_matlabs_path_preferences_or_settings_for_good_and_the_words_are_whole_words():
    """'No permanent change' (review M9): savepath, setpref, rmpref, settings and path( are as forbidden as deleting a file. They are whole words, so that fullpath( or a name that
    merely ends in 'settings' does not trip the test, and the scan sees the code with its comments cut, so a comment that names them does not either."""
    from bench import b8_matlab as b8
    for hit, word in (("savepath", r"\bsavepath\b"), ("setpref('a', 'b', 1)", r"\bsetpref\b"), ("rmpref('a')", r"\brmpref\b"), ("s = settings;", r"\bsettings\b"),
                      ("settings.matlab.general", r"\bsettings\b"), ("p = path();", r"\bpath\("), ("path(path, 'x')", r"\bpath\(")):
        assert word in b8.forbidden_in_script(hit), (hit, b8.forbidden_in_script(hit))
    for clean in ("f = fullpath(x);", "mypath(1)", "ctrl_settings = 3;", "pathway = 2;", "usersettings(1)", "x = myrmpref(2);", "filepath(3)"):
        assert b8.forbidden_in_script(clean) == [], (clean, b8.forbidden_in_script(clean))
    assert "delete(" in b8.forbidden_in_script("delete('x')") and "!" in b8.forbidden_in_script("!dir")                   # the substrings that were there before are still there
    code = b8.strip_matlab_comments("x = 1;  % savepath would be bad here\n% setpref('a','b',1)\ny = 2;")
    assert b8.forbidden_in_script(code) == []


def test_the_comments_of_a_matlab_script_are_cut_with_its_strings_in_mind():
    """The scan cut every line at its first % (re-review): a %s inside a string hid the rest of the line, a forbidden word included. A % inside a '...' or "..." literal starts no
    comment. A ' right after a letter, a digit, _, ), ], }, . or another ' is the transpose operator and starts no string. Strings stay in what is scanned (a forbidden word in one is
    still seen, which is the careful side); comments go, and so do block comments (%{ ... %}) and what follows a continuation (...)."""
    from bench import b8_matlab as b8
    savepath = r"\bsavepath\b"

    def trips(code: str) -> list[str]:
        return b8.forbidden_in_script(b8.strip_matlab_comments(code))
    assert trips(r"fprintf('%s\n', x); savepath") == [savepath]                                                         # the old cut hid it behind the %s
    assert trips(r"x = y'; % savepath") == []                                                                           # a transpose, then a comment
    for hit in ("s = sprintf('100%% of %s', a); savepath", 'disp("50% done"); savepath', "disp('it''s % not a comment'); savepath", 's = "it\'s"; savepath',
                "z = (a+b)'; savepath", "w = [1 2]'; savepath", "c = d{1}'; savepath", "t = x.'; savepath", "n = 3'; savepath", "v = a_'; savepath", "u = x''; savepath",
                "a = 'x'; savepath", "b = {'a' 'b'}; savepath", "disp('savepath')", "%{\nx\n%}\nsavepath"):
        assert trips(hit) == [savepath], hit
    for clean in ("z = (a+b)'; % settings", "w = [1 2]'; % savepath", "c = d{1}'; % setpref", "t = x.'; % savepath", "u = x''; % savepath", "a = 'x'; % savepath",
                  "b = {'a' 'b'}; % savepath", "%{\nsavepath\n%}\nx = 1;", "%{\n%{\nsavepath\n%}\nsetpref\n%}\nx = 1;", "x = 1 + ... savepath\ny = 2;", "% savepath", "x = 1;  % path(",
                  's = "it\'s"; % savepath'):
        assert trips(clean) == [], clean
    # the real script's own title line, with a word after its %s: the old cut would not have seen it
    text = (P / "bench" / "b8_mtpa_map.m").read_text(encoding="utf-8")
    assert trips(text) == []
    hidden = text.replace("version('-release'), numel(idA)", "version('-release'); savepath, numel(idA)", 1)
    assert hidden != text
    old_cut = "\n".join(line.split("%")[0] for line in hidden.splitlines())
    assert b8.forbidden_in_script(old_cut) == [] and trips(hidden) == [savepath]


def test_the_matlab_command_is_a_list_with_the_folder_the_log_the_thread_limit_and_the_script():
    from bench import b8_matlab as b8
    folder = Path("C:/scratch/B8")
    command = b8.matlab_command(Path("C:/MATLAB/bin/matlab.exe"), folder)
    assert isinstance(command, list) and command[0] == Path("C:/MATLAB/bin/matlab.exe")
    flat = [str(part) for part in command]
    assert flat[flat.index("-sd") + 1] == str(folder) and flat[flat.index("-logfile") + 1] == str(folder / "matlab.log")
    assert "-singleCompThread" in flat and "-noFigureWindows" in flat and flat[-2:] == ["-batch", "b8_mtpa_map"]


def test_matlab_prefdir_is_redirected_into_the_case_folder_for_the_run_and_put_back():
    """MathWorks' Support Team: 'To override the default location, create an environment variable named MATLAB_PREFDIR' (MATLAB reads it at startup). run_solver copies the
    environment of this process into the program's, so the variable is set around the call and removed afterwards, whatever happens."""
    from bench import b8_matlab as b8
    assert "MATLAB_PREFDIR" not in os.environ
    with b8.temporary_environment(MATLAB_PREFDIR="C:\\scratch\\prefdir"):
        assert os.environ["MATLAB_PREFDIR"] == "C:\\scratch\\prefdir"
    assert "MATLAB_PREFDIR" not in os.environ
    os.environ["B8_TEST_VARIABLE"] = "before"
    try:
        with pytest.raises(RuntimeError):
            with b8.temporary_environment(B8_TEST_VARIABLE="during"):
                assert os.environ["B8_TEST_VARIABLE"] == "during"
                raise RuntimeError("MATLAB failed")
        assert os.environ["B8_TEST_VARIABLE"] == "before"                                                              # an old value is put back, not deleted
    finally:
        os.environ.pop("B8_TEST_VARIABLE", None)
    assert b8.scratch_prefdir(Path("C:/scratch/B8")) == Path("C:/scratch/B8/prefdir")


# ----------------------------------------------------------------------------------------------- 4. the processes of the job, and the slot
def member(pid, name, ppid, first=0.0, last=1.0):
    from checks.matlab_check import Member
    return Member(pid, name, first, last, 1000 + pid, ppid)


def test_the_matlab_of_the_job_is_one_launcher_with_one_engine_below_it_and_the_service_hosts_are_not_matlab():
    from bench import b8_matlab as b8
    me = 4242
    members = [member(10, "matlab.exe", me), member(11, "MATLAB.exe", 10), member(12, "MathWorksServiceHost.exe", 11), member(13, "MathWorksServiceHost.exe", 11)]
    job = b8.classify_job(members, me)
    assert (job["launchers"], job["engines"], job["hosts"], job["graphics"], job["other_matlab"], job["others"]) == (1, 1, 2, 0, [], [])
    assert b8.exactly_one_matlab(job)
    two = b8.classify_job(members + [member(20, "matlab.exe", me), member(21, "MATLAB.exe", 20)], me)
    assert two["launchers"] == 2 and not b8.exactly_one_matlab(two)
    third = b8.classify_job(members + [member(30, "MATLAB.exe", 12)], me)                                               # a MATLAB started by something else in the job
    assert third["other_matlab"] == ["MATLAB.exe"] and not b8.exactly_one_matlab(third)
    assert not b8.exactly_one_matlab(b8.classify_job([], me)) and not b8.exactly_one_matlab(b8.classify_job([member(10, "matlab.exe", me)], me))   # no MATLAB seen, or no engine
    other = b8.classify_job(members + [member(40, "cmd.exe", me)], me)
    assert other["others"] == ["cmd.exe"]
    drawing = b8.classify_job(members + [member(50, "MATLABWindow.exe", 11), member(51, "MATLABWindow.exe", 11)], me)                       # the figure's window host: a helper
    assert drawing["graphics"] == 2 and drawing["other_matlab"] == [] and b8.exactly_one_matlab(drawing)


def test_the_slot_holder_is_looked_at_while_matlab_runs_and_every_look_inside_the_run_must_be_b8(monkeypatch):
    """B8 asks femm_queue.holder() every 0.2 s while MATLAB runs. The looks between the first and the last sighting of matlab.exe must all name B8's label and this pid; a look
    that finds the slot free, or another run's, inside that window is a failure of 'under the slot'."""
    import femm_queue
    from bench import b8_matlab as b8
    answers = iter([None, {"label": "bench B8 MATLAB", "pid": 99, "since": "x"}, {"label": "bench B8 MATLAB", "pid": 99, "since": "x"}])
    monkeypatch.setattr(femm_queue, "holder", lambda: next(answers))
    sampler = b8.SlotSampler(time.perf_counter(), poll_s=0.01)
    sampler.start()
    deadline = time.monotonic() + 5
    while len(sampler.samples) < 3 and time.monotonic() < deadline:
        time.sleep(0.01)
    sampler.stop()
    assert len(sampler.samples) >= 3 and sampler.samples[0][1] is None and sampler.samples[1][1] == ("bench B8 MATLAB", 99)
    samples = [(0.1, None), (1.0, ("bench B8 MATLAB", 99)), (2.0, ("bench B8 MATLAB", 99)), (3.0, ("bench B8 MATLAB", 99)), (6.0, None)]
    inside = b8.holders_in_window(samples, 0.9, 3.1)
    assert inside == {"looks": 3, "holders": {("bench B8 MATLAB", 99)}}
    assert b8.holders_in_window(samples, 0.0, 6.5)["holders"] == {None, ("bench B8 MATLAB", 99)}                         # a free slot inside the window shows
    assert b8.holders_in_window(samples, 7.0, 8.0) == {"looks": 0, "holders": set()}


def test_a_look_at_the_slot_that_fails_is_counted_and_is_not_a_holder(monkeypatch):
    """HolderSampler.errors counts the looks that raised an OSError; nothing read it (review M7), and now the checks that use the sampler report it. Such a look is kept as the
    sample ("unreadable", None), which is not the holder B8 expects, so it fails 'under the solver slot' as well."""
    import femm_queue
    from bench import b8_matlab as b8
    answers = iter([OSError("locked"), {"label": "bench B8 MATLAB", "pid": 99, "since": "x"}, OSError("locked again")])

    def holder():
        answer = next(answers, None)                                                                                      # after them: a free slot
        if isinstance(answer, Exception):
            raise answer
        return answer
    monkeypatch.setattr(femm_queue, "holder", holder)
    sampler = b8.SlotSampler(time.perf_counter(), poll_s=0.01)
    sampler.start()
    deadline = time.monotonic() + 5
    while len(sampler.samples) < 4 and time.monotonic() < deadline:
        time.sleep(0.01)
    sampler.stop()
    assert sampler.errors == 2 and [seen for _, seen in sampler.samples[:4]] == [("unreadable", None), ("bench B8 MATLAB", 99), ("unreadable", None), None]
    assert b8.holders_in_window(sampler.samples[:3], 0.0, 99.0)["holders"] == {("unreadable", None), ("bench B8 MATLAB", 99)}


# ----------------------------------------------------------------------------------------------- 5. no project file and no MATLAB preference changes
def test_a_folder_is_listed_and_hashed_and_a_change_of_content_size_or_name_is_found(tmp_path):
    from bench import b8_matlab as b8
    (tmp_path / "sub").mkdir()
    (tmp_path / "a.txt").write_text("alpha")
    (tmp_path / "sub" / "b.txt").write_text("beta")
    before = b8.tree_snapshot(tmp_path)
    assert set(before) == {"a.txt", "sub\\b.txt"} and before["a.txt"][0] == 5 and len(before["a.txt"][2]) == 64
    assert b8.tree_diff(before, b8.tree_snapshot(tmp_path)) == {"added": [], "removed": [], "changed": []}
    (tmp_path / "a.txt").write_text("ALPHA")                                                                           # same size, other content
    (tmp_path / "sub" / "c.txt").write_text("gamma")
    (tmp_path / "sub" / "b.txt").unlink()
    assert b8.tree_diff(before, b8.tree_snapshot(tmp_path)) == {"added": ["sub\\c.txt"], "removed": ["sub\\b.txt"], "changed": ["a.txt"]}


def test_key_and_certificate_files_are_listed_but_their_contents_are_not_read(tmp_path, monkeypatch):
    """MATLAB's preferences folder holds a key pair (matlabPkey.p12, thisMatlab.pem). A change is still found, from the size and the nanosecond time, but nothing here opens them."""
    import checks
    from bench import b8_matlab as b8
    key = tmp_path / "matlabPkey.p12"
    key.write_bytes(b"secret-bytes")
    (tmp_path / "plain.json").write_text("{}")
    opened = []
    real = checks.sha256
    monkeypatch.setattr(b8, "sha256", lambda path: (opened.append(Path(path).name), real(path))[1])
    before = b8.tree_snapshot(tmp_path)
    assert opened == ["plain.json"] and before["matlabPkey.p12"][2] is None and before["matlabPkey.p12"][0] == len(b"secret-bytes")
    key.write_bytes(b"other-bytes!")
    os.utime(key, ns=(key.stat().st_atime_ns, key.stat().st_mtime_ns + 5_000_000_000))
    assert b8.tree_diff(before, b8.tree_snapshot(tmp_path))["changed"] == ["matlabPkey.p12"]


def test_the_real_prefdir_is_the_roaming_mathworks_folder_of_matlabs_release(tmp_path, monkeypatch):
    """<APPDATA>\\MathWorks\\MATLAB\\<release>, the form of the path that MathWorks' prefdir page shows as the example of what prefdir returns (C:\\Users\\<name>\\AppData\\Roaming\\
    MathWorks\\MATLAB\\R2025a: an example output, not a statement of a default; review M10). On this PC a `matlab -batch "disp(prefdir)"` with no redirect printed
    C:\\Users\\<name>\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b on 11-06. The release is in VersionInfo.xml of MATLAB's root."""
    from bench import b8_matlab as b8
    for text in (b8.__doc__, b8.real_prefdir.__doc__):
        assert "the Windows default" not in text and "example" in text and "R2025a" in text and "R2025b" in text                  # the page's R2025a is an example; this PC's R2025b is what was seen
    (tmp_path / "VersionInfo.xml").write_text("<?xml version='1.0'?><MathWorks_version_info><version>25.2.0.3177638</version><release>R2025b</release></MathWorks_version_info>")
    assert b8.matlab_release(tmp_path) == "R2025b"
    assert b8.real_prefdir(tmp_path, appdata="C:\\Users\\user\\AppData\\Roaming") == Path("C:\\Users\\user\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b")
    monkeypatch.delenv("APPDATA", raising=False)
    assert b8.real_prefdir(tmp_path) is None                                                                            # no APPDATA, no answer (never a guess)
    monkeypatch.setenv("APPDATA", "C:\\Users\\user\\AppData\\Roaming")
    assert b8.real_prefdir(tmp_path) == Path("C:\\Users\\user\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b")
    assert b8.matlab_release(tmp_path / "nowhere") is None and b8.real_prefdir(tmp_path / "nowhere", appdata="C:\\x") is None


def test_the_scan_of_what_matlab_wrote_elsewhere_lists_the_files_changed_since_the_start_and_skips_what_is_excluded(tmp_path):
    from bench import b8_matlab as b8
    (tmp_path / "logs").mkdir()
    (tmp_path / "prefdir").mkdir()
    old = tmp_path / "logs" / "old.log"
    old.write_text("old")
    os.utime(old, ns=(1, 1))
    (tmp_path / "logs" / "new.log").write_text("new")
    (tmp_path / "prefdir" / "settings.xml").write_text("x")
    found = b8.changed_since([tmp_path, tmp_path / "missing"], time.time_ns() - 5_000_000_000, exclude=[tmp_path / "prefdir"])
    assert [Path(f["path"]).name for f in found] == ["new.log"] and found[0]["bytes"] == 3
    assert {Path(f["path"]).name for f in b8.changed_since([tmp_path], time.time_ns() - 5_000_000_000)} == {"new.log", "settings.xml"}


def test_b8_never_deletes_or_ends_anything_and_asks_the_pool_and_matlab_in_one_place_each_with_a_time_limit():
    from test_bench import queue_calls
    from test_check_tools import scopes_calling
    source = (P / "bench" / "b8_matlab.py").read_text(encoding="utf-8")
    found, missing = queue_calls(source)
    assert found == 1 and missing == []                                                                                # the one femm_pool (run_solver goes through run_program)
    assert scopes_calling(source, None, "run_program") == ["run"]
    assert scopes_calling(source, "femm_queue", "femm_pool") == ["run"]
    import ast
    calls = [n for n in ast.walk(ast.parse(source)) if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "run_program"]
    assert len(calls) == 1 and "timeout_s" in {k.arg for k in calls[0].keywords}                                       # MATLAB runs under a time limit
    for forbidden in ("rmtree", ".unlink(", "os.remove", "shutil.move", "taskkill"):
        assert forbidden not in source, forbidden


# ----------------------------------------------------------------------------------------------- 6. the checks
def good_context():
    """What make_checks needs, for a run in which everything is as the brief says."""
    from bench import b8_matlab as b8
    gamma = b8.gamma_grid()
    torque = 1.12 - 0.5e-4 * (gamma - 23.0) ** 2
    return {
        "angle": {"matlab_deg": 23.01, "python_deg": 23.01, "curve_difference": 3e-15, "angles": len(gamma), "torque_matlab_Nm": 1.12, "torque_python_Nm": 1.12},
        "control": {"closed_form_deg": 22.52, "matlab_deg": 22.52, "python_deg": 22.52},
        "run": {"returncode": 0, "timed_out": False, "seconds": 9.5, "answer_error": None, "figure_bytes": 90_000, "files": {"answer.mat": True, "answer.json": True, "figure.png": True}},
        "files": {"inputs": 2, "inputs_changed": [], "reference_files": 7, "reference_changed": [], "script_changed": False, "unreadable": []},
        "prefdir": {"real": "C:\\Users\\user\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b", "files": 7640, "bytes": 197_400_000, "added": [], "removed": [], "changed": [],
                    "used_by_matlab": "D:\\out\\bench\\B8\\prefdir", "expected": "D:\\out\\bench\\B8\\prefdir", "scratch_files": 8, "unreadable": []},
        "job": {"launchers": 1, "engines": 1, "hosts": 3, "other_matlab": [], "others": [], "alive_after": [], "empty_after_s": 0.02,
                "members": [("matlab.exe", 0.0, 9.4), ("MATLAB.exe", 0.1, 9.4), ("MathWorksServiceHost.exe", 0.8, 9.2)]},
        "slot": {"expected": ("bench B8 MATLAB", 4242), "looks": 40, "looks_total": 52, "holders": [("bench B8 MATLAB", 4242)], "unreadable_looks": 0},
        "threads": {"matlab": 1, "limit": 8},
        "priority": {"classes_seen": {"matlab.exe": ["below_normal"], "mathworksservicehost.exe": ["below_normal"], "matlabwindow.exe": ["below_normal"]}, "looks": 300, "raised": {}},
        "map": {"points": (6, 6), "solves": 72, "max_current_error": 2e-15, "circle_inside": True, "psi_m_mWb": 16.2, "ld_mH": 0.16, "lq_mH": 0.42, "psi_q_at_iq0_mWb": 0.01,
                "psi_q_max_mWb": 20.0, "torque_at_gamma0_Nm": 1.0, "torque_at_gamma90_Nm": 0.0},
        "b1": {"gamma_deg": 22.936, "psi_pm_mWb": 16.629, "curve": {"0.0": 0.998, "22.5": 1.129, "45.0": 0.81, "90.0": 0.0}, "peak_Nm": 1.1291, "path": "out/bench/B1/results.json"},
        "b1_reason": None,
        "map_vs_b1": {"max_difference_percent_of_peak": 0.9, "angles": [0.0, 22.5, 45.0, 90.0]},
        "spline": {"angle_deg": 22.4},
        "others": {"before": [(6912, "MathWorksServiceHost.exe"), (17872, "MathWorksCrashReporter.exe"), (19432, "MathWorksServiceHost-Monitor.exe")],
                   "after": [(6912, "MathWorksServiceHost.exe"), (17872, "MathWorksCrashReporter.exe"), (19432, "MathWorksServiceHost-Monitor.exe")]},
        "outside": {"count": 3, "bytes": 3_000_000, "by_root": {"C:\\Users\\user\\AppData\\Local\\MathWorks": 2, "C:\\Users\\user\\AppData\\Roaming\\MathWorks": 1}, "in_repo": [],
                    "files": [{"path": "C:\\Users\\user\\AppData\\Local\\MathWorks\\ServiceHost\\logs\\service-host.log", "bytes": 1_200_000},
                              {"path": "C:\\Users\\user\\AppData\\Local\\MathWorks\\ServiceHost\\config\\credentials.json", "bytes": 800},
                              {"path": "C:\\Users\\user\\AppData\\Roaming\\MathWorks\\licensing\\usage_data.txt", "bytes": 1_799_200}]},
    }


def verdicts(context):
    from bench import b8_matlab as b8
    return {c["name"]: c for c in b8.make_checks(context)}


def test_a_run_in_which_everything_is_as_the_brief_says_passes_every_gating_check():
    checks = verdicts(good_context())
    assert all(c["verdict"] == "PASS" for c in checks.values() if c.get("gating", True)), [n for n, c in checks.items() if c["verdict"] != "PASS" and c.get("gating", True)]
    names = " | ".join(checks)
    for needle in ("MTPA angle", "no project file and no MATLAB preference changes", "prefdir", "exactly one MATLAB", "slot", "empty", "threads", "priority"):
        assert needle in names, needle
    assert "outside the scratch folder" not in names                                                                           # the name said more than the check holds (MathWorks' own files change)


def test_matlab_half_a_degree_away_fails_the_brief_s_check_and_only_that_one_in_its_group():
    context = good_context()
    context["angle"] = {**context["angle"], "matlab_deg": 23.01 + 0.51}
    failed = [n for n, c in verdicts(context).items() if c["verdict"] == "FAIL" and c.get("gating", True)]
    assert failed == [n for n in failed if "MTPA angle" in n] and len(failed) == 1
    context["angle"] = {**context["angle"], "matlab_deg": 23.01 + 0.5}
    assert not [n for n, c in verdicts(context).items() if c["verdict"] == "FAIL" and c.get("gating", True)]               # exactly on the half degree passes


def test_each_way_the_run_can_touch_a_file_or_break_a_rule_turns_its_check_red():
    def failed(**changes):
        context = good_context()
        for key, value in changes.items():
            context[key] = {**context[key], **value}
        return {n for n, c in verdicts(context).items() if c["verdict"] == "FAIL" and c.get("gating", True)}
    scope = "no project file and no MATLAB preference changes"
    assert any(scope in n for n in failed(files={"reference_changed": ["Rotor.DXF"]}))
    assert any(scope in n for n in failed(files={"inputs_changed": ["flux_map.mat"]}))
    assert any(scope in n for n in failed(files={"script_changed": True}))
    assert any(scope in n for n in failed(prefdir={"changed": ["matlab.mlsettings"]}))                                         # MATLAB's real preferences folder is part of what must not change
    assert any(scope in n for n in failed(prefdir={"added": ["new.json"]}))
    assert any(scope in n for n in failed(prefdir={"removed": ["old.json"]}))
    assert any(scope in n for n in failed(prefdir={"found": False}))                                                           # not found: it cannot be shown unchanged
    assert any(scope in n for n in failed(files={"unreadable": ["reference\\Rotor.DXF"]}))                                     # a file that could not be read cannot be shown unchanged
    assert any(scope in n for n in failed(prefdir={"unreadable": ["matlab.mlsettings"]}))
    assert any("prefdir:" in n for n in failed(prefdir={"used_by_matlab": "C:\\Users\\user\\AppData\\Roaming\\MathWorks\\MATLAB\\R2025b"}))                # the redirect was not honoured
    assert any("exactly one MATLAB" in n for n in failed(job={"launchers": 2}))
    assert any("slot" in n for n in failed(slot={"holders": [("bench B8 MATLAB", 4242), None]}))
    assert any("slot" in n for n in failed(slot={"holders": [("someone else", 7)]}))
    assert any("slot" in n for n in failed(slot={"looks": 0, "holders": []}))                                                  # no look inside the run: it cannot be said
    assert any("slot" in n for n in failed(slot={"holders": [("bench B8 MATLAB", 4242), ("unreadable", None)]}))              # a look that could not be made is not the holder
    assert any("empty" in n for n in failed(job={"alive_after": ["MathWorksServiceHost.exe"]}))
    assert any("threads" in n for n in failed(threads={"matlab": 9}))
    assert any("threads" in n for n in failed(threads={"matlab": -1}))                                                         # MATLAB could not say
    assert any("priority" in n for n in failed(priority={"classes_seen": {"matlab.exe": ["normal"]}}))
    assert any("priority" in n for n in failed(priority={"classes_seen": {"matlab.exe": ["below_normal", "above_normal"]}}))
    assert any("priority" in n for n in failed(priority={"classes_seen": {}}))                                                  # nothing seen: it cannot pass
    assert any("ran" in n for n in failed(run={"returncode": 1}))
    assert any("ran" in n for n in failed(run={"timed_out": True, "returncode": None}))
    assert any("ran" in n for n in failed(run={"files": {"answer.mat": False, "answer.json": True, "figure.png": True}}))
    assert any("control" in n for n in failed(control={"matlab_deg": 20.0}))
    assert any("curve" in n for n in failed(angle={"curve_difference": 1e-6}))


def test_every_process_of_the_job_must_be_below_normal_and_the_graphics_host_gates_like_the_rest():
    """The solver slot's job holds every process in it at the class run_solver was asked for (femm_queue._Job, JOB_OBJECT_LIMIT_PRIORITY_CLASS), so the rule is the plain one, R10:
    every process of MATLAB's job below normal, MATLABWindow.exe (the figure's graphics host) among them. The check gates, says the rule and not the history, and its evidence is
    each image name's classes."""
    def checks(**changes):
        context = good_context()
        context["priority"] = {**context["priority"], **changes}
        return verdicts(context)

    def rule(all_checks):
        return [c for n, c in all_checks.items() if "every process of MATLAB's job" in n]
    good = rule(checks())
    assert len(good) == 1 and good[0]["verdict"] == "PASS" and good[0]["gating"] is True
    for image in ("matlab.exe", "mathworksservicehost.exe", "matlabwindow.exe"):
        assert f"{image}: below_normal" in good[0]["measured"], good[0]["measured"]                                             # each image name, with its classes
    for text in (good[0]["name"], good[0]["expected"]):
        assert not any(history in text for history in ("was seen at", "not isolated", "does not hold", "sets its own", "controller rules")), text
    assert "JOB_OBJECT_LIMIT_PRIORITY_CLASS" in good[0]["expected"] and "below normal" in good[0]["name"]
    # a graphics host seen above normal: the rule check FAILS, gating, with the evidence; the verdict of the whole case is FAIL
    seen = {"matlab.exe": ["below_normal"], "mathworksservicehost.exe": ["below_normal"], "matlabwindow.exe": ["above_normal", "below_normal", "normal"]}
    raised = {"MATLABWindow.exe (pid 77)": ["below_normal at 4.3 s", "normal at 4.4 s"], "MATLABWindow.exe (pid 78)": ["below_normal at 4.3 s", "above_normal at 4.5 s"]}
    after = checks(classes_seen=seen, raised=raised)
    item = rule(after)[0]
    assert item["verdict"] == "FAIL" and item["gating"] is True
    assert "matlabwindow.exe: above_normal, below_normal, normal" in item["measured"] and "pid 77" in item["measured"] and "normal at 4.4 s" in item["measured"]
    assert [n for n, c in after.items() if c["verdict"] == "FAIL" and c.get("gating", True)] == [item["name"]]                  # and the core check, which has its own processes, is not what failed
    assert rule(checks(classes_seen={}))[0]["verdict"] == "FAIL"                                                                # nothing seen: it cannot pass
    quiet = rule(checks(classes_seen={"matlab.exe": ["below_normal"], "mathworksservicehost.exe": ["below_normal"]}))[0]          # no figure window started in this run
    assert quiet["verdict"] == "PASS" and "matlabwindow.exe: not seen" in quiet["measured"]
    # the core processes are still judged by their own check, with the launcher, the engine and the service hosts in its name
    core = [c for n, c in checks().items() if "launcher" in n and "priority" in n]
    assert len(core) == 1 and core[0]["gating"] is True


def test_the_case_has_22_checks_of_which_18_gate_and_4_are_reported():
    """After the fix round: the graphics host's check gates (it was a reported FAIL), so 18 gating and 4 reported of the 22."""
    checks = list(verdicts(good_context()).values())
    assert len(checks) == 22
    gating, reported = [c for c in checks if c.get("gating", True)], [c for c in checks if not c.get("gating", True)]
    assert (len(gating), len(reported)) == (18, 4), [c["name"] for c in reported]
    assert all(c["verdict"] in ("PASS", "FAIL", "SKIPPED") for c in checks)
    names = " | ".join(c["name"] for c in reported)
    for needle in ("empties at once", "against B1's", "spline", "MathWorks' own files"):
        assert needle in names, needle


def test_the_reported_list_of_mathworks_files_names_every_path_with_its_size():
    """MATLAB writes some files that MATLAB_PREFDIR does not reach (the service host's logs, its credentials service's config, the license file, the usage log). They stay a
    reported list, and the list names every path, not a count by folder (a design decision)."""
    context = good_context()
    item = next(c for n, c in verdicts(context).items() if "MathWorks' own files" in n)
    assert item["gating"] is False and item["verdict"] == "PASS" and "reported, not gating" in item["expected"]
    for f in context["outside"]["files"]:
        assert f["path"] in item["measured"] and f"{f['bytes']:,}" in item["measured"], f
    context["outside"] = {**context["outside"], "in_repo": ["D:\\Storage\\x.txt"]}
    after = next(c for n, c in verdicts(context).items() if "MathWorks' own files" in n)
    assert after["verdict"] == "FAIL" and "D:\\Storage\\x.txt" in after["measured"]                                               # a project file that changed shows, as before


def test_the_slot_check_says_how_many_looks_could_not_be_made():
    """HolderSampler counts the looks that failed (an OSError); the evidence of 'under the solver slot' reports that count, as the watcher's failed looks are reported."""
    context = good_context()
    context["slot"] = {**context["slot"], "unreadable_looks": 2}
    item = next(c for n, c in verdicts(context).items() if "solver slot" in n)
    assert "2 looks could not be made" in item["measured"]
    assert "could not be made" not in next(c for n, c in verdicts(good_context()).items() if "solver slot" in n)["measured"]


def test_a_file_that_cannot_be_read_or_vanishes_is_said_in_the_snapshot_and_not_raised(tmp_path, monkeypatch):
    """The real preferences folder belongs to the user's own MATLAB: a file there can be locked, or can vanish between the listing and the look (review M4, Task 12's lesson 3).
    The snapshot keeps what it could see, says which file it could not read, and does not raise; one that vanished is simply not in it."""
    from bench import b8_matlab as b8
    for name, text in (("ok.txt", "fine"), ("locked.txt", "held by another program"), ("gone.txt", "deleted while it is looked at"), ("nostat.txt", "cannot even be looked at")):
        (tmp_path / name).write_text(text)
    real_sha = b8.sha256

    def sha(path):
        if Path(path).name == "locked.txt":
            raise PermissionError(13, "The process cannot access the file because it is being used by another process")
        if Path(path).name == "gone.txt":
            raise FileNotFoundError(2, "no such file")
        return real_sha(path)
    original_stat = Path.stat

    def stat(self, *args, **kwargs):
        if self.name == "nostat.txt":
            raise PermissionError(13, "access denied")
        return original_stat(self, *args, **kwargs)
    monkeypatch.setattr(b8, "sha256", sha)
    monkeypatch.setattr(Path, "stat", stat)
    snapshot = b8.tree_snapshot(tmp_path)
    assert set(snapshot) == {"ok.txt", "locked.txt", "nostat.txt"}                                                              # the vanished one is not there
    assert len(snapshot["ok.txt"][2]) == 64
    assert snapshot["locked.txt"][0] == len("held by another program") and snapshot["locked.txt"][2].startswith("unreadable: PermissionError")
    assert snapshot["nostat.txt"][2].startswith("unreadable: PermissionError")
    assert b8.unreadable_files(snapshot) == ["locked.txt", "nostat.txt"] and b8.unreadable_files({"a": (1, 2, None), "b": (1, 2, "ab" * 32)}) == []
    assert b8.tree_diff(snapshot, snapshot) == {"added": [], "removed": [], "changed": []}                                      # the same unreadable file twice is not a change; the check still names it


def test_the_class_sampler_keeps_each_process_s_classes_in_the_order_it_saw_them_and_names_the_ones_that_left_below_normal(monkeypatch):
    import itertools

    import bench
    from bench import b8_matlab as b8
    script = {101: [0x4000, 0x4000, 0x20], 102: [0x4000], 103: [None, 0x4000, 0x8000]}                                     # what each look finds; the last value stays
    calls = {pid: itertools.chain(values, itertools.repeat(values[-1])) for pid, values in script.items()}

    class Tracker:
        def members(self):
            return [member(101, "MATLABWindow.exe", 1), member(102, "MATLAB.exe", 1), member(103, "MATLABWindow.exe", 1)]
    monkeypatch.setattr(bench, "priority_class", lambda pid: next(calls[pid]))
    monkeypatch.setattr(b8, "process_created", lambda pid: 1000 + pid)                                                          # the creation time `member` carries: the same process
    sampler = b8.ClassSampler(Tracker(), time.perf_counter(), poll_s=0.005)
    sampler.start()
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and not ([c for _, c in sampler.sequences.get((103, "MATLABWindow.exe"), [])][-1:] == ["above_normal"]
                                               and [c for _, c in sampler.sequences.get((101, "MATLABWindow.exe"), [])][-1:] == ["normal"]):
        time.sleep(0.01)
    sampler.stop()
    seq = {key: [c for _, c in value] for key, value in sampler.sequences.items()}
    assert seq[(102, "MATLAB.exe")] == ["below_normal"]                                                                        # a class is kept once, when it changes
    assert seq[(101, "MATLABWindow.exe")] == ["below_normal", "normal"]
    assert seq[(103, "MATLABWindow.exe")] == ["below_normal", "above_normal"]                                                   # a look that could not read it (None) is no class
    raised = sampler.left_below_normal()
    assert set(raised) == {"MATLABWindow.exe (pid 101)", "MATLABWindow.exe (pid 103)"} and "MATLAB.exe (pid 102)" not in raised
    assert raised["MATLABWindow.exe (pid 101)"][0].startswith("below_normal at ") and raised["MATLABWindow.exe (pid 101)"][1].startswith("normal at ")


def test_the_class_sampler_skips_a_pid_that_now_belongs_to_another_process_and_one_whose_creation_time_it_cannot_check(monkeypatch):
    """The tracker keeps every member it ever saw, the ones that ended too, with the creation time it read. Windows can give the pid of an ended member to a stranger: that
    process has another creation time, and its class must not be added to the member's sequence (review M3). A member that carries no creation time cannot be checked and is
    not sampled either."""
    import bench
    from bench import b8_matlab as b8
    from checks.matlab_check import Member

    class Tracker:
        def members(self):
            return [member(101, "MATLABWindow.exe", 1), member(102, "MATLABWindow.exe", 1), Member(104, "MATLABWindow.exe", 0.0, 1.0, None, 1)]
    live = {101: 1101, 102: 9999, 104: 1104}                                                                                   # 102 is another process now; 104 never had a time
    monkeypatch.setattr(bench, "priority_class", lambda pid: 0x20)                                                              # whoever answers is at normal
    monkeypatch.setattr(b8, "process_created", lambda pid: live[pid])
    sampler = b8.ClassSampler(Tracker(), time.perf_counter(), poll_s=0.005)
    sampler.start()
    deadline = time.monotonic() + 10
    while time.monotonic() < deadline and (101, "MATLABWindow.exe") not in sampler.sequences:
        time.sleep(0.01)
    time.sleep(0.15)                                                                                                           # a few more passes over the members
    sampler.stop()
    assert set(sampler.sequences) == {(101, "MATLABWindow.exe")}, sampler.sequences
    assert [c for _, c in sampler.sequences[(101, "MATLABWindow.exe")]] == ["normal"]


def test_what_the_brief_does_not_gate_is_reported_with_its_own_verdict_and_never_hidden():
    checks = verdicts(good_context())
    reported = {n: c for n, c in checks.items() if not c.get("gating", True)}
    assert any("spline" in n for n in reported) and any("MathWorks" in n or "elsewhere" in n for n in reported)
    for item in reported.values():
        assert item["verdict"] in ("PASS", "FAIL", "SKIPPED") and ("reported, not gating" in item["expected"] or item["verdict"] == "SKIPPED")
    context = good_context()
    context["spline"] = {"angle_deg": 31.0}
    context["outside"] = {**context["outside"], "in_repo": ["D:\\Storage\\x.txt"]}
    after = verdicts(context)
    assert [n for n, c in after.items() if not c.get("gating", True) and c["verdict"] == "FAIL"]                              # a failing reported check shows as FAIL, it does not gate
    assert all(c["verdict"] == "PASS" for c in after.values() if c.get("gating", True))


def test_without_b1_s_results_the_comparisons_with_it_are_skipped_with_the_reason_the_shared_reader_gave():
    """'No B1 results in this root' and 'B1's results lack key K' are told apart (review M5): the reason is femm_points.b1_reference()'s, word for word."""
    for reason in ("no B1 results in this run's root: D:\\x\\B1\\results.json is not there (run B1 first)", "B1's results (D:\\x\\B1\\results.json) lack summary.none.psi_pm_mWb"):
        context = good_context()
        context["b1"], context["b1_reason"] = None, reason
        checks = verdicts(context)
        skipped = [c for c in checks.values() if c["verdict"] == "SKIPPED"]
        assert len(skipped) == 3 and all(not c["gating"] and c["measured"] == reason for c in skipped)
        assert all(c["verdict"] == "PASS" for c in checks.values() if c.get("gating", True))


def test_the_map_must_look_like_this_motor_before_it_is_trusted():
    """Hand checks on the map itself (physics, from the earlier results): psi_d at no current is the magnets' flux (B1: 16.63 mWb on the fine mesh), the saliency is there
    (Lq > Ld, so the reluctance torque exists), psi_q is zero with no q current (the machine is symmetric about its d axis) and so is the torque at gamma = 90 degrees."""
    def failed(**changes):
        context = good_context()
        context["map"] = {**context["map"], **changes}
        return {n for n, c in verdicts(context).items() if c["verdict"] == "FAIL" and c.get("gating", True)}
    assert not failed()
    assert any("magnet flux" in n for n in failed(psi_m_mWb=20.0))
    assert any("saliency" in n for n in failed(ld_mH=0.5, lq_mH=0.4))
    assert any("q axis" in n for n in failed(psi_q_at_iq0_mWb=3.0))
    assert any("driven" in n for n in failed(max_current_error=1e-3))
    assert any("grid" in n for n in failed(points=(4, 6)))
    assert any("grid" in n for n in failed(circle_inside=False))
    context = good_context()
    context["map_vs_b1"] = {**context["map_vs_b1"], "max_difference_percent_of_peak": 6.0}
    assert any("B1" in n and "curve" in n for n, c in verdicts(context).items() if c["verdict"] == "FAIL" and c.get("gating", True))


def test_matlabs_release_is_named_the_way_mathworks_names_it_with_the_r(tmp_path):
    """MATLAB's version('-release') answers '2025b'; MathWorks' name is 'R2025b' (review C2, the figure's title; the log line said '2025b' too). Not doubled, and an unknown release is said."""
    from bench import b8_matlab as b8
    assert b8.release_name("2025b") == "R2025b" and b8.release_name("R2025b") == "R2025b"
    assert b8.release_name(None) == "(release not reported)" and b8.release_name("") == "(release not reported)"


def test_the_scripts_comment_states_what_was_seen_and_not_a_fact_about_matlabs_own_theme():
    """T13-P2. The comment said '(R2025b's own is dark: grey text and a dark axes area on white)' as a fact about MATLAB, which did not reproduce; what was seen is that a figure takes the
    theme MATLAB runs with (under a dark graphics theme, set for one session in scratch preferences, the old script gave grey text and a dark panel), and that the light theme keeps this
    figure legible on white. And the word 'settings' is not used in a comment of a script that must never touch a setting: 'no MATLAB setting is touched'."""
    text = (P / "bench" / "b8_mtpa_map.m").read_text(encoding="utf-8")
    comment = "\n".join(line for line in text.splitlines() if line.lstrip().startswith("%"))
    assert "R2025b's own is dark" not in text and "whatever theme MATLAB itself runs with" not in text
    assert "a figure takes the theme MATLAB runs with" in comment and "grey text and a dark panel" in comment and "legible on white" in comment
    assert "no MATLAB setting is touched" in comment and not re.search(r"\bsettings\b", comment), comment
