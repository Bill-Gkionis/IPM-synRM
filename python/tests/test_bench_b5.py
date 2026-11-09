"""Tests for B5 (bench/b5_prepomax.py and bench/b5_prepomax_script.py): PrePoMax from the command line. Everything that needs no PrePoMax: the options and what the program says about
itself, settings.bin and the processes, B4's results, the case file, what the importer made, the picture, and the script that runs inside PrePoMax.

Run from P:   uv run pytest tests/test_bench_b5.py -q      (a second or two; PrePoMax is not started, no solver slot is taken)

In groups:
  1. The brief's 1 %, the case's constants, and the options with their sources.
  2. What the program prints (--version, --help), and what the help screen may lack.
  3. settings.bin (hash and solver fields) and the processes: what a run may leave, and what an unreadable process table says (it is unknown, and the check of it fails).
  4. B4's results as the reference (or a reason to skip), the case file, the importer's findings, the two readers of a JSON file, and the checks of a run that went right and of one
     whose script failed a step (a failed check that names the step, not a KeyError).
  5. The picture: PrePoMax's legend and model cut out of its viewport, without the scale bar, and its caption (no line that wraps).
  6. The script that runs inside PrePoMax: valid Python, the standard library only, the steps it takes, and what it does when imported without `pmx`, when a step works and when one fails.
"""
import ast
import json
import math
from pathlib import Path

import numpy as np
import pytest

import paths
from bench_test_tools import P

HELP_SCREEN = """PrePoMax 3.0.0.0
Copyright (C) 2026 Matej Borovnek

  -f, --file                     File name to be opened/imported.

  -g, --showGui                  (Default: Yes) Show graphical user interface:
                                 Yes | No. No runs without a window and can only
                                 be used for regeneration, for running a script
                                 or with the automation server (--rpc On or
                                 ReadOnly).

  --script                       The Python script file name to be executed
                                 after the application starts.

  --scriptArgs                   The arguments passed to the script as sys.argv.

  --save                         The .pmx file name to be saved after the script
                                 finishes.

  --onPrompt                     The reaction to a dialog that asks for a decision.

  --report                       The .json file name to be used for writing the
                                 machine readable run report.

  --journal                      The .py file name to be used for writing the
                                 session journal.

  --timeout                      (Default: 0) The maximum duration of the run in
                                 seconds.

  --rpc                          (Default: Off) Start the automation server.

  --rpcPort                      (Default: 0) The port the MCP endpoint listens on.

  --rpcTimeout                   (Default: 0) How many seconds the automation
                                 server waits for a request.

  --settings                     The settings file name to be used instead of
                                 the settings file located next to the executable.

  --noSettingsSave               (Default: false) Do not write the settings file.

  -o, --overwrite                (Default: No) Overwrite the .pmx file after
                                 regeneration: Yes | No.

  -p, --parameters               Overwrite the model parameters.

  --parametersExport             The .pmp file name to be used for exporting
                                 parameters.

  -r, --regenerate               The .pmx file name to be used for regeneration.

  --r0                           Regeneration of all history commands. See
                                 option -r.

  --r1                           Regeneration of pre-processing history commands
                                 only.

  --r2                           Regeneration of analysis history commands only.

  --r3                           Regeneration of post-processing history
                                 commands only.

  -u, --unitSystem               Unit system type to be used when importing.

  -w, --workDirectory            The directory path to be used as work
                                 directory.

  -x, --exitAfterRegeneration    (Default: Yes) Exit PrePoMax after
                                 regeneration: Yes | No.

  --help                         Display this help screen.

  --version                      Display version information.
"""


# ----------------------------------------------------------------------------------------------- 1. the brief's 1 %, the constants, the options
def test_the_briefs_one_percent_and_the_cases_constants_are_pinned():
    """The brief: PrePoMax reproduces B4's peak within 1 %. The case is B4's pressure model at 10,000 rpm on B4's main mesh. Change a constant and this goes red."""
    from bench import b4_ccx as b4
    from bench import b5_prepomax as b5
    assert b5.TOLERANCE_PERCENT == 1.0
    assert (b5.VARIANT, b5.RPM, b5.H_MIN) == ("pressure", 10_000, 0.1) == (b4.VERDICT_VARIANT, b4.SPEEDS_RPM[0], b4.H_MIN_MAIN)
    assert b5.peak_check(33.0 * 1.0099, 33.0, "x", "s")["verdict"] == "PASS" and b5.peak_check(33.0 * 1.0101, 33.0, "x", "s")["verdict"] == "FAIL"      # PrePoMax against the reference, in percent of the reference
    assert b5.peak_check(33.0 * 0.9901, 33.0, "x", "s")["verdict"] == "PASS" and b5.peak_check(33.0 * 0.9899, 33.0, "x", "s")["verdict"] == "FAIL"
    item = b5.peak_check(33.0, 33.0, "x", "the source")
    assert item["tolerance_percent"] == 1.0 and item["gating"] is True and item["unit"] == "MPa" and item["source"] == "the source"
    assert b5.peak_check(33.0, 40.0, "x", "s", gating=False, reason="because")["gating"] is False


def test_every_option_comes_with_what_it_does_and_the_page_it_comes_from():
    """The brief for the report: the options, with the documentation each comes from. The manual's table (section 3, pages 7-8) and the program's own help screen, which lists the script
    options that the web manual does not yet describe."""
    from bench import b5_prepomax as b5
    assert len(b5.OPTIONS) >= 19
    for option in b5.OPTIONS:
        assert set(option) == {"option", "does", "source"} and option["does"] and option["source"]
    text = " ".join(o["option"] for o in b5.OPTIONS)
    for name in b5.MANUAL_OPTIONS + b5.HELP_ONLY_OPTIONS:
        assert name in text or name in ("help", "version"), name
    manual = {o["option"] for o in b5.OPTIONS if "manual" in o["source"]}
    assert {"-f, --file", "-g, --showGui", "-r, --regenerate", "-w, --workDirectory", "-x, --exitAfterRegeneration"} <= manual
    assert all("PrePoMax.com --help" in o["source"] for o in b5.OPTIONS if o["option"].startswith("--script") or o["option"] in ("--report", "--journal", "--timeout", "--save"))
    assert set(b5.MANUAL_OPTIONS) == {"file", "showGui", "overwrite", "parameters", "parametersExport", "regenerate", "unitSystem", "workDirectory", "exitAfterRegeneration", "help", "version"}
    assert set(b5.HELP_ONLY_OPTIONS) >= {"script", "scriptArgs", "save", "onPrompt", "report", "journal", "timeout", "rpc", "settings", "noSettingsSave", "r0", "r1", "r2", "r3"}
    assert sorted(b5.EXIT_CODES) == list(range(8)) and b5.EXIT_CODES[0] == "completed" and b5.EXIT_CODES[3] == "script error" and b5.EXIT_CODES[4] == "timed out"


# ----------------------------------------------------------------------------------------------- 2. what the program prints
def test_the_help_screen_is_read_for_every_option_it_lists():
    from bench import b5_prepomax as b5
    listed = b5.parse_help(HELP_SCREEN)
    assert listed[:3] == ["file", "showGui", "script"] and "noSettingsSave" in listed and "r3" in listed and listed[-2:] == ["help", "version"]
    assert len(listed) == 27 and len(set(listed)) == 27
    assert b5.options_missing(listed) == {"manual": [], "help_only": []}


def test_a_help_screen_that_lacks_an_option_says_which():
    from bench import b5_prepomax as b5
    lacking = HELP_SCREEN.replace("  --script                       The Python", "  (removed)                      The Python").replace("  -w, --workDirectory", "  (removed)").replace("  --journal  ", "  (removed)  ")
    listed = b5.parse_help(lacking)
    missing = b5.options_missing(listed)
    assert missing["manual"] == ["workDirectory"] and "script" in missing["help_only"] and "journal" in missing["help_only"]
    assert b5.options_missing([]) == {"manual": list(b5.MANUAL_OPTIONS), "help_only": list(b5.HELP_ONLY_OPTIONS)}


def test_the_version_is_read_from_the_first_line_of_what_the_program_prints():
    from bench import b5_prepomax as b5
    assert b5.parse_version("PrePoMax 3.0.0.0\n") == "3.0.0.0" and b5.parse_version("\nPrePoMax 3.0.0\nCopyright") == "3.0.0"
    assert b5.parse_version("Version 2.23") is None and b5.parse_version("") is None


# ----------------------------------------------------------------------------------------------- 3. settings.bin and the processes
def test_settings_state_reads_the_solver_fields_and_the_hash_without_starting_anything(monkeypatch, tmp_path):
    """settings.bin is read as bytes, hashed, and its fields read with checks/nrbf.py; a path inside PrePoMax's folder is stored with a leading # and resolves to that folder."""
    import checks
    from bench import b5_prepomax as b5
    program = tmp_path / "PrePoMax v3.0.0" / "PrePoMax.exe"
    program.parent.mkdir()
    (program.parent / "settings.bin").write_bytes(b"settings-bytes")
    seen = {}

    enum = b5.nrbf.Obj(7, "CaeModel.SolverTypeEnum", ["value__"])           # an enum as nrbf.py reads it: an object with one member (found 2026-11-06: it is not a number, and not JSON)
    enum.values = [3]

    def fake_members(data, class_name):
        seen["class"] = class_name
        seen["data"] = data
        return {"_calculixExe": "#\\Solver\\ccx_dynamic.exe", "_calculixDefaultSolverType": enum}

    monkeypatch.setattr(b5.nrbf, "read_class_members", fake_members)
    state = b5.settings_state(program)
    assert seen == {"class": "PrePoMax.SolverSettings", "data": b"settings-bytes"}
    assert state["sha256"] == checks.sha256(program.parent / "settings.bin") and state["bytes"] == len(b"settings-bytes")
    assert state["calculix_exe"] == "#\\Solver\\ccx_dynamic.exe" and state["default_solver_type"] == 3 and type(state["default_solver_type"]) is int
    assert Path(state["solver"]) == program.parent / "Solver" / "ccx_dynamic.exe"
    json.dumps(state)                                                          # it goes into results.json
    assert b5.enum_value(3) == 3 and b5.enum_value(None) is None


def test_a_run_may_leave_no_process_of_its_own_and_must_leave_bill_s_alone():
    """Before the run a PrePoMax that somebody has open is recorded; afterwards the processes that were not there before are this run's and must be gone, and that one must still be there."""
    from bench import b5_prepomax as b5
    assert b5.processes_left({}, {}) == {"readable": True, "reason": "", "started_and_left": {}, "bills": {}}
    before = {100: "PrePoMax.exe"}
    assert b5.processes_left(before, {100: "PrePoMax.exe"}) == {"readable": True, "reason": "", "started_and_left": {}, "bills": {100: "PrePoMax.exe"}}
    assert b5.processes_left(before, {100: "PrePoMax.exe", 200: "PrePoMax.com", 201: "PrePoMax.exe"}) == {"readable": True, "reason": "", "started_and_left": {200: "PrePoMax.com", 201: "PrePoMax.exe"},
                                                                                                           "bills": {100: "PrePoMax.exe"}}
    assert b5.processes_left(before, {}) == {"readable": True, "reason": "", "started_and_left": {}, "bills": {}}                      # (if that window was closed meanwhile, that is not this run's)


def test_a_process_table_that_could_not_be_read_is_unknown_and_says_when():
    """Fix round 1, item 6: an unreadable table (None) is not an empty one. Nothing can be said about what a run left, and the reason names which of the two snapshots failed."""
    from bench import b5_prepomax as b5
    unknown_before = b5.processes_left(None, {})
    assert unknown_before["readable"] is False and "before the run" in unknown_before["reason"] and "after" not in unknown_before["reason"] and unknown_before["started_and_left"] == {}
    unknown_after = b5.processes_left({100: "PrePoMax.exe"}, None)
    assert unknown_after["readable"] is False and "after the run" in unknown_after["reason"] and "before" not in unknown_after["reason"] and unknown_after["bills"] == {}
    both = b5.processes_left(None, None)
    assert both["readable"] is False and "before and after the run" in both["reason"]
    assert all("process table" in r["reason"] for r in (unknown_before, unknown_after, both))


def test_a_snapshot_of_the_processes_is_said_in_words_and_an_unreadable_one_says_so():
    from bench import b5_prepomax as b5
    assert b5.snapshot_words({}) == "none" and b5.snapshot_words({10: "PrePoMax.exe"}) == "{10: 'PrePoMax.exe'}"
    assert b5.snapshot_words(None).startswith("unknown") and "could not be read" in b5.snapshot_words(None)
    assert b5.process_words(b5.processes_left({}, {}), "started_and_left") == "none"
    assert "unknown" in b5.process_words(b5.processes_left(None, {}), "bills") and "before the run" in b5.process_words(b5.processes_left(None, {}), "bills")


def test_the_process_snapshot_keeps_only_the_prepomax_programs(monkeypatch):
    from bench import b5_prepomax as b5
    monkeypatch.setattr(b5, "process_table", lambda: {1: (0, "System"), 10: (1, "PrePoMax.exe"), 11: (10, "ccx_static.exe"), 12: (1, "prepomax.COM"), 13: (1, "python.exe")})
    assert b5.prepomax_processes() == {10: "PrePoMax.exe", 12: "prepomax.COM"}

    def unreadable():
        raise OSError("no snapshot")
    monkeypatch.setattr(b5, "process_table", unreadable)
    assert b5.prepomax_processes() is None                                      # not {}: an unreadable table is not an empty one (fix round 1, item 6)


# ----------------------------------------------------------------------------------------------- 4. B4's results, the case file, the importer, the files
def b4_results(path, nodes=70534, peak=33.454, key="pressure_10000_h0.1"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"mesh": {"nodes": nodes}, "runs": {key: {"peak": {"value_mpa": peak}}}}), encoding="utf-8")
    return path


def test_b4_s_peak_is_the_reference_when_it_is_for_the_same_mesh_and_a_reason_to_skip_when_it_is_not(tmp_path):
    """The ruling: B5 compares with out\\bench\\B4\\results.json; if the file is not there the comparison is SKIPPED with the reason. It is also skipped when the file is for another
    mesh (not the same problem) or has no such run."""
    from bench import b5_prepomax as b5
    good = b4_results(tmp_path / "B4" / "results.json")
    got = b5.b4_reference(70534, good)
    assert got == {"ok": True, "peak": 33.454, "file": str(good), "nodes": 70534}
    assert b5.b4_reference(70000, good)["ok"] is False and "70,534" in b5.b4_reference(70000, good)["reason"] and "not the same problem" in b5.b4_reference(70000, good)["reason"]
    missing = b5.b4_reference(70534, tmp_path / "nope" / "results.json")
    assert missing["ok"] is False and "is not there" in missing["reason"] and "--only B4" in missing["reason"]
    other = b4_results(tmp_path / "B4b" / "results.json", key="bonded_10000_h0.1")
    assert b5.b4_reference(70534, other)["ok"] is False and "no pressure run" in b5.b4_reference(70534, other)["reason"]
    (tmp_path / "bad.json").write_text("not json", encoding="utf-8")
    assert b5.b4_reference(70534, tmp_path / "bad.json")["ok"] is False


def test_the_case_file_carries_everything_the_script_needs(tmp_path):
    from bench import b5_prepomax as b5
    case = b5.case_inputs(tmp_path / "prepomax", tmp_path / "prepomax" / "job.inp", 1047.1975511965977, 20.0)
    assert case["unit_system"] == "MM_TON_S_C" and case["model_space"] == "PlaneStress" and case["step"] == "Step-1" and case["material"] == "M19" and case["element_set"] == "STEEL"
    assert case["thickness"] == 20.0 and case["omega"] == pytest.approx(2 * math.pi * 10_000 / 60) and 1 <= case["cpus"] <= paths.setting("femm", "workers")
    assert case["deck"].endswith("job.inp") and case["answer"].endswith("answer.json") and case["image"].endswith("view.png") and case["job_dir"].endswith("job")
    assert 0 < case["solver_timeout_s"] < b5.PMX_OWN_TIMEOUT_S < b5.PMX_TIMEOUT_S                      # the analysis' own limit first, then PrePoMax's, then run_solver's
    json.dumps(case)


def imported(**changes):
    got = {"parts": [{"name": "Shell_part-1", "type": "Shell", "nodes": 100, "elements": 40}], "node_sets": {"EDGE0": 5, "EDGE90": 5, "LONG_magnet_0_up": 9, "END_magnet_0_up": 3},
           "auto_created_node_sets": 12, "element_sets": {"STEEL": 40}, "materials": ["M19"], "sections": [{"name": "Section-1", "thickness": None}], "steps": ["Step-1"],
           "boundary_conditions": ["Displacement_rotation-1", "Displacement_rotation-2"], "loads": {"concentrated_force": 12}, "field_outputs": ["NF-Output-1"], "history_outputs": []}
    got.update(changes)
    return got


EXPECTED = {"nodes": 100, "elements": 40, "node_sets": ["EDGE0", "EDGE90", "LONG_magnet_0_up", "END_magnet_0_up"], "loaded_nodes": 6, "boundary_conditions": 2}


def test_the_importer_is_held_to_the_deck_and_what_it_leaves_out_is_named():
    """What PrePoMax made of B4's deck, against what the deck holds. The two gaps found on 2026-11-06 (the *DLOAD CENTRIF line, the section's thickness) are named."""
    from bench import b5_prepomax as b5
    got = b5.import_findings(imported(), EXPECTED)
    assert got["mesh_ok"] and got["sets_ok"] and got["forces_ok"] and got["bcs_ok"] and len(got["gaps"]) == 2
    assert "DLOAD CENTRIF" in got["gaps"][0] and "thickness" in got["gaps"][1]
    full = b5.import_findings(imported(loads={"concentrated_force": 12, "centrifugal": 1}, sections=[{"name": "S", "thickness": 20.0}]), EXPECTED)
    assert full["gaps"] == []
    wrong = b5.import_findings(imported(parts=[{"name": "p", "type": "Shell", "nodes": 99, "elements": 40}], loads={"concentrated_force": 11}, boundary_conditions=["a"],
                                        node_sets={"EDGE0": 5}), EXPECTED)
    assert not wrong["mesh_ok"] and not wrong["forces_ok"] and not wrong["bcs_ok"] and not wrong["sets_ok"]
    assert not b5.import_findings(imported(parts=[]), EXPECTED)["mesh_ok"]


def test_a_json_file_that_is_missing_or_not_json_is_a_value_error_that_says_which(tmp_path):
    from bench import b5_prepomax as b5
    (tmp_path / "answer.json").write_text('{"a": 1}', encoding="utf-8")
    assert b5.read_json(tmp_path / "answer.json", "PrePoMax's answer") == {"a": 1}
    with pytest.raises(ValueError, match="PrePoMax's report .* is not there"):
        b5.read_json(tmp_path / "report.json", "PrePoMax's report")
    (tmp_path / "bad.json").write_text("{", encoding="utf-8")
    with pytest.raises(ValueError, match="is not JSON"):
        b5.read_json(tmp_path / "bad.json", "PrePoMax's report")


def test_the_checks_say_so_when_prepomax_gave_no_answer():
    """If PrePoMax produced no answer.json the case still reports what it can (its version, its options, settings.bin) and one failed check names the reason."""
    from bench import b5_prepomax as b5
    import femm_queue
    before = {"sha256": "a" * 64, "bytes": 10, "solver": str(paths.tool("ccx")), "default_solver_type": 3, "calculix_exe": "x", "path": "p"}
    results = {"version": {"text": "PrePoMax 3.0.0.0", "parsed": "3.0.0.0", "returncode": 0}, "help": {"missing": {"manual": [], "help_only": []}, "returncode": 0},
               "settings": {"expected_solver": str(paths.tool("ccx"))}, "error": "PrePoMax's answer C:\\x\\answer.json is not there"}
    pmx = femm_queue.SolverRun(3, "", "", 1.0)
    got = b5.make_checks(results, None, None, None, before, dict(before), pmx, "3.0.0.0", b5.MANUAL_OPTIONS + b5.HELP_ONLY_OPTIONS, False)
    by = {c["name"]: c for c in got}
    assert by["PrePoMax starts from the command line and says its version"]["verdict"] == "PASS"
    assert by["the options: the manual's table and the program's own help agree, and the script options are there"]["verdict"] == "PASS"
    assert by["settings.bin is not changed, and its solver fields are the ones Task 1 set"]["verdict"] == "PASS"
    failed = by["PrePoMax runs the script and answers"]
    assert failed["verdict"] == "FAIL" and failed["gating"] is True and "answer.json is not there" in failed["measured"] and "exit code 3" in failed["measured"]
    changed = dict(before, sha256="b" * 64, solver="elsewhere.exe")
    got = b5.make_checks(results, None, None, None, before, changed, pmx, "3.0.0.0", [], False)
    assert {c["name"]: c for c in got}["settings.bin is not changed, and its solver fields are the ones Task 1 set"]["verdict"] == "FAIL"        # the file changed and so did its solver


class FakeLog:
    """What make_checks reads of a bench Log: the record of the last program run (the priority classes the watcher saw)."""
    runs = [{"priority_classes_seen": {"ccx_static.exe": ["below_normal"], "prepomax.exe": ["below_normal"]}}]


def good_run():
    """The inputs of make_checks for a PrePoMax run that went right: (results, answer, report, findings, before, pmx)."""
    import femm_queue
    from bench import b5_prepomax as b5
    ccx = str(paths.tool("ccx"))
    before = {"sha256": "a" * 64, "bytes": 10654, "solver": ccx, "default_solver_type": 3, "calculix_exe": "x", "path": "p"}
    findings = b5.import_findings(imported(loads={"concentrated_force": 12, "centrifugal": 1}, sections=[{"name": "S", "thickness": 20.0}]), EXPECTED)
    results = {"version": {"text": "PrePoMax 3.0.0.0", "parsed": "3.0.0.0", "returncode": 0}, "help": {"missing": {"manual": [], "help_only": []}, "returncode": 0},
               "settings": {"expected_solver": ccx}, "processes": b5.processes_left({}, {}), "import": {"expected": EXPECTED, **findings},
               "direct": {"peak_mpa": 33.454, "same_deck_as_prepomax": True}, "reference": {"ok": True, "peak": 33.454, "file": "B4/results.json", "nodes": 100}, "frd_peak_mpa": 33.454}
    answer = {"steps": [{"step": "import B4's deck", "ok": True, "seconds": 1.5}, {"step": "run CalculiX", "ok": True, "seconds": 7.0}], "imported": imported(),
              "job": {"status": "OK", "executable": ccx, "number_of_cpus": 4, "argument": "Job-1"}, "results": {"peak": {"value": 33.454, "node": 5}},
              "image": {"path": __file__, "bytes": 1000}}
    report = {"invocation": {"showGui": "No", "settingsSaved": False}, "exitCode": 0, "exitStatus": "Ok", "durationSeconds": 12.0, "prompts": []}
    return results, answer, report, findings, before, femm_queue.SolverRun(0, "CalculiX Version 2.23\n", "", 12.0)


def run_checks(results, answer, report, findings, before, pmx):
    """make_checks on those inputs, by name."""
    from bench import b5_prepomax as b5
    return {c["name"]: c for c in b5.make_checks(results, answer, report, findings, before, dict(before), pmx, "3.0.0.0", list(b5.MANUAL_OPTIONS + b5.HELP_ONLY_OPTIONS), True, log=FakeLog())}


WINDOW_CHECK = "PrePoMax runs without its window, and closes cleanly"


def test_the_checks_of_a_good_run_all_pass_so_that_the_failures_below_are_the_ones_made():
    by = run_checks(*good_run())
    assert len(by) == 13 and [c["name"] for c in by.values() if c["gating"] and c["verdict"] != "PASS"] == []
    assert by[WINDOW_CHECK]["verdict"] == "PASS" and "processes left by this run: none" in by[WINDOW_CHECK]["measured"]


def test_no_process_left_fails_with_the_reason_when_the_process_table_could_not_be_read():
    """Fix round 1, item 6: a gate that cannot look must not pass. An unreadable table fails the check, with the reason; a process of the run that is left still fails it."""
    from bench import b5_prepomax as b5
    results, answer, report, findings, before, pmx = good_run()
    results["processes"] = b5.processes_left(None, {})
    item = run_checks(results, answer, report, findings, before, pmx)[WINDOW_CHECK]
    assert item["verdict"] == "FAIL" and item["gating"] is True and "could not be read" in item["measured"] and "before the run" in item["measured"]
    results["processes"] = b5.processes_left({}, {300: "PrePoMax.exe"})
    left = run_checks(results, answer, report, findings, before, pmx)[WINDOW_CHECK]
    assert left["verdict"] == "FAIL" and "300" in left["measured"]


def test_the_answer_of_a_script_that_failed_its_import_is_digested_without_a_key_error():
    """Fix round 1, item 5: answer.json is written in a `finally`, so a failed step leaves an answer without the keys of the later steps; the case digests it as findings = None."""
    from bench import b5_prepomax as b5
    assert b5.findings_of({"steps": [{"step": "import B4's deck", "ok": False}]}, EXPECTED) is None
    assert b5.findings_of({"steps": [], "imported": imported()}, EXPECTED)["mesh_ok"] is True


def test_a_script_that_failed_a_step_gives_a_fail_that_names_it_and_not_a_key_error():
    """Fix round 1, item 5: the failed step of answer.json and its error are in a failed check; the checks that need what the later steps would have made are not made."""
    results, answer, report, findings, before, pmx = good_run()
    answer = {"steps": [{"step": "import B4's deck", "ok": False, "seconds": 0.4, "error": "RuntimeError: Failed to import the file"}], "error": "RuntimeError: Failed to import the file"}
    report = dict(report, exitCode=3, exitStatus="ScriptError")
    results["import"] = {"expected": EXPECTED}
    by = run_checks(results, answer, report, None, before, pmx)                 # raised KeyError: 'imported' before the fix
    item = by["PrePoMax's script made what the checks below read"]
    assert item["verdict"] == "FAIL" and item["gating"] is True
    assert "import B4's deck" in item["measured"] and "Failed to import the file" in item["measured"] and "imported" in item["measured"] and "results" in item["measured"]
    assert by[WINDOW_CHECK]["verdict"] == "FAIL"
    assert "PrePoMax's peak von Mises stress reproduces the same deck solved by ccx directly" not in by
    assert by["PrePoMax starts from the command line and says its version"]["verdict"] == "PASS"             # what does not need the answer is still reported
    bare = run_checks(results, answer, {"exitCode": 3}, None, before, pmx)                                  # and a report without its "invocation" is a FAIL too, not a KeyError
    assert bare[WINDOW_CHECK]["verdict"] == "FAIL"


# ----------------------------------------------------------------------------------------------- 5. the picture
def synthetic_viewport(path, scale_bar=True):
    """A PrePoMax viewport in miniature (RGBA, transparent where nothing is drawn): the legend at the top left, the model in the middle, an info box at the top right, the triad at the
    bottom right and, far below the model, the scale bar."""
    from PIL import Image
    image = np.zeros((600, 1000, 4), np.uint8)
    image[30:200, 20:120] = (200, 30, 30, 255)                       # the legend: 100 x 170 at the left
    image[100:420, 400:700] = (30, 120, 220, 255)                    # the model: 300 x 320
    image[15:70, 760:990] = (250, 250, 250, 255)                     # the info box
    image[520:580, 930:990] = (0, 0, 0, 255)                         # the triad
    if scale_bar:
        image[520:540, 420:700] = (0, 0, 0, 255)                     # the scale bar, 100 rows under the model
    Image.fromarray(image).save(path)


def test_the_prepomax_picture_keeps_the_legend_and_the_model_and_leaves_out_the_scale_bar_the_info_box_and_the_triad(tmp_path):
    from PIL import Image
    from bench import b5_prepomax as b5
    synthetic_viewport(tmp_path / "view.png")
    done = b5.prepomax_picture(tmp_path / "view.png", "a caption\nin two lines", tmp_path / "figure.png")
    assert done["made"] is True and done["size"] == (1000, 600)
    assert done["legend"] == [170 - 1 + 17, 100 - 1 + 17] and done["model"] == [320 - 1 + 17, 300 - 1 + 17]          # the boxes, with 8 pixels of margin before the first row and 9 after the last
    assert Image.open(tmp_path / "figure.png").size[0] > 400
    synthetic_viewport(tmp_path / "view2.png", scale_bar=False)
    assert b5.prepomax_picture(tmp_path / "view2.png", "x", tmp_path / "figure2.png")["model"] == done["model"]                  # the bar made no difference: it was cut away either way


def test_the_caption_has_the_numbers_of_the_model_PrePoMax_solved_and_no_line_that_wraps():
    """The caption says the steel's node count (the model PrePoMax solved), not the mesh's (70,534 with the slots'), and each of its lines fits the width a title is wrapped to: a
    longer one wraps and leaves a word alone on a line (the first picture of the run had three of those)."""
    import bench
    from bench import b5_prepomax as b5
    text = b5.picture_caption(45_271, 33.454139709, 527)
    assert "45,271 nodes" in text and "70,534" not in text and "Peak 33.45 MPa at node 527" in text and "10,000 rpm" in text and "CalculiX 2.23" in text
    lines = text.split("\n")
    assert len(lines) == 4 and all(len(line) <= bench.TITLE_CHARS for line in lines)
    assert bench.wrap_title(text) == text                                       # nothing wraps


def test_a_picture_without_a_legend_and_a_model_is_not_made(tmp_path):
    from PIL import Image
    from bench import b5_prepomax as b5
    nothing = np.zeros((100, 200, 4), np.uint8)
    Image.fromarray(nothing).save(tmp_path / "empty.png")
    assert b5.prepomax_picture(tmp_path / "empty.png", "x", tmp_path / "f.png")["made"] is False and not (tmp_path / "f.png").exists()
    one = nothing.copy()
    one[10:60, 10:150] = (1, 2, 3, 255)
    Image.fromarray(one).save(tmp_path / "one.png")
    assert b5.prepomax_picture(tmp_path / "one.png", "x", tmp_path / "g.png")["made"] is False


def test_the_runs_of_a_row_of_flags_are_merged_when_they_are_close():
    from bench import b5_prepomax as b5
    flags = [False] * 5 + [True] * 3 + [False] * 10 + [True] * 2 + [False] * 50 + [True]
    assert b5._runs(flags, gap=20) == [(5, 19), (70, 70)]
    assert b5._runs(flags, gap=5) == [(5, 7), (18, 19), (70, 70)] and b5._runs([], gap=5) == []


# ----------------------------------------------------------------------------------------------- 6. the script that runs inside PrePoMax
def test_the_prepomax_side_script_is_valid_python_with_the_standard_library_only_and_never_starts_a_program():
    """It runs inside PrePoMax's own Python, where `pmx` is already there (so it is not imported here), and it imports the standard library only."""
    source = (P / "bench" / "b5_prepomax_script.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    compile(source, "b5_prepomax_script.py", "exec")
    imported_names = {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names} | {(n.module or "").split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom)}
    assert imported_names <= {"json", "os", "sys", "time", "traceback"}, imported_names
    assert "pmx" not in imported_names


def test_the_script_takes_every_step_the_case_needs_and_writes_its_answer_even_when_a_step_fails():
    """The steps, as calls of PrePoMax's API: the unit system and model space, the import, the read-back, the section with its thickness, the centrifugal load, the analysis and its run,
    the result's extremes, the picture; and the answer is written in a `finally`, so that a failure leaves its reason."""
    source = (P / "bench" / "b5_prepomax_script.py").read_text(encoding="utf-8")
    for call in ("pmx.import_file", "pmx.model.properties.unit_system", "pmx.model.properties.model_space", "pmx.model.sections.remove", "pmx.model.sections.add_solid",
                 "model_step.loads.add_centrifugal", "pmx.jobs.add", "job.run", "pmx.jobs.open_results", ".get_extremes(\"MISES\")", "pmx.view.save_image", "show_min_max(False, False)"):
        assert call in source, call
    tree = ast.parse(source)
    main = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "main")
    tries = [n for n in ast.walk(main) if isinstance(n, ast.Try)]
    assert any(t.finalbody and "json.dump" in ast.unparse(t.finalbody) for t in tries), "answer.json is written in a finally block"
    assert any(any(isinstance(h.type, ast.Name) and h.type.id == "Exception" for h in t.handlers) and "raise" in ast.unparse(t.handlers[0]) for t in tries)       # the error is kept and let through
    last = tree.body[-1]                                                        # the call is guarded, so that the file can be imported (and tested) where there is no `pmx`
    assert isinstance(last, ast.If) and ast.unparse(last.test) == "'pmx' in globals()" and ast.unparse(last.body[0]) == "main(sys.argv[1])" and not last.orelse
    assert "subprocess" not in source and "os.system" not in source


def load_script(monkeypatch):
    """The PrePoMax-side script as a module, executed here, where there is no `pmx` (PrePoMax puts it into the globals of the script it runs; nothing else of PrePoMax is needed
    to import the file)."""
    import importlib.util
    import sys
    monkeypatch.setattr(sys, "argv", ["b5_prepomax_script.py"])                  # PrePoMax's --scriptArgs puts the case file here; an import has none
    spec = importlib.util.spec_from_file_location("b5_prepomax_script_under_test", P / "bench" / "b5_prepomax_script.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_script_can_be_imported_where_there_is_no_pmx_and_then_runs_nothing(monkeypatch):
    """Fix round 1, item 10: the file ended with an unconditional main(sys.argv[1]), so importing it raised (no argv[1], no `pmx`); it now runs main only where `pmx` exists."""
    module = load_script(monkeypatch)
    assert callable(module.main) and callable(module.run_step) and "pmx" not in vars(module)


def test_a_step_of_the_script_records_how_long_it_took_and_what_went_wrong(monkeypatch):
    """Fix round 1, item 10: a test of what each step writes into answer.json (a count of the word "seconds" in the source said nothing of that)."""
    module = load_script(monkeypatch)
    answer = {"steps": []}
    assert module.run_step(answer, "a step", lambda: 42) == 42
    first = answer["steps"][0]
    assert set(first) == {"step", "ok", "seconds"} and first["step"] == "a step" and first["ok"] is True and 0 <= first["seconds"] < 5

    def fails():
        raise ValueError("boom")
    with pytest.raises(ValueError, match="boom"):
        module.run_step(answer, "a bad step", fails)
    second = answer["steps"][1]
    assert second["step"] == "a bad step" and second["ok"] is False and second["error"] == "ValueError: boom" and 0 <= second["seconds"] < 5


def test_main_writes_its_answer_when_a_step_fails_and_lets_the_error_through(monkeypatch, tmp_path):
    """The first step of the script reaches for pmx.model; here pmx is a stand-in whose model does not exist. The answer is written all the same, with the failed step and the error,
    and the error goes on (PrePoMax then exits with 3, a script error)."""
    module = load_script(monkeypatch)

    class Broken:
        @property
        def model(self):
            raise RuntimeError("no model")
    module.pmx = Broken()
    (tmp_path / "case.json").write_text(json.dumps({"answer": str(tmp_path / "answer.json"), "unit_system": "MM_TON_S_C", "model_space": "PlaneStress"}), encoding="utf-8")
    with pytest.raises(RuntimeError, match="no model"):
        module.main(str(tmp_path / "case.json"))
    answer = json.loads((tmp_path / "answer.json").read_text(encoding="utf-8"))
    assert [s["ok"] for s in answer["steps"]] == [False] and "RuntimeError: no model" in answer["steps"][0]["error"] and "no model" in answer["error"] and "imported" not in answer


def test_the_deck_b4_writes_is_the_one_prepomax_can_import(tmp_path):
    """B4's deck writer puts the boundary conditions inside the step and writes the axis of the centrifugal load as whole numbers: both found on 2026-11-06 with PrePoMax's importer
    (it read no boundary condition before the step, and could not read '0.')."""
    import motorgeo as mg
    from bench import b4_ccx as b4
    p = mg.Params(barrier="none")
    model = b4.build_model(p, 1.0)
    text = b4.write_deck(tmp_path, model, "pressure", 10_000).read_text(encoding="ascii")
    lines = text.splitlines()
    assert lines.index("*STEP") < lines.index("*BOUNDARY") < lines.index("*DLOAD")
    assert [l for l in lines if l.startswith("STEEL, CENTRIF")][0].split(", ")[3:] == ["0", "0", "0", "0", "0", "1"]
    assert any(l.startswith("*NSET, NSET=EDGE0") for l in lines) and any(l.startswith("*NSET, NSET=LONG_magnet_0_up") for l in lines) and "*ELEMENT, TYPE=CPS6, ELSET=STEEL" in lines


def test_the_case_is_registered_with_run_bench_between_b4_and_b6():
    from bench import run_bench
    assert list(run_bench.CASES)[3:6] == ["B4", "B5", "B6"] and run_bench.CASES["B5"] == "b5_prepomax"          # Task 13 added B7 to B9 after them
