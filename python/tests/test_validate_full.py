"""Tests for validate_full.py's job lists: no two jobs may ever share a file name.

Until 2026-11-04 the slice jobs were numbered 100 + k with k = 0..3 for the no-load jobs and again for the 10 A jobs, so
eight jobs had four names (slice_100 ... slice_103) and two FEMMs saved, solved and read each file at the same time.
Run that way validate_full.py crashed fkn.exe or FEMM, hung a FEMM at 100 % CPU, or exited 0 with wrong rows (a
flux linkage of exactly 0, or a 10 A value in a no-load row). Nothing here starts FEMM: the job lists are built inside
validate_full's __main__ block, so the tests read its source with `ast` (like test_hand_reference_check.py does for
hand_jobs) and evaluate the expressions exactly as they are written there.

Run from P:   uv run pytest tests/test_validate_full.py -q
"""
import ast
import inspect

import validate_full as vf


def assigned_expression(name):
    """The expression that validate_full's source assigns to `name` (anywhere, so also inside the __main__ block)."""
    for node in ast.walk(ast.parse(inspect.getsource(vf))):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                and node.targets[0].id == name):
            return node.value
    raise AssertionError(f"validate_full.py no longer assigns {name}")


def job_list(name):
    expression = ast.Expression(assigned_expression(name))
    ast.fix_missing_locations(expression)
    return eval(compile(expression, "validate_full.py", "eval"), dict(vars(vf)))


def test_no_two_slice_jobs_share_a_file_name():
    jobs = job_list("slice_jobs")
    ids = [job[0] for job in jobs]
    assert len(jobs) == 8
    assert len(set(ids)) == len(ids), f"slice job ids {ids}: two jobs would write the same slice_NNN.fem and .ans"


def test_no_two_hand_jobs_share_a_file_name():
    ids = [job[0] for job in job_list("hand_jobs")]
    assert len(ids) == 48 and len(set(ids)) == 48


def test_the_workers_name_their_files_after_the_job_id_and_nothing_else():
    """The ids are only unique files if the workers build the path from the id alone."""
    assert 'f"slice_{k:03d}.fem"' in inspect.getsource(vf.slice_worker)
    assert 'f"hand_{k:03d}.fem"' in inspect.getsource(vf.hand_worker)


def test_the_slice_jobs_still_cover_every_angle_at_both_currents_exactly_once():
    """Giving the jobs their own files must not change what is computed: the same 48 (angle, current) points."""
    wanted = sorted((float(th), i) for i in (0.0, vf.I_LOAD) for th in vf.THETAS)
    solved = sorted((float(th), i) for _, thetas, i in job_list("slice_jobs") for th in thetas)
    assert solved == wanted and len(solved) == 48


def test_the_main_block_checks_that_job_ids_are_unique():
    """The guard that stops the mistake from coming back unnoticed."""
    source = inspect.getsource(vf)
    assert "two jobs would share one file name" in source
