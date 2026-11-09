"""Tests for validation/hand_reference_check.py: does it make the same three changes as validate_full.hand_worker?

The check is only as good as its copy of those changes, so there are two kinds of test here.

1. Unit tests apply the changes to a tiny model and look at each one. The real files (the .FEM and a hand_000.fem that
   validate_full.py made) are not needed.
2. Tie tests read validate_full.hand_worker itself (its source, parsed with `ast`, so comments and docstrings do not
   count) and fail the day it stops doing one of the three things, does a fourth thing to the copy, or turns by another
   formula. They are what keeps the second copy of the recipe, in hand_reference_check, honest: if one fails, change
   apply_hand_worker to follow hand_worker, then change the expected lines below.

Run from P:   uv run pytest tests/test_hand_reference_check.py -q
"""
import ast
import inspect
import math
import textwrap

import pytest

import hand_reference_check as hrc
import validate_full as vf
from fem_text_tools import make_fem, write


def tiny_model():
    """A parsed model (the shape same_model.parse_fem returns) with one node and one label in each group."""
    fields = lambda h_c: {"fields": {"H_c": h_c}, "rows": []}
    return {
        "header": {}, "order": {}, "crlf": True,
        "props": {"BlockProps": {"Air": fields(0.0), "N45": fields(1030272.0)},
                  "CircuitProps": {n: {"fields": {"TotalAmps_re": 1.0}, "rows": []} for n in ("Phase A", "Phase B", "Phase C")}},
        "nodes": [{"x": 10.0, "y": 0.0, "prop": "", "group": 1, "lfs": -1.0},          # rotor: turns
                  {"x": 50.0, "y": 0.0, "prop": "", "group": 0, "lfs": -1.0}],         # stator: stays
        "segments": [], "arcs": [],
        "labels": [{"x": 5.0, "y": 0.0, "material": "N45", "group": 1, "magdir": 90.0},        # a magnet in the rotor
                   {"x": 7.0, "y": 0.0, "material": "Air", "group": 1, "magdir": 0.0},         # rotor air
                   {"x": hrc.AIR_GAP_LABEL[0], "y": hrc.AIR_GAP_LABEL[1], "material": "Air", "group": 1, "magdir": 0.0},
                   {"x": 55.0, "y": 0.0, "material": "Air", "group": 0, "magdir": 0.0}],       # stator air
    }


@pytest.fixture
def axes_65(monkeypatch):
    """validate_full.fem_axes() reads the .FEM; give it the numbers of the real one (d-axis 0.5, phase A 65.025 deg).

    Also puts validate_full.FEM back after a test that calls main(), which points it at the test's own files."""
    monkeypatch.setattr(vf, "fem_axes", lambda: (0.5, 65.025))
    monkeypatch.setattr(vf, "FEM", vf.FEM)


# ----------------------------------------------------------------------------- 1. the changes, on a tiny model
def test_the_job_turns_the_rotor_and_its_magnets_only(axes_65):
    model = tiny_model()
    out = hrc.apply_hand_worker(model, theta_e=0.0, i_a=0.0)
    move = 65.025 - 0.5
    a = math.radians(move)
    rotor, stator = out["nodes"]
    assert (rotor["x"], rotor["y"]) == pytest.approx((10 * math.cos(a), 10 * math.sin(a)))      # counter-clockwise
    assert (stator["x"], stator["y"]) == (50.0, 0.0)
    magnet, air, gap, stator_air = out["labels"]
    assert (magnet["x"], magnet["y"]) == pytest.approx((5 * math.cos(a), 5 * math.sin(a)))
    assert magnet["magdir"] == pytest.approx(90.0 + move)           # a magnet's direction turns with it
    assert air["magdir"] == 0.0                                      # other labels have no direction to turn
    assert stator_air["x"] == 55.0
    assert model["nodes"][0]["x"] == 10.0, "the argument must be left as it was"


def test_the_air_gap_label_leaves_group_1_and_stays_where_it_is(axes_65):
    out = hrc.apply_hand_worker(tiny_model(), theta_e=0.0, i_a=0.0)
    gap = out["labels"][2]
    assert gap["group"] == 0
    assert (gap["x"], gap["y"]) == hrc.AIR_GAP_LABEL                 # not turned: it left the group before the turn


def test_the_job_sets_the_three_currents_and_the_angle_depends_on_theta_e(axes_65):
    out = hrc.apply_hand_worker(tiny_model(), theta_e=60.0, i_a=10.0)
    amps = {n: b["fields"]["TotalAmps_re"] for n, b in out["props"]["CircuitProps"].items()}
    assert amps == {"Phase A": 10.0, "Phase B": -5.0, "Phase C": -5.0}
    move = 65.025 + 60.0 / vf.P.pole_pairs - 0.5                      # 4 pole pairs: 60 deg electrical = 15 deg mechanical
    assert out["labels"][0]["magdir"] == pytest.approx(90.0 + move)


def test_job_parameters_follow_validate_fulls_list():
    assert hrc.job_parameters(0) == (0.0, 0.0)                        # no load, theta_e = 0
    assert hrc.job_parameters(1) == (15.0, 0.0)
    assert hrc.job_parameters(24) == (0.0, vf.I_LOAD)                 # the second half is at load
    assert hrc.job_parameters(47) == (345.0, vf.I_LOAD)
    for bad in (-1, 48):
        with pytest.raises(ValueError, match="hand jobs 0 to 47"):
            hrc.job_parameters(bad)


def test_a_reference_without_what_the_job_changes_is_reported_not_crashed_on(axes_65):
    no_label = tiny_model()
    no_label["labels"][2]["x"] += 1.0                                # the air-gap label is not where the job looks
    with pytest.raises(ValueError, match="no label at the air-gap point"):
        hrc.apply_hand_worker(no_label, 0.0, 0.0)
    no_circuit = tiny_model()
    del no_circuit["props"]["CircuitProps"]["Phase C"]
    with pytest.raises(ValueError, match='no circuit called "Phase C"'):
        hrc.apply_hand_worker(no_circuit, 0.0, 0.0)


def test_the_command_refuses_what_it_cannot_read_and_says_when_the_reference_is_not_the_model(tmp_path, capsys, axes_65):
    good = write(tmp_path, "good.fem", make_fem())
    empty = tmp_path / "empty.fem"
    empty.write_bytes(b"")
    assert hrc.main([str(empty), "--reference", str(good)]) == 2
    assert "could not read HAND" in capsys.readouterr().err
    assert hrc.main([str(good), "--reference", str(empty)]) == 2
    assert "could not read REFERENCE" in capsys.readouterr().err
    with pytest.raises(SystemExit) as stop:
        hrc.main([str(good), "--job", "48", "--reference", str(good)])
    assert stop.value.code == 2 and "hand jobs 0 to 47" in capsys.readouterr().err
    # make_fem's model has no label at the air-gap point: that is a reference hand_worker was not written for
    assert hrc.main([str(good), "--reference", str(good)]) == 1
    assert "no label at the air-gap point" in capsys.readouterr().out


# ----------------------------------------------------------------------------- 2. tied to validate_full.hand_worker
def source_tree(thing):
    return ast.parse(textwrap.dedent(inspect.getsource(thing)))


def femm_calls(function):
    """Every femm.<name>(...) call in `function`, in the order they are written, as (name, arguments as source text)."""
    found = []
    for node in ast.walk(source_tree(function)):
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and isinstance(node.func.value, ast.Name) and node.func.value.id == "femm"):
            found.append((node.lineno, node.col_offset, node.func.attr, ", ".join(ast.unparse(a) for a in node.args)))
    return [(name, args) for _, _, name, args in sorted(found)]


def assigned(tree, name):
    """The expression assigned to the plain name `name` (the first one, wherever it sits)."""
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return node.value
    raise AssertionError(f"{name} is not assigned anywhere in the source")


def same_code(node, text):
    """True if the parsed `node` is the same expression as the source `text` (spacing and brackets do not matter)."""
    return ast.dump(node) == ast.dump(ast.parse(text, mode="eval").body)


# The calls of hand_worker that change the copy, in order: what hand_reference_check.apply_hand_worker has to repeat.
THE_THREE_CHANGES = [
    ("mi_selectlabel", "*AIR_GAP_LABEL"),    # 1. the air-gap label ...
    ("mi_setgroup", "0"),                    #    ... goes to group 0
    ("mi_selectgroup", "1"),                 # 2. the rotor (group 1) ...
    ("mi_moverotate", "0, 0, move"),         #    ... turns about the origin by `move`
    ("mi_modifycircprop", "ph, 1, i"),       # 3. each circuit gets its current
]
CALLS_THAT_LEAVE_THE_COPY_AS_IT_IS = {"mi_clearselected", "mi_saveas", "mi_analyze", "mi_loadsolution"}


def test_hand_worker_still_makes_the_three_changes_and_in_this_order():
    calls = femm_calls(vf.hand_worker)
    changes = [c for c in calls if c[0] in {name for name, _ in THE_THREE_CHANGES}]
    assert changes == THE_THREE_CHANGES, (
        "validate_full.hand_worker no longer does what hand_reference_check.apply_hand_worker copies. Make the same "
        f"change there (and in THE_THREE_CHANGES). It now calls {changes}")
    assert str(hrc.ROTOR_GROUP) == dict(THE_THREE_CHANGES)["mi_selectgroup"]


def test_hand_worker_does_not_change_the_copy_in_a_fourth_way():
    names = {name for name, _ in femm_calls(vf.hand_worker) if name.startswith("mi_")}
    unknown = names - {name for name, _ in THE_THREE_CHANGES} - CALLS_THAT_LEAVE_THE_COPY_AS_IT_IS
    assert not unknown, (f"validate_full.hand_worker now calls {sorted(unknown)}. If that changes the model, "
                         "hand_reference_check.apply_hand_worker has to make the same change; then list the call here")


def test_hand_worker_still_turns_by_the_same_angle_and_sets_the_same_currents():
    tree = source_tree(vf.hand_worker)
    assert same_code(assigned(tree, "move"), "(a0 + theta_e / P.pole_pairs) - d0")
    loops = [node.iter for node in ast.walk(tree) if isinstance(node, ast.For)]
    assert any(same_code(it, '(("Phase A", i_a), ("Phase B", -i_a / 2), ("Phase C", -i_a / 2))') for it in loops)


def test_the_list_of_hand_jobs_is_still_the_one_job_parameters_assumes():
    """The jobs are built in validate_full's __main__ block; job K of the check must be job K of the study."""
    jobs = assigned(ast.parse(inspect.getsource(vf)), "hand_jobs")
    assert same_code(jobs, "[(k, th, i) for k, (th, i) in enumerate((th, i) for i in (0.0, I_LOAD) for th in THETAS)]")


def test_the_air_gap_label_is_one_constant_that_validate_full_and_the_check_share():
    assert vf.AIR_GAP_LABEL == (45.398271307113383, 0.39618471162615876)      # the point the 10-28 copies were made with
    assert hrc.AIR_GAP_LABEL is vf.AIR_GAP_LABEL
