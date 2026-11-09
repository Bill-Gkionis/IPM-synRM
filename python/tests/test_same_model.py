"""Tests for validation/same_model.py: the tool that says whether two FEMM files hold the same model.

A "same model" verdict is only worth something if the tool also says "different" when the models differ, and refuses
what it cannot read. So most of these tests build a tiny model as text (no FEMM needed), change one thing at a time,
and check that the tool names it; that things which are not the model (how a number is written, the order of storage,
the comment) do not count; and that text which is not a model, a file that stops in the middle, a section the tool does
not know and a tolerance of zero are all dealt with out loud. Two tests (marked `femm`) start a real FEMM 5.0, save
models with it and read them back, so the tool is also tried on what FEMM really writes: they need no file that git
does not hold.

Run from P:   uv run pytest tests/test_same_model.py -q
"""
import re
from pathlib import Path

import pytest

import femm
import femm_session
import femm_test_tools as tools
import motorgeo as mg
import same_model as sm
from fem_text_tools import make_fem, write

OUT = Path(__file__).resolve().parent.parent / "out" / "tests" / "same_model"


def verdict(tmp_path, old_text, new_text, rtol=1e-9):
    old, new = write(tmp_path, "old.fem", old_text), write(tmp_path, "new.fem", new_text)
    return sm.compare(sm.parse_fem(old), sm.parse_fem(new), rtol)


def steel(h):
    """The default materials with the last B-H value of the steel set to h."""
    return [("Air", {"Mu_x": 1, "Sigma": 0}, []),
            ("M-19 Steel", {"Mu_x": 4416, "Sigma": 1.9}, [(0, 0), (0.05, 15.120714), (1.5, h)])]


# ----------------------------------------------------------------------------- what is the same, and what is not
def test_the_same_file_is_the_same_model(tmp_path):
    result = verdict(tmp_path, make_fem(), make_fem())
    assert result["differences"] == [] and result["formatting"] == [] and result["warnings"] == []
    assert result["compared"]["nodes"] == 3 and result["compared"]["arcs"] == 1
    assert result["compared"]["labels (holes included)"] == 3          # two labels and the hole
    assert result["compared"]["B-H points"] == 3 and result["compared"]["materials"] == 2


def test_how_a_number_is_written_is_formatting_not_the_model(tmp_path):
    new = make_fem(precision="1e-008", comment="Add comments here.", eol="\n",
                   materials=[("Air", {"Mu_x": 1.0, "Sigma": 0.0}, []),
                              ("M-19 Steel", {"Mu_x": 4416, "Sigma": 1.8999999999999999},
                               [(0, 0), (0.050000000000000003, 15.120714), (1.5, 1108.3255690000001)])])
    result = verdict(tmp_path, make_fem(), new)
    assert result["differences"] == []
    text = "\n".join(result["formatting"])
    assert "[Precision]" in text and "[Comment]" in text and "line endings" in text


def test_a_value_femm_computes_itself_is_listed_but_is_not_the_model(tmp_path):
    # FEMM 4.2 stores the side length it worked out for an arc (0.34), FEMM 5.0 stores 0 = "let the mesher choose"
    new = make_fem(arcs=[(1, 2, 90.0, 0.5, 2, 0, 0, 0.0)])
    result = verdict(tmp_path, make_fem(arcs=[(1, 2, 90.0, 0.5, 2, 0, 0, 0.34)]), new)
    assert result["differences"] == []
    assert len(result["formatting"]) == 1 and "mySideLength" in result["formatting"][0]
    assert "0.34 in OLD and 0 in NEW" in result["formatting"][0]


def test_the_same_things_in_another_order_are_the_same_model(tmp_path):
    # nodes reversed (segment and arc numbers follow), labels swapped, materials listed the other way round
    new = make_fem(nodes=[(0.0, 10.0, 0, 0), (10.0, 0.0, 0, 0), (0.0, 0.0, 0, 0)],
                   segments=[(2, 1, -1, 1, 0, 0)], arcs=[(1, 0, 90.0, 0.5, 2, 0, 0, 0.5)],
                   labels=[(6.0, 1.0, 2, -1, 1, 90, 2, 10, 0), (3.0, 3.0, 1, 1.5, 0, 0, 1, 1, 0)],
                   materials=[("M-19 Steel", {"Mu_x": 4416, "Sigma": 1.9}, [(0, 0), (0.05, 15.120714), (1.5, 1108.325569)]),
                              ("Air", {"Mu_x": 1, "Sigma": 0}, [])])
    result = verdict(tmp_path, make_fem(), new)
    assert result["differences"] == []
    text = "\n".join(result["formatting"])
    assert "nodes: stored in a different order" in text and "materials: listed in a different order" in text


@pytest.mark.parametrize("what, new, expect", [
    ("a node moved by 1 micrometre", dict(nodes=[(0.0, 0.0, 0, 0), (10.001, 0.0, 0, 0), (0.0, 10.0, 0, 0)]), "node"),
    ("a node in another group", dict(nodes=[(0.0, 0.0, 0, 1), (10.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)]), "group"),
    ("a segment on another boundary", dict(segments=[(0, 1, -1, 2, 0, 0)]), "boundary"),
    ("a segment with another mesh size", dict(segments=[(0, 1, 0.5, 1, 0, 0)]), "maxside"),
    ("an arc with another sweep", dict(arcs=[(1, 2, 89.0, 0.5, 2, 0, 0, 0.5)]), "sweep"),
    ("an arc run the other way", dict(arcs=[(2, 1, 90.0, 0.5, 2, 0, 0, 0.5)]), "arc"),
    ("a label with another material", dict(labels=[(3.0, 3.0, 1, 1.5, 0, 0, 1, 1, 0), (6.0, 1.0, 1, -1, 1, 90, 2, 10, 0)]), "material"),
    ("a label with another magnetisation direction", dict(labels=[(3.0, 3.0, 2, 1.5, 0, 0, 1, 1, 0), (6.0, 1.0, 1, -1, 1, 85, 2, 10, 0)]), "magdir"),
    ("a label with another number of turns", dict(labels=[(3.0, 3.0, 2, 1.5, 0, 0, 1, 1, 0), (6.0, 1.0, 1, -1, 1, 90, 2, -10, 0)]), "turns"),
    ("a label moved", dict(labels=[(3.0, 3.0, 2, 1.5, 0, 0, 1, 1, 0), (6.5, 1.0, 1, -1, 1, 90, 2, 10, 0)]), "label"),
    ("a label with a magnetisation function", dict(labels=[(3.0, 3.0, 2, 1.5, 0, 0, 1, 1, 0),
                                                             (6.0, 1.0, 1, -1, 1, 90, 2, 10, 0, '"cos(x) + 1"')]), "magdirfctn"),
    ("a node with a local size limit", dict(nodes=[(0.0, 0.0, 0, 0, 0.2), (10.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)]), "lfs"),
    ("a hole in another group", dict(holes=[(1.0, 8.0, 1)]), "group"),
    ("a hole missing (holes are labels without a material)", dict(holes=[]), "label"),
    ("a circuit with another current", dict(circuits=[("A", 6.5)]), "TotalAmps_re"),
    ("a circuit renamed", dict(circuits=[("B", 6.0)]), "circuit"),
    ("one B-H value changed in the 8th digit", dict(materials=steel(1108.3257)), "B-H point 2"),
    ("a B-H point missing", dict(materials=[
        ("Air", {"Mu_x": 1, "Sigma": 0}, []),
        ("M-19 Steel", {"Mu_x": 4416, "Sigma": 1.9}, [(0, 0), (0.05, 15.120714)])]), "B-H points"),
    ("another lamination fill", dict(materials=[
        ("Air", {"Mu_x": 1, "Sigma": 0}, []),
        ("M-19 Steel", {"Mu_x": 4416, "Sigma": 1.9, "LamFill": 0.955}, [(0, 0), (0.05, 15.120714), (1.5, 1108.325569)])]), "LamFill"),
    ("a material missing", dict(materials=[("Air", {"Mu_x": 1, "Sigma": 0}, [])], labels=[(3.0, 3.0, 0, 1.5, 0, 0, 1, 1, 0), (6.0, 1.0, 1, -1, 1, 90, 2, 10, 0)]), "M-19 Steel"),
])
def test_a_change_in_the_model_is_named(tmp_path, what, new, expect):
    result = verdict(tmp_path, make_fem(), make_fem(**new))
    assert result["differences"], f"{what}: the tool saw no difference"
    assert any(expect in line for line in result["differences"]), (what, result["differences"])


def test_one_thing_is_called_by_its_own_word_not_by_a_chopped_plural(tmp_path):
    """The words in the lines are 'boundary' and 'point property', never the chopped plurals 'boundarie' and 'propertie'."""
    old = make_fem(pointprops=[("P1", 2)], nodes=[(0.0, 0.0, 1, 0), (10.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)])
    new = make_fem(pointprops=[("P1", 3)], nodes=[(0.0, 0.0, 1, 0), (10.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)],
                   bdry=[("A=0", 1.0), ("pbc1", 4.0)], circuits=[("A", 7.0)], materials=steel(1108.325569)[:1] + [
                       ("M-19 Steel", {"Mu_x": 5000, "Sigma": 1.9}, [(0, 0), (0.05, 15.120714), (1.5, 1108.325569)])])
    lines = verdict(tmp_path, old, new)["differences"]
    assert 'point property "P1": I_re is 2 in OLD and 3 in NEW' in lines
    assert 'boundary "A=0": BdryType is 0 in OLD and 1 in NEW' in lines
    assert 'boundary "pbc0" is only in OLD' in lines and 'boundary "pbc1" is only in NEW' in lines
    assert 'circuit "A": TotalAmps_re is 6 in OLD and 7 in NEW' in lines
    assert 'material "M-19 Steel": Mu_x is 4416 in OLD and 5000 in NEW' in lines
    assert not re.search(r"boundarie\b|propertie\b", "\n".join(lines)), lines
    # and a plural count is a plural: the 'compared' line says boundaries and point properties
    counts = verdict(tmp_path, old, old)["compared"]
    assert counts["boundaries"] == 2 and counts["point properties"] == 1


def test_the_problem_definition_is_compared(tmp_path):
    for old_line, new_line, name in (("[Depth]       =  20", "[Depth]       =  25", "Depth"),
                                     ("[LengthUnits] =  millimeters", "[LengthUnits] =  meters", "LengthUnits"),
                                     ("[Precision]   =  1e-08", "[Precision]   =  1e-06", "Precision")):
        text = make_fem()
        assert old_line in text
        result = verdict(tmp_path, text, text.replace(old_line, new_line))
        assert any(f"[{name}]" in line for line in result["differences"]), name


# ----------------------------------------------------------------------------- tolerances
def test_tolerance_decides_what_counts_as_the_same_number(tmp_path):
    tiny = make_fem(nodes=[(0.0, 0.0, 0, 0), (10.0 * (1 + 1e-12), 0.0, 0, 0), (3e-15, 10.0, 0, 0)])
    assert verdict(tmp_path, make_fem(), tiny)["differences"] == []                 # rounding noise
    assert verdict(tmp_path, make_fem(), tiny, rtol=1e-14)["differences"] != []       # unless asked to be strict


def test_rtol_zero_means_exactly_equal_and_a_file_is_equal_to_itself(tmp_path, capsys):
    """With rtol = 0 the position search used to find nothing (its distance limit is exclusive): 810 'differences'."""
    a = write(tmp_path, "a.fem", make_fem())
    model = sm.parse_fem(a)
    result = sm.compare(model, model, rtol=0.0)
    assert result["differences"] == [] and result["formatting"] == []
    assert sm.main([str(a), str(a), "--rtol", "0"]) == 0
    assert "MODEL DIFFERENCES: 0" in capsys.readouterr().out
    # exactly equal means exactly equal: 1e-12 of a B-H value is a difference at rtol 0 and none at the default
    one, two = make_fem(materials=steel(1108.325569)), make_fem(materials=steel(1108.325569 * (1 + 1e-12)))
    assert verdict(tmp_path, one, two)["differences"] == []
    assert any("B-H point 2" in d for d in verdict(tmp_path, one, two, rtol=0.0)["differences"])


def test_positions_are_matched_within_a_floor_even_at_rtol_zero(tmp_path):
    """Positions are matched, not compared: nearer than 1e-12 units is one point, whatever --rtol says; 1e-9 is not."""
    near = make_fem(nodes=[(0.0, 0.0, 0, 0), (10.0 + 1e-13, 0.0, 0, 0), (0.0, 10.0, 0, 0)])
    far = make_fem(nodes=[(0.0, 0.0, 0, 0), (10.0 + 1e-9, 0.0, 0, 0), (0.0, 10.0, 0, 0)])
    assert verdict(tmp_path, make_fem(), near, rtol=0.0)["differences"] == []
    assert verdict(tmp_path, make_fem(), far, rtol=0.0)["differences"] != []


def test_there_is_no_absolute_floor_on_ordinary_numbers(tmp_path):
    """0 and 1e-20 are different numbers in the header, a property, a B-H table and a current (only positions have a floor)."""
    base = make_fem()
    header = verdict(tmp_path, base, base.replace("[Frequency]   =  0", "[Frequency]   =  1e-20"))
    assert any("[Frequency]" in d for d in header["differences"])
    current = verdict(tmp_path, make_fem(circuits=[("A", 0.0)]), make_fem(circuits=[("A", 1e-20)]))
    assert any("TotalAmps_re" in d for d in current["differences"])
    table = verdict(tmp_path, make_fem(materials=steel(0.0)), make_fem(materials=steel(1e-20)))
    assert any("B-H point 2" in d for d in table["differences"])
    field = verdict(tmp_path, make_fem(materials=[("Air", {"Mu_x": 1, "Sigma": 0}, [])] + steel(1.0)[1:]),
                    make_fem(materials=[("Air", {"Mu_x": 1, "Sigma": 1e-20}, [])] + steel(1.0)[1:]))
    assert any("Sigma" in d for d in field["differences"])


@pytest.mark.parametrize("bad", ["-1", "-1e-9", "nan"])
def test_a_negative_or_nan_rtol_is_refused(tmp_path, capsys, bad):
    a = write(tmp_path, "a.fem", make_fem())
    with pytest.raises(SystemExit) as stop:
        sm.main([str(a), str(a), f"--rtol={bad}"])
    assert stop.value.code == 2
    assert "--rtol must be zero or positive" in capsys.readouterr().err
    model = sm.parse_fem(a)
    with pytest.raises(ValueError, match="rtol must be zero or positive"):
        sm.compare(model, model, rtol=float(bad))


def test_two_nodes_of_old_on_one_node_of_new_are_reported(tmp_path):
    two_on_one = make_fem(nodes=[(0.0, 0.0, 0, 0), (0.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)], segments=[], arcs=[])
    result = verdict(tmp_path, two_on_one, make_fem(nodes=[(0.0, 0.0, 0, 0), (10.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)],
                                                     segments=[], arcs=[]))
    assert any("lie on the one" in line for line in result["differences"])


# ----------------------------------------------------------------------------- what is not a model, and what is not understood
def test_text_that_is_not_a_fem_file_is_refused_not_called_the_same(tmp_path, capsys):
    """Text against an empty .fem used to print an empty 'compared:' line, 'the same model' and exit 0."""
    notes = tmp_path / "notes.txt"
    notes.write_text("Dear reader,\nthis is not a FEMM file.\n[Section] = 3\n", encoding="utf-8")
    empty = tmp_path / "empty.fem"
    empty.write_bytes(b"")
    good = write(tmp_path, "good.fem", make_fem())
    for first, second, who in ((notes, empty, "OLD"), (good, empty, "NEW"), (empty, good, "OLD"), (notes, good, "OLD"),
                               (good, notes, "NEW")):
        assert sm.main([str(first), str(second)]) == 2
        captured = capsys.readouterr()
        assert captured.out == "", "nothing may be printed as a result"
        assert f"could not read {who}" in captured.err
        assert "no nodes, no block labels and no materials" in captured.err and "not a FEMM .fem file" in captured.err
    with pytest.raises(ValueError, match="no nodes, no block labels and no materials"):
        sm.parse_fem(empty)
    with pytest.raises(OSError):
        sm.parse_fem(tmp_path / "missing.fem")


def test_a_header_only_file_is_an_empty_model_and_is_refused(tmp_path):
    header_only = "\r\n".join(make_fem().split("\r\n")[:9]) + "\r\n[PointProps]   = 0\r\n[BlockProps]  = 0\r\n[NumPoints] = 0\r\n"
    path = write(tmp_path, "header.fem", header_only)
    with pytest.raises(ValueError, match="empty model"):
        sm.parse_fem(path)


def test_a_file_that_stops_or_breaks_in_the_middle_is_refused_and_says_where(tmp_path, capsys):
    text = make_fem()
    nodes_at = text.index("[NumPoints]")
    cases = {
        "ends inside a material": (text[:text.index("<BHPoints> = 3")], "ends inside block 2 of [BlockProps]"),
        "ends before a block": (text[:text.index("[CircuitProps]") + len("[CircuitProps]  = 1\r\n")],
                                "ends before block 1 of the 1 that [CircuitProps] announces"),
        "ends among the nodes": (text[:nodes_at + len("[NumPoints] = 3\r\n0.0\t0.0\t0\t0\r\n")],
                                 "[NumPoints] announces 3 lines and the file ends after 1"),
        "a count that is too big": (text.replace("[NumPoints] = 3", "[NumPoints] = 5"),
                                    "line 4 of [NumPoints] is not a line of at least 4 numbers: '[NumSegments] = 1'"),
        "text where a number belongs": (text.replace("10.0\t0.0\t0\t0", "10.0\tabc\t0\t0"), "could not convert string to float"),
        "a line with too few numbers": (text.replace("0.0\t10.0\t0\t0", "0.0\t10.0"), "line 3 of [NumPoints]"),
        "an edge to a node that is not there": (make_fem(segments=[(0, 7, -1, 1, 0, 0)]),
                                                "an edge names node 0 or 7, and the file has 3 nodes"),
    }
    good = write(tmp_path, "good.fem", make_fem())
    for name, (broken_text, why) in cases.items():
        broken = write(tmp_path, "broken.fem", broken_text)
        with pytest.raises(ValueError, match="is not a FEMM .fem file that this tool can read") as error:
            sm.parse_fem(broken)
        assert why in str(error.value), (name, str(error.value))
        assert sm.main([str(good), str(broken)]) == 2, name
        captured = capsys.readouterr()
        assert captured.out == "" and "could not read NEW" in captured.err and why in captured.err, name


def test_a_section_the_tool_does_not_know_is_a_warning_not_silence(tmp_path, capsys):
    odd = ("[NumFoo] = 2\r\n1\t2\r\n3\t4\r\n"
           "[ConductorProps] = 1\r\n  <BeginConductor>\r\n    <ConductorName> = \"C1\"\r\n  <EndConductor>\r\n")
    both = verdict(tmp_path, make_fem(extra=odd), make_fem(extra=odd))
    assert both["differences"] == []                      # the same unknown header lines on both sides
    assert len(both["warnings"]) == 4                     # two sections, in each of the two files
    assert both["warnings"][0].startswith("OLD: [NumFoo] is a section this tool does not know")
    assert {w.split(":")[0] for w in both["warnings"]} == {"OLD", "NEW"}
    assert not verdict(tmp_path, make_fem(), make_fem())["warnings"]
    # in one file only, the section also shows as a difference of the problem definition (its count line is not the same)
    one = verdict(tmp_path, make_fem(), make_fem(extra=odd))
    assert "problem definition: [NumFoo] is only in NEW" in one["differences"] and len(one["warnings"]) == 2
    # the command says it, and the exit status is still the verdict on what it could read
    old, new = write(tmp_path, "old.fem", make_fem(extra=odd)), write(tmp_path, "new.fem", make_fem(extra=odd))
    assert sm.main([str(old), str(new)]) == 0
    out = capsys.readouterr().out
    assert "WARNINGS (parts of a file this tool did not compare): 4" in out
    assert "[NumFoo]" in out and "[ConductorProps]" in out
    assert "RESULT: the same model (as far as the parts it understands go: see the WARNINGS)" in out
    assert sm.main([str(old), str(write(tmp_path, "plain.fem", make_fem()))]) == 1
    assert "RESULT: NOT the same model" in capsys.readouterr().out


def test_the_solution_part_of_an_ans_file_is_cut_off_before_anything_in_it_is_read(tmp_path):
    """The old code skipped it only because none of its lines happened to look like a section."""
    solution = ("[Solution]\r\n3\r\n0\t0\t1\t0\r\n"
                "[NumPoints] = 1\r\n99\t99\t0\t0\r\n"           # would replace the nodes, if it were read
                "[Precision] = 5\r\n"                            # would change the problem definition
                "[NumFoo] = 1\r\n")                              # would be a warning
    ans = write(tmp_path, "model.ans", make_fem() + solution)
    model = sm.parse_fem(ans)
    assert model["header"]["Precision"] == "1e-08" and len(model["nodes"]) == 3 and model["warnings"] == []
    assert verdict(tmp_path, make_fem(), make_fem() + solution)["differences"] == []


def test_command_line_exit_status_and_words(tmp_path, capsys):
    same_a, same_b = write(tmp_path, "a.fem", make_fem()), write(tmp_path, "b.fem", make_fem(precision="1e-008"))
    other = write(tmp_path, "c.fem", make_fem(circuits=[("A", 7.0)]))
    assert sm.main([str(same_a), str(same_b)]) == 0
    out = capsys.readouterr().out
    assert "RESULT: the same model" in out and "MODEL DIFFERENCES: 0" in out and "WARNINGS" not in out
    assert "NOT THE MODEL (formatting, notes, order, values FEMM computes): 1" in out
    assert sm.main([str(same_a), str(other)]) == 1
    out = capsys.readouterr().out
    assert "RESULT: NOT the same model" in out and "TotalAmps_re is 6 in OLD and 7 in NEW" in out
    assert sm.main([str(same_a), str(tmp_path / "missing.fem")]) == 2


# ----------------------------------------------------------------------------- what FEMM 5.0 really writes
@pytest.mark.femm
def test_two_models_that_femm_5_saved_differ_in_exactly_the_one_thing_changed():
    """learn/03's magnet, saved twice by the installed FEMM, with and without the magnet's coercivity."""
    OUT.mkdir(parents=True, exist_ok=True)
    with_magnet, without = OUT / "magnet_hc.fem", OUT / "magnet_no_hc.fem"
    femm_session.open_femm()
    try:
        for path, hc in ((with_magnet, 1_030_272), (without, 0)):
            tools.draw_magnet(hc)
            femm.mi_saveas(str(path))
            femm.mi_close()
    finally:
        femm_session.close_femm()
    a, b = sm.parse_fem(with_magnet), sm.parse_fem(without)
    assert (len(a["nodes"]), len(a["segments"]), len(a["arcs"])) == (8, 8, 0)
    assert [lb["material"] for lb in a["labels"]] == ["N45", "Air"]
    assert sm.compare(a, a)["differences"] == []
    result = sm.compare(a, b)
    assert result["differences"] == ['material "N45": H_c is 1030272 in OLD and 0 in NEW']
    assert result["formatting"] == [] and result["warnings"] == []


@pytest.mark.femm
def test_the_default_motor_as_femm_5_saves_it():
    """The same model written twice is the same file; the fine mesh differs from the coarse in [DoSmartMesh] and nothing else."""
    OUT.mkdir(parents=True, exist_ok=True)
    first, again, fine = OUT / "motor_coarse_1.fem", OUT / "motor_coarse_2.fem", OUT / "motor_fine.fem"
    p = mg.Params()
    regions = mg.build(p)
    femm_session.open_femm()
    try:
        for path, is_fine in ((first, False), (again, False), (fine, True)):
            mg.to_femm(p, regions, str(path), fine=is_fine)
            femm.mi_close()
    finally:
        femm_session.close_femm()
    a, b, c = sm.parse_fem(first), sm.parse_fem(again), sm.parse_fem(fine)
    assert (len(a["nodes"]), len(a["segments"]), len(a["arcs"]), len(a["labels"])) == (167, 136, 62, 40)
    twice = sm.compare(a, b)
    assert twice["differences"] == [] and twice["formatting"] == [] and twice["warnings"] == []
    assert sm.compare(a, c)["differences"] == ["problem definition: [DoSmartMesh] is 0 in OLD and 1 in NEW"]
