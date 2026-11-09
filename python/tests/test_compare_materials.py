"""Tests for validation/compare_materials.py: it must say "equal" only when every number is, and name what moved.

Run from P:   uv run pytest tests/test_compare_materials.py -q
"""
import json

import numpy as np
import pytest
from PIL import Image

import compare_materials as cmp_


def write_json(folder, name, data):
    folder.mkdir(exist_ok=True)
    (folder / name).write_text(json.dumps(data), encoding="utf-8")


def run(tmp_path, old, new, *extra):
    write_json(tmp_path / "old", "x.json", old)
    write_json(tmp_path / "new", "x.json", new)
    return cmp_.main([str(tmp_path / "old"), str(tmp_path / "new"), *extra])


DATA = {"a": 1.5, "b": {"c": [1, 2, 3], "d": "text"}, "e": None, "f": True}


def test_equal_files(tmp_path, capsys):
    assert run(tmp_path, DATA, json.loads(json.dumps(DATA))) == 0
    out = capsys.readouterr().out
    assert "x.json: equal (7 values)" in out and "0 differences" in out       # a, c's three, d, e and f


def test_a_moved_number_is_named_with_its_size(tmp_path, capsys):
    new = json.loads(json.dumps(DATA))
    new["b"]["c"][1] = 2.5
    assert run(tmp_path, DATA, new) == 1
    out = capsys.readouterr().out
    assert "/b/c[1] is 2 in OLD and 2.5 in NEW  (+25 %)" in out


@pytest.mark.parametrize("change", [lambda d: d.pop("a"), lambda d: d.update(z=1), lambda d: d["b"]["c"].append(4),
                                    lambda d: d["b"].update(d="other"), lambda d: d.update(f=False),
                                    lambda d: d.update(e=0), lambda d: d.update(a="1.5")])
def test_missing_extra_and_changed_things_are_differences(tmp_path, change):
    new = json.loads(json.dumps(DATA))
    change(new)
    assert run(tmp_path, DATA, new) == 1


def test_tolerance(tmp_path):
    assert run(tmp_path, {"a": 1.0}, {"a": 1.0 + 1e-10}) == 1
    assert run(tmp_path, {"a": 1.0}, {"a": 1.0 + 1e-10}, "--rtol", "1e-9") == 0


def test_a_file_missing_from_new_is_a_difference(tmp_path):
    write_json(tmp_path / "old", "only_old.json", {"a": 1})
    (tmp_path / "new").mkdir()
    assert cmp_.main([str(tmp_path / "old"), str(tmp_path / "new")]) == 1


def test_csv_compares_numbers_and_treats_two_blanks_as_equal(tmp_path):
    for folder, last in (("old", "3"), ("new", "3")):
        (tmp_path / folder).mkdir(exist_ok=True)
        (tmp_path / folder / "t.csv").write_text(f"B_T,mu_r\n0,nan\n1,{last}\n", encoding="utf-8")
    assert cmp_.main([str(tmp_path / "old"), str(tmp_path / "new")]) == 0
    (tmp_path / "new" / "t.csv").write_text("B_T,mu_r\n0,nan\n1,3.5\n", encoding="utf-8")
    assert cmp_.main([str(tmp_path / "old"), str(tmp_path / "new")]) == 1


def test_figures_are_compared_pixel_by_pixel_only_when_asked(tmp_path):
    for folder in ("old", "new"):
        write_json(tmp_path / folder, "x.json", {"a": 1})                  # something to compare besides the figure
    pixels = np.zeros((4, 4, 3), dtype=np.uint8)
    Image.fromarray(pixels).save(tmp_path / "old" / "f.png")
    changed = pixels.copy()
    changed[1, 1] = (255, 0, 0)
    Image.fromarray(changed).save(tmp_path / "new" / "f.png")
    assert cmp_.main([str(tmp_path / "old"), str(tmp_path / "new")]) == 0                 # not asked: not looked at
    assert cmp_.main([str(tmp_path / "old"), str(tmp_path / "new"), "--figures"]) == 1
    Image.fromarray(pixels).save(tmp_path / "new" / "f.png")
    assert cmp_.main([str(tmp_path / "old"), str(tmp_path / "new"), "--figures"]) == 0


def two_folders(tmp_path):
    for folder in ("old", "new"):
        (tmp_path / folder).mkdir(exist_ok=True)
    return tmp_path / "old", tmp_path / "new"


def test_two_folders_of_pictures_can_be_compared_on_their_own_when_asked(tmp_path, capsys):
    """The report compares FEMM 5.0's 16 figures with the repository's docs figures folder, which holds pictures and nothing else."""
    old, new = two_folders(tmp_path)
    pixels = np.zeros((4, 4, 3), dtype=np.uint8)
    for name in ("a.png", "b.png"):
        Image.fromarray(pixels).save(old / name)
        Image.fromarray(pixels).save(new / name)
    moved = pixels.copy()
    moved[0, 0] = (9, 9, 9)
    Image.fromarray(moved).save(new / "b.png")
    assert cmp_.main([str(old), str(new)]) == 2                                      # not asked for: nothing to compare
    assert "no JSON or CSV files in OLD" in capsys.readouterr().err
    assert cmp_.main([str(old), str(new), "--figures"]) == 1
    out = capsys.readouterr().out
    assert "a.png: equal (16 values)" in out and "1 of 16 pixels differ" in out and "2 files, 32 values compared, 1 differences" in out
    Image.fromarray(pixels).save(new / "b.png")
    assert cmp_.main([str(old), str(new), "--figures"]) == 0
    empty_old, empty_new = tmp_path / "e1", tmp_path / "e2"
    empty_old.mkdir()
    empty_new.mkdir()
    assert cmp_.main([str(empty_old), str(empty_new), "--figures"]) == 2
    assert "no JSON, CSV or PNG files in OLD" in capsys.readouterr().err


def test_a_csv_row_with_another_number_of_columns_is_a_difference(tmp_path, capsys):
    """Rows and columns used to be zipped, so a short row passed with its missing cells unseen."""
    old, new = two_folders(tmp_path)
    (old / "t.csv").write_text("a,b,c\n1,2,3\n4,5,6\n", encoding="utf-8")
    for body, expect in (("1,2,3\n4,5\n", "row 2 has 3 columns in OLD and 2 in NEW"),
                         ("1,2,3\n4,5,6,7\n", "row 2 has 3 columns in OLD and 4 in NEW"),
                         ("1,2,3\n", "2 in OLD and 1 in NEW")):
        (new / "t.csv").write_text("a,b,c\n" + body, encoding="utf-8")
        assert cmp_.main([str(old), str(new)]) == 1, body
        assert expect in capsys.readouterr().out
    (new / "t.csv").write_text("a,b\n1,2\n4,5\n", encoding="utf-8")                        # a column fewer all through
    assert cmp_.main([str(old), str(new)]) == 1


def test_a_csv_with_text_cells_is_compared_and_does_not_crash(tmp_path):
    old, new = two_folders(tmp_path)
    (old / "t.csv").write_text("name,v\nsteel,1\n", encoding="utf-8")
    (new / "t.csv").write_text("name,v\nsteel,1\n", encoding="utf-8")
    assert cmp_.main([str(old), str(new)]) == 0
    (new / "t.csv").write_text("name,v\niron,1\n", encoding="utf-8")
    assert cmp_.main([str(old), str(new)]) == 1


def test_a_file_that_only_new_has_is_reported(tmp_path, capsys):
    """Only OLD's files used to be looked at, so an extra file in NEW went unnoticed."""
    old, new = two_folders(tmp_path)
    write_json(old, "x.json", {"a": 1})
    write_json(new, "x.json", {"a": 1})
    write_json(new, "y.json", {"a": 1})
    (new / "z.csv").write_text("a\n1\n", encoding="utf-8")
    Image.fromarray(np.zeros((2, 2, 3), dtype=np.uint8)).save(new / "p.png")
    assert cmp_.main([str(old), str(new)]) == 1
    out = capsys.readouterr().out
    assert "y.json: only in NEW" in out and "z.csv: only in NEW" in out
    assert "p.png" not in out, "figures are looked at only when asked"
    assert "3 files, 1 values compared, 2 differences" in out
    assert cmp_.main([str(old), str(new), "--figures"]) == 1
    assert "p.png: only in NEW" in capsys.readouterr().out


def test_nothing_to_compare_is_an_error_not_a_pass(tmp_path, capsys):
    old, new = two_folders(tmp_path)                                                       # two empty folders
    assert cmp_.main([str(old), str(new)]) == 2
    assert "no JSON or CSV files in OLD" in capsys.readouterr().err
    assert cmp_.main([str(tmp_path / "nope"), str(new)]) == 2
    assert "OLD is not a folder" in capsys.readouterr().err
    assert cmp_.main([str(old), str(tmp_path / "nope")]) == 2
    assert "NEW is not a folder" in capsys.readouterr().err


def test_a_negative_rtol_is_refused(tmp_path, capsys):
    old, new = two_folders(tmp_path)
    with pytest.raises(SystemExit) as stop:
        cmp_.main([str(old), str(new), "--rtol=-1"])
    assert stop.value.code == 2 and "--rtol must be zero or positive" in capsys.readouterr().err
