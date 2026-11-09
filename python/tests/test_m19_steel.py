"""Tests for the "M-19 Steel" that motorgeo defines itself (motorgeo.M19_STEEL_FEMM42).

Why the model defines it: FEMM 4.2's material library had a block called "M-19 Steel" (AK Steel's 29-gauge DC curve,
47 points, lamination fill 0.98). FEMM 5.0's library has no such block; its "M-19 14 mil" is a different curve
(61 points, fill 0.955). To get the same motor on both FEMMs, to_femm() no longer asks a library for the steel: it
defines it from a copy of FEMM 4.2's block, which is the constant these tests check.

  test_constant_is_femm42s_block        the copy equals FEMM 4.2's own library file, field by field (skipped when
                                        FEMM 4.2 is no longer installed)
  test_constant_has_the_documented_values  the numbers worth remembering, so the copy stays checked without FEMM 4.2
  test_saved_fem_defines_m19_steel      FEMM 5.0 really takes the definition: the .fem that to_femm() saves holds
                                        the same block (needs FEMM 5.0)

Files go to P\\out\\tests\\ (an ordinary folder of P, which git ignores). Run from P:   uv run pytest tests/test_m19_steel.py -q
"""
import re
from pathlib import Path

import pytest

import femm_session
import motorgeo as mg

OUT = Path(__file__).resolve().parent.parent / "out" / "tests"
LIBRARY_42 = Path(r"D:\Apps\femm42\bin\matlib.dat")

# Every field of a FEMM material block, in the order FEMM writes them, with its meaning and unit.
FIELDS = ("BlockName", "Mu_x", "Mu_y", "H_c", "H_cAngle", "J_re", "J_im", "Sigma", "d_lam", "Phi_h", "Phi_hx",
          "Phi_hy", "LamType", "LamFill", "NStrands", "WireD", "BHPoints")


def read_block(path, name):
    """The material block `name` of a FEMM library or .fem file, as (fields, bh_points).

    Both file kinds write a block the same way (indentation aside):

        <BeginBlock>
        <BlockName> = "M-19 Steel"
        <Mu_x> = 4416
        ...
        <BHPoints> = 47
        0<TAB>0                      <- one "B<TAB>H" line per point
        0.05<TAB>15.12
        <EndBlock>

    Numbers come back as floats (exactly the double that the 17 digits in the file stand for), the name as text.
    """
    text = Path(path).read_text(encoding="latin-1")
    for body in re.findall(r"<BeginBlock>(.*?)<EndBlock>", text, flags=re.S):
        fields, bh = {}, []
        for line in body.splitlines():
            line = line.strip()
            tagged = re.match(r"<(\w+)>\s*=\s*(.*)$", line)
            if tagged:
                key, value = tagged.group(1), tagged.group(2).strip()
                fields[key] = value.strip('"') if value.startswith('"') else float(value)
            elif line:
                b, h = line.split()
                bh.append((float(b), float(h)))
        if fields.get("BlockName") == name:
            return fields, bh
    raise AssertionError(f'no material block called "{name}" in {path}')


def assert_block_is_the_constant(fields, bh):
    """The block read from a file equals mg.M19_STEEL_FEMM42: every field, and every (B, H) pair, exactly."""
    const = mg.M19_STEEL_FEMM42
    assert set(const) == set(fields) | {"BH"}, "the constant must hold exactly the block's fields plus BH"
    for key, value in fields.items():
        assert const[key] == value, f"{key}: the constant has {const[key]!r}, the file has {value!r}"
    assert len(bh) == fields["BHPoints"] == len(const["BH"]) == 47
    for k, ((b, h), (cb, ch)) in enumerate(zip(bh, const["BH"])):
        assert (b, h) == (cb, ch), f"B-H point {k}: the file has ({b!r}, {h!r}), the constant ({cb!r}, {ch!r})"


@pytest.mark.skipif(not LIBRARY_42.is_file(), reason="FEMM 4.2 is not installed (D:\\Apps\\femm42\\bin\\matlib.dat)")
def test_constant_is_femm42s_block():
    fields, bh = read_block(LIBRARY_42, "M-19 Steel")
    assert_block_is_the_constant(fields, bh)


def test_constant_has_the_documented_values():
    m = mg.M19_STEEL_FEMM42
    assert m["BlockName"] == "M-19 Steel"
    assert (m["Mu_x"], m["Mu_y"]) == (4416, 4416)
    assert m["Sigma"] == 1.9                        # MS/m
    assert m["d_lam"] == 0.635                      # mm
    assert m["LamFill"] == 0.98
    assert m["LamType"] == 0                        # laminated in plane
    assert m["BHPoints"] == len(m["BH"]) == 47
    b, h = zip(*m["BH"])
    assert (b[0], h[0]) == (0.0, 0.0)               # starts at the origin
    assert (b[-1], h[-1]) == (2.3, 234024.751347)   # and ends deep in saturation
    assert all(b2 > b1 and h2 > h1 for b1, b2, h1, h2 in zip(b, b[1:], h, h[1:])), "B and H must both rise"


@pytest.mark.femm
def test_saved_fem_defines_m19_steel():
    p = mg.Params()
    fem = OUT / "m19" / "ipm.fem"
    fem.parent.mkdir(parents=True, exist_ok=True)
    fem.unlink(missing_ok=True)
    femm_session.open_femm()
    try:
        mg.to_femm(p, mg.build(p), str(fem))           # the coarse mesh (fine=False)
    finally:
        femm_session.close_femm()
    fields, bh = read_block(fem, "M-19 Steel")
    assert_block_is_the_constant(fields, bh)
