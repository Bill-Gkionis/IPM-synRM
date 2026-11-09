"""Tests that femfile.py still reads what FEMM 5.0 writes.

femfile reads FEMM's text files without FEMM: the geometry in a .fem and the mesh + vector potential in a .ans.
It was written against FEMM 4.2's files. These tests check it against FEMM 5.0's:

  test_b_per_triangle_matches_femm   B computed from the .ans by femfile.element_b equals what FEMM itself reports
                                     at the centre of the same triangles (smoothing switched off)
  test_fem_reads_back                the default motor, as FEMM 5.0 saved it, has as many nodes, segments, arcs and
                                     labels as motorgeo drew

Both start a real FEMM (see femm_session.py). Files go to P\\out\\tests\\. Run from P:
    uv run pytest tests/test_femfile_femm5.py -q
"""
from pathlib import Path

import numpy as np
import pytest

import femm
import femm_session
import femfile
import femm_test_tools as tools
import motorgeo as mg

OUT = Path(__file__).resolve().parent.parent / "out" / "tests"


@pytest.mark.femm
def test_b_per_triangle_matches_femm():
    path = OUT / "b_per_triangle" / "magnet.fem"
    femm_session.open_femm()
    try:
        tools.solve_magnet(path, 1_030_272)             # learn/03's magnet; the .ans stays next to the .fem
        femm.mo_smooth("off")                           # B constant per triangle, as femfile computes it
        ans = femfile.read_ans(path.with_suffix(".ans"))
        bx, by = femfile.element_b(ans["xy"], ans["A"], ans["tri"])
        centres = ans["xy"][ans["tri"]].mean(axis=1)    # one (x, y) per triangle, in mm like the model
        # 20 triangles spread over the range of field strength, weakest to strongest, ignoring the ones with less than
        # 1 mT (two near-zeros agreeing proves little). A centre lies at least a third of the triangle's height from
        # every edge, so FEMM finds the same triangle there.
        strength = np.hypot(bx, by)
        usable = np.flatnonzero(strength > 1e-3)
        assert len(usable) > 1000, "too few triangles to pick from: the .ans was not read properly"
        ranked = usable[np.argsort(strength[usable])]
        picks = ranked[np.linspace(0, len(ranked) - 1, 20).astype(int)]
        assert strength[picks[-1]] == strength.max(), "the strongest triangle of the mesh must be among the 20"
        worst = 0.0
        for k in picks:
            femm_bx, femm_by = femm.mo_getb(float(centres[k, 0]), float(centres[k, 1]))
            worst = max(worst, abs(femm_bx - bx[k]), abs(femm_by - by[k]))
        assert worst < 1e-6, f"femfile and FEMM disagree by {worst:.3g} T on a triangle's B"
    finally:
        femm_session.close_femm()


@pytest.mark.femm
def test_fem_reads_back():
    p = mg.Params()
    regions = mg.build(p)
    fem = OUT / "reads_back" / "ipm.fem"
    fem.parent.mkdir(parents=True, exist_ok=True)
    fem.unlink(missing_ok=True)
    femm_session.open_femm()
    try:
        mg.to_femm(p, regions, str(fem))
    finally:
        femm_session.close_femm()
    # what motorgeo wrote: every corner once, every straight edge, every arc, one label per region
    edges = mg.unique_edges(regions)
    wrote = {"nodes": len({mg.key(q) for e, _ in edges for q in (e.a, e.b)}),
             "segments": sum(isinstance(e, mg.Line) for e, _ in edges),
             "arcs": sum(isinstance(e, mg.Arc) for e, _ in edges),
             "labels": len(regions)}
    read = femfile.read_fem(fem)      # the one "<No Mesh>" region is saved as a hole; femfile counts it among the labels
    got = {"nodes": len(read["points"]), "segments": len(read["segments"]),
           "arcs": len(read["arcs"]), "labels": len(read["labels"])}
    assert got == wrote
    assert read["materials"] == ["Air", "M-19 Steel", "N45", "Coil"]
    assert read["circuits"] == ["A", "B", "C"]
