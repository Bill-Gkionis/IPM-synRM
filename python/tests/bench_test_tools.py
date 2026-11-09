"""What the bench's test files share. There is no test in here (the name does not start with test_).

Four files test the bench: tests/test_bench.py (the shared part: bench/__init__.py and run_bench.py, and the scans that read every
module of bench/) and one file per case, tests/test_bench_b1.py, test_bench_b2.py and test_bench_b3.py. Tasks 12 and 13 add their own files
for B4-B9 in the same way. What more than one of them needs is here:

  * `use_config` writes a config.toml of the test's own, so that nothing is written under the real out folder.
  * `the_test_run_is_not_lowered` is an autouse fixture a file imports when one of its tests calls run_bench.main(): that call lowers its
    own process to below normal (R10), and a test must not lower the whole pytest run. It puts `lower_priority_is_off` in the place of
    femm_session.lower_priority; a test that asserts the place is held by THAT function (test_bench.py has one) cannot be fooled by what ran
    before it in a whole-file run, as a test that reads the priority class from the operating system can.
  * `tiny_mesh` writes a two-rectangle mesh with Gmsh, small enough to know every area (B2 and B3 read it back).
  * `P` is the folder above tests/.
"""
from pathlib import Path

import pytest

import paths

P = Path(__file__).resolve().parent.parent


def lower_priority_is_off():
    """What femm_session.lower_priority is in a test of the bench: a call that lowers nothing."""


@pytest.fixture(autouse=True)
def the_test_run_is_not_lowered(monkeypatch):
    """run_bench.main() lowers its own process to below normal (R10). A test that calls it must not lower the whole pytest run, which would slow the
    timing-sensitive tests of the other files that follow: the call is replaced by `lower_priority_is_off` in every test of a file that imports this fixture."""
    import femm_session
    monkeypatch.setattr(femm_session, "lower_priority", lower_priority_is_off)


def use_config(monkeypatch, tmp_path) -> Path:
    """A config.toml of the test's own, so that nothing is written under the real out\\. Returns its out folder."""
    out = tmp_path / "out"
    config = tmp_path / "config.toml"
    config.write_text(f"[tools]\n[paths]\nout = '{out}'\n[femm]\nworkers = 8\npriority = 'below_normal'\n", encoding="utf-8")
    monkeypatch.setattr(paths, "CONFIG", config)
    return out


def tiny_mesh(path):
    """Two rectangles of the unit square (0.4 and 0.6 wide), one physical group each, and the bottom edge as a curve group."""
    import gmsh

    import checks
    with checks.intact_path():
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("tiny")
            g = gmsh.model.geo
            pts = [g.addPoint(x, y, 0, 0.1) for x, y in ((0, 0), (0.4, 0), (1, 0), (1, 1), (0.4, 1), (0, 1))]
            b1, b2 = g.addLine(pts[0], pts[1]), g.addLine(pts[1], pts[2])
            r, top2, mid, top1, left = (g.addLine(pts[2], pts[3]), g.addLine(pts[3], pts[4]), g.addLine(pts[1], pts[4]),
                                        g.addLine(pts[4], pts[5]), g.addLine(pts[5], pts[0]))
            left_surface = g.addPlaneSurface([g.addCurveLoop([b1, mid, top1, left])])
            right_surface = g.addPlaneSurface([g.addCurveLoop([b2, r, top2, -mid])])
            g.synchronize()
            gmsh.model.addPhysicalGroup(2, [left_surface], name="left")
            gmsh.model.addPhysicalGroup(2, [right_surface], name="right")
            gmsh.model.addPhysicalGroup(1, [b1, b2], name="bottom")
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.model.mesh.generate(2)
            gmsh.write(str(path))
        finally:
            gmsh.finalize()
