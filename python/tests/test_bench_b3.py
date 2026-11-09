"""Tests for B3 (bench/b3_elmer.py): Gmsh -> Elmer. Everything that needs no ElmerSolver; the mesh is read and the input is written.

Run from P:   uv run pytest tests/test_bench_b3.py -q      (a few seconds; no FEMM, no Elmer, no solver slot is taken)

In groups:
  1. The solver input (every body, material, boundary and keyword), the lamination fill, the two programs' printouts, the result file, the flux
     linkage, the Park transform, the gap profile, the comparison of two waveforms, the matching of nodes before Elmer is asked.
  2. The verdict: the laminated run carries it and the plain run is reported beside it (gating or not), the three node counts, the limits that
     are the bench's own and the brief's two 3 %, the headline never without the finer-mesh trend.
"""
import math

import numpy as np
import pytest

import motorgeo as mg
from bench_test_tools import P, tiny_mesh
from conftest import process_path


# ----------------------------------------------------------------------------------------------- 1. Elmer, the input, the logs and the numbers
GRID_LOG = """Elmergrid creating and manipulating meshes:
-------------------------------------------
Scaling mesh with vector [0.001 0.001 0.001]
Reordering 8819 knots and 16914 elements in 3-dimensions.
Merging nodes close (1e-09) to one another.
8819 original nodes merged to 8550 new nodes.
Merging of nodes is complete.
Calculate Mesh Pieces in mesh:
There is a single piece in the mesh, so the mesh is conforming.

Elmergrid saving data with method 2:
-------------------------------------
Saving mesh in ElmerSolver format to directory case.
Saving 8550 coordinates to mesh.nodes.
Saving 16914 element topologies to mesh.elements.
Saving boundary elements to mesh.boundary.
"""

SOLVER_LOG = """MAIN: Version: 26.2-devel (Rev: 11c118114, Compiled: 2026-09-17)
MAIN:  Running one task without MPI parallelization.
MAIN:  Running in parallel with 4 threads per task.
MAIN:  MUMPS library linked in.
MagnetoDynamics2D: Solving equation for magnetic vector potential
ComputeChange: NS (ITER=1) (NRM,RELC): ( 0.14445670E-02  2.0000000     ) :: mgdyn2d
ComputeChange: NS (ITER=2) (NRM,RELC): ( 0.20773421E-02 0.35933642     ) :: mgdyn2d
ComputeChange: NS (ITER=3) (NRM,RELC): ( 0.22986646E-02 0.10115272     ) :: mgdyn2d
ComputeChange: NS (ITER=11) (NRM,RELC): ( 0.23626453E-02 0.33519879E-06 ) :: mgdyn2d
ComputeChange: SS (ITER=1) (NRM,RELC): ( 0.23626453E-02  2.0000000     ) :: mgdyn2d
MagnetoDynamicsCalcFields: Computing postprocessed fields
WARNING:: GetPermittivity: Permittivity not defined in material, defaulting to that of vacuum
ComputeChange: NS (ITER=1) (NRM,RELC): ( 0.48306357      2.0000000     ) :: mgdynpost
ComputeChange: NS (ITER=2) (NRM,RELC): ( 0.48335514     0.60340967E-03 ) :: mgdynpost
MagnetoDynamicsCalcFields: ElectroMagnetic Field Energy:    3.708959E+01
MAIN: *** Elmer Solver: ALL DONE ***
SOLVER TOTAL TIME(CPU,REAL):         0.98        1.40
"""


def real_names(tmp_path):
    """The physical groups of the real path-1 mesh (35 regions, 5 boundary curves) and the regions they came from."""
    from bench import b2_gmsh, b3_elmer
    b2_gmsh.make_path1_mesh(tmp_path / "path1.msh")
    return b3_elmer.physical_names(tmp_path / "path1.msh"), mg.build(mg.Params(barrier="none"), gap="sliding")


def test_the_physical_groups_of_a_mesh_are_read_as_name_to_dimension_and_tag(tmp_path):
    from bench import b3_elmer
    tiny_mesh(tmp_path / "tiny.msh")
    names = b3_elmer.physical_names(tmp_path / "tiny.msh")
    assert set(names) == {"left", "right", "bottom"}
    assert names["left"][0] == 2 and names["right"][0] == 2 and names["bottom"][0] == 1
    assert len({tag for dim, tag in names.values() if dim == 2}) == 2


def test_the_laminated_steel_is_the_b_h_curve_with_the_fill_factor_in_the_way_femm_applies_it():
    """FEMM's M-19 Steel is the 47 points of its B-H curve AND a lamination fill factor of 0.98. FEMM 5.0 applies the fill as B_eff(H) = f B(H) + (1 - f) mu0 H
    (found on 2026-11-06: FEMM with fill 1.0 and the curve changed by that formula gives, to every printed digit, FEMM with fill 0.98), so an Elmer model of the
    same steel needs the same curve. With fill 1 the curve is the 47 points exactly."""
    from bench import b3_elmer
    plain = b3_elmer.steel_table(1.0)
    assert plain == [(float(b), float(h)) for b, h in mg.M19_STEEL_FEMM42["BH"]] and len(plain) == 47 == mg.M19_STEEL_FEMM42["BHPoints"]
    assert mg.M19_STEEL_FEMM42["LamFill"] == 0.98
    laminated = b3_elmer.steel_table(mg.M19_STEEL_FEMM42["LamFill"])
    assert [h for _, h in laminated] == [h for _, h in plain] and laminated[0] == (0.0, 0.0)
    mu0 = 4e-7 * math.pi
    row = next(r for r in laminated if r[1] == 1108.325569)                  # the point at B = 1.5 T
    assert row[0] == pytest.approx(0.98 * 1.5 + 0.02 * mu0 * 1108.325569, rel=1e-12)
    assert all(a[0] < b[0] for a, b in zip(laminated, laminated[1:]))        # still a curve Elmer can interpolate: B grows with every point
    assert all(l[0] < p[0] for l, p in zip(laminated[1:], plain[1:]))        # and every point is below the unlaminated one (the steel is 98 % of the stack)


def test_the_solver_input_names_every_body_and_boundary_by_the_tags_of_the_mesh_and_carries_the_physics(tmp_path):
    from bench import b3_elmer
    names, regions = real_names(tmp_path)
    text = b3_elmer.sif_text(names, regions, fill=0.98)
    assert text.isascii() and "\t" not in text and "\r" not in text
    assert text.count("\nBody ") == 35 and text.count("\nMaterial ") == 2 + 4
    assert text.count("Relative Permittivity = 1.0") == 6                      # every material has one: ElmerSolver's field solver asks (and warns if not)
    for r in regions:                                                          # every region is a body that targets its own physical group
        assert f'Name = "{r.name}"' in text and f"Target Bodies(1) = {names[r.name][1]}" in text
    assert "Body Force" not in text and "Current Density" not in text           # no load: no current anywhere
    assert text.count("Real Cubic") == 1 and text.count('H-B Curve = Variable "dummy"') == 1
    table = text.split("Real Cubic")[1].split("End")[0].strip().splitlines()
    assert len(table) == 47 and all(len(line.split()) == 2 for line in table)
    # N45 as FEMM has it: mu_r 1.05 and the coercivity as the magnetisation, along each magnet's own direction
    magnets = [r for r in regions if r.name.startswith("magnet_")]
    assert len(magnets) == 4 and text.count("Relative Permeability = 1.05") == 4
    for r in magnets:
        assert f"Magnetization 1 = Real {1030272.0 * math.cos(math.radians(r.magdir))!r}" in text
        assert f"Magnetization 2 = Real {1030272.0 * math.sin(math.radians(r.magdir))!r}" in text
    # the boundaries: A = 0, and the periodic pair with its node-to-node (conforming) condition
    assert f"Target Boundaries(1) = {names['A0'][1]}\n  Az = Real 0.0" in text
    assert "Conforming BC = Integer 3" in text and "Radial Projector = Logical True" in text and "Anti Radial" not in text   # two poles per sector: periodic
    assert f"Target Boundaries(1) = {names['edge_0'][1]}" in text and f"Target Boundaries(1) = {names['edge_sector'][1]}" in text
    assert "Apply Conforming BCs = Logical True" in text
    # the solver: MagnetoDynamics2D, its fields, the nonlinear loop, the result file
    for line in ('Procedure = "MagnetoDynamics2D" "MagnetoDynamics2D"', 'Procedure = "MagnetoDynamics" "MagnetoDynamicsCalcFields"',
                 'Procedure = "ResultOutputSolve" "ResultOutputSolver"', "Vtu Format = Logical True", "Ascii Output = Logical True", "Save Geometry Ids = Logical True",
                 "Nonlinear System Max Iterations = 60", "Nonlinear System Convergence Tolerance = 1e-06", "CHECK KEYWORDS Warn", 'Mesh DB "." "case"'):
        assert line in text, line


def test_the_plain_curve_variant_differs_from_the_laminated_one_in_the_table_only(tmp_path):
    from bench import b3_elmer
    names, regions = real_names(tmp_path)
    laminated, plain = b3_elmer.sif_text(names, regions, fill=0.98), b3_elmer.sif_text(names, regions, fill=1.0)
    assert laminated != plain
    cut = lambda text: (text.split("Real Cubic")[0], text.split("Real Cubic")[1].split("End", 1)[1])
    assert cut(laminated) == cut(plain)
    first = plain.split("Real Cubic")[1].split("End")[0].strip().splitlines()[1].split()
    assert (float(first[0]), float(first[1])) == (0.05, 15.120714)


def test_the_grid_converter_s_account_of_the_merge_is_read_from_its_log():
    from bench import b3_elmer
    grid = b3_elmer.parse_grid_log(GRID_LOG)
    assert grid == {"nodes_read": 8819, "nodes_merged_to": 8550, "nodes_merged": 269, "conforming": True, "pieces": 1, "coordinates_saved": 8550,
                    "elements_saved": 16914}
    broken = b3_elmer.parse_grid_log(GRID_LOG.replace("There is a single piece in the mesh, so the mesh is conforming.", "There are 2 pieces in the mesh."))
    assert broken["conforming"] is False and broken["pieces"] == 2
    assert b3_elmer.parse_grid_log("nothing useful")["nodes_read"] is None


def test_the_solver_s_nonlinear_iterations_and_convergence_are_read_from_its_log():
    from bench import b3_elmer
    s = b3_elmer.parse_solver_log(SOLVER_LOG)
    assert s["iterations"] == 4 and s["last_change"] == pytest.approx(3.3519879e-07) and s["changes"][:2] == [2.0, 0.35933642]       # only the MagnetoDynamics2D lines
    assert s["finished"] is True and s["threads"] == 4 and s["version"] == "26.2-devel" and s["revision"] == "11c118114"
    assert s["energy_J"] == pytest.approx(37.08959) and s["solver_seconds_real"] == pytest.approx(1.40)
    assert s["warnings"] == ["GetPermittivity: Permittivity not defined in material, defaulting to that of vacuum"]
    unfinished = b3_elmer.parse_solver_log(SOLVER_LOG.replace("MAIN: *** Elmer Solver: ALL DONE ***\n", "ERROR:: ConformingNodePerm: This is not conforming!\n"))
    assert unfinished["finished"] is False and unfinished["errors"] == ["ConformingNodePerm: This is not conforming!"]


VTU = """<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
  <UnstructuredGrid>
    <Piece NumberOfPoints="4" NumberOfCells="3">
      <PointData>
        <DataArray type="Float64" Name="az" NumberOfComponents="1" format="ascii">
 0.0 1.0E-005 2.0E-005 3.0E-005
        </DataArray>
        <DataArray type="Float64" Name="magnetic flux density" NumberOfComponents="3" format="ascii">
 0.0 0.1 0.0 0.2 0.0 0.0 0.0 0.3 0.0 0.4 0.4 0.0
        </DataArray>
      </PointData>
      <CellData>
        <DataArray type="Int32" Name="GeometryIds" format="ascii">
 7 8 137
        </DataArray>
      </CellData>
      <Points>
        <DataArray type="Float64" NumberOfComponents="3" format="ascii">
 0.0 0.0 0.0 0.001 0.0 0.0 0.0 0.001 0.0 0.001 0.001 0.0
        </DataArray>
      </Points>
      <Cells>
        <DataArray type="Int32" Name="connectivity" format="ascii">
 0 1 2 1 3 2 0 1
        </DataArray>
        <DataArray type="Int32" Name="offsets" format="ascii">
 3 6 8
        </DataArray>
        <DataArray type="Int32" Name="types" FORMAT="ascii">
 5 5 3
        </DataArray>
      </Cells>
    </Piece>
  </UnstructuredGrid>
</VTKFile>
"""


def test_the_result_file_is_read_into_nodes_triangles_fields_and_body_ids(tmp_path):
    from bench import b3_elmer
    (tmp_path / "case_t0001.vtu").write_text(VTU, encoding="utf-8")
    v = b3_elmer.read_vtu(tmp_path / "case_t0001.vtu")
    assert v["points_mm"].shape == (4, 2) and v["points_mm"][3] == pytest.approx([1.0, 1.0])                  # metres in the file, millimetres here
    assert v["az"] == pytest.approx([0.0, 1e-5, 2e-5, 3e-5]) and v["b"].shape == (4, 2) and v["b"][3] == pytest.approx([0.4, 0.4])
    assert v["triangles"].tolist() == [[0, 1, 2], [1, 3, 2]] and v["body_ids"].tolist() == [7, 8]             # the line element (type 3) is not a triangle
    assert v["lines"].tolist() == [[0, 1]] and v["line_ids"].tolist() == [137]


def test_the_flux_linkage_of_a_phase_is_its_turns_times_the_mean_potential_of_each_coil_side_times_the_stack(tmp_path):
    """The way FEMM computes a circuit's flux linkage: for every block of the circuit, turns x (the integral of A over the block / its area) x the depth."""
    from bench import b3_elmer

    class R:                                                                 # the two coil sides of one phase, and a side of another
        def __init__(self, name, circuit, turns, material="Coil"):
            self.name, self.circuit, self.turns, self.material = name, circuit, turns, material
    regions = [R("coil_a_plus", "A", 10), R("coil_a_minus", "A", -10), R("coil_b", "B", 10), R("rotor_steel", "", 0, "M-19 Steel")]
    names = {"coil_a_plus": (2, 1), "coil_a_minus": (2, 2), "coil_b": (2, 3), "rotor_steel": (2, 4)}
    # two triangles per coil side, equal areas; A = 2e-5 + 1e-5 on the + side (mean 2e-5), 1e-5 on the - side, 5e-6 on B
    points = np.array([[0, 0], [1, 0], [0, 1], [1, 1], [2, 0], [2, 1], [3, 0], [3, 1], [4, 0], [4, 1]], float)
    triangles = np.array([[0, 1, 2], [1, 3, 2], [1, 4, 3], [4, 5, 3], [4, 6, 5], [6, 7, 5], [6, 8, 7], [8, 9, 7]])
    az = np.array([2e-5, 2e-5, 2e-5, 2e-5, 1e-5, 1e-5, 1e-5, 1e-5, 5e-6, 5e-6])
    body_ids = np.array([1, 1, 2, 2, 3, 3, 4, 4])
    vtu = {"points_mm": points, "az": az, "triangles": triangles, "body_ids": body_ids}
    psi = b3_elmer.flux_linkage(vtu, names, regions, stack_mm=20.0, sectors=4)
    mean_plus = (az[triangles[[0, 1]]].mean(axis=1) * 0.5).sum() / 1.0
    mean_minus = (az[triangles[[2, 3]]].mean(axis=1) * 0.5).sum() / 1.0
    assert psi["A"] == pytest.approx(4 * 0.020 * (10 * mean_plus - 10 * mean_minus), rel=1e-9)
    assert psi["C"] == 0.0 and set(psi) == {"A", "B", "C"}


def test_the_park_transform_gives_the_d_axis_flux_of_a_balanced_set_at_the_rotor_angle():
    from bench import b3_elmer
    theta = 37.0
    psi = [0.0166 * math.cos(math.radians(theta + s)) for s in (0, -120, 120)]
    d, q = b3_elmer.park(psi, theta)
    assert d == pytest.approx(0.0166) and q == pytest.approx(0.0, abs=1e-12)
    d, q = b3_elmer.park(psi, theta + 10)                                      # the same fluxes seen from a frame turned by 10 degrees
    assert math.hypot(d, q) == pytest.approx(0.0166) and q != 0


def test_the_gap_field_is_sampled_along_the_circle_between_the_nodes_on_it():
    from bench import b3_elmer
    angles = np.arange(0.0, 90.001, 1.5)                                       # nodes of the mid-gap circle, 1.5 degrees apart
    r = 45.4
    points = np.array([[r * math.cos(math.radians(a)), r * math.sin(math.radians(a))] for a in angles])
    b = np.array([[math.cos(math.radians(a)) * (0.3 + 0.1 * a / 90), math.sin(math.radians(a)) * (0.3 + 0.1 * a / 90)] for a in angles])      # radial: B_r = 0.3 + 0.1 a/90
    junk = np.array([[10.0, 10.0], [r, 5.0]])                                  # nodes that are not on the circle are ignored
    vtu = {"points_mm": np.vstack([junk, points[::-1]]), "b": np.vstack([[9.0, 9.0], [9.0, 9.0], b[::-1]])}          # in any order
    sample = np.array([0.0, 0.25, 10.0, 44.9, 90.0])
    br, bt = b3_elmer.gap_profile(vtu, r, sample)
    assert br == pytest.approx(0.3 + 0.1 * sample / 90, abs=1e-4) and bt == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError):
        b3_elmer.gap_profile({"points_mm": junk, "b": np.array([[1.0, 0.0], [1.0, 0.0]])}, r, sample)          # no node on that circle


def test_two_waveforms_are_compared_by_their_fundamentals_their_normalised_harmonics_and_their_shape():
    from bench import b3_elmer
    theta = np.arange(0, 90, 0.25)
    femm = sum(a * np.cos(2 * np.pi * n * theta / 90) for n, a in ((1, 0.32), (3, 0.065), (5, 0.078), (9, 0.044)))
    same = b3_elmer.compare_waveforms(femm, femm * 1.0)
    assert same["fundamental_ratio"] == pytest.approx(1.0) and same["rms_percent_of_fundamental"] == pytest.approx(0.0, abs=1e-9) and same["alike"] is True
    scaled = b3_elmer.compare_waveforms(femm, femm * 0.97)                     # an amplitude error alone is not a different pattern
    assert scaled["fundamental_percent"] == pytest.approx(-3.0) and scaled["alike"] is True
    assert scaled["top_orders_femm"] == scaled["top_orders_other"] == [5, 3, 9]
    other = femm + 0.05 * np.cos(2 * np.pi * 5 * theta / 90 + 1.0)             # a 5th harmonic 5 points of the fundamental off: another pattern
    different = b3_elmer.compare_waveforms(femm, other)
    assert different["alike"] is False and max(abs(row["difference_points"]) for row in different["harmonics"]) > 3.0
    assert [row["order"] for row in same["harmonics"]] == list(range(2, 21))
    # each of the three limits alone can fail the pattern (the limits are the bench's own, so each is pinned by a case that only it catches):
    one_harmonic = femm + 0.04 * 0.32 * np.cos(2 * np.pi * 7 * theta / 90)        # a 7th harmonic of 4 points: small in rms, the largest three unchanged
    only_points = b3_elmer.compare_waveforms(femm, one_harmonic)
    assert only_points["rms_percent_of_fundamental"] < b3_elmer.PATTERN_RMS_PERCENT and only_points["top_orders_other"] == only_points["top_orders_femm"]
    assert only_points["alike"] is False and b3_elmer.PATTERN_POINTS == 3.0
    femm2 = femm + 0.043 * np.cos(2 * np.pi * 13 * theta / 90)                    # a 13th harmonic just below the 9th (0.044): the largest three are 5, 3, 9
    other2 = femm + (-0.002) * np.cos(2 * np.pi * 9 * theta / 90) + 0.045 * np.cos(2 * np.pi * 13 * theta / 90)     # ... and just above it: 5, 3, 13
    only_order = b3_elmer.compare_waveforms(femm2, other2)
    assert only_order["top_orders_femm"] == [5, 3, 9] and only_order["top_orders_other"] == [5, 3, 13]
    assert max(abs(row["difference_points"]) for row in only_order["harmonics"]) < b3_elmer.PATTERN_POINTS and only_order["rms_percent_of_fundamental"] < b3_elmer.PATTERN_RMS_PERCENT
    assert only_order["alike"] is False                                            # every harmonic within 3 points, the shapes within 10 %: only the order of the largest differs
    assert b3_elmer.PATTERN_RMS_PERCENT == 10.0 and b3_elmer.PATTERN_TOP == 3


def test_the_two_periodic_edges_and_the_two_sides_of_the_gap_are_checked_for_matching_nodes_before_elmer_is_asked():
    from bench import b3_elmer
    radii = np.array([7.5, 9.0, 20.0, 45.0, 64.0])
    ok = b3_elmer.nodes_match(radii, radii[::-1] + 1e-13, 1e-9)
    assert ok["match"] is True and ok["count"] == (5, 5) and ok["max_difference"] < 1e-9
    assert b3_elmer.nodes_match(radii, radii[:-1], 1e-9)["match"] is False                       # a node more on one edge
    shifted = radii.copy()
    shifted[2] += 0.01
    bad = b3_elmer.nodes_match(radii, shifted, 1e-9)
    assert bad["match"] is False and bad["max_difference"] == pytest.approx(0.01)


def test_the_real_mesh_of_path_1_has_matching_periodic_edges_and_a_gap_that_merges(tmp_path):
    from bench import b2_gmsh, b3_elmer
    b2_gmsh.make_path1_mesh(tmp_path / "path1.msh")
    facts = b3_elmer.interface_facts(tmp_path / "path1.msh")
    assert facts["periodic_edges"]["match"] is True and facts["periodic_edges"]["count"] == (56, 56)
    assert facts["gap_interface"]["match"] is True and facts["gap_interface"]["count"] == (269, 269)
    assert facts["gap_interface"]["max_difference"] < 1e-9


def test_interface_facts_leaves_the_path_of_this_process_whole(tmp_path, path_is_put_back):
    """The process's own PATH, read from Windows (os.environ is never cut, so it proves nothing). Self-check: take `with intact_path():` off
    interface_facts and this fails."""
    from bench import b2_gmsh, b3_elmer
    b2_gmsh.make_path1_mesh(tmp_path / "path1.msh")
    assert process_path() == path_is_put_back
    b3_elmer.interface_facts(tmp_path / "path1.msh")
    assert process_path() == path_is_put_back


def test_a_mesh_elmer_cannot_use_is_a_failed_check_that_says_so_and_not_an_exception_while_it_is_worded():
    """When the periodic edges of a mesh do not match (Elmer would say "This is not conforming!") some numbers do not exist: the check has to be worded
    without them."""
    from bench import b3_elmer
    theta = np.arange(0, 90, 0.25)
    wave = 0.32 * np.cos(2 * np.pi * theta / 90) + 0.07 * np.cos(2 * np.pi * 3 * theta / 90) + 0.07 * np.cos(2 * np.pi * 5 * theta / 90)
    comparison = b3_elmer.compare_waveforms(wave, wave)

    def run(finished):
        return {"grid": {"conforming": None, "nodes_merged": None, "nodes_read": None, "pieces": None},
                "solver": {"finished": finished, "errors": ["ConformingNodePerm: This is not conforming!"], "last_change": None, "iterations": 0, "timed_out": False,
                           "solver_seconds_real": None, "threads": 4, "version": "26.2-devel", "revision": "11c118114", "warnings": []},
                "numbers": {"periodic_relative_difference": None, "az_max": 0.0, "comparison": comparison, "psi_pm_mWb": 16.5}}
    facts = {"periodic_edges": {"match": False, "count": (92, 91), "max_difference": None}, "gap_interface": {"match": True, "count": (269, 269), "max_difference": 1e-13}}
    checks = b3_elmer.make_checks({"laminated": run(False), "plain": run(False)}, facts, {"psi_pm_mWb": 16.6}, 16.5, [{"size_max_mm": 0.7, "periodic_edges_match": False}],
                                  MADE)
    verdicts = {c["name"]: c["verdict"] for c in checks}
    assert verdicts["mesh: the periodic edges and the two sides of the gap carry matching nodes, and ElmerGrid makes one conforming mesh"] == "FAIL"
    assert verdicts["laminated: ElmerSolver converges"] == "FAIL" and verdicts["plain: the periodic condition holds in the solution"] == "FAIL"
    assert verdicts["laminated: the agreement holds on finer meshes"] == "FAIL"
    assert any("n/a" in c["measured"] for c in checks if isinstance(c["measured"], str))


# ----------------------------------------------------------------------------------------------- 2. the verdict and the review's fixes (Task 11, fix round 1)
def synthetic_run(ratio):
    """One Elmer run as run() hands it to make_checks, with a B_r waveform that is FEMM's times `ratio` (1.0 agrees, 0.9 is 10 % low)."""
    theta = np.arange(0, 90, 0.25)
    wave = 0.32 * np.cos(2 * np.pi * theta / 90) + 0.07 * np.cos(2 * np.pi * 3 * theta / 90) + 0.07 * np.cos(2 * np.pi * 5 * theta / 90)
    from bench import b3_elmer
    return {"grid": {"conforming": True, "nodes_merged": 269, "nodes_read": 8819, "nodes_merged_to": 8550, "pieces": 1},
            "solver": {"finished": True, "errors": [], "last_change": 2.1e-7, "iterations": 11, "timed_out": False, "solver_seconds_real": 1.4, "threads": 4,
                       "version": "26.2-devel", "revision": "11c118114", "warnings": []},
            "numbers": {"periodic_relative_difference": 0.0, "az_max": 0.01, "comparison": b3_elmer.compare_waveforms(wave, wave * ratio), "psi_pm_mWb": 16.6 * ratio}}


MADE = {"nodes_in_file": 8819, "to_gmsh_nodes": 8906, "nodes_without_element": 87, "boundary_curves": {"A0": 4, "edge_0": 8}, "size_max": 1.5}      # B2's mesh as make_path1_mesh says


def synthetic_checks(laminated_ratio, plain_ratio, made=None):
    from bench import b3_elmer
    facts = {"periodic_edges": {"match": True, "count": (56, 56), "max_difference": 0.0}, "gap_interface": {"match": True, "count": (269, 269), "max_difference": 1e-13}}
    sensitivity = [{"size_max_mm": 0.5, "periodic_edges_match": True, "gap_interface_match": True, "nodes_in_file": 21867, "elmergrid_read": 21867, "elmer_nodes": 21598,
                    "merged_in_the_gap": 269, "iterations": 11, "fundamental_percent": -1.08,
                    "psi_pm_mWb": 16.35, "psi_percent": -1.48, "rms_percent_of_fundamental": 1.2, "pattern_alike": True}]
    runs = {"laminated": synthetic_run(laminated_ratio), "plain": synthetic_run(plain_ratio)}
    return b3_elmer.make_checks(runs, facts, {"psi_pm_mWb": 16.6}, 16.6, sensitivity, MADE if made is None else made)


def test_the_laminated_run_carries_the_verdict_and_the_plain_run_is_reported_beside_it_and_does_not_gate():
    """I2 of the review of Task 11 (the ledger's ruling): the plain 47-point run solves another steel than FEMM's, so its misses are reported and are not
    the case's verdict. Here the plain run misses everything (10 % low) and the laminated one agrees: every plain check still FAILS, and none gates."""
    import bench
    from bench import b3_elmer
    assert b3_elmer.VERDICT_RUN == "laminated"
    checks = synthetic_checks(0.997, 0.90)
    plain = [c for c in checks if c["name"].startswith("plain:")]
    assert len(plain) == 6 and all(c["gating"] is False for c in plain)                          # converges, periodic, fundamental, pattern, two psi
    assert [c["verdict"] for c in plain if "converges" in c["name"] or "periodic" in c["name"]] == ["PASS", "PASS"]
    assert sum(1 for c in plain if c["verdict"] == "FAIL") == 4                                    # the numbers it misses are still there, as FAIL
    for c in plain:
        assert bench.describe(c)[0].startswith("reported, not gating (the verdict is the laminated run's)"), c["name"]
    shift = next(c for c in checks if c["name"].startswith("the lamination fill"))              # it is about the plain run too
    assert shift["gating"] is False and shift["expected"].startswith("reported, not gating (the verdict is the laminated run's)")
    verdict = [c for c in checks if c not in plain and c is not shift]
    assert verdict and all(c["gating"] is True and c["verdict"] == "PASS" for c in verdict), [c["name"] for c in verdict if c["verdict"] != "PASS"]
    assert {c["name"].split(":")[0] for c in verdict} == {"laminated", "mesh"}


def test_a_laminated_miss_still_gates_the_verdict():
    checks = synthetic_checks(0.90, 1.0)
    failed = [c for c in checks if c["verdict"] == "FAIL"]
    assert failed and all(c["gating"] is True for c in failed if c["name"].startswith("laminated:"))
    assert any(c["name"].startswith("laminated: B_r") and c["verdict"] == "FAIL" and c["gating"] is True for c in checks)


def test_the_three_node_counts_are_named_each_with_what_it_counts_and_b3_agrees_with_b2_on_the_mesh():
    """M5 of the review of Task 11: "8,906 nodes in the file" was wrong. The saved file has 8,819 (what ElmerGrid reads); 8,906 is the number to_gmsh returns, which
    includes the 87 arc centres that no element uses; Elmer has 8,550 after ElmerGrid merged the 269 nodes of the gap's two sides."""
    from bench import b3_elmer
    grid = b3_elmer.parse_grid_log(GRID_LOG)
    counts = b3_elmer.node_counts(MADE, grid)
    assert (counts["nodes_in_file"], counts["to_gmsh_nodes"], counts["nodes_without_element"]) == (8819, 8906, 87)
    assert (counts["elmergrid_read"], counts["merged_in_the_gap"], counts["elmer_nodes"]) == (8819, 269, 8550)
    assert counts["elmergrid_read_b2s_file"] is True
    assert set(counts["meaning"]) == {"nodes_in_file", "to_gmsh_nodes", "nodes_without_element", "elmergrid_read", "merged_in_the_gap", "elmer_nodes"}
    assert counts["nodes_in_file"] - counts["merged_in_the_gap"] == counts["elmer_nodes"]
    assert b3_elmer.node_counts(dict(MADE, nodes_in_file=8800), grid)["elmergrid_read_b2s_file"] is False        # B3 read another file than B2 wrote
    checks = synthetic_checks(1.0, 1.0)
    mesh = next(c for c in checks if c["name"].startswith("mesh:"))
    assert mesh["verdict"] == "PASS" and mesh["gating"] is True
    for text in ("8,819 nodes in the file", "8,906 counted by to_gmsh", "87", "8,550", "269"):
        assert text in mesh["measured"], text
    other = next(c for c in synthetic_checks(1.0, 1.0, made=dict(MADE, nodes_in_file=8800)) if c["name"].startswith("mesh:"))
    assert other["verdict"] == "FAIL"                                                                    # ElmerGrid did not read the file B2 made


def test_the_lamination_shift_check_says_that_its_limit_is_the_benchs_own():
    import bench
    shift = next(c for c in synthetic_checks(1.0, 0.9) if c["name"].startswith("the lamination fill"))
    assert "(" + bench.OWN_LIMIT + ")" in shift["expected"]


def test_the_briefs_two_3_percent_limits_are_pinned_and_are_the_ones_the_comparisons_use():
    """M10 of the review of Task 11: B1's 0.5 % was pinned to the brief; B3's two 3 % limits (the air-gap flux density's fundamental, psi_PM) were not."""
    from bench import b3_elmer
    assert b3_elmer.FUNDAMENTAL_TOLERANCE_PERCENT == 3.0 and b3_elmer.PSI_TOLERANCE_PERCENT == 3.0
    used = {}
    for c in synthetic_checks(1.0, 1.0):
        if "tolerance_percent" in c:
            used.setdefault(c["name"].split(": ", 1)[1].split(",")[0].split(" from ")[0].split(" against ")[0], set()).add(c["tolerance_percent"])
    assert used and set().union(*used.values()) == {3.0}
    assert sum(1 for c in synthetic_checks(1.0, 1.0) if "tolerance_percent" in c) == 6          # fundamental, psi, psi at the position: laminated and plain


def test_the_headline_agreement_is_never_quoted_without_the_finer_mesh_trend():
    """M6 of the review of Task 11: "a third of a percent" is the figure on the default mesh (size_max 1.5). Beside it the module says that the difference grows
    on finer meshes (-1.08 % in the B_r fundamental and -1.48 % in psi_PM at size_max 0.5), names the one known contributor (the N45 area, 0.172 % short at every mesh
    size) and says that the rest is not isolated."""
    text = (P / "bench" / "b3_elmer.py").read_text(encoding="utf-8")
    places = [i for i in range(len(text)) if text.startswith("third of a percent", i)]
    assert places, "the headline is no longer quoted: update this test with the text"
    for at in places:
        near = text[max(0, at - 900):at + 1200]
        for needle in ("1.08 %", "1.48 %", "0.172 %", "not isolated", "default mesh"):
            assert needle in near, f"{needle!r} is not beside 'third of a percent' (character {at})"


# ----------------------------------------------------------------------------------------------- 3. the re-review's minors (Task 11, fix round 2)
def test_the_finer_mesh_counts_are_elmers_after_the_merge_and_say_so_beside_the_file_count():
    """P3 of the re-review of Task 11 (M5 again, for the finer meshes): B2's table counts the nodes of the FILE (10,711 at size_max 1.0), the Elmer run counts what is
    left after ElmerGrid merged the gap's two sides (10,442). The finer-mesh check names both, and so does each row of results.json."""
    finer = next(c for c in synthetic_checks(1.0, 1.0) if c["name"].startswith("laminated: the agreement holds on finer meshes"))
    assert "21,598 nodes in Elmer after ElmerGrid's merge (21,867 in the file, 269 merged across the gap)" in finer["measured"]


def test_a_row_of_the_finer_mesh_table_names_the_file_count_and_elmers_count_apart():
    from bench import b3_elmer
    edge = {"periodic_edges": {"match": True, "count": (56, 56), "max_difference": 0.0}, "gap_interface": {"match": True, "count": (269, 269), "max_difference": 0.0}}
    elmer_run = {"grid": {"nodes_read": 10711, "nodes_merged": 269, "nodes_merged_to": 10442}, "solver": {"iterations": 11}}
    numbers = {"comparison": {"fundamental_percent": -0.84, "rms_percent_of_fundamental": 4.5, "alike": True}, "psi_pm_mWb": 16.443}
    row = b3_elmer.sensitivity_row(1.0, 10711, edge, elmer_run, numbers, 16.629)
    assert (row["nodes_in_file"], row["elmergrid_read"], row["merged_in_the_gap"], row["elmer_nodes"]) == (10711, 10711, 269, 10442) and "nodes" not in row
    assert row["fundamental_percent"] == -0.84 and row["psi_percent"] == pytest.approx(100 * (16.443 / 16.629 - 1)) and row["pattern_alike"] is True
    skipped = b3_elmer.sensitivity_row(0.7, 12000, dict(edge, periodic_edges={"match": False, "count": (92, 91), "max_difference": None}))
    assert skipped == {"size_max_mm": 0.7, "nodes_in_file": 12000, "periodic_edges_match": False, "gap_interface_match": True}      # a mesh Elmer cannot use is not run


def test_the_convergence_limit_and_the_1e_9_limits_say_that_they_are_the_benchs_own():
    """P4 of the re-review of Task 11: the brief says "the solver converges, and the iteration count is reported" and gives no number: the 1e-6 and the 60 iterations are
    the bench's. So are the 1e-9 of the periodic condition and of the nodes that coincide."""
    import bench
    own = "(" + bench.OWN_LIMIT + ")"
    checks = synthetic_checks(1.0, 1.0)
    assert own in next(c for c in checks if c["name"] == "laminated: ElmerSolver converges")["expected"]
    assert own in next(c for c in checks if c["name"] == "laminated: the periodic condition holds in the solution")["expected"]
    assert own in next(c for c in checks if c["name"].startswith("mesh:"))["expected"]


def test_the_cited_manual_sentence_is_quoted_whole_and_the_ruled_wording_is_used():
    """P5 and P6 of the re-review of Task 11: the ElmerSolver Manual 9.3 sentence is quoted without the cut that was not marked; the rotor-surface chords are "at most
    they offset it on the default mesh" (the ruled words, not "are not the cause"); the question's 0.3 % says on which mesh."""
    text = (P / "bench" / "b3_elmer.py").read_text(encoding="utf-8")
    assert "a direct node-to-node mapping, or edge-to-edge mapping on the two boundaries" in text
    assert "are not the cause" not in text and "at most they offset it on the default mesh" in text
    assert "the laminated one 0.3 % (on the default mesh)" in text
