"""Tests for B2 (bench/b2_gmsh.py): Python -> Gmsh. Everything that needs no solver; Gmsh itself runs (it is a library, and the meshes are tiny).

Run from P:   uv run pytest tests/test_bench_b2.py -q      (a few seconds; no FEMM, no Elmer, no solver slot is taken)

In groups:
  1. The mesh read back (groups, areas, quality, the worst element), the PATH of the process after each function that starts Gmsh, the chord error,
     the CAD drawings (faces, bands removed in memory, meshed with exactly the faces' areas, the rotation against motorgeo's model), and the real
     drawings.
  2. The review's fixes: the mesh edges on a fillet counted from the mesh's nodes (path 1 and path 2) and the chord error held to 1 %, the node counts
     of path 1's mesh, the air-band check computed from what was found, every LINE, ARC and CIRCLE an edge of a face, path 2's mesh against
     motorgeo's areas and its slot positions, the limits that are the bench's own and the brief's two (0.1 % and 0.5 %). The real drawings and meshes
     are made once for the group (a few seconds).
  3. Fix round 2: path 2 reads SynRM rotor.DXF, the CAD's drawing that HAS the air band (Rotor.DXF never had it), and removes the band in memory. The band is
     counted from motorgeo and each removed band held against motorgeo's; the dropped edges are shown to be exactly the bands' edges (a band that shares a wall with
     a slot keeps it); after the removal the rotor IS Rotor.DXF's, with a limit derived from the arc flattening; with the removal switched off those checks fail; the
     three drawings are hashed before and after; bands.png and the words that name the drawing.
"""
import math

import numpy as np
import pytest

import motorgeo as mg
import paths
from bench_test_tools import tiny_mesh
from conftest import process_path


# ----------------------------------------------------------------------------------------------- 1. Gmsh, the meshes and the CAD drawings
def write_synthetic_cad(tmp_path, bands=True):
    """A rotor and a stator drawn like the CAD exports, small enough to know every area: the rotor is a disk of radius 10 with a shaft hole of
    radius 2, two magnet slots (rectangles 4 x 1) and, if `bands`, two air bands (rectangles 6 x 1.5); the stator is a ring (radius 20) whose bore
    (radius 11) is three tooth tops (arcs) and three slots (two radial lines and an arc each, 20 degrees wide, down to radius 15)."""
    import ezdxf
    rotor, stator = ezdxf.new("R2000"), ezdxf.new("R2000")
    msp = rotor.modelspace()
    msp.add_circle((0, 0), 2.0)
    msp.add_circle((0, 0), 10.0)

    def rectangle(space, x0, x1, y0, y1):
        for a, b in (((x0, y0), (x1, y0)), ((x1, y0), (x1, y1)), ((x1, y1), (x0, y1)), ((x0, y1), (x0, y0))):
            space.add_line(a, b)

    rectangle(msp, 4, 8, 2.5, 3.5)
    rectangle(msp, -8, -4, -3.5, -2.5)
    if bands:
        rectangle(msp, -3, 3, 5.5, 7)
        rectangle(msp, -3, 3, -7, -5.5)
    msp.add_mtext("SOLIDWORKS Educational Product.").dxf.insert = (-10, -12)          # text in the drawing: no geometry
    s = stator.modelspace()
    s.add_circle((0, 0), 20.0)
    for a0, a1 in ((10, 110), (130, 230), (250, 350)):
        s.add_arc((0, 0), 11.0, a0, a1)
    for a0, a1 in ((110, 130), (230, 250), (350, 10)):
        for a in (a0, a1):
            c, d = math.cos(math.radians(a)), math.sin(math.radians(a))
            s.add_line((11 * c, 11 * d), (15 * c, 15 * d))
        s.add_arc((0, 0), 15.0, a0, a1)
    rotor.saveas(tmp_path / "rotor.dxf")
    stator.saveas(tmp_path / "stator.dxf")
    return tmp_path / "rotor.dxf", tmp_path / "stator.dxf"


def test_the_expected_physical_groups_are_the_regions_of_the_build_and_the_five_boundaries():
    from bench import b2_gmsh
    regions = mg.build(mg.Params(barrier="none"), gap="sliding")
    names = b2_gmsh.expected_groups(regions)
    assert len(names) == 35 + 5 == len(set(names))
    assert {"rotor_steel", "stator_steel", "gap_rotor", "gap_stator", "magnet_0_up", "magnet_1_dn", "opening_8", "coil_8_1"} <= set(names)
    assert set(b2_gmsh.BOUNDARY_NAMES) == {"A0", "edge_0", "edge_sector", "sliding_rotor", "sliding_stator"}
    assert set(b2_gmsh.BOUNDARY_NAMES) <= set(names)
    assert not any(n.startswith("barrier_") for n in names)                  # R11: no band, so no group for it


def test_a_mesh_is_read_back_group_by_group_with_its_areas_its_quality_and_the_worst_element(tmp_path):
    from bench import b2_gmsh
    tiny_mesh(tmp_path / "tiny.msh")
    info = b2_gmsh.inspect_mesh(tmp_path / "tiny.msh")
    groups = {g["name"]: g for g in info["groups"]}
    assert set(groups) == {"left", "right", "bottom"}
    assert groups["left"]["dim"] == 2 and groups["left"]["area_mm2"] == pytest.approx(0.4, rel=1e-9)
    assert groups["right"]["area_mm2"] == pytest.approx(0.6, rel=1e-9)
    assert groups["bottom"]["dim"] == 1 and groups["bottom"]["length_mm"] == pytest.approx(1.0, rel=1e-9) and groups["bottom"]["elements"] >= 10
    assert info["triangles"] == groups["left"]["elements"] + groups["right"]["elements"]
    assert info["nodes"] > 100 and info["total_area_mm2"] == pytest.approx(1.0, rel=1e-9)
    q = info["quality"]
    assert 0.2 < q["min_sicn"] <= q["mean_sicn"] <= 1.0
    worst = q["worst"]
    assert worst["group"] in ("left", "right") and 0 <= worst["x_mm"] <= 1 and 0 <= worst["y_mm"] <= 1
    assert worst["r_mm"] == pytest.approx(math.hypot(worst["x_mm"], worst["y_mm"]))


# What Gmsh does to the PATH, and the four bench functions that start Gmsh. The PATH that matters is the process's own, which Windows keeps
# apart from os.environ (Python's copy is never cut, so os.environ["PATH"] == before proves nothing): process_path() in conftest.py reads it
# from Windows. Self-check: take `with intact_path():` off the function and its test fails (gmsh 4.15.2: 1,864 characters became 313).
def test_inspect_mesh_leaves_the_path_of_this_process_whole(tmp_path, path_is_put_back):
    from bench import b2_gmsh
    tiny_mesh(tmp_path / "tiny.msh")
    assert process_path() == path_is_put_back                           # the start is whole: a test that cut it earlier would show here
    b2_gmsh.inspect_mesh(tmp_path / "tiny.msh")
    assert process_path() == path_is_put_back


def test_make_path1_mesh_leaves_the_path_of_this_process_whole(tmp_path, path_is_put_back):
    from bench import b2_gmsh
    assert process_path() == path_is_put_back
    b2_gmsh.make_path1_mesh(tmp_path / "path1.msh")
    assert process_path() == path_is_put_back


def test_mesh_faces_leaves_the_path_of_this_process_whole(tmp_path, path_is_put_back):
    from bench import b2_gmsh
    rotor, stator = write_synthetic_cad(tmp_path, bands=True)
    faces = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor), b2_gmsh.read_dxf(stator), slot_size=(4.0, 1.0))
    assert process_path() == path_is_put_back
    b2_gmsh.mesh_faces(faces, tmp_path / "cad.msh", size_max=2.0)
    assert process_path() == path_is_put_back


def test_material_areas_add_the_regions_of_each_material_and_geometry_areas_come_from_the_visible_polygons():
    from bench import b2_gmsh
    p = mg.Params(barrier="none")
    regions = mg.build(p, gap="sliding")
    geometry = b2_gmsh.geometry_areas(regions)
    assert set(geometry) == {"M-19 Steel", "N45", "Air", "Coil"}
    assert geometry["N45"] == pytest.approx(4 * 16.3483, abs=2e-3)                    # four slots of the sector, 16.35 mm2 each (README: 261.57 for 16)
    assert geometry["M-19 Steel"] == pytest.approx(1480.854 + 1082.364, abs=1e-2)
    assert sum(geometry.values()) == pytest.approx(math.pi / 4 * (p.r_so ** 2 - p.r_shaft ** 2), rel=1e-4)         # the sector, whole
    visible = mg.visible(regions)
    exact = {r.name: visible[r.name].area for r in regions}
    assert b2_gmsh.material_areas(exact, regions) == pytest.approx(geometry)


def test_the_chord_error_of_a_fillet_meshed_in_two_pieces_is_what_the_magnet_area_is_short_by():
    """The README's magnet slots have 0.3 mm fillets, and to_gmsh meshes a fillet with 2 straight edges at mesh size 0.5: the mesh covers a little less
    than the drawing. The formula is the area between an arc and its n chords; here it is checked against shapely."""
    from shapely.geometry import Polygon
    from bench import b2_gmsh
    r, sweep, n = 0.3, 90.0, 2
    pie = lambda k: Polygon([(0, 0)] + [(r * math.cos(math.radians(sweep * i / k)), r * math.sin(math.radians(sweep * i / k))) for i in range(k + 1)])
    assert b2_gmsh.chord_deficit(r, sweep, n) == pytest.approx(pie(100000).area - pie(n).area, rel=1e-6)
    assert b2_gmsh.chord_deficit(r, sweep, 1) > b2_gmsh.chord_deficit(r, sweep, 2) > b2_gmsh.chord_deficit(r, sweep, 4) > 0
    assert b2_gmsh.chord_deficit(r, sweep, 2) == pytest.approx(0.045 * (math.pi / 2 - 2 * math.sin(math.pi / 4)), rel=1e-12)
    assert b2_gmsh.chord_deficit(r, sweep, 2) * 4 == pytest.approx(0.0282, abs=2e-4)                 # four fillets per slot


def test_the_dxf_is_read_into_edges_and_counts_and_the_text_in_it_is_not_geometry(tmp_path):
    from bench import b2_gmsh
    rotor_path, stator_path = write_synthetic_cad(tmp_path)
    rotor = b2_gmsh.read_dxf(rotor_path)
    assert rotor["counts"] == {"CIRCLE": 2, "LINE": 16, "MTEXT": 1}
    assert sorted(c[2] for c in rotor["circles"]) == [2.0, 10.0] and len(rotor["lines"]) == 18        # 16 lines + 2 circles (flattened)
    stator = b2_gmsh.read_dxf(stator_path)
    assert stator["counts"] == {"CIRCLE": 1, "ARC": 6, "LINE": 6}
    centred = [a for a in stator["arcs"] if abs(a[0]) < 1e-9 and abs(a[1]) < 1e-9]
    assert len(centred) == 6 and min(a[2] for a in centred) == 11.0


def test_the_slot_mouths_close_the_bore_between_the_tooth_tops(tmp_path):
    from bench import b2_gmsh
    _, stator_path = write_synthetic_cad(tmp_path)
    mouths = b2_gmsh.bore_mouths(b2_gmsh.read_dxf(stator_path))
    assert len(mouths) == 3
    for m in mouths:                                                  # each is a 20-degree arc of radius 11: its length is 11 * 20 degrees in radians
        assert m.length == pytest.approx(11 * math.radians(20), rel=1e-4)
        assert all(abs(math.hypot(x, y) - 11) < 1e-3 for x, y in m.coords)
    assert b2_gmsh.bore_mouths(b2_gmsh.read_dxf(tmp_path / "rotor.dxf")) == []       # a drawing with no bore arcs has no mouths


def test_the_faces_of_a_drawing_with_bands_come_out_with_the_bands_removed_and_every_area_known(tmp_path):
    from bench import b2_gmsh
    rotor, stator = write_synthetic_cad(tmp_path, bands=True)
    faces = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor), b2_gmsh.read_dxf(stator), slot_size=(4.0, 1.0))
    assert faces["bands_removed"] == 2 and len(faces["magnets"]) == 2 and len(faces["slots"]) == 3
    assert faces["radii"] == {"shaft": 2.0, "rotor": 10.0, "bore": 11.0, "stator_outer": 20.0}
    rel = lambda polygon, exact: polygon.area / exact - 1
    assert rel(faces["rotor_steel"], math.pi * (100 - 4) - 2 * 4.0) == pytest.approx(0, abs=2e-4)     # the bands are steel now: only the magnets are holes
    assert all(abs(rel(m, 4.0)) < 1e-9 for m in faces["magnets"])
    slot = (20 / 360) * math.pi * (15 ** 2 - 11 ** 2)
    assert all(abs(rel(s, slot)) < 2e-4 for s in faces["slots"])
    assert rel(faces["gap"], math.pi * (121 - 100)) == pytest.approx(0, abs=2e-4)
    assert rel(faces["stator_steel"], math.pi * (400 - 121) - 3 * slot) == pytest.approx(0, abs=2e-4)
    assert faces["shaft"].area == pytest.approx(math.pi * 4, rel=2e-4)
    assert len(faces["rotor_steel"].interiors) == 3 and len(faces["stator_steel"].interiors) == 1        # shaft + 2 magnets; the stator's inner contour


def test_a_drawing_with_no_bands_loses_nothing_and_says_that_it_removed_none(tmp_path):
    from bench import b2_gmsh
    rotor, stator = write_synthetic_cad(tmp_path, bands=False)
    faces = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor), b2_gmsh.read_dxf(stator), slot_size=(4.0, 1.0))
    assert faces["bands_removed"] == 0 and len(faces["magnets"]) == 2
    assert faces["rotor_steel"].area == pytest.approx(math.pi * 96 - 8.0, rel=2e-4)


def test_the_faces_are_meshed_into_groups_that_have_exactly_the_areas_of_the_faces(tmp_path):
    from bench import b2_gmsh
    rotor, stator = write_synthetic_cad(tmp_path, bands=True)
    faces = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor), b2_gmsh.read_dxf(stator), slot_size=(4.0, 1.0))
    info = b2_gmsh.mesh_faces(faces, tmp_path / "cad.msh", size_max=2.0)
    assert info["nodes"] > 200
    got = b2_gmsh.inspect_mesh(tmp_path / "cad.msh")
    names = {g["name"] for g in got["groups"]}
    assert names == set(b2_gmsh.cad_group_names(faces))
    assert {"rotor_steel", "stator_steel", "gap", "magnet_1", "magnet_2", "slot_1", "slot_3", "shaft", "outer"} <= names
    area = {g["name"]: g["area_mm2"] for g in got["groups"] if g["dim"] == 2}
    assert area["rotor_steel"] == pytest.approx(faces["rotor_steel"].area, rel=1e-9)
    assert area["magnet_1"] == pytest.approx(faces["magnets"][0].area, rel=1e-9)
    assert area["gap"] == pytest.approx(faces["gap"].area, rel=1e-9) and area["stator_steel"] == pytest.approx(faces["stator_steel"].area, rel=1e-9)
    length = {g["name"]: g["length_mm"] for g in got["groups"] if g["dim"] == 1}
    assert length["outer"] == pytest.approx(2 * math.pi * 20, rel=1e-3) and length["shaft"] == pytest.approx(2 * math.pi * 2, rel=1e-3)
    assert got["quality"]["min_sicn"] > b2_gmsh.MIN_SICN_OK                  # the bench's own floor for "no degenerate element"


def test_the_rotation_between_a_drawing_and_the_model_is_found_modulo_the_symmetry_of_the_pattern():
    from bench import b2_gmsh
    ref = np.array([[40 * math.cos(math.radians(a)), 40 * math.sin(math.radians(a))] for k in range(8) for a in (45 * k + 13.85, 45 * k + 31.15)])
    for offset in (22.5, 7.3):
        c, s = math.cos(math.radians(-offset)), math.sin(math.radians(-offset))
        cad = ref @ np.array([[c, -s], [s, c]]).T                      # the drawing is the model turned by -offset
        found, residual = b2_gmsh.cad_offset(cad, ref, period_deg=45.0)
        assert found == pytest.approx(offset, abs=1e-6) and residual < 1e-6
    shifted = ref + np.array([0.5, 0.0])                                # a drawing that is not the same pattern: the residual says so
    assert b2_gmsh.cad_offset(shifted, ref, period_deg=45.0)[1] > 0.1


def test_rotor_dxf_the_earlier_export_has_no_band_and_is_the_motor_with_the_areas_of_motorgeo():
    """What the bench found in reference\\Rotor.DXF (the CAD's export of 2026-01-13) and reference\\Stator_refined.DXF (read-only copies of the CAD exports): the rotor
    drawing Rotor.DXF has the two circles and sixteen magnet slots and NO air band (build_faces looks for one by the shape of each pocket and finds none), and every
    area is motorgeo's no-band area to better than 0.01 % (0.1 % for the slots). Since fix round 2 path 2 reads SynRM rotor.DXF, which has the band; Rotor.DXF stays
    (config.toml: dxf_rotor) as the independent check of the removal, and this test keeps what it is."""
    from bench import b2_gmsh
    p = mg.Params(barrier="none")
    rotor, stator = b2_gmsh.read_dxf(paths.path("dxf_rotor")), b2_gmsh.read_dxf(paths.path("dxf_stator"))
    assert rotor["file"] == "Rotor.DXF" and stator["file"] == "Stator_refined.DXF"
    assert rotor["counts"] == {"LINE": 64, "ARC": 64, "CIRCLE": 2, "MTEXT": 2} and stator["counts"] == {"LINE": 216, "ARC": 72, "CIRCLE": 1, "MTEXT": 2}
    faces = b2_gmsh.build_faces(rotor, stator, slot_size=(p.slot_w, p.slot_t))
    assert faces["bands_removed"] == 0 and faces["bands"] == [] and faces["band_line_indices"] == [] and len(faces["magnets"]) == 16 and len(faces["slots"]) == 36
    assert faces["radii"] == pytest.approx({"shaft": 7.5, "rotor": 45.0, "bore": 45.8, "stator_outer": 64.0})
    visible = mg.visible(mg.build(p, gap="sliding"))
    sector = lambda name: p.n_sectors * visible[name].area
    assert faces["rotor_steel"].area == pytest.approx(sector("rotor_steel"), rel=1e-4)
    assert faces["stator_steel"].area == pytest.approx(sector("stator_steel"), rel=1e-4)
    assert sum(m.area for m in faces["magnets"]) == pytest.approx(sum(sector(n) for n in visible if n.startswith("magnet_")), rel=1e-4)
    assert sum(s.area for s in faces["slots"]) == pytest.approx(sum(sector(n) for n in visible if n.startswith(("opening_", "coil_"))), rel=1e-3)


def test_synrm_rotor_dxf_has_the_air_band_16_pockets_that_are_removed_in_memory_and_the_rest_is_the_motor():
    """SynRM rotor.DXF (the CAD's export of 2026-01-22, config.toml: dxf_rotor_with_band) has what Rotor.DXF does not: sixteen air-band pockets, one beside each magnet slot.
    Read with the bench's own build_faces it gives 16 bands found and removed and exactly their 96 edges dropped (80 LINE + 16 ARC: a band's outline in motorgeo is 5 LINE
    + 1 ARC), none of them shared with a magnet slot; what is left has the areas of motorgeo's no-band motor, the same as Rotor.DXF's (the test above)."""
    from collections import Counter

    from bench import b2_gmsh
    p = mg.Params(barrier="none")
    rotor, stator = b2_gmsh.read_dxf(paths.path("dxf_rotor_with_band")), b2_gmsh.read_dxf(paths.path("dxf_stator"))
    assert rotor["file"] == "SynRM rotor.DXF" and rotor["counts"] == {"LINE": 144, "ARC": 80, "CIRCLE": 2, "MTEXT": 2}
    assert len(rotor["kinds"]) == len(rotor["lines"]) == 226 and Counter(rotor["kinds"]) == {"LINE": 144, "ARC": 80, "CIRCLE": 2}
    faces = b2_gmsh.build_faces(rotor, stator, slot_size=(p.slot_w, p.slot_t))
    assert faces["bands_removed"] == 16 and len(faces["bands"]) == 16 and len(faces["magnets"]) == 16 and len(faces["slots"]) == 36
    assert faces["band_lines_dropped"] == 96 == len(faces["band_line_indices"]) and faces["band_edge_kinds"] == {"LINE": 80, "ARC": 16}
    assert faces["band_edges_shared_with_slots"] == 0                                    # a web of steel lies between each band and its slot: no edge is shared
    assert faces["radii"] == pytest.approx({"shaft": 7.5, "rotor": 45.0, "bore": 45.8, "stator_outer": 64.0})
    visible = mg.visible(mg.build(p, gap="sliding"))
    sector = lambda name: p.n_sectors * visible[name].area
    assert faces["rotor_steel"].area == pytest.approx(sector("rotor_steel"), rel=1e-4)
    assert faces["stator_steel"].area == pytest.approx(sector("stator_steel"), rel=1e-4)
    assert sum(m.area for m in faces["magnets"]) == pytest.approx(sum(sector(n) for n in visible if n.startswith("magnet_")), rel=1e-4)
    assert sum(s.area for s in faces["slots"]) == pytest.approx(sum(sector(n) for n in visible if n.startswith(("opening_", "coil_"))), rel=1e-3)
    assert faces["steel_before_removal_mm2"] == pytest.approx(5539.54, abs=0.01)         # the steel as drawn, with the 16 band pockets cut out of it ...
    assert faces["rotor_steel"].area - faces["steel_before_removal_mm2"] == pytest.approx(sum(b.area for b in faces["bands"]), abs=1e-6)       # ... and what the removal adds
    assert sum(b.area for b in faces["bands"]) == pytest.approx(383.86, abs=0.01)         # 16 x 23.99 mm2 (95.97 mm2 in the 90-degree sector)


# ----------------------------------------------------------------------------------------------- 2. the review's fixes (Task 11, fix round 1): fillets counted on the mesh, node counts, the band, every entity, path 2's mesh, the labels and pins
# The real products of the two paths, made once for the tests below (a few seconds each: Gmsh runs in this process).
@pytest.fixture(scope="module")
def real_path1(tmp_path_factory):
    from types import SimpleNamespace

    from bench import b2_gmsh
    folder = tmp_path_factory.mktemp("path1")
    p = mg.Params(barrier="none")
    regions = mg.build(p, gap="sliding")
    made = b2_gmsh.make_path1_mesh(folder / "path1.msh")
    info = b2_gmsh.inspect_mesh(folder / "path1.msh")
    chords = b2_gmsh.fillet_chords(b2_gmsh.mesh_nodes(folder / "path1.msh"), b2_gmsh.model_fillets(regions, p))
    return SimpleNamespace(p=p, regions=regions, made=made, info=info, report=b2_gmsh.path1_report(info, regions), chords=chords,
                           table=b2_gmsh.refinement_table(folder, regions))


@pytest.fixture(scope="module")
def real_path2(tmp_path_factory):
    """Path 2 as run() does it: SynRM rotor.DXF (the one with the air band, which is removed in memory) and Stator_refined.DXF, the faces made twice (the geometry
    and the coarser reading the mesh is made from), the mesh, and Rotor.DXF's faces beside them for the comparison."""
    from types import SimpleNamespace

    from bench import b2_gmsh
    folder = tmp_path_factory.mktemp("path2")
    p = mg.Params(barrier="none")
    regions = mg.build(p, gap="sliding")
    slot = (p.slot_w, p.slot_t)
    rotor_path, earlier_path, stator_path = paths.path("dxf_rotor_with_band"), paths.path("dxf_rotor"), paths.path("dxf_stator")
    rotor, earlier_rotor, stator = b2_gmsh.read_dxf(rotor_path), b2_gmsh.read_dxf(earlier_path), b2_gmsh.read_dxf(stator_path)
    exact = b2_gmsh.build_faces(rotor, stator, slot_size=slot)
    earlier = b2_gmsh.build_faces(earlier_rotor, stator, slot_size=slot)
    meshed = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor_path, b2_gmsh.MESH_SAGITTA_MM, b2_gmsh.MIN_SEGMENT_MM),
                                 b2_gmsh.read_dxf(stator_path, b2_gmsh.MESH_SAGITTA_MM, b2_gmsh.MIN_SEGMENT_MM), slot_size=slot)
    b2_gmsh.mesh_faces(meshed, folder / "path2.msh")
    info = b2_gmsh.inspect_mesh(folder / "path2.msh")
    chords = b2_gmsh.fillet_chords(b2_gmsh.mesh_nodes(folder / "path2.msh"), b2_gmsh.drawing_fillets(rotor, p))
    return SimpleNamespace(p=p, regions=regions, rotor=rotor, stator=stator, earlier_rotor=earlier_rotor, exact=exact, earlier=earlier, meshed=meshed, info=info,
                           chords=chords, offsets=b2_gmsh.cad_offsets(exact, p, regions), folder=folder)


def path1_checks_of(r):
    from bench import b2_gmsh
    return b2_gmsh.path1_checks(r.info, r.report, r.table, r.regions, r.p, r.chords)


def path2_checks_of(r, offsets=None):
    from bench import b2_gmsh
    return b2_gmsh.path2_checks(r.exact, r.meshed, r.info, r.offsets if offsets is None else offsets, r.p, r.regions, r.chords, r.rotor, r.stator, r.earlier,
                                r.earlier_rotor)


def test_the_nodes_of_a_mesh_file_are_read_from_its_text(tmp_path):
    from bench import b2_gmsh
    tiny_mesh(tmp_path / "tiny.msh")
    nodes = b2_gmsh.mesh_nodes(tmp_path / "tiny.msh")
    assert nodes.shape == (b2_gmsh.inspect_mesh(tmp_path / "tiny.msh")["nodes"], 2)
    assert nodes.min() >= 0.0 and nodes.max() <= 1.0 and {(0.0, 0.0), (1.0, 1.0)} <= {tuple(row) for row in nodes.round(12)}


def test_the_mesh_edges_on_an_arc_are_counted_from_the_nodes_that_lie_on_it():
    from bench import b2_gmsh
    on = lambda centre, r, deg: (centre[0] + r * math.cos(math.radians(deg)), centre[1] + r * math.sin(math.radians(deg)))
    arcs = [(1.0, 2.0, 0.3, 90.0, 90.0), (5.0, 5.0, 0.3, 350.0, 90.0)]                  # a quarter circle, and one that runs through 0 degrees
    nodes = [on((1, 2), 0.3, 90), on((1, 2), 0.3, 135), on((1, 2), 0.3, 180),               # 3 nodes on the first: 2 edges
             on((1, 2), 0.3, 0), on((1, 2), 0.3, 270), on((1, 2), 0.299, 120), on((1, 2), 0.35, 150),      # the same circle outside the sweep; near it, not on it
             on((5, 5), 0.3, 350), on((5, 5), 0.3, 20), on((5, 5), 0.3, 45), on((5, 5), 0.3, 80), on((5, 5), 0.3, 90)]   # 4 nodes on the second: 3 edges
    assert b2_gmsh.fillet_chords(np.array(nodes), arcs) == [2, 3]
    assert b2_gmsh.fillet_chords(np.array(nodes[3:7]), arcs) == [0, 0]                      # an arc with no node on it has none, not -1


def test_the_fillets_of_the_model_are_the_arcs_of_the_slot_fillet_radius_whatever_it_is():
    """M2 of the review of Task 11: the check had a typed 0.3. The radius is p.slot_fillet; the sweeps are the arcs' own."""
    from bench import b2_gmsh
    for fillet in (0.3, 0.25):
        p = mg.Params(barrier="none", slot_fillet=fillet)
        arcs = b2_gmsh.model_fillets(mg.build(p, gap="sliding"), p)
        assert len(arcs) == 16 and {round(a[2], 9) for a in arcs} == {fillet} and {round(a[4], 6) for a in arcs} == {90.0}      # 4 magnets x 4 fillets in the sector


def test_path_1_names_the_nodes_of_the_file_the_ones_to_gmsh_counted_and_the_ones_that_no_element_uses(real_path1):
    """M5 of the review of Task 11: motorgeo.to_gmsh returns the nodes of its Gmsh session, which include the centre of every arc (87 here: a point of the
    geometry, no element); the saved file has fewer (it is what ElmerGrid reads). Both are said, each with what it counts."""
    made, info = real_path1.made, real_path1.info
    assert set(made) >= {"nodes_in_file", "to_gmsh_nodes", "nodes_without_element", "boundary_curves", "size_max"} and "nodes" not in made
    assert made["nodes_in_file"] == info["nodes"]                                     # B2's inspection of the file and the function that wrote it agree
    assert made["to_gmsh_nodes"] == made["nodes_in_file"] + made["nodes_without_element"] and made["nodes_without_element"] > 0
    arc_edges = sum(1 for r in real_path1.regions if r.material != "<No Mesh>" for e in r.edges if isinstance(e, mg.Arc))
    assert made["nodes_without_element"] == arc_edges                                 # one centre point for each arc edge of the meshed regions


def test_the_fillet_check_counts_the_edges_on_the_mesh_and_holds_the_chord_error_to_1_percent(real_path1):
    """M2 of the review of Task 11: the radius is p.slot_fillet, the number of mesh edges on a fillet is counted on the mesh (and has to be the same on all
    16), and the window is 1 %, not 10 %: the real agreement is about 0.003 %."""
    import copy

    from bench import b2_gmsh
    r = real_path1
    assert r.chords == [2] * 16
    fillet = next(c for c in path1_checks_of(r) if "fillet" in c["name"])
    assert fillet["verdict"] == "PASS" and "radius 0.3 mm" in fillet["expected"] and "2 mesh edges" in fillet["expected"] and "within 1 %" in fillet["expected"]
    per_slot = 4 * b2_gmsh.chord_deficit(r.p.slot_fillet, 90.0, 2)
    for shift, verdict in ((0.005, "PASS"), (0.02, "FAIL")):                          # the mesh area off by 0.5 % and by 2 % of the predicted deficit
        report = copy.deepcopy(r.report)
        report["mesh_area_mm2"]["N45"] += shift * per_slot * 4
        got = next(c for c in b2_gmsh.path1_checks(r.info, report, r.table, r.regions, r.p, r.chords) if "fillet" in c["name"])
        assert got["verdict"] == verdict, (shift, got["measured"])
    wrong_count = next(c for c in b2_gmsh.path1_checks(r.info, r.report, r.table, r.regions, r.p, [3] * 16) if "fillet" in c["name"])
    assert wrong_count["verdict"] == "FAIL"                                           # a prediction made with 3 edges to a fillet does not fit the mesh
    mixed = next(c for c in b2_gmsh.path1_checks(r.info, r.report, r.table, r.regions, r.p, [2] * 15 + [3]) if "fillet" in c["name"])
    assert mixed["verdict"] == "FAIL" and "not the same" in mixed["measured"]


def test_the_air_band_check_says_what_was_found_and_removed_from_the_drawing_and_not_what_a_comment_says(tmp_path):
    """M4 of the review of Task 11: the text hard-coded "this drawing has no band (config.toml's comment ...)". That was true of Rotor.DXF and false of SynRM rotor.DXF,
    which has sixteen, and it went stale when the comment was corrected. The text is computed from bands_removed and names the drawing; the check is called what it asserts."""
    from bench import b2_gmsh
    p = mg.Params(barrier="none")
    for bands, wanted in ((True, 2), (False, 0)):
        folder = tmp_path / f"bands_{bands}"
        folder.mkdir()
        rotor, stator = write_synthetic_cad(folder, bands=bands)
        faces = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor), b2_gmsh.read_dxf(stator), slot_size=(4.0, 1.0))
        item = b2_gmsh.band_check(faces, p)
        assert item["name"] == "path 2: no air band remains in the rotor"
        assert f"{wanted} band(s) found and removed" in item["measured"] and "config.toml" not in item["measured"] + item["expected"] + item["source"]
        assert item["verdict"] == "FAIL" and "2 pockets left, all magnet slots (4 x 1 mm)" in item["measured"]      # (the synthetic drawing has 2 magnet slots, not the motor's 16)


def test_the_air_band_check_passes_on_the_real_drawing_and_says_all_16_pockets_are_magnet_slots(real_path2):
    item = next(c for c in path2_checks_of(real_path2) if c["name"] == "path 2: no air band remains in the rotor")
    assert item["verdict"] == "PASS" and "16 band(s) found and removed from SynRM rotor.DXF" in item["measured"] and "16 pockets left, all magnet slots" in item["measured"]
    assert "config.toml" not in item["measured"] + item["expected"] + item["source"]


def test_every_line_arc_and_circle_of_a_drawing_must_become_an_edge_of_a_face_and_what_is_ignored_is_listed(tmp_path):
    """M9 of the review of Task 11: read_dxf only counted entity types it does not use. The check: every LINE, ARC and CIRCLE is an edge of a face (or was
    dropped with an air band on purpose), and the rest (the text) is named."""
    import ezdxf

    from bench import b2_gmsh
    (tmp_path / "with").mkdir()
    rotor_path, stator_path = write_synthetic_cad(tmp_path / "with", bands=True)
    rotor, stator = b2_gmsh.read_dxf(rotor_path), b2_gmsh.read_dxf(stator_path)
    faces = b2_gmsh.build_faces(rotor, stator, slot_size=(4.0, 1.0))
    assert faces["band_lines_dropped"] == 8                                             # the two bands' four lines each
    item = b2_gmsh.entity_check(faces, rotor, stator)
    assert item["verdict"] == "PASS" and "MTEXT" in item["measured"] and "8 dropped with the 2 air band(s)" in item["measured"]
    assert "rotor 16 LINE + 2 CIRCLE" in item["measured"] and "stator 6 LINE + 6 ARC + 1 CIRCLE" in item["measured"]
    (tmp_path / "without").mkdir()
    rotor_path, stator_path = write_synthetic_cad(tmp_path / "without", bands=False)
    rotor, stator = b2_gmsh.read_dxf(rotor_path), b2_gmsh.read_dxf(stator_path)
    faces = b2_gmsh.build_faces(rotor, stator, slot_size=(4.0, 1.0))
    assert faces["band_lines_dropped"] == 0 and b2_gmsh.entity_check(faces, rotor, stator)["verdict"] == "PASS"
    doc = ezdxf.readfile(str(rotor_path))                                               # a line that belongs to no face: in the steel, a dangling edge
    doc.modelspace().add_line((-9.0, 0.0), (-9.5, 0.0))
    doc.saveas(str(rotor_path))
    rotor = b2_gmsh.read_dxf(rotor_path)
    stray = b2_gmsh.entity_check(b2_gmsh.build_faces(rotor, stator, slot_size=(4.0, 1.0)), rotor, stator)
    assert stray["verdict"] == "FAIL" and "1 not an edge of any face" in stray["measured"]


def test_rotor_dxf_and_the_stator_are_consumed_whole_and_only_the_text_is_ignored(real_path2):
    """Rotor.DXF has no band, so nothing is dropped: all 130 entities of the rotor are edges of a face (SynRM rotor.DXF's accounting is in the test further down)."""
    from bench import b2_gmsh
    r = real_path2
    item = b2_gmsh.entity_check(r.earlier, r.earlier_rotor, r.stator)
    assert item["verdict"] == "PASS" and "0 not an edge of any face" in item["measured"] and "0 dropped" in item["measured"]
    assert "rotor 64 LINE + 64 ARC + 2 CIRCLE (Rotor.DXF)" in item["measured"] and "stator 216 LINE + 72 ARC + 1 CIRCLE (Stator_refined.DXF)" in item["measured"]
    assert "ignored" in item["measured"] and "MTEXT" in item["measured"]


def test_path_2_holds_the_mesh_gmsh_wrote_against_motorgeo_at_the_briefs_0_5_percent(real_path2):
    """M7 of the review of Task 11: the areas were checked on the drawing's faces; the mesh Gmsh wrote (magnets -0.172 %, steel -0.005 % / +0.007 %, slots
    -0.048 %) was neither shown nor gated. Now it is, the magnets' shortfall is explained by the fillets' chords as in path 1, and the stator's slot residual
    is a check of its own."""
    checks = path2_checks_of(real_path2)
    mesh = [c for c in checks if "the mesh Gmsh wrote against motorgeo's" in c["name"]]
    assert len(mesh) == 4 and all(c["tolerance_percent"] == 0.5 and c["verdict"] == "PASS" for c in mesh)
    magnets = next(c for c in mesh if "magnet" in c["name"])
    assert magnets["difference_percent"] == pytest.approx(-0.172, abs=0.002)
    chord = next(c for c in checks if "chord error" in c["name"])
    assert chord["verdict"] == "PASS" and "64 fillets" in chord["expected"] and "2 mesh edges" in chord["expected"] and "within 1 %" in chord["expected"]
    stator = next(c for c in checks if c["name"] == "path 2: the drawing's stator slots are motorgeo's, turned")
    assert stator["verdict"] == "PASS"
    rotor = next(c for c in checks if c["name"] == "path 2: the drawing's magnet slots are motorgeo's, turned")
    assert rotor["verdict"] == "PASS" and "stator" not in rotor["expected"]


def test_the_stator_slot_residual_decides_a_verdict_in_millimetres_at_the_bore(real_path2):
    import copy

    r = real_path2
    offsets = copy.deepcopy(r.offsets)
    bore = r.exact["radii"]["bore"]
    verdicts = {}
    for label, mm in (("wide", 0.1), ("narrow", 0.04)):                               # arc length at the bore: twice the 0.05 mm, and under it
        offsets["stator"]["max_distance_deg"] = math.degrees(mm / bore)
        checks = path2_checks_of(r, offsets)
        verdicts[label] = (next(c for c in checks if "stator slots are motorgeo's" in c["name"])["verdict"],
                           next(c for c in checks if "magnet slots are motorgeo's" in c["name"])["verdict"])
    assert verdicts == {"wide": ("FAIL", "PASS"), "narrow": ("PASS", "PASS")}        # and the rotor's own number is another check


def test_every_limit_that_b2_chose_says_so_and_the_briefs_do_not(real_path1, real_path2):
    """M1 of the review of Task 11: path 1's and path 2's SICN floors, the 0.05 mm of position, the chord check's window are the bench's own; the 0.1 % of
    path 1 and the 0.5 % of the DXF path are the brief's."""
    import bench
    own = "(" + bench.OWN_LIMIT + ")"
    one, two = path1_checks_of(real_path1), path2_checks_of(real_path2)
    assert own in next(c for c in one if c["name"].startswith("path 1: element quality"))["expected"]
    assert own in next(c for c in two if c["name"].startswith("path 2: element quality"))["expected"]
    assert own in next(c for c in one if "fillet" in c["name"])["expected"] and own in next(c for c in two if "chord error" in c["name"])["expected"]
    assert own in next(c for c in two if c["name"].startswith("path 2: the meshed areas are the faces' areas"))["expected"]         # P4 of the re-review: the 1e-6
    positions = [c for c in two if "are motorgeo's, turned" in c["name"]]
    assert len(positions) == 2 and all(own in c["expected"] for c in positions)
    for c in one + two:
        if "tolerance_percent" in c:
            assert bench.OWN_LIMIT not in c["source"], c["name"]                        # the brief's tolerances are not labelled as the bench's


def test_the_briefs_area_tolerances_are_the_ones_the_area_checks_use(real_path1, real_path2):
    """M10 of the review of Task 11: B1's 0.5 % was pinned to the brief; these two were not. Path 1: each material's area within 0.1 %; the DXF path: the steel and
    magnet areas within 0.5 %."""
    from bench import b2_gmsh
    assert b2_gmsh.AREA_TOLERANCE_PERCENT == 0.1 and b2_gmsh.CAD_AREA_TOLERANCE_PERCENT == 0.5
    one = [c for c in path1_checks_of(real_path1) if "tolerance_percent" in c]
    assert len(one) == 4 and {c["tolerance_percent"] for c in one} == {0.1}
    two = [c for c in path2_checks_of(real_path2) if "tolerance_percent" in c]
    assert len(two) == 8 and {c["tolerance_percent"] for c in two} == {0.5}


def test_the_figure_of_the_two_renders_has_no_blank_band_under_or_beside_them(tmp_path):
    """M11 of the review of Task 11: the pictures go into the project notes, so none may carry a blank margin. The two renders are cropped to their ink, put side by side
    and the figure is saved with a tight box, so what is left around them is the caption above and a small even margin."""
    import matplotlib.image as image

    from bench import b2_gmsh
    yy, xx = np.mgrid[0:700, 0:700]
    disc = np.ones((700, 700, 3))
    disc[(yy - 350) ** 2 + (xx - 350) ** 2 < 300 ** 2] = (0.1, 0.5, 0.2)                    # a model in the middle of Gmsh's white, square window
    quarter = np.ones((700, 700, 3))
    quarter[(yy - 20) ** 2 + (xx - 20) ** 2 < 650 ** 2] = (0.9, 0.8, 0.1)                   # and a quarter of one, in a corner
    image.imsave(str(tmp_path / "left.png"), quarter)
    image.imsave(str(tmp_path / "right.png"), disc)
    b2_gmsh.figure(tmp_path, "left.png", "right.png", ("Path 1: the sector\n8,819 nodes", "Path 2: the whole machine\n34,415 nodes"), tmp_path / "figure.png")
    picture = image.imread(str(tmp_path / "figure.png"))
    ink = (picture[..., :3] < 0.97)
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    height, width = picture.shape[:2]
    assert (height - 1 - rows[-1]) < 0.04 * height, f"{height - 1 - rows[-1]} blank rows under the picture of {height}"
    assert rows[0] < 0.04 * height and cols[0] < 0.04 * width and (width - 1 - cols[-1]) < 0.04 * width


# ----------------------------------------------------------------------------------------------- 3. SynRM rotor.DXF: the air band is found and removed in memory (Task 11, fix round 2)
# The brief named Rotor.DXF, which has no band; the drawing that has it is SynRM rotor.DXF, the CAD's later export (config.toml: dxf_rotor_with_band). A design decision of
# 2026-11-05: the band is removed in Python, not in the file. Path 2 reads SynRM rotor.DXF; Rotor.DXF stays as the independent check of the removal.
def write_shared_edge_cad(tmp_path, overlap=False):
    """A rotor (a disk of radius 10, a shaft hole of radius 2) whose first magnet slot (4 x 1) has an air band directly under it, touching it. overlap=False: the
    band is a 4 x 2 rectangle and the slot's bottom wall is ITS top wall too, one LINE serving both. overlap=True: the band is 6 x 1.5 with a long top edge of its own
    that lies along the slot's bottom wall (two LINEs on one another, as two exports of the same wall can be). A second slot, far away, has no band. The stator is
    write_synthetic_cad's. Returns (rotor path, stator path)."""
    import ezdxf
    _, stator_path = write_synthetic_cad(tmp_path, bands=False)
    doc = ezdxf.new("R2000")
    msp = doc.modelspace()
    msp.add_circle((0, 0), 2.0)
    msp.add_circle((0, 0), 10.0)
    lines = [((4, 3.5), (8, 3.5)), ((4, 2.5), (4, 3.5)), ((8, 2.5), (8, 3.5)), ((4, 2.5), (8, 2.5))]                   # the slot: its top, its two sides, its bottom wall
    if overlap:
        lines += [((3, 2.5), (9, 2.5)), ((3, 1), (3, 2.5)), ((9, 1), (9, 2.5)), ((3, 1), (9, 1))]                      # the band, with its own long top edge
    else:
        lines += [((4, 0.5), (4, 2.5)), ((8, 0.5), (8, 2.5)), ((4, 0.5), (8, 0.5))]                                    # the band, under the slot: three walls of its own
    lines += [((-8, -3.5), (-4, -3.5)), ((-4, -3.5), (-4, -2.5)), ((-4, -2.5), (-8, -2.5)), ((-8, -2.5), (-8, -3.5))]  # the second slot
    for a, b in lines:
        msp.add_line(a, b)
    msp.add_mtext("SOLIDWORKS Educational Product.").dxf.insert = (-10, -12)
    doc.saveas(tmp_path / "rotor_shared.dxf")
    return tmp_path / "rotor_shared.dxf", stator_path


def test_a_band_that_shares_an_edge_with_a_slot_keeps_the_shared_edge_and_the_slot_stays_closed(tmp_path):
    """F5 of fix round 2: build_faces drops the edges of a band that no magnet slot shares. The wall a band and a slot have in common is the slot's: it stays, and the slot
    stays a closed face; the steel closes over the band all the same. Both ways a CAD can draw it: one LINE serving both, or two lying along each other."""
    from bench import b2_gmsh
    for overlap, dropped in ((False, 3), (True, 4)):
        folder = tmp_path / f"overlap_{overlap}"
        folder.mkdir()
        rotor_path, stator_path = write_shared_edge_cad(folder, overlap)
        rotor, stator = b2_gmsh.read_dxf(rotor_path), b2_gmsh.read_dxf(stator_path)
        faces = b2_gmsh.build_faces(rotor, stator, slot_size=(4.0, 1.0))
        assert faces["bands_removed"] == 1 and len(faces["magnets"]) == 2, overlap
        assert all(abs(m.area - 4.0) < 1e-9 for m in faces["magnets"]), overlap                           # the slot is still a closed 4 x 1 face
        assert faces["band_lines_dropped"] == dropped and faces["band_edges_shared_with_slots"] == 1, overlap          # the band's own walls go, the shared wall stays
        assert faces["rotor_steel"].area == pytest.approx(math.pi * 96 - 8.0, rel=2e-4), overlap                       # the steel closes over the band
        assert faces["rotor_steel"].area - faces["steel_before_removal_mm2"] == pytest.approx(faces["bands"][0].area, abs=1e-9), overlap
        item = b2_gmsh.entity_check(faces, rotor, stator)
        assert item["verdict"] == "PASS" and "1 band edge(s) shared with a magnet slot kept" in item["measured"], (overlap, item["measured"])


def test_read_dxf_says_which_drawing_and_which_kind_of_entity_each_polyline_came_from(tmp_path):
    from collections import Counter

    from bench import b2_gmsh
    rotor_path, stator_path = write_synthetic_cad(tmp_path)
    rotor = b2_gmsh.read_dxf(rotor_path)
    assert rotor["file"] == "rotor.dxf" and rotor["path"] == str(rotor_path)
    assert len(rotor["kinds"]) == len(rotor["lines"]) == 18 and Counter(rotor["kinds"]) == {"CIRCLE": 2, "LINE": 16}
    for line, kind in zip(rotor["lines"], rotor["kinds"]):
        assert (len(line.coords) == 2) == (kind == "LINE") and (math.dist(line.coords[0], line.coords[-1]) < 1e-9 if kind == "CIRCLE" else True)      # a circle comes back to its start
    faces = b2_gmsh.build_faces(rotor, b2_gmsh.read_dxf(stator_path), slot_size=(4.0, 1.0))
    assert faces["files"] == {"rotor": "rotor.dxf", "stator": "stator.dxf"}


def test_motorgeo_says_how_many_air_bands_the_machine_has_and_what_one_of_them_is():
    """F2 of fix round 2: the number of bands to expect is not a 16 typed in the bench. It is the band outlines of Params(barrier="air") in the sector (4) times the sectors
    of the machine (4), and one band is what motorgeo draws: 5 LINE and 1 ARC, 16.705 x 2.531 mm, 23.99 mm2, its centroid at radius 37.9 mm, 1.46 mm from its slot."""
    from bench import b2_gmsh
    p = mg.Params(barrier="none")
    model = b2_gmsh.motorgeo_bands(p)
    assert (model["per_sector"], model["sectors"], model["count"]) == (4, 4, 16) and model["sector_deg"] == 90.0
    assert model["outline"] == {"LINE": 5, "ARC": 1} and model["edges"] == {"LINE": 80, "ARC": 16}
    assert len(model["bands"]) == 16 and len(model["slots"]) == 16
    one = model["bands"][0]
    assert one.area == pytest.approx(23.99, abs=0.01)                                  # 4 x 23.99 = 95.97 mm2 in the sector: the band B1 found between "none" and "steel"
    assert b2_gmsh.rectangle_sides(one) == pytest.approx((16.705, 2.531), abs=0.001)
    assert math.hypot(one.centroid.x, one.centroid.y) == pytest.approx(37.899, abs=0.001)
    assert min(one.distance(s) for s in model["slots"]) == pytest.approx(1.460, abs=0.001)         # a web of steel between a band and its slot
    assert b2_gmsh.motorgeo_bands(mg.Params(barrier="air"))["count"] == 16 and b2_gmsh.motorgeo_bands(mg.Params(barrier="steel"))["count"] == 16


def test_the_band_removal_check_holds_each_removed_band_against_motorgeos_band(real_path2):
    import bench
    item = next(c for c in path2_checks_of(real_path2) if c["name"].startswith("path 2: the air bands of"))
    assert item["name"] == "path 2: the air bands of SynRM rotor.DXF are found and removed, and each is motorgeo's band"
    assert item["verdict"] == "PASS" and item["gating"] is True
    measured = item["measured"]
    assert "16 found and removed from SynRM rotor.DXF" in measured and "4 band outlines in the 90-degree sector x 4 sectors = 16" in measured
    assert "96 edges dropped = 80 LINE + 16 ARC" in measured and "5 LINE + 1 ARC" in measured
    assert "(" + bench.OWN_LIMIT + ")" in item["expected"]


def test_the_band_removal_check_fails_on_a_band_that_is_missing_displaced_or_the_wrong_size_and_on_edges_that_are_not_a_band_s(real_path2):
    from shapely import affinity

    from bench import b2_gmsh
    r = real_path2
    turn = r.offsets["rotor"]["turn_deg_mod_45"]
    assert b2_gmsh.band_removal_check(r.exact, r.p, turn)["verdict"] == "PASS"
    bands = r.exact["bands"]
    cases = {"one band missing": dict(r.exact, bands=bands[:-1], bands_removed=15),
             "a band displaced by 0.2 mm": dict(r.exact, bands=[affinity.translate(bands[0], 0.2, 0.0), *bands[1:]]),
             "a band 3 % too big": dict(r.exact, bands=[affinity.scale(bands[0], 1.03, 1.03), *bands[1:]]),
             "a magnet slot taken for a band": dict(r.exact, bands=[r.exact["magnets"][0], *bands[1:]]),
             "dropped edges that are not a band's outline": dict(r.exact, band_edge_kinds={"LINE": 80, "ARC": 15})}
    for what, faces in cases.items():
        assert b2_gmsh.band_removal_check(faces, r.p, turn)["verdict"] == "FAIL", what
    assert b2_gmsh.band_removal_check(r.exact, r.p, turn + 5.0)["verdict"] == "FAIL"             # the bands sit where the slots' turn puts them, not anywhere
    assert b2_gmsh.band_removal_check(dict(r.exact, bands=[], bands_removed=0), r.p, turn)["verdict"] == "FAIL"        # none found: said, not a crash


def test_after_the_removal_the_band_drawing_is_the_earlier_export_rotor_dxf(real_path2):
    """F3 of fix round 2: the independent check. Rotor.DXF (the CAD's export of 2026-01-13) never had the band. Once build_faces has removed the bands of SynRM rotor.DXF,
    the two rotors are the same: not just close, the same faces, so the limits below are the noise floor and nothing near them is used up."""
    import bench
    from bench import b2_gmsh
    r = real_path2
    item = next(c for c in path2_checks_of(r) if c["name"] == "path 2: after the removal, SynRM rotor.DXF is Rotor.DXF")
    assert item["verdict"] == "PASS" and item["gating"] is True
    assert "(" + bench.OWN_LIMIT + ")" in item["expected"] and "sagitta" in item["expected"] and "0 band(s) found in Rotor.DXF" in item["measured"]
    assert r.earlier["bands_removed"] == 0 and r.earlier["rotor_steel"].area == pytest.approx(r.exact["rotor_steel"].area, abs=1e-6)
    assert r.exact["rotor_steel"].symmetric_difference(r.earlier["rotor_steel"]).area < 1e-9
    assert all(a.symmetric_difference(b).area < 1e-9 for a, b in zip(r.exact["magnets"], r.earlier["magnets"]))
    limit = b2_gmsh.flattening_area_limit(r.earlier_rotor)
    assert limit == pytest.approx(2 / 3 * b2_gmsh.SAGITTA_MM * 360.03, rel=1e-4) and limit == pytest.approx(0.0600, abs=1e-4)          # (2/3) x the sagitta x the arc length
    assert b2_gmsh.flattening_area_limit(dict(r.earlier_rotor, sagitta=2 * b2_gmsh.SAGITTA_MM)) == pytest.approx(2 * limit)


def test_the_comparison_with_the_earlier_export_fails_on_a_different_steel_slot_circle_or_balance(real_path2):
    from shapely import affinity

    from bench import b2_gmsh
    r = real_path2
    make = lambda **kw: b2_gmsh.earlier_export_check(kw.get("faces", r.exact), kw.get("earlier", r.earlier), kw.get("rotor", r.rotor), kw.get("earlier_rotor", r.earlier_rotor))
    assert make()["verdict"] == "PASS"
    assert make(earlier=dict(r.earlier, rotor_steel=r.earlier["rotor_steel"].buffer(0.01)))["verdict"] == "FAIL"                   # 0.01 mm of steel all round
    assert make(earlier=dict(r.earlier, magnets=[affinity.translate(m, 0.1, 0.0) for m in r.earlier["magnets"]]))["verdict"] == "FAIL"        # every slot 0.1 mm away
    assert make(earlier=dict(r.earlier, magnets=r.earlier["magnets"][:-1]))["verdict"] == "FAIL"                                   # a slot missing
    circles = [(cx, cy, radius + 0.1) for cx, cy, radius in r.earlier_rotor["circles"]]
    assert make(earlier_rotor=dict(r.earlier_rotor, circles=circles))["verdict"] == "FAIL"                                          # the circles 0.1 mm bigger
    assert make(faces=dict(r.exact, steel_before_removal_mm2=r.exact["steel_before_removal_mm2"] + 1.0))["verdict"] == "FAIL"     # the steel added is not the bands' area


def test_the_entity_check_counts_the_dropped_edges_and_shows_they_are_exactly_the_bands_edges(real_path2):
    """F2 of fix round 2: every LINE, ARC and CIRCLE of SynRM rotor.DXF is an edge of a face, except the 96 edges of the 16 bands that were dropped on purpose; those
    are exactly the bands' outlines: all of them lie on a removed band's outline, none on a magnet slot's, none is left over (it would dangle in the steel)."""
    from bench import b2_gmsh
    r = real_path2
    item = next(c for c in path2_checks_of(r) if c["name"].startswith("path 2: every LINE, ARC and CIRCLE"))
    measured = item["measured"]
    assert item["verdict"] == "PASS"
    assert "rotor 144 LINE + 80 ARC + 2 CIRCLE (SynRM rotor.DXF)" in measured and "stator 216 LINE + 72 ARC + 1 CIRCLE (Stator_refined.DXF)" in measured
    assert "96 dropped with the 16 air band(s) (80 LINE + 16 ARC)" in measured and "0 not an edge of any face" in measured
    assert "96 lie on a removed band's outline, 0 on a magnet slot's" in measured and "0 band edge(s) shared with a magnet slot kept" in measured
    assert "ignored" in measured and "MTEXT" in measured
    kept = r.exact["band_line_indices"]
    left_over = b2_gmsh.entity_check(dict(r.exact, band_line_indices=kept[:-1], band_lines_dropped=95), r.rotor, r.stator)       # one band edge not dropped: it dangles
    assert left_over["verdict"] == "FAIL" and "1 not an edge of any face" in left_over["measured"]
    zone = r.exact["magnets"][0].boundary.buffer(b2_gmsh.EDGE_ZONE_MM)
    slot_edge = next(i for i, line in enumerate(r.rotor["lines"]) if zone.covers(line))                                     # an edge of the first magnet slot
    too_many = b2_gmsh.entity_check(dict(r.exact, band_line_indices=sorted([*kept, slot_edge]), band_lines_dropped=97), r.rotor, r.stator)
    assert too_many["verdict"] == "FAIL" and "1 on a magnet slot's" in too_many["measured"]


def test_the_three_drawings_are_hashed_and_a_change_between_two_hashings_is_a_failed_check(tmp_path):
    """F6 of fix round 2: the drawings in reference\\ are read only. B2 hashes the three before it reads them and again at its end; a script that wrote into one would be a
    FAILED check, not a surprise."""
    import hashlib

    from bench import b2_gmsh
    a, b = tmp_path / "a.dxf", tmp_path / "b.dxf"
    a.write_bytes(b"alpha")
    b.write_bytes(b"beta")
    before = b2_gmsh.file_hashes([a, b])
    assert before == {"a.dxf": hashlib.sha256(b"alpha").hexdigest(), "b.dxf": hashlib.sha256(b"beta").hexdigest()}
    same = b2_gmsh.read_only_check(before, b2_gmsh.file_hashes([a, b]))
    assert same["verdict"] == "PASS" and same["gating"] is True and "unchanged" in same["measured"]
    b.write_bytes(b"beta, written to")
    changed = b2_gmsh.read_only_check(before, b2_gmsh.file_hashes([a, b]))
    assert changed["verdict"] == "FAIL" and "b.dxf" in changed["measured"] and "changed" in changed["measured"]
    assert b2_gmsh.read_only_check(before, {"a.dxf": before["a.dxf"]})["verdict"] == "FAIL"           # a file that is gone is not "unchanged"


def test_the_three_reference_drawings_are_the_files_the_controller_committed():
    """The copies in reference\\ are byte for byte what the CAD exported (stored -text, so no line-ending change can touch them): SynRM rotor.DXF (2026-01-22, the one with
    the band) 4e794cab..., Rotor.DXF (2026-01-13, no band) 388e567b..., Stator_refined.DXF (2026-01-13) 560262fe... If one of these fails, a drawing was changed or replaced."""
    from bench import b2_gmsh
    hashes = b2_gmsh.file_hashes([paths.path("dxf_rotor_with_band"), paths.path("dxf_rotor"), paths.path("dxf_stator")])
    assert hashes["SynRM rotor.DXF"] == "4e794cabc1fd0a3c9bfa1f5cec68c8e66066329f734613091a10e088cf3d63e1"
    assert hashes["Rotor.DXF"].startswith("388e567b2127cd09") and hashes["Stator_refined.DXF"].startswith("560262fe20964880")


@pytest.fixture(scope="module")
def removal_off(tmp_path_factory, real_path2):
    """SynRM rotor.DXF read with the removal SWITCHED OFF: _is_magnet_slot says yes to every pocket, so nothing is taken for a band and the sixteen band pockets stay in the
    rotor, as 'magnet slots' (32 of them), in the faces, the mesh and everything made from them."""
    from types import SimpleNamespace

    from bench import b2_gmsh
    r = real_path2
    folder = tmp_path_factory.mktemp("path2_off")
    slot = (r.p.slot_w, r.p.slot_t)
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(b2_gmsh, "_is_magnet_slot", lambda *args, **kwargs: True)
        exact = b2_gmsh.build_faces(r.rotor, r.stator, slot_size=slot)
        meshed = b2_gmsh.build_faces(b2_gmsh.read_dxf(paths.path("dxf_rotor_with_band"), b2_gmsh.MESH_SAGITTA_MM, b2_gmsh.MIN_SEGMENT_MM),
                                     b2_gmsh.read_dxf(paths.path("dxf_stator"), b2_gmsh.MESH_SAGITTA_MM, b2_gmsh.MIN_SEGMENT_MM), slot_size=slot)
    b2_gmsh.mesh_faces(meshed, folder / "path2.msh")
    chords = b2_gmsh.fillet_chords(b2_gmsh.mesh_nodes(folder / "path2.msh"), b2_gmsh.drawing_fillets(r.rotor, r.p))
    return SimpleNamespace(**{**vars(r), "exact": exact, "meshed": meshed, "info": b2_gmsh.inspect_mesh(folder / "path2.msh"), "chords": chords,
                              "offsets": b2_gmsh.cad_offsets(exact, r.p, r.regions)})


def test_with_the_removal_switched_off_the_band_check_the_comparison_with_rotor_dxf_and_the_area_checks_all_fail(real_path2, removal_off):
    """F5 of fix round 2, the RED made permanent: with the removal off nothing is found, 32 pockets stay in the rotor, and the checks that stand for the removal say so.
    With it on (real_path2) every one of them passes."""
    assert removal_off.exact["bands_removed"] == 0 and len(removal_off.exact["magnets"]) == 32
    on, off = path2_checks_of(real_path2), path2_checks_of(removal_off)
    for start in ("path 2: the air bands of", "path 2: no air band remains", "path 2: after the removal", "path 2: rotor steel area, the CAD drawing",
                  "path 2: magnet slots (N45) area, the CAD drawing", "path 2: rotor steel area, the mesh Gmsh wrote", "path 2: magnet slots (N45) area, the mesh Gmsh wrote"):
        assert next(c for c in on if c["name"].startswith(start))["verdict"] == "PASS", start
        assert next(c for c in off if c["name"].startswith(start))["verdict"] == "FAIL", start
    assert next(c for c in off if c["name"].startswith("path 2: no air band remains"))["measured"].count("only 16 of them magnet slots") == 1


def test_the_picture_of_the_removal_shows_the_steel_the_bands_and_the_slots_and_has_no_blank_margin(tmp_path):
    """F4 of fix round 2: bands.png, the band drawing's rotor faces from Python, the removed bands in a second colour and the magnet slots in a third."""
    import matplotlib.colors as colors
    import matplotlib.image as image

    from bench import b2_gmsh
    rotor_path, stator_path = write_synthetic_cad(tmp_path, bands=True)
    faces = b2_gmsh.build_faces(b2_gmsh.read_dxf(rotor_path), b2_gmsh.read_dxf(stator_path), slot_size=(4.0, 1.0))
    b2_gmsh.fig_bands(faces, "rotor.dxf, read by Python: the rotor's faces\n2 air-band pockets found and removed in memory (orange)\n2 magnet slots stay (blue)", tmp_path / "bands.png")
    picture = image.imread(str(tmp_path / "bands.png"))
    rgb = picture[..., :3]
    near = lambda colour: int((np.abs(rgb - np.array(colors.to_rgb(colour))).max(axis=2) < 0.03).sum())
    assert near(b2_gmsh.BAND_COLOURS["band"]) > 5000 and near(b2_gmsh.BAND_COLOURS["slot"]) > 3000 and near(b2_gmsh.BAND_COLOURS["steel"]) > 50000
    ink = (rgb < 0.97)
    rows, cols = np.where(ink.any(axis=2).any(axis=1))[0], np.where(ink.any(axis=2).any(axis=0))[0]
    height, width = picture.shape[:2]
    assert (height - 1 - rows[-1]) < 0.04 * height and rows[0] < 0.04 * height and cols[0] < 0.04 * width and (width - 1 - cols[-1]) < 0.04 * width
    assert width <= 1.1 * 9 * 150                                                                     # a long title does not widen it (bench.wrap_title)


def test_the_figure_the_log_and_the_record_name_the_drawing_that_was_read_and_what_was_removed(real_path2):
    """F4 of fix round 2: wherever the case says which drawing it read, it says SynRM rotor.DXF, how many bands were found and removed, and that Rotor.DXF agrees."""
    from bench import b2_gmsh
    r = real_path2
    caption = b2_gmsh.path2_caption(r.info, r.exact, r.rotor, r.stator, r.earlier, r.earlier_rotor)
    assert "SynRM rotor.DXF" in caption and "Stator_refined.DXF" in caption and "16 air bands found and removed in memory" in caption
    assert "agrees with Rotor.DXF" in caption and "0 bands found there" in caption                        # T11-P7: and that the CAD's earlier export agrees
    assert all(len(line) <= 80 for line in caption.split("\n"))                                         # it fits over half of the picture
    steel_differs = dict(r.earlier, rotor_steel=r.earlier["rotor_steel"].buffer(0.01))                    # an earlier export 0.01 mm bigger all round: the caption must not claim agreement
    differs = b2_gmsh.path2_caption(r.info, r.exact, r.rotor, r.stator, steel_differs, r.earlier_rotor)
    assert "DIFFERS from Rotor.DXF" in differs and "agrees" not in differs and all(len(line) <= 80 for line in differs.split("\n"))
    title = b2_gmsh.bands_title(r.exact, r.earlier, r.rotor, r.earlier_rotor)                            # the title of bands.png
    assert "SynRM rotor.DXF" in title and "16 air-band pockets found and removed in memory (orange)" in title and "16 magnet slots stay (blue)" in title
    assert "Rotor.DXF has 5,923.40" in title and "5,539.54 -> 5,923.40" in title and all(len(line) <= 92 for line in title.split("\n"))
    line = b2_gmsh.path2_log_line(r.info, r.exact, r.earlier, r.earlier_rotor, r.rotor, r.stator, 3.2)
    assert "SynRM rotor.DXF" in line and "16 air band(s) found and removed" in line and "Rotor.DXF" in line and "agrees" in line
    files = b2_gmsh.path2_files(r.rotor, r.earlier_rotor, r.stator, {"SynRM rotor.DXF": "aa", "Rotor.DXF": "bb", "Stator_refined.DXF": "cc"})
    assert files["rotor"].endswith("SynRM rotor.DXF") and files["rotor_earlier_export"].endswith("Rotor.DXF") and files["stator"].endswith("Stator_refined.DXF")
    assert files["names"] == {"rotor": "SynRM rotor.DXF", "rotor_earlier_export": "Rotor.DXF", "stator": "Stator_refined.DXF"} and files["sha256"]["Rotor.DXF"] == "bb"


def test_the_module_docstring_names_the_drawing_path_2_reads_and_the_one_that_checks_the_removal():
    from bench import b2_gmsh
    doc = b2_gmsh.__doc__
    assert "SynRM rotor.DXF" in doc and "Rotor.DXF" in doc and "bands.png" in doc and "16 air bands" in doc
    assert "which is why the worst triangle is in the gap." in doc                                      # P2 of the re-review: the full stop before the next sentence
    assert "the 32 pockets of SynRM rotor.DXF" not in doc and "the 16 magnet slots of SynRM rotor.DXF" in doc       # T11-P9: the final faces have no band pocket left
    assert "(32 in SynRM rotor.DXF: 16 magnet slots and 16 air bands)" in doc and "before path 2 reads them" in doc
    assert "path 1 with its table of mesh sizes 6.7 s" in doc and "6.9 s" not in doc
    assert "NO band (0 found, 0 removed)" not in doc                                                   # the old words, which named no drawing


def test_the_edge_zone_is_the_benchs_own_limit_and_the_check_says_so(real_path2):
    """T11-P12. EDGE_ZONE_MM (1e-4 mm), how close to a face's boundary an entity has to lie to count as an edge of it, is the bench's own choice and the brief gives none, so the
    check that uses it says so with the OWN_LIMIT label and names the number, like every other limit of B2."""
    from bench import b2_gmsh
    r = real_path2
    item = b2_gmsh.entity_check(r.earlier, r.earlier_rotor, r.stator)
    assert b2_gmsh.OWN_LIMIT in item["expected"] and f"{b2_gmsh.EDGE_ZONE_MM:g} mm" in item["expected"], item["expected"]
