"""Tests for B4 (bench/b4_ccx.py): Gmsh -> CalculiX. Everything that needs no solver: the hand estimate, the mesh, the deck, the result readers.

Run from P:   uv run pytest tests/test_bench_b4.py -q      (a few seconds; Gmsh runs in this process, no ccx, no cgx, no solver slot is taken)

In groups:
  1. The inputs: the note's materials, the speeds, the brief's 1 % pin, the way a number is written into a deck (CalculiX reads 20 characters of a field).
  2. The hand estimate of 04 Structural 6.1, written as a function: it must give the note's own numbers for the note's own rotor (with the air band) before its
     answer for the band-free rotor means anything.
  3. The mesh (Gmsh, second-order triangles): the areas against motorgeo's, the orientation, the shared slot walls.
  4. The magnets as a load: the slot's own axes, the force on the walls the magnet bears on, the nodal forces; and the deck of both models.
  5. Reading the answer: reactions from the .dat file, von Mises, the place of the peak, the force through a cut, the speed ratio.
  6. Fix round 1: the numbers of the real decks (at most 20 characters and 12 digits, for all three variants), the soft-filler diagnostic (its deck, its decision rule, its two
     reported checks), the peak of the slot bodies, and the words of the two pictures.
"""
import math
import re

import numpy as np
import pytest

import motorgeo as mg
from conftest import process_path

NOTE = "04-structural.md, 6.1 and 6.4"


def band_free():
    return mg.Params(barrier="none")


def with_band():
    return mg.Params()                  # the rotor the note's hand estimate was made for: barrier="air"


# ----------------------------------------------------------------------------------------------- 1. the inputs
def test_the_materials_are_the_note_s_values_and_the_ring_check_s_steel():
    """6.4 of the note: E 200 GPa, nu 0.3 and 7,650 kg/m3 for M-19 are CalculiX's inputs; N45 is 160 GPa, 0.24 and 7,600 kg/m3 (6.1, and the formula of the bridges)."""
    from bench import b4_ccx as b4
    from checks import ccx_check
    assert (b4.M19.youngs_mpa, b4.M19.poisson, b4.M19.density_kg_m3) == (200e3, 0.3, 7650.0)
    assert (b4.N45.youngs_mpa, b4.N45.poisson, b4.N45.density_kg_m3) == (160e3, 0.24, 7600.0)
    assert b4.M19.density_tonne_mm3 == pytest.approx(7.65e-9, rel=1e-12)                 # mm, N, tonne, s: 1 kg/m3 = 1e-12 tonne/mm3
    assert (ccx_check.E, ccx_check.NU, ccx_check.RHO) == (b4.M19.youngs_mpa, b4.M19.poisson, b4.M19.density_tonne_mm3)      # Task 4's ring is the same steel


def test_the_speeds_and_the_briefs_one_percent_are_pinned():
    """The brief: 10,000 and 12,000 rpm, and the two stresses scale by (12/10)^2 = 1.44 within 1 %. Change the constant and this goes red."""
    from bench import b4_ccx as b4
    assert b4.SPEEDS_RPM == (10_000, 12_000)
    assert b4.SCALE_TOLERANCE_PERCENT == 1.0
    assert b4.omega(10_000) == pytest.approx(1047.1975512, rel=1e-9)
    item = b4.scale_check("pressure", 30.0, 30.0 * 1.44)
    assert item["verdict"] == "PASS" and item["tolerance_percent"] == 1.0 and item["gating"] is True
    assert b4.scale_check("pressure", 30.0, 30.0 * 1.44 * 1.0099)["verdict"] == "PASS"        # inside 1 %
    assert b4.scale_check("pressure", 30.0, 30.0 * 1.44 * 1.0101)["verdict"] == "FAIL"        # outside it
    assert b4.scale_check("pressure", 30.0, 30.0 * 1.44 * 0.9901)["verdict"] == "PASS"
    assert b4.scale_check("pressure", 30.0, 30.0 * 1.44 * 0.9901, gating=False, reason="the verdict is the pressure model's")["gating"] is False


def test_a_number_in_a_deck_is_never_longer_than_the_twenty_characters_calculix_reads():
    """Found on 2026-11-06 with the first bonded deck: a density written as 6.627128124470999e-09 (21 characters) was cut to 6.627e-0 by ccx, a billion times too
    heavy, and the solve ended with 5.8e8 MPa and no warning. Every number a deck holds goes through deck_number."""
    from bench import b4_ccx as b4
    rho = 7.6e-9 * (9.08 * 1.57) / 16.34834244413349
    assert len(repr(rho)) > 20                                          # the number that broke the first deck
    for value in (rho, -rho, 1.096622711e6, 200000.0, 0.3, -42.806123456789, 1e-300, -1.2345678901234e-123):
        text = b4.deck_number(value)
        assert len(text) <= 20, text
        assert float(text) == pytest.approx(value, rel=1e-11)


# ----------------------------------------------------------------------------------------------- 2. the hand estimate
def test_the_hand_estimate_gives_the_note_s_numbers_for_the_note_s_own_rotor():
    """Self-check: 04 Structural 6.1 works the bridges out for the rotor WITH its air band: the pole's outer steel 162 mm2 at 40.8 mm, its two magnets 28.5 mm2 at
    39.9 mm, 29 g per pole over 20 mm, a pull of 65 kN/m (1.3 kN over 20 mm) at 10,000 rpm, the three strips 1.4 + 1.4 + 1.5 mm, a nominal stress of 15 MPa. The function that
    does the same for the band-free rotor must first reproduce these."""
    from bench import b4_ccx as b4
    h = b4.hand_estimate(with_band(), 10_000)
    assert h["cap_area_mm2"] == pytest.approx(162.0, abs=0.5) and h["cap_centroid_mm"] == pytest.approx(40.8, abs=0.05)
    assert h["magnet_area_mm2"] == pytest.approx(28.5, abs=0.05) and h["magnet_centroid_mm"] == pytest.approx(39.9, abs=0.05)
    assert h["mass_per_pole_g"] == pytest.approx(29.0, abs=0.5)
    assert h["pull_kn_per_m"] == pytest.approx(65.0, abs=0.5) and h["pull_n_over_stack"] == pytest.approx(1300.0, abs=10.0)
    assert h["strips_mm"]["centre_rib"] == pytest.approx(1.5, abs=1e-6) and h["strips_mm"]["outer_bridge"] == pytest.approx(1.4, abs=1e-6)
    assert h["strips_mm"]["total"] == pytest.approx(4.3, abs=1e-6)
    assert h["nominal_mpa"] == pytest.approx(15.0, abs=0.5)
    assert h["magnet_share_of_pull"] == pytest.approx(0.15, abs=0.01)         # "the magnets give 15 % of it"


def test_the_hand_estimate_for_the_band_free_rotor_has_other_strips_and_says_so():
    """Without the band there is no 1.5 mm rib between two bands. The cap is the steel in the V's mouth, held by the rib between the two magnet slots (3.30 mm) and the two
    bridges between a slot's outer end and the rim (1.56 mm each), so the same method gives 6.6 MPa, not 15. This is what the case holds the finite-element solve against."""
    from bench import b4_ccx as b4
    h = b4.hand_estimate(band_free(), 10_000)
    assert h["cap_area_mm2"] == pytest.approx(93.7, abs=0.2) and h["cap_centroid_mm"] == pytest.approx(42.0, abs=0.05)
    assert h["strips_mm"]["centre_rib"] == pytest.approx(3.298, abs=0.002) and h["strips_mm"]["outer_bridge"] == pytest.approx(1.5585, abs=0.002)
    assert h["strips_mm"]["total"] == pytest.approx(3.298 + 2 * 1.5585, abs=0.005)
    assert h["pull_kn_per_m"] == pytest.approx(42.5, abs=0.2) and h["nominal_mpa"] == pytest.approx(6.62, abs=0.05)
    assert h["magnet_share_of_pull"] == pytest.approx(0.22, abs=0.01)


def test_the_hand_estimate_grows_with_the_square_of_the_speed():
    from bench import b4_ccx as b4
    low, high = b4.hand_estimate(band_free(), 10_000), b4.hand_estimate(band_free(), 12_000)
    assert high["pull_kn_per_m"] / low["pull_kn_per_m"] == pytest.approx(1.44, rel=1e-9)
    assert high["nominal_mpa"] / low["nominal_mpa"] == pytest.approx(1.44, rel=1e-9)
    assert low["cap_area_mm2"] == high["cap_area_mm2"]


def test_the_two_poles_of_the_sector_are_mirror_images_so_they_have_the_same_cap():
    """Pole 1 (67.5 degrees) is pole 0 (22.5) mirrored about the 45-degree line: the same cap, the same strips."""
    from bench import b4_ccx as b4
    g0, g1 = b4.strip_geometry(band_free(), 0), b4.strip_geometry(band_free(), 1)
    assert g0["axis_deg"] == pytest.approx(22.5) and g1["axis_deg"] == pytest.approx(67.5)
    assert g0["cap"].area == pytest.approx(g1["cap"].area, rel=1e-4)
    for key in ("rib", "bridge_low", "bridge_high"):
        assert g0["widths"][key] == pytest.approx(g1["widths"][key], rel=1e-4)
        assert g0["cuts"][key].length == pytest.approx(g1["cuts"][key].length, rel=1e-3)


# ----------------------------------------------------------------------------------------------- 3. the mesh
@pytest.fixture(scope="module")
def coarse():
    from bench import b4_ccx as b4
    p = band_free()
    return p, mg.build(p, gap="sliding"), b4.mesh_rotor(p, mg.build(p, gap="sliding"), h_min=0.5, h_max=2.0)


def test_the_mesh_is_second_order_triangles_with_positive_areas_and_the_areas_of_motorgeo(coarse):
    """Self-check: the area of every 6-node triangle is integrated over its curved shape (6 Gauss points), so the sum over the steel is motorgeo's steel to a few
    parts in 10^5, and every triangle has a positive Jacobian (CalculiX refuses the others)."""
    from bench import b4_ccx as b4
    p, regions, mesh = coarse
    visible = mg.visible(regions)
    assert all(len(t) == 6 for t in mesh.steel) and mesh.steel
    steel = sum(b4.t6_area(mesh.nodes, t) for t in mesh.steel)
    assert steel == pytest.approx(visible["rotor_steel"].area, rel=5e-4)
    assert set(mesh.magnets) == {f"magnet_{k}_{s}" for k in (0, 1) for s in ("up", "dn")}
    for name, triangles in mesh.magnets.items():
        assert sum(b4.t6_area(mesh.nodes, t) for t in triangles) == pytest.approx(visible[name].area, rel=5e-4)
    assert min(b4.t6_area(mesh.nodes, t) for t in mesh.steel) > 0
    assert min(b4.t6_area(mesh.nodes, t) for ts in mesh.magnets.values() for t in ts) > 0


def test_a_magnet_shares_the_nodes_of_its_slot_walls_with_the_steel_and_the_slot_outlines_are_listed(coarse):
    from bench import b4_ccx as b4
    p, regions, mesh = coarse
    steel_nodes = {n for t in mesh.steel for n in t}
    for name, triangles in mesh.magnets.items():
        outline = mesh.boundary[f"slot:{name}"]
        assert outline and all(len(edge) == 3 for edge in outline)
        wall_nodes = {n for edge in outline for n in edge}
        assert wall_nodes <= steel_nodes                                  # every node of the slot wall belongs to a steel triangle
        assert wall_nodes <= {n for t in triangles for n in t}            # and to a magnet triangle: one node, two bodies
    assert mesh.quality["min_sicn"] > 0.1 and mesh.quality["triangles"] == len(mesh.steel) + sum(len(v) for v in mesh.magnets.values())


def test_the_mesh_is_finer_at_the_slots_and_a_smaller_h_min_gives_more_elements(coarse):
    from bench import b4_ccx as b4
    p, regions, mesh = coarse
    finer = b4.mesh_rotor(p, regions, h_min=0.3, h_max=2.0)
    assert len(finer.steel) > 1.5 * len(mesh.steel)
    edge_lengths = [math.dist(finer.nodes[e[0]], finer.nodes[e[1]]) for e in finer.boundary["slot:magnet_0_up"]]
    assert max(edge_lengths) < 0.45                                       # along the slot walls the elements are about h_min long


def test_mesh_rotor_leaves_the_path_of_this_process_whole(path_is_put_back):
    """Gmsh cuts the PATH of the process it runs in (Task 4, bench/__init__.py). The mesher is wrapped in `with intact_path():`; the process's own PATH, read from Windows,
    is the same afterwards. Self-check: take the `with` off in b4_ccx.mesh_rotor and this fails (and so does the scan in tests/test_bench.py)."""
    from bench import b4_ccx as b4
    p = band_free()
    before = process_path()
    b4.mesh_rotor(p, mg.build(p, gap="sliding"), h_min=0.8, h_max=2.0)
    assert process_path() == before == path_is_put_back


def test_the_area_of_a_curved_triangle_is_integrated_not_taken_from_its_corners():
    """A triangle (0,0), (1,0), (0,1) whose edge 1-2 is a parabola 0.1 below its chord: the area is 1/2 + (2/3) 0.1 = 0.5667 (Archimedes), not the corners' 1/2."""
    from bench import b4_ccx as b4
    nodes = {1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.0, 1.0), 4: (0.5, -0.1), 5: (0.5, 0.5), 6: (0.0, 0.5)}
    assert b4.t6_area(nodes, (1, 2, 3, 4, 5, 6)) == pytest.approx(0.5 + 2 / 3 * 0.1, rel=1e-12)
    nodes[4] = (0.5, 0.0)
    assert b4.t6_area(nodes, (1, 2, 3, 4, 5, 6)) == pytest.approx(0.5, rel=1e-12)
    cx, cy = b4.t6_first_moments(nodes, (1, 2, 3, 4, 5, 6))             # the first moments of a straight triangle: area x centroid
    assert cx == pytest.approx(0.5 / 3, rel=1e-12) and cy == pytest.approx(0.5 / 3, rel=1e-12)


# ----------------------------------------------------------------------------------------------- 4. the magnets as a load, and the deck
def toy_mesh():
    """Two 6-node triangles that share an edge: the steel (corners (0,0), (2,0), (0,2)) and a magnet (corners (2,0), (2,2), (0,2)). The shared edge 2-5-3 is the slot wall."""
    from bench import b4_ccx as b4
    nodes = {1: (0.0, 0.0), 2: (2.0, 0.0), 3: (0.0, 2.0), 4: (1.0, 0.0), 5: (1.0, 1.0), 6: (0.0, 1.0), 7: (2.0, 2.0), 8: (2.0, 1.0), 9: (1.0, 2.0)}
    return b4.Mesh(nodes=nodes, steel=[(1, 2, 3, 4, 5, 6)], magnets={"magnet_0_up": [(2, 7, 3, 8, 9, 5)]}, boundary={"slot:magnet_0_up": [(2, 3, 5)]})


def test_the_magnet_s_centrifugal_force_is_mass_times_omega_squared_times_the_radius_of_its_centre():
    """The force of one magnet at 10,000 rpm over the 20 mm stack: m w^2 r with m = 7,600 kg/m3 x (9.08 x 1.57 mm2) x 20 mm (the note's magnet, not the slot's volume), r the
    distance of the slot's centre from the axis, along that radius. Its parts along the slot's own axes (u along the slot, n across it) are worked out here from motorgeo's
    parameters, not from the slot's polygon, which is what slot_frames reads."""
    from bench import b4_ccx as b4
    p = band_free()
    frames = b4.slot_frames(p, mg.build(p, gap="sliding"), 10_000)
    assert set(frames) == {"magnet_0_up", "magnet_0_dn", "magnet_1_up", "magnet_1_dn"}
    w2 = (2 * math.pi * 10_000 / 60) ** 2
    mass = 7600.0 * (p.mag_w * 1e-3) * (p.mag_t * 1e-3) * (p.stack * 1e-3)
    for name, frame in frames.items():
        k, side = int(name.split("_")[1]), (1 if name.endswith("up") else -1)
        theta = 22.5 + 45 * k
        c = mg.rot((p.mag_x, side * p.mag_y), theta)
        u = mg.rot((math.cos(math.radians(p.v_half_deg)), side * math.sin(math.radians(p.v_half_deg))), theta)
        n = mg.rot((math.sin(math.radians(p.v_half_deg)), -side * math.cos(math.radians(p.v_half_deg))), theta)
        force = np.array(c) * 1e-3 * mass * w2
        assert frame.mass_kg == pytest.approx(mass, rel=1e-9)
        assert np.hypot(*frame.force_n) == pytest.approx(np.hypot(*force), rel=1e-4) and np.hypot(*force) == pytest.approx(94.73, abs=0.02)
        assert frame.centre == pytest.approx(c, abs=2e-3)
        assert frame.f_across == pytest.approx(abs(force @ np.array(n)), rel=1e-3) and frame.f_along == pytest.approx(abs(force @ np.array(u)), rel=1e-3)
        assert frame.f_across == pytest.approx(73.98, abs=0.05) and frame.f_along == pytest.approx(59.17, abs=0.05)
        assert frame.f_across > 0 and frame.f_along > 0                   # the magnet presses on the long wall on its n side and on the end wall on its u side
        assert np.dot(frame.n, frame.u) == pytest.approx(0.0, abs=1e-9) and np.hypot(*frame.n) == pytest.approx(1.0) and np.hypot(*frame.u) == pytest.approx(1.0)
        assert frame.half_length == pytest.approx(p.slot_w / 2, abs=1e-3) and frame.half_thickness == pytest.approx(p.slot_t / 2, abs=1e-3)


def test_the_walls_the_magnet_bears_on_get_nodal_forces_that_add_up_to_its_force_and_the_force_is_conserved(coarse):
    """The long wall on the n side takes f_across along n, the end wall on the u side takes f_along along u (a frictionless magnet in its slot). Over a whole mesh the nodal
    forces of the four slots add up to the four magnets' force vectors: what the bonded model's body force would be (mass-equivalent), to the last digit."""
    from bench import b4_ccx as b4
    p, regions, mesh = coarse
    frames = b4.slot_frames(p, regions, 10_000)
    loads = b4.slot_loads(p, mesh, frames)
    total = np.zeros(2)
    for name, frame in frames.items():
        walls = loads[name]
        assert walls["long"]["edges"] and walls["end"]["edges"]
        long_len = sum(math.dist(mesh.nodes[e[0]], mesh.nodes[e[1]]) for e in walls["long"]["edges"])
        end_len = sum(math.dist(mesh.nodes[e[0]], mesh.nodes[e[1]]) for e in walls["end"]["edges"])
        assert long_len == pytest.approx(p.slot_w - 2 * p.slot_fillet, abs=2e-3) and end_len == pytest.approx(p.slot_t - 2 * p.slot_fillet, abs=2e-3)
        assert walls["long"]["pressure_mpa"] == pytest.approx(frame.f_across / (long_len * p.stack), rel=1e-12)
        assert walls["end"]["pressure_mpa"] == pytest.approx(frame.f_along / (end_len * p.stack), rel=1e-12)
        total += np.array(frame.force_n)
    summed = np.zeros(2)
    for node, (fx, fy) in b4.nodal_forces(loads).items():
        summed += (fx, fy)
    assert summed == pytest.approx(total, rel=1e-9)


def test_nodal_forces_of_a_uniform_load_on_a_straight_three_node_edge_are_one_sixth_four_sixths_one_sixth():
    """The consistent loads of a uniform line load on a quadratic edge: the ends carry 1/6 of the force each, the mid node 4/6."""
    from bench import b4_ccx as b4
    nodes = {1: (0.0, 0.0), 2: (3.0, 0.0), 3: (1.5, 0.0)}
    got = b4.edge_nodal_forces(nodes, [(1, 2, 3)], direction=(0.0, 2.0), total_force=12.0)       # the direction is normalised: 12 N along +y
    assert got[1] == pytest.approx((0.0, 2.0)) and got[2] == pytest.approx((0.0, 2.0)) and got[3] == pytest.approx((0.0, 8.0))
    two = b4.edge_nodal_forces({1: (0.0, 0.0), 2: (1.0, 0.0), 3: (0.5, 0.0), 4: (3.0, 0.0), 5: (2.0, 0.0)}, [(1, 2, 3), (2, 4, 5)], (1.0, 0.0), 8.0)
    assert sum(f[0] for f in two.values()) == pytest.approx(8.0) and two[2][0] == pytest.approx(8.0 * (1 / 6 * 1 / 3 + 1 / 6 * 2 / 3))      # the shared node adds both edges


def test_the_pressure_deck_has_the_steel_only_the_nodal_forces_and_the_cards_ccx_needs():
    from bench import b4_ccx as b4
    mesh = toy_mesh()
    text = b4.deck_text(mesh, "pressure", 10_000, forces={5: (1.5, -2.5)}, node_sets={"LONG_magnet_0_up": [2, 3, 5]})
    lines = text.splitlines()
    assert "*ELEMENT, TYPE=CPS6, ELSET=STEEL" in lines and not any("ELSET=MAGNET" in l for l in lines)
    assert any(l.startswith("*SOLID SECTION, ELSET=STEEL, MATERIAL=M19") for l in lines) and "20" in lines[lines.index([l for l in lines if l.startswith("*SOLID SECTION")][0]) + 1]
    assert "*EQUATION" not in lines                                       # no ties in this model
    assert "*CLOAD" in lines and "5, 1, 1.5" in lines and "5, 2, -2.5" in lines
    assert "*NSET, NSET=LONG_magnet_0_up" in lines
    assert len([l for l in lines if re.fullmatch(r"\d+, [-\d.e+]+, [-\d.e+]+, 0\.0", l)]) == 6     # the steel's six nodes only: the magnet's own nodes 7, 8, 9 are not written
    assert "*NSET, NSET=EDGE0" in lines and "*NSET, NSET=EDGE90" in lines
    after = lines[lines.index("*BOUNDARY") + 1:]
    assert after[0] == "EDGE0, 2, 2" and after[1] == "EDGE90, 1, 1"
    assert lines.index("*STEP") < lines.index("*STATIC") < lines.index("*BOUNDARY") < lines.index("*DLOAD")       # the boundary conditions are the step's own: PrePoMax reads them there
    w2 = (2 * math.pi * 10_000 / 60) ** 2
    dload = [l for l in lines if l.startswith("STEEL, CENTRIF")]
    assert len(dload) == 1 and float(dload[0].split(",")[2]) == pytest.approx(w2, rel=1e-9)
    assert [t.strip() for t in dload[0].split(",")[3:]] == ["0", "0", "0", "0", "0", "1"]            # whole numbers, as PrePoMax writes them: it could not read "0." (found 2026-11-06)
    assert "*NODE PRINT, NSET=EDGE0, TOTALS=ONLY" in lines and "*EL FILE" in lines and lines[-1] == "*END STEP"
    numeric = [l for l in lines[lines.index("*HEADING") + 2:] if not l.startswith("*")]          # everything after the free-text heading line
    assert numeric and all(len(token.strip()) <= 20 for l in numeric for token in l.split(","))


def test_the_bonded_deck_ties_a_duplicate_of_every_wall_node_to_the_original_and_gives_the_magnet_the_mass_of_the_real_magnet():
    """The magnet is a second body that fills the slot: it gets its own copy of every node on the slot wall, tied to the steel's node (two equal-DOF equations per
    pair), so that CalculiX's nodal stress of the steel is the steel's alone (a shared node would average the two bodies). Its density is scaled so that the body has
    the mass of the real magnet (9.08 x 1.57 mm), not of the whole slot."""
    from bench import b4_ccx as b4
    mesh = toy_mesh()
    text = b4.deck_text(mesh, "bonded", 12_000, magnet_density_tonne_mm3=6.627128124470999e-9)
    lines = text.splitlines()
    assert "*ELEMENT, TYPE=CPS6, ELSET=MAGNET" in lines
    equations = lines[lines.index("*EQUATION") + 1:lines.index("*STEP")]
    assert len(equations) == 2 * 2 * 3 and equations[0] == "2"                   # 3 wall nodes (2, 3, 5) x 2 DOFs x (the count line + the terms)
    terms = [l for l in equations if l != "2"]
    assert len(terms) == 6 and all(re.fullmatch(r"\d+, [12], 1, \d+, [12], -1", t) for t in terms)
    originals = {int(t.split(",")[0]) for t in terms}
    duplicates = {int(t.split(",")[3]) for t in terms}
    assert originals == {2, 3, 5} and duplicates == {10, 11, 12} and not duplicates & set(mesh.nodes)       # new numbers above the largest
    magnet_row = [l for l in lines[lines.index("*ELEMENT, TYPE=CPS6, ELSET=MAGNET") + 1:] if not l.startswith("*")][0]
    assert magnet_row.split(", ")[1:] == ["10", "7", "11", "8", "9", "12"]       # the magnet's corner (2), (7), (3) and mid nodes, with the wall nodes replaced
    at = lines.index("*MATERIAL, NAME=N45")
    assert lines[at + 1:at + 5] == ["*ELASTIC", "160000, 0.24", "*DENSITY", "6.62712812447e-09"]       # E, nu, and the density at 12 digits (the 21-character one broke the first deck)
    assert "MAGNET, CENTRIF" in text and "*CLOAD" not in lines
    w2 = (2 * math.pi * 12_000 / 60) ** 2
    assert float([l for l in lines if l.startswith("MAGNET, CENTRIF")][0].split(",")[2]) == pytest.approx(w2, rel=1e-9)


def test_the_deck_total_load_is_what_the_equilibrium_check_compares_with_the_reaction():
    """external_load: the body force of every element, rho w^2 t (integral of x dA, integral of y dA) over the curved elements, plus the applied nodal forces. For the toy
    steel triangle (area 2, centroid (2/3, 2/3)) at 10,000 rpm it is rho w^2 t 2 (2/3)."""
    from bench import b4_ccx as b4
    mesh = toy_mesh()
    w2 = (2 * math.pi * 10_000 / 60) ** 2
    expected = 7.65e-9 * w2 * 20.0 * 2 * (2 / 3)
    fx, fy = b4.external_load(mesh, "pressure", 10_000, forces={5: (1.0, 2.0)})
    assert fx == pytest.approx(expected + 1.0, rel=1e-12) and fy == pytest.approx(expected + 2.0, rel=1e-12)
    bx, by = b4.external_load(mesh, "bonded", 10_000, magnet_density_tonne_mm3=5e-9)             # the magnet triangle (2,0), (2,2), (0,2): area 2, centroid (4/3, 4/3)
    assert bx == pytest.approx(7.65e-9 * w2 * 20.0 * 2 * (2 / 3) + 5e-9 * w2 * 20.0 * 2 * (4 / 3), rel=1e-12) and by == pytest.approx(bx, rel=1e-12)


# ----------------------------------------------------------------------------------------------- 5. reading the answer
DAT = """

                        S T E P       1


                                INCREMENT     1


 total force (fx,fy,fz) for set EDGE0 and time  0.1000000E+01

        4.006189E+01 -5.035066E+03 -5.426645E-13

 total force (fx,fy,fz) for set EDGE90 and time  0.1000000E+01

       -5.035063E+03  4.089694E+01 -2.175885E-13
"""


def test_the_reactions_are_read_from_the_dat_file_and_the_constrained_components_are_the_ones_that_count():
    """ccx prints the total force of a node set. Only the DOF that is constrained is a reaction: y on the x-axis edge, x on the other (the other component is what the
    body force puts on those nodes, not a reaction)."""
    from bench import b4_ccx as b4
    got = b4.parse_reactions(DAT)
    assert got["EDGE0"] == pytest.approx((40.06189, -5035.066, -5.426645e-13)) and got["EDGE90"] == pytest.approx((-5035.063, 40.89694, -2.175885e-13))
    assert b4.reaction_vector(got) == pytest.approx((-5035.063, -5035.066))      # (RF_x of EDGE90, RF_y of EDGE0): minus the external load
    assert b4.parse_reactions("nothing here") == {}


def test_von_mises_of_the_stress_states_that_have_a_known_answer():
    from bench import b4_ccx as b4
    assert b4.von_mises([100.0, 0, 0, 0, 0, 0]) == pytest.approx(100.0)                 # uniaxial: the stress itself
    assert b4.von_mises([0, 0, 0, 40.0, 0, 0]) == pytest.approx(math.sqrt(3) * 40.0)    # pure shear: sqrt(3) tau
    assert b4.von_mises([50.0, 50.0, 50.0, 0, 0, 0]) == pytest.approx(0.0, abs=1e-9)    # hydrostatic: none
    assert b4.von_mises([30.0, -30.0, 0, 0, 0, 0]) == pytest.approx(math.sqrt(30 ** 2 + 30 * 30 + 30 ** 2))      # plane stress, s1 = -s2: sqrt(3) s


def test_the_place_of_a_peak_is_named_by_the_strip_it_is_in():
    """The brief: the peak von Mises sits in a bridge or a rib, and is reported with its place. Points whose place is known: the outer corner of a slot (the farthest point
    of a slot from the axis, in the bridge to the rim), the narrowest point between the two slots of a pole (the centre rib), the shaft's bore, and a point in the bulk."""
    from bench import b4_ccx as b4
    p = band_free()
    far = b4.strip_geometry(p, 0)["cuts"]["bridge_low"]
    x, y = far.coords[0]                                                        # the inner end of the cut: the slot's farthest point
    place = b4.classify_place(p, x, y)
    assert place["feature"] == "outer bridge" and place["slot"] == "magnet_0_dn" and place["pole"] == 0 and place["in_a_strip"] is True
    assert place["r_mm"] == pytest.approx(43.44, abs=0.01) and place["angle_deg"] == pytest.approx(9.83, abs=0.05)
    cx, cy = b4.strip_geometry(p, 1)["cuts"]["rib"].interpolate(0.5, normalized=True).coords[0]
    assert b4.classify_place(p, cx, cy)["feature"] == "centre rib" and b4.classify_place(p, cx, cy)["pole"] == 1
    bore = b4.classify_place(p, 7.5 * math.cos(math.radians(60)), 7.5 * math.sin(math.radians(60)))
    assert bore["feature"] == "shaft bore" and bore["in_a_strip"] is False
    bulk = b4.classify_place(p, 25.0 * math.cos(math.radians(45)), 25.0 * math.sin(math.radians(45)))
    assert bulk["feature"] == "bulk" and bulk["in_a_strip"] is False


def test_the_force_through_a_cut_of_a_known_stress_field_is_stress_times_length_times_thickness():
    """A plain field: sigma_xx = 10 MPa everywhere on a mesh of the square (0,0)-(4,4). A vertical cut of length 3 carries 10 x 3 x t along x (the normal x); the same cut
    turned to a slope of 45 degrees carries t sigma (cut length) cos(45) of force along x with the normal n = (cos, sin)."""
    from bench import b4_ccx as b4
    xs, ys = np.meshgrid(np.linspace(0, 4, 9), np.linspace(0, 4, 9))
    nodes = {i + 1: (float(x), float(y)) for i, (x, y) in enumerate(zip(xs.ravel(), ys.ravel()))}
    triangles = []
    index = lambda i, j: j * 9 + i + 1
    for i in range(8):
        for j in range(8):
            triangles += [(index(i, j), index(i + 1, j), index(i + 1, j + 1)), (index(i, j), index(i + 1, j + 1), index(i, j + 1))]
    stress = {n: [10.0, 0.0, 0.0, 0.0, 0.0, 0.0] for n in nodes}
    from shapely.geometry import LineString
    force = b4.cut_force(nodes, triangles, stress, LineString([(2.0, 0.5), (2.0, 3.5)]), normal=(1.0, 0.0), thickness=20.0)
    assert force == pytest.approx((10.0 * 3.0 * 20.0, 0.0), rel=1e-6)
    inclined = b4.cut_force(nodes, triangles, stress, LineString([(1.0, 1.0), (3.0, 3.0)]), normal=(math.cos(math.radians(45)), -math.sin(math.radians(45))), thickness=20.0)
    assert inclined == pytest.approx((10.0 * 20.0 * 2 * math.sqrt(2) * math.cos(math.radians(45)), 0.0), rel=1e-6)
    shear = {n: [0.0, 0.0, 0.0, 5.0, 0.0, 0.0] for n in nodes}
    f = b4.cut_force(nodes, triangles, shear, LineString([(2.0, 0.5), (2.0, 3.5)]), normal=(1.0, 0.0), thickness=20.0)
    assert f == pytest.approx((0.0, 5.0 * 3.0 * 20.0), rel=1e-6)          # a shear stress on a vertical cut pulls along y


def test_the_peak_of_the_steel_ignores_the_nodes_of_the_magnet_bodies():
    """In the bonded deck the magnet's copies of the wall nodes carry the magnet's stress; the steel's peak is searched among the steel's nodes."""
    from bench import b4_ccx as b4
    mesh = toy_mesh()
    stress = {1: [1.0, 0, 0, 0, 0, 0], 2: [2.0, 0, 0, 0, 0, 0], 3: [3.0, 0, 0, 0, 0, 0], 4: [4.0, 0, 0, 0, 0, 0], 5: [5.0, 0, 0, 0, 0, 0], 6: [6.0, 0, 0, 0, 0, 0],
              7: [90.0, 0, 0, 0, 0, 0], 10: [80.0, 0, 0, 0, 0, 0]}
    value, node = b4.peak_of(stress, {n for t in mesh.steel for n in t})
    assert (value, node) == (6.0, 6)
    assert b4.peak_of({}, {1, 2}) == (None, None)


def test_the_case_is_registered_with_run_bench_after_b3():
    from bench import run_bench
    assert list(run_bench.CASES)[:4] == ["B1", "B2", "B3", "B4"] and run_bench.CASES["B4"] == "b4_ccx"


# ----------------------------------------------------------------------------------------------- 6. the run: ccx and cgx through run_program, the picture, the words about a place
def frd_text(stress: dict) -> str:
    """A .frd file with only the stress block, in CalculiX's fixed-width layout (the reader of Task 4 takes the node from columns 4-13 and six values of 12 characters)."""
    lines = [" -4  STRESS      6    1"]
    for node, values in stress.items():
        lines.append(" -1" + f"{node:10d}" + "".join(f"{v:12.5E}" for v in values))
    lines.append(" -3")
    return "\n".join(lines) + "\n"


def fake_model(coarse):
    from bench import b4_ccx as b4
    p, regions, mesh = coarse
    return {"p": p, "regions": regions, "mesh": mesh, "magnet_density": b4.mass_equivalent_density(p, regions)}


def test_the_total_load_of_the_real_sector_at_10000_rpm_is_what_ccx_confirmed(coarse):
    """ccx's reactions on the 0.1 mm mesh were 5035.066 N in each direction (and the same for the two magnet models); the load worked out from any mesh is the same to a few
    parts in 10^5, because the body force is integrated exactly over the curved elements and the magnets' force does not depend on the mesh."""
    from bench import b4_ccx as b4
    model = fake_model(coarse)
    at = b4.loads_at(model, 10_000)
    pressure = b4.external_load(model["mesh"], "pressure", 10_000, at["nodal"], model["p"].stack)
    bonded = b4.external_load(model["mesh"], "bonded", 10_000, None, model["p"].stack, magnet_density_tonne_mm3=model["magnet_density"])
    assert pressure == pytest.approx((5035.07, 5035.07), abs=0.6) and bonded == pytest.approx(pressure, rel=1e-5)
    assert math.hypot(*pressure) == pytest.approx(7120.64, abs=0.8)


def test_the_bonded_deck_of_the_real_mesh_has_two_equations_for_every_node_of_a_slot_wall_and_the_pressure_deck_none(coarse):
    from bench import b4_ccx as b4
    model = fake_model(coarse)
    mesh = model["mesh"]
    wall_nodes = {n for name, rows in mesh.boundary.items() if name.startswith("slot:") for row in rows for n in row}
    bonded = b4.deck_text(mesh, "bonded", 10_000, magnet_density_tonne_mm3=model["magnet_density"]).splitlines()
    block = bonded[bonded.index("*EQUATION") + 1:bonded.index("*STEP")]
    assert len(block) == 2 * 2 * len(wall_nodes) and block.count("2") == 2 * len(wall_nodes)
    steel_nodes = {n for t in mesh.steel for n in t}
    pressure = b4.deck_text(mesh, "pressure", 10_000, forces=b4.nodal_forces(b4.loads_at(model, 10_000)["loads"]))
    lines = pressure.splitlines()
    assert "*EQUATION" not in lines and len([l for l in lines if re.fullmatch(r"\d+, \S+, \S+, 0\.0", l)]) == len(steel_nodes)       # only the steel's nodes: the slots are holes
    assert b4.steel_node_count(mesh) == len(steel_nodes) < len(mesh.nodes)                  # the node count the pictures and the checks say for the pressure model: not the mesh's
    assert b4.mesh_counts(mesh)["nodes_steel"] == len(steel_nodes) and b4.mesh_counts(mesh)["nodes"] == len(mesh.nodes)


def test_solve_runs_ccx_by_its_full_path_through_run_program_with_a_timeout_and_reads_the_answer_back(monkeypatch, tmp_path, coarse):
    """The deck is written, `ccx -i job` is started by the path paths.tool gives (never by name), under run_program (the solver slot, below normal) with a timeout and at most
    8 threads, and job.frd and job.dat are read back. Here ccx is a stand-in that writes both."""
    from pathlib import Path

    import femm_queue
    import paths
    from bench import b4_ccx as b4
    model = fake_model(coarse)
    steel_nodes = sorted({n for t in model["mesh"].steel for n in t})
    seen = {}

    def fake_program(log, cmd, *, label, cwd, timeout_s, threads, name):
        seen.update(cmd=[str(c) for c in cmd], label=label, cwd=Path(cwd), timeout_s=timeout_s, threads=threads, name=name)
        (Path(cwd) / "job.frd").write_text(frd_text({n: [float(n % 7), 0.0, 0.0, 0.0, 0.0, 0.0] for n in steel_nodes}), encoding="ascii")
        (Path(cwd) / "job.dat").write_text(DAT, encoding="ascii")
        return femm_queue.SolverRun(0, "CalculiX Version 2.23, Copyright(C) 1998-2025 Guido Dhondt\nJob finished\n", "", 1.5)

    monkeypatch.setattr(b4, "run_program", fake_program)
    out = b4.solve(model, "pressure", 10_000, tmp_path / "run", None, "bench B4 test")
    assert seen["cmd"] == [str(paths.tool("ccx")), "-i", "job"] and seen["name"] == "ccx" and seen["cwd"] == tmp_path / "run"
    assert 0 < seen["timeout_s"] <= 3600 and 1 <= seen["threads"] <= paths.setting("femm", "workers") and seen["threads"] == b4.CCX_THREADS
    assert (tmp_path / "run" / "job.inp").is_file() and out["finished"] and out["version"] == "2.23" and out["returncode"] == 0
    assert set(out["stress"]) == set(steel_nodes) and out["reactions"]["EDGE0"][1] == pytest.approx(-5035.066)


@pytest.mark.parametrize("returncode, stdout, frd", [(1, "Job finished\n", True), (0, "*ERROR in the step: nothing was solved\n", True), (None, "", False), (0, "Job finished\n", False)])
def test_a_solve_that_did_not_finish_is_said_so_and_has_no_stress(monkeypatch, tmp_path, coarse, returncode, stdout, frd):
    """A non-zero exit, no 'Job finished' line, a timeout (the return code is None) and a missing .frd are all "not finished": nothing is read, and the check says so."""
    from pathlib import Path

    import femm_queue
    from bench import b4_ccx as b4
    model = fake_model(coarse)

    def fake_program(log, cmd, *, label, cwd, timeout_s, threads, name):
        if frd:                                                                      # the files are there: it is the exit code, the line or the timeout that says "not finished"
            (Path(cwd) / "job.frd").write_text(frd_text({1: [1.0] * 6}), encoding="ascii")
            (Path(cwd) / "job.dat").write_text(DAT, encoding="ascii")
        return femm_queue.SolverRun(returncode, stdout, "", 0.2)

    monkeypatch.setattr(b4, "run_program", fake_program)
    out = b4.solve(model, "pressure", 10_000, tmp_path / "run", None, "x")
    assert out["finished"] is False and out["stress"] is None and out["reactions"] == {}
    assert b4.analyse(model, out, 10_000) == {"finished": False}


def test_render_cgx_writes_the_script_that_makes_the_hardcopy_and_reads_cgx_s_maximum(monkeypatch, tmp_path):
    """cgx is started by its full path with -b and the script; the script reads the frd, picks the stress (dataset 2) and its von Mises (entity 7), plots faces with values, turns the
    view along z and makes the hardcopy as TGA (cgx's PNG needs ImageMagick, which is not installed). cgx does not come back from hcpy, so it is ended by the timeout and the
    frames are what counts."""
    from pathlib import Path

    import femm_queue
    import paths
    from bench import b4_ccx as b4
    seen = {}

    def fake_program(log, cmd, *, label, cwd, timeout_s, threads, name):
        seen.update(cmd=[str(c) for c in cmd], timeout_s=timeout_s, threads=threads)
        for i in range(seen.get("frames", 3)):
            (Path(cwd) / f"{i}__.tga").write_bytes(b"x" * 5000)
        return femm_queue.SolverRun(None, "2 L  101 1.000000   STRESS   Mises\n max:3.345405e+01 at node:527\n min:8.537002e-01 at node:34354\ncreate hcpy_1.png\n", "", 15.0)

    monkeypatch.setattr(b4, "run_program", fake_program)
    got = b4.render_cgx(tmp_path, None)
    assert seen["cmd"] == [str(paths.tool("cgx")), "-b", "shot.fbd"] and seen["timeout_s"] == b4.CGX_TIMEOUT_S and seen["threads"] == 1
    assert (tmp_path / "shot.fbd").read_text(encoding="ascii").split("\n")[:-1] == ["read job.frd", "ds 2 e 7", "plot fv all", "rot z 180", "hcpy tga", "quit"]
    assert got["frames_made"] and got["ended_by_timeout"] and got["max_von_mises"] == pytest.approx(33.45405)
    seen["frames"] = 1                                                                   # cgx wrote one frame only
    for name in tmp_path.glob("*.tga"):
        name.unlink()
    assert b4.render_cgx(tmp_path, None)["frames_made"] is False


def synthetic_frames(folder, model_colours):
    """cgx's frames in miniature: a 784 x 672 window with a colour bar of five swatches at the left, a 600 x 600 model frame whose coloured stripes are `model_colours`."""
    from PIL import Image
    window = np.full((672, 784, 3), 255, np.uint8)
    for i, colour in enumerate([(191, 0, 0), (255, 191, 0), (0, 128, 0), (0, 128, 153), (128, 0, 255)]):
        window[150 + 60 * i:200 + 60 * i, 8:26] = colour
    window[6:607, 172] = 0
    window[6:607, 774] = 0
    model = np.full((600, 600, 3), 255, np.uint8)
    for i, colour in enumerate(model_colours):
        model[100:400, 100 + 60 * i:160 + 60 * i] = colour
    model[300:303, :] = 0                                                       # a black outline
    Image.fromarray(window).save(folder / "0__.tga")
    Image.fromarray(model).save(folder / "1__.tga")
    Image.fromarray(np.full((60, 60, 3), 255, np.uint8)).save(folder / "2__.tga")


def test_the_cgx_picture_divides_out_cgx_s_lighting_only_after_checking_that_the_colours_are_the_legend_s(tmp_path):
    """cgx lights the model frame: every colour of it is the legend's colour times 0.298. The composition checks that (every colour of the frame is a legend colour times the factor,
    within 2.5 levels) before it divides the factor out; a frame whose colours are not the legend's is refused, not brightened."""
    from PIL import Image

    from bench import b4_ccx as b4
    legend = [(191, 0, 0), (255, 191, 0), (0, 128, 0), (0, 128, 153), (128, 0, 255)]
    lit = [tuple(round(0.298 * c) for c in colour) for colour in legend]
    synthetic_frames(tmp_path, lit)
    done = b4.cgx_picture(tmp_path, "a caption\nof two lines", tmp_path / "figure.png")
    assert done["made"] is True and done["colours_match"] == pytest.approx(1.0) and done["legend_colours"] == 5
    assert (tmp_path / "figure.png").stat().st_size > 5000 and Image.open(tmp_path / "figure.png").size[0] > 500
    other = tmp_path / "other"
    other.mkdir()
    synthetic_frames(other, [(30, 200, 40), (90, 90, 10), (10, 10, 120), (200, 20, 20), (5, 60, 60)])           # not the legend's colours at all
    refused = b4.cgx_picture(other, "x", other / "figure.png")
    assert refused["made"] is False and refused["colours_match"] < 0.9 and not (other / "figure.png").exists()


def test_the_place_on_a_slot_is_said_in_the_slot_s_own_axes(coarse):
    """Points of the slot's outline given in its axes (s along the slot, t across it; the loaded side is +t, the outer end +s) and what the case calls them."""
    from bench import b4_ccx as b4
    p, regions, mesh = coarse
    frame = b4.slot_frames(p, regions, 10_000)["magnet_0_dn"]
    at = lambda s, t: tuple(np.array(frame.centre) + s * np.array(frame.u) + t * np.array(frame.n))
    hw, ht = frame.half_length, frame.half_thickness
    assert b4.slot_detail(p, frame, *at(hw - 0.1, ht - 0.1)) == "corner fillet at the outer end, on the loaded side"
    assert b4.slot_detail(p, frame, *at(-hw + 0.1, -ht + 0.1)) == "corner fillet at the inner end, on the other side"
    assert b4.slot_detail(p, frame, *at(hw, 0.0)) == "outer end wall"
    assert b4.slot_detail(p, frame, *at(0.0, ht)) == "long wall, loaded side" and b4.slot_detail(p, frame, *at(0.0, -ht)) == "long wall, other side"
    assert b4.slot_detail(p, frame, *frame.centre) == "away from the slot's outline"


def test_the_distance_to_a_segment_is_the_distance_to_its_nearest_point():
    from bench import b4_ccx as b4
    points = np.array([[0.0, 1.0], [2.0, 0.0], [-1.0, 0.0], [5.0, 5.0], [1.0, 0.0]])
    got = b4.point_segment_distance(points, np.array([0.0, 0.0]), np.array([2.0, 0.0]))
    assert got == pytest.approx([1.0, 0.0, 1.0, math.hypot(3.0, 5.0), 0.0])


# ----------------------------------------------------------------------------------------------- 6. fix round 1: the soft filler, the decks' numbers, the pictures' words
def significant_digits(token: str) -> int:
    """The significant digits of a number as text: the mantissa's digits without the sign, the point, and the zeros at its ends ('200000.0' has 1, '6.62712812459e-09' has 12)."""
    mantissa = token.strip().lstrip("+-").lower().split("e")[0].replace(".", "")
    return len(mantissa.strip("0"))


NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")
NAME = re.compile(r"[A-Za-z_]\w*")


def tokens_of(deck: str) -> tuple[list[str], list[str]]:
    """(numbers, strange): the tokens of a deck's data lines (the lines after the free-text heading that are not keyword lines) that are numbers, and those that are neither a number nor
    a name (a set or a material: EDGE0, STEEL). The text of a numpy float, 'np.float64(2.6081)', is strange: a scan that only looked at what parses as a number would not see it."""
    numbers, strange = [], []
    lines = deck.splitlines()
    for line in lines[lines.index("*HEADING") + 2:]:
        if line.startswith("*"):
            continue
        for token in (t.strip() for t in line.split(",")):
            if NUMBER.fullmatch(token):
                numbers.append(token)
            elif not NAME.fullmatch(token):
                strange.append(token)
    return numbers, strange


def test_every_number_of_the_decks_of_the_real_mesh_is_at_most_twenty_characters_and_twelve_digits(coarse):
    """Fix round 1, item 3: the scan of the 20 characters ccx reads ran only on a toy mesh with short numbers, so a bypass of deck_number in the node, force or DLOAD lines went
    unseen. It runs on the real mesh's decks of the three variants, and also holds every number to the 12 digits of deck_number (a repr has 16 or 17): the DLOAD's number is
    18 characters when written by repr, under the 20, and only the digits show it."""
    from bench import b4_ccx as b4
    model = fake_model(coarse)
    at = b4.loads_at(model, 10_000)
    for variant in b4.ALL_VARIANTS:
        deck = b4.deck_text(model["mesh"], variant, 10_000, forces=at["nodal"], node_sets=b4.wall_sets(at["loads"]), thickness=model["p"].stack,
                            magnet_density_tonne_mm3=model["magnet_density"])
        numbers, strange = tokens_of(deck)
        assert len(numbers) > 10_000, variant                                         # the nodes' coordinates, the forces, the properties: the scan has something to look at
        too_long = [t for t in numbers if len(t) > 20]
        too_exact = [t for t in numbers if significant_digits(t) > 12]
        assert not too_long and not too_exact and not strange, (variant, too_long[:3], too_exact[:3], strange[:3])
    assert significant_digits("200000.0") == 1 and significant_digits("6.62712812459e-09") == 12 and significant_digits("-0.12345678901234567") == 17 and significant_digits("0.0") == 0


def test_the_soft_deck_is_the_bonded_deck_with_a_negligible_stiffness_and_the_same_mass_and_ties(coarse):
    """The diagnostic: the slot bodies keep their mass (the density of the bonded body) and their ties to the steel, and their Young's modulus is cut to SOFT_FILLER_YOUNGS_MPA, which
    is a negligible share (at most 0.1 %) of the steel's."""
    from bench import b4_ccx as b4
    assert b4.SOFT == "soft" and b4.ALL_VARIANTS == ("pressure", "bonded", "soft") and b4.TIED_VARIANTS == ("bonded", "soft") and b4.VARIANTS == ("pressure", "bonded")
    assert 0 < b4.SOFT_FILLER_YOUNGS_MPA <= 0.001 * b4.M19.youngs_mpa                      # at most 0.1 % of the steel's (the ruling)
    model = fake_model(coarse)
    args = dict(magnet_density_tonne_mm3=model["magnet_density"])
    bonded = b4.deck_text(model["mesh"], "bonded", 10_000, **args).splitlines()
    soft = b4.deck_text(model["mesh"], "soft", 10_000, **args).splitlines()
    assert soft.count("*EQUATION") == 1 and soft[soft.index("*EQUATION"):soft.index("*STEP")] == bonded[bonded.index("*EQUATION"):bonded.index("*STEP")]      # the same ties
    elastic = lambda lines, name: lines[lines.index(f"*MATERIAL, NAME={name}") + 2]
    assert elastic(bonded, "N45") == "160000, 0.24" and elastic(soft, "FILLER") == f"{b4.SOFT_FILLER_YOUNGS_MPA:g}, 0.24"
    density = lambda lines, name: lines[lines.index(f"*MATERIAL, NAME={name}") + 4]
    assert density(soft, "FILLER") == density(bonded, "N45") == b4.deck_number(model["magnet_density"])                      # the mass is kept
    assert "*SOLID SECTION, ELSET=MAGNET, MATERIAL=FILLER" in soft and "*SOLID SECTION, ELSET=MAGNET, MATERIAL=N45" in bonded
    assert f"MAGNET, CENTRIF, {b4.deck_number(b4.omega(10_000) ** 2)}, 0, 0, 0, 0, 0, 1" in soft and "*CLOAD" not in soft          # the body force acts on the filler, and no pressure is applied
    again = b4.deck_text(model["mesh"], "soft", 10_000, soft_youngs_mpa=2.0, **args).splitlines()
    assert elastic(again, "FILLER") == "2, 0.24"                                                                         # the modulus can be asked for (the scratch runs that showed the limit)
    with pytest.raises(ValueError, match="variant must be one of"):
        b4.deck_text(model["mesh"], "stiff", 10_000)
    assert b4.external_load(model["mesh"], "soft", 10_000, magnet_density_tonne_mm3=model["magnet_density"]) == b4.external_load(model["mesh"], "bonded", 10_000,
                                                                                                                                  magnet_density_tonne_mm3=model["magnet_density"])


def test_the_peak_of_the_slot_bodies_counts_every_node_of_the_bodies_and_none_of_the_steel():
    """Fix round 1, item 7: the bonded magnet bodies' peak was taken over the copies of the wall nodes only (the numbers above the mesh's); the bodies' interior nodes count too."""
    from bench import b4_ccx as b4
    mesh = toy_mesh()                                                                        # the steel uses nodes 1-6; the magnet's own nodes are 7, 8, 9 and the wall nodes 2, 3, 5
    stress = {1: (50.0,) * 3 + (0.0,) * 3, 5: (60.0,) * 3 + (0.0,) * 3,                         # steel nodes (5 is on the wall): large, and not the bodies'
              7: (1.0, 0.0, 0.0, 0.0, 0.0, 0.0), 8: (3.0, 0.0, 0.0, 0.0, 0.0, 0.0), 9: (2.0, 0.0, 0.0, 0.0, 0.0, 0.0),            # interior nodes of the body (above all: 3 MPa, node 8)
              101: (2.5, 0.0, 0.0, 0.0, 0.0, 0.0), 102: (0.5, 0.0, 0.0, 0.0, 0.0, 0.0)}           # copies of wall nodes
    assert b4.body_peak(mesh, stress) == pytest.approx(3.0)                                  # an interior node: the copies alone said 2.5
    stress[101] = (4.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert b4.body_peak(mesh, stress) == pytest.approx(4.0)


def test_the_soft_model_s_number_decides_between_the_two_hypotheses_when_it_is_near_one_of_the_models():
    """position = (soft - bonded) / (pressure - bonded): 1 is the pressure model's number (the lost stiffness), 0 the bonded model's (the load's path). It decides within SOFT_DECIDES of
    either; between, it does not. The run of 2026-11-06: 14.908 (pressure), 14.082 (bonded), 14.867 (soft): 95 % of the way."""
    from bench import b4_ccx as b4
    assert b4.SOFT_DECIDES == 0.25
    got = b4.soft_position(14.9085, 14.0821, 14.8674)
    assert got["position"] == pytest.approx(0.9504, abs=1e-3) and got["decided"] is True and got["hypothesis"] == "stiffness"
    low = b4.soft_position(14.9085, 14.0821, 14.15)
    assert low["decided"] is True and low["hypothesis"] == "load path" and low["position"] == pytest.approx(0.0823, abs=1e-3)
    middle = b4.soft_position(14.9085, 14.0821, 14.50)
    assert middle["decided"] is False and middle["hypothesis"] is None and 0.25 < middle["position"] < 0.75
    assert b4.soft_position(14.0, 14.0, 14.5) == {"position": None, "decided": False, "hypothesis": None}                    # nothing to decide
    edge = b4.soft_position(2.0, 1.0, 1.75)                                                   # exactly 3/4 of the way: within a quarter of the pressure model's
    assert edge["decided"] is True and edge["hypothesis"] == "stiffness"
    assert b4.soft_position(1.0, 2.0, 1.1)["hypothesis"] == "stiffness"                      # a pressure model below the bonded one: the same reading


def analysis_of(bore, bridges, load=7120.635, residual=0.00056, lame=14.0977):
    """The numbers of a B4 analysis that the soft-filler checks read."""
    return {"finished": True, "bore": {"hoop_mean_mpa": bore, "lame_mpa": lame}, "bridge_peak_mpa": bridges, "total_load_n": load,
            "equilibrium": {"residual_percent": residual, "reaction_n": (-5035.0 * load / 7120.635, -5035.0 * load / 7120.635)}}


def soft_run(returncode=0, seconds=20.5):
    return {"returncode": returncode, "timed_out": returncode is None, "finished": returncode == 0, "seconds": seconds, "version": "2.23"}


def test_the_soft_filler_checks_are_reported_and_say_which_hypothesis_the_number_supports():
    """Fix round 1: two checks, neither gating. They pass when ccx solved it cleanly (equilibrium, the pressure model's load) and when the number decides; the text states both
    hypotheses, the number, and what share of the bridges' excess is the stiffness."""
    from bench import b4_ccx as b4
    pressure, bonded = analysis_of(14.9085, 33.454), analysis_of(14.0821, 4.253)
    clean, cause = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(), "analysis": analysis_of(14.8674, 30.294)})
    assert clean["gating"] is False and cause["gating"] is False and clean["verdict"] == "PASS" and cause["verdict"] == "PASS"
    assert "ccx solves it cleanly" in clean["name"] and "+0.00000 %" in clean["measured"] and "0.00056 %" in clean["measured"]
    assert "(1)" in cause["expected"] and "(2)" in cause["expected"] and "25%" in cause["expected"] and "the bench's own limit" in cause["expected"]
    assert "95% of the way" in cause["measured"] and "hypothesis 1" in cause["measured"] and "lost hoop stiffness" in cause["measured"] and "89%" in cause["measured"]
    assert "the way the load reaches the steel barely does" in cause["measured"] and "reaches the steel does not" not in cause["measured"]          # T12-P3: 95 %, not 100
    assert "pressure 14.91 MPa (+5.75 % against Lame's 14.10)" in cause["measured"] and "bonded 14.08 (-0.11 %)" in cause["measured"] and "soft 14.87 (+5.46 %)" in cause["measured"]
    low = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(), "analysis": analysis_of(14.15, 5.0)})[1]
    assert low["verdict"] == "PASS" and "hypothesis 2" in low["measured"] and "load reaches the steel" in low["measured"]
    middle = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(), "analysis": analysis_of(14.50, 20.0)})[1]
    assert middle["verdict"] == "FAIL" and middle["gating"] is False and "does not decide" in middle["measured"] and "both contribute" in middle["measured"]      # a FAIL, with its number
    assert "14.50" in middle["measured"]


def test_the_soft_filler_checks_fail_with_the_reason_when_ccx_does_not_solve_it_cleanly_or_at_all():
    from bench import b4_ccx as b4
    pressure, bonded = analysis_of(14.9085, 33.454), analysis_of(14.0821, 4.253)
    clean = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(), "analysis": analysis_of(14.8674, 30.294, residual=0.2)})[0]
    assert clean["verdict"] == "FAIL" and "0.20000 %" in clean["measured"]                    # the reactions do not balance the load
    other = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(), "analysis": analysis_of(14.8674, 30.294, load=7200.0)})[0]
    assert other["verdict"] == "FAIL" and "+1.1" in other["measured"]                         # not the pressure model's load
    crashed = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(returncode=1), "analysis": {"finished": False}})
    assert [c["verdict"] for c in crashed] == ["FAIL", "FAIL"] and all(c["gating"] is False for c in crashed) and "ccx did not finish (exit 1" in crashed[0]["measured"]
    assert "no number" in crashed[1]["measured"]
    unsolved = b4.soft_filler_checks(pressure, bonded, None)
    assert [c["verdict"] for c in unsolved] == ["FAIL", "FAIL"] and "was not solved" in unsolved[0]["measured"]


def test_the_pictures_say_what_python_did_to_cgx_s_frames_and_what_the_white_lines_are_and_no_line_wraps():
    """Fix round 1, items 4 and C1: figure.png's caption names the 0.298 and the recomposed frames (as B5's names its crop); stress_full.png's title says that the white lines are the cuts
    of both poles; every line of both fits the width a title is wrapped to."""
    import bench
    from bench import b4_ccx as b4
    peak = {"value_mpa": 33.454, "place": {"feature": "outer bridge", "slot": "magnet_0_dn"}}
    caption, title = b4.figure_caption(peak, 45_271), b4.stress_full_title(peak, 45_271)
    assert "Python recomposed cgx's frames" in caption and "0.298" in caption and "45,271 nodes" in caption and "outer bridge of magnet_0_dn" in caption and "33.5 MPa" in caption
    assert "of each of the sector's two poles" in title and "pole 0" not in title and "45,271 nodes" in title and "the rib and the two bridges" in title
    for text in (caption, title):
        assert all(len(line) <= bench.TITLE_CHARS for line in text.split("\n")) and bench.wrap_title(text) == text


def test_the_soft_models_cause_check_decides_only_with_a_clean_solve_and_a_number_between_the_two_models():
    """T12-P5. The cause check is reported, never gating, but it must not print a hypothesis for a solve that went wrong. It decides only when the clean-solve check passed AND the soft
    model's bore number lies between the bonded model's (0) and the pressure model's (1), the ends included. Otherwise it fails and says why it cannot decide."""
    from bench import b4_ccx as b4
    pressure, bonded = analysis_of(14.9085, 33.454), analysis_of(14.0821, 4.253)

    def both(analysis):
        clean, cause = b4.soft_filler_checks(pressure, bonded, {"run": soft_run(), "analysis": analysis})
        return clean, cause

    for wild, where in ((40.0, "above"), (10.0, "below")):                                    # a diverged solve: far from both models, on either side
        clean, cause = both(analysis_of(wild, 30.294))
        assert clean["verdict"] == "PASS" and cause["verdict"] == "FAIL" and cause["gating"] is False, where
        assert "outside the range" in cause["measured"] and "does not decide" in cause["measured"], where
        assert "(hypothesis 1)" not in cause["measured"] and "(hypothesis 2)" not in cause["measured"], where
    above, below = b4.soft_position(14.9085, 14.0821, 40.0), b4.soft_position(14.9085, 14.0821, 10.0)
    assert above["position"] > 1 and below["position"] < 0
    assert above["decided"] is False and above["hypothesis"] is None and below["decided"] is False and below["hypothesis"] is None
    clean, cause = both(analysis_of(14.8674, 30.294, residual=0.2))                          # the number looks like 95 %, but the reactions do not balance the load
    assert clean["verdict"] == "FAIL" and cause["verdict"] == "FAIL" and "clean-solve check" in cause["measured"] and "does not decide" in cause["measured"]
    for edge, hypothesis in ((14.9085, "hypothesis 1"), (14.0821, "hypothesis 2")):          # exactly the pressure model's number, exactly the bonded model's: inside the range
        clean, cause = both(analysis_of(edge, 30.294))
        assert clean["verdict"] == "PASS" and cause["verdict"] == "PASS" and hypothesis in cause["measured"], edge
