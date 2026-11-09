"""B2: Python -> Gmsh (the API) -> .msh, on the motor's own geometry: motorgeo's sector, and the CAD drawings.

What it shows, and why
    Gmsh is the mesh maker for everything that is not FEMM: Elmer and CalculiX read its meshes. B2 proves the two ways the pipeline will make one,
    and that what comes out is the motor and not an approximation of it that nobody checked:
      Path 1  motorgeo.to_gmsh meshes the no-band 90-degree sector (R11: motorgeo.Params(barrier="none"), the sliding air gap that Elmer wants).
              This is the mesh B3 hands to Elmer.
      Path 2  The CAD's own drawings: the rotor as the CAD last exported it, reference\\SynRM rotor.DXF (2026-01-22), which HAS the air band, and the
              stator, reference\\Stator_refined.DXF (byte-for-byte copies of the CAD exports, read-only), through ezdxf (reads the drawing) -> shapely
              (finds the faces, and removes the air band in memory) -> Gmsh (meshes them). This is the road from the drawing to a mesh with no motorgeo
              in between, and the check that motorgeo's model IS the drawing. The brief named reference\\Rotor.DXF, the CAD's earlier export
              (2026-01-13): it never had the band, so there was nothing in it to remove. A design decision (2026-11-05) is that the band is removed in
              Python and the files are left alone, so path 2 reads the drawing that has the band, and Rotor.DXF stays as the independent check that the
              removal took the band and nothing else.
    A mesh is sound when every region of the motor is in it as a physical group, the elements are not degenerate, and the areas are the
    areas of the drawing. Each of the three is a check below, with its source.

Method
    Path 1   motorgeo.build(Params(barrier="none"), gap="sliding") -> motorgeo.to_gmsh(..., size_max=1.5) -> out\\bench\\B2\\path1.msh (MSH 2.2,
             what ElmerGrid reads). to_gmsh is wrapped in checks.intact_path(): gmsh.initialize() cuts this process's PATH short (Task 4's
             finding), and a program started by name afterwards would not be found. The file is read back in a second Gmsh session
             (inspect_mesh): every physical group with its element count and its area (or length), the area of each material against
             motorgeo's (the shapely polygons of the visible regions), and the quality of every triangle (Gmsh's minSICN, the signed inverse condition
             number: 1 for an equilateral triangle, 0 for a flat one, negative for an inverted one) with the place of the worst. A table of mesh sizes
             (size_max 1.5 ... 0.25 mm) shows how the areas converge. Two node counts are told apart: the saved file holds 8,819 nodes (what ElmerGrid
             reads and inspect_mesh counts), and motorgeo.to_gmsh RETURNS 8,906, because its Gmsh session also has the centre point of each of the
             87 arc edges, geometry that no element uses and the file does not hold (make_path1_mesh names both). The shortfall of the N45 area is
             explained, not just shown: the mesh edges lying on each fillet are COUNTED from the mesh's own nodes (fillet_chords), the fillets' radius
             is the parameter p.slot_fillet, and the shortfall has to be the chord error of those fillets within 1 %.
    Path 2   ezdxf reads the three files (SynRM rotor.DXF, Rotor.DXF, Stator_refined.DXF: LINE, ARC and CIRCLE entities; the "SOLIDWORKS Educational
             Product" text is not geometry); every arc is flattened to a polyline within 0.25 micrometre (the geometry: areas and positions are measured on
             this one). The mesh is made from a second, coarser reading (chords within 5 micrometres, no segment under 0.4 mm), with its element sizes from a
             size field (0.3 mm at the air gap, 0.5 mm at the magnet pockets, growing to 1.5 mm): the corners of a finely flattened circle would otherwise fix
             the mesh size of the whole rotor. The stator's bore is drawn as 36 tooth-top arcs; the slot openings are closed by arcs of the same radius
             (bore_mouths), so that the air gap and each slot become faces of their own. shapely polygonizes the rotor's edges alone first, to find its pockets
             (32 in SynRM rotor.DXF: 16 magnet slots and 16 air bands), and then, with the band edges dropped, all the edges together (unary_union, polygonize)
             into faces: the rotor steel, the 16 magnet slots of SynRM rotor.DXF, the shaft, the gap, the stator steel and the 36 slots. THE AIR BAND IS REMOVED IN
             MEMORY: SynRM rotor.DXF has 16 air-band pockets, one beside each magnet slot; a pocket whose smallest bounding rectangle is not the magnet slot's
             (slot_w x slot_t of motorgeo.Params) is a band, and its edges are dropped before the faces are made, so the steel closes over it (a wall that a
             band and a slot share is the slot's, and stays). The files are never written: their SHA-256 is taken before path 2 reads them and after everything B2
             did, and must be the same. The faces go to Gmsh's built-in kernel with their shared edges shared (no boolean operation, so no tolerance to
             choose), one physical group per face, the two A = 0 circles as curve groups, and the mesh is written as path2.msh.
             The removal is checked, not assumed. The number of bands to expect is counted from motorgeo (4 band outlines of Params(barrier="air") in the
             sector x 4 sectors = 16) and each band removed is held against motorgeo's: its rectangle, area, centroid radius, place after the slots' turn and
             distance to its slot. The edges dropped are shown to be exactly the bands' (96 = 16 x motorgeo's outline of 5 LINE + 1 ARC), and every other
             LINE, ARC and CIRCLE of the drawings must be an edge of a face (the text in them is listed as ignored). And after the removal SynRM rotor.DXF IS
             Rotor.DXF: the same slots, circles and steel, to a limit derived from the arc flattening, and the steel the removal added is the area of the
             bands. bands.png draws the rotor's faces, the removed bands and the magnet slots each in a colour of their own.
             The drawing is the whole cross-section (360 degrees), motorgeo's model a 90-degree sector: the areas are compared with
             4 x the sector's. The pole pattern and the slot pattern of the drawing are also compared with motorgeo's, as a rotation (the CAD
             draws its rotor and stator a little turned against motorgeo's frame; each is a check, the stator's in millimetres of arc at the bore).
             The mesh Gmsh WROTE is held against motorgeo's areas too (rotor steel, magnets, stator steel, slots, within the brief's 0.5 %), and the
             magnets' shortfall is explained by their fillets' chords as in path 1.
    The pictures   bands.png (above), and Gmsh's own window renders both meshes (a .geo with Draw; Print "...png"; Exit;, run through femm_queue.run_solver
             like every program), and figure.png puts the two renders side by side.

Run (from the Python to FEMM folder; the Gmsh library takes no solver slot, the Gmsh app that draws the picture takes it for a few seconds)
    uv run python bench/run_bench.py --only B2
Writes out\\bench\\B2\\: results.json, path1.msh, path2.msh, path1.png, path2.png, figure.png, bands.png, run.log, runs.json.

Numbers to expect (2026-11-06, Gmsh 4.15.2; the whole case is 11 s: path 1 with its table of mesh sizes 6.7 s, path 2 2.0 s, the pictures the rest)
    path 1   8,819 nodes, 16,914 triangles, all 40 physical groups; minimum SICN 0.251 in gap_rotor at (28.34, 35.30) mm (radius 45.27, angle 51.2 degrees),
             mean 0.965, 40 triangles below 0.3. Mesh area against motorgeo's: M-19 +0.037 %, N45 -0.172 %, air +0.319 %, coil -0.237 % (three of them
             outside the brief's 0.1 %: the chord error of curved boundaries; at size_max 0.25 mm: +0.010, -0.172, +0.001, -0.026 %). The rotor's surface is
             meshed with 1.47 mm edges (the steel's size) against 0.27 mm at mid-gap, which is why the worst triangle is in the gap.
             The 16 fillets of the sector are each cut by 2 mesh edges (counted on the mesh); the shortfall, 0.1127 mm2, is the chord error to 0.003 %.
    path 2   SynRM rotor.DXF: 144 lines + 80 arcs + 2 circles, stator 216 + 72 + 1 (the 4 MTEXT are ignored); 16 air bands found and removed, as motorgeo counts
             them (4 x 4); 96 edges dropped = 80 LINE + 16 ARC, none shared with a magnet slot, and every other entity is an edge of a face; each band is
             motorgeo's to 0.0001 mm (16.705 x 2.531 mm, 23.9915 mm2, centroid at radius 37.899 mm, 1.460 mm from its slot, after the slots' turn of 22.5 degrees);
             the steel goes from 5,539.54 to 5,923.40 mm2 (+383.86 mm2, the 16 bands) and Rotor.DXF (0 bands found) has the same steel to 0.0000 mm2 (the limit
             is 0.060 mm2) and the same 16 slots, shaft and outer circle. 16 magnet slots of 16.3478 mm2; 36 slots; areas of the drawing's faces against
             motorgeo's 4 x sector: rotor steel -0.0003 %, magnets -0.0018 %, stator steel -0.0079 %, slots +0.0174 %; the same against the MESH Gmsh wrote:
             rotor steel -0.005 %, magnets -0.172 %, stator steel +0.007 %, slots -0.048 % (the magnets' 0.4509 mm2 is the chord error of 64 fillets with 2
             edges each, 0.4510 mm2); the drawing is motorgeo's model turned by 22.5 degrees (rotor) and 5.0 degrees (stator), positions equal to 0.0000 mm.
             34,410 nodes, 68,229 triangles, 57 groups, minimum SICN 0.514

Self-check (the answers follow).
  1. The magnet slots have 0.3 mm fillets and the mesh covers 0.17 % less magnet area than the drawing, at every mesh size. Why does a finer
     size_max not help?
  2. The brief expects the CAD drawing's band to be removed in memory. How does the code know a pocket is a band and not a magnet?
  3. Why are the two faces of a shared wall (the rotor's surface, seen from the steel and from the gap) made from the same Gmsh line, and not
     meshed separately and glued?
  4. The brief names Rotor.DXF for path 2, and Rotor.DXF has no air band. Why does path 2 read SynRM rotor.DXF, and what is Rotor.DXF for?
Answers.
  1. The mesh size of a region is set by the size at its corner points (0.5 mm for a magnet) and size_max only caps it from above. A fillet
     of radius 0.3 mm is a quarter circle 0.47 mm long, which Gmsh meshes with two straight edges: each chord cuts a sliver of
     r^2/2 (theta - n sin(theta/n)) = 0.0070 mm2 off the circle (n = 2, theta = 90 degrees). A magnet slot has four fillets: 0.028 mm2 of
     its 16.35 mm2, which is the 0.17 %. Only a smaller size at the magnet's points would change it, and motorgeo.py is not changed by the bench.
  2. By its shape, not its position: the smallest rectangle that holds the pocket is 9.28 x 1.77 mm for the magnet slot (the parameters of
     motorgeo.Params) and 16.7 x 2.53 mm for the band of SynRM rotor.DXF (23.99 mm2 against 16.35, and the same as motorgeo's band), and the test allows
     0.05 mm. A rule that depended on where the magnets are would break as soon as the rotor is turned.
  3. Two separately meshed edges would have different nodes, and a field solver would then see a crack (a jump of the unknown) between the two
     sides unless a mortar condition glued them. With one line in both faces the nodes are shared and the mesh is conforming by construction.
  4. A rule for removing a band that has never met a band proves nothing. Rotor.DXF (the CAD's export of 2026-01-13) has no band: path 2 found 0 in it,
     and the rule had never run on a real one. SynRM rotor.DXF, the CAD's later export (2026-01-22), has sixteen, and the design decision is to remove the band in
     Python and leave the files as they are, so path 2 reads the drawing that has it. Rotor.DXF is what is left to compare with: the two exports are the same
     rotor, one with the band and one without, so after the removal they must be the same faces. That is an answer that does not come from the rule it tests.
"""
from __future__ import annotations

import dataclasses
import hashlib
import math
import time
from collections import Counter
from pathlib import Path

import numpy as np
from shapely import affinity
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import polygonize, unary_union
from shapely.prepared import prep

import motorgeo as mg
import paths
from checks import intact_path

from . import OWN_LIMIT, Log, check, compare, params, run_program, wrap_title

CASE = "B2"
TITLE = "Python -> Gmsh (the API) -> .msh: motorgeo's sector and the CAD drawings"

BOUNDARY_NAMES = ("A0", "edge_0", "edge_sector", "sliding_rotor", "sliding_stator")   # the curve groups motorgeo.to_gmsh makes (Elmer's boundaries)
SIZE_MAX = 1.5                                       # motorgeo.to_gmsh's default element size cap [mm]: the mesh B3 hands to Elmer
MESH_SIZES = (1.5, 1.0, 0.7, 0.5, 0.35, 0.25)        # the convergence table [mm]
AREA_TOLERANCE_PERCENT = 0.1                         # the brief: each material's area within 0.1 % of motorgeo's (path 1)
CAD_AREA_TOLERANCE_PERCENT = 0.5                     # the brief: the DXF path's steel and magnet areas within 0.5 % of motorgeo's no-band areas
MIN_SICN_OK = 0.1                                    # the bench's floor for "no degenerate element" (the brief only asks that the quality is stated)
CAD_POSITION_TOLERANCE_MM = 0.05                     # the drawing's pockets and slots against motorgeo's, after the best rotation (the CAD gives 0.01 mm)
SAGITTA_MM = 2.5e-4                                  # how far a flattened arc may leave the true arc: 0.25 micrometre (0.38 degrees at the rotor's surface)
MESH_SAGITTA_MM = 5e-3                               # the flattening the CAD MESH is made from: 5 micrometres (areas are still true to 2e-4)
MIN_SEGMENT_MM = 0.4                                 # ... and no polyline segment is shorter than this (a fillet of 0.47 mm gets 2 chords, a slot bottom 21)
POCKET_TOLERANCE_MM = 0.05                           # a pocket is a magnet slot if its smallest bounding rectangle is slot_w x slot_t within this
CHORD_WINDOW_PERCENT = 1.0                           # the measured shortfall of the magnets against the chord error of their fillets: within this share of it (the real agreement is 0.003 %)
GEOMETRY_KINDS = ("LINE", "ARC", "CIRCLE")           # the DXF entities that are geometry; anything else in a drawing (its text) is only counted
EDGE_ZONE_MM = 1e-4                                  # an entity is on the edge of a face if its whole polyline is within this of the faces' boundaries
BAND_COLOURS = {"steel": "#c4c9d1", "shaft": "#8d939c", "band": "#e0782f", "slot": "#1f6feb"}      # bands.png: the steel as drawn, the shaft, the removed bands, the magnet slots
SHOT_GEO = """// Gmsh's own window renders the mesh and saves it as a PNG, then closes.
Merge "{msh}";
General.GraphicsWidth = 1100;
General.GraphicsHeight = 1100;
General.Axes = 0;
General.SmallAxes = 0;
General.BackgroundGradient = 0;
General.Color.Background = {{255, 255, 255}};
General.Color.Foreground = {{40, 40, 40}};
General.Color.Text = {{20, 20, 20}};
Mesh.SurfaceFaces = 1;
Mesh.SurfaceEdges = 1;
Mesh.Lines = 1;
Mesh.LineWidth = 2;
Mesh.ColorCarousel = 2;
Print.Width = 1100;
Print.Height = 1100;
Draw;
Print "{png}";
Exit;
"""


# ============================================================================================== the groups and the areas
def expected_groups(regions) -> list[str]:
    """The physical groups a mesh of these regions must have: one for each region that is meshed, and the five boundary curve groups."""
    return [r.name for r in regions if r.material != "<No Mesh>"] + list(BOUNDARY_NAMES)


def geometry_areas(regions) -> dict[str, float]:
    """{material: area [mm2]} of what each material shows once later regions are painted over earlier ones: motorgeo's own areas (shapely)."""
    visible = mg.visible(regions)
    out: dict[str, float] = {}
    for r in regions:
        if r.material != "<No Mesh>":
            out[r.material] = out.get(r.material, 0.0) + visible[r.name].area
    return out


def material_areas(group_areas: dict[str, float], regions) -> dict[str, float]:
    """{material: area} from the area of each region's physical group: the regions of a material added up."""
    out: dict[str, float] = {}
    for r in regions:
        if r.material != "<No Mesh>":
            out[r.material] = out.get(r.material, 0.0) + group_areas[r.name]
    return out


def chord_deficit(radius: float, sweep_deg: float, segments: int) -> float:
    """The area between a circular arc of `sweep_deg` degrees and the `segments` straight chords that a mesh puts in its place:
    r^2 / 2 (theta - n sin(theta / n)). A mesh edge cuts this much off a convex arc, and adds it to the other side."""
    theta = math.radians(sweep_deg)
    return radius ** 2 / 2 * (theta - segments * math.sin(theta / segments))


def mesh_nodes(path) -> np.ndarray:
    """The (x, y) of every node of an MSH 2.2 file, in mm, read from the text of its $Nodes section (the file is plain text: no Gmsh needed)."""
    text = Path(path).read_text(encoding="ascii", errors="replace")
    block = text.split("$Nodes")[1].split("$EndNodes")[0].strip().splitlines()
    return np.array([line.split()[1:3] for line in block[1:]], float)


def fillet_chords(nodes, arcs, tolerance: float = 1e-6) -> list[int]:
    """How many mesh edges lie on each arc, counted from the mesh's own nodes. An arc is (cx, cy, r, start_deg, sweep_deg), counter-clockwise; a node is on it if its
    distance from the centre is r within `tolerance` and its angle is inside the sweep; the edges are the nodes on the arc less one (0 if none is on it)."""
    xy = np.asarray(nodes, float)
    out = []
    for cx, cy, r, start, sweep in arcs:
        dx, dy = xy[:, 0] - cx, xy[:, 1] - cy
        on_circle = np.abs(np.hypot(dx, dy) - r) < tolerance
        turned = (np.degrees(np.arctan2(dy, dx)) - start) % 360.0                    # the angle from the start of the arc, counter-clockwise
        slack = math.degrees(2 * tolerance / r)                                      # the same tolerance as an angle, so that a node at either end counts
        inside = (turned <= sweep + slack) | (turned >= 360.0 - slack)
        out.append(max(int((on_circle & inside).sum()) - 1, 0))
    return out


def model_fillets(regions, p: mg.Params) -> list[tuple]:
    """The corner fillets of motorgeo's magnet slots as arcs (cx, cy, r, start_deg, sweep_deg), counter-clockwise: the arc edges of the magnet regions whose
    radius is p.slot_fillet (the radius comes from the parameters, not from a number typed here)."""
    out = []
    for region in regions:
        if region.name.startswith("magnet_"):
            for edge in region.edges:
                if isinstance(edge, mg.Arc) and abs(edge.radius - p.slot_fillet) < 1e-6:
                    arc = edge.canonical()
                    out.append((arc.c[0], arc.c[1], arc.radius, arc.angles()[0], arc.sweep()))
    return out


def drawing_fillets(dxf: dict, p: mg.Params) -> list[tuple]:
    """The same from a drawing: its ARC entities of radius p.slot_fillet that are not centred on the axis (the bore's tooth tops are)."""
    return [(cx, cy, r, start, (end - start) % 360.0) for cx, cy, r, start, end in dxf["arcs"] if abs(r - p.slot_fillet) < 1e-6 and math.hypot(cx, cy) > 1e-9]


def chord_check(name: str, fillets: list[tuple], counts: list[int], deficit_mm2: float, area_mm2: float, source: str) -> dict:
    """Is the magnets' shortfall the chord error of their fillets? `fillets` are the arcs, `counts` the mesh edges counted on each (fillet_chords), `deficit_mm2`
    what the mesh lacks of the magnets' area, `area_mm2` that area. The prediction is the sum over the fillets of r^2/2 (theta - n sin(theta/n)) with n the COUNT on
    the mesh, so it needs one n for all of them; the shortfall has to be within CHORD_WINDOW_PERCENT of it."""
    uniform = len(counts) > 0 and len(set(counts)) == 1
    n = counts[0] if uniform else None
    radii, sweeps = sorted({round(f[2], 9) for f in fillets}), sorted({round(f[4], 6) for f in fillets})
    shape = f"radius {radii[0]:g} mm, {sweeps[0]:g} degrees" if len(radii) == 1 and len(sweeps) == 1 else "radii and sweeps as listed"
    predicted = sum(chord_deficit(f[2], f[4], n) for f in fillets) if uniform else None
    ok = uniform and abs(deficit_mm2 - predicted) <= CHORD_WINDOW_PERCENT / 100 * predicted
    expected = (f"the {len(fillets)} fillets ({shape}) are cut by {n if uniform else 'n'} mesh edges each, counted on the mesh; the shortfall is their chord error, the sum of "
                f"r^2/2 (theta - n sin(theta/n)) = {'n/a' if predicted is None else f'{predicted:.4f}'} mm2, within {CHORD_WINDOW_PERCENT:g} % of it ({OWN_LIMIT})")
    seen = ", ".join(f"{k} fillets with {v} edges" for v, k in sorted(Counter(counts).items()))
    measured = (f"{deficit_mm2:.4f} mm2 ({100 * deficit_mm2 / area_mm2:.3f} % of {area_mm2:.3f} mm2); mesh edges per fillet, counted on the mesh: {seen}"
                + ("" if uniform else "; the counts are not the same on every fillet, so no single n predicts the shortfall"))
    return check(name, ok, expected, measured, source)


# ============================================================================================== reading a mesh back
def inspect_mesh(path) -> dict:
    """Read an MSH file in a Gmsh session of its own and describe it: nodes, triangles, every physical group (dimension, name, elements, and the area of
    a surface group or the length of a curve group), the total area, and the quality of every triangle (Gmsh's minSICN) with the place of the worst.

    Areas and lengths are computed from the node coordinates with numpy, not asked of Gmsh's CAD kernel (a mesh file has none). The whole read is
    inside checks.intact_path(): gmsh.initialize() cuts this process's PATH."""
    import gmsh
    with intact_path():
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.open(str(path))
            tags, coords, _ = gmsh.model.mesh.getNodes()
            xyz = np.asarray(coords, float).reshape(-1, 3)
            index = np.full(int(max(tags)) + 1, -1, dtype=np.int64)
            index[np.asarray(tags, dtype=np.int64)] = np.arange(len(tags))
            groups: list[dict] = []
            all_q: list[np.ndarray] = []
            worst = {"min": 2.0}
            triangles = 0
            total_area = 0.0
            for dim, tag in sorted(gmsh.model.getPhysicalGroups()):
                name = gmsh.model.getPhysicalName(dim, tag)
                elements = 0
                measure = 0.0
                for entity in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                    types, element_tags, node_tags = gmsh.model.mesh.getElements(dim, entity)
                    for _, etags, ntags in zip(types, element_tags, node_tags):
                        per = dim + 1
                        corners = index[np.asarray(ntags, dtype=np.int64)].reshape(-1, per)
                        elements += len(etags)
                        if dim == 1:
                            measure += float(np.hypot(*(xyz[corners[:, 1], :2] - xyz[corners[:, 0], :2]).T).sum())
                            continue
                        a, b, c = xyz[corners[:, 0], :2], xyz[corners[:, 1], :2], xyz[corners[:, 2], :2]
                        measure += float((0.5 * np.abs((b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (c[:, 0] - a[:, 0]) * (b[:, 1] - a[:, 1]))).sum())
                        q = np.asarray(gmsh.model.mesh.getElementQualities([int(t) for t in etags], "minSICN"), float)
                        all_q.append(q)
                        triangles += len(etags)
                        k = int(np.argmin(q))
                        if q[k] < worst["min"]:
                            centre = (a[k] + b[k] + c[k]) / 3
                            worst = {"min": float(q[k]), "group": name, "element": int(etags[k]), "x_mm": float(centre[0]), "y_mm": float(centre[1]),
                                     "r_mm": float(math.hypot(*centre)), "angle_deg": float(math.degrees(math.atan2(centre[1], centre[0])) % 360)}
                entry = {"dim": dim, "tag": tag, "name": name, "elements": elements}
                if dim == 2:
                    entry["area_mm2"] = measure
                    total_area += measure
                else:
                    entry["length_mm"] = measure
                groups.append(entry)
            q_all = np.concatenate(all_q) if all_q else np.array([1.0])
            quality = {"min_sicn": float(q_all.min()), "mean_sicn": float(q_all.mean()), "worst": {k: v for k, v in worst.items() if k != "min"},
                       "below_0_2": int((q_all < 0.2).sum()), "below_0_3": int((q_all < 0.3).sum()), "inverted": int((q_all <= 0).sum())}
            return {"path": str(path), "nodes": int(len(tags)), "triangles": triangles, "groups": groups, "total_area_mm2": total_area, "quality": quality}
        finally:
            gmsh.finalize()


# ============================================================================================== path 1: motorgeo.to_gmsh
def make_path1_mesh(path, size_max: float = SIZE_MAX) -> dict:
    """motorgeo.to_gmsh on the no-band sector (the sliding air gap), written as `path` (MSH 2.2).

    Returns {"nodes_in_file", "to_gmsh_nodes", "nodes_without_element", "boundary_curves", "size_max"}. Two node counts, because they are two things:
      nodes_in_file          the nodes the saved file holds: what ElmerGrid reads ("N original nodes") and inspect_mesh counts;
      to_gmsh_nodes          the number motorgeo.to_gmsh returns, the nodes of its Gmsh session, which also has the centre point of every arc edge of the regions;
      nodes_without_element  the difference: those centre points, geometry that no element uses (87 here), which Gmsh does not save."""
    p = mg.Params(barrier="none")                        # R11: the band is left out, in Python
    regions = mg.build(p, gap="sliding")
    with intact_path():                                  # gmsh.initialize() cuts the PATH of this process
        nodes, curves = mg.to_gmsh(p, regions, str(path), size_max=size_max)
    in_file = len(mesh_nodes(path))
    return {"nodes_in_file": in_file, "to_gmsh_nodes": nodes, "nodes_without_element": nodes - in_file, "boundary_curves": curves, "size_max": size_max}


def path1_report(info: dict, regions) -> dict:
    """What path 1's mesh says against motorgeo's regions: the groups present, missing and unexpected, the area of each material against motorgeo's."""
    named = {g["name"]: g for g in info["groups"]}
    wanted = expected_groups(regions)
    present = [n for n in wanted if n in named and named[n]["elements"] > 0]
    areas = {g["name"]: g["area_mm2"] for g in info["groups"] if g["dim"] == 2}
    mesh = material_areas(areas, [r for r in regions if r.name in areas])
    geometry = geometry_areas(regions)
    return {"groups_expected": len(wanted), "groups_present": len(present), "missing": sorted(set(wanted) - set(present)),
            "unexpected": sorted(set(named) - set(wanted)), "mesh_area_mm2": mesh, "geometry_area_mm2": geometry,
            "difference_percent": {m: 100 * (mesh[m] / geometry[m] - 1) for m in geometry}}


def refinement_table(folder: Path, regions, sizes=MESH_SIZES) -> list[dict]:
    """The mesh at several size_max: nodes, triangles, worst quality and each material's area against motorgeo's. The meshes are not kept."""
    rows = []
    for size in sizes:
        scratch = Path(folder) / f"refine_{size:g}.msh"
        started = time.perf_counter()
        make_path1_mesh(scratch, size)
        seconds = time.perf_counter() - started
        info = inspect_mesh(scratch)
        scratch.unlink()
        rep = path1_report(info, regions)
        rows.append({"size_max_mm": size, "nodes": info["nodes"], "triangles": info["triangles"], "min_sicn": info["quality"]["min_sicn"],
                     "difference_percent": rep["difference_percent"], "seconds": round(seconds, 2)})
    return rows


# ============================================================================================== path 2: the CAD drawings
def arc_points(radius: float, start_deg: float, end_deg: float, sagitta: float = SAGITTA_MM, centre=(0.0, 0.0),
               min_length: float = 0.0) -> list[tuple[float, float]]:
    """Points along a counter-clockwise arc from start_deg to end_deg (a DXF arc), both ends included, never further than `sagitta` from the arc.
    With min_length > 0 the segments are not shorter than that (but there are at least two): a small arc then gets fewer, longer chords."""
    sweep = (end_deg - start_deg) % 360.0
    sweep = 360.0 if sweep == 0 else sweep
    step = 2 * math.acos(max(0.0, 1 - sagitta / radius))
    n = max(2, math.ceil(math.radians(sweep) / step))
    if min_length > 0:
        n = max(2, min(n, math.ceil(math.radians(sweep) * radius / min_length)))
    return [(centre[0] + radius * math.cos(math.radians(start_deg + sweep * i / n)), centre[1] + radius * math.sin(math.radians(start_deg + sweep * i / n)))
            for i in range(n + 1)]


def read_dxf(path, sagitta: float = SAGITTA_MM, min_length: float = 0.0) -> dict:
    """The geometry of a DXF drawing: {"lines": [LineString ...], "kinds": ["LINE" | "ARC" | "CIRCLE", one for each line], "circles": [(cx, cy, r)],
    "arcs": [(cx, cy, r, start_deg, end_deg)], "counts": {type: n}, "sagitta", "min_length", "file", "path"}. LINE, ARC and CIRCLE are geometry (arcs and circles
    are flattened to polylines within `sagitta`, with segments of at least `min_length` if that is given); anything else (the "SOLIDWORKS Educational Product"
    text) is only counted. `kinds` says which entity each polyline came from, `file` and `path` which drawing it is (the words and the record of a case name the
    drawing it read). The file is only read."""
    import ezdxf
    doc = ezdxf.readfile(str(path))
    counts: Counter = Counter()
    lines, kinds, circles, arcs = [], [], [], []
    for entity in doc.modelspace():
        kind = entity.dxftype()
        counts[kind] += 1
        if kind == "LINE":
            lines.append(LineString([(entity.dxf.start.x, entity.dxf.start.y), (entity.dxf.end.x, entity.dxf.end.y)]))
        elif kind == "ARC":
            centre, radius = (entity.dxf.center.x, entity.dxf.center.y), entity.dxf.radius
            arcs.append((centre[0], centre[1], radius, entity.dxf.start_angle, entity.dxf.end_angle))
            lines.append(LineString(arc_points(radius, entity.dxf.start_angle, entity.dxf.end_angle, sagitta, centre, min_length)))
        elif kind == "CIRCLE":
            centre, radius = (entity.dxf.center.x, entity.dxf.center.y), entity.dxf.radius
            circles.append((centre[0], centre[1], radius))
            lines.append(LineString(arc_points(radius, 0.0, 360.0, sagitta, centre, min_length)))
        else:
            continue
        kinds.append(kind)
    return {"lines": lines, "kinds": kinds, "circles": circles, "arcs": arcs, "counts": dict(counts), "sagitta": sagitta, "min_length": min_length,
            "file": Path(path).name, "path": str(path)}


def bore_mouths(dxf: dict) -> list[LineString]:
    """The arcs that close the slot openings at the bore: the tooth tops are arcs centred on the axis with the smallest radius, and between the end
    of one and the start of the next is a slot's mouth. A drawing with no such arcs has none. They are flattened like the drawing's own arcs."""
    centred = [a for a in dxf["arcs"] if abs(a[0]) < 1e-9 and abs(a[1]) < 1e-9]
    if not centred:
        return []
    r_bore = min(a[2] for a in centred)
    tops = sorted((a for a in centred if abs(a[2] - r_bore) < 1e-9), key=lambda a: a[3] % 360.0)
    mouths = []
    for k, top in enumerate(tops):
        following = tops[(k + 1) % len(tops)]
        gap = (following[3] - top[4]) % 360.0
        if 1e-6 < gap < 360.0 - 1e-6:
            mouths.append(LineString(arc_points(r_bore, top[4], following[3], dxf.get("sagitta", SAGITTA_MM), min_length=dxf.get("min_length", 0.0))))
    return mouths


def _snap(line: LineString) -> LineString:
    """Every vertex rounded to 1 nanometre (motorgeo.key): the 1e-12 gaps between entities close, and the same corner is the same point everywhere."""
    return LineString([mg.key(xy) for xy in line.coords])


def _faces(lines: list[LineString]) -> list[Polygon]:
    """The closed faces that the edges enclose (noded together, then polygonized): a face that has holes keeps them."""
    return list(polygonize(unary_union(lines)))


def rectangle_sides(face: Polygon) -> tuple[float, float]:
    """(long side, short side) of the smallest rectangle that holds the polygon [mm]: how a pocket of the rotor is told apart by its shape."""
    corners = list(face.minimum_rotated_rectangle.exterior.coords)
    long_side, short_side = sorted((math.dist(corners[0], corners[1]), math.dist(corners[1], corners[2])), reverse=True)
    return long_side, short_side


def _is_magnet_slot(face: Polygon, slot_size: tuple[float, float], tolerance: float) -> bool:
    """A pocket is a magnet slot if the smallest rectangle that holds it is slot_size (long side, short side) within `tolerance`."""
    sides = rectangle_sides(face)
    wanted = sorted(slot_size, reverse=True)
    return abs(sides[0] - wanted[0]) <= tolerance and abs(sides[1] - wanted[1]) <= tolerance


def _angle(face: Polygon) -> float:
    c = face.centroid
    return math.degrees(math.atan2(c.y, c.x)) % 360.0


def build_faces(rotor: dict, stator: dict, slot_size: tuple[float, float], tolerance: float = POCKET_TOLERANCE_MM) -> dict:
    """The faces of the two drawings, with the air band removed in memory. Returns {"rotor_steel", "shaft", "magnets", "stator_steel", "slots", "gap", "radii",
    "bands_removed", "band_lines_dropped", "slot_size", "bands", "band_line_indices", "band_edge_kinds", "band_edges_shared_with_slots", "steel_before_removal_mm2",
    "steel_as_drawn", "files"} (shapely polygons; magnets, bands and slots in the order of their angle):
      bands, bands_removed          the pockets that were taken for an air band (polygons, as drawn), and how many
      band_line_indices             which entities of the rotor drawing were dropped with them (indices into rotor["lines"]), band_lines_dropped how many,
                                    band_edge_kinds how many of each kind (LINE, ARC) when the drawing says
      band_edges_shared_with_slots  the edges that lie on a band AND on a magnet slot: they are the slot's, they stay
      steel_before_removal_mm2      the area of the rotor steel as drawn, with the band pockets cut out of it (steel_as_drawn, the polygon); rotor_steel is it
                                    after the removal, and the difference is the area of the bands
      slot_size, files              the magnet slot's size the pockets were judged by, and the names of the drawings these faces come from

    The rotor's edges are polygonized first, alone: the faces inside the rotor that are not the steel and not the shaft are pockets, and a pocket that is not a
    magnet slot (its smallest bounding rectangle is not slot_size) is a band. The edges of a band that no magnet slot shares are dropped, and everything (rotor,
    stator, the slot mouths) is polygonized together: the steel closes over the band. An edge a band and a slot have in common is kept, so the slot stays closed.
    Then the faces are sorted by where they are: inside the rotor circle, between it and the bore, outside the bore. The drawings are not written."""
    r_shaft, r_rotor = min(c[2] for c in rotor["circles"]), max(c[2] for c in rotor["circles"])
    centred = [a[2] for a in stator["arcs"] if abs(a[0]) < 1e-9 and abs(a[1]) < 1e-9]
    if not centred:
        raise ValueError("the stator drawing has no arcs centred on the axis (the bore): the air gap and the slots cannot be told apart")
    r_bore, r_outer = min(centred), max(c[2] for c in stator["circles"])
    rotor_lines = [_snap(l) for l in rotor["lines"]]
    first = _faces(rotor_lines)
    steel = max(first, key=lambda f: f.area)
    pockets = [f for f in first if f is not steel and not f.contains(Point(0, 0))]
    bands = sorted((f for f in pockets if not _is_magnet_slot(f, slot_size, tolerance)), key=_angle)
    magnet_slots = [f for f in pockets if f not in bands]
    dropped: list[int] = []
    shared = 0
    if bands:
        zones = lambda polygons: [prep(polygon.boundary.buffer(EDGE_ZONE_MM)) for polygon in polygons]
        on = lambda zone_list, line: any(zone.covers(line) for zone in zone_list)
        band_zones, slot_zones = zones(bands), zones(magnet_slots)
        for i, line in enumerate(rotor_lines):
            on_band, on_slot = on(band_zones, line), on(slot_zones, line)
            if on_band and on_slot:
                shared += 1                                        # a wall the band and a slot have in common is the slot's: it stays
            elif on_band:
                dropped.append(i)
        gone = set(dropped)
        rotor_lines = [line for i, line in enumerate(rotor_lines) if i not in gone]
    everything = rotor_lines + [_snap(l) for l in stator["lines"]] + [_snap(m) for m in bore_mouths(stator)]
    faces = _faces(everything)
    radius = lambda f: math.hypot(f.representative_point().x, f.representative_point().y)
    inside = [f for f in faces if radius(f) < r_rotor]
    shaft = next(f for f in inside if f.contains(Point(0, 0)))
    rotor_faces = [f for f in inside if f is not shaft]
    rotor_steel = max(rotor_faces, key=lambda f: f.area)
    between = [f for f in faces if r_rotor <= radius(f) < r_bore]
    outside = [f for f in faces if radius(f) >= r_bore]
    stator_steel = max(outside, key=lambda f: f.area)
    if len(between) != 1:
        raise ValueError(f"expected one air-gap face between radius {r_rotor} and {r_bore}, found {len(between)}")
    return {"rotor_steel": rotor_steel, "shaft": shaft, "magnets": sorted((f for f in rotor_faces if f is not rotor_steel), key=_angle),
            "stator_steel": stator_steel, "slots": sorted((f for f in outside if f is not stator_steel), key=_angle), "gap": between[0],
            "radii": {"shaft": r_shaft, "rotor": r_rotor, "bore": r_bore, "stator_outer": r_outer}, "bands_removed": len(bands), "band_lines_dropped": len(dropped),
            "slot_size": slot_size, "bands": bands, "band_line_indices": dropped,
            "band_edge_kinds": dict(Counter(rotor["kinds"][i] for i in dropped)) if "kinds" in rotor else {}, "band_edges_shared_with_slots": shared,
            "steel_before_removal_mm2": steel.area, "steel_as_drawn": steel, "files": {"rotor": rotor.get("file", "the rotor drawing"), "stator": stator.get("file", "the stator drawing")}}


def cad_group_names(faces: dict) -> list[str]:
    """The physical groups of the CAD mesh, in the order they are made: the faces, then the two A = 0 circles."""
    return (["rotor_steel"] + [f"magnet_{i + 1}" for i in range(len(faces["magnets"]))] + ["stator_steel"]
            + [f"slot_{i + 1}" for i in range(len(faces["slots"]))] + ["gap", "shaft", "outer"])


def mesh_faces(faces: dict, path, size_max: float = SIZE_MAX) -> dict:
    """Mesh the faces with Gmsh's built-in kernel and write `path` (MSH 2.2). One physical group per face (cad_group_names), and the curves on the
    shaft's circle and on the outer circle as the groups "shaft" and "outer". Two faces that meet along an edge use the same Gmsh line, so the mesh
    is conforming by construction.

    The element size comes from a size field, not from sizes at the points: the corners of a flattened circle are 0.3 to 1.3 mm apart, and sizes
    taken from them would keep the whole rotor at that size. The field is the smaller of two thresholds on the distance to the faces' edges: the air gap
    (0.3 mm at the gap, growing to size_max 10 mm away) and the magnet pockets (0.5 mm, growing to size_max 6 mm away). Returns {"nodes", "lines", "points"}."""
    import gmsh
    names = cad_group_names(faces)
    polygons = [("rotor_steel", faces["rotor_steel"])] + [(f"magnet_{i + 1}", f) for i, f in enumerate(faces["magnets"])] + [("stator_steel", faces["stator_steel"])] \
        + [(f"slot_{i + 1}", f) for i, f in enumerate(faces["slots"])] + [("gap", faces["gap"])]
    r_shaft, r_outer = faces["radii"]["shaft"], faces["radii"]["stator_outer"]
    with intact_path():
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("cad")
            geo = gmsh.model.geo
            points: dict[tuple, int] = {}
            curves: dict[tuple, int] = {}
            ends: dict[int, tuple] = {}
            rings: dict[str, list[int]] = {}

            def point(xy) -> int:
                key = mg.key(xy)
                if key not in points:
                    points[key] = geo.addPoint(key[0], key[1], 0.0, size_max)
                return points[key]

            def curve(a, b) -> int:                        # the signed tag of the line from a to b: one line for both of the faces that share it
                ka, kb = mg.key(a), mg.key(b)
                if (ka, kb) in curves:
                    return curves[(ka, kb)]
                if (kb, ka) in curves:
                    return -curves[(kb, ka)]
                curves[(ka, kb)] = tag = geo.addLine(point(a), point(b))
                ends[tag] = (ka, kb)
                return tag

            def loop(name, ring) -> int:
                coords = list(ring.coords)
                tags = [curve(a, b) for a, b in zip(coords[:-1], coords[1:])]
                rings.setdefault(name, []).extend(abs(t) for t in tags)
                return geo.addCurveLoop(tags)

            surfaces = {name: geo.addPlaneSurface([loop(name, poly.exterior), *[loop(name, h) for h in poly.interiors]]) for name, poly in polygons}
            geo.synchronize()
            for name in names:
                if name in surfaces:
                    gmsh.model.addPhysicalGroup(2, [surfaces[name]], name=name)
            on = lambda tag, radius: all(abs(math.hypot(*xy) - radius) < 1e-3 for xy in ends[tag])
            for name, radius in (("shaft", r_shaft), ("outer", r_outer)):
                gmsh.model.addPhysicalGroup(1, [t for t in ends if on(t, radius)], name=name)
            field = gmsh.model.mesh.field

            def threshold(curve_tags, size_min, dist_min, dist_max) -> int:
                distance = field.add("Distance")
                field.setNumbers(distance, "CurvesList", curve_tags)
                field.setNumber(distance, "Sampling", 50)
                limit = field.add("Threshold")
                for key, value in (("InField", distance), ("SizeMin", size_min), ("SizeMax", size_max), ("DistMin", dist_min), ("DistMax", dist_max)):
                    field.setNumber(limit, key, value)
                return limit

            gap_curves = sorted(set(rings["gap"]))
            pocket_curves = sorted({t for name in rings if name.startswith("magnet_") for t in rings[name]})
            smallest = field.add("Min")
            field.setNumbers(smallest, "FieldsList", [threshold(gap_curves, min(0.3, size_max), 1.0, 10.0), threshold(pocket_curves, min(0.5, size_max), 0.5, 6.0)])
            field.setAsBackgroundMesh(smallest)
            for option in ("Mesh.MeshSizeExtendFromBoundary", "Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature"):
                gmsh.option.setNumber(option, 0)
            gmsh.option.setNumber("Mesh.MeshSizeMax", size_max)
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            gmsh.model.mesh.generate(2)
            gmsh.write(str(path))
            nodes = len(gmsh.model.mesh.getNodes()[0])
        finally:
            gmsh.finalize()
    return {"nodes": nodes, "lines": len(curves), "points": len(points)}


def cad_offset(cad_points, ref_points, period_deg: float) -> tuple[float, float]:
    """The rotation, in degrees, that carries the drawing's points onto the model's, and the largest distance that is left (mm).

    The patterns repeat every `period_deg` (45 for the rotor's sixteen slots, 10 for the stator's 36), so the rotation is only known modulo that:
    it is searched over [0, period) in 0.01-degree steps, and the best one is refined by the least-squares rotation of the matched pairs."""
    cad, ref = np.asarray(cad_points, float), np.asarray(ref_points, float)
    turn = lambda pts, deg: pts @ np.array([[math.cos(math.radians(deg)), -math.sin(math.radians(deg))],
                                           [math.sin(math.radians(deg)), math.cos(math.radians(deg))]]).T
    nearest = lambda moved: np.linalg.norm(moved[:, None, :] - ref[None, :, :], axis=2)
    best = min(np.arange(0.0, period_deg, 0.01), key=lambda deg: nearest(turn(cad, deg)).min(axis=1).max())
    matched = ref[nearest(turn(cad, best)).argmin(axis=1)]
    base = turn(cad, best)
    cross = float((base[:, 0] * matched[:, 1] - base[:, 1] * matched[:, 0]).sum())
    dot = float((base * matched).sum())
    offset = (best + math.degrees(math.atan2(cross, dot))) % period_deg
    residual = float(nearest(turn(cad, offset)).min(axis=1).max())
    return float(offset), residual


# ============================================================================================== the case
def render(folder: Path, msh_name: str, png_name: str, log: Log) -> dict:
    """Gmsh's own window renders `msh_name` into `png_name` (a .geo that Merges the mesh, Draws, Prints and Exits), through run_solver. Returns the
    record of the render: whether the PNG is there, how long it took."""
    folder = Path(folder)
    geo = folder / (Path(png_name).stem + "_shot.geo")
    geo.write_text(SHOT_GEO.format(msh=msh_name, png=png_name), encoding="utf-8", newline="\n")
    started = time.perf_counter()
    result = run_program(log, [paths.tool("gmsh"), geo.name], label=f"bench B2 Gmsh render {png_name}", cwd=folder, timeout_s=120, threads=1,
                         name=Path(png_name).stem + "_gmsh")
    return {"png": png_name, "made": (folder / png_name).is_file(), "returncode": result.returncode, "timed_out": result.timed_out,
            "seconds": round(time.perf_counter() - started, 2)}


def _crop_white(image: np.ndarray, pad: int = 14) -> np.ndarray:
    """The picture without its white margin (Gmsh centres the model in a square window)."""
    ink = (image[..., :3] < 0.97).any(axis=2)
    rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
    return image[max(rows[0] - pad, 0):rows[-1] + pad, max(cols[0] - pad, 0):cols[-1] + pad]


def figure(folder: Path, left: str, right: str, captions: tuple[str, str], target: Path) -> None:
    """The two Gmsh renders side by side: path 1's sector and path 2's whole machine, each cropped to the model. The figure is saved with a tight box: the
    pictures go into the project notes, and a blank band under them (the room the captions were given) is cropped away."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    images = [_crop_white(plt.imread(str(Path(folder) / name))) for name in (left, right)]
    widths = [im.shape[1] for im in images]
    fig, axes = plt.subplots(1, 2, figsize=(15, 15 * max(im.shape[0] / sum(widths) for im in images) + 1.2), gridspec_kw={"width_ratios": widths})
    for ax, image, caption in zip(axes, images, captions):
        ax.imshow(image)
        ax.axis("off")
        ax.set_title(caption, fontsize=10.5, loc="left", color="#1d1d1f")
    fig.tight_layout()
    fig.savefig(target, dpi=110, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def run(folder: Path, log: Log) -> dict:
    """Run B2 in `folder` (already empty). Returns the results dict whose "checks" list every verdict."""
    folder = Path(folder)
    p = mg.Params(barrier="none")
    regions = mg.build(p, gap="sliding")
    started = time.perf_counter()
    checks: list[dict] = []
    results: dict = {"case": CASE, "title": TITLE, "gmsh_version": None, "params": params(p)}
    import gmsh
    results["gmsh_version"] = gmsh.__version__
    # ---- path 1
    t = time.perf_counter()
    made = make_path1_mesh(folder / "path1.msh")
    info1 = inspect_mesh(folder / "path1.msh")
    report1 = path1_report(info1, regions)
    table = refinement_table(folder, regions)
    chords1 = fillet_chords(mesh_nodes(folder / "path1.msh"), model_fillets(regions, p))
    path1_seconds = time.perf_counter() - t
    log.say(f"B2 path 1: motorgeo.to_gmsh, {made['nodes_in_file']:,} nodes in the saved file ({made['to_gmsh_nodes']:,} counted by to_gmsh, {made['nodes_without_element']} of them "
            f"arc centres that no element uses), {info1['triangles']:,} triangles, {len(info1['groups'])} physical groups, "
            f"min SICN {info1['quality']['min_sicn']:.3f} in {info1['quality']['worst']['group']} at "
            f"({info1['quality']['worst']['x_mm']:.2f}, {info1['quality']['worst']['y_mm']:.2f}) mm; {path1_seconds:.1f} s with the size table")
    results["path1"] = {"mesh": made, "inspection": info1, "report": report1, "refinement": table, "fillet_chords": chords1, "seconds": round(path1_seconds, 2)}
    checks += path1_checks(info1, report1, table, regions, p, chords1)
    # ---- path 2: SynRM rotor.DXF (the drawing that has the air band: it is removed here, in memory) and the stator; Rotor.DXF (no band) is read for the comparison
    t = time.perf_counter()
    rotor_path, earlier_path, stator_path = paths.path("dxf_rotor_with_band"), paths.path("dxf_rotor"), paths.path("dxf_stator")
    hashes_before = file_hashes([rotor_path, earlier_path, stator_path])                                       # the drawings are read only: hashed before and after
    slot = (p.slot_w, p.slot_t)
    rotor, earlier_rotor, stator = read_dxf(rotor_path), read_dxf(earlier_path), read_dxf(stator_path)         # the geometry: arcs flattened within 0.25 micrometre
    exact = build_faces(rotor, stator, slot_size=slot)
    earlier = build_faces(earlier_rotor, stator, slot_size=slot)
    meshed = build_faces(read_dxf(rotor_path, MESH_SAGITTA_MM, MIN_SEGMENT_MM), read_dxf(stator_path, MESH_SAGITTA_MM, MIN_SEGMENT_MM),
                         slot_size=slot)                                                                       # the mesh: coarser chords, so that the mesh can be coarser
    made2 = mesh_faces(meshed, folder / "path2.msh")
    info2 = inspect_mesh(folder / "path2.msh")
    chords2 = fillet_chords(mesh_nodes(folder / "path2.msh"), drawing_fillets(rotor, p))
    path2_seconds = time.perf_counter() - t
    log.say(path2_log_line(info2, exact, earlier, earlier_rotor, rotor, stator, path2_seconds))
    offsets = cad_offsets(exact, p, regions)
    model_bands = motorgeo_bands(p)
    results["path2"] = {"files": path2_files(rotor, earlier_rotor, stator, hashes_before),
                        "entities": {"rotor": rotor["counts"], "rotor_earlier_export": earlier_rotor["counts"], "stator": stator["counts"],
                                     "ignored": {"rotor": {k: n for k, n in rotor["counts"].items() if k not in GEOMETRY_KINDS},
                                                 "stator": {k: n for k, n in stator["counts"].items() if k not in GEOMETRY_KINDS}},
                                     "band_lines_dropped": exact["band_lines_dropped"]},
                        "bands": {"drawing": rotor["file"], "found_and_removed": exact["bands_removed"], "expected_from_motorgeo": model_bands["count"],
                                  "edges_dropped": exact["band_lines_dropped"], "edge_kinds_dropped": exact["band_edge_kinds"], "edge_kinds_expected": model_bands["edges"],
                                  "edges_shared_with_slots_kept": exact["band_edges_shared_with_slots"], "area_mm2_each": [round(b.area, 4) for b in exact["bands"]],
                                  "steel_before_removal_mm2": exact["steel_before_removal_mm2"], "steel_after_removal_mm2": exact["rotor_steel"].area,
                                  "steel_added_mm2": exact["rotor_steel"].area - exact["steel_before_removal_mm2"], "sum_of_band_areas_mm2": sum(b.area for b in exact["bands"])},
                        "earlier_export": {"drawing": earlier_rotor["file"], "bands_found": earlier["bands_removed"], "steel_mm2": earlier["rotor_steel"].area,
                                           "steel_difference_mm2": exact["rotor_steel"].area - earlier["rotor_steel"].area, "limit_mm2": flattening_area_limit(earlier_rotor),
                                           "faces": face_summary(earlier)},
                        "radii_mm": exact["radii"], "bands_removed": exact["bands_removed"],
                        "faces": face_summary(exact), "faces_as_meshed": face_summary(meshed),
                        "flattening": {"geometry_sagitta_mm": SAGITTA_MM, "mesh_sagitta_mm": MESH_SAGITTA_MM, "mesh_min_segment_mm": MIN_SEGMENT_MM},
                        "mesh": made2, "inspection": info2, "fillet_chords": chords2, "offsets": offsets, "seconds": round(path2_seconds, 2)}
    checks += path2_checks(exact, meshed, info2, offsets, p, regions, chords2, rotor, stator, earlier, earlier_rotor)
    # ---- the pictures: the removal, from Python (bands.png), and Gmsh's own renders of the meshes (figure.png)
    fig_bands(exact, bands_title(exact, earlier, rotor, earlier_rotor), folder / "bands.png")
    renders = [render(folder, "path1.msh", "path1.png", log), render(folder, "path2.msh", "path2.png", log)]
    results["renders"] = renders
    if all(r["made"] for r in renders):
        figure(folder, "path1.png", "path2.png",
               (f"Path 1: motorgeo.to_gmsh, the band-free 90-degree sector\n{info1['nodes']:,} nodes, {info1['triangles']:,} triangles, minimum SICN {info1['quality']['min_sicn']:.2f}",
                path2_caption(info2, exact, rotor, stator, earlier, earlier_rotor)), folder / "figure.png")
    else:
        checks.append(check("Gmsh's own render of the meshes", False, "Gmsh draws both meshes into PNG files", f"renders: {renders}", "gmsh.exe with Draw; Print"))
    checks.append(read_only_check(hashes_before, file_hashes([rotor_path, earlier_path, stator_path])))        # after everything this case did
    results["checks"] = checks
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    return results


def face_summary(faces: dict) -> dict:
    """The sizes of the faces for results.json: counts and areas."""
    return {"rotor_steel_mm2": faces["rotor_steel"].area, "rotor_steel_holes": len(faces["rotor_steel"].interiors), "shaft_mm2": faces["shaft"].area,
            "magnets": len(faces["magnets"]), "magnet_mm2_each": [round(m.area, 4) for m in faces["magnets"]], "stator_steel_mm2": faces["stator_steel"].area,
            "slots": len(faces["slots"]), "slot_mm2_each": [round(s.area, 4) for s in faces["slots"]], "gap_mm2": faces["gap"].area}


def cad_offsets(faces: dict, p: mg.Params, regions) -> dict:
    """How the drawing sits against motorgeo's model: the rotation of its sixteen magnet slots onto motorgeo's (modulo the 45-degree pole pitch) and of
    its 36 stator slots onto motorgeo's (modulo the 10-degree slot pitch), and the largest distance that is left."""
    magnets = [r.label for r in regions if r.name.startswith("magnet_")]
    model_magnets = [mg.rot(pt, i * p.sector_deg) for i in range(p.n_sectors) for pt in magnets]
    drawn_magnets = [(f.centroid.x, f.centroid.y) for f in faces["magnets"]]
    rotor_deg, rotor_residual = cad_offset(drawn_magnets, model_magnets, 360.0 / (2 * p.pole_pairs))
    unit = lambda angles: [(math.cos(math.radians(a)), math.sin(math.radians(a))) for a in angles]     # slots by their angle: points on the unit circle
    stator_deg, stator_residual = cad_offset(unit([_angle(f) for f in faces["slots"]]), unit([mg.slot_angle(p, k) for k in range(p.slots)]), 360.0 / p.slots)
    return {"rotor": {"turn_deg_mod_45": rotor_deg, "max_distance_mm": rotor_residual, "points": len(drawn_magnets)},
            "stator": {"turn_deg_mod_10": stator_deg, "max_distance_deg": math.degrees(stator_residual), "points": len(faces["slots"])}}


def motorgeo_bands(p: mg.Params) -> dict:
    """motorgeo's air band, the one the CAD still draws behind each magnet, turned into the whole machine: how many to expect and what one is.

    The count is not typed here. It is the number of `barrier_` regions of motorgeo.build(Params(barrier="air")) in the sector (4: two poles x two magnets) times the
    sectors in the machine (4), so it follows the parameters. Returns {"per_sector", "sectors", "count", "sector_deg", "outline", "edges", "bands", "slots"}:
      outline   the edges of one band's outline by kind ({"LINE": 5, "ARC": 1}), edges the same for all the bands of the machine ({"LINE": 80, "ARC": 16});
      bands     the 16 band polygons, slots the 16 magnet-slot polygons, each turned by whole sectors (motorgeo's frame, not the drawing's)."""
    air = dataclasses.replace(p, barrier="air")
    regions = mg.build(air, gap="sliding")
    visible = mg.visible(regions)
    outlines = [r for r in regions if r.name.startswith("barrier_")]
    magnets = [r for r in regions if r.name.startswith("magnet_")]
    turn = lambda polygon, k: affinity.rotate(polygon, k * air.sector_deg, origin=(0.0, 0.0))
    per_band = [Counter("ARC" if isinstance(e, mg.Arc) else "LINE" for e in r.edges) for r in outlines]
    edges: Counter = Counter()
    for counter in per_band:
        edges.update({kind: n * air.n_sectors for kind, n in counter.items()})
    return {"per_sector": len(outlines), "sectors": air.n_sectors, "count": len(outlines) * air.n_sectors, "sector_deg": air.sector_deg,
            "outline": dict(per_band[0]) if per_band else {}, "edges": dict(edges),
            "bands": [turn(visible[r.name], k) for k in range(air.n_sectors) for r in outlines],
            "slots": [turn(visible[r.name], k) for k in range(air.n_sectors) for r in magnets]}


def flattening_area_limit(dxf: dict) -> float:
    """How far two readings of the same drawing's steel can differ in area [mm2], from the flattening of its arcs: (2/3) x the sagitta x the length of its arcs and circles.

    An arc is read as chords that never leave it by more than the sagitta s. A chord of length c cuts a circular segment off the arc; for a segment as flat as these
    its area is 2/3 of the rectangle c x s. The chords of an arc of length L together cut off at most (2/3) s L, and two readings of the same arcs, each within s, differ
    by at most that, summed over every arc and circle of the drawing. With s = 0.25 micrometre and the 360 mm of arc and circle in Rotor.DXF it is 0.060 mm2, a
    relative 1e-5 of the steel. It is the bench's own limit (the brief gives none for one drawing against another)."""
    length = sum(r * math.radians((end - start) % 360.0 or 360.0) for _, _, r, start, end in dxf["arcs"]) + sum(2 * math.pi * r for _, _, r in dxf["circles"])
    return 2.0 / 3.0 * dxf["sagitta"] * length


def band_check(faces: dict, p: mg.Params) -> dict:
    """Does the rotor of the drawing come out with no air band, only the magnet slots? The words are computed from what build_faces found (how many bands it
    removed, what the pockets left are), never from what any comment says about the drawing. The number of magnet slots to expect is motorgeo's (motorgeo_bands)."""
    width, thickness = faces["slot_size"]
    slots = [_is_magnet_slot(m, faces["slot_size"], POCKET_TOLERANCE_MM) for m in faces["magnets"]]
    wanted = len(motorgeo_bands(p)["slots"])
    left = "all magnet slots" if all(slots) else f"only {sum(slots)} of them magnet slots"
    return check("path 2: no air band remains in the rotor", len(slots) == wanted and all(slots),
                 f"every pocket of the rotor that is not a magnet slot is merged into the steel; what is left are the {wanted} magnet slots (motorgeo's count)",
                 f"{faces['bands_removed']} band(s) found and removed from {faces['files']['rotor']}; {len(slots)} pockets left, {left} ({width:g} x {thickness:g} mm)",
                 "build_faces(): the smallest bounding rectangle of each pocket")


def band_removal_check(faces: dict, p: mg.Params, turn_deg: float) -> dict:
    """Is what build_faces removed motorgeo's air band, all of it and nothing else? `turn_deg` is the rotation that carries the drawing onto motorgeo's frame, found
    from the magnet slots (cad_offsets), so that a band is held to where it should sit.

    Against motorgeo (motorgeo_bands): the number found is the number to expect; each removed band, matched to the motorgeo band its turned centroid is nearest to
    (all different), has its smallest bounding rectangle within POCKET_TOLERANCE_MM and its area within CAD_AREA_TOLERANCE_PERCENT, its centroid radius, its place and its
    distance to the nearest magnet slot (a web of steel: beside a slot) within CAD_POSITION_TOLERANCE_MM; and the entities dropped with the bands are, by kind, the outlines
    of that many motorgeo bands (5 LINE + 1 ARC each)."""
    model = motorgeo_bands(p)
    drawn, slots, model_slots = faces["bands"], faces["magnets"], model["slots"]
    kinds = lambda counts: " + ".join(f"{counts[k]} {k}" for k in GEOMETRY_KINDS if counts.get(k)) or "none"
    worst = {"rectangle_mm": 0.0, "area_percent": 0.0, "radius_mm": 0.0, "position_mm": 0.0, "slot_distance_mm": 0.0}
    paired, example = [], None
    for band in drawn:
        turned = affinity.rotate(band, turn_deg, origin=(0.0, 0.0))
        j = int(np.argmin([math.hypot(turned.centroid.x - m.centroid.x, turned.centroid.y - m.centroid.y) for m in model["bands"]]))
        paired.append(j)
        m = model["bands"][j]
        sides, m_sides = rectangle_sides(band), rectangle_sides(m)
        radius, m_radius = math.hypot(band.centroid.x, band.centroid.y), math.hypot(m.centroid.x, m.centroid.y)
        slot_distance = min((band.distance(s) for s in slots), default=float("inf"))
        m_slot_distance = min(m.distance(s) for s in model_slots)
        worst["rectangle_mm"] = max(worst["rectangle_mm"], abs(sides[0] - m_sides[0]), abs(sides[1] - m_sides[1]))
        worst["area_percent"] = max(worst["area_percent"], abs(100 * (band.area / m.area - 1)))
        worst["radius_mm"] = max(worst["radius_mm"], abs(radius - m_radius))
        worst["position_mm"] = max(worst["position_mm"], math.hypot(turned.centroid.x - m.centroid.x, turned.centroid.y - m.centroid.y))
        worst["slot_distance_mm"] = max(worst["slot_distance_mm"], abs(slot_distance - m_slot_distance) if math.isfinite(slot_distance) else float("inf"))
        example = example or (sides, band.area, radius, slot_distance, m_sides, m.area, m_radius, m_slot_distance)
    edges_dropped = dict(faces["band_edge_kinds"])
    ok = (len(drawn) == model["count"] and len(set(paired)) == len(drawn) and edges_dropped == model["edges"]
          and worst["rectangle_mm"] <= POCKET_TOLERANCE_MM and worst["area_percent"] <= CAD_AREA_TOLERANCE_PERCENT
          and max(worst["radius_mm"], worst["position_mm"], worst["slot_distance_mm"]) <= CAD_POSITION_TOLERANCE_MM)
    one = (f"; one band: {example[0][0]:.3f} x {example[0][1]:.3f} mm, {example[1]:.4f} mm2, centroid at radius {example[2]:.3f} mm, {example[3]:.3f} mm from the nearest slot "
           f"(motorgeo's: {example[4][0]:.3f} x {example[4][1]:.3f} mm, {example[5]:.4f} mm2, radius {example[6]:.3f} mm, {example[7]:.3f} mm)") if example else ""
    return check(f"path 2: the air bands of {faces['files']['rotor']} are found and removed, and each is motorgeo's band", ok,
                 f"the machine has {model['count']} air bands, counted from motorgeo ({model['per_sector']} band outlines of Params(barrier=\"air\") in the "
                 f"{model['sector_deg']:g}-degree sector x {model['sectors']} sectors); each pocket of the rotor that is not a magnet slot is removed, and each is motorgeo's band: "
                 f"its rectangle within {POCKET_TOLERANCE_MM:g} mm and its area within {CAD_AREA_TOLERANCE_PERCENT:g} % of motorgeo's, its centroid radius, its place after the "
                 f"slots' turn and its distance to the nearest slot within {CAD_POSITION_TOLERANCE_MM:g} mm; the edges dropped are the outlines of those bands ({OWN_LIMIT})",
                 f"{len(drawn)} found and removed from {faces['files']['rotor']} (motorgeo's Params(barrier=\"air\"): {model['per_sector']} band outlines in the "
                 f"{model['sector_deg']:g}-degree sector x {model['sectors']} sectors = {model['count']}); {faces['band_lines_dropped']} edges dropped = {kinds(edges_dropped)} "
                 f"(motorgeo's band outline is {kinds(model['outline'])}, so {model['count']} bands are {kinds(model['edges'])}); against motorgeo's band, the largest "
                 f"difference over the {len(drawn)}: rectangle {worst['rectangle_mm']:.4f} mm, area {worst['area_percent']:.4f} %, centroid radius {worst['radius_mm']:.4f} mm, "
                 f"place {worst['position_mm']:.4f} mm after the slots' turn of {turn_deg:.3f} degrees, distance to the nearest slot {worst['slot_distance_mm']:.4f} mm" + one,
                 "build_faces(): the pockets that are not magnet slots; motorgeo.build(Params(barrier=\"air\")): its barrier regions, turned into the whole machine")


def entity_check(faces: dict, rotor: dict, stator: dict) -> dict:
    """Did every LINE, ARC and CIRCLE of the two drawings become an edge of a face? read_dxf keeps those three kinds as geometry and only COUNTS the others (the
    text in the drawing); this holds the geometry to account: each entity's polyline has to lie within EDGE_ZONE_MM of the faces' boundaries, except the edges that
    build_faces dropped on purpose with an air band. Those must be EXACTLY the edges that are not on a face (as sets of entities, not counts: one band edge left in the
    steel would dangle, and one dropped slot edge would be on a face), all of them on a removed band's outline, none on a magnet slot's. A wall a band and a slot share
    is the slot's and is kept: it is counted as such. What is ignored is listed."""
    polygons = [faces["rotor_steel"], faces["shaft"], faces["stator_steel"], faces["gap"], *faces["magnets"], *faces["slots"]]
    zone = prep(unary_union([polygon.boundary for polygon in polygons]).buffer(EDGE_ZONE_MM))
    slot_zone = prep(unary_union([m.boundary for m in faces["magnets"]]).buffer(EDGE_ZONE_MM))
    band_zone = prep(unary_union([b.boundary for b in faces["bands"]]).buffer(EDGE_ZONE_MM)) if faces["bands"] else None
    kinds = lambda counts: " + ".join(f"{counts[k]} {k}" for k in GEOMETRY_KINDS if counts.get(k))
    ignored = lambda dxf: ", ".join(f"{n} {k}" for k, n in dxf["counts"].items() if k not in GEOMETRY_KINDS) or "nothing"
    parts, ok = [], True
    for label, dxf, dropped in (("rotor", rotor, set(faces["band_line_indices"])), ("stator", stator, set())):
        lines = dxf["lines"]
        entities = sum(dxf["counts"].get(k, 0) for k in GEOMETRY_KINDS)
        stray = {i for i, line in enumerate(lines) if not zone.covers(line)}
        ok = ok and len(lines) == entities and stray == dropped
        unaccounted = stray - dropped                              # on no face, and not dropped on purpose: a band edge left in the steel, or any other stray line
        text = f"{label} {kinds(dxf['counts'])} ({dxf.get('file', 'the drawing')}): {len(lines) - len(stray)} are edges of a face"
        if label == "rotor":
            on_outline = sum(1 for i in dropped if band_zone is not None and band_zone.covers(lines[i]))
            on_slot = sum(1 for i in dropped if slot_zone.covers(lines[i]))
            ok = ok and on_outline == len(dropped) and on_slot == 0
            by_kind = Counter(dxf["kinds"][i] for i in dropped) if "kinds" in dxf else Counter()
            text += (f", {len(dropped)} dropped with the {faces['bands_removed']} air band(s)" + (f" ({kinds(by_kind)})" if dropped else "") + f", {len(unaccounted)} not an edge of any face and not dropped; "
                     f"of the dropped, {on_outline} lie on a removed band's outline, {on_slot} on a magnet slot's; {faces['band_edges_shared_with_slots']} band edge(s) shared "
                     f"with a magnet slot kept")
        else:
            text += f", {len(unaccounted)} not an edge of any face"
        parts.append(text)
    return check("path 2: every LINE, ARC and CIRCLE of the two drawings is an edge of a face", ok,
                 f"every geometry entity of both drawings is on the boundary of a face (within {EDGE_ZONE_MM:g} mm: {OWN_LIMIT}), except the edges dropped on purpose with an air band, "
                 "which are exactly the edges that are on no face, all on a removed band's outline and none on a magnet slot's; the text is ignored",
                 "; ".join(parts) + f"; ignored (not geometry): rotor {ignored(rotor)}, stator {ignored(stator)}",
                 "read_dxf(), build_faces(): the faces' boundaries")


def earlier_export_check(faces: dict, earlier: dict, rotor: dict, earlier_rotor: dict) -> dict:
    """After the removal, is the band drawing's rotor the earlier export's? `faces` are build_faces' of the band drawing (SynRM rotor.DXF, `rotor`), `earlier` those of
    the CAD's earlier export of the same rotor (Rotor.DXF, `earlier_rotor`), which never had the band: an independent answer to 'was the band removed, and nothing else?'.

    The same magnet slots at the same places (every centroid within CAD_POSITION_TOLERANCE_MM, the total area within the flattening limit), the same shaft and outer
    circle (centre and radius within CAD_POSITION_TOLERANCE_MM), the steel area equal within flattening_area_limit(), and the balance: the steel the removal added is the
    area of the bands removed, which is what Rotor.DXF's steel has over the band drawing's as drawn."""
    name, earlier_name = rotor.get("file", "the band drawing"), earlier_rotor.get("file", "the earlier export")
    slots, earlier_slots = faces["magnets"], earlier["magnets"]
    slot_shift = max((math.hypot(a.centroid.x - b.centroid.x, a.centroid.y - b.centroid.y) for a, b in zip(slots, earlier_slots)), default=0.0)
    slot_area = abs(sum(a.area for a in slots) - sum(b.area for b in earlier_slots))
    circles, earlier_circles = sorted(rotor["circles"], key=lambda c: c[2]), sorted(earlier_rotor["circles"], key=lambda c: c[2])
    circle_shift = max((max(abs(x - y) for x, y in zip(a, b)) for a, b in zip(circles, earlier_circles)), default=0.0)
    limit = flattening_area_limit(earlier_rotor)
    steel, earlier_steel = faces["rotor_steel"].area, earlier["rotor_steel"].area
    bands_area = sum(b.area for b in faces["bands"])
    added = steel - faces["steel_before_removal_mm2"]
    gap = earlier_steel - faces["steel_before_removal_mm2"]
    arc_length = 1.5 * limit / earlier_rotor["sagitta"]
    ok = (len(slots) == len(earlier_slots) and len(circles) == len(earlier_circles) and slot_shift <= CAD_POSITION_TOLERANCE_MM and circle_shift <= CAD_POSITION_TOLERANCE_MM
          and slot_area <= limit and abs(steel - earlier_steel) <= limit and abs(added - bands_area) <= limit and abs(gap - bands_area) <= limit)
    return check(f"path 2: after the removal, {name} is {earlier_name}", ok,
                 f"{earlier_name} is the CAD's earlier export of the same rotor and has no band: the same {len(earlier_slots)} magnet slots at the same places (every centroid within "
                 f"{CAD_POSITION_TOLERANCE_MM:g} mm), the same shaft and outer circle (within {CAD_POSITION_TOLERANCE_MM:g} mm), the slots' area and the steel area equal within "
                 f"{limit:.4f} mm2 = 2/3 x the sagitta {earlier_rotor['sagitta'] * 1e3:g} micrometre x the {arc_length:.1f} mm of arc and circle that bound them (the most two "
                 f"flattenings of the same arcs can differ by), and the steel the removal added = the area of the bands removed ({OWN_LIMIT})",
                 f"{len(slots)} slots, the largest centroid distance {slot_shift:.4f} mm, their area differs by {slot_area:.4f} mm2; circles (shaft, outer) differ by at most "
                 f"{circle_shift:.4f} mm; steel {steel:.4f} against {earlier_steel:.4f} mm2, difference {steel - earlier_steel:+.4f} mm2 of {limit:.4f}; the removal added "
                 f"{added:.4f} mm2 to the steel, the {len(faces['bands'])} bands removed are {bands_area:.4f} mm2, and {earlier_name}'s steel has {gap:.4f} mm2 more than "
                 f"{name}'s as drawn; {earlier['bands_removed']} band(s) found in {earlier_name}",
                 f"build_faces() of both drawings; flattening_area_limit() ({SAGITTA_MM * 1e3:g} micrometre sagitta)")


def file_hashes(files) -> dict[str, str]:
    """{file name: SHA-256 of its bytes} for each of `files`: the proof that a read-only drawing was not written."""
    return {Path(f).name: hashlib.sha256(Path(f).read_bytes()).hexdigest() for f in files}


def read_only_check(before: dict, after: dict) -> dict:
    """The drawings are read only: the hashes taken before B2 read them and after everything it did are the same, every drawing is still there."""
    changed = sorted(name for name in before if after.get(name) != before[name])
    return check("path 2: the CAD drawings are read only: byte for byte what they were before the run", not changed,
                 "the SHA-256 of each drawing is the same after the run as before it (the bench only reads reference\\)",
                 ("unchanged: " + ", ".join(f"{name} {digest[:12]}" for name, digest in before.items())) if not changed else "changed or gone: " + ", ".join(changed),
                 "hashlib, before path 2 and after its pictures")


def fig_bands(faces: dict, title: str, target) -> None:
    """bands.png: the band drawing's rotor faces from Python, so that the removal in memory can be seen at a glance. The steel as drawn (grey, with a hole for every pocket
    and the shaft), the removed air bands in a second colour (orange, where the band pockets are: the steel closes over them), the magnet slots that stay in a third
    (blue). `title` may have several lines (bench.wrap_title wraps them to the figure). Saved with a tight box, so there is no blank margin."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, PathPatch
    from matplotlib.path import Path as OutlinePath
    from shapely.geometry.polygon import orient

    def patch(polygon, colour):
        polygon = orient(polygon, 1.0)                           # the outline counter-clockwise, the holes clockwise: the holes are cut out
        rings = [np.asarray(ring.coords) for ring in (polygon.exterior, *polygon.interiors)]
        codes = np.concatenate([[OutlinePath.MOVETO] + [OutlinePath.LINETO] * (len(ring) - 2) + [OutlinePath.CLOSEPOLY] for ring in rings])
        return PathPatch(OutlinePath(np.concatenate(rings), codes), facecolor=colour, edgecolor="#4b5563", linewidth=0.4)

    fig, ax = plt.subplots(figsize=(9, 9))
    for polygon, colour in ([(faces["steel_as_drawn"], BAND_COLOURS["steel"]), (faces["shaft"], BAND_COLOURS["shaft"])]
                            + [(b, BAND_COLOURS["band"]) for b in faces["bands"]] + [(m, BAND_COLOURS["slot"]) for m in faces["magnets"]]):
        ax.add_patch(patch(polygon, colour))
    r = faces["radii"]["rotor"] * 1.02
    ax.set_xlim(-r, r)
    ax.set_ylim(-r, r)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.legend(handles=[Patch(facecolor=BAND_COLOURS["steel"], edgecolor="#4b5563", label="steel, as drawn (a hole for every pocket)"),
                       Patch(facecolor=BAND_COLOURS["band"], edgecolor="#4b5563", label=f"air band: {len(faces['bands'])} pockets, removed in memory"),
                       Patch(facecolor=BAND_COLOURS["slot"], edgecolor="#4b5563", label=f"magnet slot: {len(faces['magnets'])}, kept"),
                       Patch(facecolor=BAND_COLOURS["shaft"], edgecolor="#4b5563", label="shaft")], loc="lower center", bbox_to_anchor=(0.5, -0.06), ncol=2, frameon=False, fontsize=9)
    ax.set_title(wrap_title(title), fontsize=10.5, color="#1d1d1f", loc="left")
    fig.savefig(target, dpi=150, bbox_inches="tight", pad_inches=0.15)
    plt.close(fig)


def bands_title(faces: dict, earlier: dict, rotor: dict, earlier_rotor: dict) -> str:
    """The title of bands.png: the drawing read, what was removed from it, and the steel before and after against the earlier export's. Each line fits the figure."""
    added = faces["rotor_steel"].area - faces["steel_before_removal_mm2"]
    return (f"{rotor['file']} (the CAD's export of 2026-01-22), read by Python: the rotor's faces\n"
            f"{faces['bands_removed']} air-band pockets found and removed in memory (orange): {added:.2f} mm2 of steel added\n"
            f"{len(faces['magnets'])} magnet slots stay (blue); steel {faces['steel_before_removal_mm2']:,.2f} -> {faces['rotor_steel'].area:,.2f} mm2, "
            f"{earlier_rotor['file']} has {earlier['rotor_steel'].area:,.2f}")


def agrees_with_earlier(faces: dict, earlier: dict, earlier_rotor: dict) -> tuple[bool, float, float]:
    """(whether the band drawing's rotor steel, after the removal, is the earlier export's within the flattening limit; the difference in mm2; the limit in mm2)."""
    difference = faces["rotor_steel"].area - earlier["rotor_steel"].area
    limit = flattening_area_limit(earlier_rotor)
    return abs(difference) <= limit, difference, limit


def path2_caption(info: dict, faces: dict, rotor: dict, stator: dict, earlier: dict, earlier_rotor: dict) -> str:
    """The caption of path 2's half of figure.png (every line at most 80 characters: it stands over half of the picture): which drawings it read, how many air bands were found
    and removed, whether the result agrees with the CAD's earlier export (which never had the band), and the mesh."""
    agrees, _, _ = agrees_with_earlier(faces, earlier, earlier_rotor)
    return (f"Path 2: {rotor['file']} and {stator['file']}, the whole cross-section\n"
            f"ezdxf -> shapely -> Gmsh; {faces['bands_removed']} air bands found and removed in memory\n"
            f"the result {'agrees with' if agrees else 'DIFFERS from'} {earlier_rotor['file']} ({earlier['bands_removed']} bands found there)\n"
            f"{info['nodes']:,} nodes, {info['triangles']:,} triangles, minimum SICN {info['quality']['min_sicn']:.2f}")


def path2_log_line(info: dict, faces: dict, earlier: dict, earlier_rotor: dict, rotor: dict, stator: dict, seconds: float) -> str:
    """The line run.log gets for path 2: the drawing read, the bands found and removed, what the removal did to the steel, whether the earlier export agrees, the mesh."""
    agrees, difference, limit = agrees_with_earlier(faces, earlier, earlier_rotor)
    return (f"B2 path 2: {rotor['file']} + {stator['file']}: {faces['bands_removed']} air band(s) found and removed in memory ({faces['band_lines_dropped']} edges dropped, "
            f"{faces['rotor_steel'].area - faces['steel_before_removal_mm2']:+.2f} mm2 of steel added); the result {'agrees' if agrees else 'DIFFERS'} with "
            f"{earlier_rotor['file']} ({earlier['bands_removed']} band(s) found there; steel {difference:+.4f} mm2 apart, limit {limit:.4f}); {info['nodes']:,} nodes, "
            f"{info['triangles']:,} triangles, {len(info['groups'])} physical groups, min SICN {info['quality']['min_sicn']:.3f}; {seconds:.1f} s")


def path2_files(rotor: dict, earlier_rotor: dict, stator: dict, hashes: dict) -> dict:
    """The drawings path 2 used, for results.json: the paths, the names, and the SHA-256 of each (taken before B2 read them)."""
    return {"rotor": rotor["path"], "rotor_earlier_export": earlier_rotor["path"], "stator": stator["path"],
            "names": {"rotor": rotor["file"], "rotor_earlier_export": earlier_rotor["file"], "stator": stator["file"]}, "sha256": dict(hashes)}


def path1_checks(info: dict, report: dict, table: list[dict], regions, p: mg.Params, chords: list[int]) -> list[dict]:
    """Path 1's checks: every group is there, the quality is stated and has no inverted element, each material's area against motorgeo's, and the
    explanation of the ones that miss (the chord error of the curved boundaries, which falls with the mesh size). `chords` are the mesh edges counted on each
    fillet of the magnet slots (fillet_chords on model_fillets): the explanation is held to what the mesh has, not to a number assumed here."""
    out = [check("path 1: every physical group is present", not report["missing"] and not report["unexpected"] and report["groups_present"] == report["groups_expected"],
                 f"all {report['groups_expected']} (35 regions + 5 boundary curves) and nothing else",
                 f"{report['groups_present']} present; missing {report['missing']}; unexpected {report['unexpected']}", "motorgeo.build: the regions' names; to_gmsh's boundary groups")]
    q = info["quality"]
    w = q["worst"]
    out.append(check("path 1: element quality (minimum SICN) and the worst element's place", q["min_sicn"] > MIN_SICN_OK and q["inverted"] == 0,
                     f"stated; no inverted element and minimum SICN above {MIN_SICN_OK:g} ({OWN_LIMIT}); the brief asks only that the quality is stated",
                     f"minimum SICN {q['min_sicn']:.3f} (mean {q['mean_sicn']:.3f}) in group {w['group']}, element {w['element']}, at ({w['x_mm']:.3f}, {w['y_mm']:.3f}) mm "
                     f"= radius {w['r_mm']:.3f} mm, angle {w['angle_deg']:.2f} deg; {q['below_0_3']} of {info['triangles']:,} triangles below 0.3",
                     "Gmsh minSICN (1 = equilateral, 0 = flat, < 0 = inverted)"))
    geometry = report["geometry_area_mm2"]
    for material in ("M-19 Steel", "N45", "Air", "Coil"):
        out.append(compare(f"path 1: {material} area, mesh against motorgeo", report["mesh_area_mm2"][material], geometry[material], AREA_TOLERANCE_PERCENT,
                           "motorgeo.visible: the shapely polygons of the regions (the brief: within 0.1 %)", "mm2"))
    # the explanations (extra lines, not a way round the numbers above)
    finest = table[-1]
    d = finest["difference_percent"]
    out.append(check("path 1: the area differences are the chord error of curved boundaries, and fall with the mesh size",
                     abs(d["Air"]) < AREA_TOLERANCE_PERCENT and abs(d["Coil"]) < AREA_TOLERANCE_PERCENT and abs(d["M-19 Steel"]) < AREA_TOLERANCE_PERCENT,
                     f"steel, air and coil within {AREA_TOLERANCE_PERCENT:g} % at the smallest size_max of the table ({finest['size_max_mm']:g} mm)",
                     "; ".join(f"size_max {r['size_max_mm']:g}: " + ", ".join(f"{m} {v:+.3f} %" for m, v in r["difference_percent"].items()) + f" ({r['nodes']:,} nodes)" for r in table),
                     "the same function at smaller mesh sizes"))
    out.append(chord_check("path 1: the N45 area is short by the chord error of the fillets (the mesh edges on each counted on the mesh)", model_fillets(regions, p), chords,
                           geometry["N45"] - report["mesh_area_mm2"]["N45"], geometry["N45"],
                           "chord_deficit(), checked against shapely in tests/test_bench_b2.py; the radius is p.slot_fillet; at every size_max of the table the shortfall is "
                           + ", ".join(f"{r['difference_percent']['N45']:+.3f} %" for r in table)))
    return out


def path2_checks(faces: dict, meshed: dict, info: dict, offsets: dict, p: mg.Params, regions, chords: list[int], rotor: dict, stator: dict, earlier: dict,
                 earlier_rotor: dict) -> list[dict]:
    """Path 2's checks: the groups, the quality, the air band (found and removed: motorgeo's band, counted from motorgeo; none left), every entity of the drawings,
    the comparison with the earlier export, the areas of the drawing (`faces`, arcs flattened within 0.25 micrometre) against motorgeo's (4 x its sector), the mesh Gmsh
    wrote against the same (and against the faces it was made from, `meshed`, with coarser chords), the magnets' shortfall as the chord error of their fillets (`chords`:
    the mesh edges counted on each fillet of the rotor drawing), the pattern.
    `rotor` is SynRM rotor.DXF and `stator` Stator_refined.DXF as read_dxf gave them (`faces` and `meshed` are made from them); `earlier_rotor` is Rotor.DXF, the CAD's
    earlier export of the same rotor, which has no band, and `earlier` its faces: the independent check of the removal."""
    out = []
    wanted = cad_group_names(meshed)
    named = {g["name"]: g for g in info["groups"]}
    present = [n for n in wanted if n in named and named[n]["elements"] > 0]
    out.append(check("path 2: every physical group is present", len(present) == len(wanted) and set(named) == set(wanted),
                     f"{len(wanted)}: the rotor steel, {len(faces['magnets'])} magnet slots, the stator steel, {len(faces['slots'])} slots, the gap, and the two A = 0 circles",
                     f"{len(present)} present; missing {sorted(set(wanted) - set(present))}; unexpected {sorted(set(named) - set(wanted))}", "build_faces(): one group per face"))
    q = info["quality"]
    w = q["worst"]
    out.append(check("path 2: element quality (minimum SICN) and the worst element's place", q["min_sicn"] > MIN_SICN_OK and q["inverted"] == 0,
                     f"stated; no inverted element and minimum SICN above {MIN_SICN_OK:g} ({OWN_LIMIT})",
                     f"minimum SICN {q['min_sicn']:.3f} (mean {q['mean_sicn']:.3f}) in group {w['group']}, element {w['element']}, at ({w['x_mm']:.3f}, {w['y_mm']:.3f}) mm "
                     f"= radius {w['r_mm']:.3f} mm, angle {w['angle_deg']:.2f} deg; {q['below_0_3']} of {info['triangles']:,} triangles below 0.3", "Gmsh minSICN"))
    out.append(band_removal_check(faces, p, offsets["rotor"]["turn_deg_mod_45"]))
    out.append(band_check(faces, p))
    out.append(entity_check(faces, rotor, stator))
    out.append(earlier_export_check(faces, earlier, rotor, earlier_rotor))
    visible = mg.visible(regions)
    sector = lambda names: p.n_sectors * sum(visible[n].area for n in names)
    pairs = [("rotor steel", faces["rotor_steel"].area, sector(["rotor_steel"])),
             ("magnet slots (N45)", sum(m.area for m in faces["magnets"]), sector([n for n in visible if n.startswith("magnet_")])),
             ("stator steel", faces["stator_steel"].area, sector(["stator_steel"]))]
    drawn_in = {"rotor steel": f"{faces['files']['rotor']}, its air band removed in memory", "magnet slots (N45)": f"{faces['files']['rotor']}, its air band removed in memory",
                "stator steel": faces["files"]["stator"]}
    for name, area, reference in pairs:
        out.append(compare(f"path 2: {name} area, the CAD drawing against motorgeo's no-band areas", area, reference, CAD_AREA_TOLERANCE_PERCENT,
                           f"motorgeo's sector x {p.n_sectors} (the brief: within 0.5 %, or explained); the drawing is {drawn_in[name]}", "mm2"))
    slots_cad = sum(s.area for s in faces["slots"])
    slots_geo = sector([n for n in visible if n.startswith(("opening_", "coil_"))])
    out.append(compare("path 2: slot area (opening and winding), the CAD drawing against motorgeo's", slots_cad, slots_geo, CAD_AREA_TOLERANCE_PERCENT,
                       f"motorgeo's openings and coil layers x {p.n_sectors}; the drawing is {faces['files']['stator']}", "mm2"))
    area = {g["name"]: g["area_mm2"] for g in info["groups"] if g["dim"] == 2}
    off = max(abs(area["rotor_steel"] / meshed["rotor_steel"].area - 1), abs(area["stator_steel"] / meshed["stator_steel"].area - 1), abs(area["gap"] / meshed["gap"].area - 1),
              *[abs(area[f"magnet_{i + 1}"] / m.area - 1) for i, m in enumerate(meshed["magnets"])], *[abs(area[f"slot_{i + 1}"] / s.area - 1) for i, s in enumerate(meshed["slots"])])
    out.append(check("path 2: the meshed areas are the faces' areas (nothing lost between shapely and Gmsh)", off < 1e-6,
                     f"every group's mesh area equals its face's area to a relative 1e-6 (the faces' edges are the mesh's edges) ({OWN_LIMIT})", f"largest relative difference {off:.2e}",
                     "shapely polygon areas"))
    # the mesh Gmsh wrote, held against motorgeo's areas at the brief's 0.5 %: the same four numbers as the drawing's faces above, from the groups of path2.msh
    mesh_magnets = sum(area[f"magnet_{i + 1}"] for i in range(len(meshed["magnets"])))
    mesh_slots = sum(area[f"slot_{i + 1}"] for i in range(len(meshed["slots"])))
    for name, mesh_area, reference, note in (("rotor steel", area["rotor_steel"], pairs[0][2], ""),
                                             ("magnet slots (N45)", mesh_magnets, pairs[1][2], "; the magnets are short by the chords of their fillets (the next check)"),
                                             ("stator steel", area["stator_steel"], pairs[2][2], ""), ("slot area (opening and winding)", mesh_slots, slots_geo, "")):
        label = name if name.startswith("slot") else f"{name} area"
        out.append(compare(f"path 2: {label}, the mesh Gmsh wrote against motorgeo's", mesh_area, reference, CAD_AREA_TOLERANCE_PERCENT,
                           f"motorgeo's sector x {p.n_sectors}; the mesh is made from the drawing flattened within {MESH_SAGITTA_MM * 1e3:g} micrometres with no segment under "
                           f"{MIN_SEGMENT_MM:g} mm (the brief: within 0.5 %){note}", "mm2"))
    out.append(chord_check("path 2: the magnet area of the mesh is short by the chord error of the fillets (the mesh edges on each counted on the mesh)",
                           drawing_fillets(rotor, p), chords, pairs[1][2] - mesh_magnets, pairs[1][2],
                           "chord_deficit(), the same as path 1's (the shortfall is against motorgeo's N45 area, as there), on the rotor drawing's fillet arcs and the nodes of path2.msh"))
    r, s = offsets["rotor"], offsets["stator"]
    out.append(check("path 2: the drawing's magnet slots are motorgeo's, turned", r["max_distance_mm"] < CAD_POSITION_TOLERANCE_MM,
                     f"the sixteen slot centres coincide with motorgeo's after one rotation, within {CAD_POSITION_TOLERANCE_MM:g} mm ({OWN_LIMIT}); the CAD gives 0.01 mm",
                     f"turned by {r['turn_deg_mod_45']:.3f} degrees (modulo the 45-degree pole pitch); largest distance left {r['max_distance_mm']:.4f} mm",
                     "cad_offset(): least-squares rotation of matched points"))
    bore = faces["radii"]["bore"]
    stator_mm = math.radians(s["max_distance_deg"]) * bore
    out.append(check("path 2: the drawing's stator slots are motorgeo's, turned", stator_mm < CAD_POSITION_TOLERANCE_MM,
                     f"the {s['points']} slots' angles coincide with motorgeo's after one rotation, within {CAD_POSITION_TOLERANCE_MM:g} mm of arc at the bore radius "
                     f"({bore:g} mm) ({OWN_LIMIT})",
                     f"turned by {s['turn_deg_mod_10']:.3f} degrees (modulo the 10-degree slot pitch); largest angular distance left {s['max_distance_deg']:.4f} degrees "
                     f"= {stator_mm:.4f} mm at the bore", "cad_offset(): least-squares rotation of matched points (the slots by their angle)"))
    return out
