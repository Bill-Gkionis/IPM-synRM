"""B4: Gmsh -> CalculiX 2.23 (ccx, plane stress, second-order triangles): the band-free rotor lamination spinning at 10,000 and 12,000 rpm.

What it shows, and why
    The rotor's speed limit is expected to be its bridges (04 Structural, 4.3): thin strips of steel that hold each pole's outer steel and its magnets against the
    centrifugal pull. Task 4 proved CalculiX on a ring, where Lame gives the answer. B4 sends the MOTOR through it: the 90-degree sector of motorgeo.Params(barrier="none"),
    meshed by Gmsh, solved as plane stress, and the peak von Mises stress looked for in the steel. It is sound when that peak sits in a bridge or a rib (and the case says
    where), when two speeds give stresses in the ratio (12/10)^2 = 1.44, and when the bridges' nominal stress agrees with the note's hand estimate by the estimate's own method.

Method (the files and the protocol)
    1. The mesh. motorgeo.build gives the rotor's regions (the steel ring and the four magnet slots of the sector, band-free). Gmsh (the Python API, in this process, inside
       `with intact_path():`) fragments them into a steel body with four holes and four slot bodies that share every wall curve, and meshes them in second-order
       triangles (6 nodes, `Mesh.ElementOrder = 2`). The size is a field: h_min at the slot walls (the bridges are 1.56 mm wide, the corners are fillets of 0.3 mm) growing to
       2 mm in 8 mm. Python takes the nodes and the elements out of Gmsh as arrays (rotor.msh is also saved, with the groups, for a look in Gmsh); there is no converter in
       between: the deck is written from those arrays.
    2. The model: the 90-degree sector with its two symmetry edges (a node on the x axis may not move in y, one on the y axis may not move in x; the sector is mirror
       symmetric about both because they are q axes), no shaft (a free bore, as the note's disc), plane stress (CPS6) with the thickness of the stack, 20 mm: stress under a body
       force does not depend on the thickness, the forces in N do. The 20 mm is motorgeo's `stack`, the depth of the 2D model (FEMM's): one row of magnets, 19.14 mm long (the 2D
       work runs per row, and the stack's length is set after the optimization: the project plan's Orchestration note, section 4). The starting CAD model has 3 rows, a 60 mm stack,
       which is B6's: the forces in N here would be 3 times larger at 60 mm, and the stresses and the pressures on the walls would not change. Units: mm, N, tonne, s, so stress is
       in MPa and density in tonne/mm3. M-19: E 200 GPa, nu 0.3, 7,650 kg/m3 (the note, 6.4); N45: 160 GPa, 0.24, 7,600 kg/m3 (6.1 and its bridge formula). The load is *DLOAD
       CENTRIF about the z axis.
       Why the 90-degree sector and not one pole (the brief allows either): its two symmetry edges lie on the axes, so the supports are the plain "y fixed on this edge, x fixed on
       that one" (a pole's 45-degree edge would need a rotated constraint); it is the domain of B1 to B3's field models; and the two poles in it are mirror images that Gmsh meshes
       differently, so the four slots' peaks agreeing (they do, within 1.3 %) is a check on the mesh at no cost.
    3. The magnets are mass on the slot walls, and the model does not say how, so B4 runs both ways (and a third time, as a diagnostic) and shows what the choice changes:
         pressure  The slot is a hole. The magnet's whole centrifugal force F = m w^2 r (m the real magnet, 9.08 x 1.57 mm x 20 mm of N45, r the distance of its centre from the
                   axis, along that radius) is split along the slot's own axes, u along the slot and n across it. The part across (F.n) is a uniform pressure on the long wall the
                   magnet leans on, the part along (F.u) a uniform pressure on its outer end wall: what a magnet resting frictionless in its slot does. The pressures are written as
                   nodal forces (the consistent 1/6, 4/6, 1/6 of a quadratic edge), *CLOAD. This is the conservative model and the case's verdict.
         bonded    The slot is filled by a second body (E 160 GPa, nu 0.24) whose density is scaled so that it has the real magnet's mass, not the slot's 14.7 % more volume of
                   N45, glued to the steel all along the wall. The glue is `*EQUATION`: the body has its own copy of every wall node, tied to the steel's node in both DOFs, so that
                   CalculiX's nodal stress (an average over the elements that share a node) is the steel's alone; a shared node would mix the two materials exactly where the peak is.
         soft      A diagnostic, reported and never gating, at 10,000 rpm: the bonded deck with the slot bodies' Young's modulus cut to 20 MPa (0.01 % of the steel's, under the 0.1 % that
                   counts as a negligible share of the hoop stiffness) and their mass kept. It has the bonded model's load path (each body hangs on the walls through the ties) and
                   the pressure model's lost stiffness (the slots carry no hoop stress). So it says which of the two raises the pressure model's bore stress above Lame's disc: if the
                   soft model's bore is near the pressure model's, the slots' lost stiffness does; if it is near the bonded model's, the load's path does.
       Same steel mesh and same total load in all three (checked to six digits); only the way the load reaches the steel, and the stiffness of what fills the slots, differ.
    4. ccx_static.exe -i job through bench.run_program (the solver slot, below normal, 4 threads, a timeout). It writes job.frd (displacements, nodal stresses) and job.dat (the
       reaction forces of the two symmetry edges, `*NODE PRINT, TOTALS=ONLY`). Python reads them back: checks.ccx_check.read_frd_stress for the stress, the .dat for the reactions.
    5. The numbers. The peak von Mises over the steel's nodes (and its place: the nearest slot, the feature it sits in), the same at the other speed, the bore's hoop stress against
       Lame's disc, the sum of the reactions against the load worked out from the mesh (equilibrium), and the force that the three strips of a pole carry, by integrating the stress
       over three cuts. A cut is the narrowest line across each strip: the rib between a pole's two slots, and the line from each slot's outermost point to the rim.
    6. The hand estimate of the note (6.1): a pole's outer steel (the cap) with its magnets hangs on the three strips; pull = w^2 (rho m A r), nominal stress = pull / (the three
       widths). The note made it for the rotor WITH its air band (162 mm2, strips 1.5 + 1.4 + 1.4 mm, 15 MPa at 10,000 rpm); estimate() does the same for either rotor, and
       reproduces the note's 15 MPa on the rotor with the band before it is applied to the band-free one. Without the band the strips are 3.30 + 1.56 + 1.56 mm and the cap
       is the 94 mm2 of steel in the mouth of the V: 6.6 MPa. B4 holds the finite-element force through the strips against that.
    7. The picture: cgx's own hardcopy (`read`, `ds`, `plot fv`, `hcpy tga`: the von Mises of the verdict run on the sector). cgx writes TGA frames and then calls ImageMagick's
       mogrify and composite, which this PC does not have, so Python puts the frames together; cgx multiplies its colours by a constant 0.298 in the model frame (lighting), which
       Python divides out (checked: every colour of the frame is a colour of the legend times 0.298). cgx does not come back from `hcpy` to its script, so run_solver ends it
       at the timeout once the frames are there. If the frames are not made, the picture is plotted from the .frd in Python.

Run (from the Python to FEMM folder; the solver slot is taken for each ccx run, a minute or two in all)
    uv run python bench/run_bench.py --only B4
Writes out\\bench\\B4\\: results.json, figure.png (cgx's von Mises of the verdict run), stress_full.png (the whole rotor, from Python), rotor.msh, run.log, runs.json, and one
folder per solve (pressure_10000_h0.1\\, bonded_12000_h0.1\\, soft_10000_h0.1\\, pressure_10000_h0.2\\, ...) with job.inp, job.frd, job.dat and ccx's own printout.

Numbers to expect (this PC, 2026-11-06, CalculiX 2.23, Gmsh 4.15.2; the whole case takes 135 s: the seven ccx solves 104 s (6, 6, 20, 20 and 20 s on the 0.1 mm mesh, 2 s on the 0.2 mm
mesh and 29 s on the 0.05 mm one), cgx's picture 15 s, the rest Gmsh's three meshes and the readers and analysis of seven solves)
    mesh        h_min 0.1 mm: 70,534 nodes (mid-side nodes included; 45,271 of them in the steel's 22,130 triangles, which is all the pressure model has: its slots are voids; the
                bonded model adds 13,069 triangles in the four slots), minimum SICN 0.669, the areas within 4e-4 % of motorgeo's; the steel's nodes at h_min 0.2 mm: 17,424;
                at 0.05 mm: 132,924
    loads       one magnet (2.167 g over the 20 mm stack) at 10,000 rpm: F = 94.73 N, of which 73.98 N across the slot (0.426 MPa on the 8.68 mm long wall) and 59.17 N along it (2.53 MPa
                on the 1.17 mm end wall); the bonded body's density is 6,627 kg/m3 (7,600 x 14.26 / 16.35); the totals of the three models agree to 6 digits (7,120.635 N for the sector)
    pressure    the peak is 33.45 MPa at 10,000 rpm and 48.17 at 12,000 (ratio 1.4400), at (42.81, 7.36) mm: the outer bridge of magnet_0_dn, 1.56 mm from the rim, in the corner
                fillet at the slot's outer end on the loaded side; the four slots' bridges 33.02-33.45 MPa; the centre rib's corners 28.4-28.5 MPa; by mesh 31.82 (h_min 0.2 mm), 33.45
                (0.1), 33.85 (0.05). The bore's mean hoop stress (13 nodes): 14.91 MPa (a plain disc's, Lame's: 14.10, so +5.75 %)
    bonded      the bore's mean hoop stress is 14.08 MPa (Lame's 14.10, -0.11 %); the largest von Mises of the whole steel is also at the bore, 13.99 MPa (there the nodal radial stress
                is 0.02-0.29 MPa, so the von Mises is a little under the hoop stress, 14.02-14.15); the bridges only 4.25 MPa (7.9 times less than the pressure model's), the slot bodies
                3.77 MPa at any of their nodes (N45's tensile strength is 78 MPa)
    soft        (the diagnostic: the slot bodies' E 20 MPa, their mass kept) the bore's mean hoop stress is 14.87 MPa (+5.46 %): 95 % of the way from the bonded model's to the pressure
                model's, so the open slots' lost hoop stiffness raises the pressure model's bore stress and the way the load reaches the steel barely does (up to 5 % of the gap, about
                0.04 of 0.83 MPa, can still be the load's path); the bridges' zones 30.3 MPa: 89 % of the way, so of the difference between the pressure model's bridges and the bonded
                model's (33.45 - 4.25 MPa; their ratio is the factor 7.9) about nine tenths is the stiffness of what fills the slots and one tenth the load's path; the same total load
                (7,120.635 N) and equilibrium to 6e-4 %. E = 2 and 200 MPa (scratch runs of 2026-11-06, not part of the bench) give a bore of 14.87 and 14.86 MPa and
                bridges of 30.3 and 30.0: the soft limit is reached
    strips      the rib and the two bridges of a pole carry 787.8 N and 785.4 N (poles 0 and 1) over the stack; the hand estimate says 847.5 N (-7.1 %, -7.3 %: all of the magnets'
                pull is put on the cap), and 788.3 N with the end walls' part left off (-0.07 %, -0.37 %). The rib carries 78 % of it and each bridge 11 %. Nominal stress F / 6.41 mm:
                6.13 MPa (hand 6.61, on the same basis: the pull along the pole's axis; the note's own formula adds the forces' magnitudes, 849.6 N and 6.62 MPa), so the peak/nominal
                factor is 5.5 (hand 5.1). The note's 15 MPa is the rotor with the band (cap 162 mm2, strips 4.3 mm); its own function gives 15.11
                MPa on that rotor
    equilibrium the reactions balance the load to 6e-4 %; cgx's own maximum is 33.454 MPa, ours 33.454

Self-check (the answers follow).
  1. Why can one solve at 10,000 rpm give the stress at 12,000?
  2. The slot is filled with a body in the bonded model and the peak falls from 34 to 4 MPa. Which model is right?
  3. Why does the case write the magnet's force as pressures on two walls, and not all of it on the long wall?
  4. A density of 6.627128124470999e-09 went into a deck and the stress came out 5.8e8 MPa. Why, and where does the case prevent it?
  5. The pressure model's bore stress is 5.75 % above Lame's disc and the bonded model's is not. What makes the difference, and how does the case know?
Answers.
  1. Linear elasticity, and a load that is rho w^2 r at every point: the displacement and every stress are proportional to w^2, so the stress at 12,000 rpm is (12/10)^2 = 1.44
     times that at 10,000. It holds while nothing else loads the rotor (no press fit, no heat, no magnet that lifts off its slot): the note says so (8, question 3).
  2. Neither is the truth; they bound it. In the bonded model a stiff body (E 160 GPa) fills the hole and takes part of the load: the stress concentration of a hole is not there.
     A real magnet is brittle, glued with an adhesive 50 times softer than steel (3 GPa, the note's table), and its tensile strength is 78 MPa: it cannot be counted on to
     carry the bridges, so the pressure model, which treats it as a mass that leans on the wall, is the one a design has to meet. B4 reports both.
  3. A magnet whose centre is 8.7 degrees off the pole axis, in a slot at 60 degrees to it, is pushed along a line 51 degrees off the slot's own axis: 78 % of the force across the
     slot and 62 % along it (the sine and the cosine of 51 degrees). Across is carried by the long wall; along must go somewhere: the outer end wall. With all of it on the long
     wall the end wall would be unloaded and the long wall would carry 28 % more than the force across the slot, in a wrong direction.
  4. ccx reads 20 characters of a number: the 21-character density was cut to 6.627e-0, a billion times too heavy, and nothing warns. deck_number() writes every number of a deck
     with 12 significant digits (at most 19 characters) and a test pins it.
  5. Barely the way the load reaches the steel: the soft model has the bonded model's load path (the body hangs on the walls through ties) and the pressure model's lost stiffness
     (E = 20 MPa carries no hoop stress), and its bore is 95 % of the way to the pressure model's; the last 5 % of the gap (about 0.04 of 0.83 MPa) can still be the load's path.
     It is the slots, nearly all of it: the ring that carries the hoop stress is broken at every slot, so the rest of the steel carries more of it. Only a third model with one
     of the two differences could say so; a guess from the pressure and the bonded models could not.
"""
from __future__ import annotations

import dataclasses
import functools
import math
import re
import shutil
import time
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Point
from shapely.ops import nearest_points, unary_union

import motorgeo as mg
import paths
from checks import ccx_check, intact_path

from . import OWN_LIMIT, Log, check, compare, legend_and_model_figure, params, run_program, skipped, wrap_title

CASE = "B4"
TITLE = "Gmsh -> CalculiX 2.23 (plane stress, CPS6): the band-free rotor lamination under centrifugal load, 10,000 and 12,000 rpm"

SPEEDS_RPM = (10_000, 12_000)
VARIANTS = ("pressure", "bonded")        # the two models of the verdict and of its other bound, each at both speeds
SOFT = "soft"                            # the diagnostic (reported, never gating): the bonded deck with the slot bodies' stiffness made negligible and their mass kept; 10,000 rpm only
ALL_VARIANTS = VARIANTS + (SOFT,)
TIED_VARIANTS = ("bonded", SOFT)         # the variants that fill the slots with a body tied to the steel (a copy of every wall node, two *EQUATION lines each)
SOFT_FILLER_YOUNGS_MPA = 20.0            # the soft body's E: 0.01 % of the steel's 200 GPa, ten times under the 0.1 % that still counts as a negligible share of the hoop stiffness
VERDICT_VARIANT = "pressure"             # the conservative model: the magnet is a mass that leans on the wall, the slot is a hole
SCALE_TOLERANCE_PERCENT = 1.0            # the brief: the two speeds' stresses scale by (12/10)^2 = 1.44 within 1 %
EXPECTED_RATIO = (SPEEDS_RPM[1] / SPEEDS_RPM[0]) ** 2
NOTE_HAND_NOMINAL_MPA = 15.0             # 04 Structural 6.1 and 6.4: the bridges' nominal stress at 10,000 rpm, for the rotor WITH its air band
NOTE_TOLERANCE_PERCENT = 3.5             # the note prints "15 MPa" (two digits): the value is anywhere in 14.5-15.5, which is 3.3 % either way; 3.5 % is that rounded up

# The limits the brief does not give: the bench's own, each said where it is used (OWN_LIMIT)
EQUILIBRIUM_TOLERANCE_PERCENT = 0.05     # the reactions against the load worked out from the mesh (agreement is 1e-4 % or better: exact quadrature on both sides)
SAME_LOAD_TOLERANCE_PERCENT = 0.01       # the two magnet models carry the same total load
HAND_TOLERANCE_PERCENT = 10.0            # the force through the strips against the hand estimate (which puts all of the magnets' pull on the cap: 7 % too much, see below)
REFINED_TOLERANCE_PERCENT = 2.0          # ... and against the hand estimate with the magnets' end-wall part left off the cap (integration of a smoothed field over three cuts)
MESH_TOLERANCE_PERCENT = 2.0             # the peak on the main mesh against the finest: chosen after the exploratory runs (5.1 % from 0.2 to 0.1 mm, 1.2 % from 0.1 to 0.05)
SYMMETRY_TOLERANCE_PERCENT = 3.0         # the four slots' bridge peaks (mirror images of each other, meshed differently)
LAME_TOLERANCE_PERCENT = 1.0             # a gate on the units, the speed and the load: the bonded model's bore hoop stress against Lame's plain disc (it is not that disc: by mass alone about 0.6 % from it);
#                                          chosen after the first run (-0.11 %), not before; the pressure model, whose slots are voids, is reported only
SOFT_DECIDES = 0.25                      # the soft model's bore decides between the two hypotheses when it lies within this share of the gap between the bonded and the pressure model's, from one of them
CGX_TOLERANCE_PERCENT = 0.01             # our von Mises peak against the maximum cgx prints for the same .frd (the same numbers, two programs)
MIN_SICN_OK = 0.1                        # no degenerate element (B2's floor)
AREA_TOLERANCE_PERCENT = 0.05            # the mesh's steel and magnets against motorgeo's areas
STRIP_ZONE_MM = 2.0                      # "in a bridge or rib": within this of the cut across a strip (the bridge is 1.56 mm wide, the fillet 0.3)

H_MIN_MAIN, H_MIN_COARSE, H_MIN_FINE = 0.1, 0.2, 0.05     # mm: the element size at the slot walls (main mesh, and the two meshes of the convergence)
H_MAX, DIST_MIN, DIST_MAX = 2.0, 0.5, 8.0                 # mm: the size far from the slots, and the distances over which one grows into the other
CCX_THREADS = 4
CCX_TIMEOUT_S = 900
CGX_TIMEOUT_S = 15
CGX_LIGHT = 0.298                         # cgx's lighting multiplies every colour of the model frame by this (measured against its own legend)
SECTION_SAMPLES = 100                     # points along a cut where the stress is read


# ============================================================================================== the inputs
@dataclasses.dataclass(frozen=True)
class Material:
    """The structural constants of a material, in the units of the note (MPa, kg/m3); `density_tonne_mm3` is the same density in the units of the deck."""
    name: str
    youngs_mpa: float
    poisson: float
    density_kg_m3: float
    source: str

    @property
    def density_tonne_mm3(self) -> float:
        return self.density_kg_m3 / 1e12           # 1 kg/m3 = 1e-12 tonne/mm3


M19 = Material("M-19", 200e3, 0.3, 7650.0, "04 Structural 6.4: E 200 GPa, nu 0.3 (a handbook value), 7,650 kg/m3 (the maker's sheet)")
N45 = Material("N45", 160e3, 0.24, 7600.0, "04 Structural 6.1: E 160 GPa, nu 0.24, 7,500-7,600 kg/m3; 7,600 is the value of the note's bridge formula")


def omega(rpm: float) -> float:
    """Angular speed [rad/s] of a speed in rpm."""
    return 2 * math.pi * rpm / 60


def deck_number(value: float) -> str:
    """A number as it is written into a deck: 12 significant digits, at most 19 characters. CalculiX reads 20 characters of a field and cuts the rest: a density of
    6.627128124470999e-09 (21 characters) became 6.627e-0 and gave 5.8e8 MPa with no warning. Every number a deck holds goes through here."""
    return f"{value:.12g}"


# ============================================================================================== the hand estimate (04 Structural 6.1)
@functools.lru_cache(maxsize=8)
def _visible(p: mg.Params) -> dict:
    return mg.visible(mg.build(p, gap="sliding"))


def strip_geometry(p: mg.Params, pole: int = 0) -> dict:
    """The three strips that hold the outer steel of a pole (its "cap"), and the cap, for the rotor of `p`.

    The note's method (materials/structural.py F): the two holes of the pole (the two air bands when the rotor has them, otherwise the two magnet slots) are separated by the
    centre rib, and each hole is separated from the rim by a bridge. Cut the steel across the three strips and the cap comes free.
      rib          the shortest line between the two holes (the narrowest steel between them), its length the rib's width
      bridge_low   from the farthest point of the lower hole (the hole at the smaller angle) along the radius to the rim, its length the bridge's width; bridge_high likewise
    Returns {"axis_deg", "holes" (the two region names), "cap" (a shapely polygon), "cuts" ({"rib", "bridge_low", "bridge_high"}: the cut lines, in steel from edge to edge),
    "widths" (their lengths, mm), "pole"}."""
    vis = _visible(p)
    axis = 180 / p.pole_pairs * (pole + 0.5)
    names = (f"barrier_{pole}_dn", f"barrier_{pole}_up") if p.barrier == "air" else (f"magnet_{pole}_dn", f"magnet_{pole}_up")
    low, high = vis[names[0]], vis[names[1]]
    near_low, near_high = nearest_points(low, high)
    along = np.array([near_high.x - near_low.x, near_high.y - near_low.y])
    rib = float(np.hypot(*along))
    along /= rib
    cuts = {"rib": LineString([(near_low.x, near_low.y), (near_high.x, near_high.y)])}
    extended = [LineString([(near_low.x - 0.3 * along[0], near_low.y - 0.3 * along[1]), (near_high.x + 0.3 * along[0], near_high.y + 0.3 * along[1])])]
    widths = {"rib": rib}
    for key, hole in (("bridge_low", low), ("bridge_high", high)):
        points = np.array(hole.exterior.coords)
        radius = np.hypot(points[:, 0], points[:, 1])
        r_max = float(radius.max())
        angles = np.arctan2(points[radius > r_max - 0.01][:, 1], points[radius > r_max - 0.01][:, 0])
        middle = (angles.min() + angles.max()) / 2
        direction = np.array([math.cos(middle), math.sin(middle)])
        cuts[key] = LineString([tuple(r_max * direction), tuple(p.r_ro * direction)])
        extended.append(LineString([tuple((r_max - 0.5) * direction), tuple((p.r_ro + 0.5) * direction)]))
        widths[key] = p.r_ro - r_max
    pieces = vis["rotor_steel"].difference(unary_union([line.buffer(0.005) for line in extended]))
    candidates = [g for g in getattr(pieces, "geoms", [pieces]) if abs(math.degrees(math.atan2(g.centroid.y, g.centroid.x)) - axis) < 3]
    cap = max(candidates, key=lambda g: math.hypot(g.centroid.x, g.centroid.y))
    return {"axis_deg": axis, "holes": names, "cap": cap, "cuts": cuts, "widths": widths, "pole": pole}


def hand_estimate(p: mg.Params, rpm: float, steel: Material = M19, magnet: Material = N45) -> dict:
    """The note's hand estimate of the bridges, per metre of stack: the pull of one pole's cap and its two magnets, and the nominal stress if the three strips shared it evenly.

    pull = w^2 (rho_steel A_cap r_cap + rho_magnet A_magnets r_magnet); nominal = pull / (rib + bridge + bridge). The magnets are the real ones, 2 x 9.08 x 1.57 mm2 per pole,
    at the distance of the slot's centre. Returns the cap's area and the distance of its centroid, the magnets', the three widths, the pull in kN/m and over the stack, the nominal
    stress, the mass of one pole and the magnets' share of the pull."""
    g = strip_geometry(p, 0)
    cap = g["cap"]
    area_cap, r_cap = cap.area * 1e-6, math.hypot(cap.centroid.x, cap.centroid.y) * 1e-3
    slot = _visible(p)["magnet_0_dn"]
    area_magnet, r_magnet = 2 * p.mag_w * p.mag_t * 1e-6, math.hypot(slot.centroid.x, slot.centroid.y) * 1e-3
    moment = steel.density_kg_m3 * area_cap * r_cap + magnet.density_kg_m3 * area_magnet * r_magnet
    widths = g["widths"]
    total = widths["rib"] + widths["bridge_low"] + widths["bridge_high"]
    pull = omega(rpm) ** 2 * moment
    return {"rpm": rpm, "cap_area_mm2": cap.area, "cap_centroid_mm": r_cap * 1e3, "magnet_area_mm2": area_magnet * 1e6, "magnet_centroid_mm": r_magnet * 1e3,
            "strips_mm": {"centre_rib": widths["rib"], "outer_bridge": (widths["bridge_low"] + widths["bridge_high"]) / 2, "total": total},
            "pull_kn_per_m": pull / 1e3, "pull_n_over_stack": pull * p.stack * 1e-3, "nominal_mpa": pull / (total * 1e-3) / 1e6,
            "mass_per_pole_g": (steel.density_kg_m3 * area_cap + magnet.density_kg_m3 * area_magnet) * p.stack * 1e-3 * 1e3,
            "magnet_share_of_pull": magnet.density_kg_m3 * area_magnet * r_magnet / moment}


# ============================================================================================== the mesh
@dataclasses.dataclass
class Mesh:
    """What Gmsh made, as plain data: nodes {id: (x, y)} in mm; steel = the 6-node triangles of the steel body; magnets = {slot name: its triangles}; boundary = {name: 3-node line
    elements (end, end, mid)} for the slot outlines ("slot:<name>") and the two symmetry edges ("edge_0", "edge_90"). Every triangle has a positive Jacobian (corner order
    counter-clockwise, then the mid nodes of the edges 1-2, 2-3, 3-1: the numbering of Gmsh's type 9 and of CalculiX's CPS6)."""
    nodes: dict
    steel: list
    magnets: dict
    boundary: dict
    h_min: float = 0.0
    h_max: float = 0.0
    quality: dict = dataclasses.field(default_factory=dict)
    msh: str = ""


def _signed_area(nodes: dict, tri: tuple) -> float:
    (x1, y1), (x2, y2), (x3, y3) = nodes[tri[0]], nodes[tri[1]], nodes[tri[2]]
    return 0.5 * ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))


def _counter_clockwise(nodes: dict, tri: tuple) -> tuple:
    """The triangle with its corners counter-clockwise (corners 2 and 3 exchanged, and with them the mid nodes of the edges 1-2 and 3-1, if not)."""
    return tri if _signed_area(nodes, tri) > 0 else (tri[0], tri[2], tri[1], tri[5], tri[4], tri[3])


def mesh_rotor(p: mg.Params, regions, h_min: float, h_max: float = H_MAX, dist_min: float = DIST_MIN, dist_max: float = DIST_MAX, msh_path=None) -> Mesh:
    """Mesh the rotor of the sector (the steel and the four magnet slots of `regions`) in second-order triangles. The size is h_min within dist_min of a slot wall and h_max
    beyond dist_max, linear between (a Distance and a Threshold field). `msh_path`, if given, gets rotor.msh (MSH 2.2, with the physical groups) for a look in Gmsh."""
    import gmsh
    rotor = [r for r in regions if r.side == "rotor" and r.material in ("M-19 Steel", "N45")]
    with intact_path():                                    # gmsh.initialize() cuts the PATH of this process (bench/__init__.py, decision 6)
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.model.add("rotor")
            occ = gmsh.model.occ
            surface = {}
            for i, r in enumerate(rotor):
                points: dict = {}

                def point(xy, points=points):
                    if mg.key(xy) not in points:
                        points[mg.key(xy)] = occ.addPoint(xy[0], xy[1], 0)
                    return points[mg.key(xy)]

                curves = []
                for e in r.edges:
                    if isinstance(e, mg.Line):
                        curves.append(occ.addLine(point(e.a), point(e.b)))
                    else:                                  # OCC arcs must be < 180 degrees: pieces of at most 90
                        n = math.ceil(e.sweep() / 90)
                        a0 = e.angles()[0]
                        sweep = e.sweep() * (1 if e.ccw else -1)
                        knots = [e.a] + [mg.lin(e.c, mg.pol(1, a0 + sweep * j / n), e.radius) for j in range(1, n)] + [e.b]
                        centre = occ.addPoint(*e.c, 0)
                        curves += [occ.addCircleArc(point(q0), centre, point(q1)) for q0, q1 in zip(knots[:-1], knots[1:])]
                surface[i] = occ.addPlaneSurface([occ.addCurveLoop(curves)])
            order = list(surface)
            _, owners = occ.fragment([(2, surface[order[0]])], [(2, surface[i]) for i in order[1:]])
            owner = {}
            for i, dimtags in zip(order, owners):          # later regions win: a magnet is painted over the steel
                for _, tag in dimtags:
                    owner[tag] = i
            occ.synchronize()
            groups: dict = {}
            for tag, i in owner.items():
                groups.setdefault(rotor[i].name, []).append(tag)
            slot_curves = {name: sorted({abs(c) for _, c in gmsh.model.getBoundary([(2, t) for t in tags], oriented=False)})
                           for name, tags in groups.items() if name.startswith("magnet")}
            edge_curves = {"edge_0": [], "edge_90": []}
            for _, c in gmsh.model.getEntities(1):
                x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(1, c)
                if max(abs(y0), abs(y1)) < 1e-9:
                    edge_curves["edge_0"].append(c)
                elif max(abs(x0), abs(x1)) < 1e-9:
                    edge_curves["edge_90"].append(c)
            field = gmsh.model.mesh.field
            distance = field.add("Distance")
            field.setNumbers(distance, "CurvesList", sorted({c for curves in slot_curves.values() for c in curves}))
            field.setNumber(distance, "Sampling", 100)
            threshold = field.add("Threshold")
            field.setNumber(threshold, "InField", distance)
            field.setNumber(threshold, "SizeMin", h_min)
            field.setNumber(threshold, "SizeMax", h_max)
            field.setNumber(threshold, "DistMin", dist_min)
            field.setNumber(threshold, "DistMax", dist_max)
            field.setAsBackgroundMesh(threshold)
            for option in ("Mesh.MeshSizeFromPoints", "Mesh.MeshSizeFromCurvature", "Mesh.MeshSizeExtendFromBoundary"):
                gmsh.option.setNumber(option, 0)
            gmsh.option.setNumber("Mesh.ElementOrder", 2)
            gmsh.option.setNumber("Mesh.SecondOrderLinear", 0)       # the mid nodes of curved edges lie on the curve
            gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
            for name, tags in groups.items():
                gmsh.model.addPhysicalGroup(2, tags, name=name)
            for name, curves in {**{f"slot:{n}": c for n, c in slot_curves.items()}, **edge_curves}.items():
                if curves:
                    gmsh.model.addPhysicalGroup(1, curves, name=name)
            gmsh.model.mesh.generate(2)
            tags, coords, _ = gmsh.model.mesh.getNodes()
            xyz = np.asarray(coords).reshape(-1, 3)
            xy = {int(t): (float(x), float(y)) for t, (x, y, _) in zip(tags, xyz)}
            triangles = {}
            for name, surface_tags in groups.items():
                rows = []
                for tag in surface_tags:
                    types, _, connectivity = gmsh.model.mesh.getElements(2, tag)
                    for kind, nodes in zip(types, connectivity):
                        assert kind == 9, f"expected 6-node triangles (type 9), got type {kind}"
                        rows += [_counter_clockwise(xy, tuple(int(n) for n in row)) for row in np.asarray(nodes).reshape(-1, 6)]
                triangles[name] = rows
            lines = {}
            for _, c in gmsh.model.getEntities(1):
                types, _, connectivity = gmsh.model.mesh.getElements(1, c)
                lines[c] = [tuple(int(n) for n in row) for kind, nodes in zip(types, connectivity) for row in np.asarray(nodes).reshape(-1, 3) if kind == 8]
            boundary = {f"slot:{name}": [row for c in curves for row in lines[c]] for name, curves in slot_curves.items()}
            boundary.update({name: [row for c in curves for row in lines[c]] for name, curves in edge_curves.items()})
            element_tags, element_nodes = gmsh.model.mesh.getElementsByType(9)
            sicn = np.asarray(gmsh.model.mesh.getElementQualities(element_tags, "minSICN"))
            corner = np.asarray(element_nodes).reshape(-1, 6)[:, :3]
            worst = int(np.argmin(sicn))
            where = np.mean([xy[int(n)] for n in corner[worst]], axis=0)
            quality = {"triangles": int(len(element_tags)), "min_sicn": float(sicn.min()), "mean_sicn": float(sicn.mean()),
                       "worst": {"x_mm": float(where[0]), "y_mm": float(where[1])}}
            if msh_path is not None:
                gmsh.write(str(msh_path))
        finally:
            gmsh.finalize()
    steel = triangles["rotor_steel"]
    magnets = {name: rows for name, rows in triangles.items() if name.startswith("magnet")}
    used = {n for t in steel for n in t} | {n for rows in magnets.values() for t in rows for n in t}
    used |= {n for rows in boundary.values() for row in rows for n in row}
    return Mesh(nodes={n: xy[n] for n in sorted(used)}, steel=steel, magnets=magnets, boundary=boundary, h_min=h_min, h_max=h_max, quality=quality,
                msh=str(msh_path) if msh_path is not None else "")


# 6-point Gauss rule of degree 4 on a triangle (Dunavant): exact for the area and the first moments of a quadratic (curved) 6-node triangle
_A1, _A2 = 0.445948490915965, 0.091576213509771
_W1, _W2 = 0.223381589678011, 0.109951743655322
_GAUSS = [((_A1, _A1, 1 - 2 * _A1), _W1), ((_A1, 1 - 2 * _A1, _A1), _W1), ((1 - 2 * _A1, _A1, _A1), _W1),
          ((_A2, _A2, 1 - 2 * _A2), _W2), ((_A2, 1 - 2 * _A2, _A2), _W2), ((1 - 2 * _A2, _A2, _A2), _W2)]


def _shape(l1: float, l2: float, l3: float):
    """The six shape functions of a quadratic triangle at the area coordinates (l1, l2, l3), and their derivatives with respect to (l1, l2), l3 = 1 - l1 - l2."""
    n = np.array([l1 * (2 * l1 - 1), l2 * (2 * l2 - 1), l3 * (2 * l3 - 1), 4 * l1 * l2, 4 * l2 * l3, 4 * l3 * l1])
    d1 = np.array([4 * l1 - 1, 0.0, -(4 * l3 - 1), 4 * l2, -4 * l2, 4 * l3 - 4 * l1])
    d2 = np.array([0.0, 4 * l2 - 1, -(4 * l3 - 1), 4 * l1, 4 * l3 - 4 * l2, -4 * l1])
    return n, d1, d2


_RULE = [(_shape(*point), weight) for point, weight in _GAUSS]


def t6_integrals(nodes: dict, tri: tuple) -> tuple[float, float, float]:
    """(area, integral of x dA, integral of y dA) of a 6-node triangle over its curved shape (the mid nodes are where they are), 6 Gauss points."""
    xy = np.array([nodes[n] for n in tri])
    area = mx = my = 0.0
    for (n, d1, d2), weight in _RULE:
        j11, j12, j21, j22 = d1 @ xy[:, 0], d1 @ xy[:, 1], d2 @ xy[:, 0], d2 @ xy[:, 1]
        weighted = 0.5 * weight * (j11 * j22 - j12 * j21)
        area += weighted
        mx += weighted * (n @ xy[:, 0])
        my += weighted * (n @ xy[:, 1])
    return area, mx, my


def t6_area(nodes: dict, tri: tuple) -> float:
    """The area of a 6-node triangle, integrated over its curved shape."""
    return t6_integrals(nodes, tri)[0]


def t6_first_moments(nodes: dict, tri: tuple) -> tuple[float, float]:
    """(integral of x dA, integral of y dA) of a 6-node triangle: its area times its centroid."""
    return t6_integrals(nodes, tri)[1:]


def t6_totals(nodes: dict, triangles) -> tuple[float, float, float]:
    """(area, integral of x dA, integral of y dA) summed over many triangles: one pass, vectorised."""
    triangles = list(triangles)
    if not triangles:
        return 0.0, 0.0, 0.0
    xy = np.array([[nodes[n] for n in tri] for tri in triangles])
    area = mx = my = 0.0
    for (n, d1, d2), weight in _RULE:
        det = (d1 @ xy[:, :, 0].T) * (d2 @ xy[:, :, 1].T) - (d1 @ xy[:, :, 1].T) * (d2 @ xy[:, :, 0].T)
        weighted = 0.5 * weight * det
        area += weighted.sum()
        mx += (weighted * (xy[:, :, 0] @ n)).sum()
        my += (weighted * (xy[:, :, 1] @ n)).sum()
    return float(area), float(mx), float(my)


# ============================================================================================== the magnets as a load
@dataclasses.dataclass(frozen=True)
class SlotFrame:
    """One magnet slot with the magnet's centrifugal force in the slot's own axes.
    centre            the slot's centre [mm] (the magnet sits in the middle of it)
    u, n              unit vectors along the slot and across it, each turned the way the force pushes (so f_along and f_across are positive)
    half_length, half_thickness   half the slot's sides [mm]: 4.64 and 0.885
    mass_kg           the real magnet's mass: 9.08 x 1.57 mm x the stack x the density of N45
    force_n           the centrifugal force on it [N] at the speed asked: m w^2 r, along the radius through the slot's centre
    f_across, f_along its parts along n (across the slot: the long wall) and along u (along the slot: the outer end wall)"""
    name: str
    pole: int
    centre: tuple
    u: tuple
    n: tuple
    half_length: float
    half_thickness: float
    mass_kg: float
    force_n: tuple
    f_across: float
    f_along: float


def slot_frames(p: mg.Params, regions, rpm: float, magnet: Material = N45) -> dict:
    """{slot name: SlotFrame} of the magnet slots of the regions, at `rpm`. The slot's axes come from its own outline (the bounding rectangle of least area of its polygon, whose
    sides are the slot's 9.28 and 1.77 mm: the fillets do not change it), not from the parameters that drew it."""
    frames = {}
    for r in regions:
        if not r.name.startswith("magnet_"):
            continue
        polygon = r.polygon()
        corners = list(polygon.minimum_rotated_rectangle.exterior.coords)[:4]
        sides = [(math.dist(corners[i], corners[(i + 1) % 4]), np.array(corners[(i + 1) % 4]) - np.array(corners[i])) for i in range(4)]
        long_side = max(sides, key=lambda s: s[0])
        short_side = min(sides, key=lambda s: s[0])
        u = long_side[1] / np.linalg.norm(long_side[1])
        n = np.array([-u[1], u[0]])
        centre = np.array([polygon.centroid.x, polygon.centroid.y])
        mass = magnet.density_kg_m3 * (p.mag_w * 1e-3) * (p.mag_t * 1e-3) * (p.stack * 1e-3)
        force = mass * omega(rpm) ** 2 * centre * 1e-3
        if force @ u < 0:
            u = -u
        if force @ n < 0:
            n = -n
        frames[r.name] = SlotFrame(name=r.name, pole=int(r.name.split("_")[1]), centre=tuple(centre), u=tuple(u), n=tuple(n), half_length=long_side[0] / 2,
                                   half_thickness=short_side[0] / 2, mass_kg=mass, force_n=tuple(force), f_across=float(force @ n), f_along=float(force @ u))
    return frames


def edge_nodal_forces(nodes: dict, edges, direction, total_force: float) -> dict:
    """The consistent nodal forces of a uniform line load on quadratic edges: `total_force` [N] along `direction` (normalised) spread evenly over the length of the edges, each
    3-node edge (end, end, mid) giving 1/6, 1/6 and 4/6 of its share. {node: (fx, fy)}; a node that two edges share gets both."""
    d = np.array(direction, float)
    d /= np.hypot(*d)
    lengths = [math.dist(nodes[e[0]], nodes[e[1]]) for e in edges]
    per_mm = total_force / sum(lengths)
    out: dict = {}
    for edge, length in zip(edges, lengths):
        for node, weight in ((edge[0], 1 / 6), (edge[1], 1 / 6), (edge[2], 4 / 6)):
            f = d * per_mm * length * weight
            old = out.get(node, (0.0, 0.0))
            out[node] = (old[0] + f[0], old[1] + f[1])
    return out


def slot_loads(p: mg.Params, mesh: Mesh, frames: dict) -> dict:
    """The pressure model's loads: for every slot the two walls the magnet bears on, each with its edges, length, pressure and nodal forces.
    {slot: {"long": wall, "end": wall}}; a wall = {"edges", "length_mm", "force_n", "direction", "pressure_mpa" (force / (length x stack)), "nodal"}.
    The long wall is the flat part (between the fillets) of the long side on the n side of the slot, the end wall the flat part of the short side on its u side; an edge belongs
    to a wall when all three of its nodes lie on it."""
    tol = 1e-6
    out = {}
    for name, frame in frames.items():
        c, u, n = np.array(frame.centre), np.array(frame.u), np.array(frame.n)
        walls = {"long": [], "end": []}
        for edge in mesh.boundary[f"slot:{name}"]:
            local = [((np.array(mesh.nodes[k]) - c) @ u, (np.array(mesh.nodes[k]) - c) @ n) for k in edge]
            if all(abs(t - frame.half_thickness) < tol and abs(s) <= frame.half_length - p.slot_fillet + tol for s, t in local):
                walls["long"].append(edge)
            elif all(abs(s - frame.half_length) < tol and abs(t) <= frame.half_thickness - p.slot_fillet + tol for s, t in local):
                walls["end"].append(edge)
        entry = {}
        for kind, direction, force in (("long", n, frame.f_across), ("end", u, frame.f_along)):
            edges = walls[kind]
            length = sum(math.dist(mesh.nodes[e[0]], mesh.nodes[e[1]]) for e in edges)
            entry[kind] = {"edges": edges, "length_mm": length, "force_n": force, "direction": tuple(direction), "pressure_mpa": force / (length * p.stack),
                           "nodal": edge_nodal_forces(mesh.nodes, edges, direction, force)}
        out[name] = entry
    return out


def nodal_forces(loads: dict) -> dict:
    """All the walls' nodal forces together: {node: (fx, fy)}."""
    total: dict = {}
    for walls in loads.values():
        for wall in walls.values():
            for node, (fx, fy) in wall["nodal"].items():
                old = total.get(node, (0.0, 0.0))
                total[node] = (old[0] + fx, old[1] + fy)
    return total


def wall_sets(loads: dict) -> dict:
    """The node sets of the loaded walls as the deck names them: LONG_<slot> and END_<slot> (PrePoMax makes its pressure surfaces from them)."""
    out = {}
    for name, walls in loads.items():
        for kind, wall in walls.items():
            out[f"{kind.upper()}_{name}"] = sorted({n for e in wall["edges"] for n in e})
    return out


def mass_equivalent_density(p: mg.Params, regions, magnet: Material = N45) -> float:
    """The density [tonne/mm3] of the body that fills a slot in the bonded model: N45's, scaled so that the body has the mass of the REAL magnet (9.08 x 1.57 mm), not of the
    slot (9.28 x 1.77 mm with its fillets, 14.7 % more volume)."""
    slot_area = float(np.mean([r.polygon().area for r in regions if r.name.startswith("magnet_")]))
    return magnet.density_tonne_mm3 * (p.mag_w * p.mag_t) / slot_area


# ============================================================================================== the deck
def deck_text(mesh: Mesh, variant: str, rpm: float, forces: dict | None = None, node_sets: dict | None = None, thickness: float = 20.0,
              steel: Material = M19, magnet: Material = N45, magnet_density_tonne_mm3: float | None = None, soft_youngs_mpa: float | None = None) -> str:
    """The CalculiX input deck of one model at one speed, as text (ASCII, mm, N, tonne, s).

    variant "pressure": the steel only, `forces` ({node: (fx, fy)}, the pressures on the slot walls as nodal forces) as *CLOAD.
    variant "bonded":   the steel and the slot bodies; each slot body gets its own copy of every node on the slot wall (numbered above the largest), tied to the steel's node by two
                        *EQUATION lines each; the magnet's density is `magnet_density_tonne_mm3` (the mass of the real magnet spread over the slot).
    variant "soft":     the bonded deck with the slot bodies' Young's modulus cut to `soft_youngs_mpa` (SOFT_FILLER_YOUNGS_MPA if not given) and their mass kept: the diagnostic that
                        separates the slots' lost stiffness from the way the load reaches the steel.
    All: *DLOAD CENTRIF about the z axis, y fixed on the x-axis edge and x on the y-axis edge, the step's output (displacements, element stresses, the reactions of the two edges)
    and, when given, a node set for each of `node_sets` ({name: node ids}: the loaded walls). Every number goes through deck_number."""
    if variant not in ALL_VARIANTS:
        raise ValueError(f"variant must be one of {ALL_VARIANTS}, not {variant!r}")
    tied = variant in TIED_VARIANTS
    w2 = omega(rpm) ** 2
    steel_nodes = {n for t in mesh.steel for n in t}
    wall = {n for name, rows in mesh.boundary.items() if name.startswith("slot:") for row in rows for n in row}
    next_id = max(mesh.nodes) + 1
    duplicate: dict = {}
    magnet_rows: list = []
    if tied:
        for name in sorted(mesh.magnets):
            for t in mesh.magnets[name]:
                row = []
                for n in t:
                    if n in wall:
                        if (name, n) not in duplicate:
                            duplicate[(name, n)] = next_id
                            next_id += 1
                        row.append(duplicate[(name, n)])
                    else:
                        row.append(n)
                magnet_rows.append(tuple(row))
    used = set(steel_nodes) | ({n for t in magnet_rows for n in t} - set(duplicate.values()) if tied else set())
    coordinates = {n: mesh.nodes[n] for n in used}
    coordinates.update({d: mesh.nodes[n] for (_, n), d in duplicate.items()})
    out = ["*HEADING", f"B4 rotor sector, plane stress CPS6, {variant} magnets, {rpm:g} rpm", "*NODE"]
    out += [f"{n}, {deck_number(coordinates[n][0])}, {deck_number(coordinates[n][1])}, 0.0" for n in sorted(coordinates)]
    out.append("*ELEMENT, TYPE=CPS6, ELSET=STEEL")
    element = 0
    for t in mesh.steel:
        element += 1
        out.append(f"{element}, " + ", ".join(map(str, t)))
    if tied:
        out.append("*ELEMENT, TYPE=CPS6, ELSET=MAGNET")
        for t in magnet_rows:
            element += 1
            out.append(f"{element}, " + ", ".join(map(str, t)))

    def node_set(name: str, ids) -> list[str]:
        ids = sorted(ids)
        return [f"*NSET, NSET={name}"] + [", ".join(map(str, ids[k:k + 10])) for k in range(0, len(ids), 10)]

    out += node_set("EDGE0", [n for n in steel_nodes if abs(mesh.nodes[n][1]) < 1e-9])
    out += node_set("EDGE90", [n for n in steel_nodes if abs(mesh.nodes[n][0]) < 1e-9])
    for name, ids in (node_sets or {}).items():
        out += node_set(name, ids)
    out += ["*MATERIAL, NAME=M19", "*ELASTIC", f"{deck_number(steel.youngs_mpa)}, {deck_number(steel.poisson)}", "*DENSITY", deck_number(steel.density_tonne_mm3)]
    body_name = "FILLER" if variant == SOFT else "N45"
    if tied:
        density = magnet.density_tonne_mm3 if magnet_density_tonne_mm3 is None else magnet_density_tonne_mm3
        body_youngs = (SOFT_FILLER_YOUNGS_MPA if soft_youngs_mpa is None else soft_youngs_mpa) if variant == SOFT else magnet.youngs_mpa
        out += [f"*MATERIAL, NAME={body_name}", "*ELASTIC", f"{deck_number(body_youngs)}, {deck_number(magnet.poisson)}", "*DENSITY", deck_number(density)]
    out += ["*SOLID SECTION, ELSET=STEEL, MATERIAL=M19", deck_number(thickness)]
    if tied:
        out += [f"*SOLID SECTION, ELSET=MAGNET, MATERIAL={body_name}", deck_number(thickness), "*EQUATION"]
        for (_, n), d in sorted(duplicate.items(), key=lambda kv: kv[1]):
            for dof in (1, 2):
                out += ["2", f"{n}, {dof}, {deck_number(1.0)}, {d}, {dof}, {deck_number(-1.0)}"]
    # The boundary conditions are the step's own and the axis is written as whole numbers (0, 0, 0, 0, 0, 1), as PrePoMax writes it: PrePoMax can then import this deck whole
    # (it reads *BOUNDARY inside a step, and could not read "0." in the *DLOAD line); for ccx both ways are the same.
    out += ["*STEP", "*STATIC", "*BOUNDARY", "EDGE0, 2, 2", "EDGE90, 1, 1", "*DLOAD", f"STEEL, CENTRIF, {deck_number(w2)}, 0, 0, 0, 0, 0, 1"]
    if tied:
        out.append(f"MAGNET, CENTRIF, {deck_number(w2)}, 0, 0, 0, 0, 0, 1")
    if forces and not tied:
        out.append("*CLOAD")
        for node, (fx, fy) in sorted(forces.items()):
            out += [f"{node}, 1, {deck_number(fx)}", f"{node}, 2, {deck_number(fy)}"]
    out += ["*NODE FILE", "U", "*EL FILE", "S", "*NODE PRINT, NSET=EDGE0, TOTALS=ONLY", "RF", "*NODE PRINT, NSET=EDGE90, TOTALS=ONLY", "RF", "*END STEP"]
    return "\n".join(out) + "\n"


def external_load(mesh: Mesh, variant: str, rpm: float, forces: dict | None = None, thickness: float = 20.0, steel: Material = M19,
                  magnet_density_tonne_mm3: float | None = None, magnet: Material = N45) -> tuple[float, float]:
    """The total external load [N] on the sector, worked out from the mesh: the body force rho w^2 t (integral of x dA, integral of y dA) of the steel (and, bonded or soft, of the
    slot bodies, whose density is the same: the stiffness does not enter a load) over the curved elements, plus the applied nodal forces `forces`. The reactions of the symmetry
    edges must balance it."""
    w2 = omega(rpm) ** 2
    _, mx, my = t6_totals(mesh.nodes, mesh.steel)
    fx, fy = steel.density_tonne_mm3 * w2 * thickness * mx, steel.density_tonne_mm3 * w2 * thickness * my
    if variant in TIED_VARIANTS:
        density = magnet.density_tonne_mm3 if magnet_density_tonne_mm3 is None else magnet_density_tonne_mm3
        _, mx, my = t6_totals(mesh.nodes, [t for rows in mesh.magnets.values() for t in rows])
        fx += density * w2 * thickness * mx
        fy += density * w2 * thickness * my
    for _, (gx, gy) in (forces or {}).items():
        fx += gx
        fy += gy
    return fx, fy


# ============================================================================================== reading the answer
REACTION = re.compile(r"total force \(fx,fy,fz\) for set (\w+) and time\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)")


def parse_reactions(text: str) -> dict:
    """{node set: (fx, fy, fz)} of the `*NODE PRINT, TOTALS=ONLY` blocks of a .dat file; {} if there are none."""
    return {m.group(1): (float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in REACTION.finditer(text)}


def reaction_vector(reactions: dict) -> tuple[float, float]:
    """(RF_x of EDGE90, RF_y of EDGE0): the constrained components, the only ones that are reactions (the other component of a set's total is what the body force puts on its
    nodes). They are minus the external load."""
    return reactions["EDGE90"][0], reactions["EDGE0"][1]


def von_mises(components) -> float:
    """The von Mises stress of (sxx, syy, szz, sxy, syz, szx), the order of CalculiX's .frd."""
    sx, sy, sz, txy, tyz, tzx = components
    return math.sqrt(0.5 * ((sx - sy) ** 2 + (sy - sz) ** 2 + (sz - sx) ** 2) + 3 * (txy ** 2 + tyz ** 2 + tzx ** 2))


def peak_of(stress: dict, nodes) -> tuple:
    """(the largest von Mises stress, its node) among `nodes` of the {node: components} the .frd gave; (None, None) if none of them has a value."""
    best, where = None, None
    for n in nodes:
        if n in stress:
            v = von_mises(stress[n])
            if best is None or v > best:
                best, where = v, n
    return best, where


def body_peak(mesh: Mesh, stress: dict) -> float:
    """The largest von Mises stress [MPa] over the nodes of the slot bodies of a bonded or soft run: every node that no steel triangle uses, that is the bodies' interior nodes and the
    copies of the wall nodes (which carry numbers above the mesh's). The steel's nodes, the wall nodes included, are not counted."""
    steel = {n for t in mesh.steel for n in t}
    return max(von_mises(s) for n, s in stress.items() if n not in steel)


@functools.lru_cache(maxsize=8)
def _strips(p: mg.Params) -> tuple:
    """The strip geometries of every pole of the sector: ((pole, geometry), ...)."""
    return tuple((k, strip_geometry(p, k)) for k in range(2 * p.pole_pairs // p.n_sectors))


def classify_place(p: mg.Params, x: float, y: float) -> dict:
    """Where a point of the rotor is, in the words of the brief: {"x_mm", "y_mm", "r_mm", "angle_deg", "feature", "in_a_strip", "pole", "slot", "strip", "distance_to_strip_mm"}.
    feature is "shaft bore" (on the bore), "outer bridge" or "centre rib" (within STRIP_ZONE_MM of the cut across that strip), "slot wall" (on a slot but not near a strip), or
    "bulk". pole and slot are the nearest strip's pole and the nearest slot's name."""
    point = Point(x, y)
    r = math.hypot(x, y)
    nearest = min(((point.distance(line), k, key) for k, g in _strips(p) for key, line in g["cuts"].items()), key=lambda t: t[0])
    vis = _visible(p)
    slot = min((name for name in vis if name.startswith("magnet_")), key=lambda n: point.distance(vis[n]))
    in_a_strip = nearest[0] <= STRIP_ZONE_MM and r > p.r_shaft + 0.01
    if r <= p.r_shaft + 0.01:
        feature = "shaft bore"
    elif in_a_strip:
        feature = "centre rib" if nearest[2] == "rib" else "outer bridge"
    elif point.distance(vis[slot]) < 0.05:
        feature = "slot wall"
    else:
        feature = "bulk"
    return {"x_mm": x, "y_mm": y, "r_mm": r, "angle_deg": math.degrees(math.atan2(y, x)), "feature": feature, "in_a_strip": bool(in_a_strip), "pole": nearest[1],
            "slot": slot, "strip": nearest[2], "distance_to_strip_mm": nearest[0]}


def cut_force(nodes: dict, triangles, stress: dict, cut: LineString, normal, thickness: float, samples: int = SECTION_SAMPLES) -> tuple[float, float]:
    """The force [N] that the material on the far side of a cut exerts on the near side: thickness x the integral along the cut of (sigma n), n the unit `normal` pointing out
    of the near side. The stress is the nodal stress of the .frd interpolated linearly over the corner triangles of the mesh, read at the middle of `samples` equal
    pieces of the cut (a point that falls just outside the mesh, on a curved wall, takes the nearest valid one)."""
    import matplotlib.tri as mtri
    ids = sorted(nodes)
    index = {n: i for i, n in enumerate(ids)}
    xs, ys = np.array([nodes[n][0] for n in ids]), np.array([nodes[n][1] for n in ids])
    corners = np.array([[index[n] for n in t[:3]] for t in triangles])
    triangulation = mtri.Triangulation(xs, ys, corners)
    component = {}
    for name, k in (("sxx", 0), ("syy", 1), ("sxy", 3)):
        component[name] = mtri.LinearTriInterpolator(triangulation, np.array([stress[n][k] if n in stress else 0.0 for n in ids]))
    positions = (np.arange(samples) + 0.5) / samples
    px = np.array([cut.interpolate(f, normalized=True).x for f in positions])
    py = np.array([cut.interpolate(f, normalized=True).y for f in positions])
    values = {name: np.ma.filled(f(px, py), np.nan) for name, f in component.items()}
    good = np.isfinite(values["sxx"]) & np.isfinite(values["syy"]) & np.isfinite(values["sxy"])
    if not good.any():
        raise ValueError("the cut lies outside the mesh")
    nx, ny = np.array(normal, float) / np.hypot(*normal)
    tx = np.array([np.interp(i, np.flatnonzero(good), values["sxx"][good]) * nx + np.interp(i, np.flatnonzero(good), values["sxy"][good]) * ny for i in range(samples)])
    ty = np.array([np.interp(i, np.flatnonzero(good), values["sxy"][good]) * nx + np.interp(i, np.flatnonzero(good), values["syy"][good]) * ny for i in range(samples)])
    piece = cut.length / samples
    return float(thickness * tx.sum() * piece), float(thickness * ty.sum() * piece)


# ============================================================================================== the checks that are just a number against a number
def scale_check(label: str, low: float, high: float, **kw) -> dict:
    """The brief's: the stress at 12,000 rpm is (12/10)^2 = 1.44 times the stress at 10,000, within 1 %."""
    return compare(f"{label}: the peak at 12,000 rpm is (12/10)^2 times the peak at 10,000 rpm", high, EXPECTED_RATIO * low, SCALE_TOLERANCE_PERCENT,
                   "the brief: (12/10)^2 = 1.44 times the 10,000 rpm peak (linear elasticity: the load is rho w^2 r)", "MPa", **kw)


def slot_detail(p: mg.Params, frame: SlotFrame, x: float, y: float) -> str:
    """Where on a slot's outline a point is, in words, from its position in the slot's own axes (s along the slot, t across it, both measured from the centre, turned the way the
    magnet's force pushes): "corner fillet at the outer end, on the loaded side", "outer end wall", "long wall, loaded side", ... The loaded long wall and the outer end wall
    are the two the magnet bears on; the outer end is the end nearer the rim."""
    c, u, n = np.array(frame.centre), np.array(frame.u), np.array(frame.n)
    s, t = (np.array([x, y]) - c) @ u, (np.array([x, y]) - c) @ n
    near_end = abs(s) > frame.half_length - p.slot_fillet - 0.05
    near_wall = abs(t) > frame.half_thickness - p.slot_fillet - 0.05
    end = "outer end" if s > 0 else "inner end"
    side = "on the loaded side" if t > 0 else "on the other side"
    if near_end and near_wall:
        return f"corner fillet at the {end}, {side}"
    if near_end:
        return f"{end} wall"
    if near_wall:
        return f"long wall, {'loaded side' if t > 0 else 'other side'}"
    return "away from the slot's outline"


def point_segment_distance(points: np.ndarray, a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """The distance of each of `points` (n, 2) from the segment a-b."""
    ab = b - a
    t = np.clip(((points - a) @ ab) / (ab @ ab), 0.0, 1.0)
    return np.hypot(*(points - (a + t[:, None] * ab)).T)


# ============================================================================================== one model: build, solve, read back
def build_model(p: mg.Params, h_min: float, msh_path=None) -> dict:
    """The inputs of one mesh size: motorgeo's regions, the mesh, and the density of the bonded body. B5 builds its own copy of the case through here, in its own folder."""
    regions = mg.build(p, gap="sliding")
    mesh = mesh_rotor(p, regions, h_min, msh_path=msh_path)
    return {"p": p, "regions": regions, "mesh": mesh, "magnet_density": mass_equivalent_density(p, regions)}


def loads_at(model: dict, rpm: float) -> dict:
    """The magnets' frames, wall loads and nodal forces at `rpm` (the forces grow with the speed squared)."""
    frames = slot_frames(model["p"], model["regions"], rpm)
    loads = slot_loads(model["p"], model["mesh"], frames)
    return {"frames": frames, "loads": loads, "nodal": nodal_forces(loads)}


def write_deck(folder, model: dict, variant: str, rpm: float) -> Path:
    """job.inp of one model at one speed in `folder`; returns its path."""
    at = loads_at(model, rpm)
    text = deck_text(model["mesh"], variant, rpm, forces=at["nodal"], node_sets=wall_sets(at["loads"]), thickness=model["p"].stack,
                     magnet_density_tonne_mm3=model["magnet_density"])
    path = Path(folder) / "job.inp"
    Path(folder).mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="ascii", newline="\n")
    return path


def solve(model: dict, variant: str, rpm: float, folder: Path, log: Log, label: str) -> dict:
    """Write the deck and run ccx on it through run_program (the solver slot, below normal, CCX_THREADS threads, a timeout); read the .frd and .dat back.
    Returns {"folder", "variant", "rpm", "h_min", "returncode", "timed_out", "finished", "seconds", "version", "stress" ({node: components} or None), "reactions"}."""
    folder = Path(folder)
    write_deck(folder, model, variant, rpm)
    started = time.perf_counter()
    result = run_program(log, [paths.tool("ccx"), "-i", "job"], label=label, cwd=folder, timeout_s=CCX_TIMEOUT_S, threads=CCX_THREADS, name="ccx")
    banner = re.search(r"CalculiX Version\s+(\d+\.\d+)", result.stdout)
    finished = (not result.timed_out) and result.returncode == 0 and "Job finished" in result.stdout and (folder / "job.frd").is_file() and (folder / "job.dat").is_file()
    out = {"folder": str(folder), "variant": variant, "rpm": rpm, "h_min": model["mesh"].h_min, "returncode": result.returncode, "timed_out": result.timed_out,
           "finished": finished, "seconds": round(time.perf_counter() - started, 2), "version": banner.group(1) if banner else None, "stress": None, "reactions": {}}
    if finished:
        out["stress"] = ccx_check.read_frd_stress(folder / "job.frd")
        out["reactions"] = parse_reactions((folder / "job.dat").read_text(encoding="ascii", errors="replace"))
    return out


def bridge_zone_peaks(p: mg.Params, mesh: Mesh, stress: dict, nodes) -> dict:
    """The largest von Mises stress of the steel in each strip's zone (within STRIP_ZONE_MM of its cut): {"magnet_0_dn": {...}, ..., "rib_0": {...}}, each with the value, the
    node and its place. The bridge of a pole's lower slot is bridge_low, of its upper slot bridge_high."""
    ids = [n for n in nodes if n in stress]
    xy = np.array([mesh.nodes[n] for n in ids])
    mises = np.array([von_mises(stress[n]) for n in ids])
    out = {}
    for k, g in _strips(p):
        for key, name in (("bridge_low", f"magnet_{k}_dn"), ("bridge_high", f"magnet_{k}_up"), ("rib", f"rib_{k}")):
            a, b = (np.array(c) for c in (g["cuts"][key].coords[0], g["cuts"][key].coords[-1]))
            inside = point_segment_distance(xy, a, b) <= STRIP_ZONE_MM
            i = int(np.argmax(np.where(inside, mises, -1.0)))
            out[name] = {"value": float(mises[i]), "node": ids[i], "x_mm": float(xy[i][0]), "y_mm": float(xy[i][1])}
    return out


def pole_sections(p: mg.Params, mesh: Mesh, stress: dict, rpm: float, at: dict) -> list[dict]:
    """The force through the three strips of each pole, from the stress: the cut across each strip, the normal pointing out of the cap, F = thickness x the integral of (sigma n).
    The three forces add up to what the rest of the rotor holds the cap with. `axial_pull_n` is minus their sum along the pole's axis (what the strips carry against the
    centrifugal pull, positive), `transverse_n` their sum across it (zero by symmetry). Also each strip's own share, and the two hand-estimate pulls (the note's, and the
    note's with the magnets' end-wall part left off the cap) in N over the stack."""
    out = []
    for k, g in _strips(p):
        axis = np.array([math.cos(math.radians(g["axis_deg"])), math.sin(math.radians(g["axis_deg"]))])
        cap = g["cap"]
        total, shares = np.zeros(2), {}
        for key, cut in g["cuts"].items():
            a, b = (np.array(c) for c in (cut.coords[0], cut.coords[-1]))
            tangent = (b - a) / np.linalg.norm(b - a)
            normal = np.array([-tangent[1], tangent[0]])
            mid = np.array(cut.interpolate(0.5, normalized=True).coords[0])
            if normal @ (mid - np.array([cap.centroid.x, cap.centroid.y])) < 0:
                normal = -normal                              # pointing out of the cap
            force = np.array(cut_force(mesh.nodes, mesh.steel, stress, cut, normal, p.stack))
            total += force
            shares[key] = {"force_n": tuple(force), "axial_pull_n": float(-(force @ axis)), "width_mm": g["widths"][key]}
        steel_term = M19.density_tonne_mm3 * omega(rpm) ** 2 * p.stack * cap.area * (np.array([cap.centroid.x, cap.centroid.y]) @ axis)
        slots = [name for name in at["frames"] if at["frames"][name].pole == k]
        hand = steel_term + sum(float(np.array(at["frames"][name].force_n) @ axis) for name in slots)
        on_cap = 0.0
        for name in slots:
            for wall in at["loads"][name].values():
                points = np.array([mesh.nodes[n] for e in wall["edges"] for n in e[:2]])
                if cap.distance(Point(*points.mean(axis=0))) < 1e-3:
                    on_cap += wall["force_n"] * float(np.array(wall["direction"]) @ axis)
        carried = float(-(total @ axis))
        widths = sum(g["widths"].values())
        out.append({"pole": k, "axis_deg": g["axis_deg"], "strips": shares, "axial_pull_n": carried, "transverse_n": float(total @ np.array([-axis[1], axis[0]])),
                    "nominal_mpa": carried / p.stack / widths, "hand_pull_n": float(hand), "refined_pull_n": float(steel_term + on_cap), "cap_area_mm2": cap.area,
                    "widths_mm": g["widths"]})
    return out


def analyse(model: dict, run: dict, rpm: float) -> dict:
    """The numbers of one solve: the peak of the steel and its place, the peak in each strip's zone, the bore's hoop stress against Lame, the equilibrium (the reactions against the
    load worked out from the mesh), the strips' forces. A solve that did not finish has only {"finished": False}."""
    if not run["finished"]:
        return {"finished": False}
    p, mesh, variant = model["p"], model["mesh"], run["variant"]
    at = loads_at(model, rpm)
    stress = run["stress"]
    steel_nodes = {n for t in mesh.steel for n in t}
    value, node = peak_of(stress, steel_nodes)
    x, y = mesh.nodes[node]
    place = classify_place(p, x, y)
    place["slot_detail"] = slot_detail(p, at["frames"][place["slot"]], x, y) if place["feature"] in ("outer bridge", "centre rib", "slot wall") else None
    hoop = ccx_check.bore_hoop(mesh.nodes, stress)
    lame = ccx_check.lame_hoop(ccx_check.BORE) * (rpm / 10_000) ** 2
    fx, fy = external_load(mesh, variant, rpm, at["nodal"] if variant == "pressure" else None, p.stack, magnet_density_tonne_mm3=model["magnet_density"])
    rx, ry = reaction_vector(run["reactions"])
    zones = bridge_zone_peaks(p, mesh, stress, steel_nodes)
    bridge_values = [v["value"] for k, v in zones.items() if k.startswith("magnet_")]
    return {"finished": True, "peak": {"value_mpa": value, "node": node, "place": place}, "zone_peaks": zones, "bridge_peak_mpa": max(bridge_values),
            "bridge_spread_percent": 100 * (max(bridge_values) - min(bridge_values)) / (sum(bridge_values) / len(bridge_values)),
            "bore": {"hoop_mean_mpa": sum(hoop) / len(hoop), "nodes": len(hoop), "lame_mpa": lame}, "sections": pole_sections(p, mesh, stress, rpm, at),
            "equilibrium": {"external_n": (fx, fy), "reaction_n": (rx, ry), "residual_n": math.hypot(fx + rx, fy + ry), "residual_percent": 100 * math.hypot(fx + rx, fy + ry) / math.hypot(fx, fy)},
            "total_load_n": math.hypot(fx, fy)}


# ============================================================================================== the pictures
def render_cgx(folder: Path, log: Log) -> dict:
    """cgx's own hardcopy of the von Mises stress of the solve in `folder` (job.frd): `read`, `ds 2 e 7` (dataset 2 is the stress, entity 7 its von Mises), `plot fv all` (faces
    with values), `rot z 180` (the view along z), `hcpy tga`. cgx writes three TGA frames (0__.tga the legend and the frame, 1__.tga the model, 2__.tga the axes) and then calls
    mogrify and composite (ImageMagick), which are not installed: Python composes the frames. cgx does not return to its script after `hcpy`, so it is ended by run_solver's
    timeout once the frames are there; the record says so."""
    folder = Path(folder)
    (folder / "shot.fbd").write_text("read job.frd\nds 2 e 7\nplot fv all\nrot z 180\nhcpy tga\nquit\n", encoding="ascii", newline="\n")
    started = time.perf_counter()
    result = run_program(log, [paths.tool("cgx"), "-b", "shot.fbd"], label="bench B4 cgx hardcopy", cwd=folder, timeout_s=CGX_TIMEOUT_S, threads=1, name="cgx")
    frames = [folder / f"{i}__.tga" for i in range(3)]
    return {"frames_made": all(f.is_file() and f.stat().st_size > 1000 for f in frames[:2]), "returncode": result.returncode, "ended_by_timeout": result.timed_out,
            "seconds": round(time.perf_counter() - started, 2), "max_von_mises": (lambda m: float(m.group(1)) if m else None)(re.search(r"max:(\S+) at node", result.stdout)),
            "script": "read job.frd; ds 2 e 7; plot fv all; rot z 180; hcpy tga; quit"}


def cgx_picture(folder: Path, caption: str, target: Path) -> dict:
    """Put cgx's frames together: the model frame (its colours divided by cgx's lighting factor, after checking that they are the legend's colours times it) beside the legend,
    under `caption`. Returns {"made", "colours_match" (the share of the model frame's colours that are a legend colour times CGX_LIGHT), "legend_colours"}."""
    from PIL import Image
    folder = Path(folder)
    window = np.asarray(Image.open(folder / "0__.tga").convert("RGB"))
    model = np.asarray(Image.open(folder / "1__.tga").convert("RGB")).astype(float)
    swatches = window[:, :30].reshape(-1, 3)                                   # the colour bar is at the left of the window
    legend = np.unique(swatches[(swatches.min(axis=1) < 250) & (swatches.max(axis=1) > 20)], axis=0).astype(float)
    colours = np.unique(model.reshape(-1, 3).astype(int), axis=0)
    colours = colours[(colours.max(axis=1) > 8) & (colours.min(axis=1) < 250)]     # not the background, not the black of the outline and the text
    distance = np.abs(colours[:, None, :] - CGX_LIGHT * legend[None, :, :]).max(axis=2).min(axis=1) if len(legend) and len(colours) else np.array([])
    share = float((distance <= 2.5).mean()) if len(distance) else 0.0
    if share < 0.9:
        return {"made": False, "colours_match": share, "legend_colours": int(len(legend))}
    lit = np.clip(model / CGX_LIGHT, 0, 255).astype(np.uint8)
    ink = lambda image: np.where((image < 250).any(axis=2))
    rows, cols = ink(lit)
    lit = lit[max(rows.min() - 12, 0):rows.max() + 12, max(cols.min() - 12, 0):cols.max() + 12]
    bar_rows, bar_cols = ink(window[:, :172])
    bar = window[max(bar_rows.min() - 6, 0):bar_rows.max() + 6, :172]
    legend_and_model_figure(bar, lit, caption, target)
    return {"made": True, "colours_match": share, "legend_colours": int(len(legend))}


def fig_stress_full(model: dict, stress: dict, title: str, target: Path, marks=()) -> None:
    """The whole rotor from the sector (four copies turned by 90 degrees), the steel's von Mises stress on every triangle (turbo, 0 to its largest), the slots as holes, the
    three cuts of each pole as thin lines, and the largest value of each slot's bridge marked. The picture of the case where cgx's render is not available, and a second look
    when it is."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    mesh, p = model["mesh"], model["p"]
    ids = sorted({n for t in mesh.steel for n in t if n in stress})
    index = {n: i for i, n in enumerate(ids)}
    x, y = np.array([mesh.nodes[n][0] for n in ids]), np.array([mesh.nodes[n][1] for n in ids])
    values = np.array([von_mises(stress[n]) for n in ids])
    corners = np.array([[index[n] for n in t[:3]] for t in mesh.steel])
    fig, ax = plt.subplots(figsize=(9, 9))
    image = None
    for k in range(4):
        c, s = math.cos(math.radians(90 * k)), math.sin(math.radians(90 * k))
        image = ax.tripcolor(mtri.Triangulation(c * x - s * y, s * x + c * y, corners), values, shading="gouraud", cmap="turbo", vmin=0, vmax=float(values.max()), rasterized=True)
        if k == 0:
            for _, g in _strips(p):
                for line in g["cuts"].values():
                    ax.plot(*line.xy, color="white", lw=0.9)
            for px, py in marks:
                ax.plot(px, py, "o", mfc="none", mec="white", ms=9, mew=1.2)
    fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02).set_label("von Mises stress [MPa]")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(wrap_title(title), fontsize=10.5, color="#1d1d1f", loc="left")
    fig.savefig(target, dpi=150, bbox_inches="tight")
    plt.close(fig)


def figure_caption(peak: dict, nodes: int) -> str:
    """The caption of figure.png (cgx's render): five lines, none longer than bench.TITLE_CHARS (a longer one wraps and leaves a word alone on a line). It says whose picture it is, what it
    shows (`peak` is the verdict run's {"value_mpa", "place"}, `nodes` the model's), which way cgx's view runs, and what Python did to cgx's frames."""
    return (f"CalculiX 2.23: the band-free rotor's 90-degree sector at 10,000 rpm, cgx's own render\n"
            f"von Mises stress [MPa] of the steel, magnets as pressure on the slot walls, {nodes:,} nodes\n"
            f"peak {peak['value_mpa']:.1f} MPa in the {peak['place']['feature']} of {peak['place']['slot']}\n"
            f"cgx's view along z: the x axis runs up the right edge, the bore is at the lower right\n"
            f"Python recomposed cgx's frames, dividing the model's colours by cgx's lighting factor {CGX_LIGHT}")


def stress_full_title(peak: dict, nodes: int) -> str:
    """The title of stress_full.png: five lines, none longer than bench.TITLE_CHARS. The white lines are the cuts of both poles of the sector (the rib and the two bridges of each)."""
    return (f"CalculiX 2.23: the band-free rotor at 10,000 rpm (the 90-degree sector repeated four times)\n"
            f"von Mises stress of the steel, magnets as pressure on the slot walls\n"
            f"peak {peak['value_mpa']:.1f} MPa; {nodes:,} nodes in the sector (h_min {H_MIN_MAIN} mm)\n"
            f"circles: the largest value in each slot's bridge\n"
            f"white lines: the three cuts (the rib and the two bridges) of each of the sector's two poles")


# ============================================================================================== the case
def steel_node_count(mesh: Mesh) -> int:
    """The nodes of the steel's triangles: the nodes of the pressure model (its slots are voids), where `len(mesh.nodes)` also counts the slots' own, which only the bonded model uses."""
    return len({n for t in mesh.steel for n in t})


def mesh_counts(mesh: Mesh) -> dict:
    return {"h_min_mm": mesh.h_min, "h_max_mm": mesh.h_max, "nodes": len(mesh.nodes), "nodes_steel": steel_node_count(mesh), "triangles_steel": len(mesh.steel),
            "triangles_magnets": sum(len(v) for v in mesh.magnets.values()), "min_sicn": mesh.quality["min_sicn"], "mean_sicn": mesh.quality["mean_sicn"],
            "worst_element_at_mm": mesh.quality["worst"]}


def run(folder: Path, log: Log) -> dict:
    """Run B4 in `folder` (already empty). Returns the results dict whose "checks" list every verdict."""
    folder = Path(folder)
    started = time.perf_counter()
    p = mg.Params(barrier="none")
    note_rotor = mg.Params()                                   # the rotor of the note's hand estimate: with its air band (NOT a bench geometry: only the estimate's function is applied to it)
    hand_note = hand_estimate(note_rotor, 10_000)
    hand = {rpm: hand_estimate(p, rpm) for rpm in SPEEDS_RPM}
    log.say(f"B4 hand estimate, the note's rotor with its band: cap {hand_note['cap_area_mm2']:.1f} mm2 at {hand_note['cap_centroid_mm']:.1f} mm, pull {hand_note['pull_kn_per_m']:.1f} kN/m, "
            f"strips {hand_note['strips_mm']['total']:.1f} mm, nominal {hand_note['nominal_mpa']:.2f} MPa (the note: 15); the band-free rotor: cap {hand[10_000]['cap_area_mm2']:.1f} mm2, "
            f"strips {hand[10_000]['strips_mm']['total']:.2f} mm, nominal {hand[10_000]['nominal_mpa']:.2f} MPa")
    # ---- the main mesh, and the four solves on it
    main = build_model(p, H_MIN_MAIN, msh_path=folder / "rotor.msh")
    q = main["mesh"].quality
    log.say(f"B4 mesh h_min {H_MIN_MAIN} mm: {len(main['mesh'].nodes):,} nodes ({steel_node_count(main['mesh']):,} in the steel), {len(main['mesh'].steel):,} steel + "
            f"{sum(len(v) for v in main['mesh'].magnets.values()):,} magnet 6-node triangles, "
            f"min SICN {q['min_sicn']:.3f} at ({q['worst']['x_mm']:.2f}, {q['worst']['y_mm']:.2f}) mm")
    runs: dict = {}
    analysis: dict = {}
    for variant in VARIANTS:
        for rpm in SPEEDS_RPM:
            key = f"{variant}_{rpm}_h{H_MIN_MAIN:g}"
            runs[key] = solve(main, variant, rpm, folder / key, log, f"bench B4 ccx {variant} {rpm} rpm h {H_MIN_MAIN:g}")
            analysis[key] = analyse(main, runs[key], rpm)
            a = analysis[key]
            if a["finished"]:
                log.say(f"B4 {key}: ccx {runs[key]['version']} {runs[key]['seconds']:.1f} s; peak {a['peak']['value_mpa']:.2f} MPa at ({a['peak']['place']['x_mm']:.2f}, {a['peak']['place']['y_mm']:.2f}) mm: "
                        f"{a['peak']['place']['feature']}; bridges {a['bridge_peak_mpa']:.2f} MPa; bore {a['bore']['hoop_mean_mpa']:.3f} MPa (Lame {a['bore']['lame_mpa']:.3f}); "
                        f"equilibrium {a['equilibrium']['residual_percent']:.5f} %")
            else:
                log.say(f"B4 {key}: ccx did not finish (exit {runs[key]['returncode']}, timed out {runs[key]['timed_out']}); see {folder / key}")
    # ---- the diagnostic (reported, never gating): the bonded deck with the slot bodies made soft, at 10,000 rpm
    soft_key = f"{SOFT}_{SPEEDS_RPM[0]}_h{H_MIN_MAIN:g}"
    soft_run = solve(main, SOFT, SPEEDS_RPM[0], folder / soft_key, log, f"bench B4 ccx {SOFT} {SPEEDS_RPM[0]} rpm h {H_MIN_MAIN:g} (the slot bodies' E {SOFT_FILLER_YOUNGS_MPA:g} MPa)")
    soft_analysis = analyse(main, soft_run, SPEEDS_RPM[0])
    if soft_analysis["finished"]:
        log.say(f"B4 {soft_key}: ccx {soft_run['version']} {soft_run['seconds']:.1f} s; bore {soft_analysis['bore']['hoop_mean_mpa']:.3f} MPa (Lame {soft_analysis['bore']['lame_mpa']:.3f}); "
                f"bridges {soft_analysis['bridge_peak_mpa']:.2f} MPa; equilibrium {soft_analysis['equilibrium']['residual_percent']:.5f} %")
    else:
        log.say(f"B4 {soft_key}: ccx did not finish (exit {soft_run['returncode']}, timed out {soft_run['timed_out']}); see {folder / soft_key}")
    # ---- the mesh sizes: the verdict model at 10,000 rpm on a coarser and a finer mesh
    levels = {H_MIN_MAIN: {"key": f"{VERDICT_VARIANT}_10000_h{H_MIN_MAIN:g}", "nodes": steel_node_count(main["mesh"]), "triangles": len(main["mesh"].steel)}}
    for h in (H_MIN_COARSE, H_MIN_FINE):
        other = build_model(p, h)
        key = f"{VERDICT_VARIANT}_10000_h{h:g}"
        runs[key] = solve(other, VERDICT_VARIANT, SPEEDS_RPM[0], folder / key, log, f"bench B4 ccx {VERDICT_VARIANT} 10000 rpm h {h:g}")
        analysis[key] = analyse(other, runs[key], SPEEDS_RPM[0])
        levels[h] = {"key": key, "nodes": steel_node_count(other["mesh"]), "triangles": len(other["mesh"].steel)}
        if analysis[key]["finished"]:
            log.say(f"B4 {key}: {levels[h]['nodes']:,} nodes, ccx {runs[key]['seconds']:.1f} s, peak {analysis[key]['peak']['value_mpa']:.2f} MPa; bridges {analysis[key]['bridge_peak_mpa']:.2f} MPa")
        runs[key]["stress"] = None                              # the large dictionary of the finer meshes is not kept
    # ---- the pictures: cgx's own render of the verdict run, and the whole rotor from Python
    verdict_key = f"{VERDICT_VARIANT}_10000_h{H_MIN_MAIN:g}"
    verdict = analysis[verdict_key]
    figure_ok = False
    cgx = {"frames_made": False}
    picture = {"made": False}
    if verdict["finished"]:
        peak = verdict["peak"]
        nodes = steel_node_count(main["mesh"])
        cgx = render_cgx(folder / verdict_key, log)
        if cgx["frames_made"]:
            picture = cgx_picture(folder / verdict_key, figure_caption(peak, nodes), folder / "figure.png")
        marks = [(v["x_mm"], v["y_mm"]) for k, v in verdict["zone_peaks"].items() if k.startswith("magnet_")]
        fig_stress_full(main, runs[verdict_key]["stress"], stress_full_title(peak, nodes), folder / "stress_full.png", marks)
        if not picture["made"]:                                 # cgx's frames are not there (or not its colours): the picture is the one plotted from the .frd
            shutil.copyfile(folder / "stress_full.png", folder / "figure.png")
        figure_ok = (folder / "figure.png").is_file() and (folder / "figure.png").stat().st_size > 20_000
    results = {"case": CASE, "title": TITLE, "params": params(p), "verdict_variant": VERDICT_VARIANT, "thickness_mm": p.stack,
               "materials": {m.name: {"youngs_mpa": m.youngs_mpa, "poisson": m.poisson, "density_kg_m3": m.density_kg_m3, "source": m.source} for m in (M19, N45)},
               "bonded_magnet_density_tonne_mm3": main["magnet_density"], "bonded_magnet_density_kg_m3": main["magnet_density"] * 1e12,
               "mesh": {**mesh_counts(main["mesh"]), "msh": "rotor.msh", "levels": {str(h): v for h, v in levels.items()}},
               "hand_estimate": {"note_rotor_with_band_10000": hand_note, **{f"band_free_{rpm}": h for rpm, h in hand.items()}},
               "frames": {rpm: {n: dataclasses.asdict(f) for n, f in slot_frames(p, main["regions"], rpm).items()} for rpm in SPEEDS_RPM},
               "walls": {n: {k: {x: w[x] for x in ("length_mm", "force_n", "direction", "pressure_mpa")} for k, w in walls.items()} for n, walls in loads_at(main, SPEEDS_RPM[0])["loads"].items()},
               "soft_filler": {"youngs_mpa": SOFT_FILLER_YOUNGS_MPA, "share_of_steel_percent": 100 * SOFT_FILLER_YOUNGS_MPA / M19.youngs_mpa, "decides_within": SOFT_DECIDES},
               "runs": {k: {**{x: v for x, v in r.items() if x not in ("stress", "reactions")}, "reactions": r["reactions"], **{**analysis, soft_key: soft_analysis}[k]}
                        for k, r in {**runs, soft_key: soft_run}.items()},
               "pictures": {"cgx": cgx, "composed": picture}}
    results["checks"] = make_checks(p, main, runs, analysis, hand_note, hand, levels, figure_ok, picture, cgx.get("max_von_mises"), soft={"run": soft_run, "analysis": soft_analysis})
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    return results


# ============================================================================================== the verdict
def _percent(a: float, b: float) -> float:
    return 100.0 * (a - b) / b


def soft_position(pressure: float, bonded: float, soft: float) -> dict:
    """Where the soft model's number lies between the bonded model's (0) and the pressure model's (1): position = (soft - bonded) / (pressure - bonded). It decides between the two
    hypotheses (see soft_filler_checks) when it is within SOFT_DECIDES of 1 ("stiffness": the slots' lost hoop stiffness) or of 0 ("load path": the way the load reaches the steel),
    and only when it lies between the two models at all: 0 <= position <= 1, the ends included (T12-P5). A soft number outside that range is not between the models, and no
    hypothesis explains it: a solve that went wrong (a diverged one) looks like that. {"position", "decided", "hypothesis"}; no position (None) and no decision if the pressure and
    the bonded model have the same number."""
    gap = pressure - bonded
    if abs(gap) < 1e-9:
        return {"position": None, "decided": False, "hypothesis": None}
    position = (soft - bonded) / gap
    hypothesis = None
    if 0.0 <= position <= 1.0:
        hypothesis = "stiffness" if position >= 1 - SOFT_DECIDES else "load path" if position <= SOFT_DECIDES else None
    return {"position": position, "decided": hypothesis is not None, "hypothesis": hypothesis}


def soft_filler_checks(pressure: dict, bonded: dict, soft: dict | None) -> list[dict]:
    """The diagnostic's two checks, reported and never gating (the verdict is the pressure model's). `pressure` and `bonded` are the analyses of those models at 10,000 rpm, `soft` is
    {"run", "analysis"} of the soft model at the same speed (None: it was not solved). 1. ccx solves the soft model cleanly: it finishes, its reactions balance its load, and its load
    is the pressure model's (a body's stiffness does not change a load). 2. What raises the pressure model's bore stress above Lame's plain disc: the number against two hypotheses."""
    clean_name = "soft filler (a diagnostic): ccx solves it cleanly, with the same total load as the pressure model"
    cause_name = "soft filler (a diagnostic): what raises the pressure model's bore stress above Lame's disc"
    why = "a diagnostic: the verdict is the pressure model's"
    share = 100 * SOFT_FILLER_YOUNGS_MPA / M19.youngs_mpa
    clean_expected = (f"the bonded deck with the slot bodies' Young's modulus cut to {SOFT_FILLER_YOUNGS_MPA:g} MPa ({share:.2f} % of the steel's: under the 0.1 % that counts as a negligible share "
                      f"of the hoop stiffness) and their mass kept: ccx ends normally, its reactions balance the load within {EQUILIBRIUM_TOLERANCE_PERCENT} %, and its total load and reaction are "
                      f"the pressure model's within {SAME_LOAD_TOLERANCE_PERCENT} % ({OWN_LIMIT})")
    cause_expected = ("two hypotheses for what raises the pressure model's bore stress above Lame's plain disc: (1) its open slots take away the ring's hoop stiffness: then the soft model, whose slot bodies have the "
                      "magnets' mass but no stiffness to speak of, has a bore stress near the pressure model's; (2) the way the load reaches the steel (two pressures on two walls, and not "
                      "through a body tied all along them) raises it: then the soft model's bore is near the bonded model's. The number decides when it lies within "
                      f"{SOFT_DECIDES:.0%} of the gap between the bonded and the pressure model's from one of them ({OWN_LIMIT})")
    source = "the soft deck on the same mesh; the .frd and .dat of ccx"
    if soft is None or not soft["analysis"]["finished"]:
        run = (soft or {}).get("run") or {}
        reason = f"ccx did not finish (exit {run.get('returncode')}, timed out {run.get('timed_out')})" if soft else "the soft model was not solved"
        return [check(clean_name, False, clean_expected, reason, source, gating=False, reason=why), check(cause_name, False, cause_expected, "no number: " + reason, source, gating=False, reason=why)]
    run, a = soft["run"], soft["analysis"]
    same_load = _percent(a["total_load_n"], pressure["total_load_n"])
    same_reaction = _percent(math.hypot(*a["equilibrium"]["reaction_n"]), math.hypot(*pressure["equilibrium"]["reaction_n"]))
    residual = a["equilibrium"]["residual_percent"]
    clean_ok = residual <= EQUILIBRIUM_TOLERANCE_PERCENT and abs(same_load) <= SAME_LOAD_TOLERANCE_PERCENT and abs(same_reaction) <= SAME_LOAD_TOLERANCE_PERCENT
    out = [check(clean_name, clean_ok,
                 clean_expected, f"{run['seconds']:.1f} s (ccx {run['version']}, exit {run['returncode']}); the reactions balance the load to {residual:.5f} %; total load {a['total_load_n']:.3f} N "
                 f"against the pressure model's {pressure['total_load_n']:.3f} N ({same_load:+.5f} %), the reactions {same_reaction:+.5f} %", source, gating=False, reason=why)]
    bore = {k: m["bore"] for k, m in (("pressure", pressure), ("bonded", bonded), ("soft", a))}
    where = soft_position(bore["pressure"]["hoop_mean_mpa"], bore["bonded"]["hoop_mean_mpa"], bore["soft"]["hoop_mean_mpa"])
    bridges = soft_position(pressure["bridge_peak_mpa"], bonded["bridge_peak_mpa"], a["bridge_peak_mpa"])
    versus = {k: _percent(b["hoop_mean_mpa"], b["lame_mpa"]) for k, b in bore.items()}
    inside = where["position"] is not None and 0 <= where["position"] <= 1
    if where["position"] is None:
        verdict = "the pressure and the bonded model have the same bore stress: there is nothing to decide"
    elif not clean_ok:
        verdict = (f"{where['position']:.0%} of the way from the bonded to the pressure model's, but the clean-solve check above did not pass, so this number is not trusted and "
                   f"does not decide")
    elif not inside:
        verdict = (f"{where['position']:.0%} of the way from the bonded to the pressure model's: outside the range between the two models (0% to 100%), so the soft model's bore is "
                   f"not between them and the number does not decide (a solve that went wrong looks like this)")
    elif where["hypothesis"] == "stiffness":
        verdict = (f"{where['position']:.0%} of the way from the bonded to the pressure model's: the open slots' lost hoop stiffness raises the pressure model's bore stress (hypothesis 1), "
                   f"and the way the load reaches the steel barely does")
    elif where["hypothesis"] == "load path":
        verdict = (f"{where['position']:.0%} of the way from the bonded to the pressure model's: the way the magnets' load reaches the steel raises the pressure model's bore stress "
                   f"(hypothesis 2), and the slots' lost hoop stiffness does not")
    else:
        verdict = (f"{where['position']:.0%} of the way from the bonded to the pressure model's, more than {SOFT_DECIDES:.0%} of the gap from both: the number does not decide "
                   f"(the lost stiffness and the load's path both contribute)")
    shown = "" if bridges["position"] is None else (f"; the bridges' zones, largest von Mises: pressure {pressure['bridge_peak_mpa']:.2f}, bonded {bonded['bridge_peak_mpa']:.2f}, soft "
                                                  f"{a['bridge_peak_mpa']:.2f} MPa")
    if bridges["position"] is not None and 0 <= bridges["position"] <= 1:
        shown += (f", so of the pressure model's excess over the bonded model's the soft model puts {bridges['position']:.0%} on the stiffness and "
                  f"{1 - bridges['position']:.0%} on the load's path")
    elif bridges["position"] is not None:
        shown += ", which lies outside the range between the other two, so no share is given"
    out.append(check(cause_name, clean_ok and where["decided"], cause_expected,
                     f"bore mean hoop stress at 10,000 rpm: pressure {bore['pressure']['hoop_mean_mpa']:.2f} MPa ({versus['pressure']:+.2f} % against Lame's {bore['pressure']['lame_mpa']:.2f}), "
                     f"bonded {bore['bonded']['hoop_mean_mpa']:.2f} ({versus['bonded']:+.2f} %), soft {bore['soft']['hoop_mean_mpa']:.2f} ({versus['soft']:+.2f} %): the soft model is " + verdict + shown,
                     source, gating=False, reason=why))
    return out


def make_checks(p: mg.Params, main: dict, runs: dict, analysis: dict, hand_note: dict, hand: dict, levels: dict, figure_ok: bool, picture: dict, cgx_peak: float | None = None,
                soft: dict | None = None) -> list[dict]:
    """Every number of B4 against its expectation, with the source of the expectation. The pressure model carries the verdict; the bonded model is reported beside it
    (not gating), because it is the other bound of the same question: what the choice of the magnets' model changes. `runs` and `analysis` are the six solves that gate; `soft`
    ({"run", "analysis"}: the soft-filler diagnostic) has its own reported checks and is not among them."""
    out: list[dict] = []
    mesh = main["mesh"]
    vis = mg.visible(main["regions"])
    area_steel = sum(t6_area(mesh.nodes, t) for t in mesh.steel)
    area_diff = {"rotor_steel": _percent(area_steel, vis["rotor_steel"].area)}
    for name, rows in mesh.magnets.items():
        area_diff[name] = _percent(sum(t6_area(mesh.nodes, t) for t in rows), vis[name].area)
    worst_area = max(abs(v) for v in area_diff.values())
    out.append(check("mesh: second-order triangles, the areas of motorgeo's rotor, no degenerate element", worst_area <= AREA_TOLERANCE_PERCENT and mesh.quality["min_sicn"] >= MIN_SICN_OK,
                     f"6-node triangles (Gmsh type 9 = CalculiX CPS6) with positive Jacobians; the steel's and each slot's area within {AREA_TOLERANCE_PERCENT} % of motorgeo's "
                     f"and the minimum SICN at least {MIN_SICN_OK} ({OWN_LIMIT}); the brief asks that the model and mesh sizes are stated",
                     f"h_min {mesh.h_min} mm, h_max {mesh.h_max} mm: {len(mesh.nodes):,} nodes (those the elements use, mid-side nodes included; {steel_node_count(mesh):,} of them in the steel, the pressure model's), {len(mesh.steel):,} steel triangles and "
                     f"{sum(len(v) for v in mesh.magnets.values()):,} slot triangles (the bonded model's); areas {', '.join(f'{k} {v:+.4f} %' for k, v in area_diff.items())}; "
                     f"min SICN {mesh.quality['min_sicn']:.3f} at ({mesh.quality['worst']['x_mm']:.2f}, {mesh.quality['worst']['y_mm']:.2f}) mm, mean {mesh.quality['mean_sicn']:.3f}",
                     "motorgeo.visible (shapely) against the 6-point Gauss integral over the curved triangles; Gmsh's minSICN"))
    finished = [r for r in runs.values() if r["finished"]]
    out.append(check("ccx: every solve ends normally", len(finished) == len(runs),
                     "exit code 0, 'Job finished', the .frd file written, inside its timeout, for every solve",
                     f"{len(finished)} of {len(runs)} solves: " + ", ".join(f"{k} {r['seconds']:.1f} s (ccx {r['version']}, exit {r['returncode']})" for k, r in runs.items()),
                     "ccx_static.exe through run_solver"))
    if len(finished) != len(runs):
        return out
    residuals = {k: a["equilibrium"]["residual_percent"] for k, a in analysis.items()}
    out.append(check("equilibrium: the reactions of the two symmetry edges balance the load worked out from the mesh", max(residuals.values()) <= EQUILIBRIUM_TOLERANCE_PERCENT,
                     f"for every solve the constrained reactions (x on the y-axis edge, y on the x-axis edge) equal minus the body force of the elements plus the applied forces, within "
                     f"{EQUILIBRIUM_TOLERANCE_PERCENT} % ({OWN_LIMIT})",
                     "residual of the load vector: " + ", ".join(f"{k} {v:.5f} %" for k, v in residuals.items()),
                     "the .dat file's reactions against rho w^2 t x (integral of r dA) over the curved elements, and the applied nodal forces"))
    same = {rpm: _percent(analysis[f'bonded_{rpm}_h{H_MIN_MAIN:g}']["total_load_n"], analysis[f'pressure_{rpm}_h{H_MIN_MAIN:g}']["total_load_n"]) for rpm in SPEEDS_RPM}
    react = {rpm: _percent(math.hypot(*analysis[f'bonded_{rpm}_h{H_MIN_MAIN:g}']['equilibrium']['reaction_n']), math.hypot(*analysis[f'pressure_{rpm}_h{H_MIN_MAIN:g}']['equilibrium']['reaction_n']))
             for rpm in SPEEDS_RPM}
    out.append(check("the two magnet models carry the same total load", max(abs(v) for v in list(same.values()) + list(react.values())) <= SAME_LOAD_TOLERANCE_PERCENT,
                     f"the bonded body (the real magnet's mass spread over the slot) and the pressures on the walls are the same force: the totals, and the reactions ccx gives, agree "
                     f"within {SAME_LOAD_TOLERANCE_PERCENT} % ({OWN_LIMIT}), so the two models differ only in how the load reaches the steel",
                     ", ".join(f"{rpm:,} rpm: total {analysis[f'pressure_{rpm}_h{H_MIN_MAIN:g}']['total_load_n']:.3f} N (pressure) against {analysis[f'bonded_{rpm}_h{H_MIN_MAIN:g}']['total_load_n']:.3f} N "
                               f"(bonded), {same[rpm]:+.5f} %; the reactions {react[rpm]:+.5f} %" for rpm in SPEEDS_RPM),
                     "the loads of the two decks, and ccx's reactions"))
    for variant in VARIANTS:
        low, high = analysis[f"{variant}_10000_h{H_MIN_MAIN:g}"], analysis[f"{variant}_12000_h{H_MIN_MAIN:g}"]
        extra = {} if variant == VERDICT_VARIANT else {"gating": False, "reason": f"the verdict is the {VERDICT_VARIANT} model's"}
        out.append(scale_check(variant, low["peak"]["value_mpa"], high["peak"]["value_mpa"], **extra))
        if variant == VERDICT_VARIANT:
            out.append(scale_check(f"{variant} (the stress of the bridges' zones)", low["bridge_peak_mpa"], high["bridge_peak_mpa"]))
    verdict = analysis[f"{VERDICT_VARIANT}_10000_h{H_MIN_MAIN:g}"]
    place = verdict["peak"]["place"]
    out.append(check(f"{VERDICT_VARIANT}: the peak von Mises sits in a bridge or a rib, and its place is reported",
                     place["in_a_strip"] and place["feature"] in ("outer bridge", "centre rib"),
                     f"the largest von Mises stress of the steel is within {STRIP_ZONE_MM} mm ({OWN_LIMIT}) of the cut across a bridge or the centre rib, and its place is named",
                     f"{verdict['peak']['value_mpa']:.2f} MPa at ({place['x_mm']:.3f}, {place['y_mm']:.3f}) mm, r {place['r_mm']:.2f} mm ({p.r_ro - place['r_mm']:.2f} mm from the rim), "
                     f"{place['angle_deg']:.2f} deg: the {place['feature']} of {place['slot']} (pole {place['pole']}), {place['distance_to_strip_mm']:.2f} mm from its cut; "
                     f"{place['slot_detail']}; the largest value in each slot's bridge: "
                     + ", ".join(f"{k} {v['value']:.2f}" for k, v in verdict['zone_peaks'].items() if k.startswith('magnet_')) + f" MPa; the rib: "
                     + ", ".join(f"{k} {v['value']:.2f}" for k, v in verdict['zone_peaks'].items() if k.startswith('rib')) + " MPa",
                     "the brief; ccx's nodal von Mises stress"))
    bonded = analysis[f"bonded_10000_h{H_MIN_MAIN:g}"]
    magnet_bodies = body_peak(mesh, runs[f"bonded_10000_h{H_MIN_MAIN:g}"]["stress"])            # over every node of the slot bodies: the interior ones and the copies of the wall nodes
    out.append(check("bonded: the choice of the magnets' model changes the stress in the bridges", bonded["bridge_peak_mpa"] < verdict["bridge_peak_mpa"],
                     "the bonded body fills the slot and takes part of the load, so the stress in the bridges is below the pressure model's, which is the upper bound; the verdict is "
                     "the upper bound's (a magnet is brittle and its glue is soft, the note 4.2 and 6.1)",
                     f"at 10,000 rpm the largest stress of the steel in the bridges' zones: {verdict['bridge_peak_mpa']:.2f} MPa (pressure), {bonded['bridge_peak_mpa']:.2f} MPa (bonded), "
                     f"a factor {verdict['bridge_peak_mpa'] / bonded['bridge_peak_mpa']:.1f}; the peak of the whole steel: {verdict['peak']['value_mpa']:.2f} MPa in the "
                     f"{place['feature']} (pressure), {bonded['peak']['value_mpa']:.2f} MPa in the {bonded['peak']['place']['feature']} at ({bonded['peak']['place']['x_mm']:.2f}, "
                     f"{bonded['peak']['place']['y_mm']:.2f}) mm (bonded); the largest von Mises stress at any node of the bonded magnet bodies is {magnet_bodies:.2f} MPa (N45's tensile strength is 78 MPa)",
                     "the two decks, the same mesh and the same total load", gating=False, reason=f"the verdict is the {VERDICT_VARIANT} model's: this is the other bound"))
    # the hand estimate
    out.append(compare("hand estimate on the note's own rotor (with the air band): the bridges' nominal stress at 10,000 rpm", hand_note["nominal_mpa"], NOTE_HAND_NOMINAL_MPA,
                       NOTE_TOLERANCE_PERCENT, f"04 Structural 6.1 and 6.4: 15 MPa (the note prints two digits, so its 15 is anywhere in 14.5-15.5, 3.3 % either way; the bench allows {NOTE_TOLERANCE_PERCENT} %, "
                       f"that rounded up ({OWN_LIMIT})); "
                       f"the function that is applied to the band-free rotor below, applied first to the rotor the note worked on", "MPa"))
    sections = verdict["sections"]
    hand10 = hand[10_000]
    fe_pull = sum(s["axial_pull_n"] for s in sections) / len(sections)
    diffs = [_percent(s["axial_pull_n"], s["hand_pull_n"]) for s in sections]
    refined = [_percent(s["axial_pull_n"], s["refined_pull_n"]) for s in sections]
    widths = hand10["strips_mm"]["total"]
    nominal_fe = fe_pull / p.stack / widths
    hand_axis = sum(s["hand_pull_n"] for s in sections) / len(sections)                  # the hand pull along the pole's axis: the component the strips balance
    nominal_hand = hand_axis / p.stack / widths
    out.append(check("the force the three strips of a pole carry agrees with the hand estimate of the same method", all(abs(d) <= HAND_TOLERANCE_PERCENT for d in diffs),
                     f"the force through the rib and the two bridges (the stress integrated over the three cuts) is the pull of the pole's cap and its magnets, as the note works it out "
                     f"(w^2 sum of rho A r over the cap's steel and the two magnets): within {HAND_TOLERANCE_PERCENT} % ({OWN_LIMIT}: the note's method puts all of the magnets' pull on the "
                     f"cap, and about a third of it bears on the slots' end walls, which are not the cap's)",
                     f"pole 0 and 1: the strips carry {sections[0]['axial_pull_n']:.1f} N and {sections[1]['axial_pull_n']:.1f} N over the {p.stack:g} mm stack against {sections[0]['hand_pull_n']:.1f} N "
                     f"and {sections[1]['hand_pull_n']:.1f} N ({diffs[0]:+.2f} %, {diffs[1]:+.2f} %), the hand pull along the pole's axis (the note's formula with the magnets' forces taken "
                     f"along it; the note's own sum of the magnitudes is {hand10['pull_n_over_stack']:.1f} N, {_percent(hand10['pull_n_over_stack'], hand_axis):.2f} % more); "
                     f"the transverse force, zero by symmetry: {sections[0]['transverse_n']:.2f} N and {sections[1]['transverse_n']:.2f} N; "
                     f"nominal stress F/({widths:.2f} mm of strips), both on the axis basis: {nominal_fe:.2f} MPa (finite elements) against {nominal_hand:.2f} MPa (hand), so the peak/nominal "
                     f"factor is {verdict['peak']['value_mpa'] / nominal_fe:.2f} (finite elements) or {verdict['peak']['value_mpa'] / nominal_hand:.2f} (hand); the share of the three strips "
                     f"(rib, bridge, bridge, % of the total): " + ", ".join(f"{100 * s['axial_pull_n'] / sections[0]['axial_pull_n']:.0f}" for s in sections[0]['strips'].values()),
                     "04 Structural 6.1, the estimate applied to the band-free rotor (strip_geometry, hand_estimate); the .frd's stress integrated over the cuts"))
    out.append(check("... and with the magnets' end-wall load left off the cap (the model's own free body)", all(abs(d) <= REFINED_TOLERANCE_PERCENT for d in refined),
                     f"the same force against the hand estimate whose magnets' load is only what the walls that belong to the cap carry (the long wall; the end wall is on the far side "
                     f"of the bridge): within {REFINED_TOLERANCE_PERCENT} % ({OWN_LIMIT}: the error of integrating a smoothed nodal field over three cuts)",
                     f"{', '.join(f'pole {s['pole']}: {s['axial_pull_n']:.1f} N against {s['refined_pull_n']:.1f} N ({d:+.2f} %)' for s, d in zip(sections, refined))}",
                     "the free body of the cap: its own centrifugal force and the pressure on the long walls"))
    out.append(compare("the note's 15 MPa against this rotor's nominal stress (finite elements)", nominal_fe, NOTE_HAND_NOMINAL_MPA, NOTE_TOLERANCE_PERCENT,
                       f"04 Structural 6.1: 15 MPa at 10,000 rpm, a number for the rotor WITH its air band (its strips are 1.5 + 1.4 + 1.4 = 4.3 mm and its cap 162 mm2); "
                       f"without the band the strips are {widths:.2f} mm and the cap {hand10['cap_area_mm2']:.0f} mm2, and the note's formula (the magnitudes of the forces added) gives "
                       f"{hand10['nominal_mpa']:.2f} MPa", "MPa",
                       gating=False, reason="the note's 15 MPa is the rotor with the band: the like-for-like comparison is the check above"))
    # the bore: Lame's plain disc is the answer of the model whose slots are filled (the bonded one); the pressure model's slots are voids
    bore_bonded, bore_pressure = bonded["bore"], verdict["bore"]
    out.append(compare("bonded: the mean hoop stress at the shaft's bore against Lame's plain disc (a gate on the units, the speed and the load)", bore_bonded["hoop_mean_mpa"],
                       bore_bonded["lame_mpa"], LAME_TOLERANCE_PERCENT,
                       f"the note 5 and Task 4's ring: the hoop stress at the bore of a plain steel disc (radius 7.5 to 45 mm, 7,650 kg/m3, nu 0.3), 14.10 MPa at 10,000 rpm. The bonded model is "
                       f"not that disc (its slot bodies are lighter, 6,627 kg/m3, and softer, 160 GPa, and the ring is broken at the slots): by mass alone it is about 0.6 % from it (a rough 1D "
                       f"axisymmetric estimate in the task review, 2026-11-06: 14.01 MPa), so this check gates the units, the speed and the load and is not a verification to 0.1 %. The limit, "
                       f"{LAME_TOLERANCE_PERCENT} %, was chosen after the first run, which gave -0.11 % ({OWN_LIMIT}); before it the check was set at 5 % for the pressure model, whose slots are "
                       f"voids and which misses it (reported below): it was re-aimed at the bonded model and the limit was not widened", "MPa"))
    out.append(compare("pressure: the mean hoop stress at the shaft's bore against Lame's plain disc", bore_pressure["hoop_mean_mpa"], bore_pressure["lame_mpa"], LAME_TOLERANCE_PERCENT,
                       f"the same plain disc; the pressure model's slots are voids, so it is not that disc and its bore carries more of the load than Lame's disc would: "
                       f"{_percent(bore_pressure['hoop_mean_mpa'], bore_pressure['lame_mpa']):+.2f} % (what raises it is the question of the soft-filler diagnostic below)", "MPa",
                       gating=False, reason="a model with voids is not Lame's plain disc: the bonded model carries this check"))
    out += soft_filler_checks(verdict, bonded, soft)
    if cgx_peak is not None:
        out.append(compare("the peak von Mises stress of the steel against cgx's own", verdict["peak"]["value_mpa"], cgx_peak, CGX_TOLERANCE_PERCENT,
                           f"cgx computes the von Mises stress of the .frd's components with its own code and prints the maximum of the faces it draws; the two agree to {CGX_TOLERANCE_PERCENT} % "
                           f"({OWN_LIMIT}: they read the same numbers)", "MPa"))
    else:
        out.append(skipped("the peak von Mises stress of the steel against cgx's own", f"agreement within {CGX_TOLERANCE_PERCENT} %", "cgx printed no maximum (its frames were not made)",
                           "cgx_STATIC.exe"))
    # the mesh
    peaks = {h: analysis[v["key"]]["peak"]["value_mpa"] for h, v in levels.items() if analysis[v["key"]]["finished"]}
    finest, main_h = min(peaks), H_MIN_MAIN
    drift = _percent(peaks[main_h], peaks[finest])
    out.append(check(f"{VERDICT_VARIANT}: the peak does not depend on the mesh", abs(drift) <= MESH_TOLERANCE_PERCENT and len(peaks) == len(levels),
                     f"the peak on the main mesh (h_min {main_h} mm) within {MESH_TOLERANCE_PERCENT} % of the peak on the finest ({OWN_LIMIT}: chosen after the exploratory runs and before this one; "
                     f"the fillets are 0.3 mm, so the peak is finite and converges)",
                     ", ".join(f"h_min {h} mm: {levels[h]['nodes']:,} nodes, {peaks[h]:.2f} MPa" for h in sorted(peaks, reverse=True)) + f"; main against finest {drift:+.2f} %",
                     "the same model on three meshes"))
    spread = verdict["bridge_spread_percent"]
    out.append(check(f"{VERDICT_VARIANT}: the four slots' bridges have the same peak (they are mirror images)", spread <= SYMMETRY_TOLERANCE_PERCENT,
                     f"the largest stress in the bridge zones of the four slots (poles 0 and 1, upper and lower slot: mirror images of each other, meshed differently) spread by at most "
                     f"{SYMMETRY_TOLERANCE_PERCENT} % of their mean ({OWN_LIMIT})",
                     ", ".join(f"{k} {v['value']:.2f}" for k, v in verdict['zone_peaks'].items() if k.startswith('magnet_')) + f" MPa: spread {spread:.2f} %",
                     "the symmetry of the rotor: a spread is the mesh's"))
    out.append(check("the picture: figure.png is saved", figure_ok, "out\\bench\\B4\\figure.png exists and is not empty", f"saved: {figure_ok}", "bench/b4_ccx.py"))
    out.append(check("the picture is cgx's own render of the von Mises stress", bool(picture.get("made")),
                     "figure.png is cgx's hardcopy of the verdict run (the frames cgx writes, put together by Python after checking that their colours are the legend's); "
                     "when cgx's frames are not made the picture is plotted from the .frd in Python, and this says so",
                     f"cgx's frames made and composed: {bool(picture.get('made'))}; {picture.get('colours_match', 0):.0%} of the model frame's colours are a legend colour times cgx's lighting "
                     f"factor {CGX_LIGHT}", "cgx_STATIC.exe -b (read, ds, plot fv, hcpy tga)", gating=False, reason="the picture is saved either way; this only says whose it is"))
    return out
