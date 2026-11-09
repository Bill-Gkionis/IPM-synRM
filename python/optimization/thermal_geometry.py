"""optimization.thermal_geometry: one slot pitch of the stator drawn as heat-conduction regions (pure Python, shapely; no solver).

What it is, and why
    The thermal network (thermal.py) is checked against a 2D finite-element model of a slice of the stator. This file draws that slice once, so that
    FEMM's heat-flow solver and CalculiX's heat transfer solve the SAME regions. It builds, in millimetres, the regions of one slot pitch (10 degrees:
    the slot in the middle, half a tooth on each side; the two cut lines run through the tooth centres, where a symmetric problem has no heat flow
    across, which is what an untouched boundary means in both programs):
        winding   the block of copper, enamel and varnish inside the liner, one material of the equivalent conductivity (thermal.py, step 3)
        liner     the slot paper along the two flanks and the round bottom (not the top: the slot's top wall faces the opening)
        teeth     the stator steel from the bore down to just below the slot bottoms
        yoke      the ring behind the slots, from r_yb + margin to the outside radius (the margin keeps the two iron regions apart where the round
                  slot bottom touches r_yb)
        contact   a thin layer on the stator's outside radius: its conductivity is thickness / (contact resistance per unit area)
        housing   the aluminium wall, to the outer radius, where the cooling's heat transfer coefficient acts
    The slot opening (the air funnel from the bore to the slot's top) is a hole: its walls are untouched, so no heat crosses them.

Method
    motorgeo draws the slot (motorgeo._slot): a top chord S1-S2 at r_shoe, two radial flanks at +-slot_half_deg down to r_wall, then a round bottom
    (a 185 degree arc tangent to both flanks, centre at r_wall / cos(a_s) on the slot's axis, radius r_wall tan(a_s)). The same outline is built here
    in polar coordinates. The liner is the band of thickness t along the flanks and the bottom: shapely's negative buffer of the slot polygon (extended
    upward first, so that the top edge does not recede), clipped to the slot. Everything else is boolean operations on polygons whose arcs are
    polylines of 0.5 degree segments (a chord error of 4e-5 of the radius: 0.003 mm at the outer radius, and none of it matters here).

Run
    uv run python -m optimization.thermal_geometry        draws the slice, prints the areas against motorgeo's, and saves a picture

Numbers to expect (the optimized design, liner 0.25 mm; areas of one slot pitch): slot 72.436 mm^2 (motorgeo: 2 x 36.218 = 72.436), the winding block 65.079 and the liner 7.357
(thermal.py's closed form t x perimeter says 7.455: the liner's two top corners are not covered), the teeth 55.827 (in three blocks: shoe 5.64, fin 32.56, root 17.62), the yoke
43.920, the contact layer 1.118, the housing 56.914; the steel 99.747 mm^2 against motorgeo's 99.749 (3,591 mm^2 / 36). The regions overlap by less than 1e-13 mm^2 and leave
nothing uncovered but the opening (2.204 mm^2). The yoke starts 0.1 mm outside r_yb (4.16 mm thick at its thinnest), so that the round slot bottom, which touches r_yb, does
not pinch the two iron regions into a cusp the mesher cannot resolve. The teeth region is one connected piece (the half teeth join under the slot bottom); a boolean operation
leaves two zero-area slivers where the slot's arc touches the yoke's inner radius, and _drop_slivers removes pieces under 0.001 mm^2.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from shapely.geometry import Polygon
from shapely.ops import unary_union

import motorgeo as mg

ARC_STEP_DEG = 0.5


def pol(r: float, deg: float) -> tuple[float, float]:
    a = math.radians(deg)
    return r * math.cos(a), r * math.sin(a)


def arc_points(r: float, deg0: float, deg1: float, step: float = ARC_STEP_DEG, cx: float = 0.0, cy: float = 0.0) -> list[tuple[float, float]]:
    """Points along an arc of radius r about (cx, cy), from angle deg0 to deg1 (either direction), the end points included."""
    n = max(2, int(math.ceil(abs(deg1 - deg0) / step)) + 1)
    return [(cx + r * math.cos(math.radians(a)), cy + r * math.sin(math.radians(a))) for a in np.linspace(deg0, deg1, n)]


def ring_sector(r0: float, r1: float, deg0: float, deg1: float) -> Polygon:
    """The annular sector between radii r0 < r1 and angles deg0 < deg1."""
    outer = arc_points(r1, deg0, deg1)
    inner = arc_points(r0, deg1, deg0)
    return Polygon(outer + inner)


def _drop_slivers(geom, min_area_mm2: float = 1e-3):
    """A boolean operation can leave zero-area slivers where two curves touch (the round slot bottom touches the yoke's inner radius): drop the pieces
    smaller than min_area_mm2, which no mesh could resolve and no heat could notice."""
    parts = [g for g in getattr(geom, "geoms", [geom]) if g.area > min_area_mm2]
    return parts[0] if len(parts) == 1 else unary_union(parts)


@dataclass
class Slice:
    """The regions of one slot pitch [mm], as shapely polygons, and the numbers that define them."""
    regions: dict
    pitch_deg: float
    r_si: float
    r_so: float
    r_ho: float
    r_yb: float
    t_liner: float
    t_contact: float
    slot: Polygon
    opening: Polygon
    boundary_outer: list                # vertices of the outer arc, from -pitch/2 to +pitch/2


def build_slice(p: mg.Params, t_liner_mm: float = 0.25, t_contact_mm: float = 0.1, r_housing_out_mm: float = 69.0, yoke_margin_mm: float = 0.1,
                split_teeth: bool = False) -> Slice:
    """One slot pitch of the stator, centred on the slot (angles -5 ... +5 degrees for 36 slots), as regions in millimetres."""
    pitch = 360.0 / p.slots
    half = pitch / 2
    a_s = math.radians(p.slot_half_deg)
    sd = p.slot_half_deg
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    r_b = p.r_wall * math.tan(a_s)
    c_b = p.r_wall / math.cos(a_s)
    r_yb = c_b + r_b
    s1, s2 = pol(p.r_shoe, -sd), pol(p.r_shoe, sd)
    w1 = pol(p.r_wall, -sd)
    w2 = pol(p.r_wall, sd)
    arc = arc_points(r_b, -(90.0 + sd), 90.0 + sd, cx=c_b)             # from the tangent point on the -sd flank, round the bottom, to the +sd flank
    slot = Polygon([s1, w1] + arc[1:-1] + [w2, s2])
    o1, o2, o3, o4 = pol(p.r_si, -a_o), pol(p.r_si + p.h_open, -a_o), pol(p.r_si + p.h_open, a_o), pol(p.r_si, a_o)
    opening = Polygon([o1, o2, s1, s2, o3, o4] + arc_points(p.r_si, a_o, -a_o)[1:-1])
    # the liner: the offset band along the flanks and the bottom. Extend the slot upward by 3 mm first so that the offset keeps the top edge where it is.
    up = (-3.0, 0.0)
    extension = Polygon([s1, s2, (s2[0] + up[0], s2[1]), (s1[0] + up[0], s1[1])])
    inner = unary_union([slot, extension]).buffer(-t_liner_mm, join_style=2).intersection(slot)
    liner = slot.difference(inner)
    sector = ring_sector(p.r_si, p.r_so, -half, half)
    steel = sector.difference(slot).difference(opening)
    yoke = steel.intersection(ring_sector(r_yb + yoke_margin_mm, p.r_so, -half, half))
    teeth = steel.difference(yoke)
    contact = ring_sector(p.r_so, p.r_so + t_contact_mm, -half, half)
    housing = ring_sector(p.r_so + t_contact_mm, r_housing_out_mm, -half, half)
    teeth = _drop_slivers(teeth)
    regions = {"winding": inner, "liner": liner, "teeth": teeth, "yoke": _drop_slivers(yoke), "contact": contact, "housing": housing}
    if split_teeth:
        # the same steel in three blocks: the shoe (below the slot's top, r < r_shoe), the fin alongside the flank (r_shoe ... r_wall: the part the network's
        # tooth node stands for) and the root (r > r_wall, where the tooth widens into the yoke). Same material and heat source in all three.
        del regions["teeth"]
        regions["teeth_shoe"] = _drop_slivers(teeth.intersection(ring_sector(p.r_si - 1.0, p.r_shoe, -half, half)))
        regions["teeth_fin"] = _drop_slivers(teeth.intersection(ring_sector(p.r_shoe, p.r_wall, -half, half)))
        regions["teeth_root"] = _drop_slivers(teeth.intersection(ring_sector(p.r_wall, r_yb + yoke_margin_mm + 1.0, -half, half)))
    return Slice(regions=regions, pitch_deg=pitch, r_si=p.r_si, r_so=p.r_so, r_ho=r_housing_out_mm, r_yb=r_yb, t_liner=t_liner_mm,
                 t_contact=t_contact_mm, slot=slot, opening=opening, boundary_outer=arc_points(r_housing_out_mm, -half, half))


def check_slice(sl: Slice, p: mg.Params) -> dict:
    """Areas against motorgeo's, and that the regions tile the slice without overlap or gap."""
    a = {k: v.area for k, v in sl.regions.items()}
    if "teeth_fin" in a:
        a["teeth"] = a["teeth_shoe"] + a["teeth_fin"] + a["teeth_root"]
    total = sum(a.values())
    union = unary_union(list(sl.regions.values()))
    sector = ring_sector(p.r_si, sl.r_ho, -sl.pitch_deg / 2, sl.pitch_deg / 2)
    overlap = sum(v.area for v in sl.regions.values()) - union.area
    gap = sector.difference(union).area - sl.opening.area
    regions = mg.build(p)
    metrics, _ = mg.check(p, regions)
    return {"areas_mm2": a, "slot_mm2": sl.slot.area, "slot_motorgeo_mm2": 2 * metrics["coil_layer_area_mm2"],
            "steel_slice_mm2": a["teeth"] + a["yoke"], "steel_motorgeo_mm2_per_pitch": p.n_sectors * mg.visible(regions)["stator_steel"].area / p.slots,
            "overlap_mm2": overlap, "uncovered_mm2": gap, "opening_mm2": sl.opening.area}


def main() -> int:
    import json
    from optimization import thermal
    ov = thermal.final_overrides()
    p = thermal.design_params(ov)
    sl = build_slice(p)
    info = check_slice(sl, p)
    print(json.dumps(info, indent=1))
    g = thermal.build_geom(p)
    print(f"thermal.py: slot {g.a_slot * 1e6:.3f} mm2, liner closed form {0.25e-3 * g.perim * 1e6:.3f} mm2, yoke ring {g.a_yoke * 1e6 / p.slots:.3f} mm2 a pitch")
    ok = abs(info["slot_mm2"] - info["slot_motorgeo_mm2"]) < 0.02 and abs(info["overlap_mm2"]) < 1e-6 and abs(info["uncovered_mm2"]) < 1e-3
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 6))
    colours = {"winding": "#eb6834", "liner": "#898781", "teeth": "#2a78d6", "yoke": "#1baf7a", "contact": "#0b0b0b", "housing": "#c3c2b7"}
    for name, poly in sl.regions.items():
        for part in getattr(poly, "geoms", [poly]):
            x, y = part.exterior.xy
            ax.fill(x, y, color=colours[name], alpha=0.8, label=name, linewidth=0.2)
    ax.set_aspect("equal")
    ax.set_xlim(44, 70)
    ax.set_ylim(-6.5, 6.5)
    handles, labels = ax.get_legend_handles_labels()
    uniq = dict(zip(labels, handles))
    ax.legend(uniq.values(), uniq.keys(), fontsize=8)
    out = thermal.paths.path("out") / "optimization" / "thermal"
    out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / "slice_geometry.png", dpi=110)
    print("saved", out / "slice_geometry.png", "OK" if ok else "CHECK")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
