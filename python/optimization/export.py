"""export: a candidate of the study as a DXF of the whole motor, clean enough to open in CAD, checked by reading it back.

What it is, and why not motorgeo.to_dxf
    The study works on one 90-degree sector (two poles, nine slots), because the motor repeats four times. CAD wants the whole
    lamination. motorgeo.to_dxf(full=True) turns the sector four times, but it draws every edge of FEMM's model: the sector's
    radial seams through the steel, and the two construction arcs inside the air gap that carry FEMM's air-gap element. Neither is
    a real edge of the motor. This writer draws only what the motor has:
      ROTOR       the rotor lamination: the shaft bore and the outer circle, and the 16 magnet slots (lines and corner fillets)
      MAGNETS     the 16 magnets themselves, the bought size (motorgeo's mag_w x mag_t), centred in their slots
      MAGNET_DIR  a short arrow across each magnet: its magnetisation, north poles outward
      STATOR      the stator lamination: the outer circle and the bore, the 36 tooth faces and the 36 slot outlines
      COILS       the 72 coil layers (two per slot), each a closed outline
      LABELS      the phase of each coil layer (A+, C-, ...) as text
    Units: millimetres. The motor's centre is the origin; the d-axis of the first pole is at 22.5 degrees, as in motorgeo.

How each layer is made (from motorgeo's regions, so the drawing is the model)
    rotor slots   every edge of the regions magnet_<k>_<up|dn>, turned through 0, 90, 180 and 270 degrees
    magnets       the rectangle of +-mag_w/2 along the magnet and +-mag_t/2 across it around the slot's centre, mirrored and placed
                  as motorgeo._pole_slot places the slot
    stator        for each slot, the edges of its opening and its two coil layers that belong to only ONE of the three regions
                  (the two edges they share are inside the slot), less the mouth arc on the bore; and between two slots, the
                  tooth face: the bore arc from one opening's edge to the next
    circles       r_shaft and r_ro (ROTOR), r_so (STATOR), drawn as full circles
The check after writing (compare): the file is read back with ezdxf, each layer's lines and arcs are polygonized (shapely), and the
faces' areas are held against the model's: the rotor lamination against 4 x the visible area of rotor_steel (the ring less the slots),
the stator lamination against 4 x the visible area of stator_steel, every coil face against its region, every magnet against
mag_w x mag_t. The counts are 16 slots, 16 magnets, 72 coil faces. A preview PNG is drawn from the file read back, not from the model.

Run (from the Python to FEMM folder; no FEMM)
    uv run python -m optimization.export --overrides "{}" --out out\\optimization\\export\\reference.dxf
Expect on the reference: every area within 0.05 % of the model's (the read-back flattens arcs to 1 micrometre of sagitta), the counts
16, 16 and 72, and a preview next to the file.

Self-check (the answers follow).
  1. Why are the edges a slot's opening shares with its first coil layer left out of the STATOR layer?
  2. Why is the check made on the file read back, and not on the edges before writing?
Answers.
  1. That edge is the line between the opening's air and the copper: inside the slot, not on the steel. The steel's outline is the
     slot's outer boundary, which is made of the edges that belong to only one of the slot's three regions.
  2. Writing is where things go wrong unseen: an arc written clockwise instead of counter-clockwise, a layer name misspelt, a gap of a
     nanometre that CAD will not close. The file read back is what CAD will see.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from shapely.geometry import LineString, Point
from shapely.ops import polygonize, unary_union

import motorgeo as mg
from optimization import design

LAYERS = {"ROTOR": 8, "MAGNETS": 1, "MAGNET_DIR": 1, "STATOR": 8, "COILS": 30, "LABELS": 7}   # name: AutoCAD colour index
ARROW_MM = 1.2


def _ekey(e) -> tuple:
    """One edge's identity, whatever its direction (as motorgeo.unique_edges)."""
    return (frozenset((mg.key(e.a), mg.key(e.b))), mg.key(e.c) if isinstance(e, mg.Arc) else None)


def _magnet_rects(p: mg.Params) -> list[tuple[list, float, tuple]]:
    """The 16 magnets: (four corners, magnetisation angle [deg], centre), placed as motorgeo._pole_slot places their slots."""
    b = math.radians(p.v_half_deg)
    u, n = (math.cos(b), math.sin(b)), (math.sin(b), -math.cos(b))
    at = lambda s, t: (p.mag_x + s * u[0] + t * n[0], p.mag_y + s * u[1] + t * n[1])
    hw, ht = p.mag_w / 2, p.mag_t / 2
    out = []
    pole_pitch = 180 / p.pole_pairs
    for k in range(2 * p.pole_pairs):
        theta_d, north = pole_pitch * (k + 0.5), k % 2 == 0
        for side in (1, -1):
            g = lambda q: mg.rot((q[0], side * q[1]), theta_d)
            corners = [g(at(s, t)) for s, t in ((-hw, ht), (hw, ht), (hw, -ht), (-hw, -ht))]
            magdir = (theta_d + side * (p.v_half_deg - 90) + (0 if north else 180)) % 360
            out.append((corners, magdir, g(at(0, 0))))
    return out


def edges_by_layer(p: mg.Params) -> tuple[dict, list]:
    """{layer: [motorgeo Line/Arc, ...]} for the whole motor (no circles), and the [(text, point)] labels."""
    regions = {r.name: r for r in mg.build(p)}
    layers: dict[str, dict] = {name: {} for name in LAYERS}
    labels = []

    def add(layer, e):
        layers[layer].setdefault(_ekey(e), e)

    slots_per_sector = round(p.slots / p.n_sectors)
    for s in range(p.n_sectors):
        deg = s * p.sector_deg
        for r in regions.values():
            if r.name.startswith("magnet_"):
                for e in r.edges:
                    add("ROTOR", e.moved(deg))
        for j in range(slots_per_sector):
            parts = [regions[f"opening_{j}"], regions[f"coil_{j}_0"], regions[f"coil_{j}_1"]]
            count: dict = {}
            for r in parts:
                for e in r.edges:
                    count.setdefault(_ekey(e), [0, e])[0] += 1
            for k, (c, e) in count.items():
                mouth = isinstance(e, mg.Arc) and abs(e.radius - p.r_si) < 1e-6
                if c == 1 and not mouth:
                    add("STATOR", e.moved(deg))
            for r in parts[1:]:
                for e in r.edges:
                    add("COILS", e.moved(deg))
                sign = "+" if r.turns > 0 else "-"
                labels.append((f"{r.circuit}{sign}", mg.rot(r.label, deg)))
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    phis = [(k + 0.5) * 360 / p.slots for k in range(p.slots)]
    for k, phi in enumerate(phis):
        nxt = phis[(k + 1) % p.slots] + (360 if k + 1 == p.slots else 0)
        add("STATOR", mg.Arc(mg.pol(p.r_si, phi + a_o), mg.pol(p.r_si, nxt - a_o), True))
    for corners, _, _ in _magnet_rects(p):
        for a, b in zip(corners, corners[1:] + corners[:1]):
            add("MAGNETS", mg.Line(a, b))
    return {name: list(d.values()) for name, d in layers.items()}, labels


def write(overrides: dict, path: Path) -> mg.Params:
    """Write the whole motor of this candidate (Params overrides, or optimizer variables such as bridge_mm) to `path` (DXF R2010)."""
    import ezdxf
    from ezdxf.enums import TextEntityAlignment
    p = design.params(design.decode(overrides))
    by_layer, labels = edges_by_layer(p)
    doc = ezdxf.new("R2010", units=ezdxf.units.MM)
    msp = doc.modelspace()
    for name, colour in LAYERS.items():
        doc.layers.add(name, color=colour)
    for layer, edges in by_layer.items():
        for e in edges:
            if isinstance(e, mg.Line):
                msp.add_line(e.a, e.b, dxfattribs={"layer": layer})
            else:
                e = e.canonical()
                a0, a1 = e.angles()
                msp.add_arc(e.c, e.radius, a0, a1, dxfattribs={"layer": layer})
    for r, layer in ((p.r_shaft, "ROTOR"), (p.r_ro, "ROTOR"), (p.r_so, "STATOR")):
        msp.add_circle((0.0, 0.0), r, dxfattribs={"layer": layer})
    for _, magdir, centre in _magnet_rects(p):
        d = mg.pol(1, magdir)
        tail, tip = mg.lin(centre, d, -ARROW_MM / 2), mg.lin(centre, d, ARROW_MM / 2)
        msp.add_line(tail, tip, dxfattribs={"layer": "MAGNET_DIR"})
        for turn in (150, -150):
            msp.add_line(tip, mg.lin(tip, mg.pol(1, magdir + turn), 0.35), dxfattribs={"layer": "MAGNET_DIR"})
    for text, pt in labels:
        msp.add_text(text, dxfattribs={"layer": "LABELS", "height": 0.9}).set_placement(pt, align=TextEntityAlignment.MIDDLE_CENTER)
    path.parent.mkdir(parents=True, exist_ok=True)
    doc.saveas(path)
    return p


def _faces(path: Path, layer: str, sagitta: float = 0.001) -> list:
    """The faces that a layer's lines, arcs and circles of the file at `path` enclose (shapely polygonize)."""
    import ezdxf
    lines = []
    for e in ezdxf.readfile(path).modelspace().query(f'*[layer=="{layer}"]'):
        if e.dxftype() == "LINE":
            lines.append(LineString([e.dxf.start.xy, e.dxf.end.xy]))
        elif e.dxftype() in ("ARC", "CIRCLE"):
            pts = [v.xy for v in e.flattening(sagitta)]
            if e.dxftype() == "CIRCLE":
                pts.append(pts[0])
            lines.append(LineString(pts))
    snap = lambda ls: LineString([mg.key(xy) for xy in ls.coords])
    return list(polygonize(unary_union([snap(ls) for ls in lines])))


def compare(p: mg.Params, path: Path) -> dict:
    """Read the file back and hold its faces against the model (the module docstring). Returns {"checks": [...], "ok": bool}."""
    regions = mg.build(p)
    vis = mg.visible(regions)
    rotor, stator, coils, magnets = (_faces(path, layer) for layer in ("ROTOR", "STATOR", "COILS", "MAGNETS"))
    checks = []

    def check(name, got, want, tol=5e-4):
        checks.append({"name": name, "file": got, "model": want, "rel": got / want - 1, "ok": abs(got / want - 1) < tol})

    inside_rotor = mg.pol(p.r_shaft + 1.0, 0.0)                       # steel on a q-axis, near the shaft
    inside_stator = mg.pol(p.r_so - 0.5, 0.0)
    lam_r = next(f for f in rotor if f.contains(Point(inside_rotor)))
    lam_s = next(f for f in stator if f.contains(Point(inside_stator)))
    check("rotor lamination area [mm^2]", lam_r.area, p.n_sectors * vis["rotor_steel"].area)
    check("stator lamination area [mm^2]", lam_s.area, p.n_sectors * vis["stator_steel"].area)
    want_coil = np.mean([r.polygon().area for r in regions if r.name.startswith("coil_")])
    coil_faces = [f for f in coils]
    check("coil face area, mean [mm^2]", float(np.mean([f.area for f in coil_faces])), float(want_coil))
    check("magnet area, mean [mm^2]", float(np.mean([f.area for f in magnets])), p.mag_w * p.mag_t)
    slots = [f for f in rotor if f.area < 0.5 * lam_r.area and not f.contains(Point(0.0, 0.0))]
    counts = {"magnet slots": len(slots), "magnets": len(magnets), "coil faces": len(coil_faces), "stator laminations": 1}
    want_counts = {"magnet slots": 4 * p.pole_pairs, "magnets": 4 * p.pole_pairs, "coil faces": 2 * p.slots, "stator laminations": 1}
    for k in counts:
        checks.append({"name": f"count of {k}", "file": counts[k], "model": want_counts[k], "rel": 0.0, "ok": counts[k] == want_counts[k]})
    return {"checks": checks, "ok": all(c["ok"] for c in checks)}


def preview(path: Path, png: Path, title: str) -> None:
    """A picture of the file read back: the laminations grey, the magnets red, the coils by phase, the labels left out."""
    import ezdxf
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    phase_col = {"A": "#d08a2e", "B": "#4c9f70", "C": "#3d7ab8"}
    fig, axes = plt.subplots(1, 2, figsize=(13, 6.6), gridspec_kw={"width_ratios": [1, 1]})
    doc = ezdxf.readfile(path)
    where = lambda t: t.dxf.align_point if t.dxf.hasattr("align_point") else t.dxf.insert
    labels = [(t.dxf.text, (where(t).x, where(t).y)) for t in doc.modelspace().query("TEXT")]
    for ax, zoom in zip(axes, (None, 22.5)):
        for f in _faces(path, "STATOR"):
            if f.area > 1000:
                ax.fill(*f.exterior.xy, color="#c9cfd6", lw=0)
                for h in f.interiors:
                    ax.fill(*h.xy, color="white", lw=0)
        for f in _faces(path, "ROTOR"):
            if f.area > 1000:
                ax.fill(*f.exterior.xy, color="#aeb6bf", lw=0)
                for h in f.interiors:
                    ax.fill(*h.xy, color="white", lw=0)
        for f in _faces(path, "COILS"):
            c = f.representative_point()
            phase = min(labels, key=lambda t: math.dist(t[1], (c.x, c.y)))[0][0]
            ax.fill(*f.exterior.xy, color=phase_col.get(phase, "#999999"), lw=0.2, ec="white")
        for f in _faces(path, "MAGNETS"):
            ax.fill(*f.exterior.xy, color="#c0392b", lw=0)
        for e in doc.modelspace().query('*[layer=="ROTOR" | layer=="STATOR"]'):
            pts = list(e.flattening(0.01)) if e.dxftype() != "LINE" else [e.dxf.start, e.dxf.end]
            ax.plot([v.x for v in pts], [v.y for v in pts], color="#2b2f33", lw=0.35)   # x and y only: a Vec3 would add z as a third series
        ax.set_aspect("equal")
        ax.axis("off")
        if zoom is not None:
            c = mg.pol(41.0, zoom)
            ax.set_xlim(c[0] - 15, c[0] + 15)
            ax.set_ylim(c[1] - 15, c[1] + 15)
    axes[0].set_title("the whole motor, read back from the DXF", fontsize=10, loc="left")
    axes[1].set_title("one pole, close up (magnets red, coils by phase: A orange, B green, C blue)", fontsize=10, loc="left")
    fig.suptitle(title, fontsize=11, x=0.01, ha="left")
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--overrides", default="{}", help="the candidate: Params overrides or optimizer variables, as JSON")
    ap.add_argument("--out", required=True, help="the DXF to write")
    ap.add_argument("--title", default="", help="the preview's title")
    args = ap.parse_args(argv)
    path = Path(args.out)
    p = write(json.loads(args.overrides), path)
    result = compare(p, path)
    for c in result["checks"]:
        print(f"{'PASS' if c['ok'] else 'FAIL'}  {c['name']}: file {c['file']:.4f}, model {c['model']:.4f} ({100 * c['rel']:+.4f} %)")
    png = path.with_suffix(".png")
    preview(path, png, args.title or f"{path.name}: {json.loads(args.overrides) or 'the reference'}")
    print(f"wrote {path} and {png}")
    return 0 if result["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
