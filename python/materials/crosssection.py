"""crosssection.py: the materials map of the whole motor, drawn from the same geometry FEMM solves.

motorgeo builds one 90-degree sector (2 poles, 9 slots) with every region tagged by its material.
Rotating that sector four times gives the full cross-section. The shaft and the housing aren't in the
FEMM model (the shaft is a boundary, the housing isn't drawn), so they're sketched around it and
labelled as such.

Run from the Python to FEMM folder (motorgeo needs shapely, which is in this project's environment):
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # keeps the figure out of the docs figures folder
    uv run python materials/crosssection.py
Output: materials-cross-section.png (also copied into the folder MOTOR_FIGURES_DIR names, and into the repository's docs figures folder
only if it says `docs`: see common.save()).
"""
from __future__ import annotations

import math
import sys

import common as cm

sys.path.insert(0, str(cm.PROC))
import motorgeo as mg  # noqa: E402
from shapely.affinity import rotate  # noqa: E402

plt = cm.style()
from matplotlib.patches import Circle, FancyArrowPatch  # noqa: E402

p = mg.Params()
regions = mg.build(p)
vis = mg.visible(regions)
colour = {"M-19 Steel": cm.CAT[0], "Coil": cm.CAT[1], "N45": cm.CAT[2], "Air": cm.SURFACE, "<No Mesh>": cm.SURFACE}

fig, ax = plt.subplots(figsize=(8.6, 7.2))
# the housing (not modelled): a ring around the stator, drawn first so the stator sits on top
ax.add_patch(Circle((0, 0), 70, facecolor=cm.GRID, edgecolor=cm.AXIS, lw=0.8, ls="--"))
ax.add_patch(Circle((0, 0), p.r_so, facecolor=cm.CAT[0], edgecolor="none"))     # steel underneath: hides the
ax.add_patch(Circle((0, 0), p.r_si, facecolor=cm.SURFACE, edgecolor="none"))    # seams between the sectors
ax.add_patch(Circle((0, 0), p.r_ro, facecolor=cm.CAT[0], edgecolor="none"))
for k in range(p.n_sectors):
    for r in regions:
        poly = vis[r.name]
        for g in getattr(poly, "geoms", [poly]):
            if g.is_empty:
                continue
            g = rotate(g, k * p.sector_deg, origin=(0, 0))
            seamless = r.material == "M-19 Steel"            # no outline, so the four sectors join invisibly
            ax.fill(*g.exterior.xy, color=colour[r.material], lw=0 if seamless else 0.25,
                    ec=colour[r.material] if seamless else cm.SURFACE)
            for h in g.interiors:
                ax.fill(*h.xy, color=cm.SURFACE, lw=0)
        if r.material == "N45":                              # magnetisation arrows
            c = mg.rot(r.label, k * p.sector_deg)
            d = mg.pol(1.9, r.magdir + k * p.sector_deg)
            ax.add_patch(FancyArrowPatch((c[0] - d[0], c[1] - d[1]), (c[0] + d[0], c[1] + d[1]),
                                         arrowstyle="-|>", mutation_scale=7, color=cm.INK, lw=0.8))
# the shaft (not in FEMM: the model's inner boundary)
ax.add_patch(Circle((0, 0), p.r_shaft, facecolor=cm.AXIS, edgecolor=cm.MUTED, lw=0.8))


def callout(text, target, xytext):
    ax.annotate(text, xy=target, xytext=xytext, fontsize=9, color=cm.INK2, va="center",
                arrowprops=dict(arrowstyle="-", color=cm.MUTED, lw=0.7, shrinkA=3, shrinkB=0))


def nearest(pick, angle):
    """The label point of the region copy (any of the four sectors) whose angle is closest to `angle`."""
    best = None
    for k in range(p.n_sectors):
        for r in regions:
            if pick(r):
                pt = mg.rot(r.label, k * p.sector_deg)
                err = abs((math.degrees(math.atan2(pt[1], pt[0])) - angle + 180) % 360 - 180)
                if best is None or err < best[0]:
                    best = (err, pt)
    return best[1]


callout("stator: M-19 laminations", mg.pol(60, 40), (78, 52))
callout("copper coils, two layers per slot", nearest(lambda r: r.material == "Coil", 26), (78, 36))
callout("air gap, 0.8 mm", mg.pol(45.4, 13), (78, 20))
callout("N45 magnets: two per pole, in a V\n(arrows: the way each one points)", nearest(lambda r: r.material == "N45", 2), (78, 2))
callout("air band behind each magnet", nearest(lambda r: r.name.startswith("barrier"), -14), (78, -16))
callout("rotor: M-19 laminations", mg.pol(26, -38), (78, -34))
callout("shaft: steel (not in the FEMM model)", mg.pol(4.5, -50), (78, -52))
callout("housing: aluminium (not modelled yet)", mg.pol(67, -52), (78, -66))
ax.set_xlim(-72, 140); ax.set_ylim(-74, 74)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title("The motor's cross-section, coloured by material (36 slots, 8 poles, 128 mm across)", loc="left")
cm.save(fig, "materials-cross-section.png"); plt.close(fig)
