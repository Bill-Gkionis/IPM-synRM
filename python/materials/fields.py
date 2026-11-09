"""fields.py: what the saved FEMM solutions say about the materials at work.

The air-band study (barrier_study.py, 2026-10-28) kept six solved 90-degree sectors in out/barrier/fields:
no load, the q-axis current alone, and 10 A at the best current angle (MTPA), each with and without the
air band behind the magnets. femfile.py reads them without FEMM. This script uses them twice:

  1. the picture for the notes: flux lines (lines of constant vector potential A) over the materials;
  2. how each magnet is loaded: B along its own magnetisation direction, triangle by triangle, so the
     note can say not only where the magnets work on average but how far down their weakest spots go.

Run from the Python to FEMM folder, after barrier_study.py has written out/barrier/fields/*.ans:
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # keeps the figure out of the notes' folder
    uv run python materials/fields.py
Outputs: out/materials/fields.json and one figure (also copied into the folder MOTOR_FIGURES_DIR names, and into the notes' folder
only if it says `vault`: see common.save()).
"""
from __future__ import annotations

import math
import sys

import numpy as np

import common as cm

sys.path.insert(0, str(cm.PROC))
import femfile  # noqa: E402

N = {}
plt = cm.style()
FIELDS = cm.PROC / "out" / "barrier" / "fields"


def load(name):
    d = femfile.read_ans(FIELDS / f"{name}.ans")
    bx, by = femfile.element_b(d["xy"], d["A"], d["tri"])
    x, y = d["xy"][d["tri"], 0], d["xy"][d["tri"], 1]
    area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    mat = np.array([d["labels"][k]["material"] for k in d["label"]])
    magdir = np.radians([d["labels"][k]["magdir"] for k in d["label"]])
    return d, bx, by, area, mat, magdir


def weighted_percentile(v, w, q):
    i = np.argsort(v)
    c = np.cumsum(w[i]) / w.sum()
    return float(np.interp(q / 100, c, v[i]))


# ----------------------------------------------------------------------------- 2. how the magnets are loaded
for case in ["noload_with_band", "noload_no_band", "mtpa_with_band", "mtpa_no_band"]:
    d, bx, by, area, mat, magdir = load(case)
    m = mat == "N45"
    b_par = bx[m] * np.cos(magdir[m]) + by[m] * np.sin(magdir[m])     # B along the magnet's own direction
    b_abs = np.hypot(bx[m], by[m])
    N[case] = {"magnet_area_mm2": float(area[m].sum()),
               "B_parallel_median_T": weighted_percentile(b_par, area[m], 50),
               "B_parallel_p5_T": weighted_percentile(b_par, area[m], 5),
               "B_parallel_p1_T": weighted_percentile(b_par, area[m], 1),
               "B_parallel_min_T": float(b_par.min()),
               "B_abs_median_T": weighted_percentile(b_abs, area[m], 50),
               "steel_B_p99_T": weighted_percentile(np.hypot(bx, by)[mat == "M-19 Steel"], area[mat == "M-19 Steel"], 99)}
    # the weakest triangle: how big it is and where, and how much of the magnets' area lies below each B
    # (magnetic.py reads these against the knee at each temperature)
    k = int(np.argmin(b_par))
    cxy = d["xy"][d["tri"][m][k]].mean(axis=0)
    N[case]["weakest_element"] = {"B_parallel_T": float(b_par[k]), "area_mm2": float(area[m][k]),
                                  "share_of_magnet_area_pct": float(100 * area[m][k] / area[m].sum()),
                                  "centroid_mm": [float(cxy[0]), float(cxy[1])], "n_magnet_elements": int(m.sum())}
    N[case]["area_pct_below_B"] = {f"{t:.2f}": float(100 * area[m][b_par < t].sum() / area[m].sum())
                                   for t in np.round(np.arange(0.50, 1.001, 0.01), 2)}

# ----------------------------------------------------------------------------- 1. the picture
fills = {"M-19 Steel": "#cde2fb", "N45": cm.CAT[2], "Coil": "#f6c3ad", "Air": cm.SURFACE, "<No Mesh>": cm.SURFACE}
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.4), gridspec_kw={"width_ratios": [110, 72], "wspace": 0.04})
for ax, case, title in ((axes[0], "noload_with_band", "magnets alone (no current)"),
                        (axes[1], "mtpa_with_band", "magnets + 10 A in the coils (best angle)")):
    d, bx, by, area, mat, magdir = load(case)
    xy, tri = d["xy"], d["tri"]
    for name, col in fills.items():
        sel = mat == name
        if sel.any():
            ax.tripcolor(xy[:, 0], xy[:, 1], tri[sel], facecolors=np.zeros(sel.sum()),
                         cmap=__import__("matplotlib").colors.ListedColormap([col]), edgecolors="none")
    levels = np.linspace(d["A"].min(), d["A"].max(), 34)[1:-1]
    ax.tricontour(xy[:, 0], xy[:, 1], tri, d["A"], levels=levels, colors=cm.INK, linewidths=0.55,
                  negative_linestyles="solid")
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-2, 84); ax.set_ylim(-2, 70)
    ax.set_title(title, fontsize=10.5, loc="left")
    if case != "noload_with_band":
        continue
    # label targets taken from the solution itself, so they land on the right parts
    lab = d["labels"]
    r_of = lambda p: math.hypot(*p)
    magnet = next(l["xy"] for l in lab if l["material"] == "N45")
    band = next(l["xy"] for l in lab if l["material"] == "Air" and r_of(l["xy"]) < 43.5)
    coil = next(l["xy"] for l in lab if l["material"] == "Coil")
    cx, cy = d["xy"][tri].mean(axis=1).T                            # triangle centroids
    rr, bb = np.hypot(cx, cy), np.hypot(bx, by)
    near_surface = (mat == "M-19 Steel") & (rr > 42.5) & (rr < 44.6)
    k_bridge = np.flatnonzero(near_surface)[np.argmax(bb[near_surface])]
    bridge = (cx[k_bridge], cy[k_bridge])
    zone = near_surface & (np.hypot(cx - bridge[0], cy - bridge[1]) < 1.0)   # the bridge around its peak
    N["bridge_zone_B_T"] = {"p50": weighted_percentile(bb[zone], area[zone], 50),
                            "p90": weighted_percentile(bb[zone], area[zone], 90)}
    th_coil = math.degrees(math.atan2(coil[1], coil[0]))
    tooth = (50 * math.cos(math.radians(th_coil + 5)), 50 * math.sin(math.radians(th_coil + 5)))
    gap = (45.4 * math.cos(math.radians(62)), 45.4 * math.sin(math.radians(62)))
    ax.set_xlim(-2, 108)                       # room on the right for the labels
    marks = [("bridge: saturated steel, about 2 T", bridge), ("air gap, 0.8 mm", gap),
             ("air band behind the magnet", band), ("magnet (N45)", magnet),
             ("stator tooth (M-19)", tooth), ("slot: copper coil", coil)]
    marks.sort(key=lambda m: -math.atan2(m[1][1], m[1][0]))      # top to bottom in the order of their angles
    for (text, target), y_text in zip(marks, np.linspace(64, 12, len(marks))):
        ax.annotate(text, xy=target, xytext=(70, y_text), fontsize=8.5, color=cm.INK2, va="center",
                    arrowprops=dict(arrowstyle="-", color=cm.MUTED, lw=0.7, shrinkA=2, shrinkB=1))
fig.suptitle("Flux lines in a quarter of the motor: the steel guides them, the air resists them",
             x=0.02, ha="left", fontsize=11.5, fontweight="bold", color=cm.INK)
cm.save(fig, "materials-flux-lines.png"); plt.close(fig)

cm.write_numbers("fields", N)
for case, v in [(k, v) for k, v in N.items() if "magnet_area_mm2" in v]:
    print(f"{case:>18}: magnets {v['magnet_area_mm2']:.1f} mm^2, B|| median {v['B_parallel_median_T']:.3f} T, "
          f"p5 {v['B_parallel_p5_T']:.3f}, p1 {v['B_parallel_p1_T']:.3f}, min {v['B_parallel_min_T']:.3f}; "
          f"|B| median {v['B_abs_median_T']:.3f}; steel p99 {v['steel_B_p99_T']:.2f} T")
