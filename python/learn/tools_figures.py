"""Figures for the project note on the tools: the band-free motor, as each tool sees it.

1. motorgeo (Python): the whole cross-section of Params(barrier="none"), drawn from the geometry kernel.
2. Gmsh: the mesh of the same rotor's 90-degree sector, written by motorgeo.to_gmsh. Gmsh's own window renders it
   afterwards, through shot.geo, which this script writes beside the mesh.
3. FEMM 5.0: the field FEMM computed on 2026-11-05 for the design with the band filled with steel
   (out/flux_explorer/g00_steel.npz: the vector potential A on every node). B = curl A per triangle; the flux lines
   are contours of A. No solver runs here: it only reads saved results.

Run, from the project folder: uv run python learn/tools_figures.py [--gmsh]
It writes into out/tools_fig/ (git-ignored). With --gmsh it also starts Gmsh (the one config.toml names) on
shot.geo: a window opens for a moment, renders the mesh, saves "tools - Gmsh mesh.png" and closes. Without --gmsh it
prints that command instead. The notes' copies are made by hand: no script writes into the project notes.
"""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))  # the project folder, wherever this is run from
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import pandas as pd
from shapely import affinity

import motorgeo as mg
import paths

OUT = paths.P / "out" / "tools_fig"
INK = "#1d1d1f"
FILL = {"M-19 Steel": "#b9c2cc", "N45": "#d1495b", "Air": "#ffffff", "Coil": "#e8a33d", "<No Mesh>": "#ffffff"}
PHASE = {"A": "#e8a33d", "B": "#4c9f70", "C": "#3d7ab8"}

# Gmsh reads this: Merge and Print are relative to the .geo file's own folder, so it sits beside sector.msh.
SHOT_GEO = """// Gmsh's own window renders the band-free sector's mesh and saves it as a PNG, then closes.
Merge "sector.msh";
General.GraphicsWidth = 1400;
General.GraphicsHeight = 1100;
General.Axes = 0;
General.SmallAxes = 0;
General.BackgroundGradient = 0;
General.Color.Background = {255, 255, 255};
General.Color.Foreground = {40, 40, 40};
General.Color.Text = {20, 20, 20};
Mesh.SurfaceFaces = 1;
Mesh.SurfaceEdges = 1;
Mesh.Lines = 1;
Mesh.LineWidth = 2;
Mesh.ColorCarousel = 2;
Print.Width = 1400;
Print.Height = 1100;
Draw;
Print "tools - Gmsh mesh.png";
Exit;
"""


def polys(g):
    return [x for x in getattr(g, "geoms", [g]) if not x.is_empty]


def fig_geometry(path):
    p = mg.Params(barrier="none")
    regions = mg.build(p)
    vis = mg.visible(regions)
    fig, ax = plt.subplots(figsize=(9, 9))
    for k in range(p.n_sectors):                                   # 4 copies of the 90-degree sector
        rot = k * p.sector_deg
        for r in regions:
            if r.material == "<No Mesh>":
                continue
            col = PHASE.get(r.circuit, FILL[r.material])
            for g in polys(affinity.rotate(vis[r.name], rot, origin=(0, 0))):
                ax.fill(*g.exterior.xy, color=col, ec="#4a4a4a", lw=0.35)
                for h in g.interiors:
                    ax.fill(*h.xy, color="white", ec="#4a4a4a", lw=0.35)
            if r.material == "N45":
                lx, ly = mg.rot(r.label, rot)
                a = math.radians(r.magdir + rot)
                ax.annotate("", xy=(lx + 2.0 * math.cos(a), ly + 2.0 * math.sin(a)),
                            xytext=(lx - 2.0 * math.cos(a), ly - 2.0 * math.sin(a)),
                            arrowprops=dict(arrowstyle="-|>", color=INK, lw=0.9, mutation_scale=7))
    ax.add_patch(plt.Circle((0, 0), p.r_shaft, color="#ececec", ec="#4a4a4a", lw=0.35))
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in
               (FILL["M-19 Steel"], FILL["N45"], PHASE["A"], PHASE["B"], PHASE["C"])]
    ax.legend(handles, ["M-19 steel (rotor, stator)", "N45 magnets (arrows: magnetisation)",
                        "phase A coils", "phase B coils", "phase C coils"],
              loc="upper right", fontsize=8.5, frameon=False, bbox_to_anchor=(1.02, 1.0))
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("motorgeo (Python): the band-free motor, Params(barrier=\"none\")\n"
                 f"36 slots, 8 poles, stator {2 * p.r_so:g} mm, rotor {2 * p.r_ro:g} mm, gap {p.airgap:g} mm",
                 fontsize=11, color=INK, loc="left")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)


def mesh_for_gmsh(path):
    p = mg.Params(barrier="none")
    n, bnd = mg.to_gmsh(p, mg.build(p, gap="sliding"), str(path), size_max=1.5)   # gmsh.write takes a str
    return n, bnd


def fig_femm(path):
    d = np.load(paths.P / "out" / "flux_explorer" / "g00_steel.npz")
    xy, tri, A = d["xy"].astype(float), d["tri"], d["A"]
    st = pd.read_csv(paths.P / "out" / "flux_explorer" / "states.csv")
    st = st[(st.gid == 0) & (st.band == "steel") & (st.field)].reset_index(drop=True)
    pick = st[(st.I == 10) & (st.gamma == 30) & (st.theta_e == 0)].index[0]   # 10 A, current angle 30 deg, 0 deg el
    row = st.loc[pick]
    a = A[pick].astype(float)                                       # Wb/m, one value per node
    x, y = xy[:, 0], xy[:, 1]
    # B per triangle from the linear interpolation of A: Bx = dA/dy, By = -dA/dx (A is the z component)
    i0, i1, i2 = tri[:, 0], tri[:, 1], tri[:, 2]
    x0, x1, x2, y0, y1, y2 = x[i0], x[i1], x[i2], y[i0], y[i1], y[i2]
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    dAdx = ((a[i1] - a[i0]) * (y2 - y0) - (a[i2] - a[i0]) * (y1 - y0)) / det
    dAdy = ((a[i2] - a[i0]) * (x1 - x0) - (a[i1] - a[i0]) * (x2 - x0)) / det
    B = np.hypot(dAdy, -dAdx) * 1e3                                 # mm coordinates: x1000 gives tesla
    fig, ax = plt.subplots(figsize=(9, 9))
    for k in range(4):                                              # the field repeats every 90 deg (periodic)
        c, s = math.cos(math.radians(90 * k)), math.sin(math.radians(90 * k))
        X, Y = c * x - s * y, s * x + c * y
        T = mtri.Triangulation(X, Y, tri)
        pc = ax.tripcolor(T, facecolors=B, cmap="turbo", vmin=0, vmax=2.0, shading="flat", rasterized=True)
        ax.tricontour(T, a, levels=np.linspace(a.min(), a.max(), 26)[1:-1], colors="#111111", linewidths=0.4)
    cb = fig.colorbar(pc, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("|B| [T]")
    ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("FEMM 5.0: the band-free rotor (its band filled with steel), 10 A at a 30° current angle\n"
                 f"flux density |B| and flux lines; coarse mesh, {len(x):,} nodes per sector; at this rotor position "
                 f"T = {row.T_age:.2f} N·m, ψd = {1e3 * row.psi_d:.1f} mWb, ψq = {1e3 * row.psi_q:.1f} mWb",
                 fontsize=10.5, color=INK, loc="left")
    fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
    return row


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    fig_geometry(OUT / "tools - motorgeo cross-section.png")
    n, bnd = mesh_for_gmsh(OUT / "sector.msh")
    print(f"gmsh mesh: {n} nodes, boundary groups {bnd}")
    row = fig_femm(OUT / "tools - FEMM flux.png")
    print(f"FEMM state: I={row.I} A gamma={row.gamma} theta_e={row.theta_e} rotor_deg={row.rotor_deg:.3f} "
          f"T_age={row.T_age:.4f} psi_d={row.psi_d:.3f} psi_q={row.psi_q:.3f}")
    (OUT / "shot.geo").write_text(SHOT_GEO, encoding="utf-8", newline="\n")
    cmd = [str(paths.tool("gmsh")), str(OUT / "shot.geo")]
    if "--gmsh" in sys.argv[1:]:
        subprocess.run(cmd, check=True, timeout=120)               # renders, saves the PNG, exits by itself
        print(f"Gmsh render: {OUT / 'tools - Gmsh mesh.png'}")
    else:
        print("Gmsh render (run it, or pass --gmsh):", subprocess.list2cmdline(cmd))
    print(f"written to {OUT}")
