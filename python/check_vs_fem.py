"""Does motorgeo reproduce the hand-drawn FEMM model? Region areas, the winding, and an overlay plot. No FEMM needed.

The .FEM was saved at a different rotor position than motorgeo draws: its pole-0 d-axis is at 0.5 deg and a slot
centre at 0 deg, while motorgeo puts pole 0 at 22.5 deg and slot 0 at 5 deg. So rotor and stator are shifted separately.

Which .FEM it reads: the reference model of config.toml ([paths] reference_fem, the copy in this folder's reference\\),
or the file given as the first argument:   uv run python check_vs_fem.py [model.FEM]
"""
import math
import os
import sys

import numpy as np
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import polygonize, unary_union

import femfile
import motorgeo as mg
import paths

FEM = sys.argv[1] if len(sys.argv) > 1 else str(paths.path("reference_fem"))
SHIFT = {"rotor": 0.5 - 22.5, "stator": 0.0 - 5.0}      # add to a motorgeo angle to land on the .FEM


def mmf_axes(conductors, pole_pairs):
    """Phase axes [mech deg] from (angle_deg, phase, signed_turns) conductors, as in mg.phase_axes()."""
    th = np.radians(np.arange(0, 360, 0.05))
    axes = {}
    for ph in "ABC":
        F = np.zeros_like(th)
        for a, q, n in conductors:
            if q == ph:
                F -= n * (th > math.radians(a % 360))
        c = ((F - F.mean()) * np.exp(1j * pole_pairs * th)).sum()
        axes[ph] = (math.degrees(np.angle(c)) / pole_pairs) % (360 / pole_pairs)
    seq = round(((axes["B"] - axes["A"]) * pole_pairs) % 360)
    return axes, seq


if __name__ == "__main__":
    os.makedirs("out", exist_ok=True)
    fem = femfile.read_fem(FEM)
    snap = lambda pts: LineString([mg.key(xy) for xy in pts])         # close the 1e-12 gaps between edges
    faces = list(polygonize(unary_union([snap(pl) for pl in femfile.edges_as_polylines(fem, 0.25)])))
    p = mg.Params()
    regions = mg.build(p)
    vis = mg.visible(regions)

    print(f".FEM: {len(fem['points'])} nodes, {len(faces)} faces | motorgeo: {len(regions)} regions (one 90 deg sector)")
    worst = (None, 0.0)
    for r in regions:
        if r.name.startswith(("gap_", "opening_", "coil_")):   # air gap: one region in the .FEM, three here (AGE).
            continue                                             # slots: split differently, compared whole below
        probe = Point(mg.rot(r.label, SHIFT[r.side]))
        face = next((f for f in faces if f.contains(probe)), None)
        whole = r.name in ("rotor_steel", "stator_steel")              # one face spans the full 360 deg
        k_area = vis[r.name].area * (p.n_sectors if whole else 1)
        f_area = face.area if face else float("nan")
        err = f_area / k_area - 1
        if abs(err) > abs(worst[1]):
            worst = (r.name, err)
        if r.name.startswith(("rotor", "stator", "magnet_0_up", "barrier_0_up")):
            print(f"  {r.name:14s} motorgeo {k_area:9.3f}  .FEM {f_area:9.3f} mm2  ({err * 100:+.3f}%)")
    print(f"  worst of the {sum(1 for r in regions if r.name.startswith(('rotor', 'stator', 'magnet', 'barrier')))} "
          f"rotor/magnet/band/stator-steel regions: {worst[0]} {worst[1] * 100:+.3f}%")
    for k in range(round(p.slots / p.n_sectors)):
        mine_slot = unary_union([vis[n] for n in (f"opening_{k}", f"coil_{k}_0", f"coil_{k}_1")])
        shifted = Polygon([mg.rot(xy, SHIFT["stator"]) for xy in mine_slot.exterior.coords])
        parts = sorted((f.area for f in faces if shifted.contains(f.representative_point())), reverse=True)
        if k == 0 or abs(sum(parts) / mine_slot.area - 1) > 1e-3:
            print(f"  slot {k}: motorgeo {mine_slot.area:.3f} mm2 = opening {vis[f'opening_{k}'].area:.3f} + 2 layers of "
                  f"{vis[f'coil_{k}_0'].area:.3f} | .FEM {sum(parts):.3f} mm2 = faces {', '.join(f'{a:.3f}' for a in parts)}")

    # the winding: same machine? motorgeo's own table vs the coil labels in the .FEM
    mine = [(math.degrees(math.atan2(r.label[1], r.label[0])), r.circuit, r.turns) for r in regions if r.circuit]
    ax_m, seq_m = mmf_axes(mine, p.pole_pairs)
    name = {"Phase A": "A", "Phase B": "B", "Phase C": "C"}
    pitch = 360 / p.slots                          # .FEM labels sit anywhere in a slot: snap them to its centre line
    theirs = [(round(math.degrees(math.atan2(lb["xy"][1], lb["xy"][0])) / pitch) * pitch, name[lb["circuit"]], lb["turns"])
              for lb in fem["labels"] if lb["circuit"]]
    ax_f, seq_f = mmf_axes(theirs, p.pole_pairs)
    print(f"winding  motorgeo: {len(mine)} coil sides in the sector, turns {mine[0][2]:+d}, A->B {seq_m} deg elec, "
          f"axes {', '.join(f'{k} {v:.1f}' for k, v in ax_m.items())} mech")
    print(f"         .FEM:     {len(theirs)} coil sides in 360 deg, turns +/-{abs(theirs[0][2])}, A->B {seq_f} deg elec, "
          f"axes {', '.join(f'{k} {v:.1f}' for k, v in ax_f.items())} mech")
    # magnets: does each .FEM magnet point where its pole says it should? (the rule motorgeo uses)
    pole_pitch = 180 / p.pole_pairs
    for lb in fem["labels"]:
        if lb["material"] != "N45":
            continue
        ang = math.degrees(math.atan2(lb["xy"][1], lb["xy"][0]))
        k = round((ang - (22.5 + SHIFT["rotor"])) / pole_pitch)                # nearest pole
        d = 22.5 + SHIFT["rotor"] + k * pole_pitch
        side = 1 if ang > d else -1
        want = (d + side * (p.v_half_deg - 90) + (0 if k % 2 == 0 else 180)) % 360
        if abs((lb["magdir"] - want + 180) % 360 - 180) > 0.01:
            print(f"  magnet at {ang % 360:6.1f} deg (pole d-axis {d % 360:.1f}): magdir {lb['magdir'] % 360:.1f}, "
                  f"the pattern says {want:.1f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 2, figsize=(13, 6.5))
    for ax, side, lim in ((axs[0], "rotor", (27, 46, -15, 17)), (axs[1], "stator", (44, 66, -8, 14))):
        for pl in femfile.edges_as_polylines(fem, 0.25):
            xy = np.array(pl)
            ax.plot(xy[:, 0], xy[:, 1], color="#c43c35", lw=3.2, alpha=0.45, solid_capstyle="round")
        for e, r in mg.unique_edges(regions):
            if r.side != side or r.name.startswith("gap_"):
                continue
            pts = [mg.rot(q, SHIFT[side]) for q in e.points() + [e.b]]
            xy = np.array(pts)
            ax.plot(xy[:, 0], xy[:, 1], color="#1f2937", lw=0.8)
        ax.set_xlim(*lim[:2]); ax.set_ylim(*lim[2:]); ax.set_aspect("equal")
        ax.set_title(f"{side}: .FEM (thick red) vs motorgeo (thin black)", fontsize=10)
        ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    fig.tight_layout()
    fig.savefig("out/check_vs_fem.png", dpi=150)
    print("wrote out/check_vs_fem.png")
