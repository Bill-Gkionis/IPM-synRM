"""Figure for the project note on motorgeo: the same recipe, three sets of numbers, every tag re-made with the shapes.

Run, from the learn folder: uv run python geometry_figure.py <output.png>   (the output folder must exist; the script
imports motorgeo.py through sys.path.insert(0, ".."), which is relative to the folder it is run from)
"""
import math
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, "..")
import motorgeo as mg

INK, MUTED, GRID = "#0b0b0b", "#898781", "#e1e0d9"
FILL = {"M-19 Steel": "#d9dee4", "N45": "#d1495b", "Air": "#ffffff"}

cases = [("1. reference design\nParams()", {}),
         ("2. magnets 1 mm closer to the d-axis\nParams(mag_y=5.0)", {"mag_y": 5.0}),
         ("3. thinner web and band\nParams(web=1.2, bar_t=1.2)", {"web": 1.2, "bar_t": 1.2}),
         ("4. flatter V\nParams(v_half_deg=55)", {"v_half_deg": 55})]
fig, axs = plt.subplots(1, 4, figsize=(13, 7.2))
for ax, (title, x) in zip(axs, cases):
    p = mg.Params(**x)
    to_pole = lambda q: mg.rot(q, -22.5)                       # turn pole 0 so its d-axis is horizontal
    try:
        regions = mg.build(p)
    except mg.GeometryError as e:
        ax.text(0.5, 0.5, f"build() refused:\n{e}\n\nno FEMM run,\nscored infeasible", ha="center", va="center",
                transform=ax.transAxes, fontsize=10, color="#d03b3b")
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.set_xlim(30.5, 46); ax.set_ylim(-16.5, 16.5); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_color("#c3c2b7")
        continue
    metrics, problems = mg.check(p, regions)
    vis = mg.visible(regions)
    for r in regions:
        if r.side != "rotor" or r.name.startswith("gap_"):
            continue
        for g in getattr(vis[r.name], "geoms", [vis[r.name]]):
            if g.is_empty:
                continue
            xy = [to_pole(q) for q in g.exterior.coords]
            ax.fill(*zip(*xy), color=FILL[r.material], ec="#52514e", lw=0.5)
            for h in g.interiors:
                ax.fill(*zip(*[to_pole(q) for q in h.coords]), color="#ffffff", ec="#52514e", lw=0.5)
        if r.name.endswith(("0_up", "0_dn")):
            lx, ly = to_pole(r.label)
            ax.plot(lx, ly, "o", ms=4, color=INK)                # the label point FEMM gets
            if r.material == "N45":
                a = math.radians(r.magdir - 22.5)
                ax.annotate("", xy=(lx + 2.4 * math.cos(a), ly + 2.4 * math.sin(a)), xytext=(lx, ly),
                            arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.2))
                if r.name.endswith("up"):
                    ax.text(0.03, 0.07, f"upper magnet: N45, label at ({r.label[0]:.2f}, {r.label[1]:.2f})\n"
                            f"magdir = 22.5 + ({p.v_half_deg:g} - 90) = {r.magdir % 360:.1f}°",
                            transform=ax.transAxes, fontsize=7.5, color=INK, va="bottom",
                            bbox=dict(fc="#fcfcfb", ec="none", alpha=0.9, pad=1.5))
            elif r.name.endswith("up"):
                ax.text(lx - 0.9, ly, "Air" if r.material == "Air" else r.material, fontsize=8, color=INK,
                        ha="right", va="center")
    ax.axhline(0, color=MUTED, lw=0.8, ls="-.")
    ax.text(43.6, 0.25, "d-axis", fontsize=8, color=MUTED, ha="right")
    verdict = "check(): OK" if not problems else "check(): " + "; ".join(problems)
    ax.text(0.02, 0.02, verdict, transform=ax.transAxes, fontsize=8.5,
            color="#006300" if not problems else "#d03b3b", va="bottom")
    ax.set_title(title, fontsize=10, color=INK, loc="left")
    ax.set_xlim(30.5, 46); ax.set_ylim(-16.5, 16.5); ax.set_aspect("equal")
    ax.tick_params(colors=MUTED, labelsize=8)
    for s in ax.spines.values():
        s.set_color("#c3c2b7")
axs[0].set_ylabel("mm across the pole", color="#52514e", fontsize=9)
for ax in axs:
    ax.set_xlabel("mm from the shaft centre", color="#52514e", fontsize=9)
fig.suptitle("One recipe, new numbers: the shapes, the label points (dots) and the magnet directions (arrows)\n"
             "are re-made together by build(). Pole 0 turned so its d-axis is horizontal; magdir is the angle FEMM gets.",
             fontsize=10.5, color=INK, x=0.01, ha="left")
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=150, facecolor="#fcfcfb")
print("wrote", sys.argv[1])
