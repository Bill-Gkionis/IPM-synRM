"""Geometry-only checks, no FEMM needed: build, shapely checks, winding, DXF round trip, Gmsh, rejections."""
import os
import time

import motorgeo as mg

os.makedirs("out", exist_ok=True)
p = mg.Params()

t = time.perf_counter(); R = mg.build(p); t_build = time.perf_counter() - t
t = time.perf_counter(); m, bad = mg.check(p, R); t_check = time.perf_counter() - t
print(f"{len(R)} regions, build {t_build * 1e3:.1f} ms, check {t_check * 1e3:.1f} ms")
for k, v in m.items():
    print(f"  {k:22s} {v:9.3f}")
print("violations:", bad)

axes, kw1 = mg.phase_axes(p)
print(f"phase axes (mech deg): {axes}, kw1 = {kw1:.4f}")
w = mg.winding(p)
sign = lambda s: "+" if s > 0 else "-"
print("slots 0-8 top/bottom:", [f"{w[k][0][0]}{sign(w[k][0][1])}/{w[k][1][0]}{sign(w[k][1][1])}" for k in range(9)])

mg.plot(p, R, "out/sector.png")
mg.to_dxf(p, R, "out/ipm_full.dxf")
rows = mg.compare_dxf(R, "out/ipm_full.dxf")
worst = max(rows, key=lambda r: abs(r[2] / r[1] - 1))
print(f"DXF round trip: {len(rows)} regions, worst {worst[0]}: kernel {worst[1]:.4f}, dxf {worst[2]:.4f} mm2 "
      f"({(worst[2] / worst[1] - 1) * 100:+.3f}%)")

mg.to_solidworks_equations(p, "out/equations.txt")
t = time.perf_counter()
n, bnd = mg.to_gmsh(p, mg.build(p, gap="sliding"), "out/ipm_sector.msh")
print(f"gmsh: {n} nodes in {time.perf_counter() - t:.1f} s, boundary groups {bnd}")

# infeasible candidates are rejected before any solver runs
try:
    mg.build(mg.Params(v_half_deg=55))
except mg.GeometryError as e:
    print("v_half_deg 55 rejected:", e)
p2 = mg.Params(mag_x=38.4)
_, bad2 = mg.check(p2, mg.build(p2))
print("mag_x 38.4 violations:", bad2)
