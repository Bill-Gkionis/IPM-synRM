"""The one-screen tour of motorgeo: build the baseline rotor, check it, write every output, solve a torque sweep in FEMM.

Draws the sector plot, the full DXF for SolidWorks, the SolidWorks equations and the Gmsh mesh for Elmer into out/,
then asks one FEMM of its own (own port, below-normal priority, see femm_session.py) for 30 rotor positions over one
60 deg electrical ripple period at i_q = 10 A, and prints the mean torque and the peak-to-peak ripple.
Run: uv run python example.py   (from this folder; 11 s with FEMM 5.0)
Expect with FEMM 5.0 (coarse mesh, 2026-11-04): "mean 0.729 N.m, ripple 0.158 N.m peak-to-peak". The mean is close to
the torque that check_femm.py gets from the flux linkages (T_dq 0.7236): the ripple averages out over the period.
"""
import os
import femm
import femm_session
import motorgeo as mg

os.makedirs("out", exist_ok=True)
p = mg.Params()                                   # baseline = the hand-drawn FEMM model
regions = mg.build(p)

metrics, problems = mg.check(p, regions)          # pre-reject before any solver runs
if problems:
    raise SystemExit(problems)

mg.plot(p, regions, "out/sector.png")
mg.to_dxf(p, regions, "out/ipm_full.dxf")         # SolidWorks: open or insert into a sketch
mg.to_solidworks_equations(p, "out/equations.txt")  # or drive SolidWorks' own parametric model
mg.to_gmsh(p, mg.build(p, gap="sliding"), "out/ipm_sector.msh")   # for the Elmer tier

fem = os.path.abspath("out/ipm.fem")
femm_session.open_femm()                          # a hidden FEMM of our own: own port, below-normal (desktop first)
try:
    mg.to_femm(p, regions, fem)
    femm.mi_close()
    model = mg.FemmModel(p, fem)
    sweep = [model.solve(th, 0, 10) for th in range(0, 60, 2)]  # one ripple period, 60 deg electrical
finally:
    femm_session.close_femm()

T = [r["T_age"] for r in sweep]
print(f"mean {sum(T) / len(T):.3f} N.m, ripple {max(T) - min(T):.3f} N.m peak-to-peak")
