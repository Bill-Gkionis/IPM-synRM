"""structural.py: every number and figure of the note "04 Structural".

  A. a bond is a spring                -> stiffness and melting point rise together (coursework data)
  B. density from the crystal          -> iron, copper, and silicon steel from their unit cells
  C. a spinning ring, and a disc       -> hoop stress = density x (rim speed)^2; Lame's disc with our bore
  D. growing with heat                 -> what a 100 K rise does to the fits between parts
  E. the property table                -> what CalculiX and ROSS will need (S5, S9)
  F. the bridges, by hand              -> the pull of one pole's outer steel and magnets, over the steel that holds it

Run from the Python to FEMM folder (F reads the rotor's geometry through motorgeo; shapely is in this project's
environment, so a plain uv run is enough):
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # keeps the figures out of the notes' folder
    uv run python materials/structural.py
Outputs: out/materials/structural.json and two figures (also copied into the folder MOTOR_FIGURES_DIR names, and into the notes' folder
only if it says `vault`: see common.save()).
"""
from __future__ import annotations

import math
import sys

import numpy as np

import common as cm

sys.path.insert(0, str(cm.PROC))
import motorgeo as mg  # noqa: E402
from shapely.geometry import LineString  # noqa: E402
from shapely.ops import nearest_points, unary_union  # noqa: E402

N = {}
plt = cm.style()

# ============================================================================ A. a bond is a spring
# Coursework, lecture 3 P35-P36: melting point [deg C] and Young's modulus E [GPa] of some metals.
bonds = {"Fe": (1535, 210), "Ni": (1453, 210), "Cu": (1083, 110), "Al": (660, 70),
         "Zn": (420, 84), "Sn": (232, 43), "Pb": (327, 14)}
# The CRC Handbook's values (84th ed., as tabulated on the elastic-properties data page), where the coursework's differ.
checked = {"Fe": 211, "Ni": 200, "Cu": 130, "Al": 70, "Zn": 108, "Sn": 50, "Pb": 16}
Tm = np.array([v[0] for v in bonds.values()]) + 273.15
E = np.array([v[1] for v in bonds.values()])
slope, icpt = np.polyfit(Tm, E, 1)
a_fe = 0.28665e-9                                   # m, iron's cell edge
N["A_bonds"] = {"course_data": bonds, "fit_GPa_per_K": slope, "fit_intercept_GPa": icpt,
                "correlation_r": float(np.corrcoef(Tm, E)[0, 1]),
                "correlation_r_checked_values": float(np.corrcoef(Tm, list(checked.values()))[0, 1]),
                "bond_stiffness_N_per_m_k_eq_E_a": 200e9 * a_fe,      # cubic lattice of springs: E ~ k / a
                "M15_yield_strain_pct": 100 * 358 / 200e3}

fig, ax = plt.subplots(figsize=(6.4, 4.0))
ax.plot(Tm - 273.15, E, "o", color=cm.CAT[0], ms=7)
for (name, (tm, e)) in bonds.items():
    cm.label_end(ax, tm, e, name, dx=-7 if name == "Ni" else 7, ha="right" if name == "Ni" else "left")
tt = np.linspace(200, 1600, 10) + 273.15
ax.plot(tt - 273.15, slope * tt + icpt, color=cm.MUTED, lw=1.0, ls="--")
ax.set_xlabel("melting point [°C]"); ax.set_ylabel("Young's modulus E [GPa]")
ax.set_xlim(150, 1700); ax.set_ylim(0, 240)
ax.set_title("Stronger bonds: stiffer and harder to melt")
cm.save(fig, "materials-stiffness-melting.png"); plt.close(fig)

# ============================================================================ B. density from the crystal
def density(atoms_per_cell, molar_g, a_nm):
    """rho = (atoms in the cell x mass of one atom) / volume of the cell [g/cm^3]."""
    return atoms_per_cell * molar_g / cm.NA / (a_nm * 1e-7) ** 3

x_si = 0.0617                                              # 3.2 wt% Si as an atomic fraction: ASSUMED (no M-19 chemistry published)
m_fesi = (1 - x_si) * 55.845 + x_si * 28.086
N["B_density"] = {
    "iron_bcc": density(2, 55.845, 0.28665), "copper_fcc": density(4, 63.546, 0.36149),
    "course_copper_fcc_a0.362": density(4, 63.55, 0.362),
    "FeSi_3.2pct_same_lattice": density(2, m_fesi, 0.28665),
    "FeSi_3.2pct_lattice_0.2862": density(2, m_fesi, 0.2862),     # ASSUMED: silicon shrinks the cell a little
    "datasheet_M19": 7.65, "datasheet_N45": [7.5, 7.6],
    "NREL_29ga_M19_sheets_weighed_kg_per_m3": 7450,                # NREL TP-5400-63887, table 3: sheet samples, coating included
}

# ============================================================================ C. a spinning ring, and a disc
# A thin ring of density rho spinning with rim speed v: each half needs a pull towards the centre, which the
# tension across the two cuts supplies. Balance on a half ring gives the hoop stress sigma = rho * v^2.
RHO, r_rotor, r_bore, nu = 7650.0, 45.0e-3, 7.5e-3, 0.3
rpm = np.array([3000, 6000, 10000, 15000, 20000])
w = 2 * math.pi * rpm / 60
v = w * r_rotor
sigma_ring = RHO * v ** 2 / 1e6
N["C_ring"] = {f"{int(n)}rpm": {"omega_rad_per_s": ww, "rim_speed_m_per_s": vv, "hoop_MPa": s, "accel_g": ww ** 2 * r_rotor / 9.81}
               for n, ww, vv, s in zip(rpm, w, v, sigma_ring)}

# Lame: a disc of outer radius b with a hole of radius a, plane stress. The hoop stress is highest at the bore:
#   at the bore  sigma = (3 + nu)/4 rho w^2 [b^2 + (1 - nu)/(3 + nu) a^2]
#   at the rim   sigma = rho w^2 [(3 + nu) a^2 + (1 - nu) b^2] / 4
def lame(w_, a=r_bore, b=r_rotor):
    bore = (3 + nu) / 4 * RHO * w_ ** 2 * (b ** 2 + (1 - nu) / (3 + nu) * a ** 2)
    rim = RHO * w_ ** 2 * ((3 + nu) * a ** 2 + (1 - nu) * b ** 2) / 4
    return bore / 1e6, rim / 1e6
N["C_disc_lame_nu0.3"] = {f"{int(n)}rpm": {"bore_MPa": lame(ww)[0], "rim_MPa": lame(ww)[1]} for n, ww in zip(rpm, w)}
N["C_disc_note"] = "the FEA note checked CalculiX against this disc at 10,000 rpm: bore 14.106 MPa (0.06 %), rim within 0.27 %"

# The limits the optimizer will use until M-19's own data are found (note 04, section 6.4):
YIELD_PLACEHOLDER = 300.0             # MPa: the low end of similar 0.35 mm grades; M-15 (the neighbour grade) has 358
FATIGUE_PLACEHOLDER = 0.4 * 490.0     # MPa: 0.4 x M-15's tensile strength, polished, fully reversed (punched edges: less)
GOODMAN_PEAK = 2 / (1 / FATIGUE_PLACEHOLDER + 1 / 490.0)   # zero-to-peak cycle: amplitude = mean = peak / 2
N["C_limits_MPa"] = {"yield_placeholder": YIELD_PLACEHOLDER, "fatigue_placeholder_before_edges": FATIGUE_PLACEHOLDER,
                     "goodman_peak_zero_to_peak": GOODMAN_PEAK, "M15_yield": 358, "M15_tensile": 490}

fig, ax = plt.subplots(figsize=(6.8, 4.0))
nn = np.linspace(0, 22000, 200)
s_ring = RHO * (2 * math.pi * nn / 60 * r_rotor) ** 2 / 1e6
ax.plot(nn, s_ring, color=cm.CAT[0])
cm.label_end(ax, nn[-1], s_ring[-1], "thin ring, ρv²")
for kf, ls in ((5, "--"), (10, ":")):
    ax.plot(nn, kf * s_ring, color=cm.CAT[0], ls=ls, lw=1.3)
    n_lab = min(nn[-1], 60 / (2 * math.pi * r_rotor) * math.sqrt(455e6 / (kf * RHO)))   # where it leaves the plot
    s_lab = kf * RHO * (2 * math.pi * n_lab / 60 * r_rotor) ** 2 / 1e6
    cm.label_end(ax, n_lab, s_lab, f"x {kf} (an illustrative concentration)")
ax.axhline(YIELD_PLACEHOLDER, color=cm.CRITICAL, lw=1.0)
ax.text(300, YIELD_PLACEHOLDER + 8, f"yield, placeholder for M-19: {YIELD_PLACEHOLDER:.0f} MPa", fontsize=8.5, color=cm.INK2)
ax.axhline(FATIGUE_PLACEHOLDER, color=cm.CRITICAL, lw=1.0, ls="--")
ax.text(300, FATIGUE_PLACEHOLDER + 8, f"fatigue, polished: about {round(FATIGUE_PLACEHOLDER, -1):.0f} MPa", fontsize=8.5, color=cm.INK2)
ax.set_xlim(0, 30000); ax.set_ylim(0, 480)
ax.spines["bottom"].set_bounds(0, 22000); ax.set_xticks([0, 5000, 10000, 15000, 20000])
ax.set_xlabel("speed [rpm]"); ax.set_ylabel("stress [MPa]")
ax.set_title("Stress grows with the square of speed (rotor rim at 45 mm)")
cm.save(fig, "materials-spin-stress.png"); plt.close(fig)

# ============================================================================ D. growing with heat
cte = {"M-19": 12.0e-6,                     # Tata Steel NO35-1750H / NO20-1300H data sheets, 0-100 C (similar grades)
       "copper": 17.0e-6,                   # IEC 60028:1925; C11000 (Aurubis): 17.6e-6
       "aluminium 6061": 22.6e-6,           # NIST cryogenic materials data: slope of the 6061-T6 fit at 293 K
       "N45 along the magnetisation": 7e-6, "N45 across": -1e-6}                 # Arnold N45, 20-200 C
dT = 100.0
D_stator = 128.0
N["D_expansion"] = {
    "housing_minus_stator_mm_on_128mm_per_100K": (cte["aluminium 6061"] - cte["M-19"]) * dT * D_stator,
    "magnet_across_vs_steel_um_on_9.08mm_per_100K": (cte["M-19"] - cte["N45 across"]) * dT * 9.08e3,
    "copper_vs_steel_um_per_100mm_per_100K": (cte["copper"] - cte["M-19"]) * dT * 100e3,
    "copper_vs_steel_um_on_our_20mm_stack_per_100K": (cte["copper"] - cte["M-19"]) * dT * 20e3,
}

# ============================================================================ E. the property table (for CalculiX and ROSS)
N["E_table"] = {   # sources: the makers' sheets named in note 04, section 6.1. E [GPa], nu, density [kg/m^3], strengths [MPa], CTE [1/K]
    "M-19": {"E": "185-220 in similar 0.35 mm grades, lowest along the rolling direction", "nu": "0.29 (handbook)",
             "rho": 7650, "rho_sheets_weighed_NREL": 7450,
             "yield": f"placeholder {YIELD_PLACEHOLDER:.0f}: M-15 358; similar 0.35 mm grades 300-450",
             "fatigue": f"placeholder about {FATIGUE_PLACEHOLDER:.0f} before the punched-edge reduction", "cte": 12.0e-6},
    "N45": {"E": 160, "nu": 0.24, "rho": "7500-7600", "compressive": "850-1050", "tensile": 78,
            "flexural": "245-285 (Eclipse); 180-270 (Alliance)", "cte_par_perp": [7e-6, -1e-6]},
    "copper, C11000 annealed": {"E": 117, "nu": 0.34, "rho": 8900, "yield": 69, "tensile": "179-262",
                                 "elongation_pct": 35, "cte": 17.0e-6},           # Aurubis C11000
    "aluminium 6061-T6": {"E": 70, "nu": "0.33 (handbook)", "rho": 2700, "yield_min": 240, "tensile_min": 290, "cte": 22.6e-6},
    "shaft, C45E (AISI 1045), quenched and tempered, d <= 16 mm": {"E": 212, "nu": "0.29 (handbook)", "rho": 7850,
                                                                   "yield_min": 490, "tensile": "700-850", "k": 41.7,
                                                                   "cte": 12.5e-6},
    "bearings, AISI 52100": {"E": "190-212", "rho": 7830, "k": "20-33", "cte": 12.0e-6},
    "magnet adhesive, STRUCTALIT 5859 (one example)": {"E": 3.0, "Tg_C": 145, "lap_shear_MPa_23C_150C": [18, 15],
                                                      "cte": 62e-6, "k": 0.2},
}

# ============================================================================ F. the bridges, by hand
# Each pole's outer steel (the "cap", outside the V of air bands) and its two magnets hang on three narrow
# strips of steel: the centre rib between the two air bands, and the two outer bridges between each band's
# end and the rotor's surface. Cut the steel across those three strips and the cap comes free: its pull at
# speed is (mass x omega^2 x centroid radius). Spread over the three strips' widths, that is their average
# ("nominal") stress, if they shared it evenly in plain tension. They don't: the outer bridges also bend, and
# the corners concentrate the stress. The real peak needs the finite-element solve (stage S5).
p = mg.Params()
vis = mg.visible(mg.build(p))
steel = vis["rotor_steel"]
band_dn, band_up = vis["barrier_0_dn"], vis["barrier_0_up"]          # the two air bands of pole 0 (d-axis at 22.5 deg)
pa, pb = nearest_points(band_dn, band_up)                            # the centre rib: the narrowest steel between them
u = np.array([pb.x - pa.x, pb.y - pa.y]); rib_w = float(np.hypot(*u)); u /= rib_w
cuts = [LineString([(pa.x - 0.3 * u[0], pa.y - 0.3 * u[1]), (pb.x + 0.3 * u[0], pb.y + 0.3 * u[1])])]
for band in (band_dn, band_up):                                      # each outer bridge: radially, at the band's end
    ends = [math.atan2(y, x) for x, y in band.exterior.coords if math.hypot(x, y) > p.r_ro - p.bridge - 0.01]
    am = (min(ends) + max(ends)) / 2
    cuts.append(LineString([(43.0 * math.cos(am), 43.0 * math.sin(am)), (45.5 * math.cos(am), 45.5 * math.sin(am))]))
pieces = steel.difference(unary_union([c.buffer(0.005) for c in cuts]))
cap = max((g for g in getattr(pieces, "geoms", [pieces])
           if abs(math.degrees(math.atan2(g.centroid.y, g.centroid.x)) - 22.5) < 3), key=lambda g: math.hypot(g.centroid.x, g.centroid.y))
A_cap, r_cap = cap.area * 1e-6, math.hypot(cap.centroid.x, cap.centroid.y) * 1e-3
mag = vis["magnet_0_dn"]
A_mag, r_mag = 2 * p.mag_w * p.mag_t * 1e-6, math.hypot(mag.centroid.x, mag.centroid.y) * 1e-3   # the magnets themselves
RHO_N45 = 7600.0                                                     # Arnold N45
w_strips = (2 * p.bridge + p.rib) * 1e-3                              # m: the three strips side by side
mr = RHO * A_cap * r_cap + RHO_N45 * A_mag * r_mag                    # kg (per metre of stack) x m
N["F_bridges"] = {
    "cap_steel_mm2": A_cap * 1e6, "cap_centroid_mm": r_cap * 1e3, "magnets_mm2": A_mag * 1e6, "magnets_centroid_mm": r_mag * 1e3,
    "mass_per_pole_g_per_20mm": (RHO * A_cap + RHO_N45 * A_mag) * 0.020 * 1e3, "magnet_share_of_pull": RHO_N45 * A_mag * r_mag / mr,
    "strips_mm": {"outer_bridge": p.bridge, "centre_rib": p.rib, "rib_measured": rib_w, "total": w_strips * 1e3},
    "per_speed": {f"{int(n)}rpm": {"pull_kN_per_m": mr * ww ** 2 / 1e3, "pull_N_over_20mm": mr * ww ** 2 * 0.020,
                                   "nominal_MPa": mr * ww ** 2 / w_strips / 1e6, "thin_ring_MPa": RHO * (ww * r_rotor) ** 2 / 1e6,
                                   "factor_to_fatigue_reversed": FATIGUE_PLACEHOLDER / (mr * ww ** 2 / w_strips / 1e6),
                                   "factor_to_goodman_peak": GOODMAN_PEAK / (mr * ww ** 2 / w_strips / 1e6)}
                  for n, ww in zip(rpm, w)},
    "note": "nominal = the three strips sharing the pull evenly in tension; bending and the corners raise the peak (S5)",
}

cm.write_numbers("structural", N)
print(f"E vs melting point: {slope * 1000:.1f} GPa per 1000 K, r = {N['A_bonds']['correlation_r']:.3f} "
      f"(checked values: {N['A_bonds']['correlation_r_checked_values']:.3f}); bond k = E a = {N['A_bonds']['bond_stiffness_N_per_m_k_eq_E_a']:.1f} N/m; "
      f"M-15 yield strain {N['A_bonds']['M15_yield_strain_pct']:.3f} %")
print("density: " + ", ".join(f"{k} {v}" if isinstance(v, (list, int)) else f"{k} {v:.3f}" for k, v in N["B_density"].items()))
for k, vv in N["C_ring"].items():
    d = N["C_disc_lame_nu0.3"][k]
    print(f"  {k}: omega {vv['omega_rad_per_s']:.1f} rad/s, rim {vv['rim_speed_m_per_s']:.1f} m/s, ring {vv['hoop_MPa']:.2f} MPa, "
          f"{vv['accel_g']:,.0f} g; disc: bore {d['bore_MPa']:.3f}, rim {d['rim_MPa']:.3f} MPa")
print("expansion: " + ", ".join(f"{k} {v:.3f}" for k, v in N["D_expansion"].items()))
f_ = N["F_bridges"]
print(f"bridges: cap {f_['cap_steel_mm2']:.1f} mm^2 at {f_['cap_centroid_mm']:.2f} mm, magnets {f_['magnets_mm2']:.2f} mm^2 at "
      f"{f_['magnets_centroid_mm']:.2f} mm; {f_['mass_per_pole_g_per_20mm']:.1f} g per pole over 20 mm; magnets give "
      f"{100 * f_['magnet_share_of_pull']:.0f} % of the pull; strips {f_['strips_mm']}")
for k, vv in f_["per_speed"].items():
    print(f"  {k}: pull {vv['pull_kN_per_m']:.1f} kN/m ({vv['pull_N_over_20mm']:.0f} N over 20 mm), nominal {vv['nominal_MPa']:.1f} MPa "
          f"(thin ring {vv['thin_ring_MPa']:.1f}); factor to fatigue {vv['factor_to_fatigue_reversed']:.1f}, to Goodman's peak {vv['factor_to_goodman_peak']:.1f}")
