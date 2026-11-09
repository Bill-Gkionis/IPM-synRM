"""thermal.py: every number and figure of the note "03 Thermal".

  A. temperature is motion                -> kT at the motor's temperatures
  B. who carries heat                     -> Wiedemann-Franz: the electrons' share, against measured k
  C. how fast heat spreads                -> heat capacity, diffusivity, the time to cross 10 mm
  D. layers                               -> a laminated stack and a wound slot as series/parallel layers
  E. where the heat comes from            -> copper and iron loss per kilogram (from notes 01 and 02)
  F. what heat does to every property     -> one chart of relative change against temperature
  G. what heat does to the insulation     -> life halves every ~10 K (Arrhenius)

Run from the Python to FEMM folder:
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # keeps the figures out of the docs figures folder
    uv run python materials/thermal.py
Outputs: out/materials/thermal.json and three figures (also copied into the folder MOTOR_FIGURES_DIR names, and into the repository's docs figures folder
only if it says `docs`: see common.save()).
Values marked [CHECK report A/B] are being cross-checked online and may still change.
"""
from __future__ import annotations

import math

import numpy as np

import common as cm

N = {}
plt = cm.style()

# ============================================================================ A. temperature is motion
N["A_kT_meV"] = {f"{T}C": cm.KB * (T + 273.15) / cm.E * 1e3 for T in (20, 120, 155)}

# ============================================================================ B. who carries heat
# k_e = L0 * sigma * T: the heat the conduction electrons carry, from the electrical conductivity alone.
# 'measured' is the total, electrons + lattice vibrations (phonons).
mats = {   # name: (sigma [S/m] at 20 C, measured k [W/(m K)], source of k)
    "copper": (58.0e6, 391.0, "report B, B14: C11000 (ETP) at 20 C, Aurubis; pure Cu 401 at 300 K"),
    "aluminium 6061-T6": (0.396 * 58e6, 154.0, "report A, A26: 39.6 % IACS; k from NIST's fit at 293 K"),
    "iron (pure)": (1 / 9.71e-8, 80.0, "report B, B5 for rho; k = 80 W/(m K) for pure iron at 300 K (Incropera table A.1)"),
    "M-19 (in the sheet's plane)": (1.9e6, 21.9, "report A, A19: NREL, measured on 29-gauge M19"),
    "N45 (across the magnetisation)": (1 / 130e-8, 5.8 * 1.163, "Arnold N45: 5.8 kcal/(m h C)"),
}
T20 = 293.15
N["B_wiedemann_franz"] = {}
for name, (sig, k_meas, src) in mats.items():
    k_e = cm.LORENZ0 * sig * T20
    N["B_wiedemann_franz"][name] = {"k_electrons": k_e, "k_measured": k_meas, "electron_share": min(k_e / k_meas, 9.99),
                                     "source": src}

fig, ax = plt.subplots(figsize=(7.2, 3.9))
names = list(mats)
y = np.arange(len(names))[::-1]
k_e = [N["B_wiedemann_franz"][n]["k_electrons"] for n in names]
k_m = [N["B_wiedemann_franz"][n]["k_measured"] for n in names]
ax.barh(y + 0.18, k_m, 0.34, color=cm.CAT[0], label="measured (datasheet)")
ax.barh(y - 0.18, k_e, 0.34, color=cm.CAT[1], label="electrons alone (Wiedemann-Franz)")
for yy, a, b in zip(y, k_m, k_e):
    ax.text(a * 1.08, yy + 0.18, f"{a:.0f}", va="center", fontsize=8.5, color=cm.INK2)
    ax.text(b * 1.08, yy - 0.18, f"{b:.0f}", va="center", fontsize=8.5, color=cm.INK2)
ax.set_yticks(y, names); ax.set_xscale("log"); ax.set_xlim(1, 2000)
ax.grid(axis="y", visible=False)
ax.set_xlabel("thermal conductivity k [W/(m K)] (log)")
ax.set_title("Heat rides on the free electrons, mostly")
ax.legend(loc="lower right")
cm.save(fig, "materials-wiedemann-franz.png"); plt.close(fig)

# ============================================================================ C. how fast heat spreads
heat = {   # density [kg/m^3], specific heat [J/(kg K)], k [W/(m K)]
    "copper": (8900, 394, 391),                           # report B, B14: C11000 at 20 C (Aurubis)
    "M-19 (in plane)": (7650, 463, 21.9),                 # report A, A19: NREL measured (29 ga): 21.9 W/mK, 463 J/kgK
    "N45": (7600, 0.11 * 4186.8, 6.2),                    # Arnold: 0.11 cal/(g C), 5.3-5.8 kcal/(m h C)
    "aluminium 6061": (2700, 943, 154),                   # report A, A26: NIST fits at 293 K
    "polyimide film (Kapton HN)": (1420, 1090, 0.12),     # report B, B20/B26: k, c (DuPont); density 1.42 g/cm3 (DuPont)
    "air (still)": (1.20, 1007, 0.0251),                  # report B, B25: 20 C, 1 atm
}
L = 0.010
N["C_diffusion"] = {}
for name, (rho, c, k) in heat.items():
    a = k / (rho * c)
    N["C_diffusion"][name] = {"rho_c_MJ_per_m3K": rho * c / 1e6, "diffusivity_mm2_per_s": a * 1e6, "seconds_for_10mm": L ** 2 / a}

# ============================================================================ D. layers
# A stack of sheets: along the sheets the layers conduct side by side (parallel); across them, one after
# another (series): 1/k = sum(fraction_i / k_i). A thin poor layer dominates the series sum.
k_steel, k_coat = 21.9, 0.5                               # steel: NREL measured; coating ~0.5 W/mK (report B, B23)
t_sheet, t_coat = 0.356e-3, 2e-6                          # 29 gauge, about 1 um of coating on each face
f_coat = t_coat / (t_sheet + t_coat)
k_par = (1 - f_coat) * k_steel + f_coat * k_coat
k_ser = 1 / ((1 - f_coat) / k_steel + f_coat / k_coat)
k_through_measured = {"NREL M19, pressed 138 kPa": 1.68, "NREL M19, pressed 552 kPa": 1.92,
                      "35A300 loose": 2.9, "35A300 welded": 4.34, "35A300 glued": 10.3}   # report A A19; report B B23
contact_mm2K_per_W = 193.0                                 # NREL: per lamination interface at 138 kPa
# A wound slot, crossed by the heat: copper (fill f_cu) and insulation (the rest) in series, the crudest model.
k_cu, k_ins = 391.0, 0.25                                  # insulation: between varnish 0.22 and enamel/liner (report B)
slot = {f"fill_{f}": 1 / (f / k_cu + (1 - f) / k_ins) for f in (0.35, 0.45, 0.55)}
slot_measured = {"random-wound, measured": [0.4, 0.7], "compact, well impregnated": [1.4, 1.7]}   # report B, B24
N["D_layers"] = {"stack_in_plane": k_par, "stack_through": k_ser, "stack_through_measured": k_through_measured,
                 "slot_measured_range": slot_measured,
                 "contact_resistance_mm2K_per_W": contact_mm2K_per_W,
                 # what the measured contact resistance alone does: one interface per 0.356 mm sheet
                 "stack_through_with_contact": 0.356e-3 / (0.356e-3 / k_steel + contact_mm2K_per_W * 1e-6),
                 "slot_series_model": slot}

# ============================================================================ E. where the heat comes from
J = 6e6
rho_cu20, alpha20 = 1 / 58e6, 0.00393
cu_loss = {f"{T}C": rho_cu20 * (1 + alpha20 * (T - 20)) * J ** 2 / 8900 for T in (20, 120)}
iron_limit = 1.55 / 0.45359237                           # M-19 29 gauge, max at 1.5 T, 60 Hz (DI-MAX table 5)
N["E_sources_W_per_kg"] = {"copper_6A_per_mm2": cu_loss, "M19_29ga_1.5T_60Hz_max": iron_limit,
                           "copper_over_iron_20C": cu_loss["20C"] / iron_limit}

# ============================================================================ F. what heat does to every property
T = np.linspace(20, 180, 161)
curves = {
    "copper resistance": 100 * alpha20 * (T - 20),
    "M-19 resistivity": 100 * 0.0014 * (T - 20),         # measured on 2.5 % Si steel, 40-100 C (report B, B29)
    "N45 remanence Br": 100 * -0.0012 * (T - 20),         # Arnold, measured 20-80 C, extended
    "N45 coercivity HcJ": 100 * -0.0062 * (T - 20),
}
N["F_relative_change_pct_at"] = {name: {f"{t}C": float(np.interp(t, T, v)) for t in (80, 120, 155)} for name, v in curves.items()}

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for (name, v), col in zip(curves.items(), cm.CAT):
    t_end = 120 if "N45" in name else 160          # Arnold plots N45 to 120 C; the coefficients hold to 80 C
    ax.plot(T[T <= 80], v[T <= 80], color=col)
    seg = (T >= 80) & (T <= t_end)
    ax.plot(T[seg], v[seg], color=col, ls=":" if "N45" in name else "-")
    v_end = float(np.interp(t_end, T, v))
    cm.label_end(ax, t_end, v_end, f"{name}: {v_end:+.0f} % at {t_end} °C")
ax.axhline(0, color=cm.AXIS, lw=0.8)
ax.text(24, 70, "N45's coefficients are measured from 20 to 80 °C;\ndotted: carried on to 120 °C", fontsize=8.5,
        color=cm.INK2, va="top")
ax.set_xlim(20, 235); ax.set_ylim(-75, 75)
ax.spines["bottom"].set_bounds(20, 160); ax.set_xticks([20, 40, 60, 80, 100, 120, 140, 160])
ax.set_xlabel("temperature [°C]"); ax.set_ylabel("change from 20 °C [%]")
ax.set_title("Heat moves every property the motor relies on")
cm.save(fig, "materials-temperature-effects.png"); plt.close(fig)

# ============================================================================ G. what heat does to the insulation
# Arrhenius: a chemical ageing rate grows as exp(-Ea / kT). Near a class's temperature this looks like
# "life halves every HIC kelvin". IEC 60216's temperature index TI is the temperature for 20,000 hours.
HIC = 10.0                                                # the usual rule; Arrhenius gives 7-11 K (below)
# Arrhenius: life ~ exp(Ea/kT). The temperature step that halves it is HIC = ln2 * k T^2 / Ea (report B, B18).
N["G_HIC_K_at_155C"] = {f"Ea_{ea}eV": math.log(2) * cm.KB * (155 + 273.15) ** 2 / (ea * cm.E) for ea in (1.0, 1.5)}
classes = {"B (130 °C)": 130, "F (155 °C)": 155, "H (180 °C)": 180}
Ths = np.linspace(100, 200, 201)
fig, ax = plt.subplots(figsize=(7.2, 4.2))
N["G_insulation_life_h"] = {}
for (name, TI), col in zip(classes.items(), cm.CAT):
    life = 20000 * 2 ** ((TI - Ths) / HIC)
    ax.semilogy(Ths, life, color=col)
    ax.plot(TI, 20000, "o", color=col, ms=6)
    cm.label_end(ax, TI, 20000, f"class {name}", dx=6, dy=10)
    N["G_insulation_life_h"][name] = {f"{t}C": 20000 * 2 ** ((TI - t) / HIC) for t in (120, 130, 145, 155, 165, 180)}
ax.axhline(20000, color=cm.AXIS, lw=0.8)
ax.text(102, 2000, "grey line: 20,000 h, the life a class is\ndefined to reach at its own temperature", fontsize=8.5,
        color=cm.INK2, va="top")                                   # bottom left: the one corner no line crosses
ax.set_xlim(100, 200); ax.set_ylim(500, 3e6)
ax.set_xlabel("hot-spot temperature [°C]"); ax.set_ylabel("expected life [h] (log)")
ax.set_title(f"Insulation life halves every {HIC:.0f} K")
cm.save(fig, "materials-insulation-life.png"); plt.close(fig)

# ============================================================================ H. one slot, by hand
# How hot is the copper above the tooth around it? motorgeo's slot: radial walls at +-2.57 deg from
# r_shoe = 46.801 to r_wall = 55.744 mm, then a round bottom (the same formula as motorgeo._slot).
a_s = math.radians(2.57)
r_shoe, r_wall = 46.801, 55.744
R_b = r_wall * math.tan(a_s)
slot_mm2 = 0.5 * (r_wall ** 2 - r_shoe ** 2) * math.sin(2 * a_s) + R_b ** 2 / 2 * (math.pi + 2 * a_s + math.sin(2 * a_s))
fill = 0.45
q_per_m = rho_cu20 * (1 + alpha20 * 100) * J ** 2 * fill * slot_mm2 * 1e-6      # W per metre of slot, copper at 120 C
perimeter_m = (2 * (r_wall - r_shoe) + math.pi * R_b) * 1e-3                    # walls + round bottom
flux_W_m2 = q_per_m / perimeter_m
t_liner, k_liner = 0.25e-3, 0.139                                                # Nomex 410, 0.25 mm, at 150 C (report B, B20)
dT_liner = flux_W_m2 * t_liner / k_liner
half_width = r_wall * math.sin(a_s) * 1e-3                                       # half the slot's width at its wall end
q_vol = q_per_m / (slot_mm2 * 1e-6)
k_slot = N["D_layers"]["slot_series_model"]["fill_0.45"]
dT_inside = q_vol * half_width ** 2 / (2 * k_slot)                               # a slab heated inside, cooled at both faces
N["H_one_slot"] = {"slot_area_mm2": slot_mm2, "copper_area_mm2": fill * slot_mm2, "ampere_turns_rms": J * fill * slot_mm2 * 1e-6,
                   "loss_W_per_m_at_120C": q_per_m, "wall_flux_W_per_m2": flux_W_m2, "dT_liner_K": dT_liner,
                   "dT_inside_K": dT_inside, "half_width_mm": half_width * 1e3, "k_slot_used": k_slot}

# ============================================================================ I. one-line sums the note works in the text
# Dulong-Petit: each atom stores about 3 k_B per kelvin, so rho*c = 3 k_B n, with n the atoms per cubic metre.
atoms = {"copper": (8960, 63.546), "iron": (7874, 55.845), "aluminium": (2700, 26.982)}       # kg/m^3, g/mol
dp = {name: {"atoms_per_m3": rho / (M * 1e-3) * cm.NA, "dulong_petit_rho_c_MJ_per_m3K": 3 * cm.KB * rho / (M * 1e-3) * cm.NA / 1e6}
      for name, (rho, M) in atoms.items()}
M_fu = (2 * 144.242 + 14 * 55.845 + 10.81) / 17                    # Nd2Fe14B: the mean mass of its 17 atoms
dp["N45 (Nd2Fe14B)"] = {"atoms_per_m3": 7600 / (M_fu * 1e-3) * cm.NA,
                         "dulong_petit_rho_c_MJ_per_m3K": 3 * cm.KB * 7600 / (M_fu * 1e-3) * cm.NA / 1e6,
                         "measured_rho_c_MJ_per_m3K": N["C_diffusion"]["N45"]["rho_c_MJ_per_m3K"]}
h_ = N["H_one_slot"]
C_slot = h_["copper_area_mm2"] * 1e-6 * N["C_diffusion"]["copper"]["rho_c_MJ_per_m3K"] * 1e6      # J/(m K): copper per metre of slot
R_slot = (h_["dT_liner_K"] + h_["dT_inside_K"]) / h_["loss_W_per_m_at_120C"]                       # K m/W: copper to tooth
a_pi = N["C_diffusion"]["polyimide film (Kapton HN)"]["diffusivity_mm2_per_s"] * 1e-6
steel_mm2 = math.pi / 4 * (128 ** 2 - 91.6 ** 2) - 36 * h_["slot_area_mm2"]                       # the stator's steel section
m_iron = steel_mm2 * 1e-6 * 0.020 * 7650
N["I_worked_in_text"] = {
    "dulong_petit": dp,
    "slot_C_J_per_mK": C_slot, "slot_R_K_m_per_W": R_slot, "slot_time_constant_s": C_slot * R_slot,
    "polyimide_seconds_to_cross": {f"{t * 1e3:g}mm": t ** 2 / a_pi for t in (25e-6, 0.25e-3, 10e-3)},
    "copper_feedback": {"one_over_alpha_K": 1 / alpha20, "rise_100K_from_cold_estimate_K": 100 / (1 + alpha20 * 100),
                        "ratio_at_120C": 1 + alpha20 * 100},
    "stator_iron": {"steel_mm2": steel_mm2, "mass_kg_20mm": m_iron, "W_at_3.42_W_per_kg": m_iron * 3.42,
                    "W_at_3.19_W_per_kg": m_iron * 3.19},
    "housing_m2": {"stack_side": math.pi * 0.128 * 0.020, "assumed": 0.07, "ratio": 0.07 / (math.pi * 0.128 * 0.020)},
    "lorenz_copper_from_391_W_per_mK": 391 / (58e6 * 293.15),
}

# ============================================================================ write and print
cm.write_numbers("thermal", N)
h = N["H_one_slot"]
print(f"slot: {h['slot_area_mm2']:.1f} mm2, copper {h['copper_area_mm2']:.1f} mm2 ({h['ampere_turns_rms']:.0f} A-turns rms), "
      f"{h['loss_W_per_m_at_120C']:.1f} W/m at 120C, wall flux {h['wall_flux_W_per_m2']:.0f} W/m2, "
      f"dT liner {h['dT_liner_K']:.2f} K, inside {h['dT_inside_K']:.2f} K")
print("kT [meV]: " + ", ".join(f"{k} {v:.1f}" for k, v in N["A_kT_meV"].items()))
for n, v in N["B_wiedemann_franz"].items():
    print(f"  W-F {n}: electrons {v['k_electrons']:.1f} vs measured {v['k_measured']:.1f} W/mK (share {v['electron_share']:.2f})")
for n, v in N["C_diffusion"].items():
    print(f"  {n}: rho c {v['rho_c_MJ_per_m3K']:.2f} MJ/m3K, a {v['diffusivity_mm2_per_s']:.3g} mm2/s, 10 mm in {v['seconds_for_10mm']:.3g} s")
d = N["D_layers"]
print(f"stack: in plane {d['stack_in_plane']:.1f}, through {d['stack_through']:.2f} W/mK; slot series: "
      + ", ".join(f"{k} {v:.2f}" for k, v in d["slot_series_model"].items()))
e = N["E_sources_W_per_kg"]
print(f"copper at 6 A/mm2: {e['copper_6A_per_mm2']['20C']:.1f} (20C), {e['copper_6A_per_mm2']['120C']:.1f} W/kg (120C); "
      f"M-19 max {e['M19_29ga_1.5T_60Hz_max']:.2f} W/kg; ratio {e['copper_over_iron_20C']:.1f}")
for n, v in N["F_relative_change_pct_at"].items():
    print(f"  {n}: " + ", ".join(f"{k} {x:+.1f} %" for k, x in v.items()))
for n, v in N["G_insulation_life_h"].items():
    print(f"  life class {n}: " + ", ".join(f"{k} {x:,.0f} h" for k, x in v.items()))
iw = N["I_worked_in_text"]
print("Dulong-Petit rho c: " + ", ".join(f"{k} {v['dulong_petit_rho_c_MJ_per_m3K']:.2f} (n {v['atoms_per_m3']:.2e})" for k, v in iw["dulong_petit"].items()))
print(f"slot: C {iw['slot_C_J_per_mK']:.1f} J/(m K) x R {iw['slot_R_K_m_per_W']:.3f} K m/W = {iw['slot_time_constant_s']:.1f} s; polyimide: "
      + ", ".join(f"{k} {v:.3g} s" for k, v in iw["polyimide_seconds_to_cross"].items()))
print(f"feedback: 1/alpha {iw['copper_feedback']['one_over_alpha_K']:.0f} K, 100 K final from {iw['copper_feedback']['rise_100K_from_cold_estimate_K']:.1f} K; "
      f"stator iron {iw['stator_iron']['mass_kg_20mm']:.3f} kg, {iw['stator_iron']['W_at_3.42_W_per_kg']:.2f} W at 3.42 W/kg; "
      f"housing {iw['housing_m2']['ratio']:.1f}x the stack's side; Lorenz {iw['lorenz_copper_from_391_W_per_mK']:.3e}")
