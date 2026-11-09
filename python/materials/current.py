"""current.py: every number and figure of the note "02 Current".

The note follows charge through the motor's materials, from one electron upwards:

  A. one electron in copper            -> collisions every tau, drift speed, sigma = n e^2 tau / m
  B. a winding's resistance and loss   -> R = rho l / A, loss density rho J^2
  C. heat in the lattice               -> rho(T), the 234.5 rule, Matthiessen's rule
  D. alloying                          -> Nordheim's rule: why silicon steel resists eddy currents
  E. changing fields                   -> skin depth in copper, steel and magnet
  F. eddy currents in the laminations  -> the classical loss against the steel maker's totals
  G. eddy currents in the magnets      -> why segmenting helps (a scaling, not an absolute)
  H. the insulation                    -> no free electrons; breakdown; a void sees eps_r times the field

Run from the Python to FEMM folder:
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # keeps the figures out of the docs figures folder
    uv run python materials/current.py
Outputs: out/materials/current.json and two figures (also copied into the folder MOTOR_FIGURES_DIR names, and into the repository's docs figures folder
only if it says `docs`: see common.save()).
"""
from __future__ import annotations

import math

import numpy as np

import common as cm

N = {}
plt = cm.style()

# ============================================================================ A. one electron in copper
rho_m_cu, M_cu = 8960.0, 63.546e-3           # kg/m^3, kg/mol: pure copper, for counting its atoms
rho_m_wind = 8900.0                           # kg/m^3: the winding's commercial copper (IEC 60028: 8.89 g/cm^3; C11000)
n_cu = rho_m_cu / M_cu * cm.NA                # one conduction electron per atom
rho_cu20 = 1 / 58e6                           # Ohm m at 20 deg C, IEC 60028 annealed copper [CHECK report B, B13]
alpha20 = 0.00393                             # 1/K at 20 deg C, IEC 60028              [CHECK report B, B13]
tau = cm.ME / (n_cu * cm.E ** 2 * rho_cu20)   # mean time between collisions (Drude)
v_F = 1.57e6                                  # m/s, copper's Fermi speed               [CHECK report B, B15]
mfp = v_F * tau                               # how far an electron flies between collisions
spacing = 0.3615e-9 / math.sqrt(2)            # nearest-neighbour distance in fcc copper (a / sqrt 2)
J = 6e6                                       # A/m^2: 6 A/mm^2 rms, an example (the hub's 100 A-turns at 45 % fill)
v_drift = J / (n_cu * cm.E)
# the coursework's Exercise 2 (L8 P59-P60): 10 A in a 1.6 mm copper wire, n = 8.4e22 per cm^3
J_ex2 = 10 / (math.pi * 0.8e-3 ** 2)
N["A_copper_electron"] = {
    "n_per_m3": n_cu, "rho_20C_Ohm_m": rho_cu20, "tau_s": tau, "v_F_m_per_s": v_F, "mean_free_path_nm": mfp * 1e9,
    "mean_free_path_in_atom_spacings": mfp / spacing, "J_example_A_per_mm2": J / 1e6, "v_drift_mm_per_s": v_drift * 1e3,
    "drift_over_fermi": v_drift / v_F, "minutes_to_drift_1m": 1 / v_drift / 60, "E_field_V_per_m": rho_cu20 * J,
    "ex2_J_A_per_mm2": J_ex2 / 1e6, "ex2_v_drift_m_per_s": J_ex2 / (8.4e28 * cm.E),
}

# ============================================================================ B. a winding's resistance and loss
p_loss20 = rho_cu20 * J ** 2                  # W/m^3: the heat each cubic metre of copper makes
N["B_copper_loss"] = {"loss_W_per_cm3_20C": p_loss20 / 1e6, "loss_W_per_kg_20C": p_loss20 / rho_m_wind}

# ============================================================================ C. heat in the lattice
def rho_cu(T):
    return rho_cu20 * (1 + alpha20 * (T - 20))
alpha0 = alpha20 / (1 - 20 * alpha20)          # the same line referred to 0 deg C (the coursework's T0 = 273 K)
temps = [20, 75, 100, 120, 155, 180]
N["C_copper_hot"] = {"alpha0_per_K": alpha0, "inferred_zero_C": -1 / alpha0,
                     "rho_ratio": {f"{T}C": rho_cu(T) / rho_cu20 for T in temps},
                     "loss_W_per_kg": {f"{T}C": rho_cu(T) * J ** 2 / rho_m_wind for T in temps}}

# ============================================================================ D. alloying: silicon steel
rho_fe20 = 9.71e-8                              # Ohm m, pure iron at 20 deg C                 [CHECK report B, B5]
alpha_fe = 0.00651                              # 1/K at 20 C (report B, B5: HyperPhysics/Giancoli; the coursework: 0.0065)
sigma_m19 = cm.read_fem_blocks()["M-19 Steel"]["Sigma"] * 1e6
rho_m19 = 1 / sigma_m19                          # what FEMM uses for M-19
x_si = 0.0617                                    # atomic fraction for 3.2 wt% Si (magnetic.py) [CHECK report A]
rho_res = rho_m19 - rho_fe20                     # Matthiessen: alloy = pure-iron lattice part + impurity part
beta = rho_res / (x_si * (1 - x_si))             # Nordheim
dRdT = alpha_fe * rho_fe20                       # the lattice part rises as fast as in pure iron
N["D_silicon_steel"] = {
    "rho_iron_uOhm_cm": rho_fe20 * 1e8, "rho_M19_femm_uOhm_cm": rho_m19 * 1e8, "ratio_M19_over_iron": rho_m19 / rho_fe20,
    "rho_residual_uOhm_cm": rho_res * 1e8, "nordheim_beta_Ohm_m": beta,
    "course_beta_CuNi": 2.59e-6, "course_beta_FeCr": 4.44e-6,
    "M19_relative_coeff_pct_per_K": 100 * dRdT / rho_m19, "iron_relative_coeff_pct_per_K": 100 * alpha_fe, "M19_measured_similar_pct_per_K": 0.14,
}

# ============================================================================ E. changing fields: skin depth
def skin(rho, f, mur):
    """delta = sqrt(2 rho / (omega mu0 mur)): the depth where an AC field (and its current) falls to 1/e."""
    return math.sqrt(2 * rho / (2 * math.pi * f * cm.MU0 * mur))

rho_n45 = 130e-8                                 # Ohm m ACROSS the magnetisation (Arnold N45), where the eddy
                                                 # currents flow in a 2D model; along it Arnold gives 150e-8
f1 = 4 * 3000 / 60                               # 200 Hz: 4 pole pairs at 3000 rpm (an example until S1 sets the speed)
f_slot = 9 * f1                                  # the rotor passes 36 slots per turn = 9 per electrical cycle: 1.8 kHz
f_slot2 = 18 * f1                                # the next slot-harmonic pair (17th and 19th): 3.6 kHz
f_pwm = 20e3                                     # an example switching frequency
freqs = {"50Hz": 50, "f1_200Hz": f1, "slot_1800Hz": f_slot, "slot2_3600Hz": f_slot2, "pwm_20kHz": f_pwm}
cases = {"copper_20C": (rho_cu20, 1.0), "copper_120C": (rho_cu(120), 1.0),
         "M19_mur_4405_at_1.3T": (rho_m19, 4405), "M19_mur_1000": (rho_m19, 1000), "N45": (rho_n45, 1.05)}
N["E_skin_depth_mm"] = {c: {k: skin(r, f, mu) * 1e3 for k, f in freqs.items()} for c, (r, mu) in cases.items()}
N["E_sizes_mm"] = {"wire_26AWG": 0.405, "lam_29ga": 0.356, "lam_26ga": 0.470, "lam_24ga": 0.635,
                   "magnet_thick": 1.57, "magnet_wide": 9.08}

fig, ax = plt.subplots(figsize=(7.4, 4.6))
ff = np.logspace(1, 5, 200)
lines = [("copper, 20 °C", rho_cu20, 1.0, cm.CAT[1], "-"), ("N45, across its magnetisation", rho_n45, 1.05, cm.CAT[2], "-"),
         ("M-19, mu_r 1000", rho_m19, 1000, cm.CAT[0], "--"), ("M-19, mu_r 4400 (1.3 T)", rho_m19, 4405, cm.CAT[0], "-")]
for name, r_, mu, col, ls in lines:
    d = [skin(r_, f, mu) * 1e3 for f in ff]
    ax.loglog(ff, d, color=col, ls=ls)
    cm.label_end(ax, ff[-1], d[-1], name)
# each size's label sits where no curve passes (checked on the rendered figure)
for size, txt, x_lab, above in ((0.356, "29 gauge: 0.356 mm", 2.5e3, False), (0.635, "24 gauge (FEMM): 0.635 mm", 1.2e3, False),
                                (9.08, "magnet width: 9.08 mm", 6e3, True)):
    ax.axhline(size, color=cm.AXIS, lw=0.8)
    ax.text(x_lab, size * (1.06 if above else 0.94), txt, fontsize=8, color=cm.INK2, va="bottom" if above else "top")
for f, txt in ((f1, "200 Hz\n(f1 at 3000 rpm)"), (f_slot, "1.8 kHz\n(slot harmonic)"), (f_pwm, "20 kHz\n(PWM)")):
    ax.axvline(f, color=cm.GRID, lw=1.0)
    ax.text(f * 1.05, 0.012, txt, fontsize=8, color=cm.INK2)
ax.set_xlim(10, 1e5); ax.set_ylim(0.01, 100)
ax.set_xlabel("frequency [Hz] (log)"); ax.set_ylabel("skin depth [mm] (log)")
ax.set_title("How deep an alternating field gets in")
cm.save(fig, "materials-skin-depth.png"); plt.close(fig)

# ============================================================================ F. eddy currents in the laminations
# A sheet of thickness d in a field B sin(wt) along it: the loop currents give, per cubic metre,
# P = pi^2 sigma d^2 f^2 B^2 / 6 (valid while d is well under the skin depth).
RHO_STEEL, LB = 7650.0, 0.45359237
def eddy_W_per_kg(d, f, B):
    return math.pi ** 2 * sigma_m19 * d ** 2 * f ** 2 * B ** 2 / 6 / RHO_STEEL
gauges = {"29ga": 0.356e-3, "26ga": 0.470e-3, "24ga": 0.635e-3}
limit_60 = {"29ga": 1.55 / LB, "26ga": 1.65 / LB, "24ga": 2.00 / LB}      # DI-MAX table 5, max at 1.5 T, 60 Hz
N["F_lamination_eddy"] = {g: {"classical_60Hz_W_per_kg": eddy_W_per_kg(d, 60, 1.5), "maker_total_limit_W_per_kg": limit_60[g],
                              "share_of_limit": eddy_W_per_kg(d, 60, 1.5) / limit_60[g],
                              "classical_400Hz_W_per_kg": eddy_W_per_kg(d, 400, 1.5),
                              "d_over_skin_400Hz_mur4405": d / skin(rho_m19, 400, 4405)}
                          for g, d in gauges.items()}

fig, ax = plt.subplots(figsize=(6.6, 3.8))
names = list(gauges)
x = np.arange(len(names))
eddy = [N["F_lamination_eddy"][g]["classical_60Hz_W_per_kg"] for g in names]
rest = [limit_60[g] - e for g, e in zip(names, eddy)]
ax.bar(x, eddy, 0.55, color=cm.CAT[1], label="classical eddy current (computed)")
ax.bar(x, rest, 0.55, bottom=eddy, color=cm.CAT[0], label="the rest: hysteresis + excess", edgecolor=cm.SURFACE, linewidth=2)
for i, g in enumerate(names):
    ax.text(i, limit_60[g] + 0.08, f"{limit_60[g]:.2f} W/kg", ha="center", fontsize=9, color=cm.INK)
    ax.text(i, eddy[i] / 2, f"{100 * eddy[i] / limit_60[g]:.0f} %", ha="center", va="center", fontsize=9, color="#ffffff")
ax.set_xticks(x, [f"{g.replace('ga', ' gauge')}\n{gauges[g] * 1e3:.3f} mm" for g in names])
ax.set_ylabel("W/kg at 1.5 T, 60 Hz"); ax.set_ylim(0, 5.2); ax.grid(axis="x", visible=False)
ax.set_title("M-19's guaranteed maximum loss, and the eddy-current part of it")
ax.legend(loc="upper left")
cm.save(fig, "materials-lamination-eddy.png"); plt.close(fig)

# --- the maker's typical loss, separated. AK Steel Di-Max M-19, 29 gauge, as-sheared, typical total core
# loss in W/kg (EMERF "Lamination Steels", 3rd ed., excerpt for MIT 6.685, 2008, p. 5).
typ = {1.0: {50: 1.09, 60: 1.36, 100: 2.61, 150: 4.50, 200: 6.74, 300: 12.2, 400: 18.9, 600: 35.7,
             1000: 80.0, 1500: 158, 2000: 257},
       1.5: {50: 2.56, 60: 3.19, 100: 6.11, 150: 10.5, 200: 15.8, 300: 28.7, 400: 44.4}}
# Energy per cycle, P/f, against f, in Bertotti's three parts:
#   hysteresis  a        the same every cycle
#   classical   b * f    the loop currents of section F; b computed from sigma and d, NOT fitted
#   excess      c * f^0.5  the currents around the moving domain walls
# With b fixed, a and c are a straight-line least-squares fit of (P/f - b f) against sqrt(f).
# Points up to 400 Hz: there the field still fills a 0.356 mm sheet (d/delta under 1.3, loss error under 0.5 %).
sep = {}
for Bp, pts in typ.items():
    f_ = np.array([f for f in pts if f <= 400.0]); e_ = np.array([pts[f] / f for f in f_])
    b = math.pi ** 2 * sigma_m19 * 0.356e-3 ** 2 * Bp ** 2 / (6 * RHO_STEEL)          # J/kg per Hz, from sigma and d
    (a, c), *_ = np.linalg.lstsq(np.column_stack([np.ones_like(f_), np.sqrt(f_)]), e_ - b * f_, rcond=None)
    fit = a + b * f_ + c * np.sqrt(f_)
    parts = lambda f: {"hysteresis": a * f, "classical": b * f ** 2, "excess": c * f ** 1.5}     # W/kg
    sep[f"{Bp}T"] = {"hysteresis_J_per_kg_cycle": a, "hysteresis_J_per_m3_cycle": a * RHO_STEEL,
                     "classical_J_per_kg_Hz": b, "excess_J_per_kg_sqrtHz": c,
                     "fit_points_Hz": [int(f) for f in f_], "fit_max_error_pct": float(100 * np.max(np.abs(fit / e_ - 1))),
                     "at_60Hz_W_per_kg": parts(60.0), "at_400Hz_W_per_kg": parts(400.0),
                     "classical_share_of_speed_part_60Hz": b * 60 / (b * 60 + c * 60 ** 0.5),
                     "classical_share_of_speed_part_400Hz": b * 400 / (b * 400 + c * 400 ** 0.5)}
s10, s15 = sep["1.0T"]["hysteresis_J_per_kg_cycle"], sep["1.5T"]["hysteresis_J_per_kg_cycle"]
N["F2_loss_separation_29ga_typical"] = {**sep, "steinmetz_n_1.0_to_1.5T": math.log(s15 / s10) / math.log(1.5)}

fig, ax = plt.subplots(figsize=(7.2, 4.4))
for (Bp, pts), col in zip(typ.items(), cm.CAT):
    f_all = np.array([f for f in pts if f <= 400.0]); e_all = np.array([pts[f] / f for f in f_all])
    ax.plot(f_all, e_all * 1e3, "o", color=col, ms=6, zorder=3)
    s = sep[f"{Bp}T"]
    a, b, c = s["hysteresis_J_per_kg_cycle"], s["classical_J_per_kg_Hz"], s["excess_J_per_kg_sqrtHz"]
    ff_ = np.linspace(0, 420, 200)
    ax.plot(ff_, (a + b * ff_ + c * np.sqrt(ff_)) * 1e3, color=col, lw=1.4)
    ax.plot(ff_, (a + b * ff_) * 1e3, color=col, lw=1.2, ls=":")
    ax.plot(ff_, np.full_like(ff_, a * 1e3), color=col, lw=1.0, ls="--")
    above = Bp > 1.2                         # 1.5 T's label sits above its dashed line, 1.0 T's below (free space)
    ax.text(110 if above else 300, a * 1e3 + (1.2 if above else -1.2), f"{Bp} T: hysteresis {a * 1e3:.1f} mJ/kg",
            fontsize=9, color=cm.INK2, va="bottom" if above else "top")
ax.text(10, 122, "solid: hysteresis + classical + excess, fitted to the dots\n"
                 "dotted: hysteresis + classical only (classical computed from σ and the thickness)\n"
                 "dashed: hysteresis alone, the same every cycle", fontsize=8.5, color=cm.INK2, va="top")
ax.set_xlim(0, 420); ax.set_ylim(0, 125)
ax.set_xlabel("frequency [Hz]"); ax.set_ylabel("energy lost per cycle [mJ/kg]")
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([], [], color=col, marker="o", ms=6, lw=1.4, label=f"{Bp} T: the maker's typical loss")
                   for Bp, col in zip(sorted(typ, reverse=True), cm.CAT[1::-1])],
          loc="center right", bbox_to_anchor=(1.0, 0.45))
ax.set_title("Separating M-19's loss: the part per cycle, and the parts that grow with speed")
cm.save(fig, "materials-loss-separation.png"); plt.close(fig)

# --- concept figure: the same changing flux, a solid block against a laminated one
from matplotlib.patches import FancyArrowPatch, Rectangle

def loop(ax, x0, y0, w, h, col):
    """A closed current loop drawn as a rectangle with an arrow on its top edge (counter-clockwise)."""
    ax.add_patch(Rectangle((x0, y0), w, h, fill=False, ec=col, lw=1.4))
    ax.add_patch(FancyArrowPatch((x0 + 0.62 * w, y0 + h), (x0 + 0.38 * w, y0 + h),
                                 arrowstyle="-|>", mutation_scale=9, color=col, lw=1.4))

fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.9))
for ax, sheets, title in ((axes[0], 1, "solid steel: one big loop"),
                          (axes[1], 6, "six insulated sheets: thin loops")):
    width, gap = 6.0, 0.12
    w_sheet = (width - gap * (sheets - 1)) / sheets
    for k in range(sheets):
        x = k * (w_sheet + gap)
        ax.add_patch(Rectangle((x, 0), w_sheet, 6, color="#cde2fb", lw=0))
        for s in (0.18, 0.36) if sheets > 1 else (0.12, 0.27, 0.42):
            loop(ax, x + s * w_sheet, s * (6 if sheets == 1 else 1.2), w_sheet * (1 - 2 * s),
                 6 - 2 * s * (6 if sheets == 1 else 1.2), cm.CAT[1])
    ax.text(6.35, 5.6, "B, changing,\ninto the page", fontsize=8.5, color=cm.INK2, va="top")
    ax.set_xlim(-0.3, 8.6); ax.set_ylim(-0.6, 6.6); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(title, fontsize=10.5, loc="left")
fig.suptitle("Eddy currents: cutting the steel into sheets confines each loop to one thin sheet",
             x=0.02, ha="left", fontsize=11.5, fontweight="bold", color=cm.INK)
cm.save(fig, "materials-eddy-loops.png"); plt.close(fig)

# ============================================================================ G. eddy currents in the magnets
# Same law for any conducting block crossed by a changing field: loss per volume grows with the square of
# the width the currents can circle in. Cut a magnet into k pieces along that width: each piece's loss per
# volume falls by k^2, so the magnet's total falls by k^2 too (while each piece stays under a skin depth).
N["G_magnet_segments"] = {f"{k}_pieces": 1 / k ** 2 for k in (1, 2, 3, 4, 6)}
N["G_magnet_skin_mm"] = {k: skin(rho_n45, f, 1.05) * 1e3 for k, f in freqs.items()}

# ============================================================================ H. the insulation
kT_eV = cm.KB * 300 / cm.E
Eg = 5.0                                          # eV: an example wide gap (the coursework: insulators have Eg > 3 eV)
fraction = math.exp(-Eg / (2 * kT_eV))           # intrinsic carriers ~ exp(-Eg / 2kT)
# Two enamelled wires touching, with a thin air wedge between them: three layers in series, like three
# capacitors. The same D crosses all three, so each takes a voltage proportional to thickness / eps_r.
# The air gap's share times the applied voltage must reach Paschen's minimum (about 327 V) to spark.
eps_enamel, paschen_V, gap = 4.2, 327.0, 10e-6           # polyamide-imide eps_r 3.8-4.2 (report B, B11, B22)
def pdiv(t_enamel_one_side):
    """A first look: one gap, picked near Paschen's minimum, held at Paschen's minimum voltage."""
    share = gap / (gap + 2 * t_enamel_one_side / eps_enamel)
    return paschen_V / share, share

# The wedge holds every gap from zero up. Each gap g sparks at Paschen's V_b(g); it takes a share s(g) of
# the voltage. The pair starts to discharge at the lowest V_b(g) / s(g) over all gaps.
# Townsend's form of Paschen's curve for air: V_b = B p g / (ln(A p g) - ln(ln(1 + 1/gamma))), with
# A = 15 /(cm Torr) = 11.25 /(Pa m), B = 365 V/(cm Torr) = 273.75 V/(Pa m), gamma set so the minimum is 327 V.
A_P, B_P, P_ATM = 11.25, 273.75, 101325.0
LN_G = paschen_V / (math.e * B_P / A_P)                  # = ln(1 + 1/gamma), from V_min = e (B/A) ln(1 + 1/gamma)
def v_paschen(g):
    den = math.log(A_P * P_ATM * g) - math.log(LN_G)
    return B_P * P_ATM * g / den if den > 0.05 else float("inf")
def pdiv_worst(t_enamel_one_side):
    gs = np.linspace(2e-6, 200e-6, 19801)
    v = [v_paschen(g) / (g / (g + 2 * t_enamel_one_side / eps_enamel)) for g in gs]
    i = int(np.argmin(v))
    return v[i], gs[i]
N["H_insulation"] = {"kT_300K_eV": kT_eV, "Eg_example_eV": Eg, "carrier_fraction": fraction,
                     "void_field_factor_eps_r": {"enamel_PAI_4.2": 4.2, "nomex_410_0.25mm_2.7": 2.7},
                     "wire_0.40mm_IEC60317": {"grade1": {"increase_mm": 0.021, "min_breakdown_V_rms": 2300},
                                              "grade2": {"increase_mm": 0.040, "min_breakdown_V_rms": 4400},
                                              "grade3": {"increase_mm": 0.060, "min_breakdown_V_rms": 6600}},
                     "PDIV_estimate_V": {"check_40um_PAI (measured 900-1000 V)": pdiv(40e-6)[0],
                                         "0.40mm_grade1_10.5um": pdiv(0.0105e-3)[0],
                                         "0.40mm_grade2_20um": pdiv(0.020e-3)[0],
                                         "0.40mm_grade3_30um": pdiv(0.030e-3)[0]},
                     "air_share_40um": pdiv(40e-6)[1],
                     "PDIV_worst_gap_V_peak": {k: {"V": pdiv_worst(t)[0], "gap_um": pdiv_worst(t)[1] * 1e6}
                                               for k, t in {"check_40um_PAI (measured 900-1000 V)": 40e-6,
                                                            "0.40mm_grade1_10.5um": 0.0105e-3,
                                                            "0.40mm_grade2_20um": 0.020e-3,
                                                            "0.40mm_grade3_30um": 0.030e-3}.items()},
                     "paschen_min_gap_um": 1e6 * math.e * LN_G / (A_P * P_ATM),
                     "M19_TCR_measured_pct_per_K": 0.14}       # report B, B29: 2.5 % Si steel, 40-100 C

# ============================================================================ I. closed-form checks quoted in the note
# Sections 4.6, 4.7 and 5 of note 02: textbook solutions of section 5's equations, evaluated here.
import cmath

def jn(n, z, terms=60):
    """Bessel function J_n of a complex argument, by its power series (fine for |z| < 10)."""
    s_, term = 0j, (z / 2) ** n / math.factorial(n)
    for m_ in range(terms):
        s_ += term
        term *= -(z / 2) ** 2 / ((m_ + 1) * (m_ + 1 + n))
    return s_

def wire_ac_ratio(radius, rho, f):
    """R_ac / R_dc of a lone round wire: real part of (ka/2) J0(ka) / J1(ka), with k = (1 - j) / delta."""
    x = (1 - 1j) / skin(rho, f, 1.0) * radius
    return (x / 2 * jn(0, x) / jn(1, x)).real

def sheet_factor(xi):
    """Exact eddy loss over the classical formula, same average flux, for a sheet xi skin depths thick."""
    return 3 / xi * (math.sinh(xi) - math.sin(xi)) / (math.cosh(xi) - math.cos(xi))

def plate_factor(w, l, k, terms=199):
    """A thin w x l plate in a uniform AC field normal to it: the eddy loss per volume after cutting w into k strips,
    over the uncut plate's. The currents' stream function solves a Poisson equation on the rectangle; its series
    (odd m, n) gives the loss per area as a sum of 1 / (m^2 n^2 (m^2/a^2 + n^2/b^2)), the constant cancelling."""
    def F(a, b):
        return sum(1 / (m_ ** 2 * n_ ** 2 * (m_ ** 2 / a ** 2 + n_ ** 2 / b ** 2))
                   for m_ in range(1, terms + 1, 2) for n_ in range(1, terms + 1, 2))
    return F(w / k, l) / F(w, l)

r_wire = 0.405e-3 / 2                      # 26 AWG, bare
a_mag, b_mag = 9.08 / 2, 19.14 / 2         # mm: half the magnet's width and half its length (9.08 x 19.14 mm)
N["I_closed_form"] = {
    "wire_26AWG_Rac_over_Rdc_20C": {f"{int(f)}Hz": wire_ac_ratio(r_wire, rho_cu20, f) for f in (200, 1800, 20e3, 100e3)},
    "wire_skin_depth_in_radii_20C": {f"{int(f)}Hz": skin(rho_cu20, f, 1.0) / r_wire for f in (200, 20e3)},
    "sheet_factor_exact_over_classical": {f"xi_{xi}": sheet_factor(xi) for xi in (1.0, 1.6, 2.0, 3.0)},
    "sheet_24ga_200Hz_mur4405_xi": 0.635e-3 / skin(rho_m19, 200, 4405),
    "sheet_centre_field_over_surface_xi_1.6": 1 / abs(cmath.cosh((1 + 1j) * 1.6 / 2)),
    "magnet_elliptical_loop_factor": {f"{k}_pieces": (a_mag ** 2 + b_mag ** 2) / (a_mag ** 2 + k ** 2 * b_mag ** 2)
                                      for k in (2, 3, 4, 6)},
    "magnet_rectangular_plate_factor": {f"{k}_pieces": plate_factor(9.08, 19.14, k) for k in (2, 3, 4, 6)},
}

# ============================================================================ write and print
cm.write_numbers("current", N)
a = N["A_copper_electron"]
print(f"copper: n = {n_cu:.4e}/m^3, tau = {tau:.3e} s, mean free path {a['mean_free_path_nm']:.1f} nm "
      f"(~{a['mean_free_path_in_atom_spacings']:.0f} atom spacings)")
print(f"  at 6 A/mm^2: drift {a['v_drift_mm_per_s']:.3f} mm/s ({a['drift_over_fermi']:.1e} of the Fermi speed), "
      f"1 m in {a['minutes_to_drift_1m']:.0f} min, E = {a['E_field_V_per_m']:.3f} V/m")
print(f"  coursework Ex2: J = {a['ex2_J_A_per_mm2']:.3f} A/mm^2, drift {a['ex2_v_drift_m_per_s']:.3e} m/s")
print(f"  loss at 6 A/mm^2, 20 C: {N['B_copper_loss']['loss_W_per_cm3_20C']:.3f} W/cm^3 = {N['B_copper_loss']['loss_W_per_kg_20C']:.1f} W/kg")
print(f"  alpha0 = {alpha0:.5f}/K -> zero of the line at {-1 / alpha0:.1f} C; rho ratio: "
      + ", ".join(f"{k} {v:.3f}" for k, v in N["C_copper_hot"]["rho_ratio"].items()))
d = N["D_silicon_steel"]
print(f"M-19: rho {d['rho_M19_femm_uOhm_cm']:.1f} uOhm cm = {d['ratio_M19_over_iron']:.1f}x iron; residual {d['rho_residual_uOhm_cm']:.1f}; "
      f"beta {beta:.2e} Ohm m; dR/R {d['M19_relative_coeff_pct_per_K']:.3f} %/K vs iron {d['iron_relative_coeff_pct_per_K']:.3f}")
for c, v in N["E_skin_depth_mm"].items():
    print(f"  skin {c}: " + ", ".join(f"{k} {x:.3f} mm" for k, x in v.items()))
for g, v in N["F_lamination_eddy"].items():
    print(f"  {g}: classical {v['classical_60Hz_W_per_kg']:.3f} of {v['maker_total_limit_W_per_kg']:.2f} W/kg "
          f"({100 * v['share_of_limit']:.0f} %); 400 Hz {v['classical_400Hz_W_per_kg']:.1f} W/kg, d/delta {v['d_over_skin_400Hz_mur4405']:.2f}")
print(f"insulation: kT = {kT_eV:.4f} eV, exp(-Eg/2kT) for 5 eV = {fraction:.1e}")
for Bp, s in N["F2_loss_separation_29ga_typical"].items():
    if not isinstance(s, dict):
        continue
    p60 = s["at_60Hz_W_per_kg"]; tot = sum(p60.values())
    print(f"  loss separation {Bp}: hysteresis {s['hysteresis_J_per_kg_cycle'] * 1e3:.2f} mJ/kg/cycle = {s['hysteresis_J_per_m3_cycle']:.0f} J/m^3; "
          f"fit within {s['fit_max_error_pct']:.2f} %; 60 Hz: " + ", ".join(f"{k} {v:.3f} W/kg ({100 * v / tot:.0f} %)" for k, v in p60.items())
          + f" = {tot:.2f}; classical share of the speed part {100 * s['classical_share_of_speed_part_60Hz']:.0f} % at 60 Hz, "
          f"{100 * s['classical_share_of_speed_part_400Hz']:.0f} % at 400 Hz")
print(f"  Steinmetz n (hysteresis, 1.0 to 1.5 T) = {N['F2_loss_separation_29ga_typical']['steinmetz_n_1.0_to_1.5T']:.2f}")
print("PDIV, one gap of 10 um: " + ", ".join(f"{k} {v:.0f} V" for k, v in N["H_insulation"]["PDIV_estimate_V"].items()))
print("PDIV, worst gap: " + ", ".join(f"{k} {v['V']:.0f} V at {v['gap_um']:.1f} um" for k, v in N["H_insulation"]["PDIV_worst_gap_V_peak"].items())
      + f"; Paschen's minimum at {N['H_insulation']['paschen_min_gap_um']:.1f} um")
print(f"copper loss per kg at 6 A/mm^2 (8,900 kg/m^3): " + ", ".join(f"{k} {v:.1f}" for k, v in N["C_copper_hot"]["loss_W_per_kg"].items()))
print("skin depth N45: " + ", ".join(f"{k} {v:.2f} mm" for k, v in N["G_magnet_skin_mm"].items()))
ic = N["I_closed_form"]
print("wire Rac/Rdc: " + ", ".join(f"{k} {v:.5f}" for k, v in ic["wire_26AWG_Rac_over_Rdc_20C"].items())
      + "; delta in radii: " + ", ".join(f"{k} {v:.1f}" for k, v in ic["wire_skin_depth_in_radii_20C"].items()))
print("sheet factor: " + ", ".join(f"{k} {v:.3f}" for k, v in ic["sheet_factor_exact_over_classical"].items())
      + f"; 24 ga at 200 Hz: xi {ic['sheet_24ga_200Hz_mur4405_xi']:.2f}; centre field at xi 1.6: {ic['sheet_centre_field_over_surface_xi_1.6']:.2f}")
print("magnet pieces (elliptical loops): " + ", ".join(f"{k} 1/{1 / v:.1f}" for k, v in ic["magnet_elliptical_loop_factor"].items()))
print("magnet pieces (rectangular plate): " + ", ".join(f"{k} 1/{1 / v:.1f}" for k, v in N["I_closed_form"]["magnet_rectangular_plate_factor"].items()))
