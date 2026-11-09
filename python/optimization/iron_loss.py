"""optimization.iron_loss: the iron loss of the optimized motor. A Bertotti model fitted to the M-19 datasheet, turned into three
frequency-independent sums per region from FEMM's flux density, and the fast function iron_loss(i_d, i_q, rpm, stack_mm).

What it is, and why
    The torque-speed map of the motor (phase 2) needs the loss of the steel at every operating point: at a current (i_d, i_q) and a speed.
    A FEMM solve is magnetostatic: its field does not depend on the speed. The loss does, through the frequency f of the flux in the steel.
    So the work is split where the physics splits. FEMM is asked once, per current, what the flux density does over ONE electrical
    period in every triangle of the steel (optimization/maps_femm.py, task "loss"). From that this file builds, for each of three
    regions (the stator's teeth, the stator's yoke and the rotor), three numbers that do not contain the speed at all,
        S_h [W/Hz]    S_e [W/Hz^2]    S_x [W/Hz^1.5]
    and the loss at ANY speed is then a sum of three powers of the electrical frequency f = p n / 60 (p = 4 pole pairs, n in rpm):
        P_fe(f) = S_h f + S_e f^2 + S_x f^1.5          [W, for the whole machine at the stack length of the model]
    The three names are Simscape's own: the FEM-Parameterized PMSM block takes exactly these coefficients (k_hr, k_Jr, k_er for the rotor and
    k_hs, k_Js, k_es for the stator, as tables over (i_d, i_q); W/Hz, W/Hz^2, W/Hz^1.5), so the maps feed the plant unchanged
    (matlab/maps/, optimization/maps_matlab.py).

The model, built one step from the one before
    1. A sheet of steel in a sinusoidal flux density of peak B and frequency f loses, per kilogram (Bertotti, 1988), three parts:
           P(f, B) = k_h f B^alpha  +  k_e f^2 B^2  +  k_x f^1.5 B^1.5
         hysteresis    the area of the loop, once a cycle: k_h f B^alpha (alpha near 2);
         classical     the currents the changing flux drives round a sheet of thickness d and conductivity sigma. It needs no fit: integrate
                       E = -z dB/dt over the sheet, loss = (sigma d^2 / 12) <(dB/dt)^2> per m^3; for a sinusoid <(dB/dt)^2> = (2 pi f B)^2 / 2, so
                           k_e = pi^2 sigma d^2 / (6 rho_m)        [W/kg per Hz^2 T^2]   (rho_m = 7650 kg/m^3: the maker's density)
         excess        the currents around the moving domain walls, the part the first two do not explain: k_x f^1.5 B^1.5.
    2. The datasheet (AK Steel Di-Max M-19, 29 gauge, 0.356 mm, as-sheared; EMERF's Lamination Steels for MIT 6.685, page 5) prints the total loss in
       W/kg at 13 flux densities (0.1 ... 1.7 T) and up to 11 frequencies (50 ... 2000 Hz): 112 points. k_e is computed from the gauge and FEMM's
       conductivity (1.9 MS/m); k_h, alpha and k_x are fitted by least squares on the LOGARITHM of the loss, so that every point counts by its
       relative error (the table spans four decades). THE FIT IS AT THE DATASHEET'S OWN GAUGE, 0.356 mm (29 gauge): a loss model fitted to one sheet is the
       loss of that sheet.
       THE LAMINATION THICKNESS OF THE STUDY IS 0.3 mm, the thickness in the CAD, and it is the HEADLINE of every loss number here.
       Only the classical term knows the thickness (k_e is proportional to d^2, step 1), so it is scaled: k_e(0.3 mm) = k_e(0.356 mm) (0.3 / 0.356)^2 =
       0.710 k_e. The hysteresis and excess terms are NOT scaled: their dependence on the thickness is weaker (the loop of a sheet hardly depends on it, the
       excess term does a little, through the domain-wall spacing and the surface), and a datasheet of ONE gauge cannot pin it down (the maker's other gauges,
       26 and 24, are not in the excerpt this fit uses). So the 0.3 mm loss is the 0.356 mm loss with only the classical part cut to 0.71: if the real
       thinner sheet also loses a little less hysteresis and excess, the number is a little high, never low from this cause. iron_loss() takes
       lamination_mm = 0.3 by default and reports the datasheet's 0.356 mm on request; maps_report.py prints both.
       FEMM's own M-19 block lists a lamination thickness of 0.635 mm (24 gauge, a leftover of the library entry, Materials guide note 01 section 6.4). It does NOT
       affect any magnetostatic number in these files (it enters only FEMM's AC eddy-current solver, which is not used): the flux maps, the torque and the flux
       density waveforms would be the same with any value there. The thickness of the loss model comes from the fit and the scaling above, never from FEMM's block.
    3. A real flux density is not a sinusoid. Each triangle of the mesh has a flux density VECTOR B(t) = (Bx, By) over the period, and the three
       parts are generalised to any periodic waveform, in the way that reduces to step 1 for a sinusoid (tested in self_test()):
         classical     (sigma d^2 / 12) <|dB/dt|^2> exactly, for any waveform: with the Fourier series B = sum_k Re(c_k e^{i k w t}), the mean square
                       of the derivative is sum_k (k w)^2 |B_k|^2 / 2 (Parseval), so the harmonics count with k^2 (this is why the slot harmonics
                       matter at speed);
         hysteresis    every harmonic is a loop of its own, at its own frequency: k_h f sum_k k |B_k|^alpha, where |B_k| is the vector amplitude of
                       harmonic k (sqrt of the two components' squares). A steady part (the magnets' flux in a rotor bridge) is no loop and costs
                       nothing; an alternating field along any direction costs what the datasheet says; a circular rotating field costs
                       2^(alpha/2) times an alternating one (about twice: the two components add);
         excess        Bertotti's time-domain form, k_x / C_x <|dB/dt|^1.5> with C_x = (2 pi)^1.5 <|cos|^1.5> = 8.763, which for a sinusoid is exactly
                       k_x f^1.5 B^1.5.
       All three scale with f in the powers 1, 2, 1.5 because the waveform in the period keeps its shape when the speed changes: that is what makes
       S_h, S_e and S_x constants. Only the harmonic amplitudes (36 per triangle) and one time-domain average per triangle are stored by FEMM's side.
    4. FEMM's B is the mean over the block, and the steel carries it in the fraction `fill` of the block (FEMM's lamination fill of M-19, 0.98), so the
       steel's own flux density is B / fill and its mass is rho_m fill V. (Checked on 2026-11-09: FEMM treats the fill as steel and insulation in parallel,
       B_eff(H) = f B_steel(H) + (1 - f) mu0 H, so the steel's B is (B - (1 - f) mu0 H) / f: B / f to within 0.3 per cent at the most saturated triangle.)
       The sums are made per triangle (area A, stack length L, n_sectors = 4 copies of
       the 90 degree sector of the model) and added per region.
       The regions: rotor = every steel triangle inside the air gap; teeth = the stator's steel out to the deepest point of the slots; yoke = the rest
       (the ring behind the slots).

What is not in it (each is stated, not hidden)
    * the magnets' eddy loss and the copper's AC loss (other files); mechanical losses; the end regions of a real stack;
    * the punching and stacking of real laminations (an as-sheared Epstein strip is the datasheet's sample; cut edges add a few tens of per cent);
    * the extra hysteresis of a minor loop on a steady bias, and the rotational loss beyond the 2^(alpha/2) rule;
    * temperature (the steel's resistivity rises about 0.14 per cent per K, which lowers the eddy part; the datasheet is at room temperature);
    * the skin effect above the datasheet's 2 kHz: the exact eddy loss of a sheet over the classical one is (3 / xi) (sinh xi - sin xi) / (cosh xi - cos xi)
      with xi = d / delta; it is 0.97 at 2 kHz (mu_r about 2,000) and about 0.7 at 20 kHz, so the classical term of the slot harmonics at high speed is
      a few tens of per cent high (reported by maps_report.py as the share of S_e that sits above 2 kHz at the top speed).

The checks (the limits were written into the code after a first look at how well a Bertotti fit follows the table, about 7 % rms, and before the cross-validation or any FEMM number of
this file was run; they were not changed afterwards)
    * the fit: the relative error on the 112 points, and a cross-validation in which a whole flux density (or a whole frequency) is left out, the model
      is refitted on the rest and asked about the left-out points. Pass: over the INTERIOR folds (flux densities 0.2 ... 1.65 T, frequencies 60 ... 1500 Hz)
      the root-mean-square error is at most TOL_RMS = 10 % and the largest at most TOL_MAX = 20 %. The edge folds (0.1 T, 1.7 T, 50 Hz, 2 kHz) and
      "train up to 400 Hz, predict 600 ... 2000 Hz" ask the model to extrapolate: they are printed, not judged.
    * the free k_e: with k_e fitted too, the fit's k_e over the classical one tells whether the physics of step 2 describes the data (it should be near 1).
    * self_test(): the three operators against the closed forms (a sinusoid along x and along an oblique line, a circular rotating field, a steady bias
      on a ripple, two harmonics), to rounding.
    * the single 1.8 T entry of the table (4.40 W/kg, with no frequency heading) is not used: at 1.7 T the 50 Hz value is 3.29 and the step to 1.8 T would be
      +34 %, against +5 % per 0.05 T below, so its frequency is in doubt.

Run (from the Python to FEMM folder; no FEMM, a second)
    uv run python -m optimization.iron_loss            fit, cross-validate, self-test; writes out\\optimization\\maps\\iron_fit.json and iron_fit_points.csv
Then maps_run.py (the FEMM part), maps_report.py (builds the per-region sums with this file's build_iron_map and the checks on them).

Use
    from optimization import iron_loss
    iron_loss.iron_loss(i_d, i_q, rpm, stack_mm=20.0, turns=10, paths=1, lamination_mm=0.3)
        -> {"teeth": W, "yoke": W, "rotor": W, "total": W}        (arrays in, arrays out; a scalar in gives floats)
    i_d, i_q  the phase current in d-q [A peak, amplitude-invariant] of the winding given by turns (per coil side) and paths (parallel);
              the sums depend on the ampere-turns only, so a winding with s = turns / (10 paths) times the placeholder's series turns is the
              placeholder's map at s times the current (the same rule as the flux maps)
    lamination_mm  0.3 (the CAD's, the headline) or 0.356 (the datasheet's gauge): only the classical sums S_e change, by (lamination_mm / 0.356)^2
    The i_q axis is symmetric (the loss does not know the sign of i_q: the mirror symmetry of the rotor, checked in maps_report.py), so |i_q| is used.
    Outside the solved region (i_d from -3 to +1 times the rated current, |i_q| up to 3 times it) a ValueError, unless extrapolate="clip".
    The speed range of the study ends at the rotor's safe speed, 16,330 rpm (the stress check; a design decision: "about 16,000 rpm"): f = 1,089 Hz there.

Numbers to expect (this PC, 2026-11-09)
    the fit (python -m optimization.iron_loss, about a second): k_h = 0.02005 W/kg/(Hz T^a), a = 1.8842, k_x = 2.940e-4 W/kg/(Hz^1.5 T^1.5), k_e = 5.178e-5 (computed); relative error on the 112
    points rms 6.5 %, largest 15.4 %; with k_e free it comes out at 1.032 times the classical one; the shares at 1 T are 58 / 30 / 12 % (hysteresis / classical / excess) at 200 Hz and
    25 / 64 / 11 % at 1 kHz; the interior cross-validation folds (189 left-out points) rms 6.9 %, largest 16.2 %: PASS (limits 10 % and 20 %); train up to 400 Hz and predict 600 ... 2000 Hz
    (an extrapolation, not judged): rms 10.4 %, largest 20.7 %. The self-test passes (hysteresis and classical to 1e-15, excess to 7e-5).
    the loss of the optimized motor (maps_report.py loss; 20 mm stack; CAD lamination 0.3 mm, the datasheet's 0.356 mm in brackets), no load / at the rated current on the MTPA path / at 3 times it:
        3,000 rpm   6.2 (7.1) W / 9.0 (10.2) W / 13.3 (15.1) W         16,330 rpm   88.5 (113.0) W / 128.7 (164.5) W / 189.3 (241.9) W
    at the rated point and 3,000 rpm: teeth 5.4 W, yoke 3.4 W, rotor 0.20 W (steel masses 0.295, 0.244 and 0.888 kg). The interface: about 35 microseconds a scalar call, 0.8 microseconds a
    point in a vector call of 2,000.

Self-check (the answers follow).
  1. Why is k_e computed and not fitted, and what does the free fit say about that choice?
  2. A rotor bridge sits at 2 T from the magnets, with a ripple of 0.05 T at 1.8 kHz. What does the model charge for it, and why not k_h f 2^alpha?
  3. The speed doubles. Which of S_h, S_e, S_x is unchanged, and by what factor does each part of the loss grow?
  4. The lamination goes from the datasheet's 0.356 mm to the CAD's 0.3 mm. Which sum changes, by what factor, and which two are left alone on purpose?
Answers.
  1. The classical loss follows from Faraday's law and two measured numbers (sigma, d), so fitting it too would only let it absorb the errors of the other
     two terms. The free fit gives k_e at about 1.03 times the classical one: the data agree with the physics.
  2. Hysteresis k_h f_k A^alpha for the ripple's own amplitude at its own frequency (the steady 2 T is no loop), classical k_e (f_k A)^2 and excess
     for its own rate of change. Charging the loop of the 2 T bias at the fundamental's frequency would be charging for a cycle that is never made.
  3. All three S are unchanged (they are the shape of the waveform and the steel's constants). The hysteresis part grows by 2, the excess part by 2^1.5 = 2.83
     and the classical part by 4.
  4. S_e only, by (0.3 / 0.356)^2 = 0.710 (the classical loss goes as the thickness squared: step 1). S_h and S_x are left alone because their dependence on the
     thickness is weaker and one gauge of datasheet cannot pin it: the 0.3 mm figure is therefore, if anything, a little high.
"""
from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.optimize import least_squares

import motorgeo as mg
import paths

# ---------------------------------------------------------------------------------------------- the datasheet
# AK Steel Di-Max M-19, fully processed, 0.014 in (0.36 mm, 29 gauge), as-sheared, typical total core loss [W/kg], Epstein frame, sinusoidal flux.
# Source: EMERF, Lamination Steels 3rd ed., excerpt for MIT 6.685 (2008), page 5 "Total Core Loss - Pc vs. B", the W/kg columns (the W/lb columns
# x 2.204). The file: D:\Storage\Motor Design\Motor library\Datasheets\AK Steel - Di-Max M-19 29 gauge data (EMERF for MIT 6.685, 2008).pdf.
# Transcribed from the page's text (the table is printed with the flux density in gauss as the row's last number). The rows stop at the frequency
# where the maker's data stop. 1 T = 10,000 G.
DATASHEET_F_HZ = (50, 60, 100, 150, 200, 300, 400, 600, 1000, 1500, 2000)
DATASHEET_W_PER_KG = {                                  # B peak [T] -> loss at the first len(row) frequencies of DATASHEET_F_HZ
    0.10: (0.0176, 0.0198, 0.0375, 0.0639, 0.0926, 0.163, 0.247, 0.452, 1.02, 1.98, 3.20),
    0.20: (0.0683, 0.0860, 0.159, 0.262, 0.381, 0.661, 0.994, 1.79, 3.94, 7.43, 11.7),
    0.40: (0.240, 0.295, 0.555, 0.934, 1.37, 2.39, 3.60, 6.52, 14.0, 26.1, 40.8),
    0.70: (0.602, 0.749, 1.43, 2.44, 3.61, 6.44, 9.81, 18.0, 39.1, 74.3, 119),
    1.00: (1.09, 1.36, 2.61, 4.50, 6.74, 12.2, 18.9, 35.7, 80.0, 158, 257),
    1.20: (1.51, 1.89, 3.63, 6.30, 9.46, 17.3, 26.9, 51.8, 120, 240, 395),
    1.30: (1.79, 2.23, 4.28, 7.41, 11.2, 20.3, 31.8, 61.3, 143, 291),
    1.40: (2.14, 2.66, 5.09, 8.82, 13.2, 24.1, 37.5),
    1.50: (2.56, 3.19, 6.11, 10.5, 15.8, 28.7, 44.4),
    1.55: (2.77, 3.44, 6.59, 11.4, 17.0, 30.7, 47.6),
    1.60: (2.96, 3.67, 7.01, 12.0, 18.0),
    1.65: (3.13, 3.89, 7.44, 12.8, 19.1),
    1.70: (3.29, 4.08, 7.80, 13.4, 20.1),
}
# not used: B = 1.8 T prints one pair, 2.00 W/lb = 4.40 W/kg, with no frequency heading (see the module docstring).

RHO_STEEL = 7650.0                                      # kg/m^3, the maker's density for M-19 (7.65 g/cc, page 3 of the datasheet)
SIGMA = mg.M19_STEEL_FEMM42["Sigma"] * 1e6              # S/m: FEMM's M-19 (1.9 MS/m, 52.6 micro-ohm cm)
D_LAM = 0.356e-3                                        # m: 29 gauge, the datasheet's own gauge: the fit is made at this thickness
D_LAM_CAD = 0.300e-3                                    # m: the lamination of the CAD, the study's headline thickness (classical term scaled to it)
LAM_FILL = mg.M19_STEEL_FEMM42["LamFill"]               # 0.98: the steel's share of FEMM's block; the steel carries B / fill
K_E = math.pi ** 2 * SIGMA * D_LAM ** 2 / (6.0 * RHO_STEEL)          # [W/kg per Hz^2 T^2], the classical eddy coefficient
C_X = (2.0 * math.pi) ** 1.5 * math.gamma(1.25) / (math.sqrt(math.pi) * math.gamma(1.75))   # (2 pi)^1.5 <|cos|^1.5> = 8.7633 (Bertotti's constant)
POLE_PAIRS = 4
N_SECTORS = 4                                           # the model is the 90 degree sector: four of them make the machine
STACK_MM_MODEL = 20.0                                   # the model's depth [mm]

TOL_RMS, TOL_MAX = 0.10, 0.20                           # the declared pass limits of the interior cross-validation folds (relative error)
INTERIOR_B = (0.20, 0.40, 0.70, 1.00, 1.20, 1.30, 1.40, 1.50, 1.55, 1.60, 1.65)
INTERIOR_F = (60, 100, 150, 200, 300, 400, 600, 1000, 1500)

OUT = paths.path("out") / "optimization" / "maps"                  # the fit is the steel's, not the design's: iron_fit.json lives here
FIT_JSON = OUT / "iron_fit.json"
MAP_DIR = OUT / "chosen"                                            # the optimized design's maps (maps_report.py writes iron_map.npz here)
IRON_MAP = MAP_DIR / "iron_map.npz"
REGIONS = ("teeth", "yoke", "rotor")
TERMS = ("h", "e", "x")


def datasheet_points() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(f [Hz], B peak [T], loss [W/kg]) of every printed point, as three flat arrays (112 points)."""
    rows = [(f, b, p) for b, losses in DATASHEET_W_PER_KG.items() for f, p in zip(DATASHEET_F_HZ, losses)]
    f, b, p = (np.array(x, float) for x in zip(*rows))
    return f, b, p


# ---------------------------------------------------------------------------------------------- the model and its fit
def bertotti(f, B, kh: float, alpha: float, kx: float, ke: float = K_E):
    """The loss of a sheet in a sinusoidal flux density: hysteresis + classical + excess [W/kg]; f in Hz, B the peak in T."""
    f, B = np.asarray(f, float), np.asarray(B, float)
    return kh * f * B ** alpha + ke * f ** 2 * B ** 2 + kx * f ** 1.5 * B ** 1.5


def fit_bertotti(f, B, P, free_ke: bool = False) -> dict:
    """Least squares on log(model / data), so every point counts by its relative error. k_e is the classical one unless free_ke.
    Returns {"kh", "alpha", "kx", "ke"} (and the fit's own max and rms relative error on the points it was given)."""
    f, B, P = (np.asarray(v, float) for v in (f, B, P))

    def unpack(x):
        return (x[0], x[1], x[2], x[3]) if free_ke else (x[0], x[1], x[2], K_E)

    def residual(x):
        kh, alpha, kx, ke = unpack(x)
        return np.log(bertotti(f, B, kh, alpha, kx, ke) / P)

    x0 = [0.01, 2.0, 1e-3] + ([K_E] if free_ke else [])
    lo = [0.0, 1.0, 0.0] + ([0.0] if free_ke else [])
    hi = [1.0, 3.5, 1.0] + ([1e-3] if free_ke else [])
    scale = [0.01, 1.0, 1e-3] + ([1e-5] if free_ke else [])
    sol = least_squares(residual, x0, bounds=(lo, hi), x_scale=scale)
    kh, alpha, kx, ke = unpack(sol.x)
    err = bertotti(f, B, kh, alpha, kx, ke) / P - 1.0
    return {"kh": float(kh), "alpha": float(alpha), "kx": float(kx), "ke": float(ke),
            "rms": float(np.sqrt(np.mean(err ** 2))), "max": float(np.abs(err).max()), "n": int(len(P))}


def predict(fit: dict, f, B):
    return bertotti(f, B, fit["kh"], fit["alpha"], fit["kx"], fit["ke"])


def cross_validate() -> dict:
    """Leave one flux density (or one frequency) out, refit on the rest, ask about the left-out points (relative error of the model).
    Returns the pooled statistics of the interior folds (the pass/fail ones), and every fold as a row."""
    f, B, P = datasheet_points()
    rows, interior_err = [], []

    def one(name, value, test, interior):
        train = ~test
        fit = fit_bertotti(f[train], B[train], P[train])
        err = predict(fit, f[test], B[test]) / P[test] - 1.0
        rows.append({"fold": name, "left_out": value, "n": int(test.sum()), "interior": interior, "rms": float(np.sqrt(np.mean(err ** 2))),
                     "max": float(np.abs(err).max()), "kh": fit["kh"], "alpha": fit["alpha"], "kx": fit["kx"]})
        if interior:
            interior_err.extend(err.tolist())

    for b in DATASHEET_W_PER_KG:
        one("B", b, np.isclose(B, b), any(np.isclose(b, i) for i in INTERIOR_B))
    for fr in DATASHEET_F_HZ:
        one("f", fr, f == fr, fr in INTERIOR_F)
    # extrapolation in frequency: train up to 400 Hz, predict above (information, not judged)
    train = f <= 400
    fit = fit_bertotti(f[train], B[train], P[train])
    test = f >= 600
    err = predict(fit, f[test], B[test]) / P[test] - 1.0
    rows.append({"fold": "train f <= 400 Hz", "left_out": "f >= 600 Hz", "n": int(test.sum()), "interior": False, "rms": float(np.sqrt(np.mean(err ** 2))),
                 "max": float(np.abs(err).max()), "kh": fit["kh"], "alpha": fit["alpha"], "kx": fit["kx"]})
    e = np.array(interior_err)
    stats = {"n_interior_points_tested": int(len(e)), "rms": float(np.sqrt(np.mean(e ** 2))), "max": float(np.abs(e).max())}
    stats["pass"] = bool(stats["rms"] <= TOL_RMS and stats["max"] <= TOL_MAX)
    return {"interior": stats, "folds": rows}


def eddy_skin_factor(d: float, f, mu_r, rho: float = 1.0 / SIGMA):
    """Exact eddy loss of a sheet over the classical one, (3 / xi) (sinh xi - sin xi) / (cosh xi - cos xi), xi = d / delta,
    delta = sqrt(2 rho / (2 pi f mu0 mu_r)): 1 for a thin sheet, below 1 once the field no longer fills it."""
    mu0 = 4e-7 * math.pi
    delta = np.sqrt(2.0 * rho / (2.0 * math.pi * np.asarray(f, float) * mu0 * mu_r))
    xi = d / delta
    return 3.0 / xi * (np.sinh(xi) - np.sin(xi)) / (np.cosh(xi) - np.cos(xi))


# ---------------------------------------------------------------------------------------------- from waveforms to the three terms
def waveform_terms(bx, by) -> tuple[np.ndarray, np.ndarray]:
    """What the three loss terms need of the flux density waveform of each element, from its N samples over ONE period.

    bx, by   (N, n_elements): the two components of B [T] at N equally spaced instants of one period (N even)
    Returns  Ak     (n_elements, N/2): the vector amplitude of harmonic k = 1 ... N/2, sqrt(Bx_k^2 + By_k^2), with each component's peak amplitude
                    2 |c_k| (the Nyquist harmonic k = N/2, which N samples cannot tell the phase of, counts |c_k|, as a cosine)
             xdot15 (n_elements,): the period mean of |dB/dtau|^1.5, tau the time in periods (B in T): |dB/dtau| from the spectral derivative
                    (every harmonic times 2 pi i k; the Nyquist one has none), the length of the vector (dBx/dtau, dBy/dtau)
    The classical term needs only Ak (sum_k k^2 Ak^2), the hysteresis term only Ak (sum_k k Ak^alpha), the excess term only xdot15."""
    bx, by = np.asarray(bx, float), np.asarray(by, float)
    n = bx.shape[0]
    if n % 2 or n < 4:
        raise ValueError(f"waveform_terms: N = {n} samples, an even number of at least 4 is needed")
    cx, cy = np.fft.rfft(bx, axis=0) / n, np.fft.rfft(by, axis=0) / n          # (n/2 + 1, n_elements)
    power = np.abs(cx) ** 2 + np.abs(cy) ** 2
    amp2 = 4.0 * power
    amp2[-1] = power[-1]                                                         # the Nyquist term appears once
    Ak = np.sqrt(amp2[1:]).T
    k = np.arange(n // 2 + 1, dtype=float)[:, None]
    dx, dy = 2j * math.pi * k * cx, 2j * math.pi * k * cy
    dx[-1] = 0.0
    dy[-1] = 0.0
    rx, ry = np.fft.irfft(dx * n, n=n, axis=0), np.fft.irfft(dy * n, n=n, axis=0)
    xdot15 = np.mean(np.hypot(rx, ry) ** 1.5, axis=0)
    return Ak, xdot15


def element_sums(Ak, xdot15, area_mm2, fit: dict, *, depth_mm: float = STACK_MM_MODEL, sectors: int = N_SECTORS, fill: float = LAM_FILL) -> dict:
    """The three frequency-independent sums of every element, for the whole machine at the model's stack length.

    Ak, xdot15  waveform_terms() of the elements' flux density waveforms (FEMM's B, the block average)
    area_mm2    the triangles' areas; depth_mm the model's depth; sectors how many copies make the machine
    Returns {"h": W/Hz, "e": W/Hz^2, "x": W/Hz^1.5} per element (arrays), the steel's own B being B / fill and its mass rho_m fill V:
        S_h = k_h m sum_k k (A_k / fill)^alpha        S_e = k_e m sum_k k^2 (A_k / fill)^2        S_x = (k_x / C_x) m <|dB/dtau|^1.5> / fill^1.5"""
    Ak = np.asarray(Ak, float) / fill
    k = np.arange(1, Ak.shape[1] + 1, dtype=float)
    mass = RHO_STEEL * fill * np.asarray(area_mm2, float) * 1e-6 * depth_mm * 1e-3 * sectors            # kg of steel in the machine, per element
    return {"h": fit["kh"] * mass * (Ak ** fit["alpha"] * k).sum(axis=1),
            "e": fit["ke"] * mass * (Ak ** 2 * k ** 2).sum(axis=1),
            "x": fit["kx"] / C_X * mass * np.asarray(xdot15, float) / fill ** 1.5}


def loss_per_kg(bx, by, f: float, fit: dict) -> dict:
    """The loss in W/kg of the steel for a waveform sampled over one period at frequency f, the same operators as element_sums() with unit mass
    and no lamination fill (the closed-form tests of self_test() use it). Returns the three parts and their sum."""
    Ak, xdot15 = waveform_terms(bx, by)
    s = element_sums(Ak, xdot15, area_mm2=[1.0], fit=fit, depth_mm=1.0, sectors=1, fill=1.0)
    mass = RHO_STEEL * 1e-6 * 1e-3                                                # the kg that element_sums() weighed for area 1 mm^2, depth 1 mm
    part = {name: float(s[name][0] / mass) * f ** power for name, power in (("h", 1.0), ("e", 2.0), ("x", 1.5))}
    part["total"] = part["h"] + part["e"] + part["x"]
    return part


SELF_TEST_TOL = {"h": 1e-9, "e": 1e-9, "x": 1e-3}        # hysteresis and classical are exact for the sampled harmonics; the excess term is a mean of
                                                          # |cos|^1.5 over N samples, a quadrature of a function with a kink at every zero crossing (error ~ h^2.5)


def self_test(fit: dict | None = None, verbose: bool = True) -> bool:
    """The three operators against closed forms (no FEMM). Returns True if every case agrees: the hysteresis and classical parts to rounding (1e-9
    relative) and the excess part to 1e-3 (it is a mean of |dB/dt|^1.5 over 72 samples: the closed forms are the continuous-time values)."""
    fit = fit or {"kh": 0.02, "alpha": 1.9, "kx": 4e-4, "ke": K_E}
    n, f = 72, 173.0
    tau = np.arange(n) / n
    ok = True

    def case(name, got, want):
        nonlocal ok
        rel = {key: abs(got[key] - want[key]) / max(abs(want[key]), 1e-30) for key in want}
        good = all(rel[key] <= SELF_TEST_TOL[key] for key in want)
        ok &= good
        if verbose:
            print(f"  {'PASS' if good else 'FAIL'}  {name}: " + ", ".join(f"{key} {got[key]:.6g} (closed form {want[key]:.6g}, {rel[key]:.1e})" for key in want))

    s = np.sin(2 * math.pi * tau)
    B0 = 1.3
    want = {"h": fit["kh"] * f * B0 ** fit["alpha"], "e": fit["ke"] * f ** 2 * B0 ** 2, "x": fit["kx"] * f ** 1.5 * B0 ** 1.5}
    case("sinusoid along x", loss_per_kg(B0 * s[:, None], 0 * s[:, None], f, fit), want)
    th = 0.7
    case("the same sinusoid along an oblique line", loss_per_kg(B0 * math.cos(th) * s[:, None], B0 * math.sin(th) * s[:, None], f, fit), want)
    c = np.cos(2 * math.pi * tau)
    want_rot = {"h": want["h"] * 2 ** (fit["alpha"] / 2), "e": 2 * want["e"], "x": fit["kx"] / C_X * (2 * math.pi * B0) ** 1.5 * f ** 1.5}
    case("circular rotating field of radius B0", loss_per_kg(B0 * c[:, None], B0 * s[:, None], f, fit), want_rot)
    bias, ripple = 1.8, 0.07
    want_bias = {"h": fit["kh"] * f * ripple ** fit["alpha"], "e": fit["ke"] * f ** 2 * ripple ** 2, "x": fit["kx"] * f ** 1.5 * ripple ** 1.5}
    case("a steady 1.8 T under a 0.07 T ripple costs the ripple only", loss_per_kg((bias + ripple * s)[:, None], 0 * s[:, None], f, fit), want_bias)
    b1, b3 = 1.1, 0.12
    w2 = b1 * s + b3 * np.sin(3 * 2 * math.pi * tau + 0.4)
    got_two = loss_per_kg(w2[:, None], 0 * w2[:, None], f, fit)
    # the excess term of a non-sinusoid is the mean of |dB/dt|^1.5: its "closed form" here is a brute-force mean on a fine grid
    fine = np.linspace(0, 1, 20000, endpoint=False)
    w2f = b1 * np.sin(2 * math.pi * fine) + b3 * np.sin(3 * 2 * math.pi * fine + 0.4)
    brute = float(np.mean(np.abs(np.gradient(w2f, fine)) ** 1.5)) * fit["kx"] / C_X * f ** 1.5
    want_two = {"h": fit["kh"] * f * (b1 ** fit["alpha"] + 3 * b3 ** fit["alpha"]), "e": fit["ke"] * f ** 2 * (b1 ** 2 + 9 * b3 ** 2), "x": brute}
    case("two harmonics (1st and 3rd; excess against a fine-grid mean)", got_two, want_two)
    return bool(ok)


# ---------------------------------------------------------------------------------------------- the fitted model of the day
def load_fit(path: Path | None = None) -> dict:
    """The fit that main() wrote (iron_fit.json): {"kh", "alpha", "kx", "ke", ...}."""
    with open(path or FIT_JSON, encoding="utf-8") as f:
        return json.load(f)["fit"]


# ---------------------------------------------------------------------------------------------- the map: sums per region over (i_d, i_q)
def point_sums(point_npz: Path, mesh_npz: Path, fit: dict) -> np.ndarray:
    """The (3 regions, 3 terms) sums [W/Hz, W/Hz^2, W/Hz^1.5] of one loss point (a point_*.npz of maps_run's loss stage) on the mesh of mesh.npz, at the datasheet
    gauge and the model's depth: element_sums() per steel triangle, added per region (rows REGIONS, columns TERMS)."""
    mesh, z = np.load(mesh_npz), np.load(point_npz)
    region = mesh["region"]
    sums = element_sums(z["Ak"], z["xdot15"], mesh["area_mm2"], fit, depth_mm=float(mesh["depth_mm"]), sectors=int(mesh["sectors"]))
    return np.array([[sums[t][region == r].sum() for t in TERMS] for r in range(len(REGIONS))])


def build_iron_map(loss_dir: Path, fit: dict, out_npz: Path, rated_current: float) -> dict:
    """Read every point of the loss sweep (loss_dir/point_*.npz and mesh.npz from maps_femm's "loss" task), turn each into its nine sums
    (3 regions x S_h, S_e, S_x) with element_sums(), and write them on the (i_d, i_q) grid to out_npz. Returns a summary dict.
    A point file holds i_d, i_q, Ak (n_steel, N/2), xdot15 (n_steel,). The grid is whatever the files make: its axes are the sorted distinct i_d and
    i_q (a grid with a missing point is refused)."""
    mesh = np.load(loss_dir / "mesh.npz")
    files = sorted(loss_dir.glob("point_*.npz"))
    if not files:
        raise FileNotFoundError(f"no point_*.npz in {loss_dir}")
    pts, exact_d, exact_q = {}, {}, {}
    for path in files:
        z = np.load(path)
        key = (round(float(z["i_d"]), 9), round(float(z["i_q"]), 9))              # nodes are matched at 1e-9 A, but the axes keep the exact values
        exact_d[key[0]], exact_q[key[1]] = float(z["i_d"]), float(z["i_q"])
        pts[key] = point_sums(path, loss_dir / "mesh.npz", fit)
    kd, kq = sorted(exact_d), sorted(exact_q)
    ids, iqs = np.array([exact_d[k] for k in kd]), np.array([exact_q[k] for k in kq])
    if len(pts) != len(ids) * len(iqs):
        raise RuntimeError(f"the loss sweep is not a full grid: {len(pts)} points for {len(ids)} x {len(iqs)}")
    S = np.zeros((len(REGIONS), len(TERMS), len(ids), len(iqs)))
    for (i_d, i_q), s in pts.items():
        S[:, :, kd.index(i_d), kq.index(i_q)] = s
    np.savez(out_npz, id_A=ids, iq_A=iqs, S=S, rated_current_A=rated_current, regions=np.array(REGIONS), terms=np.array(TERMS),
             kh=fit["kh"], alpha=fit["alpha"], kx=fit["kx"], ke=fit["ke"], fill=LAM_FILL, rho_steel=RHO_STEEL, depth_mm=float(mesh["depth_mm"]),
             sectors=int(mesh["sectors"]), pole_pairs=POLE_PAIRS, fit_gauge_mm=D_LAM * 1e3)         # S_e is stored at the datasheet gauge (the fit's)
    return {"points": len(pts), "i_d": ids.tolist(), "i_q": iqs.tolist(), "file": str(out_npz)}


def eddy_scale(lamination_mm: float, fit_gauge_mm: float = D_LAM * 1e3) -> float:
    """The factor of the classical sums S_e when the lamination is lamination_mm instead of the gauge the fit was made at: (d / d_fit)^2
    (0.7101 for the CAD's 0.3 mm against the datasheet's 0.356 mm). Hysteresis and excess are left unscaled (module docstring)."""
    if not lamination_mm > 0:
        raise ValueError(f"lamination_mm = {lamination_mm!r}: a positive thickness is needed")
    return (lamination_mm / fit_gauge_mm) ** 2


class IronLossMap:
    """The nine sums as functions of the current: bicubic splines over (i_d, i_q) of the placeholder winding (10 turns a coil side, one path),
    the stack length of the model (20 mm) and the datasheet gauge (0.356 mm). __call__ turns them into watts at a speed, for any winding, stack and
    lamination thickness (default the CAD's 0.3 mm: S_e scaled by (0.3 / 0.356)^2)."""

    def __init__(self, npz: Path | None = None):
        z = np.load(npz or IRON_MAP)
        self.id_A, self.iq_A = z["id_A"].astype(float), z["iq_A"].astype(float)
        self.S = z["S"].astype(float)                          # (3, 3, n_d, n_q), iq >= 0
        self.rated = float(z["rated_current_A"])
        self.depth_mm = float(z["depth_mm"])
        self.pole_pairs = int(z["pole_pairs"])
        self.fit_gauge_mm = float(z["fit_gauge_mm"]) if "fit_gauge_mm" in z.files else D_LAM * 1e3
        iq_full = np.concatenate([-self.iq_A[:0:-1], self.iq_A])           # mirror the axis: the loss is even in i_q, so the spline is too
        self._spl = {}
        for r in range(3):
            for t in range(3):
                full = np.concatenate([self.S[r, t][:, :0:-1], self.S[r, t]], axis=1)
                self._spl[r, t] = RectBivariateSpline(self.id_A, iq_full, full, kx=3, ky=3, s=0)
        self.range_d = (self.id_A[0], self.id_A[-1])
        self.range_q = self.iq_A[-1]

    def sums(self, i_d, i_q, *, turns: int = 10, paths: int = 1, extrapolate: str = "raise") -> np.ndarray:
        """(3 regions, 3 terms, ...) in W/Hz^power for the model's stack: S(i) of the winding = S_placeholder(s i), s = turns / (10 paths)."""
        s = turns / (10.0 * paths)
        i_d, i_q = np.broadcast_arrays(np.asarray(i_d, float) * s, np.abs(np.asarray(i_q, float)) * s)
        eps = 1e-9 * self.rated
        outside = (i_d < self.range_d[0] - eps) | (i_d > self.range_d[1] + eps) | (i_q > self.range_q + eps)
        if outside.any():
            if extrapolate == "raise":
                raise ValueError(f"iron_loss: {int(outside.sum())} point(s) outside the solved region (i_d {self.range_d[0] / s:.2f} ... {self.range_d[1] / s:.2f} A, "
                                 f"|i_q| up to {self.range_q / s:.2f} A for this winding); pass extrapolate='clip' to clamp them")
            i_d, i_q = np.clip(i_d, *self.range_d), np.clip(i_q, 0.0, self.range_q)
        out = np.empty((3, 3) + i_d.shape)
        for (r, t), spl in self._spl.items():
            out[r, t] = np.maximum(spl.ev(i_d.ravel(), i_q.ravel()).reshape(i_d.shape), 0.0)
        return out

    def __call__(self, i_d, i_q, rpm, stack_mm: float = STACK_MM_MODEL, *, turns: int = 10, paths: int = 1, lamination_mm: float = D_LAM_CAD * 1e3,
                 extrapolate: str = "raise") -> dict:
        S = self.sums(i_d, i_q, turns=turns, paths=paths, extrapolate=extrapolate)
        f = self.pole_pairs * np.asarray(rpm, float) / 60.0
        scale = stack_mm / self.depth_mm
        k_e = eddy_scale(lamination_mm, self.fit_gauge_mm)
        region = {name: scale * (S[r, 0] * f + k_e * S[r, 1] * f ** 2 + S[r, 2] * f ** 1.5) for r, name in enumerate(REGIONS)}
        region["total"] = region["teeth"] + region["yoke"] + region["rotor"]
        if all(np.ndim(v) == 0 for v in (i_d, i_q, rpm)):
            return {k: float(v) for k, v in region.items()}
        return region


_MAP: IronLossMap | None = None


def iron_loss(i_d, i_q, rpm, stack_mm: float = STACK_MM_MODEL, *, turns: int = 10, paths: int = 1, lamination_mm: float = D_LAM_CAD * 1e3,
              extrapolate: str = "raise") -> dict:
    """Iron loss [W] of teeth, yoke, rotor and their total at the current (i_d, i_q) [A peak] and speed rpm, for a stack of stack_mm, the winding of
    `turns` per coil side in `paths` parallel paths and laminations of lamination_mm (the CAD's 0.3 mm by default; the datasheet's 0.356 mm on request): see
    the module docstring. Loads out/optimization/maps/chosen/iron_map.npz on first use."""
    global _MAP
    if _MAP is None:
        _MAP = IronLossMap()
    return _MAP(i_d, i_q, rpm, stack_mm, turns=turns, paths=paths, lamination_mm=lamination_mm, extrapolate=extrapolate)


# ---------------------------------------------------------------------------------------------- the command line: fit, check, write
def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    f, B, P = datasheet_points()
    print(f"the datasheet: {len(P)} points, B {B.min()} ... {B.max()} T, f {f.min():.0f} ... {f.max():.0f} Hz, gauge {D_LAM * 1e3:.3f} mm, "
          f"sigma {SIGMA / 1e6:.2f} MS/m, rho_m {RHO_STEEL:.0f} kg/m^3")
    print(f"classical coefficient k_e = pi^2 sigma d^2 / (6 rho_m) = {K_E:.5e} W/kg per Hz^2 T^2; Bertotti's constant C_x = {C_X:.4f}")
    fit = fit_bertotti(f, B, P)
    fit_free = fit_bertotti(f, B, P, free_ke=True)
    err = predict(fit, f, B) / P - 1.0
    print(f"fit (k_e classical): k_h = {fit['kh']:.5f} W/kg/(Hz T^alpha), alpha = {fit['alpha']:.4f}, k_x = {fit['kx']:.5e} W/kg/(Hz^1.5 T^1.5)")
    print(f"   relative error on the {len(P)} fitted points: rms {100 * fit['rms']:.2f} %, largest {100 * fit['max']:.1f} %")
    print(f"fit with k_e free: k_e = {fit_free['ke']:.5e} = {fit_free['ke'] / K_E:.3f} x the classical one "
          f"(k_h {fit_free['kh']:.5f}, alpha {fit_free['alpha']:.4f}, k_x {fit_free['kx']:.5e}; rms {100 * fit_free['rms']:.2f} %)")
    for fr, b in ((60, 1.5), (200, 1.0), (1000, 1.0), (2000, 1.0)):
        parts = np.array([fit["kh"] * fr * b ** fit["alpha"], K_E * fr ** 2 * b ** 2, fit["kx"] * fr ** 1.5 * b ** 1.5])
        print(f"   shares at {fr} Hz, {b} T: hysteresis {100 * parts[0] / parts.sum():.0f} %, classical {100 * parts[1] / parts.sum():.0f} %, excess {100 * parts[2] / parts.sum():.0f} %")
    cv = cross_validate()
    s = cv["interior"]
    print(f"cross-validation, interior folds ({s['n_interior_points_tested']} left-out points): rms {100 * s['rms']:.2f} % (limit {100 * TOL_RMS:.0f} %), "
          f"largest {100 * s['max']:.1f} % (limit {100 * TOL_MAX:.0f} %): {'PASS' if s['pass'] else 'FAIL'}")
    for row in cv["folds"]:
        kind = "interior" if row["interior"] else "edge/extrapolation (not judged)"
        print(f"   leave out {row['fold']:>18s} = {str(row['left_out']):>12s} ({row['n']:2d} points): rms {100 * row['rms']:5.2f} %, largest {100 * row['max']:5.1f} %   [{kind}]")
    skin = eddy_skin_factor(D_LAM, [2000.0, 20000.0], mu_r=[2200.0, 1000.0])
    print(f"skin effect on the classical term (exact / classical): {skin[0]:.3f} at 2 kHz (mu_r 2,200), {skin[1]:.3f} at 20 kHz (mu_r 1,000)")
    k_cad = eddy_scale(D_LAM_CAD * 1e3)
    print(f"the study's lamination is the CAD's {D_LAM_CAD * 1e3:.3f} mm: the classical term is scaled by ({D_LAM_CAD * 1e3:.3f} / {D_LAM * 1e3:.3f})^2 = {k_cad:.4f} "
          f"(k_e {K_E:.4e} -> {K_E * k_cad:.4e}); hysteresis and excess are left unscaled (their thickness dependence is weaker and one gauge cannot pin it)")
    for fr, b in ((200, 1.0), (1000, 1.0)):
        p356 = bertotti(fr, b, fit["kh"], fit["alpha"], fit["kx"], K_E)
        p300 = bertotti(fr, b, fit["kh"], fit["alpha"], fit["kx"], K_E * k_cad)
        print(f"   a sinusoid of {b} T at {fr} Hz: {float(p356):.2f} W/kg at 0.356 mm, {float(p300):.2f} W/kg at 0.300 mm ({100 * (float(p300) / float(p356) - 1):+.1f} %)")
    print("self-test of the waveform operators:")
    ok = self_test(fit)
    with open(FIT_JSON, "w", encoding="utf-8") as fh:
        json.dump({"fit": fit, "fit_free_ke": fit_free, "cross_validation": cv, "k_e_classical": K_E, "C_x": C_X, "gauge_mm": D_LAM * 1e3, "sigma_S_per_m": SIGMA,
                   "lamination_study_mm": D_LAM_CAD * 1e3, "eddy_scale_to_study_lamination": k_cad, "k_e_study_lamination": K_E * k_cad,
                   "rho_steel": RHO_STEEL, "lam_fill": LAM_FILL, "tolerances": {"rms": TOL_RMS, "max": TOL_MAX}, "self_test_pass": ok,
                   "source": "AK Steel Di-Max M-19 29 gauge, EMERF Lamination Steels 3rd ed. (MIT 6.685 excerpt, 2008), p. 5, W/kg"}, fh, indent=1)
    with open(OUT / "iron_fit_points.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["f_Hz", "B_T", "datasheet_W_per_kg", "model_W_per_kg", "relative_error"])
        for fi, bi, pi, ei in zip(f, B, P, err):
            w.writerow([fi, bi, pi, pi * (1 + ei), ei])
    print(f"wrote {FIT_JSON} and {OUT / 'iron_fit_points.csv'}")
    return 0 if (s["pass"] and ok) else 1


if __name__ == "__main__":
    sys.exit(main())
