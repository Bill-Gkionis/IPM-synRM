"""optimization.ac_pwm: the ripple current of a two-level three-phase inverter with space-vector PWM, from the inverter's volt-seconds to the phase current's rms, its peak-to-
peak and its harmonics. No FEMM; the winding's AC resistance at the switching frequency (optimization/ac_slot.py) turns the rms into watts.

What it is, and why
    A winding is not driven by a sinusoid: the inverter switches its three poles between the two DC rails (0 and V_dc) at f_sw, so the phase voltage is a train of pulses whose
    AVERAGE is the sinusoid the controller asked for, and the current has a ripple on top of its sinusoid. The ripple is small against the rated current, but it lives at
    f_sw (10 or 20 kHz, 10 to 100 times the fundamental), where the strands' AC resistance is several times the DC one (ac_slot.py: R_ac / R_dc = 4.3 at 20 kHz for 0.315 mm
    strands, 4.2 to 4.4 over three random layouts of the slot). So the ripple's loss is the ripple's rms squared times a resistance that grows as the square of the frequency, and the ripple itself grows as the winding's
    inductance falls: a winding of few turns has little inductance (L goes as the series turns squared) and a large ripple. This module makes that chain a number.

The method (each step built from the one before)
    1. One switching period T = 1 / f_sw, the three pole voltages v_x0 (x = a, b, c) take only the values 0 and V_dc. With the load's neutral isolated and the load balanced, the
       phase voltage across the phase inductance is v_xn = v_x0 - v_n0, v_n0 = (v_a0 + v_b0 + v_c0) / 3 (the neutral's voltage, the same common-mode term for all three).
    2. Over one period the controller wants the sinusoidal reference v*_x. Its average over the period is the average of v_xn: with the duty cycle d_x of pole x (the fraction of
       the period the pole is high), <v_x0> = d_x V_dc and <v_xn> = V_dc (d_x - (d_a + d_b + d_c) / 3) = v*_x. The load's EMF is the slowly varying reference.
    3. Symmetric ('centred') pulses, what a triangular carrier gives and the seven-segment space-vector sequence 000-X-XY-111-XY-X-000: pole x is high for |tau - 1/2| < d_x / 2
       (tau = t / T). The duty cycles of SVPWM are the sinusoid plus the common-mode term that centres the three references, d_x = 1/2 + (v*_x - (max v* + min v*) / 2) / V_dc: the
       reference's peak is m V_dc / sqrt(3), m the MODULATION INDEX (m = 1 is the largest sinusoid SVPWM can make, V_dc / sqrt(3); the line to line voltage is then V_dc).
    4. The ripple: L di_x/dt = v_xn - v*_x on top of the sinusoid, so i_x(tau) = (V_dc T / L) integral of (v_xn - v*_x) / V_dc, with the constant chosen so that the ripple
       averages to zero over the period (the sinusoid carries the rest). v_xn - v*_x is piecewise constant (at most seven pieces: six breakpoints), so i_x is piecewise linear:
       its mean square, its peak-to-peak and its Fourier coefficients have exact closed forms (ripple_period: the segment integral (i_k^2 + i_k i_{k+1} + i_{k+1}^2) dt / 3, and
       for the coefficients the second derivative of a piecewise linear function is a sum of deltas, so c_h = -sum_j (s_j - s_{j-1}) exp(-2 pi i h tau_j) / (2 pi h)^2).
    5. Average over the fundamental cycle: the duty cycles follow theta = the angle of the reference vector. Over theta the mean square of the ripple is
           I_ripple,rms^2 = (V_dc / (L f_sw))^2 h(m)^2,         h(m) = (1/2pi) integral over theta of the period's mean square
       in units of V_dc T / L: h depends on m ONLY (and not on the load's power factor: the ripple is set by the volt-second errors, and the EMF is the reference). The
       fraction of the ripple's power at the carrier's harmonics h = 1, 2, 3, ... is the same average of 2 |c_h|^2 over the period's mean square; it decides which AC
       resistance applies (R_ac / R_dc grows about as f^2), so a spectrum split across f_sw, 2 f_sw, 3 f_sw ... weighs the loss up.
    6. The tables of h(m), the peak-to-peak and the harmonic split are made once (python -m optimization.ac_pwm) and saved (out/optimization/copper_ac/pwm/ripple_table.json); a
       fast function (ripple_rms) interpolates them, in microseconds. L is the winding's phase inductance in the stationary frame as the ripple sees it: for a machine with
       L_d < L_q the saliency modulates it with the rotor position (L between L_d and L_q); the rated-current values of the 2-D model are L_d = 1.390 mH and L_q = 2.112 mH at 10
       turns a coil side, scaled with the series turns squared. Taking L_d gives the largest ripple (a bound). Not included: the end-winding leakage (the 2-D FEMM has none: it
       lowers the ripple), the incremental inductance of saturated iron (it is lower than the apparent one: it raises the ripple), dead time, device drops, the phase resistance
       (R << w L at f_sw: no effect), the common-mode and bearing currents.

Checks (python -m optimization.ac_pwm runs them, and exits 0 if all pass)
    * Parseval: the sum of the harmonics' power equals the period's mean square (closed forms against each other).
    * One leg alone (a half bridge into an inductor at the DC midpoint, the textbook case): the peak-to-peak ripple is V_dc d (1 - d) / (L f_sw), exactly the function's result for that load.
    * An independent time-domain simulation: the three poles are driven by comparing the CONTINUOUS reference (natural sampling) with a triangular carrier, the load's ODE is stepped
      on a fine grid over a whole fundamental cycle at f_sw / f_e = 360, and the ripple is what is left after the period's moving average is taken off. It must agree with h(m).
    * The limit m -> 0: the voltage errors vanish, so the ripple does (h(0) = 0).
    The method follows the standard analysis of the ripple of carrier-based and space-vector PWM; for the literature see D. Jiang, F. Wang, 'Current-Ripple Prediction for Three-Phase
    PWM Converters', IEEE Trans. Ind. Appl. 50(1), 2014 (DOI 10.1109/TIA.2013.2270224), G. Grandi, J. Loncarski, 'Evaluation of Current Ripple Amplitude in Three-Phase PWM Voltage
    Source Inverters' (Univ. of Bologna), and D. G. Holmes, T. A. Lipo, 'Pulse Width Modulation for Power Converters', IEEE Press / Wiley 2003. Those three were found by a web search and
    not read in full: CHECK BEFORE CITING (the formulas of this module are derived and checked here, not copied from them).

Run
    uv run python -m optimization.ac_pwm        (no FEMM; about 5 s) prints the checks and the table, writes out/optimization/copper_ac/pwm/ripple_table.json and fig_pwm_ripple.png.

Numbers to expect (the checks print them; units of the ripple: V_dc / (L f_sw)): h(m) is 0 at m = 0 and rises steadily to 0.02515 at m = 1; the largest peak-to-peak over the cycle
    is 0.0477 at m = 0.2, 0.1000 at m = 0.6 and exactly 1/6 at m = 1. Where the ripple's power sits: at m = 0.2 78 % of it is at 2 f_sw and 1.5 % at f_sw; at m = 0.6 it is 27 % at f_sw
    and 64 % at 2 f_sw; at m = 1 it is 90 % at f_sw. For the example winding (26 turns a coil side, 4 paths: 78 series turns a phase, L = 1.390 mH x (78 / 120)^2 = 0.587 mH) at
    V_dc = 48 V and f_sw = 20 kHz the unit is 4.09 A, so the ripple is 0.103 A rms (0.68 A peak-to-peak at most) at m = 1.

Self-check (the answers follow).
  1. Why does the ripple vanish as the reference goes to zero, although the pole voltages still jump between 0 and V_dc?
  2. A winding of N turns and a parallel paths is rewound with 2 N turns and 2 a paths (the same current). What happens to the ripple? And with 2 N turns and the same a?
  3. The ripple current goes as 1 / f_sw and the strands' proximity loss as f^2 while the strands are thin against the skin depth. Why does a higher switching frequency NOT lower the
     proximity part of the ripple loss, only its DC-like part?
Answers.
  1. At m = 0 the controller alternates the two zero vectors, 000 and 111, for equal times: all three poles move together, v_n0 follows them, and v_xn = v_x0 - v_n0 is zero: no
     voltage across the phase inductance, no ripple. The ripple is the volt-second error of the difference between the poles, not the poles' swing.
  2. The inductance goes as the series turns squared, n = 12 N / a, and the ripple as 1 / L. 2 N and 2 a keep n, so L and the ripple do not change. 2 N alone quadruples L: the ripple falls
     to a quarter (the winding's current at the same ampere-turns halves, so the ripple relative to it halves; the winding needs twice the volts, which the corner speed pays for).
  3. The ripple's rms squared falls as 1 / f_sw^2 and the proximity part of the resistance (R_ac - R_dc) rises as f_sw^2: the product does not depend on f_sw. Only the DC-like part,
     I_ripple^2 R_dc, falls as 1 / f_sw^2. (At d / delta near 1 and above the resistance rises more slowly than f^2, and a higher f_sw does lower the whole loss.)
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import numpy as np

import paths

N_HARM = 8                                   # harmonics of the carrier kept in the spectrum split (the rest is returned as 'tail')


# ---------------------------------------------------------------------------------------------- the modulator and one period's exact ripple
def duty_cycles(theta, m):
    """The pole duty cycles (n, 3) of symmetric space-vector PWM for the reference angle theta [rad] (array) at modulation index m: d_x = 1/2 + (v_x - (max v + min v) / 2) with
    v_x = (m / sqrt 3) cos(theta - 2 pi k / 3) in units of V_dc (so that m = 1 reaches 0 and 1 at the peaks: V_dc / sqrt 3 is the largest sinusoid)."""
    theta = np.atleast_1d(np.asarray(theta, dtype=float))
    v = (m / math.sqrt(3.0)) * np.cos(theta[:, None] - 2.0 * math.pi * np.arange(3)[None, :] / 3.0)
    return 0.5 + v - 0.5 * (v.max(axis=1) + v.min(axis=1))[:, None]


def ripple_period(d: np.ndarray, phase: int = 0, n_harm: int = N_HARM) -> dict:
    """Exact ripple of phase `phase` over ONE switching period for the duty cycles d (n, 3) with centred pulses, in units of V_dc T / L (so the inverter's V_dc = 1, T = 1, L = 1).
    Returns {'ms' (n,): the mean square of the ripple, 'pp' (n,): its peak-to-peak, 'harm' (n, n_harm): the power 2 |c_h|^2 of the sinusoid pair at h = 1 ... n_harm}.
    Method (module docstring, steps 3 and 4): the six edges of the three pulses cut the period into seven pieces in each of which the voltage across the phase inductance is constant;
    the current is then piecewise linear."""
    d = np.atleast_2d(np.asarray(d, dtype=float))
    n = len(d)
    edges = np.concatenate([(1.0 - d) / 2.0, (1.0 + d) / 2.0], axis=1)           # (n, 6): where the pulses start and stop
    tau = np.concatenate([np.zeros((n, 1)), np.sort(edges, axis=1), np.ones((n, 1))], axis=1)     # (n, 8) the breakpoints
    dtau = np.diff(tau, axis=1)                                                  # (n, 7) the pieces
    mid = 0.5 * (tau[:, :-1] + tau[:, 1:])                                       # (n, 7)
    high = ((1.0 - d)[:, None, :] / 2.0 < mid[:, :, None]) & (mid[:, :, None] < (1.0 + d)[:, None, :] / 2.0)       # (n, 7, 3): is pole x high in piece k?
    v_x0 = high.astype(float)
    v_xn = v_x0 - v_x0.mean(axis=2, keepdims=True)                               # phase voltage over the neutral, per unit of V_dc
    ref = d - d.mean(axis=1, keepdims=True)                                      # the average of v_xn over the period: the reference
    u = (v_xn - ref[:, None, :])[:, :, phase]                                    # (n, 7): the voltage error in each piece: the slope of the current
    i0 = np.concatenate([np.zeros((n, 1)), np.cumsum(u * dtau, axis=1)], axis=1)    # (n, 8) the current at the breakpoints, up to a constant
    mean0 = np.sum(dtau * 0.5 * (i0[:, :-1] + i0[:, 1:]), axis=1)
    i = i0 - mean0[:, None]                                                      # the ripple: zero average over the period
    ms = np.sum(dtau * (i[:, :-1] ** 2 + i[:, :-1] * i[:, 1:] + i[:, 1:] ** 2) / 3.0, axis=1)
    pp = i.max(axis=1) - i.min(axis=1)
    # Fourier coefficients of a continuous piecewise linear function from the jumps of its slope at the breakpoints
    jumps = np.concatenate([(u[:, 0] - u[:, -1])[:, None], np.diff(u, axis=1)], axis=1)     # (n, 7): the slope's jump at tau_0 = 0, tau_1, ... tau_6
    tj = tau[:, :-1]                                                              # (n, 7): the jump positions
    h = np.arange(1, n_harm + 1)
    phase_term = np.exp(-2j * math.pi * h[None, None, :] * tj[:, :, None])         # (n, 7, H)
    c = -np.sum(jumps[:, :, None] * phase_term, axis=1) / (2.0 * math.pi * h[None, :]) ** 2    # (n, H)
    return {"ms": ms, "pp": pp, "harm": 2.0 * np.abs(c) ** 2}


def ripple_vs_m(m_grid, n_theta: int = 3600, n_harm: int = N_HARM) -> dict:
    """The ripple over a fundamental cycle for each modulation index of m_grid: {'m', 'h_rms' (the rms of the ripple in units of V_dc / (L f_sw)), 'pp_max' (the largest
    peak-to-peak over the cycle, same unit), 'pp_mean', 'frac' (len(m), n_harm): the share of the ripple's power at each carrier harmonic, 'tail': the share above n_harm}."""
    theta = 2.0 * math.pi * (np.arange(n_theta) + 0.5) / n_theta                   # period centres: no sample sits on a sector edge of the modulator
    out = {"m": [], "h_rms": [], "pp_max": [], "pp_mean": [], "frac": [], "tail": []}
    for m in np.atleast_1d(m_grid):
        r = ripple_period(duty_cycles(theta, float(m)), 0, n_harm)
        ms_mean = float(np.mean(r["ms"]))
        harm = r["harm"].mean(axis=0)
        out["m"].append(float(m))
        out["h_rms"].append(math.sqrt(ms_mean))
        out["pp_max"].append(float(r["pp"].max()))
        out["pp_mean"].append(float(r["pp"].mean()))
        out["frac"].append((harm / ms_mean).tolist() if ms_mean > 0 else [0.0] * n_harm)
        out["tail"].append(float(1.0 - harm.sum() / ms_mean) if ms_mean > 0 else 0.0)
    return out


# ---------------------------------------------------------------------------------------------- independent time-domain simulation (the check)
def simulate_rms(m: float, ratio: int = 360, steps: int = 4000) -> float:
    """The ripple's rms [in units of V_dc / (L f_sw)] from a brute-force simulation: the three poles follow the comparison of the CONTINUOUS min-max-injected reference with a
    triangular carrier (natural sampling, not the frozen duty cycles of ripple_period), the phase-a current is the integral of (v_an - reference) on a grid of `steps` points a
    switching period over one whole fundamental cycle of `ratio` switching periods, and the ripple is the current minus its moving average over one period."""
    n = ratio * steps
    t = (np.arange(n) + 0.5) / steps                                               # in switching periods
    theta = 2.0 * math.pi * t / ratio
    v = (m / math.sqrt(3.0)) * np.cos(theta[:, None] - 2.0 * math.pi * np.arange(3)[None, :] / 3.0)
    dref = 0.5 + v - 0.5 * (v.max(axis=1) + v.min(axis=1))[:, None]               # the reference duty of each pole at each instant
    carrier = np.abs(2.0 * (t % 1.0) - 1.0)                                       # triangular, 1 at the period's edges, 0 at its centre
    high = dref > carrier[:, None]
    v_x0 = high.astype(float)
    v_an = v_x0[:, 0] - v_x0.mean(axis=1)
    ref_a = v[:, 0]                                                                # the reference phase voltage (the common-mode term cancels in v_xn)
    i = np.cumsum(v_an - ref_a) / steps                                            # integral over time in switching periods (V_dc T / L units)
    kernel = np.ones(steps) / steps
    pad = np.concatenate([i[-steps:], i, i[:steps]])                               # a periodic moving average over one switching period
    avg = np.convolve(pad, kernel, mode="same")[steps:-steps]
    ripple = i - avg
    return float(np.sqrt(np.mean(ripple ** 2)))


# ---------------------------------------------------------------------------------------------- the fast function used by copper_loss.py
_TABLE: dict | None = None
TABLE_PATH = paths.path("out") / "optimization" / "copper_ac" / "pwm" / "ripple_table.json"


def load_table(path: str | Path | None = None) -> dict:
    """The ripple table (written by main()); cached. The packaged copy for copper_loss.py is copper_loss_tables.json."""
    global _TABLE
    if path is not None or _TABLE is None:
        with open(path or TABLE_PATH, encoding="utf-8") as handle:
            tab = json.load(handle)
        tab["_m"] = np.array(tab["m"])
        tab["_h"] = np.array(tab["h_rms"])
        tab["_frac"] = np.array(tab["frac"])
        if path is not None:
            return tab
        _TABLE = tab
    return _TABLE


def modulation_index(v_peak: float, v_dc: float) -> float:
    """m = v_peak / (V_dc / sqrt 3): the reference phase voltage's PEAK over the largest sinusoid space-vector PWM makes (m = 1)."""
    return v_peak * math.sqrt(3.0) / v_dc


def ripple_rms(v_dc: float, f_sw: float, L: float, m: float, table: dict | None = None) -> float:
    """The phase current's ripple rms [A] at the DC bus v_dc [V], switching frequency f_sw [Hz], phase inductance L [H] and modulation index m (0 ... 1, linear range): h(m) V_dc / (L f_sw)
    with h interpolated from the table. m above 1 (overmodulation) is not covered: the value at m = 1 is returned (the caller checks m)."""
    tab = table or load_table()
    return float(np.interp(min(max(m, 0.0), 1.0), tab["_m"], tab["_h"]) * v_dc / (L * f_sw))


# ---------------------------------------------------------------------------------------------- checks, the table and the picture
def _check(name: str, value: float, expect: float, tol: float, results: list, relative: bool = True) -> None:
    err = abs(value - expect) / (abs(expect) if relative and expect != 0 else 1.0)
    ok = err <= tol
    results.append(ok)
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {value:.8g} (expected {expect:.8g}, error {err:.2e}, tolerance {tol:g})")


def self_checks() -> bool:
    out: list[bool] = []
    theta = 2.0 * math.pi * (np.arange(720) + 0.5) / 720
    for m in (0.3, 0.7, 1.0):
        r = ripple_period(duty_cycles(theta, m), 0, 200)
        _check(f"Parseval at m = {m}: sum of the harmonics' power over the period's mean square", float(r["harm"].sum(axis=1).mean()), float(r["ms"].mean()), 2e-6, out)
    # one leg alone: v_xn = v_x0 - V_dc / 2 against the DC midpoint: build it by hand with the same piecewise integration
    d = 0.3
    tau = np.array([0.0, (1 - d) / 2, (1 + d) / 2, 1.0])
    u = np.array([-d, 1.0 - d, -d])                                               # v - mean: the pole voltage minus its average d
    i0 = np.concatenate([[0.0], np.cumsum(u * np.diff(tau))])
    pp_leg = i0.max() - i0.min()
    _check("one leg, duty 0.3: peak-to-peak against V_dc d (1 - d) / (L f_sw)", float(pp_leg), d * (1 - d), 1e-12, out)
    _check("h(0) = 0 (no ripple without a reference)", ripple_vs_m([0.0], 360)["h_rms"][0], 0.0, 1e-12, out, relative=False)
    for m in (0.3, 0.6, 0.9, 1.0):
        analytic = ripple_vs_m([m], 720)["h_rms"][0]
        sim = simulate_rms(m)
        _check(f"time-domain simulation at m = {m}", sim, analytic, 2e-3, out)
    return all(out)


def figure(tab: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    surface, ink, ink2, grid, axis = "#fcfcfb", "#0b0b0b", "#52514e", "#e1e0d9", "#c3c2b7"
    blue, orange, aqua, gray = "#2a78d6", "#eb6834", "#1baf7a", "#898781"
    plt.rcParams["font.family"] = ["Segoe UI", "DejaVu Sans"]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.6), facecolor=surface)
    m = np.array(tab["m"])
    ax = axes[0]
    ax.plot(m, np.array(tab["h_rms"]) * 1e3, color=blue, lw=2.0, label="rms of the ripple (the h(m) of the text)")
    ax.plot(m, np.array(tab["pp_max"]) * 1e3, color=orange, lw=2.0, label="largest peak-to-peak over the cycle")
    ax.set_xlabel("modulation index m  (1 = the largest sinusoid of space-vector PWM, V_dc / sqrt 3 peak)", fontsize=8.5, color=ink2)
    ax.set_ylabel("phase-current ripple, in thousandths of V_dc / (L f_sw)", fontsize=8.5, color=ink2)
    ax.set_title("Ripple of a three-phase two-level inverter, space-vector PWM", fontsize=10, color=ink, loc="left", fontweight="bold")
    ax = axes[1]
    frac = np.array(tab["frac"])
    for h, col in zip((1, 2, 3, 4), (blue, orange, aqua, gray)):
        ax.plot(m, 100 * frac[:, h - 1], color=col, lw=2.0, label=f"at {h} x f_sw")
    ax.plot(m, 100 * np.array(tab["tail"]), color=ink2, lw=1.0, ls=":", label="above 8 x f_sw")
    ax.set_xlabel("modulation index m", fontsize=8.5, color=ink2)
    ax.set_ylabel("share of the ripple's power [%]", fontsize=8.5, color=ink2)
    ax.set_title("Where the ripple's power sits in frequency", fontsize=10, color=ink, loc="left", fontweight="bold")
    for ax in axes:
        ax.set_facecolor(surface)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(axis)
        ax.tick_params(colors=ink2, labelsize=8.5)
        ax.grid(True, color=grid, lw=0.75)
        ax.set_axisbelow(True)
        leg = ax.legend(frameon=False, fontsize=8.5, loc="best")
        for t in leg.get_texts():
            t.set_color(ink2)
    fig.text(0.012, 0.01, "Exact piecewise-linear ripple of one phase, centred pulses, balanced load with an isolated neutral, averaged over a fundamental cycle (3,600 periods). "
             "Source: optimization/ac_pwm.py.", fontsize=7.5, color=ink2)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(path, dpi=130, facecolor=surface)
    plt.close(fig)


def main() -> int:
    t0 = time.perf_counter()
    print("self-checks:")
    ok = self_checks()
    m_grid = np.round(np.arange(0.0, 1.0001, 0.02), 4)
    tab = ripple_vs_m(m_grid, 3600, N_HARM)
    tab["note"] = ("h_rms: ripple rms in units of V_dc/(L f_sw); pp_max: the largest peak-to-peak over the fundamental cycle, same unit; frac[m][h-1]: share of the ripple's power at "
                   "h x f_sw; tail: above N_HARM. Two-level, symmetric SVPWM, isolated neutral, balanced inductive load. Made by optimization/ac_pwm.py.")
    tab["n_theta"], tab["n_harm"] = 3600, N_HARM
    TABLE_PATH.parent.mkdir(parents=True, exist_ok=True)
    TABLE_PATH.write_text(json.dumps(tab, indent=1), encoding="utf-8")
    figure(tab, TABLE_PATH.parent / "fig_pwm_ripple.png")
    h = np.array(tab["h_rms"])
    k = int(np.argmax(h))
    print(f"\nh(m): maximum {h[k]:.5f} at m = {tab['m'][k]:.2f}; at m = 1: {h[-1]:.5f}; largest peak-to-peak over the cycle {max(tab['pp_max']):.5f} (units V_dc / (L f_sw)) at m = {tab['m'][int(np.argmax(tab['pp_max']))]:.2f}")
    for mm in (0.2, 0.4, 0.6, 0.8, 1.0):
        j = int(round(mm / 0.02))
        print(f"  m = {mm:.1f}: h = {h[j]:.5f}, pp max {tab['pp_max'][j]:.5f}, power at f_sw {100*tab['frac'][j][0]:.1f} %, 2 f_sw {100*tab['frac'][j][1]:.1f} %, 3 f_sw {100*tab['frac'][j][2]:.1f} %, "
              f"above 8 f_sw {100*tab['tail'][j]:.2f} %")
    print(f"written {TABLE_PATH} ({time.perf_counter() - t0:.1f} s)")
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
