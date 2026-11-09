"""optimization.copper_loss: the copper loss of any winding of the optimized stator, at any operating point, copper temperature and stack length: one fast function for the optimizer.

What it is, and why
    The winding study (optimization/winding.py) priced a winding by its DC loss at 100 degC and a low-frequency proximity formula at 3,000 rpm, and found the AC loss negligible. This module
    prices it properly, from the standstill to the stress limit (16,330 rpm, 1,090 Hz) and with the inverter's PWM ripple, in six parts, each with its own physics and each checked against FEMM
    in the studies of optimization/ac_*.py (their docstrings carry the numbers):
      p_dc_slot     the DC loss in the slots' active length                       1.5 R_slot i^2
      p_dc_end      the DC loss in the end turns                                  1.5 R_end i^2
      p_skin        the skin effect of the fundamental current in its own strand   (P_dc,slot + P_dc,end) (R_ac / R_dc of a lone wire - 1)
      p_prox        the proximity loss: the machine's own field (rotor and armature, its harmonics) acting on the strands of the slots, from FEMM magnetostatic solves of the sector
                    at the operating points (optimization/ac_machine.py), with the lone wire's exact high-frequency factor
      p_pwm         the loss of the inverter's PWM ripple current: its rms (optimization/ac_pwm.py) times the AC resistance of the slot at the switching frequency and its harmonics, which
                    the slot model of FEMM gives strand by strand (optimization/ac_slot.py) as the extra loss K (d k / delta^2)^2 over DC, with the iron's screening
      p_circ        the circulating currents between the strands of a turn when they are connected in parallel and are not transposed ('arrangement' of the winding): zero for a
                    transposed winding (the winding study's assumption), and for each placement of the strands the loss the FEMM slot model finds, scaled to any winding
    all in watts for the stack length given. The currents are PEAK phase currents (the amplitude of the fundamental) everywhere, as in winding.py: the three phases dissipate
    3 R (i / sqrt 2)^2 = 1.5 R i^2.

The call
    losses(winding, op, t_copper_c=100.0, stack_mm=20.0) -> dict      (watts; the full list of keys is in the function's docstring)
    winding   {'N': turns per coil side, 'a': parallel paths (1, 2 or 4), 'n_s': strands in hand, 'd_mm': bare strand diameter, 'd_over_mm': overall (enamelled) diameter,
               optional 'arrangement' ('transposed' default, 'row', 'cluster', 'random', 'line'), 'end_model' ('diamond' default or 'arc'), 'end_ext_mm', 'end_angle_deg'}
    op        {'i_peak': peak phase current [A], 'gamma_deg': current angle from the q axis toward -d, 'rpm': mechanical speed, 'f_sw': switching frequency [Hz],
               'v_dc': DC bus [V], optional 'ripple_layers' ('same' (default: the worse), 'diff' or 'mixed')}
    t_copper_c   the copper temperature [degC]: the resistivity follows it (rho_20 (1 + alpha_20 (T - 20))) in every part, with no fixed 100 degC; the skin depth, and so the strands'
                 AC resistance, follows it too
    stack_mm     the stack length [mm]
    losses_arrays(...) is the same function on arrays (numpy broadcasting), for an optimizer that evaluates thousands of windings at once.

The end turns (the improvement on the first study's stand-in, and what is still assumed)
    The first study took an end turn as a half circle on the coil's span at the slot's mean radius, plus a 2 mm straight leg at each end of it (60.4 mm). This module draws it as what a
    double-layer lap coil is, a DIAMOND: from the exit of one coil side (in the radial middle of its layer) the conductor runs straight out of the stack for e_ext (3 mm by default: past the
    lamination edge and the slot liner, to where the bend starts), then in two straight legs rising at the nose angle alpha (30 degrees by default) to the apex and down to the other coil side
    (in the middle of the OTHER layer: the two coil sides of a coil are in different layers, at radii 50.08 and 56.30 mm, so the in-plane chord between them is
    c = sqrt(r_0^2 + r_1^2 - 2 r_0 r_1 cos(40 deg)) = 36.85 mm for the pitch of 4 slots); one end turn is
        l_end = 2 e_ext + c / cos(alpha) = 6 + 42.55 = 48.55 mm,       the turn  l_mt = 2 L_stack + 2 l_end = 137.1 mm for the 20 mm stack.
    That is 20 % shorter than the first study's 61.15 mm (the half circle is the longest sensible end turn: its rise is half the chord). The rise (10.6 mm at 30 degrees) is the design's choice
    and not known until the coil is formed: the nose angle and the extension are knobs of the winding ('end_angle_deg', 'end_ext_mm'); at 20 degrees l_end is 45.2 mm and at 40 degrees 54.1 mm.
    end_model = 'arc' reproduces the first study's number (the regression test of tests/test_copper_loss.py). Not modelled: the end turns' own AC loss (the field outside the iron is weak and
    not computed here), the radial build of the end winding's bundles, the connections between coils, and the leakage of the end winding in the phase inductance.

What is, and is not, in the model (limits that a user must know)
    * The machine's field comes from the magnetostatic solves of the sector: the eddy currents of the rotor's and stator's laminations and of the magnets, which screen the high harmonics (the slot-
      passing field at 9 x f_e is 9.8 kHz at 16,330 rpm), are not in it: p_prox at high speed is an upper bound. The slot's own field at the switching frequency (p_pwm) is solved with laminated stator
      steel (FEMM's complex permeability of 0.356 mm sheets) and its loss changes by -11 % ... +16 % with the steel model (sens stage: lamination thickness, permeability, rotor surface).
    * The fit of the machine's field (ac_machine.fit_split: the no-load part plus AT^2 times a slowly varying total shared among the harmonic groups) holds over the sampled region (0 ... 3 x the
      rated ampere-turns, the current angle 0 ... 70 deg), with leave-one-out errors for the total of 5.7 % in the median and -19.5 % / +24.4 % at worst; outside the region the call clips to it
      and says so ('extrapolated' in the result).
    * The ripple uses the d-axis inductance (the smaller: the larger ripple) scaled with the series turns squared and the stack; the modulation index is from the winding's own voltage need
      (resistance plus the flux linkage of the operating point); above 1 the inverter is overmodulated, the ripple is taken at 1 and the result says so ('overmodulated').
    * The layers of the slot are 'same' phase (both layers of one phase: 6 of the 9 slots of the sector) or 'diff' (two phases 60 degrees apart: the other 3, with 0.82 of the extra loss) or
      'mixed' (the mean over the 9 slots); the ripple's loss defaults to the worse, 'same'.

Run
    uv run python -m optimization.copper_loss          prints the interface, the check against the first study and the speed per call (the tables are loaded from copper_loss_tables.json)
    uv run python -m optimization.copper_loss --build  makes copper_loss_tables.json from the studies' outputs (out/optimization/copper_ac/...): run it after the studies
    uv run pytest tests/test_copper_loss.py -q         the tests (no FEMM)

Numbers to expect (printed by the run, 2 s): the DC loss agrees with winding.py to 1e-15 with end_model = 'arc' (40.841848 W for the first study's point); the winding study's example (26 turns x 8 x
0.315 mm, 4 paths) at 118.6 ampere-turns, 3,000 rpm, 20 kHz, 48 V, 100 degC, 20 mm stack: p_dc_slot 14.16 W, p_dc_end 34.37 W (diamond end turns), p_prox 0.0070 W, p_skin 1e-6 W,
p_pwm 0.0080 W (ripple 0.094 A rms, modulation index 0.916), p_circ 0 (transposed), total 48.54 W; the speed is 189 microseconds a call of losses() and 11 microseconds a winding in an
array of 100,000 (losses_arrays), on this PC.

Self-check (the answers follow).
  1. Why is p_skin so small at 1,090 Hz, and p_prox not zero at standstill-like speeds?
  2. A winding of the same N / a, ampere-turns and fill is made of strands twice as thick. Which parts change, by how much?
  3. What does 'arrangement' do to a winding with n_s = 1?
Answers.
  1. The skin effect of a thin wire is (d / delta)^4 / 768: 1e-5 at d / delta = 0.3 and 5e-7 at the 0.137 of W1 at 1,090 Hz (the tables of copper_loss_study.py: 0.00002 W of 48.7 W). The
     proximity loss goes as the speed squared, so it is not zero at low speed, only small: it is exactly zero only at standstill (f_e = 0).
  2. The DC parts do not change (same fill, same copper area); p_prox and p_skin's field-driven part grow as d^2 (four times); p_pwm grows with the slot's extra loss, K (d k / delta^2)^2, also as
     d^2 at low d / delta; p_circ does not depend on d at a given fill (it depends on the number of strands in hand and the field's variation across the slot).
  3. Nothing: a single strand has nobody to share a turn's current with, so there is no circulating current ((1 - 1 / n_s) = 0).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.special import ive

TABLES_PATH = Path(__file__).with_name("copper_loss_tables.json")
RHO_20, ALPHA_20, MU0 = 1.7241e-8, 0.00393, 4e-7 * math.pi
COILS_PER_PHASE = 12                       # 36 slots, one coil a slot, three phases
POLE_PAIRS = 4
COIL_SIDES = 72                            # 36 slots x 2 layers
SERIES_TURNS_REF = 120                     # the FEMM model: 10 turns a coil side x 12 coils, one path
STACK_REF_MM = 20.0                        # the FEMM model's depth
ARRANGEMENTS = ("transposed", "row", "cluster", "random", "line")
_T: dict | None = None


# ---------------------------------------------------------------------------------------------- the tables (data from the studies)
def tables() -> dict:
    """The packaged data (copper_loss_tables.json), loaded once and turned into arrays."""
    global _T
    if _T is None:
        if not TABLES_PATH.exists():
            raise FileNotFoundError(f"{TABLES_PATH} is missing: run `uv run python -m optimization.copper_loss --build` after the studies (optimization/ac_*.py)")
        raw = json.loads(TABLES_PATH.read_text(encoding="utf-8"))
        t = dict(raw)
        m = raw["machine"]
        for which in ("G", "H"):                               # the field's harmonic groups: no-load part, total per AT^2 (bilinear), shares (linear)
            tb = m[which]
            t[which] = {"v0": np.array(tb["v0"]), "g_coef": np.array(tb["g_coef"]), "a_coef": np.array(tb["a_coef"]), "k_eff": np.array(tb["k_eff"])}
        t["psi_coef"] = np.array(m["psi_T_coef"])              # (6, 3): psi_d, psi_q, T at 120 series turns, 20 mm: quadratic in (i_d, i_q) of the 10-turn model
        p = raw["pwm"]
        t["pwm_m"], t["pwm_h"], t["pwm_frac"], t["pwm_tail"] = np.array(p["m"]), np.array(p["h_rms"]), np.array(p["frac"]), np.array(p["tail"])
        _T = t
    return _T


def resistivity(t_c):
    """Copper's resistivity [Ohm m] at t_c degC."""
    return RHO_20 * (1.0 + ALPHA_20 * (np.asarray(t_c, dtype=float) - 20.0))


def _q(d_m, f_hz, rho):
    """(gamma, z, q = I1(z)/I0(z)) for a wire of diameter d_m at f_hz (as optimization/ac_theory.py, with the exponentially scaled Bessel functions)."""
    delta = np.sqrt(rho / (math.pi * np.maximum(f_hz, 1e-12) * MU0))
    gamma = (1.0 + 1.0j) / delta
    z = gamma * d_m / 2.0
    return gamma, z, ive(1, z) / ive(0, z)


def skin_ratio(d_m, f_hz, rho):
    """R_ac / R_dc of a lone round wire: Re[(z / 2) / q]."""
    _, z, q = _q(d_m, f_hz, rho)
    return np.where(np.asarray(f_hz) > 0, np.real(0.5 * z / q), 1.0)


def prox_ratio(d_m, f_hz, rho):
    """The exact proximity loss of a lone wire in a uniform transverse field over the low-frequency formula (ac_theory.proximity_ratio)."""
    gamma, z, q = _q(d_m, f_hz, rho)
    a = d_m / 2.0
    delta = np.sqrt(rho / (math.pi * np.maximum(f_hz, 1e-12) * MU0))
    k2 = 2.0 / delta ** 2
    return np.where(np.asarray(f_hz) > 0, 16.0 * np.imag(gamma * np.conj(q)) / (a ** 3 * k2 ** 2), 1.0)


# ---------------------------------------------------------------------------------------------- the end turns
def end_turn_mm(model, e_ext_mm, angle_deg, geo: dict):
    """The length of one end turn [mm]: the 'diamond' (module docstring) 2 e_ext + c / cos(alpha), or the first study's 'arc' (pi / 2) span + 2 e (e = 2 mm unless e_ext is given)."""
    if model == "arc":
        return 0.5 * math.pi * geo["span_mid_mm"] + 2.0 * e_ext_mm
    return 2.0 * e_ext_mm + geo["chord_mm"] / np.cos(np.radians(angle_deg))


# ---------------------------------------------------------------------------------------------- the slot's AC resistance (from the FEMM slot model)
def _k_slot(layers: str, t: dict) -> float:
    """The slot constant K [m^2] for the phase relation of the layers: 'same' (both layers of one phase: 6 of the 9 slots), 'diff' (two phases 60 degrees apart: the other 3) or 'mixed' (the mean
    over the 9 slots)."""
    if layers == "same":
        return t["K_same_m2"]
    if layers == "diff":
        return t["K_diff_m2"]
    if layers == "mixed":
        return t["w_same"] * t["K_same_m2"] + (1.0 - t["w_same"]) * t["K_diff_m2"]
    raise ValueError(f"layers {layers!r}: 'same', 'diff' or 'mixed'")


def slot_extra_ratio(d_m, fill, f_hz, rho, layers: str = "same"):
    """(R_ac - R_dc) / R_dc of a slot of the optimized stator whose strands (diameter d_m, copper fill `fill` of the coil layer) all carry the same current at frequency f_hz: the loss of the strands
    in the slot's own leakage field, over their DC loss. The law of optimization/ac_slot_study.py:
        K (d k / delta^2)^2  x  the lone wire's exact proximity factor  x  c_iron(f)
    with K [m^2] the slot constant of the static field ('same' phase in both layers, or 'mixed': a weighted mean with the slots of two phases), delta^2 = rho / (pi f mu0), and c_iron(f) the
    ratio of FEMM's loss to the static-field law, a function of the frequency alone (the steel's laminations screening the field). Skin effect is not in it (skin_ratio)."""
    t = tables()["slot"]
    k = _k_slot(layers, t)
    f = np.maximum(np.asarray(f_hz, dtype=float), 1e-9)
    delta2 = rho / (math.pi * f * MU0)
    c_iron = np.interp(np.log(f), np.log(t["c_f"]), t["c"])
    return k * (d_m * fill / delta2) ** 2 * prox_ratio(d_m, f, rho) * c_iron


def slot_r_ratio(d_m, fill, f_hz, rho, layers: str = "same"):
    """R_ac / R_dc of the slot at f_hz: 1 + the lone strand's skin effect + slot_extra_ratio."""
    return skin_ratio(d_m, f_hz, rho) + slot_extra_ratio(d_m, fill, f_hz, rho, layers)


# ---------------------------------------------------------------------------------------------- the model
def _basis(i_d, i_q):
    """The quadratic basis (1, i_d, i_q, i_d^2, i_d i_q, i_q^2) of the flux-linkage fit, as the last axis."""
    return np.stack([np.ones_like(i_d), i_d, i_q, i_d ** 2, i_d * i_q, i_q ** 2], axis=-1)


def _groups(which: str, at, gamma_deg, t: dict):
    """The machine's field statistic of the harmonic groups at (AT, gamma): 'G' the mean square of dB/d(theta_e), 'H' the same for the potential's spread, shape (..., 7). The fit of
    optimization/ac_machine.py (fit_split): V_g = V_g(no load) + AT^2 g_tot(AT, gamma) a_g(AT, gamma), g_tot bilinear in (AT, gamma), the shares a_g linear (clipped at 0, summing to 1)."""
    tb = t[which]
    g = np.radians(gamma_deg)
    g_tot = tb["g_coef"][0] + tb["g_coef"][1] * at + tb["g_coef"][2] * g + tb["g_coef"][3] * at * g
    share = np.maximum(tb["a_coef"][0] + tb["a_coef"][1] * at[..., None] + tb["a_coef"][2] * g[..., None], 0.0)
    share = share / np.sum(share, axis=-1, keepdims=True)
    return tb["v0"] + (at ** 2 * g_tot)[..., None] * share


def losses_arrays(N, a, n_s, d_mm, d_over_mm, i_peak, gamma_deg, rpm, f_sw, v_dc, t_copper_c=100.0, stack_mm=20.0, *, arrangement: str = "transposed", end_model: str = "diamond",
                  end_ext_mm=None, end_angle_deg=30.0, ripple_layers: str = "same", l_scale=1.0) -> dict:
    """The six losses [W] of windings given as arrays (numpy broadcasting: every numeric argument may be a number or an array, all of one shape or broadcastable). The strings
    (arrangement, end_model, ripple_layers) are the same for every winding of the call. See the module docstring for the model and losses() for the keys; the values here are arrays."""
    t = tables()
    geo = t["geometry"]
    N, a, n_s, d_mm, d_over_mm, i, gam, rpm, f_sw, v_dc, T, L_mm = np.broadcast_arrays(*(np.asarray(v, dtype=float) for v in (N, a, n_s, d_mm, d_over_mm, i_peak, gamma_deg, rpm, f_sw, v_dc,
                                                                                                                            t_copper_c, stack_mm)))
    if arrangement not in ARRANGEMENTS:
        raise ValueError(f"arrangement {arrangement!r}: one of {ARRANGEMENTS}")
    rho = resistivity(T)
    sigma = 1.0 / rho
    d = d_mm * 1e-3
    S = math.pi * d ** 2 / 4.0
    L = L_mm * 1e-3
    n_ser = COILS_PER_PHASE * N / a                                       # series turns of a path (= of a phase): 12 N / a
    e_ext = (2.0 if end_model == "arc" else 3.0) if end_ext_mm is None else end_ext_mm
    l_end = end_turn_mm(end_model, e_ext, end_angle_deg, geo) * 1e-3     # one end turn [m]
    r_slot = rho * 2.0 * L * n_ser / (a * n_s * S)                        # the phase's DC resistance in the slots [Ohm]
    r_end = rho * 2.0 * l_end * n_ser / (a * n_s * S)
    r_dc = r_slot + r_end
    fill = N * n_s * S / (geo["A_layer_mm2"] * 1e-6)
    p_dc_slot, p_dc_end = 1.5 * r_slot * i ** 2, 1.5 * r_end * i ** 2
    f_e = POLE_PAIRS * rpm / 60.0
    w_e = 2.0 * math.pi * f_e
    # ---- skin effect of the fundamental current in its own strand
    p_skin = (p_dc_slot + p_dc_end) * (skin_ratio(d, f_e, rho) - 1.0)
    # ---- the machine's field: G_k at the operating point, in amperes of the 10-turn model
    m = t["machine"]
    at = N * i / a                                                       # ampere-turns a coil side, peak
    rg = m["ranges"]
    extrapolated = (at > rg["at_max"] * 1.001) | (gam > rg["gamma_max"] + 0.5) | (gam < rg["gamma_min_loaded"] - 0.5) & (at > 0.2 * m["at_rated"])
    at_c = np.minimum(at, rg["at_max"])
    gam_c = np.clip(gam, rg["gamma_min_loaded"], rg["gamma_max"])
    g_grp = _groups("G", at_c, gam_c, t)                                  # (..., 7): the mean square of dB/d(theta_e) of each harmonic group at the operating point
    k_eff = t["G"]["k_eff"]
    pr = prox_ratio(d[..., None], f_e[..., None] * k_eff, rho[..., None])
    p_prox = COIL_SIDES * N * n_s * L * (math.pi * sigma * d ** 4 * w_e ** 2 / 64.0) * np.sum(g_grp * pr, axis=-1)
    i_d_ref, i_q_ref = -(at_c / m["i_ref_turns"]) * np.sin(np.radians(gam_c)), (at_c / m["i_ref_turns"]) * np.cos(np.radians(gam_c))
    basis = _basis(i_d_ref, i_q_ref)
    # ---- the voltage the winding needs, and the inverter's ripple
    psi = (basis @ t["psi_coef"]) * (n_ser / SERIES_TURNS_REF * L_mm / STACK_REF_MM)[..., None]
    i_d, i_q = -i * np.sin(np.radians(gam)), i * np.cos(np.radians(gam))
    v_peak = np.hypot(r_dc * i_d - w_e * psi[..., 1], r_dc * i_q + w_e * psi[..., 0])
    mod = v_peak * math.sqrt(3.0) / np.maximum(v_dc, 1e-9)
    overmodulated = mod > 1.0
    mod_c = np.minimum(mod, 1.0)
    l_ph = m["L_d_ref_H"] * (n_ser / SERIES_TURNS_REF) ** 2 * (L_mm / STACK_REF_MM) * l_scale
    on = f_sw > 0
    i_rip = np.where(on, np.interp(mod_c, t["pwm_m"], t["pwm_h"]) * v_dc / (l_ph * np.maximum(f_sw, 1e-9)), 0.0)
    frac = np.stack([np.interp(mod_c, t["pwm_m"], t["pwm_frac"][:, h]) for h in range(t["pwm_frac"].shape[1])] + [np.interp(mod_c, t["pwm_m"], t["pwm_tail"])], axis=-1)   # (..., 9)
    harm = np.arange(1, frac.shape[-1] + 1, dtype=float)
    f_h = np.maximum(f_sw[..., None] * harm, 1.0)
    fac = slot_r_ratio(d[..., None], fill[..., None], f_h, rho[..., None], ripple_layers)
    skin_h = skin_ratio(d[..., None], f_h, rho[..., None])
    p_pwm = np.where(on, 3.0 * i_rip ** 2 * np.sum(frac * (r_slot[..., None] * fac + r_end[..., None] * skin_h), axis=-1), 0.0)
    # ---- circulating currents between the strands of a turn
    p_circ = _circulating(arrangement, n_s, d, fill, at_c, gam_c, w_e, f_e, rho, sigma, L, l_end, t, i_rip, frac, f_sw, r_slot, on, ripple_layers)
    total = p_dc_slot + p_dc_end + p_prox + p_skin + p_pwm + p_circ
    return {"p_dc_slot": p_dc_slot, "p_dc_end": p_dc_end, "p_prox": p_prox, "p_skin": p_skin, "p_pwm": p_pwm, "p_circ": p_circ, "p_total": total,
            "r_dc_ohm": r_dc, "r_slot_ohm": r_slot, "r_end_ohm": r_end, "fill": fill, "i_rms": i / math.sqrt(2.0), "f_e_hz": f_e,
            "d_over_delta_fund": d / np.sqrt(rho / (math.pi * np.maximum(f_e, 1e-9) * MU0)), "ripple_rms_a": i_rip, "modulation_index": mod, "l_phase_h": l_ph,
            "l_end_mm": l_end * 1e3 + 0.0 * i, "extrapolated": extrapolated, "overmodulated": overmodulated}


def _kappa(arrangement: str, n_s, tc: dict):
    """The share of the layer's spread of A that the strands of one turn see (ac_slot_study.build_circ_model): interpolated in ln n_s over the measured turns of n_s strands, 0 at n_s = 1."""
    ns = np.array([1.0] + list(map(float, tc["ns"])))
    k = np.array([0.0] + list(map(float, tc["kappa"][arrangement])))
    return np.interp(np.log(np.maximum(n_s, 1.0)), np.log(ns), k)


def _rolloff(arrangement: str, f, d, n_s, rho, tc: dict):
    """The factor by which the circulating currents' own field cuts the loss below the low-frequency law: r = 1 / (1 + (u / u_c)^p), with u the frequency (line, random, row) or the bundle size
    n_s (d / delta)^2 (cluster), the constants fitted to FEMM's points (ac_slot_study.build_circ_model)."""
    fit = tc["rolloff_fit"][arrangement]
    if arrangement == "cluster":
        u = n_s * d ** 2 * math.pi * f * MU0 / rho
    else:
        u = f
    return 1.0 / (1.0 + (u / fit["u_c"]) ** fit["p"])


def _circulating(arrangement, n_s, d, fill, at_c, gam_c, w_e, f_e, rho, sigma, L, l_end, t, i_rip, frac, f_sw, r_slot, on, ripple_layers):
    """The circulating-current loss [W]. Zero for a transposed winding (every strand sees every position equally: equal sharing) and for n_s = 1. For the other placements of the strands in a
    turn the strands share a voltage and exchange the current that the difference of the vector potential A between their positions drives (optimization/ac_slot_study.py, circ stage):
      in the 2-D slot of length L the loss of one layer's strands is (N n_s S) sigma w^2 kappa V_A L at low frequency (S the strand's area, V_A the peak-phasor spread of A over the layer, kappa
      the share of it that strands of one turn see: 1 - 1 / n_s for random placement), times a roll-off r (the exchanged currents' own field opposes the EMF); N n_s S = fill A_layer;
      the real turn is longer than the active length (two sides of length L and two end turns), the loop's resistance grows with it and the EMF does not: the loss is the 2-D one times
      L / L_turn, L_turn = 2 L + 2 l_end, summed over the 72 coil sides.
    Two parts: (1) the machine's own field, A's spread per harmonic group H_g(AT, gamma) from the sector solves (ac_machine.py), sigma w_e^2 kappa sum_g H_g r_g; (2) the PWM ripple's field,
    (2 pi f_h sigma fill A_layer)^2 kappa v_A r(f_h) times the ripple's DC loss, at the carrier's harmonics f_h (the slot's own field, from the static solve)."""
    if arrangement == "transposed":
        return np.zeros(np.shape(f_e))
    tc = t["circ"]
    a_layer = t["geometry"]["A_layer_mm2"] * 1e-6
    kappa = _kappa(arrangement, n_s, tc)
    l_turn = 2.0 * L + 2.0 * l_end
    h_grp = _groups("H", at_c, gam_c, t)
    f_g = f_e[..., None] * t["H"]["k_eff"]
    r_g = _rolloff(arrangement, np.maximum(f_g, 1e-9), d[..., None], n_s[..., None], rho[..., None], tc)
    p_fund = COIL_SIDES * fill * a_layer * sigma * w_e ** 2 * kappa * L * np.sum(h_grp * r_g, axis=-1) * L / l_turn
    harm = np.arange(1, frac.shape[-1] + 1, dtype=float)
    f_h = np.maximum(f_sw[..., None] * harm, 1.0)
    q = (2.0 * math.pi * f_h * sigma[..., None] * (fill * a_layer)[..., None]) ** 2 * (kappa * tc["v_A"])[..., None] * _rolloff(arrangement, f_h, d[..., None], n_s[..., None], rho[..., None], tc)
    q = q * _k_slot(ripple_layers, t["slot"]) / t["slot"]["K_same_m2"]            # the spread of A scales as the slot constant K (same-phase layers are the reference)
    p_pwm = np.where(on, 3.0 * i_rip ** 2 * r_slot * np.sum(frac * q, axis=-1) * L / l_turn, 0.0)
    return p_fund + p_pwm


def losses(winding: dict, op: dict, t_copper_c: float = 100.0, stack_mm: float = 20.0) -> dict:
    """The copper loss of one winding at one operating point [W for the stack length given], as a dict of floats:
      p_dc_slot, p_dc_end   DC loss in the slots' active length, and in the end turns (1.5 R i^2 with the copper's resistivity at t_copper_c)
      p_prox                proximity loss of the machine's field (rotor + armature harmonics) on the strands of the slots
      p_skin                skin effect of the fundamental current (lone-strand Bessel result times the DC loss)
      p_pwm                 loss of the inverter's PWM ripple current (ripple rms at f_sw and the harmonics of the carrier, times the slot's AC resistance there)
      p_circ                circulating currents between the strands of a turn (0 for the default 'transposed' arrangement)
      p_total               their sum
    and the numbers behind them: r_dc_ohm (phase DC resistance, slots + end turns), r_slot_ohm, r_end_ohm, fill (copper fill of a coil layer), i_rms, f_e_hz (fundamental),
    d_over_delta_fund, ripple_rms_a (the phase current's ripple rms), modulation_index, l_phase_h (the inductance the ripple sees), l_end_mm (one end turn), extrapolated (the operating
    point lies outside the FEMM points of the machine's field: the field is clipped to them), overmodulated (the inverter cannot make the winding's voltage: ripple taken at index 1).
    winding = {'N', 'a', 'n_s', 'd_mm', 'd_over_mm', optional 'arrangement', 'end_model', 'end_ext_mm', 'end_angle_deg'};
    op = {'i_peak', 'gamma_deg', 'rpm', 'f_sw', 'v_dc', optional 'ripple_layers', 'l_scale'}. Currents are peak. See the module docstring."""
    r = losses_arrays(winding["N"], winding["a"], winding["n_s"], winding["d_mm"], winding["d_over_mm"], op["i_peak"], op["gamma_deg"], op["rpm"], op["f_sw"], op["v_dc"], t_copper_c, stack_mm,
                      arrangement=winding.get("arrangement", "transposed"), end_model=winding.get("end_model", "diamond"), end_ext_mm=winding.get("end_ext_mm"),
                      end_angle_deg=winding.get("end_angle_deg", 30.0), ripple_layers=op.get("ripple_layers", "same"), l_scale=op.get("l_scale", 1.0))
    return {k: (bool(v) if np.asarray(v).dtype == bool else float(v)) for k, v in r.items()}


# ---------------------------------------------------------------------------------------------- the command line: the interface at work, the check against the first study, the speed
def main(argv=None) -> int:
    import argparse
    import time
    ap = argparse.ArgumentParser(description="copper_loss: print the interface at work, check it against the first study, time it; --build makes the tables")
    ap.add_argument("--build", action="store_true", help="make copper_loss_tables.json from the studies' outputs")
    args = ap.parse_args(argv)
    if args.build:
        from optimization import copper_loss_build
        copper_loss_build.build()
        return 0
    w1 = {"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367}
    op = {"i_peak": 118.6 * 4 / 26, "gamma_deg": 22.47, "rpm": 3000.0, "f_sw": 20000.0, "v_dc": 48.0}
    r = losses(w1, op, 100.0, 20.0)
    print("losses(winding, op, t_copper_c, stack_mm), the winding study's example (26 turns x 8 x 0.315 mm, 4 paths) at 118.6 ampere-turns, 3,000 rpm, 20 kHz, 48 V, 100 degC, 20 mm [W]:")
    for k in ("p_dc_slot", "p_dc_end", "p_prox", "p_skin", "p_pwm", "p_circ", "p_total"):
        print(f"   {k:10s} {r[k]:12.6f}")
    print("   behind them: " + ", ".join(f"{k} {r[k]:.4g}" for k in ("r_dc_ohm", "fill", "i_rms", "f_e_hz", "ripple_rms_a", "modulation_index", "l_phase_h", "l_end_mm")))
    # the check against the first study: the DC loss with its end-turn stand-in is winding.py's R i^2
    from optimization import winding
    import motorgeo as mg
    from optimization import design
    ov = json.loads((Path(__file__).resolve().parents[1] / "out" / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]
    g = winding.end_turn_geometry(design.params(ov))
    arc = losses(dict(w1, end_model="arc", end_ext_mm=2.0), dict(op, i_peak=100.0 * 4 / 26, rpm=0.0), 100.0, 20.0)
    rho = winding.resistivity(100.0)
    r_first = rho * g["l_mt_mm"] * 1e-3 * 12 * 26 / (4 ** 2 * 8 * math.pi * (0.315e-3) ** 2 / 4)
    p_first = 1.5 * r_first * (100.0 * 4 / 26) ** 2
    print(f"\nthe first study's DC loss at 100 ampere-turns (its end-turn stand-in, l_mt = {g['l_mt_mm']:.2f} mm): {p_first:.6f} W; here with end_model='arc': {arc['p_dc_slot'] + arc['p_dc_end']:.6f} W "
          f"(relative difference {(arc['p_dc_slot'] + arc['p_dc_end']) / p_first - 1:+.1e}); the diamond end turn ({r['l_end_mm']:.2f} mm against {1e3 * 0.5 * (g['l_mt_mm'] * 1e-3 - 2 * 0.020):.2f} mm) gives "
          f"{(r['p_dc_slot'] + r['p_dc_end']) * (100.0 / 118.6) ** 2:.4f} W at the same 100 ampere-turns")
    n = 3000
    t0 = time.perf_counter()
    for _ in range(n):
        losses(w1, op, 100.0, 20.0)
    per_call = (time.perf_counter() - t0) / n
    big = 100000
    nn = np.random.default_rng(0).integers(2, 40, big).astype(float)
    t0 = time.perf_counter()
    losses_arrays(nn, 4.0, 8.0, 0.315, 0.367, 118.6 * 4.0 / nn, 22.47, 3000.0, 20000.0, 48.0, 100.0, 20.0)
    per_winding = (time.perf_counter() - t0) / big
    print(f"\nspeed: {1e6 * per_call:.0f} microseconds a call of losses(); {1e6 * per_winding:.2f} microseconds a winding in an array of {big:,} (losses_arrays)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
