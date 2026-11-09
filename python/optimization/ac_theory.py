"""optimization.ac_theory: the exact AC results of round copper wires that the rest of the copper-loss work is measured against (no FEMM, no plotting).

What it is, and why
    The winding study (optimization/winding.py) treated the AC loss of the winding with one low-frequency formula and showed it was tiny. This module holds the
    exact answers that formula is the first term of, so that every FEMM time-harmonic result of the copper-loss work (optimization/ac_calibration.py,
    ac_slot.py) is compared with a known answer before it is trusted, and so that the fast model (optimization/copper_loss.py) can use the corrections at
    high frequency. Four closed-form results, each built from the one before:

    1. Skin effect. One round wire of radius a, carrying a current I of angular frequency w, alone in space. Inside the wire the vector potential A (the
       z-component; B = curl A) and the current density J = -j w sigma A (+ a constant) obey the diffusion equation
           laplacian(A) = j w mu0 sigma A = gamma^2 A,         gamma = (1 + j) / delta,   delta = sqrt(2 / (w mu0 sigma)) = sqrt(rho / (pi f mu0))
       (delta is the skin depth: 5.36 mm at 200 Hz and 0.54 mm at 20 kHz for copper at 100 degC). For a current that is the same all round the wire the solution
       depends on the radius r only, laplacian(A) = A'' + A'/r, which is Bessel's modified equation of order 0: J(r) = J_a I_0(gamma r) / I_0(gamma a) with I_0 the modified
       Bessel function and J_a the current density at the surface. The total current is the integral of J over the cross section, I = 2 pi J_a a I_1(gamma a) / (gamma I_0(gamma a))
       (using integral of r I_0(gamma r) dr = r I_1(gamma r) / gamma). The power is the real part of the voltage times the current; the voltage along the wire is
       E(a) = J_a / sigma, so the impedance per unit length is Z' = E(a) / I = (gamma / (2 pi a sigma)) I_0(gamma a) / I_1(gamma a). Dividing by the DC
       resistance per unit length, R'_dc = 1 / (sigma pi a^2):
           R_ac / R_dc  =  Re[ (gamma a / 2) I_0(gamma a) / I_1(gamma a) ]                                                       (skin_ratio)
       It is 1 + (d / delta)^4 / 768 + ... at low frequency (d = 2a; the usual 1 + (a / delta)^4 / 48), 1.0013 at d / delta = 1, 1.0205 at 2, 1.0973 at 3, and
       d / (4 delta) + 1/4 at high frequency (the current lives in a layer one skin depth thick). skin_ratio_kelvin() gets the same number from the Kelvin functions ber and bei
       (the textbook way: gamma a = x (1 + j) / sqrt 2 with x = sqrt 2 a / delta), an independent route used as a check.
    2. Proximity effect: no net current, a uniform transverse field. The same wire carries no net current but sits in a uniform field B (peak, along y; A =
       -B x far away). The field is now not round: the eddy current J = -j w sigma A depends on the angle as cos(phi). Inside, A = C I_1(gamma r) cos(phi);
       outside, A = (-B r + D / r) cos(phi). Matching A and dA/dr at r = a (the two unknowns C and D) gives C = -2 B / (gamma I_0(gamma a)). The power per
       unit length is P' = integral of |J|^2 / (2 sigma) over the disc = (pi sigma w^2 / 2) |C|^2 integral from 0 to a of |I_1(gamma r)|^2 r dr, and the integral
       has a closed form. With f = I_1(gamma r) and its conjugate g, Bessel's equation gives d/dr [r (f' g - f g')] = r (gamma^2 - conj gamma^2) f g, and
       gamma^2 - conj gamma^2 = 2 j w mu0 sigma (Lommel's integral), so the integral of r |f|^2 from 0 to a is a Im[f'(a) conj f(a)] / (w mu0 sigma), and with
       q = I_1(gamma a) / I_0(gamma a):
           P'  =  (2 pi a B^2 / (mu0^2 sigma)) Im[ gamma conj(q) ]                                                               (proximity_exact)
       At low frequency q = (gamma a / 2)(1 - (gamma a)^2 / 8 + ...) and this becomes
           P'  =  pi sigma w^2 B^2 d^4 / 128                                                                                    (proximity_lowfreq)
       the formula of the winding study (B peak; the time average of the power; d = 2a). proximity_ratio() is exact over low-frequency: 0.993 at d / delta = 1,
       0.898 at 2, 0.638 at 3. The wire is alone in the field: its neighbours' own eddy fields are NOT in this result (that is what the slot model adds).
    3. The current density in both cases is returned by skin_profile() and proximity_current_density(), for the pictures of learn/proximity.
    4. Two wires in parallel in a transverse field (the circulating current). Two identical wires, centres s apart across the field, joined in parallel at both
       ends, no net current: each has E_z = -j w A + V with the same V, A = -B x differs by B s between them, so with V = 0 (net current zero by symmetry) the
       currents are +-j w sigma S B s / 2 (S the wire's cross-section area), and the circulating-current loss of the pair, at low frequency, is
           P'_circ  =  w^2 sigma S B^2 s^2 / 4                                                                                  (two_wire_circulating)
       on top of the two wires' own proximity loss: (4 s^2 / d^2) times the proximity loss of the pair for two touching wires (s = d): circulating currents are
       4 times the eddy loss when the wires touch, and (s / d)^2 times more when they are s apart.

Copper at a temperature
    rho(T) = rho_20 (1 + alpha_20 (T - 20)), rho_20 = 1.7241e-8 Ohm m (annealed copper, IEC 60028: 1/58 Ohm mm^2/m), alpha_20 = 0.00393 /K (the same two numbers and the
    same line as optimization/winding.py and the project notes' materials guide, 02-current.md section 4.3): 2.2662e-8 Ohm m at 100 degC. Peak values everywhere
    (B and currents are peak amplitudes of sinusoids; powers are time averages), the convention of FEMM's time-harmonic solver.

Run
    uv run python -m optimization.ac_theory       (from the Python to FEMM folder; no FEMM, about a second) prints the self-checks and exits 0 if all pass.

Numbers to expect (copper at 100 degC; the checks print them)
    d / delta        0.1        0.5        1          2          3          5
    R_ac / R_dc      1.0000001  1.0000814  1.0013007  1.0204924  1.0973282  1.5050394   (skin_ratio, alone in space)
    prox ratio       0.9999993  0.9995526  0.9928916  0.8976411  0.6384099  0.2067688   (proximity_ratio, exact over low-frequency)
    The two closed forms agree with the numerical integral of the same Bessel functions (proximity_ratio_quadrature) and, for the skin effect, with the Kelvin-function
    form, to better than 1e-10.

Self-check (the answers follow).
  1. At 200 Hz (d / delta = 0.12 for a 0.63 mm wire) the skin factor is 1.00000. Does that mean a slot of such wires has no AC loss at 200 Hz?
  2. Why does one strand's proximity loss go as d^4, while the proximity loss of a whole winding of a fixed copper area goes only as d^2?
  3. Two touching wires in parallel in a transverse field: how much is their circulating-current loss, in units of the eddy loss of the same two wires?
Answers.
  1. No. The skin factor is only the wire's own current crowding against itself. The proximity loss (the slot's and the neighbours' field acting on it) and the
     circulating currents between strands come on top: they depend on the winding's field and layout, not on the lone wire.
  2. In a field changing as dB/dt, the electric field at distance y from the wire's centre is E_z = -y dB/dt, so the loss density sigma E^2 goes as y^2 and its
     integral over the disc is sigma (dB/dt)^2 pi d^4 / 64: d^4. A given copper area holds a number of strands that goes as 1 / d^2, so the winding's total goes as d^2
     (and its DC loss, for a given copper area and current, does not depend on d at all).
  3. 4 (s = d): w^2 sigma S B^2 d^2 / 4 with S = pi d^2 / 4 is pi sigma w^2 B^2 d^4 / 16, and the pair's eddy loss is 2 pi sigma w^2 B^2 d^4 / 128 = pi sigma w^2 B^2 d^4 / 64:
     the ratio is 4, so the pair loses five times its eddy loss when the two are one circuit.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.special import berp, ber, bei, beip, ive

MU0 = 4e-7 * math.pi                  # H/m
RHO_20 = 1.7241e-8                    # Ohm m: annealed copper at 20 degC (IEC 60028)
ALPHA_20 = 0.00393                    # per K: copper's temperature coefficient of resistance at 20 degC (IEC 60028)


def resistivity(t_c: float = 100.0):
    """Copper's resistivity [Ohm m] at t_c degC: rho_20 (1 + alpha_20 (t_c - 20)); 2.2662e-8 at 100 degC. Works on arrays."""
    return RHO_20 * (1.0 + ALPHA_20 * (np.asarray(t_c, dtype=float) - 20.0))


def skin_depth(f_hz, rho):
    """delta = sqrt(rho / (pi f mu0)) [m]: the depth at which an alternating current in copper has fallen to 1 / e."""
    return np.sqrt(np.asarray(rho, dtype=float) / (math.pi * np.asarray(f_hz, dtype=float) * MU0))


def _bessel_q(d_m, f_hz, rho):
    """(gamma, z, q) for a wire of diameter d_m at f_hz: gamma = (1 + j) / delta, z = gamma a, q = I_1(z) / I_0(z) (the exponentially scaled Bessel functions are used,
    so that z of 100 or more does not overflow: their common factor exp(|Re z|) cancels in the ratio)."""
    d_m, f_hz, rho = np.broadcast_arrays(np.asarray(d_m, dtype=float), np.asarray(f_hz, dtype=float), np.asarray(rho, dtype=float))
    delta = skin_depth(f_hz, rho)
    gamma = (1.0 + 1.0j) / delta
    z = gamma * d_m / 2.0
    return gamma, z, ive(1, z) / ive(0, z)


def skin_ratio(d_m, f_hz, rho):
    """R_ac / R_dc of one round wire of diameter d_m [m] alone in space, at f_hz, resistivity rho: Re[(z / 2) I_0(z) / I_1(z)] with z = (1 + j) d / (2 delta).
    1 + (d / delta)^4 / 768 at low frequency; see the module docstring."""
    _, z, q = _bessel_q(d_m, f_hz, rho)
    return np.real(0.5 * z / q)


def skin_ratio_kelvin(d_m, f_hz, rho):
    """The same R_ac / R_dc from the Kelvin functions (an independent route): with x = sqrt(2) a / delta,
    R_ac / R_dc = (x / 2) (ber x bei' x - ber' x bei x) / (ber'^2 x + bei'^2 x). Valid for x up to about 20 (ber and bei overflow beyond)."""
    d_m, f_hz, rho = np.broadcast_arrays(np.asarray(d_m, dtype=float), np.asarray(f_hz, dtype=float), np.asarray(rho, dtype=float))
    x = math.sqrt(2.0) * (d_m / 2.0) / skin_depth(f_hz, rho)
    return 0.5 * x * (ber(x) * beip(x) - berp(x) * bei(x)) / (berp(x) ** 2 + beip(x) ** 2)


def proximity_lowfreq(d_m, f_hz, b_peak, rho):
    """The low-frequency proximity loss of one round wire with no net current in a uniform transverse field, per unit length [W/m]:
    pi sigma w^2 B^2 d^4 / 128 (B the peak flux density, the power a time average): the first term of proximity_exact()."""
    d_m, f_hz, b_peak, rho = (np.asarray(v, dtype=float) for v in (d_m, f_hz, b_peak, rho))
    w = 2.0 * math.pi * f_hz
    return math.pi * w ** 2 * b_peak ** 2 * d_m ** 4 / (128.0 * rho)


def proximity_ratio(d_m, f_hz, rho):
    """The exact proximity loss of a round wire in a uniform transverse field over the low-frequency formula: 16 Im[gamma conj(q)] / (a^3 k^4), k^2 = 2 / delta^2.
    1 at low frequency; 0.993 at d / delta = 1, 0.898 at 2, 0.581 at 3 (the wire shields itself, so the formula overestimates there)."""
    d_m = np.asarray(d_m, dtype=float)
    gamma, _, q = _bessel_q(d_m, f_hz, rho)
    a = d_m / 2.0
    k2 = 2.0 / skin_depth(f_hz, rho) ** 2
    return 16.0 * np.imag(gamma * np.conj(q)) / (a ** 3 * k2 ** 2)


def proximity_exact(d_m, f_hz, b_peak, rho):
    """The exact proximity loss per unit length [W/m]: proximity_lowfreq x proximity_ratio."""
    return proximity_lowfreq(d_m, f_hz, b_peak, rho) * proximity_ratio(d_m, f_hz, rho)


def proximity_ratio_quadrature(d_m: float, f_hz: float, rho: float, n: int = 400) -> float:
    """The same ratio by brute force (an independent route): (pi sigma w^2 / 2) |C|^2 integral of |I_1(gamma r)|^2 r dr, with C = -2 B / (gamma I_0(gamma a)), over
    the low-frequency value pi sigma w^2 B^2 a^4 / 8, the integral by Gauss-Legendre points. Scalars only."""
    from scipy.special import iv
    a = d_m / 2.0
    delta = float(skin_depth(f_hz, rho))
    g = (1.0 + 1.0j) / delta
    x, w = np.polynomial.legendre.leggauss(n)
    r, wr = 0.5 * a * (x + 1.0), 0.5 * a * w
    integral = float(np.sum(np.abs(iv(1, g * r)) ** 2 * r * wr))
    return 16.0 * integral / (a ** 4 * abs(g) ** 2 * abs(iv(0, g * a)) ** 2)


def skin_profile(r_over_a, d_m, f_hz, rho):
    """|J(r)| / |J(0)|... returns the complex J(r) / (I / (pi a^2)): the current density of one wire alone at radius r, in units of the DC density of the same current.
    J(r) / J_dc = (gamma a / 2) I_0(gamma r) / I_1(gamma a). At low frequency it is 1 everywhere; at high frequency it piles up at r = a (the surface)."""
    from scipy.special import iv
    d_m = float(d_m)
    a = d_m / 2.0
    g = (1.0 + 1.0j) / float(skin_depth(f_hz, rho))
    r = np.asarray(r_over_a, dtype=float) * a
    return (g * a / 2.0) * iv(0, g * r) / iv(1, g * a)


def proximity_current_density(x_over_a, y_over_a, d_m, f_hz, b_peak, rho):
    """The eddy current density J_z [A/m^2, complex peak] of a round wire with no net current in a uniform field B (peak) along y, at the point (x, y) (in units of the
    radius a, origin at the wire's centre): J = -j w sigma C I_1(gamma r) cos(phi) with C = -2 B / (gamma I_0(gamma a)), A = -B x far away. Points outside the wire give 0."""
    from scipy.special import iv
    a = float(d_m) / 2.0
    delta = float(skin_depth(f_hz, rho))
    g = (1.0 + 1.0j) / delta
    sigma = 1.0 / float(rho)
    w = 2.0 * math.pi * float(f_hz)
    x, y = np.asarray(x_over_a, dtype=float) * a, np.asarray(y_over_a, dtype=float) * a
    r = np.hypot(x, y)
    cosphi = np.divide(x, r, out=np.zeros_like(r), where=r > 0)
    c = -2.0 * b_peak / (g * iv(0, g * a))
    j = -1j * w * sigma * c * iv(1, g * r) * cosphi
    return np.where(r <= a, j, 0.0)


def two_wire_circulating(d_m, s_m, f_hz_or_w, b_peak, rho, omega: bool = False):
    """The low-frequency circulating-current loss per unit length [W/m] of TWO identical round wires (diameter d_m) in parallel, centres s_m apart across a uniform
    transverse field b_peak (no net current): w^2 sigma S B^2 s^2 / 4 with S = pi d^2 / 4. The second argument is f in Hz, or w in rad/s with omega=True."""
    w = f_hz_or_w if omega else 2.0 * math.pi * f_hz_or_w
    sigma = 1.0 / np.asarray(rho, dtype=float)
    s_area = math.pi * np.asarray(d_m, dtype=float) ** 2 / 4.0
    return w ** 2 * sigma * s_area * np.asarray(b_peak, dtype=float) ** 2 * np.asarray(s_m, dtype=float) ** 2 / 4.0


def _check(name: str, value: float, expect: float, tol: float, results: list) -> None:
    ok = abs(value - expect) <= tol * max(1.0, abs(expect))
    results.append(ok)
    print(f"  {'ok  ' if ok else 'FAIL'} {name}: {value:.10g} (expected {expect:.10g}, tolerance {tol:g})")


def self_checks() -> bool:
    """The known answers of this module, against each other and against the first study's own proximity_ratio()."""
    rho = float(resistivity(100.0))
    out: list[bool] = []
    d = 0.63e-3
    print("copper at 100 degC: rho = %.5e Ohm m, delta(200 Hz) = %.3f mm, delta(20 kHz) = %.3f mm" % (rho, 1e3 * skin_depth(200.0, rho), 1e3 * skin_depth(2e4, rho)))
    for ratio in (0.1, 0.5, 1.0, 2.0, 3.0, 5.0):
        f = rho / (math.pi * MU0 * (d / ratio) ** 2)
        s_b, s_k = float(skin_ratio(d, f, rho)), float(skin_ratio_kelvin(d, f, rho))
        p_c, p_q = float(proximity_ratio(d, f, rho)), proximity_ratio_quadrature(d, f, rho)
        print(f"  d/delta {ratio:4.1f}: f = {f:10.1f} Hz   R_ac/R_dc {s_b:.8f}   proximity ratio {p_c:.8f}")
        _check("skin_ratio vs Kelvin functions", s_b, s_k, 1e-10, out)
        _check("proximity_ratio closed form vs quadrature", p_c, p_q, 1e-9, out)
    f = rho / (math.pi * MU0 * (d / 0.1) ** 2)
    _check("skin_ratio at d/delta = 0.1 vs 1 + (d/delta)^4/768", float(skin_ratio(d, f, rho)), 1.0 + 0.1 ** 4 / 768.0, 1e-9, out)
    f5 = rho / (math.pi * MU0 * (d / 20.0) ** 2)
    _check("skin_ratio at d/delta = 20 vs d/(4 delta) + 1/4", float(skin_ratio(d, f5, rho)), 20.0 / 4.0 + 0.25, 5e-3, out)
    _check("proximity_ratio at d/delta = 0.1 (close to 1)", float(proximity_ratio(d, f, rho)), 1.0, 1e-5, out)
    try:                                              # the first study's own function, written separately, for the same numbers
        from optimization import winding
        for ratio in (1.0, 2.0, 3.0):
            f = rho / (math.pi * MU0 * (d / ratio) ** 2)
            _check(f"winding.proximity_ratio at d/delta {ratio:g}", float(proximity_ratio(d, f, rho)), winding.proximity_ratio(d, f, rho), 1e-9, out)
    except ImportError as error:
        print("  (the first study's module could not be imported: %s)" % error)
        out.append(False)
    p_pair = float(two_wire_circulating(d, d, 1000.0, 1.0, rho))
    p_eddy_pair = 2.0 * float(proximity_lowfreq(d, 1000.0, 1.0, rho))
    _check("touching wires: circulating loss over the pair's eddy loss (4)", p_pair / p_eddy_pair, 4.0, 1e-9, out)
    return all(out)


if __name__ == "__main__":
    import sys
    ok = self_checks()
    print("all checks pass" if ok else "SOME CHECKS FAILED")
    sys.exit(0 if ok else 1)
