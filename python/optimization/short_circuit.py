"""optimization.short_circuit: what the motor does at speed if its inverter stops controlling it, and what an active short circuit costs.

What it is, and why
    A permanent-magnet motor makes a voltage whenever it turns, whether the inverter controls it or not. Its no-load flux linkage psi_m gives the
    phase voltage w_e psi_m (peak), so the line-to-line peak is sqrt 3 w_e psi_m. Above the speed where that equals the DC bus, the inverter's
    diodes conduct by themselves: if the switches stop (a fault, a trip), the motor drives current back into the bus and charges it, whatever the
    controller wants. That is "uncontrolled generation". The usual protection is to close the three low-side switches together, an ACTIVE SHORT
    CIRCUIT: the terminal voltage becomes zero and the current is set by the motor alone. This file computes that current and its braking torque
    from the flux maps, and the speeds where the magnets' voltage passes the bus.

The equations, each one built from the one before (d is the magnet axis; currents and voltages are PEAK phase values, amplitude-invariant d-q)
    1. In steady state, in the frame that turns with the rotor (motor_map.py, step 2):
           v_d = R i_d - w_e psi_q(i_d, i_q)          v_q = R i_q + w_e psi_d(i_d, i_q)
    2. The short circuit sets v_d = v_q = 0, so the current solves
           R i_d = w_e psi_q          R i_q = -w_e psi_d
       Two equations, two unknowns, with the flux maps for psi (cross-saturation included): scipy's fsolve.
    3. At high speed R << w_e L, so psi_d -> 0 and psi_q -> 0: the current turns almost all toward -d and tends to the CHARACTERISTIC current,
       the d current that cancels the magnets' flux, psi_d(i_d, 0) = 0 (brentq on the maps). With constant inductances it is psi_m / L_d.
    4. The torque is a brake, and all of the shaft power goes into the copper:
           T = 1.5 p (psi_d i_q - psi_q i_d) = -1.5 R |i|^2 / w_m          (substitute step 2: check (a) below)
       So the braking torque is small at high speed (the loss is fixed, the speed large) and largest at a few hundred rpm.
    5. The diodes start to conduct at the speed where sqrt 3 p w_m psi_m = V_dc (no load, R ignored).

The checks (run with the main block)
    (a) the torque from the maps equals -P_cu / w_m at every speed (energy, step 4), to 1e-9 relative
    (b) the residual of step 2 at the solution, below 1e-6 V
    (c) at the top speed the current is within 1 % of the characteristic current (step 3)

The inputs: the flux maps (maps.FluxMaps: FEMM 5.0, fine mesh, the optimized motor), the winding study's headline winding (9 turns a coil side,
4 paths, 60 mm) and its DC resistance at 100 C (19.261 mOhm, optimization/winding_study.py). Any winding is an argument.

What it does not do: the transient at the instant of the short. The current overshoots the steady value for a few cycles, by up to about
twice, depending on when the switches close; that peak, against the magnets' demagnetisation (step EM-6) and the switches, needs a transient
model. And the case where the battery is disconnected (the bus then rides up toward the line-to-line peak).

Call:  uv run python -m optimization.short_circuit          (prints the table, about 1 s)
       steady_short_circuit(fm, rpm, turns, paths, stack_mm, r_ohm) -> dict ;  characteristic_current(fm, turns, paths, stack_mm) -> A
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq, fsolve

from optimization.maps import FluxMaps

HEADLINE = {"turns": 9, "paths": 4, "stack_mm": 60.0, "r_ohm": 19.261e-3}       # optimization/winding_study.py, 60 mm, R at 100 C
SPEEDS = (300, 500, 1000, 2000, 3000, 5000, 10000, 16330)                         # rpm; 16,330 is the rotor's safe speed (Rotor stress)
BUSES = (48.0, 60.0)


def winding_range(fm: FluxMaps, turns: float, paths: float) -> tuple[float, float, float]:
    """The solved region of the maps in this winding's amps: (lowest i_d, highest i_d, largest |i_q|)."""
    s = (turns / paths) / 10.0                                                     # maps.winding_scale: the placeholder is 10 turns, one path
    return fm.range_d[0] / s, fm.range_d[1] / s, fm.range_q / s


def no_load_flux(fm: FluxMaps, turns: float, paths: float, stack_mm: float) -> float:
    """psi_m [Wb, peak]: the flux linkage of the magnets alone, at zero current."""
    return fm.psi_dq(0.0, 0.0, turns=turns, paths=paths, stack_mm=stack_mm)[0]


def conduction_speed(psi_m: float, pole_pairs: int, v_dc: float) -> float:
    """The speed [rpm] where the no-load line-to-line peak sqrt 3 w_e psi_m equals the bus (step 5)."""
    return v_dc / (np.sqrt(3.0) * psi_m * pole_pairs) * 60.0 / (2.0 * np.pi)


def characteristic_current(fm: FluxMaps, turns: float, paths: float, stack_mm: float) -> float:
    """The d current [A peak, negative] that cancels the magnets' flux: psi_d(i_d, 0) = 0 (step 3). NaN if it lies beyond the maps."""
    lo = winding_range(fm, turns, paths)[0] * (1 - 1e-6)
    f = lambda i: fm.psi_dq(i, 0.0, turns=turns, paths=paths, stack_mm=stack_mm)[0]
    if f(lo) > 0:
        return float("nan")
    return brentq(f, lo, 0.0, xtol=1e-9)


def steady_short_circuit(fm: FluxMaps, rpm: float, turns: float, paths: float, stack_mm: float, r_ohm: float, guess=None) -> dict:
    """The steady short-circuit current at one speed (step 2), its braking torque and copper loss."""
    w_e = rpm * 2.0 * np.pi / 60.0 * fm.pole_pairs
    w_m = rpm * 2.0 * np.pi / 60.0

    def residual(x):
        psi_d, psi_q = fm.psi_dq(x[0], x[1], turns=turns, paths=paths, stack_mm=stack_mm, fill="clip")
        return [r_ohm * x[0] - w_e * psi_q, r_ohm * x[1] + w_e * psi_d]

    x0 = guess if guess is not None else [0.8 * winding_range(fm, turns, paths)[0], 0.0]
    sol, _, ier, msg = fsolve(residual, x0, full_output=True, xtol=1e-12)
    i_d, i_q = float(sol[0]), float(sol[1])
    lo, hi, qmax = winding_range(fm, turns, paths)
    p_cu = 1.5 * r_ohm * (i_d ** 2 + i_q ** 2)                                     # three phases, peak values: 1.5 R i^2
    return {
        "rpm": rpm, "i_d": i_d, "i_q": i_q, "i": float(np.hypot(i_d, i_q)),
        "torque": fm.torque_dq(i_d, i_q, turns=turns, paths=paths, stack_mm=stack_mm, fill="clip"),
        "p_cu": p_cu, "torque_from_loss": -p_cu / w_m,
        "residual_v": float(np.max(np.abs(residual(sol)))),
        "inside_maps": bool(lo <= i_d <= hi and abs(i_q) <= qmax), "converged": ier == 1, "message": msg,
    }


def table(fm: FluxMaps | None = None, turns=HEADLINE["turns"], paths=HEADLINE["paths"], stack_mm=HEADLINE["stack_mm"],
          r_ohm=HEADLINE["r_ohm"], speeds=SPEEDS) -> list[dict]:
    """The short circuit at each speed, from low to high, each solve starting from the one before."""
    fm = fm or FluxMaps.load()
    rows, guess = [], None
    for rpm in speeds:
        row = steady_short_circuit(fm, rpm, turns, paths, stack_mm, r_ohm, guess)
        guess = [row["i_d"], row["i_q"]]
        rows.append(row)
    return rows


def main() -> None:
    fm = FluxMaps.load()
    h = HEADLINE
    psi_m = no_load_flux(fm, h["turns"], h["paths"], h["stack_mm"])
    lo, hi, qmax = winding_range(fm, h["turns"], h["paths"])
    i_ch = characteristic_current(fm, h["turns"], h["paths"], h["stack_mm"])
    w_top = SPEEDS[-1] * 2 * np.pi / 60 * fm.pole_pairs
    print(f"winding {h['turns']} / {h['paths']}, {h['stack_mm']:g} mm, R {h['r_ohm'] * 1e3:.3f} mOhm; pole pairs {fm.pole_pairs}")
    print(f"no-load flux linkage {psi_m * 1e3:.3f} mWb; the maps for this winding: i_d {lo:.1f} ... {hi:.1f} A, |i_q| up to {qmax:.1f} A")
    for v in BUSES:
        print(f"  {v:.0f} V bus: the magnets' line-to-line peak passes it above {conduction_speed(psi_m, fm.pole_pairs, v):,.0f} rpm")
    print(f"  at {SPEEDS[-1]:,} rpm: {w_top * psi_m:.1f} V a phase, {np.sqrt(3) * w_top * psi_m:.1f} V line to line (peak)")
    print(f"characteristic current (psi_d(i_d, 0) = 0): {i_ch:.1f} A")
    print("steady active short circuit:  rpm   i_d [A]  i_q [A]   |i| [A]  torque [N.m]  copper loss [W]")
    rows = table(fm)
    for r in rows:
        flag = "" if r["inside_maps"] and r["converged"] else "  (outside the maps or not converged)"
        print(f"                            {r['rpm']:6d}  {r['i_d']:7.1f}  {r['i_q']:7.1f}  {r['i']:8.1f}  {r['torque']:11.3f}  {r['p_cu']:14.1f}{flag}")
    ok_a = all(abs(r["torque"] / r["torque_from_loss"] - 1) < 1e-9 for r in rows)
    ok_b = all(r["residual_v"] < 1e-6 for r in rows)
    ok_c = abs(rows[-1]["i"] / abs(i_ch) - 1) < 0.01
    print(f"checks: (a) torque = -P_cu / w_m {'PASS' if ok_a else 'FAIL'}; (b) residual below 1e-6 V {'PASS' if ok_b else 'FAIL'}; "
          f"(c) top-speed current within 1 % of the characteristic current {'PASS' if ok_c else 'FAIL'}")


if __name__ == "__main__":
    main()
