"""optimization.maps: the optimized motor's flux-linkage maps and its iron loss, as callable functions: psi_d, psi_q, torque and ripple at any current
(i_d, i_q), for any winding and stack; the inductances and the MTPA path derived from them; and iron_loss(i_d, i_q, rpm, stack_mm).

What it is, and why
    A motor's torque and voltage are not constants times the current: the steel saturates, and the two axes share steel, so the d flux depends on the q current
    (cross-saturation). What describes the motor everywhere is the pair of maps psi_d(i_d, i_q) and psi_q(i_d, i_q). FEMM computed them once on a grid of the
    current plane (maps_run.py, flux stage; maps_report.py averages each point over the 15 rotor positions of one ripple period and writes flux_map.npz). This file
    puts a smooth interpolation over that grid and the symmetry of the rotor, and gives the maps to the torque-speed map (phase 2), the MTPA tables and the
    Simulink export (maps_matlab.py) through one interface.
    The iron loss of the same motor is iron_loss.py's map; its function is re-exported here (maps.iron_loss) so that both live at one address.

The conventions (every number is FEMM 5.0 on the FINE mesh; a 2D model 20 mm deep; the 10-turn placeholder winding; the optimized design of final.json)
    d is the magnet axis, so L_d < L_q (the saliency L_q / L_d is about 1.5). Currents are A PEAK in the amplitude-invariant d-q frame (the length of (i_d, i_q) is the
    phase peak current); flux linkages are Wb (peak, per phase); torque is N.m, positive for i_q > 0; the angle of the current is gamma from the q axis toward -d:
    i_d = -I sin(gamma), i_q = I cos(gamma), which is the convention of design.py, of the MTPA note and of Simscape's "current advance angle".
    The mirror symmetry of the rotor about the d axis would give, for the position-averaged maps, psi_d(i_d, -i_q) = psi_d, psi_q(i_d, -i_q) = -psi_q, T(i_d, -i_q) = -T.
    It was the plan to stop the grid at i_q = 0 and mirror it; the check (maps_report.py, 4.) failed its declared limits: over the whole plane the asymmetry has a median of
    0.2 % (psi_q) and 0.2 % (T) up to the rated current and 0.4 ... 0.5 % at 2 ... 3 times it, with single points at 1 ... 2 % in psi_q and up to 3 % in T (the winding of 36 slots
    in 8 poles is not mirror-symmetric about the d axis, and saturation makes it show; the ripple differs by up to 16 %),
    so FEMM solved the negative half too and the table is TWO-SIDED: nothing is mirrored (a table with only i_q >= 0 is still accepted and mirrored, for old files).
    The grid spans i_d from -3 to +1 times the rated current I_R = 11.8627 A (peak, at 10 turns a coil side) and i_q from -3 to +3 times it, a step of I_R / 4.
    Asking for a current outside it raises ValueError (extrapolating a saturating machine is how wrong torque maps are made); fill="clip" clamps instead.

The winding and the stack (the "scaling rule", which the winding study verified bit for bit: 5, 10 and 20 turns at 20, 10 and 5 A gave identical torque and field)
    FEMM multiplies the current of a block by its turns, so the field depends on the AMPERE-TURNS only. A winding with N turns a coil side in `paths` parallel
    paths has series turns n = 12 N / paths a phase and carries a phase current i; its ampere-turns per coil side are N i / paths. Against the placeholder
    (10 turns, one path, 120 series turns) the scale is  s = N / (10 paths):  its current i is the placeholder's s i, its flux linkage s times the placeholder's:
        psi(i_d, i_q) = s psi_0(s i_d, s i_q)        T(i_d, i_q) = T_0(s i_d, s i_q)        L = s^2 L_0
    A stack of stack_mm instead of the model's 20 mm multiplies psi, T, L and the ripple by stack_mm / 20 (the 2D model has no end effects; a real motor has end
    windings and end leakage on top of that, which no 2D map holds).

The interpolation
    A bicubic spline through the grid (scipy RectBivariateSpline, s = 0) on the mirrored table, so that its derivatives, the inductances, are smooth and the odd /
    even symmetry holds exactly (psi_q and T are zero at i_q = 0). method="linear" gives the bilinear interpolation (the MATLAB side's interp2 'linear', the block's own
    default): the same grid, not the same values between the nodes. The accuracy of the cubic spline between the nodes is tested against 12 off-grid FEMM points
    in maps_report.py (the numbers are in its report and at the end of this docstring).

The calls (arrays in, arrays out; a scalar in gives a float)
    fm = maps.FluxMaps.load()                               # out\\optimization\\maps\\flux_map.npz (maps_report.py writes it)
    fm.psi_dq(i_d, i_q, turns=10, paths=1, stack_mm=20.0)   -> (psi_d, psi_q) [Wb]
    fm.torque(i_d, i_q, ...)                                -> mean torque from FEMM's air-gap integral [N.m] (the study's torque)
    fm.torque_dq(i_d, i_q, ...)                             -> 1.5 p (psi_d i_q - psi_q i_d) from the maps [N.m] (what a controller computes)
    fm.ripple(i_d, i_q, ...)                                -> peak-to-peak torque ripple over the 15 positions, Fourier-rebuilt (design.fourier_pkpk) [N.m]
    fm.inductances(i_d, i_q, ...)                           -> dict of the apparent L_d, L_q ((psi - psi_PM0) / i) and the incremental L_dd, L_qq and the cross terms
                                                               L_dq = d psi_d / d i_q, L_qd = d psi_q / d i_d [H]
    fm.mtpa(I, ...)                                         -> (gamma [deg], i_d, i_q, torque) of the best angle at the current magnitude I [A peak]
    maps.iron_loss(i_d, i_q, rpm, stack_mm=20.0, turns=10, paths=1, lamination_mm=0.3)   -> {"teeth", "yoke", "rotor", "total"} [W]
    Every winding argument is `turns` (per coil side, default 10) and `paths` (parallel, default 1).

Speed of the calls (this PC, one thread, Python 3.13, numpy 2.5.4, scipy 1.18.1; measured by `python -m optimization.maps`, 2026-11-09)
    one scalar call: psi_dq 15 us, torque 10 us, torque_dq 20 us, ripple 11 us, inductances 40 us, iron_loss 36 us, mtpa 117 us (a 361-angle search for one current);
    2,000 points in one call: psi_dq 0.18 us a point, torque 0.09, torque_dq 0.18, inductances 0.46, iron_loss 0.76.

Accuracy (maps_report.py flux, the numbers of 2026-11-09)
    against 12 FEMM points on no node of the grid: psi within 0.014 % of |psi|, torque within 0.67 mN.m (a bilinear interpolation of the same grid: 0.41 % in psi; a grid of twice the step: psi
    0.15 % there and 0.35 % at the nodes it did not use, which would have missed the declared 0.25 %: the reason for the step of a quarter of the rated current). The d-q torque agrees with
    FEMM's air-gap torque to 0.006 % at the median and 0.085 % at the worst of the 425 grid points. Reciprocity of the cross terms: median 0.015 %, largest 0.58 % of sqrt(L_dd L_qq).

Run
    uv run python -m optimization.maps            loads the maps, prints a few values (the rated point, the saliency, the MTPA at 1, 2, 3 times rated) and times each call

Self-check (the answers follow).
  1. Why can the grid stop at i_q = 0, and what would break if the rotor were not symmetric about the d axis?
  2. A winding with 20 turns in one path replaces the 10-turn one. At 5.93 A (half the rated current of the placeholder) what is its flux linkage?
  3. The apparent inductance (psi_d - psi_PM0) / i_d and the incremental d psi_d / d i_d differ. Which one belongs in v = R i + L di/dt?
Answers.
  1. For the position-averaged maps the mirror about the d axis maps (i_d, i_q) to (i_d, -i_q) with psi_q and T changing sign, so the negative half is the positive
     half turned over. A rotor without that symmetry (magnets of unequal strength, a skew) would need its own negative half.
  2. s = 20 / 10 = 2: the current 5.93 A of the 20-turn winding is the placeholder's 11.86 A, and its flux linkage is twice that point's placeholder value.
  3. The incremental one: di/dt multiplies d psi / d i = L_inc in the voltage equation; the apparent one is only the flux per ampere on the way there.
"""
from __future__ import annotations

import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline

import paths
from optimization import design
from optimization.iron_loss import IronLossMap, iron_loss  # noqa: F401  (re-exported: maps.iron_loss)

OUT = paths.path("out") / "optimization" / "maps" / "chosen"          # the optimized design's maps (maps_report.py writes flux_map.npz and iron_map.npz here)
FLUX_MAP = OUT / "flux_map.npz"
POLE_PAIRS = 4
DEPTH_MM = 20.0
PLACEHOLDER_TURNS = 10


def winding_scale(turns: float = PLACEHOLDER_TURNS, paths: float = 1) -> float:
    """s = N / (10 paths): the current of the winding is the placeholder's s times it, the flux linkage s times the placeholder's."""
    if not (turns > 0 and paths > 0):
        raise ValueError(f"turns = {turns!r}, paths = {paths!r}: both must be positive")
    return turns / (PLACEHOLDER_TURNS * paths)


class FluxMaps:
    """psi_d, psi_q, torque and ripple over the current plane, from FEMM's grid. See the module docstring."""

    def __init__(self, z):
        self.id_A = np.asarray(z["id_A"], float)
        iq = np.asarray(z["iq_A"], float)
        self.rated = float(z["rated_current_A"])
        self.rated_gamma = float(z["rated_gamma_deg"])
        self.psi_pm0 = float(z["psi_pm0_Wb"])                           # the no-load flux linkage (mean over the positions), the apparent inductance's reference
        self.pole_pairs = int(z["pole_pairs"])
        self.depth_mm = float(z["depth_mm"])
        self.grid = {k: np.asarray(z[k], float) for k in ("psi_d", "psi_q", "T_age", "T_wst", "T_dq", "ripple_pp")}
        if iq[0] < 0:
            # a two-sided table: FEMM solved both halves of the plane (maps_run's flux and fluxneg stages), nothing is mirrored
            self.iq_full, self.mirrored = iq, False
        else:
            # a half table (i_q >= 0): the negative half is the mirror image (psi_d, ripple even in i_q; psi_q, the torques odd)
            parity = {"psi_d": +1, "psi_q": -1, "T_age": -1, "T_wst": -1, "T_dq": -1, "ripple_pp": +1}
            self.iq_full, self.mirrored = np.concatenate([-iq[:0:-1], iq]), True
            self.grid = {k: np.concatenate([parity[k] * g[:, :0:-1], g], axis=1) for k, g in self.grid.items()}
        self.iq_A = self.iq_full[self.iq_full >= 0]
        self._cubic, self._linear = {}, {}
        for k, g in self.grid.items():
            self._cubic[k] = RectBivariateSpline(self.id_A, self.iq_full, g, kx=3, ky=3, s=0)
            self._linear[k] = RectBivariateSpline(self.id_A, self.iq_full, g, kx=1, ky=1, s=0)
        self.range_d = (float(self.id_A[0]), float(self.id_A[-1]))
        self.range_q = float(np.abs(self.iq_full).max())

    @classmethod
    def load(cls, path: Path | None = None) -> "FluxMaps":
        z = np.load(path or FLUX_MAP)
        return cls(z)

    # ---------------------------------------------------------------------------------------------- the placeholder's current, from a winding's
    def _equivalent(self, i_d, i_q, turns, paths, fill):
        s = winding_scale(turns, paths)
        i_d, i_q = np.broadcast_arrays(np.asarray(i_d, float) * s, np.asarray(i_q, float) * s)
        eps = 1e-9 * self.rated
        outside = (i_d < self.range_d[0] - eps) | (i_d > self.range_d[1] + eps) | (np.abs(i_q) > self.range_q + eps)
        if outside.any():
            if fill == "raise":
                raise ValueError(f"maps: {int(outside.sum())} point(s) outside the solved region (i_d {self.range_d[0] / s:.2f} ... {self.range_d[1] / s:.2f} A, "
                                 f"|i_q| up to {self.range_q / s:.2f} A for this winding); pass fill='clip' to clamp them")
            i_d, i_q = np.clip(i_d, *self.range_d), np.clip(i_q, -self.range_q, self.range_q)
        return s, i_d, i_q

    def _eval(self, key, i_d0, i_q0, method, dx=0, dy=0):
        if method not in ("cubic", "linear"):
            raise ValueError(f"method = {method!r}: 'cubic' or 'linear'")
        spl = (self._cubic if method == "cubic" else self._linear)[key]
        return spl.ev(i_d0.ravel(), i_q0.ravel(), dx=dx, dy=dy).reshape(i_d0.shape)

    @staticmethod
    def _out(value, *inputs):
        return float(value) if all(np.ndim(v) == 0 for v in inputs) else value

    # ---------------------------------------------------------------------------------------------- the maps
    def psi_dq(self, i_d, i_q, *, turns: float = PLACEHOLDER_TURNS, paths: float = 1, stack_mm: float | None = None, method: str = "cubic", fill: str = "raise"):
        """(psi_d, psi_q) [Wb] at the current (i_d, i_q) [A peak] of the winding (turns a coil side, paths in parallel) and stack (default the model's 20 mm)."""
        s, a, b = self._equivalent(i_d, i_q, turns, paths, fill)
        k = (stack_mm or self.depth_mm) / self.depth_mm
        pd, pq = (self._eval(n, a, b, method) for n in ("psi_d", "psi_q"))
        return self._out(s * k * pd, i_d, i_q), self._out(s * k * pq, i_d, i_q)

    def torque(self, i_d, i_q, *, turns: float = PLACEHOLDER_TURNS, paths: float = 1, stack_mm: float | None = None, method: str = "cubic", fill: str = "raise"):
        """Mean torque [N.m] from FEMM's air-gap integral (the study's torque), averaged over the 15 positions of a ripple period."""
        _, a, b = self._equivalent(i_d, i_q, turns, paths, fill)
        k = (stack_mm or self.depth_mm) / self.depth_mm
        return self._out(k * self._eval("T_age", a, b, method), i_d, i_q)

    def torque_dq(self, i_d, i_q, *, turns: float = PLACEHOLDER_TURNS, paths: float = 1, stack_mm: float | None = None, method: str = "cubic", fill: str = "raise"):
        """1.5 p (psi_d i_q - psi_q i_d) from the interpolated flux maps [N.m]: the torque a controller or the Simulink block computes."""
        psi_d, psi_q = self.psi_dq(i_d, i_q, turns=turns, paths=paths, stack_mm=stack_mm, method=method, fill=fill)
        return self._out(1.5 * self.pole_pairs * (np.asarray(psi_d) * np.asarray(i_q, float) - np.asarray(psi_q) * np.asarray(i_d, float)), i_d, i_q)

    def ripple(self, i_d, i_q, *, turns: float = PLACEHOLDER_TURNS, paths: float = 1, stack_mm: float | None = None, method: str = "cubic", fill: str = "raise"):
        """Peak-to-peak torque ripple [N.m] over the 15 positions of one ripple period (design.fourier_pkpk of FEMM's air-gap torque)."""
        _, a, b = self._equivalent(i_d, i_q, turns, paths, fill)
        k = (stack_mm or self.depth_mm) / self.depth_mm
        return self._out(k * np.maximum(self._eval("ripple_pp", a, b, method), 0.0), i_d, i_q)

    # ---------------------------------------------------------------------------------------------- the inductances
    def inductances(self, i_d, i_q, *, turns: float = PLACEHOLDER_TURNS, paths: float = 1, stack_mm: float | None = None, fill: str = "raise") -> dict:
        """The apparent and incremental inductances [H] at (i_d, i_q) (cubic spline).
        apparent    L_d = (psi_d - psi_PM0) / i_d, L_q = psi_q / i_q, psi_PM0 the no-load flux linkage; at a zero current the limit (the incremental value)
        incremental L_dd = d psi_d / d i_d, L_qq = d psi_q / d i_q and the cross terms L_dq = d psi_d / d i_q, L_qd = d psi_q / d i_d.
        The cross terms of a conservative machine are equal (Maxwell's reciprocity): their difference is a check of the maps."""
        s, a, b = self._equivalent(i_d, i_q, turns, paths, fill)
        k = (stack_mm or self.depth_mm) / self.depth_mm
        sp = self._cubic
        d = lambda key, dx, dy: sp[key].ev(a.ravel(), b.ravel(), dx=dx, dy=dy).reshape(a.shape)
        L_dd, L_dq = d("psi_d", 1, 0), d("psi_d", 0, 1)
        L_qd, L_qq = d("psi_q", 1, 0), d("psi_q", 0, 1)
        psi_d, psi_q = d("psi_d", 0, 0), d("psi_q", 0, 0)
        with np.errstate(divide="ignore", invalid="ignore"):
            Ld_app = np.where(np.abs(a) > 1e-9, (psi_d - self.psi_pm0) / a, L_dd)
            Lq_app = np.where(np.abs(b) > 1e-9, psi_q / b, L_qq)
        scale = s * s * k                                                # L = s^2 L_0, and the stack's factor
        out = {"L_d_app": Ld_app, "L_q_app": Lq_app, "L_dd": L_dd, "L_qq": L_qq, "L_dq": L_dq, "L_qd": L_qd}
        return {key: self._out(scale * v, i_d, i_q) for key, v in out.items()}

    # ---------------------------------------------------------------------------------------------- the MTPA path
    def mtpa(self, current, *, turns: float = PLACEHOLDER_TURNS, paths: float = 1, stack_mm: float | None = None, method: str = "cubic", source: str = "dq",
             step_deg: float = 0.25) -> tuple:
        """The best current angle at each magnitude `current` [A peak of the winding]: (gamma [deg], i_d, i_q, torque). The angle is searched over 0 ... 90 degrees on a
        grid of step_deg and refined by the parabola through the best sample and its two neighbours (accurate to a few thousandths of a degree). source="dq" maximises
        1.5 p (psi_d i_q - psi_q i_d) of the maps (what the controller sees); "age" the air-gap torque map."""
        I = np.atleast_1d(np.asarray(current, float))
        gam = np.arange(0.0, 90.0 + step_deg / 2, step_deg)
        g = np.radians(gam)
        i_d = -I[:, None] * np.sin(g)[None, :]
        i_q = I[:, None] * np.cos(g)[None, :]
        if source == "dq":
            T = self.torque_dq(i_d, i_q, turns=turns, paths=paths, stack_mm=stack_mm, method=method)
        elif source == "age":
            T = self.torque(i_d, i_q, turns=turns, paths=paths, stack_mm=stack_mm, method=method)
        else:
            raise ValueError("source is 'dq' or 'age'")
        k = np.clip(np.argmax(T, axis=1), 1, len(gam) - 2)
        rows = np.arange(len(I))
        y0, y1, y2 = T[rows, k - 1], T[rows, k], T[rows, k + 1]
        denom = y0 - 2 * y1 + y2
        shift = np.where(np.abs(denom) > 0, 0.5 * (y0 - y2) / np.where(denom == 0, 1, denom), 0.0)
        best = gam[k] + np.clip(shift, -1, 1) * step_deg
        tb = y1 - 0.25 * (y0 - y2) * np.clip(shift, -1, 1)
        gr = np.radians(best)
        out = (best, -I * np.sin(gr), I * np.cos(gr), tb)
        return tuple(v if np.ndim(current) else float(v[0]) for v in out)

    def mtpa_angle(self, current, **kw):
        """The MTPA angle [deg] at the current magnitude (see mtpa)."""
        return self.mtpa(current, **kw)[0]


# ---------------------------------------------------------------------------------------------- the command line: a look and a timing
def _time(fn, n: int, *args, **kw) -> float:
    t0 = time.perf_counter()
    for _ in range(n):
        fn(*args, **kw)
    return (time.perf_counter() - t0) / n


def main() -> int:
    fm = FluxMaps.load()
    I_R = fm.rated
    i_d, i_q = design.currents(I_R, fm.rated_gamma)
    print(f"flux maps: i_d {fm.range_d[0]:.2f} ... {fm.range_d[1]:.2f} A, |i_q| up to {fm.range_q:.2f} A ({len(fm.id_A)} x {len(fm.iq_full)} grid points, "
          f"{'mirrored from i_q >= 0' if fm.mirrored else 'both halves solved'}), rated {I_R:.4f} A")
    pd, pq = fm.psi_dq(i_d, i_q)
    L = fm.inductances(i_d, i_q)
    print(f"the rated point ({i_d:.3f}, {i_q:.3f}) A: psi_d {1e3 * pd:.3f} mWb, psi_q {1e3 * pq:.3f} mWb, T {fm.torque(i_d, i_q):.4f} N.m (from the flux maps "
          f"{fm.torque_dq(i_d, i_q):.4f}), ripple {fm.ripple(i_d, i_q):.4f} N.m, L_d,app {1e3 * L['L_d_app']:.3f} mH, L_q,app {1e3 * L['L_q_app']:.3f} mH")
    print("   incremental: " + ", ".join(f"{k} {1e3 * v:.3f} mH" for k, v in L.items() if k.startswith("L_") and not k.endswith("app")))
    for mult in (1, 2, 3):
        g, a, b, T = fm.mtpa(mult * I_R)
        print(f"MTPA at {mult} x rated ({mult * I_R:.2f} A): gamma {g:.2f} deg, i_d {a:.2f} A, i_q {b:.2f} A, T {T:.3f} N.m")
    n = 2000
    vec = np.linspace(-3 * I_R, I_R, n), np.linspace(0, 3 * I_R, n)
    print("time per call, scalar (microseconds):")
    for name, fn in (("psi_dq", fm.psi_dq), ("torque", fm.torque), ("torque_dq", fm.torque_dq), ("ripple", fm.ripple), ("inductances", fm.inductances)):
        print(f"   {name:12s} {1e6 * _time(fn, 400, i_d, i_q):8.1f}")
    print(f"   {'mtpa':12s} {1e6 * _time(fm.mtpa, 40, 2 * I_R):8.1f}   (one current, a 361-angle search)")
    print(f"time per point, {n} points in one call (microseconds per point):")
    for name, fn in (("psi_dq", fm.psi_dq), ("torque", fm.torque), ("torque_dq", fm.torque_dq), ("inductances", fm.inductances)):
        print(f"   {name:12s} {1e6 * _time(fn, 20, *vec) / n:8.2f}")
    try:
        r = iron_loss(i_d, i_q, 3000.0)
        print(f"iron loss at the rated point, 3000 rpm, 0.3 mm: " + ", ".join(f"{k} {v:.2f} W" for k, v in r.items()))
        print(f"   {'iron_loss':12s} {1e6 * _time(iron_loss, 400, i_d, i_q, 3000.0):8.1f} us scalar; "
              f"{1e6 * _time(iron_loss, 20, vec[0], vec[1], 5000.0) / n:8.2f} us per point in one call of {n}")
    except FileNotFoundError:
        print("(iron_map.npz is not there yet: no iron loss)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
