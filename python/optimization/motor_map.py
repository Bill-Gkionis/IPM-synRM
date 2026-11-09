"""optimization.motor_map: the torque-speed envelope and the loss and efficiency map of the optimized motor, for any winding, stack, bus voltage and current limit.

What it is, and why
    The flux maps (maps.py) say what the motor's steel does at every current vector, and the iron-loss map says what it dissipates. This engine turns them into what a buyer asks of a motor:
    how much torque at each speed, inside the current the winding and its cooling can carry and the voltage the DC bus can make; what that costs in loss; and the efficiency. Every
    input is an argument (the winding, the stack, the bus, the current limit as a number, an array or a function of the speed, the copper model, the mechanical loss), so the final run
    swaps the chosen winding, the full copper loss and the thermal network's limits by changing inputs, not code. The thermal side is motor_map_thermal.py, the stack table's masses
    and inertia motor_map_stack.py, the run with its checks motor_map_run.py and the figures motor_map_figures.py.

The equations, each one built from the one before (d is the magnet axis; currents and voltages are PEAK values of the phase, "amplitude-invariant", so the length of (i_d, i_q) is the
phase's peak current and three phases dissipate 1.5 R i^2)
    1. The flux linkage of a winding is a map psi(i_d, i_q): psi_d falls as the d current cancels the magnets' flux, psi_q rises with i_q and saturates, and each axis feels the other's
       current (cross-saturation). The winding scales it (maps.py): psi = s k psi_0(s i), T = k T_0(s i), s = N / (10 a), k = stack / 20 mm.
    2. Faraday's law, in steady state. The phase voltage is v = R i + d psi / dt. In the frame that turns with the rotor, at the electrical angular speed w_e = p w_m, the flux linkage
       does not change in time, and the turning of the axes adds the speed voltage:
           v_d = R i_d - w_e psi_q          v_q = R i_q + w_e psi_d
       At low speed the resistive drop R i is the voltage; at high speed w_e |psi| is, and it grows with the speed.
    3. Power. The terminal power is P_in = 1.5 (v_d i_d + v_q i_q) = 1.5 R i^2 + w_m T, with the air-gap torque T = 1.5 p (psi_d i_q - psi_q i_d): copper loss plus the
       electromagnetic power. (The torque of the maps' own air-gap integral, T_age, agrees with this one to 0.5 % at worst: check (c2).) R is the winding's DC resistance; the part of
       the copper model's loss that the voltage equation does not see (FullCopper's skin, proximity, PWM-ripple and circulating-current loss, "ac"; zero for DcCopper) is dissipated
       in the winding and drawn from the terminals as well, so the input power is P_in = 1.5 (v_d i_d + v_q i_q) + P_cu,ac.
    4. The limits. |i| <= I_max (what the winding and the cooling carry) and |v| <= V_max = k V_dc / sqrt 3. Space-vector modulation can make any phase-voltage vector inside the
       hexagon whose inscribed circle has the radius V_dc / sqrt 3 (peak, a phase) without distorting the current; k < 1 leaves room for the inverter's dead time, the drops of the switches
       and cable, a regulator that does not saturate and a bus that sags (Bus: 0.95 here, an assumption).
    5. The three regions. At low speed |v| is small and the best point is MTPA at the current limit. |v| at that point grows with the speed until it reaches V_max, at the CORNER speed,
       the root of (psi_d^2 + psi_q^2) w^2 + 2 R (i_q psi_d - i_d psi_q) w + R^2 i^2 - V_max^2 = 0 (the winding study's formula). Above it the current must turn toward -d (FIELD WEAKENING):
       the best point is where the current circle meets the voltage curve |v| = V_max, whose size shrinks as 1 / speed around the point of zero flux linkage. When the curve lies inside the
       circle more current no longer helps, and the best point is the one of the curve where the torque is largest: MTPV, where the gradient of the torque is parallel to the gradient of |v|^2
       (check (e), with the maps' incremental inductances). Where no point of the circle holds the voltage the torque ends.
    6. The shaft. The air-gap torque acts on the rotor; the iron loss and the mechanical loss (bearings, windage, a shaft fan) brake it:
           T_shaft = T - (P_fe + P_mech) / w_m          P_out = T_shaft w_m          efficiency = P_out / P_in          power factor = P_in / (1.5 |v| |i|)
       (the loss torques are written without the division, so they hold at standstill: the hysteresis drag is a torque that does not vanish with the speed). The balance
       P_in = P_out + P_cu + P_fe + P_mech is check (c1). The power factor is the fundamental's (1.5 (v_d i_d + v_q i_q) over 1.5 |v| |i|), with no harmonics and no iron-loss branch.
    7. Two controls. "MTPA control" gives a requested torque the smallest current that makes it inside the voltage limit (below the corner the maximum torque per ampere, above it the point on
       the voltage curve). "Minimum-loss control" gives it the current vector of least total loss. Their difference is the loss that minimum-loss control saves.

The solver (what the code does, and why it is not a textbook one)
    The current plane is covered by a grid of rays (gamma from the q axis toward -d, 0 ... 90 degrees) and circles (the radius I), with the flux linkages, the torque and the iron-loss coefficients
    of every node computed once (Scan); a speed adds the resistance, the voltage limit and the loss torques. MAXIMUM TORQUE: the best feasible node, the best of the current circle itself (3,601
    angles: a feasible lens at the limit can be thinner than a ray spacing), then a zoom of exact evaluations (12 levels, each five times narrower). When the best point sits on the voltage
    curve inside the circle (MTPV) the zoom stops short, and a second zoom runs ALONG the curve: the torque is flat along it, so a rectangular lattice picks the node nearest the curve whatever
    its place on it (this was found by check (e), which failed on the first run with the lattice alone: sin(angle) 0.2 against its limit 1e-3; the polish brought it to 1e-8, and the
    limit stayed as it was). LEVEL SETS (a target torque): the curve of constant torque is crossed by the rays and by the circles; deep in field weakening it runs along a ray, so only the
    circles see the thin feasible lens, and near MTPA it is tangent to a circle, so only the rays see it; every crossing is evaluated exactly, the best feasible one is projected onto the exact
    level set and two zooms refine it. The solvers are checked against brute force (m, k).

The inputs and where they are swapped (the run's own inputs are the head of motor_map_run.py)
    Winding(turns N, paths a, strands n_s, wire d, overall d_over): the winding study's windings are motor_map_stack.OPTIMIZED_WINDINGS (9 / 4, 7 x 0.56 mm on 60 mm; 26 / 4, 8 x 0.315 mm on
    20 mm ...); PROVISIONAL_WINDING is the first study's example, 26 turns, 4 paths, 8 x 0.315 mm, which the demonstration below uses; scaled(f) is the family of f x the turns on the same copper.
    Bus(v_dc, k).  Machine(winding, stack_mm, copper=, mech=, lamination_mm=0.3).  The current limit and the copper temperature of envelope() are a number, an array or a function of the speed.
    copper: DcCopper() (the first study's DC loss and its 61.15 mm end turn: R = rho(T) l_mt n / (a n_s S), l_mt = 2 L + 2 l_end) or FullCopper(arrangement="cluster") (copper_loss.py: skin, proximity, PWM ripple,
    circulating currents, its diamond end turn): the same three methods, so the swap is a changed argument; its AC part is charged to the input power (the voltage equation sees only R). mech: Mechanical(fan=a function of the speed): bearings and windage as thermal.py has them, plus the
    shaft fan's power (thermal.shaft_fan_power) as a load on the shaft. The thermal limits come from motor_map_thermal.py.

The calls
    m = Machine(winding, stack_mm)           m.point(i_d, i_q, rpm, bus, t_c) -> every number at one or many operating points (the single source of the outputs and the checks)
    sc = Scan(m, i_max)                      sc.max_torque(ctx, i_lim), sc.level_set(ctx, i_lim, targets), with ctx = sc.ctx(rpm, bus, t_c, basis)  (basis "shaft" or "air_gap")
    envelope(m, bus, rpms, i_limit, t_c=...) -> the envelope: arrays over the speed, and the corner speed, its torque and the speed where the torque ends in .meta
    efficiency_map(m, bus, rpms, fractions, i_limit)   -> efficiency, loss and operating point over (speed, fraction of the envelope's torque) for the two controls, and "saving"
    uv run python -m optimization.motor_map          a short demonstration

Numbers to expect (this PC, 2026-11-09; `uv run python -m optimization.motor_map`): the example winding on 60 mm (R 0.1717 Ohm at 100 C, fill 0.448, the maps ending at 54.8 A) at 6 A/mm2 rms,
    21.16 A peak, on a 48 V bus with k = 0.95 (V_max 26.33 V): corner speed 944.0 rpm at 6.26 N.m (the shaft torque; 6.33 N.m at the air gap), the torque ends at 4,982 rpm, 6.28 N.m at standstill
    (the hysteresis drag takes 0.05); efficiency 73.5 % at 500 rpm and 84 ... 85 % between 1,000 and 2,500 rpm. A sweep of 36 speeds takes under a second. The whole run, its checks
    and their limits: motor_map_run.py.

What is not in it: the magnets' eddy loss, the inverter's loss and its current ripple, the end-winding leakage (not in the 2D maps), demagnetisation at large negative i_d (the maps use
N45 at its nominal remanence), the AC copper loss on the DC model, any current beyond 3 times the placeholder's rated current (the maps end there: asking for more raises an error or the
limit is clipped and flagged), and the harmonics of the torque (the maps hold its mean over a ripple period).

Self-check (the answers follow).
  1. Why does the same winding on a twice longer stack have half the corner speed and the same power at the corner?
  2. At 16,330 rpm the open-circuit voltage of the example winding on 60 mm is several times the bus. What lets the motor run there at all, and what limits how far it can go?
  3. Minimum-loss control saves almost nothing above the corner. Why?
Answers.
  1. The flux linkage goes with the length, so the speed at which w_e |psi| reaches V_max goes as 1 / length; the torque at a given current goes as the length, so torque x speed at the
     corner does not change.
  2. Field weakening: a current along -d cancels the magnets' flux until the flux linkage, and so the speed voltage, is back under V_max. How far it can go is the current that can be
     spared for that: the voltage curve shrinks around the point where the flux linkage is zero, at i_d of about -psi_PM / L_d, and the machine can run only while that curve touches the current circle.
  3. Above the corner the voltage curve leaves one current vector for a torque (the circle and the curve meet in a point, or the torque is at the top of the curve), so there is nothing to choose.
"""




from __future__ import annotations

import json
import math
from dataclasses import dataclass, replace

import numpy as np
from scipy.interpolate import RectBivariateSpline

import paths
from optimization import design
from optimization import winding as winding_study
from optimization.iron_loss import IronLossMap, eddy_scale
from optimization.maps import FluxMaps

PP = 4                                  # pole pairs
TOP_RPM = 16_330.0                      # the rotor's safe speed (the stress check; a design decision of 2026-11-09)
COILS_PER_PHASE = 12                    # 36 slots, one coil a slot, three phases
SQRT3 = math.sqrt(3.0)
C_SPEED = PP / (2.0 * math.pi)          # the electrical frequency is C_SPEED x the mechanical angular speed [Hz per rad/s]
ROW_MM = 20.0                           # one row of magnets and its share of the stack
V_TOL = 1e-5                            # a constraint is "active" when the solution sits within this relative distance of it
I_TOL = 1e-4                            # (the current limit: a point within 1e-4 of the circle is on it: the zoom may stop a hair inside)
SOLVER_LEVELS = 12                      # zoom levels of the maximum-torque refinement (each one shrinks the search window five-fold: 12 levels reach 1e-9 of the first window)
LEVELSET_LEVELS = 6                     # zoom levels of the level-set refinement (21 trial points a level, so the window shrinks five-fold too)
TOP_K = 6                               # candidates per target that the level set projects onto the exact level set
V_SLACK = 1e-9                          # a level-set solution may sit ON the voltage boundary: |v|^2 within this relative slack of V_max^2 counts as inside (check (o) holds it to 1e-8)

_CACHE: dict = {}


def _cached(key, make):
    if key not in _CACHE:
        _CACHE[key] = make()
    return _CACHE[key]


def load_maps() -> FluxMaps:
    """The flux maps of the optimized design (out/optimization/maps/chosen/flux_map.npz), loaded once."""
    return _cached("maps", FluxMaps.load)


def load_iron() -> IronLossMap:
    """The iron-loss map of the optimized design (iron_map.npz), loaded once."""
    return _cached("iron", IronLossMap)


def final_overrides() -> dict:
    """The optimized design: the "overrides" of out/optimization/final/final.json (motorgeo.Params fields on top of the band-free rotor)."""
    def read():
        with open(paths.path("out") / "optimization" / "final" / "final.json", encoding="utf-8") as f:
            return json.load(f)["overrides"]
    return _cached("overrides", read)


def design_params(overrides: dict | None = None):
    """motorgeo.Params of the design (default: the optimized one)."""
    ov = final_overrides() if overrides is None else overrides
    return _cached(("params", json.dumps(ov, sort_keys=True)), lambda: design.params(ov))


def layer_area_mm2(p) -> float:
    """One coil layer of the slot [mm^2] (motorgeo.check's coil_layer_area_mm2: 36.218 for the optimized design)."""
    import motorgeo as mg
    return _cached(("layer", id(p)), lambda: float(mg.check(p, mg.build(p))[0]["coil_layer_area_mm2"]))


def end_turn_stand_in(p) -> dict:
    """The first study's end-turn geometry (winding.end_turn_geometry): the half circle on the coil's span plus two legs, 61.15 mm for the optimized design."""
    return _cached(("endturn", id(p)), lambda: winding_study.end_turn_geometry(p))


# ====================================================================================================== the inputs
@dataclass(frozen=True)
class Winding:
    """A winding of the optimized stator: turns a coil side N, parallel paths a, strands in hand n_s of a bare diameter d (mm) and overall diameter d_over.
    Series turns a phase n = 12 N / a (a double-layer winding has one coil a slot, 12 coils a phase). turns and wire may be fractional: the "what if" family of
    scaled() keeps the copper fill but is not a real winding (the final run takes real ones from the winding study's enumeration)."""
    turns: float
    paths: int = 1
    strands: float = 1.0
    wire_mm: float = 0.0
    wire_over_mm: float = 0.0

    @property
    def series_turns(self) -> float:
        return COILS_PER_PHASE * self.turns / self.paths

    @property
    def strand_area_m2(self) -> float:
        return math.pi * (self.wire_mm * 1e-3) ** 2 / 4.0

    @property
    def copper_area_m2(self) -> float:
        """The copper of one turn (n_s strands) [m^2]."""
        return self.strands * self.strand_area_m2

    def fill(self, layer_area_mm2: float) -> float:
        """The copper fill of a coil layer, N n_s S / A_layer."""
        return self.turns * self.copper_area_m2 * 1e6 / layer_area_mm2

    def current_from_density(self, j_rms_a_mm2: float) -> float:
        """The phase current [A peak] at which the copper carries j_rms [A/mm^2 rms]: i = sqrt 2 J a n_s S (every strand carries i / (a n_s))."""
        return math.sqrt(2.0) * j_rms_a_mm2 * 1e6 * self.paths * self.copper_area_m2

    def density_from_current(self, i_peak: float) -> float:
        """The copper's current density [A/mm^2 rms] at the phase current i_peak [A]."""
        return i_peak / (math.sqrt(2.0) * 1e6 * self.paths * self.copper_area_m2)

    def scaled(self, f: float) -> "Winding":
        """The same slot filling with f times the turns: N f turns of wire of 1 / f the copper area (so N n_s S, the fill and the copper in the slot are unchanged).
        The family the stack table uses to hold the corner speed; its current limit (a fixed current density) goes as 1 / f."""
        r = math.sqrt(f)
        return replace(self, turns=self.turns * f, wire_mm=self.wire_mm / r, wire_over_mm=self.wire_over_mm / r)

    def label(self) -> str:
        return f"{self.turns:g} turns, {self.paths} paths, {self.strands:g} x {self.wire_mm:g} mm"


PROVISIONAL_WINDING = Winding(26, 4, 8, 0.315, 0.367)       # the winding study's example (Optimization.md 9.5): 78 series turns, fill 0.448


@dataclass(frozen=True)
class Bus:
    """The DC bus and the share of the space-vector modulation's linear limit the controller may use: V_max = k V_dc / sqrt 3 (the largest phase voltage, peak, the
    inverter can make without overmodulating). k = 0.95 here: 1.7 % for the inverter's dead time (about 0.5 us at 20 kHz), 0.3 % for the conduction drops of the switches
    and the cable, and the rest, 3 %, a reserve so that the current regulators never saturate and for the bus sagging under load. An ASSUMPTION, an input like any other."""
    v_dc: float = 48.0
    k: float = 0.95
    limit: float | None = None          # an explicit voltage limit [V] instead of k V_dc / sqrt 3: the thermal rating takes the limit away (1e9 V) and keeps v_dc for the PWM ripple

    @property
    def v_max(self) -> float:
        return self.limit if self.limit is not None else self.k * self.v_dc / SQRT3


@dataclass(frozen=True)
class Mechanical:
    """The mechanical loss, a rough estimate (bearings and windage), with the formulas and constants of thermal.py so that the network and the efficiency use one number.
    bearings: M = 0.5 mu F d for each of two deep-groove balls (mu 0.0015, F 60 N, bore 15 mm), P = 2 M w.
    windage: P = C_f pi rho w^3 r^4 L, the rotor's skin friction in the air gap (Vrancik's form), C_f the larger of the laminar Couette value 2 / Re_d and the vortex values of
    Bilgen and Boulos (0.515 (d/r)^0.3 / Re_d^0.5 up to Re_d = 1e4, 0.0325 (d/r)^0.3 / Re_d^0.2 above), Re_d = w r d / nu, air at the film temperature; windage_factor 1.2-2 for a
    slotted stator. The constants are textbook values, not looked up again here: expect a factor of two either way, which at the top speed is 5 to 10 % of the loss."""
    bearing_mu: float = 0.0015
    bearing_load_n: float = 60.0
    bearing_bore_m: float = 0.015
    windage_factor: float = 1.0
    film_c: float = 70.0
    fan: object = None                  # a shaft-mounted cooling fan: a function rpm -> shaft power [W] (None: no fan on the shaft). It is a load like the windage

    def bearing_torque(self) -> float:
        return 2.0 * 0.5 * self.bearing_mu * self.bearing_load_n * self.bearing_bore_m

    def air(self) -> tuple[float, float]:
        """(density, kinematic viscosity) of dry air at the film temperature (ideal gas, Sutherland's law)."""
        t = self.film_c + 273.15
        rho = 101325.0 / (287.05 * t)
        return rho, 1.458e-6 * t ** 1.5 / (t + 110.4) / rho

    def fan_power(self, rpm):
        """The shaft power [W] of the fan at the speed rpm (arrays allowed); zero without a fan."""
        r = np.asarray(rpm, float)
        if self.fan is None:
            return np.zeros_like(r)
        try:
            out = np.asarray(self.fan(r), float)
            if out.shape == r.shape:
                return out
        except (TypeError, ValueError):
            pass
        return np.array([float(self.fan(x)) for x in r.ravel()]).reshape(r.shape)

    def power(self, rpm, r_ro_m: float, gap_m: float, length_m: float) -> dict:
        """{"bearing", "windage", "fan", "total"} [W] at the speed rpm (arrays allowed)."""
        w = 2.0 * math.pi * np.asarray(rpm, float) / 60.0
        rho, nu = self.air()
        re = w * r_ro_m * gap_m / nu
        with np.errstate(divide="ignore", invalid="ignore"):
            cf_lam = np.where(re > 0, 2.0 / np.where(re > 0, re, 1.0), 0.0)
            cf_mid = np.maximum(cf_lam, 0.515 * (gap_m / r_ro_m) ** 0.3 / np.sqrt(np.where(re > 0, re, 1.0)))
            cf_hi = np.maximum(cf_lam, 0.0325 * (gap_m / r_ro_m) ** 0.3 / np.where(re > 0, re, 1.0) ** 0.2)
        cf = np.where(re < 500.0, cf_lam, np.where(re < 1.0e4, cf_mid, cf_hi))
        wind = np.where(re > 0, self.windage_factor * cf * math.pi * rho * w ** 3 * r_ro_m ** 4 * length_m, 0.0)
        bear = self.bearing_torque() * w
        fan = self.fan_power(rpm)
        return {"bearing": bear, "windage": wind, "fan": fan, "total": bear + wind + fan}


# ====================================================================================================== the copper models (the swap point)
class DcCopper:
    """The copper loss as DC loss with the first study's end turn: P = 1.5 R i^2 (three phases of peak current i), R = rho(T) l_mt n / (a n_s S), n the series turns,
    the mean turn l_mt = 2 L + 2 l_end, l_end = (pi / 2) span + 2 e = 61.15 mm for the optimized design (the winding study's stand-in: winding.end_turn_geometry). This is what the
    study priced and what the thermal network assumes by default. Skin, proximity, PWM-ripple and circulating-current loss are not in it ("ac" is empty);
    FullCopper below adds them. Both have the same three methods, which is what the engine calls."""
    name = "DC (first study's end turn)"

    def __init__(self, l_end_mm: float | None = None):
        self.l_end_mm = l_end_mm

    def end_mm(self, m: "Machine") -> float:
        return self.l_end_mm if self.l_end_mm is not None else end_turn_stand_in(m.p)["l_end_mm"]

    def mean_turn_m(self, m: "Machine") -> float:
        return 2.0 * m.L + 2.0 * self.end_mm(m) * 1e-3

    def resistance(self, m: "Machine", t_c):
        """The phase's DC resistance [Ohm] with the copper at t_c [degC]."""
        return winding_study.resistivity(t_c) * self.mean_turn_m(m) * m.w.series_turns / (m.w.paths * m.w.copper_area_m2)

    def slot_fraction(self, m: "Machine") -> float:
        """The share of the copper (and of its loss) that lies in the slots: L / (L + l_end), as the thermal network splits it."""
        return m.L / (m.L + self.end_mm(m) * 1e-3)

    def loss(self, m: "Machine", i_d, i_q, rpm, t_c, v_dc) -> dict:
        p = 1.5 * self.resistance(m, t_c) * (np.asarray(i_d, float) ** 2 + np.asarray(i_q, float) ** 2)
        return {"dc": p, "ac": np.zeros_like(p), "total": p}


class FullCopper:
    """The copper loss of copper_loss.py (the other module): DC in the slots and the end turns (its diamond end turn, 48.55 mm, unless end_model="arc"), skin effect,
    proximity loss of the machine's own field, the loss of the inverter's PWM ripple current and, for strands that are not transposed, circulating currents. Needs the
    switching frequency (the ripple) and the bus voltage (the modulation index). "ac" is everything that is not DC. See copper_loss.py for what each part assumes."""
    name = "full (copper_loss.py)"

    def __init__(self, f_sw: float = 20_000.0, end_model: str = "diamond", arrangement: str = "transposed", ripple_layers: str = "same",
                 end_ext_mm: float | None = None, end_angle_deg: float = 30.0):
        self.f_sw, self.end_model, self.arrangement, self.ripple_layers = f_sw, end_model, arrangement, ripple_layers
        self.end_ext_mm, self.end_angle_deg = end_ext_mm, end_angle_deg

    def _call(self, m: "Machine", i_d, i_q, rpm, t_c, v_dc, f_sw):
        from optimization import copper_loss
        i_d, i_q = np.asarray(i_d, float), np.asarray(i_q, float)
        i = np.hypot(i_d, i_q)
        gam = np.degrees(np.arctan2(-i_d, i_q))
        w = m.w
        return copper_loss.losses_arrays(w.turns, w.paths, w.strands, w.wire_mm, w.wire_over_mm, i, gam, rpm, f_sw, v_dc, t_c, m.stack_mm,
                                         arrangement=self.arrangement, end_model=self.end_model, end_ext_mm=self.end_ext_mm, end_angle_deg=self.end_angle_deg,
                                         ripple_layers=self.ripple_layers)

    def resistance(self, m: "Machine", t_c):
        return self._call(m, 0.0, 1.0, 0.0, t_c, 48.0, 0.0)["r_dc_ohm"]

    def slot_fraction(self, m: "Machine") -> float:
        r = self._call(m, 0.0, 1.0, 0.0, 100.0, 48.0, 0.0)
        return float(r["r_slot_ohm"] / r["r_dc_ohm"])

    def loss(self, m: "Machine", i_d, i_q, rpm, t_c, v_dc) -> dict:
        r = self._call(m, i_d, i_q, rpm, t_c, v_dc, self.f_sw)
        dc = r["p_dc_slot"] + r["p_dc_end"]
        return {"dc": dc, "ac": r["p_total"] - dc, "total": r["p_total"]}


# ====================================================================================================== the machine
class Machine:
    """A winding of the optimized stator on a stack, on the flux maps: everything that is a function of the current vector (i_d, i_q) [A peak] and the speed.
    The scaling rule of maps.py: with s = N / (10 a) and k = stack / 20 mm, psi = s k psi_0(s i), T = k T_0(s i), the iron loss k P_0(s i); the maps end at 3 times the
    rated current of the 10-turn placeholder, which for this winding is `i_cap` amperes."""

    def __init__(self, w: Winding, stack_mm: float = 60.0, *, overrides: dict | None = None, copper=None, mech: Mechanical | None = None, lamination_mm: float = 0.3):
        self.w, self.stack_mm, self.L = w, float(stack_mm), float(stack_mm) * 1e-3
        self.maps, self.iron_map = load_maps(), load_iron()
        self.p = design_params(overrides)
        self.s = w.turns / (10.0 * w.paths)
        self.k = self.stack_mm / self.maps.depth_mm
        self.lamination_mm = lamination_mm
        self.k_e = eddy_scale(lamination_mm, self.iron_map.fit_gauge_mm)
        self.copper = copper or DcCopper()
        self.mech = mech or Mechanical()
        self.r_ro, self.gap = self.p.r_ro * 1e-3, self.p.airgap * 1e-3
        self.layer_mm2 = layer_area_mm2(self.p)
        self.fill = w.fill(self.layer_mm2)
        self.i_cap = min(-self.maps.range_d[0], self.maps.range_q) / self.s          # the largest current the maps cover [A peak]
        self._iron = _cached(("irontotal", id(self.iron_map)), lambda: self._total_iron_splines())

    def _total_iron_splines(self):
        """Cubic splines of the three loss sums added over teeth, yoke and rotor (a spline is linear in its data, so these equal the sum of the regions' splines)."""
        im = self.iron_map
        S = im.S.sum(axis=0)
        iq_full = np.concatenate([-im.iq_A[:0:-1], im.iq_A])
        return [RectBivariateSpline(im.id_A, iq_full, np.concatenate([S[t][:, :0:-1], S[t]], axis=1), kx=3, ky=3, s=0) for t in range(3)]

    # ---- the electromagnetics
    def flux(self, i_d, i_q):
        """(psi_d, psi_q) [Wb peak] at the current (i_d, i_q) [A peak] of this winding and stack."""
        return self.maps.psi_dq(i_d, i_q, turns=self.w.turns, paths=self.w.paths, stack_mm=self.stack_mm)

    def torque_dq(self, i_d, i_q):
        psi_d, psi_q = self.flux(i_d, i_q)
        return 1.5 * PP * (psi_d * np.asarray(i_q, float) - psi_q * np.asarray(i_d, float))

    def torque_age(self, i_d, i_q):
        """FEMM's air-gap torque (the mean over the 15 positions) from the maps: the second route to the torque."""
        return self.maps.torque(i_d, i_q, turns=self.w.turns, paths=self.w.paths, stack_mm=self.stack_mm)

    def resistance(self, t_c):
        return self.copper.resistance(self, t_c)

    def back_emf_flux(self) -> float:
        """The no-load flux linkage psi_PM [Wb] of this winding and stack (the maps' psi_d at zero current)."""
        return float(self.flux(0.0, 0.0)[0])

    # ---- the losses
    def iron_coeffs(self, i_d, i_q):
        """(c_h, c_e, c_x) with P_fe = c_h f + c_e f^2 + c_x f^1.5 [W], f the electrical frequency [Hz]: the three regions added, the winding, the stack and the CAD's lamination in."""
        a = np.asarray(i_d, float) * self.s
        b = np.abs(np.asarray(i_q, float)) * self.s
        c = np.stack([np.maximum(spl.ev(a.ravel(), b.ravel()).reshape(a.shape), 0.0) for spl in self._iron])
        return self.k * np.stack([c[0], self.k_e * c[1], c[2]])

    def iron(self, i_d, i_q, rpm) -> dict:
        """{"teeth", "yoke", "rotor", "total"} [W]: iron_loss.IronLossMap for this winding, stack and lamination."""
        return self.iron_map(i_d, i_q, rpm, self.stack_mm, turns=self.w.turns, paths=self.w.paths, lamination_mm=self.lamination_mm)

    def mech_power(self, rpm) -> dict:
        return self.mech.power(rpm, self.r_ro, self.gap, self.L)

    def mech_torque(self, w_m):
        """The mechanical loss as a braking torque [N.m] at the angular speed w_m [rad/s]: bearings M, windage P / w (its limit at standstill is zero)."""
        w = np.asarray(w_m, float)
        p = self.mech.power(w * 60.0 / (2.0 * math.pi), self.r_ro, self.gap, self.L)
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.mech.bearing_torque() + np.where(w > 1e-9, (p["windage"] + p["fan"]) / np.where(w > 1e-9, w, 1.0), 0.0)

    def iron_torque(self, i_d, i_q, w_m):
        """The iron loss as a braking torque [N.m]: P_fe / w_m written without the division, so that it holds at standstill (the hysteresis drag): c_h C + c_e C^2 w + c_x C^1.5 sqrt(w)."""
        c = self.iron_coeffs(i_d, i_q)
        w = np.asarray(w_m, float)
        return c[0] * C_SPEED + c[1] * C_SPEED ** 2 * w + c[2] * C_SPEED ** 1.5 * np.sqrt(w)

    def total_loss(self, i_d, i_q, rpm, v_dc, t_c):
        """Copper + iron + mechanical [W] at the current vector and speed."""
        w = 2.0 * math.pi * np.asarray(rpm, float) / 60.0
        c = self.iron_coeffs(i_d, i_q)
        f = C_SPEED * w
        fe = c[0] * f + c[1] * f ** 2 + c[2] * f ** 1.5
        return self.copper.loss(self, i_d, i_q, rpm, t_c, v_dc)["total"] + fe + self.mech_power(rpm)["total"]

    # ---- one operating point, completely: the single source of every output and every check
    def point(self, i_d, i_q, rpm, bus: Bus, t_c=100.0) -> dict:
        """Everything about the steady state at the current vector (i_d, i_q) [A peak], the speed rpm and the copper at t_c (arrays broadcast together):
        voltages, input power, the losses, the shaft torque and power, efficiency and power factor. The conventions are in the module docstring."""
        i_d, i_q, rpm, t_c = np.broadcast_arrays(*(np.asarray(v, float) for v in (i_d, i_q, rpm, t_c)))
        w_m = 2.0 * math.pi * rpm / 60.0
        w_e = PP * w_m
        psi_d, psi_q = self.flux(i_d, i_q)
        R = self.resistance(t_c)
        v_d, v_q = R * i_d - w_e * psi_q, R * i_q + w_e * psi_d
        T_dq = 1.5 * PP * (psi_d * i_q - psi_q * i_d)
        p_fund = 1.5 * (v_d * i_d + v_q * i_q)                                     # what the voltage equation (the DC resistance) sees
        cu = self.copper.loss(self, i_d, i_q, rpm, t_c, bus.v_dc)
        p_in = p_fund + cu["ac"]                                                    # plus the copper loss it does not see: dissipated in the winding, drawn from the terminals
        fe = self.iron(i_d, i_q, rpm)
        me = self.mech_power(rpm)
        t_loss = self.iron_torque(i_d, i_q, w_m) + self.mech_torque(w_m)
        T_shaft = T_dq - t_loss
        p_out = T_shaft * w_m
        i_abs, v_abs = np.hypot(i_d, i_q), np.hypot(v_d, v_q)
        with np.errstate(divide="ignore", invalid="ignore"):
            eta = np.where(p_in > 0, p_out / np.where(p_in > 0, p_in, 1.0), np.nan)
            pf = np.where(i_abs * v_abs > 0, p_fund / (1.5 * np.where(i_abs * v_abs > 0, i_abs * v_abs, 1.0)), np.nan)
        return {"rpm": rpm, "i_d": i_d, "i_q": i_q, "I": i_abs, "gamma": np.degrees(np.arctan2(-i_d, i_q)), "psi_d": psi_d, "psi_q": psi_q, "v_d": v_d, "v_q": v_q, "v_abs": v_abs,
                "T_em": T_dq, "T_age": self.torque_age(i_d, i_q), "T_shaft": T_shaft, "P_em": T_dq * w_m, "P_out": p_out, "P_in": p_in,
                "P_cu": cu["total"], "P_cu_dc": cu["dc"], "P_cu_ac": cu["ac"], "P_fe": fe["total"], "P_fe_teeth": fe["teeth"], "P_fe_yoke": fe["yoke"], "P_fe_rotor": fe["rotor"],
                "P_bearing": me["bearing"], "P_windage": me["windage"], "P_fan": me["fan"], "P_mech": me["total"],
                "P_loss": cu["total"] + fe["total"] + me["total"], "eta": eta, "PF": pf, "R": R, "t_cu": t_c}

    def balance_error(self, pt: dict) -> np.ndarray:
        """|P_in - (P_out + P_cu + P_fe + P_mech)| / P_in: the energy balance of check (c), per point."""
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.abs(pt["P_in"] - (pt["P_out"] + pt["P_cu"] + pt["P_fe"] + pt["P_mech"])) / np.where(pt["P_in"] > 0, pt["P_in"], np.nan)


# ====================================================================================================== the solver
@dataclass
class Ctx:
    """Everything the solver needs at one speed: the speed, the winding's resistance, the squared voltage limit, the torque basis and the loss torques' speed factors."""
    rpm: float
    w_m: float
    w_e: float
    R: float
    v2: float
    v_dc: float
    t_c: float
    basis: str
    a_h: float
    a_e: float
    a_x: float
    t_mech: float


class Scan:
    """The (ray, radius) grid of a machine in the current plane: i_d = -I sin(gamma), i_q = I cos(gamma), gamma 0 ... 90 degrees (a ray every dgam) and I 0 ... i_max, with everything
    that does not depend on the speed computed once: the flux linkages, the air-gap torque and the iron-loss coefficients. A speed adds the resistance, the voltage limit and the loss
    torques (ctx). The solvers read the grid, pick the best point, and refine it by zooming (a window of 21 x 21 exact evaluations around the best, five times smaller each level)."""

    def __init__(self, m: Machine, i_max: float | None = None, n_gamma: int = 361, n_radius: int = 241):
        self.m = m
        self.i_max = min(i_max if i_max is not None else m.i_cap, m.i_cap)
        self.nG, self.nR = n_gamma, n_radius
        self.gam = np.linspace(0.0, 90.0, n_gamma)
        self.dgam = 90.0 / (n_gamma - 1)
        self.I = np.linspace(0.0, self.i_max, n_radius)
        self.dI = self.i_max / (n_radius - 1)
        g = np.radians(self.gam)[:, None]
        self.i_d, self.i_q = -self.I[None, :] * np.sin(g), self.I[None, :] * np.cos(g)
        self.psi_d, self.psi_q = m.flux(self.i_d, self.i_q)
        self.T_em = 1.5 * PP * (self.psi_d * self.i_q - self.psi_q * self.i_d)
        self.c = m.iron_coeffs(self.i_d, self.i_q)

    # ---- one speed
    def ctx(self, rpm: float, bus: Bus, t_c: float = 100.0, basis: str = "shaft") -> Ctx:
        """The solver's context at the speed rpm. basis "shaft": the torque is the shaft torque (the air-gap torque less the iron and mechanical loss torques), the quantity a load
        sees; "air_gap": the air-gap torque of the d-q model, for the checks against the MTPA table of maps.py."""
        if basis not in ("shaft", "air_gap"):
            raise ValueError("basis: 'shaft' or 'air_gap'")
        w_m = 2.0 * math.pi * rpm / 60.0
        shaft = basis == "shaft"
        return Ctx(rpm=rpm, w_m=w_m, w_e=PP * w_m, R=float(self.m.resistance(t_c)), v2=bus.v_max ** 2, v_dc=bus.v_dc, t_c=t_c, basis=basis,
                   a_h=C_SPEED if shaft else 0.0, a_e=C_SPEED ** 2 * w_m if shaft else 0.0, a_x=C_SPEED ** 1.5 * math.sqrt(w_m) if shaft else 0.0,
                   t_mech=float(self.m.mech_torque(w_m)) if shaft else 0.0)

    def grid(self, c: Ctx):
        """(T, v^2) over the grid at this speed."""
        v_d = c.R * self.i_d - c.w_e * self.psi_q
        v_q = c.R * self.i_q + c.w_e * self.psi_d
        T = self.T_em - (self.c[0] * c.a_h + self.c[1] * c.a_e + self.c[2] * c.a_x) - c.t_mech if c.basis == "shaft" else self.T_em
        return T, v_d * v_d + v_q * v_q

    def state(self, I, gam_deg, c: Ctx):
        """(T, v^2) at exact (I, gamma) arrays: the maps' splines, not the grid."""
        g = np.radians(gam_deg)
        i_d, i_q = -I * np.sin(g), I * np.cos(g)
        psi_d, psi_q = self.m.flux(i_d, i_q)
        v_d, v_q = c.R * i_d - c.w_e * psi_q, c.R * i_q + c.w_e * psi_d
        T = 1.5 * PP * (psi_d * i_q - psi_q * i_d)
        if c.basis == "shaft":
            k = self.m.iron_coeffs(i_d, i_q)
            T = T - (k[0] * c.a_h + k[1] * c.a_e + k[2] * c.a_x) - c.t_mech
        return T, v_d * v_d + v_q * v_q

    # ---- the largest torque inside the current circle and the voltage limit
    def max_torque(self, c: Ctx, i_lim: float, grid=None, levels: int = SOLVER_LEVELS, n: int = 21):
        """The point of largest torque with |i| <= i_lim and |v| <= V_max: the best feasible point of the grid, the best of the boundary ray |i| = i_lim, then a zoom.
        Returns None when no point of the plane is feasible (the voltage cannot be held with this current), else a dict with I, gam, T, v2, i_d, i_q and the two flags
        v_active (|v| within 1e-5 of V_max) and i_active (|i| = i_lim)."""
        T, vv = self.grid(c) if grid is None else grid
        i_lim = min(i_lim, self.i_max)
        kmax = min(int(math.floor(i_lim / self.dI + 1e-9)), self.nR - 1)
        Tw = np.where(vv[:, :kmax + 1] <= c.v2, T[:, :kmax + 1], -np.inf)
        j = int(np.argmax(Tw))
        ig, ir = divmod(j, kmax + 1)
        best = (float(self.I[ir]), float(self.gam[ig]), float(Tw[ig, ir]))
        gb = np.linspace(0.0, 90.0, 3601)                                   # the circle |i| = i_lim itself, finely: a thin feasible lens at the limit would slip between grid nodes
        Tb, vb = self.state(np.full(gb.shape, i_lim), gb, c)
        Tb = np.where(vb <= c.v2, Tb, -np.inf)
        jb = int(np.argmax(Tb))
        if Tb[jb] > best[2]:
            best = (float(i_lim), float(gb[jb]), float(Tb[jb]))
        if not np.isfinite(best[2]):
            return None
        I0, g0, T0 = best
        half_i, half_g = 2.0 * self.dI, 2.0 * self.dgam
        lin = np.linspace(-1.0, 1.0, n)
        for _ in range(levels):
            Ig, Gg = np.meshgrid(np.clip(I0 + half_i * lin, 0.0, i_lim), np.clip(g0 + half_g * lin, 0.0, 90.0))
            Tt, vt = self.state(Ig, Gg, c)
            Tt = np.where(vt <= c.v2, Tt, -np.inf)
            j = int(np.argmax(Tt))
            if Tt.flat[j] > T0:
                I0, g0, T0 = float(Ig.flat[j]), float(Gg.flat[j]), float(Tt.flat[j])
            half_i *= 4.0 / (n - 1)
            half_g *= 4.0 / (n - 1)
        Tf, vf = self.state(np.array(I0), np.array(g0), c)
        if np.isfinite(c.v2) and vf >= c.v2 * (1.0 - 1e-3) and I0 < i_lim * (1.0 - I_TOL):
            I0, g0 = self._along_voltage_boundary(c, i_lim, I0, g0)                    # on the voltage boundary and inside the current circle: the MTPV case
            Tf, vf = self.state(np.array(I0), np.array(g0), c)
        g = math.radians(g0)
        return {"I": I0, "gam": g0, "T": float(Tf), "v2": float(vf), "i_d": -I0 * math.sin(g), "i_q": I0 * math.cos(g),
                "v_active": bool(vf >= c.v2 * (1.0 - V_TOL)), "i_active": bool(I0 >= i_lim * (1.0 - I_TOL))}

    # ---- the voltage boundary
    def _radius_v(self, gam, I0, i_lim, c: Ctx, iters: int = 12):
        """The radius at which |v|^2 = V_max^2 on each ray (arrays of one shape), by the secant method from the guess I0: the point of the voltage boundary on that ray. Returns (I, converged)."""
        def f(I):
            return self.state(I, gam, c)[1] - c.v2
        Ia = np.clip(I0, 0.0, i_lim)
        Ib = np.clip(Ia * (1.0 + 1e-3) + 1e-6, 0.0, i_lim)
        fa, fb = f(Ia), f(Ib)
        for _ in range(iters):
            den = fb - fa
            good = np.abs(den) > 1e-300
            step = np.where(good, fb * (Ib - Ia) / np.where(good, den, 1.0), 0.0)
            Ic = np.clip(Ib - step, 0.0, i_lim)
            Ia, fa, Ib, fb = Ib, fb, Ic, f(Ic)
        return Ib, np.abs(fb) <= 1e-9 * c.v2

    def _along_voltage_boundary(self, c: Ctx, i_lim: float, I0: float, g0: float, levels: int = 14, n: int = 21):
        """The best point ON the voltage boundary near (I0, g0): where the torque is largest along the curve |v| = V_max (the MTPV point), by a zoom in the angle with the radius
        solved on the boundary for every trial angle. The 2-D zoom of max_torque cannot do this alone: along the boundary the torque is flat, so the lattice point nearest the boundary
        wins whatever its position along it, and the window then closes on a point that is not the maximum."""
        lin = np.linspace(-1.0, 1.0, n)
        T0 = float(self.state(np.array(I0), np.array(g0), c)[0])
        half = 2.0
        for _ in range(levels):
            gv = np.clip(g0 + half * lin, 0.0, 90.0)
            Iv, conv = self._radius_v(gv, np.full(n, I0), i_lim, c)
            Tv, vv = self.state(Iv, gv, c)
            ok = conv & (vv <= c.v2 * (1.0 + 1e-12)) & (Iv <= i_lim)
            Tv = np.where(ok, Tv, -np.inf)
            j = int(np.argmax(Tv))
            if Tv[j] > T0:
                I0, g0, T0 = float(Iv[j]), float(gv[j]), float(Tv[j])
            half *= 4.0 / (n - 1)
        return I0, g0

    # ---- the radius that makes a given torque on a ray
    def _radius(self, gam, target, I0, i_lim, c: Ctx, iters: int = 9):
        """The radius I at which T(I, gam) = target on each ray, by the secant method from the guess I0 (arrays of one shape). Returns (I, converged)."""
        def f(I):
            return self.state(I, gam, c)[0] - target
        Ia = np.clip(I0, 0.0, i_lim)
        Ib = np.clip(Ia * (1.0 + 1e-3) + 1e-6, 0.0, i_lim)
        fa, fb = f(Ia), f(Ib)
        for _ in range(iters):
            den = fb - fa
            step = np.where(np.abs(den) > 1e-300, fb * (Ib - Ia) / np.where(np.abs(den) > 1e-300, den, 1.0), 0.0)
            Ic = np.clip(Ib - step, 0.0, i_lim)
            Ia, fa, Ib, fb = Ib, fb, Ic, f(Ic)
        return Ib, np.abs(fb) <= 1e-9 * (1.0 + np.abs(target))

    def objective(self, name: str, I, gam, c: Ctx):
        """The quantity a level-set solve minimises at (I, gamma): "loss" (copper + iron + mechanical [W]) or "current" (|i| [A])."""
        if name == "current":
            return np.asarray(I, float)
        g = np.radians(gam)
        return self.m.total_loss(-I * np.sin(g), I * np.cos(g), c.rpm, c.v_dc, c.t_c)

    def _angle(self, I, target, g0, c: Ctx, iters: int = 9):
        """The angle gamma at which T(I, gamma) = target on each circle |i| = I, by the secant method from the guess g0 (arrays of one shape). Returns (gamma, converged)."""
        def f(g):
            return self.state(I, g, c)[0] - target
        ga = np.clip(g0, 0.0, 90.0)
        gb = np.clip(ga + 0.01, 0.0, 90.0)
        fa, fb = f(ga), f(gb)
        for _ in range(iters):
            den = fb - fa
            good = np.abs(den) > 1e-300
            step = np.where(good, fb * (gb - ga) / np.where(good, den, 1.0), 0.0)
            gc = np.clip(gb - step, 0.0, 90.0)
            ga, fa, gb, fb = gb, fb, gc, f(gc)
        return gb, np.abs(fb) <= 1e-9 * (1.0 + np.abs(target))

    def _boundary_crossings(self, c: Ctx, i_lim: float, Iw, vw, targets):
        """Where the curve of each target torque meets the voltage boundary |v| = V_max. On every ray the feasible part is an interval (the inside of the voltage curve), so the boundary
        has two points on it: the ENTRY (the first one from the origin) and the EXIT (the last one before the current limit); the level curve of the best solution often leaves the feasible
        region through the exit, deep in field weakening. Along each boundary curve the torque is smooth in the angle, so the sign changes of (torque - target) between rays are found
        however thin the feasible part of the level curve is; each crossing is then solved exactly (the secant method in the angle, the radius on the boundary solved on every trial ray).
        Returns (I, gamma, ok), each of shape (4, number of targets): the first and the last crossing of every target on each of the two curves."""
        nG, nR, nT = self.nG, vw.shape[1], len(targets)
        f = vw - c.v2
        rows = np.arange(nG)
        cols = np.arange(nT)
        outI, outG, outOk = [], [], []
        for kind in ("entry", "exit"):
            if kind == "entry":
                tr = (f[:, :-1] > 0.0) & (f[:, 1:] <= 0.0)
                k = tr.argmax(axis=1)
            else:
                tr = (f[:, :-1] <= 0.0) & (f[:, 1:] > 0.0)
                k = (nR - 2) - tr[:, ::-1].argmax(axis=1)
            has = tr.any(axis=1)
            fa, fb = f[rows, k], f[rows, k + 1]
            with np.errstate(divide="ignore", invalid="ignore"):
                guess = Iw[k] + np.clip(np.nan_to_num(fa / (fa - fb)), 0.0, 1.0) * (Iw[k + 1] - Iw[k])
            I_b, conv = self._radius_v(self.gam, guess, i_lim, c)
            T_b, _ = self.state(I_b, self.gam, c)
            good = has & conv
            S = T_b[:, None] - targets[None, :]
            pair = good[:-1, None] & good[1:, None] & (S[:-1] * S[1:] < 0.0)
            anyb = pair.any(axis=0)
            for idx in (pair.argmax(axis=0), (nG - 2) - pair[::-1].argmax(axis=0)):
                s0, s1 = S[idx, cols], S[idx + 1, cols]
                with np.errstate(divide="ignore", invalid="ignore"):
                    g0 = self.gam[idx] + np.clip(np.nan_to_num(s0 / (s0 - s1)), 0.0, 1.0) * self.dgam
                I0 = np.where(anyb, np.interp(g0, self.gam, I_b), 0.0)

                def h(g, I_guess):
                    Ix, cv = self._radius_v(g, I_guess, i_lim, c)
                    return self.state(Ix, g, c)[0] - targets, Ix, cv
                ga, gb = g0, np.clip(g0 + 0.01, 0.0, 90.0)
                fa_, Ia, _ = h(ga, I0)
                fb_, Ib_, cb = h(gb, Ia)
                for _ in range(8):
                    den = fb_ - fa_
                    good_d = np.abs(den) > 1e-300
                    step = np.where(good_d, fb_ * (gb - ga) / np.where(good_d, den, 1.0), 0.0)
                    ga, fa_ = gb, fb_
                    gb = np.clip(gb - step, 0.0, 90.0)
                    fb_, Ib_, cb = h(gb, Ib_)
                ok = anyb & cb & (np.abs(fb_) <= 1e-9 * (1.0 + np.abs(targets))) & (Ib_ <= i_lim * (1.0 + 1e-12))
                outI.append(np.where(ok, Ib_, 0.0))
                outG.append(np.where(ok, gb, 0.0))
                outOk.append(ok)
        return np.stack(outI), np.stack(outG), np.stack(outOk)

    def level_set(self, c: Ctx, i_lim: float, targets, objectives=("loss", "current"), grid=None, levels: int = LEVELSET_LEVELS, n: int = 21):
        """For each target torque (shaft or air gap, as the context's basis says) the feasible point of that torque that minimises each objective: minimum copper + iron +
        mechanical loss ("loss") or minimum current ("current": the MTPA control, which below the corner is the maximum-torque-per-ampere point and above it the
        smallest current that holds the voltage).
        Method. The curve of constant torque in the current plane is crossed by the grid's rays (a crossing on every ray: the torque rises with the radius) and by the grid's
        circles (at most two crossings on each circle: the torque is a hump in the angle). Both kinds are needed: near the MTPA point the curve meets the rays squarely and is
        tangent to a circle; deep in field weakening it runs almost along a ray, so a ray grid misses the thin feasible lens at the current limit while the circles cross it
        again and again. Every crossing is evaluated exactly, the best feasible one is picked, and two zooms refine it (in the angle with the radius re-solved on every ray,
        then in the radius with the angle re-solved on every circle), each keeping the best feasible point it has seen.
        Returns {objective: {"I", "gam", "found"}} with arrays of the length of `targets`."""
        T, vv = self.grid(c) if grid is None else grid
        i_lim = min(i_lim, self.i_max)
        kmax = min(int(math.floor(i_lim / self.dI + 1e-9)), self.nR - 1)
        Tw, vw, Iw = T[:, :kmax + 1], vv[:, :kmax + 1], self.I[:kmax + 1]
        if i_lim > Iw[-1] + 1e-12:
            Tb, vb = self.state(np.full(self.nG, i_lim), self.gam, c)
            Tw, vw, Iw = np.concatenate([Tw, Tb[:, None]], axis=1), np.concatenate([vw, vb[:, None]], axis=1), np.append(Iw, i_lim)
        nR, nT = Tw.shape[1], len(targets)
        targets = np.asarray(targets, float)
        # ---- the crossings with the rays
        Tm = np.maximum.accumulate(Tw, axis=1)
        cnt = (Tm[:, :, None] < targets[None, None, :]).sum(axis=1)                    # the first radius index whose running maximum reaches the target
        reach = (cnt >= 1) & (cnt < nR)
        k = np.clip(cnt, 1, nR - 1)
        rows = np.arange(self.nG)[:, None]
        T_lo, T_hi = Tw[rows, k - 1], Tw[rows, k]
        with np.errstate(divide="ignore", invalid="ignore"):
            frac = np.clip(np.nan_to_num((targets[None, :] - T_lo) / (T_hi - T_lo)), 0.0, 1.0)
        I_r = Iw[k - 1] + frac * (Iw[k] - Iw[k - 1])
        G_r = np.broadcast_to(self.gam[:, None], I_r.shape)
        # ---- the crossings with the circles: the first and the last sign change of (torque - target) along the angle, on every radius
        g_first, g_last, any_c = [], [], []
        for t0 in range(0, nT, 24):
            S = Tw[:, :, None] - targets[None, None, t0:t0 + 24]                          # (rays, radii, targets)
            chg = S[:-1] * S[1:] < 0.0
            any_c.append(chg.any(axis=0))
            for ii, store in ((chg.argmax(axis=0), g_first), ((self.nG - 2) - chg[::-1].argmax(axis=0), g_last)):
                s0 = np.take_along_axis(S, ii[None], axis=0)[0]
                s1 = np.take_along_axis(S, (ii + 1)[None], axis=0)[0]
                with np.errstate(divide="ignore", invalid="ignore"):
                    fr = np.clip(np.nan_to_num(s0 / (s0 - s1)), 0.0, 1.0)
                store.append(self.gam[ii] + fr * self.dgam)
        g_first, g_last, any_c = np.concatenate(g_first, axis=1), np.concatenate(g_last, axis=1), np.concatenate(any_c, axis=1)
        I_k = np.broadcast_to(Iw[:, None], g_first.shape)
        cand_I = np.concatenate([I_r, I_k, I_k], axis=0)
        cand_G = np.concatenate([G_r, g_first, g_last], axis=0)
        cand_ok = np.concatenate([reach, any_c, any_c], axis=0)
        if np.isfinite(c.v2):                                                           # the crossings with the voltage boundary itself (module docstring: thin lenses)
            bI, bG, bok = self._boundary_crossings(c, i_lim, Iw, vw, targets)
            cand_I, cand_G, cand_ok = np.concatenate([cand_I, bI], axis=0), np.concatenate([cand_G, bG], axis=0), np.concatenate([cand_ok, bok], axis=0)
        T_cd, v_cd = self.state(cand_I, cand_G, c)
        feas = cand_ok & (v_cd <= c.v2 * (1.0 + V_SLACK)) & (np.abs(T_cd - targets[None, :]) <= 1e-3 * (1.0 + np.abs(targets[None, :])))
        out = {}
        cols = np.arange(nT)
        lin = np.linspace(-1.0, 1.0, n)
        tg = np.broadcast_to(targets[:, None], (nT, n))
        for name in objectives:
            val = np.where(feas, self.objective(name, cand_I, cand_G, c), np.inf)
            # the crossings are interpolated, so their torque is a little off the target: project the best few (TOP_K) of each target onto the exact level set, along their ray and
            # along their circle, keep whichever stays feasible and take the best of those, so that everything the zooms compare it with is exact, a torque error cannot buy a
            # lower loss, and one candidate that refuses to project does not lose the target
            order = np.argsort(val, axis=0)[:TOP_K]
            Ig, Gg, Vg = (np.take_along_axis(a, order, axis=0) for a in (cand_I, cand_G, val))
            tk = np.broadcast_to(targets[None, :], Ig.shape)
            fin = np.isfinite(Vg)
            I1, c1 = self._radius(Gg, tk, Ig, i_lim, c)
            ok1 = fin & c1 & (self.state(I1, Gg, c)[1] <= c.v2 * (1.0 + V_SLACK))
            g2, c2 = self._angle(Ig, tk, Gg, c)
            ok2 = fin & c2 & (self.state(Ig, g2, c)[1] <= c.v2 * (1.0 + V_SLACK))
            allO = np.concatenate([np.where(ok1, self.objective(name, I1, Gg, c), np.inf), np.where(ok2, self.objective(name, Ig, g2, c), np.inf)], axis=0)
            allI, allG = np.concatenate([I1, Ig], axis=0), np.concatenate([Gg, g2], axis=0)
            j = np.argmin(allO, axis=0)
            best = allO[j, cols]
            found = np.isfinite(best)
            I_b, g_b = np.where(found, allI[j, cols], 0.0), np.where(found, allG[j, cols], 0.0)
            for kind in ("ray", "circle"):
                half = np.full(nT, self.dgam if kind == "ray" else self.dI)
                for _ in range(levels):
                    if kind == "ray":
                        gv = np.clip(g_b[:, None] + half[:, None] * lin[None, :], 0.0, 90.0)
                        Iv, conv = self._radius(gv, tg, np.broadcast_to(I_b[:, None], gv.shape), i_lim, c)
                    else:
                        Iv = np.clip(I_b[:, None] + half[:, None] * lin[None, :], 0.0, i_lim)
                        gv, conv = self._angle(Iv, tg, np.broadcast_to(g_b[:, None], Iv.shape), c)
                    _, vvv = self.state(Iv, gv, c)
                    ok = conv & (vvv <= c.v2 * (1.0 + V_SLACK)) & found[:, None]
                    cand = np.where(ok, self.objective(name, Iv, gv, c), np.inf)
                    j = np.argmin(cand, axis=1)
                    better = cand[cols, j] < best
                    g_b = np.where(better, gv[cols, j], g_b)
                    I_b = np.where(better, Iv[cols, j], I_b)
                    best = np.where(better, cand[cols, j], best)
                    half = half * 4.0 / (n - 1)
            out[name] = {"I": I_b, "gam": g_b, "found": found}
        return out


# ====================================================================================================== the envelope
def _per_speed(x, rpms):
    """A number, an array of the length of rpms, or a function of the speed, as an array over the speeds."""
    if callable(x):
        return np.array([float(x(r)) for r in rpms])
    a = np.asarray(x, float)
    return np.full(len(rpms), float(a)) if a.ndim == 0 else a


@dataclass
class Envelope:
    """The result of envelope(): `d` the arrays over the speed (see keys), `meta` the scalars. env["T_shaft"] etc."""
    d: dict
    meta: dict

    def __getitem__(self, key):
        return self.d[key]


def speed_grid(top: float = TOP_RPM, step: float = 100.0) -> np.ndarray:
    """0, step, 2 step ... and the top speed itself."""
    r = np.arange(0.0, top, step)
    return np.append(r, top)


def envelope(m: Machine, bus: Bus, rpms, i_limit, *, t_c=100.0, basis: str = "shaft", scan: Scan | None = None, locate: bool = True, **scan_kw) -> Envelope:
    """The torque-speed envelope: at each speed the largest torque with |i| <= i_limit(rpm) and |v| <= V_max, on the maps (MTPA, then field weakening, then MTPV).
    i_limit: a number, an array over the speeds or a function of the speed [A peak]; t_c likewise [degC of the copper]. Returns an Envelope whose arrays are everything
    Machine.point gives at the solution plus the solver's flags; speeds without a positive shaft torque carry torque 0 and region "off". With locate=True the corner speed (the
    last speed of the full current with the voltage not binding) and the speed where the torque falls to zero are found by bisection and put in meta."""
    rpms = np.asarray(rpms, float)
    lims, tcs = _per_speed(i_limit, rpms), _per_speed(t_c, rpms)
    sc = scan or Scan(m, i_max=float(lims.max()), **scan_kw)
    n = len(rpms)
    I, gam, T, v_act, i_act, feas = (np.zeros(n) for _ in range(6))
    for k, (rpm, lim, tc) in enumerate(zip(rpms, lims, tcs)):
        sol = sc.max_torque(sc.ctx(rpm, bus, tc, basis), lim)
        if sol is not None:
            I[k], gam[k], T[k], v_act[k], i_act[k], feas[k] = sol["I"], sol["gam"], sol["T"], sol["v_active"], sol["i_active"], 1.0
    g = np.radians(gam)
    pt = m.point(-I * np.sin(g), I * np.cos(g), rpms, bus, tcs)
    positive = (feas > 0) & (pt["T_shaft"] > 0 if basis == "shaft" else pt["T_em"] > 0)
    region = np.where(~positive, "off", np.where(v_act > 0, np.where(i_act > 0, "FW", "MTPV"), "MTPA"))
    d = {k: (np.where(positive, v, np.nan) if k not in ("rpm", "t_cu", "R") else v) for k, v in pt.items()}     # a speed without torque keeps its speed and no operating point
    d.update({"i_limit": lims, "v_active": v_act > 0, "i_active": i_act > 0, "feasible": feas > 0, "region": region, "T": np.where(positive, pt["T_shaft" if basis == "shaft" else "T_em"], 0.0)})
    meta = {"basis": basis, "v_max": bus.v_max, "v_dc": bus.v_dc, "k": bus.k, "stack_mm": m.stack_mm, "winding": m.w.label(), "scan": (sc.nG, sc.nR)}
    env = Envelope(d, meta)
    if locate:
        env.meta.update(_locate(m, bus, sc, rpms, lims, tcs, basis, positive, v_act > 0, I, gam, T))
    return env


def _locate(m, bus, sc, rpms, lims, tcs, basis, positive, v_act, I, gam, T) -> dict:
    """The corner speed (the highest speed whose best point is not voltage-limited) and the speed where the torque falls to zero, refined by bisection between the grid's speeds."""
    out = {"corner_rpm": None, "corner_T": None, "zero_rpm": None}
    lim_of = lambda r: float(np.interp(r, rpms, lims))
    tc_of = lambda r: float(np.interp(r, rpms, tcs))

    def solve(r):
        return sc.max_torque(sc.ctx(r, bus, tc_of(r), basis), lim_of(r))

    def solve_free(r):                                  # the maximum with no voltage limit: its voltage is what decides whether the limit binds
        return sc.max_torque(replace(sc.ctx(r, bus, tc_of(r), basis), v2=np.inf), lim_of(r))
    free = np.nonzero(positive & ~v_act)[0]
    if len(free) and free.max() + 1 < len(rpms) and v_act[free.max() + 1]:
        lo, hi = rpms[free.max()], rpms[free.max() + 1]
        for _ in range(45):
            mid = 0.5 * (lo + hi)
            if solve_free(mid)["v2"] <= bus.v_max ** 2:
                lo = mid
            else:
                hi = mid
        out["corner_rpm"], out["corner_T"] = lo, solve_free(lo)["T"]
    elif len(free) and free.max() + 1 == len(rpms):
        out["corner_rpm"] = None                       # the voltage never binds: no corner below the top speed
        out["corner_T"] = float(T[free.max()])
    pos = np.nonzero(positive)[0]
    if len(pos) and pos.max() + 1 < len(rpms):
        lo, hi = rpms[pos.max()], rpms[pos.max() + 1]
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            s = solve(mid)
            t_ok = s is not None and s["T"] > 0
            if t_ok:
                lo = mid
            else:
                hi = mid
        out["zero_rpm"] = lo
    elif len(pos):
        out["zero_rpm"] = float(rpms[pos.max()])
    return out


# ====================================================================================================== the efficiency map
def efficiency_map(m: Machine, bus: Bus, rpms, fractions, i_limit, *, t_c=100.0, scan: Scan | None = None, env: Envelope | None = None, **scan_kw) -> dict:
    """The efficiency over the (speed, shaft torque) plane inside the envelope of this machine, bus and limits, at the torques fraction x T_env(speed) for each fraction (T_env
    from the envelope at these very speeds, so that a fraction of 1 is the envelope's own point), for two controls: "mtpa" (the smallest current that makes the torque within the
    voltage limit: MTPA below the corner, field weakening above) and "minloss" (the current vector that minimises the total loss). i_limit and t_c: a number, an array over the
    speeds or a function of the speed. Returns {"rpm", "T_env", "frac", "T", "i_limit"} (the last three n_speed x n_frac) and, for each control, the solution (i_d, i_q, found)
    and Machine.point's arrays at it (n_speed x n_frac), and "saving" = loss(mtpa) - loss(minloss) [W]."""
    rpms, fractions = np.asarray(rpms, float), np.asarray(fractions, float)
    lim, tcs = _per_speed(i_limit, rpms), _per_speed(t_c, rpms)
    sc = scan or Scan(m, i_max=float(lim.max()), **scan_kw)
    if env is None:
        env = envelope(m, bus, rpms, lim, t_c=tcs, scan=sc, locate=False)
    t_env = np.asarray(env["T"], float)
    nS, nF = len(rpms), len(fractions)
    sol = {k: {"i_d": np.full((nS, nF), np.nan), "i_q": np.full((nS, nF), np.nan), "found": np.zeros((nS, nF), bool)} for k in ("mtpa", "minloss")}
    T = np.outer(t_env, fractions)
    for j, (rpm, tl, tc) in enumerate(zip(rpms, lim, tcs)):
        if t_env[j] <= 0:
            continue
        c = sc.ctx(rpm, bus, tc, "shaft")
        res = sc.level_set(c, tl, T[j], objectives=("loss", "current"))
        for key, name in (("minloss", "loss"), ("mtpa", "current")):
            r = res[name]
            g = np.radians(r["gam"])
            i_d, i_q, f = -r["I"] * np.sin(g), r["I"] * np.cos(g), r["found"].copy()
            top = fractions >= 1.0 - 1e-12                       # the envelope's own point is the only point of the top torque
            if top.any():
                i_d, i_q, f = np.where(top, env["i_d"][j], i_d), np.where(top, env["i_q"][j], i_q), f | top
            sol[key]["i_d"][j], sol[key]["i_q"][j], sol[key]["found"][j] = i_d, i_q, f
    out = {"rpm": rpms, "T_env": t_env, "frac": fractions, "T": T, "i_limit": lim}
    rr = np.broadcast_to(rpms[:, None], (nS, nF))
    tt = np.broadcast_to(tcs[:, None], (nS, nF))
    for key in ("mtpa", "minloss"):
        s_ = sol[key]
        f = s_["found"]
        pt = m.point(np.where(f, s_["i_d"], 0.0), np.where(f, s_["i_q"], 0.0), rr, bus, tt)
        out[key] = {"found": f, "i_d": s_["i_d"], "i_q": s_["i_q"], **{k: np.where(f, v, np.nan) for k, v in pt.items() if k not in ("rpm", "i_d", "i_q", "t_cu")}}
    out["saving"] = out["mtpa"]["P_loss"] - out["minloss"]["P_loss"]
    return out


# ====================================================================================================== the command line: the interface at work
def main() -> int:
    import time
    m = Machine(PROVISIONAL_WINDING, 60.0)
    bus = Bus(48.0)
    print(f"{m.w.label()}, {m.stack_mm:g} mm, fill {m.fill:.3f}, R(100 C) {m.resistance(100.0):.4f} Ohm, bus {bus.v_dc:g} V, V_max {bus.v_max:.3f} V, maps end at {m.i_cap:.1f} A")
    i_lim = m.w.current_from_density(6.0)
    t0 = time.perf_counter()
    env = envelope(m, bus, speed_grid(step=500.0), i_lim)
    dt = time.perf_counter() - t0
    print(f"envelope at {i_lim:.2f} A (6 A/mm2 rms): corner {env.meta['corner_rpm']} rpm, torque {env.meta['corner_T']} N.m, zero torque at {env.meta['zero_rpm']} rpm ({dt:.2f} s)")
    for r, t, reg, e in zip(env["rpm"], env["T"], env["region"], env["eta"]):
        print(f"   {r:8.0f} rpm  T {t:7.3f} N.m  {reg:5s}  eta {100 * e:6.2f} %")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
