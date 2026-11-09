"""optimization.thermal_machine: the thermal network driven by the motor's own models: the flux maps (torque at the MTPA point), the iron-loss map (loss by region
at that point and speed) and the copper-loss model (a function of the copper temperature), so that the continuous limit of a cooling is found with the real losses.

What it is, and why
    thermal.py takes the losses as inputs. The motor has three models of its own losses and torque: maps.py (the flux linkage and torque at any current),
    iron_loss.py (the iron loss of teeth, yoke and rotor at any current and speed) and copper_loss.py (the DC and AC copper loss at any current, speed, copper
    temperature and stack). This file puts them behind one object, Machine, in the unit the network uses (ampere-turns per coil side, peak), and solves the question
    that is asked of a cooling: how much current can the motor carry for ever, at this speed, before a temperature limit is reached, and what torque, mechanical
    power, losses and efficiency is that?

The units (the scaling rule of maps.py: the field depends on the ampere-turns only)
    The maps and the iron-loss map are those of the 10-turn, one-path placeholder winding of the FEM model. A winding of N turns a coil side in a paths in parallel carries
    the phase current i = AT a / N for AT ampere-turns a coil side (peak); in the placeholder's terms the current is I0 = AT / 10. So
        AT [A-turns, peak, per coil side]  ->  I0 = AT / 10  ->  maps.mtpa(I0) = (gamma, i_d, i_q, torque)  ->  iron_loss(i_d, i_q, rpm, stack)
    gamma is the current angle from the q axis toward -d, the convention of the copper-loss model too. The maps cover |I0| up to 3 x the rated 11.8627 A, that is up to
    AT = 355.9 (3 x 118.6); beyond it nothing is extrapolated: a continuous limit that would lie beyond is reported as limited by the maps.
    The network's current density is J_rms = AT / (sqrt 2 x the copper area of one coil side); the copper area comes from the winding's own fill.

The copper loss as the network's copper function
    copper_loss.losses(winding, op, t_copper_c, stack_mm)['p_total'] (the DC loss in the slots and the end turns at the copper temperature, plus skin, proximity, PWM-ripple
    and circulating-current loss) is given to the network as a function of the copper temperature; the slot's share is the module's own (the DC loss and its skin
    share of the end turns are the end turns'; the rest is in the slots). The network is given the same diamond end turn as the copper-loss model,
    Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6), and the winding's fill, so that the DC losses of the two agree to rounding. The example winding
    {N: 26, a: 4, n_s: 8, d: 0.315 mm, arrangement: cluster} is a stand-in for the winding the optimizer will choose; its AC loss is small (proximity 0.4 % of the DC
    loss at 16,330 rpm, PWM ripple 0.02 %), so the ampere-turn tables hold for another winding of the same fill.

The continuous limit
    continuous_limit(model, machine) bisects on the ampere-turns (through the current density): at each trial the MTPA point, the iron loss there and the copper
    function are put into thermal.ThermalModel.sustainable(losses_at=...), which finds the largest current at which neither the winding's hot spot (class F, 155 C) nor
    the magnets (80 C, or 150 C for the other grade) is passed. Every loss that depends on the current follows it: the iron loss grows with the MTPA current (the
    flux density in the steel), the copper loss with the square. At standstill the iron loss, windage and bearing loss are zero.
    LABEL: the operating point is MTPA at every speed. Above the corner speed the real motor weakens the field (less flux, less iron loss, a different current
    angle) and the motor map is the place for that; the limits here are the thermal limits AT MTPA CURRENTS, exact where the speed is below the corner and an
    upper bound on the iron loss (so a lower bound on the current) above it.

Use
    from optimization import thermal, thermal_machine as tm
    mach = tm.Machine()
    model = tm.model_for(thermal.Cooling.shaft_fan(6.0, 3000.0), stack_mm=60.0, speed_rpm=3000.0, machine=mach)
    cont = tm.continuous_limit(model, mach)                      # dict: AT_pk, J_rms, torque, losses, temperatures, binding, ...
    short = tm.short_term(model, mach, cont, factor=2.0)         # torque and time to a limit at 2 x the continuous current, from cold
    res = tm.operating_point(model, mach, 120.0, gamma_deg=40.0) # the steady state at any other point (a field-weakened one): res.margin_magnet, res.margin_winding

Run
    uv run python -m optimization.thermal_machine         one worked example and its speed (a continuous limit is about 0.1 s)

Numbers to expect (this PC, 2026-11-09; the 60 mm stack at 3,000 rpm, N45 80 C, class F, 40 C ambient; the example winding, fill 0.4476)
    natural convection      48.6 A-turns peak (2.12 A/mm2 rms), MTPA torque 2.15 N.m, T x omega 674 W (an estimate), copper 12.2 W, iron 20.9 W, efficiency 95.2 %, hot spot 81.6 C, magnets 80.0 C;
                            at 2x the current: 4.40 N.m, the magnets reach 80 C after 3,563 s from cold
    fixed 6 m/s forced air  157.8 A-turns (6.88), 7.30 N.m, 2,292 W, copper 131.8 W, iron 30.8 W, 92.9 %, hot spot 91.8 C; at 2x: 13.99 N.m, the winding reaches 155 C after 466 s
    shaft fan, 6 m/s        158.1 A-turns (6.90), 7.31 N.m: 0.2 % from the fixed fan (the same 6 m/s at this speed)
    liquid jacket           350.5 A-turns (15.29), 15.10 N.m, 4,744 W, copper 714.8 W, 86.3 %, hot spot 134.3 C; 2x is 701 A-turns, beyond the maps (up to 355.9): the winding reaches 155 C after 25 s
    A continuous limit takes 0.1 to 0.2 s (about 40 steady states, each with about 15 calls of the copper-loss model at 189 us); Machine() loads in a few ms.

Self-check (the answers follow).
  1. The maps' current is 11.8627 A at 10 turns in one path. What phase current does the example winding (26 turns in 4 paths) carry at the same ampere-turns?
  2. Why does the iron loss belong inside the bisection and not in front of it?
  3. At 10,000 rpm the natural-convection limit is reported with J = 0 and 'feasible' False. What does that say about the motor, and what does it not say?
Answers.
  1. The ampere-turns a coil side are N i / a, so i = AT a / N = 118.6 x 4 / 26 = 18.25 A (peak); the placeholder's 11.86 A is the same field (maps.winding_scale = 26 / 40 = 0.65).
  2. The MTPA flux density in the steel grows with the current, so the iron loss at 100 A-turns is not the iron loss at 50; a loss fixed before the search would find a limit for the wrong loss.
  3. That the iron, windage and bearing heat alone (116 W of iron at 60 mm and 10,000 rpm) already passes the magnets' 80 C in still-air convection: no current at all can be carried
     continuously with that cooling at that speed on MTPA-level flux. It does not say the motor cannot run there: a motor weakens its field above the corner speed and has less iron loss.
"""
from __future__ import annotations

import math
import time

import numpy as np

from optimization import copper_loss, maps
from optimization import thermal
from optimization.iron_loss import iron_loss
from optimization.thermal import (CLASS_F_C, MAGNET_LIMIT_C, RATED_RPM, Assumptions, Cooling, Losses, ThermalModel, WindingSpec)

SQRT2 = math.sqrt(2.0)
EXAMPLE_WINDING = {"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367, "arrangement": "cluster"}
END_TURN_MM, END_OVERHANG_MM = 48.55, 13.6           # the copper-loss model's diamond end turn (one end turn, and how far it sticks out of the stack)
F_SW_HZ, V_DC = 20000.0, 48.0                        # the copper-loss model's example operating point (the PWM ripple's inputs)
PLACEHOLDER_TURNS = 10                               # the maps' winding: 10 turns a coil side in one path


def assumptions(**kw) -> Assumptions:
    """The network's assumptions with the copper-loss model's end turn (and any other assumed coefficient overridden by keyword)."""
    return Assumptions(end_turn_mm=END_TURN_MM, end_overhang_mm=END_OVERHANG_MM, **kw)


class Machine:
    """The flux maps, the iron-loss map and the copper-loss model of the optimized motor, in ampere-turns per coil side (peak)."""

    def __init__(self, winding: dict | None = None, f_sw_hz: float = F_SW_HZ, v_dc: float = V_DC):
        self.fm = maps.FluxMaps.load()
        self.winding = dict(EXAMPLE_WINDING if winding is None else winding)
        self.f_sw, self.v_dc = f_sw_hz, v_dc
        self.at_rated = PLACEHOLDER_TURNS * self.fm.rated                # 118.6 A-turns: the maps' rated current
        self.at_max = 3.0 * self.at_rated                                # the largest current the maps cover (|I0| <= 3 I_R)
        self.fill = copper_loss.losses(self.winding, {"i_peak": 1.0, "gamma_deg": 0.0, "rpm": 0.0, "f_sw": 0.0, "v_dc": v_dc}, 20.0, 20.0)["fill"]

    # ------------------------------------------------------------------------------------------------ torque and iron loss at the MTPA point
    def point(self, at_pk: float, rpm: float, stack_mm: float, clip: bool = False, gamma_deg: float | None = None) -> dict:
        """The MTPA point (or, with gamma_deg, the point at that current angle from the q axis toward -d, 0 ... 90 degrees: a field-weakened point) at at_pk ampere-turns a
        coil side (peak) for the stack: {"gamma_deg", "i_d", "i_q" (placeholder terms, A peak), "torque" [N.m] (None beyond the maps; the d-q torque of the maps at an angle
        given), "iron_teeth", "iron_yoke", "iron_rotor" [W], "beyond_maps"}. Beyond the maps' range a ValueError, or with clip=True the point at the range's edge (torque None)."""
        i0 = at_pk / PLACEHOLDER_TURNS
        beyond = i0 > 3.0 * self.fm.rated * (1.0 + 1e-9)
        if beyond and not clip:
            raise ValueError(f"{at_pk:.1f} A-turns is beyond the flux maps (up to {self.at_max:.1f})")
        i0 = min(i0, 3.0 * self.fm.rated)
        if gamma_deg is None:
            gamma, i_d, i_q, torque = self.fm.mtpa(i0, stack_mm=stack_mm)
        else:
            gamma = float(gamma_deg)
            i_d, i_q = -i0 * math.sin(math.radians(gamma)), i0 * math.cos(math.radians(gamma))
            torque = float(self.fm.torque_dq(i_d, i_q, stack_mm=stack_mm))
        il = iron_loss(i_d, i_q, float(rpm), stack_mm=stack_mm)
        return {"gamma_deg": float(gamma), "i_d": float(i_d), "i_q": float(i_q), "torque": None if beyond else float(torque), "iron_teeth": float(il["teeth"]),
                "iron_yoke": float(il["yoke"]), "iron_rotor": float(il["rotor"]), "beyond_maps": bool(beyond)}

    # ------------------------------------------------------------------------------------------------ the copper loss as a function of the copper temperature
    def copper(self, at_pk: float, gamma_deg: float, rpm: float, stack_mm: float):
        """(f, share, r100): f(t_copper_c) -> the total copper loss [W] of the whole winding at that copper temperature (copper_loss.losses' p_total), share the part of
        it in the slots (the rest is the end turns' DC loss with its skin share), r100 the module's full answer at 100 C."""
        w = self.winding
        op = {"i_peak": at_pk * w["a"] / w["N"], "gamma_deg": gamma_deg, "rpm": float(rpm), "f_sw": self.f_sw, "v_dc": self.v_dc}

        def f(t_c: float) -> float:
            return copper_loss.losses(w, op, t_c, stack_mm)["p_total"]
        r100 = copper_loss.losses(w, op, 100.0, stack_mm)
        dc = r100["p_dc_slot"] + r100["p_dc_end"]
        if r100["p_total"] <= 0.0 or dc <= 0.0:
            return f, None, r100
        end = r100["p_dc_end"] * (1.0 + r100["p_skin"] / dc)
        return f, 1.0 - end / r100["p_total"], r100


def model_for(cooling: Cooling, stack_mm: float, speed_rpm: float, machine: Machine | None = None, assume: Assumptions | None = None) -> ThermalModel:
    """A ThermalModel of the optimized design with the winding's fill and the copper-loss model's end turn."""
    machine = machine or Machine()
    return ThermalModel(None, stack_mm, cooling, speed_rpm, WindingSpec(machine.fill), assume or assumptions())


def _a_side_mm2(model: ThermalModel) -> float:
    return model.ampere_turns(1.0)["a_cu_side_mm2"]


def losses_at_factory(model: ThermalModel, machine: Machine, magnet_loss_w: float = 0.0):
    """losses_at(j_rms) -> thermal.Losses at the current density j_rms [A/mm^2] for ThermalModel.sustainable: the MTPA point of the ampere-turns, its iron loss and
    the copper-loss model as a function of the copper temperature."""
    a_side = _a_side_mm2(model)

    def losses_at(j: float) -> Losses:
        at = SQRT2 * j * a_side
        pt = machine.point(min(at, machine.at_max), model.speed_rpm, model.stack_mm)
        f, share, _ = machine.copper(min(at, machine.at_max), pt["gamma_deg"], model.speed_rpm, model.stack_mm)
        return Losses(copper=f, copper_slot_fraction=share, iron_teeth=pt["iron_teeth"], iron_yoke=pt["iron_yoke"], iron_rotor=pt["iron_rotor"], magnet=magnet_loss_w)
    return losses_at


def operating_point(model: ThermalModel, machine: Machine, at_pk: float, gamma_deg: float | None = None, *, wind_limit_c: float = CLASS_F_C,
                    magnet_limit_c: float = MAGNET_LIMIT_C, magnet_loss_w: float = 0.0, tol: float = 1e-4):
    """The network's steady state at ONE operating point of the motor: at_pk ampere-turns a coil side (peak) at the current angle gamma_deg (None: MTPA) at the model's speed
    and stack, with the iron loss of that point, the copper-loss model's copper loss at the copper temperature found, and the model's windage, bearings and cooling. Returns
    the thermal.SteadyResult (temperatures, hot spot, margins to the limits: positive is inside, binding), so that the motor map can ask of any point of its torque-speed plane,
    field weakening included, whether it can be held; model.heat_flows(result) says where the heat goes."""
    pt = machine.point(at_pk, model.speed_rpm, model.stack_mm, gamma_deg=gamma_deg)
    f, share, _ = machine.copper(at_pk, pt["gamma_deg"], model.speed_rpm, model.stack_mm)
    return model.steady_state(Losses(copper=f, copper_slot_fraction=share, iron_teeth=pt["iron_teeth"], iron_yoke=pt["iron_yoke"], iron_rotor=pt["iron_rotor"],
                                     magnet=magnet_loss_w), wind_limit_c=wind_limit_c, magnet_limit_c=magnet_limit_c, tol=tol)


def continuous_limit(model: ThermalModel, machine: Machine, *, wind_limit_c: float = CLASS_F_C, magnet_limit_c: float = MAGNET_LIMIT_C,
                     magnet_loss_w: float = 0.0) -> dict:
    """The continuous (thermal) limit of the model's cooling, stack and speed at MTPA currents, with the real iron, copper, windage and bearing losses. Returns
      AT_pk, AT_rms (per coil side), J_rms [A/mm^2], feasible (False: a limit is passed even with no current: the iron, windage and bearing heat alone),
      limited_by_maps (the thermal limit lies beyond the flux maps' range: AT is the range's edge), binding ('winding' or 'magnet'), J_winding_limit, J_magnet_limit,
      gamma_deg, i_d, i_q, torque [N.m], p_em_w (T x omega: the air-gap power, an ESTIMATE from the maps), the losses (copper_slot, copper_end, iron_teeth, iron_yoke,
      iron_rotor, windage, bearing, fan, magnet), p_loss_in_motor (the heat made in the motor), efficiency (an estimate: see below), the temperatures (hot_spot,
      hot_spot_slot, hot_spot_end, magnet, T per node), the heat flows (heat) and the SteadyResult (result).
    Efficiency = (P_em - windage - bearings - shaft fan) / (P_em + copper + iron + magnet loss + an external fan): the iron loss is taken on the electrical side (it
    is supplied through the stator, as in the usual definition), the mechanical losses come off the shaft; the forced_air preset's fan is external (its power is added
    to the input), a shaft fan's is on the shaft. P_em = T x omega at the MTPA torque (the flux maps' d-q torque), at the speed of the model."""
    a_side = _a_side_mm2(model)
    j_max = machine.at_max / (SQRT2 * a_side)
    out = model.sustainable(wind_limit_c, magnet_limit_c, j_hi=j_max, losses_at=losses_at_factory(model, machine, magnet_loss_w))
    return describe(model, machine, out, magnet_loss_w, j_max)


def describe(model: ThermalModel, machine: Machine, out: dict, magnet_loss_w: float = 0.0, j_max: float | None = None) -> dict:
    """Turn ThermalModel.sustainable's answer into the full description of the operating point (see continuous_limit)."""
    a_side = _a_side_mm2(model)
    j = out["J_rms"]
    at = SQRT2 * j * a_side
    at_c = min(at, machine.at_max)
    pt = machine.point(at_c, model.speed_rpm, model.stack_mm)
    # the bisection's last solve stops at 1e-4 K; the state reported (and its heat flows, which balance to rounding) is solved again to 1e-10 K
    res = model.steady_state(losses_at_factory(model, machine, magnet_loss_w)(j), wind_limit_c=out["wind_limit_c"], magnet_limit_c=out["magnet_limit_c"],
                             t_guess=np.array([out["T"][n] for n in thermal.NODES]), max_iter=200, tol=1e-10)
    src = res.sources
    omega = 2.0 * math.pi * model.speed_rpm / 60.0
    torque = pt["torque"]
    p_em = torque * omega
    kind = model.cooling.kind
    p_fan = model.fan_power()
    fan_on_shaft = p_fan if kind == "shaft_fan" else 0.0
    fan_external = p_fan if kind == "forced_air" else 0.0
    p_mech = src["windage"] + src["bearing"] + fan_on_shaft
    p_in = p_em + res.copper_loss + src["iron_teeth"] + src["iron_yoke"] + src["iron_rotor"] + src["magnet"] + fan_external
    eff = (p_em - p_mech) / p_in if p_in > 0 and p_em > 0 else float("nan")
    limited_by_maps = bool(j_max is not None and j >= j_max * (1.0 - 1e-6) and out["feasible"])
    heat = model.heat_flows(res)
    return {"cooling": model.cooling.kind, "label": model.cooling.label, "stack_mm": model.stack_mm, "rpm": model.speed_rpm, "magnet_limit_c": out["magnet_limit_c"],
            "wind_limit_c": out["wind_limit_c"], "J_rms": j, "AT_pk": at, "AT_rms": at / SQRT2, "feasible": out["feasible"], "limited_by_maps": limited_by_maps,
            "binding": "maps" if limited_by_maps else out["binding"], "J_winding_limit": out["J_winding_limit"], "J_magnet_limit": out["J_magnet_limit"], "gamma_deg": pt["gamma_deg"], "i_d": pt["i_d"],
            "i_q": pt["i_q"], "torque": torque, "p_em_w": p_em, "copper_slot": src["copper_slot"], "copper_end": src["copper_end"], "copper": res.copper_loss,
            "iron_teeth": src["iron_teeth"], "iron_yoke": src["iron_yoke"], "iron_rotor": src["iron_rotor"], "iron": src["iron_teeth"] + src["iron_yoke"] + src["iron_rotor"],
            "windage": src["windage"], "bearing": src["bearing"], "fan": p_fan, "fan_on_shaft": bool(kind == "shaft_fan"), "magnet_loss": src["magnet"],
            "p_loss_in_motor": float(sum(src.values())), "p_in_w": p_in, "efficiency": eff, "hot_spot": res.hot_spot, "hot_spot_slot": res.hot_spot_slot,
            "hot_spot_end": res.hot_spot_end, "magnet": res.magnet, "T": dict(res.T), "margin_winding": res.margin_winding, "margin_magnet": res.margin_magnet,
            "heat": heat, "result": res, "point": pt, "v_air": model.v_air if kind in ("forced_air", "shaft_fan") else 0.0}


def short_term(model: ThermalModel, machine: Machine, cont: dict, factor: float = 2.0, *, wind_limits=(CLASS_F_C,), magnet_limit_c: float | None = None,
               t_max_s: float = 4 * 3600.0) -> dict:
    """The short-term capability at `factor` x the continuous ampere-turns: the MTPA torque from the maps (saturation included; None beyond the maps' range, where the
    iron loss is taken at the range's edge, a lower bound) and the time from COLD (everything at the ambient) until a limit is reached with the copper-loss model's
    copper loss (as an effective current density through the DC formula: the module's AC parts are about 1 %), the iron loss of the new MTPA point and the model's
    windage and bearings. Returns {"factor", "AT_pk", "torque", "t_winding_s", "t_magnet_s", "t_first_s", "first" ('winding' | 'magnet' | None), "beyond_maps",
    "j_eff"}; a time None means no limit within t_max_s. None times if the continuous limit is infeasible (no current can be carried)."""
    magnet_limit_c = cont["magnet_limit_c"] if magnet_limit_c is None else magnet_limit_c
    if not cont["feasible"]:
        return {"factor": factor, "AT_pk": 0.0, "torque": None, "t_winding_s": None, "t_magnet_s": None, "t_first_s": None, "first": None, "beyond_maps": False, "j_eff": 0.0}
    at2 = factor * cont["AT_pk"]
    pt = machine.point(at2, model.speed_rpm, model.stack_mm, clip=True)
    f, share, r100 = machine.copper(min(at2, machine.at_max), pt["gamma_deg"], model.speed_rpm, model.stack_mm)
    dc = r100["p_dc_slot"] + r100["p_dc_end"]
    a_side = _a_side_mm2(model)
    j2 = at2 / (SQRT2 * a_side)
    # the module's loss, AC parts included, as an effective current density in rho J^2 V: AC/DC at 100 C (the ratio hardly moves with the temperature)
    j_eff = j2 * math.sqrt(r100["p_total"] / dc) if dc > 0 else j2
    iron = {"iron_teeth": pt["iron_teeth"], "iron_yoke": pt["iron_yoke"], "iron_rotor": pt["iron_rotor"]}
    o = model.overload(1.0, j_eff, start="cold", wind_limits=tuple(wind_limits), magnet_limit_c=magnet_limit_c, iron=iron, t_max_s=t_max_s)
    t_w = o[f"winding_{wind_limits[0]:.0f}C_s"]
    t_m = o[f"magnet_{magnet_limit_c:.0f}C_s"]
    times = [(t, name) for t, name in ((t_w, "winding"), (t_m, "magnet")) if t is not None]
    t_first, first = min(times) if times else (None, None)
    return {"factor": factor, "AT_pk": at2, "torque": pt["torque"], "t_winding_s": t_w, "t_magnet_s": t_m, "t_first_s": t_first, "first": first, "beyond_maps": pt["beyond_maps"],
            "j_eff": j_eff, "iron_total_w": sum(iron.values())}


def main() -> int:
    t0 = time.perf_counter()
    mach = Machine()
    t_load = time.perf_counter() - t0
    print(f"Machine(): maps, iron-loss map and copper-loss model loaded in {t_load * 1e3:.0f} ms; rated {mach.at_rated:.2f} A-turns, maps up to {mach.at_max:.2f}; "
          f"winding fill {mach.fill:.4f}")
    for name, c in (("natural", Cooling.natural()), ("fixed-speed forced air", Cooling.forced_air()), ("shaft fan, 6 m/s at 3,000 rpm", Cooling.shaft_fan(6.0, 3000.0)),
                    ("liquid jacket", Cooling.liquid())):
        m = model_for(c, 60.0, RATED_RPM, mach)
        t1 = time.perf_counter()
        cl = continuous_limit(m, mach)
        dt = time.perf_counter() - t1
        sh = short_term(m, mach, cl)
        print(f"{name:32s} 60 mm, 3,000 rpm: {cl['AT_pk']:6.1f} A-turns peak ({cl['J_rms']:5.2f} A/mm2), T {cl['torque']:6.3f} N.m, P_em {cl['p_em_w']:7.0f} W, "
              f"Pcu {cl['copper']:6.1f} W, iron {cl['iron']:5.1f} W, eff {100 * cl['efficiency']:5.1f} %, hot spot {cl['hot_spot']:5.1f} C, magnets {cl['magnet']:5.1f} C, "
              f"binding {cl['binding']}; 2x: T {sh['torque'] if sh['torque'] is None else round(sh['torque'], 2)} N.m, first limit {sh['first']} after {sh['t_first_s']} s   "
              f"[{dt * 1e3:.0f} ms]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
