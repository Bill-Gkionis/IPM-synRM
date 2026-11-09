"""optimization.motor_map_thermal: the motor map's current limits from the lumped thermal network (thermal.py): the continuous limit at every speed, and the short-term one.

What it is, and why
    The envelope of motor_map.py needs a current limit, and the limit is a temperature, not a number: the N45 magnets must stay at or under 80 C and the winding's hot spot under
    155 C (class F), and what heats them is the loss at the operating point the current makes. That loss depends on the speed (the iron loss grows as f^1.5 ... f^2, the
    windage as the cube, and the cooling itself depends on the speed: the air gap conducts better as the rotor spins, a shaft fan blows in proportion to it), so the continuous
    limit is a function of the speed, and the network has to be solved at every speed, with the loss at the operating point that limit allows.

The method (per speed)
    1. At the speed n, the current limit I and the copper at T_cu, the point of largest shaft torque inside the current circle and the voltage limit is found (Scan.max_torque):
       (i_d, i_q), and from the maps its iron loss by region (teeth, yoke, rotor).
    2. The network of thermal.py, built for this speed (its air gap, end space and windage depend on it), gets those iron losses, the copper loss as a function of the copper
       temperature (1.5 R(T) I^2, the network iterates on it), and its own bearing and windage estimates (the same formulas as motor_map.Mechanical: checked). It returns the node
       temperatures; the margin is min(155 C - hot spot, 80 C - magnets).
    3. The margin falls as I rises (the copper loss goes as I^2 and the iron loss does not fall), so the continuous limit is the root of margin(I) = 0, found by bisection on I.
       Where no current of the circle can hold the voltage the margin is +infinity (nothing to heat); where the margin is still positive at the end of the maps (3 times the rated
       current of the 10-turn placeholder), the limit is the maps' edge and is flagged: the maps, not the temperature, decide.
    4. The copper's resistance in the voltage equation follows the temperature the network gives (the mean of the slot and end-turn nodes): steps 1 to 3 are repeated with the new
       temperature until it moves less than 2 K (two or three passes).
    The torque at the limit is the envelope's torque there; where it is not positive (the loss torque exceeds the air-gap torque, or the voltage cannot be held) the motor cannot run
    continuously at that speed.

The short-term limit
    The continuous limit at the rated speed with the voltage limit taken away (the thermal rating: what the winding and the magnets can carry whatever the bus), doubled: the
    inverter's peak current. It is a constant over the speed (a hardware limit), and clipped at the maps' edge. `overload_from_cold` gives how long that current can be carried from a
    cold start (everything at the ambient) before the magnets reach 80 C or the winding 155 C, with the network's transient and the iron loss of the maps. From a hot start the
    magnets pass 80 C within seconds, so every short-term figure is from a COLD start.

What is not in it
    The magnets' eddy loss (the network has an input for it; 2D FEMM cannot give it: 0 here), the inverter's heat (outside the motor), the end-winding leakage (not in the maps),
    the demagnetisation of the magnets at large negative i_d (decision EM-6) and anything the maps do not cover (currents beyond 3 times the placeholder's rated current).

Run
    uv run python -m optimization.motor_map_thermal      the thermal rating of the winding study's winding of the 20 and the 60 mm stack in each cooling at 3,000 rpm, against the network's own sustainable()

Numbers to expect (this PC, 2026-11-09; the ampere-turns a coil side, peak, with no voltage limit; this module / the network's sustainable() with its placeholder iron loss and its DC copper
    loss; copper_loss.py's full loss with the strands in a cluster and the diamond end turn; the winding study's winding of each stack, 26 / 4 8 x 0.315 mm on 20 mm and 9 / 4 7 x 0.56 mm on
    60 mm; the fans (b) and (c) sized on each stack with the diamond end turn; the magnets bind in every case but the liquid jacket on 20 mm, which the winding binds)
    shaft fan (a) 6 m/s at 3,000 rpm   20 mm 159.8 / 160.9     60 mm 154.7 / 156.6
    shaft fan (b) 50 W at 16,330 rpm   20 mm 108.5 / 108.5     60 mm  94.6 /  93.3
    shaft fan (c) 200 W at 16,330 rpm  20 mm 108.9 / 108.9     60 mm  95.4 /  94.2
    natural convection                 20 mm  72.5 /  70.8     60 mm  47.6 /  38.9
    fixed 6 m/s external fan           20 mm 158.2 / 159.3     60 mm 154.4 / 156.3
    liquid jacket (parked)             20 mm 288.6 / 289.5     60 mm 342.1 / 344.5
    The differences are the iron loss, the fill and the AC copper loss: the network's placeholder iron loss is about 7.9 W (20 mm) and 23.8 W (60 mm) at 3,000 rpm whatever the current, and
    the maps' loss at the operating point the current allows is 10.3 W and 30.6 W (shaft fan (a): a large current), 8.6 and 24.4 W (fan (b)), 7.5 and 20.8 W (natural convection: little
    armature flux at that small current), 12.6 and 39.5 W (liquid). Where the heat the cooling can shed is small, as in natural convection on 60 mm (21 W of iron loss here against 23.8 W
    there), the allowed current moves a lot (+22 %); where it is large, it moves about 1 %. The study's own ratings of its windings are reproduced to 3e-4 (check (q2) of motor_map_run).
    A copper loss that the network and this module price alike agrees to 2e-7 (check (g)). A continuous limit takes 40 ... 80 ms a speed (a bisection on the current with the maps'
    operating point, the network and the copper temperature iterated), a thermal rating 0.1 to 0.3 s.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

import numpy as np

from optimization import motor_map as mm
from optimization import thermal

RATED_RPM = 3000.0
LEVELS_SEARCH = 7                                      # zoom levels of the operating point inside the bisection (the final envelope solves it again with the full 12)
NO_VOLTAGE_LIMIT = mm.Bus(v_dc=48.0, k=1.0, limit=1.0e9)   # a bus that never limits (the thermal rating of the machine); v_dc stays 48 V for the copper model's PWM ripple


@dataclass(frozen=True)
class ThermalConfig:
    """The thermal side of a case: the cooling (thermal.Cooling.natural() / .forced_air() / .liquid(), and the shaft fan when thermal.py has it), the network's assumptions,
    the two temperature limits and the magnets' eddy loss [W] (an input: 2D FEMM does not give it)."""
    cooling: object
    assume: thermal.Assumptions = thermal.Assumptions()
    wind_limit_c: float = thermal.CLASS_F_C
    magnet_limit_c: float = thermal.MAGNET_LIMIT_C
    magnet_loss_w: float = 0.0

    @property
    def name(self) -> str:
        return getattr(self.cooling, "kind", "?")


def network(m: mm.Machine, cfg: ThermalConfig, rpm: float):
    """The thermal network of this machine's stack and fill, at this speed and cooling (thermal.get_model caches them)."""
    return thermal.get_model(None, m.stack_mm, cfg.cooling, float(rpm), m.fill, cfg.assume)


def state(m: mm.Machine, bus: mm.Bus, cfg: ThermalConfig, rpm: float, i_d: float, i_q: float, t_cu: float):
    """The network's steady state with the losses of the maps at the operating point (i_d, i_q) and speed rpm: returns (margin [K], SteadyResult, point). The copper loss is a
    function of the temperature (m.copper.loss at T), the iron loss is the maps' by region, the bearing and windage the network's own."""
    pt = m.point(i_d, i_q, rpm, bus, t_cu)
    th = network(m, cfg, rpm)

    def copper(t):
        return float(m.copper.loss(m, i_d, i_q, rpm, t, bus.v_dc)["total"])
    r = th.steady_state(thermal.Losses(copper=copper, iron_teeth=float(pt["P_fe_teeth"]), iron_yoke=float(pt["P_fe_yoke"]), iron_rotor=float(pt["P_fe_rotor"]),
                                       magnet=cfg.magnet_loss_w), fill=m.fill, wind_limit_c=cfg.wind_limit_c, magnet_limit_c=cfg.magnet_limit_c)
    margin = -1.0e3 if (r.runaway or not r.converged) else min(r.margin_winding, r.margin_magnet)
    return margin, r, pt


def copper_temperature(m: mm.Machine, r) -> float:
    """The mean copper temperature [C] the network gives: slot and end-turn nodes weighted by the share of the copper (and of the resistance) in each."""
    f = m.copper.slot_fraction(m)
    return float(f * r.T["W"] + (1.0 - f) * r.T["E"])


def continuous_limit(sc: mm.Scan, m: mm.Machine, bus: mm.Bus, cfg: ThermalConfig, rpm: float, i_cap: float, t_cu: float, i_guess: float | None = None, basis: str = "shaft",
                     tol: float = 2e-4) -> dict:
    """The largest current limit I <= i_cap at which the network keeps the magnets and the winding inside their limits at this speed (module docstring, steps 1 to 3).
    Returns {"I", "limited" ("temperature", "maps", "none"), "sol" (the max-torque solution at I or None), "margin", "result" (SteadyResult or None), "t_cu_new"}."""
    c = sc.ctx(rpm, bus, t_cu, basis)
    grid = sc.grid(c)
    cache: dict = {}

    def margin(I: float) -> float:
        sol = sc.max_torque(c, I, grid, levels=LEVELS_SEARCH)
        if sol is None:
            cache[I] = (np.inf, None, None)
            return np.inf
        mg_, r, _ = state(m, bus, cfg, rpm, sol["i_d"], sol["i_q"], t_cu)
        cache[I] = (mg_, sol, r)
        return mg_
    lo, hi = 0.0, float(i_cap)
    if margin(hi) >= 0.0:
        I, why = hi, "maps"
    else:
        if i_guess is not None and 0.0 < i_guess < hi:                    # warm start: the limit changes slowly with the speed
            a, b = max(0.0, 0.8 * i_guess), min(hi, 1.25 * i_guess)
            if margin(a) >= 0.0:
                lo = a
                if margin(b) < 0.0:
                    hi = b
            else:
                hi = a
        while hi - lo > tol * max(i_cap, 1.0):
            mid = 0.5 * (lo + hi)
            if margin(mid) >= 0.0:
                lo = mid
            else:
                hi = mid
        I, why = lo, "temperature"
        if I not in cache:
            margin(I)
    mg_, sol, r = cache[I]
    if sol is None or mg_ < 0.0:                       # nothing holds the voltage, or even the smallest current that does (or none: the loss of the iron and the air alone) overheats
        why, I = "none", 0.0
    return {"I": I, "limited": why, "sol": sol, "margin": mg_, "result": r, "t_cu_new": copper_temperature(m, r) if r is not None else t_cu}


def continuous_limits(m: mm.Machine, bus: mm.Bus, cfg: ThermalConfig, rpms, scan: mm.Scan | None = None, t_cu0: float = 100.0, passes: int = 3, t_tol: float = 2.0,
                      i_cap: float | None = None, **scan_kw) -> dict:
    """The continuous current limit at every speed in rpms, with the copper temperature the network gives (module docstring). Returns arrays over the speeds: I (the limit [A peak]),
    t_cu, limited, hot_spot, magnet (the temperatures [C]), margin (K), binding, converged (the temperature loop), and the scan used."""
    rpms = np.asarray(rpms, float)
    cap = min(i_cap if i_cap is not None else m.i_cap, m.i_cap)
    sc = scan or mm.Scan(m, i_max=cap, **scan_kw)
    n = len(rpms)
    out = {"rpm": rpms, "I": np.zeros(n), "t_cu": np.full(n, t_cu0), "limited": np.empty(n, object), "hot_spot": np.full(n, np.nan), "magnet": np.full(n, np.nan),
           "margin": np.full(n, np.nan), "binding": np.empty(n, object), "converged": np.zeros(n, bool), "scan": sc}
    prev, prev_t = None, None
    for k, rpm in enumerate(rpms):
        t = t_cu0 if prev_t is None else prev_t                           # the copper temperature changes slowly with the speed: start from the last one
        for it in range(passes):
            res = continuous_limit(sc, m, bus, cfg, rpm, cap, t, prev)
            if abs(res["t_cu_new"] - t) < t_tol:
                out["converged"][k] = True
                break
            t = res["t_cu_new"]
        out["I"][k], out["t_cu"][k], out["limited"][k] = res["I"], t, res["limited"]
        out["margin"][k] = res["margin"]
        if res["result"] is not None:
            out["hot_spot"][k], out["magnet"][k], out["binding"][k] = res["result"].hot_spot, res["result"].magnet, res["result"].binding
        else:
            out["binding"][k] = "none"
        prev = res["I"] if res["I"] > 0 else None
        prev_t = t if res["limited"] != "none" else prev_t
    return out


def thermal_rating(m: mm.Machine, cfg: ThermalConfig, rpm: float = RATED_RPM, t_cu: float = 100.0, scan: mm.Scan | None = None) -> dict:
    """The continuous current [A peak] at the rated speed with no voltage limit (the thermal rating of the machine, whatever the bus), its ampere-turns a coil side (N i / a, peak),
    its current density, and the temperatures. The peak (short-term) current is a multiple of it."""
    sc = scan or mm.Scan(m, i_max=m.i_cap)
    res = None
    t = t_cu
    for _ in range(3):
        res = continuous_limit(sc, m, NO_VOLTAGE_LIMIT, cfg, rpm, m.i_cap, t)
        if abs(res["t_cu_new"] - t) < 1.0:
            break
        t = res["t_cu_new"]
    I = res["I"]
    out = {"I": I, "AT_pk": m.w.turns * I / m.w.paths, "J_rms": m.w.density_from_current(I), "limited": res["limited"], "t_cu": t}
    if res["result"] is not None:
        out.update({"hot_spot": res["result"].hot_spot, "magnet": res["result"].magnet, "binding": res["result"].binding, "P_cu": res["result"].copper_loss,
                    "sources": res["result"].sources})
    return out


def overload_from_cold(m: mm.Machine, bus: mm.Bus, cfg: ThermalConfig, rpm: float, i_lim: float, sc: mm.Scan | None = None, t_cu: float = 100.0, t_max_s: float = 3600.0, *,
                       iron: dict | None = None, j_rms: float | None = None) -> dict:
    """How long the current limit i_lim can be carried at this speed from a COLD start (everything at the ambient), with the iron loss of the maps at the operating point it
    allows: the network's transient (thermal.ThermalModel.transient, the step that thermal.ThermalModel.overload takes) until the magnets reach their limit and the winding its own.
    Returns the operating point, the current density and the times [s] ("magnet_s", "winding_s": None if not reached in t_max_s). iron (a dict of iron_teeth, iron_yoke, iron_rotor [W]) and
    j_rms [A/mm2] replace the maps' operating point: with the network's own placeholder iron loss and 2 x its sustainable current density this is exactly ThermalModel.overload(2, j, "cold")
    (tests/test_motor_map.py holds it to that)."""
    th = network(m, cfg, rpm)
    out: dict = {"feasible": True}
    if iron is None or j_rms is None:
        sc = sc or mm.Scan(m, i_max=min(i_lim, m.i_cap))
        sol = sc.max_torque(sc.ctx(rpm, bus, t_cu, "shaft"), min(i_lim, m.i_cap))
        if sol is None:
            return {"feasible": False}
        pt = m.point(sol["i_d"], sol["i_q"], rpm, bus, t_cu)
        iron = iron or {"iron_teeth": float(pt["P_fe_teeth"]), "iron_yoke": float(pt["P_fe_yoke"]), "iron_rotor": float(pt["P_fe_rotor"])}
        j_rms = m.w.density_from_current(sol["I"]) if j_rms is None else j_rms
        out.update({"I": sol["I"], "T_shaft": float(pt["T_shaft"]), "i_d": sol["i_d"], "i_q": sol["i_q"]})
    wl, ml = cfg.wind_limit_c, cfg.magnet_limit_c

    def done(t, T, hot):
        return hot > wl + 1.0 and (T[thermal.IDX["M"]] > ml + 0.5 or t > 3600.0)
    r = th.transient(t_max_s, j_rms, start="cold", fill=m.fill, rpm=rpm, iron=iron, dt_s=1.0, dt_max_s=10.0, stop_when=done)

    def first(x, limit):
        idx = np.nonzero(x > limit)[0]
        if len(idx) == 0:
            return None
        i = int(idx[0])
        if i == 0:
            return 0.0
        f = (limit - x[i - 1]) / (x[i] - x[i - 1])
        return float(r.t[i - 1] + f * (r.t[i] - r.t[i - 1]))
    out.update({"J_rms": j_rms, "magnet_s": first(r.magnet, ml), "winding_s": first(r.hot_spot, wl), "magnet_rise_60s": float(np.interp(60.0, r.t, r.magnet) - r.magnet[0]), "iron": iron,
                "end_time_s": float(r.t[-1])})
    return out


def main() -> int:
    from optimization import motor_map_stack as st
    asm = thermal.Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6)
    copper = mm.FullCopper(arrangement="cluster")
    print("the continuous current at 3,000 rpm with no voltage limit (the thermal rating), the winding study's winding of each stack, copper_loss.py's full loss with the diamond end turn: "
          "the network's own sustainable() against this module")
    coolings = (("shaft fan (a) 6 m/s", lambda x: thermal.Cooling.shaft_fan(6.0, 3000.0)),
                ("shaft fan (b) 50 W", lambda x: thermal.size_shaft_fan(50.0, mm.TOP_RPM, 3000.0, x, assume=asm)["cooling"]),
                ("shaft fan (c) 200 W", lambda x: thermal.size_shaft_fan(200.0, mm.TOP_RPM, 3000.0, x, assume=asm)["cooling"]),
                ("natural convection", lambda x: thermal.Cooling.natural()), ("external fan 6 m/s", lambda x: thermal.Cooling.forced_air()), ("liquid (parked)", lambda x: thermal.Cooling.liquid()))
    for name, make in coolings:
        for stack in (20.0, 60.0):
            m = mm.Machine(st.OPTIMIZED_WINDINGS[stack], stack, copper=copper)
            cooling = make(stack)
            cfg = ThermalConfig(cooling, asm)
            mine = thermal_rating(m, cfg)
            net = thermal.get_model(None, stack, cooling, RATED_RPM, m.fill, cfg.assume).sustainable()
            print(f"  {name:22s} {stack:4.0f} mm: mine {mine['I']:6.2f} A ({mine['AT_pk']:6.1f} A-turns, {mine['J_rms']:5.2f} A/mm2, {mine['limited']}, {mine.get('binding')}), "
                  f"network {net['AT_pk']:6.1f} A-turns ({net['J_rms']:5.2f} A/mm2, {net['binding']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
