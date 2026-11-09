"""optimization.design: one candidate motor in, its scores out. The black box that every sweep, screening and optimizer of the study calls.

What it is, and why there is only one
    The study asks one question thousands of times: "this geometry: how much torque, and how much ripple?". If each script answered
    it its own way (another current angle, other rotor positions, another ripple measure), two scripts could disagree about the same
    rotor and nobody could tell which one was right. So the answer is written once, here:
      decode   a candidate is a dict of motorgeo.Params overrides, on the band-free rotor Params(barrier="none") (R11): {} is the
               reference design, {"mag_x": 39.2} the reference with its magnets 0.21 mm deeper
      screen   build it and check it in Python (shapely, about 10 ms, no FEMM): a candidate that cannot be built, or breaks a limit,
               gets its reason and its constraint values, and never reaches FEMM
      score    a FEMM of a femm_pool worker (one heavy run at a time, R10) writes its model, finds its own MTPA angle at the scoring
               current, and turns the rotor through one ripple period at that angle
    and returns one row: mean torque, ripple, the angle, the flux linkages, the power factor, the constraint values.

The variables (the names are motorgeo.Params fields; the study note maps them to the notes' words)
    the rotor tier    mag_x       the magnet slot's centre along the d-axis [mm]: with mag_y, the radius of the V's magnets
                      mag_y       the slot centre's distance from the d-axis [mm]: half the spacing of the two magnets of a pole
                      v_half_deg  the magnet's long side against the d-axis [deg] (60 is a V of 120 degrees)
                      bridge_mm   (instead of mag_x, the mapping layer) the outer bridge's thickness [mm]; decode() turns it into
                                  the mag_x that leaves it, so the build guard becomes the plain bound bridge_mm >= 1.4
    the stator tier   slot_half_deg  half the angle of the radial-walled slot body [deg]: the tooth is 10 - 2 x this, in degrees
                      slot_open   the slot opening at the bore [mm]
                      h_open      the straight part of the opening, the tooth tip's thinnest height [mm]
                      r_shoe      where the tooth's shoe ends and the slot walls start [mm]
                      r_wall      where the walls end and the round slot bottom starts [mm]: the slot's depth, and so the yoke
    Fixed: the envelope (r_so, r_si, the air gap, the shaft, the stack), the magnets and their slots (slot_w, slot_t, slot_fillet:
    the magnets are bought), the winding layout (36 slots, 8 poles, a coil pitch of 4) and the bridge guard (1.4 mm, decision 17).

The constraints (G <= 0 means allowed; pymoo's convention, so every value is in mm or degrees of room taken away)
    guard_bridge   the build guard: a magnet slot's farthest point must stay 1.4 mm (Params.bridge) inside the rotor surface.
                   motorgeo raises GeometryError past it; here it is computed before building, so an optimizer is told how far out it is
    centre_rib, q_rib, outer_bridge   motorgeo.check's limits (1.0 mm of steel), tooth_min (2.0 mm) and yoke (3.0 mm), the same
    shoe_taper     the shoe must rise at least 0.2 mm between the opening's straight part and the slot walls (r_shoe vs r_si + h_open)
    open_vs_slot   the opening must be at least 0.3 degrees narrower (each side) than the slot body
    wall_length    the slot walls must be at least 1.0 mm long (r_wall vs r_shoe)
    layer_split    the cut between the two coil layers must lie on the straight walls (r_mid <= r_wall; the stator sweeps found it)
    Every region's polygon must also be a valid shape; a stator drawn inside out is refused with its region's name.
    The study's limits are stricter than motorgeo.check's own: tooth_min 3.0 mm (not 2.0) and yoke 3.5 mm (not 3.0), because the
    sweeps found the optimum sitting on whatever these are set to (working values, still to be decided). A value within TOL (1e-9 mm)
    past its limit is round-off and counts as on it.

The scoring point (a working choice, written in the study note as an open decision at OP-1)
    10 A peak with the 10-turn placeholder winding: 100 ampere-turns per coil side, the bench's and the demo's point.
    scoring="AT"  every candidate gets those 100 ampere-turns. Right for the rotor tier, where the slot does not change.
    scoring="J"   every candidate gets the reference slot's current density: I = 10 A x (its coil area / the reference's coil area).
                  It stops the cheat of "AT" (at fixed ampere-turns an optimizer gets a wide tooth's lower saturation for free while
                  the current density in the smaller slot climbs where nobody looks), but the stator sweeps found its mirror image:
                  at fixed J the torque follows the slot's area exactly, so an optimizer grows the slot to the tooth and yoke limits
                  while the copper loss (Pcu_rel, logged) grows with it, up to 1.62 times, and nothing charges for it.
    scoring="P"   every candidate gets the reference's copper loss: I = 10 A x sqrt(its coil area / the reference's). The loss is
                  the heat, which is what the cooling has to carry: the torque at a fixed loss is the motor constant's view
                  (torque per root watt), and the iron-against-copper trade of the stator then has a real optimum. The study's
                  stator and joint runs use it (the mean turn's length is taken as fixed: deeper slots lengthen the end turns a
                  little, which this ignores).

How one candidate is scored (23 solves; about 8 s on one worker with the coarse mesh)
    1. no load: two rotor positions, 0 and 30 electrical degrees. Their mean d-axis flux linkage is psi_pm, the magnets' share.
    2. its own MTPA angle: the mean torque of the same two positions (30 electrical degrees apart, they cancel the 6th, 18th and
       30th harmonics of the ripple; B7) at three current angles, 10, 25 and 40 degrees. A parabola through the three gives the
       best angle (on the reference it lands within half a degree of B1's 22.9).
    3. the ripple: 15 rotor positions over one ripple period (60 electrical degrees, 4 degrees apart) at that angle.
       The mean torque is their mean. The ripple is the peak-to-peak of the waveform rebuilt from the 15 samples by their Fourier
       series (exact up to the 42nd harmonic), not of the 15 samples themselves: a peak that falls between two samples would
       otherwise look like a smaller ripple, and an optimizer finds such loopholes and lives in them.
    The power factor is the internal one, without the resistance: cos phi = (psi_d i_q - psi_q i_d) / (|psi| |i|), which is the
    torque divided by the largest torque that flux and that current could give.

The two objectives every optimizer of the study minimises (objectives()): -T_mean [N.m] and the ripple in per cent of T_mean.

Run
    Not on its own. A script opens a pool and hands it to an Evaluator:
        with femm_queue.femm_pool("optimization: rotor sweeps", task_timeout_s=design.TASK_TIMEOUT_S) as pool:
            evaluate = design.Evaluator(pool, run="sweeps_rotor", scoring="AT")
            rows = evaluate([{}, {"mag_x": 39.2}], stage="mag_x")
    uv run python -m optimization.check_design   checks this module on the reference design.

Self-check (the answers follow).
  1. Why is the ripple not simply max(T) - min(T) of the 15 samples?
  2. The stator tier is scored with scoring="J". What would an optimizer do with the tooth width under scoring="AT"?
  3. A candidate's magnet slot comes 0.3 mm closer to the rotor surface than the guard allows. What does screen() return, and why
     is it better than letting motorgeo raise its GeometryError?
Answers.
  1. The waveform's peaks fall between samples for some rotors and on a sample for others, so the sampled peak-to-peak depends on
     where the samples happen to sit (by up to about 20 % for the 18th harmonic at five samples a period). An optimizer would learn
     to prefer rotors whose peaks hide between samples. The Fourier series rebuilds the waveform between the samples, so the peaks
     are found wherever they are.
  2. Widen the teeth: they saturate less, so the torque rises, and the slot that shrinks costs nothing because its 100 ampere-turns
     stay the same. In the real motor the current density in that smaller slot would rise, and with it the copper loss and the heat.
  3. feasible False, the reason, and guard_bridge = +0.3 (mm): the optimizer learns how far outside the allowed region it is and
     which way is back. A GeometryError says only "no", so every rejected design would look equally bad.
"""
from __future__ import annotations

import csv
import json
import math
import os
import time
import traceback
import uuid
from pathlib import Path

import numpy as np

import motorgeo as mg
import paths

# ---------------------------------------------------------------------------------------------- the fixed choices
I_REF = 10.0                                   # A peak at the scoring point (10 turns a coil side: 100 ampere-turns)
GAMMAS = (10.0, 25.0, 40.0)                    # current angles of the MTPA search [deg from the q-axis toward -d]
ANGLE_POSITIONS = (0.0, 30.0)                  # electrical degrees: the pair that cancels the 6th, 18th and 30th harmonics
RIPPLE_POSITIONS = tuple(4.0 * k for k in range(15))   # one ripple period, 60 electrical degrees, 4 degrees apart
TASK_TIMEOUT_S = 300                           # one candidate is 23 solves: about 8 s coarse, 15 s fine; a hung FEMM is ended
MOTORGEO_LIMITS = {"centre_rib": 1.0, "q_rib": 1.0, "outer_bridge": 1.0, "tooth_min": 2.0, "yoke": 3.0}   # motorgeo.check's own
LIMITS = {**MOTORGEO_LIMITS, "tooth_min": 3.0, "yoke": 3.5}   # the study's working limits (the sweeps: at 2.0 and 3.0 mm the
                                                          # optimum sits on a fragile tooth and a 3 mm yoke): a decision is pending
TOL = 1e-9                                     # mm: a constraint this far past its limit is round-off, not a violation
SHOE_TAPER_MIN, OPEN_MARGIN_DEG, WALL_MIN = 0.2, 0.3, 1.0

ROTOR_VARS = ("mag_x", "mag_y", "v_half_deg")
STATOR_VARS = ("slot_half_deg", "slot_open", "h_open", "r_shoe", "r_wall")
ALL_VARS = ROTOR_VARS + STATOR_VARS
G_NAMES = ("guard_bridge", "centre_rib", "q_rib", "outer_bridge", "tooth_min", "yoke", "shoe_taper", "open_vs_slot", "wall_length",
           "layer_split")


def params(overrides: dict) -> mg.Params:
    """The candidate as motorgeo.Params: the band-free reference with these fields changed (an unknown name is a TypeError: a bug)."""
    return mg.Params(barrier="none", **overrides)


REFERENCE = params({})


def currents(amps: float, gamma_deg: float) -> tuple[float, float]:
    """(i_d, i_q) for a current of `amps` peak at the angle gamma from the q-axis toward -d (the convention of learn/dq_torque)."""
    g = math.radians(gamma_deg)
    return -amps * math.sin(g), amps * math.cos(g)


# ---------------------------------------------------------------------------------------------- screening (Python only)
def slot_reach(p: mg.Params) -> float:
    """How far from the rotor's centre a magnet slot reaches [mm]: its farthest point, on a corner fillet. The same sum as the
    build guard in motorgeo._pole_slot (the slot centre, plus s along the magnet and t across it, plus the fillet's radius)."""
    b = math.radians(p.v_half_deg)
    u, n = (math.cos(b), math.sin(b)), (math.sin(b), -math.cos(b))
    hw, ht, f = p.slot_w / 2, p.slot_t / 2, p.slot_fillet
    centres = [(-hw + f, ht - f), (hw - f, ht - f), (hw - f, -ht + f), (-hw + f, -ht + f)]
    return max(math.hypot(p.mag_x + s * u[0] + t * n[0], p.mag_y + s * u[1] + t * n[1]) for s, t in centres) + f


def bridge_of(p: mg.Params) -> float:
    """The outer bridge the slot leaves [mm]: the rotor's radius less the slot's reach (1.558 mm on the reference)."""
    return p.r_ro - slot_reach(p)


def mag_x_for_bridge(bridge_mm: float, overrides: dict | None = None) -> float:
    """The mag_x [mm] that leaves an outer bridge of `bridge_mm`, the other fields as in `overrides`. Bisection: every corner of the
    slot sits at x > 0, so the bridge shrinks steadily as mag_x grows (60 halvings of 20..r_ro: far below a nanometre)."""
    base = {k: v for k, v in (overrides or {}).items() if k != "mag_x"}
    lo, hi = 20.0, REFERENCE.r_ro
    for _ in range(60):
        mid = (lo + hi) / 2
        if bridge_of(params({**base, "mag_x": mid})) > bridge_mm:
            lo = mid
        else:
            hi = mid
    return lo          # the side where the bridge is still wider than asked: on the guard, never 1e-15 mm past it


def decode(candidate: dict) -> dict:
    """Optimizer variables -> Params overrides: the mapping layer. Every name is a Params field except bridge_mm, the outer
    bridge's thickness [mm], which replaces mag_x. Given as bridge_mm, the build guard is just the bound bridge_mm >= 1.4: an
    optimizer can sit on it, and no candidate is lost by stepping past it. ValueError if both bridge_mm and mag_x are given."""
    c = {k: float(v) for k, v in candidate.items()}
    if "bridge_mm" in c:
        if "mag_x" in c:
            raise ValueError("give bridge_mm or mag_x, not both")
        c["mag_x"] = mag_x_for_bridge(c.pop("bridge_mm"), c)
    return c


def coil_area(metrics: dict) -> float:
    """One coil layer's area [mm^2] (motorgeo.check: coil_0_0; the two layers of a slot are equal by construction)."""
    return float(metrics["coil_layer_area_mm2"])


def _stator_room(p: mg.Params) -> dict:
    """The stator's shape constraints, computed from the numbers alone (motorgeo does not check these, and would draw nonsense).
    layer_split: motorgeo cuts the slot into two coil layers of equal area at r_mid, and that cut must lie on the straight walls
    (r_mid <= r_wall). Past it the second layer's outline crosses itself; the stator sweeps found it as "invalid shapes: coil_k_1"
    with every other value inside its limit, so an optimizer could not tell how far out it was. Same sums as motorgeo._slot."""
    a_o = math.degrees(math.asin(min(1.0, p.slot_open / 2 / p.r_si)))
    a_s = math.radians(p.slot_half_deg)
    r_b = p.r_wall * math.tan(a_s)
    area = 0.5 * (p.r_wall ** 2 - p.r_shoe ** 2) * math.sin(2 * a_s) + r_b ** 2 / 2 * (math.pi + 2 * a_s + math.sin(2 * a_s))
    r_mid = math.sqrt(max(0.0, p.r_shoe ** 2 + area / math.sin(2 * a_s))) if a_s > 0 else float("inf")
    return {"shoe_taper": SHOE_TAPER_MIN - (p.r_shoe - (p.r_si + p.h_open)),
            "open_vs_slot": OPEN_MARGIN_DEG - (p.slot_half_deg - a_o),
            "wall_length": WALL_MIN - (p.r_wall - p.r_shoe),
            "layer_split": r_mid - p.r_wall}


def screen(overrides: dict) -> dict:
    """Build and check a candidate without FEMM. Returns {"feasible", "reason", "G" (every name of G_NAMES), "metrics", "params"}.
    G values that could not be computed (the drawing failed first) are 0.0: unknown, so not counted against the candidate."""
    p = params(overrides)
    G = dict.fromkeys(G_NAMES, 0.0)
    G["guard_bridge"] = p.bridge - bridge_of(p)
    G.update(_stator_room(p))
    out = {"feasible": False, "reason": "", "G": G, "metrics": {}, "params": p}
    pre = [f"{k} {v:+.4f}" for k, v in G.items() if v > TOL]
    if pre:
        out["reason"] = "before drawing: " + ", ".join(pre)
        return out
    try:
        regions = mg.build(p)
    except mg.GeometryError as error:
        out["reason"] = f"GeometryError: {error}"
        return out
    bad = [r.name for r in regions if not r.polygon().is_valid or r.polygon().area <= 0]
    if bad:
        out["reason"] = "invalid shapes: " + ", ".join(bad[:4]) + (" ..." if len(bad) > 4 else "")
        return out
    metrics, violations = mg.check(p, regions)
    for k, limit in LIMITS.items():
        G[k] = limit - metrics[k]
    out["metrics"] = metrics
    out["feasible"] = all(v <= TOL for v in G.values())            # the study's limits (LIMITS), which include motorgeo's
    stricter = [f"{k} = {metrics[k]:.2f} < {LIMITS[k]} (the study's limit)" for k in LIMITS if G[k] > TOL and LIMITS[k] > MOTORGEO_LIMITS[k]]
    out["reason"] = "; ".join([*violations, *stricter])
    return out


def scoring_amps(metrics: dict, scoring: str) -> float:
    """The scoring current [A peak] of a candidate: I_REF at fixed ampere-turns ("AT"); I_REF x its coil area / the reference's at
    the reference's current density ("J"); or I_REF x sqrt(its coil area / the reference's) at the reference's copper loss ("P":
    the loss is I^2 R and R goes as 1 / the coil area at fixed turns and fill, so an equal loss needs I^2 to grow with the area)."""
    if scoring == "AT":
        return I_REF
    if scoring == "J":
        return I_REF * coil_area(metrics) / REFERENCE_COIL_AREA
    if scoring == "P":
        return I_REF * math.sqrt(coil_area(metrics) / REFERENCE_COIL_AREA)
    raise ValueError(f"scoring={scoring!r}: 'AT', 'J' or 'P'")


# ---------------------------------------------------------------------------------------------- the waveform's measures
def best_angle(gammas, torques) -> float:
    """The current angle of the parabola's top through (gamma, torque) [deg], kept inside the three angles. If the three do not
    bend downward (no top), the best of the three."""
    a2, a1, _ = np.polyfit(gammas, torques, 2)
    if a2 >= 0:
        return float(gammas[int(np.argmax(torques))])
    return float(np.clip(-a1 / (2 * a2), gammas[0], gammas[-1]))


def fourier_pkpk(samples, upsample: int = 200) -> float:
    """Peak-to-peak of the periodic waveform whose one period the samples cover (equally spaced, the last one a step before the
    period ends), rebuilt between the samples by its Fourier series: y(t) = c0 + 2 Re sum_k c_k exp(2 pi i k t)."""
    y = np.asarray(samples, float)
    n = len(y)
    c = np.fft.rfft(y) / n
    t = np.arange(n * upsample) / (n * upsample)
    k = np.arange(len(c))
    weight = np.full(len(c), 2.0)
    weight[0] = 1.0
    if n % 2 == 0:
        weight[-1] = 1.0                         # an even count has a Nyquist term, which appears once
    wave = (weight[:, None] * (c[:, None] * np.exp(2j * np.pi * k[:, None] * t[None, :])).real).sum(axis=0)
    return float(wave.max() - wave.min())


# ---------------------------------------------------------------------------------------------- the FEMM side (a pool worker)
def score(task: dict) -> dict:
    """One task of the pool: draw the candidate on this worker's FEMM and score it (the steps of the module docstring).
    Returns the scores, or {"status": "error", "error": ...} if FEMM or the drawing failed: one bad candidate must not end the run."""
    import femm
    started = time.perf_counter()
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
    try:
        p = params(task["overrides"])
        regions = mg.build(p)
        mg.to_femm(p, regions, fem, fine=bool(task.get("fine")))
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        amps = float(task["amps"])
        no_load = [model.solve(th, 0.0, 0.0) for th in ANGLE_POSITIONS]
        angle_T = []
        for g in GAMMAS:
            i_d, i_q = currents(amps, g)
            angle_T.append(float(np.mean([model.solve(th, i_d, i_q)["T_age"] for th in ANGLE_POSITIONS])))
        gamma = best_angle(GAMMAS, angle_T)
        i_d, i_q = currents(amps, gamma)
        rows = [model.solve(th, i_d, i_q) for th in RIPPLE_POSITIONS]
        femm.mi_loadsolution()
        nodes = int(femm.mo_numnodes())
        femm.mo_close()
    except Exception as error:                    # a failed mesh or solve: say so in the row, keep the run going
        return {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=3),
                "seconds": time.perf_counter() - started}
    finally:
        try:
            femm.mi_close()
        except Exception:
            pass
        for ext in ((".ans",) if task.get("keep") else (".ans", ".fem")):
            try:
                os.remove(fem[:-4] + ext)
            except OSError:
                pass
    T = [r["T_age"] for r in rows]
    T_mean = float(np.mean(T))
    psi_pm = float(np.mean([r["psi_d"] for r in no_load]))
    psi_d, psi_q = float(np.mean([r["psi_d"] for r in rows])), float(np.mean([r["psi_q"] for r in rows]))
    pp = fourier_pkpk(T)
    pf = (psi_d * i_q - psi_q * i_d) / (math.hypot(psi_d, psi_q) * math.hypot(i_d, i_q)) if amps > 0 else float("nan")
    L_d = (psi_d - psi_pm) / i_d if abs(i_d) > 1e-9 else float("nan")
    L_q = psi_q / i_q if abs(i_q) > 1e-9 else float("nan")
    return {"status": "ok", "amps": amps, "gamma": gamma, "angle_T": angle_T, "T": T, "T_mean": T_mean,
            "ripple_pp": pp, "ripple_pp_samples": float(max(T) - min(T)), "ripple_pct": 100.0 * pp / T_mean if T_mean else float("nan"),
            "psi_pm": psi_pm, "psi_d": psi_d, "psi_q": psi_q, "i_d": i_d, "i_q": i_q, "L_d": L_d, "L_q": L_q,
            "saliency": L_q / L_d if L_d and not math.isnan(L_d) else float("nan"), "pf": pf,
            "T_dq_mean": float(np.mean([r["T_dq"] for r in rows])), "nodes": nodes, "fine": bool(task.get("fine")),
            "seconds": time.perf_counter() - started}


# ---------------------------------------------------------------------------------------------- the main process
REFERENCE_COIL_AREA = coil_area(screen({})["metrics"])

FIELDS = ["run", "stage", "gen", "n", *ALL_VARS, "bridge_mm", "feasible", "reason", *(f"G_{g}" for g in G_NAMES), "scoring", "fine",
          "amps", "coil_area_mm2", "J_rel", "Pcu_rel", "gamma", "T_mean", "ripple_pp", "ripple_pp_samples", "ripple_pct", "psi_pm",
          "psi_d", "psi_q", "L_d", "L_q", "saliency", "pf", "T_dq_mean", "angle_T", "T", "nodes", "seconds", "cached", "overrides"]


def out_dir(run: str) -> Path:
    """out\\optimization\\<run>\\, made if it is not there."""
    folder = paths.path("out") / "optimization" / run
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def _key(overrides: dict, scoring: str, fine: bool) -> tuple:
    return (tuple(sorted((k, round(float(v), 6)) for k, v in overrides.items())), scoring, bool(fine))


def objectives(row: dict) -> tuple[float, float]:
    """The two objectives every optimizer of the study minimises: (-mean torque [N.m], ripple [% of the mean torque])."""
    return -row["T_mean"], row["ripple_pct"]


class Evaluator:
    """Scores lists of candidates on an open femm_pool, writes every candidate as a CSV row as soon as its batch is back, and
    remembers what it has scored (an optimizer that proposes the same design twice gets the first answer, marked cached)."""

    def __init__(self, pool, run: str, scoring: str = "AT", fine: bool = False, keep: bool = False, csv_name: str = "evaluations.csv"):
        self.pool, self.run, self.scoring, self.fine, self.keep = pool, run, scoring, fine, keep
        self.dir = out_dir(run)
        self.fem_dir = self.dir / "fem"
        self.csv = self.dir / csv_name
        self.cache: dict = {}
        self.n = 0
        if not self.csv.exists():
            with open(self.csv, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(FIELDS)
        else:                                       # appending to a file written with other columns would scramble it
            with open(self.csv, newline="", encoding="utf-8") as f:
                header = next(csv.reader(f), [])
            if header != FIELDS:
                raise RuntimeError(f"{self.csv} was written with other columns (design.py has changed since): use another run name")

    def __call__(self, candidates: list[dict], stage: str = "", gen: int | str = "") -> list[dict]:
        rows: list[dict | None] = [None] * len(candidates)
        todo = []
        for i, overrides in enumerate(candidates):
            overrides = decode(overrides)                 # bridge_mm, if given, becomes mag_x
            key = _key(overrides, self.scoring, self.fine)
            if key in self.cache:
                rows[i] = {**self.cache[key], "cached": True}
                continue
            s = screen(overrides)
            row = {"overrides": overrides, "params": s["params"], "feasible": s["feasible"], "reason": s["reason"], "G": s["G"],
                   "metrics": s["metrics"], "scoring": self.scoring, "fine": self.fine, "cached": False}
            if s["feasible"]:
                row["amps"] = scoring_amps(s["metrics"], self.scoring)
                todo.append((i, key, row))
            else:
                self.cache[key] = row
            rows[i] = row
        tasks = [{"tag": uuid.uuid4().hex, "overrides": row["overrides"], "fem_dir": str(self.fem_dir), "amps": row["amps"],
                  "fine": self.fine, "keep": self.keep} for _, _, row in todo]
        results = self.pool.map(score, tasks) if tasks else []
        for (i, key, row), result in zip(todo, results):
            if isinstance(result, dict) and result.get("status") == "ok":
                row.update(result)
            else:
                said = result.get("error") if isinstance(result, dict) else repr(result)
                row.update({"feasible": False, "reason": f"FEMM: {said or result}"})
            self.cache[key] = row
        self._write(rows, stage, gen)
        return rows

    def _write(self, rows: list[dict], stage, gen) -> None:
        with open(self.csv, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            for row in rows:
                self.n += 1
                p, m = row["params"], row.get("metrics") or {}
                area = coil_area(m) if m else float("nan")
                amps = row.get("amps", float("nan"))
                j_rel = (amps / area) / (I_REF / REFERENCE_COIL_AREA) if m else float("nan")
                pcu_rel = (amps / I_REF) ** 2 * (REFERENCE_COIL_AREA / area) if m else float("nan")
                values = {"run": self.run, "stage": stage, "gen": gen, "n": self.n, **{v: getattr(p, v) for v in ALL_VARS},
                          "bridge_mm": bridge_of(p), "feasible": row["feasible"], "reason": row.get("reason", ""),
                          **{f"G_{g}": row["G"][g] for g in G_NAMES}, "scoring": row["scoring"], "fine": row["fine"], "amps": amps,
                          "coil_area_mm2": area, "J_rel": j_rel, "Pcu_rel": pcu_rel, "cached": row.get("cached", False),
                          "overrides": json.dumps(row["overrides"]), "angle_T": json.dumps(row.get("angle_T", [])),
                          "T": json.dumps(row.get("T", []))}
                for k in ("gamma", "T_mean", "ripple_pp", "ripple_pp_samples", "ripple_pct", "psi_pm", "psi_d", "psi_q", "L_d", "L_q",
                          "saliency", "pf", "T_dq_mean", "nodes", "seconds"):
                    values[k] = row.get(k, "")
                w.writerow([values.get(k, "") for k in FIELDS])


# ---------------------------------------------------------------------------------------------- the front and its knee
def nondominated(F) -> np.ndarray:
    """A boolean mask of the rows of F (minimised objectives) that no other row beats: at least as good in every objective and
    better in one. O(n^2), plenty for a few thousand rows."""
    F = np.asarray(F, float)
    keep = np.ones(len(F), bool)
    for i in range(len(F)):
        if keep[i]:
            beaten = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
            if beaten.any():
                keep[i] = False
    return keep


def knee(F, lo=None, hi=None) -> int:
    """The index of the knee of a front F (minimised objectives): the point nearest the ideal corner once every objective is
    scaled to 0..1 between lo and hi (by default the front's own best and worst). Pass the same lo and hi to compare two fronts."""
    F = np.asarray(F, float)
    lo = F.min(axis=0) if lo is None else np.asarray(lo, float)
    hi = F.max(axis=0) if hi is None else np.asarray(hi, float)
    span = np.where(hi - lo > 0, hi - lo, 1.0)
    return int(np.argmin(np.linalg.norm((F - lo) / span, axis=1)))
