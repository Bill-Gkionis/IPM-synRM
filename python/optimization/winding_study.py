"""optimization.winding_study: the winding optimizer. Choose the turns, the parallel paths, the strands in hand and the wire of the optimized motor by
ENUMERATING every winding the constraints allow, with every loss modelled (copper_loss.py), the continuous rating from the thermal network (thermal.py) and the voltage from
the flux maps (maps.py), and keep the exact front of three objectives.

What it is, and why
    The first winding study (optimization/winding.py) chose a winding on a 20 mm stack with a DC loss and a one-line AC formula, and found that the loss hardly moves along
    the front while the turns trade the inverter's current against the speed range. Since then the motor has gained what that study lacked: the copper's AC losses with the
    proximity effect, the PWM ripple and the circulating currents between strands (copper_loss.py); a lumped thermal network with the chosen cooling, a fan on the shaft
    (thermal.py); and the flux maps with the cross-saturation and the iron loss (maps.py, iron_loss.py). This study puts them behind the same question: which winding, and
    what does it imply for the inverter and the speed range?

The variables and constraints
    N     turns a coil side, 1 ... 60.                          a     parallel paths, 1, 2 or 4 (the winding's four-fold symmetry: winding.allowed_paths).
    d     the wire: the IEC 60317 R20 sizes of winding.catalogue(), 0.200 ... 1.000 mm (1.12 mm and up do not pass the slot opening, below).
    n_s   strands in hand, 1 ... n_s_cap (the cap is chosen and justified below: 8).
    arrangement   how the strands of a turn lie in the slot: 'cluster' (a compact bundle, the hand-wound default), 'row' (side by side at one depth) or 'transposed' (every
                  strand in every position: an ideal, Litz-like) are the options; 'random' and 'line' (the strands of a turn spread over the slot's depth) are the cases to
                  avoid. The arrangement changes only the circulating-current loss (copper_loss.py).
    Constraints: the copper fill N n_s (pi d^2 / 4) / A_layer <= 0.45 (A_layer = 36.218 mm^2, one layer of the optimized slot; random-wound round wire); the strand's overall
    diameter + 2 x 0.2 mm <= the slot opening (1.5 mm, so d_over <= 1.1 mm); the thermal limits at the continuous rating (hot spot <= 155 C, class F; magnets <= 80 C, N45;
    40 C ambient; the network of thermal.py with the copper loss of copper_loss.py in it); and the voltage: at the continuous torque the corner speed on the bus must be at
    least 1.1 x 3,000 = 3,300 rpm (the first study's rule: the continuous torque holds to 3,000 rpm with 10 % of the voltage to spare).
    The design decisions behind the headline case (the Case class holds them): 60 mm stack (3 rows; 20, 40 and 80 mm besides), a 48 V bus (60 V besides), the shaft fan of the
    thermal study held to 50 W at the top speed of 16,330 rpm (sizing (b); (a), 6 m/s at 3,000 rpm, besides), N45 magnets at 80 C (150 C besides), class F, the end turn of 48.55 mm.

The model, one winding in, its numbers out (each step is built from the one before)
    1. THE COMMON CONTINUOUS RATING. For a winding the current is i = AT a / N [A peak] for AT ampere-turns a coil side (the field depends on the ampere-turns only: the maps'
       scaling rule). At AT the MTPA point of the maps gives (i_d, i_q), the torque and the iron loss of teeth, yoke and rotor (iron_loss.py, 0.3 mm laminations);
       copper_loss.losses gives the winding's total copper loss there as a function of the copper temperature (DC in the slots and the end turns, skin, proximity, PWM ripple,
       circulating currents); the network (thermal.py: the stack, the cooling, 3,000 rpm, the winding's own fill) returns the temperatures, and the margin min(155 C - hot spot,
       80 C - magnets) falls as AT rises. The DC copper loss at a given AT depends on the FILL (the first study: R i^2 = rho l_mt c AT^2 / (k A_layer)), so a winding with
       less copper runs hotter at the same ampere-turns. If every winding were rated at its OWN thermal limit, an under-filled winding would be rated lower: a smaller short-term
       current and, having less torque to give, a higher corner speed, so it would sit on the front by weakness and not by merit (the study first ran that way and put the least-filled windings on its front). So every winding of a case is held to the SAME continuous torque: the reference winding's (the thermal study's example, 26 turns, 4 paths, 8 x 0.315 mm, fill 0.448,
       full copper loss) less an allowance of 5 % (RATING_MARGIN): AT_req = 0.95 AT_ref. The thermal limits become a CONSTRAINT at AT_req: the winding must keep both
       temperatures under their limits with its own fill and its own AC losses. The allowance says how far under the reference a winding's fill and AC losses may put it; the
       report's allowance table gives the sensitivity (0, 3, 5, 8, 12 %). Each winding's own rating (the root of the margin: AT_cont, a secant iteration of 3 to 5 solves) is
       reported as well, with the headroom over AT_req, and the short-term time from a cold start.
    2. THE SHORT-TERM CURRENT, which sizes the inverter: I_st = 2 AT_req a / N (twice the continuous one: the motor map's convention; the maps end at 3 x the rated
       ampere-turns, 355.9, and a larger current is clipped there).
    3. THE CORNER SPEED at the continuous torque. At fixed (i_d, i_q) the terminal voltage |v|^2 = (R i_d - w psi_q)^2 + (R i_q + w psi_d)^2 is a quadratic in the
       electrical speed w, so the speed at which it reaches V_max = V_dc / sqrt 3 is its positive root, in closed form (the first study's equation, on the new maps with their
       cross-saturation). R is the winding's DC resistance at 100 C with the diamond end turn. On 48 V and on 60 V.
    4. THE HIGH-SPEED POINT, where the AC losses differ. All windings are asked for the SAME shaft-torque demand at the same speed: T_hs = 0.3 of the continuous torque, at 10,000
       rpm (constant power from 3,000 rpm: the torque falls as 1 / speed, so 3,000 / 10,000 = 0.3), on the bus, with the smallest current that holds the voltage (MTPA below the corner,
       field weakening above it; the current may not exceed the winding's I_st). In the placeholder's frame (maps.py: psi = s k psi_0(i_0), i = i_0 / s, T = k T_0(i_0), with
       s = N / (10 a) and k = stack / 20 mm) the constant-torque curve is the same for every winding and the voltage limit is a closed form along it, so one polar grid serves all
       the windings of a case (HsCurve; its answer is within 0.2 % of a search four times finer, and the motor-map engine's own solver agrees to three decimals). The objective is
       the winding's total copper loss at that point, copper at 100 C. A winding that cannot hold the voltage at that torque within its I_st cannot give the demand and is out.
       The point is a DEMAND, not a thermal rating: the report says which of the windings the network could carry at it for
       ever (those of many turns could not: the recommended winding would reach 97.5 C at the hot spot and 90.5 C at the magnets; those of few turns could).
       The machine's field of copper_loss.py was sampled for current angles up to 70 degrees; the point lies at 6 ... 81 degrees, so for the windings of many turns (N / a of 1.5
       or more) the field is clipped to 70 (the 'extrapolated' flag of copper_loss.losses): the proximity part is approximate there, and it is 0.7 ... 2.2 % of the loss.
    5. THE OBJECTIVES (all three are a winding's own): f1 = I_st, the peak phase current at the short-term torque (A; it sizes the inverter), minimised; f2 = the corner speed
       on the bus at the continuous torque (rpm; it sets the speed range), maximised; f3 = the copper loss at the high-speed point (W), minimised. Three, because the turns trade
       f1 against f2 along one line (both go as a / N), and f3 is what the AC terms, the resistance and the field weakening add to that trade: the copper loss
       above the corner, where they matter.

The method, and why it is not an optimizer
    The first study enumerated every combination exactly (13 ms for 102,600) and graded NSGA-II and NSGA-III against it: with a tenth of the evaluations they found 90 ... 95 % of the
    hypervolume but the exact optimum of only about a third of the front's points. Here one evaluation costs about 2.5 ms (the thermal solve is most of it) and the combinations
    are 20,316 for 30 strands in hand at most (11,418 for 8): after the prune below the whole case takes 5 s, so the EXACT front is read off, as before, and verified against
    all the designs (winding.verify_front). A search would have to be graded against the enumeration anyway, and would cost the exactness for nothing. The prune: the DC rating
    of a thin-wire winding of each fill (a table of 13 fills, the secant's starting point) is an UPPER bound for any winding of that fill (AC losses only lower a rating), so
    the windings whose fill is below the floor that carries AT_req (0.392 at the headline case) cannot satisfy the thermal constraint and are not solved: 1,965 of the 20,316 are.
    The pool (--workers n) runs chunks of the enumeration in worker processes at below-normal priority under the solver slot; one process is as fast as it needs to be.

The cap on strands in hand
    Hand winding feeds all the strands of a turn from separate spools at once, through a 1.5 mm slot opening, with a needle: each strand more is another spool, another tension to
    keep even and another chance for a crossing that spreads the bundle (and the circulating-current loss with it). The first study's practical variant took 6 (its enumeration
    took the strands as ideal, up to 30). The cap is chosen from the data: cap_table() gives, for caps of 4, 6, 8, 12 and 30, the windings that satisfy every constraint, the
    smallest short-term current among them and what the cap costs in copper loss at the high-speed point. At 60 mm a cap of 6 or less cannot reach the smallest current (81.5 A, 7
    strands of 0.56 mm: the winding at the voltage rule's knee); 8 can, at a mean cost of 2.6 % in the high-speed loss against no cap (9 strands of 0.5 mm would save 2.2 % at
    that current: not worth a ninth spool); 12 leaves a mean 1.7 %, 30 none. The circulating-current model of the 'cluster' bundle was measured at 2, 4, 8, 13, 16, 26, 52
    strands (and interpolated in ln n_s between them), so a cap of 8 also keeps every winding inside its measured range. The same cap reaches the smallest current at 20, 40 and 80 mm too (80 mm needs it: 7 strands of 0.63 mm).
    8 is a choice, not a result: the cap is a field of Case, and the report gives the cost of the others.

Run (from the Python to FEMM folder; no FEMM; the pool takes the solver slot)
    uv run python -m optimization.winding_study --stage headline          the front, the shortlist and the recommended winding at 60 mm, 48 V, the shaft fan (7 s)
    uv run python -m optimization.winding_study --stage variants          the stacks, the magnet grade, the cooling, the arrangements, the cap, the allowance
    uv run python -m optimization.winding_study                           everything and the figures (about 95 s)
    --workers n    a pool of n worker processes     --no-figures    no PNGs
    uv run pytest tests/test_winding_study.py -q                          37 tests, 8 s
Writes out\\optimization\\winding_study\\: report.txt (everything the run says), front_60mm.csv, shortlist_<stack>mm.csv, results.json, fig_ws_*.png (winding_study_figures.py).

Numbers to expect (60 mm, 48 V, the shaft fan held to 50 W, N45 80 C, class F, compact bundles, at most 8 strands)
    The reference winding carries 96.57 A-turns a coil side continuously (4.374 N.m); the common rating is 91.74 A-turns, 4.145 N.m; the high-speed demand 1.244 N.m at 10,000 rpm.
    20,316 windings enumerated (11,418 within 8 strands); the network solved 1,965; the thermal limits keep 1,627 (654 within 8 strands); the voltage rule and the demand leave
    22; the exact front has 8 designs (N / a: I_st [A] / corner speed on 48 V [rpm] / copper loss at 10,000 rpm [W]):
        9 / 4  7 x 0.56    81.5 / 3,464 / 85.9        8 / 4  4 x 0.8    91.7 / 3,912 / 70.0         7 / 4  7 x 0.63   104.8 / 4,479 / 57.7        6 / 4  4 x 0.9   122.3 / 5,241 / 40.9
        5 / 4  8 x 0.71   146.8 / 6,311 / 21.7        5 / 4  5 x 0.9   146.8 / 6,312 / 21.9        4 / 4  8 x 0.8   183.5 / 7,913 / 7.3         3 / 4  8 x 0.9   244.6 / 10,576 / 5.8
    Recommended: 9 turns, 4 paths, 7 x 0.56 mm in a cluster, fill 0.428, R(100 C) 19.26 mOhm; at the common rating 40.8 A continuous and 81.5 A short-term, corner speed
    3,464 rpm (48 V) and 4,353 rpm (60 V), 85.9 W of copper at the high-speed point; its own rating is 94.6 A-turns = 42.0 A (29.7 A rms), 4.28 N.m, 84.8 C at the hot spot and 80.0 C
    at the magnets. The inverter classes of 100, 150, 200 and 300 A take 8/4 4 x 0.8, 5/4 8 x 0.71, 4/4 8 x 0.8 and 3/4 8 x 0.9; no winding satisfies the voltage rule below 81.5 A.
    What the variants change: the bus does not change the rating (60 V: 11 / 4, 6 x 0.56 at 66.7 A short-term); 150 C magnets or the larger fan raise it to about 150 ... 155
    A-turns and the winding becomes 8 / 4, 4 x 0.8 at 150 ... 155 A short-term; random or line strands add 17 ... 19 W of circulating-current loss at 10,000 rpm and take 1 % off the rating.

Limits (what the numbers do not say)
    The corner speeds use the maps' two-dimensional flux linkages, which have no end-winding leakage (the first study's caveat: they are optimistic) and R at 100 C (the winding
    sits near 84 C at its rating: the corner speed of the recommended winding would be 0.1 % higher). The circulating-current loss of copper_loss.py assumes the strands of a turn
    are joined at the ends of every turn; if they were joined only at the coil's terminals and crossed at random, 'random' would cost about N times less (and 'line' the same). The machine's field is
    clipped at 70 degrees of current angle (step 4). The recommendation depends on the continuous torque, and so on the cooling: it moves with the fan's sizing and the magnet
    grade (the variants table). The short-term current is the motor map's convention (twice the continuous), and the time it can be held from cold (about 25 minutes at the rating of
    the shortlisted windings) is the network's, not a test of the winding. Neither the inverter's dead time nor the bus ripple is modelled beyond the 10 % of the voltage rule.

Self-check (the answers follow).
  1. Two windings have the same fill and different N / a. Which of f1, f2, f3 do they share?
  2. Why does a higher bus voltage not change the continuous rating, and what does it change?
  3. The strands of a turn are connected in parallel and spread over the slot's depth. Which objective notices, and why do the compact bundles not?
  4. Why are all the windings of a case held to one continuous torque, and not each to its own thermal limit?
Answers.
  1. None exactly, but the continuous ampere-turns, and so the torque, are the same (the common rating; the DC loss depends on the fill). f1 goes as a / N, f2 as roughly a / N too
     (the voltage goes with the turns), and f3 with the field-weakening current the voltage forces, which grows with N / a.
  2. The continuous rating is a heat balance: the bus enters only through the PWM ripple loss (a fraction of a watt). The bus sets the voltage limit: the corner speed and the
     high-speed point (60 V takes the recommendation from 81.5 A to 66.7 A short-term).
  3. f3 through the copper loss (the circulating current is loss: 17 W at 10,000 rpm for the recommended winding in a random layout, a sixth of the copper loss it would then have); a compact bundle sees almost the same
     vector potential on all its strands, so no EMF difference drives a current between them.
  4. Because a winding with less copper then has less torque to give, and looks better on the other objectives (less current, more speed) for it. Held to one torque, the objectives
     are the cost of the same task, and the thermal limit is what it should be: a constraint.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import warnings
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np

import paths
from optimization import copper_loss as cl
from optimization import maps, thermal, winding
from optimization.iron_loss import iron_loss

# ---------------------------------------------------------------------------------------------- the fixed choices (each one is a ruling or an assumption)
FILL_MAX = 0.45                          # copper fill of a coil layer, bare copper area over the layer's area (random-wound round wire)
CLEARANCE_MM = 0.2                       # room left on EACH side of the strand's overall diameter in the slot opening (the needle)
SLOT_OPEN_MM = 1.5                       # the optimized design's slot opening [mm] (final.json: slot_open)
N_RANGE = (1, 60)                        # turns a coil side
PATHS = (1, 2, 4)                        # parallel paths (the winding's symmetry)
NS_CAP = 8                               # strands in hand at most (module docstring)
RATED_RPM = 3000.0                       # the speed of the continuous rating
HS_RPM = 10000.0                         # the high-speed point
TOP_RPM = 16330.0                        # the rotor's safe speed
RULE_RPM = 1.1 * RATED_RPM               # the corner speed the continuous torque must reach: 10 % of the voltage to spare at 3,000 rpm
ST_FACTOR = 2.0                          # the short-term current over the continuous one
RATING_MARGIN = 0.05                     # the allowance below the reference winding's thermal rating that a winding's own fill and AC losses may cost
HS_FRACTION = 0.3                        # the torque demand of the high-speed point over the continuous torque at 3,000 rpm: constant power, 3,000 / 10,000
T_COPPER = 100.0                         # deg C: the copper in the voltage equation and in the objectives' loss
F_SW = 20000.0                           # Hz: the inverter's switching frequency
V_DC, V_DC_ALT = 48.0, 60.0              # the bus [V]
PP = 4                                   # pole pairs
SQRT3 = math.sqrt(3.0)
END_TURN_MM, END_OVERHANG_MM = 48.55, 13.6
LAMINATION_MM = 0.3
OUT = lambda: paths.path("out") / "optimization" / "winding_study"


@dataclass(frozen=True)
class Case:
    """One case of the study: everything that is not a winding. Hashable (the contexts are cached by it)."""
    stack_mm: float = 60.0
    cooling: str = "fan_b"               # 'fan_b' (the shaft fan held to fan_budget_w at the top speed), 'fan_a' (6 m/s at 3,000 rpm), 'forced', 'natural', 'liquid'
    fan_budget_w: float = 50.0
    magnet_limit_c: float = thermal.MAGNET_LIMIT_C
    wind_limit_c: float = thermal.CLASS_F_C
    arrangement: str = "cluster"
    n_s_cap: int = NS_CAP
    v_dc: float = V_DC
    f_sw: float = F_SW
    hs_rpm: float = HS_RPM
    hs_fraction: float = HS_FRACTION
    rating_margin: float = RATING_MARGIN

    def label(self) -> str:
        return (f"{self.stack_mm:g} mm, {self.cooling}{'' if self.cooling != 'fan_b' else f' {self.fan_budget_w:g} W'}, magnets {self.magnet_limit_c:g} C, "
                f"{self.arrangement}, <= {self.n_s_cap} strands, {self.v_dc:g} V")


def assumptions() -> thermal.Assumptions:
    """The network's assumptions with the copper-loss model's diamond end turn (a design decision: one end turn in both models)."""
    return thermal.Assumptions(end_turn_mm=END_TURN_MM, end_overhang_mm=END_OVERHANG_MM)


def make_cooling(case: Case) -> thermal.Cooling:
    """The cooling of a case. 'fan_b' is the shaft fan whose power at 16,330 rpm is held to case.fan_budget_w (the thermal study's sizing (b), 50 W); 'fan_a' the one
    sized for the rated speed (6 m/s in the fin channels at 3,000 rpm); the others are thermal.py's presets."""
    if case.cooling == "fan_b":
        return thermal.size_shaft_fan(case.fan_budget_w, stack_mm=case.stack_mm, assume=assumptions())["cooling"]
    if case.cooling == "fan_a":
        return thermal.Cooling.shaft_fan(6.0, RATED_RPM)
    if case.cooling in ("forced", "forced_air"):
        return thermal.Cooling.forced_air()
    if case.cooling == "natural":
        return thermal.Cooling.natural()
    if case.cooling == "liquid":
        return thermal.Cooling.liquid()
    raise ValueError(f"cooling {case.cooling!r}: fan_b, fan_a, forced, natural or liquid")


# ---------------------------------------------------------------------------------------------- the candidates
def wires() -> dict:
    """The wire catalogue's sizes that pass the slot opening: {"d_mm", "d_over_mm"} (arrays), with d_over + 2 x 0.2 mm <= 1.5 mm."""
    cat = winding.catalogue()
    keep = cat["d_over_mm"] + 2.0 * CLEARANCE_MM <= SLOT_OPEN_MM + 1e-9
    return {"d_mm": np.asarray(cat["d_mm"])[keep], "d_over_mm": np.asarray(cat["d_over_mm"])[keep]}


def layer_area_mm2() -> float:
    """One coil layer of the optimized slot [mm^2] (copper_loss's tables: 36.218)."""
    return float(cl.tables()["geometry"]["A_layer_mm2"])


def candidates(n_s_cap: int = NS_CAP, fill_max: float = FILL_MAX) -> dict:
    """Every allowed winding as arrays (one entry each): N, a, n_s, wire index di, d_mm, d_over_mm, fill. Allowed: N 1 ... 60, a in (1, 2, 4), n_s 1 ... n_s_cap, the wires
    that pass the opening, fill <= fill_max."""
    w = wires()
    n, di, ns, a = np.meshgrid(np.arange(N_RANGE[0], N_RANGE[1] + 1), np.arange(len(w["d_mm"])), np.arange(1, n_s_cap + 1), np.asarray(PATHS), indexing="ij")
    n, di, ns, a = n.ravel(), di.ravel(), ns.ravel(), a.ravel()
    fill = n * ns * math.pi * w["d_mm"][di] ** 2 / 4.0 / layer_area_mm2()
    ok = fill <= fill_max + 1e-12
    n, di, ns, a, fill = n[ok], di[ok], ns[ok], a[ok], fill[ok]
    return {"N": n.astype(int), "a": a.astype(int), "n_s": ns.astype(int), "di": di.astype(int), "d_mm": w["d_mm"][di], "d_over_mm": w["d_over_mm"][di], "fill": fill}


def winding_dict(c: dict, i: int, arrangement: str | None = None) -> dict:
    """Candidate i of candidates() as copper_loss.losses' winding dict."""
    w = {"N": int(c["N"][i]), "a": int(c["a"][i]), "n_s": int(c["n_s"][i]), "d_mm": float(c["d_mm"][i]), "d_over_mm": float(c["d_over_mm"][i])}
    if arrangement:
        w["arrangement"] = arrangement
    return w


# ---------------------------------------------------------------------------------------------- the machine's own models, built once per process
_CONTEXTS: dict = {}


class Context:
    """What one case needs again and again, built once per process: the flux maps, the network of the case's stack, cooling and the rated speed, the table of the DC rating
    against the fill (the secant's starting point), and the constants."""

    def __init__(self, case: Case):
        self.case = case
        self.fm = maps.FluxMaps.load()
        self.assume = assumptions()
        self.cooling = make_cooling(case)
        self.tm = thermal.ThermalModel(None, case.stack_mm, self.cooling, RATED_RPM, thermal.WindingSpec(0.45), self.assume)
        self.stack = case.stack_mm
        self.k = case.stack_mm / self.fm.depth_mm
        self.at_max = 3.0 * 10.0 * self.fm.rated                          # the maps' range: 3 x the rated ampere-turns a coil side (355.9)
        self.a_layer = layer_area_mm2()
        self.rho100 = float(cl.resistivity(T_COPPER))
        self.l_mt = 2.0 * case.stack_mm * 1e-3 + 2.0 * END_TURN_MM * 1e-3          # the mean turn [m]: the copper-loss model's diamond end turn
        self._guess = None
        self._hs: dict = {}

    # ---- the DC rating against the fill: the rating of a thin-wire (so, nearly AC-free) winding of each fill, in log-log; the root-finder's starting point and an UPPER bound for any winding
    def _table(self):
        if self._guess is None:
            fills = np.array([0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.38, 0.4, 0.42, 0.44, 0.45])
            rows = []
            for f in fills:
                n_s = max(1, int(round(f * self.a_layer / (10 * math.pi * 0.2 ** 2 / 4.0))))
                w = {"N": 10, "a": 1, "n_s": n_s, "d_mm": math.sqrt(4.0 * f * self.a_layer / (10 * n_s * math.pi)), "d_over_mm": 0.239}     # a mathematical wire of exactly this fill
                r = continuous_rating(self, w, "transposed", 100.0, exact_fill=f)
                rows.append(r["AT"])
            self._guess = (np.log(fills), np.log(np.maximum(rows, 1e-3)))
        return self._guess

    def at_guess(self, fill: float) -> float:
        lf, la = self._table()
        return float(np.exp(np.interp(math.log(max(fill, 1e-4)), lf, la)))

    def fill_floor(self, at: float, tolerance: float = 0.01) -> float:
        """The smallest fill whose DC rating reaches (1 - tolerance) x at: no winding of a smaller fill can carry `at` (AC losses and the end turns' share only lower a rating below
        its DC value, the arrangement 'transposed' having none), so it need not be solved. 1.0 if no fill in the table can."""
        lf, la = self._table()
        target = math.log((1.0 - tolerance) * at)
        if target <= la[0]:
            return 0.0
        if target > la[-1]:
            return 1.0
        return float(np.exp(np.interp(target, la, lf)))


def get_context(case: Case) -> Context:
    key = case
    if key not in _CONTEXTS:
        if len(_CONTEXTS) > 12:
            _CONTEXTS.clear()
        _CONTEXTS[key] = Context(case)
    return _CONTEXTS[key]


# ---------------------------------------------------------------------------------------------- the continuous rating of one winding
def copper_at(ctx: Context, w: dict, arrangement: str, at: float, rpm: float, gamma: float, temps=(40.0, 100.0, 160.0), v_dc: float | None = None) -> dict:
    """copper_loss.losses_arrays of the winding at AT ampere-turns a coil side (peak), the current angle gamma and the speed, for an ARRAY of copper temperatures in one call."""
    n, a = w["N"], w["a"]
    return cl.losses_arrays(n, a, w["n_s"], w["d_mm"], w["d_over_mm"], at * a / n, gamma, rpm, ctx.case.f_sw, ctx.case.v_dc if v_dc is None else v_dc,
                            np.asarray(temps, float), ctx.stack, arrangement=arrangement)


def margin_at(ctx: Context, w: dict, arrangement: str, fill: float, at: float, t_guess=None) -> tuple[float, "thermal.SteadyResult", dict]:
    """The thermal margin min(wind limit - hot spot, magnet limit - magnets) [K] at AT ampere-turns, the network's SteadyResult and the operating point."""
    case = ctx.case
    n, a = w["N"], w["a"]
    at = min(at, ctx.at_max)
    gamma, i_d, i_q, torque = ctx.fm.mtpa(at * a / n, turns=n, paths=a, stack_mm=ctx.stack)
    iron = iron_loss(i_d, i_q, RATED_RPM, ctx.stack, turns=n, paths=a, lamination_mm=LAMINATION_MM)
    r = copper_at(ctx, w, arrangement, at, RATED_RPM, float(gamma))
    temps = np.array([40.0, 100.0, 160.0])
    coef = np.polyfit(temps, np.asarray(r["p_total"], float), 2)
    p100 = float(r["p_total"][1])
    dc100 = float(r["p_dc_slot"][1] + r["p_dc_end"][1])
    # the share of the copper loss in the slots: the end turns carry their DC loss and its skin share; proximity, the PWM ripple's and the circulating loss are in the slots
    end = float(r["p_dc_end"][1]) * (1.0 + float(r["p_skin"][1]) / dc100) if dc100 > 0 else 0.0
    share = 1.0 - end / p100 if p100 > 0 else None
    res = ctx.tm.steady_state(thermal.Losses(copper=lambda t: float(np.polyval(coef, t)), copper_slot_fraction=share, iron_teeth=iron["teeth"], iron_yoke=iron["yoke"],
                                             iron_rotor=iron["rotor"]), fill=fill, wind_limit_c=case.wind_limit_c, magnet_limit_c=case.magnet_limit_c, t_guess=t_guess)
    m = -1.0e3 if (res.runaway or not res.converged) else min(res.margin_winding, res.margin_magnet)
    return m, res, {"gamma": float(gamma), "i_d": float(i_d), "i_q": float(i_q), "torque": float(torque), "iron": iron, "share": share}


def continuous_rating(ctx: Context, w: dict, arrangement: str, guess: float | None = None, tol_k: float = 2e-3, exact_fill: float | None = None) -> dict:
    """The continuous ampere-turns a coil side (peak) of the winding at the rated speed: the root of the thermal margin (module docstring, step 1), by a secant iteration.
    Returns {"AT", "feasible" (False: a limit is passed with no current at all), "limited_by_maps" (the limit lies beyond the maps' 355.9), "binding", "hot_spot",
    "magnet", "iterations", "result" (the network's), "point" (the MTPA point)}."""
    n, a = w["N"], w["a"]
    fill = exact_fill if exact_fill is not None else n * w["n_s"] * math.pi * w["d_mm"] ** 2 / 4.0 / ctx.a_layer
    x0 = guess if guess is not None else (ctx.at_guess(fill) if ctx._guess is not None or exact_fill is None else 100.0)
    x0 = min(max(x0, 0.05), ctx.at_max)
    t_guess = [None]
    calls = [0]
    RUNAWAY = -900.0                                                        # a margin at or below this is the network's thermal runaway: no slope to trust

    def f(x: float):
        calls[0] += 1
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            m, res, pt = margin_at(ctx, w, arrangement, fill, x, t_guess[0])
        if not res.runaway and np.isfinite(m):
            t_guess[0] = np.array([res.T[k] for k in thermal.NODES])
        return m, res, pt
    # ---- 1. a bracket: margin > 0 at `lo` (cool enough), < 0 at `hi` (too hot)
    x, (fx, rx, px) = x0, f(x0)
    pts = {x0: (fx, rx, px)}
    if fx > 0:
        lo, hi, step = x0, None, 1.0 + 0.04
        while hi is None:
            if lo >= ctx.at_max * (1 - 1e-12):
                return _rating(ctx.at_max, pts[lo][1], pts[lo][2], True, True, calls[0])
            xn = min(lo * step, ctx.at_max)
            pts[xn] = f(xn)
            if pts[xn][0] > 0:
                lo, step = xn, step * 1.6
            else:
                hi = xn
    else:
        hi, lo, step = x0, None, 1.0 - 0.04
        while lo is None:
            if hi < 0.02:
                r2 = pts[hi]
                return {"AT": 0.0, "feasible": False, "limited_by_maps": False, "binding": r2[1].binding, "hot_spot": float(r2[1].hot_spot), "magnet": float(r2[1].magnet),
                        "iterations": calls[0], "result": r2[1], "point": r2[2]}
            xn = hi * step
            pts[xn] = f(xn)
            if pts[xn][0] <= 0:
                hi, step = xn, step * 0.6
            else:
                lo = xn
    # ---- 2. the root: Illinois (regula falsi) with a bisection step wherever the high side is a runaway
    f_lo, f_hi = pts[lo][0], pts[hi][0]
    side = 0
    best = (lo, pts[lo])
    for _ in range(40):
        if abs(pts[lo][0]) < tol_k or hi - lo < 1e-7 * hi:
            break
        if f_hi <= RUNAWAY:
            xn = 0.5 * (lo + hi)
        else:
            xn = (lo * f_hi - hi * f_lo) / (f_hi - f_lo)
            if not (lo < xn < hi):
                xn = 0.5 * (lo + hi)
        fn = f(xn)
        pts[xn] = fn
        if fn[0] > 0:
            lo, f_lo = xn, fn[0]
            if side == 1 and f_hi > RUNAWAY:
                f_hi *= 0.5
            side = 1
        else:
            hi, f_hi = xn, fn[0]
            if side == -1:
                f_lo *= 0.5
            side = -1
        best = (lo, pts[lo])
        if abs(fn[0]) < tol_k:
            best = (xn, fn)
            break
    xb, (fb, rb, pb) = best
    return _rating(xb, rb, pb, False, True, calls[0])


def _rating(x, res, pt, capped, feasible, calls) -> dict:
    return {"AT": float(x), "feasible": feasible, "limited_by_maps": bool(capped), "binding": res.binding, "hot_spot": float(res.hot_spot), "magnet": float(res.magnet),
            "iterations": calls, "result": res, "point": pt}


# ---------------------------------------------------------------------------------------------- the closed forms (arrays over windings)
def series_resistance(ctx: Context, n, a, n_s, d_mm, t_c: float = T_COPPER):
    """The phase's DC resistance [Ohm] at t_c with the diamond end turn: rho l_mt (12 N / a) / (a n_s pi d^2 / 4) (the copper-loss model's, in closed form)."""
    rho = float(cl.resistivity(t_c))
    return rho * ctx.l_mt * 12.0 * np.asarray(n, float) / np.asarray(a, float) / (np.asarray(a, float) * np.asarray(n_s, float) * math.pi * (np.asarray(d_mm, float) * 1e-3) ** 2 / 4.0)


def corner_speed(psi_d, psi_q, i_d, i_q, r, v_max):
    """The mechanical speed [rpm] at which the terminal voltage of a winding at the fixed current (i_d, i_q) reaches v_max (arrays): the positive root of
    (psi_d^2 + psi_q^2) w^2 + 2 R (i_q psi_d - i_d psi_q) w + R^2 (i_d^2 + i_q^2) - v_max^2 = 0; 0 where the resistive drop alone is over the limit."""
    qa = psi_d ** 2 + psi_q ** 2
    qb = 2.0 * r * (i_q * psi_d - i_d * psi_q)
    qc = r ** 2 * (i_d ** 2 + i_q ** 2) - v_max ** 2
    ok = qc < 0.0
    w = np.where(ok, (-qb + np.sqrt(np.where(ok, qb ** 2 - 4.0 * qa * qc, 0.0))) / (2.0 * qa), 0.0)
    return w * 60.0 / (2.0 * math.pi * PP)


def voltage(psi_d, psi_q, i_d, i_q, r, rpm):
    """The peak phase voltage [V] at the current (i_d, i_q) and the speed (arrays)."""
    w = PP * 2.0 * math.pi * rpm / 60.0
    return np.hypot(r * i_d - w * psi_q, r * i_q + w * psi_d)


class HsCurve:
    """The constant-torque curve of the high-speed demand in the placeholder's frame (module docstring, step 4): for each current angle gamma the smallest placeholder current
    i_0 whose torque k T_0(i_0) is the demand (the first crossing along the ray), with the flux linkages there. A winding (s, R) then asks, along the curve, for the smallest
    current at which |v| <= V_max and |i| <= its limit: solve()."""

    def __init__(self, fm: "maps.FluxMaps", stack_mm: float, torque_nm: float, n_gamma: int = 721, n_radius: int = 480):
        self.k = stack_mm / fm.depth_mm
        self.target = torque_nm
        gam = np.linspace(0.0, 90.0, n_gamma)
        i0 = np.linspace(0.0, min(-fm.range_d[0], fm.range_q), n_radius)
        g = np.radians(gam)[:, None]
        d0, q0 = -i0[None, :] * np.sin(g), i0[None, :] * np.cos(g)
        pd, pq = fm.psi_dq(d0, q0)
        t0 = 1.5 * PP * (pd * q0 - pq * d0)
        want = torque_nm / self.k
        reach = (t0 >= want).any(axis=1)
        first = np.argmax(t0 >= want, axis=1)
        rows = np.flatnonzero(reach & (first > 0))
        f = (want - t0[rows, first[rows] - 1]) / (t0[rows, first[rows]] - t0[rows, first[rows] - 1])
        i_c = i0[first[rows] - 1] + f * (i0[first[rows]] - i0[first[rows] - 1])
        th = np.radians(gam[rows])
        self.gamma = gam[rows]
        self.i0 = i_c
        self.d0, self.q0 = -i_c * np.sin(th), i_c * np.cos(th)
        self.psi_d0, self.psi_q0 = fm.psi_dq(self.d0, self.q0)

    def solve(self, s, r, rpm: float, v_max: float, i_lim, chunk: int = 4000) -> dict:
        """For arrays s (N / (10 a)), r (the winding's R at the operating temperature) and i_lim (A peak): the smallest feasible current. Returns arrays I [A], gamma [deg],
        feasible (bool); infeasible entries carry NaN."""
        s, r, i_lim = np.broadcast_arrays(np.asarray(s, float), np.asarray(r, float), np.asarray(i_lim, float))
        w_e = PP * 2.0 * math.pi * rpm / 60.0
        n = len(s)
        out_i, out_g = np.full(n, np.nan), np.full(n, np.nan)
        if len(self.gamma) == 0:
            return {"I": out_i, "gamma": out_g, "feasible": np.zeros(n, bool)}
        for lo in range(0, n, chunk):
            sl = slice(lo, min(lo + chunk, n))
            ss, rr, ll = s[sl][:, None], r[sl][:, None], i_lim[sl][:, None]
            i_d, i_q = self.d0[None, :] / ss, self.q0[None, :] / ss
            p_d, p_q = ss * self.k * self.psi_d0[None, :], ss * self.k * self.psi_q0[None, :]
            v = np.hypot(rr * i_d - w_e * p_q, rr * i_q + w_e * p_d)
            cur = self.i0[None, :] / ss
            ok = (v <= v_max) & (cur <= ll)
            has = ok.any(axis=1)
            cur_ok = np.where(ok, cur, np.inf)
            j = np.argmin(cur_ok, axis=1)
            rows = np.arange(ok.shape[0])
            ii, gg = cur_ok[rows, j], self.gamma[j]
            # refine: the feasible point of smallest current lies at the voltage boundary when the neighbour towards the smaller current is infeasible
            for step in (-1, 1):
                jn = np.clip(j + step, 0, len(self.gamma) - 1)
                bad = has & ~ok[rows, jn] & (cur[rows, jn] < cur[rows, j]) & (jn != j) & (cur[rows, jn] <= ll[:, 0])
                dv_j, dv_n = v[rows, j] - v_max, v[rows, jn] - v_max
                with np.errstate(divide="ignore", invalid="ignore"):
                    f = np.where(bad & (dv_n != dv_j), (0.0 - dv_j) / (dv_n - dv_j), 0.0)
                f = np.clip(f, 0.0, 1.0)
                ii = np.where(bad, cur[rows, j] + f * (cur[rows, jn] - cur[rows, j]), ii)
                gg = np.where(bad, self.gamma[j] + f * (self.gamma[jn] - self.gamma[j]), gg)
            out_i[sl], out_g[sl] = np.where(has, ii, np.nan), np.where(has, gg, np.nan)
        return {"I": out_i, "gamma": out_g, "feasible": np.isfinite(out_i)}


REFERENCE = {"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367}          # the thermal study's example winding (fill 0.4476): the rating's reference


def requirement(ctx: Context) -> dict:
    """The common continuous rating every winding is held to (module docstring, step 1): the reference winding's own thermal rating at the rated speed (the full losses,
    ctx.case.arrangement), less the case's allowance for a winding's own fill and AC losses: AT_req = (1 - rating_margin) AT_ref. Returns {"AT_ref", "T_ref", "AT" (AT_req),
    "T" (the MTPA torque at AT_req), "gamma", "limited_by_maps", "reference"}. Cached in the context."""
    if "req" not in ctx._hs:
        r = continuous_rating(ctx, REFERENCE, ctx.case.arrangement)
        at_req = (1.0 - ctx.case.rating_margin) * r["AT"]
        g, i0d, i0q, t_req = ctx.fm.mtpa(at_req / 10.0, stack_mm=ctx.stack)
        ctx._hs["req"] = {"AT_ref": r["AT"], "T_ref": float(r["point"]["torque"]), "AT": float(at_req), "T": float(t_req), "gamma": float(g), "limited_by_maps": r["limited_by_maps"],
                          "reference": dict(REFERENCE)}
    return ctx._hs["req"]


def hs_curve(ctx: Context) -> HsCurve:
    """The constant-torque curve of the high-speed demand: hs_fraction x the torque of the common rating (the same for every winding)."""
    if "curve" not in ctx._hs:
        ctx._hs["curve"] = HsCurve(ctx.fm, ctx.stack, ctx.case.hs_fraction * requirement(ctx)["T"])
    return ctx._hs["curve"]


# ---------------------------------------------------------------------------------------------- the thermal check of every winding at the common rating
def check_rows(task: dict) -> list[dict]:
    """The pool's function: the network at the common rating for each winding of a chunk, with the winding's full copper loss. task = {"case", "arrangement", "rows": [(N, a, n_s,
    d_mm, d_over_mm), ...]}. Returns for each: margin (the smaller of the two margins, K; negative: a limit is passed), margin_winding, margin_magnet, hot_spot, magnet, binding,
    p_cu (the copper loss at the solution), iterations, runaway."""
    case = task["case"]
    ctx = get_context(case)
    at = requirement(ctx)["AT"]
    guess = None
    out = []
    for (n, a, n_s, d, d_over) in task["rows"]:
        w = {"N": int(n), "a": int(a), "n_s": int(n_s), "d_mm": float(d), "d_over_mm": float(d_over)}
        fill = w["N"] * w["n_s"] * math.pi * w["d_mm"] ** 2 / 4.0 / ctx.a_layer
        try:
            with np.errstate(all="ignore"), warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                m, res, _ = margin_at(ctx, w, task["arrangement"], fill, at, guess)
            if not res.runaway and np.isfinite(m):
                guess = np.array([res.T[k] for k in thermal.NODES])
            out.append({"margin": float(m), "margin_winding": float(res.margin_winding), "margin_magnet": float(res.margin_magnet), "hot_spot": float(res.hot_spot),
                        "magnet": float(res.magnet), "binding": res.binding, "p_cu": float(res.copper_loss), "iterations": int(res.iterations), "runaway": bool(res.runaway)})
        except Exception as error:                                                         # a failed candidate must not end a run
            out.append({"margin": -1.0e3, "margin_winding": float("nan"), "margin_magnet": float("nan"), "hot_spot": float("nan"), "magnet": float("nan"),
                        "binding": "error: " + str(error)[:80], "p_cu": float("nan"), "iterations": 0, "runaway": True})
    return out


def run_check(case: Case, cand: dict, arrangement: str | None = None, workers: int = 1, chunk: int = 500) -> dict:
    """check_rows for every candidate: arrays margin, margin_winding, margin_magnet, hot_spot, magnet, binding, p_cu, iterations, runaway, solved. The windings whose fill is below
    the DC rating's floor (Context.fill_floor: no smaller fill can carry the common rating) are not solved: margin -1000, binding 'pruned', solved False. workers > 1: a pool of
    worker processes at below-normal priority under the solver slot (femm_queue, no FEMM)."""
    arr = arrangement or case.arrangement
    ctx = get_context(case)
    floor = ctx.fill_floor(requirement(ctx)["AT"])
    todo = np.flatnonzero(cand["fill"] >= floor)
    n = len(cand["N"])
    rows = list(zip(cand["N"][todo].tolist(), cand["a"][todo].tolist(), cand["n_s"][todo].tolist(), cand["d_mm"][todo].tolist(), cand["d_over_mm"][todo].tolist()))
    tasks = [{"case": case, "arrangement": arr, "rows": rows[i:i + chunk]} for i in range(0, len(rows), chunk)]
    if workers > 1 and len(tasks) > 1:
        import femm_queue
        with femm_queue.femm_pool(f"winding study {case.label()}", workers=workers, task_timeout_s=900, femm=False) as pool:
            parts = pool.map(check_rows, tasks)
    else:
        parts = [check_rows(t) for t in tasks]
    flat = [r for part in parts for r in part]
    if len(flat) != len(rows) or any(not isinstance(r, dict) or "margin" not in r for r in flat):
        raise RuntimeError("a worker did not return its chunk: " + str([r for r in flat if not (isinstance(r, dict) and "margin" in r)][:2]))
    out = {"margin": np.full(n, -1.0e3), "margin_winding": np.full(n, np.nan), "margin_magnet": np.full(n, np.nan), "hot_spot": np.full(n, np.nan), "magnet": np.full(n, np.nan),
           "binding": np.full(n, "pruned", dtype=object), "p_cu": np.full(n, np.nan), "iterations": np.zeros(n, int), "runaway": np.zeros(n, bool), "solved": np.zeros(n, bool)}
    for key in ("margin", "margin_winding", "margin_magnet", "hot_spot", "magnet", "binding", "p_cu", "iterations", "runaway"):
        if flat:
            out[key][todo] = np.array([r[key] for r in flat], dtype=out[key].dtype)
    out["solved"][todo] = True
    out["fill_floor"] = floor
    return out


# ---------------------------------------------------------------------------------------------- the objectives
def objectives(case: Case, cand: dict, check: dict, arrangement: str | None = None) -> dict:
    """The three objectives, the constraints and the quantities behind them, for every candidate (arrays), at the common rating (module docstring, steps 1 to 5). `check` is
    run_check's answer. Keys: thermal_ok, AT, I_cont, I_st, AT_st, T_cont, T_st, R100, corner48, corner60, I_hs, gamma_hs, hs_feasible, P_hs, P_hs_dc, hs_extrapolated, T_hs, fill."""
    ctx = get_context(case)
    arr = arrangement or case.arrangement
    req = requirement(ctx)
    n, a, n_s, d, d_over = (cand[k] for k in ("N", "a", "n_s", "d_mm", "d_over_mm"))
    s = n / (10.0 * a)                                                      # the winding against the placeholder (maps.py): psi = s psi_0(i_0), i = i_0 / s, i_0 = AT / 10
    at = np.full(len(n), req["AT"])
    gam, i0d, i0q, t_cont = ctx.fm.mtpa(at / 10.0, stack_mm=ctx.stack)                       # the MTPA point at the rating in the placeholder's terms (any winding's)
    p0d, p0q = ctx.fm.psi_dq(i0d, i0q, stack_mm=ctx.stack)
    psi_d, psi_q, i_d, i_q = s * p0d, s * p0q, i0d / s, i0q / s
    i_cont = at * a / n
    r100 = series_resistance(ctx, n, a, n_s, d)
    corner_bus = corner_speed(psi_d, psi_q, i_d, i_q, r100, case.v_dc / SQRT3)
    corner48 = corner_speed(psi_d, psi_q, i_d, i_q, r100, V_DC / SQRT3)
    corner60 = corner_speed(psi_d, psi_q, i_d, i_q, r100, V_DC_ALT / SQRT3)
    at_st = np.minimum(ST_FACTOR * at, ctx.at_max)
    i_st = at_st * a / n
    t_st = ctx.fm.mtpa(at_st / 10.0, stack_mm=ctx.stack)[3]
    curve = hs_curve(ctx)
    sol = curve.solve(s, r100, case.hs_rpm, case.v_dc / SQRT3, i_st)
    thermal_ok = (check["margin"] >= 0.0) & ~check["runaway"]
    ok = sol["feasible"]
    i_hs = np.where(ok, sol["I"], 0.0)
    g_hs = np.where(ok, sol["gamma"], 0.0)
    lh = cl.losses_arrays(n, a, n_s, d, d_over, i_hs, g_hs, case.hs_rpm, case.f_sw, case.v_dc, T_COPPER, ctx.stack, arrangement=arr)
    p_hs = np.where(ok, lh["p_total"], np.inf)
    return {"thermal_ok": thermal_ok, "AT": at, "I_cont": i_cont, "I_st": i_st, "AT_st": at_st, "T_cont": t_cont, "T_st": t_st, "gamma_cont": gam, "R100": r100, "corner_bus": corner_bus,
            "corner48": corner48, "corner60": corner60, "I_hs": i_hs, "gamma_hs": g_hs, "hs_feasible": ok, "P_hs": p_hs, "P_hs_dc": np.where(ok, lh["p_dc_slot"] + lh["p_dc_end"], np.inf),
            "hs_extrapolated": np.asarray(lh["extrapolated"]) & ok, "T_hs": curve.target, "fill": cand["fill"]}


def allowed(res: dict, cand: dict, rule: bool = True, cap: int | None = None, hs: bool = True) -> np.ndarray:
    """The mask of the windings that satisfy the constraints: the thermal limits at the common rating, the high-speed demand holdable (hs=True), the voltage rule (rule=True: the
    corner speed on the case's bus at least 3,300 rpm), and at most `cap` strands in hand."""
    ok = res["thermal_ok"].copy()
    if hs:
        ok &= res["hs_feasible"]
    if rule:
        ok &= res["corner_bus"] >= RULE_RPM
    if cap is not None:
        ok &= cand["n_s"] <= cap
    return ok


def front(res: dict, cand: dict, rule: bool = True, cap: int | None = None, hs: bool = True) -> np.ndarray:
    """The indices of the exact front of (I_st minimised, corner speed on 48 V maximised, copper loss at the high-speed point minimised) among the allowed windings (allowed())."""
    idx = np.flatnonzero(allowed(res, cand, rule, cap, hs))
    if not len(idx):
        return idx
    f = np.stack([res["I_st"][idx], -res["corner_bus"][idx], res["P_hs"][idx] if hs else np.zeros(len(idx))], axis=1)
    keep = winding.pareto_mask(winding.round_sig(f, 9))
    return idx[keep]


# ---------------------------------------------------------------------------------------------- one case, run
@dataclass
class CaseRun:
    """A case run over every allowed winding: the candidates, the thermal check at the common rating and the objectives (arrays), and the time it took."""
    case: Case
    cand: dict
    check: dict
    obj: dict
    seconds: float
    arrangement: str

    def allowed(self, rule: bool = True, cap: int | None = None, hs: bool = True) -> np.ndarray:
        return allowed(self.obj, self.cand, rule, self.case.n_s_cap if cap is None else cap, hs)

    def front(self, rule: bool = True, cap: int | None = None, hs: bool = True) -> np.ndarray:
        return front(self.obj, self.cand, rule, self.case.n_s_cap if cap is None else cap, hs)


def run_case(case: Case, workers: int = 1, enum_cap: int = 30, arrangement: str | None = None) -> CaseRun:
    """Enumerate every winding with up to enum_cap strands (so that the cap can be varied afterwards), check each against the thermal limits at the common rating with its full
    copper loss, and compute the objectives."""
    t0 = time.perf_counter()
    arr = arrangement or case.arrangement
    cand = candidates(max(enum_cap, case.n_s_cap))
    chk = run_check(case, cand, arr, workers=workers)
    obj = objectives(case, cand, chk, arr)
    return CaseRun(case, cand, chk, obj, time.perf_counter() - t0, arr)


# ---------------------------------------------------------------------------------------------- one winding, in full
def _copper_poly(ctx: Context, w: dict, arrangement: str, at: float, rpm: float, gamma: float):
    """(f(T), share, r100): the winding's total copper loss as a function of the copper temperature (a quadratic through 40, 100 and 160 C), the share of it in the slots, and the
    module's answer at 100 C, at the ampere-turns `at`, speed and current angle."""
    r = copper_at(ctx, w, arrangement, at, rpm, gamma)
    coef = np.polyfit(np.array([40.0, 100.0, 160.0]), np.asarray(r["p_total"], float), 2)
    p100 = float(r["p_total"][1])
    dc = float(r["p_dc_slot"][1] + r["p_dc_end"][1])
    end = float(r["p_dc_end"][1]) * (1.0 + float(r["p_skin"][1]) / dc) if dc > 0 else 0.0
    share = 1.0 - end / p100 if p100 > 0 else None
    return (lambda t: float(np.polyval(coef, t))), share, {k: float(np.asarray(v)[1]) if np.ndim(v) else float(v) for k, v in r.items() if k not in ("extrapolated", "overmodulated")}


def manufacturing(w: dict) -> dict:
    """What it takes to wind it: the strands in hand, the wire, how many strands go abreast through the opening (with the 0.2 mm of the needle on each side), the room to spare, and
    the stiffness of the wire against the 0.315 mm of the first study's example (a wire's bending stiffness goes as d^4)."""
    room = SLOT_OPEN_MM - 2.0 * CLEARANCE_MM
    return {"n_s": w["n_s"], "d_mm": w["d_mm"], "d_over_mm": w["d_over_mm"], "strands_abreast": int(room // w["d_over_mm"]), "opening_room_mm": room - w["d_over_mm"],
            "stiffness_vs_0315": (w["d_mm"] / 0.315) ** 4, "bundle_stiffness_vs_example": w["n_s"] * (w["d_mm"] / 0.315) ** 4 / 8.0,
            "strands_per_coil_side": w["N"] * w["n_s"], "turns_per_coil": w["N"], "coils_per_path": 12 // w["a"]}


def detail(case: Case, w: dict, arrangement: str | None = None) -> dict:
    """Everything about one winding (a flat dict of numbers): its own continuous rating from the network; the losses at it at the network's copper temperature (DC in the slots and
    the end turns, skin, proximity, PWM, circulating); R; the currents and the corner speeds on 48 and 60 V; the common-rating objectives (what the front used); the high-speed point
    with its losses and the network's temperatures there; and the manufacturing note."""
    ctx = get_context(case)
    arr = arrangement or case.arrangement
    n, a = w["N"], w["a"]
    fill = n * w["n_s"] * math.pi * w["d_mm"] ** 2 / 4.0 / ctx.a_layer
    s = n / (10.0 * a)
    req = requirement(ctx)
    ww = dict(w, arrangement=arr)
    # ---- the winding's own continuous rating and the state of the network there
    rt = continuous_rating(ctx, w, arr)
    at_c, res, pt = rt["AT"], rt["result"], rt["point"]
    t_cu = float(pt["share"] * res.T["W"] + (1.0 - pt["share"]) * res.T["E"]) if pt["share"] is not None else float(res.T["W"])
    i_c = at_c * a / n
    op = {"i_peak": i_c, "gamma_deg": pt["gamma"], "rpm": RATED_RPM, "f_sw": case.f_sw, "v_dc": case.v_dc}
    lc = cl.losses(ww, op, t_cu, ctx.stack)
    r100 = float(series_resistance(ctx, n, a, w["n_s"], w["d_mm"]))
    out = {"N": n, "a": a, "n_s": w["n_s"], "d_mm": w["d_mm"], "d_over_mm": w["d_over_mm"], "arrangement": arr, "fill": fill, "N_over_a": n / a, "series_turns": 12.0 * n / a,
           "R100_ohm": r100, "R_hot_ohm": float(series_resistance(ctx, n, a, w["n_s"], w["d_mm"], t_cu)), "t_cu_c": t_cu, "AT_cont": at_c, "AT_req": req["AT"],
           "AT_headroom_pct": 100.0 * (at_c / req["AT"] - 1.0), "limited_by_maps": rt["limited_by_maps"], "binding": rt["binding"], "I_cont": i_c, "I_cont_rms": i_c / math.sqrt(2.0),
           "gamma_cont": pt["gamma"], "T_cont": pt["torque"], "hot_spot_c": rt["hot_spot"], "hot_spot_slot_c": float(res.hot_spot_slot), "hot_spot_end_c": float(res.hot_spot_end),
           "magnet_c": rt["magnet"]}
    for k in ("p_dc_slot", "p_dc_end", "p_prox", "p_skin", "p_pwm", "p_circ", "p_total"):
        out[f"cont_{k}"] = float(lc[k])
    out["cont_iron_w"] = float(pt["iron"]["total"])
    out["cont_ripple_a_rms"], out["cont_modulation_index"] = float(lc["ripple_rms_a"]), float(lc["modulation_index"])
    # ---- the short-term current and the corner speeds at the winding's own rating
    at_st = min(ST_FACTOR * at_c, ctx.at_max)
    t_st = ctx.fm.mtpa(at_st / 10.0, stack_mm=ctx.stack)[3]
    i0d, i0q = ctx.fm.mtpa(at_c / 10.0, stack_mm=ctx.stack)[1:3]
    p0d, p0q = ctx.fm.psi_dq(i0d, i0q, stack_mm=ctx.stack)
    out.update({"AT_st": at_st, "I_st": at_st * a / n, "T_st": float(t_st), "short_term_clipped": bool(ST_FACTOR * at_c > ctx.at_max),
                "corner48_rpm": float(corner_speed(s * p0d, s * p0q, i0d / s, i0q / s, r100, V_DC / SQRT3)),
                "corner60_rpm": float(corner_speed(s * p0d, s * p0q, i0d / s, i0q / s, r100, V_DC_ALT / SQRT3))})
    # ---- how long the short-term current can be carried from a cold start at the rated speed (the network's transient with the copper loss as an effective current density)
    g_st, i_dst, i_qst = (float(v) for v in ctx.fm.mtpa(at_st * a / n, turns=n, paths=a, stack_mm=ctx.stack)[:3])
    lst = cl.losses(ww, {"i_peak": at_st * a / n, "gamma_deg": g_st, "rpm": RATED_RPM, "f_sw": case.f_sw, "v_dc": case.v_dc}, T_COPPER, ctx.stack)
    j_st = at_st / (math.sqrt(2.0) * ctx.tm.ampere_turns(1.0, fill)["a_cu_side_mm2"])
    dc_st = float(lst["p_dc_slot"] + lst["p_dc_end"])
    ist_iron = iron_loss(i_dst, i_qst, RATED_RPM, ctx.stack, turns=n, paths=a, lamination_mm=LAMINATION_MM)
    with np.errstate(all="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ov = ctx.tm.overload(1.0, j_st * math.sqrt(float(lst["p_total"]) / dc_st) if dc_st > 0 else j_st, start="cold", wind_limits=(case.wind_limit_c,), magnet_limit_c=case.magnet_limit_c,
                             fill=fill, rpm=RATED_RPM, iron={"iron_teeth": ist_iron["teeth"], "iron_yoke": ist_iron["yoke"], "iron_rotor": ist_iron["rotor"]}, t_max_s=3600.0)
    t_w, t_m = ov[f"winding_{case.wind_limit_c:.0f}C_s"], ov[f"magnet_{case.magnet_limit_c:.0f}C_s"]
    times = [(t, name) for t, name in ((t_w, "winding"), (t_m, "magnet")) if t is not None]
    out.update({"st_time_s": min(times)[0] if times else None, "st_first": min(times)[1] if times else None, "st_cu_loss_w": float(lst["p_total"]), "st_iron_w": float(ist_iron["total"])})
    # ---- the common-rating numbers (what the front used) and the high-speed point
    one = {"N": np.array([n]), "a": np.array([a]), "n_s": np.array([w["n_s"]]), "d_mm": np.array([w["d_mm"]]), "d_over_mm": np.array([w["d_over_mm"]]), "fill": np.array([fill])}
    chk = check_rows({"case": case, "arrangement": arr, "rows": [(n, a, w["n_s"], w["d_mm"], w["d_over_mm"])]})[0]
    obj = objectives(case, one, {k: np.array([v]) for k, v in chk.items()}, arr)
    out.update({"common_I_st": float(obj["I_st"][0]), "common_I_cont": float(obj["I_cont"][0]), "common_T_cont": float(obj["T_cont"][0]), "common_T_st": float(obj["T_st"][0]),
                "common_p_cu": chk["p_cu"], "common_corner_bus": float(obj["corner_bus"][0]), "common_corner48": float(obj["corner48"][0]), "common_corner60": float(obj["corner60"][0]), "common_margin_k": chk["margin"],
                "common_hot_spot_c": chk["hot_spot"], "common_magnet_c": chk["magnet"], "hs_feasible": bool(obj["hs_feasible"][0]), "hs_T_nm": float(obj["T_hs"]),
                "hs_I": float(obj["I_hs"][0]), "hs_gamma": float(obj["gamma_hs"][0]), "hs_p_total": float(obj["P_hs"][0]), "hs_extrapolated": bool(obj["hs_extrapolated"][0])})
    nan = float("nan")
    for k in ("p_dc_slot", "p_dc_end", "p_prox", "p_skin", "p_pwm", "p_circ", "p_total", "modulation_index", "iron_w", "hot_spot_c", "magnet_c", "margin_k", "fan_w"):
        out.setdefault(f"hs_{k}", nan)
    out.setdefault("hs_continuous", False)
    if out["hs_feasible"]:
        at_hs = out["hs_I"] * n / a
        lh = cl.losses(ww, {"i_peak": out["hs_I"], "gamma_deg": out["hs_gamma"], "rpm": case.hs_rpm, "f_sw": case.f_sw, "v_dc": case.v_dc}, T_COPPER, ctx.stack)
        for k in ("p_dc_slot", "p_dc_end", "p_prox", "p_skin", "p_pwm", "p_circ", "p_total"):
            out[f"hs_{k}"] = float(lh[k])
        out["hs_modulation_index"] = float(lh["modulation_index"])
        gam = math.radians(out["hs_gamma"])
        iron = iron_loss(-out["hs_I"] * math.sin(gam), out["hs_I"] * math.cos(gam), case.hs_rpm, ctx.stack, turns=n, paths=a, lamination_mm=LAMINATION_MM)
        f_hs, share_hs, _ = _copper_poly(ctx, w, arr, at_hs, case.hs_rpm, out["hs_gamma"])
        tm_hs = thermal.ThermalModel(None, ctx.stack, ctx.cooling, case.hs_rpm, thermal.WindingSpec(fill), ctx.assume)
        with np.errstate(all="ignore"), warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            rh = tm_hs.steady_state(thermal.Losses(copper=f_hs, copper_slot_fraction=share_hs, iron_teeth=iron["teeth"], iron_yoke=iron["yoke"], iron_rotor=iron["rotor"]), fill=fill,
                                    wind_limit_c=case.wind_limit_c, magnet_limit_c=case.magnet_limit_c)
        mh = min(rh.margin_winding, rh.margin_magnet)
        out.update({"hs_iron_w": float(iron["total"]), "hs_hot_spot_c": float(rh.hot_spot), "hs_magnet_c": float(rh.magnet), "hs_margin_k": float(mh),
                    "hs_continuous": bool(mh >= 0 and not rh.runaway), "hs_fan_w": float(tm_hs.fan_power())})
    out.update({f"make_{k}": v for k, v in manufacturing(w).items()})
    return out


# ---------------------------------------------------------------------------------------------- choosing
CLASSES = (15, 20, 30, 50, 75, 100, 150, 200, 300)          # inverter classes: the peak (short-term) phase current [A]


def _tie_key(run: CaseRun, j: int) -> tuple:
    c = run.cand
    return (int(c["n_s"][j]), int(c["a"][j]), int(c["N"][j]))


def pick_for_class(run: CaseRun, i_limit: float, rule: bool = True, cap: int | None = None, hs: bool = True) -> int | None:
    """The winding for an inverter of i_limit A peak: among the allowed windings (allowed(); rule=False drops the voltage rule) whose short-term current (the common rating's: objective f1;
    a winding's own rating, which detail() gives, is up to a few percent above it) is within the limit, the one
    with the highest corner speed on the case's bus; among windings within 0.2 % of it the lowest copper loss at the high-speed point, then the fewest strands in hand, then the
    fewest paths. None if no winding fits."""
    ok = np.flatnonzero(run.allowed(rule, cap, hs))
    ok = ok[run.obj["I_st"][ok] <= i_limit + 1e-9]
    if not len(ok):
        return None
    best = run.obj["corner_bus"][ok].max()
    near = ok[run.obj["corner_bus"][ok] >= (1.0 - 2e-3) * best]
    p = run.obj["P_hs"][near]
    near = near[p <= p.min() * (1.0 + 1e-3)]
    return int(min(near, key=lambda j: _tie_key(run, j)))


def recommend(run: CaseRun, cap: int | None = None) -> int | None:
    """The first study's rule on the new problem: the winding of smallest short-term current among those that satisfy every constraint (the voltage rule: the continuous torque to
    3,300 rpm); among windings of that current the lowest copper loss at the high-speed point, then the fewest strands in hand. None if none does."""
    ok = np.flatnonzero(run.allowed(True, cap))
    if not len(ok):
        return None
    i_min = run.obj["I_st"][ok].min()
    near = ok[run.obj["I_st"][ok] <= (1.0 + 1e-6) * i_min]
    p = run.obj["P_hs"][near]
    near = near[p <= p.min() * (1.0 + 1e-3)]
    return int(min(near, key=lambda j: _tie_key(run, j)))


# ---------------------------------------------------------------------------------------------- the speed range, from the motor-map engine
def speed_range(case: Case, w: dict, arrangement: str | None = None, rpms=(0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 8000.0, 10000.0, 12000.0, 14000.0, TOP_RPM)) -> dict | None:
    """The torque-speed envelope of the winding on the case's bus from optimization/motor_map.py (the other module's voltage-limit and field-weakening solver, on the same maps):
    the d-q torque at the largest torque inside the voltage limit (V_dc / sqrt 3) and the current limit, for the winding's continuous current and its short-term one, at each speed
    (copper at 100 C, the full copper loss in the engine's own accounting), and the corner speed it finds. None if the engine cannot be loaded."""
    try:
        from optimization import motor_map as mm
    except Exception:                                                                       # the engine is another module's, in progress: its absence must not stop this study
        return None
    ctx = get_context(case)
    arr = arrangement or case.arrangement
    d = detail(case, w, arr)
    mw = mm.Winding(float(w["N"]), int(w["a"]), float(w["n_s"]), float(w["d_mm"]), float(w["d_over_mm"]))
    m = mm.Machine(mw, ctx.stack, copper=mm.FullCopper(f_sw=case.f_sw, arrangement=arr))
    bus = mm.Bus(case.v_dc, k=1.0)
    out = {"rpm": np.asarray(rpms, float), "I_cont": d["I_cont"], "I_st": d["I_st"]}
    sc = mm.Scan(m, i_max=min(d["I_st"], m.i_cap))
    for key, lim in (("cont", d["I_cont"]), ("st", d["I_st"])):
        env = mm.envelope(m, bus, out["rpm"], min(lim, m.i_cap), basis="air_gap", scan=sc, locate=True)
        out[f"T_{key}"] = np.asarray(env["T"], float)
        out[f"region_{key}"] = [str(r) for r in env["region"]]
        out[f"corner_{key}"] = env.meta.get("corner_rpm")
        out[f"zero_{key}"] = env.meta.get("zero_rpm")
    return out


# ---------------------------------------------------------------------------------------------- reports
class Report:
    """What a stage says: printed and kept, so that the run's text goes to report.txt."""

    def __init__(self):
        self.lines: list[str] = []

    def say(self, text: str = "") -> None:
        print(text, flush=True)
        self.lines.append(text)


def _fmt(v, nd: int = 1) -> str:
    if v is None:
        return "-"
    if isinstance(v, (bool, np.bool_)):
        return "yes" if v else "no"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, float) and not np.isfinite(v):
        return "-" if np.isnan(v) else "inf"
    return f"{v:.{nd}f}"


def front_table(rep: Report, run: CaseRun, idx: np.ndarray, title: str) -> None:
    """The windings of idx (indices of the run's candidates), by short-term current, with the numbers the front used and the network's temperatures at the common rating."""
    c, o, k = run.cand, run.obj, run.check
    rep.say(title)
    rep.say("    N  a  n_s    d [mm]   fill   R100 [mOhm]  I_cont [A]  I_st [A]  corner 48 V  corner 60 V  P_hs [W]  I_hs [A]  gamma_hs  hot spot  magnets  margin [K]")
    for j in idx[np.argsort(o["I_st"][idx], kind="stable")]:
        rep.say(f"  {c['N'][j]:3d} {c['a'][j]:2d} {c['n_s'][j]:4d}   {c['d_mm'][j]:6.3f}  {c['fill'][j]:6.3f}  {1e3 * o['R100'][j]:10.2f}  {o['I_cont'][j]:10.1f}  {o['I_st'][j]:8.1f}  {o['corner48'][j]:11.0f}  "
                f"{o['corner60'][j]:11.0f}  {o['P_hs'][j]:8.1f}  {o['I_hs'][j]:8.1f}  {o['gamma_hs'][j]:8.1f}  {k['hot_spot'][j]:8.1f}  {k['magnet'][j]:7.1f}  {k['margin'][j]:10.2f}")


def winding_block(rep: Report, d: dict, title: str) -> None:
    """One winding in full, as the shortlist wants it."""
    rep.say(title)
    rep.say(f"   winding        N {d['N']} turns a coil side, a {d['a']} paths, {d['n_s']} strands of {d['d_mm']:g} mm ({d['d_over_mm']:g} mm over the enamel) in a {d['arrangement']} bundle; fill {d['fill']:.3f}; "
            f"{d['series_turns']:g} series turns a phase; R(100 C) {1e3 * d['R100_ohm']:.2f} mOhm")
    rep.say(f"   rating         at the common rating ({d['AT_req']:.1f} A-turns a coil side, {d['common_T_cont']:.2f} N.m: what every winding is held to and the class current is read at): continuous "
            f"{d['common_I_cont']:.1f} A peak, short-term {d['common_I_st']:.1f} A peak ({d['common_T_st']:.2f} N.m)")
    rep.say(f"                  its own thermal rating {d['AT_cont']:.1f} A-turns = {d['I_cont']:.1f} A peak ({d['I_cont_rms']:.1f} A rms), {d['T_cont']:.2f} N.m at {d['gamma_cont']:.1f} deg "
            f"(headroom {d['AT_headroom_pct']:+.1f} % over the common rating, {d['binding']} binding); short-term {d['I_st']:.1f} A peak, {d['T_st']:.2f} N.m, "
            f"{_fmt(d['st_time_s'], 0)} s from cold until the {d['st_first']} limit")
    rep.say(f"   speed range    corner speed at the continuous torque: at the common rating {d['common_corner48']:.0f} rpm on 48 V, {d['common_corner60']:.0f} rpm on 60 V; "
            f"at its own rating {d['corner48_rpm']:.0f} and {d['corner60_rpm']:.0f} rpm")
    rep.say(f"   temperatures   at 3,000 rpm: at the common rating hot spot {d['common_hot_spot_c']:.1f} C, magnets {d['common_magnet_c']:.1f} C; at its own rating hot spot {d['hot_spot_c']:.1f} C "
            f"(slot {d['hot_spot_slot_c']:.1f}, end turn {d['hot_spot_end_c']:.1f}), magnets {d['magnet_c']:.1f} C, copper {d['t_cu_c']:.1f} C")
    rep.say(f"   copper loss    3,000 rpm at the winding's own continuous rating [W]: DC slots {d['cont_p_dc_slot']:.2f}, DC end turns {d['cont_p_dc_end']:.2f}, proximity {d['cont_p_prox']:.3f}, skin {d['cont_p_skin']:.5f}, "
            f"PWM {d['cont_p_pwm']:.3f}, circulating {d['cont_p_circ']:.3f}; total {d['cont_p_total']:.2f} (iron {d['cont_iron_w']:.1f})")
    if d["hs_feasible"]:
        rep.say(f"                  {HS_RPM:,.0f} rpm, {d['hs_T_nm']:.2f} N.m, {d['hs_I']:.1f} A at {d['hs_gamma']:.1f} deg (copper 100 C) [W]: DC slots {d['hs_p_dc_slot']:.2f}, DC end turns {d['hs_p_dc_end']:.2f}, "
                f"proximity {d['hs_p_prox']:.3f}, skin {d['hs_p_skin']:.5f}, PWM {d['hs_p_pwm']:.3f}, circulating {d['hs_p_circ']:.3f}; total {d['hs_p_total']:.2f} (iron {d['hs_iron_w']:.1f}); "
                f"the network there: hot spot {d['hs_hot_spot_c']:.1f} C, magnets {d['hs_magnet_c']:.1f} C ({'sustainable' if d['hs_continuous'] else 'NOT sustainable: a short-term demand'})")
    else:
        rep.say(f"                  the {HS_RPM:,.0f} rpm demand cannot be held on this bus within the short-term current")
    rep.say(f"   to wind it     {d['make_n_s']} strands in hand of {d['make_d_mm']:g} mm, {d['make_strands_abreast']} abreast through the opening ({d['make_opening_room_mm']:.2f} mm to spare); "
            f"a strand is {d['make_stiffness_vs_0315']:.1f} x as stiff as the 0.315 mm of the first study's example, the bundle {d['make_bundle_stiffness_vs_example']:.1f} x as stiff as that example's 8 x 0.315 mm; "
            f"{d['make_strands_per_coil_side']} strands in a coil side, {d['make_turns_per_coil']} turns a coil, {d['make_coils_per_path']} coils in a path")


def arrangement_costs(case: Case, w: dict, arrangements=("transposed", "row", "cluster", "random", "line")) -> list[dict]:
    """The winding in each arrangement: its own continuous rating, the copper loss at 3,000 rpm and at the high-speed point, the network's temperatures, and what each costs against
    the case's arrangement."""
    rows = [detail(case, w, a) for a in arrangements]
    base = next(r for r in rows if r["arrangement"] == case.arrangement) if any(r["arrangement"] == case.arrangement for r in rows) else rows[0]
    for r in rows:
        r["d_AT_pct"] = 100.0 * (r["AT_cont"] / base["AT_cont"] - 1.0)
        r["d_cont_loss_pct"] = 100.0 * (r["cont_p_total"] / base["cont_p_total"] - 1.0)
        r["d_hs_loss_pct"] = 100.0 * (r["hs_p_total"] / base["hs_p_total"] - 1.0) if np.isfinite(base["hs_p_total"]) and np.isfinite(r["hs_p_total"]) else float("nan")
    return rows


def arrangement_table(rep: Report, rows: list[dict]) -> None:
    rep.say("   arrangement   AT_cont  d AT [%]  copper 3,000 rpm [W] (circulating)  d [%]   copper at the high-speed point [W] (circulating)  d [%]   hot spot / magnets there [C]")
    for r in rows:
        rep.say(f"   {r['arrangement']:11s}  {r['AT_cont']:7.2f}  {r['d_AT_pct']:+7.2f}  {r['cont_p_total']:12.2f} ({r['cont_p_circ']:6.3f})  {r['d_cont_loss_pct']:+6.2f}   "
                f"{r['hs_p_total']:14.2f} ({r['hs_p_circ']:7.3f})  {r['d_hs_loss_pct']:+7.2f}      {r['hs_hot_spot_c']:6.1f} / {r['hs_magnet_c']:6.1f}")


def classes_for(run: CaseRun) -> list[float]:
    """The inverter classes worth showing for the case: the ladder CLASSES from just below the smallest current that satisfies every constraint up to 4 x it."""
    ok = np.flatnonzero(run.allowed(True))
    if not len(ok):
        return [c for c in CLASSES]
    lo, hi = run.obj["I_st"][ok].min(), run.obj["I_st"][ok].max()
    return [c for c in CLASSES if c >= 0.5 * lo and c <= 1.35 * hi] or [c for c in CLASSES]


def shortlist(rep: Report, run: CaseRun, classes=None) -> list[dict]:
    """The best winding of each inverter class (pick_for_class), in full (detail), with what the classes below the voltage rule would give."""
    classes = classes or classes_for(run)
    rows = []
    seen = set()
    for lim in classes:
        j = pick_for_class(run, lim)
        if j is None:
            jr = pick_for_class(run, lim, rule=False, cap=None, hs=False)
            alt = None
            if jr is not None:
                alt = f"best without the voltage rule: N {run.cand['N'][jr]}, a {run.cand['a'][jr]}, {run.cand['n_s'][jr]} x {run.cand['d_mm'][jr]:g} mm, corner speed {run.obj['corner_bus'][jr]:.0f} rpm"
            rep.say(f"  class {lim:g} A: no winding satisfies the constraints (the voltage rule needs {run.obj['I_st'][np.flatnonzero(run.allowed(True))].min():.1f} A or more); {alt or 'no winding at all'}")
            continue
        key = (int(run.cand["N"][j]), int(run.cand["a"][j]), int(run.cand["n_s"][j]), float(run.cand["d_mm"][j]))
        if key in seen:
            rep.say(f"  class {lim:g} A: the same winding as the class below")
            continue
        seen.add(key)
        d = detail(run.case, winding_dict(run.cand, j), run.arrangement)
        d["class_A"] = lim
        winding_block(rep, d, f"  CLASS {lim:g} A peak (short-term):")
        rows.append(d)
    return rows


def flatten(d: dict) -> dict:
    return {k: (round(v, 6) if isinstance(v, float) else v) for k, v in d.items()}


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(flatten(r) for r in rows)


def stage_headline(rep: Report, out: Path, workers: int = 1) -> dict:
    """The headline case (60 mm, 48 V, the shaft fan held to 50 W at the top speed, N45 at 80 C, compact bundles, 8 strands at most): the setup, the front, the shortlist, the
    recommended winding with what each arrangement would cost, the phase-1 rules applied, and the speed range."""
    case = Case()
    run = run_case(case, workers)
    ctx = get_context(case)
    req = requirement(ctx)
    c, o, k = run.cand, run.obj, run.check
    n_all, n_cap = len(c["N"]), int((c["n_s"] <= case.n_s_cap).sum())
    in_cap = c["n_s"] <= case.n_s_cap
    rep.say("=" * 190)
    rep.say(f"THE HEADLINE CASE: {case.label()}")
    rep.say(f"  cooling {ctx.cooling.label}; the shaft fan's power {thermal.shaft_fan_power(RATED_RPM, ctx.cooling, ctx.stack):.1f} W at 3,000 rpm, {thermal.shaft_fan_power(TOP_RPM, ctx.cooling, ctx.stack):.1f} W at {TOP_RPM:,.0f} rpm")
    rep.say(f"  the reference winding (26 turns, 4 paths, 8 x 0.315 mm, fill 0.448) carries {req['AT_ref']:.2f} A-turns a coil side continuously at 3,000 rpm ({req['T_ref']:.3f} N.m, class F / N45 80 C, "
            f"full copper loss); the common rating is {100 * (1 - case.rating_margin):.0f} % of it: {req['AT']:.2f} A-turns, {req['T']:.3f} N.m; short-term {ST_FACTOR:g} x; "
            f"the high-speed demand {case.hs_fraction:g} x {req['T']:.3f} = {case.hs_fraction * req['T']:.3f} N.m at {case.hs_rpm:,.0f} rpm")
    rep.say(f"  windings enumerated {n_all:,} (up to 30 strands; {n_cap:,} with at most {case.n_s_cap}); the thermal network solved {int(k['solved'].sum()):,} of them (the others have a fill below "
            f"{k['fill_floor']:.3f}, which no winding can carry the rating with); thermal limits kept {int(o['thermal_ok'].sum()):,}; of those within {case.n_s_cap} strands {int((o['thermal_ok'] & in_cap).sum()):,}; "
            f"the voltage rule (corner speed >= {RULE_RPM:,.0f} rpm on 48 V) then leaves {int((o['thermal_ok'] & in_cap & (o['corner_bus'] >= RULE_RPM)).sum())}; "
            f"the high-speed demand {int(run.allowed(True).sum())}; the exact front has {len(run.front())} designs. {run.seconds:.1f} s.")
    allowed_idx = np.flatnonzero(run.allowed(True))
    front_idx = run.front(True)
    rep.say(f"  the front is verified against all {len(allowed_idx)} allowed designs: "
            f"{winding.verify_front(np.stack([o['I_st'][allowed_idx], -o['corner_bus'][allowed_idx], o['P_hs'][allowed_idx]], axis=1), np.flatnonzero(np.isin(allowed_idx, front_idx)))}")
    rep.say()
    front_table(rep, run, allowed_idx, "ALL ALLOWED WINDINGS (the constraints satisfied), by short-term current (currents, corner speeds and temperatures at the common rating); the front follows")
    rep.say()
    front_table(rep, run, front_idx, "* THE EXACT FRONT (I_st minimised, corner speed on 48 V maximised, copper loss at the high-speed point minimised)")
    rep.say()
    # ---- the shortlist
    rep.say("THE SHORTLIST: the best winding of each inverter class (the highest corner speed within the class, then the lowest copper loss at the high-speed point, then the fewest strands in hand)")
    rows = shortlist(rep, run)
    write_csv(out / "shortlist_60mm.csv", rows)
    # ---- the recommendation
    j = recommend(run)
    rec = None
    if j is not None:
        w = winding_dict(run.cand, j)
        rec = detail(case, w, case.arrangement)
        rep.say()
        winding_block(rep, rec, "THE RECOMMENDED HEADLINE WINDING (the smallest inverter current that still holds the continuous torque to 3,300 rpm on 48 V; the lowest copper loss at the high-speed point)")
        rep.say("  what each arrangement would cost this winding:")
        arrangement_table(rep, arrangement_costs(case, w))
    return {"case": case, "run": run, "rows": rows, "rec": rec, "rec_index": j, "req": req}


def cap_table(rep: Report, run: CaseRun, caps=(4, 6, 8, 12, 30)) -> list[dict]:
    """What the cap on strands in hand costs: for each cap the number of allowed windings, the smallest short-term current that satisfies every constraint, and the best winding at that
    current (strands, wire, copper loss at the high-speed point), with the loss against the cap-30 optimum at the same N / a."""
    rows = []
    c, o = run.cand, run.obj
    best30 = {}
    ok30 = np.flatnonzero(run.allowed(True, 30))
    for j in ok30:
        key = round(c["N"][j] / c["a"][j], 6)
        if key not in best30 or o["P_hs"][j] < o["P_hs"][best30[key]]:
            best30[key] = j
    rep.say("  cap   allowed  front | smallest I_st [A]  best winding there (N a n_s x d)         fill    P_hs [W]  vs the cap-30 optimum at its N/a  |  mean P_hs excess over the cap-30 optimum, all N/a that exist at the cap")
    for cap in caps:
        idx = np.flatnonzero(run.allowed(True, cap))
        row = {"cap": cap, "allowed": len(idx), "front": len(run.front(True, cap))}
        if len(idx):
            i_min = o["I_st"][idx].min()
            near = idx[o["I_st"][idx] <= (1 + 1e-6) * i_min]
            j = near[np.argmin(o["P_hs"][near])]
            key = round(c["N"][j] / c["a"][j], 6)
            ex = []
            for noa in {round(c["N"][i] / c["a"][i], 6) for i in idx}:
                sel = idx[np.isclose(c["N"][idx] / c["a"][idx], noa)]
                ex.append(100.0 * (o["P_hs"][sel].min() / o["P_hs"][best30[noa]] - 1.0))
            row.update({"i_st_min": float(i_min), "N": int(c["N"][j]), "a": int(c["a"][j]), "n_s": int(c["n_s"][j]), "d_mm": float(c["d_mm"][j]), "fill": float(c["fill"][j]),
                        "P_hs": float(o["P_hs"][j]), "vs_best30_pct": 100.0 * (o["P_hs"][j] / o["P_hs"][best30[key]] - 1.0), "mean_excess_pct": float(np.mean(ex)), "max_excess_pct": float(np.max(ex))})
            rep.say(f"  {cap:3d}  {len(idx):7d}  {row['front']:5d} | {i_min:16.1f}  {c['N'][j]:3d} {c['a'][j]} {c['n_s'][j]:3d} x {c['d_mm'][j]:5.3f}                  {c['fill'][j]:6.3f}  {o['P_hs'][j]:8.2f}  "
                    f"{row['vs_best30_pct']:+7.2f} %                          | {row['mean_excess_pct']:+6.2f} % (largest {row['max_excess_pct']:+6.2f} %)")
        else:
            rep.say(f"  {cap:3d}  {len(idx):7d}  {row['front']:5d} | no winding satisfies the constraints")
        rows.append(row)
    return rows


def compact_header(rep: Report) -> None:
    rep.say("    N  a  n_s    d [mm]  arr.      fill  R100 [mOhm]  | at the common rating: I_cont  I_st [A]  corner 48 V  corner 60 V  P 3,000 rpm [W]  P high-speed [W]  hot spot / magnets [C]  | own rating: AT_cont  I_cont  I_st [A]")


def compact_row(rep: Report, d: dict) -> None:
    rep.say(f"  {d['N']:3d} {d['a']:2d} {d['n_s']:4d}   {d['d_mm']:6.3f}  {d['arrangement'][:8]:8s} {d['fill']:5.3f}  {1e3 * d['R100_ohm']:10.2f}  | {d['common_I_cont']:30.1f}  {d['common_I_st']:8.1f}  "
            f"{d['common_corner48']:11.0f}  {d['common_corner60']:11.0f}  {d['common_p_cu']:14.2f}  {d['hs_p_total']:15.2f}  {d['common_hot_spot_c']:10.1f} / {d['common_magnet_c']:5.1f}       | "
            f"{d['AT_cont']:8.1f}  {d['I_cont']:6.1f}  {d['I_st']:8.1f}")


def case_summary(rep: Report, name: str, run: CaseRun) -> dict:
    """One line about a case: its rating, the allowed windings, the recommendation. Returns the numbers."""
    ctx = get_context(run.case)
    req = requirement(ctx)
    j = recommend(run)
    c, o = run.cand, run.obj
    row = {"case": name, "AT_ref": req["AT_ref"], "AT_req": req["AT"], "T_req": req["T"], "allowed": int(run.allowed(True).sum()), "front": len(run.front(True)), "seconds": run.seconds}
    s = (f"  {name:28s} AT_ref {req['AT_ref']:6.1f}  common rating {req['AT']:6.1f} A-turns ({req['T']:5.2f} N.m)  allowed {row['allowed']:3d}  front {row['front']:2d}")
    if j is not None:
        d = detail(run.case, winding_dict(c, j), run.arrangement)
        row.update({"N": d["N"], "a": d["a"], "n_s": d["n_s"], "d_mm": d["d_mm"], "fill": d["fill"], "I_cont": d["I_cont"], "I_st": d["I_st"], "corner48": d["corner48_rpm"],
                    "corner60": d["corner60_rpm"], "AT_cont": d["AT_cont"], "P_cont": d["cont_p_total"], "P_hs": d["hs_p_total"], "T_cont": d["T_cont"], "st_time_s": d["st_time_s"],
                    "I_cont_common": d["common_I_cont"], "I_st_common": d["common_I_st"], "corner48_common": d["common_corner48"], "corner60_common": d["common_corner60"]})
        s += (f" | recommended N {d['N']} a {d['a']} {d['n_s']} x {d['d_mm']:g} mm fill {d['fill']:.3f}: at the common rating {d['common_I_cont']:.1f} A continuous, {d['common_I_st']:.1f} A short-term, "
              f"corner {d['common_corner48']:.0f} / {d['common_corner60']:.0f} rpm (48 / 60 V); its own rating {d['AT_cont']:.1f} A-turns = {d['I_cont']:.1f} A, {d['T_cont']:.2f} N.m; P_hs {d['hs_p_total']:.1f} W")
    else:
        idx = np.flatnonzero(run.allowed(False, None, False))
        if len(idx):
            b = idx[np.argmax(o["corner_bus"][idx])]
            s += f" | NO winding satisfies the voltage rule; the best reaches {o['corner_bus'][b]:.0f} rpm (N {c['N'][b]}, a {c['a'][b]}, {c['n_s'][b]} x {c['d_mm'][b]:g} mm, {o['I_st'][b]:.1f} A short-term)"
    rep.say(s)
    return row


def class_table(rep: Report, run: CaseRun, title: str) -> list[dict]:
    """The best winding of each inverter class of a run (pick_for_class: the class limit is on the short-term current of the common rating, the objective f1), one compact row each; a
    class below the voltage rule says so and names the best winding without it. Returns the detail rows (with 'class_A')."""
    rep.say(title)
    compact_header(rep)
    rows = []
    for lim in classes_for(run):
        j = pick_for_class(run, lim)
        if j is None:
            jr = pick_for_class(run, lim, rule=False, hs=False)
            note = "" if jr is None else (f" (the best without the rule: N {run.cand['N'][jr]}, a {run.cand['a'][jr]}, {run.cand['n_s'][jr]} x {run.cand['d_mm'][jr]:g} mm, "
                                          f"corner {run.obj['corner_bus'][jr]:.0f} rpm)")
            rep.say(f"  class {lim:g} A: below the voltage rule{note}")
            continue
        d = detail(run.case, winding_dict(run.cand, j), run.arrangement)
        d["class_A"] = lim
        compact_row(rep, d)
        rows.append(d)
    return rows


def rule_table(rep: Report, run: CaseRun, margins=(0.0, 0.05, 0.10, 0.15, 0.20, 0.30)) -> list[dict]:
    """The voltage rule's margin: for each margin of the corner speed at the continuous torque over 3,000 rpm (10 % is the rule), the winding that satisfies every other constraint with the
    smallest short-term current (then the lowest copper loss at the high-speed point, the fewest strands, the fewest paths), and what it costs."""
    c, o = run.cand, run.obj
    base = np.flatnonzero(run.allowed(False))
    rows = []
    rep.say("  margin  corner >= [rpm]    N  a  n_s    d [mm]   I_st [A]  corner [rpm]  P_hs [W]  windings")
    for m in margins:
        thr = RATED_RPM * (1.0 + m)
        ok = base[o["corner_bus"][base] >= thr]
        if not len(ok):
            rep.say(f"  {100 * m:4.0f} %  {thr:14.0f}    no winding")
            continue
        near = ok[o["I_st"][ok] <= (1.0 + 1e-6) * o["I_st"][ok].min()]
        near = near[o["P_hs"][near] <= o["P_hs"][near].min() * (1.0 + 1e-3)]
        j = int(min(near, key=lambda k: _tie_key(run, k)))
        rows.append({"margin_pct": 100.0 * m, "rule_rpm": thr, "N": int(c["N"][j]), "a": int(c["a"][j]), "n_s": int(c["n_s"][j]), "d_mm": float(c["d_mm"][j]), "I_st": float(o["I_st"][j]),
                     "corner": float(o["corner_bus"][j]), "P_hs": float(o["P_hs"][j]), "windings": int(len(ok))})
        rep.say(f"  {100 * m:4.0f} %  {thr:14.0f}   {c['N'][j]:3d} {c['a'][j]:2d} {c['n_s'][j]:4d}   {c['d_mm'][j]:6.3f}  {o['I_st'][j]:9.1f}  {o['corner_bus'][j]:12.0f}  {o['P_hs'][j]:8.1f}  {len(ok):8d}")
    return rows


def stage_variants(rep: Report, out: Path, workers: int = 1, headline: dict | None = None) -> dict:
    """The cap on strands in hand, the rating allowance, the stacks, the bus, the magnet grade, the fan's sizing and the arrangement, each against the headline case."""
    h = headline or stage_headline(Report(), out, workers)
    base: Case = h["case"]
    res: dict = {}
    # ---- 1. the cap on strands in hand, on every stack
    rep.say("=" * 190)
    rep.say("THE CAP ON STRANDS IN HAND: for each cap the number of allowed windings, the smallest short-term current that satisfies every constraint and the best winding there")
    res["cap"] = {}
    stacks = {60.0: h["run"]}
    for stack in (20.0, 40.0, 80.0):
        stacks[stack] = run_case(replace(base, stack_mm=stack), workers)
    for stack in (20.0, 40.0, 60.0, 80.0):
        rep.say(f" {stack:g} mm stack:")
        res["cap"][stack] = cap_table(rep, stacks[stack])
    # ---- 2. the stacks
    rep.say()
    rep.say("=" * 190)
    rep.say("THE STACKS (48 V, the shaft fan sized for the stack, N45 80 C, compact bundles, at most 8 strands): the recommended winding and the shortlist of each")
    res["stacks"], res["stack_shortlists"] = {}, {}
    res["fronts"] = {}
    for stack in (20.0, 40.0, 60.0, 80.0):
        run = stacks[stack]
        res["stacks"][stack] = case_summary(rep, f"{stack:g} mm", run)
        res["fronts"][stack] = [dict(winding_dict(run.cand, int(j)), I_st=float(run.obj["I_st"][j]), corner_bus=float(run.obj["corner_bus"][j]), P_hs=float(run.obj["P_hs"][j]))
                                for j in run.front(True)]
    for stack in (20.0, 40.0, 60.0, 80.0):
        run = stacks[stack]
        rep.say()
        rows = class_table(rep, run, f" {stack:g} mm: the best winding of each inverter class")
        for d in rows:
            d["stack_mm"] = stack
        res["stack_shortlists"][stack] = rows
        write_csv(out / f"shortlist_{stack:g}mm.csv", rows)
    # ---- 3. the rating allowance
    rep.say()
    rep.say("=" * 190)
    rep.say("THE RATING ALLOWANCE (60 mm): how far below the reference winding's thermal rating a winding may fall (its own fill and AC losses)")
    res["allowance"] = []
    for m in (0.0, 0.03, 0.05, 0.08, 0.12):
        run = h["run"] if abs(m - base.rating_margin) < 1e-12 else run_case(replace(base, rating_margin=m), workers)
        res["allowance"].append(case_summary(rep, f"allowance {100 * m:.0f} %", run))
    # ---- 3b. the voltage rule's margin
    rep.say()
    rep.say("=" * 190)
    rep.say("THE VOLTAGE RULE'S MARGIN (60 mm, 48 V): the winding of smallest short-term current whose corner speed at the continuous torque is at least 3,000 rpm plus the margin (10 % is the rule)")
    res["rule"] = rule_table(rep, h["run"])
    # ---- 4. the bus, the magnet grade, the fan
    rep.say()
    rep.say("=" * 190)
    rep.say("THE BUS, THE MAGNET GRADE AND THE FAN (60 mm)")
    res["variants"] = {}
    for name, case in (("48 V (headline)", base), ("60 V bus", replace(base, v_dc=V_DC_ALT)), ("150 C magnets", replace(base, magnet_limit_c=150.0)),
                       ("fan (a): 6 m/s at 3,000 rpm", replace(base, cooling="fan_a")), ("fan (a) and 150 C magnets", replace(base, cooling="fan_a", magnet_limit_c=150.0))):
        run = h["run"] if case == base else run_case(case, workers)
        res["variants"][name] = case_summary(rep, name, run)
        if name == "60 V bus":
            run60 = run
    rep.say()
    res["shortlist_60V"] = class_table(rep, run60, " 60 mm on the 60 V bus: the best winding of each inverter class")
    write_csv(out / "shortlist_60mm_60V.csv", res["shortlist_60V"])
    # ---- 5. the arrangement, whole runs
    rep.say()
    rep.say("=" * 190)
    rep.say("THE ARRANGEMENT OVER THE WHOLE ENUMERATION (60 mm): the recommended winding when every strand placement is assumed")
    res["arrangement"] = {}
    for arr in ("transposed", "row", "cluster", "random", "line"):
        case = replace(base, arrangement=arr)
        run = h["run"] if arr == base.arrangement else run_case(case, workers)
        res["arrangement"][arr] = case_summary(rep, f"arrangement {arr}", run)
    return res


def stage_envelopes(rep: Report, h: dict) -> list[dict]:
    """The torque-speed envelope on the case's bus of the recommended winding and of each shortlisted one, from motor_map.py's engine (speed_range): the d-q torque at the largest torque
    inside the voltage limit and the current limit, at the winding's own continuous and short-term currents. [] if the engine cannot be loaded."""
    case = h["case"]
    wins = ([("recommended", h["rec"])] if h["rec"] is not None else []) + [(f"class {d['class_A']:g} A", d) for d in h["rows"]]
    out: list[dict] = []
    for name, d in wins:
        w = {"N": d["N"], "a": d["a"], "n_s": d["n_s"], "d_mm": d["d_mm"], "d_over_mm": d["d_over_mm"]}
        sr = speed_range(case, w, case.arrangement)
        if sr is None:
            return []
        sr.update({"name": name, "N": d["N"], "a": d["a"], "n_s": d["n_s"], "d_mm": d["d_mm"]})
        out.append(sr)
    if not out:
        return out
    rep.say(f"THE SPEED RANGE ON {case.v_dc:g} V OF THE RECOMMENDED WINDING AND OF THE SHORTLIST (motor_map.py's solver: the d-q torque at the largest torque inside the voltage limit V_dc / sqrt 3 "
            f"and the winding's current limit, copper at 100 C, the winding's own continuous and short-term currents) [N.m]")
    rep.say("    speed [rpm]" + " " * 44 + " ".join(f"{r:6.0f}" for r in out[0]["rpm"]))
    for e in out:
        for key, lab in (("cont", "continuous"), ("st", "short-term")):
            rep.say(f"    {e['name']:11s} N {e['N']:2d} a {e['a']} {e['n_s']} x {e['d_mm']:5.3f}  {lab:10s} {e['I_' + key]:6.1f} A  " + " ".join(f"{t:6.2f}" for t in e["T_" + key])
                    + f"   corner {_fmt(e['corner_' + key], 0)} rpm, no torque above {_fmt(e['zero_' + key], 0)} rpm")
    return out


def rules_applied(rep: Report, h: dict) -> dict:
    """The phase-1 rules against this study's numbers: (1) strands of 0.4 mm or less toward 16,000 rpm, (2) compact bundles, (3) N / a high enough for the ripple, (4) thin strands near
    the opening."""
    case, run, rec = h["case"], h["run"], h["rec"]
    ctx = get_context(case)
    c, o, k = run.cand, run.obj, run.check
    out: dict = {}
    rep.say("=" * 190)
    rep.say("THE PHASE-1 RULES APPLIED")
    # ---- 1. fine strands
    j0 = h["rec_index"]
    noa = c["N"][j0] / c["a"][j0]
    fine = np.flatnonzero(run.allowed(True, 30) & (c["d_mm"] <= 0.4 + 1e-9) & np.isclose(c["N"] / c["a"], noa))
    rep.say(f"  1. STRANDS OF 0.4 MM OR LESS toward 16,000 rpm. The recommended winding has N / a = {noa:g}; the best thinner winding that the constraints allow at that N / a (any number of strands):")
    top_hs = {}
    for name, j in (("recommended", j0), ("<= 0.4 mm", int(fine[np.argmin(o["P_hs"][fine])]) if len(fine) else None)):
        if j is None:
            rep.say("     no winding of 0.4 mm or less satisfies the constraints at that N / a")
            continue
        w = winding_dict(c, j)
        d = detail(case, w, case.arrangement)
        gam = d["hs_gamma"]
        lt = cl.losses(dict(w, arrangement=case.arrangement), {"i_peak": d["hs_I"], "gamma_deg": gam, "rpm": TOP_RPM, "f_sw": case.f_sw, "v_dc": case.v_dc}, T_COPPER, ctx.stack)
        rep.say(f"     {name:12s} {w['n_s']} x {w['d_mm']:g} mm (N {w['N']}, a {w['a']}), fill {d['fill']:.3f}: copper at {HS_RPM:,.0f} rpm {d['hs_p_total']:.2f} W = DC {d['hs_p_dc_slot'] + d['hs_p_dc_end']:.2f} + "
                f"proximity {d['hs_p_prox']:.3f} + circulating {d['hs_p_circ']:.3f} + PWM {d['hs_p_pwm']:.3f}; at the top speed {TOP_RPM:,.0f} rpm and the same current: proximity {float(lt['p_prox']):.3f}, "
                f"circulating {float(lt['p_circ']):.3f}, ({100 * float(lt['p_prox'] + lt['p_circ']) / float(lt['p_dc_slot'] + lt['p_dc_end']):.2f} % of the DC loss)")
        top_hs[name] = (d, lt)
    out["fine"] = {k: v[0] for k, v in top_hs.items()}
    # ---- 2. compact bundles: the table printed with the recommendation; here the share
    rows = arrangement_costs(case, winding_dict(c, j0))
    out["arrangements"] = rows
    rand = next(r for r in rows if r["arrangement"] == "random")
    clus = next(r for r in rows if r["arrangement"] == "cluster")
    rep.say(f"  2. COMPACT BUNDLES. Strands of a turn spread over the slot's depth (random, line) add {rand['hs_p_circ']:.1f} W of circulating-current loss at {HS_RPM:,.0f} rpm ({rand['d_hs_loss_pct']:+.1f} % of the "
            f"copper loss there, {rand['d_AT_pct']:+.1f} % of the continuous rating, magnets {rand['hs_magnet_c'] - clus['hs_magnet_c']:+.1f} K hotter); a compact bundle (cluster) or a row adds {clus['hs_p_circ']:.2f} W "
            f"({100 * clus['hs_p_circ'] / clus['hs_p_total']:.2f} % of it). The rule holds, with the strands of each turn kept together.")
    # ---- 3. the ripple
    rep.say("  3. N / a HIGH ENOUGH FOR THE RIPPLE (20 kHz, 48 V): the PWM ripple's loss on the windings of the shortlist")
    rep.say("       N  a  N/a   ripple at 3,000 rpm [A rms]  PWM loss at 3,000 rpm [W] (share of the copper loss)   PWM loss at the high-speed point [W] (share)")
    out["ripple"] = []
    for d in h["rows"] + ([h["rec"]] if h["rec"] is not None else []):
        s3 = 100 * d["cont_p_pwm"] / d["cont_p_total"]
        s4 = 100 * d["hs_p_pwm"] / d["hs_p_total"] if d["hs_feasible"] else float("nan")
        rep.say(f"     {d['N']:3d} {d['a']:2d} {d['N_over_a']:5.2f}   {d['cont_ripple_a_rms']:26.3f}  {d['cont_p_pwm']:22.3f} ({s3:5.2f} %)  {d['hs_p_pwm']:32.3f} ({s4:5.2f} %)")
        out["ripple"].append({"N": d["N"], "a": d["a"], "ripple_a_rms": d["cont_ripple_a_rms"], "p_pwm": d["cont_p_pwm"], "share_pct": s3, "hs_p_pwm": d["hs_p_pwm"], "hs_share_pct": s4})
    # ---- 4. thin strands near the opening
    d = h["rec"]
    saving = 0.9 * d["hs_p_prox"]
    rep.say(f"  4. THIN STRANDS NEAR THE OPENING. The machine's field makes about 90 % of the proximity loss in the layer next to the opening (phase 1), but a coil has one side in each layer, so one wire size serves both: "
            f"the wire cannot differ by layer. The most it could save here is 90 % of the recommended winding's proximity loss, {saving:.2f} W of {d['hs_p_total']:.1f} W ({100 * saving / d['hs_p_total']:.2f} %) at "
            f"{HS_RPM:,.0f} rpm: not worth a mixed winding.")
    out["thin_near_opening_saving_w"] = saving
    return out


def write_json(path: Path, obj) -> None:
    def conv(v):
        if isinstance(v, (np.floating, float)):
            return None if not np.isfinite(v) else float(v)
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.bool_,)):
            return bool(v)
        if isinstance(v, np.ndarray):
            return [conv(x) for x in v.tolist()]
        if isinstance(v, dict):
            return {str(k): conv(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [conv(x) for x in v]
        if isinstance(v, (Case,)):
            return str(v)
        return v
    path.write_text(json.dumps(conv(obj), indent=1), encoding="utf-8")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="the winding optimizer: the front, the shortlist and the recommended winding of the optimized motor")
    ap.add_argument("--stage", default="all", choices=("headline", "variants", "all"))
    ap.add_argument("--workers", type=int, default=1, help="worker processes (below-normal priority, under the solver slot); 1 runs in this process")
    ap.add_argument("--no-figures", action="store_true")
    args = ap.parse_args(argv)
    out = OUT()
    out.mkdir(parents=True, exist_ok=True)
    rep = Report()
    t0 = time.perf_counter()
    h = stage_headline(rep, out, args.workers)
    results = {"headline": {"front": [winding_dict(h["run"].cand, int(j)) for j in h["run"].front()], "recommended": h["rec"], "requirement": h["req"], "shortlist": h["rows"]}}
    rep.say()
    rep.say("THE CAP ON STRANDS IN HAND AT 60 MM, AND WHAT EACH CAP COSTS")
    results["headline"]["cap"] = cap_table(rep, h["run"])
    rep.say()
    results["headline"]["envelopes"] = stage_envelopes(rep, h)
    rep.say()
    results["rules"] = rules_applied(rep, h)
    if args.stage in ("variants", "all"):
        rep.say()
        results["variants"] = stage_variants(rep, out, args.workers, h)
    rep.say()
    rep.say(f"done in {time.perf_counter() - t0:.1f} s")
    (out / "report.txt").write_text("\n".join(rep.lines) + "\n", encoding="utf-8")
    write_csv(out / "front_60mm.csv", [{"N": int(h["run"].cand["N"][j]), "a": int(h["run"].cand["a"][j]), "n_s": int(h["run"].cand["n_s"][j]), "d_mm": float(h["run"].cand["d_mm"][j]),
                                         "fill": float(h["run"].cand["fill"][j]), "I_st": float(h["run"].obj["I_st"][j]), "corner48": float(h["run"].obj["corner48"][j]),
                                         "P_hs": float(h["run"].obj["P_hs"][j])} for j in h["run"].front()])
    write_json(out / "results.json", results)
    if not args.no_figures:
        try:
            from optimization import winding_study_figures
            winding_study_figures.draw_all(h, results, out)
        except ImportError:
            print("(winding_study_figures.py is not there yet: no figures)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
