"""B7: Python -> SALib -> FEMM (femm_pool): Morris screening of four rotor dimensions of the band-free motor, on the coarse-mesh torque at 10 A and a fixed current angle.

What it shows, and why
    Before an optimizer spends hundreds of FEMM runs on a motor's dimensions, a screening asks which of them move the answer at all. Morris's method (SALib) does it cheaply:
    it walks the parameters one at a time along a few random trajectories, and for each parameter averages how much the output changes per step. B7 sends the motor through it:
    four rotor dimensions of motorgeo.Params (barrier="none", R11), the output the mean torque at 10 A and the MTPA current angle, FEMM 5.0 on the coarse mesh through femm_pool.
    The case is sound when two runs give identical indices (the same seed, the same sample, the same tool: nothing may differ in the last digit), when SALib's wiring reproduces
    a known answer, and when the ranking of mu* is the ranking of the motor's physics, argued below in one line each.

The four dimensions, and the ranges (a finding as much as a choice)
    PARAMETERS = slot_t, slot_w, mag_y, v_half_deg: the magnet slot's thickness (1.77 mm) and width (9.28 mm; the slot is filled with N45, so these are the magnet), how far the two
    magnets of a V sit from the d-axis (mag_y, 6 mm) and the V's half-angle (60 degrees). Why these: with barrier="none" the rotor's other fields draw nothing (rib, bridge, web,
    bar_t and bar_apex belong to the band or only guard the build; a test shows the drawing does not change when they do), and mag_x, the magnet's depth, has almost no room outward (+3 % is rejected)
    (see below). Each one changes the drawing, so each has an index that means something.
    The ranges are +-RANGE_PERCENT of the design value, the same for all four, so mu* is a like-for-like ranking (the change of torque for a +-1.5 % change of a dimension).
    Why so narrow: the design sits 0.158 mm from motorgeo's guard of the outer bridge (bridge = 1.4 mm; the slot's outer corner is at 1.558 mm from the rotor surface), and three of the
    four dimensions push the slot's outer end toward the surface. At the corners of a box, where they all push together, the room is used up. Every level of the Morris grid (4 levels
    in 4 parameters: 256 designs, 16 corners among them) is built and checked by motorgeo BEFORE any FEMM run, and the case tried a ladder of ranges until none was rejected: see
    `ranges_tried` in results.json (+-5 %: 48 of 256 rejected; +-3 %: 16; +-2 %: 1, a corner; +-1.5 %: 0). A rejected sample would break a Morris trajectory, so none is allowed.
    (This is a finding about the design, not a choice: a wider screening needs the guard relaxed. `bridge` only guards the build in the slot layout, and the drawn rotor is identical with
    bridge = 1.0 mm, motorgeo.check's own limit for the outer bridge; with it the same four parameters pass at +-5 % (0 of 256 rejected; +-7 %: 1). The case keeps R11's
    Params(barrier="none") as it is and does not use it; a later task can rule on relaxing the guard.)

Method (the files and the protocol)
    1. SALib 1.6.0 (Morris, read from the package metadata: it has no __version__): morris.sample(problem, r, num_levels=4, seed) gives r (k + 1) = 6 x 5 = 30 samples in r
       trajectories of k + 1 points, each one step of 2/3 of a range from the last in one parameter. morris.analyze(problem, X, Y, seed) gives mu (the mean elementary effect), mu*
       (the mean of its size), sigma (its spread) and the confidence interval of mu* by bootstrap (100 resamples, 95 %): the same seed in both, so the bootstrap repeats.
       The defaults and the definitions are SALib's own, read in its installed source (site-packages\\SALib\\, opened 2026-11-06; SALIB_SOURCES holds them and a test reads each line
       back): 4 levels (sample/morris/morris.py:27 and analyze/morris.py:24), 100 resamples (analyze/morris.py:20) at 95 % (analyze/morris.py:21), the step 2/3 of the unit range,
       num_levels / (2 (num_levels - 1)) (util/__init__.py:399), an elementary effect as the change of the output over that step (analyze/morris.py:415), mu as the mean of the
       elementary effects (analyze/morris.py:211), mu* as the mean of their absolute values (analyze/morris.py:212) and sigma as their standard deviation (analyze/morris.py:213).
    2. A sample is a design: femm_points.solve_points builds it (R11), has motorgeo.check refuse it if it is rejected, and has a FEMM of the pool write the model (coarse mesh) and
       solve it at two rotor positions, theta_e = 0 and 30 electrical degrees, with 10 A at GAMMA_DEG = 22.9 degrees (id = -10 sin gamma, iq = 10 cos gamma). The output of a sample is
       the MEAN of the two air-gap-element torques. Why two positions: the torque of this machine ripples with the rotor position (in an exploration of 11-06, coarse mesh, at the design
       and this current: the torque at theta_e = 0 was 9 % off the mean of six positions 10 degrees apart, B1's); the pair 30 electrical degrees apart cancels the 6th, 18th and 30th
       harmonics, among them the cogging torque's fundamental (period 20 degrees), and the mean of the pair was 0.27 % from the six-position mean there. So 30 samples x 2 positions = 60 FEMM solves per run, at most about 60 as ruled.
    3. Two runs, two femm_pool sessions, one after the other, each with its own call of sample() with the same seed. The samples, the torques and every index must be EXACTLY equal.
       If they are not, the case FAILS and says by how much; the check is never turned into a tolerance. (FEMM 5.0 gave identical numbers in B1 to 1e-13.)
    4. Cross-checks that need no other sampling: SALib on a LINEAR function through this very sample (every elementary effect is then the coefficient times the range, so mu* = |c| x
       range exactly: the wiring of columns and ranges); and the one-at-a-time swings, 9 designs (the centre and each end of each range) in the first session: for a monotonic response
       mu* is about |swing| (the change of torque across the whole range), so the two rankings must agree. The design's own mean torque is held against B1's MTPA torque on the fine mesh.
    5. The picture (the ruling): two panels. Left, mu* per parameter as horizontal bars, the largest on top, each with SALib's confidence interval as an error bar. Right, the Morris plane,
       sigma against mu*, each parameter a labelled point, with the line sigma = mu* drawn: below it a parameter acts mostly alone and linearly, above it nonlinearly or together
       with others. In N.m, the torque's units. matplotlib only.

Run (from the Python to FEMM folder; the solver slot is taken for two pool sessions of about 20 s each)
    uv run python bench/run_bench.py --only B7
Writes out\\bench\\B7\\: results.json (every number and every check), figure.png, fem\\ (the models, small), run.log, runs.json.

Numbers to expect (this PC, 2026-11-06, FEMM 5.0, SALib 1.6.0; four real runs gave the same numbers; the case takes 43 to 48 s: 35 s of it the screening of the ladder and the grid)
    ranges      +-5 %: 48 of 256 grid designs rejected, +-3 %: 16, +-2 %: 1 (a corner: slot_t, slot_w and mag_y high, v_half_deg low), +-1.5 %: 0 (chosen); 0 of the 30 samples rejected
    size        r (k + 1) = 6 x 5 = 30 samples, 60 solves a run (2.8 s on 8 workers), and 18 more for the one-at-a-time designs (1.3 s), in the first session only: 60 + 60 + 18 = 138
                solves in the case; two sessions of 8 workers; the models have 8,850 to 9,150 nodes
    indices     slot_w   mu* 0.1261 N.m (+-0.0017)  mu +0.1261  sigma 0.0024   one-at-a-time swing +0.1251 N.m   (psi_d +2.27 mWb across its range, psi_q -0.04 mWb)
                mag_y    mu* 0.0414 N.m (+-0.0009)  mu -0.0414  sigma 0.0011   swing -0.0421 N.m                  (psi_d -0.76 mWb, psi_q +0.03 mWb)
                slot_t   mu* 0.0186 N.m (+-0.0013)  mu +0.0186  sigma 0.0019   swing +0.0183 N.m                  (psi_d +0.34 mWb, psi_q 0.00 mWb)
                v_half_deg mu* 0.0074 N.m (+-0.0014)  mu -0.0074  sigma 0.0020   swing -0.0070 N.m                (psi_d -0.11 mWb, psi_q -0.03 mWb)
                sigma is 1.9 % (slot_w), 2.8 % (mag_y), 10 % (slot_t) and 26.5 % (v_half_deg) of its mu* (sigma_over_mu_star in results.json): all four sit below the line sigma = mu*, so
                each dimension acts mostly alone and in one direction; the smallest effect, v_half_deg, has the largest share
    two runs    the 30 torques, the sample and every index are identical to the last bit (FEMM 5.0 is deterministic here)
    wiring      SALib on a linear function through this sample: mu* equals |c| x range to 5e-15 of the largest
    the design  mean torque 1.12493 N.m (coarse mesh, 0 and 30 electrical degrees) against B1's 1.12914 N.m (fine mesh, six positions): -0.37 %; a single rotor position was 9 % off

Self-check (the answers follow).
  1. A parameter has mu* = 0.0075 N.m and sigma = 0.0006 N.m. Is it a parameter to fix at its design value?
  2. Why is the output a MEAN over two rotor positions, and not the torque at one?
  3. The two runs agree to 1e-12 but not exactly. Why is that a FAIL here and not a pass?
Answers.
  1. Probably, for this operating point and this range: its mean effect is the smallest of the four and it acts almost linearly (sigma is a tenth of mu*, below the line sigma = mu*).
     Fixing it costs a torque error of at most mu* when it moves across its whole range. Whether it matters for the cogging or the ripple is another output: B7 asks only the mean torque.
  2. The torque of an IPM machine ripples with the rotor position, and the ripple itself depends on the dimensions (the slot openings meet the magnets differently). At one position
     the index of a dimension is its effect on the mean torque AND on the ripple's phase; the mean over two positions 30 electrical degrees apart takes most of the ripple out.
  3. The brief asks for identical indices from a deterministic tool. A difference in the last digits means something in the chain is not deterministic (a thread, a file read in another
     order); that is worth knowing, and a tolerance would hide it. The case reports the difference, so the reader can judge it.
"""
from __future__ import annotations

import importlib.metadata
import itertools
import time
from pathlib import Path

import numpy as np

import femm_queue
import motorgeo as mg

from . import OWN_LIMIT, Log, check, compare, params, skipped, wrap_title
from . import femm_points
from .b1_femm import currents, pool_record
from .femm_points import POSITIONS_EL, TASK_TIMEOUT_S

CASE = "B7"
TITLE = "Python -> SALib -> FEMM (femm_pool): Morris screening of four rotor dimensions on the coarse-mesh torque"

# ---- the fixed choices (the brief and the rulings)
PARAMETERS = ("slot_t", "slot_w", "mag_y", "v_half_deg")          # the columns of the sample, in this order: k = 4
LEVELS = 4                                                         # Morris's grid: p = 4 levels (SALib's default; it wants an even number)
TRAJECTORIES = 6                                                   # r: 6 x (4 + 1) = 30 samples
SEED = 1                                                           # in sample() AND in analyze() (its bootstrap)
GAMMA_DEG = 22.9
GAMMA_SOURCE = ("B1's MTPA angle on the fine mesh: out\\bench\\B1\\results.json of 11-06, summary.none.mtpa.bench_method.gamma_deg = 22.936 degrees "
                "(its best 1-degree sample 23.0); a design decision types it as 22.9. At gamma = 0 only the magnet torque counts and the saliency parameters would be hidden")
AMPS = 10.0                                                        # the brief: 10 A (peak, 10 turns a coil side: 100 A-turns)
SOLVE_BUDGET = 60                                                  # the ruling: at most about 60 FEMM solves per run, on the coarse mesh (POSITIONS_EL and TASK_TIMEOUT_S: femm_points)
RANGE_PERCENT = 1.5                                                # +- of the design value, for every parameter: the widest rung of the ladder that motorgeo accepts
RANGE_LADDER_PERCENT = (5.0, 3.0, 2.0, 1.5)                        # the rungs tried, widest first
CONF_LEVEL, RESAMPLES = 0.95, 100                                  # the bootstrap of mu*'s confidence interval (SALib's default resamples)
# Where SALib's defaults and definitions come from: SALib 1.6.0's own source, as installed in this project's environment (<project>\.venv\Lib\site-packages\SALib\), opened on
# 2026-11-06: (what, the file under SALib\, the line, the text on that line). tests/test_bench_b7.py reads each line back, so a SALib that moved one or changed a default shows.
SALIB_SOURCES = (
    ("4 levels (the sampler's default)", "sample/morris/morris.py", 27, "num_levels: int = 4"),
    ("4 levels (the analysis's default)", "analyze/morris.py", 24, "num_levels: int = 4"),
    ("100 resamples (the bootstrap)", "analyze/morris.py", 20, "num_resamples: int = 100"),
    ("95 % confidence", "analyze/morris.py", 21, "conf_level: float = 0.95"),
    ("the step 2/3 of the unit range", "util/__init__.py", 399, "return num_levels / (2.0 * (num_levels - 1))"),
    ("an elementary effect is the change of the output over that step", "analyze/morris.py", 415, "elementary_effects = difference / delta"),
    ("mu is the mean of the elementary effects", "analyze/morris.py", 211, "mu = np.average(elementary_effects, 1)"),
    ("mu* is the mean of their absolute values", "analyze/morris.py", 212, "mu_star = np.average(np.abs(elementary_effects), 1)"),
    ("sigma is their standard deviation (n - 1)", "analyze/morris.py", 213, "sigma = np.std(elementary_effects, axis=1, ddof=1)"),
)
B1_TOLERANCE_PERCENT = 2.0                                         # the design's coarse two-position torque against B1's fine six-position MTPA torque: the bench's own limit
KNOWN_COEFFICIENTS = (4.0, -2.5, 0.0, 10.0)                        # the linear test function of the SALib wiring check (N.m per unit of each parameter)

PHYSICS = {
    "slot_w": "the magnet's length (the slot is filled with N45): it moves the d-axis flux psi_d by much more than its share of magnet area, since the slot's outer end also nears the "
              "bridge that leaks the magnet's flux (the two are not separated here); psi_q does not feel it",
    "mag_y": "how far the two magnets of a V sit from the d-axis: moved apart they widen the V's centre rib, which would leak more of their flux, so psi_d and the torque fall (the sign "
             "of mu says so; the mechanism is argued, not isolated); psi_q does not feel it",
    "slot_t": "the magnet's thickness: more of it is more magnetising force in series with the air gap, so psi_d and the torque rise; it is also a wider barrier to the d-axis flux "
              "(not isolated), but psi_q hardly moves",
    "v_half_deg": "the V's half-angle: turning the magnets away from the d-axis (a flatter V) lets them push a little less flux through the pole, a small loss of psi_d (argued, not isolated)",
}


# ============================================================================================== the design, the ranges, the grid (no solver)
def samples() -> int:
    """r (k + 1): the samples of one Morris run."""
    return TRAJECTORIES * (len(PARAMETERS) + 1)


def solves_per_run() -> int:
    """The FEMM solves of one Morris run: every sample at every rotor position."""
    return samples() * len(POSITIONS_EL)


def oat_solves() -> int:
    """The FEMM solves of the one-at-a-time designs (the centre and each end of each range, in the first session only): every design at every rotor position."""
    return (1 + 2 * len(PARAMETERS)) * len(POSITIONS_EL)


def solves_in_the_case() -> int:
    """All the FEMM solves of the case: the Morris solves of each of the two sessions, and the one-at-a-time ones of the first."""
    return 2 * solves_per_run() + oat_solves()


def bounds(percent: float = RANGE_PERCENT) -> list[list[float]]:
    """[low, high] of each parameter: the design value less and plus `percent` per cent of it."""
    base = mg.Params(barrier="none")
    return [[getattr(base, name) * (1 - percent / 100), getattr(base, name) * (1 + percent / 100)] for name in PARAMETERS]


def problem(percent: float = RANGE_PERCENT) -> dict:
    """SALib's problem definition."""
    return {"num_vars": len(PARAMETERS), "names": list(PARAMETERS), "bounds": bounds(percent)}


def grid_overrides(percent: float = RANGE_PERCENT, corners_only: bool = False) -> list[dict]:
    """Every design of the Morris grid: the LEVELS levels of each parameter (or only its two ends, the corners of the box) as {parameter: value}, 256 of them for 4 levels in 4."""
    levels = np.linspace(-1.0, 1.0, LEVELS)
    choices = levels[[0, -1]] if corners_only else levels
    base = mg.Params(barrier="none")
    return [{name: getattr(base, name) * (1 + percent / 100 * c) for name, c in zip(PARAMETERS, combo)} for combo in itertools.product(choices, repeat=len(PARAMETERS))]


def rejection(overrides: dict) -> str | None:
    """Why motorgeo rejects the design Params(barrier="none", **overrides): None if it accepts it, else the reason (build() refuses it, or check() reports a violation)."""
    try:
        p = mg.Params(barrier="none", **overrides)
        _, violations = mg.check(p, mg.build(p))
    except mg.GeometryError as error:
        return f"GeometryError: {error}"
    return "; ".join(violations) if violations else None


def screen(designs: list[dict]) -> dict:
    """motorgeo's verdict on each design: {"points", "count" (rejected), "rejected": [(overrides, reason)]}."""
    rejected = [(d, reason) for d in designs if (reason := rejection(d)) is not None]
    return {"points": len(designs), "count": len(rejected), "rejected": rejected}


def ranges_tried(ladder=RANGE_LADDER_PERCENT) -> list[dict]:
    """The ladder of ranges the case tried, each on the whole grid: how many of its designs motorgeo rejects. {"percent", "points", "rejected", "first_reason"}."""
    out = []
    for percent in ladder:
        screened = screen(grid_overrides(percent))
        out.append({"percent": percent, "points": screened["points"], "rejected": screened["count"],
                    "first_reason": screened["rejected"][0][1] if screened["rejected"] else None})
    return out


# ============================================================================================== the sample and the analysis (SALib)
def sample(seed: int = SEED, trajectories: int = TRAJECTORIES, percent: float = RANGE_PERCENT) -> np.ndarray:
    """SALib's Morris sample: (trajectories x (k + 1), k), one column per parameter in PARAMETERS' order, in the parameters' own units."""
    from SALib.sample import morris as morris_sample
    return morris_sample.sample(problem(percent), trajectories, num_levels=LEVELS, seed=seed)


def overrides_of(row) -> dict:
    """A row of the sample as the overrides of motorgeo.Params: {parameter: value}, in the order of PARAMETERS."""
    return {name: float(value) for name, value in zip(PARAMETERS, row)}


def analyse(X, Y, seed: int = SEED) -> dict:
    """SALib's Morris analysis as plain lists: {"names", "mu", "mu_star", "sigma", "mu_star_conf"} (the torque's units; mu_star_conf is the 95 % confidence interval of mu*)."""
    from SALib.analyze import morris as morris_analyze
    result = morris_analyze.analyze(problem(), np.asarray(X, float), np.asarray(Y, float), conf_level=CONF_LEVEL, num_resamples=RESAMPLES, print_to_console=False,
                                    num_levels=LEVELS, seed=seed)
    plain = {key: [float(v) for v in np.ma.filled(np.ma.asarray(result[key], dtype=float), np.nan)] for key in ("mu", "mu_star", "sigma", "mu_star_conf")}
    return {"names": [str(n) for n in result["names"]], **plain}


def known_answer(X) -> dict:
    """SALib on a linear function through this very sample: every elementary effect of c.x is c_i x range_i, so mu*_i = |c_i| range_i, mu_i = c_i range_i, sigma = 0.
    {"coefficients", "expected_mu_star", "mu_star", "max_relative_error", "sigma_max", "ok"}; the error is scaled by the largest mu* (the one parameter with c = 0 has no scale)."""
    coef = np.array(KNOWN_COEFFICIENTS, float)
    ranges = np.array([high - low for low, high in problem()["bounds"]])
    result = analyse(X, np.asarray(X, float) @ coef)
    expected_star, expected_mu = np.abs(coef) * ranges, coef * ranges
    scale = float(expected_star.max())
    error = max(float(np.max(np.abs(np.array(result["mu_star"]) - expected_star))), float(np.max(np.abs(np.array(result["mu"]) - expected_mu)))) / scale
    sigma_max = float(np.max(np.abs(result["sigma"]))) / scale
    return {"coefficients": list(KNOWN_COEFFICIENTS), "expected_mu_star": expected_star.tolist(), "mu_star": result["mu_star"], "max_relative_error": error,
            "sigma_max": sigma_max, "ok": bool(error < 1e-9 and sigma_max < 1e-9 and result["names"] == list(PARAMETERS))}


def exactly_equal(a, b) -> bool:
    """a and b are the same to the last bit: arrays, lists and numbers by value, dicts key by key (the same keys), text by text. NaN is not equal to anything: a NaN is a failure."""
    if isinstance(a, dict) or isinstance(b, dict):
        return isinstance(a, dict) and isinstance(b, dict) and set(a) == set(b) and all(exactly_equal(a[k], b[k]) for k in a)
    if isinstance(a, str) or isinstance(b, str):
        return a == b
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)) and any(isinstance(v, str) for v in (*a, *b)):
        return list(a) == list(b)
    first, second = np.asarray(a, float), np.asarray(b, float)
    return first.shape == second.shape and bool(np.array_equal(first, second))


def largest_difference(a, b) -> float:
    """The largest absolute difference between two runs' numbers (arrays, or dicts of lists; the names are skipped); infinity if their shapes differ. 0.0 when they are equal."""
    if isinstance(a, dict):
        differences = [largest_difference(a[k], b[k]) for k in a if k in b and k != "names"]
        return max(differences) if differences else 0.0
    first, second = np.asarray(a, float), np.asarray(b, float)
    if first.shape != second.shape:
        return float("inf")
    return float(np.max(np.abs(first - second))) if first.size else 0.0


def ranking(analysis: dict) -> list[str]:
    """The parameters' names, the one with the largest mu* first (equal values keep their order)."""
    order = np.argsort(-np.asarray(analysis["mu_star"], float), kind="stable")
    return [analysis["names"][i] for i in order]


# ============================================================================================== the one-at-a-time designs and their swings
def oat_designs(percent: float = RANGE_PERCENT) -> list[dict]:
    """The design itself and each parameter alone at each end of its range: [{"name", "param", "side", "overrides"}], 1 + 2 k = 9 of them."""
    base = mg.Params(barrier="none")
    designs = [{"name": "centre", "param": None, "side": None, "overrides": {}}]
    for name in PARAMETERS:
        for side, factor in (("low", 1 - percent / 100), ("high", 1 + percent / 100)):
            designs.append({"name": f"{name} {side}", "param": name, "side": side, "overrides": {name: getattr(base, name) * factor}})
    return designs


def swings(torque_by_design: dict) -> dict:
    """{parameter: the torque at the high end of its range less the torque at the low end} (signed, N.m), from {design name: torque}: the change across the whole range."""
    return {name: torque_by_design[f"{name} high"] - torque_by_design[f"{name} low"] for name in PARAMETERS}


def swing_ranking(swing: dict) -> list[str]:
    """The parameters by the SIZE of their swing, the largest first."""
    return sorted(swing, key=lambda name: -abs(swing[name]))


def ranking_lines(analysis: dict, swing: dict, flux: dict | None = None) -> list[str]:
    """The printed ranking, one line for each parameter, the largest mu* first: its mu*, the confidence interval, mu (signed), sigma, the one-at-a-time swing and the physics.
    `flux` ({parameter: (change of psi_d, change of psi_q across the range, in mWb)}) says where the effect enters."""
    lines = []
    for place, name in enumerate(ranking(analysis), start=1):
        i = analysis["names"].index(name)
        where = ""
        if flux and name in flux:
            where = f" [across the range psi_d {flux[name][0]:+.3f} mWb, psi_q {flux[name][1]:+.3f} mWb]"
        lines.append(f"{place}. {name}: mu* {analysis['mu_star'][i]:.4f} N.m (+-{analysis['mu_star_conf'][i]:.4f}), mu {analysis['mu'][i]:.4f}, sigma {analysis['sigma'][i]:.4f}; "
                     f"one-at-a-time swing {swing[name]:+.4f} N.m{where} - {PHYSICS[name]}")
    return lines


# ============================================================================================== the FEMM side
def current_point(theta_e: float) -> dict:
    """The point to solve at the rotor position theta_e: 10 A at GAMMA_DEG (b1_femm.currents: id = -I sin(gamma), iq = I cos(gamma))."""
    i_d, i_q = currents(GAMMA_DEG, AMPS)
    return {"theta_e": float(theta_e), "id": i_d, "iq": i_q}


def make_task(task_id: int, tag: str, overrides: dict, fem_dir) -> dict:
    """A task for femm_points.solve_points: this design at every rotor position of POSITIONS_EL."""
    return {"id": task_id, "tag": tag, "overrides": overrides, "fem_dir": str(fem_dir), "points": [current_point(t) for t in POSITIONS_EL]}


def mean_torque(result: dict) -> float:
    """The output of a sample: the mean over its rotor positions of the torque from the air-gap element [N.m]."""
    return float(np.mean([row["T_age"] for row in result["rows"]]))


def mean_flux(result: dict, key: str) -> float:
    """The mean over the rotor positions of psi_d or psi_q [Wb]."""
    return float(np.mean([row[key] for row in result["rows"]]))


def run_session(number: int, X, designs: list[dict], fem_dir, log: Log) -> dict:
    """One femm_pool session: the Morris sample's tasks (and, when `designs` are given, the one-at-a-time ones) solved on up to 8 FEMMs. Both runs come through here."""
    label = f"bench B7 FEMM (run {number})"
    tasks = [make_task(i, f"run{number}_s{i:03d}", overrides_of(row), fem_dir) for i, row in enumerate(X)]
    extra = [make_task(100 + i, f"run{number}_{d['name'].replace(' ', '_')}", d["overrides"], fem_dir) for i, d in enumerate(designs)]
    started = time.perf_counter()
    with log.tee():
        with femm_queue.femm_pool(label, task_timeout_s=TASK_TIMEOUT_S) as pool:
            sample_results = femm_points.run_tasks(pool, tasks, f"B7 run {number}: the Morris sample", log)
            oat_results = femm_points.run_tasks(pool, extra, f"B7 run {number}: the one-at-a-time designs", log) if extra else []
            record = pool_record(label, pool, sample_results + oat_results, time.perf_counter() - started, task_timeout_s=TASK_TIMEOUT_S)
    log.add_run(record)
    return {"samples": sample_results, "oat": oat_results, "Y": np.array([mean_torque(r) for r in sample_results]), "seconds": time.perf_counter() - started, "record": record}


# ============================================================================================== the picture
def draw_morris(analysis: dict, title: str):
    """The two Morris views as a matplotlib figure (not saved): left, mu* per parameter as horizontal bars, the largest on top, with SALib's confidence interval of mu* as an error
    bar; right, sigma against mu* with a labelled point for each parameter and the line sigma = mu*. In N.m."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    names = list(analysis["names"])
    mu_star, sigma, conf = (np.asarray(analysis[k], float) for k in ("mu_star", "sigma", "mu_star_conf"))
    order = np.argsort(-mu_star, kind="stable")
    text = wrap_title(title)
    lines = text.count("\n") + 1
    height = 4.6 + 0.2 * lines
    fig = plt.figure(figsize=(11, height))
    grid = fig.add_gridspec(1, 2, left=0.115, right=0.985, bottom=0.135, top=1 - (0.2 * lines + 0.3) / height, wspace=0.27)
    left, right = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
    fig.text(0.005, 1 - 0.05 / height, text, ha="left", va="top", fontsize=10.5, color="#1d1d1f")
    # left: mu* with its confidence interval, the largest on top
    y = np.arange(len(names))[::-1]
    left.barh(y, mu_star[order], xerr=conf[order], height=0.62, color="#2f6fb5", ecolor="#222222", capsize=4)
    left.set_yticks(y)
    left.set_yticklabels([names[i] for i in order])
    for row, i in zip(y, order):
        left.text(mu_star[i] + conf[i] + 0.012 * float(mu_star.max()), row, f"{mu_star[i]:.4f}", va="center", fontsize=8.5, color="#333333")
    left.set_xlim(0, 1.22 * float((mu_star + conf).max()))
    left.set_xlabel("mu*  [N.m]   (mean size of the change of torque across a parameter's whole range)")
    left.set_title("how much each dimension moves the torque (bars: 95 % confidence interval)", fontsize=9.5, loc="left")
    left.grid(axis="x", alpha=0.3)
    # right: the Morris plane, sigma against mu*, with the line sigma = mu*
    top = 1.18 * float(max(mu_star.max(), sigma.max()))
    right.plot([0, top], [0, top], "--", color="#888888", lw=1.2)
    right.scatter(mu_star, sigma, s=46, color="#c0392b", zorder=3)
    ascending = [int(i) for i in np.argsort(mu_star, kind="stable")]
    groups: list[list[int]] = []                                              # points whose labels would overlap (neighbours closer than 11 % of the axis) form one group
    for i in ascending:
        if groups and mu_star[i] - mu_star[groups[-1][-1]] < 0.11 * top:
            groups[-1].append(i)
        else:
            groups.append([i])
    for group in groups:                                                      # a label a little above its point; in a group the lowest point's label is the highest, on a leader line
        for place, i in enumerate(group):
            lift = len(group) - 1 - place
            right.annotate(names[i], (mu_star[i], sigma[i]), xytext=(4, 16 + 15 * lift), textcoords="offset points", fontsize=9.5, ha="left",
                           arrowprops={"arrowstyle": "-", "color": "#999999", "lw": 0.8} if lift else None)
    right.set(xlim=(0, top), ylim=(0, top), xlabel="mu*  [N.m]", ylabel="sigma  [N.m]")
    right.text(0.97, 0.30, "below the line sigma = mu*:\nacts mostly alone, linearly", transform=right.transAxes, ha="right", va="center", fontsize=8.5, color="#555555")
    right.text(0.04, 0.96, "above it: nonlinear,\nor together with others", transform=right.transAxes, ha="left", va="top", fontsize=8.5, color="#555555")
    right.set_title("the Morris plane: sigma against mu*", fontsize=9.5, loc="left")
    right.grid(alpha=0.3)
    return fig


def fig_morris(analysis: dict, title: str, target) -> None:
    """draw_morris, saved with a tight box (so a long title cannot leave the picture blank at the sides)."""
    import matplotlib.pyplot as plt
    fig = draw_morris(analysis, title)
    fig.savefig(target, dpi=130, bbox_inches="tight", pad_inches=0.12)
    plt.close(fig)


def figure_title(version: str, percent: float = RANGE_PERCENT) -> str:
    """The title of figure.png: what was varied, the output, the number of solves, the levels and the seed (the ruling). Each line is at most bench.TITLE_CHARS long, so that
    wrap_title leaves it as written."""
    return (f"Morris screening (SALib {version}) of four rotor dimensions of the band-free motor:\n"
            f"{', '.join(PARAMETERS)}, each +-{percent:g} % of its design value\n"
            f"output: mean torque [N.m] at {AMPS:g} A, gamma = {GAMMA_DEG:g} deg, rotor at {' and '.join(f'{t:g}' for t in POSITIONS_EL)} electrical degrees;\n"
            f"FEMM 5.0, coarse mesh; r x (k + 1) = {TRAJECTORIES} x {len(PARAMETERS) + 1} = {samples()} samples = {solves_per_run()} solves a run; "
            f"{LEVELS} levels, seed {SEED}")


def sigma_over_mu_star(analysis: dict) -> dict[str, float]:
    """{parameter: sigma / mu*} in the order of the analysis: how large the spread of a parameter's elementary effects is against their mean size (0 for a linear response)."""
    return {name: float(s) / float(m) for name, s, m in zip(analysis["names"], analysis["sigma"], analysis["mu_star"])}


# ============================================================================================== the case
def run(folder: Path, log: Log) -> dict:
    """Run B7 in `folder` (already empty). Returns the results dict, whose "checks" list every verdict."""
    folder = Path(folder)
    started = time.perf_counter()
    version = importlib.metadata.version("SALib")
    # ---- no solver: the ranges, every level of the grid, the sample
    t = time.perf_counter()
    tried = ranges_tried()
    grid = grid_overrides()
    grid_screen = screen(grid)
    X1 = sample()
    sample_screen = screen([overrides_of(row) for row in X1])
    log.say(f"B7: ranges tried (rejected of {grid_screen['points']} grid designs): " + "; ".join(f"+-{r['percent']:g} %: {r['rejected']}" for r in tried)
            + f"; chosen +-{RANGE_PERCENT:g} %: {grid_screen['count']} rejected, {sample_screen['count']} of {len(X1)} samples rejected ({time.perf_counter() - t:.1f} s)")
    designs = oat_designs()
    # ---- two femm_pool sessions, each with its own call of sample() with the same seed
    first = run_session(1, X1, designs, folder / "fem" / "run1", log)
    X2 = sample()
    second = run_session(2, X2, [], folder / "fem" / "run2", log)
    analysis1, analysis2 = analyse(X1, first["Y"]), analyse(X2, second["Y"])
    torque_by_design = {d["name"]: mean_torque(r) for d, r in zip(designs, first["oat"])}
    swing = swings(torque_by_design)
    by_design = {d["name"]: result for d, result in zip(designs, first["oat"])}
    flux = {name: tuple(1e3 * (mean_flux(by_design[f"{name} high"], key) - mean_flux(by_design[f"{name} low"], key)) for key in ("psi_d", "psi_q")) for name in PARAMETERS}
    known = known_answer(X1)
    lines = ranking_lines(analysis1, swing, flux)
    for line in lines:
        log.say("B7 " + line)
    figure = folder / "figure.png"
    fig_morris(analysis1, figure_title(version), figure)
    nodes = [r["nodes"] for r in first["samples"] + first["oat"] if r["nodes"]]
    centre = torque_by_design["centre"]
    b1, b1_reason = femm_points.b1_reference()
    solves = {"run1": sum(len(r["rows"]) for r in first["samples"]), "run2": sum(len(r["rows"]) for r in second["samples"]), "oat": sum(len(r["rows"]) for r in first["oat"])}
    context = {"X1": X1, "X2": X2, "Y1": first["Y"], "Y2": second["Y"], "analysis1": analysis1, "analysis2": analysis2, "rejected_grid": grid_screen,
               "grid_points": grid_screen["points"], "rejected_samples": sample_screen, "solves": solves["run1"], "solves_run2": solves["run2"], "oat_solves": solves["oat"],
               "swings": swing, "known": known, "figure": figure, "centre_torque": centre, "b1": b1, "b1_reason": b1_reason, "mesh_nodes": (min(nodes), max(nodes)) if nodes else (0, 0)}
    results = {"case": CASE, "title": TITLE, "salib_version": version, "params": params(mg.Params(barrier="none")), "parameters": list(PARAMETERS),
               "range_percent": RANGE_PERCENT, "bounds": bounds(), "design_values": {n: getattr(mg.Params(barrier="none"), n) for n in PARAMETERS},
               "ranges_tried": tried, "grid": {"points": grid_screen["points"], "rejected": grid_screen["count"], "corners": 2 ** len(PARAMETERS)},
               "operating_point": {"gamma_deg": GAMMA_DEG, "gamma_source": GAMMA_SOURCE, "amps": AMPS, "positions_el_deg": list(POSITIONS_EL)},
               "morris": {"levels": LEVELS, "trajectories": TRAJECTORIES, "samples": samples(), "solves_per_run": solves_per_run(), "seed": SEED, "resamples": RESAMPLES,
                          "conf_level": CONF_LEVEL, "X": X1, "Y_run1": first["Y"], "Y_run2": second["Y"], "analysis_run1": analysis1, "analysis_run2": analysis2,
                          "ranking": ranking(analysis1), "largest_difference_Y": largest_difference(first["Y"], second["Y"]),
                          "sigma_over_mu_star": sigma_over_mu_star(analysis1), "salib_sources": [list(s) for s in SALIB_SOURCES]},
               "solves": {"morris_run1": solves["run1"], "morris_run2": solves["run2"], "one_at_a_time": solves["oat"], "in_the_case": sum(solves.values())},
               "one_at_a_time": {"torque": torque_by_design, "swing": swing, "swing_ranking": swing_ranking(swing), "psi_d_psi_q_swing_mWb": flux},
               "known_answer": known, "physics": PHYSICS, "ranking_lines": lines, "mesh": {"nodes_min_max": list(context["mesh_nodes"])}, "b1_reference": context["b1"],
               "b1_reference_missing": context["b1_reason"],
               "pools": [first["record"], second["record"]], "seconds": {"run1": round(first["seconds"], 1), "run2": round(second["seconds"], 1)}}
    results["checks"] = make_checks(context)
    results["seconds_total"] = round(time.perf_counter() - started, 1)
    return results


def make_checks(context: dict) -> list[dict]:
    """Every number of B7 against its expectation, with the source of the expectation (see good_context in tests/test_bench_b7.py for what `context` holds)."""
    out: list[dict] = []
    grid, samples_screen = context["rejected_grid"], context["rejected_samples"]
    out.append(check("ranges: every level of the Morris grid passes motorgeo's checks (build and check)", grid["count"] == 0,
                     f"0 of {context['grid_points']} designs rejected: all {LEVELS} levels of the {len(PARAMETERS)} parameters, the {2 ** len(PARAMETERS)} corners among them, checked before any FEMM run "
                     f"(the rule for B7: the ranges are chosen so that no design of the grid is rejected, and the whole grid is checked before any FEMM run)",
                     f"{grid['count']} of {grid['points']} rejected" + (f"; first: {grid['rejected'][0][1]}" if grid["rejected"] else f" at +-{RANGE_PERCENT:g} %"),
                     "motorgeo.build and motorgeo.check on every grid design"))
    out.append(check("rejected samples: none of the Morris sample is rejected by motorgeo", samples_screen["count"] == 0,
                     f"0 rejected (the brief: 'report rejected samples'; a rejected sample breaks a trajectory, so the ranges were narrowed until there was none)",
                     f"{samples_screen['count']} of {samples_screen['points']} rejected", "motorgeo.build and motorgeo.check on each row of the sample"))
    total = context["solves"] + context["solves_run2"] + context["oat_solves"]
    out.append(check("solves: the Morris run is r x (k + 1) x the rotor positions, at most about 60 solves a session; the case also solves the one-at-a-time designs",
                     context["solves"] == solves_per_run() <= SOLVE_BUDGET and context["solves_run2"] == solves_per_run() and context["oat_solves"] == oat_solves(),
                     f"{TRAJECTORIES} x {len(PARAMETERS) + 1} x {len(POSITIONS_EL)} = {solves_per_run()} Morris solves a session, not more than {SOLVE_BUDGET} (the ruling), in each of the two "
                     f"sessions; and, in the first only, {1 + 2 * len(PARAMETERS)} one-at-a-time designs x {len(POSITIONS_EL)} positions = {oat_solves()} more (the physics cross-check of the "
                     f"ranking, unrequested and disclosed): {solves_in_the_case()} solves in the case",
                     f"{context['solves']} Morris solves in session 1, {context['solves_run2']} in session 2, plus {context['oat_solves']} one-at-a-time in session 1: {total} in the case",
                     "the rows the pool returned"))
    out.append(check("SALib on a known answer: a linear function through this very sample gives mu* = |c| x range", context["known"]["ok"],
                     f"mu* and mu equal |c| x range and c x range to 1e-9 of the largest, sigma 0, the names in the order of the sample's columns (Morris's definition: every elementary "
                     f"effect of a linear function is c_i x range_i; {OWN_LIMIT})",
                     f"largest error {context['known']['max_relative_error']:.2e}, largest sigma {context['known']['sigma_max']:.2e} (of the largest mu*)", "SALib 1.6.0 morris.analyze"))
    out.append(check("two runs, the same seed: the samples are identical", exactly_equal(context["X1"], context["X2"]),
                     "morris.sample called twice with the same seed gives the same array to the last bit", "identical" if exactly_equal(context["X1"], context["X2"])
                     else f"DIFFERENT: the largest difference is {largest_difference(context['X1'], context['X2']):.3e}", "SALib 1.6.0 morris.sample, seed " + str(SEED)))
    same_y = exactly_equal(context["Y1"], context["Y2"])
    out.append(check("two runs, the same seed: the torques are identical", same_y,
                     f"the {len(context['Y1'])} mean torques of the two femm_pool sessions are equal to the last bit (FEMM 5.0 was deterministic in B1 to 1e-13; the brief: 'identical')",
                     "identical" if same_y else f"DIFFERENT: the largest difference is {largest_difference(context['Y1'], context['Y2']):.3e} N.m", "FEMM 5.0 through femm_pool, twice"))
    same_a = exactly_equal(context["analysis1"], context["analysis2"])
    out.append(check("two runs, the same seed: the indices are identical (mu, mu*, sigma and the confidence intervals)", same_a,
                     "the analysis of the second run equals the first's in every number (the seed is in analyze() too, so its bootstrap repeats); the brief: 'identical indices'",
                     "identical" if same_a else f"DIFFERENT: the largest difference is {largest_difference(context['analysis1'], context['analysis2']):.3e} N.m", "SALib 1.6.0 morris.analyze"))
    y = np.asarray(context["Y1"], float)
    out.append(check("FEMM answered every sample with a finite, positive torque", bool(np.isfinite(y).all() and (y > 0).all()),
                     "the motor makes torque at 10 A and the MTPA angle in every design of the box", f"{len(y)} torques between {y.min():.4f} and {y.max():.4f} N.m",
                     "FEMM 5.0 through femm_pool"))
    by_mu, by_swing = ranking(context["analysis1"]), swing_ranking(context["swings"])
    out.append(check("the mu* ranking is the ranking of the one-at-a-time swings from the design", by_mu == by_swing,
                     f"mu* (the mean size of a step's effect over {TRAJECTORIES} random trajectories) orders the parameters as the change of torque across each range does at the design "
                     f"(a monotonic response gives mu* about |swing|) ({OWN_LIMIT})",
                     f"mu*: {', '.join(by_mu)}; swings: {', '.join(by_swing)} (" + ", ".join(f"{n} {context['swings'][n]:+.4f}" for n in by_swing) + " N.m)",
                     "9 more designs in the first session: the centre and each end of each range"))
    reference = context["b1"]
    if reference is None:
        out.append(skipped("the design's mean torque against B1's MTPA torque (fine mesh)", f"within {B1_TOLERANCE_PERCENT:g} % ({OWN_LIMIT})", context["b1_reason"], "B1"))
    else:
        out.append(compare("the design's mean torque (coarse mesh, two rotor positions) against B1's MTPA torque (fine mesh, six positions)", context["centre_torque"],
                           reference["best_sample_T_Nm"], B1_TOLERANCE_PERCENT, f"{reference['path']}: summary.none.mtpa.bench_method.best_sample_T_Nm at "
                           f"{reference['best_sample_gamma_deg']:g} degrees; {B1_TOLERANCE_PERCENT:g} % is {OWN_LIMIT} (a coarse mesh and two positions against a fine mesh and six)", "N.m"))
    figure = context["figure"]
    saved = figure is not None and Path(figure).is_file() and Path(figure).stat().st_size > 20_000
    out.append(check("the picture: figure.png is saved", bool(saved), "out\\bench\\B7\\figure.png: the two Morris views (mu* with its confidence interval, sigma against mu*)",
                     f"saved: {bool(saved)}" + (f", {Path(figure).stat().st_size:,} bytes" if saved else ""), "bench/b7_salib.py"))
    return out
