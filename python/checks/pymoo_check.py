"""pymoo_check: NSGA-II finds the Pareto front of ZDT1, and the distance to the true front (IGD) is measured.

What it does, and why
    The design study is a multi-objective optimization: more torque costs more ripple, so the answer is not one design but a
    front of best trade-offs (a Pareto front). pymoo's NSGA-II breeds designs generation after generation and keeps the
    non-dominated ones. ZDT1 is the standard test problem for that: two objectives, 30 variables in [0, 1],
        f1 = x1,    f2 = g (1 - sqrt(f1 / g)),    g = 1 + 9 (x2 + ... + x30) / 29,
    whose true front is known exactly: g = 1, that is f2 = 1 - sqrt(f1), from (0, 1) to (1, 0).
    The check measures how close the optimizer's front gets to it with the IGD (inverted generational distance): take 100 points
    spread over the TRUE front and, for each, the distance to the nearest point the optimizer found; the IGD is their mean. It is
    small only if the found points are on the front AND spread along all of it.

Method
    NSGA2(pop_size=40), 250 generations, seed 1 (so that the number is the same on every run), pymoo's own ZDT1 and its own
    IGD and true front (problem.pareto_front(): 100 points). It runs in-process: 10,000 evaluations of a formula, under a second.
    Numbers to expect (2026-11-06, pymoo 0.6.2): IGD 0.0134 for seed 1, 40 non-dominated points, 0.4 s; the limit is 0.05.

The ruling that sets the generations at 250 (2026-11-06)
    The brief said "NSGA-II on ZDT1, 40 x 100 ... the IGD is below 0.05". Those two numbers contradict each other. With a
    population of 40 and 100 generations the IGD is 0.0688 for seed 1, and over seeds 1 to 100 it lies between 0.0514 and 0.1846
    (median 0.0837): not one seed is below 0.05. The reason is that 40 individuals in 100 generations have not yet converged on a
    30-variable problem (the mean g is 1.14; the front has g = 1) and 40 points cannot cover 100 reference points. The evidence
    (pymoo 0.6.2, IGD against problem.pareto_front(), 100 seeds each):
        population  40 x 100 generations   IGD 0.0514 to 0.1846 over 100 seeds   the brief: 0 of 100 seeds below 0.05
        population 100 x 100 generations   IGD 0.0127 to 0.0294 over 100 seeds   100 of 100 below 0.05
        population  40 x 250 generations   IGD 0.0115 to 0.0153 over 100 seeds   100 of 100 below 0.05 (the ruling)
    A design decision of 2026-11-06 set 40 x 250 and left the limit at 0.05: more budget keeps the meaning of the limit, while a
    looser one would pass a run that has not converged. 40 x 250 and 100 x 100 cost the same 10,000 evaluations; 40 x 250 keeps
    the brief's population and has the wider margin (the worst of 100 seeds 0.0153, against 0.0294).

Self-check (the answers follow).
  1. The IGD of a front that is exactly on the true front but has only one point, at (0.5, 0.29), would be what, roughly?
  2. Why is a fixed seed used when the result of an optimizer is random?
  3. Would the generational distance (the mean distance of the FOUND points to the front) see what the IGD sees?
Answers.
  1. About 0.32 (measured 0.3236): the reference points near the two ends of the true front, (0, 1) and (1, 0), are 0.87 and 0.58
     away from that point, and the mean over all 100 reference points is large although the one point is exactly on the front.
     A perfect but short front is penalised.
  2. So that the row says the same thing every time it is run. The seed does not make the case easier: the 100-seed spread above
     is the honest statement of what NSGA-II gives for these settings.
  3. Only half of it: GD says whether the found points are near the front, not whether they cover all of it. Forty points huddled
     near (0.5, 0.29) have a perfect GD and a bad IGD.
"""
from __future__ import annotations

import importlib.metadata
import json
import time
from pathlib import Path

from . import FAIL, PASS, Row, fresh

POPULATION = 40
GENERATIONS = 250                       # the brief said 100; a design decision of 2026-11-06 set 250 (see the docstring: 100 never reaches LIMIT)
SEED = 1
LIMIT = 0.05


def run(folder: Path):
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.indicators.igd import IGD
    from pymoo.optimize import minimize
    from pymoo.problems import get_problem

    fresh(folder)
    problem = get_problem("zdt1")
    front = problem.pareto_front()
    case = f"NSGA-II on ZDT1, {POPULATION} individuals x {GENERATIONS} generations, seed {SEED}"
    expected = f"IGD to the true front below {LIMIT}"
    started = time.perf_counter()
    result = minimize(problem, NSGA2(pop_size=POPULATION), ("n_gen", GENERATIONS), seed=SEED, verbose=False)
    seconds = round(time.perf_counter() - started, 2)
    igd = float(IGD(front)(result.F))
    rows = "".join(f"{a:.9g},{b:.9g}\n" for a, b in sorted(map(tuple, result.F)))
    (folder / "front.csv").write_text("f1,f2\n" + rows, encoding="utf-8")        # the found front, for a look at it
    (folder / "result.json").write_text(json.dumps({
        "pymoo_version": importlib.metadata.version("pymoo"), "population": POPULATION, "generations": GENERATIONS, "seed": SEED,
        "igd": igd, "limit": LIMIT, "non_dominated_points": len(result.F)}, indent=1), encoding="utf-8")
    measured = f"pymoo {importlib.metadata.version('pymoo')}: IGD {igd:.4f} ({len(result.F)} non-dominated points, {result.algorithm.evaluator.n_eval:,} evaluations)"
    yield Row("pymoo", case, expected, measured, PASS if igd < LIMIT else FAIL, seconds)
