"""salib_check: SALib finds the Sobol indices of the Ishigami function, a function whose indices are known exactly.

What it does, and why
    Sensitivity analysis asks "which inputs does the output depend on?". Sobol's first-order index S_i is the share of the
    output's variance (its spread) that input i explains on its own. The study of the motor will use it to find which
    dimensions (a barrier angle, a magnet width) matter before an optimizer spends FEMM runs on them. SALib computes it from a
    sample of the function; the check is that on a function with a known answer it finds the known answer.

The Ishigami function (the standard test of sensitivity tools)
    f(x1, x2, x3) = sin(x1) + a sin^2(x2) + b x3^4 sin(x1),   each x uniform on [-pi, pi],   a = 7, b = 0.1
    Its variance and the parts of it have closed forms (integrate over the cube):
        V   = 1/2 + a^2/8 + b pi^4/5 + b^2 pi^8/18   = 13.845
        V1  = (1 + b pi^4/5)^2 / 2                   = 4.346     (what x1 explains alone)
        V2  = a^2/8                                  = 6.125     (x2 alone)
        V3  = 0                                      (x3 only acts together with x1: the term b x3^4 sin(x1))
    so S1 = V1/V = 0.3139, S2 = V2/V = 0.4424, S3 = 0. The remaining 24 % of the variance is the x1-x3 interaction, which the
    TOTAL index ST1 = 0.5576 and ST3 = 0.2437 show and the first-order indices cannot.

Method
    SALib 1.6.0 (it has no __version__: the version is read from the package metadata). sobol.sample(problem, N = 1024, seed) draws
    a scrambled Sobol sequence of N (2 D + 2) = 8,192 points (D = 3 inputs, second-order indices on); the function is evaluated
    here in numpy (SALib's own copy of Ishigami is compared with it at every point, so a bug in either shows); sobol.analyze
    estimates the indices (with bootstrap confidence intervals, 100 resamples). The estimate has a sampling error, which is why
    the check allows 0.05; the seed is fixed so that the number is the same every time.

Numbers to expect (2026-11-06): S1 = 0.3156, 0.4383, 0.0014 for seed 1 (errors +0.002, -0.004, +0.001), 8,192 evaluations,
    0.1 s. Over 20 seeds the worst error of any index was 0.024, so the limit of 0.05 is about twice the worst case.

Self-check (the answers follow).
  1. S1 + S2 + S3 = 0.7563. Where is the rest of the variance?
  2. Why is the sample 8,192 points when N = 1024?
  3. The estimate of S3 is 0.0014 and not 0. Is SALib wrong?
Answers.
  1. In the interaction of x1 and x3 (the term b x3^4 sin(x1) is not a sum of a function of x1 and a function of x3): 24.4 % of the
     variance. Second-order indices (SALib computes them too) would show S13 = 0.244.
  2. The method needs N points for each of the A and B matrices, and N more for each input with one column swapped, in both
     directions for second-order indices: N (2 D + 2) = 1024 x 8.
  3. No. A Monte Carlo estimate from a finite sample has an error of this size; the exact value is 0 and the estimate would be
     off by +-0.02 for other seeds. The bootstrap confidence interval of S3 (0.0014 +- 0.057) contains 0.
"""
from __future__ import annotations

import importlib.metadata
import json
import math
import time
from pathlib import Path

from . import FAIL, PASS, Row, fresh

A, B = 7.0, 0.1
N = 1024
SEED = 1
TOLERANCE = 0.05


def ishigami_first_order(a: float, b: float) -> tuple[float, float, float]:
    """The exact first-order Sobol indices (S1, S2, S3) of the Ishigami function."""
    v1 = (1 + b * math.pi ** 4 / 5) ** 2 / 2
    v2 = a ** 2 / 8
    variance = 0.5 + a ** 2 / 8 + b * math.pi ** 4 / 5 + b ** 2 * math.pi ** 8 / 18
    return v1 / variance, v2 / variance, 0.0


def run(folder: Path):
    import numpy as np
    from SALib.analyze import sobol as analyze
    from SALib.sample import sobol as sample
    from SALib.test_functions import Ishigami

    fresh(folder)
    problem = {"num_vars": 3, "names": ["x1", "x2", "x3"], "bounds": [[-math.pi, math.pi]] * 3}
    exact = ishigami_first_order(A, B)
    case = f"Sobol indices of the Ishigami function (a = {A:g}, b = {B:g}), N = {N}"
    expected = f"first-order S1, S2, S3 within {TOLERANCE} of {exact[0]:.3f}, {exact[1]:.3f}, {exact[2]:.3f}"
    started = time.perf_counter()
    points = sample.sample(problem, N, seed=SEED)
    f = np.sin(points[:, 0]) + A * np.sin(points[:, 1]) ** 2 + B * points[:, 2] ** 4 * np.sin(points[:, 0])
    ours_vs_theirs = float(np.max(np.abs(f - Ishigami.evaluate(points, A, B))))
    result = analyze.analyze(problem, f, seed=SEED)
    seconds = round(time.perf_counter() - started, 2)
    s1 = [float(v) for v in result["S1"]]
    errors = [s - e for s, e in zip(s1, exact)]
    version = importlib.metadata.version("SALib")
    (folder / "result.json").write_text(json.dumps({
        "salib_version": version, "N": N, "seed": SEED, "evaluations": len(points), "S1": s1, "S1_conf": [float(v) for v in result["S1_conf"]],
        "exact_S1": list(exact), "errors": errors, "tolerance": TOLERANCE}, indent=1), encoding="utf-8")
    measured = (f"SALib {version}: S1 = {s1[0]:.3f}, {s1[1]:.3f}, {s1[2]:.3f} (errors {errors[0]:+.3f}, {errors[1]:+.3f}, {errors[2]:+.3f}), "
                f"{len(points):,} evaluations, seed {SEED}")
    ok = all(abs(e) < TOLERANCE for e in errors) and ours_vs_theirs < 1e-12
    if ours_vs_theirs >= 1e-12:
        measured += f"; SALib's own Ishigami differs from the formula by {ours_vs_theirs:.2e}"
    yield Row("SALib", case, expected, measured, PASS if ok else FAIL, seconds)
