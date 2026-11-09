"""check_design: the study's black box (optimization/design.py) held against numbers it must reproduce, before anybody trusts it.

What it checks, and against what
    Python only (no FEMM):
      1. the bridge: design.bridge_of() on the reference gives the bench's 1.558 mm (B7: the slot's outer corner, 0.158 mm from the
         guard), and slot_reach() is the sum motorgeo's own build guard makes: at the edge, 0.16 mm more mag_x is still built by both
         (the corner moves outward by less than mag_x does, since it sits off the d-axis), and 0.17 mm more is refused by both
      2. the Fourier peak-to-peak: a sine of the 18th harmonic sampled 15 times over 60 electrical degrees has a peak-to-peak of
         2 x its amplitude wherever the samples fall (within 1e-4: the rebuilt waveform is read on a grid of 3,000 points); the
         samples' own max - min is smaller
      3. screen() refuses what motorgeo refuses, with a reason and a positive constraint value, and passes the reference
    FEMM (a pool of 3 workers, coarse and fine mesh):
      4. the reference at 10 A: its own MTPA angle within 1 degree of B1's 22.9 (fine mesh) and its mean torque within 1 % of B1's
         1.12914 N.m on the fine mesh, and within 1 % of B7's 1.12493 N.m (coarse, two positions at 22.9 degrees) on the coarse one
      5. scoring="J" gives the reference exactly the scoring="AT" current (its coil area is the reference's)
      6. the torque from the flux linkages (T_dq) and from the air-gap element agree within 2 % on the mean

Run (from the Python to FEMM folder; the solver slot is taken for about a minute)
    uv run python -m optimization.check_design
Writes out\\optimization\\check_design\\evaluations.csv and prints every check with PASS or FAIL. Exit status 0 if all pass.
"""
from __future__ import annotations

import math
import sys

import numpy as np

import femm_queue
import motorgeo as mg
from optimization import design

B1_FINE_T, B1_GAMMA, B7_COARSE_T = 1.12914, 22.936, 1.12493


def main() -> int:
    results = []

    def check(name, ok, said):
        results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {said}")

    ref = design.REFERENCE
    b = design.bridge_of(ref)
    check("1. the reference's outer bridge is B7's 1.558 mm", abs(b - 1.558) < 0.001, f"{b:.4f} mm, guard {ref.bridge} mm")
    for more, both_build in ((0.16, True), (0.17, False)):
        try:
            mg.build(design.params({"mag_x": ref.mag_x + more}))
            built = True
        except mg.GeometryError:
            built = False
        s = design.screen({"mag_x": ref.mag_x + more})
        agree = built == both_build and s["feasible"] == both_build and (s["G"]["guard_bridge"] <= 0) == both_build
        check(f"1. {more} mm more mag_x: motorgeo and screen() both {'build' if both_build else 'refuse'} it", agree,
              f"motorgeo builds it: {built}; screen feasible: {s['feasible']} ({s['reason'] or 'no reason'}); guard_bridge {s['G']['guard_bridge']:+.4f} mm")

    t = np.arange(15) / 15
    worst = min(np.ptp(np.sin(2 * math.pi * 3 * t + phase)) for phase in np.linspace(0, 2 * math.pi, 50))
    exact = max(abs(design.fourier_pkpk(np.sin(2 * math.pi * 3 * t + phase)) - 2.0) for phase in np.linspace(0, 2 * math.pi, 50))
    check("2. the Fourier peak-to-peak of the 18th harmonic is 2 x its amplitude wherever the samples fall", exact < 1e-4,
          f"largest error {exact:.2e}; the samples' own max - min falls to {worst:.3f}")

    s_ref = design.screen({})
    check("3. screen() passes the reference", s_ref["feasible"] and all(v <= 0 for v in s_ref["G"].values()),
          "G: " + ", ".join(f"{k} {v:+.3f}" for k, v in s_ref["G"].items()))
    s_bad = design.screen({"r_shoe": ref.r_si + 0.1})
    check("3. screen() refuses a shoe below the opening's straight part", (not s_bad["feasible"]) and s_bad["G"]["shoe_taper"] > 0, s_bad["reason"])

    with femm_queue.femm_pool("optimization: check_design", workers=3, task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        coarse = design.Evaluator(pool, run="check_design", scoring="AT")([{}], stage="coarse AT")[0]
        coarse_j = design.Evaluator(pool, run="check_design", scoring="J")([{}], stage="coarse J")[0]
        fine = design.Evaluator(pool, run="check_design", scoring="AT", fine=True)([{}], stage="fine AT")[0]
    for name, row in (("coarse", coarse), ("fine", fine)):
        print(f"      {name}: gamma {row['gamma']:.2f} deg, T_mean {row['T_mean']:.5f} N.m, ripple {row['ripple_pp']:.4f} N.m "
              f"({row['ripple_pct']:.2f} %; samples {row['ripple_pp_samples']:.4f}), psi_pm {1e3 * row['psi_pm']:.3f} mWb, "
              f"L_d {1e3 * row['L_d']:.3f} mH, L_q {1e3 * row['L_q']:.3f} mH, pf {row['pf']:.3f}, {row['nodes']} nodes, {row['seconds']:.1f} s")
    check("4. the fine mesh's own MTPA angle is B1's within 1 degree", abs(fine["gamma"] - B1_GAMMA) < 1.0, f"{fine['gamma']:.2f} against {B1_GAMMA}")
    check("4. the fine mesh's mean torque is B1's within 1 %", abs(fine["T_mean"] / B1_FINE_T - 1) < 0.01,
          f"{fine['T_mean']:.5f} against {B1_FINE_T} N.m ({100 * (fine['T_mean'] / B1_FINE_T - 1):+.2f} %)")
    check("4. the coarse mesh's mean torque is B7's within 1 %", abs(coarse["T_mean"] / B7_COARSE_T - 1) < 0.01,
          f"{coarse['T_mean']:.5f} against {B7_COARSE_T} N.m ({100 * (coarse['T_mean'] / B7_COARSE_T - 1):+.2f} %)")
    check("5. scoring J gives the reference the scoring AT current", abs(coarse_j["amps"] - coarse["amps"]) < 1e-9,
          f"{coarse_j['amps']:.6f} A and {coarse['amps']:.6f} A")
    p_amps = design.scoring_amps(s_ref["metrics"], "P")
    check("5. scoring P gives the reference the scoring AT current too", abs(p_amps - coarse["amps"]) < 1e-9, f"{p_amps:.6f} A")
    on_guard = [design.screen(design.decode({"bridge_mm": 1.4, "mag_y": y, "v_half_deg": v})) for y in (5.2, 5.6, 6.0, 6.4) for v in (50, 60, 70, 80)]
    worst = max(s["G"]["guard_bridge"] for s in on_guard)
    check("1. bridge_mm = 1.4 exactly is never refused by the guard (16 designs on the guard)", worst <= 0, f"largest guard_bridge {worst:+.2e} mm")
    check("6. the torque from the flux linkages agrees with the air-gap element's", abs(fine["T_dq_mean"] / fine["T_mean"] - 1) < 0.02,
          f"T_dq {fine['T_dq_mean']:.5f}, T_age {fine['T_mean']:.5f} N.m")
    print(f"{sum(results)} of {len(results)} checks pass")
    return 0 if all(results) else 1


if __name__ == "__main__":
    sys.exit(main())
