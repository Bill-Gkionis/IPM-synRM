"""The hand-drawn rotor with and without the air band behind each magnet. Same magnets, same mesh: only the band's material changes.

Phase 1: no-load sweep, cogging sweep, torque vs current angle at 4 current levels.
Phase 2: torque ripple at each variant's own MTPA angle, plus field snapshots (.ans kept for plots).
Writes out/barrier/*.csv and out/barrier/fields/*.ans. Plots: barrier_plots.py.
About 59 min measured on 10 workers (fine mesh, FEMM 4.2); 8 workers will take somewhat longer. (FEMM 5.0: about 2 min
on 8 workers; the numbers are in the README, Expected numbers.)
Run: uv run python barrier_study.py [--coarse]   (from this folder; one FEMM-using script at a time)
"""
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import femm_session
import motorgeo as mg

OUT = os.path.abspath("out/barrier")
VARIANTS = {"with band": "air", "no band": "steel"}   # the band is drawn in both, so both get the identical mesh
WORKERS = 8                             # 10 made the desktop stutter on 10-28; 8 leaves room for other work on the PC
FINE = "--coarse" not in sys.argv
I_LEVELS = (5, 10, 15, 20)              # A peak; x 10 turns per coil side = 50 ... 200 A-turns
I_RATED = 10                            # 100 A-turns peak per coil side, about 6 A/mm2 rms at 45 % slot fill
GAMMAS = np.arange(0, 90.1, 7.5)        # current angle, measured from +q toward -d (0 = all q current)
MEAN_POS = range(0, 60, 10)             # 6 rotor positions whose average cancels the 6th/12th/18th/24th ripple
GAP_DEG = list(np.arange(0, 90, 0.25))  # air-gap samples, stator frame, one 90 deg sector


def worker(chunk):
    """One FEMM, one private .fem/.ans per chunk, every job of the chunk on the same variant."""
    import femm
    cid, variant, jobs = chunk
    p = mg.Params(barrier=VARIANTS[variant])
    fem = os.path.join(OUT, "fem", f"c{cid:03d}.fem")
    rows = []
    femm_session.open_femm()                  # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(p, mg.build(p), fem, fine=FINE)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        for job in jobs:
            t = time.perf_counter()
            ans = os.path.join(OUT, "fields", job["ans"] + ".ans") if job.get("ans") else None
            r = model.solve(job["theta_e"], job["id"], job["iq"], gap_deg=GAP_DEG if ans else (), keep_ans=ans)
            rows.append({**job, **r, "variant": variant, "seconds": time.perf_counter() - t, "pid": os.getpid()})
    finally:
        femm_session.close_femm()
    return rows


def run(jobs, label, size=12):
    """Split jobs into chunks per variant and solve them on the pool. Returns a DataFrame."""
    chunks, cid = [], 0
    for v in VARIANTS:
        mine = [j for j in jobs if j["variant"] == v]
        for k in range(0, len(mine), size):
            chunks.append((cid, v, mine[k:k + size]))
            cid += 1
    t = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part in pool.map(worker, chunks) for r in part]
    df = pd.DataFrame(rows).drop(columns=["gap_br", "gap_bt"])
    print(f"{label}: {len(rows)} solves in {time.perf_counter() - t:.0f} s on {WORKERS} workers "
          f"({df.seconds.mean():.1f} s per solve, {df.pid.nunique()} worker processes)", flush=True)
    for r in rows:                                    # air-gap samples go to their own file
        if r.get("ans"):
            pd.DataFrame({"deg": GAP_DEG, "br": r["gap_br"], "bt": r["gap_bt"]}).to_csv(
                os.path.join(OUT, f"gap_{r['ans']}.csv"), index=False)
    return df


def mtpa(df_gamma, current):
    """Current angle of maximum mean torque, from a parabola through the best sampled angle and its neighbours."""
    g = df_gamma[df_gamma.I == current].groupby("gamma").T_age.mean()
    k = int(np.argmax(g.values))
    k = min(max(k, 1), len(g) - 2)
    a, b, c = np.polyfit(g.index[k - 1:k + 2], g.values[k - 1:k + 2], 2)
    return float(-b / (2 * a))


if __name__ == "__main__":
    for d in ("fem", "fields"):
        os.makedirs(os.path.join(OUT, d), exist_ok=True)
    tag = lambda v: v.replace(" ", "_")

    jobs = []
    for v in VARIANTS:
        jobs += [{"study": "noload", "variant": v, "theta_e": float(th), "id": 0.0, "iq": 0.0,
                  "ans": f"noload_{tag(v)}" if th == 0 else ""} for th in range(0, 360, 3)]
        jobs += [{"study": "cogging", "variant": v, "theta_e": float(th), "id": 0.0, "iq": 0.0}
                 for th in np.arange(0, 20.01, 0.5)]
        jobs += [{"study": "gamma", "variant": v, "theta_e": float(th), "I": float(i), "gamma": float(g),
                  "id": -i * math.sin(math.radians(g)), "iq": i * math.cos(math.radians(g))}
                 for i in I_LEVELS for g in GAMMAS for th in MEAN_POS]
    df1 = run(jobs, "phase 1")
    df1.to_csv(os.path.join(OUT, "phase1.csv"), index=False)

    jobs = []
    for v in VARIANTS:
        g = mtpa(df1[(df1.study == "gamma") & (df1.variant == v)], I_RATED)
        i_d, i_q = -I_RATED * math.sin(math.radians(g)), I_RATED * math.cos(math.radians(g))
        print(f"  {v}: MTPA at {I_RATED} A is gamma = {g:.1f} deg (id {i_d:.2f}, iq {i_q:.2f} A)")
        jobs += [{"study": "ripple", "variant": v, "theta_e": float(th), "I": float(I_RATED), "gamma": g,
                  "id": i_d, "iq": i_q, "ans": f"mtpa_{tag(v)}" if th == 0 else ""} for th in range(0, 61)]
        jobs += [{"study": "snapshot", "variant": v, "theta_e": 0.0, "I": float(I_RATED), "gamma": 0.0,
                  "id": 0.0, "iq": float(I_RATED), "ans": f"qaxis_{tag(v)}"}]
    df2 = run(jobs, "phase 2", size=8)
    df2.to_csv(os.path.join(OUT, "phase2.csv"), index=False)
