"""Reference point for the air-band study: the magnet INSIDE the air band (the 10-22 prototype, Params(pocket="open")).

Same tests as barrier_study.py's MTPA sweep at 10 A, fine mesh: no-load flux linkage, torque vs current angle.
Writes out/barrier/reference_open.csv. About 7 min on 8 workers (FEMM 4.2; FEMM 5.0: about 13 s; the numbers are in the
README, Expected numbers).
Run: uv run python reference_open.py   (from this folder)
"""
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor

import pandas as pd

import femm_session
import motorgeo as mg
from barrier_study import GAMMAS, MEAN_POS, OUT

WORKERS = 8


def worker(chunk):
    import femm
    cid, jobs = chunk
    p = mg.Params(pocket="open")
    fem = os.path.join(OUT, "fem", f"ref{cid:03d}.fem")
    rows = []
    femm_session.open_femm()                  # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(p, mg.build(p), fem, fine=True)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        for job in jobs:
            rows.append({**job, **model.solve(job["theta_e"], job["id"], job["iq"])})
    finally:
        femm_session.close_femm()
    return rows


if __name__ == "__main__":
    jobs = [{"study": "noload", "theta_e": float(th), "id": 0.0, "iq": 0.0} for th in range(0, 360, 15)]
    jobs += [{"study": "gamma", "theta_e": float(th), "I": 10.0, "gamma": float(g),
              "id": -10 * math.sin(math.radians(g)), "iq": 10 * math.cos(math.radians(g))}
             for g in GAMMAS for th in MEAN_POS]
    chunks = [(k, jobs[k::WORKERS]) for k in range(WORKERS)]
    t = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part in pool.map(worker, chunks) for r in part]
    df = pd.DataFrame(rows).drop(columns=["gap_br", "gap_bt"])
    df.to_csv(os.path.join(OUT, "reference_open.csv"), index=False)
    m = mg.check(mg.Params(pocket="open"), mg.build(mg.Params(pocket="open")))[0]
    g = df[df.study == "gamma"].groupby("gamma").T_age.mean()
    print(f"{len(df)} solves in {time.perf_counter() - t:.0f} s | magnet area {m['magnet_area_mm2']:.1f} mm2 | "
          f"no-load psi_d {df[df.study == 'noload'].psi_d.mean() * 1e3:.2f} mWb | "
          f"best mean torque at 10 A {g.max():.4f} N.m at gamma {g.idxmax():.1f} deg")
