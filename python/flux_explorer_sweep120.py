"""The 120-deg-electrical end position for flux_explorer_sweep.py's field states.

The picture repeats every 120 deg electrical (rotated by 30 deg mech), but blending from 90 to 120 inside one period
needs the 120 position on the same mesh, so it is solved explicitly here: 9 field states per rotor and band.
Writes out/flux_explorer/g{gid:02d}_{band}_p120.npz and states_p120.csv. About 2 min on 8 workers (FEMM 4.2; FEMM 5.0:
about 20 s; the numbers are in the README, Expected numbers).
Run: uv run python flux_explorer_sweep120.py   (from this folder, after flux_explorer_sweep.py)
"""
import math
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import femfile
import femm_session
import motorgeo as mg
from flux_explorer_sweep import BANDS, OUT, WORKERS, field_states, geometries


def worker(task):
    import femm
    gid, x, band = task
    p = mg.Params(barrier=band, **x)
    fem = os.path.join(OUT, "fem", f"g{gid:02d}_{band}_p120.fem")
    rows, A, xy = [], [], None
    femm_session.open_femm()                  # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(p, mg.build(p), fem)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        for i, g in field_states():
            r = model.solve(120.0, -i * math.sin(math.radians(g)), i * math.cos(math.radians(g)))
            rows.append({"gid": gid, "band": band, "I": i, "gamma": g, "theta_e": 120.0, "field": True,
                         **{k: r[k] for k in ("T_age", "psi_d", "psi_q", "psi_a", "psi_b", "psi_c", "rotor_deg",
                                              "psi_retries")}})
            sol = femfile.read_ans(fem[:-4] + ".ans")
            A.append(sol["A"].astype(np.float32))
            xy = sol["xy"].astype(np.float32)          # kept to prove it's the same mesh as the main sweep
    finally:
        femm_session.close_femm()
    np.savez_compressed(os.path.join(OUT, f"g{gid:02d}_{band}_p120.npz"), A=np.array(A), xy=xy)
    return rows


if __name__ == "__main__":
    tasks = [(gid, x, band) for gid, x in geometries() for band in BANDS]
    t = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part in pool.map(worker, tasks) for r in part]
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "states_p120.csv"), index=False)
    print(f"{len(rows)} solves in {time.perf_counter() - t:.0f} s")
