"""Precompute the FEMM results for the interactive flux page (flux_explorer_build.py makes the HTML from them).

One slider at a time from the first-pass design: 4 sliders x 5 stops = 17 rotors, each with the band as air and as steel.
Per rotor and band (coarse mesh, 60 solves):
  fields  : 0 A, plus 10 and 20 A at current angle 0/30/60/90 deg = 9 current states, each at 4 rotor positions
            (0/30/60/90 deg electrical; the whole picture repeats every 120 deg electrical, rotated by 30 deg mech)
            -> the vector potential A on every mesh node (= the flux lines) + torque and flux linkages
  scalars : 10 and 20 A at current angle 15/45/75 deg, same 4 positions -> torque and flux linkages only
Writes out/flux_explorer/g{gid:02d}_{band}.npz and states.csv. About 12 min on 8 workers (FEMM 4.2; FEMM 5.0: about 2 min;
the numbers are in the README, Expected numbers).
Run: uv run python flux_explorer_sweep.py   (from this folder; then flux_explorer_sweep120.py, then flux_explorer_build.py)
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

OUT = os.path.abspath("out/flux_explorer")
WORKERS = 8
BASE = {"mag_y": 6.0, "mag_x": 39.4136, "v_half_deg": 60.0, "web": 1.5}
SLIDERS = {"mag_y": [5.0, 5.5, 6.0, 6.3, 6.6],                 # all buildable (scanned 10-28: 4.9 .. 6.6)
           "mag_x": [39.08, 39.25, 39.4136, 39.5, 39.57],      # 39.075 .. 39.575
           "v_half_deg": [57.0, 58.5, 60.0, 61.0, 61.75],      # 56.75 .. 61.75
           "web": [1.05, 1.3, 1.5, 1.65, 1.8]}                 # 1.05 .. 1.8
BANDS = ("air", "steel")
FIELD_I, FIELD_G, POS, SCALAR_G = (10.0, 20.0), (0.0, 30.0, 60.0, 90.0), (0.0, 30.0, 60.0, 90.0), (15.0, 45.0, 75.0)


def geometries():
    """gid 0 = the first-pass design, then the 4 other stops of each slider (the rest stay at the first-pass design)."""
    out = [(0, dict(BASE))]
    for name, stops in SLIDERS.items():
        for v in stops:
            if v != BASE[name]:
                out.append((len(out), {**BASE, name: v}))
    return out


def field_states():
    """(I, gamma) pairs whose fields are kept, in the order the page expects."""
    return [(0.0, 0.0)] + [(i, g) for i in FIELD_I for g in FIELD_G]


def worker(task):
    import femm
    gid, x, band = task
    p = mg.Params(barrier=band, **x)
    fem = os.path.join(OUT, "fem", f"g{gid:02d}_{band}.fem")
    rows, A, xy, tri = [], [], None, None
    femm_session.open_femm()                              # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(p, mg.build(p), fem)                   # coarse mesh: the field is for looking at
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        cur = lambda i, g: (-i * math.sin(math.radians(g)), i * math.cos(math.radians(g)))
        jobs = [(i, g, th, True) for i, g in field_states() for th in POS]
        jobs += [(i, g, th, False) for i in FIELD_I for g in SCALAR_G for th in POS]
        for i, g, th, keep in jobs:
            r = model.solve(th, *cur(i, g))
            rows.append({"gid": gid, "band": band, "I": i, "gamma": g, "theta_e": th, "field": keep,
                         **{k: r[k] for k in ("T_age", "psi_d", "psi_q", "psi_a", "psi_b", "psi_c", "rotor_deg",
                                              "psi_retries")}})
            if keep:
                sol = femfile.read_ans(fem[:-4] + ".ans")
                A.append(sol["A"].astype(np.float32))
                if xy is None:
                    xy, tri = sol["xy"].astype(np.float32), sol["tri"].astype(np.int32)
    finally:
        femm_session.close_femm()
    np.savez_compressed(os.path.join(OUT, f"g{gid:02d}_{band}.npz"), xy=xy, tri=tri, A=np.array(A))
    return rows


if __name__ == "__main__":
    os.makedirs(os.path.join(OUT, "fem"), exist_ok=True)
    tasks = [(gid, x, band) for gid, x in geometries() for band in BANDS]
    t = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part in pool.map(worker, tasks) for r in part]
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "states.csv"), index=False)
    print(f"{len(tasks)} models, {len(df)} solves in {time.perf_counter() - t:.0f} s on {WORKERS} workers; "
          f"flux-linkage retries: {int(df.psi_retries.sum())}")
