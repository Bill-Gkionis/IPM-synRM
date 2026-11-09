"""Torque over one full on-load ripple period (0-60 deg electrical), coarse vs fine mesh, on parallel FEMM workers.

Writes out/torque_mesh.png and prints mean, ripple and the 6th/12th/18th/24th harmonics.
The fine-mesh pass takes about 1.5 minutes on 6 workers (FEMM 4.2). (FEMM 5.0: about 7 s for the whole script; the
numbers are in the README, Expected numbers.)
Run: uv run python mesh_study.py   (from this folder; ~2 min in all with FEMM 4.2)
"""
import os
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import femm_session
import motorgeo as mg

WORKERS = 6


def worker(job):
    import femm
    wid, fine, thetas = job
    p = mg.Params()
    fem = os.path.abspath(f"out/w{wid}_{int(fine)}.fem")    # each worker owns its .fem and .ans
    femm_session.open_femm()                                # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(p, mg.build(p), fem, fine=fine)
        femm.mi_close()
        m = mg.FemmModel(p, fem)
        return [m.solve(th, 0, 10) for th in thetas]
    finally:
        femm_session.close_femm()


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs("out", exist_ok=True)
    thetas = list(range(0, 61, 2))
    res = {}
    for fine in (False, True):
        t = time.perf_counter()
        jobs = [(i, fine, thetas[i::WORKERS]) for i in range(WORKERS)]
        with ProcessPoolExecutor(WORKERS) as ex:
            res[fine] = sorted((r for part in ex.map(worker, jobs) for r in part), key=lambda r: r["theta_e"])
        T = [r["T_age"] for r in res[fine]]
        h = np.abs(np.fft.rfft(T[:-1])) / (len(T) - 1) * 2       # 0..58 deg in 2-deg steps = one 60-deg period
        print(f"{'fine' if fine else 'coarse'}: {len(T)} solves in {time.perf_counter() - t:.1f} s on {WORKERS} workers"
              f" | mean {np.mean(T):.4f}, p-p {max(T) - min(T):.4f} N.m | T(0) {T[0]:.4f}, T(60) {T[-1]:.4f}"
              f" | max |T_age - T_wst| {max(abs(r['T_age'] - r['T_wst']) for r in res[fine]):.4f}")
        print("   harmonics per electrical cycle:", ", ".join(f"{6 * k}th {h[k]:.4f}" for k in range(1, 5)))

    fig, ax = plt.subplots(figsize=(7, 3.6))
    for fine, style in ((True, "-o"), (False, "--s")):
        ax.plot([r["theta_e"] for r in res[fine]], [r["T_age"] for r in res[fine]], style, ms=3.5,
                label="fine mesh (smartmesh on, ~24k nodes)" if fine else "coarse mesh (smartmesh off, ~8k nodes)")
    ax.set_xlabel("rotor position, electrical degrees (full on-load ripple period = 60)")
    ax.set_ylabel("torque [N·m]")
    ax.set_title("i_q = 10 A, i_d = 0, 10 turns/coil side (placeholder winding)", fontsize=10)
    ax.grid(alpha=0.3); ax.legend(fontsize=8); fig.tight_layout()
    fig.savefig("out/torque_mesh.png", dpi=160)
