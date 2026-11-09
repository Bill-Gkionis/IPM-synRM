"""How much faster do more FEMM workers make a batch? Same 48 coarse-mesh solves on 1, 2, 4, 6, 8, 10 workers.

Each worker: its own process, its own FEMM, its own .fem/.ans. Writes out/scaling/scaling.csv and a figure.
Run: uv run python worker_scaling.py [figure.png]   (from this folder; about 4 min with FEMM 4.2, about 33 s with FEMM 5.0:
the numbers are in the README, Expected numbers; nothing else should be running)
"""
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

import femm_session
import motorgeo as mg

OUT = os.path.abspath("out/scaling")
THETAS = list(range(0, 48))                       # 48 rotor positions, i_q = 10 A


def worker(args):
    import femm
    wid, n_workers, thetas = args
    p = mg.Params()
    fem = os.path.join(OUT, f"n{n_workers}_w{wid}.fem")        # private file per worker
    femm_session.open_femm()                  # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(p, mg.build(p), fem)
        femm.mi_close()
        m = mg.FemmModel(p, fem)
        return [(th, m.solve(th, 0, 10)["T_age"], os.getpid()) for th in thetas]
    finally:
        femm_session.close_femm()


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    rows, ref = [], None
    for n in (1, 2, 4, 6, 8, 10):
        t = time.perf_counter()
        with ProcessPoolExecutor(n) as pool:
            res = [r for part in pool.map(worker, [(w, n, THETAS[w::n]) for w in range(n)]) for r in part]
        wall = time.perf_counter() - t
        T = np.array([r[1] for r in sorted(res)])
        ref = T if ref is None else ref
        rows.append((n, wall, len({r[2] for r in res}), float(np.abs(T - ref).max())))
        print(f"{n:2d} workers: {wall:6.1f} s, speed-up {rows[0][1] / wall:4.2f}x, {rows[-1][2]} processes, "
              f"max |T - T(1 worker)| = {rows[-1][3]:.2e} N.m", flush=True)
    with open(os.path.join(OUT, "scaling.csv"), "w") as f:
        f.write("workers,wall_s,processes,max_diff_Nm\n" + "".join(f"{a},{b:.2f},{c},{d:.3e}\n" for a, b, c, d in rows))

    if len(sys.argv) > 1:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        n = np.array([r[0] for r in rows]); s = rows[0][1] / np.array([r[1] for r in rows])
        fig, ax = plt.subplots(figsize=(6.4, 3.8), facecolor="#fcfcfb")
        ax.set_facecolor("#fcfcfb")
        ax.plot(n, n, color="#c3c2b7", lw=1.2, label="perfect (n workers = n times faster)")
        ax.plot(n, s, "-o", color="#2a78d6", lw=1.8, ms=5, label="measured")
        for a, b in zip(n, s):
            ax.annotate(f"{b:.1f}x", (a, b), textcoords="offset points", xytext=(0, -14), ha="center",
                        fontsize=8, color="#52514e")
        ax.set_xlabel("FEMM workers (i7-12700K: 8 fast + 4 efficient cores)", color="#52514e")
        ax.set_ylabel("speed-up vs 1 worker", color="#52514e")
        ax.set_xticks(n); ax.grid(color="#e1e0d9", lw=0.6); ax.set_axisbelow(True)
        for sp in ax.spines.values():
            sp.set_color("#c3c2b7")
        ax.tick_params(colors="#898781")
        ax.legend(frameon=False, loc="upper left", fontsize=8.5)
        ax.set_title(f"48 coarse-mesh solves: {rows[0][1]:.0f} s on 1 worker, {min(r[1] for r in rows):.0f} s at best",
                     loc="left", fontsize=10, color="#0b0b0b")
        fig.tight_layout()
        fig.savefig(sys.argv[1], dpi=150)
        print("wrote", sys.argv[1])
