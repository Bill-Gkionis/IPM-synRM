"""MTPA on our band-free rotor, from 5 to 40 A: the FEMM runs behind the project note MTPA.

The question it answers: does the best current angle keep sliding toward 45 deg as the current grows, as the
constant-parameter formula says, or does saturation bend it back? So, magnets on, fine mesh, one FEMM under the solver slot:
  - the no-load flux at 6 rotor positions (psi_PM);
  - at 5, 10, 15, ... 40 A peak: the current angle gamma from 0 to 60 deg in 5 deg steps, plus 90 deg (all on -d, for
    L_d), each at the 6 rotor positions whose average cancels the 6th, 12th and 18th ripple harmonics.
It reuses learn/dq_torque/common.py: the same rotor (Params(barrier="steel")), the same positions and helpers.

Run, from the code folder:  uv run python learn/mtpa/run.py        (678 solves, about 7 minutes)
Then:                       uv run python learn/mtpa/report.py     (the note's numbers and its lookup table)
                            uv run python learn/mtpa/figures.py    (the note's two figures)
Writes out/learn/mtpa/: solves.csv and mtpa.fem.
"""
import sys
import time
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / "dq_torque"))   # common.py, shared with the d-q note (after this folder, so its own report.py wins)
import numpy as np
import pandas as pd

import common as dq
import femm_queue
import femm_session
import motorgeo as mg
import paths

OUT = paths.path("out") / "learn" / "mtpa"
CURRENTS = tuple(float(i) for i in range(5, 41, 5))           # A peak
ANGLES = tuple(float(g) for g in np.arange(0, 60.1, 5)) + (90.0,)

if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    p = dq.params()
    fem = str(OUT / "mtpa.fem")
    rows, t0 = [], time.perf_counter()
    with femm_queue.solver_slot("learn/mtpa: the angle sweep to 40 A"):
        femm_session.open_femm()                               # this process's own FEMM: own port, below-normal priority
        try:
            mg.to_femm(p, mg.build(p), fem, fine=True)
            import femm
            femm.mi_close()
            model = mg.FemmModel(p, fem)
            rows += [{"I": 0.0, "gamma": np.nan, **model.solve(float(th), 0.0, 0.0)} for th in dq.MEAN_POS]
            for I in CURRENTS:
                for g in ANGLES:
                    i_d, i_q = dq.currents(I, g)
                    rows += [{"I": I, "gamma": g, **model.solve(float(th), i_d, i_q)} for th in dq.MEAN_POS]
                print(f"  {I:g} A done", flush=True)
        finally:
            femm_session.close_femm()
    df = pd.DataFrame(rows).drop(columns=["gap_br", "gap_bt"])
    df.to_csv(OUT / "solves.csv", index=False)
    print(f"{len(df)} solves in {time.perf_counter() - t0:.0f} s; flux read-back retries: {int(df.psi_retries.sum())}")
