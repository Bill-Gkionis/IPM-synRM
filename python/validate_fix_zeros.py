"""Re-solve the validation slice points whose FEMM flux linkage came back as exactly 0.0 (2026-10-28 anomaly).

In validate_full.py, 11 of the 24 no-current slice solves returned psi_a = 0.000 while their torque was right.
The same sequence re-run alone returned correct values, which made the cause look environmental. It was not: until
2026-11-04 validate_full.py gave its no-load job and its 10 A job of the same angles ONE file name, so two FEMMs saved,
solved and read the same slice_10k.fem and .ans at once, and one of them could be handed a half-written file or the other
job's answer (a flux linkage of exactly 0 with the right torque is one way that shows; a FEMM that crashes or hangs is
another). validate_full.py now gives every job files of its own, and its rows come out identical from run to run, so
there is nothing left to re-solve: with FEMM 5.0 this prints "0 rows to re-solve" and writes validate.csv back with the
same numbers. It stays as the check that no slice row is exactly 0.0, and as the repair for a table that an older
validate_full.py made. It re-solves just those rows, one fresh FEMM per point, and patches validate.csv.
Run: uv run python validate_fix_zeros.py   (from this folder, after validate_full.py; ~2 min with FEMM 4.2; a second
with FEMM 5.0, where there is nothing to re-solve)
"""
import os

import pandas as pd

import femm_session
import motorgeo as mg

OUT = os.path.abspath("out/validate")

if __name__ == "__main__":
    import femm
    path = os.path.join(OUT, "validate.csv")
    df = pd.read_csv(path)
    bad = df[(df.model == "motorgeo 90") & (df.psi_a == 0.0)]
    print(f"{len(bad)} rows to re-solve: theta_e = {sorted(bad.theta_e)}")
    P = mg.Params()
    for i, row in bad.iterrows():
        fem = os.path.join(OUT, f"fix_{int(row.theta_e):03d}.fem")
        femm_session.open_femm()                  # a fresh FEMM of our own for every point: own port, below-normal
        try:
            mg.to_femm(P, mg.build(P), fem, fine=True)
            femm.mi_close()
            r = mg.FemmModel(P, fem).solve_abc(row.theta_e, row.I, -row.I / 2, -row.I / 2)
        finally:
            femm_session.close_femm()
        print(f"  theta_e {row.theta_e:5.1f}: psi_a {r['psi_a'] * 1e3:8.4f} mWb (was 0), T {r['T_age']:+.4f} (was {row['T']:+.4f})")
        df.loc[i, ["psi_a", "T"]] = r["psi_a"], r["T_age"]
    df.to_csv(path, index=False)
    print("patched", path)
