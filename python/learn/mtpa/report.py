"""The numbers of the project note MTPA, from run.py's results (out/learn/mtpa/solves.csv).

For each current: FEMM's best angle (the top of a parabola through the best sampled angle and its neighbours), the
torque there, the apparent L_d and L_q, and the constant-parameter formula's angle,

    sin(gamma) = ( -psi_PM + sqrt(psi_PM^2 + 8 dL^2 I^2) ) / (4 dL I),     dL = L_q - L_d,

once with the 10 A constants (what is obtained by measuring the motor once) and once with that current's own L_d and L_q.
Then what the 10 A formula's angle costs in torque, and the lookup table: torque request -> (i_d, i_q).
It writes summary.json beside the results (figures.py reads it).

Run, from the code folder:  uv run python learn/mtpa/report.py
"""
import json
import math
import sys
from pathlib import Path

sys.path.append(str(Path(__file__).resolve().parents[1] / "dq_torque"))   # after this folder, so its own report.py wins
import numpy as np
import pandas as pd

import common as dq
import paths

OUT = paths.path("out") / "learn" / "mtpa"
PP = dq.params().pole_pairs
LUT_TORQUES = (0.5, 1.0, 1.5, 2.0, 3.0, 4.0)             # N.m, the rows of the lookup table printed below


def formula_angle(psi_pm, dL, I):
    """The MTPA angle [deg] of the constant-parameter d-q model (dT/dgamma = 0 solved for sin gamma)."""
    s = (-psi_pm + math.sqrt(psi_pm ** 2 + 8 * dL ** 2 * I ** 2)) / (4 * dL * I)
    return math.degrees(math.asin(s))


def peak(gm):
    """FEMM's best angle and torque: a parabola through the best sampled angle and its two neighbours."""
    g, T = gm.index.values, gm.T_age.values
    k = min(max(int(np.argmax(T)), 1), len(T) - 2)
    a2, a1, a0 = np.polyfit(g[k - 1:k + 2], T[k - 1:k + 2], 2)
    best = -a1 / (2 * a2)
    return best, a0 + a1 * best + a2 * best * best


if __name__ == "__main__":
    df = pd.read_csv(OUT / "solves.csv")
    psi_pm = df[df.I == 0].psi_d.mean()
    rows = {}
    for I, d in df[df.I > 0].groupby("I"):
        gm = d.groupby("gamma")[["T_age", "psi_d", "psi_q", "id", "iq"]].mean()
        Ld = (gm.loc[90.0].psi_d - psi_pm) / gm.loc[90.0].id
        Lq = gm.loc[0.0].psi_q / gm.loc[0.0].iq
        best, T = peak(gm[gm.index <= 60])
        rows[I] = {"gamma": best, "T": T, "Ld": Ld, "Lq": Lq, "curve": gm}
    c10 = rows[10.0]
    dL10 = c10["Lq"] - c10["Ld"]

    table = {}
    for I, r in rows.items():
        g10 = formula_angle(psi_pm, dL10, I)
        g_own = formula_angle(psi_pm, r["Lq"] - r["Ld"], I)
        curve = r["curve"][r["curve"].index <= 60]
        T_at_g10 = float(np.interp(g10, curve.index.values, curve.T_age.values))   # FEMM's torque at the formula's angle
        T_at_own = float(np.interp(g_own, curve.index.values, curve.T_age.values))
        T_pm = 1.5 * PP * psi_pm * I * math.cos(math.radians(r["gamma"]))
        table[I] = {"FEMM angle [deg]": r["gamma"], "formula, 10 A constants [deg]": g10,
                    "formula, own constants [deg]": g_own, "torque at FEMM's angle [N.m]": r["T"],
                    "torque at the 10 A formula's angle [N.m]": T_at_g10,
                    "lost by the formula's angle [%]": 100 * (1 - T_at_g10 / r["T"]),
                    "lost by the own-constants angle [%]": 100 * (1 - T_at_own / r["T"]),
                    "L_d [mH]": r["Ld"] * 1e3, "L_q [mH]": r["Lq"] * 1e3, "L_q - L_d [mH]": (r["Lq"] - r["Ld"]) * 1e3,
                    "reluctance share [%]": 100 * (r["T"] - T_pm) / r["T"]}
    t = pd.DataFrame(table)
    pd.set_option("display.width", 170)
    print(f"psi_PM = {psi_pm * 1e3:.3f} mWb (no load); FEMM 5.0, fine mesh, band-free rotor, magnets on\n")
    print("1. the best angle against current\n" + t.round(3).to_string())

    # 2. the lookup table: torque request -> current and angle, by interpolation along FEMM's MTPA points
    I_ = np.array([0.0] + list(rows)); T_ = np.array([0.0] + [r["T"] for r in rows.values()])
    G_ = np.array([0.0] + [r["gamma"] for r in rows.values()])
    lut = {}
    for Tr in LUT_TORQUES:
        if Tr > T_.max():
            continue
        I_r = float(np.interp(Tr, T_, I_)); g_r = float(np.interp(Tr, T_, G_))
        i_d, i_q = dq.currents(I_r, g_r)
        lut[Tr] = {"current [A]": I_r, "angle [deg]": g_r, "i_d [A]": i_d, "i_q [A]": i_q}
    print("\n2. the lookup table (linear interpolation between FEMM's MTPA points)\n" + pd.DataFrame(lut).round(2).to_string())

    S = {"psi_pm": psi_pm, "dL10": dL10,
         "rows": {f"{I:g}": {k: float(v) for k, v in vals.items()} for I, vals in table.items()},
         "points": [{"I": float(I), "gamma": r["gamma"], "T": r["T"]} for I, r in rows.items()],
         "lut": {f"{k:g}": v for k, v in lut.items()}}
    with open(OUT / "summary.json", "w", encoding="utf-8") as f:
        json.dump(S, f, indent=1)
    print(f"\nwrote {OUT / 'summary.json'}")
