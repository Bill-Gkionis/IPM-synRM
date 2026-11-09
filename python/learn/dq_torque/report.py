"""The numbers of the project note Torque in d-q - field, current and saliency, from run.py's results.

It solves nothing: it reads out/learn/dq_torque/solves.csv, probes.csv and gap.csv, prints seven small tables, and writes
summary.json beside them (figures.py reads it). Every quantity uses the Air band study's own method (barrier_plots.py):
  psi_PM = psi_d at no load, averaged over the 6 rotor positions;
  L_d = (psi_d at gamma 90 - psi_PM) / i_d,   L_q = psi_q at gamma 0 / i_q,   both averaged over the same 6 positions;
  the MTPA angle = the top of a parabola through the best sampled angle and its two neighbours.

Run, from the code folder:  uv run python learn/dq_torque/report.py
"""
import json
import math

import numpy as np
import pandas as pd

import common as c

PP = c.params().pole_pairs


def by_gamma(df, variant, I=c.I_RATED):
    """Mean over the 6 rotor positions, one row per current angle."""
    d = df[(df.variant == variant) & (df.I == I) & df.part.isin(["gamma", "inductance"])]
    return d.groupby("gamma")[["T_age", "T_wst", "T_dq", "psi_d", "psi_q", "id", "iq"]].mean()


def at(gm, g):
    return gm.iloc[int(np.argmin(np.abs(gm.index.values - g)))]


def psi_pm(df, variant):
    return df[(df.variant == variant) & (df.part == "noload")].psi_d.mean()


def inductances(df, variant, I):
    gm, pm = by_gamma(df, variant, I), psi_pm(df, variant)
    Ld = (at(gm, 90).psi_d - pm) / at(gm, 90).id
    Lq = at(gm, 0).psi_q / at(gm, 0).iq
    return Ld, Lq


def mtpa(gm):
    k = min(max(int(np.argmax(gm.T_age.values)), 1), len(gm) - 2)
    a2, a1, a0 = np.polyfit(gm.index.values[k - 1:k + 2], gm.T_age.values[k - 1:k + 2], 2)
    g = -a1 / (2 * a2)
    return g, a0 + a1 * g + a2 * g * g


if __name__ == "__main__":
    df = pd.read_csv(c.OUT / "solves.csv")
    probes = pd.read_csv(c.OUT / "probes.csv")
    gap = pd.read_csv(c.OUT / "gap.csv")
    pd.set_option("display.width", 170)
    S = {}

    # 1. magnets on against magnets off, at 10 A
    t1 = {}
    for v in c.VARIANTS:
        gm = by_gamma(df, v)
        Ld, Lq = inductances(df, v, c.I_RATED)
        g, T = mtpa(gm)
        t1[v] = {"psi_PM [mWb]": psi_pm(df, v) * 1e3, "L_d [mH]": Ld * 1e3, "L_q [mH]": Lq * 1e3, "saliency L_q/L_d": Lq / Ld,
                 "MTPA angle [deg]": g, "torque at MTPA [N.m]": T,
                 "torque, all q (gamma 0) [N.m]": at(gm, 0).T_age, "torque, all -d (gamma 90) [N.m]": at(gm, 90).T_age,
                 "psi_d with all-q current [mWb]": at(gm, 0).psi_d * 1e3, "psi_q with all-d current [mWb]": at(gm, 90).psi_q * 1e3}
        S[v] = {"psi_pm": psi_pm(df, v), "Ld": Ld, "Lq": Lq, "mtpa_gamma": g, "T_mtpa": T}
    print("1. magnets on and off, 10 A peak, FEMM 5.0 fine mesh\n" + pd.DataFrame(t1).round(4).to_string())

    # 2. the three terms at the MTPA angle, from psi_PM, L_d and L_q (the d-q model), against FEMM's torque
    on = S["on"]
    i_d, i_q = c.currents(c.I_RATED, on["mtpa_gamma"])
    psi_d, psi_q = on["psi_pm"] + on["Ld"] * i_d, on["Lq"] * i_q
    terms = c.three_terms(PP, on["psi_pm"], psi_d, psi_q, i_d, i_q)
    print(f"\n2. the three terms at MTPA ({on['mtpa_gamma']:.1f} deg: i_d {i_d:.3f} A, i_q {i_q:.3f} A), magnets on")
    for k, v in terms.items():
        print(f"   {k:22s} {v:+.4f} N.m")
    print(f"   {'sum (the d-q model)':22s} {sum(terms.values()):+.4f} N.m    FEMM at the same angle: {on['T_mtpa']:+.4f} N.m")
    S["terms_at_mtpa"] = terms

    # 3. two torque routes inside FEMM (the air-gap element against the d-q formula on FEMM's own flux linkages)
    gm = by_gamma(df, "on")
    r = at(gm, 22.5)
    print(f"\n3. two routes at gamma 22.5 deg, magnets on: AGE {r.T_age:.4f}, stress tensor {r.T_wst:.4f}, "
          f"1.5 p (psi_d i_q - psi_q i_d) {r.T_dq:.4f} N.m")

    # 4. L_d, L_q and the saliency against current
    t4 = {}
    for v in c.VARIANTS:
        for I in c.CURRENTS:
            Ld, Lq = inductances(df, v, I)
            t4[(v, I)] = {"L_d [mH]": Ld * 1e3, "L_q [mH]": Lq * 1e3, "saliency": Lq / Ld}
    t4 = pd.DataFrame(t4).T
    print("\n4. against current (A peak)\n" + t4.round(4).to_string())
    S["saliency"] = {v: {f"{I:g}": float(t4.loc[(v, I), "saliency"]) for I in c.CURRENTS} for v in c.VARIANTS}
    S["inductance"] = {v: {f"{I:g}": [float(t4.loc[(v, I), "L_d [mH]"]), float(t4.loc[(v, I), "L_q [mH]"])]
                           for I in c.CURRENTS} for v in c.VARIANTS}

    # 4b. the axis check, magnets off at the highest current (the only state with saliency and no magnet flux)
    I_hi = max(c.CURRENTS)
    gm_hi = by_gamma(df, "off", I_hi)
    peak = gm_hi.T_age.abs().max()
    t0, t90 = at(gm_hi, 0).T_age, at(gm_hi, 90).T_age
    cross = max(abs(at(gm_hi, 0).psi_d / at(gm_hi, 0).psi_q), abs(at(gm_hi, 90).psi_q / at(gm_hi, 90).psi_d))
    print(f"\n4b. axes, magnets off at {I_hi:g} A: torque peak {peak:.4f} N.m at {gm_hi.T_age.idxmax():g} deg; "
          f"at 0 deg {t0:.4f}, at 90 deg {t90:.4f} N.m ({100 * max(abs(t0), abs(t90)) / peak:.2f} % of the peak); "
          f"a current on one axis makes flux on the other of at most {100 * cross:.3f} % of its own")
    S["axis_check"] = {"I": I_hi, "T_peak": peak, "gamma_peak": float(gm_hi.T_age.idxmax()), "T_0": t0, "T_90": t90,
                       "cross_pct": 100 * cross}

    # 5. |B| in pole 0's bridges and rib
    print("\n5. |B| [T] at the middle of pole 0's bridges and rib, rotor at 0 deg\n"
          + probes.pivot_table(index="where", columns=["variant", "state"], values="B_T").round(3).to_string())
    S["probes"] = probes.to_dict(orient="records")

    # 6. the hand route to the magnet torque (sum of B I l over the winding) and the net flux at MTPA
    p = c.params()
    B1 = c.harmonic_amplitudes(gap.br.values)[1]
    r_gap, length = (p.r_ro + p.airgap / 2) / 1e3, p.stack / 1e3
    import motorgeo as mg
    _, kw1 = mg.phase_axes(p)
    Ns = c.series_turns(p)
    T_hand = 3 * length * r_gap * kw1 * Ns * B1 * i_q
    psi_from_gap = kw1 * Ns * 2 * B1 * r_gap * length / p.pole_pairs     # flux per pole 2 B1 r l / p, linked k_w N_s times
    ang_psi = math.degrees(math.atan2(psi_q, psi_d))
    ang_i = math.degrees(math.atan2(i_q, i_d))
    pf = math.cos(math.radians(ang_psi + 90 - ang_i))
    print(f"\n6. air-gap fundamental at no load B1 = {B1:.4f} T (r {r_gap * 1e3:.1f} mm, l {length * 1e3:.0f} mm, "
          f"k_w {kw1:.4f}, N_s {Ns})")
    print(f"   3 l r k_w N_s B1 i_q = {T_hand:.4f} N.m    against the magnet term {terms['magnets x i_q']:.4f} N.m")
    print(f"   k_w N_s 2 B1 r l / p = {psi_from_gap * 1e3:.3f} mWb    against psi_PM read off the winding "
          f"{on['psi_pm'] * 1e3:.3f} mWb")
    print(f"   net flux at MTPA: psi_d {psi_d * 1e3:.2f}, psi_q {psi_q * 1e3:.2f} mWb, {ang_psi:.1f} deg from d; "
          f"current at {ang_i:.1f} deg; power factor {pf:.3f}")
    S.update({"B1": B1, "T_hand": T_hand, "psi_from_gap": psi_from_gap, "psi_angle_deg": ang_psi, "pf": pf,
              "kw1": kw1, "Ns": Ns})

    with open(c.OUT / "summary.json", "w", encoding="utf-8") as f:
        json.dump(S, f, indent=1, default=float)
    print(f"\nwrote {c.OUT / 'summary.json'}")
