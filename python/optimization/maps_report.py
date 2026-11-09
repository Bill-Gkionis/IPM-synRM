"""optimization.maps_report: the analysis of phase 1: it builds the flux map from FEMM's rows, runs every check, and prints the tables the report quotes.

What it is, and why
    maps_run.py leaves FEMM's raw results in out\\optimization\\maps\\chosen\\ (rows_<stage>.csv, frozen\\, loss\\). Nothing in them is a result yet. This file turns
    them into the maps and says how far they can be trusted, section by section (python -m optimization.maps_report SECTION ...; "all" runs them in this order):
      flux      the flux map: the mean of every grid point over its positions (psi_d, psi_q, the three torques, the Fourier peak-to-peak ripple) -> flux_map.npz;
                then the checks on it: the regression against compare.csv, how many positions the mean needs (the 60-position reference, and 6 / 15 against 30 at
                every grid point), the d-q torque against FEMM's, the mirror symmetry, the interpolation error at 12 off-grid points and of a coarser grid, and Maxwell's
                reciprocity of the cross terms
      cross     cross-saturation: the apparent and incremental inductances and the cross terms at the points that matter, and how far the cross terms move the MTPA
                angle and the torque at 1, 2 and 3 times the rated current, against the constant-parameter formulas and against the old rotor's finding (learn\\mtpa)
      frozen    FEMM's frozen permeability: the magnets' flux and the current's flux told apart at six points, against the apparent and incremental inductances
      loss      the iron loss: the per-region sums of the loss grid -> iron_map.npz, the checks (positions, interpolation, mirror, period) and the loss tables at both lamination
                thicknesses (0.3 mm, the headline, and the datasheet's 0.356 mm)
    Every pass limit below was written before the numbers existed and is not changed afterwards: a FAIL is reported with its number, not argued away.

The maths the tables use, each step built from the one before
    d-q torque               T = 1.5 p (psi_d i_q - psi_q i_d)        (amplitude-invariant d-q: 3/2 from the three phases, p pole pairs)
    apparent inductance      L_d,app = (psi_d - psi_PM0) / i_d,  L_q,app = psi_q / i_q          (psi_PM0: the no-load flux linkage)
    incremental inductance   L_dd = d psi_d / d i_d,  L_qq = d psi_q / d i_q,  cross terms L_dq = d psi_d / d i_q,  L_qd = d psi_q / d i_d
    reciprocity              psi_d = d Wc / d i_d and psi_q = d Wc / d i_q (Wc the coenergy, up to one constant factor), so L_dq = L_qd
    constant-parameter MTPA  T(gamma) = 1.5 p [psi_PM I cos(gamma) + (L_d - L_q) I^2 sin(gamma) cos(gamma)]; dT/dgamma = 0 gives
                             sin(gamma) = ( -psi_PM + sqrt(psi_PM^2 + 8 dL^2 I^2) ) / (4 dL I),  dL = L_q - L_d        (learn\\mtpa)
    the same without cross terms (the "decoupled" model): psi_d(i_d, i_q) replaced by psi_d(i_d, 0) and psi_q by psi_q(0, i_q), each axis saturating on its own current only;
                             its MTPA angle against the full maps' angle is what the cross terms move
    "torque lost"            the full maps' torque at a model's angle, against the maps' best at that current: what trusting the model's angle costs

Run (from the Python to FEMM folder; no FEMM, seconds)
    uv run python -m optimization.maps_report flux            after maps_run's regression, kcheck, mirror, flux, heldout, flux2 and fluxneg stages
    uv run python -m optimization.maps_report cross frozen    needs flux_map.npz (frozen also the frozen stage)
    uv run python -m optimization.maps_report loss            after maps_run's loss and losscheck stages
Writes report_<section>.json beside the rows, flux_map.npz, iron_map.npz; prints everything.

Numbers to expect (2026-11-09; the printouts are saved as report_<section>.txt beside the rows)
    flux     1. regression 2e-14 relative (limit 1e-9): PASS. 2. positions: the means are good to 1e-5 of |psi| and 0.3 mN.m with 15, to 2e-4 and 8 mN.m with 6; the ripple at 15 positions misses its
             declared limit at the weak-ripple point (4.5 mN.m against 3.96) and reads the rated point's 2.4 % low, so the maps use 30 positions (2 degrees apart): FAIL at 15, PASS at 20 and 30.
             3. d-q torque against the air-gap torque: median 0.006 %, largest 0.085 %: PASS. 4. mirror symmetry at the 7 check points: FAIL (psi_q 0.52 ... 0.76 % of |psi| at 3 of them, limit
             0.5 %), over the whole plane the median asymmetry is 0.2 ... 0.7 %: the negative half was solved. 5. interpolation, 12 FEMM points on no node: 0.014 % of |psi|, 0.67 mN.m: PASS;
             the grid of twice the step would miss the limit at the nodes it did not see (0.35 % against 0.25 %). 6. reciprocity of the cross terms: median 0.015 %, largest 0.58 %.
    cross    at 1, 2, 3 times the rated current the maps' MTPA angle is 23.5, 37.4, 43.3 deg (the torque 1.81, 3.65, 5.09 N.m); the constant-parameter formula with the rated point's constants gives
             17.3, 26.2, 30.9 deg (0.8, 3.3, 3.6 % of the torque lost); the model with saturation but no cross terms 17.4, 10.9, 0 deg (0.7, 14.6, 31 % lost).
    frozen   closure of the frozen method 1e-7 relative: PASS; the magnets' flux in the loaded steel is 1.05, 1.15, 1.21 times the no-load flux along the MTPA path.
    loss     14. period PASS; 15. 72 against 144 positions within 1.1 % of a region's loss: PASS; 16. negative i_q within 1.8 %: PASS, off-grid points: teeth 2.3 %, yoke 4.1 %, rotor 5.5 % at worst
             (limit 5 % per region): FAIL for the rotor at one of three points, the total within 2.8 %; 17b. the full loss is above the fundamental-only one in the stator (x1.2 ... 3.0): PASS.

Self-check (the answers follow).
  1. The d-q torque from the flux maps and FEMM's air-gap torque agree to a fraction of a per cent on the position mean. Why would they not agree at one position?
  2. Why does the decoupled model (no cross terms) tell how much the cross terms matter, and why is it better than comparing with constant inductances?
  3. A check passes at 15 positions. What must still be true before the same 15 are used at every grid point?
Answers.
  1. The d-q formula is the energy balance averaged over the cycle; at one position the instantaneous torque has terms from d psi / d theta (the ripple) that the formula does
     not contain. They average to zero over a period.
  2. Constants differ from the maps in two ways at once (saturation along the axis and cross-coupling between them); the decoupled model keeps the saturation and removes only the
     cross-coupling, so the difference is the cross terms' effect and nothing else.
  3. That the positions chosen resolve the harmonics that matter at the WORST point, not only at the point they were tested at: the 60-position reference is run at points from no
     load to the corner of the grid, and 6, 15 and 30 positions are compared at every point.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import motorgeo as mg
import paths
from optimization import design, iron_loss, maps, maps_femm, maps_run

OUT = maps_run.OUT
I_REGRESSION = 1e-9                  # relative, against compare.csv: the same code path must give the same numbers
K_TOL_PSI, K_TOL_T, K_TOL_T_ABS = 1e-3, 2e-3, 2e-3         # a K-position mean against the 60-position one: psi relative to |psi|; T relative to |T| plus 2 mN.m
K_TOL_RIPPLE, K_TOL_RIPPLE_ABS = 0.05, 2e-3                   # the Fourier peak-to-peak: 5 % plus 2 mN.m
TOL_DQ_REL, TOL_DQ_ABS = 5e-3, 0.01                           # |T_dq - T_age| at a grid point: 0.5 % of T plus 0.01 N.m
TOL_MIRROR_PSI, TOL_MIRROR_T_REL, TOL_MIRROR_T_ABS = 5e-3, 1e-2, 0.01   # the mirrored point: psi within 0.5 % of |psi|; T within 1 % plus 0.01 N.m
TOL_INTERP_PSI, TOL_INTERP_T_REL, TOL_INTERP_T_ABS = 2.5e-3, 5e-3, 0.01   # the interpolation at an off-grid FEMM point: 0.25 % of |psi|; T 0.5 % plus 0.01 N.m
TOL_LOSS_INTERP, TOL_LOSS_MIRROR, TOL_LOSS_N = 0.05, 0.05, 0.05   # the iron loss of a region at 3000 and 16,330 rpm: interpolation, mirror and 72-against-144-position differences
TOL_PERIOD_PSI, TOL_PERIOD_T_REL, TOL_PERIOD_T_ABS = 1e-3, 2e-3, 2e-3   # the mean over the full 360 degrees against the mean over 60: psi of |psi|; T as the K check's
K_LIST = (2, 3, 4, 5, 6, 10, 12, 15, 20, 30)                  # numbers of equally spaced positions per 60 degrees that divide 60
PP = 4
RPM_RATED, RPM_TOP = 3000.0, 16330.0                          # the study's rated speed and the rotor's safe top speed (stress check; a design decision: about 16,000 rpm)
LAMINATIONS_MM = (0.3, 0.356)                                 # the CAD's (headline) and the datasheet's


# ============================================================================================== the rows
def load_rows(stage: str) -> pd.DataFrame:
    """rows_<stage>.csv, duplicates (a job repeated after a crash) dropped, sorted."""
    df = pd.read_csv(OUT / f"rows_{stage}.csv")
    return df.drop_duplicates(subset=["key", "i_d", "i_q", "theta_e"]).sort_values(["key", "theta_e"]).reset_index(drop=True)


def point_summary(g: pd.DataFrame) -> dict:
    """The mean over the positions of one point and its torque ripple (Fourier-rebuilt peak-to-peak and the samples' own)."""
    g = g.sort_values("theta_e")
    T = g.T_age.values
    return {"i_d": float(g.i_d.iloc[0]), "i_q": float(g.i_q.iloc[0]), "n": int(len(g)), "psi_d": float(g.psi_d.mean()), "psi_q": float(g.psi_q.mean()),
            "T_age": float(T.mean()), "T_wst": float(g.T_wst.mean()), "T_dq": float(g.T_dq.mean()),
            "ripple_pp": design.fourier_pkpk(T) if len(T) >= 8 else float("nan"), "ripple_samples": float(T.max() - T.min()),
            "psi_d_pp": float(g.psi_d.max() - g.psi_d.min()), "psi_q_pp": float(g.psi_q.max() - g.psi_q.min()), "retries": int(g.psi_retries.sum())}


class Report:
    def __init__(self):
        self.ctx = maps_run.Context()
        self.I_R = self.ctx.I_R
        self.out: dict = {}

    def say(self, text: str = "") -> None:
        print(text, flush=True)

    def verdict(self, ok: bool) -> str:
        return "PASS" if ok else "FAIL"

    def save(self, name: str) -> None:
        with open(OUT / f"report_{name}.json", "w", encoding="utf-8") as f:
            json.dump(self.out.get(name, {}), f, indent=1, default=float)

    # ============================================================================================== flux
    def build_flux_map(self) -> dict:
        pos = pd.concat([load_rows("flux"), load_rows("flux2")]).drop_duplicates(subset=["key", "i_d", "i_q", "theta_e"])
        neg = load_rows("fluxneg")
        summ = [point_summary(g) for _, g in pos.groupby("key")] + [point_summary(g) for _, g in neg.groupby("key")]
        n_pos, n_neg = len(maps_run.FLUX_ID_STEPS) * len(maps_run.FLUX_IQ_STEPS), len(maps_run.FLUX_ID_STEPS) * (len(maps_run.FLUX_IQ_STEPS) - 1)
        want_n = len(maps_run.POSITIONS_FLUX) + len(maps_run.POSITIONS_FLUX2)
        if len(summ) != n_pos + n_neg or any(s["n"] != (want_n if s["i_q"] >= -1e-9 else len(maps_run.POSITIONS_FLUX)) for s in summ):
            raise RuntimeError(f"the flux stages are not complete: {len(summ)} points of {n_pos + n_neg}")
        exact_d, exact_q = {round(s["i_d"], 9): s["i_d"] for s in summ}, {round(s["i_q"], 9): s["i_q"] for s in summ}        # nodes are matched at 1e-9 A, but kept exact
        ids = np.array([exact_d[k] for k in sorted(exact_d)])
        iqs = np.array([exact_q[k] for k in sorted(exact_q)])
        grid = {k: np.full((len(ids), len(iqs)), np.nan) for k in ("psi_d", "psi_q", "T_age", "T_wst", "T_dq", "ripple_pp", "ripple_samples", "psi_d_pp", "psi_q_pp")}
        retries = 0
        for s in summ:
            a, b = sorted(exact_d).index(round(s["i_d"], 9)), sorted(exact_q).index(round(s["i_q"], 9))
            for k in grid:
                grid[k][a, b] = s[k]
            retries += s["retries"]
        a0, b0 = int(np.argmin(np.abs(ids))), int(np.argmin(np.abs(iqs)))
        psi_pm0 = grid["psi_d"][a0, b0]
        np.savez(OUT / "flux_map.npz", id_A=ids, iq_A=iqs, **grid, rated_current_A=self.I_R, rated_gamma_deg=self.ctx.gamma_R, psi_pm0_Wb=psi_pm0, pole_pairs=PP,
                 depth_mm=20.0, positions_per_point_iq_ge_0=want_n, positions_per_point_iq_lt_0=len(maps_run.POSITIONS_FLUX))
        self.say(f"flux map: {len(ids)} x {len(iqs)} points (i_d {ids[0]:.2f} ... {ids[-1]:.2f} A, i_q {iqs[0]:.2f} ... {iqs[-1]:.2f} A: both halves solved by FEMM), "
                 f"{want_n} positions per point for i_q >= 0 and {len(maps_run.POSITIONS_FLUX)} for i_q < 0, flux read-back retries {retries}; no-load psi_PM {1e3 * psi_pm0:.4f} mWb; "
                 "written flux_map.npz")
        self.say(f"   torque {np.nanmin(grid['T_age']):.3f} ... {np.nanmax(grid['T_age']):.3f} N.m, ripple (Fourier pk-pk) {np.nanmin(grid['ripple_pp']):.3f} ... "
                 f"{np.nanmax(grid['ripple_pp']):.3f} N.m")
        return {"ids": ids, "iqs": iqs, "grid": grid, "psi_pm0": psi_pm0, "retries": retries}

    def check_regression(self) -> bool:
        import csv
        df = load_rows("regression")
        with open(paths.path("out") / "optimization" / "final" / "compare.csv", encoding="utf-8") as f:
            cmp_ = {r["design"]: r for r in csv.DictReader(f)}
        ok_all, rows = True, []
        for label, name in (("chosen", "chosen (fine knee)"), ("reference", "reference")):
            nl = df[df.key == f"{label}-noload"]
            mt = df[df.key == f"{label}-mtpa"].sort_values("theta_e")
            psi_pm = nl.psi_d.mean()
            i_d, i_q = mt.i_d.iloc[0], mt.i_q.iloc[0]
            mine = {"T_mean": mt.T_age.mean(), "ripple_pp": design.fourier_pkpk(mt.T_age.values), "psi_pm_mWb": 1e3 * psi_pm,
                    "L_d_mH": 1e3 * (mt.psi_d.mean() - psi_pm) / i_d, "L_q_mH": 1e3 * mt.psi_q.mean() / i_q}
            for k, v in mine.items():
                want = float(cmp_[name][k])
                rel = abs(v / want - 1)
                ok_all &= rel <= I_REGRESSION
                rows.append({"design": label, "quantity": k, "this run": v, "compare.csv": want, "relative difference": rel})
        t = pd.DataFrame(rows)
        self.say("1. regression: this pipeline against out\\optimization\\final\\compare.csv (the old reference rotor at 10 A and the chosen design at 11.86 A, "
                 f"fine mesh, limit {I_REGRESSION:g} relative): {self.verdict(ok_all)}")
        self.say(t.to_string(index=False, float_format=lambda v: f"{v:.10g}"))
        self.out.setdefault("flux", {})["regression"] = {"pass": bool(ok_all), "max_relative_difference": float(t["relative difference"].max())}
        return ok_all

    def check_positions(self) -> dict:
        """How many positions does the mean need? The 60-position reference (1 degree apart) at the 8 points of the kcheck stage, sub-sampled to each K."""
        df = load_rows("kcheck")
        table = []
        for key, g in df.groupby("key"):
            g = g.sort_values("theta_e")
            T60, pd60, pq60 = g.T_age.values, g.psi_d.values, g.psi_q.values
            ref = {"psi_d": pd60.mean(), "psi_q": pq60.mean(), "T": T60.mean()}
            psi_abs = math.hypot(ref["psi_d"], ref["psi_q"])
            rip60 = design.fourier_pkpk(T60)
            for K in K_LIST:
                idx = np.arange(0, 60, 60 // K)
                e_psi = max(abs(pd60[idx].mean() - ref["psi_d"]), abs(pq60[idx].mean() - ref["psi_q"])) / psi_abs
                e_T = abs(T60[idx].mean() - ref["T"])
                row = {"point": key, "i_d": float(g.i_d.iloc[0]), "i_q": float(g.i_q.iloc[0]), "K": K, "psi_err": float(e_psi), "T_err_Nm": float(e_T), "T": float(ref["T"]),
                       "ok_mean": bool(e_psi <= K_TOL_PSI and e_T <= K_TOL_T * abs(ref["T"]) + K_TOL_T_ABS)}
                if K >= 15:
                    rip = design.fourier_pkpk(T60[idx])
                    row["ripple_err"] = float(abs(rip - rip60))
                    row["ripple60"] = float(rip60)
                    row["ripple_K"] = float(rip)
                    row["ok_ripple"] = bool(abs(rip - rip60) <= K_TOL_RIPPLE * rip60 + K_TOL_RIPPLE_ABS)
                table.append(row)
        t = pd.DataFrame(table)
        k_min = next((K for K in K_LIST if t[t.K == K].ok_mean.all()), None)
        k_rip = next((K for K in K_LIST if K >= 15 and t[t.K == K].ok_ripple.all()), None)
        self.say("2. positions: the mean of K equally spaced positions per 60 degrees against the 60-position mean, at 8 points from no load to (-3, 3) x rated "
                 f"(limits: psi {K_TOL_PSI:.0e} of |psi|, T {K_TOL_T:.0e} of T + {1e3 * K_TOL_T_ABS:.0f} mN.m; ripple {100 * K_TOL_RIPPLE:.0f} % + {1e3 * K_TOL_RIPPLE_ABS:.0f} mN.m)")
        for K in K_LIST:
            s = t[t.K == K]
            line = (f"   K = {K:2d}: worst psi error {s.psi_err.max():.2e} of |psi|, worst T error {1e3 * s.T_err_Nm.max():7.3f} mN.m, all means within limits: {self.verdict(bool(s.ok_mean.all()))}")
            if K >= 15:
                line += f"; worst ripple error {1e3 * s.ripple_err.max():.3f} mN.m: {self.verdict(bool(s.ok_ripple.all()))}"
            self.say(line)
        self.say(f"   the fewest positions that give every mean within its limit: K = {k_min}; the fewest that also give the ripple: K = {k_rip} "
                 "(15 is the project's standard, design.score's; the rated point's ripple at K = 15 reads "
                 f"{1e3 * float(t[(t.K == 15) & (t.point == 'k1')].ripple_K.iloc[0]):.2f} mN.m against {1e3 * float(t[(t.K == 15) & (t.point == 'k1')].ripple60.iloc[0]):.2f} at K = 60)")
        self.out.setdefault("flux", {})["positions"] = {"K_min_mean": k_min, "K_min_ripple": k_rip, "table": table}
        return {"k_min": k_min, "k_ripple": k_rip, "table": t}

    def check_positions_grid(self) -> None:
        """6, 15 and 30 positions at EVERY positive grid point: the 30 positions are the flux and flux2 stages together (2 degrees apart); the 15 are the flux stage's
        (4 degrees apart), the 6 are 10 degrees apart (the MTPA study's), the 5 are 12 degrees apart."""
        both = pd.concat([load_rows("flux"), load_rows("flux2")]).drop_duplicates(subset=["key", "i_d", "i_q", "theta_e"])
        errs = {K: {"psi": [], "T": [], "ripple": [], "ripple_ref": []} for K in (3, 5, 6, 10, 15)}
        for _, g in both.groupby("key"):
            g = g.sort_values("theta_e")
            th = g.theta_e.values.astype(int)
            T30, pd30, pq30 = g.T_age.values, g.psi_d.values, g.psi_q.values
            psi_abs = math.hypot(pd30.mean(), pq30.mean())
            rip30 = design.fourier_pkpk(T30)
            for K in errs:
                sel = (th % (60 // K)) == 0
                errs[K]["psi"].append(max(abs(pd30[sel].mean() - pd30.mean()), abs(pq30[sel].mean() - pq30.mean())) / psi_abs)
                errs[K]["T"].append(abs(T30[sel].mean() - T30.mean()))
                if K >= 15:
                    errs[K]["ripple"].append(design.fourier_pkpk(T30[sel]) - rip30)
                    errs[K]["ripple_ref"].append(rip30)
        self.say(f"2b. positions at every one of the {len(errs[6]['psi'])} positive grid points: the mean of K positions against the mean of 30 (the maps' own)")
        out = {}
        for K, e in errs.items():
            e_psi, e_T = np.array(e["psi"]), np.array(e["T"])
            line = (f"   K = {K:2d} ({60 // K:2d} deg apart): psi error median {np.median(e_psi):.1e}, largest {e_psi.max():.1e} of |psi|; "
                    f"T error median {1e3 * np.median(e_T):.2f}, largest {1e3 * e_T.max():.2f} mN.m")
            out[K] = {"psi_max": float(e_psi.max()), "psi_median": float(np.median(e_psi)), "T_max_Nm": float(e_T.max()), "T_median_Nm": float(np.median(e_T))}
            if K >= 15:
                r, ref = np.array(e["ripple"]), np.array(e["ripple_ref"])
                line += (f"; ripple (Fourier pk-pk) minus the 30-position one: median {1e3 * np.median(r):+.2f}, range {1e3 * r.min():+.2f} ... {1e3 * r.max():+.2f} mN.m "
                         f"(relative: median {100 * np.median(r / ref):+.2f} %, range {100 * (r / ref).min():+.1f} ... {100 * (r / ref).max():+.1f} %)")
                out[K].update({"ripple_diff_median_Nm": float(np.median(r)), "ripple_diff_min_Nm": float(r.min()), "ripple_diff_max_Nm": float(r.max())})
            self.say(line)
        self.out.setdefault("flux", {})["positions_grid"] = out

    def check_dq_torque(self, built: dict) -> bool:
        g = built["grid"]
        dT = g["T_dq"] - g["T_age"]
        lim = TOL_DQ_REL * np.abs(g["T_age"]) + TOL_DQ_ABS
        ok = bool((np.abs(dT) <= lim).all())
        k = np.unravel_index(np.argmax(np.abs(dT) / lim), dT.shape)
        big = np.abs(g["T_age"]) > 0.3
        rel = np.abs(dT[big]) / np.abs(g["T_age"][big])
        wst = g["T_wst"] - g["T_age"]
        self.say(f"3. d-q torque 1.5 p (psi_d i_q - psi_q i_d) against FEMM's air-gap torque, on the position mean of all {dT.size} grid points "
                 f"(limit {100 * TOL_DQ_REL:g} % of T + {TOL_DQ_ABS} N.m): {self.verdict(ok)}")
        self.say(f"   largest |difference| {np.abs(dT).max():.4f} N.m; where |T| > 0.3 N.m: median {100 * np.median(rel):.3f} %, 95th percentile {100 * np.quantile(rel, 0.95):.3f} %, "
                 f"largest {100 * rel.max():.3f} %; the worst point against its limit is (i_d {built['ids'][k[0]]:.2f}, i_q {built['iqs'][k[1]]:.2f}) A with {dT[k]:+.4f} N.m "
                 f"of T {g['T_age'][k]:.3f}")
        self.say(f"   (the weighted-stress-tensor torque differs from the air-gap torque by at most {np.abs(wst).max():.4f} N.m, median {np.median(np.abs(wst)):.5f})")
        self.out.setdefault("flux", {})["dq_torque"] = {"pass": ok, "max_abs_N_m": float(np.abs(dT).max()), "median_rel": float(np.median(rel)), "p95_rel": float(np.quantile(rel, 0.95)),
                                                         "max_rel": float(rel.max()), "mean_sign_of_T_dq_minus_T_age": float(np.sign(dT[big]).mean()),
                                                         "max_T_wst_minus_T_age": float(np.abs(wst).max())}
        return ok

    def check_mirror(self, built: dict) -> bool:
        plus, minus = load_rows("kcheck"), load_rows("mirror")
        rows, ok_all = [], True
        for km, gm in minus.groupby("key"):
            kp = "k" + km[1:]
            gp = plus[plus.key == kp].sort_values("theta_e")
            gm = gm.sort_values("theta_e")
            pp, pm = point_summary(gp), point_summary(gm)
            scale = math.hypot(pp["psi_d"], pp["psi_q"])
            e_d, e_q = abs(pm["psi_d"] - pp["psi_d"]) / scale, abs(pm["psi_q"] + pp["psi_q"]) / scale
            e_T = abs(pm["T_age"] + pp["T_age"])
            ok = e_d <= TOL_MIRROR_PSI and e_q <= TOL_MIRROR_PSI and e_T <= TOL_MIRROR_T_REL * abs(pp["T_age"]) + TOL_MIRROR_T_ABS
            ok_all &= ok
            Tp, Tm = gp.T_age.values, gm.T_age.values
            best = min(range(60), key=lambda c: float(np.sqrt(np.mean((Tm + Tp[(c - np.arange(60)) % 60]) ** 2))))
            rms_flip = float(np.sqrt(np.mean((Tm + Tp[(best - np.arange(60)) % 60]) ** 2)))
            rows.append({"i_d": pp["i_d"], "i_q": pp["i_q"], "psi_d_err": e_d, "psi_q_err": e_q, "T_plus": pp["T_age"], "T_minus": pm["T_age"], "T_sum": e_T, "ok": bool(ok),
                         "ripple_plus_mNm": 1e3 * pp["ripple_pp"], "ripple_minus_mNm": 1e3 * pm["ripple_pp"], "best_flip_c": best, "rms_T_flip_mNm": 1e3 * rms_flip})
        t = pd.DataFrame(rows)
        self.say(f"4. mirror symmetry psi_d(i_d, -i_q) = psi_d, psi_q(i_d, -i_q) = -psi_q, T(-i_q) = -T on the position mean (60 positions), at {len(t)} points "
                 f"(limits: psi {100 * TOL_MIRROR_PSI:g} % of |psi|; T {100 * TOL_MIRROR_T_REL:g} % + {TOL_MIRROR_T_ABS} N.m): {self.verdict(ok_all)}")
        self.say(t.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        # the whole plane, from the two solved halves
        g, ids, iqs = built["grid"], built["ids"], built["iqs"]
        jp = [j for j, v in enumerate(iqs) if v > 1e-9]
        jm = [int(np.argmin(np.abs(iqs + iqs[j]))) for j in jp]
        mag = np.hypot(ids[:, None], iqs[jp][None, :]) / self.I_R
        psi_scale = np.hypot(g["psi_d"][:, jp], g["psi_q"][:, jp])
        e_d = np.abs(g["psi_d"][:, jm] - g["psi_d"][:, jp]) / psi_scale
        e_q = np.abs(g["psi_q"][:, jm] + g["psi_q"][:, jp]) / psi_scale
        e_T = np.abs(g["T_age"][:, jm] + g["T_age"][:, jp])
        e_T_rel = e_T / np.maximum(np.abs(g["T_age"][:, jp]), 0.05)
        self.say("   the same symmetry over the whole plane, from the two solved halves (error relative to |psi| for psi_d and psi_q, absolute and relative for the torque), by current magnitude:")
        bands = {}
        for lo, hi in ((0, 1), (1, 2), (2, 3), (3, 9)):
            m = (mag > lo) & (mag <= hi)
            if m.any():
                bands[f"{lo}-{hi}"] = {"n": int(m.sum()), "psi_d_max": float(e_d[m].max()), "psi_d_median": float(np.median(e_d[m])), "psi_q_max": float(e_q[m].max()),
                                       "psi_q_median": float(np.median(e_q[m])), "T_max_Nm": float(e_T[m].max()), "T_median_Nm": float(np.median(e_T[m])),
                                       "T_rel_max": float(e_T_rel[m].max()), "T_rel_median": float(np.median(e_T_rel[m]))}
                b = bands[f"{lo}-{hi}"]
                self.say(f"   |i| {lo} ... {hi} x rated ({b['n']:3d} points): psi_d median {100 * b['psi_d_median']:.3f} %, largest {100 * b['psi_d_max']:.3f} %; psi_q median "
                         f"{100 * b['psi_q_median']:.3f} %, largest {100 * b['psi_q_max']:.3f} %; T median {1e3 * b['T_median_Nm']:.1f} mN.m ({100 * b['T_rel_median']:.2f} %), "
                         f"largest {1e3 * b['T_max_Nm']:.1f} mN.m ({100 * b['T_rel_max']:.2f} %)")
        self.say("   (the sign of the asymmetry: the generating half has the larger psi_d and the larger |T| at heavy load: T_minus is more negative than T_plus is positive)")
        self.out.setdefault("flux", {})["mirror"] = {"pass_7_points": bool(ok_all), "table": rows, "plane_bands": bands}
        return ok_all

    def fm_from(self, built: dict, ids=None, iqs=None) -> maps.FluxMaps:
        """A FluxMaps from the built grid (optionally a coarser one: every index of ids / iqs)."""
        sel_a = np.arange(len(built["ids"])) if ids is None else np.asarray(ids)
        sel_b = np.arange(len(built["iqs"])) if iqs is None else np.asarray(iqs)
        z = {"id_A": built["ids"][sel_a], "iq_A": built["iqs"][sel_b], "rated_current_A": self.I_R, "rated_gamma_deg": self.ctx.gamma_R, "psi_pm0_Wb": built["psi_pm0"],
             "pole_pairs": PP, "depth_mm": 20.0}
        for k in ("psi_d", "psi_q", "T_age", "T_wst", "T_dq", "ripple_pp"):
            z[k] = built["grid"][k][np.ix_(sel_a, sel_b)]
        return maps.FluxMaps(z)

    def check_interpolation(self, built: dict) -> dict:
        held = load_rows("heldout")
        rows = []
        fm = self.fm_from(built)
        n_q = len(built["iqs"])
        coarse = self.fm_from(built, ids=np.arange(0, len(built["ids"]), 2), iqs=np.arange(0, n_q, 2))
        for key, g in held.groupby("key"):
            s = point_summary(g)
            row = {"key": key, "i_d": s["i_d"], "i_q": s["i_q"], "T_femm": s["T_age"]}
            scale = math.hypot(s["psi_d"], s["psi_q"])
            for name, m in (("fine", fm), ("coarse", coarse)):
                pd_, pq_ = m.psi_dq(s["i_d"], s["i_q"])
                row[f"{name}_psi_err"] = max(abs(pd_ - s["psi_d"]), abs(pq_ - s["psi_q"])) / scale
                row[f"{name}_T_err"] = abs(m.torque(s["i_d"], s["i_q"]) - s["T_age"])
                row[f"{name}_ripple_err"] = abs(m.ripple(s["i_d"], s["i_q"]) - s["ripple_pp"])
            row["fine_ok"] = bool(row["fine_psi_err"] <= TOL_INTERP_PSI and row["fine_T_err"] <= TOL_INTERP_T_REL * abs(s["T_age"]) + TOL_INTERP_T_ABS)
            row["linear_psi_err"] = max(abs(v - w) for v, w in zip(fm.psi_dq(s["i_d"], s["i_q"], method="linear"), (s["psi_d"], s["psi_q"]))) / scale
            rows.append(row)
        t = pd.DataFrame(rows)
        ok = bool(t.fine_ok.all())
        self.say(f"5. interpolation: the bicubic spline of the grid (step rated/4) against {len(t)} FEMM points that are on no node (15 positions each; limits: psi "
                 f"{100 * TOL_INTERP_PSI:g} % of |psi|, T {100 * TOL_INTERP_T_REL:g} % + {TOL_INTERP_T_ABS} N.m): {self.verdict(ok)}")
        self.say(t[["i_d", "i_q", "T_femm", "fine_psi_err", "fine_T_err", "fine_ripple_err", "linear_psi_err", "coarse_psi_err", "coarse_T_err"]].to_string(index=False, float_format=lambda v: f"{v:.3g}"))
        self.say(f"   worst: fine grid psi {100 * t.fine_psi_err.max():.3f} % of |psi|, T {1e3 * t.fine_T_err.max():.2f} mN.m, ripple {1e3 * t.fine_ripple_err.max():.2f} mN.m (the 12 points' own "
                 f"ripple is the 15-position estimate); bilinear on the same grid psi {100 * t.linear_psi_err.max():.3f} %; a grid of step rated/2 (every second node): "
                 f"psi {100 * t.coarse_psi_err.max():.3f} %, T {1e3 * t.coarse_T_err.max():.2f} mN.m")
        # the coarse grid against the nodes it did not see: more exact test points
        ia, ib = np.arange(0, len(built["ids"]), 2), np.arange(0, n_q, 2)
        mask = np.ones((len(built["ids"]), n_q), bool)
        mask[np.ix_(ia, ib)] = False
        errs, errT = [], []
        for a, b in zip(*np.nonzero(mask)):
            i_d, i_q = built["ids"][a], built["iqs"][b]
            pd_, pq_ = coarse.psi_dq(i_d, i_q)
            sc = math.hypot(built["grid"]["psi_d"][a, b], built["grid"]["psi_q"][a, b])
            errs.append(max(abs(pd_ - built["grid"]["psi_d"][a, b]), abs(pq_ - built["grid"]["psi_q"][a, b])) / sc)
            errT.append(abs(coarse.torque(i_d, i_q) - built["grid"]["T_age"][a, b]))
        errs, errT = np.array(errs), np.array(errT)
        self.say(f"   the coarser grid (step rated/2) at the {len(errs)} grid nodes it did not use: psi error median {100 * np.median(errs):.3f} %, 95th percentile "
                 f"{100 * np.quantile(errs, 0.95):.3f} %, largest {100 * errs.max():.3f} % of |psi| (limit {100 * TOL_INTERP_PSI:g} %: {self.verdict(bool(errs.max() <= TOL_INTERP_PSI))}); "
                 f"T error median {1e3 * np.median(errT):.2f}, largest {1e3 * errT.max():.2f} mN.m")
        self.out.setdefault("flux", {})["interpolation"] = {"pass_fine": ok, "table": rows, "coarse_nodes_psi_err_max": float(errs.max()), "coarse_nodes_psi_err_p95": float(np.quantile(errs, 0.95)),
                                                              "coarse_nodes_T_err_max": float(errT.max()), "coarse_pass": bool(errs.max() <= TOL_INTERP_PSI)}
        return {"ok": ok, "table": t}

    def check_reciprocity(self, built: dict) -> None:
        fm = self.fm_from(built)
        a = np.linspace(built["ids"][0], built["ids"][-1], 41)
        b = np.linspace(0, built["iqs"][-1], 31)
        A, B = np.meshgrid(a, b, indexing="ij")
        L = fm.inductances(A, B)
        d = L["L_dq"] - L["L_qd"]
        ref = np.sqrt(L["L_dd"] * L["L_qq"])
        self.say("6. reciprocity of the cross terms (Maxwell): d psi_d / d i_q against d psi_q / d i_d from the spline, on a 41 x 31 grid of the current plane")
        self.say(f"   difference relative to sqrt(L_dd L_qq): median {100 * np.median(np.abs(d) / ref):.3f} %, 95th percentile {100 * np.quantile(np.abs(d) / ref, 0.95):.3f} %, "
                 f"largest {100 * (np.abs(d) / ref).max():.3f} %; the cross terms themselves reach {1e3 * np.abs(L['L_dq']).max():.3f} mH ({1e3 * np.abs(L['L_qd']).max():.3f} mH the other way)")
        self.out.setdefault("flux", {})["reciprocity"] = {"median_rel": float(np.median(np.abs(d) / ref)), "p95_rel": float(np.quantile(np.abs(d) / ref, 0.95)),
                                                           "max_rel": float((np.abs(d) / ref).max()), "L_dq_max_mH": float(1e3 * np.abs(L["L_dq"]).max())}

    def run_flux(self) -> None:
        self.say("=" * 30 + " flux map and its checks " + "=" * 30)
        pd.set_option("display.width", 220)
        pd.set_option("display.max_columns", 30)
        built = self.build_flux_map()
        self.check_regression()
        self.check_positions()
        self.check_positions_grid()
        self.check_dq_torque(built)
        self.check_mirror(built)
        self.check_interpolation(built)
        self.check_reciprocity(built)
        self.save("flux")

    # ============================================================================================== cross-saturation and the MTPA shift
    @staticmethod
    def formula_angle(psi_pm: float, dL: float, I: float) -> float:
        """The MTPA angle [deg] of the constant-parameter d-q model (learn\\mtpa's formula): dT/dgamma = 0 solved for sin(gamma)."""
        s = (-psi_pm + math.sqrt(psi_pm ** 2 + 8 * dL ** 2 * I ** 2)) / (4 * dL * I)
        return math.degrees(math.asin(min(1.0, s)))

    @staticmethod
    def model_torque(psi_pm: float, Ld: float, Lq: float, I: float, gamma_deg: float) -> float:
        """The constant-parameter model's own torque at (I, gamma) [N.m]: 1.5 p (psi_PM i_q + (L_d - L_q) i_d i_q)."""
        g = math.radians(gamma_deg)
        i_d, i_q = -I * math.sin(g), I * math.cos(g)
        return 1.5 * PP * (psi_pm * i_q + (Ld - Lq) * i_d * i_q)

    def inductance_table(self, fm: maps.FluxMaps) -> pd.DataFrame:
        I_R = fm.rated
        pts = []
        for mult in (1, 2, 3):
            g, i_d, i_q, T = fm.mtpa(mult * I_R)
            pts.append((f"MTPA, {mult} x rated", i_d, i_q))
        for mult in (1, 2, 3):
            pts.append((f"q axis, {mult} x rated", 0.0, mult * I_R))
        for mult in (1, 2, 3):
            pts.append((f"-d axis, {mult} x rated", -mult * I_R, 0.0))
        pts.append(("no load", 0.0, 0.0))
        rows = []
        for name, i_d, i_q in pts:
            psi_d, psi_q = fm.psi_dq(i_d, i_q)
            L = fm.inductances(i_d, i_q)
            rows.append({"point": name, "i_d [A]": i_d, "i_q [A]": i_q, "psi_d [mWb]": 1e3 * psi_d, "psi_q [mWb]": 1e3 * psi_q, "T [N.m]": fm.torque(i_d, i_q),
                         "L_d,app [mH]": 1e3 * L["L_d_app"], "L_q,app [mH]": 1e3 * L["L_q_app"], "L_dd [mH]": 1e3 * L["L_dd"], "L_qq [mH]": 1e3 * L["L_qq"],
                         "L_dq [mH]": 1e3 * L["L_dq"], "L_qd [mH]": 1e3 * L["L_qd"]})
        return pd.DataFrame(rows)

    def run_cross(self) -> None:
        self.say("=" * 30 + " cross-saturation and the MTPA shift " + "=" * 30)
        fm = maps.FluxMaps.load()
        I_R = fm.rated
        psi_pm0 = fm.psi_pm0
        pd.set_option("display.width", 250)
        pd.set_option("display.max_columns", 30)
        tab = self.inductance_table(fm)
        self.say("7. the inductances at the points that matter (FEMM 5.0, fine mesh; psi_PM0 = no-load flux linkage, so L_app = (psi - psi_PM0)/i for d and psi/i for q; "
                 "incremental = partial derivatives of the maps; the cross terms L_dq = d psi_d / d i_q and L_qd = d psi_q / d i_d)")
        self.say(tab.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        small = 0.25 * I_R
        Lq0 = fm.inductances(0.0, small)["L_qq"]
        Ld0 = fm.inductances(-small, 0.0)["L_dd"]
        self.say(f"   on the axes the incremental inductances fall from L_qq {1e3 * Lq0:.3f} mH (at i_q = rated/4) to "
                 + ", ".join(f"{1e3 * fm.inductances(0.0, m * I_R)['L_qq']:.3f} ({m} x)" for m in (1, 2, 3))
                 + f" mH, and L_dd from {1e3 * Ld0:.3f} mH (at i_d = -rated/4) to "
                 + ", ".join(f"{1e3 * fm.inductances(-m * I_R, 0.0)['L_dd']:.3f} ({m} x)" for m in (1, 2, 3)) + " mH")
        g1, d1, q1, T1 = fm.mtpa(I_R)
        L1 = fm.inductances(d1, q1)
        Ld_c, Lq_c = L1["L_d_app"], L1["L_q_app"]
        self.say(f"   the constants of the rated point (the way compare.csv gives them): psi_PM {1e3 * psi_pm0:.3f} mWb, L_d {1e3 * Ld_c:.4f} mH, L_q {1e3 * Lq_c:.4f} mH "
                 f"(saliency {Lq_c / Ld_c:.3f}), MTPA angle {g1:.3f} deg, T {T1:.4f} N.m")

        def T_decoupled(i_d, i_q):
            psi_d_axis = fm.psi_dq(i_d, np.zeros_like(np.asarray(i_q, float)))[0]
            psi_q_axis = fm.psi_dq(np.zeros_like(np.asarray(i_d, float)), i_q)[1]
            return 1.5 * PP * (np.asarray(psi_d_axis) * np.asarray(i_q, float) - np.asarray(psi_q_axis) * np.asarray(i_d, float))

        def argmax_gamma(Tfun, I, step=0.1):
            gam = np.arange(0.0, 90.0 + step / 2, step)
            g = np.radians(gam)
            T = Tfun(-I * np.sin(g), I * np.cos(g))
            k = int(np.clip(np.argmax(T), 1, len(gam) - 2))
            y0, y1, y2 = T[k - 1], T[k], T[k + 1]
            den = y0 - 2 * y1 + y2
            shift = 0.5 * (y0 - y2) / den if den != 0 else 0.0
            return float(gam[k] + np.clip(shift, -1, 1) * step)

        cur = lambda I, gam: (-I * math.sin(math.radians(gam)), I * math.cos(math.radians(gam)))
        # the old rotor's "10 A constants": L_d and L_q measured on the axes at the rated current, used at every current
        Lq_axR = float(fm.psi_dq(0.0, I_R)[1]) / I_R
        Ld_axR = (float(fm.psi_dq(-I_R, 0.0)[0]) - psi_pm0) / (-I_R)
        rows = []
        for mult in (1, 2, 3):
            I = mult * I_R
            g_dq, _, _, T_dq_best = fm.mtpa(I, source="dq")
            g_age, _, _, T_age_best = fm.mtpa(I, source="age")
            Tm = lambda gam: float(fm.torque_dq(*cur(I, gam)))
            Ta = lambda gam: float(fm.torque(*cur(I, gam)))
            Lq_ax = float(fm.psi_dq(0.0, I)[1]) / I
            Ld_ax = (float(fm.psi_dq(-I, 0.0)[0]) - psi_pm0) / (-I)
            Lop = fm.inductances(*cur(I, g_dq))
            Ld_op, Lq_op = Lop["L_d_app"], Lop["L_q_app"]
            models = {
                "full maps (truth)": (g_dq, None),
                "on-axis L of the rated current": (self.formula_angle(psi_pm0, Lq_axR - Ld_axR, I), (psi_pm0, Ld_axR, Lq_axR)),
                "constants of the rated point": (self.formula_angle(psi_pm0, Lq_c - Ld_c, I), (psi_pm0, Ld_c, Lq_c)),
                "each current's own on-axis L": (self.formula_angle(psi_pm0, Lq_ax - Ld_ax, I), (psi_pm0, Ld_ax, Lq_ax)),
                "L at the operating point": (self.formula_angle(psi_pm0, Lq_op - Ld_op, I), (psi_pm0, Ld_op, Lq_op)),
                "saturation, no cross terms": (argmax_gamma(T_decoupled, I), None),
            }
            for name, (gam, const) in models.items():
                if const:
                    predicted = self.model_torque(*const, I, gam)
                elif name.startswith("saturation"):
                    predicted = float(T_decoupled(*cur(I, gam)))
                else:
                    predicted = Tm(gam)
                rows.append({"I [A]": I, "x rated": mult, "model": name, "angle [deg]": gam, "shift vs maps [deg]": gam - g_dq,
                             "maps' torque at this angle [N.m]": Tm(gam), "torque lost [%]": 100 * (1 - Tm(gam) / T_dq_best),
                             "T_age lost [%]": 100 * (1 - Ta(gam) / T_age_best), "the model's own torque [N.m]": predicted,
                             "its error [%]": 100 * (predicted / Tm(gam) - 1)})
        t = pd.DataFrame(rows)
        self.say("8. the MTPA angle: the full maps against five simpler models, at 1, 2 and 3 times the rated current (the maps' torque is 1.5 p (psi_d i_q - psi_q i_d) of the interpolated "
                 "flux maps; 'torque lost' = the maps' torque at the model's angle against their best; 'its error' = the model's own torque against the maps' at that angle)")
        self.say(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        g_age_rows = [(m, float(fm.mtpa(m * I_R, source="age")[0])) for m in (1, 2, 3)]
        self.say("   the MTPA angle from FEMM's air-gap torque map instead of the flux maps: " + ", ".join(f"{m} x: {g:.2f} deg" for m, g in g_age_rows))
        old = pd.DataFrame({"x rated (old rotor: 10, 20, 30 A)": [1, 2, 3], "FEMM angle [deg]": [22.909, 31.190, 38.621], "10 A constants [deg]": [23.509, 31.671, 35.372],
                            "own on-axis constants [deg]": [23.509, 30.8, 30.0], "10 A formula lost [%]": [float("nan"), 0.20, 0.48], "own-constants formula lost [%]": [0.135, 0.12, 2.64]})
        self.say("   for comparison, learn\\mtpa on the OLD rotor (the band-free hand-drawn rotor, air band filled with steel; its rated current 10 A; MTPA note sections 4.1-4.2: "
                 "FEMM angle 14.3 / 22.9 / 27.7 / 31.2 / 35.3 / 38.6 / 41.1 / 43.3 deg at 5 ... 40 A; the 10 A-constants formula 23.5 / 31.7 / 35.4 at 10 / 20 / 30 A; the formula "
                 "fed each current's own L_d, L_q 23.5 / 30.8 / 30.0; torque lost 0.20 and 0.48 % (10 A constants) at 20 and 30 A, 0.12 and 2.64 % (own constants); "
                 "at 40 A: 43.3 against 37.5 and 25.3 deg, 1.46 and 10.07 % lost):")
        self.say(old.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        I_tab = np.linspace(0.0, 3 * I_R, 31)
        g_t, d_t, q_t, T_t = fm.mtpa(I_tab[1:])
        table = {"I_A": I_tab.tolist(), "gamma_deg": [0.0] + g_t.tolist(), "i_d_A": [0.0] + d_t.tolist(), "i_q_A": [0.0] + q_t.tolist(), "T_Nm": [0.0] + T_t.tolist()}
        self.say("9. the MTPA lookup table (the maps' path, torque from the flux maps): current, angle, i_d, i_q, torque")
        self.say(pd.DataFrame(table).iloc[::3].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        self.out["cross"] = {"inductance_table": tab.to_dict(orient="records"), "mtpa_shift": t.to_dict(orient="records"), "mtpa_table": table,
                             "constants_rated": {"psi_pm_Wb": psi_pm0, "L_d_H": Ld_c, "L_q_H": Lq_c, "gamma_deg": g1, "T_Nm": T1}}
        self.save("cross")

    # ============================================================================================== frozen permeability
    def run_frozen(self) -> None:
        self.say("=" * 30 + " FEMM's frozen permeability: the magnets' flux and the current's flux " + "=" * 30)
        fm = maps.FluxMaps.load()
        rows, close = [], []
        folder = OUT / "frozen"
        order = ["mtpa_1x", "mtpa_2x", "mtpa_3x", "q_1x", "q_2x", "q_3x"]
        cols = {n: k for k, n in enumerate(maps_femm.FROZEN_COLUMNS)}
        for label in order:
            z = np.load(folder / f"{label}.npz")
            R, delta = z["rows"], float(z["test_current"])
            i_d, i_q = float(z["i_d"]), float(z["i_q"])
            pair = lambda name: np.column_stack([R[:, cols[name + "_psi_d"]], R[:, cols[name + "_psi_q"]]])
            ref, f0, fd, fq, fr = pair("ref"), pair("F0"), pair("Fd"), pair("Fq"), pair("Fref")
            native = np.column_stack([R[:, cols["native_psi_d"]], R[:, cols["native_psi_q"]]])
            Lfp = np.stack([(fd - f0) / delta, (fq - f0) / delta], axis=2)             # (n_pos, 2, 2): column j = response to the test current on axis j
            pred = f0 + Lfp @ np.array([i_d, i_q])
            scale = np.hypot(ref[:, 0], ref[:, 1])
            close.append({"point": label, "equivalent vs native": float(np.max(np.abs(ref - native).max(axis=1) / scale)),
                          "Fref vs loaded": float(np.max(np.abs(fr - ref).max(axis=1) / scale)), "F0 + L i vs loaded": float(np.max(np.abs(pred - ref).max(axis=1) / scale))})
            psi_pm_fp, L = f0.mean(axis=0), Lfp.mean(axis=0)
            mdl = fm.psi_dq(i_d, i_q)
            Lm = fm.inductances(i_d, i_q)
            rows.append({"point": label, "i_d [A]": i_d, "i_q [A]": i_q, "psi_PM,fp d [mWb]": 1e3 * psi_pm_fp[0], "psi_PM,fp q [mWb]": 1e3 * psi_pm_fp[1],
                         "psi_PM,fp d / no-load": psi_pm_fp[0] / fm.psi_pm0, "L_dd,fp [mH]": 1e3 * L[0, 0], "L_qq,fp [mH]": 1e3 * L[1, 1],
                         "L_dq,fp [mH]": 1e3 * L[0, 1], "L_qd,fp [mH]": 1e3 * L[1, 0], "L_d,app [mH]": 1e3 * Lm["L_d_app"], "L_q,app [mH]": 1e3 * Lm["L_q_app"],
                         "L_dd,inc [mH]": 1e3 * Lm["L_dd"], "L_qq,inc [mH]": 1e3 * Lm["L_qq"], "L_dq,inc [mH]": 1e3 * Lm["L_dq"], "L_qd,inc [mH]": 1e3 * Lm["L_qd"],
                         "psi_d (maps) [mWb]": 1e3 * mdl[0], "psi_d (frozen runs, 6-position mean) [mWb]": 1e3 * ref[:, 0].mean()})
        t = pd.DataFrame(rows)
        pd.set_option("display.width", 250)
        pd.set_option("display.max_columns", 40)
        c = pd.DataFrame(close)
        ok = bool((c[["Fref vs loaded", "F0 + L i vs loaded"]].max().max()) < 1e-5 and c["equivalent vs native"].max() < 1e-9)
        self.say("10. closure of the frozen method, worst over the 6 positions of each point, relative to |psi| (limits: Fref against the loaded solution and psi_PM,fp + L_fp i against "
                 f"it 1e-5, the fill-1 equivalent model against the native one 1e-9): {self.verdict(ok)}")
        self.say(c.to_string(index=False, float_format=lambda v: f"{v:.2e}"))
        self.say("11. the split of the loaded flux linkage: psi = psi_PM,fp + L_fp i (frozen permeability, FEMM 5.0, fine mesh, mean over 6 positions) against the maps' apparent "
                 "and incremental inductances at the same point")
        self.say(t.T.to_string(float_format=lambda v: f"{v:.4g}"))
        self.out["frozen"] = {"closure_pass": ok, "closure": close, "table": rows}
        self.save("frozen")

    # ============================================================================================== the iron loss
    def run_loss(self) -> None:
        self.say("=" * 30 + " the iron loss " + "=" * 30)
        pd.set_option("display.width", 250)
        pd.set_option("display.max_columns", 40)
        fit = iron_loss.load_fit()
        loss_dir = OUT / "loss"
        mesh = np.load(loss_dir / "mesh.npz")
        summary = iron_loss.build_iron_map(loss_dir, fit, OUT / "iron_map.npz", self.I_R)
        imap = iron_loss.IronLossMap(OUT / "iron_map.npz")
        names = iron_loss.REGIONS
        self.say(f"12. the iron-loss map: {summary['points']} points ({len(summary['i_d'])} x {len(summary['i_q'])}) -> iron_map.npz; the fit: k_h {fit['kh']:.5f}, alpha {fit['alpha']:.4f}, "
                 f"k_x {fit['kx']:.4e}, k_e {fit['ke']:.4e} (datasheet gauge {iron_loss.D_LAM * 1e3:.3f} mm; the study's {iron_loss.D_LAM_CAD * 1e3:.3f} mm scales the classical term by "
                 f"{iron_loss.eddy_scale(iron_loss.D_LAM_CAD * 1e3):.4f})")
        # the mesh and the masses
        area = mesh["steel_area_mm2"]
        regions = mg.build(design.params(self.ctx.overrides))
        visible = mg.visible(regions)                    # the steel regions as they show once the slots and the magnets are painted over them (the saved mesh.npz has the unpainted ones)
        poly = [float(sum(visible[r.name].area for r in regions if r.material == "M-19 Steel" and r.side != "rotor")),
                float(sum(visible[r.name].area for r in regions if r.material == "M-19 Steel" and r.side == "rotor"))]
        mass = iron_loss.RHO_STEEL * iron_loss.LAM_FILL * area * 1e-6 * float(mesh["depth_mm"]) * 1e-3 * int(mesh["sectors"])
        self.say(f"    the mesh: {int(mesh['nodes'])} nodes, {int(mesh['triangles'])} triangles; steel area in one sector: teeth {area[0]:.1f}, yoke {area[1]:.1f}, rotor {area[2]:.1f} mm^2; "
                 f"stator {area[0] + area[1]:.1f} against its polygon {poly[0]:.1f} ({100 * ((area[0] + area[1]) / poly[0] - 1):+.3f} %), rotor {area[2]:.1f} against {poly[1]:.1f} "
                 f"({100 * (area[2] / poly[1] - 1):+.3f} %); the slots' deepest radius {float(mesh['r_split']):.3f} mm (teeth inside, yoke outside)")
        self.say(f"    steel mass of the machine (20 mm stack, fill {iron_loss.LAM_FILL}, {iron_loss.RHO_STEEL:.0f} kg/m^3): teeth {mass[0]:.3f} kg, yoke {mass[1]:.3f} kg, rotor {mass[2]:.3f} kg, "
                 f"total {mass.sum():.3f} kg")
        # the rated point and the table of the loss against speed
        i_d1, i_q1 = design.currents(self.I_R, self.ctx.gamma_R)
        fm = maps.FluxMaps.load()
        pts = [("no load", 0.0, 0.0)] + [(f"MTPA {m} x rated", *fm.mtpa(m * self.I_R)[1:3]) for m in (1, 2, 3)]
        out_tab = []
        self.say(f"13. the iron loss of the whole machine (20 mm stack, 10-turn placeholder winding, which sets the ampere-turns only), by region, at the rated speed {RPM_RATED:.0f} rpm "
                 f"(f = {PP * RPM_RATED / 60:.0f} Hz) and the top speed {RPM_TOP:.0f} rpm (f = {PP * RPM_TOP / 60:.0f} Hz); both lamination thicknesses, the CAD's 0.3 mm first (the headline)")
        for lam in LAMINATIONS_MM:
            rows = []
            for name, i_d, i_q in pts:
                for rpm in (RPM_RATED, 8000.0, RPM_TOP):
                    r = imap(i_d, i_q, rpm, lamination_mm=lam)
                    rows.append({"lamination [mm]": lam, "point": name, "rpm": rpm, "teeth [W]": r["teeth"], "yoke [W]": r["yoke"], "rotor [W]": r["rotor"], "total [W]": r["total"]})
            t = pd.DataFrame(rows)
            out_tab += rows
            self.say(t.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
        # the three terms at the rated MTPA point
        S = imap.sums(i_d1, i_q1)
        rows = []
        for r, nm in enumerate(names):
            for f_rpm in (RPM_RATED, RPM_TOP):
                f = PP * f_rpm / 60
                for lam in LAMINATIONS_MM:
                    k = iron_loss.eddy_scale(lam)
                    ph, pe, px = S[r, 0] * f, k * S[r, 1] * f ** 2, S[r, 2] * f ** 1.5
                    rows.append({"region": nm, "rpm": f_rpm, "lamination [mm]": lam, "hysteresis [W]": float(ph), "classical [W]": float(pe), "excess [W]": float(px),
                                 "total [W]": float(ph + pe + px)})
        self.say(f"    the three terms at the rated point ({i_d1:.2f}, {i_q1:.2f}) A:")
        self.say(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        self.say("    the sums S_h [W/Hz], S_e [W/Hz^2], S_x [W/Hz^1.5] per region at the rated point (datasheet gauge; S_e x 0.710 at 0.3 mm): "
                 + "; ".join(f"{nm}: {S[r, 0]:.4g}, {S[r, 1]:.4g}, {S[r, 2]:.4g}" for r, nm in enumerate(names)))
        self.out["loss"] = {"summary": summary, "table": out_tab, "terms_at_rated": rows, "mass_kg": mass.tolist(), "fit": fit}
        # the checks
        self.check_loss_period(summary)
        self.check_loss_positions(fit, mesh)
        self.check_loss_interp(imap, fit, loss_dir / "mesh.npz")
        self.check_loss_loo(imap)
        self.check_loss_fundamental(fit, mesh)
        self.check_loss_skin(fit, mesh)
        self.save("loss")

    def check_loss_fundamental(self, fit: dict, mesh) -> None:
        """A sanity check of the whole chain by an independent, simpler route: the datasheet's own formula applied to the FUNDAMENTAL of every triangle's flux density alone (each triangle a
        sinusoid of its fundamental's vector amplitude, the textbook calculation) against the waveform-based sums. The full result must be at least the fundamental-only one (the harmonics
        only add loss), by a modest factor in the stator and by a large one in the rotor, whose fundamental is a steady field (nothing at all in the sinusoidal reading)."""
        rows = []
        depth, sectors = float(mesh["depth_mm"]), int(mesh["sectors"])
        mass_e = iron_loss.RHO_STEEL * iron_loss.LAM_FILL * mesh["area_mm2"] * 1e-6 * depth * 1e-3 * sectors
        reg = mesh["region"]
        for key in ("p+0.00_0.00", "p-0.50_1.00", "p-1.00_1.00", "p-2.00_2.00"):
            path = OUT / "loss" / f"point_{key}.npz"
            if not path.exists():
                continue
            z = np.load(path)
            B1 = z["Ak"][:, 0].astype(float) / iron_loss.LAM_FILL
            S = iron_loss.point_sums(path, OUT / "loss" / "mesh.npz", fit)
            for rpm in (RPM_RATED, RPM_TOP):
                f = PP * rpm / 60
                full = self.region_watts(S, rpm, lam=iron_loss.D_LAM * 1e3)
                for r, rn in enumerate(iron_loss.REGIONS):
                    fund = float((mass_e[reg == r] * iron_loss.bertotti(f, B1[reg == r], fit["kh"], fit["alpha"], fit["kx"], fit["ke"])).sum())
                    rows.append({"point": key, "rpm": rpm, "region": rn, "fundamental only [W]": fund, "waveform-based [W]": float(full[r]), "ratio": float(full[r] / fund) if fund > 0 else float("inf")})
        t = pd.DataFrame(rows)
        stator = t[t.region != "rotor"]
        ok = bool((stator.ratio >= 1.0 - 1e-9).all())
        self.say(f"17b. the chain against the textbook route: the datasheet formula on each triangle's fundamental alone against the waveform-based loss (datasheet gauge, whole machine). The full "
                 f"loss must not be below the fundamental-only one in the stator: {self.verdict(ok)}")
        self.say(t.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        self.out.setdefault("loss", {})["fundamental_only"] = {"pass": ok, "rows": rows}

    def check_loss_period(self, summary: dict) -> None:
        """The 72 positions of every loss point cover 360 degrees: the mean over the full period against the mean over the first 60 (12 positions, 5 degrees apart),
        and against the flux grid's own mean at the same node (an independent run of the same point), and T(theta + 60) against T(theta)."""
        built = np.load(OUT / "flux_map.npz")
        ids, iqs = built["id_A"], built["iq_A"]
        rows = []
        for path in sorted((OUT / "loss").glob("point_*.npz")):
            z = np.load(path)
            n = int(z["n"])
            th = np.asarray(z["positions"])
            sel60 = th < 60 - 1e-9
            pd_, pq_, T = z["psi_d"], z["psi_q"], z["T_age"]
            a, b = int(np.argmin(np.abs(ids - float(z["i_d"])))), int(np.argmin(np.abs(iqs - float(z["i_q"]))))
            scale = math.hypot(pd_.mean(), pq_.mean())
            rows.append({"i_d": float(z["i_d"]), "i_q": float(z["i_q"]), "psi_360_vs_60": max(abs(pd_.mean() - pd_[sel60].mean()), abs(pq_.mean() - pq_[sel60].mean())) / scale,
                         "T_360_vs_60": abs(T.mean() - T[sel60].mean()), "T": float(T.mean()),
                         "psi_vs_flux_grid": max(abs(pd_.mean() - built["psi_d"][a, b]), abs(pq_.mean() - built["psi_q"][a, b])) / scale,
                         "T_vs_flux_grid": abs(T.mean() - built["T_age"][a, b]),
                         "T_shift60_rms": float(np.sqrt(np.mean((np.roll(T, -n // 6) - T) ** 2)))})
        t = pd.DataFrame(rows)
        ok = bool((t.psi_360_vs_60 <= TOL_PERIOD_PSI).all() and (t.T_360_vs_60 <= TOL_PERIOD_T_REL * t["T"].abs() + TOL_PERIOD_T_ABS).all())
        self.say(f"14. the period: over {len(t)} loss points the mean over the full 360 degrees against the mean over the first 60 (limits: psi {TOL_PERIOD_PSI:.0e} of |psi|, T "
                 f"{TOL_PERIOD_T_REL:.0e} of T + {1e3 * TOL_PERIOD_T_ABS:.0f} mN.m): {self.verdict(ok)}")
        self.say(f"    worst psi {t.psi_360_vs_60.max():.2e}, worst T {1e3 * t.T_360_vs_60.max():.2f} mN.m; the same node in the flux grid (30 positions over 60 degrees, a different run): psi "
                 f"{t.psi_vs_flux_grid.max():.2e} of |psi|, T {1e3 * t.T_vs_flux_grid.max():.2f} mN.m at most; T(theta + 60 deg) against T(theta): rms difference at most "
                 f"{1e3 * t.T_shift60_rms.max():.2f} mN.m (median {1e3 * t.T_shift60_rms.median():.2f})")
        self.out.setdefault("loss", {})["period"] = {"pass": ok, "psi_360_vs_60_max": float(t.psi_360_vs_60.max()), "T_360_vs_60_max_Nm": float(t.T_360_vs_60.max()),
                                                     "psi_vs_flux_grid_max": float(t.psi_vs_flux_grid.max()), "T_vs_flux_grid_max_Nm": float(t.T_vs_flux_grid.max()),
                                                     "T_shift60_rms_max_Nm": float(t.T_shift60_rms.max())}

    def region_watts(self, S: np.ndarray, rpm: float, lam: float = 0.3) -> np.ndarray:
        """The three regions' watts of a (3, 3) sums array at a speed and lamination."""
        f = PP * rpm / 60
        k = iron_loss.eddy_scale(lam)
        return S[:, 0] * f + k * S[:, 1] * f ** 2 + S[:, 2] * f ** 1.5

    def check_loss_positions(self, fit: dict, mesh) -> None:
        """What 72 positions over 360 degrees do to the sums, from the two nodes solved at 144, and which harmonics carry them."""
        rows, spec = [], {}
        for key in ("fine-1.00_1.00", "fine-2.00_2.00"):
            p144 = OUT / "losscheck" / f"{key}.npz"
            if not p144.exists():
                continue
            a, b = float(key[4:9]), float(key[10:])
            node72 = OUT / "loss" / f"point_p{a:+.2f}_{b:.2f}.npz"
            S144 = iron_loss.point_sums(p144, OUT / "loss" / "mesh.npz", fit)
            S72 = iron_loss.point_sums(node72, OUT / "loss" / "mesh.npz", fit)
            for term, name in enumerate(iron_loss.TERMS):
                for r, rn in enumerate(iron_loss.REGIONS):
                    rows.append({"point": key, "term": name, "region": rn, "S_72": S72[r, term], "S_144": S144[r, term], "72 vs 144": S72[r, term] / S144[r, term] - 1})
            w72, w144 = self.region_watts(S72, RPM_TOP), self.region_watts(S144, RPM_TOP)
            w72r, w144r = self.region_watts(S72, RPM_RATED), self.region_watts(S144, RPM_RATED)
            rows.append({"point": key, "term": "watts at top speed", "region": "total", "S_72": w72.sum(), "S_144": w144.sum(), "72 vs 144": w72.sum() / w144.sum() - 1})
            for rn, x, y, xr, yr in zip(iron_loss.REGIONS, w72, w144, w72r, w144r):
                rows.append({"point": key, "term": "watts at top speed", "region": rn, "S_72": x, "S_144": y, "72 vs 144": x / y - 1})
                rows.append({"point": key, "term": "watts at rated speed", "region": rn, "S_72": xr, "S_144": yr, "72 vs 144": xr / yr - 1})
            # which harmonics carry S_e (the classical sum), by region, from the 144-position run
            z = np.load(p144)
            Ak = z["Ak"].astype(float) / iron_loss.LAM_FILL
            reg = np.load(OUT / "loss" / "mesh.npz")["region"]
            area = np.load(OUT / "loss" / "mesh.npz")["area_mm2"]
            k = np.arange(1, Ak.shape[1] + 1)
            share = {}
            for r, rn in enumerate(iron_loss.REGIONS):
                e = (Ak[reg == r] ** 2 * area[reg == r][:, None]).sum(axis=0) * k ** 2
                share[rn] = {"k=1": float(e[k == 1].sum() / e.sum()), "k=2..6": float(e[(k >= 2) & (k <= 6)].sum() / e.sum()), "k=7..18": float(e[(k >= 7) & (k <= 18)].sum() / e.sum()),
                             "k=19..36": float(e[(k >= 19) & (k <= 36)].sum() / e.sum()), "k=37..72": float(e[k >= 37].sum() / e.sum()),
                             "k=9 (slot passing)": float(e[k == 9].sum() / e.sum()), "k=18": float(e[k == 18].sum() / e.sum())}
            spec[key] = share
        t = pd.DataFrame(rows)
        worst = {}
        for term in ("h", "e", "x", "watts at top speed", "watts at rated speed"):
            s = t[t.term == term]
            worst[term] = float(s["72 vs 144"].abs().max())
        sel = t[(t.term.str.startswith("watts")) & (t.region != "total")]
        ok = bool(sel["72 vs 144"].abs().max() <= TOL_LOSS_N)
        self.say(f"15. positions of the loss runs: the sums from 72 positions against 144 at two nodes ((-1, 1) and (-2, 2) x rated): the loss of each region at {RPM_RATED:.0f} and {RPM_TOP:.0f} rpm "
                 f"(limit {100 * TOL_LOSS_N:g} %): {self.verdict(ok)}")
        self.say(t[t.term.isin(["h", "e", "x"])].pivot_table(index=["point", "region"], columns="term", values="72 vs 144").to_string(float_format=lambda v: f"{100 * v:+.2f} %"))
        self.say(t[t.term.str.startswith("watts")].to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        for key, share in spec.items():
            self.say(f"    the classical sum S_e by harmonic order k of the electrical frequency, from the 144-position run at {key}:")
            self.say(pd.DataFrame(share).to_string(float_format=lambda v: f"{100 * v:.1f} %"))
        self.out.setdefault("loss", {})["positions"] = {"pass": ok, "worst_relative": worst, "rows": rows, "spectrum": spec}

    def check_loss_interp(self, imap: iron_loss.IronLossMap, fit: dict, mesh_path: Path) -> None:
        folder = OUT / "losscheck"
        rows = []
        for path in sorted(folder.glob("*.npz")):
            key = path.stem
            if key.startswith("fine"):
                continue
            z = np.load(path)
            S = iron_loss.point_sums(path, mesh_path, fit)
            i_d, i_q = float(z["i_d"]), float(z["i_q"])
            kind = "held-out (off-grid)" if key.startswith("held") else "mirror (i_q < 0)"
            for rpm in (RPM_RATED, RPM_TOP):
                want = self.region_watts(S, rpm)
                got = imap(i_d, i_q, rpm)
                for r, rn in enumerate(iron_loss.REGIONS):
                    rows.append({"kind": kind, "point": key, "i_d [A]": i_d, "i_q [A]": i_q, "rpm": rpm, "region": rn, "FEMM [W]": want[r], "map [W]": got[rn],
                                 "difference": got[rn] / want[r] - 1})
                rows.append({"kind": kind, "point": key, "i_d [A]": i_d, "i_q [A]": i_q, "rpm": rpm, "region": "total", "FEMM [W]": want.sum(), "map [W]": got["total"],
                             "difference": got["total"] / want.sum() - 1})
        t = pd.DataFrame(rows)
        held = t[(t.kind.str.startswith("held")) & (t.region != "total")]
        mirror = t[(t.kind.str.startswith("mirror")) & (t.region != "total")]
        ok_i = bool(held.difference.abs().max() <= TOL_LOSS_INTERP)
        ok_m = bool(mirror.difference.abs().max() <= TOL_LOSS_MIRROR)
        self.say(f"16. the loss map against FEMM at points it was not built from (the map uses |i_q|; loss at 0.3 mm): off-grid points: each region within {100 * TOL_LOSS_INTERP:g} %: {self.verdict(ok_i)}; "
                 f"negative i_q (the mirror assumption) within {100 * TOL_LOSS_MIRROR:g} %: {self.verdict(ok_m)}")
        self.say(t.to_string(index=False, float_format=lambda v: f"{v:.4g}"))
        self.say(f"    worst off-grid difference {100 * held.difference.abs().max():.2f} % (region), {100 * t[(t.kind.str.startswith('held')) & (t.region == 'total')].difference.abs().max():.2f} % (total); "
                 f"worst mirror difference {100 * mirror.difference.abs().max():.2f} % (region), {100 * t[(t.kind.str.startswith('mirror')) & (t.region == 'total')].difference.abs().max():.2f} % (total)")
        self.out.setdefault("loss", {})["interpolation_and_mirror"] = {"pass_interp": ok_i, "pass_mirror": ok_m, "rows": t.to_dict(orient="records")}

    def check_loss_loo(self, imap: iron_loss.IronLossMap) -> None:
        """A second look at the interpolation, from the grid alone: leave one whole interior line of the 7 x 7 grid out (a row of i_d or a column of i_q), build the map from the rest, and ask it
        for the left-out nodes, which FEMM solved. The gap it must bridge is double the grid's, so this is harder than the real use; it is a bound, with 210 test values per speed."""
        from scipy.interpolate import RectBivariateSpline
        S, ids, iqs = imap.S, imap.id_A, imap.iq_A
        errs = {RPM_RATED: [], RPM_TOP: []}
        for axis in (0, 1):
            n = S.shape[2 + axis]
            for drop in range(1, n - 1):
                keep = [i for i in range(n) if i != drop]
                ids_k, iqs_k = (ids[keep], iqs) if axis == 0 else (ids, iqs[keep])
                if iqs_k[0] != 0:
                    continue
                iq_full = np.concatenate([-iqs_k[:0:-1], iqs_k])
                pred = np.zeros_like(S[:, :, [drop], :] if axis == 0 else S[:, :, :, [drop]])
                for r in range(3):
                    for t in range(3):
                        g = S[r, t][keep, :] if axis == 0 else S[r, t][:, keep]
                        full = np.concatenate([g[:, :0:-1], g], axis=1)
                        spl = RectBivariateSpline(ids_k, iq_full, full, kx=3, ky=3, s=0)
                        A, B = np.meshgrid(ids[[drop]] if axis == 0 else ids, iqs if axis == 0 else iqs[[drop]], indexing="ij")
                        pred[r, t] = np.maximum(spl.ev(A.ravel(), B.ravel()).reshape(A.shape), 0)
                truth = S[:, :, [drop], :] if axis == 0 else S[:, :, :, [drop]]
                for rpm in errs:
                    for a in range(pred.shape[2]):
                        for b in range(pred.shape[3]):
                            errs[rpm] += list(self.region_watts(pred[:, :, a, b], rpm) / self.region_watts(truth[:, :, a, b], rpm) - 1)
        out = {}
        self.say("16b. the loss map's interpolation from the grid alone (leave one interior row or column out, rebuild, ask for it; a doubled gap, so a bound): relative error of a region's watts")
        for rpm, e in errs.items():
            e = np.abs(np.array(e))
            out[f"{rpm:.0f}"] = {"n": int(len(e)), "median": float(np.median(e)), "p95": float(np.quantile(e, 0.95)), "max": float(e.max())}
            self.say(f"    {rpm:.0f} rpm: median {100 * np.median(e):.2f} %, 95th percentile {100 * np.quantile(e, 0.95):.2f} %, largest {100 * e.max():.1f} % (n = {len(e)})")
        self.out.setdefault("loss", {})["leave_one_line_out"] = out

    def check_loss_skin(self, fit: dict, mesh) -> None:
        """How much of the classical loss sits at harmonics where the sheet's skin effect (not in the model) would cut it: from the node nearest the rated point."""
        key = "p-0.50_1.00"
        path = OUT / "loss" / f"point_{key}.npz"
        z = np.load(path)
        reg, area = mesh["region"], mesh["area_mm2"]
        Ak = z["Ak"].astype(float) / iron_loss.LAM_FILL
        k = np.arange(1, Ak.shape[1] + 1)
        f1 = PP * RPM_TOP / 60
        rows = []
        for r, rn in enumerate(iron_loss.REGIONS):
            w = (Ak[reg == r] ** 2 * area[reg == r][:, None]).sum(axis=0) * k ** 2          # the classical sum's weight of each harmonic
            for lam in LAMINATIONS_MM:
                for mu in (1000.0, 3000.0):
                    fac = iron_loss.eddy_skin_factor(lam * 1e-3, k * f1, mu)
                    rows.append({"region": rn, "lamination [mm]": lam, "mu_r": mu, "exact / classical (weighted, all harmonics)": float((w * fac).sum() / w.sum()),
                                 "share of S_e above 2 kHz": float(w[k * f1 > 2000].sum() / w.sum())})
        t = pd.DataFrame(rows)
        self.say(f"17. what the classical term leaves out: the skin effect of the sheet at {RPM_TOP:.0f} rpm (f = {f1:.0f} Hz), weighting each harmonic k f by its share of S_e at the node "
                 f"(-0.5, 1) x rated (mu_r a stated assumption, 1,000 or 3,000; the datasheet ends at 2 kHz):")
        self.say(t.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
        self.out.setdefault("loss", {})["skin"] = t.to_dict(orient="records")


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    rep = Report()
    todo = ["flux", "cross", "frozen", "loss"] if args in ([], ["all"]) else args
    for s in todo:
        if s == "flux":
            rep.run_flux()
        elif s == "cross":
            rep.run_cross()
        elif s == "frozen":
            rep.run_frozen()
        elif s == "loss":
            rep.run_loss()
        else:
            raise SystemExit(f"unknown section {s!r}: flux, cross, frozen, loss or all")
    return 0


if __name__ == "__main__":
    sys.exit(main())
