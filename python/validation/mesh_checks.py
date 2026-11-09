"""mesh_checks.py: is the difference between FEMM 4.2 and FEMM 5.0 the mesh? Solve the model again at other densities.

Why this exists. compare.py finds a few numbers that moved by more than the plan allows (a cogging torque, a small torque
harmonic, a flux linkage where the d-current has cancelled most of the magnet flux, ...). The model is the same on both FEMMs
(same_model.py proves it), and FEMM 5.0's mesher makes fewer, larger triangles from the same mesh sizes: 9,580 nodes where
FEMM 4.2 made 11,644 (coarse), 17,934 where it made 27,939 (fine). The question is whether that alone explains the moved
numbers. Two checks answer it, and neither changes the model (the mesh sizes of the labels are multiplied by one factor f,
which is what a finer mesh is; the geometry, the materials and the currents are the model's own).

  1. A MATCHED DENSITY. Find the f that gives FEMM 5.0 FEMM 4.2's node count, solve the same points again, and compare with
     FEMM 4.2's numbers. If the difference was the mesh, it shrinks. This is not tuning the model to match FEMM 4.2: nothing
     is adjusted until the numbers agree, f is chosen by the node count alone, before any torque is looked at.
  2. A LADDER. The same quantities at five or six densities. A number that is converged does not move along the ladder. A
     number that is still moving is mesh noise, and the 4.2-5.0 difference should be no bigger than that movement.

Run (from the Python to FEMM folder, one FEMM-using thing at a time, about 23 minutes with 8 workers, 22 min 45 s on 2026-11-05,
about 7,800 solves):
    uv run python validation/mesh_checks.py femm42-out              # all of it
    uv run python validation/mesh_checks.py femm42-out --only calibrate,slice,states,post42
OLD_OUT (the first argument) is FEMM 4.2's run, read only: its numbers are the reference. femm42-out is a checked copy (every
file's SHA-256 equal, femm42-out.sha256) of the `out` folder that FEMM 4.2 left in D:/dev/scripts/automation/proc/out.
Outputs go to out/validation/mesh/ (CSV files and mesh_convergence.png, which this script also copies to validation/, beside
it, so that the copy cannot go stale). Every worker opens its own FEMM (femm_session) and every job has files of its own, because
two FEMMs on one file name crash, hang or swap answers (learn/04_workers.py, parts 3 and 4: part 3 shows the failure, part 4 the
rule that prevents it).

"The report" below is the project note "4.2 vs 5.0"; its section 5 holds
the numbers this script prints, number for number, and section 2.12 the flags (F2, F3, F5, F6, F7) that it backs with checks.

The thirteen parts, in the order they run (--only takes any of them, comma-separated; the later ones read the files of the
earlier ones, so run calibrate and slice first on a clean out folder):
  calibrate  find f for each target node count (calibration.csv). Targets: FEMM 4.2's two counts, and a ladder around them.
  slice      the barrier study's no-load, cogging, torque-vs-angle and ripple points at every density (slice_ladder.csv,
             slice_summary.csv, gap_*.csv), and the four check_femm points (check_femm_ladder.csv)
  states     the flux explorer's 60 states of five rotors at coarse densities, against FEMM 4.2's states (states_ladder.csv)
  post42     FEMM 5.0's post-processor on FEMM 4.2's own six saved solutions, against what FEMM 4.2 recorded for them (post42.csv)
  post42full the same for the 48 saved solutions of the hand-drawn 360-degree model (post42_full.csv): is the full model's torque offset
             made by the post-processor, or already in the solution?
  fullmodel  the hand-drawn 360-degree model with the rotor turned by tiny extra angles: how much a re-mesh alone moves its torque, and how
             much of the difference from FEMM 4.2 is one common offset of the 24-position mean (full_model.csv)
  fulldense  the hand-drawn 360-degree model at 96 positions per electrical cycle instead of 24: is that offset the sampling (full_dense.csv)
  noise      the slice ladder's quantities again with ONE parameter moved by one unit in the last place: how much a number moves when the
             mesh changes and the physics does not (noise_ladder.csv, noise_summary.csv)
  flagged    the two flags that were argued by analogy, solved: the 20 A states of the torque-vs-angle grid (flagged_states.csv)
             and the no-band rotor's ripple harmonics (flagged_noband.csv), each on a ladder and with one bit moved
  timing     one FEMM, one job at a time: seconds per solve at four node counts (timing.csv), for the speed-up split
  candidate  the optimizer's candidate whose best current angle moved most, solved again as it was and with one bit of a parameter changed
  ulp        all 48 feasible designs of the optimizer, solved again as they were (they must equal the sweep's rows) and with one bit moved
  figure     mesh_convergence.png

What to expect (2026-11-04 and 05, this PC): the report, section "5. The node counts, and whether the mesh explains the
differences" (5.1 to 5.7), and 2.12 for the flags F2, F3, F5 and F6. The short version: FEMM 5.0 at FEMM 4.2's node count
gives the cogging torque within 0.2 % of FEMM 4.2's (at its own it is 7 % high); the mean torques and the magnet flux come into
line only at 45k to 70k nodes; the slice ladder, the one-bit noise, the saved solutions and the full model each print a table
that section 5 reproduces number for number. Which part prints which table of the report:
  slice, check_femm ladder   5.3 (the ladder), the numbers of 2.12 F1;           states    2.12 F7;       timing   4 (the speed)
  post42                     5.2 (FEMM 5.0 reads FEMM 4.2's solutions);           noise     5.5;           flagged  5.7 (and 2.12 F2, F3)
  post42full, fullmodel, fulldense   5.6 (the hand-drawn full 360-degree model) and 2.12 F5;                           candidate, ulp   2.12 F6
  figure                     the picture of 5.3 (mesh_convergence.png)

Self-check (the answers follow).
  1. FEMM 5.0 at FEMM 4.2's node count gives a cogging torque that is within 2 % of FEMM 4.2's, while at its own node count
     it is 7 % off. What did f do, and what did it not do to the model?
  2. A quantity changes by 0.4 % along the ladder and FEMM 4.2 differs from FEMM 5.0 by 0.3 %. Is that a difference between the
     programs, or the same program at two densities?
  3. Why does each worker write to its own file names, even though the jobs are independent?
  4. Moving a parameter by 1e-15 mm changes a torque by 0.04 %. Is the torque "wrong" by 0.04 %? What does it say about a difference
     of 0.2 % between two programs?
  5. The mean of 24 no-load torques of the full model is +1.69 mN.m on FEMM 4.2 and -0.73 on FEMM 5.0. Give two reasons why the
     mean of samples of a periodic waveform need not be the waveform's mean, and say which one `fulldense` tests.

Answers.
  1. f multiplied every label's mesh size by one number, so the mesher cut the same drawing into more (f < 1) or fewer triangles.
     It changed nothing in the geometry, the materials or the currents, and f was chosen by the node count alone, before any torque
     was looked at, so nothing was tuned toward FEMM 4.2's answer. That the cogging torque then agrees says that the 7 % came
     from the mesh's density, not from the model or the solver.
  2. Not decidable from those two numbers alone. If the quantity moves by 0.4 % between two densities of FEMM 5.0 itself, then 0.3 %
     is inside the movement of ONE program, and nothing says that the programs differ. The ladder is there to give that
     comparison: a difference between the programs counts only where it is bigger than the movement along the ladder (and bigger
     than the one-bit noise).
  3. A FEMM saves its model to its file before every solve, and every solver then rewrites the .ans. Two FEMMs on one name read
     each other's half-written files: they crash, hang, or one job is handed the other job's answer (the "0.0 glitch" of 10-28 was
     this). That the jobs are independent in the mathematics does not make their files independent on the disk;
     learn/04_workers.py, part 3, makes it happen on purpose and part 4 shows the private names that prevent it.
  4. No. One unit in the last place of a parameter is no change of the motor, yet FEMM meshes the slightly different drawing
     differently, so the torque moves by 0.04 %: that is the mesh's own noise at that quantity. A 0.2 % difference between two
     programs is five times it, so it is not noise, but it is not a fault of either program either: it is what the two meshes of
     the same motor differ by.
  5. (a) Aliasing: 24 positions per electrical cycle cannot tell the 24th harmonic of the torque (nor the 48th ...) from a constant,
     so the sample mean is the mean plus those harmonics, with a sign that depends on where the grid starts. (b) The waveform's
     own mean is not zero: the torque is the slope of an energy, and its mean over one cycle is 0 only if the rotor is the same
     rotor again after the cycle, which the hand-drawn rotor is not quite (two magnets point 5 degrees off the pattern); and a
     solver's energy has its own error, which differs between two meshes. `fulldense` tests (a): four 24-position grids, shifted
     against each other, give means that differ by only 0.9 mN.m at no load, and the 96-position mean differs from FEMM 5.0's own
     24-position mean by 0.4 (0.37 of the 2.42 between the programs at no load, 0.51 of 3.76 at 10 A), so the sampling is about a
     seventh of it. (b) is what is left: it cannot be separated without a third, trusted answer, and FEMM 4.2 cannot be run again.
"""
from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
import shutil
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
P = HERE.parent
for _folder in (str(P), str(HERE)):
    if _folder not in sys.path:
        sys.path.insert(0, _folder)
import femm_session  # noqa: E402
import flux_explorer_sweep as sweep  # noqa: E402
import motorgeo as mg  # noqa: E402
import paths  # noqa: E402  (config.toml: the reference model, and every other path that is not P's own)
from compare import README_CHECK_FEMM, README_NODES, harmonic_amplitudes  # noqa: E402  (FEMM 4.2's numbers live in one place)

OUT = P / "out" / "validation" / "mesh"
WORKERS = 8
PARTS = ("calibrate", "slice", "states", "post42", "post42full", "fullmodel", "fulldense", "noise", "flagged", "timing", "candidate", "ulp",
         "figure")                                              # the order they run in
# the ladder: (tag, smart mesh on?, target nodes). "matched" is FEMM 4.2's own node count (README_NODES), "natural" FEMM 5.0's (f = 1)
CONFIGS = [("c-8k", False, 8000), ("c-natural", False, None), ("c-matched", False, README_NODES["coarse"]), ("c-17k", False, 17000),
           ("c-28k", False, 28000),
           ("f-natural", True, None), ("f-matched", True, README_NODES["fine"]), ("f-45k", True, 45000), ("f-70k", True, 70000)]
STATES_ROTORS = (0, 7, 9, 10, 12)                               # the reference design, and the rotors of the flux explorer that compare.py flags
MEAN_POS = range(0, 60, 10)                                     # barrier_study's six rotor positions
GAMMA_MTPA_42 = 28.105146968413635                              # FEMM 4.2's best current angle with the band (its phase2.csv)
IQ_ID_42 = (-4.7109112213837685, 8.82084550733321)              # the (id, iq) FEMM 4.2 solved its ripple sweep at
GAP_DEG = list(np.arange(0, 90, 0.25))


def current_vector(amps: float, gamma_deg: float) -> tuple[float, float]:
    """(i_d, i_q) for a current of `amps` at the angle gamma from the q axis toward -d (0 = all q), as the studies write it."""
    return -amps * math.sin(math.radians(gamma_deg)), amps * math.cos(math.radians(gamma_deg))


# ----------------------------------------------------------------------------------------------- the model at density f
def scaled_regions(p, f):
    """The model's regions with every label's mesh size multiplied by f. Nothing else about them changes."""
    return [dataclasses.replace(r, mesh=r.mesh * f) for r in mg.build(p)]


def open_model(p, f, fine, path):
    """Inside an open FEMM: write the model with mesh sizes x f, and open it for solving. Returns the FemmModel."""
    import femm
    mg.to_femm(p, scaled_regions(p, f), str(path), fine=fine)
    femm.mi_close()
    return mg.FemmModel(p, str(path))


def count_nodes(model) -> int:
    """Solve once and ask FEMM how many nodes the mesh has."""
    import femm
    femm.mi_analyze(1)
    femm.mi_loadsolution()
    n = int(femm.mo_numnodes())
    femm.mo_close()
    return n


def calibrate_one(task):
    """One worker: find f for a target node count. Returns (tag, f, nodes)."""
    import femm
    tag, fine, target = task
    p = mg.Params()
    path = OUT / "fem" / f"cal_{tag}.fem"
    femm_session.open_femm()
    try:
        f, tried = 1.0, []
        for _ in range(10):
            n = count_nodes(open_model(p, f, fine, path))
            tried.append((f, n))
            if target is None or abs(n / target - 1) < 0.01:    # within 1 %: the mesher's count is not smooth in f anyway
                break
            slope = -2.0                                         # nodes go roughly as f^-2; later steps use what was measured
            if len(tried) >= 2 and tried[-1][1] != tried[-2][1] and tried[-1][0] != tried[-2][0]:
                measured = math.log(tried[-1][1] / tried[-2][1]) / math.log(tried[-1][0] / tried[-2][0])
                slope = measured if measured < -0.3 else slope
            f *= (target / n) ** (1 / slope)
        f, n = min(tried, key=lambda fn: abs(fn[1] / target - 1)) if target else tried[-1]
        return tag, f, n                                         # the closest one tried: f and its own node count together
    finally:
        femm_session.close_femm()


def run_calibration() -> pd.DataFrame:
    (OUT / "fem").mkdir(parents=True, exist_ok=True)
    with ProcessPoolExecutor(min(WORKERS, len(CONFIGS))) as pool:
        got = list(pool.map(calibrate_one, CONFIGS))
    rows = [{"tag": t, "smart_mesh": fine, "target_nodes": target, "f": f, "nodes": n}
            for (t, fine, target), (_, f, n) in zip(CONFIGS, got)]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "calibration.csv", index=False)
    return df


# ----------------------------------------------------------------------------------------------- the slice ladder
def slice_jobs():
    """The points of the barrier study that matter here, each one a row of FEMM 4.2's phase1.csv or phase2.csv."""
    jobs = [{"set": "noload", "theta": float(t), "id": 0.0, "iq": 0.0, "gap": t == 0} for t in range(0, 360, 15)]
    jobs += [{"set": "cogging", "theta": float(t), "id": 0.0, "iq": 0.0, "gap": False} for t in np.arange(0, 20.01, 0.5)]
    for g in (0, 30, 90):
        i_d, i_q = current_vector(10.0, g)
        jobs += [{"set": f"gamma{g}", "theta": float(t), "id": i_d, "iq": i_q, "gap": False} for t in MEAN_POS]
    jobs += [{"set": "ripple", "theta": float(t), "id": IQ_ID_42[0], "iq": IQ_ID_42[1], "gap": False} for t in range(0, 61)]
    jobs += [{"set": f"check{k}", "theta": 0.0, "id": float(a), "iq": float(b), "gap": False}
             for k, (a, b) in enumerate(((0, 0), (0, 10), (-6, 8), (6, 8)))]
    return jobs


def slice_worker(task):
    """One worker, one chunk of jobs of one density: its own FEMM, its own .fem and .ans."""
    tag, fine, f, chunk, jobs, override = task
    p = mg.Params(**override)
    path = OUT / "fem" / f"{tag}_{chunk:02d}.fem"
    rows, gaps = [], {}
    femm_session.open_femm()
    try:
        model = open_model(p, f, fine, path)
        for job in jobs:
            t = time.perf_counter()
            r = model.solve(job["theta"], job["id"], job["iq"], gap_deg=GAP_DEG if job["gap"] else ())
            rows.append({"tag": tag, "set": job["set"], "theta": job["theta"], "id": job["id"], "iq": job["iq"],
                         **{k: r[k] for k in ("psi_a", "psi_d", "psi_q", "T_age", "T_wst")}, "seconds": time.perf_counter() - t})
            if job["gap"]:
                gaps[tag] = list(r["gap_br"])
    finally:
        femm_session.close_femm()
    return rows, gaps


def run_slice(cal: pd.DataFrame):
    jobs = slice_jobs()
    tasks = []
    for _, c in cal.iterrows():
        fine = bool(c.smart_mesh)
        for k in range(0, len(jobs), 12):
            tasks.append((c.tag, fine, float(c.f), k // 12, jobs[k:k + 12], {}))
    rows, gaps = [], {}
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        for part_rows, part_gaps in pool.map(slice_worker, tasks):
            rows += part_rows
            gaps.update(part_gaps)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "slice_ladder.csv", index=False)
    for tag, br in gaps.items():
        pd.DataFrame({"deg": GAP_DEG, "br": br}).to_csv(OUT / f"gap_{tag}.csv", index=False)
    print(f"slice: {len(df)} solves at {len(cal)} densities in {time.perf_counter() - t0:.0f} s on {WORKERS} workers")
    return df


# ----------------------------------------------------------------------------------------------- the measures
def ripple_measures(T) -> dict:
    """Mean torque, peak to peak and the harmonics 6, 12, ..., 36 [mN.m] of one electrical period of torque samples.

    T is the 60 torques of positions 0 ... 59 degrees (the 61st sample, at 60 degrees, repeats the first and is not passed).
    harmonic{6k} is the amplitude of the 6k-th harmonic of the electrical angle: the 60 samples hold 1 cycle of the 6th.
    """
    T = np.asarray(T, float)
    m = {"T_mean_mtpa_Nm": float(T.mean()), "ripple_pp_Nm": float(np.ptp(T))}
    h = harmonic_amplitudes(T)
    for k in range(1, 7):
        m[f"harmonic{6 * k}_mNm"] = float(h[k] * 1e3)
    return m


def measures(d: pd.DataFrame) -> dict:
    """The numbers the comparison is about, from the solved points of one density (or from FEMM 4.2's old rows).

    d needs the columns set, theta, psi_a, psi_d, psi_q, T_age. Sets: noload (24 positions, 15 degrees apart), cogging (41),
    gamma0 / gamma30 / gamma90 (six positions each at 10 A), ripple (61 positions at FEMM 4.2's best angle).
    """
    m = {}
    nl = d[d["set"] == "noload"].sort_values("theta")
    m["psi_pm_mWb"] = float(harmonic_amplitudes(nl.psi_a.values)[1] * 1e3)
    m["psi_d_noload_mean_mWb"] = float(nl.psi_d.mean() * 1e3)
    cg = d[d["set"] == "cogging"]
    m["cogging_pp_mNm"] = float(np.ptp(cg.T_age) * 1e3)
    for g in (0, 30, 90):
        s = d[d["set"] == f"gamma{g}"]
        m[f"T_mean_gamma{g}_Nm"] = float(s.T_age.mean())
        m[f"psi_d_gamma{g}_mWb"] = float(s.psi_d.mean() * 1e3)
        m[f"psi_q_gamma{g}_mWb"] = float(s.psi_q.mean() * 1e3)
    m.update(ripple_measures(d[d["set"] == "ripple"].sort_values("theta").T_age.values[:-1]))
    return m


def old_rows(old: Path) -> pd.DataFrame:
    """FEMM 4.2's rows for the same points, in the same layout (set, theta, psi_a, psi_d, psi_q, T_age)."""
    p1 = pd.read_csv(old / "barrier" / "phase1.csv")
    p2 = pd.read_csv(old / "barrier" / "phase2.csv")
    w = p1[p1.variant == "with band"]
    parts = [w[(w.study == "noload") & (w.theta_e % 15 == 0)].assign(set="noload"),
             w[w.study == "cogging"].assign(set="cogging")]
    for g in (0, 30, 90):
        parts.append(w[(w.study == "gamma") & (w.I == 10) & (w.gamma == g)].assign(set=f"gamma{g}"))
    parts.append(p2[(p2.study == "ripple") & (p2.variant == "with band")].assign(set="ripple"))
    d = pd.concat(parts, ignore_index=True).rename(columns={"theta_e": "theta"})
    return d[["set", "theta", "psi_a", "psi_d", "psi_q", "T_age"]]


def gap_deviation(old_gap: pd.DataFrame, new_gap: pd.DataFrame) -> tuple[float, float]:
    """(largest, rms) of |B_new - B_old| over the 360 samples, as a percentage of the old waveform's peak."""
    d = np.abs(new_gap.br.values - old_gap.br.values) / np.abs(old_gap.br).max() * 100
    return float(d.max()), float(np.sqrt((d ** 2).mean()))


def summarise_slice(df: pd.DataFrame, cal: pd.DataFrame, old: Path) -> pd.DataFrame:
    rows = []
    ref = measures(old_rows(old))
    old_gap = pd.read_csv(old / "barrier" / "gap_noload_with_band.csv")
    rows.append({"tag": "FEMM 4.2 (fine, its own mesh)", "nodes": README_NODES["fine"], **ref, "gap_max_pct": 0.0, "gap_rms_pct": 0.0})
    for _, c in cal.iterrows():
        d = df[df.tag == c.tag]
        gap_path = OUT / f"gap_{c.tag}.csv"
        gmax, grms = gap_deviation(old_gap, pd.read_csv(gap_path)) if gap_path.is_file() else (math.nan, math.nan)
        rows.append({"tag": c.tag, "nodes": int(c.nodes), **measures(d), "gap_max_pct": gmax, "gap_rms_pct": grms})
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "slice_summary.csv", index=False)
    return s


def print_slice(s: pd.DataFrame):
    show = ["psi_pm_mWb", "psi_d_gamma90_mWb", "T_mean_gamma0_Nm", "T_mean_mtpa_Nm", "cogging_pp_mNm", "ripple_pp_Nm", "harmonic24_mNm",
            "harmonic30_mNm", "gap_max_pct"]
    print("\nthe barrier study's numbers (with band, 10 A) at several densities; the first row is FEMM 4.2")
    print(f"{'density':32s} {'nodes':>7s}  " + " ".join(f"{c[:15]:>15s}" for c in show))
    for _, r in s.iterrows():
        print(f"{r.tag:32s} {int(r.nodes):7,d}  " + " ".join(f"{r[c]:15.5g}" for c in show))
    ref = s.iloc[0]
    print("\ndifference from FEMM 4.2, in percent of FEMM 4.2's value (gap: percent of the waveform's peak)")
    print(f"{'density':32s} {'nodes':>7s}  " + " ".join(f"{c[:15]:>15s}" for c in show))
    for _, r in s.iloc[1:].iterrows():
        print(f"{r.tag:32s} {int(r.nodes):7,d}  " + " ".join(
            f"{r[c]:15.3f}" if c == "gap_max_pct" else f"{100 * (r[c] / ref[c] - 1):+15.3f}" for c in show))


# ----------------------------------------------------------------------------------------------- check_femm at every density
def summarise_check_femm(df: pd.DataFrame, cal: pd.DataFrame):
    """The four check_femm points (theta = 0) at every density, against FEMM 4.2's coarse numbers (README_CHECK_FEMM in compare.py)."""
    rows = []
    for _, c in cal.iterrows():
        d = df[df.tag == c.tag]
        row = {"tag": c.tag, "nodes": int(c.nodes)}
        for k in range(4):
            r = d[d["set"] == f"check{k}"].iloc[0]
            row[f"psi_d_{k}_mWb"], row[f"psi_q_{k}_mWb"] = r.psi_d * 1e3, r.psi_q * 1e3
            row[f"T_age_{k}"], row[f"T_wst_{k}"] = r.T_age, r.T_wst
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(OUT / "check_femm_ladder.csv", index=False)
    r42 = README_CHECK_FEMM
    print(f"\ncheck_femm's four points at several densities (FEMM 4.2, coarse mesh, {README_NODES['coarse']:,} nodes: psi_d "
          f"{r42[(0, 0)]['psi_d']}, T {r42[(0, 10)]['T_age']} / {r42[(-6, 8)]['T_age']} / {r42[(6, 8)]['T_age']})")
    print(f"{'density':14s} {'nodes':>7s} {'psi_d(0,0)':>11s} {'T(0,10)':>9s} {'T(-6,8)':>9s} {'T(+6,8)':>9s}")
    for _, r in out.iterrows():
        print(f"{r.tag:14s} {int(r.nodes):7,d} {r.psi_d_0_mWb:11.3f} {r.T_age_1:9.4f} {r.T_age_2:9.4f} {r.T_age_3:9.4f}")
    return out


# ----------------------------------------------------------------------------------------------- the flux explorer's states
def states_worker(task):
    """One worker: one rotor of the flux explorer at one coarse density: all 60 states of its air-band model."""
    tag, f, gid, x = task
    p = mg.Params(barrier="air", **x)
    path = OUT / "fem" / f"states_{tag}_g{gid:02d}.fem"
    jobs = [(i, g, th) for i, g in sweep.field_states() for th in sweep.POS]
    jobs += [(i, g, th) for i in sweep.FIELD_I for g in sweep.SCALAR_G for th in sweep.POS]
    rows = []
    femm_session.open_femm()
    try:
        model = open_model(p, f, False, path)
        nodes = count_nodes(model)
        for i, g, th in jobs:
            r = model.solve(th, *current_vector(i, g))
            rows.append({"tag": tag, "nodes": nodes, "gid": gid, "I": i, "gamma": g, "theta_e": th, "T_age": r["T_age"],
                         "psi_d": r["psi_d"], "psi_q": r["psi_q"]})
    finally:
        femm_session.close_femm()
    return rows


def run_states(cal: pd.DataFrame, old: Path):
    geoms = dict(sweep.geometries())
    configs = [c for c in cal.itertuples() if c.tag in ("c-natural", "c-matched", "c-17k", "c-28k")]
    tasks = [(c.tag, float(c.f), gid, geoms[gid]) for c in configs for gid in STATES_ROTORS]
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part in pool.map(states_worker, tasks) for r in part]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "states_ladder.csv", index=False)
    print(f"states: {len(df)} solves in {time.perf_counter() - t0:.0f} s")
    return summarise_states(df, old)


def summarise_states(df: pd.DataFrame, old: Path) -> pd.DataFrame:
    """Each density's 60 states per rotor against FEMM 4.2's states of the same rotor (the old states.csv)."""
    o = pd.read_csv(old / "flux_explorer" / "states.csv")
    o = o[o.band == "air"]
    key = ["gid", "I", "gamma", "theta_e"]
    t_scale = float(o.T_age.abs().max())
    f_scale = float(max(o.psi_d.abs().max(), o.psi_q.abs().max()))
    rows = []
    for (tag, gid), d in df.groupby(["tag", "gid"]):
        m = d.merge(o, on=key, suffixes=("_n", "_o"))
        dT = 100 * (m.T_age_n - m.T_age_o) / np.maximum(m.T_age_o.abs(), 0.1 * t_scale)
        dd = 100 * (m.psi_d_n - m.psi_d_o) / np.maximum(m.psi_d_o.abs(), 0.1 * f_scale)
        dq = 100 * (m.psi_q_n - m.psi_q_o) / np.maximum(m.psi_q_o.abs(), 0.1 * f_scale)
        rows.append({"tag": tag, "gid": gid, "nodes": int(d.nodes.iloc[0]), "n_states": len(m),
                     "T_median": dT.abs().median(), "T_max": dT.abs().max(), "psi_d_median": dd.abs().median(), "psi_d_max": dd.abs().max(),
                     "psi_q_median": dq.abs().median(), "psi_q_max": dq.abs().max(), "psi_d_mean_offset_mWb": (m.psi_d_n - m.psi_d_o).mean() * 1e3})
    s = pd.DataFrame(rows)
    s.to_csv(OUT / "states_summary.csv", index=False)
    print("\nthe flux explorer's states (air band, coarse mesh) against FEMM 4.2's, by density (|difference| in percent, the "
          "judged percentage of compare.py)")
    print(f"{'density':12s} {'rotor':>5s} {'nodes':>7s} {'T med':>7s} {'T max':>7s} {'psi_d med':>10s} {'psi_d max':>10s} {'psi_q max':>10s} "
          f"{'psi_d offset':>13s}")
    for _, r in s.sort_values(["gid", "nodes"]).iterrows():
        print(f"{r.tag:12s} {int(r.gid):5d} {int(r.nodes):7,d} {r.T_median:7.3f} {r.T_max:7.3f} {r.psi_d_median:10.3f} {r.psi_d_max:10.3f} "
              f"{r.psi_q_max:10.3f} {r.psi_d_mean_offset_mWb:+11.3f}mWb")
    return s


# ----------------------------------------------------------------------------------------------- seconds per solve
def timing_worker(task):
    """One FEMM, one job at a time: the median seconds a solve() takes at this density (nothing else of ours is running)."""
    tag, fine, f, nodes = task
    p = mg.Params()
    path = OUT / "fem" / f"timing_{tag}.fem"
    femm_session.open_femm()
    try:
        model = open_model(p, f, fine, path)
        model.solve(0.0, 0.0, 10.0)                                     # the first solve of a document is slower: not counted
        times = []
        for th in range(0, 30, 2):
            t = time.perf_counter()
            model.solve(float(th), 0.0, 10.0)
            times.append(time.perf_counter() - t)
        return {"tag": tag, "nodes": nodes, "median_s": statistics.median(times), "min_s": min(times), "max_s": max(times), "solves": len(times)}
    finally:
        femm_session.close_femm()


def run_timing(cal: pd.DataFrame):
    rows = []
    for tag in ("c-natural", "c-matched", "f-natural", "f-matched"):
        c = cal[cal.tag == tag].iloc[0]
        rows.append(timing_worker((tag, bool(c.smart_mesh), float(c.f), int(c.nodes))))      # strictly one after the other
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "timing.csv", index=False)
    print("\nseconds per solve, one FEMM alone on the machine (15 solves, the first one not counted)")
    for _, r in df.iterrows():
        print(f"  {r.tag:12s} {int(r.nodes):7,d} nodes: median {r.median_s:.3f} s (range {r.min_s:.3f} to {r.max_s:.3f})")
    return df


# ----------------------------------------------------------------------------------------------- one bit of a design
def exact_csv_rows(path: Path) -> list[dict]:
    """The rows of a CSV, every number still text. Convert with Python's float(), which is exact.

    pandas' default CSV reader is not: it can differ from the number that was written in the last bit (here in 27 of the
    192 parameters of the optimizer's 48 designs). For this model one unit in the last place of a parameter is enough to
    change FEMM's mesh, so a design read back through pandas' default reader is not quite the design that was solved. Use
    pd.read_csv(..., float_precision="round_trip") when a design is to be solved again.
    """
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def evaluate_worker(design):
    """optimize_demo.evaluate, as the optimizer's workers run it, but with its scratch files in this folder."""
    import optimize_demo as od
    od.OUT = str(OUT / "opt")
    os.makedirs(os.path.join(od.OUT, "fem"), exist_ok=True)
    return od.evaluate(design)


def run_candidate(old: Path):
    """The optimizer's generation-1 candidate whose best current angle moved most between FEMM 4.2 and FEMM 5.0, solved again.

    The best angle of a candidate is the vertex of a parabola through three mean torques (at 25, 45 and 65 degrees), so it is as
    sensitive as three torques are. The question is how sensitive that is: the candidate is solved by the optimizer's own function
    as it is (it must give the sweep's numbers to the last digit), and again with ONE parameter moved by ONE unit in the last
    place (about 1e-15 mm: no physical change at all). Whatever those do to the angle is the noise of the number.
    """
    import optimize_demo as od
    ro = [r for r in exact_csv_rows(old / "opt" / "evaluations.csv") if r["gen"] == "1"]
    rn = [r for r in exact_csv_rows(P / "out" / "opt" / "evaluations.csv") if r["gen"] == "1"]
    both = [i for i in range(len(rn)) if ro[i]["feasible"] == "True" and rn[i]["feasible"] == "True"]
    k = max(both, key=lambda i: abs(float(rn[i]["gamma"]) - float(ro[i]["gamma"])) / float(ro[i]["gamma"]))
    design = {n: float(rn[k][n]) for n in od.NAMES}
    variants = {"as solved by the sweep": design}
    for name in od.NAMES:
        for sign, toward in (("+", math.inf), ("-", -math.inf)):
            v = dict(design)
            v[name] = math.nextafter(v[name], toward)
            variants[f"{name} {sign}1 ulp"] = v
    with ProcessPoolExecutor(min(WORKERS, len(variants))) as pool:
        results = list(pool.map(evaluate_worker, variants.values()))
    rows = [{"variant": name, "gamma": r["gamma"], "T_mean": r["T_mean"], "ripple": r["ripple"]} for name, r in zip(variants, results)]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "candidate.csv", index=False)
    same = (results[0]["gamma"], results[0]["T_mean"], results[0]["ripple"]) == (float(rn[k]["gamma"]), float(rn[k]["T_mean"]), float(rn[k]["ripple"]))
    print(f"\ncandidate {k} of generation 1: best current angle {float(ro[k]['gamma']):.2f} deg with FEMM 4.2 (torque {float(ro[k]['T_mean']):.4f} N.m), "
          f"{float(rn[k]['gamma']):.2f} deg in the FEMM 5.0 sweep ({float(rn[k]['T_mean']):.4f} N.m)")
    print(f"  solved again as it is: {'identical to the sweep to the last digit' if same else 'NOT identical to the sweep'}")
    print(f"  with one parameter moved by one unit in the last place (8 variants): best angle {df.gamma[1:].min():.2f} to {df.gamma[1:].max():.2f} deg, "
          f"mean torque {df.T_mean[1:].min():.5f} to {df.T_mean[1:].max():.5f} N.m, ripple {df.ripple[1:].min():.5f} to {df.ripple[1:].max():.5f} N.m")
    return df


def run_ulp():
    """All the optimizer's feasible designs, solved again as they were and with mag_y one unit in the last place away.

    Two results. First, a design read back EXACTLY gives the sweep's row to the last digit, every time (FEMM 5.0 is repeatable; what looked
    like a difference in an earlier attempt was pandas' reader, see exact_csv_rows). Second, one bit of one parameter changes the answer:
    the noise of the optimizer's numbers (best angle, mean torque, ripple), which is what a difference between two programs is to be
    judged against.
    """
    rows = [r for r in exact_csv_rows(P / "out" / "opt" / "evaluations.csv") if r["feasible"] == "True"]
    names = ["mag_y", "mag_x", "v_half_deg", "web"]
    exact = [{n: float(r[n]) for n in names} for r in rows]
    moved = [dict(d, mag_y=math.nextafter(d["mag_y"], math.inf)) for d in exact]
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        r_exact = list(pool.map(evaluate_worker, exact))
        r_moved = list(pool.map(evaluate_worker, moved))
    same = [(a["gamma"], a["T_mean"], a["ripple"]) == (float(r["gamma"]), float(r["T_mean"]), float(r["ripple"])) for a, r in zip(r_exact, rows)]
    dg = np.array([abs(m["gamma"] - a["gamma"]) for m, a in zip(r_moved, r_exact)])
    dt = np.array([abs(m["T_mean"] / a["T_mean"] - 1) * 100 for m, a in zip(r_moved, r_exact)])
    dr = np.array([abs(m["ripple"] / a["ripple"] - 1) * 100 for m, a in zip(r_moved, r_exact)])
    df = pd.DataFrame({"gamma_sweep": [float(r["gamma"]) for r in rows], "gamma_exact": [a["gamma"] for a in r_exact], "gamma_1ulp": [m["gamma"] for m in r_moved],
                       "T_sweep": [float(r["T_mean"]) for r in rows], "T_exact": [a["T_mean"] for a in r_exact], "T_1ulp": [m["T_mean"] for m in r_moved],
                       "ripple_sweep": [float(r["ripple"]) for r in rows], "ripple_exact": [a["ripple"] for a in r_exact], "ripple_1ulp": [m["ripple"] for m in r_moved]})
    df.to_csv(OUT / "ulp.csv", index=False)
    print(f"ulp: {2 * len(rows)} evaluations in {time.perf_counter() - t0:.0f} s on {WORKERS} workers")
    print(f"  {len(rows)} feasible designs solved again exactly as they were: identical to the sweep's row (best angle, mean torque, ripple) in {sum(same)} of {len(rows)}")
    print(f"  mag_y moved by one unit in the last place: the result changed in {int(((dg > 0) | (dt > 0)).sum())} of {len(rows)} designs")
    print(f"    best angle: median change {np.median(dg):.3f} deg, 95th percentile {np.percentile(dg, 95):.3f}, largest {dg.max():.3f}")
    print(f"    mean torque: median {np.median(dt):.3f} %, largest {dt.max():.2f} %;  ripple: median {np.median(dr):.2f} %, largest {dr.max():.1f} %")
    return df


# ----------------------------------------------------------------------------------------------- the noise of the mesh itself
NOISE_VARIANTS = {"baseline": {},
                  "mag_y + 1 ulp": {"mag_y": math.nextafter(mg.Params().mag_y, math.inf)},
                  "mag_x - 1 ulp": {"mag_x": math.nextafter(mg.Params().mag_x, -math.inf)}}
NOISE_DENSITIES = ("c-natural", "f-natural", "f-matched")


def run_noise(cal: pd.DataFrame, old: Path):
    """How much does a number move when the mesh changes and the physics does not?

    One parameter of the design is moved by one unit in the last place (about 1e-15 mm), which no model of a magnet could
    call a change, and FEMM meshes the slightly different drawing differently. Every quantity of the slice ladder is worked out
    again. What it moves by is that quantity's mesh noise: a difference between FEMM 4.2 and FEMM 5.0 that is not bigger than
    this is a difference between two meshes of the same thing, not between two programs.
    """
    base = pd.read_csv(OUT / "slice_ladder.csv")
    jobs = slice_jobs()
    tasks, owner = [], {}
    for tag in NOISE_DENSITIES:
        c = cal[cal.tag == tag].iloc[0]
        for v, (name, override) in enumerate(NOISE_VARIANTS.items()):
            if not override:
                continue
            for k in range(0, len(jobs), 12):
                job_tag = f"noise-{tag}-v{v}"
                tasks.append((job_tag, bool(c.smart_mesh), float(c.f), k // 12, jobs[k:k + 12], override))
                owner[job_tag] = (tag, name)
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [r for part_rows, _ in pool.map(slice_worker, tasks) for r in part_rows]
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "noise_ladder.csv", index=False)
    print(f"noise: {len(df)} solves in {time.perf_counter() - t0:.0f} s on {WORKERS} workers")
    summary = []
    show = ["psi_pm_mWb", "psi_d_gamma90_mWb", "T_mean_gamma0_Nm", "T_mean_mtpa_Nm", "cogging_pp_mNm", "ripple_pp_Nm", "harmonic18_mNm",
            "harmonic24_mNm", "harmonic30_mNm", "harmonic36_mNm"]
    for tag in NOISE_DENSITIES:
        values = {"baseline": measures(base[base.tag == tag])}
        for job_tag, (t, name) in owner.items():
            if t == tag:
                values[name] = measures(df[df.tag == job_tag])
        for q in show:
            vals = [values[n][q] for n in NOISE_VARIANTS]
            summary.append({"density": tag, "quantity": q, **{n: values[n][q] for n in NOISE_VARIANTS},
                            "spread_pct": 100 * (max(vals) - min(vals)) / abs(values["baseline"][q])})
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "noise_summary.csv", index=False)
    ref = pd.read_csv(OUT / "slice_summary.csv").iloc[0]
    print("\nmesh noise: the same quantities when one parameter moves by one unit in the last place (a spread is max - min over the three, in "
          "percent of the baseline)")
    print(f"{'quantity':20s} " + " ".join(f"{t:>12s}" for t in NOISE_DENSITIES) + "   (FEMM 4.2 against FEMM 5.0 at its natural fine mesh)")
    fine = s[s.density == "f-natural"].set_index("quantity")
    for q in show:
        diff = 100 * (fine.loc[q, "baseline"] / ref[q] - 1)
        print(f"{q:20s} " + " ".join(f"{s[(s.density == t) & (s.quantity == q)].spread_pct.iloc[0]:11.2f}%" for t in NOISE_DENSITIES) + f"   {diff:+8.2f} %")
    return s


# ----------------------------------------------------------------------------------------------- the flags that were argued by analogy
FLAGGED_LADDER = ("f-natural", "f-matched", "f-45k", "f-70k")  # the smart-mesh densities; f-natural is the barrier study's own mesh
FLAGGED_NOISE = ("f-natural", "f-matched")                      # where the no-band rotor also gets its one-bit moves (as run_noise does)
STATES_20A = ((20.0, 90.0, 0.0), (20.0, 90.0, 40.0))            # (I, gamma, theta_e) of the two torques that compare.py flags (F2)


def ripple_current_no_band(old: Path) -> tuple[float, float]:
    """(i_d, i_q) FEMM 4.2 solved the no-band rotor's ripple sweep at (its own best current angle, 22.87 degrees), from its phase2.csv."""
    p2 = pd.read_csv(Path(old) / "barrier" / "phase2.csv", float_precision="round_trip")
    r = p2[(p2.study == "ripple") & (p2.variant == "no band")].sort_values("theta_e").iloc[0]
    return float(r["id"]), float(r["iq"])


def sweep_ripple_no_band(base: Path) -> dict:
    """ripple_measures of the no-band rotor's ripple sweep in a barrier study's phase2.csv (FEMM 4.2's, or FEMM 5.0's own)."""
    p2 = pd.read_csv(Path(base) / "barrier" / "phase2.csv", float_precision="round_trip")
    r = p2[(p2.study == "ripple") & (p2.variant == "no band")].sort_values("theta_e")
    return ripple_measures(r.T_age.values[:-1])


def sweep_state_torque(base: Path, current: float, gamma: float, theta: float) -> float:
    """The with-band torque T_age [N.m] of one (I, gamma, theta_e) state in a barrier study's phase1.csv."""
    p1 = pd.read_csv(Path(base) / "barrier" / "phase1.csv", float_precision="round_trip")
    r = p1[(p1.study == "gamma") & (p1.variant == "with band") & (p1.I == current) & (p1.gamma == gamma) & (p1.theta_e == theta)]
    return float(r.T_age.iloc[0])


def run_flagged(cal: pd.DataFrame, old: Path):
    """Solve the two flags that the report (2.12) first argued by analogy.

    F2: the two torques of the torque-vs-angle grid that are beyond 1 % are at 20 A all on the d axis (gamma 90 degrees), where the
    torque is a small remainder. They are solved on the smart-mesh ladder, as they are and with one parameter one bit away (the
    noise of the mesh). The first column must equal FEMM 5.0's own sweep: the same model on the same mesh gives the same number.

    F3: the no-band rotor's 24th and 30th torque harmonics are +37 % and +61 % against FEMM 4.2. Its ripple sweep (61 positions at
    FEMM 4.2's own best angle for it) is solved on the same ladder and with the same one-bit moves. "No band" is the barrier
    study's barrier="steel": the band is drawn and filled with steel, so both variants get the identical mesh (and the node counts
    of the with-band calibration hold for it).
    """
    i_d_nb, i_q_nb = ripple_current_no_band(old)
    states = [{"set": f"state{k}", "theta": th, "id": current_vector(i, g)[0], "iq": current_vector(i, g)[1], "gap": False}
              for k, (i, g, th) in enumerate(STATES_20A)]
    ripple = [{"set": "ripple", "theta": float(t), "id": i_d_nb, "iq": i_q_nb, "gap": False} for t in range(0, 61)]
    tasks, owner = [], {}
    for kind, jobs, base in (("states", states, {}), ("noband", ripple, {"barrier": "steel"})):
        for tag in FLAGGED_LADDER:
            c = cal[cal.tag == tag].iloc[0]
            for v, (name, override) in enumerate(NOISE_VARIANTS.items()):
                if v and kind == "noband" and tag not in FLAGGED_NOISE:
                    continue
                job_tag = f"flag-{kind}-{tag}-v{v}"
                for k in range(0, len(jobs), 12):
                    tasks.append((job_tag, bool(c.smart_mesh), float(c.f), k // 12, jobs[k:k + 12], {**base, **override}))
                owner[job_tag] = (kind, tag, name, int(c.nodes))
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        rows = [{**r, "kind": owner[r["tag"]][0], "density": owner[r["tag"]][1], "variant": owner[r["tag"]][2], "nodes": owner[r["tag"]][3]}
                for part_rows, _ in pool.map(slice_worker, tasks) for r in part_rows]
    df = pd.DataFrame(rows)
    print(f"flagged: {len(df)} solves in {time.perf_counter() - t0:.0f} s on {WORKERS} workers")

    # F2: the two 20 A states
    st = df[df.kind == "states"]
    st.to_csv(OUT / "flagged_states.csv", index=False)
    new_out = P / "out"
    old_t = [sweep_state_torque(old, *s) * 1e3 for s in STATES_20A]
    new_t = [sweep_state_torque(new_out, *s) * 1e3 for s in STATES_20A] if (new_out / "barrier" / "phase1.csv").is_file() else [math.nan] * 2
    print("\nF2: with band, 20 A on the d axis (gamma = 90 degrees): the two torques that compare.py flags, T_age in mN.m")
    print(f"  rotor position {STATES_20A[0][2]:g} degrees / {STATES_20A[1][2]:g} degrees:  FEMM 4.2 {old_t[0]:.3f} / {old_t[1]:.3f};  "
          f"FEMM 5.0's sweep {new_t[0]:.3f} / {new_t[1]:.3f}")
    print(f"  {'density':10s} {'nodes':>7s}   " + "   ".join(f"{'as solved':>9s} {'mag_y+1':>8s} {'mag_x-1':>8s} {'spread':>7s}" for _ in STATES_20A)
          + "   (one bit = one unit in the last place)")
    for tag in FLAGGED_LADDER:
        d = st[st.density == tag]
        cells = []
        for k in range(len(STATES_20A)):
            v = [float(d[(d.variant == name) & (d["set"] == f"state{k}")].T_age.iloc[0]) * 1e3 for name in NOISE_VARIANTS]
            cells.append(f"{v[0]:9.3f} {v[1]:8.3f} {v[2]:8.3f} {max(v) - min(v):7.3f}")
        print(f"  {tag:10s} {int(d.nodes.iloc[0]):7,d}   " + "   ".join(cells))

    # F3: the no-band rotor's ripple
    nb = df[df.kind == "noband"]
    summary = []
    for (density, variant), d in nb.groupby(["density", "variant"], sort=False):
        summary.append({"density": density, "variant": variant, "nodes": int(d.nodes.iloc[0]),
                        **ripple_measures(d.sort_values("theta").T_age.values[:-1])})
    s = pd.DataFrame(summary)
    s.to_csv(OUT / "flagged_noband.csv", index=False)
    r42, r50 = sweep_ripple_no_band(old), (sweep_ripple_no_band(new_out) if (new_out / "barrier" / "phase2.csv").is_file() else None)
    print(f"\nF3: the no-band rotor, ripple sweep at FEMM 4.2's best angle for it ({math.degrees(math.atan2(-i_d_nb, i_q_nb)):.2f} degrees): "
          "harmonics in mN.m, torque in N.m")
    print(f"  {'':34s} {'nodes':>7s} {'24th':>7s} {'30th':>7s} {'ripple pp':>10s} {'mean T':>8s}")
    print(f"  {'FEMM 4.2':34s} {'':>7s} {r42['harmonic24_mNm']:7.3f} {r42['harmonic30_mNm']:7.3f} {r42['ripple_pp_Nm']:10.4f} {r42['T_mean_mtpa_Nm']:8.4f}")
    if r50:
        print(f"  {'FEMM 5.0 sweep':34s} {'':>7s} {r50['harmonic24_mNm']:7.3f} {r50['harmonic30_mNm']:7.3f} {r50['ripple_pp_Nm']:10.4f} {r50['T_mean_mtpa_Nm']:8.4f}")
    for _, r in s.iterrows():
        print(f"  {r.density + ', ' + r.variant:34s} {int(r.nodes):7,d} {r.harmonic24_mNm:7.3f} {r.harmonic30_mNm:7.3f} {r.ripple_pp_Nm:10.4f} {r.T_mean_mtpa_Nm:8.4f}")
    ladder = s[s.variant == "baseline"]
    print(f"  ladder (baseline, {len(ladder)} densities): 24th {ladder.harmonic24_mNm.min():.3f} to {ladder.harmonic24_mNm.max():.3f}, "
          f"30th {ladder.harmonic30_mNm.min():.3f} to {ladder.harmonic30_mNm.max():.3f}")
    for tag in FLAGGED_NOISE:
        d = s[s.density == tag]
        print(f"  one-bit moves at {tag}: 24th {d.harmonic24_mNm.min():.3f} to {d.harmonic24_mNm.max():.3f} "
              f"({100 * (d.harmonic24_mNm.max() - d.harmonic24_mNm.min()) / d.harmonic24_mNm.iloc[0]:.0f} % of the baseline), "
              f"30th {d.harmonic30_mNm.min():.3f} to {d.harmonic30_mNm.max():.3f} "
              f"({100 * (d.harmonic30_mNm.max() - d.harmonic30_mNm.min()) / d.harmonic30_mNm.iloc[0]:.0f} %)")
    return df, s


# ----------------------------------------------------------------------------------------------- FEMM 5.0 reading FEMM 4.2's solutions
def run_post42(old: Path) -> pd.DataFrame:
    """Does FEMM 5.0 read a FEMM 4.2 solution the way FEMM 4.2 did? Open the six solved sectors that the barrier study kept
    (FEMM 4.2's own .ans files: its mesh, its solved field), let FEMM 5.0's post-processor work out the flux linkages, the two
    torques and the air-gap B from them, and compare with what FEMM 4.2 recorded for the same solution.

    A copy is opened, never the original. If the numbers agree, everything that differs between the two programs' results
    was decided before the post-processor: by the mesh and the solver. If one of them does not, that quantity is worked out
    differently by the two programs, whatever the mesh.
    """
    import femm
    src, dst = Path(old) / "barrier" / "fields", OUT / "ans42"
    dst.mkdir(parents=True, exist_ok=True)
    p1 = pd.read_csv(Path(old) / "barrier" / "phase1.csv")
    p2 = pd.read_csv(Path(old) / "barrier" / "phase2.csv")

    def recorded(case, variant):
        v = "with band" if variant == "with_band" else "no band"
        if case == "noload":
            return p1[(p1.study == "noload") & (p1.variant == v) & (p1.theta_e == 0)].iloc[0]
        if case == "qaxis":
            return p2[(p2.study == "snapshot") & (p2.variant == v)].iloc[0]
        return p2[(p2.study == "ripple") & (p2.variant == v) & (p2.theta_e == 0)].iloc[0]

    rows = []
    femm_session.open_femm()
    try:
        for case in ("noload", "qaxis", "mtpa"):
            for variant in ("with_band", "no_band"):
                name = f"{case}_{variant}.ans"
                shutil.copyfile(src / name, dst / name)
                femm.opendocument(str(dst / name))
                n = int(femm.mo_numnodes())
                psi = [4 * femm.mo_getcircuitproperties(ph)[2] for ph in "ABC"]        # 4 sectors make the machine
                t_age = femm.mo_gapintegral("AGE", 0)
                femm.mo_groupselectblock(1)
                t_wst = 4 * femm.mo_blockintegral(22)
                br = np.array([femm.mo_getgapb("AGE", float(a))[0] for a in GAP_DEG])
                femm.mo_close()
                r = recorded(case, variant)
                gap = pd.read_csv(Path(old) / "barrier" / f"gap_{case}_{variant}.csv")
                rows.append({"solution": name[:-4], "nodes": n,
                             "psi_max_rel_diff": max(abs(psi[0] / r.psi_a - 1), abs(psi[1] / r.psi_b - 1), abs(psi[2] / r.psi_c - 1)),
                             "T_wst_old": r.T_wst, "T_wst_diff": t_wst - r.T_wst, "T_age_old": r.T_age, "T_age_diff": t_age - r.T_age,
                             "gap_br_max_diff_pct_of_peak": float(np.abs(br - gap.br.values).max() / np.abs(gap.br).max() * 100)})
    finally:
        femm_session.close_femm()
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "post42.csv", index=False)
    print("\nFEMM 5.0's post-processor on FEMM 4.2's own solutions, against what FEMM 4.2 recorded for them")
    print(f"{'solution':18s} {'nodes':>7s} {'psi: largest rel. diff':>23s} {'T_wst diff [N.m]':>17s} {'T_age diff [N.m]':>17s} {'T_age diff, % of T':>19s} "
          f"{'gap B_r diff, % of peak':>24s}")
    for _, r in df.iterrows():
        pct = f"{100 * r.T_age_diff / abs(r.T_age_old):.2f}" if abs(r.T_age_old) > 0.1 else "(no load)"
        print(f"{r.solution:18s} {int(r.nodes):7,d} {r.psi_max_rel_diff:23.1e} {r.T_wst_diff:17.1e} {r.T_age_diff:17.2e} {pct:>19s} "
              f"{r.gap_br_max_diff_pct_of_peak:24.3f}")
    return df


# ----------------------------------------------------------------------------------------------- the hand-drawn 360-degree model
FULL_OFFSETS = (0.0, 0.02, 0.04)                               # extra rotor angle in electrical degrees: a new mesh, the same physics
FULL_DENSE_STEP = 3.75                                         # electrical degrees between the positions of the dense grid: 96 per cycle
FULL_LABEL, SLICE_LABEL = "hand-drawn 360", "motorgeo 90"      # the two values of the model column of validate.csv


def full_model_jobs() -> list[tuple[int, float, float]]:
    """validate_full.py's 48 solves of the hand-drawn model, in its order: (k, rotor angle [electrical degrees], phase-A current [A]).
    k numbers the files hand_000 ... hand_047."""
    import validate_full as vf
    return [(k, float(th), i) for k, (th, i) in enumerate((th, i) for i in (0.0, vf.I_LOAD) for th in vf.THETAS)]


def run_post42full(old: Path) -> pd.DataFrame:
    """Does the full model's torque offset come from the post-processor? Open FEMM 4.2's own 48 saved solutions of the hand-drawn 360-degree
    model (its mesh, its solved field) in FEMM 5.0, work out phase A's flux linkage and the stress-tensor torque again, exactly as
    validate_full.py does, and compare with what FEMM 4.2 recorded in validate.csv for the same solutions.

    The full model has no air-gap element (its only boundaries are the two outer ones), so it has no T_age: only the stress tensor
    T_wst, which is also what the slice check showed to be read identically by both programs. A copy is opened, never the original.
    If the numbers agree, the offset was decided before the post-processor: by the mesh and the solver.
    """
    import femm
    src, dst = Path(old) / "validate", OUT / "ans42_full"
    dst.mkdir(parents=True, exist_ok=True)
    ref = pd.read_csv(src / "validate.csv", float_precision="round_trip")
    ref = ref[ref.model == FULL_LABEL].set_index(["I", "theta_e"])
    rows = []
    femm_session.open_femm()
    try:
        for k, theta, current in full_model_jobs():
            name = f"hand_{k:03d}.ans"
            shutil.copyfile(src / name, dst / name)
            femm.opendocument(str(dst / name))
            circuit = femm.mo_getcircuitproperties("Phase A")          # (current, voltage, flux linkage)
            femm.mo_groupselectblock(1)                                # the rotor
            t_wst = femm.mo_blockintegral(22)
            nodes = int(femm.mo_numnodes())
            femm.mo_close()
            r = ref.loc[(current, theta)]
            rows.append({"k": k, "I": current, "theta_e": theta, "nodes": nodes, "I_read": circuit[0], "psi_a_old": r.psi_a,
                         "psi_a_diff": circuit[2] - r.psi_a, "T_old": r["T"], "T_new": t_wst, "T_diff": t_wst - r["T"]})
    finally:
        femm_session.close_femm()
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "post42_full.csv", index=False)
    print("\nFEMM 5.0's post-processor on FEMM 4.2's 48 saved solutions of the hand-drawn 360-degree model, against what FEMM 4.2 recorded for them")
    for current, d in df.groupby("I"):
        print(f"  {current:4.1f} A: {int(d.nodes.min()):,} to {int(d.nodes.max()):,} nodes; largest |flux linkage difference| "
              f"{d.psi_a_diff.abs().max():.1e} Wb ({d.psi_a_diff.abs().max() / d.psi_a_old.abs().max():.1e} of the peak); "
              f"largest |T_wst difference| {d.T_diff.abs().max():.1e} N.m; the 24-position mean torque: FEMM 4.2 recorded "
              f"{d.T_old.mean() * 1e3:+.3f}, FEMM 5.0 reading the same solutions {d.T_new.mean() * 1e3:+.3f} mN.m")
    return df


FULL_FEM = paths.path("reference_fem")                          # the hand-drawn model: config.toml's reference, the file validate_full.py reads


def full_worker(task):
    """One solve of the hand-drawn 360-degree model, as validate_full.py does it, with the rotor turned a little further.

    validate_full.hand_worker copies the model to a file of this job's own, turns the rotor, sets the currents, solves and reads
    phase A's flux linkage and the stress-tensor torque. Only one thing is changed here: the job's files go to out/validation/
    mesh/full (not into the study's out/validate). The model is the one validate_full.py reads, config.toml's reference_fem, the
    copy that is kept beside the code in P/reference (SHA-256 3AC2F628...; the reference commit explains why); it is named again
    below so that this check does not depend on validate_full's own default.
    """
    k, theta_e, i_a = task
    import validate_full as vf
    vf.OUT = str(OUT / "full")                                  # the files of this check, not the study's
    if FULL_FEM.is_file():
        vf.FEM = str(FULL_FEM)
    return vf.hand_worker((k, theta_e, i_a))


def offset_analysis(old_t, new_t, meshes) -> dict:
    """How much of the difference between two programs' torque waveforms, at the same rotor positions, is ONE common offset.

    old_t, new_t: the torque [N.m] of FEMM 4.2 and of FEMM 5.0 at the positions. meshes: FEMM 5.0's torque [N.m] on several meshes of
    the same physics (the rotor turned by tiny extra angles), positions as rows, meshes as columns. Everything is returned in mN.m.
    The noise of a position is the largest minus the smallest value over the meshes.
    """
    o, n, a = np.asarray(old_t, float) * 1e3, np.asarray(new_t, float) * 1e3, np.asarray(meshes, float) * 1e3
    d, d0 = n - o, (n - n.mean()) - (o - o.mean())              # d0: each waveform's own mean removed
    noise, noise0 = np.ptp(a, axis=1), np.ptp(a - a.mean(axis=0), axis=1)
    means = a.mean(axis=0)
    return {"mean_old": float(o.mean()), "means_new": [float(x) for x in means], "spread_of_means": float(np.ptp(means)),
            "shift": float(d.mean()), "n_negative": int((d < 0).sum()), "n_positions": len(d),
            "median_abs_diff": float(np.median(np.abs(d))), "noise_median": float(np.median(noise)), "noise_largest": float(noise.max()),
            "within_noise": int((np.abs(d) <= noise).sum()),
            "median_abs_diff_no_offset": float(np.median(np.abs(d0))), "noise_median_no_offset": float(np.median(noise0)),
            "within_noise_no_offset": int((np.abs(d0) <= noise0).sum())}


def print_full_model(df: pd.DataFrame, ref: pd.DataFrame, new: pd.DataFrame):
    """The tables of the report's section 5.6: re-meshing noise against the difference between the programs, and the common offset."""
    ref_full, new_full = ref[ref.model == FULL_LABEL], new[new.model == FULL_LABEL]
    print("\nthe hand-drawn full 360-degree model: torque noise from re-meshing (rotor turned by 0, 0.02 and 0.04 electrical degrees) against the "
          "difference between FEMM 4.2 and FEMM 5.0, in mN.m")
    for label, current in (("no load (cogging)", 0.0), ("10 A", 10.0)):
        a = df[df.I == current].pivot(index="theta_e", columns="offset_deg", values="T")
        o = ref_full[ref_full.I == current].set_index("theta_e")["T"].reindex(a.index)
        n = new_full[new_full.I == current].set_index("theta_e")["T"].reindex(a.index)
        noise = (a.max(axis=1) - a.min(axis=1)) * 1e3
        d = ((n - o).abs() * 1e3)
        same = float(np.abs(a[0.0].values - n.values).max())
        print(f"  {label:18s} noise over three meshes: median {noise.median():6.2f}, largest {noise.max():6.2f}   |FEMM 5.0 - FEMM 4.2|: median "
              f"{d.median():6.2f}, largest {d.max():6.2f}   positions where the difference is within the noise: {int((d <= noise).sum())} of {len(d)}   "
              f"(offset 0 against validate.csv: largest difference {same:.1e} N.m)")
        if current == 0.0:
            print("  cogging peak to peak over the 24 positions [mN.m]: FEMM 4.2 %.2f; FEMM 5.0 at offset 0 / 0.02 / 0.04: %s" % (
                np.ptp(o) * 1e3, " / ".join(f"{np.ptp(a[c]) * 1e3:.2f}" for c in FULL_OFFSETS)))
        s = offset_analysis(o.values, n.values, a.values)
        print(f"    mean of the 24 positions [mN.m]: FEMM 4.2 {s['mean_old']:+.2f}; FEMM 5.0 on its three meshes "
              f"{' / '.join(f'{x:+.2f}' for x in s['means_new'])} (they differ by {s['spread_of_means']:.2f})")
        print(f"    common shift (mean of FEMM 5.0 - FEMM 4.2) {s['shift']:+.2f}, negative at {s['n_negative']} of {s['n_positions']} positions: "
              f"{100 * abs(s['shift']) / s['median_abs_diff']:.0f} % of the median |difference| ({s['median_abs_diff']:.2f}), "
              f"{abs(s['shift']) / s['noise_median']:.1f} times the median noise ({s['noise_median']:.2f}), "
              f"{abs(s['shift']) / s['spread_of_means']:.0f} times the spread of the three means")
        print(f"    with each waveform's own mean removed: median |difference| {s['median_abs_diff_no_offset']:.2f} against median noise "
              f"{s['noise_median_no_offset']:.2f}; positions within the noise {s['within_noise_no_offset']} of {s['n_positions']} "
              f"(with the offset: {s['within_noise']})")
    print("  the 90-degree slice at the same 24 positions (it has no mesh asymmetry between its sectors):")
    for label, current in (("no load", 0.0), ("10 A", 10.0)):
        so = ref[(ref.model == SLICE_LABEL) & (ref.I == current)]["T"]
        sn = new[(new.model == SLICE_LABEL) & (new.I == current)]["T"]
        line = f"    {label:8s} mean torque [mN.m]: FEMM 4.2 {so.mean() * 1e3:+.3f}, FEMM 5.0 {sn.mean() * 1e3:+.3f}"
        if current == 0.0:
            line += f";  cogging peak to peak [mN.m]: FEMM 4.2 {np.ptp(so) * 1e3:.2f}, FEMM 5.0 {np.ptp(sn) * 1e3:.2f} ({100 * (np.ptp(sn) / np.ptp(so) - 1):+.1f} %)"
        print(line)


def run_fullmodel(old: Path) -> pd.DataFrame:
    """Turn the rotor of the full model by two tiny extra angles. FEMM re-meshes at every angle, so each is a different mesh
    of the same physics (0.02 electrical degrees is 0.005 mechanical degrees). The torque's spread over the three is the mesh
    noise of this model, and it is what FEMM 4.2's and FEMM 5.0's numbers must be compared with. Then split the difference between
    the programs into one common offset of the 24-position mean and what is left (offset_analysis)."""
    import validate_full as vf
    (OUT / "full").mkdir(parents=True, exist_ok=True)
    thetas = [float(t) for t in vf.THETAS]
    tasks, meta = [], []
    for k, (offset, current, theta) in enumerate((o, i, t) for o in FULL_OFFSETS for i in (0.0, vf.I_LOAD) for t in thetas):
        tasks.append((2000 + k, theta + offset, current))
        meta.append((offset, current, theta))
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        got = list(pool.map(full_worker, tasks))
    df = pd.DataFrame([{"offset_deg": o, "I": i, "theta_e": t, "psi_a": g["psi_a"], "T": g["T"], "nodes": g["nodes"]}
                       for (o, i, t), g in zip(meta, got)])
    df.to_csv(OUT / "full_model.csv", index=False)
    print(f"fullmodel: {len(df)} solves in {time.perf_counter() - t0:.0f} s on {WORKERS} workers")
    ref = pd.read_csv(Path(old) / "validate" / "validate.csv", float_precision="round_trip")
    new = pd.read_csv(P / "out" / "validate" / "validate.csv", float_precision="round_trip")
    print_full_model(df, ref, new)
    return df


def dense_analysis(torque) -> dict:
    """A torque waveform sampled at 96 positions of one electrical cycle: its mean, the means of its four interleaved 24-position
    grids (the rotor shifted by 0, 3.75, 7.5 and 11.25 degrees) and the amplitudes of some harmonics, all in the unit of `torque`.

    The mean over the 96 positions is the mean of the waveform. A 24-position grid also contains every harmonic whose number is a
    multiple of 24 (the 24th above all, then the 48th), because 24 positions cannot tell such a harmonic from a constant; a grid
    shifted by half that harmonic's period (7.5 degrees) sees it with the opposite sign, and one shifted by a quarter (3.75 or
    11.25 degrees) not at all. So the four grid means differ from one another by exactly this "sampling" part, and the spread of
    the four is a measure of how far the choice of 24 positions can move a mean.
    """
    t = np.asarray(torque, float)
    h = harmonic_amplitudes(t)
    return {"mean": float(t.mean()), "grid_means": [float(t[k::4].mean()) for k in range(4)],
            "harmonics": {k: float(h[k]) for k in (6, 12, 18, 24, 30, 36)}}


def run_fulldense(old: Path) -> pd.DataFrame:
    """The full 360-degree model with FEMM 5.0 at 96 rotor positions per electrical cycle, at no load and at 10 A.

    The question. The mean of the 24 no-load torques is +1.69 mN.m on FEMM 4.2 and -0.73 on FEMM 5.0, and the 10 A means
    differ in the same way (+2.58 against -1.18). Is that difference the sampling? 24 positions per cycle cannot tell the 24th
    harmonic of the torque from a constant, so each program's 24-position mean contains it. 96 positions can (up to the 96th).

    What is compared. FEMM 5.0's mean over the 96 positions, against its own four interleaved 24-position grids (the first of
    them, shifted by 0 degrees, is validate_full.py's grid): the difference between a grid's mean and the 96-position mean is
    what the sampling does to that mean. If it were the cause of the 2.4 mN.m between the programs, it would be of that size.
    It cannot be compared with FEMM 4.2's own sampling, which was solved on the 24-position grid only.

    A caution about the "true" mean. For a rotor with the symmetry of one electrical cycle (the 90-degree slice) the energy
    returns to itself after the cycle, so the mean torque over it is 0 at no load and at any fixed stator field. The hand-drawn
    rotor has two magnets that point 5 degrees off the pattern (check_vs_fem.py), so after one electrical cycle it is not the same
    rotor, and the mean over one cycle need not be exactly 0: a small non-zero mean is not by itself a fault of the solver.
    """
    import validate_full as vf
    (OUT / "full").mkdir(parents=True, exist_ok=True)
    thetas = [float(t) for t in np.arange(0, 360, FULL_DENSE_STEP)]
    tasks, meta = [], []
    for k, (current, theta) in enumerate((i, t) for i in (0.0, vf.I_LOAD) for t in thetas):
        tasks.append((3000 + k, theta, current))
        meta.append((current, theta))
    t0 = time.perf_counter()
    with ProcessPoolExecutor(WORKERS) as pool:
        got = list(pool.map(full_worker, tasks))
    df = pd.DataFrame([{"I": i, "theta_e": t, "psi_a": g["psi_a"], "T": g["T"], "nodes": g["nodes"]} for (i, t), g in zip(meta, got)])
    df.to_csv(OUT / "full_dense.csv", index=False)
    print(f"fulldense: {len(df)} solves in {time.perf_counter() - t0:.0f} s on {WORKERS} workers")
    ref = pd.read_csv(Path(old) / "validate" / "validate.csv", float_precision="round_trip")
    ref = ref[ref.model == FULL_LABEL]
    print(f"\nthe hand-drawn full 360-degree model on FEMM 5.0 at {len(thetas)} positions per electrical cycle ({FULL_DENSE_STEP} degrees apart), torque in mN.m")
    for label, current in (("no load", 0.0), ("10 A", vf.I_LOAD)):
        a = dense_analysis(df[df.I == current].sort_values("theta_e")["T"].values * 1e3)
        o = ref[ref.I == current]["T"].mean() * 1e3
        print(f"  {label:8s} mean over all {len(thetas)} positions {a['mean']:+.3f};  the four interleaved 24-position grids, rotor shifted by "
              f"0 / 3.75 / 7.5 / 11.25 degrees: {' / '.join(f'{x:+.3f}' for x in a['grid_means'])} (they differ by {np.ptp(a['grid_means']):.2f});  "
              f"FEMM 4.2's 24-position mean (shift 0): {o:+.3f}")
        print(f"           harmonics of the electrical angle: " + ", ".join(f"{k}th {v:.2f}" for k, v in a["harmonics"].items()))
        shift, sampling = a["grid_means"][0] - o, a["grid_means"][0] - a["mean"]
        print(f"           FEMM 5.0's 24-position mean minus FEMM 4.2's: {shift:+.2f}.  The 24-position sampling moves FEMM 5.0's mean by {sampling:+.2f} "
              f"({100 * abs(sampling / shift):.0f} % of that); what is left, {shift - sampling:+.2f}, is not the sampling")
    return df


# ----------------------------------------------------------------------------------------------- the picture
def make_figure(s: pd.DataFrame, cal: pd.DataFrame):
    """Quantities against node count: FEMM 5.0's two ladders (lines), FEMM 4.2 (stars, at its own node counts)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = [("cogging_pp_mNm", "cogging torque, peak to peak [mN.m]"), ("harmonic24_mNm", "torque harmonic 24 [mN.m]"),
              ("psi_pm_mWb", "magnet flux linkage psi_PM [mWb]"), ("T_mean_gamma0_Nm", "mean torque, 10 A on the q axis [N.m]"),
              ("psi_d_gamma90_mWb", "psi_d, 10 A on the d axis (almost cancelled) [mWb]"), ("gap_max_pct", "air-gap B: worst sample against FEMM 4.2 [% of peak]")]
    fig, axes = plt.subplots(2, 3, figsize=(13, 6.8))
    ref = s.iloc[0]
    for ax, (col, title) in zip(axes.ravel(), panels):
        for prefix, label, colour in (("c-", "FEMM 5.0, smart mesh off", "#2a78d6"), ("f-", "FEMM 5.0, smart mesh on", "#eb6834")):
            d = s[s.tag.str.startswith(prefix)].sort_values("nodes")
            ax.plot(d.nodes, d[col], "-o", ms=4, color=colour, label=label)
        if col != "gap_max_pct":
            ax.plot([README_NODES["fine"]], [ref[col]], "*", ms=14, color="#0b0b0b", label=f"FEMM 4.2 (smart mesh on, {README_NODES['fine']:,} nodes)")
        ax.axvline(README_NODES["fine"], color="#c3c2b7", lw=0.8)
        ax.set_xscale("log")
        ticks = [8e3, 12e3, 18e3, 28e3, 45e3, 70e3]
        ax.set_xticks(ticks)
        ax.set_xticklabels([f"{t / 1e3:.0f}k" for t in ticks], fontsize=8)
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
        ax.set_title(title, fontsize=9.5, loc="left")
        ax.set_xlabel("mesh nodes", fontsize=8.5)
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=8, frameon=False)
    fig.suptitle("The same model on FEMM 5.0 at several mesh densities, and FEMM 4.2 at its own", x=0.02, ha="left", fontsize=11.5)
    fig.tight_layout()
    path = OUT / "mesh_convergence.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    shutil.copyfile(path, HERE / path.name)                # the copy that is kept beside this script, so a re-run cannot leave it stale
    print(f"wrote {path} and a copy in {HERE}")


# ----------------------------------------------------------------------------------------------- the command
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Solve the model on FEMM 5.0 at other mesh densities and compare with FEMM 4.2.")
    ap.add_argument("old", help="FEMM 4.2's out folder (read only)")
    ap.add_argument("--only", default=",".join(PARTS), help="comma-separated parts; they always run in this order: " + ", ".join(PARTS))
    args = ap.parse_args(argv)
    old = Path(args.old)
    if not (old / "barrier" / "phase1.csv").is_file():
        print(f"{old} does not look like FEMM 4.2's out folder (no barrier/phase1.csv)", file=sys.stderr)
        return 2
    parts = [x.strip() for x in args.only.split(",") if x.strip()]
    unknown = set(parts) - set(PARTS)
    if unknown:
        print(f"unknown part: {', '.join(sorted(unknown))}", file=sys.stderr)
        return 2
    femm_session.lower_priority()                          # every worker inherits it, so the desktop stays quick
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "fem").mkdir(exist_ok=True)
    if "calibrate" in parts:
        cal = run_calibration()
        print(cal.to_string(index=False))
    else:
        cal = pd.read_csv(OUT / "calibration.csv")
    if "slice" in parts:
        df = run_slice(cal)
        summarise_check_femm(df, cal)
        print_slice(summarise_slice(df, cal, old))
    if "states" in parts:
        run_states(cal, old)
    if "post42" in parts:
        run_post42(old)
    if "post42full" in parts:
        run_post42full(old)
    if "fullmodel" in parts:
        run_fullmodel(old)
    if "fulldense" in parts:
        run_fulldense(old)
    if "noise" in parts:
        run_noise(cal, old)
    if "flagged" in parts:
        run_flagged(cal, old)
    if "timing" in parts:
        run_timing(cal)
    if "candidate" in parts:
        run_candidate(old)
    if "ulp" in parts:
        run_ulp()
    if "figure" in parts:
        make_figure(pd.read_csv(OUT / "slice_summary.csv"), cal)
    return 0


if __name__ == "__main__":
    sys.exit(main())
