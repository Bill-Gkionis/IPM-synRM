"""provenance.py: where the report's hand-made numbers come from. Four small analyses of files that other scripts wrote.

Why this exists. compare.py and mesh_checks.py print most of the numbers in the report (the project note "4.2 vs 5.0",
"the report" below). A few numbers in it are not printed by either, because
they are small analyses of the files those scripts left behind. A number nobody can make again is a number nobody can check, so
each of them is made here, from the files, and printed in the words the report uses:

  offset   F1 (report 2.12): the flux linkage psi_d of the torque-vs-angle grid is higher on FEMM 5.0 by one constant
           amount in every state. How big is it, and is it the same in all of the 52 states? (barrier/phase1.csv, both runs)
  gap      F4 (report 2.12): the air-gap flux density B_r of the six barrier-study cases differs most at its steep samples.
           Where is the worst sample, how steep is the waveform there, is the difference a shift in angle, and how does it
           compare with FEMM 5.0 against its own denser meshes? (barrier/gap_*.csv of both runs, and the ladder's gap files
           that mesh_checks.py `slice` writes)
  pandas   report 3: how many of the optimizer's numbers does pandas' default CSV reader read one unit in the last place off?
           (out/opt/evaluations.csv)
  picture  report 2.1: where do the pixels of check_vs_fem.png differ, by how much, and in which direction?

Run (from the Python to FEMM folder; no FEMM, about 3 seconds; run mesh_checks.py first for the last table of `gap`):
    uv run python validation/provenance.py femm42-out [--only offset,gap,pandas,picture]
femm42-out is FEMM 4.2's run (a checked copy of the `out` folder it left); FEMM 5.0's is `out`. Nothing is written.

Method, in one line each.
  offset   For every (current, current angle) of the grid: the mean over the six rotor positions of psi_d, on both runs. The
           difference new - old is the offset; its mean, its spread and how many of the 52 states it is positive in are printed.
  gap      |B_new - B_old| / peak of B_old at the 360 samples (0.25 degrees apart over the 90-degree sector). The slope is the
           central difference of FEMM 4.2's waveform. The best shift moves FEMM 5.0's waveform along the angle (a periodic cubic
           spline, from -0.3 to +0.3 degrees in steps of 0.01) and keeps the shift with the smallest rms difference.
  pandas   Every parameter of the optimizer's feasible designs is read twice: with Python's float() (exact) and with
           pandas.read_csv's default; the two are compared bit for bit.
  picture  The two images are compared pixel by pixel in their red, green, blue and alpha values.

Numbers to expect (2026-11-05, this PC), the same as in the report:
  offset   with band: +0.0392 mWb (spread 0.0041), positive in 52 of 52 states; no band: +0.0464 (0.0044), 52 of 52. With band that
           is 0.32 % of the magnet flux (12.09 mWb) and 0.11 % of the largest flux linkage of the grid (34.5 mWb).
  gap      six cases: median difference 0.13 to 0.24 % of the peak, worst sample 1.26 to 1.83 %, and the best shift in angle is
           +0.00 degrees in every case. The with-band no-load case's worst sample is at 52.5 degrees, on a steep edge of the
           waveform (0.140 T per degree; the steepest anywhere is 0.152); the other five have theirs at 51.25 degrees, on the foot
           of the same edge. The large differences lean towards the steep samples but are not confined to them: of the 90 largest
           differences 37 to 42 are among the 90 steepest samples (a random choice would give 22 or 23); the five largest of the first
           case are at 52.50 (-1.83 %), 33.50 (-1.78 %), 51.25 (+1.71 %), 7.00 (-1.59 %) and 76.50 degrees (+1.56 %, where the waveform
           is almost flat: 0.012 T per degree). A shift of 0.02 degrees at the
           steepest part would move B by 1.0 to 1.2 % of the peak. Against FEMM 5.0's own 70k-node field the largest differences
           are 2.54 % (FEMM 5.0's natural mesh), 2.04 % (FEMM 4.2) and 0.53 % (FEMM 5.0 at 45k nodes).
  pandas   27 of the 192 parameters (48 feasible designs, four parameters each).
  picture  850 of 1,901,250 pixels, all in the right half (the stator's), none in the left half; at most 43 of 255 in one colour
           channel; 415 pixels brighter, 435 darker.

Self-check (the answers follow).
  1. Two waveforms differ by 0.02 degrees in angle, and the steepest part of the waveform rises 0.16 T per degree. By how much
     do they differ there, and why is that where two meshes of the same field disagree most?
  2. pandas reads the number 39.227554466990604 from a CSV file as a different float than Python's float() does. A parameter of
     the model that is one unit in the last place different (about 1e-15 mm) is "no change at all" to the motor. Why does it
     still matter here?
Answers.
  1. 0.02 x 0.16 = 0.0032 T. Wherever the waveform is flat, a small shift in angle changes B by almost nothing; where it is steep,
     the same shift changes B by the slope times the shift. Two meshes put their triangle corners in slightly different places,
     so the field they solve is slightly displaced in effect, and the difference tends to show up where the field changes fastest
     (not only there: see the `large&steep` column of the gap table).
  2. FEMM meshes the drawing it is given, and one unit in the last place is enough to make it cut the drawing into different
     triangles (mesh_checks.py `noise` measures how much a number moves when only that happens: torque 0.04 %). A design read back
     through pandas' default reader is therefore not exactly the design that was solved, and solving it again does not give the
     CSV's numbers. Python's float(), or pandas with float_precision="round_trip", reads exactly.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
P = HERE.parent
for _folder in (str(P), str(HERE)):
    if _folder not in sys.path:
        sys.path.insert(0, _folder)
from mesh_checks import gap_deviation  # noqa: E402  (the one definition of "% of the old waveform's peak")

CASES = [(case, variant) for case in ("noload", "qaxis", "mtpa") for variant in ("with_band", "no_band")]
SHIFTS = np.round(np.arange(-30, 31) / 100, 2)               # degrees tried for the best shift: -0.30 to +0.30 in steps of 0.01 (0.0 is exact)
LADDER = ("f-natural", "f-45k", "f-70k")                     # FEMM 5.0's smart-mesh fields of mesh_checks.py `slice`; 70k is the densest


# ----------------------------------------------------------------------------------------------- F1: the offset of psi_d
def psi_d_offset(old: Path, new: Path) -> list[dict]:
    """For each variant: the constant by which FEMM 5.0's psi_d [mWb] of the torque-vs-angle grid is higher, state by state."""
    po = pd.read_csv(Path(old) / "barrier" / "phase1.csv")
    pn = pd.read_csv(Path(new) / "barrier" / "phase1.csv")
    summary = pd.read_csv(Path(old) / "barrier" / "summary.csv", index_col=0)
    rows = []
    for variant in ("with band", "no band"):
        grid = {}
        for label, df in (("old", po), ("new", pn)):
            g = df[(df.study == "gamma") & (df.variant == variant)]
            grid[label] = g.groupby(["I", "gamma"])[["psi_d", "psi_q"]].mean() * 1e3
        diff = (grid["new"].psi_d - grid["old"].psi_d)
        largest = float(max(grid["old"].psi_d.abs().max(), grid["old"].psi_q.abs().max()))
        pm = float(summary.loc["psi_pm_mWb", variant])
        rows.append({"variant": variant, "states": len(diff), "mean": float(diff.mean()), "std": float(diff.std()),
                     "positive": int((diff > 0).sum()), "pct_of_magnet_flux": 100 * float(diff.mean()) / pm, "magnet_flux": pm,
                     "pct_of_largest": 100 * float(diff.mean()) / largest, "largest": largest,
                     "psi_d_min": float(grid["old"].psi_d.min()), "psi_d_max": float(grid["old"].psi_d.max())})
    return rows


# ----------------------------------------------------------------------------------------------- F4: the air-gap flux density
def slope_per_degree(deg, b) -> np.ndarray:
    """|dB/d(angle)| [T per degree] by central differences, the waveform taken as periodic over the sector."""
    step = float(deg[1] - deg[0])
    return np.abs(np.roll(b, -1) - np.roll(b, 1)) / (2 * step)


def best_shift(deg, old_b, new_b, shifts=SHIFTS) -> float:
    """The angle [degrees] by which FEMM 5.0's waveform has to be moved along the angle to come nearest (rms) to FEMM 4.2's.

    The waveform is periodic over the sector (90 degrees), so a periodic cubic spline gives its value between samples.
    """
    deg, old_b, new_b = np.asarray(deg, float), np.asarray(old_b, float), np.asarray(new_b, float)
    period = float(deg[-1] - deg[0] + (deg[1] - deg[0]))
    spline = CubicSpline(np.append(deg, deg[0] + period), np.append(new_b, new_b[0]), bc_type="periodic")
    rms = [float(np.sqrt(np.mean((spline((deg + s - deg[0]) % period + deg[0]) - old_b) ** 2))) for s in shifts]
    return float(shifts[int(np.argmin(rms))])


def gap_case(old_gap: pd.DataFrame, new_gap: pd.DataFrame) -> dict:
    """Where the air-gap flux density of two runs differs most, how steep it is there, and whether a shift in angle explains it."""
    deg, b_old, b_new = old_gap.deg.values, old_gap.br.values, new_gap.br.values
    peak = float(np.abs(b_old).max())
    d = np.abs(b_new - b_old) / peak * 100
    slope = slope_per_degree(deg, b_old)
    worst = int(np.argmax(d))
    quarter = len(d) // 4                     # a quarter of the samples: 90 of 360. A random choice of two such sets shares a quarter of them
    overlap = len(set(np.argsort(-d)[:quarter]) & set(np.argsort(-slope)[:quarter]))
    signed = (b_new - b_old) / peak * 100
    top = [(float(deg[i]), float(signed[i]), float(slope[i])) for i in np.argsort(-d)[:5]]       # (angle, signed % of peak, slope there)
    return {"top": top, "median_pct": float(np.median(d)), "max_pct": float(d.max()), "worst_deg": float(deg[worst]), "slope_at_worst": float(slope[worst]),
            "steepest": float(slope.max()), "peak": peak, "shift_deg": best_shift(deg, b_old, b_new),
            "pct_of_peak_per_0p02_deg_at_steepest": 100 * 0.02 * float(slope.max()) / peak, "quarter": quarter, "steep_overlap": overlap}


def gap_ladder(old: Path, new_mesh_dir: Path) -> dict | None:
    """The with-band no-load field: FEMM 5.0's natural fine mesh, FEMM 4.2 and FEMM 5.0 at 45k nodes, each against FEMM 5.0 at 70k.

    (largest, rms) of |difference| as a percentage of the 70k field's peak, from mesh_checks.gap_deviation. None where the
    ladder's gap files are not there (mesh_checks.py `slice` has not run).
    """
    files = {tag: Path(new_mesh_dir) / f"gap_{tag}.csv" for tag in LADDER}
    if not all(f.is_file() for f in files.values()):
        return None
    ref = pd.read_csv(files["f-70k"])
    fields = {"FEMM 5.0, its natural mesh": pd.read_csv(files["f-natural"]), "FEMM 4.2": pd.read_csv(Path(old) / "barrier" / "gap_noload_with_band.csv"),
              "FEMM 5.0, 45k nodes": pd.read_csv(files["f-45k"])}
    return {name: gap_deviation(ref, field) for name, field in fields.items()}


# ----------------------------------------------------------------------------------------------- the pandas trap
def pandas_misreads(evaluations: Path, names=("mag_y", "mag_x", "v_half_deg", "web")) -> tuple[int, int]:
    """(parameters that pandas' default CSV reader reads one bit off, parameters in all) over the feasible designs of the optimizer."""
    with open(evaluations, newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r["feasible"] == "True"]
    default = pd.read_csv(evaluations)
    default = default[default.feasible.astype(str) == "True"].reset_index(drop=True)
    wrong = sum(float(rows[i][n]) != float(default.loc[i, n]) for i in range(len(rows)) for n in names)
    return int(wrong), len(rows) * len(names)


# ----------------------------------------------------------------------------------------------- the picture
def picture_difference(old_png: Path, new_png: Path) -> dict:
    """How two images of one size differ: pixels, halves, the largest change of one colour channel, brighter and darker."""
    from PIL import Image
    a = np.asarray(Image.open(old_png).convert("RGBA")).astype(int)
    b = np.asarray(Image.open(new_png).convert("RGBA")).astype(int)
    if a.shape != b.shape:
        return {"same_size": False}
    differ = (a != b).any(axis=2)
    half = a.shape[1] // 2
    brightness = lambda x: 0.299 * x[..., 0] + 0.587 * x[..., 1] + 0.114 * x[..., 2]
    change = brightness(b) - brightness(a)
    return {"same_size": True, "differ": int(differ.sum()), "total": int(differ.size), "left": int(differ[:, :half].sum()),
            "right": int(differ[:, half:].sum()), "max_channel": int(np.abs(a - b).max()), "brighter": int((change > 0).sum()),
            "darker": int((change < 0).sum())}


# ----------------------------------------------------------------------------------------------- the command
def print_offset(old: Path, new: Path):
    print("\nF1: psi_d of the torque-vs-angle grid, FEMM 5.0 minus FEMM 4.2, state by state (mean over six rotor positions)")
    for r in psi_d_offset(old, new):
        print(f"  {r['variant']:9s}: offset {r['mean']:+.4f} mWb (spread {r['std']:.4f}), positive in {r['positive']} of {r['states']} states; "
              f"{r['pct_of_magnet_flux']:.2f} % of the magnet flux ({r['magnet_flux']:.2f} mWb), {r['pct_of_largest']:.2f} % of the largest "
              f"flux linkage ({r['largest']:.1f} mWb); psi_d runs from {r['psi_d_min']:.1f} to {r['psi_d_max']:.1f} mWb")


def print_gap(old: Path, new: Path):
    print("\nF4: air-gap flux density B_r, FEMM 5.0 against FEMM 4.2, 360 samples over the 90-degree sector")
    print(f"  {'case':22s} {'median':>8s} {'largest':>8s} {'at':>7s} {'slope there':>12s} {'steepest':>9s} {'best shift':>11s} {'large&steep':>12s}   "
          "(percent of the peak; slopes in T per degree; large&steep: how many of the 90 largest differences are among the 90 steepest samples)")
    for case, variant in CASES:
        name = f"gap_{case}_{variant}.csv"
        r = gap_case(pd.read_csv(Path(old) / "barrier" / name), pd.read_csv(Path(new) / "barrier" / name))
        print(f"  {case + ', ' + variant.replace('_', ' '):22s} {r['median_pct']:8.3f} {r['max_pct']:8.3f} {r['worst_deg']:6.2f}d {r['slope_at_worst']:12.3f} "
              f"{r['steepest']:9.3f} {r['shift_deg']:+10.2f}d {r['steep_overlap']:8d} of {r['quarter']}   a shift of 0.02 degrees at the steepest part moves B by "
              f"{r['pct_of_peak_per_0p02_deg_at_steepest']:.2f} % of the peak")
    first = gap_case(pd.read_csv(Path(old) / "barrier" / "gap_noload_with_band.csv"), pd.read_csv(Path(new) / "barrier" / "gap_noload_with_band.csv"))
    print("  the five largest differences of the first case (angle, percent of the peak, slope there in T per degree): "
          + "; ".join(f"{a:.2f}d {pct:+.2f} % ({sl:.3f})" for a, pct, sl in first["top"]))
    ladder = gap_ladder(old, new / "validation" / "mesh")
    if ladder is None:
        print("  (the ladder's gap files are not there: run mesh_checks.py with the part `slice` first)")
        return
    print("  against FEMM 5.0's own 70k-node field (largest, rms; percent of its peak), with band, no load:")
    for name, (largest, rms) in ladder.items():
        print(f"    {name:28s} {largest:6.2f} %   rms {rms:5.2f} %")


def print_pandas(new: Path):
    wrong, total = pandas_misreads(Path(new) / "opt" / "evaluations.csv")
    print(f"\nthe pandas trap: pandas' default CSV reader reads {wrong} of the {total} parameters of the optimizer's feasible designs one bit off")


def print_picture(old: Path, new: Path):
    r = picture_difference(Path(old) / "check_vs_fem.png", Path(new) / "check_vs_fem.png")
    if not r["same_size"]:
        print("\ncheck_vs_fem.png: the two pictures have different sizes")
        return
    print(f"\ncheck_vs_fem.png: {r['differ']:,} of {r['total']:,} pixels differ; in the left half {r['left']}, in the right half {r['right']}; "
          f"at most {r['max_channel']} of 255 in one colour channel; {r['brighter']} pixels brighter, {r['darker']} darker")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Make the report's hand-made numbers again from the files of both runs.")
    ap.add_argument("old", help="FEMM 4.2's out folder (femm42-out), read only")
    ap.add_argument("new", nargs="?", default="out", help="FEMM 5.0's out folder (default: out)")
    ap.add_argument("--only", default="offset,gap,pandas,picture", help="comma-separated parts (default: all four)")
    args = ap.parse_args(argv)
    old, new = Path(args.old), Path(args.new)
    parts = [x.strip() for x in args.only.split(",") if x.strip()]
    unknown = set(parts) - {"offset", "gap", "pandas", "picture"}
    if unknown or not old.is_dir() or not new.is_dir():
        print(f"unknown part: {', '.join(sorted(unknown))}" if unknown else f"not a folder: {old if not old.is_dir() else new}", file=sys.stderr)
        return 2
    if "offset" in parts:
        print_offset(old, new)
    if "gap" in parts:
        print_gap(old, new)
    if "pandas" in parts:
        print_pandas(new)
    if "picture" in parts:
        print_picture(old, new)
    return 0


if __name__ == "__main__":
    sys.exit(main())
