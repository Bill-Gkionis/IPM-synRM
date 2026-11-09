"""m19_femm5.py: FEMM 5.0's "M-19 14 mil" against FEMM 4.2's "M-19 Steel" and against the maker's own curve.

Why. The motor model was built on FEMM 4.2's library block "M-19 Steel" (47 points of its B-H curve, lamination fill 0.98).
FEMM 5.0's library has no block of that name. Its nearest, "M-19 14 mil" (61 points, fill 0.955), is a different curve, so the
model carries FEMM 4.2's block itself (motorgeo.M19_STEEL_FEMM42, rule R5) and the two programs solve the same steel. This script
puts numbers on the question "what would change if the model used FEMM 5.0's block instead?", so that it can be decided whether
to adopt it. The decision, the numbers and what they mean are in the project note M19
("the note" below), which also holds the open question "adopt M-19 14 mil?".

What it prints, in order.
  1. The H a curve needs for a flux density B of 0.5, 1.0, 1.5 and 1.8 T, for the three curves, and the difference of each
     library curve from the maker's: materials/research/m19_dc_maker.csv, the DC curve of the steel's maker (16 points; the
     author of FEMM 4.2's block took its numbers from the same source). The maker's table has no point at 0.5 T, so there its
     value is read off its 16 points with a smooth monotone interpolation (PCHIP on ln H), good to about 1 %.
  2. check_femm's four solves (rotor position 0, currents (0, 0), (0, 10), (-6, 8) and (+6, 8) A) on both meshes, with the
     steel of the model changed, and nothing else: FEMM 4.2's curve and fill (what the model uses), FEMM 5.0's "M-19 14 mil"
     as the library has it (its curve, its fill 0.955), and the two mixed (14 mil's curve with 4.2's fill, 4.2's curve with
     14 mil's fill) to tell which of the two differences does the work. A last row gives 4.2's curve with 14 mil's
     conductivity and lamination thickness: those only matter when the field alternates, so this problem (DC) must not change.

Run (from the Python to FEMM folder; it starts one FEMM of its own, about 20 s):
    uv run python validation/m19_femm5.py
Outputs: the tables above on screen, and out/validation/m19_femm5.csv, out/validation/m19_curves.csv and
out/validation/m19_curves.png. The picture is also copied to validation/m19_curves.png, beside this script, so the copy that
is kept cannot go stale. Nothing is written to the project notes, and the model (motorgeo.py) is not changed: the steel is swapped for
the length of one run by replacing motorgeo.define_m19_steel in this process.

Numbers to expect (2026-11-04, FEMM 5.0.20261101; the note, section "The effect on check_femm's four solves", has them too).
  H [A/m] at B           0.5 T            1.0 T           1.5 T           1.8 T
    FEMM 4.2 "M-19 Steel"   51.44 (-0.48 %)   106.20 (+0.19 %)  1108.33 (-0.06 %)  9708.82 (-0.02 %)   against the maker's curve
    FEMM 5.0 "M-19 14 mil"  56.09 (+8.52 %)   118.12 (+11.43 %) 1107.38 (-0.15 %)  12085.74 (+24.45 %)
    the maker               51.69 (read between its points; it has none at 0.5 T), 106.00, 1109.00, 9711.00
  check_femm's torque T_age [N.m] and psi_d [mWb], coarse mesh (9,580 nodes), the change from the model's own steel in brackets:
                                           T (0, 10) A        T (-6, 8) A       T (+6, 8) A        psi_d, no load
    4.2 curve, fill 0.98 (the model's)     0.8087             0.9104            0.3944             11.991
    14 mil curve, fill 0.955 (as in library) 0.8388 (+3.72 %)   0.9310 (+2.26 %)  0.4218 (+6.93 %)   12.472 (+4.02 %)
    14 mil curve, fill 0.98                0.8174 (+1.08 %)   0.9162 (+0.63 %)  0.4032 (+2.23 %)   12.133 (+1.19 %)
    4.2 curve, fill 0.955                  0.8303 (+2.67 %)   0.9254 (+1.65 %)  0.4133 (+4.78 %)   12.332 (+2.85 %)
    4.2 curve, 14 mil's sigma and d_lam    unchanged, +0.00 % everywhere (bit-identical)
  The fine mesh (17,934 nodes) says the same: 14 mil as the library has it gives +3.75 %, +2.31 %, +7.00 % in the three torques and
  +4.02 % in psi_d. So the steel moves the torque by 2 to 7 % and psi_d by 4 %, more than any difference between FEMM 4.2 and
  FEMM 5.0 that compare.py finds, and the lamination fill (+2.7 % alone) does more of it than the curve (+1.1 % alone).

Self-check (the answers follow).
  1. At which flux density do the two library curves differ most from the maker's, and why does that matter more for the
     bridges and the stator teeth than for the middle of the magnets' flux path?
  2. Why can the fill factor change a DC torque at all, when it only says how much of a lamination is steel?
  3. 14 mil's conductivity and lamination thickness differ from 4.2's. Why must the last row of the second table equal the first?

Answers.
  1. At 1.8 T. There FEMM 5.0's curve needs 24 % more H than the maker's (and +11 % at 1.0 T, +8.5 % at 0.5 T), while FEMM 4.2's is
     within 0.5 % everywhere. Between 0.5 and 1.5 T the steel is far from saturation and takes almost no magnetizing force
     (mu_r above 1000), so 10 % more or less of an already small H changes little. The thin bridges, the web and the tooth tips
     are saturated, around 1.8 T and above, where the curve is flat: a small change of B needs a large change of H, so a curve
     that wants 24 % more H at the same B is a much weaker steel exactly there, and those parts decide how much of the magnets'
     flux is lost around them.
  2. Because FEMM applies it to the curve. A lamination stack is steel plus insulation, and the field strength H is the same in
     both, so at one H the stack carries fill x B_steel: a fill of 0.955 is a weaker steel at every point of the curve (the
     saturated bridges carry less flux than with 0.98), and the torque and the flux linkage follow. Whether this is the whole
     reason was not investigated: the table only shows that fill 0.98 to 0.955 changes psi_d by 2.9 %.
  3. A DC (magnetostatic) problem has no eddy currents. Conductivity and lamination thickness only enter the terms that
     describe currents induced in the steel by a changing field; with a field that does not change they are never used, so FEMM
     gives the same numbers to the last bit.
"""
from __future__ import annotations

import csv
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator

HERE = Path(__file__).resolve().parent
P = HERE.parent
for _folder in (str(P), str(P / "materials")):
    if _folder not in sys.path:
        sys.path.insert(0, _folder)
import common as materials_common  # noqa: E402  (materials/common.py: read_matlib)
import femm_session  # noqa: E402
import motorgeo as mg  # noqa: E402

MAKER = P / "materials" / "research" / "m19_dc_maker.csv"
OUT = P / "out" / "validation"
FLUX_DENSITIES = (0.5, 1.0, 1.5, 1.8)                                   # T
POINTS = ((0, 0), (0, 10), (-6, 8), (6, 8))                              # check_femm's (id, iq) at rotor position 0


# ----------------------------------------------------------------------------------------------- the curves
def maker_curve() -> np.ndarray:
    """The maker's DC curve: rows (B [T], H [A/m])."""
    with open(MAKER, newline="", encoding="utf-8") as f:
        return np.array([[float(r["B_T"]), float(r["H_A_per_m"])] for r in csv.DictReader(f)])


def library_block(name: str = "M-19 14 mil") -> dict:
    """FEMM 5.0's library block of that name, as materials/common.py reads it from matlib.dat."""
    return next(b for b in materials_common.read_matlib() if b["BlockName"] == name)


def femm42_block() -> dict:
    """FEMM 4.2's "M-19 Steel", from the copy the model carries, in the same shape as library_block()."""
    return materials_common.m19_steel_block()


def h_at(bh: np.ndarray, b: float) -> float:
    """H [A/m] a B-H table needs for flux density b: PCHIP (smooth, monotone) through ln H, so it needs no tabulated point at b."""
    bb, hh = bh[:, 0], bh[:, 1]
    keep = hh > 0                                                       # ln 0 does not exist: leave out the (0, 0) point
    return float(np.exp(PchipInterpolator(bb[keep], np.log(hh[keep]))(b)))


def curve_table(blocks: dict, maker: np.ndarray) -> pd.DataFrame:
    """H at the chosen B for every curve, and each library curve's difference from the maker's."""
    rows = []
    for b in FLUX_DENSITIES:
        h = {name: h_at(block["bh"], b) for name, block in blocks.items()}
        h["maker"] = h_at(maker, b)
        rows.append({"B_T": b, **{f"H_{k}": v for k, v in h.items()},
                     **{f"pct_{k}": 100 * (v / h["maker"] - 1) for k, v in h.items() if k != "maker"},
                     "maker_is_a_table_point": bool(np.isclose(maker[:, 0], b).any())})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------- solving with another steel
def define_steel(block: dict, fill: float | None = None, bh: np.ndarray | None = None, name: str = "M-19 Steel",
                 sigma: float | None = None, d_lam: float | None = None):
    """Inside an open FEMM: define the material `name` from a library block, optionally with another fill, curve, sigma or d_lam.

    The same two kinds of call as motorgeo.define_m19_steel(): mi_addmaterial() makes the material with its constants (FEMM wants
    them in the order of its manual: name, mu_x, mu_y, H_c, J, conductivity [MS/m], lamination thickness [mm], hysteresis lag,
    fill, lamination type, lag_x, lag_y, strands, strand diameter), then mi_addbhpoint() adds the curve's points one by one.
    """
    import femm
    curve = block["bh"] if bh is None else bh
    femm.mi_addmaterial(name, block["Mu_x"], block["Mu_y"], block["H_c"], block["J_re"], block["Sigma"] if sigma is None else sigma,
                        block["d_lam"] if d_lam is None else d_lam, block["Phi_h"], block["LamFill"] if fill is None else fill,
                        block["LamType"], block["Phi_hx"], block["Phi_hy"], block["NStrands"], block["WireD"])
    for b, h in curve:
        femm.mi_addbhpoint(name, float(b), float(h))


def solve_points(label: str, steel: callable, fine: bool, path: Path) -> list[dict]:
    """Build the default model with `steel()` as the material "M-19 Steel", solve check_femm's four points, return the numbers."""
    import femm
    p = mg.Params()
    original = mg.define_m19_steel
    mg.define_m19_steel = steel                        # to_femm() calls define_m19_steel() by name: this swaps the steel, nothing else
    try:
        mg.to_femm(p, mg.build(p), str(path), fine=fine)
    finally:
        mg.define_m19_steel = original                 # always put the model's own steel back
    femm.mi_close()
    model = mg.FemmModel(p, str(path))
    rows = []
    for i_d, i_q in POINTS:
        r = model.solve(0.0, i_d, i_q)
        rows.append({"steel": label, "mesh": "fine" if fine else "coarse", "id": i_d, "iq": i_q,
                     **{k: r[k] for k in ("psi_d", "psi_q", "T_age", "T_wst")}})
    return rows


def run_solves(b42: dict, b50: dict) -> pd.DataFrame:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "m19").mkdir(exist_ok=True)
    variants = {
        "4.2 curve, fill 0.98 (the model's)": lambda: define_steel(b42),
        "14 mil curve, fill 0.955 (library as is)": lambda: define_steel(b50),
        "14 mil curve, fill 0.98": lambda: define_steel(b50, fill=b42["LamFill"]),
        "4.2 curve, fill 0.955": lambda: define_steel(b42, fill=b50["LamFill"]),
        "4.2 curve, 14 mil's sigma and d_lam": lambda: define_steel(b42, sigma=b50["Sigma"], d_lam=b50["d_lam"]),
    }
    rows = []
    femm_session.open_femm()                           # one FEMM of our own, on its own port, at below-normal priority
    try:
        for k, (label, steel) in enumerate(variants.items()):
            for fine in (False, True):
                rows += solve_points(label, steel, fine, OUT / "m19" / f"v{k}_{'fine' if fine else 'coarse'}.fem")
    finally:
        femm_session.close_femm()
    return pd.DataFrame(rows)


def print_solves(df: pd.DataFrame):
    base_label = "4.2 curve, fill 0.98 (the model's)"
    for mesh in ("coarse", "fine"):
        print(f"\ncheck_femm's four solves, {mesh} mesh: torque T_age [N.m] and flux linkage psi_d [mWb]; "
              f"in brackets the change from the model's own steel")
        print(f"{'steel':40s} " + " ".join(f"{'(' + str(a) + ',' + str(b) + ') T':>17s}" for a, b in POINTS) + f"   {'psi_d (0,0)':>17s}")
        base = df[(df.steel == base_label) & (df.mesh == mesh)].set_index(["id", "iq"])
        for label in df.steel.unique():
            d = df[(df.steel == label) & (df.mesh == mesh)].set_index(["id", "iq"])
            cells = []
            for pt in POINTS:
                t, t0 = d.loc[pt, "T_age"], base.loc[pt, "T_age"]
                cells.append(f"{t:8.4f} ({100 * (t / t0 - 1):+6.2f}%)" if label != base_label else f"{t:8.4f}          ")
            psi, psi0 = d.loc[(0, 0), "psi_d"] * 1e3, base.loc[(0, 0), "psi_d"] * 1e3
            tail = f"{psi:8.3f} ({100 * (psi / psi0 - 1):+6.2f}%)" if label != base_label else f"{psi:8.3f}          "
            print(f"{label:40s} " + " ".join(f"{c:>17s}" for c in cells) + f"   {tail:>17s}")


def make_figure(table: pd.DataFrame, blocks: dict, maker: np.ndarray):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (a, b) = plt.subplots(1, 2, figsize=(12, 4.6))
    styles = {"femm42": ("FEMM 4.2 \"M-19 Steel\" (47 points, fill 0.98)", "#2a78d6", "-"),
              "femm50": ("FEMM 5.0 \"M-19 14 mil\" (61 points, fill 0.955)", "#eb6834", "-")}
    for key, (label, colour, ls) in styles.items():
        bh = blocks[key]["bh"]
        a.plot(bh[1:, 1], bh[1:, 0], ls, color=colour, label=label, lw=1.6)
    a.plot(maker[:, 1], maker[:, 0], "o", color="#0b0b0b", ms=4, label="maker's DC curve (16 points)")
    a.set_xscale("log")
    a.set_xlabel("H [A/m]")
    a.set_ylabel("B [T]")
    a.set_ylim(0, 2.3)
    a.set_title("The three curves", loc="left", fontsize=10)
    a.legend(fontsize=8, frameon=False, loc="lower right")
    a.grid(alpha=0.3)
    grid = np.linspace(0.2, 2.0, 91)
    for key, (label, colour, ls) in styles.items():
        bh = blocks[key]["bh"]
        b.plot(grid, [100 * (h_at(bh, g) / h_at(maker, g) - 1) for g in grid], ls, color=colour, label=label.split(" (")[0], lw=1.6)
    b.axhline(0, color="#c3c2b7", lw=0.8)
    for x in FLUX_DENSITIES:
        b.axvline(x, color="#e1e0d9", lw=0.8)
    b.set_xlabel("flux density B [T]")
    b.set_ylabel("H of the curve against the maker's H at the same B [%]")
    b.set_title("H each library curve needs, against the maker's, at the same B", loc="left", fontsize=9.5)
    b.text(0.99, 0.02, "between the maker's 16 points its curve is interpolated,\nso a wiggle of a few % is the interpolation", transform=b.transAxes,
           ha="right", va="bottom", fontsize=7.5, color="#52514e")
    b.legend(fontsize=8, frameon=False)
    b.grid(alpha=0.3)
    fig.tight_layout()
    path = OUT / "m19_curves.png"
    fig.savefig(path, dpi=130)
    plt.close(fig)
    shutil.copyfile(path, HERE / path.name)                # the copy that is kept beside this script, so a re-run cannot leave it stale
    print(f"\nwrote {path} and a copy in {HERE}")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    maker = maker_curve()
    blocks = {"femm42": femm42_block(), "femm50": library_block()}
    table = curve_table(blocks, maker)
    print("H [A/m] a curve needs for a given flux density, and how far each library curve is from the maker's")
    print(f"{'B [T]':>6s} {'FEMM 4.2':>11s} {'FEMM 5.0':>11s} {'maker':>11s} {'4.2 vs maker':>14s} {'5.0 vs maker':>14s}")
    for _, r in table.iterrows():
        print(f"{r.B_T:6.1f} {r.H_femm42:11.2f} {r.H_femm50:11.2f} {r.H_maker:11.2f} {r.pct_femm42:+13.2f}% {r.pct_femm50:+13.2f}%"
              + ("" if r.maker_is_a_table_point else "   (the maker's table has no point here: interpolated)"))
    for key in ("femm42", "femm50"):
        b = blocks[key]
        print(f"{b['BlockName']:14s}: {len(b['bh'])} points, lamination fill {b['LamFill']}, thickness {b['d_lam']} mm, conductivity {b['Sigma']} MS/m")
    table.to_csv(OUT / "m19_curves.csv", index=False)
    df = run_solves(blocks["femm42"], blocks["femm50"])
    df.to_csv(OUT / "m19_femm5.csv", index=False)
    print_solves(df)
    make_figure(table, blocks, maker)
    return 0


if __name__ == "__main__":
    sys.exit(main())
