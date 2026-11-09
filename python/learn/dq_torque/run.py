"""Torque in d-q on our band-free rotor, magnets on against magnets off: the FEMM runs behind the project note
Documentation/Torque in d-q - field, current and saliency.

What it solves, on one FEMM (its own port, below-normal priority), fine mesh, under the machine's solver slot:
  part 1  at 10 A peak, the current angle gamma swept 0 to 90 deg in 7.5 deg steps, 6 rotor positions each, plus the
          no-load flux at the same 6 positions. Magnets on, then off.
  part 2  L_d and L_q at 5, 20, 30 and 40 A (10 A comes from part 1): gamma = 0 (all q) and 90 deg (all -d),
          6 positions each. Magnets on, then off; with the magnets off at 40 A, the whole angle sweep (the axis check).
  part 3  field snapshots at rotor angle 0, kept as .ans files for the figures, with |B| read by FEMM at the middle of
          pole 0's two outer bridges and its centre rib: magnets on at no load (the air-gap field sampled too), and
          10 A on -d with the magnets on and off; and 40 A on -d with the magnets off.
"Magnets off" is the same file with N45's coercivity set to 0; mu_r stays 1.05 (common.py says why).

Run, from the code folder:  uv run python learn/dq_torque/run.py          (about 3.5 minutes, 330 solves)
Then:                       uv run python learn/dq_torque/report.py       (the note's numbers)
                            uv run python learn/dq_torque/figures.py      (the note's three figures)
Writes out/learn/dq_torque/: solves.csv, probes.csv, gap.csv, the two .fem files and fields/*.ans.
"""
import math
import shutil
import time

import femm                                  # the binding: its calls reach the FEMM femm_session.open_femm() starts
import pandas as pd

import common as c
import femm_queue
import femm_session
import motorgeo as mg


def snapshot(model, theta_e, i_d, i_q, points, ans_path, gap_deg=()):
    """Solve one state and read FEMM's |B| at `points` before the solution closes; keep the .ans as `ans_path`.

    The same steps as FemmModel.solve_abc: turn the rotor side of the air-gap element, set the three phase currents by the
    inverse Park transform, solve, load the solution. Then mo_getb at each point, and the air-gap field if asked."""
    p = model.p
    femm.mi_modifyboundprop("AGE", 10, theta_e / p.pole_pairs + model.axes["A"] - 90 / p.pole_pairs)
    for ph, t in zip("ABC", [math.radians(theta_e + s) for s in (0, -120, 120)]):
        femm.mi_modifycircprop(ph, 1, i_d * math.cos(t) - i_q * math.sin(t))
    femm.mi_analyze(1)
    femm.mi_loadsolution()
    b = {k: math.hypot(*femm.mo_getb(x, y)[:2]) for k, (x, y) in points.items()}
    gap = [femm.mo_getgapb("AGE", a) for a in gap_deg]
    femm.mo_close()
    shutil.copyfile(model.path[:-4] + ".ans", ans_path)
    return b, gap


def solve_variant(name, hc):
    """Parts 1 to 3 for one variant. Returns (solve rows, probe rows, air-gap rows)."""
    p = c.params()
    regions = mg.build(p)
    fem = str(c.OUT / f"{name}.fem")                   # one file per variant: two models never share a file name
    mg.to_femm(p, regions, fem, fine=True)             # N45 written with H_c = 1,030,272 A/m; the document stays open
    if hc != c.HC_N45:
        femm.mi_modifymaterial("N45", 3, hc)           # property 3 is H_c [A/m]; properties 1 and 2 (mu_r) are untouched
        femm.mi_saveas(fem)
    femm.mi_close()
    model = mg.FemmModel(p, fem)
    row = lambda part, I, g, r: {"variant": name, "part": part, "I": I, "gamma": g, **r}

    rows = [row("noload", 0.0, float("nan"), model.solve(float(th), 0.0, 0.0)) for th in c.MEAN_POS]
    for g in c.GAMMAS:                                                                    # part 1
        i_d, i_q = c.currents(c.I_RATED, g)
        rows += [row("gamma", c.I_RATED, float(g), model.solve(float(th), i_d, i_q)) for th in c.MEAN_POS]
    for I in c.CURRENTS:                                                                  # part 2
        if I == c.I_RATED:
            continue                                                                      # part 1 has gamma 0 and 90 at 10 A
        # magnets off at 40 A: the whole angle sweep, for the axis check (is the torque zero at 0 and 90 deg?)
        angles = c.GAMMAS if (name == "off" and I == max(c.CURRENTS)) else (0.0, 90.0)
        for g in angles:
            i_d, i_q = c.currents(I, g)
            rows += [row("inductance", I, float(g), model.solve(float(th), i_d, i_q)) for th in c.MEAN_POS]

    points = c.probe_points(p, regions)                                                   # part 3
    states = [("no load", 0.0, 0.0, 0.0)] if name == "on" else []
    states += [("10 A on -d", 10.0, -10.0, 0.0)] + ([("40 A on -d", 40.0, -40.0, 0.0)] if name == "off" else [])
    probes, gap_rows = [], []
    for state, I, i_d, i_q in states:
        tag = state.replace(" ", "_").replace("-", "minus")
        want_gap = c.GAP_DEG if state == "no load" else ()
        b, gap = snapshot(model, 0.0, i_d, i_q, points, str(c.FIELDS / f"{name}_{tag}.ans"), want_gap)
        probes += [{"variant": name, "state": state, "I": I, "where": k, "B_T": v} for k, v in b.items()]
        gap_rows += [{"deg": d, "br": br, "bt": bt} for d, (br, bt) in zip(want_gap, gap)]
    femm.mi_close()
    print(f"  magnets {name}: {len(rows) + len(states)} solves", flush=True)
    return rows, probes, gap_rows


if __name__ == "__main__":
    c.FIELDS.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    rows, probes, gap_rows = [], [], []
    with femm_queue.solver_slot("learn/dq_torque: magnets on and off"):
        femm_session.open_femm()
        try:
            for name, hc in c.VARIANTS.items():
                r, pr, gp = solve_variant(name, hc)
                rows, probes, gap_rows = rows + r, probes + pr, gap_rows + gp
        finally:
            femm_session.close_femm()
    df = pd.DataFrame(rows).drop(columns=["gap_br", "gap_bt"])
    df.to_csv(c.OUT / "solves.csv", index=False)
    pd.DataFrame(probes).to_csv(c.OUT / "probes.csv", index=False)
    pd.DataFrame(gap_rows).to_csv(c.OUT / "gap.csv", index=False)
    # The magnets-off no-load solves have a flux linkage of exactly 0.0, which makes FemmModel retry each read 3 times
    # (its guard against a wrong 0.0, see its _flux_linkage docstring). Here 0.0 is the right answer, so those retries
    # are expected: 6 solves x 3 = 18. Any other retry would be news.
    unexpected = int(df[~((df.variant == "off") & (df.part == "noload"))].psi_retries.sum())
    print(f"{len(df)} solves in {time.perf_counter() - t0:.0f} s; unexpected flux read-back retries: {unexpected}")
