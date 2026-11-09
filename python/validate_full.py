"""Is the 90 deg motorgeo slice the same machine as the full 360 deg hand-drawn .FEM?

Both get the same test: a fixed stator field (i_a = I, i_b = i_c = -I/2) while the rotor turns through one
electrical period, at no load and at load. Compared: phase-A flux linkage and torque vs the electrical angle
between the rotor's north-pole d-axis and phase A's axis. The .FEM is copied first; the original is never touched.
Writes out/validate/validate.csv, and a .fem and .ans for every job: hand_000 ... hand_047 (the full model) and
slice_100 ... slice_107 (the 90 deg slice). Every job has files of ITS OWN, and that is not a detail: a FEMM saves its model
to its file before every solve and every solver then rewrites the .ans, so two FEMMs on one name read each other's
half-written files. They crash (fkn.exe, or FEMM itself), hang at 100 % CPU, or hand a job the other job's answer.
About 27 s on 8 workers with FEMM 5.0 (the FEMM 4.2 run took 29 min).
Run: uv run python validate_full.py   (from this folder; reads the reference .FEM that config.toml names, [paths] reference_fem)

Expect with FEMM 5.0 (2026-11-04), 96 solves; every run gives these to the last digit:
  hand-drawn 360  I= 0.0 A: psi_a -11.651 ..  11.655 mWb, T -0.0207 ..  0.0191 N.m, 58376 nodes
  hand-drawn 360  I=10.0 A: psi_a   0.270 ..  25.169 mWb, T -0.9502 ..  0.9523 N.m, 58376 nodes
  motorgeo 90     I= 0.0 A: psi_a -11.540 ..  11.539 mWb, T -0.0121 ..  0.0126 N.m
  motorgeo 90     I=10.0 A: psi_a   0.156 ..  24.993 mWb, T -0.9377 ..  0.9386 N.m
So the slice is within 1.13 % of the full model on flux linkage (the largest difference at any rotor position, as a share of
the full model's peak flux linkage, the worse of no load and 10 A) and within 1.43 % on the peak torque at 10 A (the slice's
peak torque against the full model's); FEMM 4.2's two numbers were 1.15 % and 1.48 %. (validation/compare.py's
slice_vs_full() computes them from validate.csv; the four lines above show only each column's smallest and largest value.)
FEMM 4.2 made 90,611 nodes for the full model (FEMM 5.0 makes fewer from the same mesh sizes) and gave, for it, psi_a
-11.623 .. 11.622 mWb at no load and 0.324 .. 25.186 mWb with T -0.9433 .. 0.9511 N.m at 10 A.
"""
import math
import os
import shutil
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

import femfile
import femm_session
import motorgeo as mg
import paths
from check_vs_fem import mmf_axes

FEM = str(paths.path("reference_fem"))         # the full 360 deg model: the copy in this folder's reference\ (config.toml)
OUT = os.path.abspath("out/validate")
WORKERS = 8                                    # 10 made the desktop stutter on 10-28; 8 leaves room for other work on the PC
I_LOAD = 10.0                                  # A peak in phase A (x 10 turns), -5 A in B and C
THETAS = np.arange(0, 360, 15.0)               # electrical degrees, rotor d-axis minus phase-A axis
P = mg.Params()
AIR_GAP_LABEL = (45.398271307113383, 0.39618471162615876)   # mm: the .FEM's air-gap label. It sits in the rotor's group 1, so
                                                            # hand_worker moves it to group 0 (validation/hand_reference_check.py shares it)


def fem_axes():
    """Rotor d-axis (pole 0, N) and phase-A axis of the .FEM, both in mechanical degrees."""
    fem = femfile.read_fem(FEM)
    pitch = 360 / P.slots
    name = {"Phase A": "A", "Phase B": "B", "Phase C": "C"}
    cond = [(round(math.degrees(math.atan2(lb["xy"][1], lb["xy"][0])) / pitch) * pitch, name[lb["circuit"]], lb["turns"])
            for lb in fem["labels"] if lb["circuit"]]
    axes, _ = mmf_axes(cond, P.pole_pairs)
    return 0.5, axes["A"]                       # the .FEM's pole 0 d-axis sits at 0.5 deg (rib corners)


def hand_worker(job):
    """Full model: rotate the rotor group by moving its geometry (FEMM re-meshes), then solve.

    The copy of the .FEM is changed in three ways before it is solved: the air-gap label goes to group 0, the rotor
    (group 1) is turned about the origin, and the three circuits get their currents. validation/hand_reference_check.py
    makes the same three changes to the reference and compares it with a copy made here, and
    tests/test_hand_reference_check.py fails when this function changes what it does to the copy: if it is changed,
    change apply_hand_worker there too.
    """
    import femm
    k, theta_e, i_a = job
    d0, a0 = fem_axes()
    move = (a0 + theta_e / P.pole_pairs) - d0          # mechanical rotation that puts the d-axis there
    path = os.path.join(OUT, f"hand_{k:03d}.fem")
    shutil.copyfile(FEM, path)
    femm_session.open_femm()                           # this worker's own FEMM: own port, below-normal priority
    try:
        femm.opendocument(path)
        femm.mi_selectlabel(*AIR_GAP_LABEL)                              # the air-gap label sits in the rotor
        femm.mi_setgroup(0)                                              # group: keep it out of the torque set
        femm.mi_clearselected()
        femm.mi_selectgroup(1)
        femm.mi_moverotate(0, 0, move)
        femm.mi_clearselected()
        for ph, i in (("Phase A", i_a), ("Phase B", -i_a / 2), ("Phase C", -i_a / 2)):
            femm.mi_modifycircprop(ph, 1, i)
        femm.mi_saveas(path)
        t = time.perf_counter()
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        psi_a = femm.mo_getcircuitproperties("Phase A")[2]
        femm.mo_groupselectblock(1)
        torque = femm.mo_blockintegral(22)
        n = femm.mo_numnodes()
        femm.mo_close()
        return {"model": "hand-drawn 360", "theta_e": theta_e, "I": i_a, "psi_a": psi_a, "T": torque,
                "nodes": n, "seconds": time.perf_counter() - t}
    finally:
        femm_session.close_femm()


def slice_worker(job):
    """motorgeo 90 deg slice, same stator field, rotor turned with the air-gap element."""
    import femm
    k, thetas, i_a = job
    path = os.path.join(OUT, f"slice_{k:03d}.fem")
    rows = []
    femm_session.open_femm()                           # this worker's own FEMM: own port, below-normal priority
    try:
        mg.to_femm(P, mg.build(P), path, fine=True)
        femm.mi_close()
        m = mg.FemmModel(P, path)
        for th in thetas:
            t = time.perf_counter()
            r = m.solve_abc(th, i_a, -i_a / 2, -i_a / 2)
            rows.append({"model": "motorgeo 90", "theta_e": th, "I": i_a, "psi_a": r["psi_a"], "T": r["T_age"],
                         "nodes": None, "seconds": time.perf_counter() - t})
    finally:
        femm_session.close_femm()
    return rows


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    t = time.perf_counter()
    hand_jobs = [(k, th, i) for k, (th, i) in enumerate((th, i) for i in (0.0, I_LOAD) for th in THETAS)]
    # One slice job = one set of angles at one current, and its id names its files: slice_100 ... slice_107. Until
    # 2026-11-04 the id was 100 + k, the same for the no-load job and the 10 A job of the same angles, so two FEMMs
    # saved, solved and read one slice_10k.fem / .ans at the same time (every mi_analyze saves the model to its file
    # first, and every solver then rewrites the .ans). Run that way the script crashed fkn.exe or FEMM itself, left a
    # FEMM spinning at 100 % CPU, or, when it did not fail, handed a job the other job's answer: a flux linkage of
    # exactly 0 with the right torque (the "0.0 glitch" of 2026-10-28), or a 10 A value in a no-load row. Two FEMMs
    # must never share a file name: learn/04_workers.py, part 3, shows it.
    slice_jobs = [(100 + n, list(THETAS[k::4]), i)
                  for n, (i, k) in enumerate((i, k) for i in (0.0, I_LOAD) for k in range(4))]
    for jobs in (hand_jobs, slice_jobs):
        assert len({job[0] for job in jobs}) == len(jobs), "two jobs would share one file name"
    with ProcessPoolExecutor(WORKERS) as pool:
        f_hand = pool.map(hand_worker, hand_jobs)
        f_slice = pool.map(slice_worker, slice_jobs)
        rows = list(f_hand) + [r for part in f_slice for r in part]
    df = pd.DataFrame(rows).sort_values(["model", "I", "theta_e"])
    df.to_csv(os.path.join(OUT, "validate.csv"), index=False)
    print(f"{len(df)} solves in {time.perf_counter() - t:.0f} s")
    for (model, i), g in df.groupby(["model", "I"]):
        print(f"  {model:15s} I={i:4.1f} A: psi_a {g.psi_a.min() * 1e3:7.3f} .. {g.psi_a.max() * 1e3:7.3f} mWb, "
              f"T {g['T'].min():7.4f} .. {g['T'].max():7.4f} N.m, {g.seconds.mean():.1f} s/solve"
              + (f", {int(g.nodes.mean())} nodes" if g.nodes.notna().any() else ""))
