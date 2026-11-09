"""optimization.ac_machine: the machine's own field at the copper, over its operating range: the rotor's and the armature's field harmonics seen by a strand, from FEMM
magnetostatic solves of the 90 degree sector, and the proximity loss they make as a fast function of the strand, the operating point and the speed.

What it is, and why
    The first winding study measured the field the copper feels at ONE point, the rated current at its own angle, and took 3,000 rpm. The copper loss model needs it wherever the
    machine runs: from a fraction of the rated current to three times it, on the MTPA path and in field weakening (more negative i_d, the way the voltage limit is met above the corner
    speed), and from 3,000 rpm to the 16,330 rpm the stress check allows (1,090 Hz). This module solves the sector at those points, takes the field every point of the coil
    area sees over one electrical period (the rotor turns, the currents turn with it), and keeps what the loss needs:
      G_k      the mean square of dB/d(theta_e) at harmonic k of the electrical angle, over the coil area [T^2 per rad^2]: the first study's G, split by harmonic and by where in the slot;
      H_k      the same for the vector potential A: the spatial variance of A's harmonic within a coil layer (what drives circulating currents between strands of one turn);
      flux linkages psi_d, psi_q and the mean torque at every point (for the voltage the winding needs, and so the modulation index of the inverter).
    The field of a coil depends only on its AMPERE-TURNS (the first study showed that to the last digit), so the points are in amperes at the model's 10 turns a coil side and apply to every
    winding at the same ampere-turns; a winding with N turns and a paths at a peak phase current I has N I / a ampere-turns a coil side, i.e. N I / (10 a) amperes of the 10-turn model.

The method (each step from the one before)
    1. One solve = one rotor position theta_e (electrical degrees) with the d-q currents (i_d, i_q) turned into phase currents (motorgeo.FemmModel.solve): the fine mesh (smart
       mesh on, about 15.7 thousand nodes, 30 thousand triangles of the sector, 6,176 of them in the copper), FEMM 5.0, the band-free optimized rotor. The air-gap element turns the
       rotor without re-meshing, so a triangle is the same triangle at every position (checked, as in optimization/winding_femm.py).
    2. 48 positions, 7.5 degrees apart, cover ONE electrical period (the 90 degree sector holds exactly one: 4 pole pairs); B in a coil triangle is read from the nodal vector potential of the
       .ans (first-order triangles: B = (dA/dy, -dA/dx), constant in a triangle) and so is A (the mean of the three nodes).
    3. Each triangle's B(theta) and A(theta) are taken apart in their Fourier series over theta, c_k = (1/48) sum_j f(theta_j) exp(-i k theta_j), k = 1 ... 24. A sinusoid of peak B_k at
       harmonic k has |c_k| = B_k / 2. Parseval turns the mean square of the derivative into a sum:
           G = sum over triangles (area weighted) of sum_{k>=1} 2 k^2 (|cBx_k|^2 + |cBy_k|^2),     G_k the k-th term of that sum,
       exactly winding_femm.field_statistics (and checked against it for the first study's own point: the same G to the last digit).
    4. The proximity loss of a strand at that harmonic. A strand of diameter d in a field B_k cos(k w_e t) loses (low frequency) pi sigma d^4 (k w_e)^2 B_k^2 L / 128, i.e.
       pi sigma d^4 w_e^2 G_k L / 64 with G_k = k^2 B_k^2 / 2 (the mean over the strand's position); at higher frequency the exact lone-wire factor of ac_theory.proximity_ratio (d / delta at
       k f_e) multiplies it. For the winding: strands uniformly spread over the coil area (as the first study), 72 N n_s of them for 36 slots with 2 layers,
           P_prox = 72 N n_s L (pi sigma d^4 w_e^2 / 64) sum_k G_k(i_d, i_q) c_k(d, k f_e),   c_k = proximity_ratio * (the slot model's discreteness factor).
       It goes as the speed squared, and, for a given copper area, as d^2.
    5. G_k depends on the operating point. If the field were linear in the currents (the rotor's flux plus i_d times the d-axis armature field plus i_q times the q-axis one, per harmonic
       and per place) then G_k would be a quadratic polynomial in (i_d, i_q): 6 coefficients a harmonic (fit_quadratic). It is not: the quadratic form has to be 0 + a no-load field and then
       grow as the square of the ampere-turns, and with the iron's saturation it fails on the points it was not made from (leave-one-out errors of 20 ... 1700 % for the total G). What works
       is fit_split: V_g = V_g(no load) + AT^2 g_tot(AT, gamma) a_g(AT, gamma) per harmonic group g (the armature's part goes as AT^2 times a slowly varying total, shared among the harmonic
       groups in proportions that vary linearly), whose leave-one-out error for the total G is 5.7 % in the median of the nine loaded points and -19.5 % / +24.4 % at worst (field weakening at
       70 degrees, where the samples are fewest). The fit's limits are the region sampled (the 'extrapolated' flag of copper_loss.losses). The proximity loss that comes out is an UPPER BOUND for one reason this
       module does not model: the field is taken from magnetostatic solves, and at the harmonics of the slot passing (the 9th, 11th, ... of the electrical angle, 10 kHz and up at the top speed) the
       laminations' eddy currents would screen it (ac_slot_study.py: 4 % at 5 kHz, 11 % at 20 kHz for the slot's own field).

The operating points (the amperes are peak, at 10 turns a coil side; the angle gamma is from the q axis toward -d, i_d = -I sin(gamma), i_q = I cos(gamma))
    no load                      0 A  (the rotor's field alone: the magnets' slot harmonics, nothing from the armature)
    the first study's point      10.00 A at 19.657 deg (100 ampere-turns; 'reference': the check against the first study)
    MTPA at 0.5, 1, 2, 3 x rated 5.93, 11.86, 23.72, 35.58 A at the angle of maximum torque of each current (found here by a parabola through the best three of six angles: 10 ... 60 deg)
    field weakening              11.86 A and 23.72 A at 45 and 70 deg (more negative i_d than MTPA at that current)
    The rated point is the brief's 11.86 A (118.6 ampere-turns) at the MTPA angle. The brief gives 22.47 deg for the optimized design; the parabola through the six angles gives 23.02 deg here
    (the torque is flat there: by the parabola the two angles differ by about 0.01 % in torque), and the points use 23.02 deg.

Run (from the Python to FEMM folder; the solver slot is taken in two sessions, the angle search of a few seconds and the fields, 80 tasks, of about 1 minute (52 s of FEMM), and waited for if
another run holds it)
    uv run python -m optimization.ac_machine                    all stages (angles, fields, analysis, figures)
    uv run python -m optimization.ac_machine --stage analyse    from the saved raw results only (no FEMM)
Writes out/optimization/copper_ac/machine/: raw/<tag>.pkl (one file a pool task: resumable), points.json (the operating points with their angles), summary.json (every number), tables.json
(the tables copper_loss.py packages), fig_machine_field.png.

Numbers to expect (this PC, FEMM 5.0.20261101, 15.7 thousand nodes, 6,176 coil triangles; the printed summary has them all)
    the check against the first study's own point (100 ampere-turns at 19.657 deg): G = 4.8222292e-4 T^2/rad^2 here and there, relative difference 2e-15, the harmonics' 5e-15;
    MTPA angles 12.33 deg at 0.5 x rated, 23.02 at 1 x, 37.28 at 2 x and 43.37 at 3 x; the mean torque at them 0.878, 1.807, 3.645 and 5.085 N.m;
    G over the coil area (1e-4 T^2/rad^2): no load 0.522, MTPA 0.5 / 1 / 2 / 3 x rated 1.92 / 6.86 / 30.6 / 87.0, field weakening at 1 x rated 6.01 (45 deg) and 5.92 (70 deg), at 2 x rated
    29.7 and 24.2; at the rated MTPA point the top layer (next to the opening) has 12.3 against 1.39 for the bottom layer, and the fundamental carries 69 % of G;
    the no-load field is the rotor's slot harmonics: its G is 0.52e-4, 1 / 13 of the rated point's, and the fundamental is 15 % of it;
    the fits: total G leave-one-out -19.5 % ... +24.4 % (median 5.7 %), the vector potential's spread H up to 13.7 %, the flux linkage psi_d 1.2 mWb of 27.7 and psi_q 3.4 mWb of 38.0 (quadratic).

Self-check (the answers follow).
  1. The loss goes as sum_k G_k times k-dependent factors, not as G times one factor. Why does the harmonic split matter at 16,000 rpm?
  2. Why is the no-load field not zero in the copper, and why does it matter for the proximity loss?
  3. Would the points of this module apply to a winding with 30 % of the ampere-turns on a coil side?
Answers.
  1. Each harmonic k is at frequency k f_e, and the lone-wire factor (the strand shielding itself) depends on d / delta at THAT frequency: at 1,090 Hz the 9th harmonic is at 9.8 kHz and
     the 17th at 18.5 kHz, where a 0.63 mm strand has d / delta 0.82 and 1.13 and the factor is 0.997 and 0.988 (a 1.4 mm wire: 1.83 and 2.52, 0.93 and 0.78); and the circulating-current
     and ripple losses care about the frequency too. At low d / delta the factor is 1 and the split does not matter: P goes as the speed squared times the total G.
  2. The magnets' flux is concentrated in the teeth, but a part of it passes the slot opening and enters the top of the coil: the slot-passing harmonic (the 9th of the electrical
     angle: 36 slots over 4 pole pairs) is the rotor's alone. It is there at no load and at any speed, it needs no current, and it makes the strands next to the opening lose power.
  3. Yes: the field is a function of the ampere-turns only, and 30 % of the rated 118.6 is 35.6 ampere-turns, inside the sampled range (0 ... 356). Outside it the fit is not trusted.
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
import time
from pathlib import Path

import numpy as np

import femfile
import femm_queue
import motorgeo as mg
import paths
from optimization import ac_theory, design, winding_femm

I_REF_TURNS = 10                                 # turns a coil side of the FEMM model
AT_RATED = 118.6                                 # ampere-turns a coil side at the rated point (the brief's 11.86 A at 10 turns)
GAMMA_FIRST_STUDY = 19.656668507770743           # deg: the first study's own MTPA angle for the 100 ampere-turns point of the optimized design (winding/optimized/femm.json)
N_POS = 48                                       # rotor positions over one electrical period
POSITIONS = tuple(360.0 * j / N_POS for j in range(N_POS))
POS_PER_TASK = 6
ANGLE_GRID = (10.0, 20.0, 30.0, 40.0, 50.0, 60.0)
ANGLE_POSITIONS = (0.0, 30.0)                    # as design.ANGLE_POSITIONS: their mean cancels the 6th, 18th and 30th harmonics of the torque ripple
OUT = paths.path("out") / "optimization" / "copper_ac" / "machine"
N_DEPTH_BINS = 14
TASK_TIMEOUT_S = 300


def overrides() -> dict:
    """The optimized design's geometry overrides (out/optimization/final/final.json)."""
    return json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))["overrides"]


# ---------------------------------------------------------------------------------------------- the pool's tasks (they run in a worker, on its own FEMM)
def _model_files(task: dict) -> tuple[str, str]:
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
    return fem, fem[:-4] + ".ans"


def angle_task(task: dict) -> dict:
    """One pool task: the mean torque (at the two positions that cancel the ripple's low harmonics) for a list of current angles at one current amplitude.
    task = {'tag', 'overrides', 'fem_dir', 'amps' (A peak at 10 turns), 'gammas' (deg)}."""
    import femm
    import traceback
    started = time.perf_counter()
    fem, ans = _model_files(task)
    try:
        p = design.params(task["overrides"])
        mg.to_femm(p, mg.build(p), fem, fine=True)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        torque = []
        for g in task["gammas"]:
            i_d, i_q = design.currents(task["amps"], g)
            torque.append(float(np.mean([model.solve(th, i_d, i_q)["T_age"] for th in ANGLE_POSITIONS])))
    except Exception as error:
        return {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4), "seconds": time.perf_counter() - started}
    finally:
        _cleanup(femm, fem, ans)
    return {"status": "ok", "kind": "angle", "tag": task["tag"], "amps": task["amps"], "gammas": list(task["gammas"]), "torque": torque, "seconds": time.perf_counter() - started}


def field_task(task: dict) -> dict:
    """One pool task: B and A in every coil triangle of the sector at a list of rotor positions at one operating point (i_d, i_q).
    task = {'tag', 'overrides', 'fem_dir', 'i_d', 'i_q', 'positions' (electrical degrees)}. Returns arrays bx, by, a of shape (positions, coil triangles) [T, T, Wb/m] and the
    triangles' area [mm^2], centroid [mm], layer index (0 ... 17: coil_<slot>_<layer> in the order of motorgeo's regions), and for every position psi_d, psi_q, the air-gap torque.
    As winding_femm.rated_field, which it follows: the mesh must be the same at every position (checked); the triangles of a layer are found from the block labels (matched to the
    regions by the label POINT, not by position in the list)."""
    import femm
    import traceback
    started = time.perf_counter()
    fem, ans = _model_files(task)
    try:
        p = design.params(task["overrides"])
        regions = mg.build(p)
        mg.to_femm(p, regions, fem, fine=True)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        bx_all, by_all, a_all, psi_d, psi_q, torque, static = [], [], [], [], [], [], None
        for theta in task["positions"]:
            row = model.solve(theta, task["i_d"], task["i_q"])
            sol = femfile.read_ans(ans)
            bx, by = femfile.element_b(sol["xy"], sol["A"], sol["tri"])
            a_tri = sol["A"][sol["tri"]].mean(axis=1)
            if static is None:
                static = winding_femm.coil_triangles(sol, regions)
                static["tri"] = sol["tri"][static["mask"]]
                static["xy"] = sol["xy"]
                static["n_triangles"] = int(len(sol["tri"]))
                if static["outside"]:
                    raise RuntimeError(f"{static['outside']} coil triangles lie outside their layer's polygon: the label matching is wrong")
            elif not (np.array_equal(sol["tri"][static["mask"]], static["tri"]) and np.array_equal(sol["xy"], static["xy"])):
                raise RuntimeError(f"the mesh at theta_e = {theta} differs from the first position's: the elements are not comparable")
            m = static["mask"]
            bx_all.append(bx[m])
            by_all.append(by[m])
            a_all.append(a_tri[m])
            psi_d.append(row["psi_d"])
            psi_q.append(row["psi_q"])
            torque.append(row["T_age"])
        nodes = int(len(static["xy"]))
    except Exception as error:
        return {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4), "seconds": time.perf_counter() - started}
    finally:
        _cleanup(femm, fem, ans)
    return {"status": "ok", "kind": "field", "tag": task["tag"], "positions": list(task["positions"]), "bx": np.array(bx_all), "by": np.array(by_all), "a": np.array(a_all),
            "area_mm2": static["area_mm2"], "centroid": static["centroid"], "layer": static["layer"], "names": static["names"], "psi_d": psi_d, "psi_q": psi_q, "T_age": torque,
            "nodes": nodes, "triangles": static["n_triangles"], "coil_triangles": int(len(static["area_mm2"])), "seconds": time.perf_counter() - started}


def run_task(task: dict) -> dict:
    """The one function the pool maps: an 'angle' task or a 'field' task."""
    if task["kind"] == "angle":
        return angle_task(task)
    if task["kind"] == "field":
        return field_task(task)
    raise ValueError(f"unknown kind of task: {task['kind']!r}")


def _cleanup(femm, fem: str, ans: str) -> None:
    import os
    try:
        femm.mi_close()
    except Exception:
        pass
    for path in (ans, fem):
        try:
            os.remove(path)
        except OSError:
            pass


# ---------------------------------------------------------------------------------------------- the operating points
def mtpa_angle(gammas, torque) -> float:
    """The current angle of maximum mean torque [deg]: a parabola through the best of the sampled angles and its two neighbours (the angles of the grid are 10 deg apart, and the torque
    against the angle is not a parabola over the whole range, only near its top)."""
    g, t = np.asarray(gammas, float), np.asarray(torque, float)
    k = int(np.argmax(t))
    k = min(max(k, 1), len(g) - 2)
    a2, a1, _ = np.polyfit(g[k - 1:k + 2], t[k - 1:k + 2], 2)
    if a2 >= 0:
        return float(g[int(np.argmax(t))])
    return float(np.clip(-a1 / (2.0 * a2), g[k - 1], g[k + 1]))


def operating_points(mtpa: dict[float, float]) -> list[dict]:
    """The list of operating points: {'name', 'amps' (A peak at 10 turns), 'gamma' (deg), 'i_d', 'i_q', 'at' (ampere-turns a coil side), 'group'}; mtpa maps each current to its MTPA angle."""
    rated = AT_RATED / I_REF_TURNS
    pts = [{"name": "no load", "amps": 0.0, "gamma": 0.0, "group": "no load"},
           {"name": "first study (100 AT)", "amps": 10.0, "gamma": GAMMA_FIRST_STUDY, "group": "reference"}]
    for fac in (0.5, 1.0, 2.0, 3.0):
        amps = round(rated * fac, 6)
        pts.append({"name": f"MTPA {fac:g} x rated", "amps": amps, "gamma": mtpa[amps], "group": "mtpa"})
    for fac in (1.0, 2.0):
        for g in (45.0, 70.0):
            pts.append({"name": f"FW {fac:g} x rated, {g:g} deg", "amps": round(rated * fac, 6), "gamma": g, "group": "fw"})
    for q in pts:
        q["i_d"], q["i_q"] = design.currents(q["amps"], q["gamma"])
        q["at"] = q["amps"] * I_REF_TURNS
    return pts


# ---------------------------------------------------------------------------------------------- the stage that needs FEMM
def run_angles(workers: int = 8, log=print) -> dict:
    """One pool session: the MTPA angle at 0.5, 1, 2 and 3 x the rated ampere-turns. Returns {amps: gamma}, saved in OUT/mtpa_angles.json."""
    OUT.mkdir(parents=True, exist_ok=True)
    rated = AT_RATED / I_REF_TURNS
    amps_list = [round(rated * f, 6) for f in (0.5, 1.0, 2.0, 3.0)]
    tasks = [{"kind": "angle", "tag": f"cu_ang_{a:g}".replace(".", "p"), "overrides": overrides(), "fem_dir": str(OUT / "fem"), "amps": a, "gammas": list(ANGLE_GRID)} for a in amps_list]
    with femm_queue.femm_pool("copper ac: MTPA angles", workers=min(workers, len(tasks)), task_timeout_s=TASK_TIMEOUT_S) as pool:
        res = pool.map(run_task, tasks)
    bad = [(t["tag"], r) for t, r in zip(tasks, res) if not (isinstance(r, dict) and r.get("status") == "ok")]
    if bad:
        raise RuntimeError("angle tasks failed: " + "; ".join(f"{tag}: {r.get('error', r) if isinstance(r, dict) else r}" for tag, r in bad))
    out = {}
    for r in res:
        gam = mtpa_angle(r["gammas"], r["torque"])
        out[r["amps"]] = gam
        log(f"MTPA at {r['amps']:.2f} A ({r['amps'] * I_REF_TURNS:.1f} ampere-turns): {gam:.2f} deg   (torque at 10 ... 60 deg: {', '.join(f'{t:.3f}' for t in r['torque'])})")
    (OUT / "mtpa_angles.json").write_text(json.dumps({"angles": {str(k): v for k, v in out.items()}, "raw": [{k: v for k, v in r.items()} for r in res]}, indent=1), encoding="utf-8")
    return out


def run_fields(points: list[dict], workers: int = 8, log=print) -> None:
    """One pool session: the 48 rotor positions of every operating point (POS_PER_TASK a task); results saved one file a task (OUT/raw/<tag>.pkl), tasks already saved are skipped (so a
    stopped run resumes)."""
    raw = OUT / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    tasks = []
    for k, q in enumerate(points):
        for w in range(N_POS // POS_PER_TASK):
            tag = f"cu_fld_{k:02d}_{w}"
            if (raw / f"{tag}.pkl").exists():
                continue
            tasks.append({"kind": "field", "tag": tag, "overrides": overrides(), "fem_dir": str(OUT / "fem"), "i_d": q["i_d"], "i_q": q["i_q"],
                          "positions": list(POSITIONS[w * POS_PER_TASK:(w + 1) * POS_PER_TASK])})
    if not tasks:
        log("all field tasks are saved already")
        return
    t0 = time.perf_counter()
    with femm_queue.femm_pool("copper ac: machine field", workers=workers, task_timeout_s=TASK_TIMEOUT_S) as pool:
        log(f"{len(tasks)} field tasks to run ({time.perf_counter() - t0:.1f} s waiting for the slot)")
        t1 = time.perf_counter()
        res = pool.map(run_task, tasks)
    bad = [(t["tag"], r) for t, r in zip(tasks, res) if not (isinstance(r, dict) and r.get("status") == "ok")]
    for t, r in zip(tasks, res):
        if isinstance(r, dict) and r.get("status") == "ok":
            with open(raw / f"{t['tag']}.pkl", "wb") as handle:
                pickle.dump(r, handle)
    log(f"{len(tasks) - len(bad)} of {len(tasks)} field tasks done in {time.perf_counter() - t1:.1f} s of FEMM")
    if bad:
        raise RuntimeError("field tasks failed: " + "; ".join(f"{tag}: {r.get('error', r) if isinstance(r, dict) else r}" for tag, r in bad))


def load_point(k: int) -> dict:
    """The raw fields of operating point k, the tasks' pieces put together in position order: arrays (48, coil triangles)."""
    pieces = []
    for w in range(N_POS // POS_PER_TASK):
        with open(OUT / "raw" / f"cu_fld_{k:02d}_{w}.pkl", "rb") as handle:
            pieces.append(pickle.load(handle))
    first = pieces[0]
    for r in pieces[1:]:
        for name in ("area_mm2", "centroid", "layer"):
            if not np.array_equal(first[name], r[name]):
                raise RuntimeError(f"point {k}: the field tasks have different meshes ({name} differs)")
    cat = lambda name: np.concatenate([r[name] for r in pieces], axis=0)
    return {"bx": cat("bx"), "by": cat("by"), "a": cat("a"), "area_mm2": first["area_mm2"], "centroid": first["centroid"], "layer": first["layer"], "names": first["names"],
            "psi_d": np.concatenate([r["psi_d"] for r in pieces]), "psi_q": np.concatenate([r["psi_q"] for r in pieces]), "T_age": np.concatenate([r["T_age"] for r in pieces]),
            "positions": np.concatenate([r["positions"] for r in pieces]), "nodes": first["nodes"], "coil_triangles": first["coil_triangles"]}


# ---------------------------------------------------------------------------------------------- the analysis (no FEMM)
def depth_of(centroid: np.ndarray) -> np.ndarray:
    """The radius [mm] of each triangle's centroid: the depth coordinate of the slot (the opening is at the bore, 45.8 mm; the slot's coil region runs from 46.7 to about 59.8 mm)."""
    return np.hypot(centroid[:, 0], centroid[:, 1])


def harmonic_stats(raw: dict, n_depth: int = N_DEPTH_BINS) -> dict:
    """The field statistics of one operating point from its raw fields: G_k and H_k (module docstring, step 3) over the whole coil area, per coil layer (0 = next to the opening,
    1 = the bottom: motorgeo's names coil_<slot>_<layer>), and per depth bin; B^2_k the mean square of the harmonic's flux density; and the mean flux linkages and torque.
    The FFT is over the 48 positions (axis 0); k runs from 1 to 24 (index k: the Nyquist term k = 24 counts with weight 1/2, as winding_femm.field_statistics)."""
    bx, by, a = raw["bx"], raw["by"], raw["a"]
    n = bx.shape[0]
    cbx, cby, ca = (np.fft.rfft(v, axis=0) / n for v in (bx, by, a))               # (25, triangles)
    k = np.arange(n // 2 + 1, dtype=float)
    w = 2.0 * k ** 2
    w[0] = 0.0
    w[-1] = 0.5 * k[-1] ** 2
    area = raw["area_mm2"]
    layer_idx = raw["layer"]                                                         # 0 ... 17: slot * 2 + layer
    layer_in_slot = layer_idx % 2
    pg = w[:, None] * (np.abs(cbx) ** 2 + np.abs(cby) ** 2)                          # per triangle and harmonic: the contribution to <(dB/dtheta)^2>
    pb = 2.0 * (np.abs(cbx) ** 2 + np.abs(cby) ** 2)
    pb[0] = 0.0
    pa_dev = np.zeros_like(pg)                                                       # filled below: the spatial variance of A's harmonic (per layer of each slot)
    r_tri = depth_of(raw["centroid"])
    edges = np.linspace(r_tri.min(), r_tri.max() + 1e-9, n_depth + 1)
    bins = np.clip(np.digitize(r_tri, edges) - 1, 0, n_depth - 1)

    def wmean(x, sel):
        return (x[:, sel] * area[sel]).sum(axis=1) / area[sel].sum()

    sel_all = np.ones(len(area), bool)
    out = {"n_positions": int(n), "n_triangles": int(len(area)), "G_k": wmean(pg, sel_all).tolist(), "B2_k": wmean(pb, sel_all).tolist(), "G": float(wmean(pg, sel_all).sum()),
           "G_k_layer": [wmean(pg, layer_in_slot == j).tolist() for j in range(2)], "B2_k_layer": [wmean(pb, layer_in_slot == j).tolist() for j in range(2)],
           "G_layer_of_slot": [float(wmean(pg, layer_idx == j).sum()) for j in range(int(layer_idx.max()) + 1)]}
    out["G_depth"] = [float(wmean(pg, bins == b).sum()) if (bins == b).any() else 0.0 for b in range(n_depth)]
    out["G_k_depth"] = [wmean(pg, bins == b).tolist() if (bins == b).any() else [0.0] * (n // 2 + 1) for b in range(n_depth)]
    out["depth_mm"] = [float(0.5 * (edges[b] + edges[b + 1])) for b in range(n_depth)]
    out["depth_area_mm2"] = [float(area[bins == b].sum()) for b in range(n_depth)]
    # the spatial variance of A's complex harmonic inside each layer of each slot: sum over triangles of area x |c_k - mean|^2, divided by the layer's area (then 4 for the peak amplitude)
    hk = np.zeros((int(layer_idx.max()) + 1, n // 2 + 1))
    for j in range(int(layer_idx.max()) + 1):
        sel = layer_idx == j
        wsum = area[sel].sum()
        mean = (ca[:, sel] * area[sel]).sum(axis=1) / wsum
        hk[j] = (np.abs(ca[:, sel] - mean[:, None]) ** 2 * area[sel]).sum(axis=1) / wsum
    out["H_k_layer_of_slot"] = (w[None, :] * hk).tolist()                           # the mean square of d(A - mean A)/dtheta per harmonic, per layer of each slot: 2 k^2 var(c_k)
    out["H_k"] = (w[None, :] * hk).mean(axis=0).tolist()
    out["psi_d"], out["psi_q"], out["T_mean"] = float(raw["psi_d"].mean()), float(raw["psi_q"].mean()), float(raw["T_age"].mean())
    out["T_ripple_pp"] = float(raw["T_age"].max() - raw["T_age"].min())
    return out


def fit_quadratic(points: list[dict], values: np.ndarray, names: list[str] | None = None) -> dict:
    """Least-squares fit of values (n points, m outputs) to a quadratic polynomial in (i_d, i_q): 1, i_d, i_q, i_d^2, i_d i_q, i_q^2. Returns {'coef' (6, m), 'rms' (m,) of the residuals,
    'loo_max_rel' (m,): the largest relative error of any point predicted by the fit of the others (leave one out, over the points whose value is larger than 1e-3 of the maximum)}."""
    x = np.array([[1.0, q["i_d"], q["i_q"], q["i_d"] ** 2, q["i_d"] * q["i_q"], q["i_q"] ** 2] for q in points])
    y = np.asarray(values, float)
    coef, *_ = np.linalg.lstsq(x, y, rcond=None)
    resid = x @ coef - y
    loo = np.zeros_like(y)
    for i in range(len(points)):
        keep = np.arange(len(points)) != i
        c_i, *_ = np.linalg.lstsq(x[keep], y[keep], rcond=None)
        loo[i] = x[i] @ c_i - y[i]
    scale = np.maximum(np.abs(y), 1e-3 * np.abs(y).max(axis=0, keepdims=True))
    return {"coef": coef, "rms": np.sqrt((resid ** 2).mean(axis=0)), "loo_max_rel": np.abs(loo / scale).max(axis=0), "rel_loo": np.abs(loo / scale), "loo": loo}


GROUPS = tuple([k] for k in (1, 3, 5, 7, 9, 11)) + (list(range(13, 25)),)     # the harmonic groups of the model: the odd harmonics up to 11 one by one, then 13 ... 24 together
# (the field's harmonics are all odd: the even ones are zero to 1e-10, the winding and the rotor being symmetric over a pole pitch)


def group_sums(v: np.ndarray) -> np.ndarray:
    """v[..., k] for k = 0 ... 24 summed over each harmonic group: (..., 7)."""
    return np.stack([np.asarray(v)[..., g].sum(axis=-1) for g in GROUPS], axis=-1)


def group_k_eff(v: np.ndarray) -> np.ndarray:
    """The representative harmonic of each group: the mean of its k weighted by the statistic v (a group of one is its k): the frequency the group's loss is evaluated at."""
    v = np.asarray(v)
    out = []
    for g in GROUPS:
        w = np.maximum(v[..., g].mean(axis=tuple(range(v.ndim - 1))) if v.ndim > 1 else v[g], 0.0)
        out.append(float(np.sum(np.array(g) * w) / np.sum(w)) if np.sum(w) > 0 else float(np.mean(g)))
    return np.array(out)


def ghat_basis(at, gamma_deg) -> np.ndarray:
    """The basis of the field's total fit: (1, AT, gamma, AT gamma) with AT in ampere-turns a coil side and gamma in radians, as the last axis."""
    at, g = np.asarray(at, float), np.radians(np.asarray(gamma_deg, float))
    return np.stack([np.ones_like(at), at, g, at * g], axis=-1)


def share_basis(at, gamma_deg) -> np.ndarray:
    """The basis of the harmonic shares' fit: (1, AT, gamma) with AT in ampere-turns and gamma in radians."""
    at, g = np.asarray(at, float), np.radians(np.asarray(gamma_deg, float))
    return np.stack([np.ones_like(at), at, g], axis=-1)


def fit_split(points: list[dict], values: np.ndarray) -> dict:
    """Fit of a field statistic V_k(AT, gamma) (G_k or H_k; values (n points, 25 harmonics k = 0 ... 24)) as
        V_g = V_g(no load) + AT^2 g_tot(AT, gamma) a_g(AT, gamma)         for each harmonic group g (GROUPS),
        g_tot = c0 + c1 AT + c2 gamma + c3 AT gamma        (bilinear: the armature part's total per AT^2),
        a_g   = b0 + b1 AT + b2 gamma                      (linear: the group's share of it, normalised to 1 when used).
    Why this form and not a polynomial in (i_d, i_q): the armature's field is proportional to the current, so what it adds to a mean square goes as AT^2 and a polynomial in the currents would
    have to represent a function that is 0 at no load and then grows as AT^2 (faster with the iron's saturation); AT^2 g_tot has a g_tot that is nearly constant over the range
    (4.0 ... 6.8e-8) and four coefficients describe it (the leave-one-out error of the total G is 7 % typically, 25 % at worst, against 20 ... 1700 % for the quadratic form). The no-load point gives
    V_g(no load) (the rotor's field alone), exactly. Returns {'v0' (7,), 'g_coef' (4,), 'a_coef' (3, 7), 'k_eff' (7,), 'loaded', 'loo_total' (n loaded: predicted - true V_total), 'loo_total_rel',
    'loo_group_rel' (n loaded, 7: error of each group's V in per cent of the point's total), 'rms_total_rel'}."""
    at = np.array([q["at"] for q in points], float)
    gam = np.array([q["gamma"] for q in points], float)
    y = np.asarray(values, float)
    no_load = int(np.flatnonzero(at == 0.0)[0])
    loaded = np.flatnonzero(at > 0.0)
    yg = group_sums(y)                                           # (n, 7)
    v0 = yg[no_load]
    ghat = (yg[loaded] - v0[None, :]) / (at[loaded] ** 2)[:, None]       # (n loaded, 7): the armature part per AT^2
    tot = ghat.sum(axis=1)
    shares = ghat / tot[:, None]

    def solve(idx):
        xg = ghat_basis(at[loaded][idx], gam[loaded][idx])
        xa = share_basis(at[loaded][idx], gam[loaded][idx])
        cg, *_ = np.linalg.lstsq(xg, tot[idx], rcond=None)
        ca, *_ = np.linalg.lstsq(xa, shares[idx], rcond=None)
        return cg, ca

    def predict(cg, ca, i):
        gt = ghat_basis(at[loaded][i], gam[loaded][i]) @ cg
        sh = np.maximum(share_basis(at[loaded][i], gam[loaded][i]) @ ca, 0.0)
        sh = sh / sh.sum()
        return v0 + at[loaded][i] ** 2 * gt * sh

    cg, ca = solve(np.arange(len(loaded)))
    loo = np.zeros((len(loaded), len(GROUPS)))
    fit = np.zeros_like(loo)
    for i in range(len(loaded)):
        cg_i, ca_i = solve(np.arange(len(loaded)) != i)
        loo[i] = predict(cg_i, ca_i, i) - yg[loaded][i]
        fit[i] = predict(cg, ca, i) - yg[loaded][i]
    total = yg[loaded].sum(axis=1)
    return {"v0": v0, "g_coef": cg, "a_coef": ca, "k_eff": group_k_eff(y[loaded]), "loaded": loaded, "loo_total": loo.sum(axis=1), "loo_total_rel": loo.sum(axis=1) / total,
            "loo_group_rel": loo / total[:, None], "rms_total_rel": float(np.sqrt(np.mean((fit.sum(axis=1) / total) ** 2)))}


def main() -> int:
    ap = argparse.ArgumentParser(description="the machine's field at the copper over its operating range (module docstring)")
    ap.add_argument("--stage", choices=("all", "angles", "fields", "analyse"), default="all")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    if args.stage in ("all", "angles") and not (OUT / "mtpa_angles.json").exists():
        run_angles(args.workers)
    mtpa = {float(k): v for k, v in json.loads((OUT / "mtpa_angles.json").read_text(encoding="utf-8"))["angles"].items()}
    points = operating_points(mtpa)
    (OUT / "points.json").write_text(json.dumps(points, indent=1), encoding="utf-8")
    if args.stage in ("all", "fields"):
        run_fields(points, args.workers)
    if args.stage in ("all", "analyse"):
        from optimization import ac_machine_report
        ac_machine_report.report(points)
    print(f"done in {time.perf_counter() - started:.1f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
