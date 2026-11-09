"""optimization.winding_femm: the FEMM side of the winding study: the rated point, the field the copper sits in, and the proof that the turns do not change it.

What it is, and why
    The winding study (optimization/winding.py) picks the turns, the wire, the strands in hand and the parallel paths for a FIXED geometry. All of that is arithmetic
    except three things only FEMM can say, and this file asks FEMM for them in one session of the solver slot (116 solves, about 25 s on 5 workers):
      the rated point          the geometry at 100 ampere-turns a coil side (10 A x 10 turns) at its own MTPA angle, fine mesh, by design.Evaluator. Its mean torque is the
                               rated torque of the whole study; its flux linkages psi_d and psi_q are what every winding's corner speed is scaled from.
      the field of the copper  what a conductor in a slot feels at that point over one electrical period: G, the mean over the copper area of the mean square of dB/d(theta_e).
                               It is the one thing the proximity loss needs from FEMM.
      the turns check          the same geometry with 5, 10 and 20 turns at 20, 10 and 5 A. The winding model scales ONE solved map to every winding by 'flux linkage with the
                               series turns, current with 1 / turns'. That is right only if the turns do not change the field, so the study solves it instead of assuming it.

The method
    1. The rated point. design.Evaluator(pool, run='winding/<label>', scoring='AT', fine=True) scores the geometry in 23 solves: the no-load pair (0 and 30 electrical degrees), the
       mean torque of that pair at three current angles (10, 25, 40 degrees) and a parabola through them for the best angle, then the 15 rotor positions of one ripple period, 4 degrees
       apart, at it. The row gives gamma, the mean torque, the ripple, the mean psi_d and psi_q, i_d, i_q and the power factor; check() gives the layer's area.
    2. The field of the copper (rated_field, two pool tasks). A conductor feels the CHANGE of the field, not the field. A strand of diameter d in a field B(t) across its axis has
       an axial eddy current density J_z = -sigma y dB/dt (y across the strand; the net current is zero), so its loss per unit length is sigma (dB/dt)^2 integral y^2 dA =
       pi sigma d^4 (dB/dt)^2 / 64. What the winding model needs is therefore the mean of (dB/dt)^2, over the copper area (a strand lies anywhere in its layer) and over time. The
       field energy of a coil block (FEMM's block integral 2, which is the integral of B^2 / 2 mu_0 because the coil is mu_r = 1) gives only the mean of B^2: it does not know how fast B
       changes, and it counts a steady part (the magnets' stray field) that makes no eddy current. So the time series of B at every point is taken instead: the model is solved at n rotor
       positions equally spaced over ONE electrical period (the rotor turns, and the currents turn with it, FemmModel.solve; the 90 degree sector is one electrical period, and the field
       of a slot repeats after it exactly), and B is read from FEMM's own .ans. The solution file holds A at the nodes; a first-order triangle has one B, curl A (femfile.element_b), and the
       triangles of a coil layer are found from the block labels, matched to the regions by the label POINT (the file saves the <No Mesh> region as a hole, which moves every later label
       one place). Only the air-gap element turns with the rotor, so the mesh is the same at every position; that is checked, as is that every coil triangle's centre lies inside its
       layer's polygon. Each triangle's B(theta), theta = theta_e, is then taken apart in its Fourier series, c_k = (1 / n) sum_j B(theta_j) exp(-i k theta_j), and Parseval gives
           G = sum over triangles (weighted by area) of  sum_{k >= 1} 2 k^2 (|cx_k|^2 + |cy_k|^2),     [T^2 per rad^2 of electrical angle]
       the mean square of dB/d(theta) over a period (2 counts k and -k; the Nyquist term, whose phase n samples cannot tell, is taken as cos(k theta), weight k^2 / 2). With
       theta = omega t, <(dB/dt)^2> = omega^2 G. The sinusoidal reading, the 'B-hat^2 / 2' of the textbook formula, is the k = 1 term alone, G1. Two checks against FEMM itself, at the
       first position: B at the centre of 24 triangles (mo_getb with smoothing off gives the triangle's own value) and the field energy of each of the 18 coil layers (block integral 2)
       against the sum over its triangles of area x B^2 / (2 mu_0) x depth, with the layer areas (block integral 5). A third, independent one: forward differences between neighbouring samples.
       Two tasks of 24 positions each, 15 degrees apart (0, 15, ... 345: the brief's set) and 7.5, 22.5, ... 352.5: together 48 positions, 7.5 degrees apart, so that the sampling of the
       period can be judged: G is made from 12, 24 and 48 of them and the 48 are used.
    3. The turns check (turns_scaling, three pool tasks): the 15 ripple positions of design.score at the rated angle (the angle is given, not searched, so that the runs differ in nothing
       but the turns), for 5, 10 and 20 turns at 20, 10 and 5 A: mean torque, ripple (design.fourier_pkpk), mean psi_d, power factor, and the field of the copper at the first position,
       triangle by triangle.
    Pool tasks: femm_queue.femm_pool maps ONE function over its items, so run_task(task) looks at task['kind'] and calls rated_field or turns_scaling. Every task writes a model file of its
    own (design.Evaluator's rule: two FEMMs never share a name) and removes it. Nothing here starts FEMM except the pool, and nothing ends a process.

Run
    Not on its own: winding.py calls run_stage(overrides, folder, label). That opens the pool (it waits, printing 'waiting for the solver slot', while another run holds the slot), runs the
    Evaluator, then the five tasks together, and writes <folder>/femm.json (every number that comes back) and field.npz (the triangles' B at the 48 positions, so the statistics can be made
    again without FEMM). winding.py --reuse-femm reads femm.json and needs no solver slot.

Numbers to expect (this PC, 2026-11-08 / 09, FEMM 5.0, the reference geometry; the same to the last digit in three runs, so FEMM 5.0 is deterministic here)
    rated point      gamma 22.457 deg, i_d -3.8200 A, i_q 9.2416 A, T_mean 1.12880 N.m (ripple 18.45 %, pf 0.8099), psi_pm 16.638 mWb, L_d 1.1669 mH, L_q 2.1400 mH, mean psi_d 12.1803 mWb,
                     mean psi_q 19.7770 mWb; the torque from the flux linkages is 1.12868 N.m (-0.011 %); 15,859 nodes; 10.9 s
    the layer        25.737 mm^2 (motorgeo.check; the two layers of a slot are equal by construction); the triangles of the 18 layers sum to 25.725 ... 25.737 mm^2 (the round slot
                     bottom is chords in the mesh); 5,869 of the 30,610 triangles are coil triangles
    the checks       B from the .ans against mo_getb at 24 triangles: 2.3e-14 T; layer energies against block integral 2: 2.3e-14 relative; areas: 2e-15
    G (T^2/rad^2)    12 positions 5.9848e-4; 24 positions 6.4955e-4 (the shifted set 6.4947e-4); 48 positions 6.6824e-4 (used); the fundamental alone 5.4126e-4 (81 % of G);
                     forward differences 6.5118e-4 (0.974 of the spectral value). Harmonics of G: k = 1 81.0 %, 3 1.5 %, 5 6.5 %, 9 5.8 % (36 slots / 4 pole pairs = 9: the slot-passing order), 13 2.0 %, 15 1.0 %, k >= 17 1.7 %;
                     only odd k, the steady part of B is 2.7e-13 T^2 of 5.45e-4. Each doubling of the positions adds 0.34 of the last gain (+8.5 %, +2.9 %): the true G is perhaps 1.5 % higher.
                     The layers next to the bore feel G = 10.76 ... 13.12e-4, the deep ones 1.10 ... 1.12e-4: eleven times less. B_rms in the copper 23.3 mT, peak 0.124 T.
    the turns check  5, 10, 20 turns at 20, 10, 5 A: T_mean 1.1288033440 N.m, ripple 18.44840354 %, pf 0.8098925635 in all three; psi_d 6.09014479, 12.18028958, 24.36057916 mWb (times 2, 1, 1/2
                     = 12.18028958 each); the field of the copper identical, 0.0 T apart (up to 0.1092 T).
    time             116 solves: rated point 10.9 ... 11.8 s, then the five tasks 13.8 ... 14.3 s (24 positions 13.7 ... 14.3 s, 15 positions 8.1 ... 8.3 s): 25.2 ... 26.2 s of FEMM in all. The
                     slot was free for the first run (0.5 s) and held by another sweep for the last (449 s of waiting; the pool prints 'waiting for the solver slot' every 30 s).

Self-check (the answers follow).
  1. The coil blocks' field energy is one FEMM call per layer. Why is it not enough for the proximity loss?
  2. The three turns runs agree to the last digit. Why, and why is that not a sign that FEMM ignored the turns?
  3. G(24) is 2.8 % below G(48), and G(12) 9 % below G(24). Why does sampling the period more finely raise G?
Answers.
  1. It gives the mean of B^2 over the layer at one position. A conductor loses power in proportion to (dB/dt)^2: how fast the field changes, which B^2 does not say (the same B^2 can be slow
     or fast), and it counts a steady field, which makes no eddy current. Here the steady part is nil, but G (6.68e-4) is still 23 % larger than the mean of B^2 (5.45e-4), because the
     harmonics change faster in proportion to k, and their squares count with k^2: 19 % of G comes from harmonics that carry 0.65 % of B^2.
  2. FEMM multiplies each circuit's current by the block's turns, so the source in every slot is N I = 100 A.turns in all three runs, bit for bit (20 x 5 = 10 x 10 = 5 x 20): the same mesh, the
     same equations, the same A. It did not ignore the turns: the flux linkage, which is the turns times the flux, is 6.09, 12.18 and 24.36 mWb.
  3. A harmonic of order k changes k times as fast as the fundamental, so it weighs k^2 in the derivative's mean square. n samples see only k up to n / 2: the 12-position set cannot see
     k > 6, and 24 positions cannot tell k = 13 ... 23, which carry 4.8 % of G (the 48-position shares). A harmonic above the limit is not just lost: it is folded onto a lower k, where its
     weight is smaller, so G(24) is 2.8 % low, not 4.8 %.
"""
from __future__ import annotations

import json
import math
import os
import time
import traceback
from pathlib import Path

import numpy as np

import femfile
import femm_queue
import motorgeo as mg
from optimization import design

MU0 = 4e-7 * math.pi                                          # H/m
FIELD_POSITIONS_A = tuple(15.0 * j for j in range(24))        # electrical degrees: one full electrical period, 15 degrees apart (the brief's set)
FIELD_POSITIONS_B = tuple(7.5 + 15.0 * j for j in range(24))  # the same period shifted by half a step: with A it makes 48 positions 7.5 degrees apart
TURNS_CHECK = (5, 10, 20)                                     # the turns of the invariance check; the current is 100 / turns, the ampere-turns stay 100
N_PROBE = 24                                                  # triangles whose B is also read from FEMM itself, to check the .ans route


# ============================================================================================== the field of the copper (pure Python)
def coil_triangles(sol: dict, regions: list) -> dict:
    """Which triangles of a solution belong to a coil layer, and everything about them that does not change with the rotor.

    sol     femfile.read_ans(...) of a solved model of this design
    regions mg.build(...) of the same design
    Returns {"mask" (bool per triangle), "area_mm2", "centroid" (n, 2) [mm], "layer" (index into names), "names" (coil_k_j), "outside"}.
    A triangle is a layer's when its block label (a number into the file's label list) is that layer's label; the label list is matched to
    the regions by the label POINT (the file keeps every label's x y), not by position in the list, because the <No Mesh> region is saved as
    a hole and moves everything after it by one place. `outside` counts the triangles whose centroid is not inside their layer's polygon:
    0 says the matching is right."""
    import shapely
    layers = [r for r in regions if r.material == "Coil"]
    names = [r.name for r in layers]
    label_of = {}
    for j, r in enumerate(layers):
        hits = [i for i, lab in enumerate(sol["labels"]) if lab["material"] == "Coil"
                and abs(lab["xy"][0] - r.label[0]) < 1e-6 and abs(lab["xy"][1] - r.label[1]) < 1e-6]
        if len(hits) != 1:
            raise RuntimeError(f"coil layer {r.name}: {len(hits)} block labels of the .ans sit at its label point {r.label}")
        label_of[hits[0]] = j
    tri_label = sol["label"]
    layer_all = np.full(len(tri_label), -1, int)
    for i, j in label_of.items():
        layer_all[tri_label == i] = j
    mask = layer_all >= 0
    corners = sol["xy"][sol["tri"][mask]]                                   # (n, 3, 2)
    d1, d2 = corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0]
    area = 0.5 * np.abs(d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0])
    centroid = corners.mean(axis=1)
    layer = layer_all[mask]
    outside = 0
    for j, r in enumerate(layers):
        sel = layer == j
        if sel.any():
            outside += int((~shapely.contains_xy(r.polygon(), centroid[sel, 0], centroid[sel, 1])).sum())
    return {"mask": mask, "area_mm2": area, "centroid": centroid, "layer": layer, "names": names, "outside": outside}


def field_statistics(bx: np.ndarray, by: np.ndarray, area_mm2: np.ndarray, layer: np.ndarray, n_layers: int) -> dict:
    """What the copper feels: the mean over the copper area of the mean square of dB/dtheta_e, and the same cut up by harmonic and by layer.

    bx, by   (n positions, n elements): the flux density of every coil triangle at n rotor positions equally spaced over ONE electrical
             period (position j is theta_e = j 360/n degrees, the currents turning with the rotor); the field repeats after a period
    area_mm2 the triangles' areas: every triangle counts by its area, so the copper area is averaged, not the mesh
    layer    which coil layer each triangle belongs to (0 .. n_layers - 1)

    A conductor feels only the change of the field, never the field itself, so each triangle's time series is taken apart in its Fourier
    series, B(theta) = sum_k c_k exp(i k theta), with c_k = (1/n) sum_j B(theta_j) exp(-i k theta_j). The mean square of the derivative over
    a period is then Parseval's sum:  <(dB/dtheta)^2> = sum over k >= 1 of 2 k^2 (|cx_k|^2 + |cy_k|^2)  (the factor 2 counts k and -k).
    n samples resolve k = 1 .. n/2 - 1 and the Nyquist term k = n/2, which cannot tell its phase and is taken as cos(k theta):
    the mean square of that derivative is k^2 |c|^2 / 2. The constant term k = 0 (the magnets' steady stray field) drops out of a derivative.
    Returns, G in T^2 per rad^2 of electrical angle (multiply by omega_e^2 for T^2 per s^2):
      G, G_k (k = 0 .. n/2; G_0 = 0), G1 (the fundamental alone: the sinusoidal reading, B-hat^2 / 2), G_layer (per layer),
      B2_mean (the time and area mean of |B|^2), B2_dc (the part of it that is constant in time), B2_1 (the fundamental's share), B_peak
      (the largest |B| of any triangle at any position), area_total_mm2, area_layer_mm2."""
    n = bx.shape[0]
    cx, cy = np.fft.rfft(bx, axis=0) / n, np.fft.rfft(by, axis=0) / n          # (n/2 + 1, n elements)
    k = np.arange(n // 2 + 1, dtype=float)
    w = 2.0 * k ** 2
    w[0] = 0.0
    if n % 2 == 0:
        w[-1] = 0.5 * k[-1] ** 2                                                # the Nyquist term (see above)
    power = np.abs(cx) ** 2 + np.abs(cy) ** 2                                   # |cx_k|^2 + |cy_k|^2 per harmonic and element
    g_el = w[:, None] * power                                                   # each harmonic's share of <(dB/dtheta)^2> per element
    a = np.asarray(area_mm2, float)
    weight = a / a.sum()
    G_k = g_el @ weight
    per_el = g_el.sum(axis=0)
    G_layer = np.array([float(per_el[layer == j] @ a[layer == j] / a[layer == j].sum()) for j in range(n_layers)])
    b2 = (bx ** 2 + by ** 2).mean(axis=0)                                       # time mean of |B|^2 per element
    b2_1 = 2.0 * power[1]                                                       # the mean square of the fundamental part of B
    step = 2.0 * math.pi / n                                                    # an independent estimate: forward differences between neighbouring samples
    fd = ((np.roll(bx, -1, axis=0) - bx) ** 2 + (np.roll(by, -1, axis=0) - by) ** 2).mean(axis=0) / step ** 2
    return {"n_positions": n, "n_elements": int(bx.shape[1]), "G": float(G_k.sum()), "G_k": G_k.tolist(), "G1": float(G_k[1]),
            "G_forward_difference": float(fd @ weight),
            "G_layer": G_layer.tolist(), "B2_mean": float(b2 @ weight), "B2_dc": float(power[0] @ weight), "B2_1": float(b2_1 @ weight),
            "B_peak": float(np.sqrt(bx ** 2 + by ** 2).max()), "area_total_mm2": float(a.sum()),
            "area_layer_mm2": [float(a[layer == j].sum()) for j in range(n_layers)]}


# ============================================================================================== the pool's tasks (they run in a worker, on its own FEMM)
def _fem_paths(task: dict) -> tuple[str, str]:
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")                  # a file of its own per task: two FEMMs never share a name
    return fem, fem[:-4] + ".ans"


def _cleanup(femm, fem: str, ans: str) -> None:
    try:
        femm.mi_close()
    except Exception:
        pass
    for path in (ans, fem):
        try:
            os.remove(path)
        except OSError:
            pass


def _refuse_if_rejected(p, regions, task: dict) -> None:
    """motorgeo.check on the design, as bench/femm_points.solve_points does: a design motorgeo rejects is an error that names the reason, never a number."""
    _, violations = mg.check(p, regions)
    if violations:
        raise RuntimeError(f"{task['tag']}: motorgeo.check rejects this design {task['overrides']}: {'; '.join(violations)}")


def rated_field(task: dict) -> dict:
    """One pool task: the flux density in every coil triangle of the sector, at every rotor position of a list.

    task = {"tag", "overrides", "fem_dir", "i_d", "i_q", "positions"}: the geometry, the d-q currents of the rated point, the rotor positions
    [electrical degrees]. The worker draws the model on the fine mesh, solves each position with the currents turning with the rotor
    (FemmModel.solve), reads the solution FEMM wrote (femfile.read_ans: first-order triangles, so B is one number per triangle) and keeps
    the coil triangles' B. The mesh must be the same at every position (only the air-gap element turns); it is checked, not assumed.
    Two checks against FEMM itself at the first position: B at the middle of 24 triangles (FEMM's mo_getb with smoothing off gives the
    triangle's own value) and the field energy of every layer (block integral 2, = the sum over triangles of area x |B|^2 / (2 mu0) x depth
    because the coil material has mu_r = 1) with its cross-section area (block integral 5)."""
    import femm
    started = time.perf_counter()
    fem, ans = _fem_paths(task)
    try:
        p = design.params(task["overrides"])
        regions = mg.build(p)
        _refuse_if_rejected(p, regions, task)
        mg.to_femm(p, regions, fem, fine=True)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        layers = [r for r in regions if r.material == "Coil"]
        bx_all, by_all, psi_d, psi_q, torque = [], [], [], [], []
        static = None
        checks: dict = {}
        for n, theta in enumerate(task["positions"]):
            row = model.solve(theta, task["i_d"], task["i_q"])
            sol = femfile.read_ans(ans)
            bx, by = femfile.element_b(sol["xy"], sol["A"], sol["tri"])
            if static is None:
                static = coil_triangles(sol, regions)
                static["tri"] = sol["tri"][static["mask"]]
                static["xy"] = sol["xy"]
                static["n_triangles"] = int(len(sol["tri"]))
                if static["outside"]:
                    raise RuntimeError(f"{static['outside']} coil triangles lie outside their layer's polygon: the label matching is wrong")
            elif not (np.array_equal(sol["tri"][static["mask"]], static["tri"]) and np.array_equal(sol["xy"], static["xy"])):
                raise RuntimeError(f"the mesh at theta_e = {theta} differs from the one at {task['positions'][0]}: the elements are not comparable")
            m = static["mask"]
            bx_all.append(bx[m])
            by_all.append(by[m])
            psi_d.append(row["psi_d"])
            psi_q.append(row["psi_q"])
            torque.append(row["T_age"])
            if n == 0:
                checks = _check_against_femm(femm, layers, static, bx[m], by[m], p)
        nodes = int(len(static["xy"]))
    except Exception as error:                                   # one bad task must not end the run: say so in the result
        return {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4),
                "seconds": time.perf_counter() - started}
    finally:
        _cleanup(femm, fem, ans)
    return {"status": "ok", "kind": "field", "positions": list(task["positions"]), "bx": np.array(bx_all), "by": np.array(by_all),
            "area_mm2": static["area_mm2"], "centroid": static["centroid"], "layer": static["layer"], "names": static["names"],
            "psi_d": psi_d, "psi_q": psi_q, "T_age": torque, "nodes": nodes, "triangles": static["n_triangles"],
            "coil_triangles": int(len(static["area_mm2"])), "checks": checks, "seconds": time.perf_counter() - started}


def _check_against_femm(femm, layers, static: dict, bx: np.ndarray, by: np.ndarray, p) -> dict:
    """The two checks of rated_field, at the position just solved (its .fem is still open in this FEMM)."""
    femm.mi_loadsolution()
    try:
        femm.mo_smooth("off")                                    # FEMM smooths B between triangles unless told not to
        rng = np.random.default_rng(1)
        pick = rng.choice(len(bx), N_PROBE, replace=False)
        worst = 0.0
        for i in pick:
            fb = femm.mo_getb(float(static["centroid"][i, 0]), float(static["centroid"][i, 1]))
            worst = max(worst, math.hypot(float(np.real(fb[0])) - bx[i], float(np.real(fb[1])) - by[i]))
        depth_m = p.stack * 1e-3
        energy, energy_mine, area_femm, area_mine = [], [], [], []
        for j, r in enumerate(layers):
            femm.mo_selectblock(*r.label)
            energy.append(complex(femm.mo_blockintegral(2)).real)
            area_femm.append(complex(femm.mo_blockintegral(5)).real)
            femm.mo_clearblock()
            sel = static["layer"] == j
            a_m2 = static["area_mm2"][sel] * 1e-6
            energy_mine.append(float((a_m2 * (bx[sel] ** 2 + by[sel] ** 2)).sum() * depth_m / (2 * MU0)))
            area_mine.append(float(static["area_mm2"][sel].sum()))
    finally:
        femm.mo_close()
    return {"B_probe_max_diff_T": worst, "n_probe": N_PROBE, "energy_femm_J": energy, "energy_mine_J": energy_mine,
            "area_femm": area_femm, "area_mine_mm2": area_mine}


def turns_scaling(task: dict) -> dict:
    """One pool task: one geometry with `turns` turns a coil side, at the ampere-turns of the rated point and a given current angle.

    task = {"tag", "overrides" (with "turns"), "fem_dir", "amps", "gamma", "positions"}. Solves the 15 ripple positions of design.score at
    that fixed angle (no search of its own: the three runs of the check must differ in nothing but the turns) and returns the mean torque,
    the ripple (design.fourier_pkpk), the mean flux linkages, the power factor and the field of the copper at the first position
    (bx, by per coil triangle) so that the three runs' fields can be set side by side."""
    import femm
    started = time.perf_counter()
    fem, ans = _fem_paths(task)
    try:
        p = design.params(task["overrides"])
        regions = mg.build(p)
        _refuse_if_rejected(p, regions, task)
        mg.to_femm(p, regions, fem, fine=True)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        i_d, i_q = design.currents(float(task["amps"]), float(task["gamma"]))
        rows, field = [], None
        for n, theta in enumerate(task["positions"]):
            rows.append(model.solve(theta, i_d, i_q))
            if n == 0:
                sol = femfile.read_ans(ans)
                static = coil_triangles(sol, regions)
                bx, by = femfile.element_b(sol["xy"], sol["A"], sol["tri"])
                field = (static["centroid"], bx[static["mask"]], by[static["mask"]])
    except Exception as error:
        return {"status": "error", "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=4),
                "seconds": time.perf_counter() - started}
    finally:
        _cleanup(femm, fem, ans)
    T = [r["T_age"] for r in rows]
    psi_d, psi_q = float(np.mean([r["psi_d"] for r in rows])), float(np.mean([r["psi_q"] for r in rows]))
    pp = design.fourier_pkpk(T)
    pf = (psi_d * i_q - psi_q * i_d) / (math.hypot(psi_d, psi_q) * math.hypot(i_d, i_q))
    return {"status": "ok", "kind": "turns", "turns": int(p.turns), "amps": float(task["amps"]), "gamma": float(task["gamma"]), "i_d": i_d, "i_q": i_q,
            "T": T, "T_mean": float(np.mean(T)), "ripple_pp": pp, "ripple_pct": 100.0 * pp / float(np.mean(T)), "psi_d": psi_d, "psi_q": psi_q,
            "pf": pf, "T_dq_mean": float(np.mean([r["T_dq"] for r in rows])), "centroid": field[0], "bx": field[1], "by": field[2],
            "seconds": time.perf_counter() - started}


def run_task(task: dict) -> dict:
    """The one function the pool maps (a pool maps one function over its items, so the kind of task says which one to run)."""
    if task["kind"] == "field":
        return rated_field(task)
    if task["kind"] == "turns":
        return turns_scaling(task)
    raise ValueError(f"unknown kind of task: {task['kind']!r}")


# ============================================================================================== the stage: one pool session, everything FEMM
def jsonable(x):
    """numpy scalars and arrays, and Paths, as plain JSON values (a number that is not finite becomes the text 'inf', '-inf' or 'nan': JSON has no such numbers)."""
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return jsonable(x.tolist())
    if isinstance(x, (np.floating, np.integer)):
        return jsonable(x.item())
    if isinstance(x, float) and not math.isfinite(x):
        return str(x)
    if isinstance(x, Path):
        return str(x)
    return x


def analyse(field_a: dict, field_b: dict, turn_results: list[dict]) -> dict:
    """Numbers from what the pool brought back (no FEMM): the field statistics at 12, 24 and 48 rotor positions, FEMM's own checks, the turns check.

    field_a, field_b   the two rated_field results: positions 0, 15, ... 345 and 7.5, 22.5, ... 352.5 electrical degrees of the same period.
                       They must be the same mesh (the same triangles in the same order): the 48 positions are the two sets interleaved.
    The statistics are made three ways from the same solves, so the reader can see whether the period was sampled finely enough:
      12 positions  every second sample of A (30 degrees apart)
      24 positions  A alone (the brief's 15 degrees), and B alone (the same, shifted by half a step: a second opinion)
      48 positions  A and B interleaved (7.5 degrees)
    G (T^2 per rad^2) is the quantity the loss needs; the 48-position value is the one the study uses."""
    for name in ("centroid", "area_mm2", "layer"):
        if not np.array_equal(field_a[name], field_b[name]):
            raise RuntimeError(f"the two field tasks have different meshes ({name} differs): their positions cannot be interleaved")
    n_layers = len(field_a["names"])
    bx48 = np.empty((48, field_a["bx"].shape[1]))
    by48 = np.empty_like(bx48)
    bx48[0::2], bx48[1::2] = field_a["bx"], field_b["bx"]
    by48[0::2], by48[1::2] = field_a["by"], field_b["by"]
    area, layer = field_a["area_mm2"], field_a["layer"]
    stats = {"48": field_statistics(bx48, by48, area, layer, n_layers),
             "24": field_statistics(field_a["bx"], field_a["by"], area, layer, n_layers),
             "24_shifted": field_statistics(field_b["bx"], field_b["by"], area, layer, n_layers),
             "12": field_statistics(field_a["bx"][::2], field_a["by"][::2], area, layer, n_layers)}
    chk = [field_a["checks"], field_b["checks"]]
    e_f, e_m = np.array(chk[0]["energy_femm_J"]), np.array(chk[0]["energy_mine_J"])
    a_f, a_m = np.array(chk[0]["area_femm"]), np.array(chk[0]["area_mine_mm2"])
    femm_checks = {"B_probe_max_diff_T": max(c["B_probe_max_diff_T"] for c in chk), "n_probe": chk[0]["n_probe"],
                   "energy_max_rel_diff": float(np.max(np.abs(e_m / e_f - 1))), "energy_femm_J_layer0": float(e_f[0]),
                   "energy_mine_J_layer0": float(e_m[0]), "energy_total_femm_J": float(e_f.sum()), "energy_total_mine_J": float(e_m.sum()),
                   "area_ratio_femm_over_mine": float(np.mean(a_f / (a_m * 1e-6))), "area_ratio_min": float(np.min(a_f / (a_m * 1e-6))),
                   "area_ratio_max": float(np.max(a_f / (a_m * 1e-6)))}                   # FEMM's block area is in m^2, the triangles' sum in mm^2
    ref = next((r for r in turn_results if r["turns"] == 10), turn_results[0])
    turns = []
    for r in turn_results:
        if not np.array_equal(r["centroid"], ref["centroid"]):
            raise RuntimeError("the turns runs have different meshes: their fields cannot be compared element by element")
        t = {k: v for k, v in r.items() if k not in ("centroid", "bx", "by", "kind", "status")}
        t["field_max_diff_T"] = float(np.max(np.hypot(r["bx"] - ref["bx"], r["by"] - ref["by"])))
        t["field_max_T"] = float(np.max(np.hypot(r["bx"], r["by"])))
        turns.append(t)
    return {"statistics": stats, "femm_checks": femm_checks, "turns_check": turns, "turns_reference": int(ref["turns"]),
            "layers": field_a["names"], "mesh": {"nodes": field_a["nodes"], "triangles": field_a["triangles"], "coil_triangles": field_a["coil_triangles"]},
            "positions": {"A": list(field_a["positions"]), "B": list(field_b["positions"])},
            "T_age_A": field_a["T_age"], "psi_d_A": field_a["psi_d"], "psi_q_A": field_a["psi_q"]}


def run_stage(overrides: dict, folder: Path, label: str, log=print) -> dict:
    """Everything of the study that needs FEMM, in ONE session of the solver slot (the slot is machine-wide: other runs wait while this holds it).

    1. the rated point: design.Evaluator on the fine mesh (its own MTPA angle, then 15 rotor positions): gamma, mean torque, flux linkages;
    2. then, in parallel, five tasks: the field of the copper over one electrical period at 24 positions (rated_field, set A) and at 24 more,
       half a step on (set B), and the turns check (turns_scaling with 5, 10 and 20 turns at 20, 10 and 5 A, at the rated angle).
    Writes <folder>/femm.json (every number that comes back) and <folder>/field.npz (the triangles' B, so the statistics can be made again
    without FEMM) and returns the same dict."""
    folder.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    with femm_queue.femm_pool(f"winding study: {label}", workers=5, task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        waited = time.perf_counter() - t0
        t1 = time.perf_counter()
        evaluate = design.Evaluator(pool, run=f"winding/{label}", scoring="AT", fine=True)
        row = evaluate([overrides], stage="rated")[0]
        if not row.get("feasible") or "T_mean" not in row:
            raise RuntimeError(f"the geometry {overrides} cannot be scored: {row.get('reason', 'no reason given')}")
        rated_seconds = time.perf_counter() - t1
        log(f"rated point: gamma {row['gamma']:.3f} deg, T_mean {row['T_mean']:.5f} N.m, ripple {row['ripple_pct']:.2f} %, "
            f"psi_pm {1e3 * row['psi_pm']:.3f} mWb, pf {row['pf']:.4f} ({rated_seconds:.1f} s)")
        at = float(row["amps"]) * design.params(overrides).turns                      # ampere-turns a coil side: 100 at the reference
        tag = f"winding_{label}".replace("/", "_")
        fem_dir = str(folder / "fem")
        tasks = [{"kind": "field", "tag": f"{tag}_fieldA", "overrides": overrides, "fem_dir": fem_dir, "i_d": float(row["i_d"]),
                  "i_q": float(row["i_q"]), "positions": list(FIELD_POSITIONS_A)},
                 {"kind": "field", "tag": f"{tag}_fieldB", "overrides": overrides, "fem_dir": fem_dir, "i_d": float(row["i_d"]),
                  "i_q": float(row["i_q"]), "positions": list(FIELD_POSITIONS_B)}]
        for turns in TURNS_CHECK:
            tasks.append({"kind": "turns", "tag": f"{tag}_turns{turns}", "overrides": {**overrides, "turns": turns}, "fem_dir": fem_dir,
                          "amps": at / turns, "gamma": float(row["gamma"]), "positions": list(design.RIPPLE_POSITIONS)})
        t2 = time.perf_counter()
        results = pool.map(run_task, tasks)
        task_seconds = time.perf_counter() - t2
    bad = [(t["tag"], r) for t, r in zip(tasks, results) if not (isinstance(r, dict) and r.get("status") == "ok")]
    if bad:
        raise RuntimeError("FEMM tasks failed: " + "; ".join(f"{tag}: {r.get('error', r) if isinstance(r, dict) else r}" for tag, r in bad))
    field_a, field_b, turn_results = results[0], results[1], results[2:]
    p = design.params(overrides)
    rated = {k: row[k] for k in ("gamma", "T_mean", "T", "ripple_pct", "ripple_pp", "psi_pm", "psi_d", "psi_q", "L_d", "L_q", "saliency", "pf",
                                 "T_dq_mean", "i_d", "i_q", "amps", "nodes", "seconds", "fine")}
    rated["coil_layer_area_mm2"] = float(row["metrics"]["coil_layer_area_mm2"])
    rated["ampere_turns_per_coil_side"] = at
    rated["turns"] = int(p.turns)
    data = {"overrides": overrides, "rated": rated, **analyse(field_a, field_b, turn_results),
            "seconds": {"waiting_for_slot": waited, "rated": rated_seconds, "tasks": task_seconds, "total": time.perf_counter() - t0,
                        "field_tasks": [field_a["seconds"], field_b["seconds"]], "turns_tasks": [r["seconds"] for r in turn_results]},
            "solves": {"rated": 23, "field": len(FIELD_POSITIONS_A) + len(FIELD_POSITIONS_B), "turns": len(TURNS_CHECK) * len(design.RIPPLE_POSITIONS)}}
    data["solves"]["total"] = sum(data["solves"].values())
    np.savez_compressed(folder / "field.npz", bx_A=field_a["bx"], by_A=field_a["by"], bx_B=field_b["bx"], by_B=field_b["by"],
                        area_mm2=field_a["area_mm2"], centroid=field_a["centroid"], layer=field_a["layer"],
                        positions_A=np.array(field_a["positions"]), positions_B=np.array(field_b["positions"]))
    (folder / "femm.json").write_text(json.dumps(jsonable(data), indent=1), encoding="utf-8")
    return data
