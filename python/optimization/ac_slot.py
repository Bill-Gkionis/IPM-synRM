"""optimization.ac_slot: one slot of the optimized stator in FEMM, strand by strand: the winding's AC resistance, the proximity loss it has there, and the circulating
currents between strands. The FEMM side of the copper-loss work; its studies are in optimization/ac_slot_study.py.

What it is, and why
    The first winding study (optimization/winding.py) put the AC loss in one low-frequency formula fed by the field of a HOMOGENISED coil (a block of uniform current density,
    the strands never drawn) and found it negligible at 3,000 rpm. That left four things open: what the formula does at higher frequency (the PWM switching frequency,
    10 ... 20 kHz, and the machine's slot harmonics), whether it holds strand by strand (the strands near the slot opening feel a field that changes within one strand
    diameter), what the neighbouring strands add (their own eddy fields and the field of their currents), and what happens when the strands of one turn are joined in
    parallel (circulating currents). This module answers them the direct way: it draws one slot of the optimized stator with every round strand of both coil layers in it,
    and lets FEMM's time-harmonic solver find the eddy currents in each strand.

The model (build_model, in the FEMM of a pool worker)
    Geometry   the real slot outline of motorgeo (motorgeo._slot: the tooth tip, the 1.5 mm opening, the radial walls at +-3.16 deg, the round bottom), a cell of ONE SLOT
               PITCH (10 deg: from the centre line of one tooth to the centre line of the next) with the stator steel in it up to the outer radius, the 0.8 mm air gap above the
               opening, and the rotor surface at the rotor radius. The slot is drawn pointing along +y (x across the slot, y radial: the depth into the slot is y).
    Boundaries the outer radius: A = 0 (no flux leaves). The rotor's surface: FEMM's natural boundary (field lines hit it at right angles: a rotor of infinite permeability).
               The two cuts through the tooth centres: A = 0 by default ('dirichlet': the field lines run along the cut) or the natural boundary ('neumann': they cross it at right
               angles). Neither is the real neighbouring slots. The Neumann cell is ill-posed for a slot that carries a net current (the slot's ampere-turns have no way to close their
               flux through natural boundaries), and in the first prototype it gave 3.56 for R_ac / R_dc at 20 kHz where the Dirichlet cell gave 4.37: the studies use Dirichlet cuts, and
               report the cell's width as a sensitivity (ac_slot_study.py: 1, 3 and 5 pitches give 4.420, 4.397 and 4.301 at 20 kHz for the 26 x 8 x 0.315 mm winding, +0.5 % and -2.2 %
               of the 3-pitch value).
    Steel      linear: mu_r 4416 (the M-19 curve's initial slope in this project), conductivity 1.9 MS/m, lamination thickness 0.356 mm (29 gauge, the working assumption of
               the materials guide), lamination fill 0.98. FEMM's time-harmonic solver handles a laminated block as one material with a complex effective permeability (the
               one-dimensional eddy-current solution through the sheet), so the eddy currents of the laminations screen the slot field at high frequency without a mesh fine
               enough for the skin depth. The real teeth are saturated by the rotor's flux; this slot field is the leakage field of the slot's own current on top of that,
               and its permeability is the small-signal one (the sensitivity stage varies it).
    Strands    every strand a copper disc of the bare diameter d (conductivity 1 / rho(T): 44.1 MS/m at 100 degC), in the slot's two coil layers (the layers of motorgeo: equal
               areas, 36.218 mm^2 each, the coil region between the tooth-tip line and the round bottom). The strands are packed at random (pack_disks) as discs of the
               OVERALL diameter (the enamel) that do not overlap and stay inside the layer: for the 0.315 mm strand in 0.367 mm overall, 208 strands a layer is a packing
               fraction of 0.608 of the layer; the enamel itself is not drawn (it is air to FEMM: mu_r = 1, no conductivity). Each strand carries a current FEMM imposes:
               'series' circuits (one circuit a strand, the current given: equal sharing, the ideal transposed winding) or 'parallel' circuits (one circuit a TURN: the
               n_s strands of a turn share a voltage and the turn's current: whatever their positions make of the division).
    Currents   strand currents are PEAK phasors, complex numbers (FEMM's convention). The two layers carry in-phase currents of the same direction ('same phase', the case
               of 6 of the 9 slots of the sector: the double-layer winding with coil pitch 4 puts the same phase in both layers of those) or a 60 degree phase difference (the
               other 3 slots: layers of two phases).
    Losses     the power in strand k is FEMM's block integral 'resistive losses' of its block (mo_blockintegral 4: the time average of |J|^2 / (2 sigma) over the strand, times
               the depth), its current the block integral 7, its area 5. DC reference: the same currents at zero frequency, |I_k|^2 rho L / (2 A_k) with FEMM's own polygon
               area A_k (the discs are polygons of ARC_DEG segments, so the area is 0.5 % low: with the DC loss made from the same area the ratio is not affected).

The strand-by-strand comparison (the formula against FEMM)
    Three numbers for every strand at every frequency, from three solves of the same strand layout:
      FEMM       the time-harmonic loss P_k of the strand (above).
      formula    the first study's formula applied strand by strand: P_k = pi sigma w^2 d^4 L |B_k|^2 / 128 with B_k the local field in a MAGNETOSTATIC solve of the same slot: the
                 mean over the strand's cross-section of B (peak, at the strand's current amplitudes); the strand's own field averages to zero over its disc, so B_k is the
                 field of everything else (the other strands, the slot's walls, the opening). 'formula_exact' multiplies it by the exact-over-low-frequency factor of
                 ac_theory.proximity_ratio (a lone wire's shielding of itself) and adds the lone wire's skin effect (ac_theory.skin_ratio): the whole lone-wire theory.
      low-f      the exact low-frequency loss of the same magnetostatic solution, P_k = (sigma w^2 L / 2) integral over the strand of |A - A_mean|^2 dS, with A the vector
                 potential (J = I / S - j w sigma (A - A_mean) when the currents' own fields are small; the integral includes the strand's own field and the gradient of the
                 neighbours', which the uniform-field formula does not). It is what FEMM must give as the frequency goes to 0 (checked), and the difference between it and
                 'formula' is what the field's variation over a strand adds.
    'homogenised' is the same formula with B taken from the HOMOGENISED slot (FEMM's coil block of uniform current density in each layer, the first study's model), at the
    strand's centre: the difference to 'formula' is the discreteness of the strands' currents (the field of the neighbours seen at distance).

Why not trust FEMM's meshing blindly
    optimization/ac_calibration.py runs this module's own discretisation (ARC_DEG, the mesh sizes) on the two cases with exact answers (a lone wire's skin effect, a wire in a
    uniform field), and its convergence stage repeats the slot at a finer mesh. Every number of the studies carries the mesh behind it.

How a script uses it (the pool and its slot are the caller's: nothing here starts FEMM by itself)
    from optimization import ac_slot
    with femm_queue.femm_pool("copper ac: slot", workers=8, task_timeout_s=600) as pool:
        results = pool.map(ac_slot.run_task, [ac_slot.task_ac(...), ...])
    The tasks (task_ac, task_ms, task_hom) are plain dicts; run_task(task) is the one function the pool maps (it dispatches on task['kind']) and returns a plain dict with
    numpy arrays, or {'status': 'error', ...}: a bad task never ends a run. Every task writes a model file of its own (two FEMMs never share a name) and removes it.
    The module has no command line; its checks and studies are optimization/ac_slot_study.py (python -m optimization.ac_slot_study).

Numbers to expect (this PC, FEMM 5.0.20261101, the 26 turns x 8 strands x 0.315 mm winding of the winding study's example at 118.6 ampere-turns a coil side, one random layout)
    a slot of 416 strands in the 3-pitch cell is 60,523 nodes; a whole task (drawing, solving, reading 416 strands) takes about 3 s (the full study, 245 tasks on 8 workers, 106 s of FEMM);
    R_ac / R_dc of the same-phase layers is 1.00038 at 200 Hz, 1.230 at 5 kHz and 4.397 at 20 kHz for the random layout of seed 1 (the three layouts of the study give 4.17 to 4.40 at 20 kHz,
    mean 4.27), the proximity loss of the slot's own field.

Self-check (the answers follow).
  1. The strands of a layer are drawn as discs of the OVERALL diameter packed without overlap, then shrunk to the bare diameter. Why not draw the enamel too?
  2. The DC loss that R_ac / R_dc divides by is made from FEMM's area of each strand and not from pi d^2 / 4. Why?
  3. In a 'parallel' circuit the loss of the strands is more than in 'series' circuits with the same total current. What is the extra, physically?
Answers.
  1. The enamel is a non-conductor with mu_r = 1, the same as the air between the strands to FEMM; it matters only because it keeps the copper from touching, which the
     packing already enforces. Drawing it would add regions (and elements) and change nothing in the solution.
  2. The polygon FEMM meshes (ARC_DEG segments on a circle) has an area that is 0.5 % smaller than the circle's. The DC loss and the AC loss are made of the same polygon, so their
     ratio is the physics; made from the circle's area, the ratio would carry a 0.5 % error that has nothing to do with frequency.
  3. The strands of one turn sit at different places in the slot's field, so their vector potential A differs; with a common terminal voltage each carries a current that makes up
     the difference (I_k - I_mean = -j w sigma S (A_k - A_mean) at low frequency). Those circulating currents flow in loops closed through the strands' ends and are loss with no
     torque: the extra loss against equal sharing.
"""
from __future__ import annotations

import math
import os
import time
import traceback
from pathlib import Path

import numpy as np

import motorgeo as mg
from optimization import ac_theory, design

# ---------------------------------------------------------------------------------------------- the fixed choices of the model
HALF_PITCH_DEG = 5.0              # half a slot pitch (the pitch is 10 deg): a one-slot cell reaches from tooth centre to tooth centre
SLOT_PHI_DEG = 90.0               # the slot's axis: +y (x across the slot, y radial); the packing and the plots use this frame
MU0 = 4e-7 * math.pi
DEFAULT_OPTIONS = {
    "arc_deg": 10.0,              # a strand's outline: segments of at most this many degrees (36 a circle). The mesher replaces the arcs by these chords
    "area_correct": True,         # draw the strand's circle a little larger (radius x sqrt(2 pi / (n sin(2 pi / n))), 0.25 % for n = 36) so that the polygon of the mesh has the
                                  # area of the true circle: without it every strand is 0.5 % small in area, and R_ac / R_dc (which goes as d^2) 1 % low
    "mesh_cu_over_d": 1.0 / 6.0,  # the mesh size inside a strand, in strand diameters
    "mesh_air_mm": 0.15,          # the air (the opening, the gap, the spaces between strands): the mesher refines it near the strands by itself
    "mesh_fe_mm": 0.5,            # the stator steel
    "steel": {"mu_r": 4416.0, "sigma_MS": 1.9, "d_lam_mm": 0.356, "lam_fill": 0.98},
    "n_pitch": 3,                 # the cell's width in slot pitches: 1 (tooth centre to tooth centre), 3 or 5 (the neighbouring slots drawn, empty)
    "lateral": "dirichlet",       # the cuts through the tooth centres: 'dirichlet' (A = 0) or 'neumann' (natural)
    "rotor": "neumann",           # the rotor surface: 'neumann' (natural: an infinitely permeable rotor) or 'dirichlet' (A = 0)
    "precision": 1e-8,            # FEMM's relative residual of the linear solver
}
LAYERS = 2


# ---------------------------------------------------------------------------------------------- geometry (shapely, no FEMM)
def slot_points(p: mg.Params, phi: float = SLOT_PHI_DEG) -> dict:
    """The drawing points of the slot whose axis is at the angle phi (the central slot of the cell: +y, so x is across the slot and y radial, mm): the opening corners O1 ... O4,
    the tooth-tip ends S1, S2 (where the coil region starts), the wall ends W1, W2 (where the round bottom starts), the bottom's centre cb and lowest point B."""
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    sd = p.slot_half_deg
    a_s = math.radians(sd)
    r_b = p.r_wall * math.tan(a_s)
    P = lambda r, d: mg.pol(r, phi + d)
    return {"a_o": a_o, "sd": sd, "R_b": r_b, "phi": phi,
            "O1": P(p.r_si, -a_o), "O2": P(p.r_si + p.h_open, -a_o), "O3": P(p.r_si + p.h_open, a_o), "O4": P(p.r_si, a_o),
            "S1": P(p.r_shoe, -sd), "S2": P(p.r_shoe, sd), "W1": P(p.r_wall, -sd), "W2": P(p.r_wall, sd),
            "cb": P(p.r_wall / math.cos(a_s), 0), "B": P(p.r_wall / math.cos(a_s) + r_b, 0)}


def slot_regions(p: mg.Params) -> list:
    """motorgeo's regions of one slot drawn along +y: [opening, coil layer 0 (next to the opening), coil layer 1 (the bottom)]."""
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    return mg._slot(p, 0, SLOT_PHI_DEG, a_o, (("A", 1), ("A", 1)))


def layer_polygons(p: mg.Params) -> list:
    """The two coil layers as shapely polygons (mm), layer 0 next to the opening, layer 1 at the bottom: 36.218 mm^2 each for the optimized design."""
    return [r.polygon() for r in slot_regions(p) if r.name.startswith("coil")]


def pack_disks(poly, n: int, d_over: float, rng: np.random.Generator, tries: int = 25, max_iter: int = 4000) -> np.ndarray:
    """n centres (n, 2) [mm] of non-overlapping discs of diameter d_over inside the polygon, packed at random.
    Method: start from n uniformly random centres inside the polygon shrunk by the disc radius, then minimise the overlap energy E = sum over pairs closer than d_over
    of (d_over - distance)^2 plus the squared distance by which any centre lies outside the shrunk polygon, with L-BFGS (analytic gradient: the pair forces from a k-d tree, the wall forces from
    the nearest point of the shrunk polygon). E reaches 0 when a packing exists; at the packing fractions of a winding (0.6 of the layer; jamming is above 0.8) it does within some
    tens of iterations, in 0.01 s. Raises RuntimeError if E stays above 1e-10 mm^2 after `tries` random starts (the fraction asked is too high)."""
    import shapely
    from scipy.optimize import minimize
    from scipy.spatial import cKDTree
    inner = poly.buffer(-d_over / 2.0)
    if inner.is_empty:
        raise RuntimeError("the layer is thinner than one strand")
    minx, miny, maxx, maxy = inner.bounds

    def energy(v):
        xy = v.reshape(-1, 2)
        grad = np.zeros_like(xy)
        e = 0.0
        pairs = cKDTree(xy).query_pairs(d_over, output_type="ndarray")
        if len(pairs):
            i, j = pairs[:, 0], pairs[:, 1]
            dv = xy[i] - xy[j]
            dist = np.maximum(np.hypot(dv[:, 0], dv[:, 1]), 1e-12)
            over = d_over - dist
            e += float(np.sum(over ** 2))
            force = (-2.0 * over / dist)[:, None] * dv
            np.add.at(grad, i, force)
            np.add.at(grad, j, -force)
        outside = np.flatnonzero(~shapely.contains_xy(inner, xy[:, 0], xy[:, 1]))
        if len(outside):
            line = shapely.shortest_line(shapely.points(xy[outside]), inner)
            nearest = shapely.get_coordinates(line).reshape(-1, 2, 2)[:, 1]
            dv = xy[outside] - nearest
            e += float(np.sum(dv ** 2))
            grad[outside] += 2.0 * dv
        return e, grad.ravel()

    from scipy.optimize import minimize
    for _ in range(tries):
        pts: list = []
        while len(pts) < n:
            c = rng.uniform([minx, miny], [maxx, maxy], size=(4 * n, 2))
            pts.extend(c[shapely.contains_xy(inner, c[:, 0], c[:, 1])].tolist())
        res = minimize(energy, np.array(pts[:n]).ravel(), jac=True, method="L-BFGS-B", options={"maxiter": max_iter, "ftol": 1e-30, "gtol": 1e-13})
        if res.fun < 1e-10:
            return res.x.reshape(-1, 2)
    raise RuntimeError(f"cannot pack {n} discs of diameter {d_over} mm into the layer (packing fraction {n * math.pi * d_over ** 2 / 4 / poly.area:.3f})")


def strand_layout(p: mg.Params, n_per_layer: int, d_over_mm: float, seed: int) -> dict:
    """The strands of both layers: {'xy' (n, 2) mm, 'layer' (n,) 0 / 1, 'n_per_layer'}: layer 0's strands first, then layer 1's, each layer packed with its own random draw
    of the seed (so a layout is repeatable and the seed is the only thing that changes between two layouts of the same winding)."""
    rng = np.random.default_rng(seed)
    xy, layer = [], []
    for li, poly in enumerate(layer_polygons(p)):
        c = pack_disks(poly, n_per_layer, d_over_mm, rng)
        xy.append(c)
        layer.append(np.full(len(c), li))
    return {"xy": np.vstack(xy), "layer": np.concatenate(layer), "n_per_layer": int(n_per_layer)}


ARRANGEMENTS = ("random", "line", "cluster", "row")


def assign_turns(xy: np.ndarray, layer: np.ndarray, n_turns: int, n_s: int, arrangement: str, seed: int) -> np.ndarray:
    """Which turn each strand belongs to, within its layer: an int array (n,), the turns of the two layers numbered 0 ... n_turns - 1 each. Every turn has exactly n_s strands.
    Strands of one turn are the n_s wires wound together; where they sit decides what their parallel connection does (circulating currents):
      'random'   a random partition of the layer's strands into turns (the n_s strands of a turn anywhere in the layer);
      'line'     the WORST ordered arrangement: sort the layer's strands by depth (y), give turn t the strands of ranks t, t + n_turns, t + 2 n_turns, ...: every turn has one
                 strand in each n_turns-th of the depth, a line of strands from the slot bottom to the opening;
      'cluster'  compact blobs: sorted by depth, cut into bands of ceil(sqrt(n_s)) turns, each band sorted across the slot and cut into turns of n_s neighbours;
      'row'      the strands of a turn side by side at one depth: sorted by depth, consecutive n_s strands a turn (the best case for a field that points across the slot)."""
    rng = np.random.default_rng(seed)
    turn = np.full(len(xy), -1, dtype=int)
    for li in range(LAYERS):
        idx = np.flatnonzero(layer == li)
        if len(idx) != n_turns * n_s:
            raise ValueError(f"layer {li} has {len(idx)} strands, not {n_turns} x {n_s}")
        if arrangement == "random":
            order = rng.permutation(idx)
            for t in range(n_turns):
                turn[order[t * n_s:(t + 1) * n_s]] = t
        elif arrangement == "line":
            order = idx[np.argsort(xy[idx, 1], kind="stable")]
            for t in range(n_turns):
                turn[order[t::n_turns]] = t
        elif arrangement == "row":
            order = idx[np.argsort(xy[idx, 1], kind="stable")]
            for t in range(n_turns):
                turn[order[t * n_s:(t + 1) * n_s]] = t
        elif arrangement == "cluster":
            order = idx[np.argsort(xy[idx, 1], kind="stable")]
            band = max(1, math.ceil(math.sqrt(n_s))) * n_s
            t = 0
            for s0 in range(0, len(order), band):
                chunk = order[s0:s0 + band]
                chunk = chunk[np.argsort(xy[chunk, 0], kind="stable")]
                for c0 in range(0, len(chunk), n_s):
                    turn[chunk[c0:c0 + n_s]] = t
                    t += 1
            if t != n_turns:                                          # the last band is shorter: its turns still have n_s strands, so the count holds
                raise RuntimeError(f"cluster arrangement made {t} turns, not {n_turns}")
        else:
            raise ValueError(f"arrangement {arrangement!r}: one of {ARRANGEMENTS}")
    return turn


# ---------------------------------------------------------------------------------------------- the FEMM model (runs in a pool worker)
def _arc_angle_ccw(a, b, c) -> float:
    return math.degrees(math.atan2(b[1] - c[1], b[0] - c[0]) - math.atan2(a[1] - c[1], a[0] - c[0])) % 360.0


def _draw_cell(femm, p: mg.Params, o: dict, *, freq: float, rho: float | None, layered: bool, drive_field_T: float | None) -> list:
    """Draw the cell (steel, the air gap, the slot or slots) and its boundaries in the open FEMM document; returns the list of the slots' angles. The cell is o['n_pitch']
    slot pitches wide (1: one slot, tooth centre to tooth centre; 3, 5: that many slot outlines, the neighbours empty), centred on the slot at SLOT_PHI_DEG. `layered`: the
    central slot is cut into its two coil layers by lines (S1-S2 and M1-M2: the homogenised model's blocks); otherwise the slot, its opening and the gap are ONE air region
    (the strands are drawn inside it by the caller). Materials: Air, Cu (if rho is given: conductivity 1 / rho), Fe (the steel of o['steel']), Coil (mu_r 1, no conductivity).
    The steel is laminated only in a time-harmonic problem (freq > 0)."""
    steel = o["steel"]
    ac = freq > 0
    femm.newdocument(0)
    femm.mi_smartmesh(0)
    femm.mi_probdef(freq, "millimeters", "planar", o["precision"], p.stack, 30)
    femm.mi_getmaterial("Air")
    if rho is not None:
        femm.mi_addmaterial("Cu", 1, 1, 0, 0, 1.0 / rho / 1e6, 0, 0, 1, 0, 0, 0, 0, 0)
    femm.mi_addmaterial("Coil", 1, 1, 0)
    femm.mi_addmaterial("Fe", steel["mu_r"], steel["mu_r"], 0, 0, steel["sigma_MS"] if ac else 0.0, steel["d_lam_mm"] if ac else 0.0, 0,
                        steel["lam_fill"] if ac else 1.0, 0, 0, 0, 0, 0)
    femm.mi_addboundprop("A0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    if drive_field_T is not None:
        femm.mi_addboundprop("UNI", 0, -drive_field_T * 1e-3, 0, 0, 0, 0, 0, 0, 0, 0, 0)     # A = A1 x with x in mm: B_y = -A1 * 1e3 tesla
    n_pitch = int(o["n_pitch"])
    half = n_pitch * HALF_PITCH_DEG
    pitch = 2.0 * HALF_PITCH_DEG
    phi0 = SLOT_PHI_DEG
    phis = [phi0 + pitch * (j - (n_pitch - 1) // 2) for j in range(n_pitch)]
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    seg = lambda a, b: femm.mi_addsegment(*a, *b)
    cut_l = [mg.pol(r, phi0 - half) for r in (p.r_ro, p.r_si, p.r_so)]                         # the cuts through the tooth centres: rotor surface, bore, outer radius
    cut_r = [mg.pol(r, phi0 + half) for r in (p.r_ro, p.r_si, p.r_so)]
    for q in (*cut_l, *cut_r):
        femm.mi_addnode(*q)
    for phi in phis:
        pts = slot_points(p, phi)
        O1, O2, O3, O4, S1, S2, W1, W2, B, cb = (pts[k] for k in ("O1", "O2", "O3", "O4", "S1", "S2", "W1", "W2", "B", "cb"))
        for q in (O1, O2, O3, O4, S1, S2, W1, W2, B):
            femm.mi_addnode(*q)
        if layered and phi == phi0:
            regs = slot_regions(p)
            m1, m2 = regs[1].edges[1].a, regs[1].edges[1].b
            femm.mi_addnode(*m1)
            femm.mi_addnode(*m2)
            seg(O1, O2); seg(O2, S1); seg(S1, m1); seg(m1, W1); seg(W2, m2); seg(m2, S2); seg(S2, O3); seg(O3, O4)
            seg(S1, S2); seg(m1, m2)
        else:
            seg(O1, O2); seg(O2, S1); seg(S1, W1); seg(W2, S2); seg(S2, O3); seg(O3, O4)     # the outline only: the air and the strands share one region
        femm.mi_addarc(*W1, *B, _arc_angle_ccw(W1, B, cb), 5.0)                                # the round bottom: two counter-clockwise arcs about cb, each under 180 deg
        femm.mi_addarc(*B, *W2, _arc_angle_ccw(B, W2, cb), 5.0)

    def arc_between(r, a_deg, b_deg, maxseg):
        femm.mi_addarc(*mg.pol(r, a_deg), *mg.pol(r, b_deg), b_deg - a_deg, maxseg)
    bore = [phi0 - half] + [a for phi in phis for a in (phi - a_o, phi + a_o)] + [phi0 + half]      # the bore's arcs: cut ... opening | opening ... opening | opening ... cut
    for k in range(0, len(bore), 2):
        arc_between(p.r_si, bore[k], bore[k + 1], 1.0)
    arc_between(p.r_so, phi0 - half, phi0 + half, 1.0)                                         # the outer radius
    arc_between(p.r_ro, phi0 - half, phi0 + half, 1.0)                                         # the rotor's surface
    seg(cut_l[0], cut_l[1]); seg(cut_l[1], cut_l[2]); seg(cut_r[0], cut_r[1]); seg(cut_r[1], cut_r[2])
    uni = drive_field_T is not None

    def set_arc(r, prop):
        femm.mi_selectarcsegment(*mg.pol(r, phi0))
        femm.mi_setarcsegmentprop(1.0, prop, 0, 0)
        femm.mi_clearselected()

    def set_seg(a, b, prop):
        femm.mi_selectsegment((a[0] + b[0]) / 2, (a[1] + b[1]) / 2)
        femm.mi_setsegmentprop(prop, 0, 1, 0, 0)
        femm.mi_clearselected()
    set_arc(p.r_so, "UNI" if uni else "A0")
    if uni:
        set_arc(p.r_ro, "UNI")
    elif o["rotor"] == "dirichlet":
        set_arc(p.r_ro, "A0")
    if uni or o["lateral"] == "dirichlet":
        for cut in (cut_l, cut_r):
            set_seg(cut[0], cut[1], "UNI" if uni else "A0")
            set_seg(cut[1], cut[2], "UNI" if uni else "A0")
    for pt, mat, mesh in ((mg.pol(p.r_so - 1.5, phi0), "Fe", o["mesh_fe_mm"]), (mg.pol((p.r_ro + p.r_si) / 2, phi0), "Air", o["mesh_air_mm"])):
        femm.mi_addblocklabel(*pt)
        femm.mi_selectlabel(*pt)
        femm.mi_setblockprop(mat, 0, mesh, "", 0, 0, 0)
        femm.mi_clearselected()
    return phis


def build_model(femm, p: mg.Params, fem: str, *, freq: float, rho: float, d_mm: float, xy: np.ndarray, circuit: list[str], circuits: dict, options: dict,
                drive_field_T: float | None = None) -> None:
    """Draw the slot cell and the strands in the FEMM of this worker and save it as `fem`.
    freq      Hz (0: magnetostatic)
    rho       the copper's resistivity [Ohm m]
    d_mm      the strands' bare diameter
    xy        (n, 2) the strands' centres [mm]
    circuit   the circuit name of each strand; circuits: {name: (current (complex, peak A), 1 series | 0 parallel)}: every name used must be there
    drive_field_T  if given, a uniform field B (peak, along y) is applied: A = -B x on the outer boundary of the cell (the outer arc, the two cuts, the rotor surface), the
              'rotor field' drive of the strand study (the strands then usually carry circuits of 0 A)
    The strands' block labels are the points xy themselves, which is how the results are read back."""
    o = {**DEFAULT_OPTIONS, **options}
    _draw_cell(femm, p, o, freq=freq, rho=rho, layered=False, drive_field_T=drive_field_T)
    for name, (current, kind) in circuits.items():
        femm.mi_addcircprop(name, current, kind)
    r = d_mm / 2.0
    if o["area_correct"]:
        n_seg = 2 * math.ceil(180.0 / o["arc_deg"])
        r *= math.sqrt(2.0 * math.pi / (n_seg * math.sin(2.0 * math.pi / n_seg)))
    for (x, y), name in zip(xy, circuit):
        femm.mi_addnode(x + r, y)
        femm.mi_addnode(x - r, y)
        femm.mi_addarc(x + r, y, x - r, y, 180, o["arc_deg"])
        femm.mi_addarc(x - r, y, x + r, y, 180, o["arc_deg"])
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop("Cu", 0, d_mm * o["mesh_cu_over_d"], name, 0, 0, 1)
        femm.mi_clearselected()
    femm.mi_saveas(fem)


def build_homogenised(femm, p: mg.Params, fem: str, *, layer_currents: tuple, options: dict) -> None:
    """The same cell with the strands replaced by two coil blocks of uniform current density (the first study's model of the coil): layer_currents are the two layers' TOTAL
    currents (A, series circuits, turns 1). Magnetostatic only (the blocks have no conductivity). The air above the opening and the two layers are separate regions here."""
    o = {**DEFAULT_OPTIONS, **options}
    _draw_cell(femm, p, o, freq=0.0, rho=None, layered=True, drive_field_T=None)
    for i, cur in enumerate(layer_currents):
        femm.mi_addcircprop(f"L{i}", cur, 1)
    for i, reg in enumerate(slot_regions(p)[1:]):
        femm.mi_addblocklabel(*reg.label)
        femm.mi_selectlabel(*reg.label)
        femm.mi_setblockprop("Coil", 0, 0.1, f"L{i}", 0, 0, 1)
        femm.mi_clearselected()
    femm.mi_saveas(fem)


# ---------------------------------------------------------------------------------------------- reading what FEMM wrote (the .ans)
def read_ans(path: str) -> dict:
    """The mesh and the solution of an .ans (magnetostatic or time-harmonic): {'xy' (n, 2) mm, 'A' (n,) complex [Wb/m], 'tri' (m, 3), 'label' (m,), 'labels' ([x, y] of
    every block label in file order), 'freq'}. First-order triangles: A is linear in a triangle and B constant. The label of a triangle indexes the list of block labels."""
    with open(path, encoding="latin-1") as handle:
        text = handle.read()
    head, _, body = text.partition("[Solution]")
    freq = float(next(ln.split("=")[1] for ln in head.splitlines() if ln.startswith("[Frequency]")))
    lines = body.split("\n")
    nn = int(lines[1])
    cols = 4 if freq > 0 else 3
    nodes = np.array([ln.split()[:cols] for ln in lines[2:2 + nn]], dtype=float)
    ne = int(lines[2 + nn])
    el = np.array([ln.split()[:4] for ln in lines[3 + nn:3 + nn + ne]], dtype=int)
    labels, in_block = [], False
    for ln in head.splitlines():
        if ln.startswith("[NumBlockLabels]"):
            n_lab = int(ln.split("=")[1])
            in_block = True
            labels = []
            continue
        if in_block:
            if len(labels) < n_lab:
                parts = ln.split()
                labels.append((float(parts[0]), float(parts[1])))
            else:
                break
    A = nodes[:, 2] + (1j * nodes[:, 3] if freq > 0 else 0.0)
    return {"xy": nodes[:, :2], "A": A, "tri": el[:, :3], "label": el[:, 3], "labels": labels, "freq": freq}


def strand_triangles(sol: dict, xy: np.ndarray, tol: float = 1e-4) -> np.ndarray:
    """For every triangle of the solution the index of the strand it belongs to (-1 if none): the triangle's block label point is the strand's centre (found with a k-d tree:
    a label within tol mm of a strand's centre is that strand's; the two airs and the steel are far from every centre)."""
    from scipy.spatial import cKDTree
    dist, nearest = cKDTree(xy).query(np.array(sol["labels"], dtype=float))
    lab_to_strand = np.where(dist < tol, nearest, -1)
    if len(set(lab_to_strand[lab_to_strand >= 0].tolist())) != len(xy):
        raise RuntimeError("the labels of the solution do not match the strands one to one")
    return lab_to_strand[sol["label"]]


def triangle_geometry(sol: dict) -> dict:
    """Area [mm^2], centroid [mm] and the nodal values A per triangle, and B = (dA/dy, -dA/dx) [T] (complex when A is) of every triangle."""
    x, y = sol["xy"][sol["tri"], 0], sol["xy"][sol["tri"], 1]
    a = sol["A"][sol["tri"]]
    det = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0])
    dadx = ((a[:, 1] - a[:, 0]) * (y[:, 2] - y[:, 0]) - (a[:, 2] - a[:, 0]) * (y[:, 1] - y[:, 0])) / det
    dady = ((a[:, 2] - a[:, 0]) * (x[:, 1] - x[:, 0]) - (a[:, 1] - a[:, 0]) * (x[:, 2] - x[:, 0])) / det
    return {"area_mm2": 0.5 * np.abs(det), "cx": x.mean(axis=1), "cy": y.mean(axis=1), "A": a, "Bx": 1e3 * dady, "By": -1e3 * dadx}


def strand_field_stats(sol: dict, xy: np.ndarray) -> dict:
    """Per strand, from a MAGNETOSTATIC (real) solution: the area [m^2], the mean of A over the strand's cross section [Wb/m], the mean of B (Bx, By) [T], and the integral over the
    strand of (A - A_mean)^2 [Wb^2/m^2 x m^2] (the low-frequency eddy loss of the strand is sigma w^2 L / 2 times it). Exact for first-order triangles."""
    tri_strand = strand_triangles(sol, xy)
    g = triangle_geometry(sol)
    ok = tri_strand >= 0
    k = tri_strand[ok]
    n = len(xy)
    area = g["area_mm2"][ok] * 1e-6
    a = np.real(g["A"][ok])
    s = np.bincount(k, weights=area, minlength=n)
    a_mean_tri = a.mean(axis=1)
    a_int = np.bincount(k, weights=area * a_mean_tri, minlength=n)
    a_sq = np.bincount(k, weights=area / 12.0 * ((a ** 2).sum(axis=1) + a.sum(axis=1) ** 2), minlength=n)    # integral of A^2 over each triangle: T/12 (sum A_i^2 + (sum A_i)^2)
    a_bar = a_int / s
    return {"area_m2": s, "A_mean": a_bar, "Bx_mean": np.bincount(k, weights=area * np.real(g["Bx"][ok]), minlength=n) / s,
            "By_mean": np.bincount(k, weights=area * np.real(g["By"][ok]), minlength=n) / s, "A_dev_integral": a_sq - s * a_bar ** 2,
            "n_triangles": np.bincount(k, minlength=n)}


def strand_cross_integrals(sol_a: dict, sol_b: dict, xy: np.ndarray) -> np.ndarray:
    """integral over each strand of (A_a - mean A_a)(A_b - mean A_b) for two real solutions on the SAME mesh (the same drawing: the mesh only differs when the geometry does).
    With the two solutions the layers' unit-current fields, the loss of any complex pair of layer currents is made from the three integrals (a, a), (b, b), (a, b)."""
    tri_strand = strand_triangles(sol_a, xy)
    ga, gb = triangle_geometry(sol_a), triangle_geometry(sol_b)
    ok = tri_strand >= 0
    k = tri_strand[ok]
    area = ga["area_mm2"][ok] * 1e-6
    fa, fb = np.real(ga["A"][ok]), np.real(gb["A"][ok])
    n = len(xy)
    s = np.bincount(k, weights=area, minlength=n)
    ma = np.bincount(k, weights=area * fa.mean(axis=1), minlength=n) / s
    mb = np.bincount(k, weights=area * fb.mean(axis=1), minlength=n) / s
    ab = np.bincount(k, weights=area / 12.0 * ((fa * fb).sum(axis=1) + fa.sum(axis=1) * fb.sum(axis=1)), minlength=n)
    return ab - s * ma * mb


# ---------------------------------------------------------------------------------------------- the tasks (plain dicts) and the one function the pool maps
def _wire_rho(task: dict) -> float:
    return float(ac_theory.resistivity(task.get("t_c", 100.0)))


def task_ac(tag: str, fem_dir: str, overrides: dict, *, n_turns: int, n_s: int, d_mm: float, d_over_mm: float, f: float, i_layers: tuple, seed: int = 1,
            arrangement: str = "random", circuit: str = "series", t_c: float = 100.0, options: dict | None = None, drive_field_T: float | None = None,
            want_field: bool = False) -> dict:
    """A time-harmonic task: the strand-by-strand loss at frequency f.
    i_layers  the two layers' STRAND currents [A, peak, complex]: every strand of layer 0 carries i_layers[0], every strand of layer 1 i_layers[1] ('series' circuits: equal
              sharing). With circuit='parallel' the strands of one turn are ONE circuit carrying n_s times that (the turn's current) and share it as their positions make them.
    arrangement  how strands are given to turns (assign_turns): used by circuit='parallel' only
    drive_field_T  an external uniform field instead of the strands' own currents (then pass i_layers (0, 0))
    want_field  also return, for every strand triangle, the current density (complex, A/m^2 peak) for the pictures"""
    return {"kind": "ac", "tag": tag, "fem_dir": fem_dir, "overrides": overrides, "n_turns": n_turns, "n_s": n_s, "d_mm": d_mm, "d_over_mm": d_over_mm, "f": f,
            "i_layers": tuple(i_layers), "seed": seed, "arrangement": arrangement, "circuit": circuit, "t_c": t_c, "options": options or {}, "drive_field_T": drive_field_T,
            "want_field": want_field}


def task_ms(tag: str, fem_dir: str, overrides: dict, *, n_turns: int, n_s: int, d_mm: float, d_over_mm: float, seed: int = 1, t_c: float = 100.0,
            options: dict | None = None, drive_field_T: float | None = None) -> dict:
    """A magnetostatic task on the same strand layout: two solves, one with every strand of layer 0 at 1 A (layer 1 at 0) and one the other way round (or, with
    drive_field_T, one solve in the uniform field and no strand currents); per strand the area, the mean A, the mean B and the integrals that make the low-frequency loss
    (strand_field_stats, strand_cross_integrals)."""
    return {"kind": "ms", "tag": tag, "fem_dir": fem_dir, "overrides": overrides, "n_turns": n_turns, "n_s": n_s, "d_mm": d_mm, "d_over_mm": d_over_mm, "seed": seed,
            "t_c": t_c, "options": options or {}, "drive_field_T": drive_field_T}


def task_hom(tag: str, fem_dir: str, overrides: dict, *, n_turns: int, n_s: int, d_mm: float, d_over_mm: float, seed: int = 1, options: dict | None = None) -> dict:
    """The homogenised slot (two coil blocks of uniform current density, 1 A a strand: n_turns x n_s A in each layer's circuit): B at the centre of every strand of the layout
    of this seed, for layer 0 alone and for layer 1 alone (magnetostatic)."""
    return {"kind": "hom", "tag": tag, "fem_dir": fem_dir, "overrides": overrides, "n_turns": n_turns, "n_s": n_s, "d_mm": d_mm, "d_over_mm": d_over_mm, "seed": seed,
            "options": options or {}}


def _paths(task: dict) -> tuple[str, str]:
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
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


def _layout(task: dict, p: mg.Params) -> dict:
    n_layer = task["n_turns"] * task["n_s"]
    lay = strand_layout(p, n_layer, task["d_over_mm"], task["seed"])
    lay["turn"] = assign_turns(lay["xy"], lay["layer"], task["n_turns"], task["n_s"], task.get("arrangement", "random"), task["seed"] + 1000)
    return lay


def _run_ac(task: dict) -> dict:
    import femm
    started = time.perf_counter()
    fem, ans = _paths(task)
    p = design.params(task["overrides"])
    lay = _layout(task, p)
    xy, layer, turn = lay["xy"], lay["layer"], lay["turn"]
    n = len(xy)
    rho = _wire_rho(task)
    i_layers = task["i_layers"]
    if task["circuit"] == "series":
        names = [f"s{k}" for k in range(n)]
        circuits = {names[k]: (complex(i_layers[layer[k]]), 1) for k in range(n)}
    elif task["circuit"] == "parallel":
        names = [f"t{layer[k]}_{turn[k]}" for k in range(n)]
        circuits = {f"t{li}_{t}": (complex(i_layers[li]) * task["n_s"], 0) for li in range(LAYERS) for t in range(task["n_turns"])}
    else:
        raise ValueError(f"circuit {task['circuit']!r}: 'series' or 'parallel'")
    t_draw = time.perf_counter()
    build_model(femm, p, fem, freq=task["f"], rho=rho, d_mm=task["d_mm"], xy=xy, circuit=names, circuits=circuits, options=task["options"], drive_field_T=task["drive_field_T"])
    femm.mi_close()
    femm.opendocument(fem)
    t_solve = time.perf_counter()
    femm.mi_analyze(1)
    femm.mi_loadsolution()
    t_read = time.perf_counter()
    loss, cur, area = np.empty(n), np.empty(n, dtype=complex), np.empty(n)
    for k, (x, y) in enumerate(xy):
        femm.mo_selectblock(float(x), float(y))
        loss[k] = float(np.real(femm.mo_blockintegral(4)))
        cur[k] = complex(femm.mo_blockintegral(7))
        area[k] = float(np.real(femm.mo_blockintegral(5)))
        femm.mo_clearblock()
    nodes, elements = int(femm.mo_numnodes()), int(femm.mo_numelements())
    femm.mo_close()
    out = {"status": "ok", "kind": "ac", "tag": task["tag"], "f": task["f"], "xy": xy, "layer": layer, "turn": turn, "loss_W": loss, "current_A": cur, "area_m2": area,
           "rho": rho, "nodes": nodes, "elements": elements, "stack_mm": p.stack,
           "seconds": {"draw": t_solve - t_draw, "solve": t_read - t_solve, "read": time.perf_counter() - t_read, "total": time.perf_counter() - started}}
    if task["want_field"]:
        sol = read_ans(ans)
        g = triangle_geometry(sol)
        ts = strand_triangles(sol, xy)
        ok = ts >= 0
        k = ts[ok]
        a_tri = g["A"][ok]
        w = 2.0 * math.pi * task["f"]
        s_k = area
        a_mean = np.bincount(k, weights=g["area_mm2"][ok] * 1e-6 * np.real(a_tri.mean(axis=1)), minlength=n) + 1j * np.bincount(
            k, weights=g["area_mm2"][ok] * 1e-6 * np.imag(a_tri.mean(axis=1)), minlength=n)
        a_mean = a_mean / np.bincount(k, weights=g["area_mm2"][ok] * 1e-6, minlength=n)
        sigma = 1.0 / rho
        j_tri = cur[k][:, None] / s_k[k][:, None] - 1j * w * sigma * (a_tri - a_mean[k][:, None])        # J = I_k / S_k - j w sigma (A - A_mean): the nodal values, per triangle
        out["field"] = {"cx": g["cx"][ok], "cy": g["cy"][ok], "area_mm2": g["area_mm2"][ok], "strand": k, "J_nodes": j_tri,
                        "tri_xy": sol["xy"][sol["tri"][ok]], "A_mean": a_mean}
    return out


def _run_ms(task: dict) -> dict:
    import femm
    started = time.perf_counter()
    fem, ans = _paths(task)
    p = design.params(task["overrides"])
    lay = _layout(task, p)
    xy, layer = lay["xy"], lay["layer"]
    n = len(xy)
    rho = _wire_rho(task)
    names = [f"s{k}" for k in range(n)]
    drive = task["drive_field_T"]
    patterns = [(0.0, 0.0)] if drive is not None else [(1.0, 0.0), (0.0, 1.0)]
    sols = []
    for pat in patterns:
        circuits = {names[k]: (complex(pat[layer[k]]), 1) for k in range(n)}
        if not sols:
            build_model(femm, p, fem, freq=0.0, rho=rho, d_mm=task["d_mm"], xy=xy, circuit=names, circuits=circuits, options=task["options"], drive_field_T=drive)
            femm.mi_close()
            femm.opendocument(fem)
        else:
            for k in range(n):
                femm.mi_modifycircprop(names[k], 1, complex(pat[layer[k]]))
        femm.mi_analyze(1)
        sols.append(read_ans(ans))
    nodes = int(len(sols[0]["xy"]))
    femm.mi_close()
    stats = [strand_field_stats(s, xy) for s in sols]
    out = {"status": "ok", "kind": "ms", "tag": task["tag"], "xy": xy, "layer": layer, "area_m2": stats[0]["area_m2"], "rho": rho, "nodes": nodes, "stack_mm": p.stack,
           "A_mean": [s["A_mean"] for s in stats], "Bx_mean": [s["Bx_mean"] for s in stats], "By_mean": [s["By_mean"] for s in stats],
           "A_dev_integral": [s["A_dev_integral"] for s in stats]}
    if len(sols) == 2:
        if not (np.array_equal(sols[0]["tri"], sols[1]["tri"]) and np.array_equal(sols[0]["xy"], sols[1]["xy"])):
            raise RuntimeError("the two magnetostatic solves of a strand layout were meshed differently: their fields cannot be combined element by element")
        out["A_dev_cross"] = strand_cross_integrals(sols[0], sols[1], xy)
    out["seconds"] = time.perf_counter() - started
    return out


def _run_hom(task: dict) -> dict:
    import femm
    started = time.perf_counter()
    fem, ans = _paths(task)
    p = design.params(task["overrides"])
    lay = _layout(task, p)
    xy, layer = lay["xy"], lay["layer"]
    n_layer = task["n_turns"] * task["n_s"]
    bx, by = [], []
    build_homogenised(femm, p, fem, layer_currents=(float(n_layer), 0.0), options=task["options"])
    femm.mi_close()
    femm.opendocument(fem)
    for cur in ((float(n_layer), 0.0), (0.0, float(n_layer))):
        femm.mi_modifycircprop("L0", 1, cur[0])
        femm.mi_modifycircprop("L1", 1, cur[1])
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        b = np.array([np.real(femm.mo_getb(float(x), float(y))) for x, y in xy], dtype=float)
        femm.mo_close()
        bx.append(b[:, 0])
        by.append(b[:, 1])
    femm.mi_close()
    return {"status": "ok", "kind": "hom", "tag": task["tag"], "xy": xy, "layer": layer, "Bx": bx, "By": by, "seconds": time.perf_counter() - started}


def run_task(task: dict) -> dict:
    """The one function the pool maps: runs a task of kind 'ac', 'ms' or 'hom' on this worker's FEMM. An error inside is returned as {'status': 'error', ...}, never raised: one bad
    task must not end a run. The task's model files are removed whatever happens."""
    import femm
    fem, ans = _paths(task)
    started = time.perf_counter()
    try:
        return {"ac": _run_ac, "ms": _run_ms, "hom": _run_hom}[task["kind"]](task)
    except Exception as error:
        return {"status": "error", "kind": task.get("kind"), "tag": task.get("tag"), "error": f"{type(error).__name__}: {error}", "trace": traceback.format_exc(limit=5),
                "seconds": time.perf_counter() - started}
    finally:
        _cleanup(femm, fem, ans)


# ---------------------------------------------------------------------------------------------- numbers from the results (no FEMM)
def dc_loss_W(res: dict, i_layers: tuple) -> np.ndarray:
    """Each strand's DC loss [W] at the amplitudes of i_layers: |I|^2 rho L / (2 A_k), A_k FEMM's own area of that strand's polygon. For a parallel circuit pass the strand current
    that equal sharing would give (the turn's current over n_s)."""
    i = np.array([abs(complex(i_layers[li])) for li in res["layer"]])
    return i ** 2 * res["rho"] * res["stack_mm"] * 1e-3 / (2.0 * res["area_m2"])


def r_ratio(res: dict, i_layers: tuple) -> float:
    """R_ac / R_dc of the slot: the total loss of all strands over the DC loss of the same strand currents. For a complex pair of layer currents the DC loss uses the amplitudes."""
    return float(res["loss_W"].sum() / dc_loss_W(res, i_layers).sum())
