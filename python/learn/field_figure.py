"""The field figure: the band-free motor at work, drawn the way the Rotor Flux Explorer draws it.

What it makes, and why
    out/learn/field_figure/field_figure.png, one picture 1400 px wide: the 90 degree sector of the band-free motor,
    Params(barrier="steel") (the rotor with its air band filled with steel), at one working point. The working point is 10 A
    peak at the current angle that learn/mtpa found best at 10 A, 22.9 degrees from +q toward -d (learn/dq_torque/common.py
    has the conventions), with the rotor at theta_e = 0: the d axis of its north pole on phase A's axis. It is meant as the
    opening picture of the motor project's README, so it comes from a fresh FEMM solve made with the project's own code, and it
    takes the Explorer's look: |B| from near-white to deep blue on a fixed 0 ... 2.2 T scale, flux lines in thin dark ink, every
    slot layer tinted by the sign of its current and labelled with its phase, magnets outlined in red with an arrow to their north.
    The sector lies on its side in the first quadrant, as in the Gmsh picture of the tools note: the centre at the bottom left,
    one cut edge along the bottom, the other up the left side, the rim at the top right.

Run, from the code folder (about 6 seconds, of which the solves take 3.5):
    uv run python learn/field_figure.py                  solve, check, draw
    uv run python learn/field_figure.py --render-only    draw again from the saved solution, no FEMM: to try another look
    options: --theta-e DEG (rotor position, default 0), --lines N (about this many flux lines, default 36)
It reads learn/mtpa's results (out/learn/mtpa/summary.json) for the angle and for the torque it checks against; without
them it takes 22.9 degrees and skips that check. Writes into out/learn/field_figure/ (git ignores out): field_figure.png, the
picture; field_figure.ans, the solution it draws; field_figure.json, the facts of the run; field_figure_work.fem and .ans, FEMM's
working files. Their names are this script's own, because two FEMMs must never share a file name.

The method, in the order of the code
 1. The solve. One FEMM of its own, below normal priority, under the machine's solver slot (femm_queue), on the fine mesh
    that learn/mtpa used. motorgeo draws the model and FemmModel.solve puts the three phase currents in and keeps the .ans.
    The other five rotor positions of the 60 degree ripple period are solved too, only for a check: the mean torque of the six
    is the number learn/mtpa reports. A watchdog (faulthandler) ends the run after 300 s if a FEMM call hangs; the FEMM goes with
    the process, because a headless FEMM quits when its caller hangs up, and the operating system frees the slot (checked: a
    3 s timer on a deliberately hung run left no FEMM and a free slot).
 2. The three phase currents. A current of peak I at angle gamma from +q toward -d has i_d = -I sin(gamma), i_q = I cos(gamma).
    The inverse Park transform turns that pair into the phase currents for a rotor at theta_e:
        i_k = i_d cos(th_k) - i_q sin(th_k),   th_k = theta_e, theta_e - 120, theta_e + 120 degrees for k = a, b, c.
    Written out, i_k = -I sin(th_k + gamma): three sines 120 degrees apart, peak I, and gamma slides all three along together.
    Here (10 A, 22.909 degrees, theta_e = 0): i_a = -3.893 A, i_b = +9.924 A, i_c = -6.031 A, which add up to zero.
 3. |B| = |curl A| on each triangle. FEMM solves for the magnetic vector potential A. In this planar model the current and the
    A it makes point along z only, so A is one number A(x, y) per point and curl A is B = (dA/dy, -dA/dx), with |B| = |grad A|.
    The elements are first-order triangles: A is stored at the three corners and varies linearly between them, so its gradient is
    one constant per triangle, and so is B. With corners (x_i, y_i) and A_i, and det = (x1 - x0)(y2 - y0) - (x2 - x0)(y1 - y0),
        dA/dx = [(A1 - A0)(y2 - y0) - (A2 - A0)(y1 - y0)] / det,    dA/dy = [(A2 - A0)(x1 - x0) - (A1 - A0)(x2 - x0)] / det,
    and x and y in mm make B come out in T once multiplied by 1000 (femfile.element_b does it). Every triangle gets one colour from
    its |B|. FEMM's own post-processor smooths B between triangles unless told not to; with mo_smooth("off") its B at the middle
    of a triangle is the triangle's, and the script checks that the two agree (they do, to 1e-14 T).
 4. Flux lines are contours of A. B = curl(A z) is at right angles to grad A, so B runs along the lines where A is constant:
    those are the field lines. Equal spacing of the levels means equal flux between neighbours: the flux through a strip L deep
    that runs from a point P to a point Q is L (A(Q) - A(P)), whatever path the strip takes (Stokes' theorem; only the two edges of
    the strip, where A points along them, count). Lines at A = (k + 1/2) s for whole numbers k therefore enclose L s of flux between
    each pair of neighbours, anywhere along them. Here s = 5e-4 Wb/m and L = 20 mm, so every gap carries 10 microwebers. Where the
    lines crowd, the field is strong: |B| = s / spacing. (The levels sit half a step off zero, so that no line traces the rim and the
    shaft, where A = 0 is the boundary condition, and so that the gap across A = 0 is no wider than the others.)
 5. Where the rotor really is. The .ans keeps the rotor where it was drawn, whatever its angle: FEMM turns the rotor only inside
    the air-gap element. The rotor's nodes are turned here, counter-clockwise by
        rotor_deg = theta_e / 4 + (phase A's axis, 65.03 degrees) - (pole 0's drawn d axis, 22.5 degrees) = 42.53 degrees,
    and three copies of everything, 90 degrees apart, are laid side by side, because the sector's boundaries are periodic (the field
    repeats every 90 degrees, with the same A); the window cuts the copies to the sector. The thin ring in the middle of the air gap
    is not meshed (FEMM solves it by formula), so its two edges of nodes are moved to mid-gap to close it. The test that the rotor
    is turned the right way (turn_check): across that ring A hardly changes, so A at the rotor's edge nodes, turned, must equal A
    on the stator side at the same angle. It does to 0.06 % of the range of A; turned the other way it is off by 6.7 %, not turned
    by 37 %.
 6. Each slot layer's phase and sign. motorgeo.winding(p) gives, for every slot k, (phase, sign) of its inner layer (nearer the
    bore) and of its outer layer. build() gives each layer region circuit = phase and turns = sign x 10, and FEMM multiplies the
    circuit's current by the turns: the layer carries sign x 10 x i_phase amperes, and a positive current flows out of the screen
    (FEMM's convention: checked on a single wire, whose B circles it counter-clockwise). So the label "A+" says which phase and
    which way that layer is wound, and the tint says which way its current flows now: red where sign x i_phase is positive, blue where
    it is negative, deeper with more current (the Explorer scales the tint to its 20 A slider; here the full tint is this
    picture's 10 A peak, so the weakest phase still shows its sign). Two checks: the labels are read from the regions and compared
    with the winding table, and FEMM's own total current of every layer (block integral 7) is compared with sign x 10 x i_phase.
 7. The magnets. Each arrow is the model's magnetisation direction, region.magdir (degrees from +x, as drawn), turned with the rotor,
    so it points to the magnet's north; it starts a little before the middle of the slot and ends past its north face. The check:
    in the .ans (the rotor's own frame) the mean B inside each magnet points the same way as its arrow, to within 0.2 degrees.

Numbers to expect (FEMM 5.0, 2026-11-08; the solves repeat to the last digit)
    mesh            17,934 nodes, 34,758 triangles; A from -0.0092 to +0.0092 Wb/m, so 36 flux lines at 5e-4 Wb/m
    torque          T (air-gap integral) at theta_e = 0, 10, ... 50: 1.2244, 1.0815, 1.1452, 1.0367, 1.1737, 1.1033 N.m, mean
                    1.1275; learn/mtpa: 1.1275 (-0.004 %). Its rows in solves.csv agree position by position to 0.001 %. The figure's
                    own position sits 8.6 % above the mean: it is a crest of the torque ripple.
    |B|             up to 2.83 T, in the rotor's outer bridges (the 1.4 mm of steel between the end of a magnet slot and the rotor
                    surface); 2.5 % of the triangles are above the 2.2 T end of the scale, almost all within 2 mm of a magnet
                    slot, and are drawn in its deepest blue
    checks          |B| from the .ans against FEMM's: 1e-14 T; 18 slot layers agree with FEMM's total current; magnets within
                    0.2 degrees; the rotor turn: 0.06 % against 6.7 % and 37 %
    picture         1400 x 1518 px, about 0.90 MB (a 256-colour palette would make it a third of that; this keeps every colour)
"""
import argparse
import faulthandler
import json
import math
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))               # the code folder, wherever this is run from
sys.path.append(str(Path(__file__).resolve().parent / "dq_torque"))        # common.py, shared with the d-q and MTPA notes
import matplotlib
matplotlib.use("Agg")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
import shapely
from matplotlib.colors import LinearSegmentedColormap, to_rgb, to_rgba
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle, Wedge
from PIL import Image

import common as dq
import femfile
import femm_queue
import femm_session
import motorgeo as mg
import paths

# ----------------------------------------------------------------------------- 0. what is drawn, and where it goes
OUT = paths.path("out") / "learn" / "field_figure"
WORK = OUT / "field_figure_work.fem"       # FEMM's working file; its .ans (field_figure_work.ans) is overwritten by every solve
STATE = OUT / "field_figure.ans"           # the solution of the picture's own state, copied out under a name of its own
FACTS = OUT / "field_figure.json"          # what the solve found, so that --render-only can draw without FEMM
PNG = OUT / "field_figure.png"

I_PEAK = dq.I_RATED                        # 10 A peak
GAMMA_FALLBACK = 22.9                      # [deg] used only if learn/mtpa has not been run
B_MAX = 2.2                                # the colour scale of the Rotor Flux Explorer: 0 ... 2.2 T
ARROW = (-1.0, 3.3)                        # a magnet's arrow runs from here to there [mm], along its magnetisation, from its middle
TIMEOUT_S = 300                            # the solves take about 4 s; a FEMM that hangs cannot hold the slot for longer

# The Rotor Flux Explorer's colours (flux_explorer_template.html, light theme)
RAMP = ["#fbfbf9", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]      # |B| = 0 ... B_MAX, evenly spaced stops
CMAP = LinearSegmentedColormap.from_list("explorer", RAMP)
INK, INK2, MUTED, RULE = "#15181b", "#4b4f53", "#84898d", "#dcdfda"
FLUX, MAGNET = "#1e2328", "#d1495b"
POS, NEG, ZERO = "#e34948", "#2a78d6", "#ecebe6"                              # slot current: out of / into the screen
FONT = ["Segoe UI", "DejaVu Sans"]
FONT_DATA = ["Consolas", "DejaVu Sans Mono"]                                  # the phase labels, in the Explorer's monospace

# The picture's layout, in pixels (the figure is 1400 px wide at 100 dpi, so 1 pt = 100/72 px)
W, SIDE, TOP, BOTTOM, GAP, ROW = 1400, 22, 16, 18, 14, 38
PX = 72 / 100                              # points per pixel: the line widths and font sizes below are written in pixels


# ----------------------------------------------------------------------------- 1. the working point
def mtpa_point():
    """(best current angle [deg], torque [N.m]) at 10 A from learn/mtpa's report; (GAMMA_FALLBACK, None) if it was not run."""
    try:
        summary = json.loads((paths.path("out") / "learn" / "mtpa" / "summary.json").read_text(encoding="utf-8"))
        row = summary["rows"][f"{I_PEAK:g}"]
        return row["FEMM angle [deg]"], row["torque at FEMM's angle [N.m]"]
    except (OSError, KeyError, ValueError):
        return GAMMA_FALLBACK, None


def phase_currents(gamma, theta_e, I=I_PEAK):
    """(i_a, i_b, i_c) [A] for a current of peak I at angle gamma (from +q toward -d), rotor at theta_e [electrical degrees].

    Step 1, the d and q currents (common.currents): i_d = -I sin(gamma), i_q = I cos(gamma).
    Step 2, the inverse Park transform, the one FemmModel.solve applies: phase k sees the rotor angle shifted by 0, -120, +120
    degrees, and i_k = i_d cos(th_k) - i_q sin(th_k). Written out, i_k = -I sin(th_k + gamma)."""
    i_d, i_q = dq.currents(I, gamma)
    return [i_d * math.cos(th) - i_q * math.sin(th) for th in (math.radians(theta_e + s) for s in (0, -120, 120))]


# ----------------------------------------------------------------------------- 2. the solve: the only part that needs FEMM
def check_against_femm(regions, currents):
    """Ask FEMM itself two things about the solution that was just made (its .fem is still open in this FEMM).

    1. |B|: does the |B| this script computes from the saved .ans (curl A on each triangle) equal the one FEMM's post-processor
       gives at the middle of 40 triangles spread over the mesh? Returns the largest difference [T].
    2. The slot layers: FEMM's own "total current" of every coil layer (block integral 7, positive = out of the screen)
       against turns x phase current, which is what the picture colours by. Returns [(layer, phase, turns, FEMM's, mine)]."""
    import femm
    sol = femfile.read_ans(str(STATE))
    bx, by = femfile.element_b(sol["xy"], sol["A"], sol["tri"])
    mine = np.hypot(bx, by)
    middle = sol["xy"][sol["tri"]].mean(axis=1)
    pick = np.random.default_rng(1).choice(len(mine), 40, replace=False)
    amps = dict(zip("ABC", currents))
    femm.mi_loadsolution()
    femm.mo_smooth("off")                                           # FEMM's default smooths B between triangles; off, it gives the triangle's own
    try:
        theirs = [math.hypot(*femm.mo_getb(*middle[i])[:2]) for i in pick]
        layers = []
        for r in regions:
            if r.name.startswith("coil_"):
                femm.mo_selectblock(*r.label)
                total = complex(femm.mo_blockintegral(7)).real
                femm.mo_clearblock()
                layers.append((r.name, r.circuit, r.turns, total, r.turns * amps[r.circuit]))
    finally:
        femm.mo_close()
    return float(max(abs(a - mine[i]) for a, i in zip(theirs, pick))), layers


def magnet_check(regions):
    """For every magnet: the angle between the arrow (its magnetisation direction) and the mean B inside it [deg].

    Both are read in the rotor's own frame, the one the .ans keeps the rotor in. B inside a magnet points the way the magnet
    does, give or take the stray field of the rest of the motor, so the two angles agree to a fraction of a degree."""
    sol = femfile.read_ans(str(STATE))
    bx, by = femfile.element_b(sol["xy"], sol["A"], sol["tri"])
    corners = sol["xy"][sol["tri"]]
    mid = corners.mean(axis=1)
    area = 0.5 * np.abs((corners[:, 1, 0] - corners[:, 0, 0]) * (corners[:, 2, 1] - corners[:, 0, 1])
                        - (corners[:, 2, 0] - corners[:, 0, 0]) * (corners[:, 1, 1] - corners[:, 0, 1]))
    rows = []
    for r in regions:
        if r.material == "N45":
            inside = shapely.contains_xy(r.polygon(), mid[:, 0], mid[:, 1])
            direction = math.degrees(math.atan2(np.average(by[inside], weights=area[inside]),
                                                np.average(bx[inside], weights=area[inside])))
            rows.append((r.name, r.magdir, direction, (direction - r.magdir + 180) % 360 - 180))
    return rows


def solve(gamma, theta_e, torque_mtpa):
    """Solve the working point at theta_e and, for the torque check, the other five rotor positions; write the .ans and the facts."""
    import femm
    p = dq.params()
    regions = mg.build(p)
    OUT.mkdir(parents=True, exist_ok=True)
    i_d, i_q = dq.currents(I_PEAK, gamma)
    with femm_queue.solver_slot("learn/field_figure: the field of one working point"):
        faulthandler.dump_traceback_later(TIMEOUT_S, exit=True)     # the net under this run: see the module docstring, method 1
        femm_session.open_femm()                                    # this process's own FEMM: own port, below-normal priority
        try:
            t0 = time.perf_counter()
            mg.to_femm(p, regions, str(WORK), fine=True)
            femm.mi_close()
            model = mg.FemmModel(p, str(WORK))
            rows = {theta_e: model.solve(theta_e, i_d, i_q, keep_ans=str(STATE))}
            b_err, layers = check_against_femm(regions, phase_currents(gamma, theta_e))
            for th in dq.MEAN_POS:
                if th not in rows:
                    rows[th] = model.solve(float(th), i_d, i_q)
            seconds = time.perf_counter() - t0
        finally:
            femm_session.close_femm()
            faulthandler.cancel_dump_traceback_later()
    first = rows[theta_e]
    mine = phase_currents(gamma, theta_e)
    given = [first["ia"], first["ib"], first["ic"]]
    assert np.allclose(mine, given, atol=1e-9), f"the currents FEMM was given {given} are not {mine}"
    assert b_err < 1e-6, f"|B| from the .ans differs from FEMM's by {b_err} T"
    for name, phase, turns, total, expected in layers:
        assert abs(total - expected) < 1e-6 * max(1.0, abs(expected)), f"{name}: FEMM says {total} A, the picture colours {expected} A"
    magnets = magnet_check(regions)
    facts = {"I": I_PEAK, "gamma": gamma, "theta_e": theta_e, "rotor_deg": first["rotor_deg"], "i_abc": mine,
             "T_age": first["T_age"], "T_wst": first["T_wst"], "T_dq": first["T_dq"],
             "T_by_position": {f"{th:g}": rows[th]["T_age"] for th in sorted(rows)},
             "T_mean_6": float(np.mean([rows[th]["T_age"] for th in dq.MEAN_POS])), "T_mtpa": torque_mtpa,
             "b_check_max_diff_T": b_err, "slot_layers_checked": len(layers), "layers": layers,
             "magnet_B_vs_arrow_deg": {name: off for name, _, _, off in magnets}, "solve_seconds": seconds}
    FACTS.write_text(json.dumps(facts, indent=1), encoding="utf-8")
    return facts


# ----------------------------------------------------------------------------- 3. the field: from the .ans to triangles
def field_triangles(p, rotor_deg):
    """The solution cut into triangles that cover the 0 to 90 degree window, with the rotor turned to where it really is.

    What the .ans holds: nodes (x, y) in mm, the vector potential A at every node, and the triangles. The stator side sits
    where it was drawn. The rotor side also sits where it was drawn, whatever the rotor angle: FEMM turns the rotor only inside
    the air-gap element, never in the mesh. So the rotor's nodes are turned here, by rotor_deg, counter-clockwise (turn_check
    tests that this is the right turn). Turned by rotor_deg, the rotor's 90 degree wedge no longer lines up with the window, so
    three copies of everything, turned by -90, 0 and +90 degrees more, are laid side by side (the sector has periodic boundaries:
    the field repeats every 90 degrees, same A) and the triangles with a corner in the window are kept. The window cuts the rest off.
    The ring in the middle of the air gap is not meshed, and would show as a white seam: the edge nodes of the two sides, at a
    third and two thirds of the gap, are moved to the middle of it.
    Returns (matplotlib Triangulation, A at every node of it, |B| of every triangle in T)."""
    sol = femfile.read_ans(str(STATE))
    xy, A, tri = sol["xy"], sol["A"], sol["tri"]
    bx, by = femfile.element_b(xy, A, tri)                          # B = curl A on each triangle (method 3 of the module docstring)
    radius = np.hypot(xy[:, 0], xy[:, 1])
    is_rotor = radius < p.r_ro + p.airgap / 2                       # the unmeshed ring of the gap lies between the two sides
    r_mid = p.r_ro + p.airgap / 2
    for r_edge in (p.r_ro + p.airgap / 3, p.r_ro + 2 * p.airgap / 3):   # close the ring: each side's edge nodes move to mid-gap
        edge = np.abs(radius - r_edge) < 1e-3
        xy[edge] *= (r_mid / radius[edge])[:, None]
    n = len(xy)
    X, Y, T = [], [], []
    for k, copy in enumerate((-1, 0, 1)):
        turn = np.radians(np.where(is_rotor, rotor_deg, 0.0) + 90.0 * copy)
        c, s = np.cos(turn), np.sin(turn)
        X.append(c * xy[:, 0] - s * xy[:, 1])
        Y.append(s * xy[:, 0] + c * xy[:, 1])
        T.append(tri + k * n)
    X, Y, T = np.concatenate(X), np.concatenate(Y), np.vstack(T)
    corner_deg = np.degrees(np.arctan2(Y[T], X[T]))
    keep = ((corner_deg > -0.01) & (corner_deg < 90.01)).any(axis=1)
    return mtri.Triangulation(X, Y, T[keep]), np.tile(A, 3), np.tile(np.hypot(bx, by), 3)[keep]


def turn_check(p, rotor_deg):
    """How well A matches across the unmeshed ring of the gap for three ways of placing the rotor [% of the range of A].

    The test that the rotor is turned the right way (counter-clockwise by rotor_deg) and not the other way, or not at all.
    Across the ring, 0.27 mm thick, A hardly changes. So A at a rotor-side node, turned to where the rotor really is, must
    equal A on the stator side at that angle. Only the right turn makes it so, and the picture relies on it being under 1 %.
    The two wrong ones are for scale: they give tens of per cent, except that the other turn lands on the same place, modulo
    the 90 degree period, when rotor_deg is 45 degrees (theta_e = 10): there the two cannot be told apart, and need not be."""
    sol = femfile.read_ans(str(STATE))
    xy, A = sol["xy"], sol["A"]
    radius, angle = np.hypot(xy[:, 0], xy[:, 1]), np.degrees(np.arctan2(xy[:, 1], xy[:, 0]))
    rotor_side = np.abs(radius - (p.r_ro + p.airgap / 3)) < 1e-3
    stator_side = np.abs(radius - (p.r_ro + 2 * p.airgap / 3)) < 1e-3
    order = np.argsort(angle[stator_side])
    a_stator, A_stator = angle[stator_side][order], A[stator_side][order]
    out = {}
    for name, turn in (("turned by rotor_deg", rotor_deg), ("turned the other way", -rotor_deg), ("not turned", 0.0)):
        stator_A = np.interp(angle[rotor_side] + turn, a_stator, A_stator, period=p.sector_deg)   # the field repeats every 90 deg
        out[name] = 100 * float(np.sqrt(np.mean((stator_A - A[rotor_side]) ** 2)) / (A.max() - A.min()))
    return out


def flux_levels(A, lines):
    """Levels of A for the flux lines: equal steps s of a round number, about `lines` of them over the range of A.

    The levels are A = (k + 1/2) s, half a step off zero. The rim and the shaft sit at A = 0 (the boundary condition), and a line
    there would only trace their outlines; and with no level at 0 the gap across A = 0 is no wider than any other."""
    step = (A.max() - A.min()) / lines
    exp = 10 ** math.floor(math.log10(step))
    step = min(m * exp for m in (1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) if m * exp >= step * 0.9)
    ks = np.arange(math.ceil(A.min() / step - 0.5), math.floor(A.max() / step - 0.5) + 1)
    return step, (ks + 0.5) * step


# ----------------------------------------------------------------------------- 4. the picture
def rotated(pts, deg):
    t = math.radians(deg)
    c, s = math.cos(t), math.sin(t)
    pts = np.asarray(pts, dtype=float)
    return np.column_stack([c * pts[:, 0] - s * pts[:, 1], s * pts[:, 0] + c * pts[:, 1]])


def arc(r):
    t = np.radians(np.linspace(0, 90, 361))
    return r * np.cos(t), r * np.sin(t)


def slot_colour(amps):
    """The Explorer's slot tint: grey for no current, red for current out of the screen, blue for into it, deeper with more."""
    t = min(1.0, abs(amps) / I_PEAK)
    end = to_rgb(POS if amps >= 0 else NEG)
    return tuple(z + (e - z) * t for z, e in zip(to_rgb(ZERO), end))


def draw_stator(ax, p, regions, amps):
    """The slots (two layers each: tint and phase label), the openings, and the stator's grey outlines."""
    table = mg.winding(p)
    for r in regions:
        if r.name.startswith("coil_"):
            k, j = (int(v) for v in r.name.split("_")[1:])
            sign = 1 if r.turns > 0 else -1
            assert table[k][j] == (r.circuit, sign), f"{r.name}: the region and the winding table disagree"
            tint = to_rgba(slot_colour(sign * amps[r.circuit]), 0.72)          # sign x phase current: red out of the screen, blue into it
            ax.add_patch(Polygon(np.array(r.polygon().exterior.coords), closed=True, facecolor=tint, edgecolor=MUTED,
                                 linewidth=1.0 * PX, zorder=3.5))              # over the flux lines, as in the Explorer: they show through
        elif r.name.startswith("opening_"):
            ax.add_patch(Polygon(np.array(r.polygon().exterior.coords), closed=True, fill=False, edgecolor=MUTED,
                                 linewidth=1.2 * PX, zorder=4))
    for r_, lw in ((p.r_so, 1.4), (p.r_si, 1.2)):                              # the rim and the bore
        ax.plot(*arc(r_), color=MUTED, linewidth=lw * PX, zorder=4, solid_capstyle="butt")
    halo = [pe.withStroke(linewidth=3.2 * PX, foreground="white", alpha=0.75)]
    for r in regions:                                                          # the labels go over everything
        if r.name.startswith("coil_"):
            x, y = r.polygon().centroid.coords[0]
            ax.text(x, y, f"{r.circuit}{'+' if r.turns > 0 else '−'}", ha="center", va="center", fontsize=26 * PX,
                    fontfamily=FONT_DATA, fontweight="bold", color=INK, zorder=7, path_effects=halo)


def draw_rotor(ax, window, p, regions, rotor_deg):
    """The rotor's grey outlines, and every magnet: red outline and an arrow to its north, turned with the rotor."""
    for r_ in (p.r_ro, p.r_shaft):
        ax.plot(*arc(r_), color=MUTED, linewidth=1.2 * PX, zorder=4, solid_capstyle="butt")
    for r in regions:
        if r.material != "N45":
            continue
        outline = np.array(r.polygon().exterior.coords)
        for copy in (-1, 0, 1):                                                # the same three copies as the field's
            turn = rotor_deg + 90 * copy
            ax.add_patch(Polygon(rotated(outline, turn), closed=True, fill=False, edgecolor=MAGNET, linewidth=2.2 * PX,
                                 zorder=5)).set_clip_path(window)
            centre = rotated([r.label], turn)[0]
            d = math.radians(r.magdir + turn)                                  # the magnetisation direction, turned with the rotor
            u = np.array([math.cos(d), math.sin(d)])
            arrow = FancyArrowPatch(centre + ARROW[0] * u, centre + ARROW[1] * u, arrowstyle="-|>,head_length=0.55,head_width=0.22",
                                    mutation_scale=34 * PX, color=MAGNET, linewidth=2.6 * PX, shrinkA=0, shrinkB=0, zorder=6,
                                    path_effects=[pe.withStroke(linewidth=5.5 * PX, foreground="white", alpha=0.8)])
            ax.add_patch(arrow).set_clip_path(window)


def draw_legend(fig, ov, facts):
    """The strip under the sector, three rows. `ov` spans the figure in pixel units (origin bottom left), so every item sits where it is put."""
    renderer = fig.canvas.get_renderer()
    size = 21 * PX

    def put(x, y, s, color=INK2, weight="normal"):                             # write text at (x, y), return where it ends
        t = ov.text(x, y, s, ha="left", va="center", fontsize=size, fontfamily=FONT, fontweight=weight, color=color)
        return x + t.get_window_extent(renderer).width

    y1, y2, y3 = BOTTOM + 2.5 * ROW, BOTTOM + 1.5 * ROW, BOTTOM + 0.5 * ROW
    x = SIDE + 4
    x = put(x, y1, "|B|", INK2, "bold") + 12                                   # row 1: the colour bar, then the magnet
    ov.imshow(np.linspace(0, 1, 256)[None, :], extent=(x, x + 240, y1 - 9, y1 + 9), cmap=CMAP, aspect="auto", zorder=2)
    ov.add_patch(Rectangle((x, y1 - 9), 240, 18, fill=False, edgecolor=RULE, linewidth=1.2 * PX, zorder=3))
    x = put(x + 240 + 12, y1, f"0 … {B_MAX:g} T") + 44
    ov.add_patch(Rectangle((x, y1 - 7), 40, 14, fill=False, edgecolor=MAGNET, linewidth=2.2 * PX))
    ov.add_patch(FancyArrowPatch((x + 20, y1 - 4), (x + 20, y1 + 24), arrowstyle="-|>,head_length=0.55,head_width=0.22",
                                 mutation_scale=34 * PX, color=MAGNET, linewidth=2.6 * PX, shrinkA=0, shrinkB=0))
    put(x + 52, y1, "red outline: magnet, arrow = its north")

    x = put(SIDE + 4, y2, "slots:", INK2, "bold") + 12                         # row 2: the two tints
    ov.add_patch(Rectangle((x, y2 - 9), 30, 18, facecolor=to_rgba(POS, 0.72), edgecolor=MUTED, linewidth=1.0 * PX))
    x = put(x + 40, y2, "current out of screen") + 14
    x = put(x, y2, "/") + 14
    ov.add_patch(Rectangle((x, y2 - 9), 30, 18, facecolor=to_rgba(NEG, 0.72), edgecolor=MUTED, linewidth=1.0 * PX))
    x = put(x + 40, y2, "into screen") + 14
    put(x, y2, "(deeper tint = larger current)", MUTED)

    put(SIDE + 4, y3, f"band-free rotor · {facts['I']:g} A peak at {facts['gamma']:.1f}° (MTPA) · θe {facts['theta_e']:g}°", INK, "bold")   # row 3


def draw(p, regions, facts, lines):
    """Draw the figure from the saved solution and the facts of the run; returns (flux-line step, its levels, the picture's height in px)."""
    rotor_deg = facts["rotor_deg"]
    amps = dict(zip("ABC", facts["i_abc"]))
    R = p.r_so
    tri, A, B = field_triangles(p, rotor_deg)
    step, levels = flux_levels(A, lines)
    print(f"field: A from {A.min():+.4f} to {A.max():+.4f} Wb/m; |B| up to {B.max():.2f} T, "
          f"{100 * np.mean(B > B_MAX):.1f} % of the triangles above the {B_MAX:g} T end of the scale")

    pad = 0.5                                                                  # mm of white around the sector
    scale = (W - 2 * SIDE) / (R + 2 * pad)                                     # pixels per mm
    side = round((R + 2 * pad) * scale)
    H = TOP + side + GAP + 3 * ROW + BOTTOM
    fig = plt.figure(figsize=(W / 100, H / 100), dpi=100, facecolor="white")
    ax = fig.add_axes([SIDE / W, (BOTTOM + 3 * ROW + GAP) / H, side / W, side / H])
    ax.set_xlim(-pad, R + pad)
    ax.set_ylim(-pad, R + pad)
    ax.set_aspect("equal")                                                     # 1 mm is the same length along x and y
    ax.axis("off")
    window = Wedge((0, 0), R, 0, 90, transform=ax.transData)                   # the field is cut to the 90 degree sector

    # |B|: one flat colour per triangle (B is constant on a triangle), the Explorer's ramp, 0 ... 2.2 T. No antialiasing here:
    # the triangles share their edges, and antialiased they would leave hairline seams between them.
    pc = ax.tripcolor(tri, facecolors=B, cmap=CMAP, vmin=0, vmax=B_MAX, edgecolors="face", linewidth=0.5, antialiased=False, zorder=1)
    pc.set_clip_path(window)
    # the flux lines: contours of A at equal steps (negative levels solid, not dashed)
    cs = ax.tricontour(tri, A, levels=levels, colors=FLUX, linewidths=1.25 * PX, alpha=0.75, zorder=3, negative_linestyles="solid")
    cs.set_clip_path(window)
    for x0, x1, y0, y1 in ((p.r_shaft, R, 0, 0), (0, 0, p.r_shaft, R)):         # the two cut edges of the window
        ax.plot([x0, x1], [y0, y1], color=RULE, linewidth=1.2 * PX, zorder=4, solid_capstyle="butt")

    draw_stator(ax, p, regions, amps)
    draw_rotor(ax, window, p, regions, rotor_deg)

    ov = fig.add_axes([0, 0, 1, 1], facecolor="none")                          # the legend's own axes, in pixels
    ov.set_xlim(0, W)
    ov.set_ylim(0, H)
    ov.set_autoscale_on(False)                                                 # imshow would otherwise reset the limits to the bar's extent
    ov.axis("off")
    draw_legend(fig, ov, facts)
    ov.set_xlim(0, W)
    ov.set_ylim(0, H)
    fig.canvas.draw()                                                          # saved without matplotlib's alpha channel: the same pixels, 8 % fewer bytes
    Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3]).save(PNG, optimize=True, compress_level=9)
    plt.close(fig)
    return step, levels, H


# ----------------------------------------------------------------------------- 5. command line
def report(facts):
    """Print what the solve found: the torque against learn/mtpa's, and the checks made on FEMM's own answers."""
    th = facts["theta_e"]
    print(f"torque (air-gap integral, fine mesh) at theta_e {th:g}: {facts['T_age']:.4f} N.m  "
          f"(weighted stress tensor {facts['T_wst']:.4f}, d-q formula {facts['T_dq']:.4f})")
    print("torque by rotor position [N.m]  " + "  ".join(f"{k}:{v:.4f}" for k, v in facts["T_by_position"].items()))
    line = f"mean over the six positions {facts['T_mean_6']:.4f} N.m"
    if facts["T_mtpa"] is not None:
        line += f"; learn/mtpa's value at 10 A, best angle: {facts['T_mtpa']:.4f} N.m ({100 * (facts['T_mean_6'] / facts['T_mtpa'] - 1):+.3f} %)"
    print(line)
    print(f"|B| from the .ans against FEMM's own mo_getb at 40 triangles: largest difference {facts['b_check_max_diff_T']:.1e} T")
    print(f"slot layers: FEMM's total current against turns x phase current, {facts['slot_layers_checked']} layers, all agree:")
    for name, phase, turns, total, expected in facts["layers"]:
        print(f"    {name:9s} phase {phase}  {turns:+d} turns   FEMM {total:+8.3f} A   turns x i {expected:+8.3f} A   "
              f"{'out of the screen (red)' if total >= 0 else 'into the screen (blue)'}")
    print("magnets: the mean B inside each, against its arrow [deg]: " + ", ".join(f"{n} {d:+.1f}" for n, d in facts["magnet_B_vs_arrow_deg"].items()))
    print(f"the six solves, with the fine mesh built: {facts['solve_seconds']:.1f} s")


def main():
    ap = argparse.ArgumentParser(description="The field figure: the band-free motor at 10 A at its MTPA angle.")
    ap.add_argument("--render-only", action="store_true", help="draw again from the saved solution; no FEMM")
    ap.add_argument("--theta-e", type=float, default=0.0, help="rotor position [electrical degrees] (default 0)")
    ap.add_argument("--lines", type=int, default=36, help="about this many flux lines (default 36)")
    args = ap.parse_args()
    p = dq.params()
    regions = mg.build(p)
    gamma, torque_mtpa = mtpa_point()
    if args.render_only:
        facts = json.loads(FACTS.read_text(encoding="utf-8"))
    else:
        facts = solve(gamma, args.theta_e, torque_mtpa)
        report(facts)
    print(f"working point: {facts['I']:g} A peak at {facts['gamma']:.3f} deg, theta_e {facts['theta_e']:g} deg, rotor turned {facts['rotor_deg']:.2f} deg")
    print("phase currents [A]: A %+.3f  B %+.3f  C %+.3f" % tuple(facts["i_abc"]))
    mismatch = turn_check(p, facts["rotor_deg"])
    print("A across the gap, rms mismatch as % of the range of A: " + ", ".join(f"{k} {v:.2f} %" for k, v in mismatch.items()))
    assert mismatch["turned by rotor_deg"] < 1.0, "the rotor is not where the picture puts it"
    step, levels, H = draw(p, regions, facts, args.lines)
    print(f"flux lines: {len(levels)} at steps of {step:.2e} Wb/m = {step * p.stack * 1e-3 * 1e6:.2f} uWb per gap over the {p.stack:g} mm stack")
    print(f"wrote {PNG}: {W} x {H} px, {PNG.stat().st_size / 1e3:.0f} kB")


if __name__ == "__main__":
    main()
