"""stress_ccx: the optimized rotor under centrifugal load in CalculiX (bench case B4, run on the optimized geometry), its safe speeds, and what its two thin strips cost.

What it shows, and why
    The optimizer thinned the two strips that hold each pole's outer steel (its "cap") and its magnets against the centrifugal pull: the centre rib is 1.0 mm (the reference's is 3.3 mm)
    and the outer bridges are 1.4 mm (the reference's 1.56 mm). Nothing in the study looked at the stress in them: the bridge guard of 1.4 mm and the rib floor of 1.0 mm are
    manufacturing limits, not strength limits. This script puts the optimized rotor through the bench's validated structural case, B4 (bench/b4_ccx.py), and answers four questions:
      1. Does the case still give B4's number on the reference rotor? (the regression gate: 33.45 MPa within 0.5 %, or nothing else is run)
      2. What is the stress in the optimized rotor, where does it peak, and are the checks that made B4 trustworthy still met (equilibrium, the (12/10)^2 scaling, the four slots'
         symmetry, the mesh, the strips' force against the hand estimate, the bore against Lame, the magnets' own stress)?
      3. At what speeds does the peak reach the placeholder strength limits (04 Structural 6.4)? There is no ruled top speed (the supply is 48-60 V), so the speeds are the answer.
      4. How much does the peak move when the rib, the bridge or the slot's corner radius is changed? (these are the limits the optimizer's design sits on, which are still open decisions)
    It reports stress only. The electromagnetic cost of a wider rib or bridge is not in it.

Method (the maths, then the files)
    1. The load. A bit of steel of density rho at distance r from the axis, in a rotor turning at w rad/s (w = 2 pi n / 60), is held on its circle by a force rho w^2 r per unit
       volume: the centrifugal body force. CalculiX's `*DLOAD ... CENTRIF` applies exactly that, about the z axis, with w^2 as its magnitude. A magnet is a mass m on a wall: its
       whole force is F = m w^2 r_c along the radius through its centre (r_c), m the real magnet's mass (9.08 x 1.57 mm x the stack x 7,600 kg/m3). The pressure model (B4's verdict)
       lets the magnet rest loose in its slot: F is split along the slot's own axes, u along it and n across it. F.n is a uniform pressure on the long wall the magnet leans on,
       F.u a uniform pressure on its outer end wall (F.n / (wall length x stack), in MPa): what a frictionless magnet does. The pressures reach the mesh as the consistent nodal forces
       of a quadratic edge (1/6, 4/6, 1/6). The slot is a hole and carries no hoop stress, so the model is the conservative one: a glued magnet (the "bonded" model, a body tied to
       the walls) would take part of the load and cut the bridges' peak several times over (B4: 7.9 times); a glue far softer than steel and a magnet that cracks at 78 MPa in tension
       cannot be counted on for that.
    2. The model. Plane stress (CPS6, six-node triangles), the 90-degree sector with y fixed on the x axis and x fixed on the y axis (both are q axes of the rotor), a free bore, the
       stack's 20 mm. Gmsh meshes the steel finely at the slot walls (h_min) and coarsely elsewhere; units mm, N, tonne, s, so stress is in MPa. All of it is b4_ccx.py's code, imported:
       build_model, loads_at, write_deck, solve (the solver slot, below normal, 4 threads, a timeout), analyse, hand_estimate, strip_geometry. Nothing in bench/, design.py or motorgeo is
       changed. B4 is a case of the reference rotor only; every function of it takes the motorgeo.Params of the rotor, and design.params(overrides) gives the optimized one.
    3. Why one solve gives every speed. The deck is K u = f, with the stiffness K of the steel and a load f that is w^2 times a load f1 that does not depend on the speed (rho r for
       the steel, m r_c for each magnet, both fixed by the geometry). Linear elasticity then gives u = w^2 K^-1 f1 and every stress s = w^2 (C B K^-1 f1): the whole field is
       proportional to w^2. This holds while nothing else loads the rotor (no press fit, no heat, no magnet that lifts off its wall): the pressure model has none of those.
       So the peak at n is the peak at the solved speed n_ref times (n / n_ref)^2, and the speed at which the peak reaches a limit s_lim is
           n = n_ref x sqrt(s_lim / s_peak(n_ref)).
       The case checks it: the 12,000 rpm solve gives the 10,000 rpm peak times (12/10)^2 = 1.44.
    4. The limits, from 04 Structural 6.4 (placeholders until M-19's own values are found: the maker's sheet skips M-19's mechanics, so they come from its neighbour M-15 and from
       similar grades). Yield 300 MPa (the low end of similar grades, 300-450; M-15 has 358). Fatigue about 200 MPa: 0.4 x M-15's 490 MPa tensile strength = 196 MPa, polished and fully
       reversed, which a punched edge lowers (by a factor 04 Structural leaves to S5). A bridge's cycle is zero to peak, so its mean stress counts: Goodman's line gives the peak
       2 / (1/196 + 1/490) = 280 MPa. The speeds asked for, for each rotor: (a) the speed at which the peak reaches the yield; (b) yield / 1.5 and yield / 2, each taken at a 1.2 x
       overspeed, so the maximum operating speed is that speed divided by 1.2; (c) the speeds at which the peak reaches 196 MPa and 280 MPa (one zero-to-peak cycle per start and stop).
       The peak is a linear-elastic value at a 0.3 mm fillet: yielding or a notch's size effect would make the real local stress lower, so the speeds are cautious on that count; the
       punched edge (a weaker edge) and the mesh (the main mesh reads a little low, see 6) go the other way.
    5. The numbers of the case, for each rotor (the reference first: the optimized rotor is judged against what B4 trusted). Solves of the main mesh (h_min 0.1 mm): pressure and bonded
       models at 10,000 and 12,000 rpm. Equilibrium: the reactions of the two symmetry edges against the load worked out from the mesh. Where the peak sits (the strip's zone within 2 mm
       of its cut, and the corner of the slot). The three strips' force from the stress integrated over the cut across each, F = thickness x integral of (sigma n), their nominal stress
       F / (width x stack), against the hand estimate (04 Structural 6.1: the cap's and the magnets' pull shared evenly by the three widths); the peak/nominal factor is the stress
       concentration the hand estimate lacks. The bore's mean hoop stress against Lame's plain disc. The bonded model's largest stress in a magnet body against N45's tensile strength.
    6. The mesh. The verdict peak at h_min 0.2, 0.1 and 0.05 mm (the element size at the slot walls, halved each time). If the error of the peak falls as h^p, the extrapolated value is
       P* = P(0.05) + (P(0.05) - P(0.1)) / (2^p - 1), with the order p = log2((P(0.2) - P(0.1)) / (P(0.1) - P(0.05))) read from the three values (B4 found p = 2.0, as quadratic
       triangles should give). The main mesh reads (P(0.1) - P*) / P* against that.
    7. The sensitivities, all on the optimized design with one thing changed and the verdict model at 10,000 rpm (design.screen checks each geometry first and the strips it really has
       are printed). The rib: moving mag_y moves the rib twice as far (the rib is twice the height of the slot's lowest point above the d-axis, and a slot that is moved along x
       keeps that height), so mag_y is moved until the rib has the asked width (the design's 1.0, then 1.25, 1.4, 1.5, 2.0 and 2.5 mm: the brief's 1.5, 2.0, 2.5 and two more where the
       bridges take over from the rib), and mag_x is re-solved (design.mag_x_for_bridge) so that the bridge stays at 1.4 mm.
       The bridge: mag_x is re-solved for 1.4, 1.6, 1.8 and 2.0 mm; mag_y does not move, and the rib stays what it is (a slot moved along x keeps its rib). The corner radius:
       motorgeo.Params.slot_fillet (0.3 mm in the design; a design note advises 0.5 mm or more) at 0.5 mm, (a) with the slot where it is, and (b) with the rib and the bridge put back to
       the design's (mag_y and mag_x re-solved), so that the radius is the only difference; (b) is also solved at h_min 0.05 mm, because a larger radius is resolved better at the same
       mesh and the comparison at one h_min would flatter it. A 0.5 mm radius leaves the 9.08 x 1.57 mm magnet's sharp corner outside the rounded slot corner (the printed corner
       clearance): the magnet's corners would need to be chamfered or rounded by about 0.07 mm, a fact for the magnet's maker, not for the stress.
    8. The pictures (PNG, in the output folder): the optimized rotor's von Mises stress with its peak marked, the reference and the optimized rotor on one colour scale, and the peak
       against the rib and against the bridge.

Run (from the Python to FEMM folder; every ccx run takes the machine's one solver slot, below normal priority, at most 4 threads, and waits if another run holds it)
    uv run python -m optimization.stress_ccx                 the whole case: the gate, both rotors, the sensitivities, the pictures
    uv run python -m optimization.stress_ccx --gate-only     only the regression gate (one mesh, one solve, the reference rotor), written to out\\optimization\\stress\\gate\\
    uv run python -m optimization.stress_ccx --redo          no solve: the printed report and the pictures again from state.pkl, which the last full run saved (a few seconds); its own
                                                             log and results.json go to redo\\, and the full run's results.json is left as it was
Exit status 0 if the gate and every gating check pass (the solves end, equilibrium, the two speeds scale), 1 if the gate or a gating check fails, 130 after Ctrl-C.
Writes out\\optimization\\stress\\ (git ignores out\\): results.json (every number), run.log and runs.json (what the run said and every program it started, with the priority class it
was seen at), state.pkl (what --redo reads), the pictures, and solves\\<rotor or sweep>\\<run>\\ with job.inp, job.frd, job.dat and ccx's printout for every solve. The folder solves\\
is emptied at the start of a run; the other files are overwritten.

Numbers to expect (this PC, 2026-11-09, CalculiX 2.23, Gmsh 4.15.2; the whole case took 331 s here: 23 ccx solves of 2 to 31 s (268 s of it, every one seen at below-normal priority on 4
threads) and the rest Gmsh's meshes, the readers and the pictures; the same numbers came out of a run on a busier PC in 401 s, so allow 5 to 7 minutes; the gate alone is about 10 s)
    gate        the reference rotor's peak is 33.454 MPa at (42.81, 7.36) mm, +0.012 % from the brief's 33.45: B4's number to the digit (48.174 MPa at 12,000 rpm; bonded bridges 4.25 MPa;
                bore 14.908 MPa; equilibrium 0.00056 %; the three meshes 31.82, 33.45, 33.85 MPa, order 2.04, extrapolated 33.98 MPa)
    optimized   the peak is 39.061 MPa at 10,000 rpm and 56.248 MPa at 12,000 (ratio 1.4400): 16.8 % above the reference's, and in the CENTRE RIB, not in a bridge: the corner fillet at
                the inner end of magnet_1_up on its loaded side, at (13.95, 34.99) mm, 0.016 mm from the rib's cut, and on all three meshes. The bridges' zones peak at 34.70, 34.38, 34.68,
                34.52 MPa (spread 0.92 %) and the two ribs' at 38.76 and 39.06 MPa. design.screen finds the rib 1.0004 mm and the bridge 1.4000 mm (the reference's 3.2979 and 1.5585).
                The rib carries 465.1 N of the 636 N that a pole's strips carry (73.1 %; 23.25 MPa nominal over its 1.0 mm), the bridges 85.6 and 85.4 N (13.5 and 13.4 %; 3.06 and 3.05 MPa):
                the reference's rib carried 615.9 N (78.3 %; 9.34 MPa). The hand estimate gives 706.9 N and the strips carry 10.1 % and 9.9 % less (B4's limit is 10 %: pole 0 misses it by 0.12 points, because the magnets
                are 26.8 % of this pole's pull against 22.3 % and the hand method puts all of it on the cap; with the magnets' end-wall load left off the cap the free body is 636.2 N and the
                strips agree to -0.12 % and +0.08 %). Peak / nominal: 4.67 on the three strips together (hand 4.20), 1.68 on the rib's own nominal stress; the reference's 5.46 (5.06) and 12.19 on
                its bridge's. The bore's mean hoop stress is 14.941 MPa (+5.99 % on Lame's disc, 14.098) in the pressure model and 14.093 MPa (-0.03 %) in the bonded one, whose bridges peak at
                4.25 MPa (8.2 times less), its rib at 4.90 MPa and its magnet bodies at 3.77 MPa von Mises and 4.57 MPa principal (tension), against N45's 78 MPa; equilibrium 0.00052 % at worst.
    mesh        38.96, 39.06, 39.64 MPa at h_min 0.2, 0.1, 0.05 mm: the changes grow (+0.10, then +0.58 MPa), so no order can be read and no extrapolation is made; the main mesh reads 1.46 %
                under the finest (an assumed order 2 would give 39.83 MPa and 1.94 %), and every speed below is 0.7 % lower with the finest mesh's peak. The bridges' zones peak at 32.97, 34.70
                and 35.36 MPa on the same meshes (changes +1.73 and +0.66: order about 1.4), the ribs' at 38.96, 39.06 and 39.64; the two ribs, mirror images meshed differently, differ by
                about 1 % on one mesh (38.76 and 39.06 MPa on the main one).
    speeds      n = 10,000 x sqrt(limit / peak) [rpm]: the optimized rotor reaches the yield (300 MPa) at 27,713, yield / 1.5 (200) at 22,628 (maximum operating speed 18,857 with the 1.2 x
                overspeed), yield / 2 (150) at 19,596 (16,330), the fatigue placeholder (196) at 22,400 and Goodman's 280 MPa at 26,774; the reference 29,946, 24,451 (20,376), 21,175 (17,646),
                24,205 and 28,930. A tip speed of 80 m/s (16,977 rpm) gives the optimized peak 112.6 MPa; ROSS's first critical speed (22,838 rpm) 203.7 MPa.
    rib         the peak [MPa] at a rib of 1.000, 1.25, 1.4, 1.5, 2.0, 2.5 mm (bridge held at 1.40 mm): 39.06, 35.48, 34.32, 34.22, 33.77, 33.57; the first two in the rib, the others in a bridge;
                the rib's own zone: 39.06, 35.48, 33.90, 32.70, 29.64, 28.86 and the bridges' 34.70, 34.38, 34.32, 34.22, 33.77, 33.57, so the bridges take over at about 1.36 mm.
    bridge      the peak at a bridge of 1.4, 1.6, 1.8, 2.0 mm (rib held at 1.0004 mm): 39.06, 39.03, 39.50, 39.63, always in the rib; the bridges' zones fall only from 34.70 to 33.86.
    fillet      0.5 mm instead of 0.3, at the design's rib and bridge: 32.03 MPa (-18.0 %; -18.6 % on the finest mesh, 32.27 against 39.64), the bridges' zones -17.2 %; with the slot where it
                was, 29.89 MPa (-23.5 %: the rib is then 1.157 mm and the bridge 1.483 mm). The magnet's sharp corner then stands 0.066 mm inside the steel (+0.017 mm clearance at 0.3 mm).
    exit        0 (the 9 gating checks pass, the gate among them); one reported check does not: the hand estimate within 10 % (above).

Self-check (the answers follow).
  1. The optimized rotor's peak is 39.06 MPa at 10,000 rpm. At what speed does it reach 150 MPa, yield / 2, and why does one solve give the answer?
  2. The optimizer thinned the rib from 3.3 to 1.0 mm. Why did the verdict peak rise by only 17 %, from 33.5 to 39.1 MPa?
  3. The rib sweep moves mag_y and re-solves mag_x. Why does the bridge sweep not need mag_y to move, and why does the rib not change in it?
  4. The three meshes give 31.82, 33.45 and 33.85 MPa. Why is the estimate for an infinitely fine mesh 33.98 and not 33.85, and why does it not exist for the optimized rotor's 38.96, 39.06, 39.64?
  5. Why is the nominal stress of a strip F / (width x stack), and what does the peak/nominal factor tell the designer?
Answers.
  1. 10,000 x sqrt(150 / 39.06) = 19,596 rpm. Every load is proportional to w^2 (rho w^2 r on the steel, m w^2 r_c on each magnet) and the model is linear (K u = f), so the displacement and
     every stress scale with w^2 and a stress that is 3.84 times larger needs a speed sqrt(3.84) = 1.96 times larger. It holds while nothing else loads the rotor: no press fit, no heat, no
     magnet that lifts off its wall.
  2. The rib's nominal stress rose 2.5 times (9.34 to 23.25 MPa: 465 N now on 1.0 mm, 616 N then on 3.3 mm), but its peak/nominal factor fell from 3.05 (28.4 MPa in the old rib) to 1.68:
     the same 0.3 mm fillet is a blunter notch on a 1.0 mm strip (radius over width 0.30, not 0.09), so it concentrates less. The peak that was 33.5 MPa was in the bridges, 12 times their
     nominal stress; it is set by the load on the slot's end beside the fillet and the fillet itself, hardly by the bridge's width (34.7 MPa on a 1.4 mm bridge, 33.9 on a 2.0 mm one).
  3. The rib is twice the height of the slot's lowest point above the d-axis, and moving a slot along x does not change that height: so mag_x alone gives the bridge sweep its bridges and
     leaves the rib at 1.0004 mm. Moving mag_y changes the rib by twice as much (1 mm of rib for 0.5 mm of mag_y) and also changes how far the slot reaches toward the rim, so mag_x
     has to be solved again to keep the bridge at 1.4 mm.
  4. P* = P3 + (P3 - P2) / (2^p - 1) adds what the finer meshes would still change: if each halving of h cuts the change by 2^p (1.63 to 0.40 MPa is a factor 4.1, p = 2.04), the changes still to
     come add up to 0.40 / (2^p - 1) = 0.13 MPa. For the optimized rotor the changes grow (0.10, then 0.58 MPa), which is what a sequence does that has not yet reached the regime where
     the error falls as h^p (the rib is 1.0 mm wide and its fillets 0.3 mm: the coarser meshes have few elements in them); an order read from it would mean nothing. The finest mesh's
     value is then the best one, and the true peak may be a little higher still.
  5. F is in newtons and the cut's area is its width (mm) times the stack (mm), so F / (width x stack) is in N/mm^2, which is MPa: the stress that the strip would have if it carried its
     force evenly. The peak/nominal factor is the stress concentration that the hand estimate lacks (the note's 15 MPa was a nominal stress): 4.67 for the three strips together, 1.68 at the
     rib, 11 to 12 at the bridges, where the peak is set locally (the magnet's end-wall pressure beside a small fillet) and not by the tension through the strip.
"""
from __future__ import annotations

import argparse
import json
import math
import pickle
import shutil
import sys
import time
import traceback
from pathlib import Path

import numpy as np

import bench
import bench.b4_ccx as b4
import femm_queue
import motorgeo as mg
import paths
from optimization import design

CASE = "stress_ccx"
RPM_REF = 10_000                                   # the speed of the solve that the safe speeds are scaled from (B4's first speed)
SPEEDS_RPM = (10_000, 12_000)                      # B4's two speeds
H_MAIN, H_COARSE, H_FINE = b4.H_MIN_MAIN, b4.H_MIN_COARSE, b4.H_MIN_FINE      # mm: the element size at the slot walls (0.1, 0.2, 0.05)
GATE_MPA, GATE_TOLERANCE_PERCENT = 33.45, 0.5      # the brief: B4's peak on the reference rotor (pressure model, 10,000 rpm, h_min 0.1 mm), within 0.5 %
BRIDGE_GUARD_MM = 1.4                              # design's bridge guard (Params.bridge, decision 17 of the Map)
RIB_TARGETS_MM = (1.25, 1.4, 1.5, 2.0, 2.5)        # the rib sweep, after the design's own 1.0 mm (the brief's 1.5, 2.0, 2.5, and 1.25 and 1.4 to find where the bridges take over)
BRIDGE_TARGETS_MM = (1.6, 1.8, 2.0)                # the bridge sweep, after the design's own 1.4 mm
FILLET_NEW_MM = 0.5                                # a design note advises 0.5 mm or more; motorgeo's is 0.3 mm

# The strength limits (04 Structural 6.4: placeholders until M-19's own values are found)
YIELD_MPA = 300.0                                  # the low end of similar grades (300-450 MPa; M-15, the neighbouring grade, has 358)
S_ULTIMATE_MPA = 490.0                             # M-15's tensile strength
FATIGUE_MPA = 196.0                                # 0.4 x M-15's 490 MPa: "about 200", polished and fully reversed (punched edges are lower: S5 sets the reduction)
GOODMAN_MPA = 2.0 / (1.0 / FATIGUE_MPA + 1.0 / S_ULTIMATE_MPA)    # 280 MPa: the largest peak of a zero-to-peak cycle on Goodman's line (the note's 6.1)
SAFETY_FACTORS = (1.5, 2.0)                        # the brief: yield / 1.5 and yield / 2
OVERSPEED = 1.2                                    # the brief: those limits are met at 1.2 x the maximum operating speed
N45_TENSILE_MPA = 78.0                             # 04 Structural 4.2 and 6.1: tension; from Eclipse (2011) and Alliance (2010), which give it for sintered NdFeB generally; Arnold's N45 sheet gives flexure
TIP_SPEED_COMFORT_M_S = 80.0                       # the brief: a laminated 90 mm rotor is comfortable up to about this tip speed
ROSS_FIRST_CRITICAL_RPM = 22_838                   # the brief: bench B6, a 60 mm stack
STRIP_NAMES = (("rib", "centre rib"), ("bridge_low", "outer bridge, lower slot"), ("bridge_high", "outer bridge, upper slot"))


# ============================================================================================== folders and small helpers
def out_folder() -> Path:
    """out\\optimization\\stress\\ (made if it is not there)."""
    return design.out_dir("stress")


def fresh_folder(folder: Path, root: Path) -> Path:
    """Empty `folder` and make it again. Only a folder strictly inside `root` (this script's own output folder) may be emptied: the path is resolved first, so a ".." cannot lead out."""
    target, base = Path(folder).resolve(), Path(root).resolve()
    if base not in target.parents:
        raise ValueError(f"{target} is not inside {base}: refusing to empty it")
    if target.exists():
        shutil.rmtree(target)
    target.mkdir(parents=True)
    return Path(folder)


def slug(text: str) -> str:
    """A folder name from a run's tag: letters, digits and dots, the rest as underscores."""
    return "".join(c if c.isalnum() or c == "." else "_" for c in text).strip("_")


def percent(a: float, b: float) -> float:
    return 100.0 * (a - b) / b


def mean(values) -> float:
    values = list(values)
    return sum(values) / len(values)


def tip_speed_m_s(rpm: float, radius_mm: float) -> float:
    return 2 * math.pi * rpm / 60 * radius_mm * 1e-3


def rpm_at_tip_speed(v_m_s: float, radius_mm: float) -> float:
    return v_m_s / (radius_mm * 1e-3) * 60 / (2 * math.pi)


def load_final() -> dict:
    """The optimized design: out\\optimization\\final\\final.json, whose key "overrides" is design.params's input."""
    return json.loads((paths.path("out") / "optimization" / "final" / "final.json").read_text(encoding="utf-8"))


# ============================================================================================== a geometry, as the study sees it
def corner_clearance(p: mg.Params) -> float:
    """How far the real magnet's sharp corner lies inside the slot's rounded corner [mm]: the fillet's radius less the distance from the fillet's centre to the magnet's corner.
    Positive: the corner fits (0.017 mm at 0.3 mm); negative: it would stand in the steel (-0.066 mm at 0.5 mm), and the magnet's corner would have to be chamfered or rounded."""
    hw, ht, f = p.slot_w / 2, p.slot_t / 2, p.slot_fillet
    return f - math.hypot(hw - f - p.mag_w / 2, ht - f - p.mag_t / 2)


def geometry_of(ov: dict) -> dict:
    """What design.screen says about the candidate `ov`, and the strips it really has: {"feasible", "reason", "rib_mm", "bridge_mm", "q_rib_mm", "mag_x", "mag_y", "v_half_deg",
    "slot_fillet_mm", "corner_clearance_mm"}. The rib is motorgeo.check's (the narrowest steel between a pole's two slots) and the bridge is r_ro less the slot's reach. screen() leaves
    out its metrics when it refuses a candidate before drawing it, so they are taken from motorgeo.check then (a GeometryError, a slot past the guard, is for the caller to meet)."""
    p = design.params(ov)
    s = design.screen(ov)
    metrics = s["metrics"] or mg.check(p, mg.build(p))[0]
    return {"feasible": bool(s["feasible"]), "reason": s["reason"], "rib_mm": metrics["centre_rib"], "bridge_mm": design.bridge_of(p), "q_rib_mm": metrics["q_rib"], "mag_x": p.mag_x,
            "mag_y": p.mag_y, "v_half_deg": p.v_half_deg, "slot_fillet_mm": p.slot_fillet, "corner_clearance_mm": corner_clearance(p)}


def rib_of(ov: dict) -> float:
    p = design.params(ov)
    return mg.check(p, mg.build(p))[0]["centre_rib"]


def candidate(base: dict, bridge_mm: float = BRIDGE_GUARD_MM, **changes) -> dict:
    """`base` with `changes` (mag_y, slot_fillet, ...) and its mag_x re-solved so that the outer bridge is `bridge_mm`."""
    ov = {**base, **changes}
    ov["mag_x"] = design.mag_x_for_bridge(bridge_mm, ov)
    return ov


def with_rib(base: dict, rib_mm: float, bridge_mm: float = BRIDGE_GUARD_MM) -> dict:
    """`base` with mag_y moved until the centre rib is `rib_mm`, and mag_x re-solved so that the bridge is `bridge_mm`. The rib is twice the height of the slot's lowest point above the
    d-axis, and that height is mag_y less a constant, so the slope is 2: each step moves mag_y by half the rib still to gain, and the first step is nearly all of it."""
    y = base["mag_y"]
    ov = candidate(base, bridge_mm, mag_y=y)
    for _ in range(4):
        error = rib_mm - rib_of(ov)
        if abs(error) < 1e-7:
            break
        y += error / 2
        ov = candidate(base, bridge_mm, mag_y=y)
    return ov


# ============================================================================================== one solve, and what is read from it
def say_solve(log: bench.Log, label: str, key: str, run: dict, a: dict) -> None:
    if a["finished"]:
        pl = a["peak"]["place"]
        log.say(f"  {label} {key}: ccx {run['version']} {run['seconds']:.1f} s; peak {a['peak']['value_mpa']:.2f} MPa at ({pl['x_mm']:.2f}, {pl['y_mm']:.2f}) mm: {pl['feature']} of "
                f"{pl['slot']}; bridges {a['bridge_peak_mpa']:.2f} MPa; bore {a['bore']['hoop_mean_mpa']:.3f} MPa (Lame {a['bore']['lame_mpa']:.3f}); "
                f"equilibrium {a['equilibrium']['residual_percent']:.5f} %")
    else:
        log.say(f"  {label} {key}: ccx did not finish (exit {run['returncode']}, timed out {run['timed_out']}); see {run['folder']}")


def peak_record(a: dict) -> dict:
    """The verdict peak of an analysis, flat: value, node, place and the corner of the slot."""
    pk, pl = a["peak"], a["peak"]["place"]
    return {"value_mpa": pk["value_mpa"], "node": pk["node"], "x_mm": pl["x_mm"], "y_mm": pl["y_mm"], "r_mm": pl["r_mm"], "angle_deg": pl["angle_deg"], "feature": pl["feature"],
            "slot": pl["slot"], "pole": pl["pole"], "strip": pl["strip"], "distance_to_strip_mm": pl["distance_to_strip_mm"], "detail": pl.get("slot_detail")}


def describe_peak(rec: dict) -> str:
    return (f"{rec['feature']} of {rec['slot']} (pole {rec['pole']}), {rec['detail']}; at ({rec['x_mm']:.2f}, {rec['y_mm']:.2f}) mm, r {rec['r_mm']:.2f} mm, {rec['angle_deg']:.2f} deg, "
            f"{rec['distance_to_strip_mm']:.3f} mm from the strip's cut")


def zone_groups(a: dict) -> dict:
    """The largest von Mises stress in each strip's zone: {"bridges": {slot: MPa}, "ribs": {rib_k: MPa}}."""
    zones = a["zone_peaks"]
    return {"bridges": {k: v["value"] for k, v in zones.items() if k.startswith("magnet_")}, "ribs": {k: v["value"] for k, v in zones.items() if k.startswith("rib")}}


def principal_max(components) -> float:
    """The largest principal stress of a plane-stress state (sxx, syy, szz, sxy, ...): the tensile stress a brittle magnet fails by (von Mises overstates it where the state is compressive)."""
    sx, sy, _, txy, _, _ = components
    return 0.5 * (sx + sy) + math.hypot(0.5 * (sx - sy), txy)


def magnet_body_stress(mesh: b4.Mesh, stress: dict) -> dict:
    """The bonded model's stress at the nodes of the slot bodies (every node that no steel triangle uses: the interior ones and the copies of the wall nodes): the largest von Mises
    (B4's number) and the largest principal stress (the brittle criterion), in MPa."""
    steel = {n for t in mesh.steel for n in t}
    body = [s for n, s in stress.items() if n not in steel]
    return {"von_mises_mpa": max(b4.von_mises(s) for s in body), "principal_max_mpa": max(principal_max(s) for s in body), "nodes": len(body)}


def strip_table(p: mg.Params, a: dict) -> dict:
    """The three strips of the verdict run as the case sees them: per strip its width, the force through its cut (the mean of the two poles: they are mirror images), its share of the
    total and its nominal stress F / (width x stack); the total carried and its nominal stress over the three widths; the hand estimate's pull and nominal stress on the same basis
    (the pull along the pole's axis, which is what the strips balance) and by the note's own formula; and the peak/nominal factors: against the mean nominal stress (B4's), against the
    hand estimate's, and against the nominal stress of the strip the peak lies in."""
    sections = a["sections"]
    hand = b4.hand_estimate(p, RPM_REF)
    strips = {}
    for key, name in STRIP_NAMES:
        force = mean(s["strips"][key]["axial_pull_n"] for s in sections)
        width = sections[0]["strips"][key]["width_mm"]
        strips[key] = {"name": name, "width_mm": width, "force_n": force, "forces_by_pole_n": [s["strips"][key]["axial_pull_n"] for s in sections], "nominal_mpa": force / (p.stack * width)}
    carried = mean(s["axial_pull_n"] for s in sections)
    total_width = sum(v["width_mm"] for v in strips.values())
    for v in strips.values():
        v["share_percent"] = 100 * v["force_n"] / carried
    peak = a["peak"]["value_mpa"]
    nominal_mean = carried / p.stack / total_width
    hand_axis = mean(s["hand_pull_n"] for s in sections)
    nominal_hand = hand_axis / p.stack / total_width
    place = a["peak"]["place"]
    local = strips[place["strip"]]["nominal_mpa"] if place["in_a_strip"] else None
    return {"strips": strips, "carried_n": carried, "carried_by_pole_n": [s["axial_pull_n"] for s in sections], "total_width_mm": total_width, "nominal_mean_mpa": nominal_mean,
            "hand_pull_axis_n": hand_axis, "hand_pull_by_pole_n": [s["hand_pull_n"] for s in sections], "hand_pull_note_n": hand["pull_n_over_stack"],
            "refined_pull_by_pole_n": [s["refined_pull_n"] for s in sections], "hand_nominal_mpa": nominal_hand, "hand_nominal_note_mpa": hand["nominal_mpa"],
            "against_hand_percent": [percent(s["axial_pull_n"], s["hand_pull_n"]) for s in sections], "against_refined_percent": [percent(s["axial_pull_n"], s["refined_pull_n"]) for s in sections],
            "cap_area_mm2": hand["cap_area_mm2"], "magnet_share_of_pull": hand["magnet_share_of_pull"], "transverse_n": [s["transverse_n"] for s in sections],
            "peak_over_nominal_mean": peak / nominal_mean, "peak_over_nominal_hand": peak / nominal_hand, "peak_over_nominal_local": peak / local if local else None,
            "peak_strip": place["strip"] if place["in_a_strip"] else None}


def richardson(peaks: dict) -> dict:
    """The peak on the three meshes of the convergence study (h_min 0.2, 0.1, 0.05 mm: each half the one before) and what they say of an infinitely fine mesh.
    p = log2((P1 - P2) / (P2 - P3)) is the order the error falls with (it needs P1, P2, P3 to approach a limit: monotone, each change smaller than the last); P* = P3 + (P3 - P2) / (2^p - 1).
    Also P* for p = 2 (the order of quadratic triangles on a smooth field, and the one B4 read), the main mesh's (P2) departure from P*, and from P3."""
    h = sorted(peaks, reverse=True)
    p1, p2, p3 = (peaks[x] for x in h)
    d1, d2 = p2 - p1, p3 - p2
    out = {"h_mm": h, "peaks_mpa": [p1, p2, p3], "changes_mpa": [d1, d2], "main_against_finest_percent": percent(p2, p3), "observed_order": None, "extrapolated_mpa": None,
           "extrapolated_p2_mpa": p3 + d2 / (2 ** 2 - 1), "main_against_extrapolated_percent": None, "main_against_extrapolated_p2_percent": percent(p2, p3 + d2 / 3), "usable": False}
    if d1 * d2 > 0 and abs(d2) < abs(d1):
        order = math.log2(d1 / d2)
        out["observed_order"] = order
        if 0.5 <= order <= 4.0:
            out["extrapolated_mpa"] = p3 + d2 / (2 ** order - 1)
            out["main_against_extrapolated_percent"] = percent(p2, out["extrapolated_mpa"])
            out["usable"] = True
    return out


def safe_speeds(peak_mpa: float, radius_mm: float, rpm_ref: float = RPM_REF) -> list[dict]:
    """The speeds at which the verdict peak (peak_mpa at rpm_ref) reaches each limit: n = n_ref x sqrt(limit / peak). (a) the yield; (b) yield / 1.5 and yield / 2 at a 1.2 x overspeed:
    the maximum operating speed is the speed at which the peak reaches the limit, divided by 1.2; (c) the fatigue placeholder and Goodman's zero-to-peak limit. Each row: its name, limit
    [MPa], the speed [rpm], the maximum operating speed [rpm] (for (b), else None) and the tip speeds [m/s]."""
    rows = [("(a) yield", YIELD_MPA, False)]
    rows += [(f"(b) yield / {sf:g}, at {OVERSPEED:g} x overspeed", YIELD_MPA / sf, True) for sf in SAFETY_FACTORS]
    rows += [("(c) fatigue placeholder (196 MPa)", FATIGUE_MPA, False), ("(c) Goodman, zero-to-peak (280 MPa)", GOODMAN_MPA, False)]
    out = []
    for name, limit, operating in rows:
        n = rpm_ref * math.sqrt(limit / peak_mpa)
        n_op = n / OVERSPEED if operating else None
        out.append({"name": name, "limit_mpa": limit, "speed_rpm": n, "operating_rpm": n_op, "tip_speed_m_s": tip_speed_m_s(n, radius_mm),
                    "operating_tip_speed_m_s": tip_speed_m_s(n_op, radius_mm) if n_op else None})
    return out


# ============================================================================================== one rotor: the checks of B4, on the rotor of `p`
def study_rotor(label: str, p: mg.Params, folder: Path, root: Path, log: bench.Log, first: dict | None = None) -> dict:
    """B4's case on the rotor of `p`: the pressure and bonded models at 10,000 and 12,000 rpm on the main mesh, the verdict model on the coarser and the finer mesh, and everything that
    is read from them. `first`, if given, is a solve already made on the main mesh ({"model", "key", "run", "analysis"}: the regression gate's), which is not made again.
    Returns {"data": every number (JSON-able), "model", "stress": the verdict run's, "checks": [...]}. A solve that does not finish raises RuntimeError: nothing after it means anything."""
    log.say(f"--- {label} rotor: the main mesh, h_min {H_MAIN:g} mm")
    main = first["model"] if first else b4.build_model(p, H_MAIN)
    q = main["mesh"].quality
    log.say(f"  {label}: {len(main['mesh'].nodes):,} nodes ({b4.steel_node_count(main['mesh']):,} in the steel), {len(main['mesh'].steel):,} steel triangles, min SICN {q['min_sicn']:.3f}")
    runs: dict = {}
    analysis: dict = {}
    for variant in ("pressure", "bonded"):
        for rpm in SPEEDS_RPM:
            key = f"{variant}_{rpm}_h{H_MAIN:g}"
            if first and key == first["key"]:
                runs[key], analysis[key] = first["run"], first["analysis"]
                say_solve(log, label, key + " (the gate's)", runs[key], analysis[key])
                continue
            runs[key] = b4.solve(main, variant, rpm, fresh_folder(folder / key, root), log, f"stress_ccx {label} {variant} {rpm} rpm h {H_MAIN:g}")
            analysis[key] = b4.analyse(main, runs[key], rpm)
            say_solve(log, label, key, runs[key], analysis[key])
            if not analysis[key]["finished"]:
                raise RuntimeError(f"{label} {key}: ccx did not finish")
    verdict_key = f"pressure_{RPM_REF}_h{H_MAIN:g}"
    bonded_key = f"bonded_{RPM_REF}_h{H_MAIN:g}"
    bodies = magnet_body_stress(main["mesh"], runs[bonded_key]["stress"])
    # ---- the mesh: the verdict model at 10,000 rpm on the coarser and the finer mesh
    levels = {H_MAIN: {"key": verdict_key, "nodes": b4.steel_node_count(main["mesh"]), "peak_mpa": analysis[verdict_key]["peak"]["value_mpa"], "peak": peak_record(analysis[verdict_key])}}
    for h in (H_COARSE, H_FINE):
        other = b4.build_model(p, h)
        key = f"pressure_{RPM_REF}_h{h:g}"
        runs[key] = b4.solve(other, "pressure", RPM_REF, fresh_folder(folder / key, root), log, f"stress_ccx {label} pressure {RPM_REF} rpm h {h:g}")
        analysis[key] = b4.analyse(other, runs[key], RPM_REF)
        say_solve(log, label, key, runs[key], analysis[key])
        if not analysis[key]["finished"]:
            raise RuntimeError(f"{label} {key}: ccx did not finish")
        levels[h] = {"key": key, "nodes": b4.steel_node_count(other["mesh"]), "peak_mpa": analysis[key]["peak"]["value_mpa"], "peak": peak_record(analysis[key])}
        runs[key]["stress"] = None                          # the finer meshes' fields are not kept
    conv = richardson({h: v["peak_mpa"] for h, v in levels.items()})
    verdict, bonded = analysis[verdict_key], analysis[bonded_key]
    verdict_stress = runs[verdict_key]["stress"]
    runs[bonded_key]["stress"] = None
    for key in (f"pressure_{SPEEDS_RPM[1]}_h{H_MAIN:g}", f"bonded_{SPEEDS_RPM[1]}_h{H_MAIN:g}"):
        runs[key]["stress"] = None
    peak = verdict["peak"]["value_mpa"]
    table = strip_table(p, verdict)
    speeds = safe_speeds(peak, p.r_ro)
    hand10 = b4.hand_estimate(p, RPM_REF)
    loads = b4.loads_at(main, RPM_REF)["loads"]
    first_slot = sorted(loads)[0]
    data = {"label": label, "geometry": geometry_of({k: getattr(p, k) for k in design.ALL_VARS} | {"slot_fillet": p.slot_fillet}),
            "mesh": b4.mesh_counts(main["mesh"]),
            "walls": {n: {k: {x: w[x] for x in ("length_mm", "force_n", "pressure_mpa")} for k, w in walls.items()} for n, walls in loads.items()},
            "wall_example": {"slot": first_slot, **{k: {x: w[x] for x in ("length_mm", "force_n", "pressure_mpa")} for k, w in loads[first_slot].items()}},
            "hand": {"cap_area_mm2": hand10["cap_area_mm2"], "pull_n_over_stack": hand10["pull_n_over_stack"], "nominal_mpa": hand10["nominal_mpa"], "strips_mm": hand10["strips_mm"],
                     "magnet_share_of_pull": hand10["magnet_share_of_pull"]},
            "verdict_key": verdict_key, "peak": peak_record(verdict), "peak_12000_mpa": analysis[f"pressure_{SPEEDS_RPM[1]}_h{H_MAIN:g}"]["peak"]["value_mpa"],
            "bridge_peak_mpa": verdict["bridge_peak_mpa"], "bridge_spread_percent": verdict["bridge_spread_percent"], "zones": zone_groups(verdict),
            "bore": {"pressure": verdict["bore"], "bonded": bonded["bore"]},
            "bonded": {"peak": peak_record(bonded), "bridge_peak_mpa": bonded["bridge_peak_mpa"], "zones": zone_groups(bonded), "magnet_bodies": bodies},
            "strips": table, "convergence": {**conv, "levels": {str(h): v for h, v in levels.items()}}, "speeds": speeds,
            "runs": {k: {"folder": r["folder"], "seconds": r["seconds"], "finished": r["finished"], "version": r["version"], "h_min": r["h_min"], "analysis": analysis[k]} for k, r in runs.items()}}
    checks = rotor_checks(label, p, data, analysis)
    return {"data": data, "model": main, "stress": verdict_stress, "checks": checks}


def rotor_checks(label: str, p: mg.Params, d: dict, analysis: dict) -> list[dict]:
    """B4's checks, each against B4's own tolerance (imported from b4_ccx, never widened here), on the rotor `label`. The gating ones are the solver's own (the solves end, equilibrium, the
    same load in both magnet models, the two speeds scale): a miss there means the numbers cannot be trusted. The others are reported with their verdict: on the optimized rotor they
    describe the design, and a miss is a finding, not a fault."""
    out: list[dict] = []
    runs = d["runs"]
    finished = [k for k, r in runs.items() if r["finished"]]
    out.append(bench.check(f"{label}: every solve ends normally", len(finished) == len(runs), "exit 0, 'Job finished', the .frd written, inside the timeout, for every solve",
                           f"{len(finished)} of {len(runs)} solves: " + ", ".join(f"{k} {r['seconds']:.1f} s" for k, r in runs.items()), "ccx_static.exe through run_solver"))
    residuals = {k: a["equilibrium"]["residual_percent"] for k, a in analysis.items()}
    out.append(bench.check(f"{label}: equilibrium (the reactions of the symmetry edges balance the load worked out from the mesh)", max(residuals.values()) <= b4.EQUILIBRIUM_TOLERANCE_PERCENT,
                           f"within {b4.EQUILIBRIUM_TOLERANCE_PERCENT} % for every solve (B4's own limit)", ", ".join(f"{k} {v:.5f} %" for k, v in residuals.items()), "the .dat file's reactions, the mesh's load"))
    same = [percent(analysis[f"bonded_{rpm}_h{H_MAIN:g}"]["total_load_n"], analysis[f"pressure_{rpm}_h{H_MAIN:g}"]["total_load_n"]) for rpm in SPEEDS_RPM]
    out.append(bench.check(f"{label}: the two magnet models carry the same total load", max(abs(v) for v in same) <= b4.SAME_LOAD_TOLERANCE_PERCENT,
                           f"within {b4.SAME_LOAD_TOLERANCE_PERCENT} % at both speeds (B4's own limit): they differ only in how the load reaches the steel",
                           ", ".join(f"{rpm:,} rpm {v:+.5f} %" for rpm, v in zip(SPEEDS_RPM, same)), "the decks' loads"))
    for variant, gating in (("pressure", True), ("bonded", False)):
        low, high = analysis[f"{variant}_{SPEEDS_RPM[0]}_h{H_MAIN:g}"], analysis[f"{variant}_{SPEEDS_RPM[1]}_h{H_MAIN:g}"]
        extra = {} if gating else {"gating": False, "reason": "the verdict is the pressure model's"}
        out.append(b4.scale_check(f"{label} {variant}", low["peak"]["value_mpa"], high["peak"]["value_mpa"], **extra))
    out.append(bench.check(f"{label}: the four slots' bridges have the same peak (mirror images, meshed differently)", d["bridge_spread_percent"] <= b4.SYMMETRY_TOLERANCE_PERCENT,
                           f"spread within {b4.SYMMETRY_TOLERANCE_PERCENT} % of their mean (B4's own limit)",
                           ", ".join(f"{k} {v:.2f}" for k, v in d["zones"]["bridges"].items()) + f" MPa: spread {d['bridge_spread_percent']:.2f} %; the two ribs "
                           + ", ".join(f"{v:.2f}" for v in d["zones"]["ribs"].values()) + " MPa", "the rotor's symmetry", gating=False, reason="a finding on the design, not a solver fault"))
    conv = d["convergence"]
    out.append(bench.check(f"{label}: the peak does not depend on the mesh", abs(conv["main_against_finest_percent"]) <= b4.MESH_TOLERANCE_PERCENT,
                           f"the peak on the main mesh within {b4.MESH_TOLERANCE_PERCENT} % of the finest (B4's own limit)",
                           f"h_min 0.2, 0.1, 0.05 mm: " + ", ".join(f"{v:.2f}" for v in conv["peaks_mpa"]) + f" MPa; main against finest {conv['main_against_finest_percent']:+.2f} %",
                           "the same model on three meshes", gating=False, reason="a finding on the design, not a solver fault"))
    st = d["strips"]
    diffs, refined = st["against_hand_percent"], st["against_refined_percent"]
    out.append(bench.check(f"{label}: the force the three strips carry agrees with the hand estimate", all(abs(x) <= b4.HAND_TOLERANCE_PERCENT for x in diffs),
                           f"within {b4.HAND_TOLERANCE_PERCENT} % for both poles (B4's own limit; the hand method puts all of the magnets' pull on the cap)",
                           f"poles 0 and 1: {st['carried_by_pole_n'][0]:.1f} N and {st['carried_by_pole_n'][1]:.1f} N against {st['hand_pull_by_pole_n'][0]:.1f} N and "
                           f"{st['hand_pull_by_pole_n'][1]:.1f} N ({diffs[0]:+.2f} %, {diffs[1]:+.2f} %)", "04 Structural 6.1, strip_geometry, hand_estimate",
                           gating=False, reason="the hand method's error depends on the magnets' share of the pull"))
    out.append(bench.check(f"{label}: ... and with the magnets' end-wall load left off the cap", all(abs(x) <= b4.REFINED_TOLERANCE_PERCENT for x in refined),
                           f"within {b4.REFINED_TOLERANCE_PERCENT} % for both poles (B4's own limit)",
                           f"{st['carried_by_pole_n'][0]:.1f} N and {st['carried_by_pole_n'][1]:.1f} N against {st['refined_pull_by_pole_n'][0]:.1f} N and {st['refined_pull_by_pole_n'][1]:.1f} N "
                           f"({refined[0]:+.2f} %, {refined[1]:+.2f} %)", "the free body of the cap", gating=False, reason="the cut forces against the model's own free body"))
    bb = d["bore"]["bonded"]
    out.append(bench.compare(f"{label} bonded: the bore's mean hoop stress against Lame's plain disc", bb["hoop_mean_mpa"], bb["lame_mpa"], b4.LAME_TOLERANCE_PERCENT,
                             "Lame's disc, 7.5 to 45 mm (B4's own limit: a gate on units, speed and load)", "MPa", gating=False, reason="a diagnostic of B4, reported"))
    bodies = d["bonded"]["magnet_bodies"]
    out.append(bench.check(f"{label} bonded: the largest stress in a magnet body is below N45's tensile strength", bodies["principal_max_mpa"] < N45_TENSILE_MPA and bodies["von_mises_mpa"] < N45_TENSILE_MPA,
                           f"below {N45_TENSILE_MPA:g} MPa (04 Structural 4.2: Eclipse and Alliance, for sintered NdFeB)",
                           f"largest von Mises {bodies['von_mises_mpa']:.2f} MPa, largest principal (tensile) {bodies['principal_max_mpa']:.2f} MPa over {bodies['nodes']:,} body nodes",
                           "the bonded deck's nodal stress", gating=False, reason="the bonded model is the other bound; a brittle magnet is judged by its principal stress"))
    return out


# ============================================================================================== the regression gate
def gate(p_ref: mg.Params, folder: Path, root: Path, log: bench.Log) -> dict:
    """B4 on the reference rotor (pressure model, 10,000 rpm, h_min 0.1 mm): its peak is the brief's 33.45 MPa within 0.5 %. Returns {"model", "key", "run", "analysis", "check"}."""
    key = f"pressure_{RPM_REF}_h{H_MAIN:g}"
    model = b4.build_model(p_ref, H_MAIN)
    run = b4.solve(model, "pressure", RPM_REF, fresh_folder(folder / key, root), log, f"stress_ccx gate: reference pressure {RPM_REF} rpm h {H_MAIN:g}")
    a = b4.analyse(model, run, RPM_REF)
    say_solve(log, "reference (the gate)", key, run, a)
    if not a["finished"]:
        return {"model": model, "key": key, "run": run, "analysis": a, "check": bench.check("regression gate: B4 on the reference rotor", False, f"ccx ends and the peak is {GATE_MPA} MPa within "
                f"{GATE_TOLERANCE_PERCENT} %", "ccx did not finish", "bench/b4_ccx.py")}
    peak = a["peak"]["value_mpa"]
    item = bench.compare("regression gate: B4's peak on the reference rotor (pressure model, 10,000 rpm, h_min 0.1 mm)", peak, GATE_MPA, GATE_TOLERANCE_PERCENT,
                         "the brief: B4's 33.45 MPa (bench/b4_ccx.py, 2026-11-06)", "MPa")
    return {"model": model, "key": key, "run": run, "analysis": a, "check": item}


# ============================================================================================== the sensitivities
def solve_geometry(tag: str, ov: dict, folder: Path, root: Path, log: bench.Log, h_min: float = H_MAIN) -> dict:
    """One geometry of a sweep: screen it, mesh it, solve the verdict model at 10,000 rpm, read it. Returns the row: {"tag", "overrides", "geometry", "peak", "zones", ...}; a geometry that cannot
    be drawn or solved comes back with "error" and what is known."""
    row: dict = {"tag": tag, "overrides": ov, "h_min": h_min}
    try:
        row["geometry"] = geometry_of(ov)
        p = design.params(ov)
        model = b4.build_model(p, h_min)
        run = b4.solve(model, "pressure", RPM_REF, fresh_folder(folder / slug(tag), root), log, f"stress_ccx {tag} pressure {RPM_REF} rpm h {h_min:g}")
        a = b4.analyse(model, run, RPM_REF)
        say_solve(log, tag, f"h {h_min:g}", run, a)
        if not a["finished"]:
            row["error"] = "ccx did not finish"
            return row
        row.update(sweep_numbers(p, model, a, run))
        run["stress"] = None
    except Exception as error:                              # one bad geometry must not end the sweep
        row["error"] = f"{type(error).__name__}: {error}"
        log.say(f"  {tag}: {row['error']}")
        log.say(traceback.format_exc(limit=3))
    return row


def sweep_numbers(p: mg.Params, model: dict, a: dict, run: dict | None = None) -> dict:
    """What a sweep keeps of a solve: the peak and its place, the largest stress in the rib's and the bridges' zones, the rib's share of the force the strips carry, the walls' pressures."""
    zones = zone_groups(a)
    sections = a["sections"]
    carried = mean(s["axial_pull_n"] for s in sections)
    rib_force = mean(s["strips"]["rib"]["axial_pull_n"] for s in sections)
    walls = b4.loads_at(model, RPM_REF)["loads"]["magnet_0_dn"]
    return {"peak": peak_record(a), "rib_zone_mpa": max(zones["ribs"].values()), "bridge_zone_mpa": a["bridge_peak_mpa"], "zones": zones, "carried_n": carried, "rib_force_n": rib_force,
            "rib_share_percent": 100 * rib_force / carried, "equilibrium_percent": a["equilibrium"]["residual_percent"], "nodes": b4.steel_node_count(model["mesh"]),
            "walls": {k: {"length_mm": w["length_mm"], "pressure_mpa": w["pressure_mpa"]} for k, w in walls.items()}, "seconds": run["seconds"] if run else None}


def design_row(opt: dict) -> dict:
    """The optimized design as the first row of the rib and bridge sweeps: its own verdict solve, read the way every other row is."""
    d = opt["data"]
    a = d["runs"][d["verdict_key"]]["analysis"]
    p = opt["model"]["p"]
    row = {"tag": "design", "overrides": {k: getattr(p, k) for k in design.ALL_VARS}, "h_min": H_MAIN, "geometry": d["geometry"]}
    row.update(sweep_numbers(p, opt["model"], a, d["runs"][d["verdict_key"]]))
    return row


def sensitivities(opt: dict, base: dict, folder: Path, root: Path, log: bench.Log) -> dict:
    """The rib sweep, the bridge sweep and the corner radius, on the optimized design `base` (design.params's overrides)."""
    log.say("--- sensitivities on the optimized design (pressure model, 10,000 rpm)")
    first = design_row(opt)
    rib_rows, bridge_rows = [first], [first]
    for rib in RIB_TARGETS_MM:
        ov = with_rib(base, rib)
        rib_rows.append(solve_geometry(f"rib {rib:g}", ov, folder, root, log))
    for bridge in BRIDGE_TARGETS_MM:
        ov = candidate(base, bridge)
        bridge_rows.append(solve_geometry(f"bridge {bridge:g}", ov, folder, root, log))
    # the corner radius: (a) the slot where it is, (b) the rib and the bridge as in the design, and (b) on the finest mesh
    rib_design = first["geometry"]["rib_mm"]
    fillet_a = {**base, "slot_fillet": FILLET_NEW_MM}
    fillet_b = with_rib({**base, "slot_fillet": FILLET_NEW_MM}, rib_design)
    rows_a = solve_geometry(f"fillet {FILLET_NEW_MM:g} (a) slot unmoved", fillet_a, folder, root, log)
    rows_b = solve_geometry(f"fillet {FILLET_NEW_MM:g} (b) same rib and bridge", fillet_b, folder, root, log)
    rows_b_fine = solve_geometry(f"fillet {FILLET_NEW_MM:g} (b) same rib and bridge, fine", fillet_b, folder, root, log, h_min=H_FINE)
    levels = opt["data"]["convergence"]["levels"]
    fine_design = levels[f"{H_FINE:g}"]
    return {"rib": rib_rows, "bridge": bridge_rows,
            "fillet": {"design": first, "a": rows_a, "b": rows_b, "b_fine": rows_b_fine, "design_fine_peak_mpa": fine_design["peak_mpa"], "design_fine_peak": fine_design["peak"],
                       "design_extrapolated_mpa": opt["data"]["convergence"]["extrapolated_mpa"]}}


# ============================================================================================== the pictures
def rotated(x, y, degrees: float):
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    return c * np.asarray(x) - s * np.asarray(y), s * np.asarray(x) + c * np.asarray(y)


def field_of(model: dict, stress: dict) -> dict:
    """The steel's von Mises stress at its nodes, on the corner triangles of the mesh: {"x", "y", "v", "tri"}."""
    mesh = model["mesh"]
    ids = sorted({n for t in mesh.steel for n in t if n in stress})
    index = {n: i for i, n in enumerate(ids)}
    return {"x": np.array([mesh.nodes[n][0] for n in ids]), "y": np.array([mesh.nodes[n][1] for n in ids]), "v": np.array([b4.von_mises(stress[n]) for n in ids]),
            "tri": np.array([[index[n] for n in t[:3]] for t in mesh.steel])}


def draw_view(ax, model: dict, field: dict, vmax: float, view: str, peak: dict, title: str) -> None:
    """One picture of the rotor. view "full": the whole rotor, the 90-degree sector repeated four times. view "pole": the pole the peak lies in, turned so that its axis runs along +x, with the
    three cuts of its strips (the rib and the two bridges) as black lines. The slots are the grey voids (the loose magnets sit in them); the peak is the black ring."""
    import matplotlib.tri as mtri
    p = model["p"]
    slots = [poly for name, poly in mg.visible(model["regions"]).items() if name.startswith("magnet_")]
    if view == "full":
        turns, tilt, geom = (0, 90, 180, 270), 0.0, None
    else:
        geom = b4.strip_geometry(p, peak["pole"])
        turns, tilt = (0,), -geom["axis_deg"]
    for turn in turns:
        x, y = rotated(field["x"], field["y"], turn + tilt)
        ax.tripcolor(mtri.Triangulation(x, y, field["tri"]), field["v"], shading="gouraud", cmap="turbo", vmin=0, vmax=vmax, rasterized=True)
        for poly in slots:
            sx, sy = rotated(*(np.array(c) for c in poly.exterior.xy), turn + tilt)
            ax.fill(sx, sy, facecolor="#dcdcdc", edgecolor="#444444", linewidth=0.4)
    if geom is not None:
        for line in geom["cuts"].values():
            lx, ly = rotated(*(np.array(c) for c in line.xy), tilt)
            ax.plot(lx, ly, color="black", lw=1.0)
    px, py = rotated([peak["x_mm"]], [peak["y_mm"]], tilt)
    ax.plot(px, py, "o", mfc="none", mec="black", ms=15 if view == "pole" else 12, mew=1.8)
    ax.set_aspect("equal")
    if view == "full":
        ax.set_xlim(-p.r_ro - 2, p.r_ro + 2)
        ax.set_ylim(-p.r_ro - 2, p.r_ro + 2)
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
    else:
        ax.set_xlim(27, 47.5)
        ax.set_ylim(-12.5, 12.5)
        ax.set_xlabel("distance along the pole's axis [mm]")
        ax.set_ylabel("distance across the pole [mm]")
    ax.set_title(bench.wrap_title(title), fontsize=9.5, loc="left")


def peak_label(rec: dict) -> str:
    return f"peak {rec['value_mpa']:.1f} MPa in the {rec['feature']} ({rec['detail']})"


def fig_optimized(model: dict, stress: dict, d: dict, target: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    rec, p = d["peak"], model["p"]
    f = field_of(model, stress)
    fig, axes = plt.subplots(1, 2, figsize=(14.5, 7.2), constrained_layout=True)
    draw_view(axes[0], model, f, rec["value_mpa"], "full", rec, "the whole rotor (the 90-degree sector repeated four times); grey: the slots, with the magnets loose in them")
    draw_view(axes[1], model, f, rec["value_mpa"], "pole", rec, f"the pole of the peak, turned to its axis; black lines: the rib's and the two bridges' cuts; ring: the peak, {rec['value_mpa']:.1f} MPa, "
              f"{rec['detail']}; white: outside the rotor and its sector")
    fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0, rec["value_mpa"]), cmap="turbo"), ax=axes, shrink=0.8, pad=0.02).set_label("von Mises stress [MPa]")
    fig.suptitle(bench.wrap_title(f"Optimized rotor at {RPM_REF:,} rpm: von Mises stress, magnets loose (pressure model), CalculiX plane stress, h_min {H_MAIN:g} mm; centre rib "
                                  f"{d['geometry']['rib_mm']:.2f} mm, outer bridges {d['geometry']['bridge_mm']:.2f} mm"), fontsize=11)
    fig.savefig(target, dpi=140, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def fig_compare(ref: dict, opt: dict, target: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    vmax = max(ref["data"]["peak"]["value_mpa"], opt["data"]["peak"]["value_mpa"])
    fig, axes = plt.subplots(2, 2, figsize=(14.5, 14), constrained_layout=True)
    for column, (name, R) in enumerate((("Reference rotor", ref), ("Optimized rotor", opt))):
        d, model = R["data"], R["model"]
        f = field_of(model, R["stress"])
        title = f"{name}: rib {d['geometry']['rib_mm']:.2f} mm, bridges {d['geometry']['bridge_mm']:.2f} mm; {peak_label(d['peak'])}"
        draw_view(axes[0, column], model, f, vmax, "full", d["peak"], title)
        draw_view(axes[1, column], model, f, vmax, "pole", d["peak"], f"{name}, the pole of the peak, turned to its axis (peak {d['peak']['value_mpa']:.1f} MPa)")
    fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(0, vmax), cmap="turbo"), ax=axes, shrink=0.6, pad=0.02).set_label("von Mises stress [MPa], one scale for both rotors")
    fig.suptitle(bench.wrap_title(f"Reference and optimized rotor side by side at {RPM_REF:,} rpm: von Mises stress, magnets loose (pressure model), CalculiX plane stress, h_min {H_MAIN:g} mm\n"
                                  f"rings: the peaks; black lines: the three strips' cuts (the rib and the two bridges); white: outside the rotor and its sector"), fontsize=11)
    fig.savefig(target, dpi=120, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


def fig_sweep(rows: list[dict], which: str, target: Path) -> None:
    """The verdict peak against the rib (which="rib") or against the outer bridge (which="bridge"), with the largest stress of the rib's zone and of the bridges' zones, and the strip the
    verdict peak sits in as the marker's shape. The first row is the design itself."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    good = [r for r in rows if "peak" in r]
    x = [r["geometry"]["rib_mm" if which == "rib" else "bridge_mm"] for r in good]
    fig, ax = plt.subplots(figsize=(8.6, 5.6), constrained_layout=True)
    ax.plot(x, [r["peak"]["value_mpa"] for r in good], "-", color="#bbbbbb", lw=5.0, label="verdict peak (largest in the steel)", zorder=1)
    ax.plot(x, [r["bridge_zone_mpa"] for r in good], "--", color="#2a7fb8", marker="v", label="largest stress in the outer bridges' zones", zorder=3)
    ax.plot(x, [r["rib_zone_mpa"] for r in good], "--", color="#c0504d", marker="^", label="largest stress in the centre rib's zone", zorder=3)
    for xi, r in zip(x, good):
        ax.plot(xi, r["peak"]["value_mpa"], "o" if r["peak"]["feature"] == "outer bridge" else "s", color="black", ms=9, zorder=4)
        ax.annotate(f"{r['peak']['value_mpa']:.1f}", (xi, r["peak"]["value_mpa"]), textcoords="offset points", xytext=(0, 11), ha="center", fontsize=9)
    ax.plot(x[0], good[0]["peak"]["value_mpa"], "*", color="gold", mec="black", ms=17, zorder=5, label="the optimized design")
    ax.set_xlabel(f"{'centre rib' if which == 'rib' else 'outer bridge'} width [mm]" + ("  (outer bridge held at 1.40 mm)" if which == "rib" else "  (centre rib held at 1.00 mm)"))
    ax.set_ylabel(f"von Mises stress at {RPM_REF:,} rpm [MPa]")
    ax.set_ylim(0, max(r["peak"]["value_mpa"] for r in good) * 1.25)
    ax.grid(alpha=0.3)
    ax.legend(loc="lower right", fontsize=9, title="marker: circle = verdict peak in a bridge, square = in the rib", title_fontsize=8.5)
    ax.set_title(bench.wrap_title(f"Optimized rotor: stress at {RPM_REF:,} rpm against the {'centre rib' if which == 'rib' else 'outer bridge'} width, magnets loose (pressure model), "
                                  f"CalculiX, h_min {H_MAIN:g} mm"), fontsize=10, loc="left")
    fig.savefig(target, dpi=140, bbox_inches="tight", pad_inches=0.2)
    plt.close(fig)


# ============================================================================================== the printed report
def line(log: bench.Log, cells, widths) -> None:
    log.say("  " + "  ".join(str(c).ljust(w) if i == 0 else str(c).rjust(w) for i, (c, w) in enumerate(zip(cells, widths))))


def report_tables(log: bench.Log, ref: dict, opt: dict) -> None:
    """The reference against the optimized rotor, every check of the case, as the run's own words."""
    R, O = ref["data"], opt["data"]
    log.say("")
    log.say("=== the reference rotor against the optimized rotor (pressure model = the verdict, unless said; 10,000 rpm, h_min 0.1 mm, unless said)")
    w = (58, 26, 26)
    line(log, ("", "reference", "optimized"), w)

    def row(name, a, b):
        line(log, (name, a, b), w)

    g1, g2 = R["geometry"], O["geometry"]
    row("centre rib [mm] (design.screen / strip cut)", f"{g1['rib_mm']:.4f} / {R['hand']['strips_mm']['centre_rib']:.4f}", f"{g2['rib_mm']:.4f} / {O['hand']['strips_mm']['centre_rib']:.4f}")
    row("outer bridge [mm] (design.screen / strip cut)", f"{g1['bridge_mm']:.4f} / {R['hand']['strips_mm']['outer_bridge']:.4f}", f"{g2['bridge_mm']:.4f} / {O['hand']['strips_mm']['outer_bridge']:.4f}")
    row("screen: feasible, q rib [mm]", f"{g1['feasible']}, {g1['q_rib_mm']:.2f}", f"{g2['feasible']}, {g2['q_rib_mm']:.2f}")
    row("slot fillet [mm], corner clearance [mm]", f"{g1['slot_fillet_mm']:g}, {g1['corner_clearance_mm']:+.3f}", f"{g2['slot_fillet_mm']:g}, {g2['corner_clearance_mm']:+.3f}")
    row("mesh: nodes in the steel, steel triangles, min SICN", f"{R['mesh']['nodes_steel']:,}, {R['mesh']['triangles_steel']:,}, {R['mesh']['min_sicn']:.3f}",
        f"{O['mesh']['nodes_steel']:,}, {O['mesh']['triangles_steel']:,}, {O['mesh']['min_sicn']:.3f}")
    row("loads on one magnet's walls [MPa]: long wall, end wall", f"{R['wall_example']['long']['pressure_mpa']:.3f}, {R['wall_example']['end']['pressure_mpa']:.3f}",
        f"{O['wall_example']['long']['pressure_mpa']:.3f}, {O['wall_example']['end']['pressure_mpa']:.3f}")
    row("peak, pressure, 10,000 rpm [MPa]", f"{R['peak']['value_mpa']:.3f}", f"{O['peak']['value_mpa']:.3f}")
    row("peak, pressure, 12,000 rpm [MPa]; ratio to 10,000 (1.4400)", f"{R['peak_12000_mpa']:.3f}; {R['peak_12000_mpa'] / R['peak']['value_mpa']:.4f}",
        f"{O['peak_12000_mpa']:.3f}; {O['peak_12000_mpa'] / O['peak']['value_mpa']:.4f}")
    row("where: strip, slot, corner", f"{R['peak']['feature']}, {R['peak']['slot']}, {R['peak']['detail']}", f"{O['peak']['feature']}, {O['peak']['slot']}, {O['peak']['detail']}")
    row("   position (x, y) [mm], r [mm], distance to the cut [mm]", f"({R['peak']['x_mm']:.2f}, {R['peak']['y_mm']:.2f}), {R['peak']['r_mm']:.2f}, {R['peak']['distance_to_strip_mm']:.3f}",
        f"({O['peak']['x_mm']:.2f}, {O['peak']['y_mm']:.2f}), {O['peak']['r_mm']:.2f}, {O['peak']['distance_to_strip_mm']:.3f}")
    row("largest in each outer bridge's zone [MPa]", ", ".join(f"{v:.2f}" for v in R["zones"]["bridges"].values()), ", ".join(f"{v:.2f}" for v in O["zones"]["bridges"].values()))
    row("   spread of the four [%] (limit 3)", f"{R['bridge_spread_percent']:.2f}", f"{O['bridge_spread_percent']:.2f}")
    row("largest in each centre rib's zone [MPa]", ", ".join(f"{v:.2f}" for v in R["zones"]["ribs"].values()), ", ".join(f"{v:.2f}" for v in O["zones"]["ribs"].values()))
    eq = lambda d: max(a["equilibrium"]["residual_percent"] for a in (r["analysis"] for r in d["runs"].values()))
    row("equilibrium, worst of the solves [%] (limit 0.05)", f"{eq(R):.5f}", f"{eq(O):.5f}")
    a, b = R["runs"][f"bonded_{RPM_REF}_h{H_MAIN:g}"]["analysis"], O["runs"][f"bonded_{RPM_REF}_h{H_MAIN:g}"]["analysis"]
    row("bonded: peak of the whole steel [MPa], where", f"{a['peak']['value_mpa']:.2f}, {a['peak']['place']['feature']}", f"{b['peak']['value_mpa']:.2f}, {b['peak']['place']['feature']}")
    row("bonded: largest in the bridges' zones [MPa]; the pressure model's is this times", f"{a['bridge_peak_mpa']:.2f}; {R['bridge_peak_mpa'] / a['bridge_peak_mpa']:.1f}",
        f"{b['bridge_peak_mpa']:.2f}; {O['bridge_peak_mpa'] / b['bridge_peak_mpa']:.1f}")
    row("bonded: largest in the rib's zone [MPa]", ", ".join(f"{v:.2f}" for v in R["bonded"]["zones"]["ribs"].values()), ", ".join(f"{v:.2f}" for v in O["bonded"]["zones"]["ribs"].values()))
    ba, bb = R["bonded"]["magnet_bodies"], O["bonded"]["magnet_bodies"]
    row(f"bonded: largest in a magnet body [MPa]: von Mises; principal (N45 {N45_TENSILE_MPA:g})", f"{ba['von_mises_mpa']:.2f}; {ba['principal_max_mpa']:.2f}", f"{bb['von_mises_mpa']:.2f}; {bb['principal_max_mpa']:.2f}")
    row("bore, mean hoop stress [MPa], Lame's 14.10: pressure; bonded", f"{R['bore']['pressure']['hoop_mean_mpa']:.3f} ({percent(R['bore']['pressure']['hoop_mean_mpa'], R['bore']['pressure']['lame_mpa']):+.2f} %); "
        f"{R['bore']['bonded']['hoop_mean_mpa']:.3f} ({percent(R['bore']['bonded']['hoop_mean_mpa'], R['bore']['bonded']['lame_mpa']):+.2f} %)",
        f"{O['bore']['pressure']['hoop_mean_mpa']:.3f} ({percent(O['bore']['pressure']['hoop_mean_mpa'], O['bore']['pressure']['lame_mpa']):+.2f} %); "
        f"{O['bore']['bonded']['hoop_mean_mpa']:.3f} ({percent(O['bore']['bonded']['hoop_mean_mpa'], O['bore']['bonded']['lame_mpa']):+.2f} %)")
    row("mesh: peak at h_min 0.2, 0.1, 0.05 mm [MPa]", ", ".join(f"{v:.2f}" for v in R["convergence"]["peaks_mpa"]), ", ".join(f"{v:.2f}" for v in O["convergence"]["peaks_mpa"]))
    short = {"outer bridge": "bridge", "centre rib": "rib"}
    row("   where the peak sits on each mesh", ", ".join(short.get(R["convergence"]["levels"][h]["peak"]["feature"], "?") for h in ("0.2", "0.1", "0.05")),
        ", ".join(short.get(O["convergence"]["levels"][h]["peak"]["feature"], "?") for h in ("0.2", "0.1", "0.05")))
    for group, name in (("bridges", "outer bridges"), ("ribs", "centre ribs")):
        row(f"   largest in the {name}' zones at h_min 0.2 | 0.1 | 0.05 [MPa]",
            " | ".join(",".join(f"{v:.1f}" for v in zone_groups(R["runs"][f"pressure_{RPM_REF}_h{h}"]["analysis"])[group].values()) for h in ("0.2", "0.1", "0.05")),
            " | ".join(",".join(f"{v:.1f}" for v in zone_groups(O["runs"][f"pressure_{RPM_REF}_h{h}"]["analysis"])[group].values()) for h in ("0.2", "0.1", "0.05")))
    row("   main mesh against the finest [%] (limit 2)", f"{R['convergence']['main_against_finest_percent']:+.2f}", f"{O['convergence']['main_against_finest_percent']:+.2f}")
    row("   observed order; extrapolated peak [MPa]", f"{R['convergence']['observed_order']:.2f}; {R['convergence']['extrapolated_mpa']:.2f}" if R["convergence"]["usable"] else "no usable order",
        f"{O['convergence']['observed_order']:.2f}; {O['convergence']['extrapolated_mpa']:.2f}" if O["convergence"]["usable"] else "no usable order")
    row("   extrapolated for p = 2 [MPa]; main mesh against it [%]", f"{R['convergence']['extrapolated_p2_mpa']:.2f}; {R['convergence']['main_against_extrapolated_p2_percent']:+.2f}",
        f"{O['convergence']['extrapolated_p2_mpa']:.2f}; {O['convergence']['main_against_extrapolated_p2_percent']:+.2f}")
    log.say(f"  reference peak: {describe_peak(R['peak'])}")
    log.say(f"  optimized peak: {describe_peak(O['peak'])}")
    log.say("")
    log.say("=== the three strips of a pole (the mean of the two poles), pressure model, 10,000 rpm: the force through each cut, from the stress")
    w2 = (28, 14, 14, 12, 14, 14, 12, 14)
    line(log, ("", "ref width", "ref F [N]", "ref share", "ref nominal", "opt width", "opt F [N]", "opt nominal"), w2)
    for key, name in STRIP_NAMES:
        a, b = R["strips"]["strips"][key], O["strips"]["strips"][key]
        line(log, (name, f"{a['width_mm']:.3f} mm", f"{a['force_n']:.1f}", f"{a['share_percent']:.1f} %", f"{a['nominal_mpa']:.2f} MPa", f"{b['width_mm']:.3f} mm", f"{b['force_n']:.1f}",
                   f"{b['nominal_mpa']:.2f} MPa"), w2)
    ra, ob = R["strips"], O["strips"]
    log.say(f"  opt share of the three strips (rib, bridge, bridge): {', '.join(f'{v['share_percent']:.1f} %' for v in ob['strips'].values())}; "
            f"ref: {', '.join(f'{v['share_percent']:.1f} %' for v in ra['strips'].values())}")
    row("total carried (poles 0, 1) [N]", " ".join(f"{v:.1f}" for v in ra["carried_by_pole_n"]), " ".join(f"{v:.1f}" for v in ob["carried_by_pole_n"]))
    row("hand estimate, pull along the axis (poles 0, 1) [N]; the note's sum [N]", f"{' '.join(f'{v:.1f}' for v in ra['hand_pull_by_pole_n'])}; {ra['hand_pull_note_n']:.1f}",
        f"{' '.join(f'{v:.1f}' for v in ob['hand_pull_by_pole_n'])}; {ob['hand_pull_note_n']:.1f}")
    row("strips against the hand estimate [%] (limit 10); against the refined pull [%] (limit 2)", f"{ra['against_hand_percent'][0]:+.2f} {ra['against_hand_percent'][1]:+.2f}; "
        f"{ra['against_refined_percent'][0]:+.2f} {ra['against_refined_percent'][1]:+.2f}", f"{ob['against_hand_percent'][0]:+.2f} {ob['against_hand_percent'][1]:+.2f}; "
        f"{ob['against_refined_percent'][0]:+.2f} {ob['against_refined_percent'][1]:+.2f}")
    row("cap area [mm2]; the magnets' share of the pull [%]", f"{ra['cap_area_mm2']:.1f}; {100 * ra['magnet_share_of_pull']:.1f}", f"{ob['cap_area_mm2']:.1f}; {100 * ob['magnet_share_of_pull']:.1f}")
    row("nominal stress over the three widths [MPa]: finite elements; hand (axis); hand (note)", f"{ra['nominal_mean_mpa']:.2f}; {ra['hand_nominal_mpa']:.2f}; {ra['hand_nominal_note_mpa']:.2f}",
        f"{ob['nominal_mean_mpa']:.2f}; {ob['hand_nominal_mpa']:.2f}; {ob['hand_nominal_note_mpa']:.2f}")
    row("peak / nominal: against the finite elements'; the hand's; the peak's own strip", f"{ra['peak_over_nominal_mean']:.2f}; {ra['peak_over_nominal_hand']:.2f}; "
        f"{ra['peak_over_nominal_local']:.2f} ({ra['peak_strip']})", f"{ob['peak_over_nominal_mean']:.2f}; {ob['peak_over_nominal_hand']:.2f}; {ob['peak_over_nominal_local']:.2f} ({ob['peak_strip']})")
    log.say("")
    log.say("=== safe speeds: n = n_ref x sqrt(limit / peak), n_ref = 10,000 rpm, the verdict peak (pressure model, main mesh); the limits are the placeholders of 04 Structural 6.4")
    w3 = (42, 8, 11, 12, 11, 11, 11, 12, 11)
    line(log, ("limit", "MPa", "ref [rpm]", "ref max op.", "ref v [m/s]", "opt [rpm]", "opt max op.", "opt v [m/s]", ""), w3)
    for a, b in zip(R["speeds"], O["speeds"]):
        line(log, (a["name"], f"{a['limit_mpa']:.0f}", f"{a['speed_rpm']:,.0f}", f"{a['operating_rpm']:,.0f}" if a["operating_rpm"] else "-", f"{a['tip_speed_m_s']:.1f}", f"{b['speed_rpm']:,.0f}",
                   f"{b['operating_rpm']:,.0f}" if b["operating_rpm"] else "-", f"{b['tip_speed_m_s']:.1f}", ""), w3)
    r_ro = design.REFERENCE.r_ro
    comfort = rpm_at_tip_speed(TIP_SPEED_COMFORT_M_S, r_ro)
    log.say(f"  context: a tip speed of {TIP_SPEED_COMFORT_M_S:g} m/s (rotor radius {r_ro:g} mm) is {comfort:,.0f} rpm; ROSS's first critical speed is {ROSS_FIRST_CRITICAL_RPM:,} rpm (a 60 mm stack)")
    for name, d in (("reference", R), ("optimized", O)):
        at = {f"{TIP_SPEED_COMFORT_M_S:g} m/s ({comfort:,.0f} rpm)": comfort, "12,000 rpm": 12_000, "15,000 rpm": 15_000, f"ROSS critical ({ROSS_FIRST_CRITICAL_RPM:,} rpm)": ROSS_FIRST_CRITICAL_RPM}
        log.say(f"  {name}: the verdict peak at " + "; ".join(f"{k}: {d['peak']['value_mpa'] * (n / RPM_REF) ** 2:.1f} MPa" for k, n in at.items()))
    for name, d in (("reference", R), ("optimized", O)):
        conv = d["convergence"]
        for what, value in (("the finest mesh's peak", conv["peaks_mpa"][2]), ("the extrapolated peak", conv["extrapolated_mpa"])):
            if value is None:
                log.say(f"  {name}: no extrapolated peak (the order of convergence cannot be read from the three meshes)")
                continue
            ex = safe_speeds(value, r_ro)
            log.say(f"  {name}: with {what} ({value:.2f} MPa) the speeds fall by {100 * (1 - ex[0]['speed_rpm'] / d['speeds'][0]['speed_rpm']):.2f} %: yield {ex[0]['speed_rpm']:,.0f} rpm; "
                    f"yield/2 at 1.2 x: max operating {ex[2]['operating_rpm']:,.0f} rpm; Goodman {ex[4]['speed_rpm']:,.0f} rpm")


def report_sweeps(log: bench.Log, s: dict, opt_data: dict) -> None:
    def table(title: str, rows: list[dict]) -> None:
        log.say("")
        log.say(f"=== {title}")
        w = (36, 9, 9, 7, 9, 9, 10, 12, 11, 10, 10, 9, 9, 60)
        line(log, ("run", "rib [mm]", "br. [mm]", "q rib", "feasible", "peak MPa", "rib zone", "bridge zone", "carried N", "rib share", "rib nom.", "y/2 rpm", "Goodm.rpm",
                   "where the peak sits"), w)
        for r in rows:
            if "peak" not in r:
                line(log, (r["tag"], "", "", "", "", "error", "", "", "", "", "", "", "", r.get("error", "")), w)
                continue
            g, pk = r["geometry"], r["peak"]
            rib_nominal = r["rib_force_n"] / (design.REFERENCE.stack * g["rib_mm"])
            speeds = safe_speeds(pk["value_mpa"], design.REFERENCE.r_ro)
            line(log, (r["tag"], f"{g['rib_mm']:.3f}", f"{g['bridge_mm']:.3f}", f"{g['q_rib_mm']:.2f}", g["feasible"], f"{pk['value_mpa']:.2f}", f"{r['rib_zone_mpa']:.2f}", f"{r['bridge_zone_mpa']:.2f}",
                       f"{r['carried_n']:.1f}", f"{r['rib_share_percent']:.1f} %", f"{rib_nominal:.2f}", f"{speeds[2]['operating_rpm']:,.0f}", f"{speeds[4]['speed_rpm']:,.0f}",
                       f"{pk['feature']}, {pk['detail']}"), w)
            if g["reason"]:
                log.say(f"      screen says: {g['reason']}")

    log.say("")
    log.say("the sweeps' columns: rib zone and bridge zone are the largest von Mises stress [MPa] within 2 mm of the rib's cut and of the four bridges' cuts; carried N is the force the three strips of a pole "
            "carry (mean of the two poles); rib nom. is the rib's force / (its width x the stack) [MPa]; y/2 rpm is the maximum operating speed for yield / 2 at a 1.2 x overspeed, Goodm.rpm the speed at "
            "which the peak reaches Goodman's 280 MPa (both from this row's peak, scaled by the speed squared)")
    table(f"the rib: mag_y moved, mag_x re-solved for a {BRIDGE_GUARD_MM:g} mm bridge (pressure model, 10,000 rpm, h_min {H_MAIN:g} mm)", s["rib"])
    good = [r for r in s["rib"] if "peak" in r]
    for a, b in zip(good[:-1], good[1:]):                       # where the rib's own peak falls below the bridges': a linear interpolation between the two solves that bracket it
        da, db = a["rib_zone_mpa"] - a["bridge_zone_mpa"], b["rib_zone_mpa"] - b["bridge_zone_mpa"]
        if da > 0 >= db:
            x = a["geometry"]["rib_mm"] + (b["geometry"]["rib_mm"] - a["geometry"]["rib_mm"]) * da / (da - db)
            log.say(f"  the rib's zone stops being the highest at a rib of about {x:.2f} mm (linear interpolation between the solves at {a['geometry']['rib_mm']:.3f} and "
                    f"{b['geometry']['rib_mm']:.3f} mm): below it the rib sets the verdict peak, above it the bridges do")
            break
    table(f"the bridge: mag_x re-solved, mag_y as in the design, so the rib stays {opt_data['geometry']['rib_mm']:.3f} mm (pressure model, 10,000 rpm, h_min {H_MAIN:g} mm)", s["bridge"])
    f = s["fillet"]
    rows = [f["design"], f["a"], f["b"], f["b_fine"]]
    table(f"the slot's corner radius: 0.3 mm (the design), 0.5 mm with the slot unmoved (a), 0.5 mm with the rib and bridge put back (b), and (b) on the finest mesh", rows)
    for r in rows:
        if "geometry" in r and "peak" in r:
            g = r["geometry"]
            log.say(f"  {r['tag']}: fillet {g['slot_fillet_mm']:g} mm, corner clearance {g['corner_clearance_mm']:+.3f} mm, walls on a magnet: long {r['walls']['long']['length_mm']:.2f} mm at "
                    f"{r['walls']['long']['pressure_mpa']:.3f} MPa, end {r['walls']['end']['length_mm']:.2f} mm at {r['walls']['end']['pressure_mpa']:.3f} MPa; h_min {r['h_min']:g} mm; "
                    f"mag_x {g['mag_x']:.4f}, mag_y {g['mag_y']:.4f}")
    if all("peak" in f[k] for k in ("a", "b", "b_fine")):
        main_design, main_b, main_a = f["design"]["peak"]["value_mpa"], f["b"]["peak"]["value_mpa"], f["a"]["peak"]["value_mpa"]
        fine_b = f["b_fine"]["peak"]["value_mpa"]
        log.say(f"  the peak with a {FILLET_NEW_MM:g} mm fillet against the design's 0.3 mm, at the same rib and bridge (b): {percent(main_b, main_design):+.1f} % on the main mesh "
                f"({main_b:.2f} against {main_design:.2f} MPa), {percent(fine_b, f['design_fine_peak_mpa']):+.1f} % on the finest mesh, h_min {H_FINE:g} mm ({fine_b:.2f} against "
                f"{f['design_fine_peak_mpa']:.2f} MPa); with the slot unmoved (a): {percent(main_a, main_design):+.1f} % on the main mesh ({main_a:.2f} MPa)")
        log.say(f"  the largest stress in the outer bridges' zones: 0.3 mm {f['design']['bridge_zone_mpa']:.2f} MPa; {FILLET_NEW_MM:g} mm (b) {f['b']['bridge_zone_mpa']:.2f} MPa "
                f"({percent(f['b']['bridge_zone_mpa'], f['design']['bridge_zone_mpa']):+.1f} %); in the rib's zone: {f['design']['rib_zone_mpa']:.2f} and {f['b']['rib_zone_mpa']:.2f} MPa "
                f"({percent(f['b']['rib_zone_mpa'], f['design']['rib_zone_mpa']):+.1f} %)")


# ============================================================================================== the case
def run_studies(p_ref: mg.Params, p_opt: mg.Params, base: dict, solves: Path, root: Path, log: bench.Log, gated: dict) -> tuple[dict, dict, dict]:
    """Both rotors with the same checks (the reference first, from the gate's solve), then the sensitivities on the optimized design. Returns (reference, optimized, sweeps)."""
    ref = study_rotor("reference", p_ref, solves / "reference", root, log, first=gated)
    opt = study_rotor("optimized", p_opt, solves / "optimized", root, log)
    return ref, opt, sensitivities(opt, base, solves / "sweeps", root, log)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="stress_ccx", description="The optimized rotor under centrifugal load in CalculiX, its safe speeds and its sensitivities.")
    parser.add_argument("--gate-only", action="store_true", help="only the regression gate: B4 on the reference rotor")
    parser.add_argument("--redo", action="store_true", help="no solve: the printed report and the pictures again from the state.pkl that the last full run saved (its log and results go to the redo folder)")
    args = parser.parse_args(argv)
    import femm_session
    femm_session.lower_priority()                   # R10: this script's own Python (Gmsh's library, numpy, matplotlib) runs below normal too; never raised
    for stream in (sys.stdout, sys.stderr):
        reconfigure = getattr(stream, "reconfigure", None)
        if reconfigure is not None:
            reconfigure(errors="replace")
    started = time.perf_counter()
    root = out_folder()
    work = root / "gate" if args.gate_only else root          # a gate-only run leaves the full run's files alone
    work.mkdir(parents=True, exist_ok=True)
    log = bench.Log(work / "redo" if args.redo else work)
    state_path = root / "state.pkl"
    results: dict = {"case": CASE, "time": time.strftime("%Y-%m-%d %H:%M:%S"), "redo": bool(args.redo)}
    checks: list[dict] = []
    code = 0

    def save_results() -> None:
        results["checks"] = checks
        results["seconds_total"] = round(time.perf_counter() - started, 1)
        results["exit_code"] = code
        target = (work / "redo" if args.redo else work) / "results.json"        # a --redo leaves the full run's results.json as it was
        target.parent.mkdir(parents=True, exist_ok=True)
        bench.write_json(target, results)

    try:
        base = load_final()["overrides"]
        p_ref, p_opt = design.REFERENCE, design.params(base)
        log.say(f"stress_ccx: the reference rotor and the optimized rotor ({', '.join(f'{k} {v:.4f}' for k, v in base.items())}), CalculiX plane stress, loads at {RPM_REF:,} rpm")
        log.say(f"  output: {work}")
        if args.redo:
            state = pickle.loads(state_path.read_bytes())
            ref, opt, sweep, checks, results["gate"] = state["ref"], state["opt"], state["sweep"], state["checks"], state["gate"]
            results["seconds_studies"] = state["seconds_studies"]
            log.say("--redo: the report and the pictures from state.pkl; no solve is made")
        else:
            solves = fresh_folder(work / "solves", root)
            holder = femm_queue.holder()
            log.say("  the solver slot is free" if holder is None else f"  the solver slot is held by '{holder['label']}' (pid {holder['pid']}): the first ccx run waits for it")
            # ---- 1. the regression gate
            log.say("--- the regression gate: B4's method on the reference rotor (pressure model, 10,000 rpm, h_min 0.1 mm)")
            gated = gate(p_ref, solves / "reference", root, log)
            item = gated["check"]
            checks.append(item)
            results["gate"] = {"check": item, "peak": peak_record(gated["analysis"]) if gated["analysis"]["finished"] else None}
            log.say(f"  {item['verdict']}: {item['name']}: " + (f"{item['measured']:.4f} MPa against {GATE_MPA} MPa ({item['difference_percent']:+.3f} %, limit {GATE_TOLERANCE_PERCENT} %)"
                                                               if "difference_percent" in item else item["measured"]))
            if item["verdict"] != "PASS":
                log.say("the regression gate FAILED: nothing else is run (the case no longer reproduces B4 on the reference rotor)")
                code = 1
                return code
            if args.gate_only:
                log.say("--gate-only: stopped after the gate")
                return code
            # ---- 2 and 3. both rotors, the same checks; the sensitivities on the optimized design
            ref, opt, sweep = run_studies(p_ref, p_opt, base, solves, root, log, gated)
            checks += ref["checks"] + opt["checks"]
            results["seconds_studies"] = round(time.perf_counter() - started, 1)
            state_path.write_bytes(pickle.dumps({"ref": ref, "opt": opt, "sweep": sweep, "checks": checks, "gate": results["gate"], "seconds_studies": results["seconds_studies"]}, protocol=4))
        names = {"optimized": "stress_optimized_10000rpm.png", "compare": "stress_reference_vs_optimized_10000rpm.png", "rib": "peak_vs_rib.png", "bridge": "peak_vs_bridge.png"}
        results.update({"inputs": {"overrides": base, "reference": bench.params(p_ref), "optimized": bench.params(p_opt)}, "reference": ref["data"], "optimized": opt["data"], "sweeps": sweep,
                        "pictures": {k: str(root / v) for k, v in names.items()},
                        "limits": {"yield_mpa": YIELD_MPA, "fatigue_mpa": FATIGUE_MPA, "goodman_mpa": GOODMAN_MPA, "s_ultimate_mpa": S_ULTIMATE_MPA, "safety_factors": SAFETY_FACTORS,
                                   "overspeed": OVERSPEED, "n45_tensile_mpa": N45_TENSILE_MPA}})
        save_results()                                       # the numbers are kept before the report and the pictures are made
        # ---- 4. what was found, and the pictures
        report_tables(log, ref, opt)
        report_sweeps(log, sweep, opt["data"])
        fig_optimized(opt["model"], opt["stress"], opt["data"], root / names["optimized"])
        fig_compare(ref, opt, root / names["compare"])
        fig_sweep(sweep["rib"], "rib", root / names["rib"])
        fig_sweep(sweep["bridge"], "bridge", root / names["bridge"])
        log.say("")
        log.say("=== pictures: " + ", ".join(str(root / n) for n in names.values()))
        log.say("")
        log.say("=== the checks")
        for item in checks:
            expected, measured = bench.describe(item)
            word = item["verdict"] if bench.is_gating(item) else f"REPORTED ({item['verdict']})"
            log.say(f"  {word}  {item['name']}")
            log.say(f"        expected: {expected}")
            log.say(f"        measured: {measured}")
        failed = [c for c in checks if bench.is_gating(c) and c["verdict"] != "PASS"]
        reported_failed = [c for c in checks if not bench.is_gating(c) and c["verdict"] == "FAIL"]
        log.say("")
        log.say(f"stress_ccx: {'FAIL' if failed else 'PASS'}: {len(failed)} of {sum(bench.is_gating(c) for c in checks)} gating checks fail; {len(reported_failed)} of "
                f"{sum(not bench.is_gating(c) for c in checks)} reported checks fail")
        for c in reported_failed:
            log.say(f"  reported, did not pass: {c['name']}")
        code = 1 if failed else 0
        return code
    except KeyboardInterrupt:
        log.say("stress_ccx: interrupted")
        code = 130
        return code
    except Exception:                                         # keep the traceback in run.log, say so in results.json, and let it end the run with its own exit status
        code = 1
        log.say("stress_ccx: stopped by an error:\n" + traceback.format_exc())
        raise
    finally:
        log.say(f"stress_ccx: {time.perf_counter() - started:.0f} s")
        try:
            log.add_run(bench.process_record())
            log.save()
            save_results()
        except Exception as error:                            # the run's own words are already on the screen
            print(f"stress_ccx: the results could not be saved: {error}", file=sys.stderr)


if __name__ == "__main__":
    sys.exit(main())
