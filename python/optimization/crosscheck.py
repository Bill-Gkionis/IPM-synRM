"""optimization.crosscheck: the joint problem searched a second way: a model of the torque waveform tells the search where a FEMM run is worth its seconds.

What it is, and why
    The joint problem has six variables (bridge_mm, mag_y, v_half_deg, slot_open, slot_half_deg, r_wall; h_open held at 0.5 mm), two objectives to minimise (-mean
    torque and the ripple in per cent of it), a torque floor of 0.8 N.m, and at most 2,000 designs scored in FEMM (design.Evaluator, scoring "P": every design is
    given the reference's copper loss; coarse mesh). This file is the optimizer of a cross-check: the same black box and the same budget handed to a second method
    that was chosen on its own, so that two fronts which agree cannot both be accidents of one method. The method was chosen before any code was written, from what
    the sweeps and the screening had already measured (out\\optimization\\sweeps_rotor, sweeps_stator, screening; nothing was re-measured):
      smooth torque   the mean torque has a noise floor of 0.56 mN.m (0.05 %: a design a hair away is meshed again) against a spread of 0.5 to 2 N.m across the
                      box, and it moves almost monotonically: up as the outer bridge thins to its guard, as mag_y falls to the centre rib's limit, as v_half_deg grows,
                      as the slot grows. A model can learn it from a hundred designs.
      rugged ripple   the ripple has valleys, not a slope: along mag_y and along v_half_deg it dips V-shaped, 1.8 mm and about 39 degrees apart, and on the (mag_y,
                      v_half_deg) map the valleys run as parallel diagonal bands; Morris puts sigma / mu* of the ripple at 0.76 to 1.10 on all six variables (they act
                      together and nonlinearly), with mu* of 4.8 to 7.0 points for bridge_mm, mag_y, v_half_deg and slot_open and 0.8 and 0.5 for slot_half_deg and r_wall.
                      The ripple's noise floor is 0.14 points, and a central difference needs a step of 0.07 mm in mag_y and 0.8 to 1.6 degrees in v_half_deg to be
                      usable: as wide as the valleys. A gradient is no use.
      explicit limits three of the 10 constraints of design.screen are ever active inside the box (the centre rib, the tooth, the yoke), each is a formula, and together they
                      take a third of the box away (672 of 1,024 random designs are feasible); the best torque sits on them.
      cost and room   a design costs about 1.3 s of pool time (the first block: 148 designs in 192 s) and 63 ms to screen in Python; 2,000 designs is 3.5 a dimension, and
                      every design returns the whole 15-sample torque waveform, not only its two scores.
    So: a problem where a model can be learned, where gradients cannot be used, where the constraints are known exactly, and where a few hundred designs can go
    into learning. That is a surrogate-assisted search, and this one has four parts.
      1. A feasible cube. The three limits are solved for in closed form, so that every point of the unit cube [0,1]^6 is a design the screen accepts, and the
         faces of the cube are the limits themselves (the best designs sit on them).
      2. A surrogate of the waveform, not of the ripple. Five Gaussian processes predict the 15 numbers of one ripple period of the torque (the mean and the Fourier
         coefficients of the first seven harmonics), and the ripple is computed from the predicted waveform by the black box's own peak-to-peak rule.
      3. An acquisition, the expected improvement of the front. Each round the surrogate says, for about 4,000 candidate designs, how much hypervolume scoring each one
         would add to the front, allowing for what it does not know; the best 21 (kept apart) and the 3 most uncertain go to FEMM.
      4. Two looks that do not use the pool search. A model-free probe (96 designs the surrogate did not choose) that must find nothing before a run may stop, and
         the surrogate's own optima found by another optimizer (differential evolution over the whole cube) and a map of the two variables that make the valleys.
    Rejected, and why. A genetic algorithm (NSGA-II and its relatives) needs no model and would work with a budget this size, but it learns nothing from the smoothness
    and has to cover the valleys by chance; and a cross-check is worth more when the two methods would fail differently. A weighted sum of the objectives (or a decomposition into weights) cannot reach a front
    made of valleys: the dents of the front are invisible to it. Gradient and quasi-Newton methods: the noise floor and the valleys (above). A Gaussian process on the ripple
    itself, as in most Bayesian optimizers: tested on the sweeps' own data it is worse than the waveform model (leave-one-out ripple error, same kernel: 0.73 against
    0.29 points on the 216 designs of the rotor sweeps, 1.26 against 1.07 on the 84 of the Morris design). A grid or a Latin hypercube alone: no use of the smoothness.

The method, and the maths
    The cube. The tooth is the chord between two shoe corners, 2 r_shoe sin(5 deg - a), a = slot_half_deg, so tooth >= 3.0 mm is a <= a_max = 5 deg - asin(1.5 / r_shoe)
    = 3.16332 deg. The yoke is r_so - r_wall (1 + sin a) / cos a (the round slot bottom reaches r_wall / cos a + r_wall tan a), so yoke >= 3.5 mm is r_wall <= (r_so - 3.5) cos a
    / (1 + sin a) (57.25 mm at a_max). The centre rib is 2 (mag_y - h(v)), h depending on the magnet's angle v only (a rotated rounded rectangle), so rib >= 1.0 mm is
    mag_y >= ymin(v); ymin is read from the screen's own rib number at one mag_y for each v on a 0.25 degree grid (rib_table.json) and interpolated (4.04 mm at 40 degrees,
    4.85 at 60, 5.15 at 89). Checked against design.screen: the first two to 1e-14, the rib wherever the two slots do not overlap, and 144 designs of the cube (corners
    included) were all accepted. The map is then
        bridge_mm = lo + u1 (hi - lo)        mag_y = ymin(v) + u2 (7 - ymin(v))       v_half_deg = lo + u3 (hi - lo)       slot_open = lo + u4 (hi - lo)
        slot_half_deg = lo + u5 (min(hi, a_max) - lo)        r_wall = lo + u6 (min(hi, (r_so - 3.5) cos a / (1 + sin a)) - lo)
    with margins of 1e-6 mm (bridge, tooth, yoke) and 2e-4 mm (rib), so that a face is inside the limit, never past it. u = 0 and u = 1 are the limits' edges.
    The waveform. The black box samples one ripple period (60 electrical degrees) at 15 rotor positions T_n. Their discrete Fourier series c_k = (1/15) sum_n T_n
    exp(-2 pi i k n / 15), k = 0..7, has c_0 = the mean torque and c_k the 6k-th harmonic of the electrical frequency (a complex number): 15 real numbers in all
    (c_0, Re c_k, Im c_k). The ripple is the peak-to-peak of the series rebuilt between the samples, as design.fourier_pkpk does it (200 points a sample); the surrogate
    does it on 8 points a sample (0.03 points of difference at most), with one inverse FFT. The harmonics' sizes (rms over the designs of the Morris and rotor-sweep rows): the
    6th 0.02 N.m, the 12th 0.003, the 18th (the slot harmonic) 0.03, the 24th to 42nd 0.002 or less.
    The Gaussian process. The prediction of a function f from n designs, if f is drawn from a Gaussian process with the squared-exponential kernel k(x, x') = exp(-1/2 sum_d
    ((x_d - x'_d) / l_d)^2) (one length scale l_d per variable: a variable that does not matter gets a long one) and each measurement is f plus noise of variance tau^2, is
    the mean k(x, X) K^-1 y and the variance k(x, x) - k(x, X) K^-1 k(X, x) + tau^2, with K = k(X, X) + tau^2 I. The l_d and tau^2 maximise the marginal likelihood
    -1/2 y' K^-1 y - 1/2 log |K| - n/2 log 2 pi; its gradient is analytic (1/2 tr[(a a' - K^-1) dK/dt], a = K^-1 y, dK/d log l_d = K * (dx_d / l_d)^2; tested against finite
    differences) and L-BFGS-B climbs it on at most 600 designs. Five groups of outputs, each its own kernel, standardised to zero mean and unit variance: the mean torque (all
    scored designs), and, on the designs with T >= 0.55 N.m, the 6th harmonic (Re, Im), the 12th, the 18th, and the 24th to 42nd together. Separate kernels matter: one kernel
    for all the harmonics gave 0.89 points of leave-one-out ripple error on the rotor sweeps, one per group 0.44, and the squared-exponential kernel in place of Matern 5/2
    0.29. The data prefer the smoother kernel: the black box is the solution of a PDE, smooth up to mesh noise.
    The acquisition. For a candidate x the surrogate gives a mean and a standard deviation for each of the 15 numbers. Draw 48 scenarios (the same 48 standard normal vectors
    for every candidate): mean + sd eps, rebuilt into (-T, ripple). The hypervolume a scenario's point p adds to the current front is the area it dominates that the front
    does not. The front is a staircase of m points that cuts the torque axis into m + 1 strips, in each of which it dominates every ripple above a fixed level, so
        HVI(p) = sum over strips of (the part of the strip right of p_1) x max(0, level - p_2),       level_0 = the reference point's 40 %
    (checked against the difference of two hypervolumes). The mean over the 48 scenarios is the expected hypervolume improvement (EHVI): high where the candidate is likely
    to beat the front by a lot, which includes "the model does not know". The batch is greedy: take the candidate with the most EHVI, add its predicted outcome to the front
    (the next must beat what the first promises), rule out every candidate closer than 0.05 (box units) to a chosen one, repeat 21 times; then 3 picks of the candidates
    whose predicted ripple varies most across the scenarios among those expected to clear 0.85 N.m. A design beyond the reference point (T < 0.8 N.m or ripple > 40 %) adds
    nothing, so the torque floor needs no rule of its own.
    The candidates of a round. 2,048 fresh scrambled Sobol' points of the cube (look everywhere); around each of up to 60 members of the front, 24 perturbations (Gaussian
    steps of 0.01, 0.03 and 0.08 of the cube, each coordinate moved with probability 1/2, clipped to the cube, so that a member on a face keeps half its cloud on the face);
    and each member with one coordinate set to 0 or to 1. About 4,000 in all, none a design already scored. Optionally (--de-every N) the surrogate's own optima, below.
    The first block. 128 Sobol' points of the cube and 20 seeds (rotor values read from the sweeps' fronts, with the bridge on its guard, each with the reference stator and
    with a wide-slot, narrow-opening stator). The model is then fitted and judged by leave-one-out on the designs at T >= 0.8; if its R^2 for the ripple is under 0.70, a block
    of 64 more Sobol' points is scored, twice at most (it was not needed).
    The stopping rule. After 12 rounds, if the hypervolume gained less than 0.1 % over the last 8, a probe is scored: 96 designs the surrogate did not choose (48 Sobol' points
    where the torque model says T >= 0.85, 48 perturbations of random front members at 0.005 to 0.1 of the cube). If it also gains less than 0.1 %, the run has converged;
    otherwise it goes on. 0.1 % is the size at which mesh noise alone, picking the luckiest of many near-equal designs, lets the front creep (0.14 points of ripple over
    the 0.6 N.m the front spans is 0.2 % of the area).
    The surrogate's own optima (probe). model_optima(): for each of the torque levels 0.9, 1.1, ..., 1.9 N.m, differential evolution (240 designs, 250 generations, two
    seeds) over the whole cube on the model's mean prediction, ripple plus 200 points per N.m of missing torque; the best three designs of each run are scored. It asks what
    the pool search cannot: whether a better design lay between its candidates. And map_points(): v_half_deg over its range against mag_y from the rib limit to 7 mm
    (20 x 14 designs) with the other four variables at the front's high-torque corner, a picture of the valleys that needs no model.

Run (from the Python to FEMM folder; the solver slot is held while a pool is open, so another run waits)
    .venv/Scripts/python.exe -u -m optimization.crosscheck run [--budget 2000] [--batch 24] [--no-verify] [--max-minutes N] [--de-every N]
        the rounds, until the budget or the convergence test; resumable (evaluations.csv is its memory: it is read, nothing scored is scored again). A file STOP in the
        output folder ends it after the round in hand
    .venv/Scripts/python.exe -u -m optimization.crosscheck probe      the surrogate's optima and the map, in one pool session
    .venv/Scripts/python.exe -m optimization.crosscheck report        front.csv and knee.json from evaluations.csv (no FEMM; run and probe end with it)
    .venv/Scripts/python.exe -m optimization.crosscheck_figure        figure.png
    The run reported below was: run (converged at 723 designs), probe (it found 0.22 % more hypervolume in the low-torque tail), run --de-every 3 (the surrogate's optima in
    the pool every 3rd round; converged at 1,447 designs).
Writes out\\optimization\\crosscheck\\: evaluations.csv (every design, the schema of design.FIELDS; the stage column says design, infill, verify, probe-model or probe-map; the
n column counts from 1 in every session, so the order of the file, not n, is the order of scoring), progress.csv (a line a round), rib_table.json, front.csv, knee.json,
figure.png, and the logs.

Numbers to expect (this PC, 2026-11-09, FEMM 5.0, coarse mesh, scoring "P"; three pool sessions on the slot)
    size            1,447 designs reached FEMM (72 % of the budget) and every one came back: 148 first block, 792 infill (33 rounds of 24), 187 verify (two probes), 40 model
                    optima, 280 map; 0 refused by the screen (the cube), 0 FEMM failures; 1,425 of them at T >= 0.8 N.m; torque 0.48 to 1.97 N.m, ripple 3.59 to 28.2 %
    time            2,571 s of pool time (43 minutes) and 415 s of waiting for the slot, in three waits (259, 155 and 1 s); a round of 24 designs 32 to 40 s in all, a round
                    with the differential-evolution step 108 to 145 s
    first block     148 designs in 192 s; 126 admissible, best torque 1.926 N.m, least ripple 5.58 %, a front of 7, hypervolume 37.28; the model fitted on them: leave-one-out
                    error 0.0042 N.m on the torque and 0.70 points on the ripple (sd 4.98), R^2 0.980; length scales (box units) for the 18th harmonic 1.2, 0.27, 0.35, 1.4, 2.5, 4.5
                    (bridge, mag_y, v_half_deg, slot_open, slot_half_deg, r_wall): mag_y and v_half_deg make the valleys
    rounds          hypervolume 40.97 after 276 designs, 41.17 after 400, 41.21 after 500, 41.2161 after 723 (converged: 20 rounds, then 95 probe designs that gained 0.002 %);
                    the surrogate's errors on the 24 designs it had just chosen, over rounds 5 to 36: torque 0.0022 N.m (median; 0.0011 to 0.0077), ripple 0.30 points (0.15 to
                    1.27), mean bias +0.03 points (observed minus predicted, the model slightly hopeful); leave-one-out at round 36: 0.0023 N.m and 0.36 points, R^2 0.996
    probes          the model optima, 40 designs: new front members at 1.10 to 1.70 N.m with 3.64 to 4.04 % ripple (a bridge up to 1.71 mm, slot_half_deg down to 2.0 deg, r_wall
                    near 54 mm): hypervolume +0.22 %; the model's best ripple at each torque level was within 0.1 to 0.3 points of FEMM's. The 280-design map: +0.003 %; its
                    lowest ripple is 4.41 % at v_half_deg 55.5, mag_y on the rib limit (T 1.80 N.m); a second line of valleys at 9 to 10 %, ending on v_half_deg 89 at mag_y 5.58
    continuation    13 more rounds with the optima in the pool: hypervolume 41.3831 (the last 8 rounds +0.04 %), then 92 probe designs that gained 0.0005 %: 41.3833; the
                    front's least ripple went from 3.64 to 3.59 %
    the front       116 designs, 1.359 to 1.973 N.m, ripple 3.59 to 13.97 %, hypervolume 41.3833 N.m % to (0.8 N.m, 40 %): a plateau of 3.6 to 4.4 % from 1.36 to 1.80 N.m, then a
                    steep rise to 9 % at 1.82 N.m and a slope to 14 % at 1.97. Of the 116 designs 108 have the bridge on its guard (1.40 mm), 102 the opening at its narrowest (1.5 mm),
                    98 the slot at the tooth limit, 67 the magnets at 89 degrees (the end of the box), 52 mag_y on the rib limit. r_wall is the variable that stays inside its room:
                    56.4 to 56.8 mm for the 98 designs above 1.8 N.m, 0.4 to 0.9 mm short of the yoke limit (57.25 mm); the tail below 1.8 N.m goes to 54.0 mm, with the slot
                    down to 2.5 degrees and the bridge up to 1.55 mm. Found by: 101 infill, 10 verify, 4 map, 1 model-optimum design
    the knee        design.knee on the front's own range: 1.802 N.m, 4.41 % ripple, MTPA angle 22.3 degrees, power factor 0.868 (bridge 1.400, mag_y 4.707, v_half_deg 55.47,
                    slot_open 1.50, slot_half_deg 3.163, r_wall 56.63 mm): the corner of the plateau and the rise, on the valley that ends at the rib limit

Self-check (the answers follow).
  1. The cube puts the best designs exactly on the constraints' edges. Why is that better than proposing any point of the box and refusing the infeasible ones?
  2. The surrogate predicts the 15 numbers of the waveform, not the ripple. Why, and where does the ripple's V-shaped valley come from if the model never sees one?
  3. The hypervolume has gained 0.04 % in the last 8 rounds. Is the run finished?
Answers.
  1. A third of the box is infeasible, so refusing wastes a third of the proposals, and the designs that matter most lie on an edge, where a point is refused or accepted by
     round-off. In the cube a face is the edge with a margin on the safe side: nothing is refused, and the optimizer can sit on a limit.
  2. The ripple is a peak-to-peak of a sum of harmonics. Where one harmonic's amplitude passes through zero the peak-to-peak has a corner, which a smooth model smears; the
     amplitudes themselves are smooth functions of the geometry and a Gaussian process learns them well (0.29 against 0.73 points of leave-one-out error on the rotor sweeps).
     The V comes from the max-minus-min rule applied to the predicted waveform.
  3. Not yet, and the first run is the proof: it passed this test at 723 designs, and the surrogate's own optima then found 0.22 % more in the low-torque tail of the front,
     where the expected improvement of every candidate was tiny (probably because the area a low-torque design adds is small, and a few thousand random points rarely land on a
     valley floor). A model can be wrong where it has no designs and a pool can miss what lies between its candidates, so the stop needs a look that uses neither: the probe,
     the differential evolution, the map. And 0.04 % is below what mesh noise lets the front creep by luck, so one more look would not move it.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import differential_evolution, minimize
from scipy.stats import qmc

import femm_queue
import paths
from optimization import design

# ---------------------------------------------------------------------------------------------- the problem (fixed by the brief)
RUN = "crosscheck"
SCORING = "P"
HELD = {"h_open": 0.5}                                   # held for every candidate
VARS = ("bridge_mm", "mag_y", "v_half_deg", "slot_open", "slot_half_deg", "r_wall")
D = len(VARS)
REF_POINT = np.array([-0.8, 40.0])                        # hypervolume reference in the minimised form (-T [N.m], ripple [%])
T_FLOOR = 0.8                                             # the torque floor [N.m]
BUDGET = 2000                                             # most designs that may reach FEMM

# ---------------------------------------------------------------------------------------------- the method's settings
BATCH = 24                                                # designs per round (a multiple of the 8 workers)
N_SOBOL = 128                                             # the first space-filling block (a power of two: Sobol' balance)
SEED_BLOCK = 64                                           # an extra Sobol' block when the first model is not good enough
R2_NEEDED = 0.70                                          # ... "good enough": leave-one-out R^2 of the ripple at T >= T_FLOOR
N_SAMPLES = 48                                            # posterior scenarios per candidate in the expected-improvement estimate
N_EXPLORE = 3                                             # of each batch: the most uncertain ripple inside the torque region
N_GLOBAL = 2048                                           # Sobol' candidates in a round's pool
LOCAL_SCALES = (0.01, 0.03, 0.08)                         # sizes of the local clouds around the front (fraction of the unit cube)
LOCAL_PER_SCALE = 8
MAX_FRONT_CENTRES = 60                                    # front members that get a cloud
EXCLUDE = 0.05                                            # two designs of a batch are at least this far apart (normalised box, Euclid)
MAX_FIT = 600                                             # hyperparameters are fitted on at most this many rows
REFIT_GROWTH = 1.3                                        # ... again when the data has grown by this factor
TRAIN_FLOOR = {"T": 0.0, "H": 0.55}                       # rows below these torques do not train the torque model (none) / the harmonic models
UPS = 8                                                   # zero-padding of the waveform in the surrogate (the black box uses 200)
PATIENCE, TOL, MIN_ROUNDS = 8, 1.0e-3, 12                 # convergence: hypervolume gain over PATIENCE rounds below TOL, after MIN_ROUNDS
N_VERIFY = 96                                             # the model-free probe that must find nothing before the run stops
SEED = 20261109

# margins that keep a design strictly inside the feasible region (the screen's own tolerance is 1e-9)
M_BRIDGE, M_RIB, M_TOOTH, M_YOKE = 1e-6, 2e-4, 1e-6, 1e-6

NS = len(design.RIPPLE_POSITIONS)                         # 15 torque samples over one ripple period
NH = NS // 2                                              # 7 harmonics: the 6th, 12th, ..., 42nd of the electrical frequency
NC = 1 + 2 * NH                                           # 15 numbers: the mean, then Re and Im of the 7 harmonics
GROUPS = (("T", (0,)), ("k1", (1, 1 + NH)), ("k2", (2, 2 + NH)), ("k3", (3, 3 + NH)),
          ("k4-7", (4, 5, 6, 7, 4 + NH, 5 + NH, 6 + NH, 7 + NH)))
INIT_ELL = {"T": [3.9, 1.4, 1.6, 15.0, 3.0, 3.0], "k1": [2.0, 0.6, 0.8, 8.0, 2.0, 2.0], "k2": [1.5, 0.35, 0.65, 2.0, 3.0, 2.0],
            "k3": [1.7, 0.25, 0.45, 2.4, 6.0, 8.0], "k4-7": [0.8, 0.3, 0.5, 2.0, 3.0, 2.0]}
INIT_NOISE = {"T": 1e-5, "k1": 1e-3, "k2": 3e-3, "k3": 1e-3, "k4-7": 1e-1}
ELL_BOUNDS = {"T": (0.05, 40.0), "k1": (0.05, 40.0), "k2": (0.05, 40.0), "k3": (0.05, 40.0), "k4-7": (0.2, 40.0)}
NOISE_BOUNDS = {"T": (1e-7, 0.3), "k1": (1e-6, 0.5), "k2": (1e-6, 0.5), "k3": (1e-6, 0.5), "k4-7": (1e-2, 0.9)}

# ---------------------------------------------------------------------------------------------- the box and the feasible cube


def load_box() -> tuple[np.ndarray, np.ndarray]:
    """The 'joint' entry of out\\optimization\\boxes.json as (lo, hi) arrays in the order of VARS."""
    box = json.loads((paths.path("out") / "optimization" / "boxes.json").read_text(encoding="utf-8"))["joint"]
    return np.array([box[v][0] for v in VARS], float), np.array([box[v][1] for v in VARS], float)


LO, HI = load_box()


def out_folder() -> Path:
    return design.out_dir(RUN)


class Cube:
    """The feasible part of the box as a unit cube: u in [0,1]^6 -> a design that the screen always accepts.

    Why: of the 10 constraints of design.screen only three are ever active inside the box, and each is explicit.
      centre_rib   the rib between the two magnets is 2 (mag_y - h(v_half_deg)), exactly affine in mag_y: so mag_y >= ymin(v_half_deg),
                   and ymin is read off the screen's own number at one mag_y for each v_half_deg (a table, one screen per row)
      tooth_min    2 r_shoe sin(5 deg - slot_half_deg) >= 3.0 mm: slot_half_deg <= a_max (a constant)
      yoke         r_so - r_wall (1 + sin a) / cos a >= 3.5 mm: r_wall <= (r_so - 3.5) cos a / (1 + sin a), a function of slot_half_deg
    (checked against design.screen to 1e-14 for the last two, to the table's rounding for the first). Each variable is then
    "its lower end + u times the room that is left", so the corner u = 0 or 1 is exactly the constraint's edge, where the best designs
    sit, and an optimizer or a sampler can never waste a candidate on the infeasible part (a third of the box)."""

    def __init__(self, folder: Path | None = None):
        folder = folder or out_folder()
        ref = design.REFERENCE
        self.a_tooth = 5.0 - math.degrees(math.asin(design.LIMITS["tooth_min"] / (2.0 * ref.r_shoe))) - M_TOOTH
        self.r_so_yoke = ref.r_so - design.LIMITS["yoke"]
        v, ymin = self._rib_table(folder)
        self.v_grid, self.ymin_grid = v, ymin

    @staticmethod
    def _rib_table(folder: Path):
        file = folder / "rib_table.json"
        if file.exists():
            t = json.loads(file.read_text(encoding="utf-8"))
            return np.array(t["v_half_deg"]), np.array(t["mag_y_min"])
        v = np.arange(LO[2], HI[2] + 1e-9, 0.25)
        ymin = []
        for vi in v:
            s = design.screen(design.decode({"bridge_mm": 1.5, "mag_y": 6.0, "v_half_deg": float(vi), **HELD}))
            ymin.append(6.0 - (s["metrics"]["centre_rib"] - design.LIMITS["centre_rib"]) / 2.0)   # rib = rib(6) + 2 (mag_y - 6)
        folder.mkdir(parents=True, exist_ok=True)
        file.write_text(json.dumps({"v_half_deg": v.tolist(), "mag_y_min": ymin,
                                    "note": "mag_y at which the centre rib is exactly its limit, per v_half_deg (from design.screen)"}), encoding="utf-8")
        return v, np.array(ymin)

    def ymin(self, v):
        """The lowest mag_y [mm] the centre rib allows at magnet angle v (interpolated from the table, plus the margin)."""
        return np.maximum(np.interp(v, self.v_grid, self.ymin_grid) + M_RIB, LO[1])

    def r_wall_max(self, a_deg):
        """The deepest r_wall [mm] the yoke allows at slot_half_deg a (and the box's own top, 58 mm)."""
        a = np.radians(a_deg)
        return np.minimum(self.r_so_yoke * np.cos(a) / (1.0 + np.sin(a)) - M_YOKE, HI[5])

    def u_to_x(self, U):
        """Unit-cube points (N, 6) -> designs (N, 6) in the order of VARS: each variable is its lower end plus u times the room that is left."""
        U = np.atleast_2d(np.asarray(U, float))
        X = np.empty_like(U)
        X[:, 0] = LO[0] + M_BRIDGE + U[:, 0] * (HI[0] - LO[0] - M_BRIDGE)
        X[:, 2] = LO[2] + U[:, 2] * (HI[2] - LO[2])
        y0 = self.ymin(X[:, 2])
        X[:, 1] = y0 + U[:, 1] * (HI[1] - y0)
        X[:, 3] = LO[3] + U[:, 3] * (HI[3] - LO[3])
        a_hi = min(HI[4], self.a_tooth)
        X[:, 4] = LO[4] + U[:, 4] * (a_hi - LO[4])
        r_hi = self.r_wall_max(X[:, 4])
        X[:, 5] = LO[5] + U[:, 5] * (r_hi - LO[5])
        return X

    def x_to_u(self, X):
        """The inverse of u_to_x (clipped to the cube: a design a hair outside, by a margin, lands on the face)."""
        X = np.atleast_2d(np.asarray(X, float))
        U = np.empty_like(X)
        U[:, 0] = (X[:, 0] - LO[0] - M_BRIDGE) / (HI[0] - LO[0] - M_BRIDGE)
        U[:, 2] = (X[:, 2] - LO[2]) / (HI[2] - LO[2])
        y0 = self.ymin(X[:, 2])
        U[:, 1] = (X[:, 1] - y0) / (HI[1] - y0)
        U[:, 3] = (X[:, 3] - LO[3]) / (HI[3] - LO[3])
        a_hi = min(HI[4], self.a_tooth)
        U[:, 4] = (X[:, 4] - LO[4]) / (a_hi - LO[4])
        U[:, 5] = (X[:, 5] - LO[5]) / (self.r_wall_max(X[:, 4]) - LO[5])
        return np.clip(U, 0.0, 1.0)


def normalise(X):
    """The box scaled to the unit cube (the GP's inputs): (x - lo) / (hi - lo)."""
    return (np.asarray(X, float) - LO) / (HI - LO)


def candidate(x) -> dict:
    """One design as the dict design.Evaluator takes: the six variables by name, and the held h_open."""
    return {**HELD, **{v: float(x[i]) for i, v in enumerate(VARS)}}


def xkey(x) -> tuple:
    """What counts as 'the same design': the six variables rounded to 5 decimals (the Evaluator's own cache rounds to 6)."""
    return tuple(np.round(np.asarray(x, float), 5))


# ---------------------------------------------------------------------------------------------- the waveform and its objectives
def to_coeffs(T) -> np.ndarray:
    """15 torque samples (one ripple period) -> 15 numbers [mean, Re c1..c7, Im c1..c7]: the discrete Fourier series, c_k = mean(T e^{-2 pi i k n / 15})."""
    T = np.asarray(T, float)
    c = np.fft.rfft(T, axis=-1) / NS
    return np.concatenate([c[..., :1].real, c[..., 1:].real, c[..., 1:].imag], axis=-1)


def objectives_of(C, ups: int = UPS):
    """(T_mean, ripple %) of coefficient vectors C (..., 15): the waveform rebuilt from its harmonics on NS x ups points,
    exactly design.fourier_pkpk's way (ups = 200 there, 8 here: 0.05 points at most), and its peak-to-peak in per cent of the mean."""
    C = np.asarray(C, float)
    spec = np.zeros(C.shape[:-1] + (NH + 1,), complex)
    spec[..., 0] = C[..., 0]
    spec[..., 1:] = C[..., 1:1 + NH] + 1j * C[..., 1 + NH:]
    m = NS * ups
    wave = np.fft.irfft(spec * m, n=m, axis=-1)
    pp = wave.max(-1) - wave.min(-1)
    t = C[..., 0]
    return t, np.where(t > 1e-3, 100.0 * pp / np.maximum(t, 1e-3), 1e3)


# ---------------------------------------------------------------------------------------------- the hypervolume (two objectives, minimised)
def front_inside(F, ref=REF_POINT) -> np.ndarray:
    """The non-dominated points of F that lie inside the reference box, sorted by the first objective (so the second falls)."""
    F = np.asarray(F, float).reshape(-1, 2)
    F = F[(F[:, 0] < ref[0]) & (F[:, 1] < ref[1])]
    if len(F) == 0:
        return F
    F = F[design.nondominated(F)]
    return F[np.argsort(F[:, 0])]


def hv2d(F, ref=REF_POINT) -> float:
    """The area dominated by F inside the box up to ref: the sum of the strips between neighbouring points of the staircase."""
    P = front_inside(F, ref)
    if len(P) == 0:
        return 0.0
    right = np.r_[P[1:, 0], ref[0]]
    return float(((right - P[:, 0]) * (ref[1] - P[:, 1])).sum())


def hvi(points, front, ref=REF_POINT):
    """Hypervolume that each point (..., 2) would add to the sorted staircase `front` (m, 2): the area it dominates and the front does not.
    The staircase cuts the plane into m + 1 vertical strips with a fixed ripple level in each; inside a strip the point adds
    (width of the strip that lies right of the point) x (level - the point's ripple, if positive)."""
    a, b = front[:, 0], front[:, 1]
    lo = np.r_[-np.inf, a]
    hi = np.r_[a, ref[0]]
    level = np.r_[ref[1], b]
    p1, p2 = points[..., 0:1], points[..., 1:2]
    width = np.clip(hi - np.maximum(lo, p1), 0.0, None)
    height = np.clip(level - p2, 0.0, None)
    return (width * height).sum(-1)


# ---------------------------------------------------------------------------------------------- the Gaussian process

def sqdist(A, B, ell):
    """r2[i, j] = sum_d ((A_id - B_jd) / ell_d)^2 (the squared distance in units of the length scales)."""
    As, Bs = A / ell, B / ell
    r2 = (As ** 2).sum(1)[:, None] + (Bs ** 2).sum(1)[None, :] - 2.0 * As @ Bs.T
    return np.maximum(r2, 0.0)


class SEGP:
    """A Gaussian process with the squared-exponential kernel k(x, x') = exp(-1/2 sum_d ((x_d - x'_d) / ell_d)^2) (one length scale per
    variable: automatic relevance determination), unit signal variance, and noise tau^2 on the diagonal, shared by several outputs
    that have been standardised to zero mean and unit variance. The log marginal likelihood of m outputs y_j with the same K is
        -1/2 sum_j y_j' K^-1 y_j - m/2 log|K| - m n/2 log 2 pi,
    and its gradient with respect to a parameter t is 1/2 tr[(sum_j a_j a_j' - m K^-1) dK/dt] with a_j = K^-1 y_j, where
    dK/d log ell_d = K0 * (x_id - x_jd)^2 / ell_d^2 and dK/d log tau^2 = tau^2 I. L-BFGS-B climbs it from the starting values."""

    def __init__(self, ell, noise):
        self.ell = np.asarray(ell, float)
        self.noise = float(noise)
        self.X = self.c = self.alpha = None

    def _nll(self, theta, X, Y):
        d = X.shape[1]
        ell, tau2 = np.exp(theta[:d]), np.exp(theta[d])
        K0 = np.exp(-0.5 * sqdist(X, X, ell))
        n, m = Y.shape
        try:
            c = cho_factor(K0 + (tau2 + 1e-9) * np.eye(n), lower=True)
        except np.linalg.LinAlgError:
            return 1e12, np.zeros_like(theta)
        alpha = cho_solve(c, Y)
        nll = 0.5 * (Y * alpha).sum() + m * np.log(np.diag(c[0])).sum() + 0.5 * m * n * np.log(2 * np.pi)
        W = alpha @ alpha.T - m * cho_solve(c, np.eye(n))
        grad = np.empty_like(theta)
        for k in range(d):
            grad[k] = -0.5 * (W * K0 * (X[:, k][:, None] - X[:, k][None, :]) ** 2).sum() / ell[k] ** 2
        grad[d] = -0.5 * np.trace(W) * tau2
        return nll, grad

    def fit(self, X, Y, ell_bounds, noise_bounds, restarts=1, rng=None, maxiter=60):
        """Maximise the marginal likelihood from the current values (and restarts-1 jittered ones); keep the best."""
        rng = rng or np.random.default_rng(0)
        d = X.shape[1]
        bounds = [(np.log(ell_bounds[0]), np.log(ell_bounds[1]))] * d + [(np.log(noise_bounds[0]), np.log(noise_bounds[1]))]
        start = np.r_[np.log(self.ell), np.log(self.noise)]
        best = None
        for r in range(restarts):
            x0 = start if r == 0 else np.clip(start + np.r_[rng.normal(0, 0.5, d), rng.normal(0, 1.0)], [b[0] for b in bounds], [b[1] for b in bounds])
            res = minimize(self._nll, x0, args=(X, Y), jac=True, method="L-BFGS-B", bounds=bounds, options={"maxiter": maxiter})
            if best is None or res.fun < best.fun:
                best = res
        self.ell, self.noise = np.exp(best.x[:d]), float(np.exp(best.x[d]))
        return self

    def set_data(self, X, Y):
        """Factor K for these rows with the current hyperparameters (a bad factorisation is retried with more noise)."""
        self.X, self.Y = X, Y
        K0 = np.exp(-0.5 * sqdist(X, X, self.ell))
        extra = 1e-9
        for _ in range(8):
            try:
                self.c = cho_factor(K0 + (self.noise + extra) * np.eye(len(X)), lower=True)
                break
            except np.linalg.LinAlgError:
                extra *= 10.0
        else:
            raise np.linalg.LinAlgError("kernel matrix is not positive definite")
        self.alpha = cho_solve(self.c, Y)
        return self

    def predict(self, Xs):
        """Mean (N, m) and predictive variance (N,) of a new noisy observation (standardised units)."""
        Ks = np.exp(-0.5 * sqdist(Xs, self.X, self.ell))
        v = solve_triangular(self.c[0], Ks.T, lower=True)
        return Ks @ self.alpha, np.maximum(1.0 + self.noise - (v ** 2).sum(0), self.noise)

    def loo(self):
        """Leave-one-out mean (n, m) and variance (n,) for every training row, in closed form from K^-1."""
        kinv = cho_solve(self.c, np.eye(len(self.X)))
        d = np.diag(kinv)
        return self.Y - self.alpha / d[:, None], 1.0 / d


class Surrogate:
    """The torque waveform of a design as five Gaussian processes over the 15 numbers of to_coeffs: the mean torque (T), the 6th harmonic (k1),
    the 12th (k2), the 18th (k3: the slot harmonic, the biggest) and the 24th to 42nd together (k4-7). Each group has its own length scales
    and noise (one kernel for all the harmonics gave 0.89 points of leave-one-out ripple error on the rotor sweeps, one per group 0.44, the
    squared-exponential kernel 0.29), and the ripple is predicted through the same peak-to-peak rule as the black box's, so the V-shaped valleys
    of the ripple come from the rule, not from the GP."""

    def __init__(self):
        self.gp = {name: SEGP(INIT_ELL[name], INIT_NOISE[name]) for name, _ in GROUPS}
        self.n_fit = 0

    def _rows(self, name, Tm):
        return Tm >= TRAIN_FLOOR["T" if name == "T" else "H"]

    def fit(self, X, C, Tm, refit: bool, rng):
        Xn = normalise(X)
        self.mu, self.sd, self.keep = {}, {}, {}
        for name, idx in GROUPS:
            keep = self._rows(name, Tm)
            idx = list(idx)
            Yg = C[keep][:, idx]
            self.mu[name], self.sd[name] = Yg.mean(0), Yg.std(0) + 1e-12
            Z = (Yg - self.mu[name]) / self.sd[name]
            gp = self.gp[name]
            if refit:
                sel = np.arange(len(Z))
                if len(sel) > MAX_FIT:
                    sel = np.sort(rng.choice(len(Z), MAX_FIT, replace=False))
                gp.fit(Xn[keep][sel], Z[sel], ELL_BOUNDS[name], NOISE_BOUNDS[name], restarts=2 if self.n_fit == 0 else 1, rng=rng)
            gp.set_data(Xn[keep], Z)
            self.keep[name] = keep
        if refit:
            self.n_fit = len(X)

    def predict_mean(self, X, chunk=2000):
        """Only the means of the 15 numbers (no variances: a model optimizer asks for hundreds of thousands of these)."""
        Xn = normalise(X)
        mean = np.empty((len(Xn), NC))
        for name, idx in GROUPS:
            gp = self.gp[name]
            for i in range(0, len(Xn), chunk):
                Ks = np.exp(-0.5 * sqdist(Xn[i:i + chunk], gp.X, gp.ell))
                mean[i:i + chunk, list(idx)] = (Ks @ gp.alpha) * self.sd[name] + self.mu[name]
        return mean

    def describe(self) -> str:
        """The fitted length scales and noise of each group, for the log (box units: 1 is the whole range of a variable)."""
        return "; ".join(f"{name}: ell {np.round(self.gp[name].ell, 2).tolist()} noise {self.gp[name].noise:.1e}" for name, _ in GROUPS)

    def predict(self, X, chunk=2000):
        """Mean and standard deviation of the 15 numbers at designs X (N, 6), in the black box's units."""
        Xn = normalise(X)
        mean = np.empty((len(Xn), NC))
        sd = np.empty((len(Xn), NC))
        for name, idx in GROUPS:
            idx = list(idx)
            for i in range(0, len(Xn), chunk):
                m, v = self.gp[name].predict(Xn[i:i + chunk])
                mean[i:i + chunk, idx] = m * self.sd[name] + self.mu[name]
                sd[i:i + chunk, idx] = np.sqrt(v)[:, None] * self.sd[name]
        return mean, sd

    def loo_report(self, X, C, Tm, R):
        """How well the surrogate predicts a design it has not seen: the leave-one-out error of the mean torque (N.m) and of the ripple (points)
        on the rows at or above the torque floor, the ripple rebuilt from the leave-one-out harmonics and the observed mean torque."""
        pred = C.copy()
        for name, idx in GROUPS:
            m, _ = self.gp[name].loo()
            m = m * self.sd[name] + self.mu[name]
            keep = self.keep[name]
            rows = np.flatnonzero(keep)
            if name == "T":
                pred_t = np.full(len(C), np.nan)
                pred_t[rows] = m[:, 0]
            else:
                pred[np.ix_(rows, list(idx))] = m
        ok = (Tm >= T_FLOOR) & self.keep["k3"] & self.keep["T"]
        _, rp = objectives_of(pred)
        e_t = pred_t[ok] - Tm[ok]
        e_r = rp[ok] - R[ok]
        var_r = R[ok].var()
        return {"n": int(ok.sum()), "T_rmse": float(np.sqrt(np.mean(e_t ** 2))), "R_rmse": float(np.sqrt(np.mean(e_r ** 2))),
                "R_r2": float(1.0 - np.mean(e_r ** 2) / var_r) if var_r > 0 else float("nan"), "R_sd": float(np.sqrt(var_r))}


# ---------------------------------------------------------------------------------------------- the scored designs
class Dataset:
    """Every design that reached FEMM and came back with a waveform (X, C, T, R), plus the bookkeeping: which designs were tried (also those
    that failed), and how many reached FEMM (the budget)."""

    def __init__(self):
        self.X = np.empty((0, D))
        self.C = np.empty((0, NC))
        self.T = np.empty(0)
        self.R = np.empty(0)
        self.keys: set = set()
        self.reached = 0

    def add(self, x, waveform, t_mean, ripple):
        self.X = np.vstack([self.X, x])
        self.C = np.vstack([self.C, to_coeffs(waveform)])
        self.T = np.append(self.T, t_mean)
        self.R = np.append(self.R, ripple)

    def add_rows(self, rows) -> None:
        """Rows as design.Evaluator returns them."""
        for r in rows:
            p = r["params"]
            x = np.array([design.bridge_of(p), p.mag_y, p.v_half_deg, p.slot_open, p.slot_half_deg, p.r_wall])
            if r.get("cached"):
                continue
            self.keys.add(xkey(x))
            if "seconds" in r or str(r.get("reason", "")).startswith("FEMM:"):
                self.reached += 1
            if r["feasible"] and "T" in r and np.isfinite(r["T_mean"]) and np.isfinite(r["ripple_pct"]):
                self.add(x, r["T"], r["T_mean"], r["ripple_pct"])

    def load_csv(self, file: Path) -> None:
        """The rows already in evaluations.csv (a run that is resumed): same bookkeeping as add_rows."""
        df = pd.read_csv(file)
        df = df[df["cached"].astype(str) != "True"]
        for _, r in df.iterrows():
            x = np.array([r[v] for v in VARS], float)
            self.keys.add(xkey(x))
            reason = str(r["reason"]) if pd.notna(r["reason"]) else ""
            if pd.notna(r["seconds"]) or reason.startswith("FEMM:"):
                self.reached += 1
            if str(r["feasible"]) == "True" and pd.notna(r["T_mean"]) and pd.notna(r["ripple_pct"]):
                self.add(x, json.loads(r["T"]), float(r["T_mean"]), float(r["ripple_pct"]))

    @property
    def admissible(self) -> np.ndarray:
        return self.T >= T_FLOOR

    def objectives(self) -> np.ndarray:
        return np.c_[-self.T, self.R]

    def front(self) -> np.ndarray:
        """Indices (into X) of the non-dominated admissible designs."""
        idx = np.flatnonzero(self.admissible)
        if len(idx) == 0:
            return idx
        return idx[design.nondominated(self.objectives()[idx])]

    def hv(self) -> float:
        idx = np.flatnonzero(self.admissible)
        return hv2d(self.objectives()[idx]) if len(idx) else 0.0


# ---------------------------------------------------------------------------------------------- candidates: where the surrogate may look
def sobol_u(n, seed) -> np.ndarray:
    return qmc.Sobol(D, scramble=True, seed=seed).random(n)


def local_u(centres, scales, per_scale, rng) -> np.ndarray:
    """Clouds around each centre in the unit cube: Gaussian steps of each size, each coordinate moved with probability 1/2 (at least one),
    clipped to the cube, so a centre on a face has half its cloud on the face."""
    out = []
    for s in scales:
        for _ in range(per_scale):
            step = rng.normal(0.0, s, centres.shape)
            mask = rng.random(centres.shape) < 0.5
            mask[np.arange(len(centres)), rng.integers(0, D, len(centres))] = True
            out.append(np.clip(centres + step * mask, 0.0, 1.0))
    return np.vstack(out) if out else np.empty((0, D))


def face_snaps(centres) -> np.ndarray:
    """Each centre with one coordinate pushed to 0 or to 1 (the best designs of a monotone problem sit on faces)."""
    out = []
    for i in range(D):
        for face in (0.0, 1.0):
            c = centres.copy()
            c[:, i] = face
            out.append(c)
    return np.vstack(out)


def front_centres(data: Dataset, cube: Cube, limit=MAX_FRONT_CENTRES) -> np.ndarray:
    idx = data.front()
    if len(idx) == 0:                                  # no admissible design yet: the best torques
        idx = np.argsort(-data.T)[:limit]
    if len(idx) > limit:
        idx = idx[np.round(np.linspace(0, len(idx) - 1, limit)).astype(int)]
    return cube.x_to_u(data.X[idx])


def make_pool(data: Dataset, cube: Cube, rng, seed) -> np.ndarray:
    """The candidate designs of one round, as unit-cube points: a fresh Sobol' set (look everywhere) and clouds and face-snaps around the front."""
    centres = front_centres(data, cube)
    U = np.vstack([sobol_u(N_GLOBAL, seed), local_u(centres, LOCAL_SCALES, LOCAL_PER_SCALE, rng), face_snaps(centres)])
    return U


def drop_known(U, cube: Cube, data: Dataset):
    X = cube.u_to_x(U)
    fresh = np.array([xkey(x) not in data.keys for x in X])
    return U[fresh], X[fresh]


# ---------------------------------------------------------------------------------------------- the acquisition: expected hypervolume improvement
def sample_objectives(mean, sd, eps):
    """For every candidate, one scenario per row of eps: the 15 numbers drawn from their predictive normal (mean + sd eps), turned into
    (-T, ripple). Returns (N, S, 2). The same eps for all candidates (a fixed set of scenarios: nearby candidates share their luck)."""
    n, s = len(mean), len(eps)
    out = np.empty((n, s, 2))
    for i in range(0, n, 300):
        Y = mean[i:i + 300, None, :] + sd[i:i + 300, None, :] * eps[None, :, :]
        t, r = objectives_of(Y)
        out[i:i + 300, :, 0] = -t
        out[i:i + 300, :, 1] = r
    return out


def expected_hvi(scen, front, chunk=400):
    """Mean over scenarios of the hypervolume each candidate would add to `front`: the expected improvement of the front (EHVI)."""
    out = np.empty(len(scen))
    for i in range(0, len(scen), chunk):
        out[i:i + chunk] = hvi(scen[i:i + chunk], front).mean(-1)
    return out


def select_batch(mean, sd, Xn_pool, front, q_ehvi, q_explore, rng, t_gate=0.85):
    """Greedy batch: q_ehvi times, take the candidate with the largest expected improvement of the front, add its predicted outcome to
    the front (so the next one must improve on what the last promises), and rule out everything within EXCLUDE of it; then q_explore
    times, the candidate whose predicted ripple is the most uncertain among those expected to clear the torque floor with room.
    Returns (indices, EHVI of each pick)."""
    eps = rng.standard_normal((N_SAMPLES, NC))
    scen = sample_objectives(mean, sd, eps)
    t_hat, r_hat = objectives_of(mean)
    believed = np.c_[-t_hat, r_hat]
    free = np.ones(len(mean), bool)
    chosen, gains = [], []
    fr = front.copy()
    for _ in range(q_ehvi):
        e = np.where(free, expected_hvi(scen, fr), -1.0)
        i = int(np.argmax(e))
        if e[i] <= 1e-12:                              # nothing left that the model expects to improve the front
            break
        chosen.append(i)
        gains.append(float(e[i]))
        free &= np.linalg.norm(Xn_pool - Xn_pool[i], axis=1) > EXCLUDE
        fr = front_inside(np.vstack([fr, believed[i]]))
    spread = scen[:, :, 1].std(1)                     # how unsure the model is of each candidate's ripple
    ok = free & (t_hat >= t_gate)
    for _ in range(q_ehvi + q_explore - len(chosen)):  # the exploration picks, and any EHVI pick that could not be made
        if not ok.any():
            break
        i = int(np.argmax(np.where(ok, spread, -1.0)))
        chosen.append(i)
        gains.append(0.0)
        ok &= np.linalg.norm(Xn_pool - Xn_pool[i], axis=1) > EXCLUDE
    return np.array(chosen, int), np.array(gains)


# ---------------------------------------------------------------------------------------------- the run
def say(msg: str = "") -> None:
    print(msg, flush=True)


PROGRESS_COLS = ["round", "stage", "reached", "admissible", "T_best", "ripple_min", "front", "hv", "hv_gain_pct", "acq_max", "loo_T_rmse",
                 "loo_R_rmse", "loo_R2", "pred_T_rmse", "pred_R_rmse", "pred_R_bias", "round_s", "total_s"]


def prospective_errors(rows, mean_pick) -> dict | None:
    """How far the surrogate was from FEMM on the designs it had just chosen (it had not seen them): rms error of the mean torque (N.m) and of the
    ripple (points), and the ripple's mean error observed - predicted (positive: the model was too hopeful about the ripple)."""
    t_hat, r_hat = objectives_of(mean_pick)
    got = [(i, r["T_mean"], r["ripple_pct"]) for i, r in enumerate(rows) if r.get("feasible") and "T_mean" in r and np.isfinite(r["T_mean"])]
    if not got:
        return None
    i = np.array([g[0] for g in got])
    e_t = np.array([g[1] for g in got]) - t_hat[i]
    e_r = np.array([g[2] for g in got]) - r_hat[i]
    return {"T_rmse": float(np.sqrt(np.mean(e_t ** 2))), "R_rmse": float(np.sqrt(np.mean(e_r ** 2))), "R_bias": float(np.mean(e_r))}


def seed_designs(cube: Cube) -> list[np.ndarray]:
    """Designs taken from the one-at-a-time sweeps' fronts (rotor values; bridge on the guard) with two stators: the reference's and a big
    one (the widest slot the tooth and yoke limits allow, narrowest opening). They go in the first block, so the model sees the front's region
    from the start."""
    rotor = [(89.0, 5.225), (89.0, 5.55), (59.5, 4.90), (54.58, 4.90), (54.0, 4.90), (79.17, 5.875), (69.33, 6.2), (84.0, 7.0), (45.0, 4.40), (60.0, 6.0)]
    stators = [(2.0, 2.57, 55.744), (1.5, 3.1, 57.0)]
    out = []
    for v, y in rotor:
        for so, a, rw in stators:
            out.append(np.array([LO[0] + M_BRIDGE, max(y, float(cube.ymin(v))), v, so, a, rw]))
    return out


def run(args) -> None:
    folder = out_folder()
    t_start = time.perf_counter()
    cube = Cube(folder)
    say(f"cube ready: tooth limit {cube.a_tooth:.5f} deg, rib table {len(cube.v_grid)} rows")
    data = Dataset()
    csv_file = folder / "evaluations.csv"
    if csv_file.exists():
        data.load_csv(csv_file)
        say(f"resumed from {csv_file.name}: {data.reached} designs reached FEMM, {len(data.T)} scored")
    prog_file = folder / "progress.csv"
    if not prog_file.exists():
        prog_file.write_text(",".join(PROGRESS_COLS) + "\n", encoding="utf-8")
    t_wait = time.perf_counter()
    with femm_queue.femm_pool("optimization: crosscheck", task_timeout_s=design.TASK_TIMEOUT_S, poll_s=3600.0) as pool:
        waited = time.perf_counter() - t_wait
        say(f"pool open after {waited:.0f} s of waiting for the slot")
        ev = design.Evaluator(pool, run=RUN, scoring=SCORING)
        optimize(ev, cube, data, prog_file, args, t_start, waited)
    report(folder)


def log_round(prog_file: Path, row: dict) -> None:
    with open(prog_file, "a", encoding="utf-8") as f:
        f.write(",".join(str(row.get(c, "")) for c in PROGRESS_COLS) + "\n")


def record(data: Dataset, prog_file: Path, hv_hist: list, rnd: int, stage: str, t0: float, t_start: float, acq=float("nan"), loo=None, pred=None) -> None:
    """One line of progress.csv and of the log: the state of the front after a round (stage design / infill / verify / probe-...)."""
    idx = data.front()
    hv = data.hv()
    gain = 100.0 * (hv - hv_hist[-1]) / hv_hist[-1] if hv_hist and hv_hist[-1] > 0 else float("nan")
    adm = data.admissible
    row = {"round": rnd, "stage": stage, "reached": data.reached, "admissible": int(adm.sum()),
           "T_best": round(float(data.T.max()), 4) if len(data.T) else "", "ripple_min": round(float(data.R[adm].min()), 3) if adm.any() else "",
           "front": len(idx), "hv": round(hv, 4), "hv_gain_pct": round(gain, 3) if gain == gain else "", "acq_max": round(acq, 5) if acq == acq else "",
           "loo_T_rmse": round(loo["T_rmse"], 4) if loo else "", "loo_R_rmse": round(loo["R_rmse"], 3) if loo else "",
           "loo_R2": round(loo["R_r2"], 3) if loo else "",
           "pred_T_rmse": round(pred["T_rmse"], 4) if pred else "", "pred_R_rmse": round(pred["R_rmse"], 3) if pred else "",
           "pred_R_bias": round(pred["R_bias"], 3) if pred else "",
           "round_s": round(time.perf_counter() - t0, 1), "total_s": round(time.perf_counter() - t_start, 1)}
    hv_hist.append(hv)
    log_round(prog_file, row)
    say(f"[round {rnd:3d} {stage:7s}] reached {data.reached:4d}  admissible {row['admissible']:4d}  T_best {row['T_best']}  ripple_min {row['ripple_min']}  "
        f"front {row['front']:3d}  HV {hv:8.3f} ({row['hv_gain_pct']}%)  acq {row['acq_max']}  LOO T {row['loo_T_rmse']} R {row['loo_R_rmse']} (R2 {row['loo_R2']})  "
        f"picks: T err {row['pred_T_rmse']} R err {row['pred_R_rmse']} bias {row['pred_R_bias']}  {row['round_s']} s")


def optimize(ev, cube: Cube, data: Dataset, prog_file: Path, args, t_start: float, waited: float = 0.0) -> None:
    budget, batch = args.budget, args.batch
    t_open = time.perf_counter()
    prog_rows = pd.read_csv(prog_file)
    rnd = int(prog_rows["round"].max()) + 1 if len(prog_rows) else 0
    hv_hist = list(prog_rows["hv"].astype(float)) if len(prog_rows) else []
    surrogate = Surrogate()

    def snapshot(stage, t0, acq=float("nan"), loo=None, pred=None):
        record(data, prog_file, hv_hist, rnd, stage, t0, t_start, acq, loo, pred)

    def asked_to_stop() -> bool:
        """A file named STOP in the output folder, or --max-minutes of pool time used up: finish the round in hand and stop."""
        if (prog_file.parent / "STOP").exists():
            say("STOP file found: stopping after this round")
            return True
        if getattr(args, "max_minutes", None) and (time.perf_counter() - t_open) / 60.0 > args.max_minutes:
            say(f"--max-minutes {args.max_minutes} used up: stopping after this round")
            return True
        return False

    # ---- block 1: space-filling designs and the sweeps' seeds (feasible by construction)
    if data.reached == 0:
        t0 = time.perf_counter()
        U = sobol_u(N_SOBOL, SEED)
        cands = [candidate(x) for x in cube.u_to_x(U)] + [candidate(x) for x in seed_designs(cube)]
        for i in range(0, len(cands), 48):
            data.add_rows(ev(cands[i:i + 48], stage="design", gen=rnd))
        snapshot("design", t0)
        rnd += 1
    # ---- is the first model good enough? if not, more space-filling designs (at most two blocks)
    extra = 0
    rng = np.random.default_rng([SEED, 99])
    last_refit_n = 0
    while len(data.T) >= 40:
        surrogate.fit(data.X, data.C, data.T, refit=True, rng=rng)
        last_refit_n = len(data.T)
        loo = surrogate.loo_report(data.X, data.C, data.T, data.R)
        say(f"first model on {len(data.T)} designs: leave-one-out on {loo['n']} rows at T >= {T_FLOOR}: T rmse {loo['T_rmse']:.4f} N.m, "
            f"ripple rmse {loo['R_rmse']:.2f} points (sd {loo['R_sd']:.2f}), R2 {loo['R_r2']:.3f}")
        say("  length scales (box units) " + surrogate.describe())
        if loo["R_r2"] >= R2_NEEDED or extra >= 2 or data.reached + SEED_BLOCK > budget:
            break
        t0 = time.perf_counter()
        fresh = drop_known(sobol_u(SEED_BLOCK, SEED + 1 + extra), cube, data)[1]
        data.add_rows(ev([candidate(x) for x in fresh], stage="design", gen=rnd))
        snapshot("design", t0)
        rnd += 1
        extra += 1
    # ---- the rounds: fit, propose, score
    n_infill, loo = 0, None
    while data.reached < budget and not asked_to_stop():
        t0 = time.perf_counter()
        q = min(batch, budget - data.reached)
        rng = np.random.default_rng([SEED, rnd])
        refit = len(data.T) >= REFIT_GROWTH * last_refit_n
        surrogate.fit(data.X, data.C, data.T, refit=refit, rng=rng)
        if refit:
            last_refit_n = len(data.T)
            say("  refit, length scales (box units) " + surrogate.describe())
        loo = surrogate.loo_report(data.X, data.C, data.T, data.R) if (refit or n_infill % 4 == 0) else None
        pool_u = make_pool(data, cube, rng, SEED + 1000 + rnd)
        if getattr(args, "de_every", 0) and n_infill % args.de_every == 0:     # the surrogate's own best designs join the pool (see model_optima)
            pool_u = np.vstack([pool_u, np.array([u for _, u in model_optima(surrogate, cube, seeds=(1 + n_infill,), per_run=2)])])
        U, X = drop_known(pool_u, cube, data)
        mean, sd = surrogate.predict(X)
        fr = front_inside(data.objectives()[data.admissible]) if data.admissible.any() else np.empty((0, 2))
        n_explore = min(N_EXPLORE, q)
        pick, gains = select_batch(mean, sd, normalise(X), fr, q - n_explore, n_explore, rng)
        if len(pick) == 0:
            say("no candidate left to propose: stopping")
            break
        rows = ev([candidate(X[i]) for i in pick], stage="infill", gen=rnd)
        data.add_rows(rows)
        snapshot("infill", t0, acq=float(gains.max()) if len(gains) else float("nan"), loo=loo, pred=prospective_errors(rows, mean[pick]))
        rnd += 1
        n_infill += 1
        # convergence: the hypervolume has not moved for PATIENCE rounds, then a probe that does not use the model must find nothing either
        if n_infill >= MIN_ROUNDS and len(hv_hist) > PATIENCE and (hv_hist[-1] - hv_hist[-1 - PATIENCE]) <= TOL * hv_hist[-1] and not args.no_verify:
            if data.reached + N_VERIFY > budget:
                break
            t0 = time.perf_counter()
            hv_before = hv_hist[-1]
            surrogate.fit(data.X, data.C, data.T, refit=False, rng=rng)
            Xv = verification_points(data, cube, surrogate, np.random.default_rng([SEED, 7000 + rnd]))
            data.add_rows(ev([candidate(x) for x in Xv], stage="verify", gen=rnd))
            snapshot("verify", t0)
            rnd += 1
            if hv_hist[-1] - hv_before <= TOL * hv_hist[-1]:
                say(f"converged: no gain above {TOL:.1%} over {PATIENCE} rounds, and the model-free probe of {len(Xv)} designs found none either")
                break
            say("the probe found an improvement: the run goes on (the patience window now contains it, so it starts again)")
    say(f"stopped with {data.reached} designs reached FEMM ({data.reached / budget:.0%} of the budget), {time.perf_counter() - t_start:.0f} s, "
        f"{waited:.0f} s of it waiting for the slot")


def verification_points(data: Dataset, cube: Cube, surrogate: Surrogate, rng) -> np.ndarray:
    """N_VERIFY designs that the optimizer's surrogate did not choose: half fresh Sobol' points inside the torque region (the surrogate is used
    only to keep the probe where T can clear the floor), half random perturbations of front members at sizes from 0.005 to 0.1 of the cube."""
    n_half = N_VERIFY // 2
    U = sobol_u(4096, int(rng.integers(1, 2 ** 31 - 1)))
    X = cube.u_to_x(U)
    t_hat = surrogate.predict(X)[0][:, 0]
    glob = X[t_hat >= 0.85][:n_half]
    centres = front_centres(data, cube)
    pick = centres[rng.integers(0, len(centres), n_half)]
    scale = rng.choice([0.005, 0.01, 0.02, 0.05, 0.1], n_half)[:, None]
    step = rng.normal(0.0, 1.0, pick.shape) * scale
    mask = rng.random(pick.shape) < 0.5
    mask[np.arange(n_half), rng.integers(0, D, n_half)] = True
    loc = cube.u_to_x(np.clip(pick + step * mask, 0.0, 1.0))
    allx = np.vstack([glob, loc])
    keep = np.array([xkey(x) not in data.keys for x in allx])
    return allx[keep]


# ---------------------------------------------------------------------------------------------- after convergence: two looks that the rounds could not take
T_LADDER = (0.9, 1.1, 1.3, 1.5, 1.6, 1.7, 1.8, 1.9)        # torque levels at which the surrogate's lowest ripple is looked for [N.m]
MAP_NV, MAP_NY = 20, 14                                     # the map: 20 magnet angles x 14 distances from the rib limit
MAP_FIXED = {0: 0.0, 3: 0.0, 4: 1.0, 5: 0.81}               # the other four coordinates, at the front's high-torque corner: bridge on its guard,
                                                            # narrowest opening, widest slot (the tooth limit), r_wall 0.81 of its room (56.6 mm)


def model_optima(surrogate: Surrogate, cube: Cube, seeds=(1, 2), per_run=3):
    """The surrogate's own best designs: for each torque level, the lowest ripple the model predicts among designs it predicts at that torque or
    more, found by differential evolution over the whole unit cube (a population of 240 designs, 250 generations; the penalty is 200 points of
    ripple per N.m of missing torque). A different optimizer from the rounds' pool search and with no pool to limit it: the question it asks is
    whether the pool search left a better design lying between its candidates. Returns [(level, u)], at most per_run distinct designs of each run's
    final population."""
    found = []
    for level in T_LADDER:
        def cost(Ut, level=level):                          # scipy hands over the whole population as an array (D, S)
            t, r = objectives_of(surrogate.predict_mean(cube.u_to_x(Ut.T)))
            return r + 200.0 * np.maximum(0.0, level - t)
        for seed in seeds:
            res = differential_evolution(cost, [(0.0, 1.0)] * D, vectorized=True, updating="deferred", popsize=40, maxiter=250, tol=1e-9, polish=False,
                                         init="latinhypercube", strategy="rand1bin", mutation=(0.4, 1.0), recombination=0.8, seed=100 * seed + int(10 * level))
            kept = []
            for k in np.argsort(res.population_energies):
                if all(np.linalg.norm(res.population[k] - w) > 0.03 for w in kept):
                    kept.append(res.population[k])
                if len(kept) == per_run:
                    break
            found += [(level, u) for u in kept]
    return found


def map_points(cube: Cube) -> np.ndarray:
    """The ripple map of the two variables that make its valleys: v_half_deg over its whole range against mag_y from the rib limit to its top, the
    other four at the front's high-torque corner (MAP_FIXED). 20 x 14 designs, every one feasible by construction."""
    U = []
    for uv in np.linspace(0.0, 1.0, MAP_NV):
        for uy in np.linspace(0.0, 1.0, MAP_NY):
            u = np.zeros(D)
            u[1], u[2] = uy, uv
            for i, val in MAP_FIXED.items():
                u[i] = val
            U.append(u)
    return cube.u_to_x(np.array(U))


def run_probe(ev, cube: Cube, data: Dataset, prog_file: Path, t_start: float) -> None:
    """The two looks, after the rounds have converged (they use the converged data and add to it):
    A. model optima: differential evolution on the fitted surrogate at each torque level, the best designs scored in FEMM. If the pool search had
       left a better design between its candidates, the model would predict it here; the printed table sets what the model promised against what
       FEMM gave and against what the front already had at that torque.
    B. a map of the ripple valleys at the high-torque corner of the stator, 280 designs, which needs no model: if the surrogate were blind to a valley,
       the map would find it."""
    prog = pd.read_csv(prog_file)
    rnd = int(prog["round"].max()) + 1 if len(prog) else 0
    hv_hist = list(prog["hv"].astype(float)) if len(prog) else []
    rng = np.random.default_rng([SEED, 31337])
    surrogate = Surrogate()
    surrogate.fit(data.X, data.C, data.T, refit=True, rng=rng)
    loo = surrogate.loo_report(data.X, data.C, data.T, data.R)
    say(f"probe: surrogate refitted on {len(data.T)} designs (leave-one-out: T {loo['T_rmse']:.4f} N.m, ripple {loo['R_rmse']:.2f} points, R2 {loo['R_r2']:.3f})")
    say("  length scales (box units) " + surrogate.describe())
    # ---- A
    t0 = time.perf_counter()
    found = model_optima(surrogate, cube)
    X = cube.u_to_x(np.array([u for _, u in found]))
    levels = np.array([lv for lv, _ in found])
    fresh = np.array([xkey(x) not in data.keys for x in X])
    X, levels = X[fresh], levels[fresh]
    mean = surrogate.predict_mean(X)
    t_hat, r_hat = objectives_of(mean)
    front_before = front_inside(data.objectives()[data.admissible])
    rows = ev([candidate(x) for x in X], stage="probe-model", gen=rnd)
    data.add_rows(rows)
    say(f"A. model optima: {len(X)} new designs; for each torque level the model's best, FEMM's answer, and the front's ripple at that torque or more before this probe")
    for lv in T_LADDER:
        idx = [i for i in range(len(X)) if levels[i] == lv and rows[i].get("feasible") and "T_mean" in rows[i]]
        if not idx:
            continue
        i = min(idx, key=lambda k: r_hat[k])
        before = front_before[front_before[:, 0] <= -lv][:, 1]
        say(f"   T >= {lv:4.2f}: model T {t_hat[i]:.3f} ripple {r_hat[i]:5.2f} %   FEMM T {rows[i]['T_mean']:.3f} ripple {rows[i]['ripple_pct']:5.2f} %   "
            f"front before: {before.min() if len(before) else float('nan'):5.2f} %")
    record(data, prog_file, hv_hist, rnd, "probe-model", t0, t_start, pred=prospective_errors(rows, mean))
    rnd += 1
    # ---- B
    t0 = time.perf_counter()
    Xm = map_points(cube)
    Xm = Xm[np.array([xkey(x) not in data.keys for x in Xm])]
    for i in range(0, len(Xm), 48):
        data.add_rows(ev([candidate(x) for x in Xm[i:i + 48]], stage="probe-map", gen=rnd))
    say(f"B. map of {len(Xm)} designs (v_half_deg x mag_y from the rib limit up, the other four at the high-torque corner)")
    record(data, prog_file, hv_hist, rnd, "probe-map", t0, t_start)


def probe(args) -> None:
    folder = out_folder()
    t_start = time.perf_counter()
    cube = Cube(folder)
    data = Dataset()
    data.load_csv(folder / "evaluations.csv")
    say(f"probe: {data.reached} designs reached FEMM so far, {len(data.T)} scored")
    prog_file = folder / "progress.csv"
    t_wait = time.perf_counter()
    with femm_queue.femm_pool("optimization: crosscheck probe", task_timeout_s=design.TASK_TIMEOUT_S, poll_s=3600.0) as pool:
        say(f"pool open after {time.perf_counter() - t_wait:.0f} s of waiting for the slot")
        run_probe(design.Evaluator(pool, run=RUN, scoring=SCORING), cube, data, prog_file, t_start)
    report(folder)


# ---------------------------------------------------------------------------------------------- the answer
def read_scored(folder: Path) -> pd.DataFrame:
    """Every admissible design in evaluations.csv (feasible, scored, T_mean >= the floor), repeats dropped."""
    df = pd.read_csv(folder / "evaluations.csv")
    df = df[(df["feasible"].astype(str) == "True") & df["T_mean"].notna() & df["ripple_pct"].notna()]
    df = df[df["cached"].astype(str) != "True"].drop_duplicates(subset="overrides")
    return df[df["T_mean"] >= T_FLOOR]


def report(folder: Path | None = None) -> dict:
    """front.csv (the non-dominated admissible designs, every column of evaluations.csv, by falling torque), knee.json, and the numbers."""
    folder = folder or out_folder()
    adm = read_scored(folder)
    F = np.c_[-adm["T_mean"].to_numpy(), adm["ripple_pct"].to_numpy()]
    mask = design.nondominated(F)
    front = adm[mask].sort_values("T_mean", ascending=False)
    Ff = np.c_[-front["T_mean"].to_numpy(), front["ripple_pct"].to_numpy()]
    front.to_csv(folder / "front.csv", index=False)
    k = design.knee(Ff)
    row = front.iloc[k]
    allrows = pd.read_csv(folder / "evaluations.csv")
    reached = int(((allrows["cached"].astype(str) != "True") & (allrows["seconds"].notna() | allrows["reason"].astype(str).str.startswith("FEMM:"))).sum())
    knee = {"overrides": json.loads(row["overrides"]), "variables": {v: float(row[v]) for v in VARS}, "held": HELD,
            "T_mean": float(row["T_mean"]), "ripple_pct": float(row["ripple_pct"]), "gamma_deg": float(row["gamma"]), "pf": float(row["pf"]),
            "amps": float(row["amps"]), "Pcu_rel": float(row["Pcu_rel"]), "saliency": float(row["saliency"]),
            "how": "design.knee on the front's own (-T_mean, ripple_pct), each scaled to 0..1 between the front's best and worst",
            "front_size": int(len(front)), "front_T_range": [float(front["T_mean"].min()), float(front["T_mean"].max())],
            "front_ripple_range": [float(front["ripple_pct"].min()), float(front["ripple_pct"].max())],
            "hypervolume": hv2d(Ff), "hypervolume_reference": REF_POINT.tolist(), "designs_reached_femm": reached,
            "admissible_designs": int(len(adm))}
    (folder / "knee.json").write_text(json.dumps(knee, indent=2), encoding="utf-8")
    say(f"front {len(front)} designs, T {knee['front_T_range'][0]:.4f}..{knee['front_T_range'][1]:.4f} N.m, ripple "
        f"{knee['front_ripple_range'][0]:.3f}..{knee['front_ripple_range'][1]:.3f} %, hypervolume {knee['hypervolume']:.4f}; knee T {knee['T_mean']:.4f} "
        f"ripple {knee['ripple_pct']:.3f} gamma {knee['gamma_deg']:.2f} pf {knee['pf']:.4f}; {reached} designs reached FEMM")
    return knee


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m optimization.crosscheck", description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("run", help="the optimizer: opens the FEMM pool (holds the solver slot) and scores designs until the budget or the convergence test")
    p.add_argument("--budget", type=int, default=BUDGET, help="most designs that may reach FEMM (default 2000)")
    p.add_argument("--batch", type=int, default=BATCH, help="designs per round (default 24)")
    p.add_argument("--no-verify", action="store_true", help="skip the model-free probe before stopping")
    p.add_argument("--max-minutes", type=float, default=None, help="stop after this many minutes with the pool open (default: no limit)")
    p.add_argument("--de-every", type=int, default=0, help="every N rounds, add the surrogate's own optima (model_optima) to the round's candidates (default 0: never)")
    sub.add_parser("probe", help="after convergence: the surrogate's own optima (differential evolution) and a 280-design ripple map, both scored in FEMM")
    sub.add_parser("report", help="front.csv and knee.json from evaluations.csv (no FEMM)")
    args = parser.parse_args(argv)
    if args.cmd == "run":
        run(args)
    elif args.cmd == "probe":
        probe(args)
    else:
        report()
    return 0


if __name__ == "__main__":
    sys.exit(main())
