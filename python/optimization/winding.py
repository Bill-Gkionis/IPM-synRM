"""optimization.winding: the winding study. Choose the turns, the wire, the strands in hand and the parallel paths of a fixed geometry by ENUMERATING every
combination, keep the exact Pareto front, and grade NSGA-II and NSGA-III against it.

What it is, and why
    The geometry is fixed (a dict of motorgeo.Params overrides: {} is the reference); what is left to choose is the winding: N turns a coil side, one of 19 wire
    sizes, n_s strands in hand, a parallel paths. That is 60 x 19 x 30 x 3 = 102,600 combinations, and one numpy sweep evaluates all of them in 13 ms, so the
    exact Pareto front of the three objectives is not estimated, it is read off. That is the first reason for the study. The second: this is the one problem in
    the project where an evolutionary optimizer can be GRADED against the truth before it is trusted on the geometry, where nobody knows the truth. The same
    problem, the same objective function, NSGA-II and NSGA-III with a tenth of the enumeration's evaluations, 5 seeds each, then the hypervolume ratio, the IGD and
    the share of the exact front's points found exactly.
    The FEMM part (the rated point, the field of the copper, the turns check) is optimization/winding_femm.py, one session of the solver slot of about 25 s; this
    file reads its numbers (or its saved femm.json) and is Python from there on. The pictures are optimization/winding_figures.py.

The variables
    N     turns a coil side, 1 ... 60 (Params.turns is the placeholder 10).
    d     the bare wire's diameter: the 19 IEC 60317 preferred (R 20) sizes from 0.200 to 1.600 mm, each with its grade 2 maximum overall diameter. Published for 14
          of them: Elektrisola's "Technical data for enamelled copper wire based on IEC 60317" (sheet MW.IECJIS.eng/2015.2, 0.200 ... 0.500 mm) and the Multicomp
          Pro "Enamelled Copper Wire" sheet (IEC 60317-0-1:2014, Grade 2, V1.2 of 17/10/23; 0.560, 0.710, 0.800, 1.000, 1.250 mm); where they overlap they agree to the
          last digit (URLs in the constants). 0.630, 0.900, 1.120, 1.400 and 1.600 mm are in no table I could read: their overall diameter is INTERPOLATED (the increase
          over the nominal diameter, linearly between its published neighbours), and marked * in the report and wire_catalogue.csv. It only decides the slot-opening
          limit, which falls between 1.400 mm (1.503 interpolated) and 1.600 mm (1.709), so the interpolation cannot change a result.
    n_s   strands in hand, 1 ... 30: a turn is n_s round wires wound together, taken as transposed (ideal: they share the current equally, no circulating current).
    a     parallel paths, and only 1, 2 or 4. motorgeo.winding gives 12 coils a phase (double layer: one coil a slot, 36 / 3), and their emf phasors, s (exp(j p alpha_top) -
          exp(j p alpha_bottom)), are 3 distinct ones of equal size (at -30, -10 and +10 degrees) each repeated 4 times, because the winding repeats every 90 degrees
          (gcd(36 slots, 4 pole pairs) = 4: the sector the FEMM model draws). Paths in parallel must have equal emf or the difference drives a current round the loop.
          allowed_paths() counts, by brute force, the ways to split the 12 coils into a equal groups with equal emf sums: 1, 108, 0, 576, 0, 0 for a = 1, 2, 3, 4, 6, 12.
          So a is a divisor of 4. (a changes the winding's series turns and current but not its slot: every path still has N turns a coil.)
    Constraints (G <= 0 allowed): the copper fill k = N n_s (pi d^2 / 4) / A_layer <= 0.45 (A_layer = motorgeo.check's coil_layer_area_mm2, one layer: 25.737 mm^2;
    random-wound round wire); the strand's overall diameter plus 2 x 0.2 mm (room each side for the needle that winds it) <= the 2.0 mm slot opening, which removes
    the 1.600 mm wire; and R i < V_max (a winding whose resistive drop alone is more than the voltage limit cannot drive the rated current at all: no corner speed).

The model (Model.evaluate: a winding in, the three objectives out; the maths, with its assumptions)
    Rated point, from FEMM (winding_femm.py): 100 ampere-turns a coil side at the geometry's own MTPA angle (gamma = 22.457 degrees, i_d = -3.8200 A, i_q = 9.2416 A
    for the reference winding of 10 turns and one path, 120 series turns): the mean torque 1.12880 N.m is THE rated torque of the whole study, psi_d = 12.1803 mWb,
    psi_q = 19.7770 mWb (means over 15 rotor positions).
    series turns   n = c N / a, c = 12: each of the a paths has c / a coils of N turns, and a path's voltage is the phase's.
    current        i = AT a / N  [A peak], AT = 100: the coil side's ampere-turns, N turns each carrying i / a, are the 100 the map was solved at.
    flux linkages  psi' = psi n / 120 and i' = i_ref i / 10 (i_ref = the 10 A of the map): psi' i' = psi i, so the torque 1.5 p (psi_d i_q - psi_q i_d), 1.12868 N.m
                   from the flux linkages (1.12880 from the air-gap element), is the same for every one of the 102,600 combinations (checked). The scaling is right only if
                   the turns do not change the field: winding_femm solves 5, 10 and 20 turns at 20, 10 and 5 A and gets the same numbers, bit for bit.
    resistance     R = rho l_mt n / (a n_s pi d^2 / 4)  =  rho l_mt c N / (a^2 n_s A_s): n series turns in each path, a n_s strand-areas of copper per turn.
                   rho = rho_20 (1 + alpha (T - 20)), rho_20 = 1.7241e-8 Ohm m (the IEC standard for annealed copper), alpha = 0.00393 /K, T = 100 deg C: 2.2662e-8 Ohm m.
    mean turn length  l_mt = 2 L + pi span + 4 e. The 2-D model has no end windings, so this is a stand-in, written down: one end turn is a half circle on the coil's span
                   (the chord of the 4-slot pitch at the radius of the slot's two layers' centroid, 52.45 mm: span 35.88 mm, half circle 56.4 mm) plus a straight leg of
                   e = 2 mm at each of its two ends; two ends and two active sides make a turn: 2 x 20 + 2 x 60.4 = 160.7 mm for the 20 mm placeholder stack. R and the DC
                   loss are proportional to it, the AC loss does not use it (over the active length only). The corner speed depends on it through R only.
    DC loss        1.5 R i^2 (amplitude-invariant: three phases, 3 R (i / sqrt 2)^2). With k = N n_s A_s / A_layer: R i^2 = rho l_mt c AT^2 / (N n_s A_s) = rho l_mt c AT^2 /
                   (k A_layer): the turns and the paths cancel (R grows as (N / a)^2, i^2 falls as (a / N)^2) and so does the wire: the DC loss depends on the FILL only,
                   1.5 rho l_mt c AT^2 / (k A_layer) = 56.604 W at k = 0.45 (509 W at the 5 % fill of the 10-turn placeholder). Checked over all 102,600 combinations
                   (agreement 6.7e-16), and in a table of eight windings of one fill (173.9 W each) whose currents run 10 ... 40 A and corner speeds 1,796 ... 10,407 rpm.
    AC loss        a round strand of diameter d in a uniform transverse field B(t): the changing field drives an axial electric field E_z = -y dB/dt (zero net current), so
                   the loss per unit length is sigma integral of E_z^2 over the disc = sigma (dB/dt)^2 pi d^4 / 64, i.e. P' = pi sigma d^4 <(dB/dt)^2> / 64, which for a
                   sinusoid of peak B-hat is the textbook pi sigma omega^2 B-hat^2 d^4 / 128. winding_femm measures G = the mean over the copper area of the mean square of
                   dB/d(theta_e) = 6.6824e-4 T^2/rad^2 (48 rotor positions of one electrical period; the sinusoidal reading, the fundamental alone, is 5.4126e-4, 81 %);
                   with theta_e = omega t, <(dB/dt)^2> = omega^2 G, omega = 2 pi p rpm / 60 = 1256.6 rad/s at the 3000 rpm of the loss objective. Over every strand of every
                   coil side, 72 N n_s of them for L = 20 mm: P_ac = 72 N n_s L pi sigma d^4 omega^2 G / 64 = k A_layer n_sides L sigma d^2 omega^2 G / 16: for a given fill it
                   grows with d^2 only. Assumptions: the field is the 2-D slot field of the homogenised coil (the local fields of neighbouring wires are not modelled), the
                   strands are spread evenly over the layer, and a wire's own current does not matter (skin effect: d / delta = 0.26 at 200 Hz).
    validity       delta = sqrt(rho / (pi f mu_0)) = 5.36 mm at 200 Hz and 1.09 mm at 4.8 kHz. The low-frequency formula above is the first term of the exact result for a
                   round wire in a transverse field (Bessel functions, proximity_ratio()): the exact loss over it is 0.993 at d / delta = 1, 0.965 at 1.5, 0.898 at 2, so it
                   is within 1 % up to d / delta = 1 and within 5 % up to about 1.65. For the largest allowed wire (1.4 mm) d / delta = 0.26 at the fundamental and 1.01
                   at the 15th harmonic, the highest that carries 1 % of G; the harmonic-weighted ratio is 0.9993: the formula holds for every allowed wire at 3000 rpm.
                   It would stop holding (5 %) for 1.4 mm wire above about 8.0 kHz of field frequency (the 40th harmonic at 3000 rpm; the fundamental itself only at 120,000 rpm).
    corner speed   at the rated point the terminal voltage of a winding at electrical speed w is |v|^2 = (R i_d - w psi_q)^2 + (R i_q + w psi_d)^2 (steady state, amplitude-
                   invariant d-q); it reaches V_max = V_dc / sqrt 3 = 27.713 V (space-vector modulation, 48 V, the bus's worst case) at the root of
                   (psi_d^2 + psi_q^2) w^2 + 2 R (i_q psi_d - i_d psi_q) w + R^2 (i_d^2 + i_q^2) - V_max^2 = 0, w = (-B + sqrt(B^2 - 4 A C)) / 2 A; rpm = 60 w / (2 pi p).
                   R is the DC resistance (the AC loss is at most 0.15 % of the DC loss); the end-winding leakage is not in the 2-D flux linkages, so the corner speeds are optimistic.
    objectives     all minimised: f1 = P_dc + P_ac [W] at the rated point and 3000 rpm; f2 = i [A peak] (the inverter's size); f3 = -corner speed [rpm]. The turns trade
                   f2 against f3 (psi and i scale oppositely, so the corner speed goes as a / N, as the current does); at a given fill the wire and strands trade the DC loss
                   (the fill the catalogue allows) against the AC loss (the strand size) - and the second turns out to be 0.002 ... 0.03 W of 57 on the front (below).

The front, and why its objectives are compared after round_sig
    The same winding written twice, (N, a, n_s) and (2N, 2a, n_s / 2), is bit-identical in all three objectives (the integer products are formed before the
    divisions: checked). Windings that agree only in exact arithmetic (an equal fill reached with another wire) differ by about 1e-16, and compared as they come that noise would let
    one beat the other in one objective and put both on the front; so every comparison of objective vectors is made after rounding to 10 significant digits.
    pareto_mask() sorts lexicographically and sweeps (a design that beats another comes before it), and verify_front() proves the result against ALL the allowed
    designs: no front member is beaten, and every design off the front is beaten by a front member. pymoo's hypervolume of the front is checked against an independent
    slicing computation (hv3_slicing).

NSGA-II and NSGA-III, and how they are graded
    Same problem, same constraints, same objective function (Model.evaluate, so the objectives are bit-identical to the enumeration's), population 100, pymoo 0.6.2's
    NSGA2 and NSGA3, the same variation operators (pymoo's defaults for each type of variable), duplicate elimination, seeds 1 ... 5. n_gen counts the first, random
    generation: 50 generations are 5,000 evaluations (4.9 % of the 102,600 combinations), 100 are 10,000 (9.7 %); a run to 100 passes through 50 exactly as a run stopped there
    (the grades are read at both). NSGA-III uses 91 Das-Dennis reference directions (12 divisions of 3 objectives; the population is at least that many).
    Encoding, 'mixed' (the one reported first): pymoo's own mixed-variable support. N and n_s are Integer variables: SBX (eta 15, applied to 90 % of the matings, each gene
    crossed with probability 0.5) and polynomial mutation (eta 20, applied to 90 % of the offspring, then each gene with pymoo's min(0.5, 1 / n) = 0.5 for a group of two),
    the results rounded. The wire and the number of paths are Choice variables, categorical, without order: uniform crossover (each gene swaps with probability 0.5) and
    mutation by drawing a new value from the set (each gene with probability 0.5). 'index': all four are integers in one array (N, the wire's list position, n_s, the position
    of a in (1, 2, 4)), SBX and polynomial mutation as above with rounding, each gene mutated with probability 0.25, so the list positions are treated as ordered numbers.
    Both are run (--encoding).
    Grades, against the exact front's 120 distinct points: HYPERVOLUME RATIO = HV of the set over HV of the exact front, the objectives scaled 0 ... 1 between the front's
    ideal and nadir (one scaling for every run), one common reference point (1.1, 1.1, 1.1), pymoo's HV (moocore underneath); IGD = pymoo's, the mean distance from each
    exact point to the nearest point of the set, in the same scaled space; EXACT POINTS FOUND = the front's objective vectors present, to 10 digits, in the set (any winding with
    that vector counts). Two sets are graded: the final POPULATION's non-dominated feasible designs (what pymoo returns; at most 100, and the front has 120 points, so at most
    83 % can be found) and the ARCHIVE (the non-dominated designs among everything the run evaluated: what the geometry study's CSV keeps). describe_misses() says what the
    misses are: points the run evaluated but its population dropped, points whose N / a it tried with another wire or strand count (and how much loss that costs), and
    points whose N / a it never tried. Mean, sample standard deviation (n - 1), minimum and maximum over the seeds.

The example winding
    pick_example(): the front member with the smallest peak current among those whose corner speed is at least 10 % above the 3000 rpm reference (3300 rpm: the rated torque held
    to the reference speed with 10 % of the voltage to spare on the 48 V bus); windings that are one front point written several ways are settled by the fewest strands, then the
    fewest paths. It is an EXAMPLE: the stack (20 mm) is a placeholder and the turns are sized again once it is ruled. pick_for_current() gives the front member for an inverter of
    a stated peak current, which is the other way to choose.

Run (from the Python to FEMM folder; the solver slot is taken for about 25 s, and waited for if another run holds it)
    uv run python -m optimization.winding --label reference                       (or .venv/Scripts/python.exe -m optimization.winding ...)
    --overrides '{"mag_y": 5.5}' --label mag_y_5p5     another geometry: FEMM stage and all (the output folder carries the label)
    --reuse-femm                    read femm.json of the label instead of running FEMM (no solver slot); --encoding mixed|index|both; --seeds 5 --pop 100 --gens 100
    --budgets 50,100  --vdc 48  --temp 100  --rpm 3000  --no-figures
Writes out\\optimization\\winding\\<label>\\: report.txt (everything printed), results.json (every number), femm.json and field.npz (the FEMM stage), evaluations.csv (the rated point,
design.Evaluator's own), front_exact.csv, wire_catalogue.csv, enumeration_all.npz (all 102,600 combinations), nsga_runs.csv, nsga_history.csv (every generation of every run),
nsga_summary.csv, and four figures (winding_figures.py). Exit status 0 if all 13 checks pass. The seeds are fixed: runs give the same numbers (a 20-seed run of the same study, made once in a scratch label, gives NSGA-II 0.964 (0.035) and
37.8 (8.4) exact points, NSGA-III 0.917 (0.073) and 41.0 (6.4) at 100 generations: the same picture with a better estimate of the spread).

Numbers to expect (this PC, 2026-11-08 / 09, FEMM 5.0, pymoo 0.6.2; a whole run takes a minute or so with a free slot: 26 s of FEMM, then 25 ... 40 s of enumeration, optimizers and pictures,
depending on what else the PC is doing; the numbers below were the same in the runs of both days)
    enumeration  102,600 combinations in 13 ... 14 ms; 14,241 allowed (13.9 %): the fill removes 86,709, the slot opening 5,400, the standstill drop 1,620 (a combination can break
                 several); the front sweep 0.13 ... 0.15 s: 136 designs = 120 distinct objective vectors (16 copies of windings written twice); verified exact. One winding for each distinct
                 N / a (120 values from 0.25 to 60: current 1.67 ... 400 A, corner speed 92 ... 113,621 rpm), loss 56.62 ... 62.49 W, fill 0.408 ... 0.450, wires 0.200 ... 0.450 and
                 0.800 mm. The loss varies along the front only because the best fill the catalogue allows changes with N; the AC loss is 0.002 ... 0.031 W on the front (0.00001 ...
                 0.089 W over all the allowed designs), at most 0.055 % of the loss of a front design (0.145 % of any allowed design's).
    optimizers   population grades after 100 generations (10,000 evaluations; about 8,300 distinct designs), mean (standard deviation) over 5 seeds, mixed encoding:
                 NSGA-II   HV ratio 0.949 (0.032), IGD 0.0124, exact points found 43.8 (8.7) of 120 = 36.5 %   (archive: 50.8 = 42.3 %; after 50 generations: 0.931, 20.6 = 17.2 %)
                 NSGA-III  HV ratio 0.896 (0.109), IGD 0.0145, exact points found 41.0 (10.4) of 120 = 34.2 %  (archive: 51.0 = 42.5 %; after 50 generations: 0.840, 23.6 = 19.7 %)
                 the index encoding gives HV ratios 0.941 / 0.929 and 42.8 / 35.2 exact points (archives 58.6 / 53.0). A run takes 1 ... 3 s. The optimizers find the trade-off (95 % of
                 the hypervolume, 90 % for NSGA-III, whose spread over seeds is three times larger: one seed stops at 0.74) but reach the exact optimum of about a third of the 120 N / a values.
                 The misses are not at one end of the front: NSGA-II runs miss a point in 57 ... 69 % of the cases in each of five ranges of N / a and 60 ... 65 % in each range of
                 the exact loss (NSGA-III 60 ... 80 % and 55 ... 70 %), whichever range. The windings the runs tried at the N / a of a missed point are, at the median, 2 % worse in
                 loss; the worst (110 % for NSGA-II, 300 % for NSGA-III) are values of N / a such as 26.5, 12.75 or 10.75 that only one pair (N, a) makes. Each exact point is found
                 by one run or another (104 of 120 by at least one of the 5 NSGA-II runs), but only 2 by all 5.
    example      N = 31, a = 4, 6 strands of 0.280 mm, fill 0.4450, 12.90 A peak, R = 0.2292 Ohm, corner speed 3350 rpm, loss 57.245 W (DC 57.241 + AC 0.0038).

Self-check (the answers follow).
  1. Why is the exact front one winding for each N / a, and not a curve with a knee? What would make it several windings at one N / a?
  2. Two windings have the same fill and different turns. Why do they lose the same DC power, and what do the turns change?
  3. An optimizer returns a set with HV ratio 0.95 that holds 52 of the 120 exact points. Is it good?
Answers.
  1. f2 = 100 a / N and f3 = -(corner speed), which goes as a / N too, depend on N / a first of all (the resistance lowers f3 more at high N / a, where R i is a large part
     of V_max): two windings of different N / a are better in one of them and worse in the other, so none beats another; at one N / a only the lowest loss survives, because a
     higher fill lowers R and the DC loss together (the corner speed rises as R falls). Several windings per N / a would appear if the AC loss were large enough to pay for a lower fill with a finer wire: a thin-
     wire design with a slightly higher R would then be better in f1 and worse in f3. Here the AC loss is 0.1 % of the DC loss, so the only ties are the same winding written twice.
  2. R grows as (N / a)^2 and i^2 falls as (a / N)^2, so 1.5 R i^2 is 1.5 rho l_mt c AT^2 / (k A_layer): the copper that carries 100 ampere-turns does the same work in few
     thick turns or many thin ones. The turns set the winding's voltage: psi goes with N / a, so the current is 100 a / N and the corner speed goes as a / N.
  3. For the designer, mostly: 95 % of the volume the front dominates means its designs are on or near the trade-off, and most of the 68 it misses are near misses (a few per cent
     of the loss at the same current, with a few much worse). For an exact claim ('the optimum winding for this current') it is not: it found the exact optimum for fewer than half of the currents.
     Grade an optimizer by the indicators that mean something for the use, and use the exact hits only where the truth is known, as here.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

import motorgeo as mg
import paths
from optimization import design

# ---------------------------------------------------------------------------------------------- the fixed choices (each one is a ruling or an assumption)
MU0 = 4e-7 * math.pi                    # H/m
V_DC = 48.0                             # V: the worst case of the 48 ... 60 V bus
RPM_REF = 3000.0                        # the speed at which the loss objective is evaluated [rpm mechanical]
T_COPPER = 100.0                        # deg C: the winding's temperature in the resistance and in the conductivity
RHO_20 = 1.7241e-8                      # Ohm m: annealed copper at 20 deg C (the IEC standard value, 0.017241 Ohm mm^2/m)
ALPHA_20 = 0.00393                      # per K: copper's temperature coefficient of resistance at 20 deg C (IEC 60228)
FILL_MAX = 0.45                         # copper fill of a coil layer, bare copper area over layer area (random-wound round wire)
CLEARANCE_MM = 0.2                      # room left on EACH side of the strand's overall diameter in the slot opening [mm]
END_EXTENSION_MM = 2.0                  # straight conductor beyond the stack before an end turn bends, each of the 4 legs of a turn [mm]
N_RANGE = (1, 60)                       # turns per coil side
NS_RANGE = (1, 30)                      # strands in hand
VOLTAGE_MARGIN = 1.0                    # the standstill drop R i must stay below this share of the voltage limit

# IEC 60317 enamelled round copper wire, preferred (R 20) nominal conductor diameters from 0.200 to 1.600 mm, with the Grade 2 maximum overall diameter.
# Published values (nominal d -> maximum overall diameter of Grade 2, mm), and where each comes from:
#   E  Elektrisola, "Technical data for enamelled copper wire based on IEC 60317", sheet MW.IECJIS.eng/2015.2, the IEC half of the table (columns
#      Grade 2 min. / max. of "Outer Diameter"): https://www.elektrisola.com/Attachments/TechnicalDataBySize/ELEKTRISOLA_EnCuWire_IECJIS_Datasheet_eng.pdf
#      (read 2026-11-08; the sheet ends at 0.500 mm).
#   M  Multicomp Pro "Enamelled Copper Wire" data sheet V1.2 of 17/10/23, Standard IEC 60317-0-1:2014, Grade 2 (2L), column "Max. overall Diameter":
#      https://www.farnell.com/datasheets/3974480.pdf (the same numbers as Farnell's older sheets 1876784.pdf and 1942907.pdf, which say BS EN 60317-20).
#      Where the two sources overlap (0.200 ... 0.500 mm) they agree to the last digit.
PUBLISHED_G2_MAX = {0.200: (0.239, "E"), 0.224: (0.266, "E"), 0.250: (0.297, "E"), 0.280: (0.329, "E"), 0.315: (0.367, "E"), 0.355: (0.411, "E"),
                    0.400: (0.459, "E"), 0.450: (0.513, "E"), 0.500: (0.566, "E"), 0.560: (0.630, "M"), 0.710: (0.789, "M"), 0.800: (0.884, "M"),
                    1.000: (1.094, "M"), 1.250: (1.349, "M"), 1.500: (1.606, "M")}      # 1.500 (an R 40 size) is used only as an anchor for 1.400 and 1.600
R20_MM = (0.200, 0.224, 0.250, 0.280, 0.315, 0.355, 0.400, 0.450, 0.500, 0.560, 0.630, 0.710, 0.800, 0.900, 1.000, 1.120, 1.250, 1.400, 1.600)
PATHS_ALLOWED = (1, 2, 4)               # parallel paths: derived from motorgeo.winding by allowed_paths() below, written here as the answer it must give


def catalogue() -> dict:
    """The wire catalogue: {"d_mm", "d_over_mm", "published" (bool), "source" (E, M or 'interpolated')} for the 19 R 20 sizes.

    The five sizes no table I could read lists (0.630, 0.900, 1.120, 1.400, 1.600 mm) get an INTERPOLATED overall diameter: the increase of the Grade 2
    maximum overall diameter over the nominal one (published: 0.039 mm at 0.200 ... 0.099 at 1.250, 0.106 at 1.500) is interpolated linearly in d between its
    published neighbours, and at 1.600 mm extended with the slope of the last two (1.250 to 1.500). The only thing the overall diameter decides in this study is
    the slot-opening limit, and that limit falls between 1.400 mm (interpolated 1.503) and 1.600 mm (1.709): a few hundredths of a millimetre of error in either
    changes nothing. FLAGGED because it is an approximation, not a table value."""
    anchors = sorted(PUBLISHED_G2_MAX)
    inc = np.array([PUBLISHED_G2_MAX[d][0] - d for d in anchors])
    d_over, published, source = [], [], []
    for d in R20_MM:
        if d in PUBLISHED_G2_MAX:
            d_over.append(PUBLISHED_G2_MAX[d][0])
            published.append(True)
            source.append(PUBLISHED_G2_MAX[d][1])
        else:
            if d > anchors[-1]:
                slope = (inc[-1] - inc[-2]) / (anchors[-1] - anchors[-2])
                extra = inc[-1] + slope * (d - anchors[-1])
            else:
                extra = float(np.interp(d, anchors, inc))
            d_over.append(round(d + extra, 4))
            published.append(False)
            source.append("interpolated")
    return {"d_mm": np.array(R20_MM), "d_over_mm": np.array(d_over), "published": np.array(published), "source": source}


# ---------------------------------------------------------------------------------------------- the winding's symmetry: how many parallel paths are allowed
def phase_coils(p: mg.Params, phase: str = "A") -> list[dict]:
    """The coils of one phase from motorgeo.winding: {"top" (slot), "bottom" (slot), "sign", "emf" (complex)}.
    motorgeo gives, for every slot, the (phase, sign) of its top layer and of its bottom layer; the bottom layer of slot k + coil_pitch carries the return
    side of the coil whose top side is in slot k (the same phase, the opposite sign), so a coil is (k, k + pitch). A conductor of sign s at the mechanical
    angle alpha has an electromotive force proportional to s exp(j pole_pairs alpha); the coil's is the sum over its two sides."""
    table = mg.winding(p)
    out = []
    for k in range(p.slots):
        (ph, sign), _ = table[k]
        if ph != phase:
            continue
        k2 = (k + p.coil_pitch) % p.slots
        ph2, sign2 = table[k2][1]
        if ph2 != ph or sign2 != -sign:
            raise RuntimeError(f"slot {k} and slot {k2} do not make a coil: {table[k]} {table[k2]}")
        a1, a2 = (math.radians(mg.slot_angle(p, kk)) for kk in (k, k2))
        emf = sign * (np.exp(1j * p.pole_pairs * a1) - np.exp(1j * p.pole_pairs * a2))
        out.append({"top": k, "bottom": k2, "sign": sign, "emf": complex(emf)})
    return out


def equal_emf_partitions(emfs: np.ndarray, groups: int, tol: float = 1e-9) -> int:
    """How many ways the coils (their emf phasors) split into `groups` groups of equal size whose emf sums are all equal. Brute force: the groups of one
    phase are 12 coils at most, so a few thousand partitions."""
    n = len(emfs)
    if n % groups:
        return 0
    size = n // groups
    target = emfs.sum() / groups
    count = 0

    def split(remaining: list[int], left: int) -> None:
        nonlocal count
        if left == 1:
            count += int(abs(emfs[remaining].sum() - target) < tol)
            return
        first, rest = remaining[0], remaining[1:]
        for combo in itertools.combinations(rest, size - 1):        # the group that holds the first coil: fixes the order of the groups
            group = (first, *combo)
            if abs(emfs[list(group)].sum() - target) < tol:
                split([c for c in remaining if c not in group], left - 1)

    split(list(range(n)), groups)
    return count


def allowed_paths(p: mg.Params) -> dict:
    """The numbers of parallel paths a symmetric winding allows, derived from motorgeo.winding.
    Paths in parallel must have the same emf (the same sum of their coils' phasors), or the difference drives a circulating current. So `a` is allowed when
    the phase's coils split into `a` equal groups with equal emf sums. Result: {"allowed", "counts" ({a: number of such splits}), "coils", "distinct_emfs",
    "periodic_sectors", "emf_angles_deg"}. The reason is the winding's periodicity: it repeats every 360 / gcd(slots, pole pairs) degrees (4 times here: the
    90 degree sector the FEMM model draws), so the coils fall in 4 copies of the same 3 emf phasors, and only a that divides 4 can give equal sums."""
    coils = phase_coils(p)
    emf = np.array([c["emf"] for c in coils])
    counts = {a: equal_emf_partitions(emf, a) for a in range(1, len(coils) + 1) if len(coils) % a == 0}
    distinct = len({(round(abs(e), 6), round(math.degrees(np.angle(e)), 4)) for e in emf})
    return {"allowed": sorted(a for a, c in counts.items() if c > 0), "counts": counts, "coils": len(coils), "distinct_emfs": distinct,
            "periodic_sectors": p.n_sectors, "emf_angles_deg": sorted({round(math.degrees(np.angle(e)), 1) for e in emf})}


# ---------------------------------------------------------------------------------------------- copper at the winding temperature, and the skin effect
def resistivity(t_c: float = T_COPPER) -> float:
    """Copper's resistivity at t_c [Ohm m]: rho_20 (1 + alpha_20 (t_c - 20)); 2.266e-8 at 100 deg C."""
    return RHO_20 * (1.0 + ALPHA_20 * (t_c - 20.0))


def skin_depth(f_hz: float, rho: float) -> float:
    """sqrt(rho / (pi f mu_0)) [m]: where the current of a wire in an alternating field has fallen to 1/e."""
    return math.sqrt(rho / (math.pi * f_hz * MU0))


def proximity_ratio(d_m: float, f_hz: float, rho: float) -> float:
    """The exact proximity loss of a round wire in a uniform transverse field, divided by the low-frequency formula (pi sigma omega^2 B^2 d^4 / 128).
    Inside the wire (radius a) the vector potential is A = C I_1(g r) cos(phi), g = (1 + j) / delta; outside, A = (-B r + D / r) cos(phi). Matching A and
    dA/dr at r = a gives C = -2 B / (g I_0(g a)), and the power per unit length is (pi sigma omega^2 / 2) |C|^2 integral of |I_1(g r)|^2 r dr from 0 to a.
    Divided by pi sigma omega^2 B^2 a^4 / 8 (the low-frequency value: I_1(x) = x/2) it is 16 integral / (a^4 |g|^2 |I_0(g a)|^2): 1 at low frequency,
    0.993 at d / delta = 1, 0.965 at 1.5 and 0.898 at 2: the low-frequency formula OVERestimates a little where it stops being valid (the wire shields itself). The wire is
    alone in a uniform external field: the neighbouring wires' own eddy fields are not modelled (the textbook's 'external proximity' loss)."""
    from scipy.special import iv
    a = d_m / 2.0
    g = (1.0 + 1.0j) / skin_depth(f_hz, rho)
    x, w = np.polynomial.legendre.leggauss(96)
    r, wr = 0.5 * a * (x + 1.0), 0.5 * a * w
    integral = float(np.sum(np.abs(iv(1, g * r)) ** 2 * r * wr))
    return 16.0 * integral / (a ** 4 * abs(g) ** 2 * abs(iv(0, g * a)) ** 2)


# ---------------------------------------------------------------------------------------------- the model: one winding in, its three objectives out
class Model:
    """The winding model of one geometry. Built from the geometry's Params, its rated point and its field statistics (the FEMM stage), then evaluate()
    turns arrays of (N, d index, n_s, a) into the objectives and the constraints, all at once, in numpy."""

    def __init__(self, p: mg.Params, femm: dict, *, v_dc: float = V_DC, t_cu: float = T_COPPER, rpm_ref: float = RPM_REF, g_key: str = "48"):
        self.p, self.v_dc, self.t_cu, self.rpm_ref = p, v_dc, t_cu, rpm_ref
        r = femm["rated"]
        self.rated = r
        self.turns_ref, self.amps_ref, self.at = int(r["turns"]), float(r["amps"]), float(r["ampere_turns_per_coil_side"])
        self.psi_d_ref, self.psi_q_ref, self.i_d_ref, self.i_q_ref = (float(r[k]) for k in ("psi_d", "psi_q", "i_d", "i_q"))
        self.torque = float(r["T_mean"])
        self.a_layer = float(r["coil_layer_area_mm2"]) * 1e-6                 # m^2
        self.slot_open = p.slot_open                                          # mm
        self.stack = p.stack * 1e-3                                           # m
        self.coils_per_phase = p.slots // 3                                   # a double-layer winding has one coil a slot: 12 a phase
        self.n_sides = 2 * p.slots                                            # coil sides in the machine: two layers in every slot
        self.series_ref = self.coils_per_phase * self.turns_ref               # the reference's series turns a phase (one path): 120
        self.v_max = v_dc / math.sqrt(3.0)                                    # peak phase voltage at the edge of the linear range of space-vector modulation
        self.rho = resistivity(t_cu)
        self.sigma = 1.0 / self.rho
        self.f_ref = p.pole_pairs * rpm_ref / 60.0                            # electrical frequency at the reference speed [Hz]
        self.w_ref = 2.0 * math.pi * self.f_ref                               # rad/s
        stats = femm["statistics"][g_key]
        self.G, self.G_k = float(stats["G"]), np.array(stats["G_k"])          # T^2 per rad^2 of electrical angle
        self.cat = catalogue()
        self.d = self.cat["d_mm"] * 1e-3                                      # m
        self.d_over = self.cat["d_over_mm"]
        self.paths = np.array(PATHS_ALLOWED)
        g = end_turn_geometry(p)
        self.r_mean, self.span, self.l_mt = g["r_mean_mm"] * 1e-3, g["span_mm"] * 1e-3, g["l_mt_mm"] * 1e-3
        self.geometry = g

    # ---- the pieces of the maths, written once
    def evaluate(self, N, di, ns, a) -> dict:
        """The winding's numbers. N turns per coil side, di the index of the wire in the catalogue, ns strands in hand, a parallel paths (arrays, one element
        a design, or numbers). The order of the operations is fixed so that two designs that are the same winding written twice, (N, a, n_s) and (2N, 2a, n_s/2),
        come out bit for bit equal: products of integers first, then the divisions.
        series turns   12 N / a           the coils of a phase in one path, times their turns
        current        i = 100 a / N      A peak: the ampere-turns of a coil side stay 100 (N i / a), which is what the map was solved at
        flux linkages  psi' = psi (12 N / a) / 120   and   i' = i_ref (i / 10): their product, and so the torque, does not change
        resistance     R = rho L_mt (12 N / a) / (a n_s pi d^2 / 4)
        fill           k = N n_s (pi d^2 / 4) / A_layer
        DC loss        1.5 R i^2 (three phases, peak current i: 3 R (i / sqrt 2)^2)
        AC loss        72 N n_s L (pi sigma d^4 omega^2 G / 64): every strand of every coil side, over the active length only
        corner speed   the electrical w that solves (R i_d - w psi_q)^2 + (R i_q + w psi_d)^2 = V^2, V = V_dc / sqrt 3
        Returns arrays: F (n, 3) = [DC + AC loss W, current A, minus the corner speed rpm], G (n, 3) <= 0 when allowed (fill, slot opening, standstill drop),
        and the pieces: k, R, i, p_dc, p_ac, rpm, d_over, series."""
        N, ns, a = (np.asarray(v, float) for v in (N, ns, a))
        di = np.asarray(di, int)
        N, ns, a, di = np.broadcast_arrays(N, ns, a, di)
        d = self.d[di]
        a_strand = np.pi * d ** 2 / 4.0
        series = self.coils_per_phase * N / a
        i = self.at * a / N
        R = self.rho * self.l_mt * series / ((a * ns) * a_strand)
        k = (N * ns) * a_strand / self.a_layer
        p_dc = 1.5 * R * i ** 2
        p_ac = self.n_sides * (N * ns) * self.stack * (np.pi * self.sigma * d ** 4 * self.w_ref ** 2 * self.G / 64.0)
        s_psi, s_i = series / self.series_ref, i / self.amps_ref
        i_d, i_q = self.i_d_ref * s_i, self.i_q_ref * s_i
        psi_d, psi_q = self.psi_d_ref * s_psi, self.psi_q_ref * s_psi
        qa = psi_d ** 2 + psi_q ** 2
        qb = 2.0 * R * (i_q * psi_d - i_d * psi_q)
        qc = R ** 2 * (i_d ** 2 + i_q ** 2) - self.v_max ** 2
        ok = qc < 0.0
        w_c = np.where(ok, (-qb + np.sqrt(np.where(ok, qb ** 2 - 4.0 * qa * qc, 0.0))) / (2.0 * qa), 0.0)
        rpm = w_c * 60.0 / (2.0 * math.pi * self.p.pole_pairs)
        F = np.stack([p_dc + p_ac, i, -rpm], axis=-1)
        G = np.stack([k - FILL_MAX, (self.d_over[di] + 2.0 * CLEARANCE_MM - self.slot_open) / self.slot_open, (R * i - VOLTAGE_MARGIN * self.v_max) / self.v_max], axis=-1)
        return {"F": F, "G": G, "k": k, "R": R, "i": i, "p_dc": p_dc, "p_ac": p_ac, "rpm": rpm, "d_over": self.d_over[di], "series": series,
                "psi_d": psi_d, "psi_q": psi_q, "i_d": i_d, "i_q": i_q, "w_c": w_c}

    def torque_of(self, e: dict) -> np.ndarray:
        """1.5 p (psi_d i_q - psi_q i_d) of the scaled flux linkages and currents: the torque whatever the winding (it must equal the reference's)."""
        return 1.5 * self.p.pole_pairs * (e["psi_d"] * e["i_q"] - e["psi_q"] * e["i_d"])

    def dc_loss_at_fill(self, k) -> np.ndarray:
        """The DC loss as a function of the fill alone: 1.5 rho L_mt (coils a phase) AT^2 / (k A_layer): no N, no a, no d."""
        return 1.5 * self.rho * self.l_mt * self.coils_per_phase * self.at ** 2 / (np.asarray(k, float) * self.a_layer)

    def ac_validity(self) -> list[dict]:
        """For every wire: the skin depth ratio d / delta at the fundamental and at the highest harmonic that carries 1 % of G or more, and the exact-over-
        low-frequency ratio of the proximity loss weighted by the harmonics' shares of G (the factor the objective would carry if it used the exact formula)."""
        share = self.G_k / self.G_k.sum()
        harmonics = [k for k in range(1, len(self.G_k)) if share[k] >= 0.01]
        out = []
        for j, d in enumerate(self.d):
            ratios = {k: proximity_ratio(d, k * self.f_ref, self.rho) for k in range(1, len(self.G_k)) if share[k] > 0}
            factor = sum(share[k] * ratios[k] for k in ratios) / sum(share[k] for k in ratios)
            k_top = max(harmonics)
            out.append({"d_mm": d * 1e3, "d_over_mm": float(self.d_over[j]), "d_over_delta_fund": d / skin_depth(self.f_ref, self.rho), "k_top": k_top,
                        "d_over_delta_top": d / skin_depth(k_top * self.f_ref, self.rho), "ratio_fund": ratios[1], "ratio_top": ratios[k_top],
                        "ratio_weighted": factor, "valid_5pct": bool(min(ratios.values()) >= 0.95)})
        return out


def end_turn_geometry(p: mg.Params) -> dict:
    """The mean turn length of a coil, from the slot the geometry draws (the end-turn model, written down):
    r_mean   the radius of the centroid of a slot's two coil layers [mm] (shapely, from motorgeo.build)
    span     the chord of the coil pitch at that radius, 2 r_mean sin(pi pitch / slots) [mm]: the straight-line distance between the two slots of a coil
    l_end    one end turn is modelled as a half circle on the span: (pi / 2) span, plus the two legs' straight extensions, 2 e
    l_mt     2 L_stack + 2 l_end = 2 L + pi span + 4 e  [mm]   (e = END_EXTENSION_MM)
    Assumption: the 2-D model has no end windings, so this is a stand-in; a longer or shorter end turn scales R and both losses by the same factor
    (l_mt appears once in R, and the AC loss does not use it), and the corner speed through R only."""
    from shapely.ops import unary_union
    regions = mg.build(p)
    layers = [r.polygon() for r in regions if r.name in ("coil_0_0", "coil_0_1")]
    centre = unary_union(layers).centroid
    r_mean = math.hypot(centre.x, centre.y)
    span = 2.0 * r_mean * math.sin(math.pi * p.coil_pitch / p.slots)
    l_end = 0.5 * math.pi * span + 2.0 * END_EXTENSION_MM
    return {"r_mean_mm": r_mean, "span_mm": span, "l_end_mm": l_end, "l_mt_mm": 2.0 * p.stack + 2.0 * l_end, "stack_mm": p.stack}


# ---------------------------------------------------------------------------------------------- the exact front
def pareto_mask(F) -> np.ndarray:
    """A boolean mask of the rows of F (minimised objectives) that no other row beats (at least as good in every objective and better in one). Identical rows
    do not beat each other, so every copy of a non-dominated row is kept. Method: sort the rows lexicographically; a row that beats another comes before it in
    that order, so each row only has to be compared with the front found so far (a few hundred rows), not with all."""
    F = np.asarray(F, float)
    order = np.lexsort(F.T[::-1])
    keep = np.zeros(len(F), bool)
    front = np.empty((len(F), F.shape[1]))
    m = 0
    for idx in order:
        f = F[idx]
        if m and np.any(np.all(front[:m] <= f, axis=1) & np.any(front[:m] < f, axis=1)):
            continue
        keep[idx] = True
        front[m] = f
        m += 1
    return keep


def round_sig(F, digits: int = 10) -> np.ndarray:
    """F rounded to `digits` significant digits. Two windings that are the same winding in other numbers (equal fill reached with a different wire, say) get
    objective values that agree to about 1e-16 but not to the last bit; compared as they come, that rounding noise would let one of them 'beat' the other in
    one objective and put both on the front. Rounded to 10 digits the noise is gone and every difference that matters (a loss that differs in the 6th digit)
    stays. Used for every comparison of objective vectors in this module: the front, the duplicates, the grading."""
    F = np.asarray(F, float)
    with np.errstate(divide="ignore"):
        mag = np.floor(np.log10(np.abs(F)))
    mag = np.where(np.isfinite(mag), mag, 0.0)
    scale = 10.0 ** (digits - 1 - mag)
    return np.round(F * scale) / scale


def distinct_rows(F) -> tuple[np.ndarray, np.ndarray]:
    """(index of the first row of each distinct objective vector, the group number of every row), after round_sig: the same winding written twice, (N, a, n_s)
    and (2N, 2a, n_s / 2), is ONE point of the front."""
    _, first, group = np.unique(round_sig(F), axis=0, return_index=True, return_inverse=True)
    return np.sort(first), group.ravel()


def hv2(points: np.ndarray, ref: np.ndarray) -> float:
    """The area dominated by 2-D points (minimised) up to ref: sweep along the first objective, adding a rectangle each time the second one improves."""
    best, area = ref[1], 0.0
    for x, y in points[np.argsort(points[:, 0], kind="stable")]:
        if y < best:
            area += (ref[0] - x) * (best - y)
            best = y
    return float(area)


def hv3_slicing(F: np.ndarray, ref: np.ndarray) -> float:
    """The volume dominated by 3-D points (minimised) up to ref, by slicing along the third objective (exact; an independent check of pymoo's indicator)."""
    F = np.asarray(F, float)
    F = F[np.all(F < ref, axis=1)]
    if len(F) == 0:
        return 0.0
    F = F[np.argsort(F[:, 2], kind="stable")]
    z = np.append(F[:, 2], ref[2])
    return float(sum((z[i + 1] - z[i]) * hv2(F[: i + 1, :2], ref[:2]) for i in range(len(F)) if z[i + 1] > z[i]))


class Grader:
    """Grades a set of objective vectors against the exact front: hypervolume ratio, IGD and the share of the front's points found exactly.

    Objectives are scaled 0 ... 1 between the exact front's ideal (best of each) and nadir (worst of each), the same scaling for every run, so that loss [W],
    current [A] and speed [rpm] weigh alike; the hypervolume's reference point is 1.1 in every scaled objective (one common point). pymoo's HV and IGD do the
    counting (HV: moocore's algorithm, IGD: the mean distance from every front point to the nearest point of the set). The share found exactly counts the front's
    DISTINCT objective vectors that appear in the set (two windings that are the same winding written twice are one point)."""

    def __init__(self, F_front: np.ndarray, ref_value: float = 1.1):
        from pymoo.indicators.hv import HV
        from pymoo.indicators.igd import IGD
        F_front = np.asarray(F_front, float)
        self.ideal, nadir = F_front.min(axis=0), F_front.max(axis=0)
        self.span = np.where(nadir - self.ideal > 0, nadir - self.ideal, 1.0)
        self.ref = np.full(F_front.shape[1], ref_value)
        first, _ = distinct_rows(F_front)
        self.front = F_front[first]                                      # the distinct points
        self.keys = {self.key(f) for f in self.front}
        self._hv = HV(ref_point=self.ref, norm_ref_point=False)
        self.hv_exact = self.hypervolume(F_front)
        self._igd = IGD(self.scale(self.front))

    def scale(self, F) -> np.ndarray:
        return (np.asarray(F, float) - self.ideal) / self.span

    @staticmethod
    def key(f) -> tuple:
        """A vector as a hashable key: its objectives to 10 significant digits."""
        return tuple(round_sig(np.asarray(f, float)).tolist())

    def hypervolume(self, F) -> float:
        S = self.scale(F)
        S = S[np.all(S < self.ref, axis=1)]
        return float(self._hv(S)) if len(S) else 0.0

    def grade(self, F) -> dict:
        F = np.asarray(F, float).reshape(-1, len(self.ideal))
        if len(F) == 0:
            return {"n": 0, "hv_ratio": 0.0, "igd": float("nan"), "hits": 0, "hit_share": 0.0}
        found = {self.key(f) for f in F}
        hits = len(found & self.keys)
        return {"n": int(len(F)), "hv_ratio": self.hypervolume(F) / self.hv_exact, "igd": float(self._igd(self.scale(F))), "hits": hits,
                "hit_share": hits / len(self.keys)}


# ---------------------------------------------------------------------------------------------- the enumeration
def enumerate_all(model: Model) -> dict:
    """Every combination of N (1 ... 60), wire (19 sizes), n_s (1 ... 30) and a (1, 2, 4): the objectives and constraints of each, in numpy, and the exact front."""
    t0 = time.perf_counter()
    N, di, ns, path_idx = np.meshgrid(np.arange(N_RANGE[0], N_RANGE[1] + 1), np.arange(len(model.d)), np.arange(NS_RANGE[0], NS_RANGE[1] + 1),
                                np.arange(len(model.paths)), indexing="ij")
    N, di, ns, path_idx = N.ravel(), di.ravel(), ns.ravel(), path_idx.ravel()
    e = model.evaluate(N, di, ns, model.paths[path_idx])
    feasible = np.all(e["G"] <= 0.0, axis=1)
    seconds_eval = time.perf_counter() - t0
    idx = np.flatnonzero(feasible)
    t1 = time.perf_counter()
    nd = pareto_mask(round_sig(e["F"][idx]))
    front = idx[nd]
    seconds_front = time.perf_counter() - t1
    return {"N": N, "di": di, "ns": ns, "path_idx": path_idx, "e": e, "feasible": feasible, "front": front, "total": len(N), "n_feasible": int(feasible.sum()),
            "violations": {"fill": int((e["G"][:, 0] > 0).sum()), "slot": int((e["G"][:, 1] > 0).sum()), "voltage": int((e["G"][:, 2] > 0).sum())},
            "seconds_eval": seconds_eval, "seconds_front": seconds_front}


def verify_front(F_all: np.ndarray, front_local: np.ndarray) -> dict:
    """The front is exact when (a) no design of the whole feasible set beats a front member and (b) every design that is not on the front is beaten by a front
    member. (a) and (b) together say front = the set of non-dominated designs. Both are plain comparisons against ALL the designs, not against the sweep's own
    bookkeeping. F_all: the feasible designs' objectives (n, 3); front_local: indices into it. Objectives are compared after round_sig, like the front."""
    F_all = round_sig(F_all)
    Ff = F_all[front_local]
    beaten_front = 0
    for f in Ff:
        beaten_front += int(np.any(np.all(F_all <= f, axis=1) & np.any(F_all < f, axis=1)))
    on = np.zeros(len(F_all), bool)
    on[front_local] = True
    missed = 0
    for j in np.flatnonzero(~on):
        if not np.any(np.all(Ff <= F_all[j], axis=1) & np.any(Ff < F_all[j], axis=1)):
            missed += 1
    return {"front_members_beaten_by_any_design": beaten_front, "other_designs_not_beaten_by_front": missed, "ok": beaten_front == 0 and missed == 0}


# ---------------------------------------------------------------------------------------------- NSGA-II and NSGA-III on the same problem
def design_key(N, di, ns, a) -> tuple:
    """A design as four plain ints: (N, wire index, n_s, a)."""
    return (int(N), int(di), int(ns), int(a))


class ArchiveFront:
    """The non-dominated feasible designs among everything a run has evaluated so far (the run's archive), kept up to date a few designs at a time."""

    def __init__(self):
        self.F = np.empty((0, 3))
        self.X: list[tuple] = []

    def add(self, F_new: np.ndarray, X_new: list[tuple]) -> None:
        if len(F_new) == 0:
            return
        F = np.vstack([self.F, F_new])
        X = self.X + list(X_new)
        keep = pareto_mask(round_sig(F))
        self.F, self.X = F[keep], [x for x, k in zip(X, keep) if k]


def make_problem(model: Model, encoding: str):
    """The study's problem for pymoo, in one of two encodings of the same four variables; every evaluation is logged (problem.log) in order.

    encoding "mixed"  pymoo's own mixed-variable support: N and n_s are Integer(1 ... 60) and Integer(1 ... 30), the wire is a Choice among the 19 catalogue
                      sizes and the number of paths a Choice among (1, 2, 4): categorical variables, which have no order, are crossed by uniform crossover
                      and mutated by drawing a new value from the set, and the integers by SBX and polynomial mutation, rounded
    encoding "index"  all four as integers in a plain array (N, wire index, n_s, index of a in (1, 2, 4)); SBX and polynomial mutation followed by rounding
                      on every gene, so the wire's and the path number's list positions are treated as ordered numbers
    Constraints (G <= 0 allowed): the fill, the slot opening, the standstill voltage drop; objectives: the three of Model.evaluate."""
    from pymoo.core.problem import ElementwiseProblem, Problem
    from pymoo.core.variable import Choice, Integer
    n_d = len(model.d)
    if encoding == "mixed":
        class Mixed(ElementwiseProblem):
            def __init__(self):
                super().__init__(vars={"N": Integer(bounds=N_RANGE), "d": Choice(options=list(range(n_d))), "ns": Integer(bounds=NS_RANGE),
                                       "a": Choice(options=list(PATHS_ALLOWED))}, n_obj=3, n_ieq_constr=3)
                self.log: list[tuple] = []

            def _evaluate(self, X, out, *args, **kwargs):
                key = design_key(X["N"], X["d"], X["ns"], X["a"])
                self.log.append(key)
                e = model.evaluate(*key[:3], key[3])
                out["F"], out["G"] = e["F"], e["G"]
        return Mixed()

    class Index(Problem):
        def __init__(self):
            super().__init__(n_var=4, n_obj=3, n_ieq_constr=3, xl=np.array([N_RANGE[0], 0, NS_RANGE[0], 0]),
                             xu=np.array([N_RANGE[1], n_d - 1, NS_RANGE[1], len(PATHS_ALLOWED) - 1]), vtype=int)
            self.log: list[tuple] = []

        def _evaluate(self, X, out, *args, **kwargs):
            X = np.asarray(X).astype(int)
            self.log.extend(design_key(x[0], x[1], x[2], PATHS_ALLOWED[x[3]]) for x in X)
            e = model.evaluate(X[:, 0], X[:, 1], X[:, 2], np.array(PATHS_ALLOWED)[X[:, 3]])
            out["F"], out["G"] = e["F"], e["G"]
    return Index()


def make_algorithm(name: str, encoding: str, pop_size: int):
    """NSGA-II or NSGA-III of pymoo with the same variation operators (pymoo's defaults for the variable types: SBX eta 15, polynomial mutation eta 20) and
    duplicate elimination, so that what differs is what the two algorithms differ in: NSGA-II keeps the best non-dominated ranks and breaks ties by crowding
    distance (binary tournament selection); NSGA-III keeps the ranks and breaks ties by niching around 91 reference directions (Das-Dennis, 12 divisions of the
    3 objectives; population 100 is at least the number of directions, as pymoo asks)."""
    from pymoo.algorithms.moo.nsga2 import NSGA2, binary_tournament
    from pymoo.algorithms.moo.nsga3 import NSGA3, comp_by_cv_then_random
    from pymoo.operators.selection.tournament import TournamentSelection
    from pymoo.util.ref_dirs import get_reference_directions
    ref_dirs = get_reference_directions("das-dennis", 3, n_partitions=12)
    select = TournamentSelection(func_comp=binary_tournament if name == "NSGA-II" else comp_by_cv_then_random)
    if encoding == "mixed":
        from pymoo.core.mixed import MixedVariableDuplicateElimination, MixedVariableMating, MixedVariableSampling
        common = dict(sampling=MixedVariableSampling(), mating=MixedVariableMating(selection=select, eliminate_duplicates=MixedVariableDuplicateElimination()),
                      eliminate_duplicates=MixedVariableDuplicateElimination())
    else:
        from pymoo.operators.crossover.sbx import SBX
        from pymoo.operators.mutation.pm import PM
        from pymoo.operators.repair.rounding import RoundingRepair
        from pymoo.operators.sampling.rnd import IntegerRandomSampling
        common = dict(sampling=IntegerRandomSampling(), selection=select, crossover=SBX(prob=0.9, eta=15, vtype=float, repair=RoundingRepair()),
                      mutation=PM(eta=20, vtype=float, repair=RoundingRepair()), eliminate_duplicates=True)
    if name == "NSGA-II":
        return NSGA2(pop_size=pop_size, **common)
    return NSGA3(ref_dirs=ref_dirs, pop_size=pop_size, **common)


def run_nsga(model: Model, grader: Grader, name: str, encoding: str, seed: int, pop_size: int, generations: int) -> dict:
    """One run: `generations` generations of population `pop_size` (n_gen counts the first, random one: 100 generations are 100 x 100 = 10,000 evaluations).
    After every generation the run is graded twice against the exact front: its POPULATION's non-dominated feasible designs (what pymoo returns as the result;
    at most 100 designs) and its ARCHIVE (the non-dominated feasible designs among everything it ever evaluated; what the geometry study keeps in its CSV).
    Returns the per-generation records and the final sets."""
    from pymoo.core.callback import Callback
    from pymoo.optimize import minimize
    problem = make_problem(model, encoding)

    def tuple_of(x) -> tuple:
        if isinstance(x, dict):
            return design_key(x["N"], x["d"], x["ns"], x["a"])
        return design_key(x[0], x[1], x[2], PATHS_ALLOWED[int(x[3])])

    class Recorder(Callback):
        def __init__(self):
            super().__init__()
            self.records, self.seen, self.done, self.archive = [], set(), 0, ArchiveFront()
            self.evaluated_F: list[np.ndarray] = []
            self.evaluated_X: list[tuple] = []

        def notify(self, algorithm):
            new = [k for k in dict.fromkeys(problem.log[self.done:]) if k not in self.seen]
            self.done = len(problem.log)
            if new:
                self.seen.update(new)
                arr = np.array(new)
                e = model.evaluate(arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3])
                ok = np.all(e["G"] <= 0.0, axis=1)
                self.archive.add(e["F"][ok], [k for k, o in zip(new, ok) if o])
                self.evaluated_F.append(e["F"][ok])
                self.evaluated_X.extend(k for k, o in zip(new, ok) if o)
            pop = algorithm.pop
            F, CV, X = pop.get("F"), pop.get("CV")[:, 0], [tuple_of(x) for x in pop.get("X")]
            ok = CV <= 0.0
            nd = pareto_mask(round_sig(F[ok])) if ok.any() else np.zeros(0, bool)
            Fp = F[ok][nd]
            self.pop_X = [x for x, o in zip(X, ok) if o]
            self.pop_X = [x for x, k in zip(self.pop_X, nd) if k]
            self.pop_F = Fp
            self.records.append({"gen": int(algorithm.n_gen), "n_eval": int(algorithm.evaluator.n_eval), "n_distinct": len(self.seen),
                                 "pop_feasible": int(ok.sum()), "pop": grader.grade(Fp), "archive": grader.grade(self.archive.F)})

    rec = Recorder()
    t0 = time.perf_counter()
    minimize(problem, make_algorithm(name, encoding, pop_size), ("n_gen", generations), seed=seed, verbose=False, callback=rec)
    return {"algorithm": name, "encoding": encoding, "seed": seed, "records": rec.records, "pop_F": rec.pop_F, "pop_X": rec.pop_X, "archive_F": rec.archive.F,
            "archive_X": rec.archive.X, "evaluated_F": np.vstack(rec.evaluated_F) if rec.evaluated_F else np.empty((0, 3)), "evaluated_X": rec.evaluated_X,
            "seconds": time.perf_counter() - t0, "n_eval": rec.records[-1]["n_eval"], "n_distinct": rec.records[-1]["n_distinct"]}


def summarise(runs: list[dict], budgets: tuple[int, ...]) -> list[dict]:
    """Mean, sample standard deviation, minimum and maximum over the seeds of every grade, for each algorithm, set (population / archive) and budget (the
    generation at which the grade is read: the runs record every generation, and a run to generation 100 passes through generation 50 exactly as a run
    stopped there would)."""
    rows = []
    for name in dict.fromkeys(r["algorithm"] for r in runs):
        mine = [r for r in runs if r["algorithm"] == name]
        for gens in budgets:
            recs = [next(x for x in r["records"] if x["gen"] == gens) for r in mine]
            for which in ("pop", "archive"):
                row = {"algorithm": name, "set": which, "generations": gens, "seeds": len(mine), "evaluations": recs[0]["n_eval"],
                       "distinct_designs_mean": float(np.mean([x["n_distinct"] for x in recs]))}
                for metric in ("hv_ratio", "igd", "hits", "hit_share", "n"):
                    v = np.array([x[which][metric] for x in recs], float)
                    row.update({f"{metric}_mean": float(v.mean()), f"{metric}_std": float(v.std(ddof=1)) if len(v) > 1 else 0.0, f"{metric}_min": float(v.min()),
                                f"{metric}_max": float(v.max())})
                rows.append(row)
    return rows


def describe_misses(model: Model, en: dict, grader: Grader, runs: list[dict]) -> dict:
    """What the exact points an algorithm missed look like. For every run, every distinct exact point missing from the final POPULATION's result is put in
    one of three kinds: 'lost' (the run evaluated it, and the population did not keep it: the population cannot hold more than 100 designs, the front has
    more points than that, and survival dropped it), 'visited r' (the run never evaluated that exact winding but did evaluate others with the same current, i.e.
    the same N / a: it had the turns right and the wire or strands slightly off; the loss of its best such design is compared with the exact one) and 'unvisited'
    (no design with that current was ever evaluated). Also which exact points were missed by every run of an algorithm."""
    F_ex = grader.front
    keys = [Grader.key(f) for f in F_ex]
    out = {}
    for name in dict.fromkeys(r["algorithm"] for r in runs):
        mine = [r for r in runs if r["algorithm"] == name]
        kinds = {"lost": 0, "visited_r": 0, "unvisited": 0}
        excess, missed_by = [], np.zeros(len(F_ex), int)
        for r in mine:
            pop_keys = {Grader.key(f) for f in r["pop_F"]}
            arch_keys = {Grader.key(f) for f in r["archive_F"]}
            EF = r["evaluated_F"]
            for j, (key, f) in enumerate(zip(keys, F_ex)):
                if key in pop_keys:
                    continue
                missed_by[j] += 1
                if key in arch_keys:
                    kinds["lost"] += 1
                    continue
                same_i = EF[np.isclose(EF[:, 1], f[1], rtol=1e-9)] if len(EF) else EF
                if len(same_i):
                    kinds["visited_r"] += 1
                    excess.append(float(same_i[:, 0].min() / f[0] - 1.0))
                else:
                    kinds["unvisited"] += 1
        n_runs = len(mine)
        r_all = model.at / F_ex[:, 1]                                       # N / a of every exact point (its current is AT a / N)

        def share(mask):
            return float(missed_by[mask].mean() / n_runs) if mask.any() else float("nan")
        by_r = [(lo, hi, int(m.sum()), share(m)) for lo, hi in ((0.0, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0), (20.0, np.inf)) for m in [(r_all >= lo) & (r_all < hi)]]
        by_loss = [(lo, hi, int(m.sum()), share(m)) for lo, hi in ((0.0, 57.0), (57.0, 58.0), (58.0, 60.0), (60.0, np.inf)) for m in [(F_ex[:, 0] >= lo) & (F_ex[:, 0] < hi)]]
        out[name] = {"mean_missed_per_run": float(sum(kinds.values()) / n_runs), "kinds_per_run": {k: v / n_runs for k, v in kinds.items()},
                     "missed_share_by_r": by_r, "missed_share_by_loss": by_loss,
                     "loss_excess_visited_median_pct": 100.0 * float(np.median(excess)) if excess else float("nan"),
                     "loss_excess_visited_max_pct": 100.0 * float(np.max(excess)) if excess else float("nan"),
                     "missed_by_every_run": [int(j) for j in np.flatnonzero(missed_by == n_runs)], "missed_by_some_run": int((missed_by > 0).sum()),
                     "missed_counts": missed_by.tolist()}
    return out


# ---------------------------------------------------------------------------------------------- the front as a table, and an example winding from it
def design_row(model: Model, en: dict, j: int) -> dict:
    """Everything about design j of the enumeration, as plain numbers."""
    e = en["e"]
    N, a, d_mm, ns = int(en["N"][j]), int(model.paths[en["path_idx"][j]]), float(model.cat["d_mm"][en["di"][j]]), int(en["ns"][j])
    return {"N": N, "a": a, "d_mm": d_mm, "d_over_mm": float(e["d_over"][j]), "ns": ns, "r_turns_over_paths": N / a, "series_turns": float(e["series"][j]),
            "strands_per_coil_side": N * ns, "fill": float(e["k"][j]), "R_hot_ohm": float(e["R"][j]), "R_20C_ohm": float(e["R"][j] * RHO_20 / model.rho),
            "i_peak_A": float(e["i"][j]), "i_rms_A": float(e["i"][j] / math.sqrt(2.0)), "p_dc_W": float(e["p_dc"][j]), "p_ac_W": float(e["p_ac"][j]),
            "p_loss_W": float(e["F"][j, 0]), "corner_rpm": float(e["rpm"][j]), "wire_index": int(en["di"][j])}


def terminal_voltage(model: Model, row: dict, rpm: float) -> float:
    """The peak phase voltage at the rated point and `rpm` for a winding (a row of design_row): sqrt((R i_d - w psi_q)^2 + (R i_q + w psi_d)^2)."""
    e = model.evaluate(row["N"], row["wire_index"], row["ns"], row["a"])
    w = 2.0 * math.pi * model.p.pole_pairs * rpm / 60.0
    return float(math.hypot(e["R"] * e["i_d"] - w * e["psi_q"], e["R"] * e["i_q"] + w * e["psi_d"]))


def pick_example(model: Model, en: dict, margin: float = 1.10) -> dict:
    """The example winding: the front member with the smallest peak current among those whose corner speed is at least `margin` x the reference speed (the rated
    torque held up to 3300 rpm with 10 % to spare on the 48 V bus), and, among windings that are one front point written several ways, the one with the fewest
    strands in hand, then the fewest paths. Returns its row, the front point's other ways of being written, and the need."""
    e, fi = en["e"], en["front"]
    ok = fi[e["rpm"][fi] >= margin * model.rpm_ref]
    met = bool(len(ok))
    if not met:                                                              # no front winding is fast enough on this bus: show the fastest one
        ok = fi[[int(np.argmax(e["rpm"][fi]))]]
    j0 = ok[np.argmin(e["i"][ok])]
    same = fi[np.all(round_sig(e["F"][fi]) == round_sig(e["F"][j0]), axis=1)]
    j = min(same, key=lambda s: (int(en["ns"][s]), int(model.paths[en["path_idx"][s]]), int(en["N"][s])))
    return {"index": int(j), "row": design_row(model, en, int(j)), "equivalent": [design_row(model, en, int(s)) for s in same], "rpm_needed": margin * model.rpm_ref,
            "margin": margin, "met": met}


def pick_for_current(model: Model, en: dict, i_max: float) -> dict | None:
    """The front member with the highest corner speed among those whose peak current is at most i_max (an inverter's limit); among windings that are one front
    point written several ways, the fewest strands in hand, then the fewest paths."""
    e, fi = en["e"], en["front"]
    ok = fi[e["i"][fi] <= i_max + 1e-9]
    if not len(ok):
        return None
    j0 = ok[np.argmax(e["rpm"][ok])]
    same = fi[np.all(round_sig(e["F"][fi]) == round_sig(e["F"][j0]), axis=1)]
    j = min(same, key=lambda s: (int(en["ns"][s]), int(model.paths[en["path_idx"][s]]), int(en["N"][s])))
    return design_row(model, en, int(j))


def practical_variant(model: Model, en: dict, j: int, max_strands: int = 6) -> dict | None:
    """The same peak current as design j (the same N / a) with at most `max_strands` strands in hand: the feasible design of least loss, and what it costs against
    design j (loss and corner speed)."""
    e = en["e"]
    same = np.flatnonzero(en["feasible"] & np.isclose(e["i"], e["i"][j], rtol=1e-9) & (en["ns"] <= max_strands))
    if not len(same):
        return None
    k = same[np.argmin(e["F"][same, 0])]
    row = design_row(model, en, int(k))
    row["loss_excess_pct"] = 100.0 * (e["F"][k, 0] / e["F"][j, 0] - 1.0)
    row["speed_change_pct"] = 100.0 * (e["rpm"][k] / e["rpm"][j] - 1.0)
    return row


# ---------------------------------------------------------------------------------------------- reports and files
class Report:
    """Prints every line and keeps it, so that the same text goes to the screen and to report.txt."""

    def __init__(self):
        self.lines: list[str] = []

    def say(self, text: str = "") -> None:
        print(text, flush=True)
        self.lines.append(text)

    def check(self, name: str, ok: bool, said: str) -> bool:
        self.say(f"{'PASS' if ok else 'FAIL'}  {name}: {said}")
        return bool(ok)


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def flat(row: dict, prefix: str = "") -> dict:
    """A nested dict of numbers as one flat row (keys joined with an underscore)."""
    out = {}
    for k, v in row.items():
        if isinstance(v, dict):
            out.update(flat(v, f"{prefix}{k}_"))
        else:
            out[f"{prefix}{k}"] = v
    return out


# ---------------------------------------------------------------------------------------------- the study
def out_folder(label: str) -> Path:
    """out\\optimization\\winding\\<label>\\, made if it is not there."""
    folder = paths.path("out") / "optimization" / "winding" / label
    folder.mkdir(parents=True, exist_ok=True)
    return folder


def section_rated(rep: Report, p: mg.Params, femm: dict, model: Model) -> dict:
    r = femm["rated"]
    rep.say("== 1. The rated point (FEMM 5.0, fine mesh, design.Evaluator) ==")
    rep.say(f"geometry: Params(barrier='none') with overrides {femm['overrides'] or '{}'}; {p.slots} slots, {p.pole_pairs * 2} poles, stack {p.stack:g} mm, "
            f"slot opening {p.slot_open:g} mm, coil layer {r['coil_layer_area_mm2']:.3f} mm^2 (check(): coil_layer_area_mm2)")
    rep.say(f"{r['ampere_turns_per_coil_side']:g} ampere-turns a coil side ({r['amps']:g} A x {r['turns']} turns), the geometry's own MTPA angle: gamma {r['gamma']:.3f} deg, "
            f"i_d {r['i_d']:.4f} A, i_q {r['i_q']:.4f} A")
    rep.say(f"mean torque {r['T_mean']:.5f} N.m = THE RATED TORQUE of the study (15 positions; ripple {r['ripple_pct']:.2f} %), pf {r['pf']:.4f}, "
            f"{r['nodes']:,} nodes, {r['seconds']:.1f} s")
    rep.say(f"flux linkages (means over the 15 positions, {model.series_ref} series turns, one path): psi_d {1e3 * r['psi_d']:.4f} mWb, psi_q {1e3 * r['psi_q']:.4f} mWb; "
            f"psi_pm {1e3 * r['psi_pm']:.3f} mWb, L_d {1e3 * r['L_d']:.4f} mH, L_q {1e3 * r['L_q']:.4f} mH (saliency {r['saliency']:.3f})")
    t_dq = model.torque_of(model.evaluate(model.turns_ref, 0, 1, 1))
    rep.say(f"torque from the flux linkages, 1.5 p (psi_d i_q - psi_q i_d) = {float(t_dq):.5f} N.m against the air-gap element's {r['T_mean']:.5f} "
            f"({100 * (float(t_dq) / r['T_mean'] - 1):+.3f} %)")
    return {}


def section_field(rep: Report, femm: dict, model: Model) -> None:
    s = femm["statistics"]
    rep.say("")
    rep.say("== 2. The field the copper sits in (FEMM, fine mesh, the rated point, one electrical period) ==")
    m = femm["mesh"]
    rep.say(f"{len(femm['positions']['A']) + len(femm['positions']['B'])} rotor positions: 24 at 15 deg (0, 15, ... 345 electrical degrees) and 24 more half a step on; the currents "
            f"turn with the rotor (FemmModel.solve); {m['nodes']:,} nodes, {m['triangles']:,} triangles, {m['coil_triangles']:,} of them in the {len(femm['layers'])} coil layers "
            f"of the 9 slots of the sector (all 36 slots follow from the periodic sector)")
    c = femm["femm_checks"]
    rep.say(f"B is read from FEMM's .ans (one value a first-order triangle) and checked against FEMM itself at the first position: B at {c['n_probe']} triangle centres "
            f"(mo_getb, smoothing off) agrees to {c['B_probe_max_diff_T']:.1e} T; the field energy of every layer (block integral 2, mu_r = 1) agrees with "
            f"the sum over triangles of area x B^2 / (2 mu_0) x depth to {c['energy_max_rel_diff']:.1e} (relative); the layers' areas (block integral 5) to "
            f"{abs(c['area_ratio_min'] - 1):.0e}")
    rep.say("what a strand feels is dB/dt, not B, so the statistic is G = the mean over the copper area of the mean square of dB/d(theta_e) [T^2 / rad^2], from each triangle's Fourier "
            "series (Parseval), times omega_e^2 for T^2/s^2:")
    for key, label in (("12", "12 positions, 30 deg"), ("24", "24 positions, 15 deg (brief)"), ("24_shifted", "24 positions, shifted half a step"), ("48", "48 positions, 7.5 deg (used)")):
        v = s[key]
        rep.say(f"  {label:36s}  G = {v['G']:.4e}   fundamental alone G1 = {v['G1']:.4e}   ({100 * v['G1'] / v['G']:.1f} % of G)")
    v = s["48"]
    if "G_forward_difference" in v:
        rep.say(f"  an independent estimate, forward differences between neighbouring samples, (B_j+1 - B_j)^2 / dtheta^2 averaged the same way: {v['G_forward_difference']:.4e} "
                f"({v['G_forward_difference'] / v['G']:.3f} of the spectral value; differences blur the harmonics: a harmonic k is read as (2 sin(k dtheta / 2) / dtheta)^2 instead of k^2)")
    share = np.array(v["G_k"]) / v["G"]
    top = [(k, 100 * share[k]) for k in np.argsort(-share)[:6] if share[k] > 0]
    rep.say("  harmonics of G (48 positions): " + ", ".join(f"k = {k}: {x:.1f} %" for k, x in sorted(top)) + f"; k >= 17: {100 * share[17:].sum():.1f} %; only odd k (half-wave symmetry)")
    g24, g48, g12 = s["24"]["G"], s["48"]["G"], s["12"]["G"]
    gain1, gain2 = g24 / g12 - 1.0, g48 / g24 - 1.0
    ratio = gain2 / gain1
    rep.say(f"  convergence: G(12) -> G(24) -> G(48) = {g12:.4e} -> {g24:.4e} ({100 * gain1:+.1f} %) -> {g48:.4e} ({100 * gain2:+.1f} %): each doubling adds {ratio:.2f} of the "
            f"last gain, so if that goes on the true G is about {100 * gain2 * ratio / (1 - ratio):.1f} % above the 48-position value (the loss it feeds is under 0.2 % of the total: it does not matter)")
    gl = np.array(v["G_layer"])
    rep.say(f"  by layer: the layers next to the bore (coil_k_0) have G = {1e4 * gl[0::2].min():.2f} ... {1e4 * gl[0::2].max():.2f} e-4, the deep layers (coil_k_1) "
            f"{1e4 * gl[1::2].min():.2f} ... {1e4 * gl[1::2].max():.2f} e-4: the strands near the slot opening feel {gl[0::2].mean() / gl[1::2].mean():.0f} times the dB/dt^2 of the deep ones")
    rep.say(f"  time-mean |B|^2 over the copper {v['B2_mean']:.3e} T^2 (B_rms {1e3 * math.sqrt(v['B2_mean']):.1f} mT), of which constant in time {v['B2_dc']:.1e}: nothing; "
            f"peak |B| of any triangle at any position {v['B_peak']:.3f} T")


def section_turns(rep: Report, femm: dict) -> dict:
    rep.say("")
    rep.say("== 3. The turns do not change the field (FEMM, fine mesh, the same current angle) ==")
    rep.say("same geometry, 5 / 10 / 20 turns a coil side at 20 / 10 / 5 A (100 ampere-turns), gamma "
            f"{femm['turns_check'][0]['gamma']:.3f} deg, the 15 ripple positions; psi_d x ratio = psi_d x (10 / turns)")
    rows = femm["turns_check"]
    ref = next(r for r in rows if r["turns"] == femm["turns_reference"])
    rep.say("  turns   amps   T_mean [N.m]   ripple [%]   psi_d [mWb]   psi_d x ratio   pf        max |dB| in the copper vs 10 turns")
    worst = {"T": 0.0, "ripple": 0.0, "psi": 0.0, "pf": 0.0, "B": 0.0}
    for r in rows:
        ratio = femm["turns_reference"] / r["turns"]
        rep.say(f"  {r['turns']:5d}  {r['amps']:5.1f}  {r['T_mean']:.10f}  {r['ripple_pct']:.8f}  {1e3 * r['psi_d']:.8f}  {1e3 * r['psi_d'] * ratio:.8f}  {r['pf']:.10f}  "
                f"{r['field_max_diff_T']:.1e} T (field up to {r['field_max_T']:.4f} T)")
        worst["T"] = max(worst["T"], abs(r["T_mean"] / ref["T_mean"] - 1))
        worst["ripple"] = max(worst["ripple"], abs(r["ripple_pct"] / ref["ripple_pct"] - 1))
        worst["psi"] = max(worst["psi"], abs(r["psi_d"] * ratio / ref["psi_d"] - 1))
        worst["pf"] = max(worst["pf"], abs(r["pf"] / ref["pf"] - 1))
        worst["B"] = max(worst["B"], r["field_max_diff_T"])
    rep.say(f"  largest relative differences from the 10-turn run: torque {worst['T']:.1e}, ripple {worst['ripple']:.1e}, psi_d x ratio {worst['psi']:.1e}, pf {worst['pf']:.1e}; "
            f"field {worst['B']:.1e} T. Same mesh, same ampere-turns bit for bit (20 x 5 = 10 x 10 = 5 x 20): the solutions are the same numbers, and the flux linkage is the "
            f"turns times the same flux")
    return worst


def section_model(rep: Report, model: Model, p: mg.Params, pa: dict) -> dict:
    rep.say("")
    rep.say("== 4. The winding model ==")
    g = model.geometry
    rep.say(f"bus {model.v_dc:g} V -> peak phase voltage V_dc / sqrt 3 = {model.v_max:.3f} V; reference speed {model.rpm_ref:g} rpm = {model.f_ref:g} Hz electrical "
            f"(omega = {model.w_ref:.1f} rad/s, p = {p.pole_pairs}); copper at {model.t_cu:g} deg C: rho = {model.rho:.4e} Ohm m (rho_20 {RHO_20:.4e}, alpha {ALPHA_20}), "
            f"sigma = {model.sigma:.4e} S/m")
    rep.say(f"parallel paths: motorgeo.winding gives {pa['coils']} coils a phase in {pa['distinct_emfs']} distinct emf phasors (equal magnitude, at {pa['emf_angles_deg']} deg), "
            f"repeating in {pa['periodic_sectors']} sectors of 90 degrees (gcd(slots, pole pairs) = {math.gcd(p.slots, p.pole_pairs)}). Splits into equal-emf groups, by brute force: "
            + ", ".join(f"a = {a}: {c}" for a, c in pa["counts"].items()) + f" -> allowed a = {pa['allowed']} (the divisors of 4: the paths must have equal emf, and the coils come in "
            f"4 copies of 3 phasors; 3, 6 and 12 would need unequal phasors to add up equal)")
    rep.say(f"end turn model: slot centroid radius {g['r_mean_mm']:.2f} mm, coil span (the chord of {p.coil_pitch} slot pitches) {g['span_mm']:.2f} mm, each end turn a half circle on it "
            f"({0.5 * math.pi * g['span_mm']:.1f} mm) plus 2 x {END_EXTENSION_MM:g} mm of straight leg = {g['l_end_mm']:.1f} mm; mean turn length l_mt = 2 x {p.stack:g} + 2 x "
            f"{g['l_end_mm']:.1f} = {g['l_mt_mm']:.1f} mm")
    rep.say(f"constraints: fill <= {FILL_MAX} (bare copper over one layer's {model.a_layer * 1e6:.3f} mm^2); the strand's grade 2 overall diameter + 2 x {CLEARANCE_MM:g} mm of clearance "
            f"<= the {p.slot_open:g} mm slot opening (overall <= {p.slot_open - 2 * CLEARANCE_MM:g} mm); standstill drop R i < V_max (else the rated current cannot be driven at all)")
    cat = model.cat
    rep.say("wire catalogue (IEC 60317 R 20, 0.200 ... 1.600 mm; grade 2 maximum overall diameter; E = Elektrisola IEC 60317 sheet, M = Multicomp Pro IEC 60317-0-1:2014 sheet): "
            + ", ".join(f"{d:.3f}/{o:.3f}{'' if pub else '*'}" for d, o, pub in zip(cat["d_mm"], cat["d_over_mm"], cat["published"]))
            + f"   (* = interpolated, not a table value: {int((~cat['published']).sum())} of {len(cat['d_mm'])})")
    return {}


def section_losses(rep: Report, model: Model, en: dict) -> dict:
    """The numbers behind 'the DC loss does not depend on N or a' and 'what the AC loss is'."""
    e = en["e"]
    f = en["feasible"]
    rep.say("")
    rep.say("  the DC loss at a given fill does not depend on N or a. With R = rho l_mt c (N / a) / (a n_s A_s), i = AT a / N and k = N n_s A_s / A_layer (c = 12 coils a phase, "
            "AT = 100):  1.5 R i^2 = 1.5 rho l_mt c AT^2 / (k A_layer), the turns and the paths cancel.")
    const = model.dc_loss_at_fill(1.0)
    dev = np.abs(e["p_dc"] * e["k"] / const - 1.0)
    rep.say(f"  over ALL {en['total']:,} combinations: P_dc x k / (1.5 rho l_mt c AT^2 / A_layer) = 1 to within {dev.max():.1e}; the closed form gives {model.dc_loss_at_fill(FILL_MAX):.3f} W "
            f"at k = {FILL_MAX} and {model.dc_loss_at_fill(0.05):.0f} W at the 5 % fill of the placeholder winding")
    rep.say("  one fill, many windings (wire 0.200 mm, N x n_s = 120, k = "
            f"{float(120 * math.pi * 0.2e-3 ** 2 / 4 / model.a_layer):.4f}): the same loss, different current and speed")
    rep.say("    N  n_s  a   i [A]   R [Ohm]   P_dc [W]   corner [rpm]")
    for N, ns, a in ((4, 30, 1), (6, 20, 1), (6, 20, 2), (10, 12, 1), (10, 12, 4), (20, 6, 2), (15, 8, 4), (30, 4, 4)):
        x = model.evaluate(N, 0, ns, a)
        rep.say(f"  {N:3d}  {ns:3d}  {a}  {float(x['i']):6.2f}  {float(x['R']):8.4f}  {float(x['p_dc']):9.4f}   {float(x['rpm']):9.1f}")
    share = e["p_ac"][f] / e["F"][f, 0]
    rep.say(f"  AC (proximity) loss: over the {en['n_feasible']:,} allowed combinations it is {e['p_ac'][f].min():.4f} ... {e['p_ac'][f].max():.4f} W against a DC loss of "
            f"{e['p_dc'][f].min():.1f} ... {e['p_dc'][f].max():.1f} W: at most {100 * share.max():.3f} % of the total (d^2 grows it, but d <= 1.4 mm here)")
    val = model.ac_validity()
    worst_top = min(v["ratio_top"] for v in val)
    allowed = [v for v in val if v["d_over_mm"] + 2 * CLEARANCE_MM <= model.slot_open]
    k_max = len(model.G_k) - 1
    rep.say(f"  skin depth at the fundamental {1e3 * skin_depth(model.f_ref, model.rho):.2f} mm ({model.f_ref:g} Hz) and {1e3 * skin_depth(k_max * model.f_ref, model.rho):.2f} mm at "
            f"{k_max * model.f_ref / 1e3:.1f} kHz (the {k_max}th harmonic, the highest the 48 positions resolve); the largest allowed wire has d / delta = "
            f"{max(v['d_over_delta_fund'] for v in allowed):.2f} at the fundamental and {max(v['d_over_delta_top'] for v in allowed):.2f} at the {allowed[0]['k_top']}th harmonic (the "
            f"highest that carries 1 % of G). The exact formula for a round wire in a transverse field (Bessel functions, Model.ac_validity) over the low-frequency one: "
            f"{min(v['ratio_weighted'] for v in allowed):.4f} (G-weighted, the allowed wires), {min(v['ratio_top'] for v in allowed):.4f} at the highest of those harmonics; the "
            f"low-frequency formula is within 1 % for d / delta <= 1 and within 5 % up to about 1.65: here it "
            f"{'holds for every allowed wire' if all(v['valid_5pct'] for v in allowed) else 'fails for some allowed wire'}; it would stop holding (5 %) for 1.4 mm wire above about "
            f"{(1.65 / 1.4e-3) ** 2 * model.rho / (math.pi * MU0) / 1e3:.1f} kHz")
    return {"ac_share_max": float(share.max()), "validity": val, "worst_ratio_top": worst_top}


def section_enumeration(rep: Report, model: Model, en: dict) -> dict:
    e, fi = en["e"], en["front"]
    rep.say("")
    rep.say("== 5. Enumeration (numpy, every allowed combination) ==")
    n_d = len(model.d)
    rep.say(f"{len(range(N_RANGE[0], N_RANGE[1] + 1))} values of N x {n_d} wires x {len(range(NS_RANGE[0], NS_RANGE[1] + 1))} values of n_s x {len(model.paths)} values of a = "
            f"{en['total']:,} combinations, all evaluated in {1e3 * en['seconds_eval']:.0f} ms; {en['n_feasible']:,} are allowed ({100 * en['n_feasible'] / en['total']:.1f} %); "
            f"violations (a combination can break several): fill > {FILL_MAX}: {en['violations']['fill']:,}, slot opening: {en['violations']['slot']:,}, standstill drop: "
            f"{en['violations']['voltage']:,}")
    first, group = distinct_rows(e["F"][fi])
    rep.say(f"exact non-dominated set by one sweep over the {en['n_feasible']:,} allowed designs ({1e3 * en['seconds_front']:.0f} ms): {len(fi)} designs = {len(first)} distinct objective vectors "
            f"(the same winding written twice, (N, a, n_s) = (2N, 2a, n_s / 2), is one point: {len(fi) - len(first)} extra copies)")
    v = en["verify"]
    rep.say(f"verified against ALL allowed designs: front members beaten by any design: {v['front_members_beaten_by_any_design']}; designs off the front not beaten by a front member: "
            f"{v['other_designs_not_beaten_by_front']} -> {'exact' if v['ok'] else 'NOT EXACT'}")
    ff = e["F"][fi]
    r_values = sorted({round(float(en['N'][j] / model.paths[en['path_idx'][j]]), 6) for j in fi})
    ac_share = float((e["p_ac"][fi] / e["F"][fi, 0]).max())
    rpm0 = 60.0 / (2.0 * math.pi * model.p.pole_pairs) * model.v_max / np.hypot(e["psi_d"][fi], e["psi_q"][fi])      # the corner speed of a winding with no resistance
    lost = 1.0 - e["rpm"][fi] / rpm0
    one_each = "one winding each" if len(first) == len(r_values) else f"{len(first)} distinct points: some N / a have more than one"
    rep.say(f"the front has {len(r_values)} distinct values of N / a ({one_each}) from {r_values[0]:g} to {r_values[-1]:g} (current {ff[:, 1].min():.2f} ... {ff[:, 1].max():.0f} A, corner "
            f"speed {-ff[:, 2].max():.0f} ... {-ff[:, 2].min():.0f} rpm); the loss along it varies from {ff[:, 0].min():.2f} to {ff[:, 0].max():.2f} W "
            + ("only because the best fill that the catalogue allows changes with N (the AC loss is at most " if ac_share < 0.01 else "because the best fill that the catalogue allows changes "
               "with N and because of the AC loss (up to ")
            + f"{100 * ac_share:.3f} % of it): fill {e['k'][fi].min():.3f} ... {e['k'][fi].max():.3f}; the resistance lowers the corner speed below V_max / |psi| by "
            f"{100 * lost.min():.0f} % (N / a = {r_values[0]:g}) ... {100 * lost.max():.0f} % (N / a = {r_values[-1]:g}), {100 * np.median(lost):.0f} % at the median; wires on the front: "
            + ", ".join(f"{d:.3f}" for d in sorted({float(model.cat['d_mm'][en['di'][j]]) for j in fi})) + " mm")
    return {"distinct": len(first), "front": len(fi), "n_values_r": len(r_values)}


def section_nsga(rep: Report, rows: list[dict], runs: list[dict], grader: Grader, misses: dict, pop: int, budgets, total_combos: int, encoding: str) -> None:
    K = len(grader.keys)
    rep.say("")
    rep.say(f"== 6. NSGA-II and NSGA-III on the same problem ({encoding} encoding, population {pop}, {len(set(r['seed'] for r in runs))} seeds each) ==")
    rep.say(f"graded against the exact front: {K} distinct points; hypervolume ratio = HV of the set / HV of the exact front, objectives scaled 0 ... 1 between the front's ideal and nadir, "
            f"one common reference point (1.1, 1.1, 1.1), pymoo's HV (exact HV of the front: {grader.hv_exact:.6f}); IGD in the same scaled space; 'exact' = the front's objective "
            f"vector reproduced to 10 digits. A population of {pop} can hold at most {pop} of the {K} points ({100 * min(1, pop / K):.0f} %).")
    rep.say("  gens  evals (share of combinations)  algorithm   set        HV ratio mean (std) [min .. max]       IGD mean (std)         exact points found mean (std) [min .. max] of "
            f"{K}  (share)")
    for row in rows:
        rep.say(f"  {row['generations']:4d}  {row['evaluations']:6d} ({100 * row['evaluations'] / total_combos:4.1f} %; {row['distinct_designs_mean']:.0f} distinct)  {row['algorithm']:9s}   "
                f"{'population' if row['set'] == 'pop' else 'archive':10s} {row['hv_ratio_mean']:.4f} ({row['hv_ratio_std']:.4f}) [{row['hv_ratio_min']:.4f} .. {row['hv_ratio_max']:.4f}]   "
                f"{row['igd_mean']:.5f} ({row['igd_std']:.5f})   {row['hits_mean']:6.1f} ({row['hits_std']:4.1f}) [{row['hits_min']:3.0f} .. {row['hits_max']:3.0f}]  "
                f"({100 * row['hit_share_mean']:.1f} %)")
    secs = {name: np.mean([r["seconds"] for r in runs if r["algorithm"] == name]) for name in dict.fromkeys(r["algorithm"] for r in runs)}
    rep.say("  time of one run: " + ", ".join(f"{k} {v:.1f} s" for k, v in secs.items()))
    rep.say("  what the misses look like (final population, per run, after the last generation):")
    for name, m in misses.items():
        k = m["kinds_per_run"]
        rep.say(f"    {name}: {m['mean_missed_per_run']:.1f} of {K} exact points missed per run: {k['lost']:.1f} evaluated but dropped by survival (the population holds at most {pop}), "
                f"{k['visited_r']:.1f} never evaluated as such, though designs with the same N / a (the same current) were (the best of them {m['loss_excess_visited_median_pct']:.2f} % "
                f"median, {m['loss_excess_visited_max_pct']:.2f} % worst above the exact loss), {k['unvisited']:.1f} whose N / a was never tried; found by at least one of the "
                f"{len(runs) // 2} runs: {len(grader.keys) - len(m['missed_by_every_run'])} points (missed by every run: {len(m['missed_by_every_run'])}); found by every run: "
                f"{len(grader.keys) - m['missed_by_some_run']}")
        rep.say("      share of the runs that miss a point, by N / a: " + ", ".join(f"{lo:g} ... {hi:g}: {100 * sh:.0f} % ({n} points)" for lo, hi, n, sh in m["missed_share_by_r"] if n)
                + "; by the exact loss [W]: " + ", ".join(f"{lo:g} ... {hi:g}: {100 * sh:.0f} % ({n})" for lo, hi, n, sh in m["missed_share_by_loss"] if n))


def section_example(rep: Report, model: Model, en: dict, ex: dict, practical: dict | None) -> None:
    r = ex["row"]
    rep.say("")
    rep.say("== 7. An example winding from the exact front (EXAMPLE ONLY: the turns are sized again once the stack is ruled) ==")
    if not ex.get("met", True):
        rep.say(f"NO front winding reaches {ex['rpm_needed']:.0f} rpm on this bus: the fastest one is shown instead")
    rep.say(f"requirement: rated torque {model.torque:.4f} N.m held up to {model.rpm_ref:g} rpm with {100 * (ex['margin'] - 1):.0f} % to spare on the {model.v_dc:g} V bus "
            f"(corner speed >= {ex['rpm_needed']:.0f} rpm), the smallest inverter current that does it; windings that are the same front point written several ways are "
            f"settled by the fewest strands in hand, then the fewest paths")
    rep.say(f"  N = {r['N']} turns a coil side, a = {r['a']} parallel path(s), {r['ns']} strands of {r['d_mm']:.3f} mm (overall {r['d_over_mm']:.3f} mm, {model.slot_open - r['d_over_mm']:.3f} mm "
            f"narrower than the {model.slot_open:g} mm opening), {r['strands_per_coil_side']} strands in a coil side, fill {r['fill']:.4f}")
    rep.say(f"  series turns a phase 12 N / a = {r['series_turns']:g}; R = {r['R_hot_ohm']:.4f} Ohm at {model.t_cu:g} deg C ({r['R_20C_ohm']:.4f} at 20); peak current i = 100 a / N = "
            f"{r['i_peak_A']:.2f} A ({r['i_rms_A']:.2f} A rms)")
    v3000 = terminal_voltage(model, r, model.rpm_ref)
    rep.say(f"  at the rated point: corner speed {r['corner_rpm']:.0f} rpm; terminal voltage at {model.rpm_ref:g} rpm {v3000:.2f} V of {model.v_max:.2f} V "
            f"({100 * (1 - v3000 / model.v_max):.1f} % to spare); copper loss at {model.rpm_ref:g} rpm {r['p_loss_W']:.3f} W = DC {r['p_dc_W']:.3f} + AC {r['p_ac_W']:.4f} "
            )
    if len(ex["equivalent"]) > 1:
        rep.say("  the same front point, written other ways: " + "; ".join(f"N = {q['N']}, a = {q['a']}, n_s = {q['ns']}" for q in ex["equivalent"] if q != r))
    if practical and r["ns"] > 6:
        rep.say(f"  a practical variant at the same current (N / a = {practical['r_turns_over_paths']:g}), at most 6 strands in hand: N = {practical['N']}, a = {practical['a']}, "
                f"{practical['ns']} strands of {practical['d_mm']:.3f} mm, fill {practical['fill']:.4f}: loss {practical['p_loss_W']:.3f} W ({practical['loss_excess_pct']:+.2f} %), corner speed "
                f"{practical['corner_rpm']:.0f} rpm ({practical['speed_change_pct']:+.2f} %)")
    elif practical:
        rep.say(f"  (it already has {r['ns']} strands in hand, within the limit of 6 assumed for a practical variant: no variant is needed)")
    rep.say("  the same front for other inverters (the front member with the highest corner speed within a peak-current limit; equal points: fewest strands, then fewest paths):")
    rep.say("    limit [A peak]   N  a  strands x d [mm]   fill    i [A]   corner [rpm]   loss [W]")
    for limit in (10, 15, 20, 30, 50, 100):
        q = pick_for_current(model, en, limit)
        if q:
            rep.say(f"    {limit:12d}  {q['N']:3d}  {q['a']}  {q['ns']:3d} x {q['d_mm']:.3f}      {q['fill']:.4f}  {q['i_peak_A']:6.2f}  {q['corner_rpm']:11.0f}   {q['p_loss_W']:.3f}")


def parse(argv=None):
    ap = argparse.ArgumentParser(prog="python -m optimization.winding", description="The winding study: the exact front by enumeration, graded against NSGA-II and NSGA-III.")
    ap.add_argument("--label", default="reference", help="the output folder's name: out/optimization/winding/<label>/ (default: reference)")
    ap.add_argument("--overrides", default="{}", help="JSON of motorgeo.Params overrides: the geometry, e.g. '{\"mag_y\": 5.5}' (default: {}, the reference)")
    ap.add_argument("--reuse-femm", action="store_true", help="read femm.json of this label instead of running the FEMM stage (no solver slot needed)")
    ap.add_argument("--seeds", type=int, default=5, help="seeds per algorithm (default 5: 1 ... 5)")
    ap.add_argument("--pop", type=int, default=100, help="population size (default 100)")
    ap.add_argument("--gens", type=int, default=100, help="generations (default 100: 10,000 evaluations)")
    ap.add_argument("--budgets", default="50,100", help="generations at which the grades are read (default 50,100)")
    ap.add_argument("--encoding", default="both", choices=["mixed", "index", "both"], help="mixed: pymoo's Integer and Choice variables; index: plain integers; both: the two, one after the other (default both)")
    ap.add_argument("--vdc", type=float, default=V_DC, help="DC bus [V] (default 48)")
    ap.add_argument("--temp", type=float, default=T_COPPER, help="copper temperature [deg C] (default 100)")
    ap.add_argument("--rpm", type=float, default=RPM_REF, help="reference speed of the loss objective [rpm] (default 3000)")
    ap.add_argument("--no-figures", action="store_true", help="skip the PNGs")
    return ap.parse_args(argv)


def main(argv=None) -> int:
    args = parse(argv)
    try:
        sys.stdout.reconfigure(errors="replace")
    except Exception:
        pass
    t_start = time.perf_counter()
    overrides = design.decode(json.loads(args.overrides))
    folder = out_folder(args.label)
    rep = Report()
    p = design.params(overrides)
    screened = design.screen(overrides)
    if not screened["feasible"]:
        rep.say(f"the geometry {overrides} cannot be built or breaks a limit: {screened['reason']}")
        return 2
    rep.say(f"winding study, label '{args.label}', geometry overrides {overrides or '{}'}, output {folder}")
    # ---- 1-3: FEMM (or the saved numbers)
    from optimization import winding_femm
    if args.reuse_femm:
        femm = json.loads((folder / "femm.json").read_text(encoding="utf-8"))
        sec = femm.get("seconds", {})
        rep.say(f"FEMM stage: read from {folder / 'femm.json'} (--reuse-femm); it was made in one pool session with {femm.get('solves', {}).get('total', '?')} solves, "
                f"{sec.get('total', float('nan')) - sec.get('waiting_for_slot', 0.0):.1f} s of FEMM after {sec.get('waiting_for_slot', float('nan')):.1f} s of waiting for the solver slot")
    else:
        t = time.perf_counter()
        data = winding_femm.run_stage(overrides, folder, args.label, log=rep.say)
        femm = json.loads(json.dumps(winding_femm.jsonable(data)))
        rep.say(f"FEMM stage: {femm['solves']['total']} solves in one pool session ({femm['seconds']['total']:.1f} s, of which {femm['seconds']['waiting_for_slot']:.1f} s waiting for the "
                f"solver slot), {time.perf_counter() - t:.1f} s wall")
    model = Model(p, femm, v_dc=args.vdc, t_cu=args.temp, rpm_ref=args.rpm)
    pa = allowed_paths(p)
    checks: list[bool] = []
    rated = section_rated(rep, p, femm, model)
    section_field(rep, femm, model)
    worst = section_turns(rep, femm)
    section_model(rep, model, p, pa)
    # ---- 5: the enumeration, and what it shows about the losses
    en = enumerate_all(model)
    en["verify"] = verify_front(en["e"]["F"][en["feasible"]], np.searchsorted(np.flatnonzero(en["feasible"]), en["front"]))
    losses = section_losses(rep, model, en)
    enum = section_enumeration(rep, model, en)
    F_front = en["e"]["F"][en["front"]]
    grader = Grader(F_front)
    # ---- 6: NSGA-II and NSGA-III
    budgets = tuple(int(b) for b in args.budgets.split(",") if int(b) <= args.gens) or (args.gens,)
    encodings = ["mixed", "index"] if args.encoding == "both" else [args.encoding]
    all_runs: dict[str, list[dict]] = {}
    all_rows: dict[str, list[dict]] = {}
    all_miss: dict[str, dict] = {}
    for enc in encodings:
        runs = []
        t = time.perf_counter()
        for name in ("NSGA-II", "NSGA-III"):
            for seed in range(1, args.seeds + 1):
                runs.append(run_nsga(model, grader, name, enc, seed, args.pop, args.gens))
        rep.say(f"  ({enc}: {len(runs)} runs in {time.perf_counter() - t:.1f} s)")
        all_runs[enc] = runs
        all_rows[enc] = summarise(runs, budgets)
        all_miss[enc] = describe_misses(model, en, grader, runs)
        section_nsga(rep, all_rows[enc], runs, grader, all_miss[enc], args.pop, budgets, en["total"], enc)
    # ---- 7: the example
    ex = pick_example(model, en)
    practical = practical_variant(model, en, ex["index"])
    section_example(rep, model, en, ex, practical)
    # ---- checks
    rep.say("")
    rep.say("== checks ==")
    r = femm["rated"]
    if not overrides and args.vdc == V_DC:
        checks.append(rep.check("the reference's rated point is the brief's (gamma 22.46, T 1.12880, ripple 18.45 %, psi_pm 16.638 mWb, pf 0.810)",
                                abs(r["gamma"] - 22.46) < 0.01 and abs(r["T_mean"] - 1.12880) < 2e-5 and abs(r["ripple_pct"] - 18.45) < 0.01
                                and abs(1e3 * r["psi_pm"] - 16.638) < 0.001 and abs(r["pf"] - 0.810) < 0.001,
                                f"gamma {r['gamma']:.3f}, T {r['T_mean']:.5f}, ripple {r['ripple_pct']:.2f} %, psi_pm {1e3 * r['psi_pm']:.3f} mWb, pf {r['pf']:.4f}"))
    checks.append(rep.check("the turns check: torque, ripple, psi_d x ratio and pf agree to 1e-9, the field to 1e-9 T",
                            max(worst["T"], worst["ripple"], worst["psi"], worst["pf"]) < 1e-9 and worst["B"] < 1e-9,
                            f"largest relative difference {max(worst['T'], worst['ripple'], worst['psi'], worst['pf']):.1e}, field {worst['B']:.1e} T"))
    c = femm["femm_checks"]
    checks.append(rep.check("B from the .ans equals FEMM's own (mo_getb) and the layer energies equal FEMM's block integral 2",
                            c["B_probe_max_diff_T"] < 1e-9 and c["energy_max_rel_diff"] < 1e-9, f"{c['B_probe_max_diff_T']:.1e} T, {c['energy_max_rel_diff']:.1e}"))
    checks.append(rep.check("the allowed parallel paths derived from motorgeo.winding are (1, 2, 4)", tuple(pa["allowed"]) == PATHS_ALLOWED, f"{pa['allowed']}"))
    tq = model.torque_of(en["e"])
    checks.append(rep.check("the torque is the same for every one of the combinations (1.5 p (psi_d' i_q' - psi_q' i_d') is winding-independent)",
                            float(np.ptp(tq) / abs(tq.mean())) < 1e-12, f"spread {float(np.ptp(tq) / abs(tq.mean())):.1e} of {float(tq.mean()):.6f} N.m"))
    dev = np.abs(en["e"]["p_dc"] * en["e"]["k"] / model.dc_loss_at_fill(1.0) - 1.0).max()
    checks.append(rep.check("the DC loss depends only on the fill: P_dc x k is the same for all combinations", dev < 1e-12, f"largest deviation {dev:.1e}"))
    a1 = model.evaluate(3, 7, 4, 1)["F"]
    a2, a3 = model.evaluate(6, 7, 2, 2)["F"], model.evaluate(12, 7, 1, 4)["F"]
    checks.append(rep.check("the same winding written three ways, (3, a=1, n_s=4), (6, a=2, n_s=2), (12, a=4, n_s=1), gives bit-identical objectives",
                            np.array_equal(a1, a2) and np.array_equal(a1, a3), f"{a1.tolist()}"))
    rng = np.random.default_rng(7)
    pick = rng.choice(np.flatnonzero(en["feasible"]), 6, replace=False)
    volts = []
    for j in pick:
        row = design_row(model, en, int(j))
        volts.append(terminal_voltage(model, row, row["corner_rpm"]))
    checks.append(rep.check("at its corner speed every winding's terminal voltage is the limit (the quadratic's root, put back in the voltage equation)",
                            max(abs(v / model.v_max - 1) for v in volts) < 1e-9, f"6 random allowed designs: largest deviation {max(abs(v / model.v_max - 1) for v in volts):.1e} of {model.v_max:.3f} V"))
    rng2 = np.random.default_rng(3)
    same = 0
    for _ in range(100):
        F_test = rng2.integers(0, 6, size=(int(rng2.integers(1, 300)), 3)).astype(float)       # many ties and duplicates
        same += int(np.array_equal(pareto_mask(F_test), design.nondominated(F_test)))
    checks.append(rep.check("pareto_mask agrees with design.nondominated (the O(n^2) definition) on 100 random sets full of ties", same == 100, f"{same} of 100 identical"))
    checks.append(rep.check("the front is exact (verified against all allowed designs)", en["verify"]["ok"],
                            f"{en['verify']['front_members_beaten_by_any_design']} front members beaten, {en['verify']['other_designs_not_beaten_by_front']} others unbeaten"))
    hv_slice = hv3_slicing(grader.scale(F_front), grader.ref)
    checks.append(rep.check("pymoo's hypervolume of the exact front equals an independent slicing computation", abs(hv_slice / grader.hv_exact - 1) < 1e-9,
                            f"{grader.hv_exact:.9f} against {hv_slice:.9f}"))
    sample = [x for runs in all_runs.values() for r_ in runs[:2] for x in r_["archive_X"][:40]]
    if sample:
        arr = np.array(sample)
        path_idx = np.searchsorted(model.paths, arr[:, 3])
        flat_i = (((arr[:, 0] - 1) * len(model.d) + arr[:, 1]) * NS_RANGE[1] + (arr[:, 2] - 1)) * len(model.paths) + path_idx
        ok = np.array_equal(en["e"]["F"][flat_i], model.evaluate(arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3])["F"]) and bool(en["feasible"][flat_i].all())
        checks.append(rep.check("the designs the optimizers found are the enumeration's designs, objective for objective (same function, bit for bit)", ok, f"{len(arr)} archive designs compared"))
    final_hits = {(enc, name): max(r_["records"][-1]["archive"]["hits"] for r_ in runs if r_["algorithm"] == name) for enc, runs in all_runs.items()
                  for name in ("NSGA-II", "NSGA-III")}
    checks.append(rep.check("no optimizer 'finds' more exact points than there are (grading sanity)", all(h <= len(grader.keys) for h in final_hits.values()),
                            ", ".join(f"{e_}/{n_}: {h}" for (e_, n_), h in final_hits.items()) + f" of {len(grader.keys)}"))
    rep.say(f"{sum(checks)} of {len(checks)} checks pass; total {time.perf_counter() - t_start:.1f} s")
    # ---- the files
    ex_row = ex["row"]
    results = {"label": args.label, "overrides": overrides, "args": vars(args), "rated": femm["rated"], "femm_checks": femm["femm_checks"], "statistics": femm["statistics"],
               "turns_check": femm["turns_check"], "paths": pa, "catalogue": {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in model.cat.items()},
               "constants": {"v_dc": model.v_dc, "v_max": model.v_max, "rpm_ref": model.rpm_ref, "f_ref_hz": model.f_ref, "t_cu": model.t_cu, "rho": model.rho, "fill_max": FILL_MAX,
                             "clearance_mm": CLEARANCE_MM, "end_extension_mm": END_EXTENSION_MM, "slot_open_mm": model.slot_open, "a_layer_mm2": model.a_layer * 1e6,
                             "G": model.G, "end_turn": model.geometry},
               "enumeration": {"total": en["total"], "feasible": en["n_feasible"], "violations": en["violations"], "front": len(en["front"]), **enum,
                               "seconds_eval": en["seconds_eval"], "seconds_front": en["seconds_front"], "verify": en["verify"]},
               "losses": {k: v for k, v in losses.items() if k != "validity"}, "ac_validity": losses["validity"],
               "nsga": {enc: {"summary": all_rows[enc], "misses": all_miss[enc]} for enc in encodings}, "example": {"row": ex_row, "equivalent": ex["equivalent"], "practical": practical},
               "checks": {"passed": int(sum(checks)), "total": len(checks)}, "seconds_total": time.perf_counter() - t_start}
    (folder / "results.json").write_text(json.dumps(winding_femm.jsonable(results), indent=1), encoding="utf-8")
    fr = en["front"]
    first, group = distinct_rows(en["e"]["F"][fr])
    rank = {int(group[f]): k for k, f in enumerate(first)}                 # the point's number: the order of first appearance, the order the grader and the figures use
    front_rows = [{**design_row(model, en, int(j)), "point": rank[int(g)]} for j, g in zip(fr, group)]
    write_csv(folder / "front_exact.csv", front_rows)
    write_csv(folder / "wire_catalogue.csv", [{"d_mm": d, "d_over_mm": o, "published": bool(pub), "source": s} for d, o, pub, s in
                                                zip(model.cat["d_mm"], model.cat["d_over_mm"], model.cat["published"], model.cat["source"])])
    np.savez_compressed(folder / "enumeration_all.npz", N=en["N"], wire_index=en["di"], ns=en["ns"], paths_index=en["path_idx"], feasible=en["feasible"], F=en["e"]["F"], G=en["e"]["G"],
                        fill=en["e"]["k"], R=en["e"]["R"], i=en["e"]["i"], p_dc=en["e"]["p_dc"], p_ac=en["e"]["p_ac"], rpm=en["e"]["rpm"], front=en["front"])
    run_rows, hist_rows = [], []
    for enc, runs in all_runs.items():
        for r_ in runs:
            last = r_["records"][-1]
            run_rows.append({"encoding": enc, "algorithm": r_["algorithm"], "seed": r_["seed"], "generations": last["gen"], "evaluations": r_["n_eval"], "distinct_designs": r_["n_distinct"],
                             "seconds": r_["seconds"], **flat({"pop": last["pop"], "archive": last["archive"]})})
            hist_rows.extend({"encoding": enc, "algorithm": r_["algorithm"], "seed": r_["seed"], "gen": x["gen"], "n_eval": x["n_eval"], "n_distinct": x["n_distinct"],
                              **flat({"pop": x["pop"], "archive": x["archive"]})} for x in r_["records"])
    write_csv(folder / "nsga_runs.csv", run_rows)
    write_csv(folder / "nsga_history.csv", hist_rows)
    write_csv(folder / "nsga_summary.csv", [dict(enc=enc, **row) for enc in encodings for row in all_rows[enc]])
    # ---- the figures
    if not args.no_figures:
        from optimization import winding_figures
        t = time.perf_counter()
        made = winding_figures.make_all({"model": model, "en": en, "grader": grader, "runs": all_runs[encodings[0]], "example": ex, "practical": practical,
                                         "pop": args.pop, "folder": folder, "label": args.label, "encoding": encodings[0]})
        rep.say(f"figures ({time.perf_counter() - t:.1f} s): " + ", ".join(str(m) for m in made))
    (folder / "report.txt").write_text("\n".join(rep.lines) + "\n", encoding="utf-8")
    return 0 if all(checks) else 1


if __name__ == "__main__":
    sys.exit(main())
