"""Shared settings and helpers of learn/dq_torque, the FEMM runs behind the project note
Documentation/Torque in d-q - field, current and saliency.

The rotor is the band-free one: Params(barrier="steel"), the hand-drawn rotor with its air band filled with steel, the
same model (and the same mesh) as the Air band study's "no band" column. Two variants of it are solved:
  "on":  N45 as the model has it, mu_r 1.05 and coercivity H_c = 1,030,272 A/m;
  "off": the same file with one number changed, H_c = 0. Same slots, same mesh, same mu_r: magnet material with no
         magnetisation, so whatever changes between the two is the magnets' flux and nothing else.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))     # the code folder, wherever this is run from
import numpy as np
from shapely.ops import nearest_points

import motorgeo as mg
import paths

OUT = paths.path("out") / "learn" / "dq_torque"
FIELDS = OUT / "fields"

I_RATED = 10.0                            # A peak: 100 ampere-turns per coil side with the 10-turn placeholder winding
GAMMAS = np.arange(0, 90.1, 7.5)          # current angle, from +q toward -d (0 = all q current), as barrier_study.py
MEAN_POS = tuple(range(0, 60, 10))        # 6 rotor positions (electrical degrees) over the 60 deg ripple period: their
                                          # average cancels the 6th, 12th and 18th torque harmonics
CURRENTS = (5.0, 10.0, 20.0, 30.0, 40.0)  # A peak, for L_d and L_q against current
HC_N45 = 1030272                          # A/m, the coercivity motorgeo.to_femm() gives N45
VARIANTS = {"on": HC_N45, "off": 0}
GAP_DEG = list(np.arange(0, 90, 0.25))    # air-gap samples over the 90 deg sector: one electrical period (p = 4)


def params():
    return mg.Params(barrier="steel")


def series_turns(p):
    """Series turns per phase: 10 turns per coil side, 24 coil sides per phase (36 slots, two layers), two sides a turn."""
    return p.turns * (2 * p.slots // 3) // 2


def probe_points(p, regions):
    """The middle of pole 0's two outer bridges and of its centre rib, in drawn coordinates [mm].

    The rotor side of the air-gap element never moves in the mesh (FEMM turns the rotor inside the AGE), so these points
    are in pole 0's steel at every rotor angle. Pole 0's d axis is at 90/p = 22.5 deg as drawn."""
    poly = {r.name: r.polygon() for r in regions}
    pts = {}
    for tag, side in (("dn", "outer bridge, right"), ("up", "outer bridge, left")):
        x, y = max(poly[f"magnet_0_{tag}"].exterior.coords, key=lambda c: math.hypot(*c))   # slot corner nearest the surface
        r = math.hypot(x, y)
        s = (r + (p.r_ro - r) / 2) / r                                                       # halfway to the rotor surface
        pts[side] = (x * s, y * s)
    a, b = nearest_points(poly["magnet_0_up"], poly["magnet_0_dn"])                        # where the two slots come closest
    pts["centre rib"] = ((a.x + b.x) / 2, (a.y + b.y) / 2)
    return pts


def currents(I, gamma_deg):
    """(i_d, i_q) for a current of peak I at angle gamma from +q toward -d."""
    g = math.radians(gamma_deg)
    return -I * math.sin(g), I * math.cos(g)


def three_terms(pp, psi_pm, psi_d, psi_q, i_d, i_q):
    """The d-q torque T = 1.5 p (psi_d i_q - psi_q i_d), split into three products of a field and a current [N.m].

    psi_d = psi_pm + (the d flux the current adds), so 1.5 p psi_d i_q splits into the magnets' part and the rest."""
    k = 1.5 * pp
    return {"magnets x i_q": k * psi_pm * i_q,                 # the magnets' field, pushed by i_q
            "q field x i_d": -k * psi_q * i_d,                 # the iron's q field (made by i_q), pushed by i_d
            "d field of i_d x i_q": k * (psi_d - psi_pm) * i_q}  # i_d's own d field (negative i_d: it weakens d), pushed by i_q


def harmonic_amplitudes(x):
    """Peak amplitude of each harmonic of one period sampled evenly: index 1 is the fundamental."""
    X = np.fft.rfft(np.asarray(x, dtype=float))
    a = 2 * np.abs(X) / len(x)
    a[0] /= 2
    return a
