"""Pack the flux-explorer sweep into one self-contained HTML page: flux_explorer_template.html + the data, embedded.

Per rotor (17): the coarse FEMM mesh, and the vector potential A on every node for 2 bands x 9 currents x 5 rotor
positions, at 12-bit precision, gzipped and base64-encoded. Torque and flux linkages at every solved point go into a
JSON block. Run, from this folder, after flux_explorer_sweep.py and flux_explorer_sweep120.py:
    uv run python flux_explorer_build.py
Writes pages/Rotor Flux Explorer.html, next to this script.

The footer of the page says where the pictures come from: which FEMM, how many solves, which day. None of that is
typed in here. provenance() reads the version from the installed FEMM and the date from the sweeps' own tables, and
the rest (solves, nodes, turns, the steps of the sweep) is counted from the data and from the sweep's constants.
"""
import base64
import datetime
import functools
import gzip
import json
import math
import os

import numpy as np
import pandas as pd

import femm_session
import motorgeo as mg
import paths
from flux_explorer_sweep import BANDS, BASE, FIELD_G, FIELD_I, OUT, POS, SCALAR_G, SLIDERS, field_states, geometries

LABELS = {"mag_y": ("Distance from the pole centre line", "mm"),
          "mag_x": ("Height along the pole centre line", "mm"),
          "v_half_deg": ("V half-angle (each magnet turns)", "°"),
          "web": ("Steel web under the magnet", "mm")}
GAMMAS = [0, 15, 30, 45, 60, 75, 90]
HERE = os.path.dirname(os.path.abspath(__file__))


def flat(edges, chord=0.3):
    """Edges -> one list of points; arcs cut into ~0.3 mm chords (enough for drawing)."""
    pts = []
    for e in edges:
        if isinstance(e, mg.Arc):
            n = max(2, math.ceil(math.radians(e.sweep()) * e.radius / chord))
            pts += e.points(step_deg=e.sweep() / n)
        else:
            pts += e.points()
    return [[round(x, 3), round(y, 3)] for x, y in pts]


def arc_pts(r, a0, a1, step=0.5):
    return [[round(r * math.cos(math.radians(a)), 3), round(r * math.sin(math.radians(a)), 3)]
            for a in np.arange(a0, a1 + 1e-9, step)]


def provenance(out_dir):
    """What the footer says made the pictures, and when: the engine ("FEMM 5.0.20261101") and the day (YYYY-MM-DD).

    Both come from facts, not from typing. The date is the day the sweeps finished writing their tables
    (states.csv and states_p120.csv in out_dir). The version is what the installed FEMM states about itself
    (femm_session.femm_version()), and it is named only if the data cannot be older than that FEMM: tables written
    before the FEMM program file's own date came from an earlier FEMM, so the footer then says just "FEMM" instead of
    naming a version that did not make them. Build the page right after the sweeps, with the FEMM that ran them, and
    the footer names that FEMM.
    """
    written = [os.path.getmtime(os.path.join(out_dir, name)) for name in ("states.csv", "states_p120.csv")]
    day = datetime.date.fromtimestamp(max(written)).isoformat()
    try:
        version, built = femm_session.femm_version(), os.path.getmtime(femm_session.femm_program())
    except FileNotFoundError:                          # FEMM is not installed here, so its version cannot be known
        return "FEMM", day
    if version is None or min(written) < built:
        why = "the sweep tables are older than the installed FEMM" if version else "the installed FEMM states none"
        print(f"note: the footer does not name a FEMM version ({why})")
        return "FEMM", day
    return f"FEMM {version}", day


if __name__ == "__main__":
    geoms = geometries()
    p0 = mg.Params()
    states = pd.read_csv(os.path.join(OUT, "states.csv"))

    # ---- fields: load, check that every file of a rotor has the same mesh, stack in the page's order
    fields, meshes, amax = {}, {}, 0.0
    for gid, _ in geoms:
        stack, ref = [], None
        for band in BANDS:
            main = np.load(os.path.join(OUT, f"g{gid:02d}_{band}.npz"))
            p120 = np.load(os.path.join(OUT, f"g{gid:02d}_{band}_p120.npz"))
            if ref is None:
                ref = (main["xy"], main["tri"])
            assert np.array_equal(main["xy"], ref[0]) and np.array_equal(main["tri"], ref[1]), f"mesh differs, g{gid} {band}"
            assert np.array_equal(p120["xy"], ref[0]), f"120-deg mesh differs, g{gid} {band}"
            A = main["A"].reshape(len(field_states()), len(POS), -1)               # [current][position][node]
            A = np.concatenate([A, p120["A"][:, None, :]], axis=1)                 # + the 120 deg position
            stack.append(A)
        fields[gid] = np.stack(stack)                                               # [band][current][position][node]
        meshes[gid] = ref
        amax = max(amax, float(np.abs(fields[gid]).max()))
    print(f"{len(geoms)} rotors, A max {amax:.5f} Wb/m, nodes {[m[0].shape[0] for m in meshes.values()][:3]}...")

    # ---- binary blobs, packed so gzip can do its job (2026-10-28 test: 1.81 -> 0.77 MB per rotor):
    #   nodes sorted ring by ring (neighbours in the list are neighbours on the drawing), coordinates in whole micrometres,
    #   A at 12-bit precision (|B| error ~0.01-0.02 T), stored as the difference to the previous node, low bytes then
    #   high bytes. Layout: [N, M, states] uint32 | x, y uint16 (um) | triangles uint16 | low bytes | high bytes
    blobs, raw_total = [], 0
    for gid, _ in geoms:
        xy, tri = meshes[gid]
        r, ang = np.hypot(xy[:, 0], xy[:, 1]), np.arctan2(xy[:, 1], xy[:, 0])
        order = np.lexsort((ang, np.round(r / 0.25)))
        new_index = np.empty_like(order); new_index[order] = np.arange(len(order))
        A = fields[gid].reshape(-1, xy.shape[0])[:, order]
        q = np.round(A / amax * 2047).astype(np.int32)
        d = np.diff(q, axis=1, prepend=0).astype("<i2").ravel().view(np.uint8).reshape(-1, 2)
        um = np.round(xy[order] * 1000)
        assert um.min() >= 0 and um.max() < 65535, "coordinates must fit uint16 micrometres"
        raw = (np.array([xy.shape[0], tri.shape[0], A.shape[0]], "<u4").tobytes() + um.astype("<u2").tobytes()
               + new_index[tri].astype("<u2").tobytes() + d[:, 0].tobytes() + d[:, 1].tobytes())
        raw_total += len(raw)
        blobs.append((gid, base64.b64encode(gzip.compress(raw, 9)).decode()))

    # ---- outlines: rotor per rotor, stator once
    rotor = {}
    for gid, x in geoms:
        regions = mg.build(mg.Params(**x))
        rr = {"magnets": [], "bands": [], "od": arc_pts(p0.r_ro, 0, 90)}
        for r in regions:
            if r.name.startswith("magnet_"):
                rr["magnets"].append({"pts": flat(r.edges), "lab": [round(v, 3) for v in r.label], "dir": round(r.magdir, 3)})
            elif r.name.startswith("barrier_"):
                rr["bands"].append(flat(r.edges))
        rotor[gid] = rr
    coils, edges = [], [arc_pts(p0.r_si, 0, 90), arc_pts(p0.r_so, 0, 90)]
    for r in mg.build(p0):
        if r.name.startswith("coil_"):
            coils.append({"pts": flat(r.edges), "lab": [round(v, 3) for v in r.label], "phase": r.circuit,
                          "sign": 1 if r.turns > 0 else -1})
        elif r.name.startswith("opening_"):
            edges.append(flat(r.edges) + [flat(r.edges)[0]])

    # ---- solved torque and flux linkages: [I level 0/10/20][gamma 0..90 step 15][position 0/30/60/90]
    scal = {}
    for gid, _ in geoms:
        scal[gid] = []
        for band in BANDS:
            s = states[(states.gid == gid) & (states.band == band)]
            out = {}
            for name, col in (("T", "T_age"), ("psi_d", "psi_d"), ("psi_q", "psi_q")):
                lev = []
                for I in (0.0, *FIELD_I):
                    rows = []
                    for g in GAMMAS:
                        sel = s[(s.I == I) & ((s.gamma == g) if I > 0 else (s.gamma == 0))].set_index("theta_e")[col]
                        rows.append([round(float(sel[th]), 6) for th in POS])
                    lev.append(rows)
                out[name] = lev
            scal[gid].append(out)

    gid_of = {}
    for gid, x in geoms[1:]:
        key = next(k for k in SLIDERS if x[k] != BASE[k])
        gid_of[f"{key}:{SLIDERS[key].index(x[key])}"] = gid
    axes, _ = mg.phase_axes(p0)
    n_solves = len(states) + len(pd.read_csv(os.path.join(OUT, "states_p120.csv")))
    engine, day = provenance(OUT)
    step_g = functools.reduce(math.gcd, [int(g) for g in FIELD_G + SCALAR_G])      # the finest step of the current angle
    meta = {
        "sliderOrder": list(SLIDERS),
        "sliders": {k: {"label": LABELS[k][0], "unit": LABELS[k][1], "stops": [round(v, 2) for v in SLIDERS[k]],
                        "base": SLIDERS[k].index(BASE[k])} for k in SLIDERS},
        "gidOf": gid_of, "bands": len(BANDS), "nCur": len(field_states()), "nPos": len(POS) + 1,
        "aScale": amax, "aSpan": 2 * amax, "axisA": round(axes["A"], 4), "polePairs": p0.pole_pairs,
        "rSplit": p0.r_ro + p0.airgap / 2, "turns": p0.turns,
        "rotor": rotor, "stator": {"coils": coils, "edges": edges}, "scal": scal,
        "footer": [
            f"Every picture comes from {engine} solves of the first-pass rotor ({n_solves:,} solves, {day}): 90° periodic slice, "
            f"coarse mesh ({meshes[0][0].shape[0]:,} nodes), N45 magnets as drawn in the hand-drawn .FEM, M-19 steel, "
            f"{p0.turns} turns per coil side (placeholder), {p0.stack:g} mm deep.",
            f"Between solved points the page blends the fields: current every {FIELD_I[1] - FIELD_I[0]:g} A and "
            f"{FIELD_G[1] - FIELD_G[0]:g}° of current angle, rotor every {POS[1] - POS[0]:g}° electrical. "
            "The whole picture repeats every 120° electrical (rotated by 30° mechanical), so a third of a turn "
            f"was solved. Torque read-outs use solved values at {step_g}° steps of current angle; 'torque now' is sampled "
            f"every {POS[1] - POS[0]:g}° of rotor position, so it misses the fastest slot ripple.",
            f"Sliders change one thing at a time from the first-pass design. Code: {paths.P} "      # the folder this script is in (paths.py)
            "(flux_explorer_sweep.py, flux_explorer_build.py). Notes: the air band study in the project notes.",
        ],
    }
    with open(os.path.join(HERE, "flux_explorer_template.html"), encoding="utf-8") as f:
        page = f.read()
    assert page.count("__META__") == 1 and page.count("__BLOBS__") == 1, "each placeholder must appear once in the template"
    page = page.replace("__META__",json.dumps(meta, separators=(",", ":"), ensure_ascii=False))
    page = page.replace("__BLOBS__", "\n".join(f'<script type="application/octet-stream" id="g{gid}">{b}</script>'
                                               for gid, b in blobs))
    path = os.path.join(HERE, "pages", "Rotor Flux Explorer.html")     # tracked in git, as a binary file
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(page)
    print(f"wrote {path}: {os.path.getsize(path) / 1e6:.1f} MB (fields {raw_total / 1e6:.1f} MB raw, "
          f"{sum(len(b) for _, b in blobs) / 1e6:.1f} MB gzipped+base64)")
