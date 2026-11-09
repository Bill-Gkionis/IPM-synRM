"""B1: Python -> FEMM 5.0 (femm_pool) on the motor's own geometry: no-load flux linkage, MTPA torque and cogging torque of the no-band sector.

What it shows, and why
    FEMM is the solver the whole motor study stands on, and Task 4 proved it on a motor at one operating point. B1 proves it the way the
    pipeline will use it: a pool of eight FEMMs, each with its own model file, solving hundreds of rotor positions and current angles of
    the band-free motor (R11), and the numbers that come back are the three the design is judged on:
      * psi_PM, the flux the magnets put through the winding, from the no-load sweep (the fundamental of the phase-A flux linkage);
      * the MTPA torque, the largest torque the machine makes from 10 A (peak) by choosing the current angle;
      * the cogging torque, the ripple the magnets and the slots make with no current at all, over its one period.
    Each is held against the README's band-free row (barrier_study.py, FEMM 5.0, fine mesh): 16.61 mWb and 1.127 N.m, within 0.5 %.

Two geometries, because the row is not the bench's geometry
    barrier_study.py's VARIANTS are {"with band": "air", "no band": "steel"}. Its "no band" row is barrier="steel": the band is still
    DRAWN (its outline is in the drawing, so the mesh has to follow it) and filled with the same M-19 as the rotor. The bench's geometry,
    barrier="none" (R11), does not draw the band at all. Same steel, same magnets, same physics; a different mesh, because FEMM meshes
    what is drawn: barrier="steel" has four 0.5 mm blocks in the rotor steel of the sector, barrier="none" does not. same_machine()
    proves the first half (no solver) and the node counts show the second. So B1 solves BOTH:
      steel  must REPRODUCE the row (same model, same mesh, same positions as the study: here the comparison is with itself);
      none   is the bench's number: the 0.5 % test of the brief applies to it. If it misses, the case shows it is the mesh (the node
             counts, and a finer-mesh point where every block's mesh size is halved) and reports it as missed.
    A third variant, air (the design WITH its band), is solved for the cogging sweep only: the README's cogging number (25.7 mN.m, the one
    "of the design") is for this variant, and the case checks that it is (the band-free variants give about twice as much) instead of assuming it.

Method (every solve is on a private FEMM, in a pool, under the solver slot)
    Positions and angles are barrier_study.py's, so the numbers are comparable (tests/test_bench_b1.py pins them against that file):
      no load      theta_e = 0, 3, ..., 357 (120 positions over one electrical period): psi_PM = the fundamental of phase A's flux
                   linkage (an FFT of the 120 samples), mean psi_d, and the spread of |psi_dq| over the positions;
      MTPA         I = 10 A peak (x 10 turns = 100 A-turns), current angle gamma from +q towards -d, id = -I sin(gamma), iq = I cos(gamma).
                   Coarse: gamma = 0, 7.5, ..., 90, then 1-degree steps over +-4 degrees around the coarse peak. At every gamma the torque is
                   the mean over six rotor positions (theta_e = 0, 10, ..., 50), which cancels the 6th to 30th harmonics of the ripple and
                   the cogging fundamental. The peak is a parabola through the best sample and its two neighbours (the study's estimator),
                   applied to the 1-degree samples (the bench's MTPA) and, to show the row is reproduced, to the coarse samples (the row's);
      cogging      no current, theta_e = 0 ... 40 in 0.5-degree steps (81 positions). The period is 5 degrees mechanical because the slot
                   pattern (36) and the pole pattern (8) repeat together every 360/lcm(36, 8) = 5 degrees = 20 electrical degrees at 4 pole
                   pairs. Peak to peak over the first period; the second period must repeat the first;
      drawn        one no-load solve at the position the geometry is drawn in (the rotor angle of the model is 0), with 360 samples of the
                   air-gap flux density on the mid-gap circle (stator frame, 0.25 degrees apart). B3 holds Elmer against it; the case's
                   figure is this solution.
    A model is one .fem per chunk of jobs (one variant per chunk, fine mesh = FEMM's smart mesh on), with a file name of its own, so no two
    FEMMs ever share a file. Python hands a FEMM its model by saving it (motorgeo.to_femm writes the .fem through the FEMM that
    femm_session started) and reads the answer through the binding's calls (mo_getcircuitproperties for the flux linkages, mo_gapintegral
    for the torque, mo_getgapb for the gap field): FEMM's own TCP protocol, one line of text per command, on a port the pool's workers each
    asked the operating system for (never 48015, never femm.openfemm()).

Run (from the Python to FEMM folder; it takes the solver slot and about two minutes)
    uv run python bench/run_bench.py --only B1
Writes out\\bench\\B1\\: results.json (every number and every check), jobs.csv (every solve), figure.png, curves.png, run.log, runs.json.

Numbers to expect (2026-11-06, FEMM 5.0.20261101, 8 workers; the whole case is 65 s, 117 s with the finer-mesh point; every number is in results.json)
    meshes (nodes of the 90-degree sector, FEMM's fine mesh)   none 15,859 (30,610 triangles), steel 17,934 (34,758); blocks halved: 23,584 and 27,199
    steel   psi_PM 16.6073 mWb (-0.016 % against the README's 16.61), MTPA 1.12745 N.m (+0.040 % against 1.127; 1.12747 from the 1-degree steps),
            cogging 52.43 mN.m; psi_PM, MTPA and cogging equal out\\barrier\\summary.csv's no-band numbers to 0.000 % (that file is git-ignored: on a PC that
            never ran barrier_study.py these three checks are in the list all the same, SKIPPED with the reason and not gating)
    none    psi_PM 16.629 mWb (+0.114 %), MTPA 1.12914 N.m at 22.9 degrees (+0.190 %), cogging 50.31 mN.m peak to peak
    air     cogging 25.68 mN.m (the README's 25.7 is this variant: -0.08 %)
    both    the cogging wave repeats (the 5 % limit is the bench's own): the second period differs by 1.1 % (none) and 3.0 % (steel) of the peak to peak; the torque from the air-gap element and
            from the weighted stress tensor agree to 0.025 % (none) and 0.017 % (steel)
    mesh    none against steel on halved blocks: psi_PM +0.131 % -> +0.067 %, torque +0.149 % -> +0.059 %, cogging -4.04 % -> +1.09 %: the difference between
            the two geometries is the mesh
    time    641 + 108 solves in 63 s on 8 workers (stage 1: 46 tasks of 14 solves, 54 s; stage 2: 8 tasks, 9 s); the finer-mesh point 336 solves, 53 s

Self-check (the answers follow).
  1. The brief asks for ten amperes "peak", and 10 A x 10 turns is called 100 A-turns. Where does the 10 turns come from, and why is the
     torque at "10 A" not simply 100 A-turns x something?
  2. Why is the torque at each current angle the MEAN over six rotor positions, and not the torque at one?
  3. psi_PM is the fundamental of an FFT of 120 samples of one period. What would a window of 100 samples that does not cover a whole
     period do to it?
  4. barrier="steel" and barrier="none" have the same steel and the same magnets. Why can their psi_PM still differ in the third digit?
Answers.
  1. motorgeo.Params.turns = 10 per coil side, as in the hand-drawn .FEM (a placeholder, README "Placeholders"). The torque is not linear in
     the current: the steel saturates (the bridges first), and the reluctance part goes with the product id x iq. The MTPA angle moves from
     about 14 degrees at 5 A to about 31 degrees at 20 A (README, Air band study).
  2. The torque of an IPM machine ripples with the rotor position (slotting, harmonics of the winding and of the magnets): one position can
     be 10 % above or below the mean. Six positions 10 electrical degrees apart cancel the ripple's 6th, 12th, 18th, 24th and 30th
     harmonics, and the cogging's fundamental (period 20 degrees), which is why a smooth curve over gamma can be fitted at all.
  3. It would not be a fundamental. The DFT assumes the samples are one period repeated; with a partial period the energy of the
     fundamental leaks into the neighbouring frequencies, and the amplitude is off by a few percent or more.
  4. FEMM meshes the drawing. "steel" has the band's outline in it, which forces triangle edges along it and puts 0.5 mm blocks into the
     steel; "none" meshes the same steel with its own sizes. A first-order triangle mesh has a discretisation error of its own, a
     fraction of a percent here, and the two meshes make different errors. That is a mesh effect and not a physical one; the case
     shows it by halving the block sizes of both and watching the difference fall.
"""
from __future__ import annotations

import csv
import math
import os
import time
from pathlib import Path

import numpy as np

import femm_queue
import motorgeo as mg
import paths

from . import OWN_LIMIT, Log, check, compare, fig_field, params, skipped

CASE = "B1"
TITLE = "Python -> FEMM 5.0 (femm_pool): the no-band sector, flux linkage, MTPA torque, cogging"

# The yardsticks. The expectation is the README's band-free row ("Expected numbers", barrier_study.py, FEMM 5.0, fine mesh).
ROW = {"psi_pm_mWb": 16.61, "T_mtpa_Nm": 1.127}
ROW_SOURCE = "README.md, Expected numbers: barrier_study.py's band-free row (FEMM 5.0, fine mesh): 16.61 mWb, 1.127 N.m"
README_COGGING_WITH_BAND_MNM = 25.7                   # README: "The cogging torque of the design, peak to peak, is 25.7 mN.m": the design WITH its band
TOLERANCE_PERCENT = 0.5                               # the brief's: within 0.5 % of the row
SAME_MODEL_TOLERANCE_PERCENT = 0.01                   # the same model on the same FEMM, against its own earlier numbers (out\barrier\summary.csv)
TORQUE_METHODS_TOLERANCE_PERCENT = 0.1                # T from the air-gap element and from the weighted stress tensor (femm_check's own limit)
REPEAT_TOLERANCE_PERCENT_OF_PP = 5.0                  # the second cogging period against the first, as a share of the peak to peak (set before the run)

# barrier_study.py's own constants, repeated here and pinned to that file by tests/test_bench_b1.py (the study is a script: not imported)
VARIANTS = {"none": "none", "steel": "steel", "air": "air"}      # B1's name -> motorgeo's barrier: none = R11's geometry, steel = the row's
I_RATED = 10                                          # A peak; x 10 turns per coil side = 100 A-turns
GAMMAS = np.arange(0, 90.1, 7.5)                      # current angle, from +q towards -d (0 = all q current)
MEAN_POS = range(0, 60, 10)                           # six rotor positions: their mean cancels the 6th ... 30th ripple harmonics
GAP_DEG = list(np.arange(0, 90, 0.25))                # air-gap samples, stator frame, one 90-degree sector
FINE_HALF_WIDTH = 4                                   # the 1-degree steps cover the coarse peak +- 4 degrees
COGGING_STEP_EL = 0.5                                 # electrical degrees
CHUNK = 14                                            # jobs per task of the pool (one model per task)
TASK_TIMEOUT_S = 300                                  # a chunk is about 15 s of solves; a hung FEMM is ended after 5 minutes and the chunk tried again
MESH_SCALES = (0.5,)                                  # the finer-mesh point: every block's mesh size times this
CELL = {"none": "barrier=\"none\" (R11)", "steel": "barrier=\"steel\" (the row's)", "air": "barrier=\"air\" (the design)"}


# ============================================================================================== the arithmetic (no solver)
def noload_positions() -> list[float]:
    """theta_e = 0, 3, ..., 357: 120 positions over one electrical period (barrier_study.py's no-load sweep)."""
    return [float(t) for t in range(0, 360, 3)]


def coarse_gammas() -> list[float]:
    """gamma = 0, 7.5, ..., 90 degrees: 13 current angles."""
    return [float(g) for g in GAMMAS]


def mean_positions() -> list[float]:
    """The six rotor positions (electrical degrees) whose mean torque is the torque at one current angle."""
    return [float(t) for t in MEAN_POS]


def cogging_period_mech(p: mg.Params) -> float:
    """The cogging torque repeats every 360 / lcm(slots, poles) mechanical degrees: 5 for 36 slots and 8 poles (lcm 72)."""
    return 360.0 / math.lcm(p.slots, 2 * p.pole_pairs)


def cogging_period_el(p: mg.Params) -> float:
    """The same period in electrical degrees: 20 for 4 pole pairs."""
    return cogging_period_mech(p) * p.pole_pairs


def cogging_positions(periods: int = 2) -> list[float]:
    """theta_e = 0, 0.5, ..., 20 * periods (81 positions for two periods): the first period for the number, the second to see it repeat."""
    steps = int(round(periods * cogging_period_el(mg.Params()) / COGGING_STEP_EL))
    return [COGGING_STEP_EL * k for k in range(steps + 1)]


def currents(gamma_deg: float, amps: float) -> tuple[float, float]:
    """(id, iq) for a current of `amps` at the angle gamma from the q axis towards -d: id = -I sin(gamma), iq = I cos(gamma)."""
    g = math.radians(gamma_deg)
    return -amps * math.sin(g), amps * math.cos(g)


def harmonic_amplitudes(samples) -> np.ndarray:
    """Amplitudes of harmonics 0 ... N/2 of one period sampled evenly (the last point not repeated): barrier_plots.py's harmonics()."""
    c = np.fft.rfft(np.asarray(samples, float)) / len(samples)
    c[1:] *= 2
    return np.abs(c)


def fundamental(samples) -> float:
    """The amplitude of the fundamental of one sampled period (barrier_plots.py: how the study's psi_PM is made)."""
    return float(harmonic_amplitudes(samples)[1])


def parabola_peak(x, y) -> tuple[float, float]:
    """(x, y) of the top of a parabola through the best sample and its two neighbours (barrier_study.py's mtpa(): the study's estimator).
    A best sample at either end takes its neighbours from inside. At least three samples."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) < 3:
        raise ValueError("a parabola needs at least three samples")
    k = min(max(int(np.argmax(y)), 1), len(x) - 2)
    a, b, c = np.polyfit(x[k - 1:k + 2], y[k - 1:k + 2], 2)
    return float(-b / (2 * a)), float(c - b * b / (4 * a))


def refine_gammas(best: float, half_width: int = FINE_HALF_WIDTH) -> list[float]:
    """The 1-degree steps around a coarse peak: the nearest whole degree +- half_width, kept inside 0 ... 90."""
    centre = int(math.floor(best + 0.5))
    return [float(g) for g in range(centre - half_width, centre + half_width + 1) if 0 <= g <= 90]


def drawn_theta_e(p: mg.Params, axis_a: float) -> float:
    """The electrical angle at which FemmModel.solve_abc turns the air-gap element's rotor side by 0: the rotor is where it is DRAWN, and
    the rotor frame and the stator frame are one frame. rotor_deg = theta_e / pole_pairs + axis_a - theta_d0 = 0, theta_d0 = 90 / pole_pairs."""
    return p.pole_pairs * (90 / p.pole_pairs - axis_a)


def stage_one_jobs() -> list[dict]:
    """Every solve of the first stage, as plain dicts that can cross a pipe: study, variant, theta_e, id, iq (and I, gamma, gap_deg, keep).

    none and steel: the no-load sweep, the coarse current-angle sweep (13 angles x 6 positions), the cogging sweep, the drawn position
    (with the air-gap samples, and its .ans kept). air: the cogging sweep only."""
    p = mg.Params()
    drawn = drawn_theta_e(p, mg.phase_axes(p)[0]["A"])
    jobs: list[dict] = []
    for variant in ("none", "steel"):
        jobs += [{"variant": variant, "study": "noload", "theta_e": t, "id": 0.0, "iq": 0.0} for t in noload_positions()]
        for g in coarse_gammas():
            i_d, i_q = currents(g, I_RATED)
            jobs += [{"variant": variant, "study": "gamma", "theta_e": t, "I": float(I_RATED), "gamma": g, "id": i_d, "iq": i_q}
                     for t in mean_positions()]
        jobs += [{"variant": variant, "study": "cogging", "theta_e": t, "id": 0.0, "iq": 0.0} for t in cogging_positions()]
        jobs.append({"variant": variant, "study": "drawn", "theta_e": drawn, "id": 0.0, "iq": 0.0,
                     "gap_deg": [float(a) for a in GAP_DEG], "keep": True})
    jobs += [{"variant": "air", "study": "cogging", "theta_e": t, "id": 0.0, "iq": 0.0} for t in cogging_positions()]
    return jobs


def fine_jobs(variant: str, best_gamma: float, skip: set[float]) -> list[dict]:
    """The second stage for one variant: the 1-degree steps around its coarse peak (six positions each), minus the angles already solved."""
    jobs = []
    for g in refine_gammas(best_gamma):
        if g in skip:
            continue
        i_d, i_q = currents(g, I_RATED)
        jobs += [{"variant": variant, "study": "gamma", "theta_e": t, "I": float(I_RATED), "gamma": g, "id": i_d, "iq": i_q}
                 for t in mean_positions()]
    return jobs


def make_chunks(jobs: list[dict], size: int, mesh_scale: float, fem_dir, first_id: int = 0) -> list[dict]:
    """Cut the jobs into tasks of at most `size`, one variant in each (a task builds ONE model and solves its jobs on it).
    Returns [{"id", "variant", "jobs", "mesh_scale", "fem_dir"}], ids unique, counted from first_id (so that the model files of two stages never
    share a name)."""
    chunks: list[dict] = []
    for variant in dict.fromkeys(j["variant"] for j in jobs):
        mine = [j for j in jobs if j["variant"] == variant]
        for k in range(0, len(mine), size):
            chunks.append({"id": first_id + len(chunks), "variant": variant, "jobs": mine[k:k + size], "mesh_scale": mesh_scale, "fem_dir": str(fem_dir)})
    return chunks


def same_machine() -> dict:
    """The check that comes before any solve: is barrier="none" the geometry of barrier_study.py's band-free row? It is not; this says exactly what
    differs. (No solver: motorgeo's regions and shapely areas.) The two builds have the same regions except the band outlines (four in the 90-degree
    sector); the rotor steel of "none" is the rotor steel of "steel" plus the bands filled with steel, to the last square micrometre; every other
    region has the same outline, label and mesh size. What differs is only what FEMM is asked to mesh: the band outlines and their 0.5 mm blocks."""
    none, steel = mg.build(mg.Params(barrier="none")), mg.build(mg.Params(barrier="steel"))
    by_none, by_steel = {r.name: r for r in none}, {r.name: r for r in steel}
    vis_none, vis_steel = mg.visible(none), mg.visible(steel)
    bands = sorted(n for n in by_steel if n.startswith("barrier_"))
    shared = sorted(set(by_none) & set(by_steel))
    same = all(by_none[n].edges == by_steel[n].edges and by_none[n].label == by_steel[n].label and by_none[n].material == by_steel[n].material
               and by_none[n].mesh == by_steel[n].mesh and by_none[n].magdir == by_steel[n].magdir for n in shared)
    return {"only_in_steel": sorted(set(by_steel) - set(by_none)), "only_in_none": sorted(set(by_none) - set(by_steel)),
            "n_band_regions": len(bands), "band_material": {by_steel[n].material for n in bands},
            "steel_area_none_mm2": vis_none["rotor_steel"].area, "steel_area_steel_mm2": vis_steel["rotor_steel"].area,
            "band_area_mm2": sum(vis_steel[n].area for n in bands),
            "magnet_area_none_mm2": sum(vis_none[n].area for n in by_none if n.startswith("magnet_")),
            "magnet_area_steel_mm2": sum(vis_steel[n].area for n in by_steel if n.startswith("magnet_")),
            "other_regions_identical": bool(same), "mesh_size_band_mm": by_steel[bands[0]].mesh, "mesh_size_rotor_steel_mm": by_steel["rotor_steel"].mesh}


# ============================================================================================== the worker (runs in a pool worker, on its own FEMM)
def solve_chunk(chunk: dict) -> dict:
    """One task of the pool: build this chunk's model on the worker's own FEMM, solve every job of the chunk on it, return the rows.

    The model file has a name of its own (variant, mesh scale, chunk id), so no two FEMMs share a file. Rows carry what barrier_study.py's rows do
    (flux linkages, three torques, the rotor angle, the seconds) and, for a job that asked for it, 360 samples of the air-gap flux density and the
    node and element counts of the mesh. The worker's facts (pid, FEMM's pid, both priority classes) come back too: they are the log of B9."""
    import femm
    import femm_session

    from . import priority_class
    variant, scale = chunk["variant"], chunk["mesh_scale"]
    p = mg.Params(barrier=VARIANTS[variant])
    regions = mg.build(p)
    if scale != 1.0:
        for r in regions:                                  # the finer-mesh point: every block's mesh size scaled (the Region is plain data)
            r.mesh *= scale
    fem_dir = Path(chunk["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{variant}_s{scale:g}_c{chunk['id']:03d}.fem")
    started = time.perf_counter()
    mg.to_femm(p, regions, fem, fine=True)                 # fine = FEMM's smart mesh on, as in the numbers of the README
    femm.mi_close()
    model = mg.FemmModel(p, fem)
    build_s = time.perf_counter() - started
    rows, meshes = [], []
    try:
        for job in chunk["jobs"]:
            t = time.perf_counter()
            ans = str(fem_dir / f"{variant}_s{scale:g}_{job['study']}.ans") if job.get("keep") else None
            r = model.solve(job["theta_e"], job["id"], job["iq"], gap_deg=job.get("gap_deg", ()), keep_ans=ans)
            row = {"study": job["study"], "variant": variant, "mesh_scale": scale, "theta_e": job["theta_e"], "I": job.get("I", 0.0),
                   "gamma": job.get("gamma"), "id": r["id"], "iq": r["iq"], "psi_a": r["psi_a"], "psi_b": r["psi_b"], "psi_c": r["psi_c"],
                   "psi_d": r["psi_d"], "psi_q": r["psi_q"], "T_age": r["T_age"], "T_wst": r["T_wst"], "T_dq": r["T_dq"],
                   "rotor_deg": r["rotor_deg"], "psi_retries": r["psi_retries"], "seconds": time.perf_counter() - t}
            if job.get("gap_deg"):
                row["gap_deg"], row["gap_br"], row["gap_bt"] = job["gap_deg"], r["gap_br"], r["gap_bt"]
            if ans:                                        # how big is the mesh? ask FEMM, which still has the solution on disk
                femm.mi_loadsolution()
                row["nodes"], row["elements"] = int(femm.mo_numnodes()), int(femm.mo_numelements())
                femm.mo_close()
                row["ans"] = ans
                meshes.append({"variant": variant, "mesh_scale": scale, "nodes": row["nodes"], "elements": row["elements"]})
            rows.append(row)
    finally:
        femm.mi_close()                                    # the pool's FEMM lives on for the next task: leave no document open in it
    me, his = os.getpid(), femm_session.pid(femm._default)
    return {"chunk": chunk["id"], "variant": variant, "mesh_scale": scale, "rows": rows, "meshes": meshes, "build_seconds": build_s,
            "worker": {"pid": me, "class": priority_class(me), "femm_pid": his, "femm_class": priority_class(his) if his else None}}


# ============================================================================================== one stage on the pool
def run_stage(pool, jobs: list[dict], fem_dir, label: str, log: Log, mesh_scale: float = 1.0, first_id: int = 0) -> list[dict]:
    """Solve `jobs` on the pool in chunks. Returns the chunk results. A chunk that timed out twice or whose worker died is an error: a missing
    solve would silently change a mean or a fundamental, so nothing is returned instead of a table with holes."""
    chunks = make_chunks(jobs, CHUNK, mesh_scale, fem_dir, first_id)
    started = time.perf_counter()
    results = pool.map(solve_chunk, chunks)
    lost = [(c["id"], r) for c, r in zip(chunks, results) if not isinstance(r, dict) or "rows" not in r]       # (the chunk's id, what the pool said): the results go by position
    if lost:
        raise RuntimeError(f"{label}: chunk(s) {[i for i, _ in lost]} did not come back ({[said for _, said in lost]}): the pool ended their workers; see run.log")
    log.say(f"B1 {label}: {len(jobs)} solves in {len(chunks)} chunks, {time.perf_counter() - started:.1f} s")
    return results


def table(results: list[dict]) -> list[dict]:
    """The rows of every chunk, in one list."""
    return [row for r in results for row in r["rows"]]


# ============================================================================================== reading the solves
def psi_pm_of(rows: list[dict]) -> dict:
    """psi_PM and the no-load numbers of one variant from its no-load rows (any order)."""
    nl = sorted((r for r in rows if r["study"] == "noload"), key=lambda r: r["theta_e"])
    mag = np.hypot([r["psi_d"] for r in nl], [r["psi_q"] for r in nl]) * 1e3
    return {"n": len(nl), "psi_pm_mWb": fundamental([r["psi_a"] for r in nl]) * 1e3, "psi_d0_mWb": float(np.mean([r["psi_d"] for r in nl]) * 1e3),
            "psi_dq_min_mWb": float(mag.min()), "psi_dq_max_mWb": float(mag.max()), "psi_dq_mean_mWb": float(mag.mean())}


def torque_curve(rows: list[dict], i_rated: float = I_RATED) -> dict:
    """{gamma: (mean T_age, mean T_wst, n positions)} over the six rotor positions, from the current-angle rows at i_rated."""
    out: dict[float, list] = {}
    for r in rows:
        if r["study"] == "gamma" and r["I"] == float(i_rated):
            out.setdefault(r["gamma"], []).append((r["T_age"], r["T_wst"]))
    return {g: (float(np.mean([a for a, _ in v])), float(np.mean([b for _, b in v])), len(v)) for g, v in sorted(out.items())}


def coarse_peak(rows: list[dict], coarse: list[float]) -> float:
    """The coarse current angle (one of `coarse`) with the largest mean torque: where the 1-degree steps of the second stage are centred."""
    curve = torque_curve(rows)
    g_c = [g for g in curve if g in set(coarse)]
    return g_c[int(np.argmax([curve[g][0] for g in g_c]))]


def mtpa_of(rows: list[dict], coarse: list[float]) -> dict:
    """The MTPA of one variant: the row's estimator on the coarse samples, and the bench's on the 1-degree samples around the coarse peak."""
    curve = torque_curve(rows)
    g_all = list(curve)
    g_c = [g for g in g_all if g in set(coarse)]
    t_c = [curve[g][0] for g in g_c]
    row_gamma, row_t = parabola_peak(g_c, t_c)
    best_coarse = g_c[int(np.argmax(t_c))]
    fine = [g for g in g_all if g == int(g) and abs(g - best_coarse) <= FINE_HALF_WIDTH + 0.5]
    fine_gamma, fine_t = parabola_peak(fine, [curve[g][0] for g in fine])
    k = int(np.argmax([curve[g][0] for g in fine]))
    return {"row_method": {"gamma_deg": row_gamma, "T_Nm": row_t, "samples": len(g_c)},
            "bench_method": {"gamma_deg": fine_gamma, "T_Nm": fine_t, "samples": len(fine), "best_sample_gamma_deg": fine[k],
                             "best_sample_T_Nm": curve[fine[k]][0], "best_sample_T_wst_Nm": curve[fine[k]][1]},
            "best_coarse_gamma_deg": best_coarse, "curve": {str(g): {"T_age_Nm": v[0], "T_wst_Nm": v[1], "positions": v[2]} for g, v in curve.items()},
            "methods_apart_percent_at_best": 100 * abs(curve[fine[k]][0] - curve[fine[k]][1]) / abs(curve[fine[k]][0])}


def cogging_of(rows: list[dict], p: mg.Params) -> dict:
    """Peak to peak over the first period, and how well the second period repeats it."""
    cg = sorted((r for r in rows if r["study"] == "cogging"), key=lambda r: r["theta_e"])
    theta = np.array([r["theta_e"] for r in cg])
    torque = np.array([r["T_age"] for r in cg])
    period = cogging_period_el(p)
    first = theta <= period + 1e-9
    pp = float(torque[first].max() - torque[first].min())
    index = {round(t, 6): i for i, t in enumerate(theta)}
    pairs = [(i, index[round(t + period, 6)]) for i, t in enumerate(theta) if t <= period + 1e-9 and round(t + period, 6) in index]
    repeat = max(abs(torque[b] - torque[a]) for a, b in pairs)
    return {"period_el_deg": period, "period_mech_deg": cogging_period_mech(p), "positions": len(cg), "positions_first_period": int(first.sum()),
            "pp_mNm": pp * 1e3, "mean_mNm": float(torque[first].mean() * 1e3), "repeat_max_difference_mNm": float(repeat * 1e3),
            "repeat_percent_of_pp": 100 * float(repeat) / pp, "theta_e": theta.tolist(), "T_Nm": torque.tolist()}


def stored_row() -> dict | None:
    """The unrounded numbers behind the README's rows: out\\barrier\\summary.csv (git-ignored, so only on the PC that ran barrier_study.py).
    None if the file is not there. {column: {name: value}} for 'with band' and 'no band'."""
    path = paths.path("out") / "barrier" / "summary.csv"
    if not path.is_file():
        return None
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    out: dict = {"with band": {}, "no band": {}}
    for row in rows:
        for column in out:
            try:
                out[column][row[""]] = float(row[column])
            except (KeyError, ValueError):
                pass
    return out


# ============================================================================================== pictures
def read_field(ans: Path) -> tuple:
    """(nodes in mm, triangles, A at the nodes, |B| of each triangle) of FEMM's own mesh, from the .ans it wrote (femfile.read_ans): what bench.fig_field draws."""
    import femfile
    sol = femfile.read_ans(str(ans))
    xy, tri, a = sol["xy"], sol["tri"], sol["A"]
    bx, by = femfile.element_b(xy, a, tri)
    return xy, tri, a, np.hypot(bx, by)


def fig_curves(summary: dict, target: Path) -> None:
    """Three panels: the no-load flux linkage, the torque against the current angle, the cogging torque over two periods."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    colours = {"none": "#1f6feb", "steel": "#e0782f", "air": "#2a9d5c"}
    for v in ("none", "steel"):
        s = summary[v]
        axes[0].plot(s["noload"]["theta_e"], s["noload"]["psi_a_mWb"], color=colours[v], label=f"{v}: psi_PM {s['psi_pm_mWb']:.2f} mWb")
        curve = s["mtpa"]["curve"]
        g = sorted(float(k) for k in curve)
        axes[1].plot(g, [curve[str(x)]["T_age_Nm"] for x in g], "o-", ms=3.5, color=colours[v],
                     label=f"{v}: MTPA {s['mtpa']['bench_method']['T_Nm']:.4f} N.m at {s['mtpa']['bench_method']['gamma_deg']:.1f} deg")
    for v in ("none", "steel", "air"):
        c = summary[v]["cogging"]
        axes[2].plot(c["theta_e"], [t * 1e3 for t in c["T_Nm"]], color=colours[v], label=f"{v}: {c['pp_mNm']:.1f} mN.m peak to peak")
    axes[0].set(xlabel="rotor position [electrical degrees]", ylabel="phase A flux linkage [mWb]", title="No load: flux linkage of phase A")
    axes[1].set(xlabel="current angle gamma [degrees]", ylabel="mean torque [N.m]", title="10 A: torque against the current angle")
    axes[2].set(xlabel="rotor position [electrical degrees] (one period = 20)", ylabel="cogging torque [mN.m]", title="No current: cogging torque, two periods")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(target, dpi=130)
    plt.close(fig)


# ============================================================================================== the case
def write_jobs_csv(path: Path, rows: list[dict]) -> None:
    """Every solve, one line each (the field samples are left out: they are in results.json for the drawn position)."""
    keep = ["study", "variant", "mesh_scale", "theta_e", "I", "gamma", "id", "iq", "psi_a", "psi_b", "psi_c", "psi_d", "psi_q", "T_age", "T_wst",
            "T_dq", "rotor_deg", "psi_retries", "seconds"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keep, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def pool_record(label: str, pool, chunk_results: list[dict], seconds: float, task_timeout_s: float | None = None) -> dict:
    """The run's record for runs.json: how many workers, and the priority class each worker and its FEMM was SEEN at (read by the worker itself). `task_timeout_s` is the limit
    the pool was opened with (B1's own TASK_TIMEOUT_S if none is given): B7 and B8 open theirs with femm_points.TASK_TIMEOUT_S and say so."""
    from . import class_name
    workers = {r["worker"]["pid"]: r["worker"] for r in chunk_results}
    classes = lambda key: _count(class_name(w[key]) if w[key] is not None else "unreadable" for w in workers.values())
    return {"kind": "femm_pool", "label": label, "workers": pool.workers, "workers_limit": paths.setting("femm", "workers"), "threads": None,
            "workers_that_solved": len(workers), "worker_pids": sorted(workers), "femm_pids": sorted(w["femm_pid"] for w in workers.values() if w["femm_pid"]),
            "worker_priority_classes_seen": classes("class"), "femm_priority_classes_seen": classes("femm_class"),
            "task_timeout_s": TASK_TIMEOUT_S if task_timeout_s is None else task_timeout_s, "tasks": len(chunk_results), "seconds": round(seconds, 2),
            "time": time.strftime("%Y-%m-%d %H:%M:%S")}


def _count(items) -> dict:
    out: dict = {}
    for item in items:
        out[item] = out.get(item, 0) + 1
    return out


def run(folder: Path, log: Log, mesh_study: bool | None = None) -> dict:
    """Run B1 in `folder` (already empty). Returns the results dict, whose "checks" list every verdict (run_bench.py writes it as results.json).
    mesh_study: also solve the finer mesh; None takes it from `run_bench.py --mesh-study` (bench.OPTIONS). It runs by itself when a number misses."""
    import femm_session

    from . import OPTIONS
    if mesh_study is None:
        mesh_study = bool(OPTIONS.get("mesh_study"))
    folder = Path(folder)
    fem_dir = folder / "fem"
    p = mg.Params()
    started = time.perf_counter()
    facts = same_machine()
    log.say(f"B1: barrier=\"none\" is not the row's geometry. Against barrier=\"steel\": regions only in steel {facts['only_in_steel']}, only in none "
            f"{facts['only_in_none']}; the band (M-19, mesh {facts['mesh_size_band_mm']} mm) is {facts['band_area_mm2']:.3f} mm2 of the sector.")
    jobs = stage_one_jobs()
    checks: list[dict] = []
    seconds: dict = {}
    with log.tee():
        with femm_queue.femm_pool("bench B1 FEMM", task_timeout_s=TASK_TIMEOUT_S) as pool:
            t0 = time.perf_counter()
            first = run_stage(pool, jobs, fem_dir, "stage 1: no load, current angles, cogging, drawn position", log)
            seconds["stage_1"] = time.perf_counter() - t0
            rows = table(first)
            coarse = coarse_gammas()
            t0 = time.perf_counter()
            fines = []
            for v in ("none", "steel"):
                fines += fine_jobs(v, coarse_peak([r for r in rows if r["variant"] == v], coarse), set(coarse))
            second = run_stage(pool, fines, fem_dir, "stage 2: 1-degree steps around the coarse peak", log, first_id=len(first))
            seconds["stage_2"] = time.perf_counter() - t0
            rows += table(second)
            chunk_results = first + second
            pool_facts = pool_record("bench B1 FEMM", pool, chunk_results, sum(seconds.values()))
            extra: list[dict] = []
            summary, checks = analyse(rows, p, facts)
            missed = [c for c in checks if c["verdict"] == "FAIL" and c["name"].startswith("none:")]
            if missed or mesh_study:
                t0 = time.perf_counter()
                log.say("B1: " + ("the bench's number misses the row" if missed else "mesh study asked for") + ": running the finer-mesh point (every block's mesh size x "
                        f"{MESH_SCALES[0]:g}) for both variants")
                extra_jobs = []
                for v in ("none", "steel"):
                    g_ref = summary[v]["mtpa"]["bench_method"]["best_sample_gamma_deg"]
                    i_d, i_q = currents(g_ref, I_RATED)
                    extra_jobs += [j for j in stage_one_jobs() if j["variant"] == v and j["study"] in ("noload", "drawn")]
                    extra_jobs += [{"variant": v, "study": "cogging", "theta_e": t, "id": 0.0, "iq": 0.0} for t in cogging_positions(periods=1)]
                    extra_jobs += [{"variant": v, "study": "gamma", "theta_e": t, "I": float(I_RATED), "gamma": g_ref, "id": i_d, "iq": i_q}
                                   for t in mean_positions()]
                for scale in MESH_SCALES:
                    res = run_stage(pool, extra_jobs, fem_dir, f"mesh study: block sizes x {scale:g}", log, mesh_scale=scale, first_id=len(chunk_results))
                    extra += table(res)
                    chunk_results += res
                seconds["mesh_study"] = time.perf_counter() - t0
                pool_facts = pool_record("bench B1 FEMM", pool, chunk_results, sum(seconds.values()))
    log.add_run(pool_facts)
    mesh_sizes = {f"{m['variant']}_s{m['mesh_scale']:g}": {"nodes": m["nodes"], "elements": m["elements"]}
                  for r in chunk_results for m in r["meshes"]}
    if extra:
        checks += mesh_study_checks(summary, rows + extra)
    results = {"case": CASE, "title": TITLE, "femm_version": femm_session.femm_version(), "params": params(p), "expected_source": ROW_SOURCE,
               "geometry": {k: (sorted(v) if isinstance(v, set) else v) for k, v in facts.items()},
               "variants": {CELL[v]: v for v in VARIANTS}, "mesh": mesh_sizes,
               "stored_row": stored_row(), "summary": summary, "mesh_study": mesh_study_summary(summary, rows + extra) if extra else None,
               "checks": checks, "pool": pool_facts, "seconds": {k: round(v, 1) for k, v in seconds.items()},
               "seconds_total": round(time.perf_counter() - started, 1)}
    drawn = next(r for r in rows if r["study"] == "drawn" and r["variant"] == "none" and r["mesh_scale"] == 1.0)
    results["drawn_reference"] = {"variant": "none", "theta_e": drawn["theta_e"], "rotor_deg": drawn["rotor_deg"], "psi_abc_Wb": [drawn["psi_a"], drawn["psi_b"], drawn["psi_c"]],
                                  "psi_d_Wb": drawn["psi_d"], "psi_q_Wb": drawn["psi_q"], "nodes": drawn["nodes"], "gap_deg": drawn["gap_deg"],
                                  "gap_br_T": drawn["gap_br"], "gap_bt_T": drawn["gap_bt"]}
    write_jobs_csv(folder / "jobs.csv", rows + extra)
    fig_field(*read_field(Path(drawn["ans"])),                # a title of short lines (bench.fig_field wraps any line over its width, but a wrap can leave an orphan word)
              f"FEMM 5.0 (femm_pool): the band-free rotor, barrier=\"none\", no load at the drawn position\n"
              f"flux density |B| and flux lines (contours of A); fine mesh, {drawn['nodes']:,} nodes in the sector\n"
              f"psi_PM {summary['none']['psi_pm_mWb']:.2f} mWb, MTPA {summary['none']['mtpa']['bench_method']['T_Nm']:.4f} N.m, "
              f"cogging {summary['none']['cogging']['pp_mNm']:.1f} mN.m peak to peak", folder / "figure.png")
    fig_curves(summary, folder / "curves.png")
    for leftover in fem_dir.glob("*.ans"):                  # the solutions are megabytes each; the case keeps the drawn position's two and nothing else
        if "_drawn" not in leftover.name:
            leftover.unlink()
    return results


def analyse(rows: list[dict], p: mg.Params, facts: dict) -> tuple[dict, list[dict]]:
    """(summary, checks) from the rows of the two stages (mesh scale 1): every number of B1 and its verdict. No solver."""
    by_variant = {v: [r for r in rows if r["variant"] == v and r["mesh_scale"] == 1.0] for v in VARIANTS}
    summary = {v: summarise(by_variant[v], p, v) for v in ("none", "steel")}
    summary["air"] = {"cogging": cogging_of(by_variant["air"], p)}
    return summary, make_checks(summary, facts)


def summarise(rows: list[dict], p: mg.Params, variant: str) -> dict:
    """Everything B1 says about one variant (none or steel) from its rows."""
    pm = psi_pm_of(rows)
    nl = sorted((r for r in rows if r["study"] == "noload"), key=lambda r: r["theta_e"])
    nodes = next((r["nodes"] for r in rows if r["study"] == "drawn"), None)
    return {"barrier": VARIANTS[variant], "nodes": nodes, "psi_pm_mWb": pm["psi_pm_mWb"], "noload": {**pm, "theta_e": [r["theta_e"] for r in nl],
            "psi_a_mWb": [r["psi_a"] * 1e3 for r in nl]}, "mtpa": mtpa_of(rows, coarse_gammas()), "cogging": cogging_of(rows, p)}


def make_checks(summary: dict, facts: dict) -> list[dict]:
    """Every number of B1 against its expectation, with the source of the expectation. The names start with the variant (none: or steel:)."""
    stored = stored_row()
    out: list[dict] = []
    out.append(check("geometry: barrier=\"none\" is the row's machine without the band outlines",
                     facts["only_in_none"] == [] and facts["n_band_regions"] == 4 and abs(facts["steel_area_none_mm2"] - facts["steel_area_steel_mm2"] - facts["band_area_mm2"]) < 1e-6
                     and facts["other_regions_identical"],
                     f"same regions except the four band outlines of barrier=\"steel\"; rotor steel(none) = rotor steel(steel) + band, to 1e-6 mm2 ({OWN_LIMIT})",
                     f"only in steel: {facts['only_in_steel']}; steel(none) - steel(steel) - band = "
                     f"{facts['steel_area_none_mm2'] - facts['steel_area_steel_mm2'] - facts['band_area_mm2']:+.2e} mm2",
                     "motorgeo.build, shapely areas (no solver)"))
    for v in ("none", "steel"):
        s = summary[v]
        out.append(compare(f"{v}: psi_PM (fundamental of phase A, no load)", s["psi_pm_mWb"], ROW["psi_pm_mWb"], TOLERANCE_PERCENT, ROW_SOURCE, "mWb"))
        method = "bench_method" if v == "none" else "row_method"
        label = "MTPA torque at 10 A (1-degree steps)" if v == "none" else "MTPA torque at 10 A (the row's estimator on the coarse samples)"
        out.append(compare(f"{v}: {label}", s["mtpa"][method]["T_Nm"], ROW["T_mtpa_Nm"], TOLERANCE_PERCENT, ROW_SOURCE, "N.m"))
        c = s["cogging"]
        out.append(check(f"{v}: cogging torque repeats every 5 degrees mechanical", c["repeat_percent_of_pp"] <= REPEAT_TOLERANCE_PERCENT_OF_PP
                         and abs(c["period_mech_deg"] - 5.0) < 1e-9 and c["positions_first_period"] >= 11,
                         f"period 360/lcm(36, 8) = 5 degrees, at least 11 positions in it (the brief); T(theta + 20 el) = T(theta) within {REPEAT_TOLERANCE_PERCENT_OF_PP:g} % "
                         f"of the peak to peak ({OWN_LIMIT})",
                         f"period {c['period_mech_deg']:g} deg mech = {c['period_el_deg']:g} el, {c['positions_first_period']} positions in it; the second period differs by at "
                         f"most {c['repeat_max_difference_mNm']:.3f} mN.m = {c['repeat_percent_of_pp']:.2f} % of the {c['pp_mNm']:.2f} mN.m peak to peak",
                         "lcm(slots, poles), and the sweep of two periods"))
        m = s["mtpa"]["bench_method"]
        out.append(compare(f"{v}: torque from the air-gap element against the weighted stress tensor, at the MTPA sample",
                           m["best_sample_T_Nm"], m["best_sample_T_wst_Nm"], TORQUE_METHODS_TOLERANCE_PERCENT,
                           f"femm_check.py: the two methods share only the field solution, so they must agree ({OWN_LIMIT})", "N.m"))
    # the same model on the same FEMM against the study's own unrounded numbers: the three checks are in the list on every PC. out\barrier\summary.csv is
    # git-ignored, so where it is not there they are SKIPPED (never gating, with the reason) and not silently left out
    s = summary["steel"]
    for name, measured, column, unit, source in (
            ("steel: psi_PM against the study's own unrounded number", s["psi_pm_mWb"], "psi_pm_mWb", "mWb",
             f"out\\barrier\\summary.csv, column 'no band' (the same model, the same positions, the same FEMM); tolerance {SAME_MODEL_TOLERANCE_PERCENT:g} % ({OWN_LIMIT})"),
            ("steel: MTPA torque (the row's estimator) against the study's own unrounded number", s["mtpa"]["row_method"]["T_Nm"], "T_mtpa_10A", "N.m",
             f"out\\barrier\\summary.csv, column 'no band'; tolerance {SAME_MODEL_TOLERANCE_PERCENT:g} % ({OWN_LIMIT})"),
            ("steel: cogging peak to peak against the study's own number", s["cogging"]["pp_mNm"], "cogging_pp_mNm", "mN.m",
             f"out\\barrier\\summary.csv, column 'no band' (the README gives the design's number only); tolerance {SAME_MODEL_TOLERANCE_PERCENT:g} % ({OWN_LIMIT})")):
        if stored is None:
            out.append(skipped(name, f"the study's own number, within {SAME_MODEL_TOLERANCE_PERCENT:g} % of it, the same model on the same FEMM ({OWN_LIMIT})",
                               "out\\barrier\\summary.csv is not there: barrier_study.py makes it", source))
        else:
            out.append(compare(name, measured, stored["no band"][column], SAME_MODEL_TOLERANCE_PERCENT, source, unit))
    air = summary["air"]["cogging"]
    out.append(compare("air: cogging peak to peak: the README's 25.7 mN.m is this variant (the design with its band)", air["pp_mNm"], README_COGGING_WITH_BAND_MNM,
                       TOLERANCE_PERCENT, f"README.md, Expected numbers: 'The cogging torque of the design, peak to peak, is 25.7 mN.m'; tolerance {TOLERANCE_PERCENT:g} % ({OWN_LIMIT})",
                       "mN.m"))
    return out


def mesh_study_summary(summary: dict, rows: list[dict]) -> dict:
    """psi_PM and the torque at the MTPA angle for both variants, on the normal mesh and on the finer one."""
    out: dict = {}
    for v in ("none", "steel"):
        g_ref = summary[v]["mtpa"]["bench_method"]["best_sample_gamma_deg"]
        entry: dict = {"gamma_deg": g_ref}
        for scale in (1.0, *MESH_SCALES):
            mine = [r for r in rows if r["variant"] == v and r["mesh_scale"] == scale]
            t = [r["T_age"] for r in mine if r["study"] == "gamma" and r["gamma"] == g_ref]
            cog = [r["T_age"] for r in mine if r["study"] == "cogging" and r["theta_e"] <= 20.0 + 1e-9]
            entry[f"s{scale:g}"] = {"psi_pm_mWb": psi_pm_of(mine)["psi_pm_mWb"], "T_at_gamma_Nm": float(np.mean(t)), "positions": len(t),
                                    "cogging_pp_mNm": 1e3 * (max(cog) - min(cog)), "nodes": next((r["nodes"] for r in mine if r["study"] == "drawn"), None)}
        out[v] = entry
    return out


def mesh_study_checks(summary: dict, rows: list[dict]) -> list[dict]:
    """Is the difference between the two geometries a mesh effect? Then it must shrink on the finer mesh (or already be below 0.05 %): the bench's own
    criterion, since the brief only asks that a miss be shown to be the mesh."""
    study = mesh_study_summary(summary, rows)
    out = []
    for key, unit in (("psi_pm_mWb", "mWb"), ("T_at_gamma_Nm", "N.m"), ("cogging_pp_mNm", "mN.m")):
        a1, b1 = study["none"]["s1"][key], study["steel"]["s1"][key]
        a2, b2 = study["none"][f"s{MESH_SCALES[0]:g}"][key], study["steel"][f"s{MESH_SCALES[0]:g}"][key]
        d1, d2 = 100 * (a1 / b1 - 1), 100 * (a2 / b2 - 1)
        out.append(check(f"mesh study: none against steel in {key}: the gap falls when both meshes are finer", abs(d2) < abs(d1) or abs(d2) < 0.05,
                         f"the difference between the two geometries is a mesh effect, so it shrinks with the mesh size (or is below 0.05 % already) ({OWN_LIMIT})",
                         f"{unit}: normal mesh {a1:.5g} against {b1:.5g} ({d1:+.3f} %), blocks x {MESH_SCALES[0]:g}: {a2:.5g} against {b2:.5g} ({d2:+.3f} %); "
                         f"nodes none {study['none']['s1']['nodes']:,} -> {study['none'][f's{MESH_SCALES[0]:g}']['nodes']:,}, steel {study['steel']['s1']['nodes']:,} -> "
                         f"{study['steel'][f's{MESH_SCALES[0]:g}']['nodes']:,}", "FEMM, the same model on a finer mesh"))
    return out


def reference(folder: Path, log: Log, label: str) -> dict:
    """FEMM's side of B3: B1's band-free model, solved on a femm_pool the way B1 does, at the no-load rotor positions of the first stage (120 positions over
    one electrical period) and at the drawn position (with its 360 samples of the air-gap flux density). Returns {"drawn": the drawn-position row, "psi_pm_mWb":
    psi_PM as B1 defines it (the fundamental of phase A's flux linkage over the 120 positions), "rows": every row}."""
    fem_dir = Path(folder) / "femm_reference"
    jobs = [j for j in stage_one_jobs() if j["variant"] == "none" and j["study"] in ("noload", "drawn")]
    started = time.perf_counter()
    with log.tee():
        with femm_queue.femm_pool(label, task_timeout_s=TASK_TIMEOUT_S) as pool:
            results = run_stage(pool, jobs, fem_dir, "reference: the 120 no-load positions and the drawn position", log)
            record = pool_record(label, pool, results, time.perf_counter() - started)
    log.add_run(record)
    rows = table(results)
    return {"drawn": next(r for r in rows if r["study"] == "drawn"), "psi_pm_mWb": psi_pm_of([r for r in rows if r["study"] == "noload"])["psi_pm_mWb"], "rows": rows}
