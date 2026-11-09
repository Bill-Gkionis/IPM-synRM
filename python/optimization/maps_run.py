"""optimization.maps_run: the batch driver of the flux maps and the iron-loss map: it plans every FEMM job of phase 1, runs them through the solver slot in
batches that let the slot go, and saves every result as it comes, so a stopped run resumes where it stopped.

What it is, and why
    Phase 1 asks FEMM for several thousand solves of the optimized motor: the flux-linkage grid, the checks that decide the grid and the number of positions,
    the frozen-permeability points, and the flux density waveforms of the iron-loss grid. Other runs of the project (the copper-loss and thermal studies among them) use the same solver slot, and one heavy run
    at a time is a rule of this PC (R10), so this driver
      * takes the slot ONLY through femm_queue.femm_pool (at most 8 workers, below-normal priority, a per-task timeout), never around it;
      * works in BATCHES: one pool holds the slot for about --batch-min minutes (default 15) of work, then it is closed and the slot is let go; the driver
        waits --gap-s seconds (default 30: the waiting runs look at the lock about once a second) before it asks for the slot again, so a run that is
        waiting gets it first. Waiting for the slot is normal ("waiting for the solver slot" lines, every 30 s);
      * saves every finished job at once (CSV rows, one npz per point, a done-list), and skips what is saved: run the same command again and it continues;
      * never ends a process itself: only the pool ends its own workers.
    The stage lists and the check points are fixed here before the numbers they check were seen, so no check can be tuned to its result. Two stages were ADDED afterwards, and say so
    (flux2 and fluxneg: the positions check and the mirror check had failed their declared limits); the frozen stage's points are read from the finished maps.

The stages (python -m optimization.maps_run STAGE [STAGE ...]; "all" runs the planned ones in this order)
    regression   the old reference rotor and the chosen design at the points of out\\optimization\\final\\compare.csv (no-load pair, then the 15 positions at
                 the MTPA angle), the way design.score does them: the pipeline must give back the project's known numbers before it is trusted
    lamfill      does FEMM apply the 0.98 lamination fill of M-19 in this model? (decides whether the loss model divides B by the fill)
    kcheck       8 points x 60 positions (1 degree apart): the reference for how many positions the flux maps need
    mirror       the same points with i_q negated: the symmetry psi_d(i_d, -i_q) = psi_d, psi_q(i_d, -i_q) = -psi_q, T(-i_q) = -T that lets the map stop at i_q = 0
    flux         the grid: i_d from -3 to +1 times the rated current, i_q from 0 to 3 times it, a step of one quarter of it (17 x 13 points), 15 positions each
    heldout      12 points that are on no node of the grid (a seeded random draw), 15 positions each: the independent test of the interpolation
    flux2        the same grid at the 15 positions between the first 15's (with them: 30 positions, 2 degrees apart), because 15 do not resolve the ripple (below)
    fluxneg      the other half of the grid, i_q from -0.25 to -3 times the rated current (17 x 12 points), 15 positions each: solved, not mirrored, because the
                 mirror check failed its declared limits at heavy load (psi_q asymmetry up to 0.76 % of |psi| at the 7 check points and 1.9 % at single points of the plane,
                 torque up to 3 %; the report has the numbers)
    frozen       FEMM's frozen permeability at the MTPA point of 1, 2 and 3 times the rated current and on the q axis at the same currents (needs flux_map.npz)
    loss         the iron-loss grid: 7 x 7 points (i_d / rated = -3 -2 -1.5 -1 -0.5 0 1, i_q / rated = 0 ... 3 in steps of 0.5), 72 positions over 360 degrees each
    losscheck    two nodes of that grid at 144 positions (what 36, 72 and 144 positions say about the sums), three off-grid points and two negative-i_q points at 72
  flux-point jobs are one point each; the loss jobs one point each; every job has its own model files.

Run (from the Python to FEMM folder; the output folder is out\\optimization\\maps\\chosen\\)
    uv run python -m optimization.maps_run regression lamfill            the first checks (about 1 minute of FEMM)
    uv run python -m optimization.maps_run kcheck mirror heldout flux    then the map (about 15 minutes)
    uv run python -m optimization.maps_run loss losscheck                the loss grid (about 20 ... 30 minutes)
    uv run python -m optimization.maps_run status                        what is done, what is left, who holds the slot (no FEMM, no slot)
    options: --batch-min 15  --gap-s 30  --workers 8  --chunk 24  --dry-run (list the jobs, run nothing)
    The design is out\\optimization\\final\\final.json's "overrides"; the rated current is its "fine" block's amps (11.8627 A peak at 10 turns/coil side).
Files: rows_<stage>.csv (flux-type results, one row per point and position), frozen\\*.npz, loss\\point_*.npz and mesh.npz, lamfill.json, meta.json, run.log.

Numbers to expect (this PC, 2026-11-09, FEMM 5.0, 8 workers at below-normal priority, the slot free or shared in 15-minute batches): regression 10 s; lamfill 3 s; kcheck 48 s (60
positions, 38 ... 48 s a job); mirror 50 s; heldout 22 s; flux 221 jobs in 7 min (a 15-position job 9 ... 24 s on its worker; the first chunk, with the corner points at i_d = -3 x rated and
the pool starting, had one of 110 s); flux2 (the other 15 positions) 6.5 min; fluxneg (204 jobs) 6 min; frozen 6 points in 18 s; losscheck 7 jobs in 126 s (a 72-position job 63 s, a
144-position job 126 s on its worker); loss 49 jobs (72 positions each) about 12 min. 15,202 FEMM solves in all (34 + 6 + 480 + 420 + 180 + 3,315 + 3,315 + 3,060 + 216 + 3,528 +
648), 0.1 ... 0.2 s of wall time each on 8 workers.

Self-check (the answers follow).
  1. The driver closes the pool and waits 30 s between batches although it could keep it open. Why?
  2. A job times out twice. What does the driver do, and what does the next run do?
  3. Why are the check points fixed in this file and not chosen after the maps exist?
Answers.
  1. The slot is one lock for the whole PC and other runs wait for it; a run that never lets go locks them out. Closing the pool frees the slot, and the pause
     gives a waiting run (it looks about once a second) the chance to take it before this one asks again.
  2. femm_pool ends the worker and gives the item the result {"status": "timeout"}; the driver prints it, saves nothing for it and goes on. The next run finds
     the job still undone and tries it again (a job that fails every time shows up in `status`).
  3. A check whose points are picked after seeing the data can be made to pass. Fixed in advance, the check says what the maps do at points nobody tuned.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
import time
import uuid
from pathlib import Path

import numpy as np

import femm_queue
import femm_session
import paths
from optimization import design, iron_loss, maps_femm

FINAL_JSON = paths.path("out") / "optimization" / "final" / "final.json"
OUT = paths.path("out") / "optimization" / "maps" / "chosen"
POSITIONS_FLUX = tuple(design.RIPPLE_POSITIONS)                  # 15 positions, 4 degrees apart over the 60 degree ripple period (design.score's own)
POSITIONS_FINE = tuple(float(k) for k in range(60))              # 60 positions, 1 degree apart
POSITIONS_NOLOAD = tuple(design.ANGLE_POSITIONS)                 # (0, 30): design.score's pair for psi_pm
N_LOSS, N_LOSS_FINE = 72, 144                                    # positions over 360 electrical degrees
POSITIONS_FLUX2 = tuple(2.0 + 4.0 * k for k in range(15))        # the 15 positions BETWEEN those of POSITIONS_FLUX: together 30 positions, 2 degrees apart
TASK_TIMEOUT = {"flux": 300, "kcheck": 900, "frozen": 600, "loss": 1200, "losscheck": 2400, "lamfill": 300, "regression": 300, "mirror": 900, "heldout": 300,
                "flux2": 300, "fluxneg": 300}
# the points, in units of the rated current I_R (fixed before any number is seen)
FLUX_STEP = 0.25                                                 # of I_R
FLUX_ID_STEPS, FLUX_IQ_STEPS = range(-12, 5), range(0, 13)       # i_d from -3 to +1 I_R, i_q from 0 to 3 I_R
KCHECK_UNITS = ((0.0, 0.0), ("mtpa", None), (0.0, 1.0), (-1.0, 1.0), (-2.0, 2.0), (-3.0, 3.0), (1.0, 1.0), (-1.5, 3.0))
LOSS_ID_UNITS = (-3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 1.0)
LOSS_IQ_UNITS = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
LOSS_FINE_UNITS = ((-1.0, 1.0), (-2.0, 2.0))                      # nodes of the loss grid, solved again at 144 positions
LOSS_HELD_UNITS = ((-0.75, 0.75), (-1.75, 2.25), (0.5, 1.25))     # off-grid, for the interpolation error of iron_loss()
LOSS_MIRROR_UNITS = ((-1.0, -1.0), (-2.0, -2.0))                  # negative i_q
HELDOUT_SEED, HELDOUT_N = 20261109, 12
LAMFILL_FILLS = (0.98, 1.0, 0.5)
LAMFILL_UNITS = ((0.0, 0.0), (-1.0, 2.0))


# ============================================================================================== the design and the points
def read_final() -> dict:
    with open(FINAL_JSON, encoding="utf-8") as f:
        return json.load(f)


class Context:
    """The design, the rated current and the folders of this run."""

    def __init__(self):
        final = read_final()
        self.overrides = final["overrides"]
        self.I_R = float(final["fine"]["amps"])                      # A peak at the 10-turn placeholder winding
        self.gamma_R = float(final["fine"]["gamma"])                 # the MTPA angle at that current [deg]
        self.out = OUT
        self.fem_dir = OUT / "fem"
        self.out.mkdir(parents=True, exist_ok=True)

    def mtpa_point(self) -> tuple[float, float]:
        return design.currents(self.I_R, self.gamma_R)

    def units(self, u) -> tuple[float, float]:
        """A point given in units of I_R, or ("mtpa", None) for the rated MTPA point, as (i_d, i_q) in A."""
        if u[0] == "mtpa":
            return self.mtpa_point()
        return u[0] * self.I_R, u[1] * self.I_R

    def write_meta(self) -> None:
        meta = {"overrides": self.overrides, "rated_current_A": self.I_R, "rated_gamma_deg": self.gamma_R, "mesh": "fine (smart mesh on)",
                "femm_version": femm_session.femm_version(), "positions_flux_deg": list(POSITIONS_FLUX), "n_loss_positions": N_LOSS, "pole_pairs": 4,
                "winding": "placeholder: 10 turns per coil side, 12 coils per phase in series (120 series turns)", "depth_mm": 20.0,
                "currents": "A peak, amplitude-invariant d-q (the length of (i_d, i_q) is the phase peak current); d = magnet axis"}
        path = self.out / "meta.json"
        if path.exists():
            old = json.loads(path.read_text(encoding="utf-8"))
            for key in ("overrides", "rated_current_A", "rated_gamma_deg"):
                if old[key] != meta[key]:
                    raise RuntimeError(f"{path} was written for another design ({key}: {old[key]} against {meta[key]}); use another folder")
        else:
            path.write_text(json.dumps(meta, indent=1), encoding="utf-8")


def _tag() -> str:
    return uuid.uuid4().hex


def flux_job(ctx: Context, stage: str, key: str, label: str, points, positions, overrides=None) -> dict:
    return {"stage": stage, "key": key, "label": label, "fn": maps_femm.flux_task,
            "task": {"tag": _tag(), "overrides": ctx.overrides if overrides is None else overrides, "fem_dir": str(ctx.fem_dir),
                     "points": [list(p) for p in points], "positions": list(positions)}}


# ============================================================================================== the stage lists
def jobs_regression(ctx: Context) -> list[dict]:
    """The reference rotor and the chosen design, each as design.score does it: the no-load pair, then the 15 positions at the MTPA angle."""
    cmp_path = paths.path("out") / "optimization" / "final" / "compare.csv"
    with open(cmp_path, encoding="utf-8") as f:
        rows = {r["design"]: r for r in csv.DictReader(f)}
    jobs = []
    for label, overrides, name in (("chosen", ctx.overrides, "chosen (fine knee)"), ("reference", {}, "reference")):
        r = rows[name]
        i_d, i_q = design.currents(float(r["amps"]), float(r["gamma"]))
        jobs.append(flux_job(ctx, "regression", f"{label}-noload", label, [(0.0, 0.0)], POSITIONS_NOLOAD, overrides))
        jobs.append(flux_job(ctx, "regression", f"{label}-mtpa", label, [(i_d, i_q)], POSITIONS_FLUX, overrides))
    return jobs


def jobs_kcheck(ctx: Context) -> list[dict]:
    return [flux_job(ctx, "kcheck", f"k{k}", "chosen", [ctx.units(u)], POSITIONS_FINE) for k, u in enumerate(KCHECK_UNITS)]


def jobs_mirror(ctx: Context) -> list[dict]:
    jobs = []
    for k, u in enumerate(KCHECK_UNITS):
        i_d, i_q = ctx.units(u)
        if abs(i_q) > 1e-12:
            jobs.append(flux_job(ctx, "mirror", f"m{k}", "chosen", [(i_d, -i_q)], POSITIONS_FINE))
    return jobs


def flux_grid(ctx: Context) -> list[tuple[float, float]]:
    return [(a * FLUX_STEP * ctx.I_R, b * FLUX_STEP * ctx.I_R) for a in FLUX_ID_STEPS for b in FLUX_IQ_STEPS]


def jobs_flux(ctx: Context) -> list[dict]:
    return [flux_job(ctx, "flux", f"g{a}_{b}", "chosen", [(a * FLUX_STEP * ctx.I_R, b * FLUX_STEP * ctx.I_R)], POSITIONS_FLUX)
            for a in FLUX_ID_STEPS for b in FLUX_IQ_STEPS]


def jobs_flux2(ctx: Context) -> list[dict]:
    """The same grid at the 15 positions between the first 15's: with them 30 positions 2 degrees apart (the positions check found that 15 resolve the mean of
    psi and T to 1e-5 but not the ripple to the declared limit; 30 do)."""
    return [flux_job(ctx, "flux2", f"g{a}_{b}", "chosen", [(a * FLUX_STEP * ctx.I_R, b * FLUX_STEP * ctx.I_R)], POSITIONS_FLUX2)
            for a in FLUX_ID_STEPS for b in FLUX_IQ_STEPS]


def jobs_fluxneg(ctx: Context) -> list[dict]:
    """The negative-i_q half of the grid, solved, not mirrored: the mirror check found the symmetry good to 0.1 ... 0.2 % up to the rated current and 0.5 ... 1 %
    in psi_q and T at 2 ... 3 times it, so the maps carry both halves from FEMM. 15 positions each (the means are good to 1e-5 with them)."""
    return [flux_job(ctx, "fluxneg", f"n{a}_{b}", "chosen", [(a * FLUX_STEP * ctx.I_R, -b * FLUX_STEP * ctx.I_R)], POSITIONS_FLUX)
            for a in FLUX_ID_STEPS for b in range(1, 13)]


def heldout_points(ctx: Context) -> list[tuple[float, float]]:
    rng = np.random.default_rng(HELDOUT_SEED)
    return [(float(rng.uniform(-3.0, 1.0)) * ctx.I_R, float(rng.uniform(0.0, 3.0)) * ctx.I_R) for _ in range(HELDOUT_N)]


def jobs_heldout(ctx: Context) -> list[dict]:
    return [flux_job(ctx, "heldout", f"h{k}", "chosen", [p], POSITIONS_FLUX) for k, p in enumerate(heldout_points(ctx))]


def jobs_lamfill(ctx: Context) -> list[dict]:
    return [{"stage": "lamfill", "key": "lamfill", "label": "chosen", "fn": maps_femm.lamfill_task,
             "task": {"tag": _tag(), "overrides": ctx.overrides, "fem_dir": str(ctx.fem_dir), "fills": list(LAMFILL_FILLS), "theta_e": 0.0,
                      "points": [list(ctx.units(u)) for u in LAMFILL_UNITS]}}]


def frozen_points(ctx: Context) -> list[tuple[str, float, float]]:
    """The MTPA point of 1, 2 and 3 times the rated current (from the maps) and the q-axis point at each: (label, i_d, i_q). Needs flux_map.npz."""
    from optimization import maps
    fm = maps.FluxMaps.load()
    out = []
    for mult in (1, 2, 3):
        gamma = fm.mtpa_angle(mult * ctx.I_R)
        i_d, i_q = design.currents(mult * ctx.I_R, gamma)
        out.append((f"mtpa_{mult}x", i_d, i_q))
    for mult in (1, 2, 3):
        out.append((f"q_{mult}x", 0.0, mult * ctx.I_R))
    return out


def jobs_frozen(ctx: Context) -> list[dict]:
    jobs = []
    for label, i_d, i_q in frozen_points(ctx):
        jobs.append({"stage": "frozen", "key": label, "label": label, "fn": maps_femm.frozen_task,
                     "task": {"tag": _tag(), "overrides": ctx.overrides, "fem_dir": str(ctx.fem_dir), "i_d": i_d, "i_q": i_q,
                              "positions": list(range(0, 60, 10)), "test_current": 10.0}})
    return jobs


def loss_job(ctx: Context, stage: str, key: str, i_d: float, i_q: float, n: int) -> dict:
    return {"stage": stage, "key": key, "label": key, "fn": maps_femm.loss_task, "i_d": i_d, "i_q": i_q, "n": n,
            "task": {"tag": _tag(), "overrides": ctx.overrides, "fem_dir": str(ctx.fem_dir), "i_d": i_d, "i_q": i_q, "n_positions": n}}


def jobs_loss(ctx: Context) -> list[dict]:
    return [loss_job(ctx, "loss", f"p{a:+.2f}_{b:.2f}", a * ctx.I_R, b * ctx.I_R, N_LOSS) for a in LOSS_ID_UNITS for b in LOSS_IQ_UNITS]


def jobs_losscheck(ctx: Context) -> list[dict]:
    jobs = [loss_job(ctx, "losscheck", f"fine{a:+.2f}_{b:.2f}", a * ctx.I_R, b * ctx.I_R, N_LOSS_FINE) for a, b in LOSS_FINE_UNITS]
    jobs += [loss_job(ctx, "losscheck", f"held{a:+.2f}_{b:.2f}", a * ctx.I_R, b * ctx.I_R, N_LOSS) for a, b in LOSS_HELD_UNITS]
    jobs += [loss_job(ctx, "losscheck", f"mirror{a:+.2f}_{b:+.2f}", a * ctx.I_R, b * ctx.I_R, N_LOSS) for a, b in LOSS_MIRROR_UNITS]
    return jobs


STAGES = {"regression": jobs_regression, "lamfill": jobs_lamfill, "kcheck": jobs_kcheck, "mirror": jobs_mirror, "flux": jobs_flux, "heldout": jobs_heldout,
          "flux2": jobs_flux2, "fluxneg": jobs_fluxneg, "frozen": jobs_frozen, "loss": jobs_loss, "losscheck": jobs_losscheck}
ORDER = ("regression", "lamfill", "kcheck", "mirror", "flux", "heldout", "flux2", "fluxneg", "frozen", "loss", "losscheck")


# ============================================================================================== saving what comes back
FLUX_HEADER = ["stage", "label", "key", "i_d", "i_q", "theta_e", "psi_d", "psi_q", "T_age", "T_wst", "T_dq", "psi_retries"]


def _done_file(ctx: Context, stage: str) -> Path:
    return ctx.out / f"done_{stage}.txt"


def done_keys(ctx: Context, stage: str) -> set[str]:
    path = _done_file(ctx, stage)
    return set(path.read_text(encoding="utf-8").split()) if path.exists() else set()


def _mark_done(ctx: Context, stage: str, key: str) -> None:
    with open(_done_file(ctx, stage), "a", encoding="utf-8") as f:
        f.write(key + "\n")


def save_result(ctx: Context, job: dict, res: dict) -> None:
    """Write one finished job to disk, THEN mark it done (a crash in between repeats the job, and duplicates are dropped when the rows are read)."""
    stage, key = job["stage"], job["key"]
    if res["kind"] == "flux":
        path = ctx.out / f"rows_{stage}.csv"
        new = not path.exists()
        with open(path, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new:
                w.writerow(FLUX_HEADER)
            for r in res["rows"]:
                w.writerow([stage, job["label"], key, *(repr(float(v)) for v in r)])
    elif res["kind"] == "frozen":
        (ctx.out / "frozen").mkdir(exist_ok=True)
        np.savez(ctx.out / "frozen" / f"{key}.npz", rows=res["rows"], i_d=res["i_d"], i_q=res["i_q"], test_current=res["test_current"], seconds=res["seconds"])
    elif res["kind"] == "lamfill":
        (ctx.out / "lamfill.json").write_text(json.dumps({"rows": res["rows"].tolist(),
                                                          "columns": ["fill", "i_d", "i_q", "psi_d", "psi_q", "T_age", "T_dq"]}, indent=1), encoding="utf-8")
    elif res["kind"] == "loss":
        folder = ctx.out / "loss" if stage == "loss" else ctx.out / "losscheck"
        folder.mkdir(exist_ok=True)
        mesh = res["mesh"]
        mesh_path = ctx.out / "loss" / "mesh.npz"
        (ctx.out / "loss").mkdir(exist_ok=True)
        if not mesh_path.exists():
            np.savez(mesh_path, region=mesh["region"], area_mm2=mesh["area_mm2"], centroid=mesh["centroid"], r_split=mesh["r_split"], nodes=mesh["nodes"],
                     triangles=mesh["triangles"], steel_area_mm2=np.array(mesh["steel_area_mm2"]), polygon_area_mm2_stator_rotor=np.array(mesh["polygon_area_mm2_stator_rotor"]),
                     depth_mm=20.0, sectors=iron_loss.N_SECTORS)
        else:
            saved = np.load(mesh_path)
            if not (np.array_equal(saved["region"], mesh["region"]) and np.array_equal(saved["area_mm2"], mesh["area_mm2"])):
                raise RuntimeError(f"{key}: the mesh of this job differs from the one saved in {mesh_path}")
        name = f"point_{key}.npz" if stage == "loss" else f"{key}.npz"
        np.savez(folder / name, i_d=job["i_d"], i_q=job["i_q"], n=job["n"], positions=np.array(res["positions"]), psi_d=np.array(res["psi_d"]),
                 psi_q=np.array(res["psi_q"]), T_age=np.array(res["T_age"]), T_wst=np.array(res["T_wst"]), Ak=res["Ak"], xdot15=res["xdot15"], seconds=res["seconds"])
    else:
        raise ValueError(f"unknown result kind {res['kind']!r}")
    _mark_done(ctx, stage, key)


# ============================================================================================== running a stage in batches
def stamp() -> str:
    return time.strftime("%H:%M:%S")


def run_stage(ctx: Context, stage: str, args) -> dict:
    """Run every undone job of a stage through femm_pool, in batches. Returns {"done": n, "failed": [keys]}."""
    jobs = [j for j in STAGES[stage](ctx) if j["key"] not in done_keys(ctx, stage)]
    total = len(STAGES[stage](ctx))
    print(f"[{stamp()}] stage {stage}: {total - len(jobs)} of {total} jobs already done, {len(jobs)} to run", flush=True)
    failed: list[str] = []
    finished = 0
    batch_no = 0
    while jobs:
        batch_no += 1
        t_batch = time.monotonic()
        with femm_queue.femm_pool(f"maps_run: {stage} (batch {batch_no})", workers=args.workers, task_timeout_s=TASK_TIMEOUT[stage]) as pool:
            print(f"[{stamp()}] stage {stage}: the slot is ours, batch {batch_no}", flush=True)
            while jobs and (time.monotonic() - t_batch) < args.batch_min * 60:
                chunk, jobs = jobs[:args.chunk], jobs[args.chunk:]
                t0 = time.monotonic()
                results = pool.map(chunk[0]["fn"], [j["task"] for j in chunk])
                for job, res in zip(chunk, results):
                    if isinstance(res, dict) and res.get("status") == "ok":
                        save_result(ctx, job, res)
                        finished += 1
                    else:
                        said = res.get("error") or res.get("status") if isinstance(res, dict) else repr(res)
                        failed.append(job["key"])
                        print(f"[{stamp()}]   FAILED {stage} {job['key']}: {said}", flush=True)
                secs = [r["seconds"] for r in results if isinstance(r, dict) and r.get("status") == "ok"]
                print(f"[{stamp()}] stage {stage}: {finished} done, {len(jobs)} left; this chunk of {len(chunk)}: {time.monotonic() - t0:.0f} s"
                      + (f" (a job {min(secs):.0f} ... {max(secs):.0f} s on its worker)" if secs else ""), flush=True)
        if jobs:
            print(f"[{stamp()}] stage {stage}: batch {batch_no} over after {(time.monotonic() - t_batch) / 60:.1f} min; the slot is let go for {args.gap_s:g} s", flush=True)
            time.sleep(args.gap_s)
    return {"done": finished, "failed": failed}


def status(ctx: Context) -> None:
    info = femm_queue.holder()
    print("the solver slot is free" if info is None else f"the solver slot is held by {femm_queue._who(info)}")
    for stage in ORDER:
        try:
            total = len(STAGES[stage](ctx))
        except Exception as error:                                   # frozen needs flux_map.npz
            print(f"  {stage:11s} (cannot plan yet: {type(error).__name__}: {error})")
            continue
        print(f"  {stage:11s} {len(done_keys(ctx, stage)):4d} of {total:4d} jobs done")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stages", nargs="+", help="stage names, 'all', or 'status'")
    ap.add_argument("--batch-min", type=float, default=15.0, help="minutes of work per hold of the slot (default 15)")
    ap.add_argument("--gap-s", type=float, default=30.0, help="seconds the slot is let go between batches (default 30)")
    ap.add_argument("--workers", type=int, default=8, help="pool workers, at most 8 (config.toml)")
    ap.add_argument("--chunk", type=int, default=24, help="jobs handed to the pool at once (default 24)")
    ap.add_argument("--dry-run", action="store_true", help="list the jobs, run nothing")
    args = ap.parse_args(argv)
    ctx = Context()
    ctx.write_meta()
    if args.stages == ["status"]:
        status(ctx)
        return 0
    stages = list(ORDER) if args.stages == ["all"] else args.stages
    for s in stages:
        if s not in STAGES:
            ap.error(f"unknown stage {s!r}; the stages are {', '.join(ORDER)}")
    if args.dry_run:
        for s in stages:
            jobs = STAGES[s](ctx)
            print(f"{s}: {len(jobs)} jobs, {len(done_keys(ctx, s))} done")
            for j in jobs[:6]:
                t = j["task"]
                what = {"points": [[round(x, 3) for x in p] for p in t["points"]]} if "points" in t else {"i_d": round(t["i_d"], 3), "i_q": round(t["i_q"], 3)}
                for k in ("positions", "n_positions", "fills"):
                    if k in t:
                        what[k] = len(t[k]) if k == "positions" else t[k]
                print("   ", j["key"], what)
        return 0
    summary = {}
    for s in stages:
        summary[s] = run_stage(ctx, s, args)
    bad = {s: v["failed"] for s, v in summary.items() if v["failed"]}
    print(f"[{stamp()}] finished: " + ", ".join(f"{s} {v['done']} done" for s, v in summary.items()) + (f"; FAILED {bad}" if bad else ""), flush=True)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
