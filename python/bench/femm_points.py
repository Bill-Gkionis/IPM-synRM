"""femm_points: the FEMM side that B7 and B8 share: one coarse-mesh model of the band-free motor, solved at a list of (rotor position, id, iq) points, in a pool worker.

What it is, and why it is its own module
    B7 (Morris screening of four rotor dimensions) and B8 (a flux map) both need the same thing from FEMM: "build THIS design of the motor on the coarse mesh, then solve it at
    these rotor positions and currents, and give me the flux linkages and the torque". B1 does it for three fixed variants and a fine mesh; these two need a design that changes
    from task to task (B7's sample) and a mesh that is the coarse one (the ruling: coarse mesh, so that sixty solves take seconds). So the worker is here, once, and a pool worker
    imports it by its name, bench.femm_points.solve_points (a pool sends a function by name, so it has to live at the top level of a module).

The protocol (how Python passes the model in and reads the numbers out)
    A task is a plain dict: {"id", "tag", "overrides", "fem_dir", "points"}.
      overrides  the motorgeo.Params fields that differ from motorgeo.Params(barrier="none") (R11: the band is left out in Python), {} for the design itself
      points     [{"theta_e", "id", "iq"}]: the rotor position in electrical degrees and the d and q currents in amperes (peak, 10 turns a coil side as in motorgeo)
      tag        the model's file name: <fem_dir>\\<tag>.fem, so that no two FEMMs ever share a file (the lesson of 2026-10-28, motorgeo.FemmModel._flux_linkage)
    The worker builds the design (motorgeo.build), runs motorgeo.check on it and REFUSES a design it rejects (an error naming the reason, never a number), has the FEMM of this
    worker write the model (motorgeo.to_femm, fine=False: FEMM's smart mesh off, about 0.3 s a solve), opens it, and solves each point (motorgeo.FemmModel.solve). FEMM's own
    TCP protocol carries the commands, one line of text each, on the port this worker's FEMM asked the operating system for (never 48015, never femm.openfemm()).
    It returns {"id", "tag", "overrides", "rows", "nodes", "elements", "metrics", "build_seconds", "worker"}: one row per point (the flux linkages psi_d and psi_q in Wb, the
    torque from the air-gap element T_age, from the weighted stress tensor T_wst and from the flux linkages T_dq, in N.m), the size of the mesh FEMM made (read back from FEMM),
    motorgeo.check's metrics, and the worker's facts (its pid and FEMM's, and the priority class of each, read back) that B9 reads.

What else the two share (written here once)
    POSITIONS_EL    the rotor positions both solve at, TASK_TIMEOUT_S    the time limit of one task of the pool, and b1_reference(): B1's numbers from its results.json in the
    ACTIVE root, for the comparisons with the fine mesh; (the numbers, None) or (None, why they are missing: "no B1 results in this run's root" or "B1's results lack <key>").

Run
    Not on its own: B7 and B8 call run_tasks(pool, tasks, label, log) inside a femm_queue.femm_pool. `uv run pytest tests/test_bench_b7.py -q` tests what needs no FEMM.

Self-check (the answers follow).
  1. A design is rejected by motorgeo.check, and the worker raises an error. Why not return a row of NaN and let the case count them?
  2. The model file is written by the worker's FEMM and then read back by the same FEMM. Why not write one model in the main process and hand the path to all the workers?
Answers.
  1. A rejected design is not a number, and a case that carries a NaN into a mean or a Morris index poisons all of it without a word. The case screens every design before any
     FEMM run, so a rejection here means the screen and the worker disagree: that is a bug, and it should stop the run and say which design.
  2. FEMM 5.0 writes its solution next to the model (<tag>.ans) and the circuit currents are set in the open document. Two FEMMs on one file would solve each other's
     currents. A file of its own per task, in a folder of its own, costs nothing and makes it impossible.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import motorgeo as mg

from . import root

POSITIONS_EL = (0.0, 30.0)         # electrical degrees: the pair that cancels the 6th, 18th and 30th harmonics of the torque ripple (B7's choice; B8's map uses it too)
TASK_TIMEOUT_S = 120               # one task is a design at a few points (two to nine solves, a few seconds); a hung FEMM is ended after 2 minutes and the task tried again


def _dig(data, dotted: str):
    """The value at a dotted path of a JSON document ('summary.none.psi_pm_mWb'); a KeyError that names the path as far as it got if a step is not there."""
    here: list[str] = []
    for key in dotted.split("."):
        here.append(key)
        try:
            data = data[key]
        except (KeyError, TypeError, IndexError):
            raise KeyError(".".join(here)) from None
    return data


def b1_reference() -> tuple[dict | None, str | None]:
    """B1's numbers on the fine mesh, from its results.json in the ACTIVE root (bench.root() / B1 / results.json), for B7's and B8's comparisons: (the numbers, None), or (None, why).

    The numbers: path, psi_pm_mWb (summary.none.psi_pm_mWb), gamma_deg (the MTPA angle by B1's own method, summary.none.mtpa.bench_method.gamma_deg), best_sample_gamma_deg and
    best_sample_T_Nm (its best one-degree sample and the torque there), curve ({current angle: T_age_Nm}, summary.none.mtpa.curve) and peak_Nm (the largest of it).

    Why not: 'no B1 results in this run's root: <path> is not there (run B1 first)' (a run on another root has its own B1, or none), 'B1's results (<path>) are not JSON' or
    'cannot be read', 'lack <the dotted key>' (the file is there and is not what was expected) or 'have a value that is not a number'. The comparisons are then SKIPPED
    with that reason: the list of checks is the same on every PC, and a reader can tell a missing run from a damaged file."""
    path = root() / "B1" / "results.json"
    if not path.is_file():
        return None, f"no B1 results in this run's root: {path} is not there (run B1 first)"
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except OSError as error:
        return None, f"B1's results ({path}) cannot be read: {error}"
    except ValueError as error:
        return None, f"B1's results ({path}) are not JSON: {error}"
    try:
        base = "summary.none."
        _dig(data, "summary.none")                                          # so that a missing 'summary' or 'none' is named before anything under it
        psi_pm = _dig(data, base + "psi_pm_mWb")
        values = {key: _dig(data, base + "mtpa.bench_method." + key) for key in ("gamma_deg", "best_sample_gamma_deg", "best_sample_T_Nm")}
        curve_raw = _dig(data, base + "mtpa.curve")
        curve: dict[float, float] = {}
        for angle, entry in curve_raw.items():
            try:
                curve[float(angle)] = float(entry["T_age_Nm"])
            except (KeyError, TypeError):
                raise KeyError(f"{base}mtpa.curve.{angle}.T_age_Nm") from None
        if not curve:
            raise KeyError("summary.none.mtpa.curve (it is empty)")
        reference = {"path": str(path), "psi_pm_mWb": float(psi_pm), **{k: float(v) for k, v in values.items()}, "curve": curve, "peak_Nm": max(curve.values())}
    except KeyError as missing:
        return None, f"B1's results ({path}) lack {missing.args[0]}"
    except (TypeError, ValueError, AttributeError) as error:
        return None, f"B1's results ({path}) have a value that is not a number: {error}"
    return reference, None


def solve_points(task: dict) -> dict:
    """One task of the pool: build `task["overrides"]` on this worker's own FEMM, solve every point of the task on it, and return the rows (see the module docstring).

    Runs in a pool worker, on the FEMM that femm_queue started for it. Raises (the pool hands the error and the worker's traceback to the case) if motorgeo rejects the design."""
    import femm
    import femm_session

    from . import priority_class
    p = mg.Params(barrier="none", **task["overrides"])                       # R11: the band is left out, in Python
    regions = mg.build(p)                                                    # GeometryError (a ValueError) if the design cannot be drawn
    metrics, violations = mg.check(p, regions)
    if violations:
        raise RuntimeError(f"{task['tag']}: motorgeo.check rejects this design {task['overrides']}: {'; '.join(violations)}")
    fem_dir = Path(task["fem_dir"])
    fem_dir.mkdir(parents=True, exist_ok=True)
    fem = str(fem_dir / f"{task['tag']}.fem")
    started = time.perf_counter()
    mg.to_femm(p, regions, fem, fine=False)                                  # the coarse mesh: FEMM's smart mesh off
    femm.mi_close()
    model = mg.FemmModel(p, fem)
    build_seconds = time.perf_counter() - started
    rows: list[dict] = []
    nodes = elements = None
    try:
        for point in task["points"]:
            t = time.perf_counter()
            r = model.solve(point["theta_e"], point["id"], point["iq"])
            rows.append({"theta_e": point["theta_e"], "id": r["id"], "iq": r["iq"], "psi_d": r["psi_d"], "psi_q": r["psi_q"], "T_age": r["T_age"], "T_wst": r["T_wst"],
                         "T_dq": r["T_dq"], "rotor_deg": r["rotor_deg"], "psi_retries": r["psi_retries"], "seconds": time.perf_counter() - t})
        if rows:                                                             # how big is the mesh? ask FEMM, which still has the last solution on disk
            femm.mi_loadsolution()
            nodes, elements = int(femm.mo_numnodes()), int(femm.mo_numelements())
            femm.mo_close()
            ans = Path(fem[:-4] + ".ans")
            if ans.is_file():
                ans.unlink()                                                 # the solution is megabytes; the model (.fem) stays, small
    finally:
        femm.mi_close()                                                      # the pool's FEMM lives on for the next task: leave no document open in it
    me, his = os.getpid(), femm_session.pid(femm._default)
    return {"id": task["id"], "tag": task["tag"], "overrides": task["overrides"], "rows": rows, "nodes": nodes, "elements": elements, "metrics": metrics,
            "build_seconds": build_seconds,
            "worker": {"pid": me, "class": priority_class(me), "femm_pid": his, "femm_class": priority_class(his) if his else None}}


def run_tasks(pool, tasks: list[dict], label: str, log) -> list[dict]:
    """Solve `tasks` on the pool. Returns the results in the order of the tasks. A task that timed out twice or whose worker died is an error: a missing solve would silently
    change a mean or an index, so nothing is returned instead of a table with holes (the same rule as B1's run_stage)."""
    started = time.perf_counter()
    results = pool.map(solve_points, tasks)
    lost = [(t["id"], r) for t, r in zip(tasks, results) if not isinstance(r, dict) or "rows" not in r]
    if lost:
        raise RuntimeError(f"{label}: task(s) {[i for i, _ in lost]} did not come back ({[said for _, said in lost]}): the pool ended their workers; see run.log")
    log.say(f"{label}: {sum(len(t['points']) for t in tasks)} solves in {len(tasks)} tasks, {time.perf_counter() - started:.1f} s")
    return results
