"""run_grid: a full-factorial grid over a box, every point scored, the front read from all of them. The optimizer of choice when the
variables are few, the evaluation cheap and the landscape nearly separable: then evaluating everything at a set resolution is exact
at that resolution, and has nothing random to repeat.

Why a grid for the stator tier (the study note, section 5)
    The stator sweeps (out\\optimization\\sweeps_stator) found three variables that matter (slot_open, slot_half_deg, r_wall), the
    torque running monotonically into the tooth and yoke limits, and the ripple's basins set by the slot opening, with the
    interaction between slot_open and slot_half_deg only 3.3 % of the ripple's variance. A few hundred designs cover such a box at
    a resolution finer than any manufacturing tolerance would keep, at about 0.8 s of wall time each. An evolutionary search would
    spend the same budget to approximate what the grid measures outright. (The rotor tier is the opposite case: diagonal ripple
    valleys, strong interaction; it gets NSGA-II, run_moo.py.)

Run (from the Python to FEMM folder)
    uv run python -m optimization.run_grid --run stator_grid --box stator --levels "slot_open=11,slot_half_deg=9,r_wall=6" \\
        --scoring P --fix "{\\"h_open\\": 0.5}" --fix-from out\\optimization\\rotor_nsga2\\knee.json
    --levels   points per variable, evenly spaced from bound to bound
    --fix      overrides held for every point; --fix-from adds a knee.json's overrides (the other tier's chosen design)
Writes out\\optimization\\<run>\\: evaluations.csv, front.csv, knee.json, front.png (run_moo.write_front, the same books).

Self-check (the answers follow).
  1. The grid has 11 x 9 x 6 = 594 points and 60 are refused by the screen. How many FEMM runs does it cost?
  2. When would a grid be the wrong choice here?
Answers.
  1. 534 designs, 23 solves each: the refused ones cost about 10 ms of Python and no FEMM.
  2. With more variables (the count multiplies: 10 levels in 6 variables is a million designs), or with a landscape whose
     features are narrower than the grid's spacing (a valley 0.05 mm wide between points 0.1 mm apart is missed). The sweeps
     measured the features' widths first; that is what made the spacing safe to choose.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

import femm_queue
from optimization import design
from optimization.run_moo import TORQUE_FLOOR, write_front


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", required=True)
    ap.add_argument("--box", required=True)
    ap.add_argument("--levels", required=True, help='e.g. "slot_open=11,slot_half_deg=9,r_wall=6"')
    ap.add_argument("--scoring", default="P", choices=("AT", "J", "P"))
    ap.add_argument("--fix", default="{}")
    ap.add_argument("--fix-from", default="")
    args = ap.parse_args(argv)
    box = json.loads((design.out_dir("") / "boxes.json").read_text(encoding="utf-8"))[args.box]
    levels = {k: int(v) for k, v in (item.split("=") for item in args.levels.split(","))}
    names = list(box)
    fix = json.loads(args.fix)
    if args.fix_from:
        fix = {**json.loads(Path(args.fix_from).read_text(encoding="utf-8"))["overrides"], **fix}
    axes = [np.linspace(box[n][0], box[n][1], levels[n]) for n in names]
    candidates = [{**fix, **dict(zip(names, map(float, point)))} for point in itertools.product(*axes)]
    folder = design.out_dir(args.run)
    (folder / "setup.json").write_text(json.dumps({"args": vars(args), "box": box, "levels": levels, "fix": fix}, indent=2), encoding="utf-8")
    print(f"{args.run}: {len(candidates)} grid points over {names}, fixed {fix}")
    with femm_queue.femm_pool(f"optimization: {args.run}", task_timeout_s=design.TASK_TIMEOUT_S) as pool:
        rows = design.Evaluator(pool, run=args.run, scoring=args.scoring)(candidates, stage=args.run, gen=0)
    ok = [r for r in rows if r["feasible"] and "T_mean" in r]
    print(f"{len(ok)} scored, {len(rows) - len(ok)} refused or failed; {sum(r['T_mean'] < TORQUE_FLOOR for r in ok)} under the torque floor")
    knee = write_front(folder, names)
    print("knee:", json.dumps(knee))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
