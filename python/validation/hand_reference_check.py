"""hand_reference_check.py: is the hand-drawn reference .FEM still the model that validate_full.py works on?

validate_full.py compares the motorgeo slice with the hand-drawn full 360-degree model, `IPM synRM analysis.FEM`. It never touches
that file: each job copies it (hand_000.fem, hand_001.fem, ...), changes three things on purpose, saves, and solves.
This check asks: if I make those same three changes to the reference as it is today, do I get the copy that job made?
If yes, nothing else in the reference has changed since that copy was made (a re-save by another FEMM, an edit in the
GUI), and validate_full's numbers still describe the same model.

The three changes are the ones in validate_full.hand_worker (read it next to this script):
  1. the air-gap label (validate_full.AIR_GAP_LABEL, at 45.398, 0.396 mm, which sits in the rotor) is put into group 0,
     so that the torque integral over group 1 leaves it out;
  2. everything in group 1 (the rotor) is rotated about the origin by  move = (phase-A axis + theta_e / pole pairs)
     - 0.5 degrees, which puts the rotor's d-axis where the job wants it. A magnet turns with the rotor, so FEMM adds
     the same angle to the magnetisation direction of each magnet label it rotates (FEMM 4.2 did: 390.5 deg became
     455.025 deg, the difference being exactly `move`);
  3. the three circuits get the job's currents: Phase A = I, Phase B = Phase C = -I/2 (zero at no load).
Nothing else is changed. Rotating is done here with the textbook formula x' = x cos(a) - y sin(a), y' = x sin(a) + y cos(a)
(counter-clockwise), and the sign is itself checked: a wrong sign would put the 228 rotor nodes in the wrong place.

This is a second copy of what hand_worker does, so nothing but a test ties the two together. The air-gap label is one
constant that both import. The rest is held by tests/test_hand_reference_check.py, which reads hand_worker's source
(with `ast`, so comments and strings do not count) and fails if it stops doing one of the three things, does a fourth
thing to the copy, or turns by a different formula; and the list of jobs is pinned the same way. If it fails, change
apply_hand_worker below to follow hand_worker, then the test. The rule "a magnet's direction turns with it" was read
from FEMM 4.2's copies: if a copy that FEMM 5.0 made shows magdir differences on labels that are not magnets, FEMM 5.0
turns more than 4.2 did; read those lines before changing the rule.

Run (from the Python to FEMM folder):
    uv run python validation/hand_reference_check.py HAND.fem [--job K] [--reference REF.fem]
HAND.fem is one copy that validate_full.py made (job 0 by default: no load, theta_e = 0; --job K picks another job
from the same list validate_full uses). Prints the plain differences (to show what the job changed), then the
differences that are left once those changes are applied to the reference. Exit status 0 if nothing is left, 1 if
something is left (or the reference lacks what the job changes), 2 if a file cannot be read or the job number is wrong.

Numbers to expect for the FEMM 4.2 copy of 2026-10-28 (job 0) against the reference as FEMM 5.0 re-saved it on 2026-11-04:
  plain comparison: 984 differences = 3 circuits + 456 nodes (the 228 rotor nodes, each missing from one side) + 288
  segments (144) + 168 arcs (84) + 69 labels (34 rotor labels and the air-gap label's group); after the three changes:
  0 differences. Only the way [Precision] is written (1e-008 against 1e-08) remains, and that is not the model.
  The same holds for each of the 48 copies, with --job K for hand_KKK.fem.
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

P = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(P))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import validate_full as vf  # noqa: E402  (its constants and fem_axes(): the numbers the job itself used)
import same_model as sm  # noqa: E402
from validate_full import AIR_GAP_LABEL  # noqa: E402  (the label hand_worker moves to group 0: one constant, not a copy)

ROTOR_GROUP = 1                                             # hand_worker: femm.mi_selectgroup(1)


def job_parameters(job: int) -> tuple[float, float]:
    """(theta_e [deg], phase-A current [A]) of validate_full's hand job number `job` (0 to 47)."""
    jobs = [(th, i) for i in (0.0, vf.I_LOAD) for th in vf.THETAS]
    if not 0 <= job < len(jobs):
        raise ValueError(f"validate_full has hand jobs 0 to {len(jobs) - 1}, not {job}")
    return float(jobs[job][0]), float(jobs[job][1])


def apply_hand_worker(model: dict, theta_e: float, i_a: float) -> dict:
    """The reference after the three changes of validate_full.hand_worker, on the parsed model (the file is not touched).

    Returns a new model; the argument is left as it was. Raises ValueError if the reference has no air-gap label or no
    circuit called Phase A, B or C: then it is no longer the model hand_worker was written for.
    """
    try:
        d0, a0 = vf.fem_axes()                              # (0.5, phase-A axis of the reference), both in degrees
    except KeyError as exc:
        raise ValueError(f"the winding of the reference has a circuit that validate_full does not know: {exc}") from exc
    move = (a0 + theta_e / vf.P.pole_pairs) - d0
    out = {**model, "nodes": [dict(n) for n in model["nodes"]], "labels": [dict(lb) for lb in model["labels"]],
           "props": {k: {n: {"fields": dict(b["fields"]), "rows": b["rows"]} for n, b in blocks.items()}
                     for k, blocks in model["props"].items()}}
    # 1. the air-gap label: the label nearest to the point the job selects goes to group 0 (and so does not turn)
    gap = min(out["labels"], key=lambda lb: math.hypot(lb["x"] - AIR_GAP_LABEL[0], lb["y"] - AIR_GAP_LABEL[1]))
    if math.hypot(gap["x"] - AIR_GAP_LABEL[0], gap["y"] - AIR_GAP_LABEL[1]) > 1e-6:
        raise ValueError(f"the reference has no label at the air-gap point {AIR_GAP_LABEL}")
    gap["group"] = 0
    # 2. everything left in group 1 turns about the origin, and a magnet's direction turns with it
    c, s = math.cos(math.radians(move)), math.sin(math.radians(move))
    materials = out["props"].get("BlockProps", {})
    for item in (*out["nodes"], *out["labels"]):
        if item["group"] == ROTOR_GROUP:
            item["x"], item["y"] = item["x"] * c - item["y"] * s, item["x"] * s + item["y"] * c
            if item.get("material") in materials and materials[item["material"]]["fields"].get("H_c", 0.0) != 0.0:
                item["magdir"] += move
    # 3. the circuit currents
    circuits = out["props"].get("CircuitProps", {})
    for name, amps in (("Phase A", i_a), ("Phase B", -i_a / 2), ("Phase C", -i_a / 2)):
        if name not in circuits:
            raise ValueError(f'the reference has no circuit called "{name}"')
        circuits[name]["fields"]["TotalAmps_re"] = amps
    return out


def summary(differences: list[str]) -> str:
    """The count of differences by what they are about, for example 'node 448, segment 120, circuit 3'."""
    words = Counter(d.split()[0] for d in differences)
    return ", ".join(f"{word} {n}" for word, n in words.most_common()) or "none"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Does the reference .FEM still equal the copy a validate_full job made?")
    ap.add_argument("hand", help="a copy that validate_full.py made, for example out/validate/hand_000.fem")
    ap.add_argument("--job", type=int, default=0, help="which validate_full hand job made that copy (default 0)")
    ap.add_argument("--reference", default=vf.FEM, help="the hand-drawn .FEM (default: the one validate_full uses)")
    args = ap.parse_args(argv)
    vf.FEM = args.reference                                  # fem_axes() reads its winding from this file
    try:
        theta_e, i_a = job_parameters(args.job)
    except ValueError as exc:
        ap.error(str(exc))
    models = []
    for label, path in (("HAND", args.hand), ("REFERENCE", args.reference)):
        try:
            models.append(sm.parse_fem(path))
        except (OSError, ValueError) as exc:
            print(f"could not read {label}: {exc}", file=sys.stderr)
            return 2
    hand, reference = models

    plain = sm.compare(hand, reference)
    print(f"HAND       {args.hand}   (job {args.job}: theta_e = {theta_e:g} deg, phase A = {i_a:g} A)")
    print(f"REFERENCE  {args.reference}")
    print(f"\nPlain comparison (HAND is OLD, REFERENCE is NEW): {len(plain['differences'])} differences: "
          f"{summary(plain['differences'])}")

    try:
        changed = apply_hand_worker(reference, theta_e, i_a)
    except ValueError as exc:
        print(f"\nThe reference is not the model hand_worker changes: {exc}")
        return 1
    result = sm.compare(hand, changed)
    print("\nAfter making hand_worker's three changes to the reference (HAND is OLD, the changed reference is NEW):\n")
    print(sm.report(args.hand, "the reference with hand_worker's changes", result))
    return 0 if not result["differences"] else 1


if __name__ == "__main__":
    sys.exit(main())
