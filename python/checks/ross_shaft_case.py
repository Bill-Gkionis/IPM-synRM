"""ross_shaft_case.py: ROSS's half of the ROSS row of check_tools.py. It runs INSIDE ROSS's own environment, never in P's.

What it does, and why it is a separate file
    ROSS 3.1.0 pulls in about seventy packages (an AI client, Excel automation, thermodynamics libraries) that the motor's own
    environment must not carry, so it lives in envs\\ross with a virtual environment of its own (R3 of the tools plan). Python
    cannot import a package from another environment, so the work is split: ross_check.py (P's side) starts the python of that
    environment through femm_queue.run_solver and hands it this script, and this script builds the model with ROSS and writes
    its numbers to a file. P's side reads the file. They share nothing but two JSON files in the case folder:
        case.json    in:   what to build  {E, rho, G_s, diameter, length, elements, bearing_stiffness}
        result.json  out:  what ROSS said  {ross_version, wn_rad_s, ...}
    So this file imports nothing from P (it could not: P's modules are not in this environment), only the standard library
    and ROSS.

The model (a uniform shaft on two very stiff bearings at its ends)
    A steel bar of diameter d and length L, cut into `elements` equal beam elements (ROSS's ShaftElement, a Timoshenko beam: it
    includes shear and rotary inertia, which ROSS switches on by default and which the pipeline's real shafts will use). A
    bearing at each end (node 0 and node `elements`) is a spring in x and in y with stiffness k, and no damping. A spring with
    k = 1e11 N/m is, for a shaft whose own bending stiffness is of the order of 1e6 N/m, a wall: the ends cannot move sideways but
    can turn freely, so this is the "simply supported" shaft of the textbook. run_modal(speed=0) gives the natural frequencies at
    rest (no gyroscopic effect), in rad/s, in ascending order; the shaft bends in x and in y alike, so they come in equal pairs.

Run (by ross_check.py; by hand for a look, with ROSS's python)
    envs\\ross\\.venv\\Scripts\\python.exe -I checks\\ross_shaft_case.py <case folder that holds case.json>
Expect (2026-11-06): ross_version 3.1.0; for the brief's shaft (E 211 GPa, rho 7,810 kg/m^3, d 15 mm, L 300 mm) the first natural
frequency is about 340 Hz, a little under the Euler-Bernoulli formula's 340.2 Hz because the Timoshenko beam is softer.
"""
import json
import math
import sys
from pathlib import Path


def main(folder: str) -> None:
    folder = Path(folder)
    case = json.loads((folder / "case.json").read_text(encoding="utf-8"))
    import ross as rs                                           # after the case is read: a bad case.json fails fast, without ROSS's start-up

    steel = rs.Material(name="steel", rho=case["rho"], E=case["E"], G_s=case["G_s"])
    n = int(case["elements"])
    length = case["length"] / n
    shaft = [rs.ShaftElement(L=length, idl=0.0, odl=case["diameter"], material=steel) for _ in range(n)]
    k = case["bearing_stiffness"]
    bearings = [rs.BearingElement(node, kxx=k, kyy=k, cxx=0.0, cyy=0.0) for node in (0, n)]
    rotor = rs.Rotor(shaft, [], bearings)
    modal = rotor.run_modal(speed=0, num_modes=24, sparse=False)
    wn = [float(w) for w in modal.wn[:12]]
    result = {
        "ross_version": rs.__version__,
        "elements": n,
        "ndof": int(rotor.ndof),
        "mass_kg": float(rotor.m),
        "wn_rad_s": wn,
        "wn_hz": [w / (2 * math.pi) for w in wn],
    }
    (folder / "result.json").write_text(json.dumps(result, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
