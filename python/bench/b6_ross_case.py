"""b6_ross_case.py: ROSS's half of B6. It runs INSIDE ROSS's own environment (envs\\ross), never in P's.

What it does, and why it is a separate file
    ROSS 3.1.0 pulls in about seventy packages that the motor's own environment must not carry (R3 of the tools plan), so it has a virtual environment of its own, and
    Python cannot import a package from another environment. So the work is split exactly as checks/ross_check.py and checks/ross_shaft_case.py split it: B6's
    module (bench/b6_ross.py, P's side) starts the python of that environment through bench.run_program and hands it this script, and this script builds the model with
    ROSS and writes its numbers to a file. They share nothing but two JSON files in the case folder:
        case.json    in:   what to build: the shaft, the disc, the bearings, the Campbell sweep, the variants (see bench/b6_ross.py, case_inputs)
        result.json  out:  what ROSS said: version, mass, the natural frequencies at rest, the critical speeds with their whirl, the Campbell arrays, one line for each variant
    So this file imports nothing from P (it could not: P's modules are not in that environment): the standard library, numpy and ROSS only, and ROSS only inside main(), after
    case.json is read, so that a wrong case fails at once and this file can be imported anywhere. (ROSS prints a REFPROP notice when it is imported: harmless, ignored; the
    numbers come back in result.json, not in what is printed.)

The model (a Timoshenko shaft, a rigid disc at its middle, two bearings at its ends)
    The shaft is `elements` equal beam elements (ROSS's ShaftElement: it includes shear deformation, rotary inertia and the gyroscopic effect, on by default). The disc is a
    DiskElement at the middle node: a mass and two inertias, Id about a diameter and Ip about the axis (ROSS needs them for the tilting modes and for the gyroscopic effect).
    Each bearing is a spring in x and in y, k each, no damping, at the first and the last node. The variants of case.json change the bearing's stiffness k and the shaft's length
    (the span between the bearings) and rebuild everything.
What ROSS is asked
    run_modal(speed=0): the natural frequencies at rest, in rad/s (an isotropic rotor has them in pairs: x and y).
    run_critical_speed(): the speeds at which a natural frequency equals the speed (the iteration of ROSS), with the whirl of each: only a FORWARD whirl is excited by unbalance.
    run_campbell(speed_range, frequencies, "wd"): the damped natural frequencies over the speed range, one row for each speed, with whirl_values (0 forward, 1 backward).

Run (by bench/b6_ross.py; by hand for a look, with ROSS's python)
    envs\\ross\\.venv\\Scripts\\python.exe -I bench\\b6_ross_case.py <case folder that holds case.json>
Expect (2026-11-06): ross_version 3.1.0; for the motor's rotor (a 2.76 kg disc on a 15 mm shaft, 105 mm between two bearings of 3.7e7 N/m) a first natural frequency of about
381 Hz at rest and a first forward critical speed of about 22,800 rpm; about 15 s with the variants.
"""
import json
import math
import sys
from pathlib import Path


def build_rotor(rs, case: dict, stiffness: float, length: float):
    """The ROSS rotor of case.json with bearings of `stiffness` [N/m] and a shaft of `length` [m] between them."""
    shaft = case["shaft"]
    steel = rs.Material(name="shaft", rho=shaft["rho"], E=shaft["E"], G_s=shaft["G_s"])
    n = int(shaft["elements"])
    elements = [rs.ShaftElement(L=length / n, idl=0.0, odl=shaft["diameter"], material=steel) for _ in range(n)]
    disk = rs.DiskElement(int(case["disk"]["node"]), m=case["disk"]["m"], Id=case["disk"]["Id"], Ip=case["disk"]["Ip"])
    damping = case["bearing"]["c"]
    bearings = [rs.BearingElement(node, kxx=stiffness, kyy=stiffness, cxx=damping, cyy=damping) for node in (0, n)]
    return rs.Rotor(elements, [disk], bearings)


def critical_speeds(rotor, num_modes: int) -> dict:
    """ROSS's critical speeds of the rotor: {"wn_rad_s": [...], "whirl": ["Forward", "Backward", ...]} (the whirl each belongs to)."""
    found = rotor.run_critical_speed(num_modes=num_modes)
    return {"wn_rad_s": [float(w) for w in found.wn()], "whirl": [str(w) for w in found.whirl_direction]}


def one_variant(rs, case: dict, variant: dict) -> dict:
    """The first natural frequency at rest and the first forward critical speed of one variant; an error ROSS raises is kept as text, not hidden and not fatal."""
    out = {"name": variant["name"], "k": variant["k"], "length": variant["length"], "f1_hz": None, "critical_forward_rpm": None, "error": None}
    try:
        rotor = build_rotor(rs, case, variant["k"], variant["length"])
        out["f1_hz"] = float(rotor.run_modal(speed=0, num_modes=12, sparse=False).wn[0]) / (2 * math.pi)
        found = critical_speeds(rotor, case["critical"]["num_modes"])
        forward = [w for w, whirl in zip(found["wn_rad_s"], found["whirl"]) if whirl == "Forward"]
        out["critical_forward_rpm"] = min(forward) * 30 / math.pi if forward else None
    except Exception as error:                                              # a variant that ROSS cannot solve is a row with its reason
        out["error"] = f"{type(error).__name__}: {error}"
    return out


def main(folder: str) -> None:
    folder = Path(folder)
    case = json.loads((folder / "case.json").read_text(encoding="utf-8"))
    import numpy as np
    import ross as rs                                                       # after the case is read: a bad case.json fails fast, without ROSS's start-up

    base = case["variants"][0]
    rotor = build_rotor(rs, case, base["k"], base["length"])
    modal = rotor.run_modal(speed=0, num_modes=24, sparse=False)
    speeds = np.linspace(0.0, case["campbell"]["max_rpm"] * math.pi / 30, int(case["campbell"]["points"]))
    campbell = rotor.run_campbell(speeds, frequencies=int(case["campbell"]["frequencies"]), frequency_type="wd")
    result = {
        "ross_version": rs.__version__,
        "elements": int(case["shaft"]["elements"]),
        "ndof": int(rotor.ndof),
        "rotor_mass_kg": float(rotor.m),
        "wn_hz_at_rest": [float(w) / (2 * math.pi) for w in modal.wn[:12]],
        "critical": critical_speeds(rotor, case["critical"]["num_modes"]),
        "campbell": {"speed_rad_s": [float(s) for s in campbell.speed_range], "wd_rad_s": np.asarray(campbell.wd, float).tolist(),
                     "whirl_values": np.asarray(campbell.whirl_values, float).tolist()},
        "variants": [one_variant(rs, case, variant) for variant in case["variants"]],
    }
    (folder / "result.json").write_text(json.dumps(result, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main(sys.argv[1])
