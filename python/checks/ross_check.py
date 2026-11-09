"""ross_check: ROSS, in its own environment, finds a shaft's first natural frequency, and the formula agrees.

What it does, and why
    Every shaft has frequencies at which it likes to vibrate, its natural frequencies. A shaft that spins at one of them shakes
    itself apart, so the motor's rotor has to be designed to stay clear of them. ROSS is the library that computes them for a
    shaft with disks and bearings; this row checks that it is installed and gives the right answer for the one shaft whose
    answer is a formula: a uniform bar supported at both ends, which bends like a guitar string (half a sine wave).

The formula (Euler-Bernoulli beam theory, simply supported at both ends)
    A bar of length L, Young's modulus E, density rho, cross-section area A and second moment of area I vibrates in bending with
        omega_n = (n pi / L)^2 sqrt(E I / (rho A))      n = 1, 2, 3, ...
    The first mode (n = 1) in hertz:
        f1 = omega_1 / (2 pi) = pi / (2 L^2) * sqrt(E I / (rho A))
    For a round bar I = pi d^4 / 64 and A = pi d^2 / 4, so E I / (rho A) = E d^2 / (16 rho): the bar's stiffness-to-mass ratio.
    With E = 211 GPa, rho = 7,810 kg/m^3, d = 15 mm and L = 300 mm: sqrt(E d^2 / (16 rho)) = 19.49 m^2/s, and
    f1 = pi / (2 * 0.3^2) * 19.49 = 340.2 Hz.

How ROSS is run (ross_shaft_case.py has the model)
    The python of ROSS's own environment, envs\\ross\\.venv\\Scripts\\python.exe (built from paths.P, never typed), runs
    checks\\ross_shaft_case.py through femm_queue.run_solver. It reads case.json (the numbers above, 30 elements, bearings of
    1e11 N/m at both ends) and writes result.json. This module compares ROSS's first frequency with the formula.
    ROSS's beam elements are Timoshenko beams (they include shear deformation and rotary inertia, which a thick bar feels), so
    ROSS lands a little BELOW the Euler-Bernoulli formula: about 0.3 % for this bar, which has L/d = 20. The row allows 2 %.

Numbers to expect (2026-11-06): ROSS 3.1.0; the formula 340.18 Hz; ROSS 339.17 Hz (-0.30 %), twice (x and y bending are the same
    for a round bar), then 1344.7 Hz (the formula's second mode is 4 times the first, 1360.7 Hz; the Timoshenko beam is softer
    for the shorter wave); about 7 s, most of it ROSS's start-up.

Self-check (the answers follow).
  1. Why does the shaft's first frequency come twice from ROSS?
  2. The bar is cut in half (L = 150 mm). What happens to f1, and what does the test file say?
  3. Would the answer change if the bearings were 100 times softer, 1e9 N/m? (Think about the bar's own stiffness in N/m.)
Answers.
  1. The bar is round, so it bends in the x direction exactly as in the y direction: two modes with the same frequency. A rotor
     with a keyway or bearings of different stiffness in x and y would split the pair.
  2. f1 goes with 1/L^2, so it becomes four times higher, 1,360.7 Hz (tests/test_check_tools.py checks the formula that way).
  3. Hardly: the bar's own bending stiffness, 48 E I / L^3, is about 1e6 N/m, so a spring of 1e9 N/m is already a thousand times
     stiffer than the bar and behaves like a wall (measured 2026-11-06: 339.04 Hz at 1e9 N/m, 339.15 at 1e10, 339.17 at 1e11 and
     at 1e12, so the 1e11 of this row is at the end of the curve, not on its slope).
"""
from __future__ import annotations

import json
import math
import time
from pathlib import Path

import femm_queue
import paths

from . import FAIL, PASS, Row, fresh, need_file, save_log

SHAFT = {"E": 211e9, "rho": 7810.0, "G_s": 81.2e9, "diameter": 0.015, "length": 0.300, "elements": 30, "bearing_stiffness": 1e11}
TOLERANCE_PERCENT = 2.0
VERSION = "3.1.0"                           # the pin of envs\ross\pyproject.toml
TIMEOUT_S = 300


def euler_bernoulli_f1(E: float, rho: float, diameter: float, length: float) -> float:
    """f1 = pi / (2 L^2) * sqrt(E I / (rho A)) [Hz] of a round bar of this diameter [m] and length [m], supported at both ends."""
    inertia = math.pi * diameter ** 4 / 64
    area = math.pi * diameter ** 2 / 4
    return math.pi / (2 * length ** 2) * math.sqrt(E * inertia / (rho * area))


def ross_python() -> Path:
    """The python of ROSS's own environment, built from the project folder."""
    return paths.P / "envs" / "ross" / ".venv" / "Scripts" / "python.exe"


def run(folder: Path):
    python = need_file("ross", ross_python())
    fresh(folder)
    (folder / "case.json").write_text(json.dumps(SHAFT, indent=1), encoding="utf-8")
    case = f"uniform {SHAFT['diameter'] * 1e3:.0f} mm x {SHAFT['length'] * 1e3:.0f} mm steel shaft on two very stiff bearings (ROSS, its own environment)"
    formula = euler_bernoulli_f1(SHAFT["E"], SHAFT["rho"], SHAFT["diameter"], SHAFT["length"])
    expected = f"ROSS {VERSION}; f1 within {TOLERANCE_PERCENT:g} % of {formula:.1f} Hz"
    started = time.perf_counter()
    script = paths.P / "checks" / "ross_shaft_case.py"
    # -I: isolated mode, so nothing of the caller's environment (PYTHONPATH, P's own packages) can leak into ROSS's python
    result = femm_queue.run_solver([python, "-I", script, folder], label="check_tools ROSS shaft", cwd=folder, timeout_s=TIMEOUT_S, threads=2)
    seconds = round(time.perf_counter() - started, 2)
    save_log(folder, "ross", result)                    # ross.out.txt / ross.err.txt: ROSS prints a REFPROP notice at import, harmless
    answer = folder / "result.json"
    if result.timed_out or result.returncode != 0 or not answer.is_file():
        why = f"timed out after {TIMEOUT_S} s" if result.timed_out else f"exit code {result.returncode}, no result.json"
        tail = " | ".join(result.stderr.strip().splitlines()[-2:])
        yield Row("ROSS", case, expected, f"the case script failed: {why}. {tail}", FAIL, seconds)
        return
    data = json.loads(answer.read_text(encoding="utf-8"))
    f1 = data["wn_hz"][0]
    error = 100 * (f1 - formula) / formula
    pair = abs(data["wn_hz"][1] - f1) < 1e-3 * f1
    measured = (f"ROSS {data['ross_version']}: f1 = {f1:.2f} Hz ({error:+.2f} %)" + (", twice (x and y)" if pair else "") +
                f"; next {data['wn_hz'][2]:.0f} Hz; {data['elements']} elements")
    ok = data["ross_version"] == VERSION and abs(error) < TOLERANCE_PERCENT
    yield Row("ROSS", case, expected, measured, PASS if ok else FAIL, seconds)
