"""femm_check: FEMM 5.0 solves the motor at 10 A, and its two ways of computing the torque agree.

What it does, and why
    FEMM (finite element method magnetics) is the solver the whole motor study stands on. This row does not re-implement a
    test; it runs check_femm.py, which already solves the default motor at four operating points, and reads its answer. At the
    point iq = 10 A (id = 0) the torque is computed in two independent ways from the same solution:
      * T_age: from the air-gap element, FEMM's own device for the torque of a rotating machine (a ring of special elements in the
        air gap carries the field, and the torque comes from the energy of that ring);
      * T_wst: from the weighted stress tensor, an integral of the Maxwell stress over the rotor's surface.
    Two methods that share only the field solution must give the same torque when the field is right and the mesh is fine
    enough. They agree to 0.1 % (0.8087 and 0.8085 N.m on the coarse mesh: 0.025 %). If the solver, the mesh or the air-gap
    model were wrong, the two would drift apart first.

How it runs FEMM (R10)
    `python check_femm.py` is started as a child through femm_queue.run_solver, from the case folder (check_femm.py writes its
    out\\ipm.fem next to where it runs): the FEMM it opens, and the fkn.exe that FEMM starts for each solve, are inside the
    child's job, below normal priority, under a timeout, and are ended with it. It never touches port 48015 (check_femm opens a
    private FEMM on a free port). The check reads the printed line
        th_e=   0 id=  0 iq= 10 | psi_d= ... | T_age= 0.8087 T_wst= 0.8085 T_dq= ... N.m | 0.4 s
    and does not repeat check_femm.py's logic.

Numbers to expect (2026-11-06, FEMM 5.0.20261101): T_age 0.8087 and T_wst 0.8085 N.m, 0.025 % apart (T_dq 0.7236); 1.6 s.
    The printout is kept in out\\check_tools\\femm\\check_femm.out.txt.

Self-check (the answers follow).
  1. Why is the torque from the flux linkages (T_dq, also printed) not part of the comparison?
  2. Why is the check run in a folder of its own?
Answers.
  1. T_dq comes from the d and q flux linkages and sees only the average torque; at this rotor position the ripple is at its
     top, so T_dq is 10 % lower (0.7236 against 0.8087 N.m; check_femm.py's docstring explains). The other two see the torque
     at the position itself, which is why they are the pair that must agree.
  2. check_femm.py writes out\\ipm.fem and its solution beside where it runs. Another run in P (a study, a test) must not share
     those names with this check: two FEMMs on one file name crash, hang or swap answers.
"""
from __future__ import annotations

import re
import sys
import time
from pathlib import Path

import femm_queue
import paths

from . import FAIL, PASS, Row, fresh, need_file, need_tool, save_log

TOLERANCE_PERCENT = 0.1
TIMEOUT_S = 300
LINE = re.compile(r"th_e=\s*0\s+id=\s*0\s+iq=\s*10\s*\|.*?T_age=\s*(-?[\d.]+)\s+T_wst=\s*(-?[\d.]+)")


def run(folder: Path):
    need_tool("femm")
    script = need_file("check_femm.py", paths.P / "check_femm.py")
    fresh(folder)
    import femm_session
    case = "check_femm.py's 10 A solve: torque from the air-gap element and from the weighted stress tensor"
    expected = f"T_age and T_wst agree within {TOLERANCE_PERCENT} %"
    started = time.perf_counter()
    result = femm_queue.run_solver([sys.executable, script], label="check_tools FEMM (check_femm.py)", cwd=folder,
                                   timeout_s=TIMEOUT_S, threads=2)
    seconds = round(time.perf_counter() - started, 2)
    save_log(folder, "check_femm", result)
    if result.timed_out or result.returncode != 0:
        why = f"timed out after {TIMEOUT_S} s" if result.timed_out else f"exit code {result.returncode}"
        tail = " | ".join((result.stderr or result.stdout).strip().splitlines()[-2:])
        yield Row("FEMM", case, expected, f"check_femm.py failed: {why}. {tail}", FAIL, seconds)
        return
    found = LINE.search(result.stdout)
    solves = result.stdout.count("th_e=")
    if not found or solves != 4:
        yield Row("FEMM", case, expected, f"check_femm.py printed {solves} solve line(s) and no 10 A line with both torques", FAIL, seconds)
        return
    t_age, t_wst = float(found.group(1)), float(found.group(2))
    difference = 100 * abs(t_age - t_wst) / abs(t_age)
    version = femm_session.femm_version() or "?"
    measured = f"FEMM {version}: T_age {t_age:.4f}, T_wst {t_wst:.4f} N.m, {difference:.3f} % apart"
    yield Row("FEMM", case, expected, measured, PASS if difference < TOLERANCE_PERCENT else FAIL, seconds)
