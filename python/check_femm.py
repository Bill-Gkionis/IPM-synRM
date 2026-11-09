"""Physics sanity checks on the FEMM model at one rotor position.

Expect: psi_q ~ 0 at no load (d sits on the magnet), the two direct torque methods agreeing (T_age, from the air-gap
element, and T_wst, the weighted stress tensor over the rotor: within 0.1 %), and id < 0 helping torque more than
id > 0 (Ld < Lq). T_dq comes from the flux linkages and sees only the average torque, so it sits below the other two
here, where the ripple is at its top (example.py prints the mean over a ripple period: 0.729 N.m, T_dq 0.7236).
Run: uv run python check_femm.py [--fine]   (from this folder; 2 s, fine mesh 3 s with FEMM 5.0)
--fine selects the fine mesh (the numbers in the project notes); the default is the coarse loop mesh.
Two runs print the same numbers; only the seconds differ.

Numbers to expect with FEMM 5.0 (2026-11-04), coarse mesh (9,580 mesh nodes) and fine mesh (17,934 nodes):
    no load             psi_d / psi_q [mWb]   11.991 / 0.031      12.032 / 0.031
    iq = 10 A           T_age, T_wst [N.m]    0.8087, 0.8085      0.8082, 0.8082
    (id, iq) = (-6, 8)  T_age [N.m]           0.9104              0.9203
    (id, iq) = (+6, 8)  T_age [N.m]           0.3944              0.3826
FEMM 4.2 gave 11,644 and 27,939 nodes, and on the coarse mesh psi_d 11.941, T_age 0.8061, T_wst 0.8059 (iq = 10 A),
0.9110 and 0.3886 N.m. It is the same model (validation/same_model.py): FEMM 5.0 makes fewer triangles from the
same mesh sizes, and the torques move by 0.1 to 1.5 %.
"""
import os
import sys
import time

import femm
import femm_session
import motorgeo as mg

os.makedirs("out", exist_ok=True)
p = mg.Params()
fem = os.path.abspath("out/ipm.fem")

femm_session.open_femm()                          # a hidden FEMM of our own: own port, below-normal (desktop first)
try:
    t = time.perf_counter()
    mg.to_femm(p, mg.build(p), fem, fine="--fine" in sys.argv)
    print(f"to_femm {time.perf_counter() - t:.2f} s")
    femm.mi_close()
    m = mg.FemmModel(p, fem)
    for theta_e, i_d, i_q in [(0, 0, 0), (0, 0, 10), (0, -6, 8), (0, 6, 8)]:     # 10 A x 10 turns = 100 A-turns
        t = time.perf_counter()
        r = m.solve(theta_e, i_d, i_q)
        print(f"th_e={theta_e:4.0f} id={i_d:3.0f} iq={i_q:3.0f} | psi_d={r['psi_d'] * 1e3:7.3f} psi_q={r['psi_q'] * 1e3:7.3f} mWb"
              f" | T_age={r['T_age']:7.4f} T_wst={r['T_wst']:7.4f} T_dq={r['T_dq']:7.4f} N.m | {time.perf_counter() - t:.1f} s")
finally:
    femm_session.close_femm()
