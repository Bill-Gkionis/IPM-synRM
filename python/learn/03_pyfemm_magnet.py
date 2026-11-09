"""pyfemm in 30 lines: one magnet in a box of air, solved by FEMM.
Run, from the learn folder: uv run python 03_pyfemm_magnet.py   (Windows + FEMM 5.0)

Numbers to expect (FEMM 5.0.20261101, run on 2026-11-05; they are the same on every run):
    B above the magnet: Bx = -0.0002 T, By = 0.1378 T
    B inside the magnet: By = 0.1582 T
    mesh: 4872 nodes
Bx is zero in principle (the magnet is symmetric about the y axis, so above its middle the field can only point along
y); the -0.0002 T that prints is the mesh's rounding, not physics. FEMM 4.2 (run on 2026-10-28) printed Bx = +0.0002 T
and By = 0.1383 T above the magnet, By = 0.1572 T inside it and 4646 nodes: it cut the same drawing into other
triangles, so the digits differ a little.
04_workers.py draws the same magnet and reads 0.1662 T inside it: see the note there.
"""
import os
import sys

# `import femm_session` must find femm_session.py, which lives one folder up. Python looks for a module in the folders
# on sys.path, and when a script is run it puts the folder OF THE SCRIPT on that list (here learn\), not the folder
# above it. So this line adds the folder above: abspath(__file__) is this file's full path, the first dirname() takes
# the file name off (learn\), the second goes up one folder (Python to FEMM), and insert(0, ...) puts it first.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import femm
import femm_session

# Start a hidden FEMM that belongs to this script alone. Not femm.openfemm(1): that always dials port 48015, the
# number an interactive FEMM window answers on, and JOINS whatever answers there. Then closefemm() at the end would
# send quit() to that FEMM, an interactive window included. open_femm() picks a free port instead (and sets the
# priority before FEMM starts, so the desktop stays quick). It is explained in femm_session.py.
femm_session.open_femm()
try:
    femm.newdocument(0)                            # 0 = magnetics problem
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, 20, 30)       # units, 2D planar, 20 mm deep
    femm.mi_getmaterial("Air")                     # from FEMM's library
    femm.mi_addmaterial("N45", 1.05, 1.05, 1030272)                  # name, mu_x, mu_y, coercivity Hc [A/m]

    femm.mi_drawrectangle(-30, -30, 30, 30)        # the air box
    femm.mi_drawrectangle(-5, -1, 5, 1)            # the magnet, 10 x 2 mm

    for (x, y), mat, magdir in (((0, 0), "N45", 90), ((0, 20), "Air", 0)):   # one label per closed area
        femm.mi_addblocklabel(x, y)
        femm.mi_selectlabel(x, y)
        femm.mi_setblockprop(mat, 1, 0, "", magdir, 0, 0)          # material, automesh, -, circuit, magnet direction
        femm.mi_clearselected()

    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)   # no flux leaves the box
    for x, y in ((0, -30), (30, 0), (0, 30), (-30, 0)):            # pick each outer edge by a point on it
        femm.mi_selectsegment(x, y)
    femm.mi_setsegmentprop("A=0", 0, 1, 0, 0)
    femm.mi_clearselected()

    femm.mi_saveas(os.path.abspath("magnet.fem"))  # FEMM needs a file; the answer goes to magnet.ans
    femm.mi_analyze(1)                             # mesh + solve (1 = no progress window)
    femm.mi_loadsolution()
    bx, by = femm.mo_getb(0, 2)                    # flux density 1 mm above the magnet's top face
    print(f"B above the magnet: Bx = {bx:.4f} T, By = {by:.4f} T")
    print(f"B inside the magnet: By = {femm.mo_getb(0, 0)[1]:.4f} T")
    print("mesh:", femm.mo_numnodes(), "nodes")
finally:
    femm_session.close_femm()                      # always: quits this FEMM now (it would also go when the script ends)
