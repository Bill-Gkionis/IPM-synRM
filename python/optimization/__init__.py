"""optimization: the geometry study of the band-free motor (project notes: the Optimization study).

What it is
    The study takes the band-free rotor of the bench (motorgeo.Params(barrier="none"), R11) and the stator around it, and looks
    for the geometries where the torque cannot rise any more without the ripple rising too: the Pareto front of mean torque
    against torque ripple, at the scoring point, inside the fixed envelope (the outer and bore radii, the air gap, the shaft, the
    stack, and the magnets' size, since they are bought). Everything in it is driven from one black box, design.py.

The files, in the order of the study
    design.py        one candidate in, its scores out: decode, screen in Python, score in FEMM (the module every other file calls)
    check_design.py  the black box on the reference design, held against the bench's numbers (B1 and B7) before anybody trusts it
    then, one file per step, each with its own docstring: the one-at-a-time sweeps (the look), the Morris screening on the real
    box, the optimizer runs of each tier, the fine-mesh re-check of the front, the winding study and the export of the chosen
    design to a DXF of the whole motor.

Where it writes
    out\\optimization\\<step>\\ (out\\ is git-ignored, and on this PC only). Every candidate's scores become a CSV row as they come.

The rules it keeps (the tools plan)
    R10  one heavy run at a time: every FEMM goes through femm_queue.femm_pool. Nothing here starts FEMM or ends a process itself.
    R11  the band is left out, in Python: Params(barrier="none").
    The bridge guard stays at 1.4 mm (Params.bridge, decision 17 in the Map): relaxing it is a manufacturing call, not the study's.

Run every file of it from the Python to FEMM folder as a module, e.g.
    uv run python -m optimization.check_design
"""
