"""A small FEMM problem file as text, for the tests of the tools that read .fem files.

Shared by tests/test_same_model.py and tests/test_hand_reference_check.py. make_fem() writes a complete model in FEMM
4.x/5.0's layout (header, point properties, boundaries, materials with a B-H table, circuits, nodes, segments, arcs,
holes and block labels) and every part of it can be changed through an argument, so a test makes a variant of the
default model by naming the one thing that differs. No FEMM is needed to make or read these files.
"""
from pathlib import Path

BDRY = [("A=0", 0.0), ("pbc0", 4.0)]


def make_fem(*, precision="1e-08", comment="", nodes=None, segments=None, arcs=None, labels=None, holes=None,
             materials=None, circuits=None, bdry=None, pointprops=None, extra="", eol="\r\n"):
    """A small but complete FEMM 4.x/5.0 problem file as text. Every argument can be changed to make a variant.

    The model: three nodes, one segment on boundary A=0, one 90-degree arc, two labels (steel with a B-H table, and a
    magnet), one circuit, one hole. `extra` is text (lines ending in \\r\\n) put in front of the node lines.
    """
    nodes = nodes if nodes is not None else [(0.0, 0.0, 0, 0), (10.0, 0.0, 0, 0), (0.0, 10.0, 0, 0)]
    segments = segments if segments is not None else [(0, 1, -1, 1, 0, 0)]
    arcs = arcs if arcs is not None else [(1, 2, 90.0, 0.5, 2, 0, 0, 0.5)]
    labels = labels if labels is not None else [(3.0, 3.0, 2, 1.5, 0, 0, 1, 1, 0), (6.0, 1.0, 1, -1, 1, 90, 2, 10, 0)]
    holes = holes if holes is not None else [(1.0, 8.0, 3)]
    materials = materials if materials is not None else [
        ("Air", {"Mu_x": 1, "Sigma": 0}, []),
        ("M-19 Steel", {"Mu_x": 4416, "Sigma": 1.9}, [(0, 0), (0.05, 15.120714), (1.5, 1108.325569)]),
    ]
    circuits = circuits if circuits is not None else [("A", 6.0)]
    bdry = bdry if bdry is not None else BDRY
    pointprops = pointprops if pointprops is not None else []
    out = ["[Format]      =  4.0", "[Frequency]   =  0", f"[Precision]   =  {precision}", "[MinAngle]    =  30",
           "[DoSmartMesh] =  0", "[Depth]       =  20", "[LengthUnits] =  millimeters", "[ProblemType] =  planar",
           f'[Comment]     =  "{comment}"', f"[PointProps]   = {len(pointprops)}"]
    for name, amps in pointprops:
        out += ["  <BeginPoint>", f'    <PointName> = "{name}"', f"    <I_re> = {amps}", "  <EndPoint>"]
    out.append(f"[BdryProps]   = {len(bdry)}")
    for name, kind in bdry:
        out += ["  <BeginBdry>", f'    <BdryName> = "{name}"', f"    <BdryType> = {kind:g}", "  <EndBdry>"]
    out.append(f"[BlockProps]  = {len(materials)}")
    for name, fields, bh in materials:
        out += ["  <BeginBlock>", f'    <BlockName> = "{name}"']
        out += [f"    <{k}> = {v}" for k, v in fields.items()]
        out.append(f"    <BHPoints> = {len(bh)}")
        out += [f"      {b}\t{h}" for b, h in bh]
        out.append("  <EndBlock>")
    out.append(f"[CircuitProps]  = {len(circuits)}")
    for name, amps in circuits:
        out += ["  <BeginCircuit>", f'    <CircuitName> = "{name}"', f"    <TotalAmps_re> = {amps}", "  <EndCircuit>"]
    text = eol.join(out) + eol + extra
    out = [f"[NumPoints] = {len(nodes)}"]
    out += ["\t".join(map(str, n)) for n in nodes]
    out.append(f"[NumSegments] = {len(segments)}")
    out += ["\t".join(map(str, s)) for s in segments]
    out.append(f"[NumArcSegments] = {len(arcs)}")
    out += ["\t".join(map(str, a)) for a in arcs]
    out.append(f"[NumHoles] = {len(holes)}")
    out += ["\t".join(map(str, h)) for h in holes]
    out.append(f"[NumBlockLabels] = {len(labels)}")
    out += ["\t".join(map(str, lb)) for lb in labels]
    return text + eol.join(out) + eol


def write(tmp_path, name, text):
    path = tmp_path / name
    path.write_bytes(text.encode("latin-1"))
    return path
