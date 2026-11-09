"""ezdxf in 20 lines: write a DXF a CAD program can open, then read it back.
Run, from the learn folder: uv run python 02_ezdxf_basics.py   (writes tiny.dxf into the current directory)"""
import ezdxf

doc = ezdxf.new("R2010", units=ezdxf.units.MM)    # a new, empty drawing in millimetres
doc.layers.add("N45", color=1)                    # layers = named groups, e.g. one per material (1 = red)
msp = doc.modelspace()                            # the sheet we draw on
msp.add_line((0, 0), (10, 0))                                                      # a line from a to b
msp.add_arc(center=(0, 0), radius=10, start_angle=0, end_angle=90)                 # always counter-clockwise
msp.add_lwpolyline([(2, 2), (6, 2), (6, 3), (2, 3)], close=True, dxfattribs={"layer": "N45"})   # a 'magnet'
doc.saveas("tiny.dxf")

for e in ezdxf.readfile("tiny.dxf").modelspace():                                  # read it back
    if e.dxftype() == "LINE":
        print("LINE", e.dxf.start.xy, "->", e.dxf.end.xy, "layer", e.dxf.layer)
    elif e.dxftype() == "ARC":
        pts = list(e.flattening(0.01))           # an arc as points, max 0.01 mm from the true curve
        print(f"ARC r={e.dxf.radius} {e.dxf.start_angle}->{e.dxf.end_angle} deg, {len(pts)} points when flattened")
    else:
        print(e.dxftype(), "layer", e.dxf.layer, "corners", [(float(v[0]), float(v[1])) for v in e.get_points()])
