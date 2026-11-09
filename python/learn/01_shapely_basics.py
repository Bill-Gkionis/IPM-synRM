"""shapely in 20 lines: a calculator for flat shapes. Run, from the learn folder: uv run python 01_shapely_basics.py"""
from shapely.affinity import rotate
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import polygonize, unary_union

square = Polygon([(0, 0), (4, 0), (4, 4), (0, 4)])          # corners, in order
hole = Polygon([(1, 1), (3, 1), (3, 3), (1, 3)])
print("square area:", square.area)                            # 16.0
print("square with a hole cut out:", square.difference(hole).area)   # 12.0

magnet = Polygon([(6, 1), (8, 1), (8, 3), (6, 3)])
print("gap between square and magnet:", square.distance(magnet))     # 2.0 -> the 'steel bridge' between them
both = unary_union([square, magnet])
print("square + magnet is a", both.geom_type)                 # MultiPolygon: two separate pieces

print("is (2, 2) inside the square?", square.contains(Point(2, 2)))
print("magnet turned 90 deg about (0, 0):", [tuple(round(v, 3) for v in xy) for xy in rotate(magnet, 90, origin=(0, 0)).exterior.coords][:4])

# loose lines in, closed faces out: how a DXF or .fem file (lines only) becomes regions again
lines = [LineString([(0, 0), (4, 0)]), LineString([(4, 0), (4, 4)]), LineString([(4, 4), (0, 4)]),
         LineString([(0, 4), (0, 0)]), LineString([(2, 0), (2, 4)])]          # a box with a wall down the middle
faces = list(polygonize(unary_union(lines)))
print("faces found:", len(faces), "areas:", [f.area for f in faces])          # 2 faces of 8.0
