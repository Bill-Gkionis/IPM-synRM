"""motorgeo: one parametric IPM synRM cross-section, written out to FEMM, DXF, Gmsh and SolidWorks.

The geometry lives here as Python data (edges + regions). Every tool downstream is just a writer.
Baseline values come from the hand-drawn model `IPM synRM analysis.FEM` (in the IPM-synRM repository's femm/ folder; reference/IPM synRM analysis.FEM here is a byte-exact copy)
unless marked PLACEHOLDER.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from shapely.affinity import rotate
from shapely.geometry import LineString, Point, Polygon
from shapely.ops import polygonize, unary_union


# ----------------------------------------------------------------------------- 1. parameters
@dataclass(frozen=True)
class Params:
    slots: int = 36
    pole_pairs: int = 4
    # stator
    r_so: float = 64.0          # stator outer radius
    r_si: float = 45.8          # bore radius
    slot_open: float = 2.0      # slot opening width at the bore [mm]
    h_open: float = 0.3         # straight (radial) part of the opening
    r_shoe: float = 46.801      # where the tooth shoe ends and the slot walls start
    slot_half_deg: float = 2.57 # half angle of the radial-walled slot body
    r_wall: float = 55.744      # where the walls end and the round slot bottom starts
    # airgap, rotor, stack
    airgap: float = 0.8
    r_shaft: float = 7.5
    stack: float = 20.0
    # V magnets, shared by both pocket layouts
    pocket: str = "slot"        # "slot": the hand-drawn .FEM, magnet in its own slot + air band behind it. "open": magnet inside the band
    v_half_deg: float = 60.0    # magnet long side vs d-axis (V angle 120 deg)
    rib: float = 1.5            # centre rib between the two pockets
    bridge: float = 1.4         # outer bridge, pocket end to rotor surface
    # "slot" layout (pole frame: d-axis = +x, upper magnet). Values read from IPM synRM analysis.FEM
    barrier: str = "air"        # the band behind each magnet: "air", "steel" (same mesh, no barrier) or "none" (not drawn)
    mag_x: float = 39.4136      # magnet-slot centre along the d-axis (its outer corner lands on x = 42.5)
    mag_y: float = 6.0          # magnet-slot centre, distance from the d-axis
    slot_w: float = 9.28        # the .FEM fills the slot with N45: magnet 9.08 x 1.57 plus 0.1 mm clearance all round
    slot_t: float = 1.77
    slot_fillet: float = 0.3
    web: float = 1.5            # steel between magnet slot and air band
    bar_t: float = 1.5          # air band thickness
    bar_apex: float = 33.4771   # shaft-side edge of the band's short part next to the rib
    # "open" layout (the 10-22 prototype)
    apex_r: float = 33.5        # radius of the pocket corner nearest the d-axis
    mag_w: float = 9.08
    mag_t: float = 1.57
    # winding
    coil_pitch: int = 4         # in slots (pole pitch = 4.5)
    turns: int = 10             # turns per coil side, as in the .FEM (PLACEHOLDER: 26 AWG x 10 turns fills ~5% of a slot)

    @property
    def r_ro(self): return self.r_si - self.airgap
    @property
    def n_sectors(self): return math.gcd(self.slots, self.pole_pairs)   # 4 -> 90 deg sector
    @property
    def sector_deg(self): return 360 / self.n_sectors


class GeometryError(ValueError):
    """Raised for parameter sets that can't be drawn. The optimizer scores these as infeasible."""


# ----------------------------------------------------------------------------- 2. primitives
def pol(r, deg):
    a = math.radians(deg)
    return (r * math.cos(a), r * math.sin(a))

def rot(pt, deg):
    a = math.radians(deg); c, s = math.cos(a), math.sin(a)
    return (c * pt[0] - s * pt[1], s * pt[0] + c * pt[1])

def lin(p, d, s):           # p + s*d
    return (p[0] + s * d[0], p[1] + s * d[1])

def key(pt):                # node identity, 1 nm resolution
    return (round(pt[0], 6) + 0.0, round(pt[1], 6) + 0.0)


@dataclass(frozen=True)
class Line:
    a: tuple
    b: tuple
    def points(self): return [self.a]
    def moved(self, deg): return Line(rot(self.a, deg), rot(self.b, deg))


@dataclass(frozen=True)
class Arc:
    a: tuple                 # traversed from a to b ...
    b: tuple
    ccw: bool = True         # ... counter-clockwise (or clockwise) around c
    c: tuple = (0.0, 0.0)

    @property
    def radius(self): return math.dist(self.a, self.c)
    def angles(self):
        a0 = math.degrees(math.atan2(self.a[1] - self.c[1], self.a[0] - self.c[0]))
        a1 = math.degrees(math.atan2(self.b[1] - self.c[1], self.b[0] - self.c[0]))
        return a0, a1
    def sweep(self):
        a0, a1 = self.angles()
        return (a1 - a0) % 360 if self.ccw else (a0 - a1) % 360
    def points(self, step_deg=0.25):
        a0, _ = self.angles(); sw = self.sweep() * (1 if self.ccw else -1)
        n = max(2, math.ceil(abs(sw) / step_deg))
        return [lin(self.c, pol(1, a0 + sw * i / n), self.radius) for i in range(n)]
    def canonical(self):     # same arc, counter-clockwise
        return self if self.ccw else Arc(self.b, self.a, True, self.c)
    def moved(self, deg): return Arc(rot(self.a, deg), rot(self.b, deg), self.ccw, rot(self.c, deg))


@dataclass
class Region:
    """A closed loop of edges. Later regions are painted over earlier ones (magnets over rotor steel)."""
    name: str
    material: str            # "M-19 Steel", "N45", "Air", "Coil", "<No Mesh>"
    edges: list
    label: tuple             # a point inside the visible part of the region
    side: str                # "rotor" or "stator" (what moves, what doesn't)
    mesh: float = 1.0
    magdir: float = 0.0
    circuit: str = ""
    turns: int = 0

    def polygon(self):
        return Polygon([pt for e in self.edges for pt in e.points()])


def ring(r0, r1, deg0, deg1, inner_cuts=()):
    """Edges of an annular sector, counter-clockwise. inner_cuts: extra nodes on the r0 arc."""
    cuts = [deg0, *inner_cuts, deg1]
    inner = [Arc(pol(r0, b), pol(r0, a), False) for a, b in zip(cuts[:-1], cuts[1:])][::-1]
    return [Line(pol(r0, deg0), pol(r1, deg0)), Arc(pol(r1, deg0), pol(r1, deg1)),
            Line(pol(r1, deg1), pol(r0, deg1)), *inner]


# ----------------------------------------------------------------------------- 3. the geometry kernel
def build(p: Params, gap: str = "age") -> list[Region]:
    """Regions of one symmetric sector. gap='age' for FEMM, 'sliding' for Gmsh/Elmer."""
    sec, g = p.sector_deg, p.airgap
    if (2 * p.pole_pairs // p.n_sectors) % 2:
        raise GeometryError("odd poles per sector: needs anti-periodic boundaries")
    R = [Region("rotor_steel", "M-19 Steel", ring(p.r_shaft, p.r_ro, 0, sec),
                pol(p.r_shaft + 5, sec / 2), "rotor", mesh=1.5)]
    pole_pitch = 180 / p.pole_pairs
    pole = {"slot": _pole_slot, "open": _pole}[p.pocket]
    for k in range(2 * p.pole_pairs // p.n_sectors):
        R += pole(p, theta_d=pole_pitch * (k + 0.5), north=(k % 2 == 0), k=k)

    if gap == "age":    # rotor air | unmeshed band carrying FEMM's air-gap element | stator air
        r_a, r_b = p.r_ro + g / 3, p.r_ro + 2 * g / 3
        R.append(Region("gap_rotor", "Air", ring(p.r_ro, r_a, 0, sec), pol((p.r_ro + r_a) / 2, sec / 2), "rotor", g / 3))
        R.append(Region("gap_age", "<No Mesh>", ring(r_a, r_b, 0, sec), pol((r_a + r_b) / 2, sec / 2), "stator"))
    else:               # two air bands meeting at mid-gap: the sliding (mortar) interface for Elmer
        r_b = p.r_ro + g / 2
        R.append(Region("gap_rotor", "Air", ring(p.r_ro, r_b, 0, sec), pol((p.r_ro + r_b) / 2, sec / 2), "rotor", g / 3))

    slots = [(k + 0.5) * 360 / p.slots for k in range(round(p.slots / p.n_sectors))]
    a_o = math.degrees(math.asin(p.slot_open / 2 / p.r_si))
    bore_cuts = [a for phi in slots for a in (phi - a_o, phi + a_o)]
    gap_stator = ring(r_b, p.r_si, 0, sec)
    gap_stator[1:2] = [Arc(pol(p.r_si, a), pol(p.r_si, b)) for a, b in zip([0, *bore_cuts], [*bore_cuts, sec])]
    R.append(Region("gap_stator", "Air", gap_stator, pol((r_b + p.r_si) / 2, sec / 2), "stator", g / 3))
    R.append(Region("stator_steel", "M-19 Steel", ring(p.r_si, p.r_so, 0, sec, bore_cuts),
                    pol(p.r_so - 1.5, sec / 2), "stator", mesh=1.5))

    table = winding(p)
    for k, phi in enumerate(slots):
        R += _slot(p, k, phi, a_o, table[k])
    return R


def _hit(pt, u, r):
    """Distance s along unit vector u from pt to the circle of radius r."""
    pu = pt[0] * u[0] + pt[1] * u[1]
    return -pu + math.sqrt(pu * pu - (pt[0] ** 2 + pt[1] ** 2) + r * r)


def _pole(p, theta_d, north, k):
    b = math.radians(p.v_half_deg)
    u = (math.cos(b), math.sin(b))                # along the magnet, outward
    n = (math.sin(b), -math.cos(b))               # across the magnet, into the V
    c_out = (p.apex_r, p.rib / 2)                 # corner nearest the d-axis
    c_in = lin(c_out, n, -p.mag_t)
    r_lim = p.r_ro - p.bridge
    s_out, s_in = _hit(c_out, u, r_lim), _hit(c_in, u, r_lim)
    if min(s_out, s_in) < p.mag_w + 0.3:
        raise GeometryError(f"magnet reaches the bridge (pole {k})")
    m_out, m_in = lin(c_out, u, p.mag_w), lin(c_in, u, p.mag_w)
    b_out, b_in = lin(c_out, u, s_out), lin(c_in, u, s_in)

    out = []
    for side, tag in ((1, "up"), (-1, "dn")):
        f = lambda q: rot((q[0], side * q[1]), theta_d)   # mirror for the lower magnet, then place
        mag = [Line(f(c_out), f(m_out)), Line(f(m_out), f(m_in)), Line(f(m_in), f(c_in)), Line(f(c_in), f(c_out))]
        bar = [Line(f(m_out), f(b_out)), Arc(f(b_out), f(b_in), side > 0), Line(f(b_in), f(m_in)), Line(f(m_in), f(m_out))]
        magdir = theta_d + side * (p.v_half_deg - 90) + (0 if north else 180)
        mid = lambda *q: (sum(x[0] for x in q) / len(q), sum(x[1] for x in q) / len(q))
        out.append(Region(f"magnet_{k}_{tag}", "N45", mag, f(mid(c_out, m_in)), "rotor", 0.5, magdir % 360))
        out.append(Region(f"barrier_{k}_{tag}", "Air", bar, f(mid(m_out, b_out, b_in, m_in)), "rotor", 0.5))
    return out


def arc3(a, b, c):
    """The short arc from a to b around c. The direction comes from the points, so it survives mirroring."""
    return Arc(a, b, (a[0] - c[0]) * (b[1] - c[1]) - (a[1] - c[1]) * (b[0] - c[0]) > 0, c)


def _pole_slot(p, theta_d, north, k):
    """The hand-drawn layout: each magnet in its own rounded slot, a V air band behind it, a steel web in between."""
    b = math.radians(p.v_half_deg)
    u = (math.cos(b), math.sin(b))                # along the magnet, outward
    n = (math.sin(b), -math.cos(b))               # across the magnet, the way it points on a north pole
    c = (p.mag_x, p.mag_y)                        # slot centre, upper magnet, pole frame (d-axis = +x)
    at = lambda s, t: lin(lin(c, u, s), n, t)     # s along the magnet, t across it, measured from its centre
    hw, ht, f = p.slot_w / 2, p.slot_t / 2, p.slot_fillet
    fil = [((-hw, ht - f), (-hw + f, ht), (-hw + f, ht - f)),     # (start, end, centre) of each corner fillet
           ((hw - f, ht), (hw, ht - f), (hw - f, ht - f)),
           ((hw, -ht + f), (hw - f, -ht), (hw - f, -ht + f)),
           ((-hw + f, -ht), (-hw, -ht + f), (-hw + f, -ht + f))]

    r_lim = p.r_ro - p.bridge
    t_out, t_in = -ht - p.web, -ht - p.web - p.bar_t             # band edges, parallel to the magnet
    kink = lambda t, x: at((x - c[0] - t * n[0]) / u[0], t)       # band edge meets the short part at x
    end = lambda t: lin(at(0, t), u, _hit(at(0, t), u, r_lim))    # band edge meets the bridge circle
    R1, R2 = (p.bar_apex, p.rib / 2), (p.bar_apex + p.bar_t, p.rib / 2)
    K_in, K_out, E_in, E_out = kink(t_in, R1[0]), kink(t_out, R2[0]), end(t_in), end(t_out)
    if p.barrier != "none" and min(K_in[1], K_out[1]) < p.rib / 2 + 0.1:
        raise GeometryError(f"air band misses the centre rib (pole {k})")
    if p.barrier != "none" and math.dist(K_out, E_out) < 1.0:
        raise GeometryError(f"air band too short (pole {k})")
    if max(math.hypot(*at(*cc)) + f for _, _, cc in fil) > r_lim:       # farthest point is on a fillet
        raise GeometryError(f"magnet slot reaches the bridge (pole {k})")

    mid = lambda *q: (sum(x[0] for x in q) / len(q), sum(x[1] for x in q) / len(q))
    out = []
    for side, tag in ((1, "up"), (-1, "dn")):
        g = lambda q: rot((q[0], side * q[1]), theta_d)           # mirror for the lower magnet, then place
        slot = []
        for (_, e0, _), (s1, e1, c1) in zip(fil[-1:] + fil[:-1], fil):
            slot.append(Line(g(at(*e0)), g(at(*s1))))
            if f > 0:
                slot.append(arc3(g(at(*s1)), g(at(*e1)), g(at(*c1))))
        if p.barrier != "none":
            band = [Line(g(R1), g(K_in)), Line(g(K_in), g(E_in)), arc3(g(E_in), g(E_out), (0.0, 0.0)),
                    Line(g(E_out), g(K_out)), Line(g(K_out), g(R2)), Line(g(R2), g(R1))]
            out.append(Region(f"barrier_{k}_{tag}", "Air" if p.barrier == "air" else "M-19 Steel", band,
                              g(mid(K_in, K_out, E_in, E_out)), "rotor", 0.5))
        magdir = theta_d + side * (p.v_half_deg - 90) + (0 if north else 180)
        out.append(Region(f"magnet_{k}_{tag}", "N45", slot, g(c), "rotor", 0.5, magdir % 360))
    return out


def _slot(p, k, phi, a_o, layers):
    a_s = math.radians(p.slot_half_deg)
    R_b = p.r_wall * math.tan(a_s)                                   # round bottom, tangent to both walls (185 deg arc)
    area = 0.5 * (p.r_wall ** 2 - p.r_shoe ** 2) * math.sin(2 * a_s) + R_b ** 2 / 2 * (math.pi + 2 * a_s + math.sin(2 * a_s))
    r_mid = math.sqrt(p.r_shoe ** 2 + area / math.sin(2 * a_s))     # splits the coil area in two equal layers
    P = lambda r, d: pol(r, phi + d)
    O1, O2, O3, O4 = P(p.r_si, -a_o), P(p.r_si + p.h_open, -a_o), P(p.r_si + p.h_open, a_o), P(p.r_si, a_o)
    sd = p.slot_half_deg
    S1, M1, W1, S2, M2, W2 = P(p.r_shoe, -sd), P(r_mid, -sd), P(p.r_wall, -sd), P(p.r_shoe, sd), P(r_mid, sd), P(p.r_wall, sd)
    cb, B = P(p.r_wall / math.cos(a_s), 0), P(p.r_wall / math.cos(a_s) + R_b, 0)   # FEMM takes arcs up to 180 deg:
    bottom = [Arc(W1, B, True, cb), Arc(B, W2, True, cb)]                            # split it at the bottom node
    regs = [Region(f"opening_{k}", "Air", [Line(O1, O2), Line(O2, S1), Line(S1, S2), Line(S2, O3), Line(O3, O4),
                   Arc(O4, O1, False)], P((p.r_si + p.r_shoe) / 2, 0), "stator", 0.4)]
    for j, (edges, r_lab) in enumerate((
            ([Line(S1, M1), Line(M1, M2), Line(M2, S2), Line(S2, S1)], (p.r_shoe + r_mid) / 2),
            ([Line(M1, W1), *bottom, Line(W2, M2), Line(M2, M1)], (r_mid + p.r_wall) / 2))):
        phase, sign = layers[j]
        regs.append(Region(f"coil_{k}_{j}", "Coil", edges, P(r_lab, 0), "stator", 1.0,
                           circuit=phase, turns=sign * p.turns))
    return regs


# ----------------------------------------------------------------------------- 4. winding
def slot_angle(p, k): return (k + 0.5) * 360 / p.slots

def winding(p):
    """slot -> ((phase, sign) top layer, (phase, sign) bottom layer). 60-degree phase belts, lap coils."""
    belts = [("A", 1), ("C", -1), ("B", 1), ("A", -1), ("C", 1), ("B", -1)]
    top = lambda k: belts[int(((p.pole_pairs * slot_angle(p, k) + 30) % 360) // 60)]
    return {k: (top(k), (top((k - p.coil_pitch) % p.slots)[0], -top((k - p.coil_pitch) % p.slots)[1]))
            for k in range(p.slots)}

def phase_axes(p):
    """Mechanical angle where each phase's MMF fundamental points outward, plus the winding factor kw1."""
    th = np.radians(np.arange(0, 360, 0.05))
    axes, table = {}, winding(p)
    for ph in "ABC":
        F = np.zeros_like(th)   # crossing a +z conductor counter-clockwise drops the MMF by its ampere-turns
        cond = [(slot_angle(p, k), s) for k, lay in table.items() for q, s in lay if q == ph]
        for a, s in cond:
            F -= s * (th > math.radians(a))
        c = ((F - F.mean()) * np.exp(1j * p.pole_pairs * th)).sum()      # angle(c) = axis, electrical
        axes[ph] = (math.degrees(np.angle(c)) / p.pole_pairs) % (360 / p.pole_pairs)
        if ph == "A":
            kw1 = abs(sum(s * np.exp(1j * p.pole_pairs * math.radians(a)) for a, s in cond)) / len(cond)
    seq = round(((axes["B"] - axes["A"]) * p.pole_pairs) % 360)
    assert seq == 120, f"phase sequence is {seq} deg, swap B and C"
    return axes, kw1


# ----------------------------------------------------------------------------- 5. checks (before any solver runs)
def visible(regions):
    """Polygon of what each region shows once later regions are painted over it."""
    out, covered = {}, Polygon()
    for r in reversed(regions):
        poly = r.polygon()
        out[r.name] = poly.difference(covered)
        covered = covered.union(poly)
    return out

def check(p, regions, min_web=1.0, min_tooth=2.0, min_yoke=3.0):
    """Manufacturability metrics and violations. Shapely, ~10 ms, no solver."""
    poly = {r.name: r.polygon() for r in regions}
    air = {r.name for r in regions if r.material == "Air"}
    n_poles = 2 * p.pole_pairs // p.n_sectors
    pole_keys = [(k, t) for k in range(n_poles) for t in ("up", "dn")]
    pocket = {(k, t): poly[f"magnet_{k}_{t}"].union(poly[f"barrier_{k}_{t}"]) if f"barrier_{k}_{t}" in air
              else poly[f"magnet_{k}_{t}"] for k, t in pole_keys}           # a steel-filled band is no pocket
    wrap = rotate(pocket[(n_poles - 1, "up")], -p.sector_deg, origin=(0, 0))   # last pole, seen across the sector edge
    q_ribs = [pocket[(k, "up")].distance(pocket[(k + 1, "dn")]) for k in range(n_poles - 1)] + [wrap.distance(pocket[(0, "dn")])]
    corners = lambda g: [xy for part in getattr(g, "geoms", [g]) for xy in part.exterior.coords]   # Polygon or MultiPolygon
    slot = [unary_union([poly[n] for n in poly if n.startswith(("opening_", "coil_")) and n.split("_")[1] == str(k)])
            for k in range(round(p.slots / p.n_sectors))]
    m = {
        "centre_rib": min(pocket[(k, "up")].distance(pocket[(k, "dn")]) for k in range(n_poles)),
        "q_rib": min(q_ribs),
        "outer_bridge": p.r_ro - max(math.hypot(*xy) for pk in pocket.values() for xy in corners(pk)),
        "tooth_min": min(a.distance(b) for a, b in zip(slot[:-1], slot[1:])),
        "yoke": p.r_so - max(math.hypot(*xy) for s in slot for xy in s.exterior.coords),
        "magnet_area_mm2": p.n_sectors * sum(poly[n].area for n in poly if n.startswith("magnet_")),
        "coil_layer_area_mm2": poly["coil_0_0"].area,
    }
    m["magnet_mass_g"] = m["magnet_area_mm2"] * p.stack * 7.5e-3      # NdFeB ~7.5 g/cm3
    limits = {"centre_rib": min_web, "q_rib": min_web, "outer_bridge": min_web, "tooth_min": min_tooth, "yoke": min_yoke}
    if p.pocket == "slot" and p.barrier == "air":                       # steel between magnet slot and air band
        m["web"] = min(poly[f"magnet_{k}_{t}"].distance(poly[f"barrier_{k}_{t}"]) for k, t in pole_keys)
        limits["web"] = min_web
    return m, [f"{k} = {m[k]:.2f} < {v}" for k, v in limits.items() if m[k] < v]


# ----------------------------------------------------------------------------- 6. writers
def unique_edges(regions):
    """Every edge once, arcs made counter-clockwise. FEMM and DXF both want this."""
    seen = {}
    for r in regions:
        for e in r.edges:
            e = e.canonical() if isinstance(e, Arc) else e
            k = (frozenset((key(e.a), key(e.b))), key(e.c) if isinstance(e, Arc) else None)
            seen[k] = (e, r)          # later regions win, so an edge takes the layer of what's painted on top
    return list(seen.values())


# FEMM 4.2's material "M-19 Steel", copied field by field from its library and kept on purpose (ruling R5: test the move
# to FEMM 5.0 like for like, with the same steel). Where it comes from:
#   - FEMM 4.2 (21Apr2019), file bin\matlib.dat, block "M-19 Steel" in the folder "Soft Magnetic Materials / Silicon Iron";
#   - the curve is AK Steel's 29-gauge DC magnetisation curve, to within 0.2 % (project note
#     IPM synRM\docs\1-materials\01-magnetic.md, data in materials\research\m19_dc_maker.csv);
#   - FEMM 5.0's library has no "M-19 Steel". Its nearest block, "M-19 14 mil", is a different curve (61 points and a
#     lamination fill of 0.955), so asking FEMM 5.0 for the old name fails and the nearest name would change the motor.
# The keys are the field names of a block in FEMM's files (<Mu_x>, <Sigma>, ...), so tests/test_m19_steel.py can compare
# this dictionary with the library and with a saved .fem field by field. The numbers are FEMM's 17-digit text read as
# doubles and written back in their shortest exact form (0.05 is what the file shows as 0.050000000000000003: the same
# double). B-H curve: B [T] and H [A/m] of the magnetisation curve, the solver interpolates between the points.
M19_STEEL_FEMM42 = {
    "BlockName": "M-19 Steel",
    "Mu_x": 4416,                 # relative permeability along x. The B-H curve below replaces it in the solver
    "Mu_y": 4416,                 # relative permeability along y
    "H_c": 0,                     # coercivity [A/m]: 0, this is not a magnet
    "H_cAngle": 0,                # direction of the coercivity [deg]
    "J_re": 0,                    # applied current density, real part [A/mm^2]
    "J_im": 0,                    # applied current density, imaginary part [A/mm^2]
    "Sigma": 1.9,                 # electrical conductivity [MS/m]
    "d_lam": 0.635,               # lamination thickness [mm]
    "Phi_h": 0,                   # hysteresis lag angle at saturation [deg]
    "Phi_hx": 0,                  # hysteresis lag along x, for linear materials [deg]
    "Phi_hy": 0,                  # hysteresis lag along y, for linear materials [deg]
    "LamType": 0,                 # 0 = not laminated, or laminated in the plane of the problem
    "LamFill": 0.98,              # share of the lamination that is steel, the rest is insulation [-]
    "NStrands": 0,                # wire strands (0: not a wire)
    "WireD": 0,                   # strand diameter [mm] (0: not a wire)
    "BHPoints": 47,               # the number of (B, H) pairs below
    "BH": (                       # (B [T], H [A/m]), lowest first
        (0.0, 0.0),
        (0.05, 15.120714),
        (0.1, 22.718292),
        (0.15, 27.842733),
        (0.2, 31.871434),
        (0.25, 35.365044),
        (0.3, 38.600588),
        (0.35, 41.736202),
        (0.4, 44.873979),
        (0.45, 48.087807),
        (0.5, 51.437236),
        (0.55, 54.975221),
        (0.6, 58.752993),
        (0.65, 62.823644),
        (0.7, 67.245285),
        (0.75, 72.084406),
        (0.8, 77.4201),
        (0.85, 83.350021),
        (0.9, 89.999612),
        (0.95, 97.537353),
        (1.0, 106.201406),
        (1.05, 116.348464),
        (1.1, 128.547329),
        (1.15, 143.765431),
        (1.2, 163.754169),
        (1.25, 191.868158),
        (1.3, 234.833507),
        (1.35, 306.509769),
        (1.4, 435.255202),
        (1.45, 674.911968),
        (1.5, 1108.325569),
        (1.55, 1813.085468),
        (1.6, 2801.217421),
        (1.65, 4053.653117),
        (1.7, 5591.10689),
        (1.75, 7448.318413),
        (1.8, 9708.81567),
        (1.85, 12486.931615),
        (1.9, 16041.483644),
        (1.95, 21249.420624),
        (2.0, 31313.495878),
        (2.05, 53589.446877),
        (2.1, 88477.484601),
        (2.15, 124329.41054),
        (2.2, 159968.5693),
        (2.25, 197751.604272),
        (2.3, 234024.751347),
    ),
}


def define_m19_steel():
    """Define the material "M-19 Steel" in the FEMM document that is open, from M19_STEEL_FEMM42.

    Two calls do it, both ordinary FEMM calls: mi_addmaterial() creates the material with its constants, then
    mi_addbhpoint() adds the 47 points of the B-H curve one by one. FEMM wants the constants in this order
    (the order of its own manual): name, mu_x, mu_y, H_c [A/m], J [A/mm^2], conductivity [MS/m], lamination
    thickness [mm], hysteresis lag [deg], lamination fill, lamination type, lag_x, lag_y, strands, strand diameter [mm].
    mi_addmaterial() cannot set <H_cAngle> and <J_im>; both are 0 in this block, which is also FEMM's default.
    """
    import femm
    m = M19_STEEL_FEMM42
    femm.mi_addmaterial(m["BlockName"], m["Mu_x"], m["Mu_y"], m["H_c"], m["J_re"], m["Sigma"], m["d_lam"],
                        m["Phi_h"], m["LamFill"], m["LamType"], m["Phi_hx"], m["Phi_hy"], m["NStrands"], m["WireD"])
    for b, h in m["BH"]:
        femm.mi_addbhpoint(m["BlockName"], b, h)


def to_femm(p, regions, path, fine=False):
    """Write the regions as a FEMM problem file at `path`, using the FEMM that femm_session.open_femm() started.

    In FEMM's terms, step by step (each is one or a few femm.* calls):
      1. a new magnetics document: planar, millimeters, depth = the stack length, no frequency, 30 deg minimum mesh angle;
      2. the materials: Air from FEMM's library; "M-19 Steel" from define_m19_steel() (FEMM 5.0's library has none);
         N45 (mu_r 1.05, coercivity 1,030,272 A/m, as in the hand-drawn .FEM) and Coil (a plain mu_r = 1 material);
      3. three series circuits A, B, C (their currents are set later, by FemmModel.solve) and the boundaries:
         A = 0 on the shaft and the outer rim, the periodic "AGE" air gap (the rotor turns inside it without a re-mesh)
         and one periodic pair per segment on the sector's 0-degree edge;
      4. every node, then every line and arc once (arcs get a mesh limit, finer near the air gap);
      5. one block label per region with its material, mesh size, circuit, magnetisation direction, group and turns
         (group 1 = rotor bodies, whose torque is integrated; 2 = stator; 3 = the air gap);
      6. mi_saveas(path).
    fine=True switches FEMM's smart mesh on, which gives a finer mesh: the nodes of the default model are about 17,900
    against 9,600 with FEMM 5.0 (FEMM 4.2 made 27,900 against 11,600 from the same file), and the numbers of the
    project notes come from the fine mesh. At one rotor position the two meshes differ in torque by up to about 3 %
    (check_femm.py prints either one, check_femm.py --fine the other).
    """
    import femm
    femm.newdocument(0)
    femm.mi_smartmesh(1 if fine else 0)     # off: ~9.6k nodes, ~0.3 s per solve with FEMM 5.0; on: ~17.9k nodes, ~0.5 s
    femm.mi_probdef(0, "millimeters", "planar", 1e-8, p.stack, 30)
    femm.mi_getmaterial("Air")
    define_m19_steel()                                    # FEMM 4.2's steel, defined here: FEMM 5.0's library lacks it
    femm.mi_addmaterial("N45", 1.05, 1.05, 1030272)       # mu_r, mu_r, Hc [A/m] as in the hand-drawn .FEM
    femm.mi_addmaterial("Coil", 1, 1, 0)
    for ph in "ABC":
        femm.mi_addcircprop(ph, 0, 1)
    femm.mi_addboundprop("A=0", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    femm.mi_addboundprop("AGE", 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0)   # 6 = periodic air gap

    edges = unique_edges(regions)
    for pt in {key(q) for e, _ in edges for q in (e.a, e.b)}:
        femm.mi_addnode(*pt)
    for e, _ in edges:
        if isinstance(e, Line):
            femm.mi_addsegment(*e.a, *e.b)
        else:
            near_gap = p.r_ro - 0.01 < e.radius < p.r_si + 0.01
            femm.mi_addarc(*e.a, *e.b, e.sweep(), 0.5 if near_gap else min(5, e.sweep() / 4))

    def arc_prop(r, prop):    # boundary arcs: shaft, outer, the two AGE arcs
        femm.mi_selectarcsegment(*pol(r, p.sector_deg / 2))
        femm.mi_setarcsegmentprop(0.5, prop, 0, 0)
        femm.mi_clearselected()
    arc_prop(p.r_shaft, "A=0"); arc_prop(p.r_so, "A=0")
    arc_prop(p.r_ro + p.airgap / 3, "AGE"); arc_prop(p.r_ro + 2 * p.airgap / 3, "AGE")

    edge0 = sorted((math.hypot(*e.a), math.hypot(*e.b)) for e, _ in edges
                   if isinstance(e, Line) and abs(e.a[1]) < 1e-9 and abs(e.b[1]) < 1e-9)
    for i, (r1, r2) in enumerate(edge0):          # pair each segment on the 0-deg edge with its twin
        femm.mi_addboundprop(f"pbc{i}", 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0)
        for deg in (0, p.sector_deg):
            femm.mi_selectsegment(*pol((r1 + r2) / 2, deg))
        femm.mi_setsegmentprop(f"pbc{i}", 0, 1, 0, 0)
        femm.mi_clearselected()

    for r in regions:
        femm.mi_addblocklabel(*r.label)
        femm.mi_selectlabel(*r.label)
        group = 3 if r.name.startswith("gap_") else 1 if r.side == "rotor" else 2   # 1 = rotor bodies (torque)
        femm.mi_setblockprop(r.material, 0, r.mesh, r.circuit, r.magdir, group, r.turns)
        femm.mi_clearselected()
    femm.mi_saveas(path)


def to_dxf(p, regions, path, full=True):
    """Whole cross-section (or one sector), one layer per material. Opens in SolidWorks as a sketch."""
    import ezdxf
    doc = ezdxf.new("R2010", units=ezdxf.units.MM)
    msp = doc.modelspace()
    for name in ("M-19", "N45", "AIR", "COIL"):
        doc.layers.add(name)
    layer = {"M-19 Steel": "M-19", "N45": "N45", "Coil": "COIL"}
    done = set()
    for e0, r in unique_edges(regions):
        for i in range(p.n_sectors if full else 1):
            e = e0.moved(i * p.sector_deg)
            k = (frozenset((key(e.a), key(e.b))), key(e.c) if isinstance(e, Arc) else None)
            if k in done:
                continue
            done.add(k)
            att = {"layer": layer.get(r.material, "AIR")}
            if isinstance(e, Line):
                msp.add_line(e.a, e.b, dxfattribs=att)
            else:
                a0, a1 = e.angles()
                msp.add_arc(e.c, e.radius, a0, a1, dxfattribs=att)
    doc.saveas(path)


def compare_dxf(regions, path, sagitta=0.001):
    """Regression test: area of every kernel region vs the DXF face that contains its label point."""
    import ezdxf
    lines = []
    for e in ezdxf.readfile(path).modelspace():
        if e.dxftype() == "LINE":
            lines.append(LineString([e.dxf.start.xy, e.dxf.end.xy]))
        elif e.dxftype() == "ARC":
            lines.append(LineString([v.xy for v in e.flattening(sagitta)]))
    snap = lambda ls: LineString([key(xy) for xy in ls.coords])     # close the 1e-12 gaps between entities
    faces = list(polygonize(unary_union([snap(ls) for ls in lines])))
    vis, rows = visible(regions), []
    for r in regions:
        if r.material == "<No Mesh>":
            continue
        face = next((f for f in faces if f.contains(Point(r.label))), None)
        rows.append((r.name, vis[r.name].area, face.area if face else float("nan")))
    return rows


def to_gmsh(p, regions, path, size_max=1.5):
    """2D mesh with one physical group per region, rotor and stator meshed separately (sliding interface).

    The process's PATH is the same afterwards. gmsh.initialize() (4.15.2, on Windows) cuts the PATH of the running process short (1,864 characters became 313) and
    finalize() does not put it back, so a program started by name later in the same process (powershell.exe, which lives only on the PATH) would not be found. The PATH
    is saved before gmsh.initialize() and restored in a finally after finalize(), whatever happens in between, so that every caller is safe (check_geometry.py, example.py,
    learn/tools_figures.py, and the guards of check_tools and the bench, which stay)."""
    import os
    saved = os.environ.get("PATH")
    try:
        return _to_gmsh(p, regions, path, size_max)
    finally:
        if saved is not None:
            os.environ["PATH"] = saved                      # (assigning os.environ["PATH"] sets the process's own variable again)


def _to_gmsh(p, regions, path, size_max):
    import gmsh
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("ipm_sector")
    occ = gmsh.model.occ
    surf = {}
    for i, r in enumerate(regions):
        if r.material == "<No Mesh>":
            continue
        pts = {}
        def P(xy):
            if key(xy) not in pts:
                pts[key(xy)] = occ.addPoint(xy[0], xy[1], 0, r.mesh)
            return pts[key(xy)]
        curves = []
        for e in r.edges:
            if isinstance(e, Line):
                curves.append(occ.addLine(P(e.a), P(e.b)))
            else:                       # OCC arcs must be < 180 deg: split into <= 90 deg pieces
                n = math.ceil(e.sweep() / 90)
                a0 = e.angles()[0]; sw = e.sweep() * (1 if e.ccw else -1)
                knots = [e.a] + [lin(e.c, pol(1, a0 + sw * j / n), e.radius) for j in range(1, n)] + [e.b]
                c = occ.addPoint(*e.c, 0)
                curves += [occ.addCircleArc(P(q0), c, P(q1)) for q0, q1 in zip(knots[:-1], knots[1:])]
        surf[i] = occ.addPlaneSurface([occ.addCurveLoop(curves)])
    owner = {}
    for side in ("rotor", "stator"):    # fragment each side on its own, so the mid-gap arcs stay two curves
        idx = [i for i in surf if regions[i].side == side]
        _, omap = occ.fragment([(2, surf[idx[0]])], [(2, surf[i]) for i in idx[1:]])
        for i, dimtags in zip(idx, omap):
            for _, t in dimtags:
                owner[t] = i            # later regions win
    occ.synchronize()
    groups = {}
    for t, i in owner.items():
        groups.setdefault(regions[i].name, []).append(t)
    for name, tags in groups.items():
        gmsh.model.addPhysicalGroup(2, tags, name=name)
    bnd = {}                            # boundary curves Elmer needs: A=0 arcs, periodic edges, sliding arcs
    r_gap = p.r_ro + p.airgap / 2
    for _, c in gmsh.model.getEntities(1):
        ends = [gmsh.model.getValue(0, t, []) for _, t in gmsh.model.getBoundary([(1, c)], oriented=False)]
        rad = [math.hypot(x, y) for x, y, _ in ends]
        ang = [math.degrees(math.atan2(y, x)) for x, y, _ in ends]
        if gmsh.model.getType(1, c) == "Circle" and max(abs(r - rad[0]) for r in rad) < 1e-6:
            if abs(rad[0] - p.r_so) < 1e-6 or abs(rad[0] - p.r_shaft) < 1e-6:
                bnd.setdefault("A0", []).append(c)
            elif abs(rad[0] - r_gap) < 1e-6:
                surf_up = gmsh.model.getAdjacencies(1, c)[0]
                bnd.setdefault(f"sliding_{regions[owner[surf_up[0]]].side}", []).append(c)
        elif all(abs(a) < 1e-6 for a in ang):
            bnd.setdefault("edge_0", []).append(c)
        elif all(abs(a - p.sector_deg) < 1e-6 for a in ang):
            bnd.setdefault("edge_sector", []).append(c)
    for name, tags in bnd.items():
        gmsh.model.addPhysicalGroup(1, tags, name=name)
    gmsh.option.setNumber("Mesh.MeshSizeMax", size_max)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)      # what ElmerGrid reads
    gmsh.model.mesh.generate(2)
    gmsh.write(path)
    n_nodes = len(gmsh.model.mesh.getNodes()[0])
    gmsh.finalize()
    return n_nodes, {k: len(v) for k, v in bnd.items()}


def to_solidworks_equations(p, path):
    """Global variables for a SolidWorks part whose equations are linked to this file."""
    rows = {"Magnet thickness": f"{p.mag_t}mm", "Magnet short side": f"{p.mag_w}mm", "Inner diameter": f"{2 * p.r_si}mm",
            "Stator OD": f"{2 * p.r_so}mm", "Slots": p.slots, "Slot opening": f"{p.slot_open}mm",
            "Airgap": f"{p.airgap}mm", "V half angle": f"{p.v_half_deg}deg", "Centre rib": f"{p.rib}mm",
            "Outer bridge": f"{p.bridge}mm", "Apex radius": f"{p.apex_r}mm"}
    if p.pocket == "slot":
        rows.update({"Magnet centre x": f"{p.mag_x}mm", "Magnet centre y": f"{p.mag_y}mm", "Slot length": f"{p.slot_w}mm",
                     "Slot thickness": f"{p.slot_t}mm", "Slot fillet": f"{p.slot_fillet}mm", "Web": f"{p.web}mm",
                     "Barrier thickness": f"{p.bar_t}mm", "Barrier apex": f"{p.bar_apex}mm",
                     "Barrier": 0 if p.barrier == "none" else 1})
    with open(path, "w", encoding="utf-8-sig", newline="\r\n") as f:
        f.writelines(f'"{k}"= {v}\n' for k, v in rows.items())


def plot(p, regions, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colour = {"M-19 Steel": "#9aa5b1", "N45": "#d1495b", "Air": "#ffffff", "Coil": "#e8a33d", "<No Mesh>": "#ffffff"}
    phase_col = {"A": "#e8a33d", "B": "#4c9f70", "C": "#3d7ab8"}
    fig, ax = plt.subplots(figsize=(8, 8))
    vis = visible(regions)
    for r in regions:
        poly = vis[r.name]
        for g in getattr(poly, "geoms", [poly]):
            if g.is_empty:
                continue
            ax.fill(*g.exterior.xy, color=phase_col.get(r.circuit, colour[r.material]), lw=0.3, ec="#333")
            for h in g.interiors:
                ax.fill(*h.xy, color="white", lw=0.3, ec="#333")
        if r.material == "N45":
            ax.annotate("", xy=lin(r.label, pol(1, r.magdir), 2.2), xytext=lin(r.label, pol(1, r.magdir), -2.2),
                        arrowprops=dict(arrowstyle="->", lw=1.2))
        if r.circuit:
            ax.text(*r.label, f"{r.circuit}{'+' if r.turns > 0 else '-'}", ha="center", va="center", fontsize=6)
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title("motorgeo sector: 9 slots, 2 poles (arrows = magnetisation)")
    fig.savefig(path, dpi=160, bbox_inches="tight"); plt.close(fig)


# ----------------------------------------------------------------------------- 7. FEMM solves on the drawn model
# Starting FEMM (one private FEMM per process, below-normal priority) is femm_session.open_femm()'s job, not this module's.
class FemmModel:
    """Open a .fem written by to_femm and evaluate (theta_e, id, iq) points. Rotation = AGE inner angle, no remesh.

    Works on the FEMM that femm_session.open_femm() started in this process: the module-level femm.* calls reach it.
    """

    def __init__(self, p, path):
        import femm
        self.femm, self.p, self.path = femm, p, path
        self.axes, self.kw1 = phase_axes(p)
        femm.opendocument(path)

    def solve(self, theta_e, i_d, i_q, **kw):
        """One rotor position, currents given in d/q (amplitude-invariant inverse Park)."""
        a = [math.radians(theta_e + s) for s in (0, -120, 120)]
        return self.solve_abc(theta_e, *(i_d * math.cos(t) - i_q * math.sin(t) for t in a), **kw)

    def solve_abc(self, theta_e, i_a, i_b, i_c, gap_deg=(), keep_ans=None):
        """One rotor position, phase currents given directly. gap_deg: angles (stator frame) to sample the
        air-gap B at. keep_ans: copy the solution file there (for field plots)."""
        femm, p = self.femm, self.p
        theta_d0 = 90 / p.pole_pairs                    # d-axis of pole 0 as drawn (half a pole pitch)
        rotor_deg = theta_e / p.pole_pairs + self.axes["A"] - theta_d0
        femm.mi_modifyboundprop("AGE", 10, rotor_deg)   # rotate the rotor side of the air-gap element
        for ph, i in zip("ABC", (i_a, i_b, i_c)):
            femm.mi_modifycircprop(ph, 1, i)
        femm.mi_analyze(1)
        femm.mi_loadsolution()
        psi, retries = self._flux_linkage()
        th = [math.radians(theta_e + s) for s in (0, -120, 120)]
        park = lambda x: (2 / 3 * sum(v * math.cos(t) for v, t in zip(x, th)), -2 / 3 * sum(v * math.sin(t) for v, t in zip(x, th)))
        (psi_d, psi_q), (i_d, i_q) = park(psi), park((i_a, i_b, i_c))
        torque_age = femm.mo_gapintegral("AGE", 0)
        femm.mo_groupselectblock(1)
        torque_wst = p.n_sectors * femm.mo_blockintegral(22)
        gap = [femm.mo_getgapb("AGE", a) for a in gap_deg]
        femm.mo_close()
        if keep_ans:
            import shutil
            shutil.copyfile(self.path[:-4] + ".ans", keep_ans)
        return {"theta_e": theta_e, "id": i_d, "iq": i_q, "ia": i_a, "ib": i_b, "ic": i_c,
                "psi_a": psi[0], "psi_b": psi[1], "psi_c": psi[2], "psi_d": psi_d, "psi_q": psi_q,
                "T_age": torque_age, "T_wst": torque_wst, "T_dq": 1.5 * p.pole_pairs * (psi_d * i_q - psi_q * i_d),
                "rotor_deg": rotor_deg, "gap_br": [b[0] for b in gap], "gap_bt": [b[1] for b in gap],
                "psi_retries": retries}

    def _flux_linkage(self):
        """Phase flux linkages, whole machine, with a safety net against one known wrong answer: exactly 0.0.

        What happened on 2026-10-28: in validate_full.py, 11 of the 24 no-current slice rows had a flux linkage of exactly
        0.0 while their torque was right, and re-running them alone gave correct values, so it looked like FEMM
        misreading under load. It was not a FEMM or CPU problem (found 2026-11-04, and reproduced with FEMM 5.0). The
        no-load job and the 10 A job of the same angles had been given ONE pair of file names, so two FEMMs saved,
        solved and read the same .fem and .ans at the same moment, and one of them could be handed the other job's
        answer. validate_full.py now gives every job files of its own, and this loop has not fired since: the psi_retries
        column of every study CSV that has one (3,518 solves) is 0 with FEMM 5.0, and none returned 0.0.

        The loop stays because a wrong flux linkage that raises no error is the worst kind of result, and checking is
        cheap: if a value is exactly 0.0, read again, then reload the solution, then re-solve (at most 3 retries). It
        cannot repair a shared file (reading the other job's file again gives the other job's answer), so the real defence
        is the rule that two FEMMs never share a file name (learn/04_workers.py, parts 3 and 4: part 3 makes the accident
        happen, part 4 shows the private names that prevent it). Returns the values and how many retries it took."""
        femm, p = self.femm, self.p
        for retries in range(4):
            psi = [p.n_sectors * femm.mo_getcircuitproperties(ph)[2] for ph in "ABC"]
            if all(x != 0.0 for x in psi) or retries == 3:
                break
            if retries >= 1:                       # 1st retry: just read again; 2nd: reload; 3rd: re-solve
                femm.mo_close()
                if retries == 2:
                    femm.mi_analyze(1)
                femm.mi_loadsolution()
        return psi, retries
