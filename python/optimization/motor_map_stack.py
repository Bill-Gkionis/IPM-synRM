"""optimization.motor_map_stack: what a stack length changes besides the electromagnetics: the masses, the rotor's inertia, the magnets to buy, and the winding of each stack.

What it is, and why
    The motor map is asked for 1 to 4 rows of magnets (a 20, 40, 60 or 80 mm stack). The maps scale the electromagnetics with the stack (psi, torque and iron loss go with the
    length, motor_map.Machine); this file adds the rest of the table a buyer reads: the steel, the magnets and the copper, the rotor's moment of inertia, the number of magnets to
    buy, and the ways to wind a longer stack.

The windings of a stack
    "same"  the same winding on every stack (the base winding: the headline winding of the winding study). A longer stack links more flux (psi goes with the length), so the same
            winding reaches the voltage limit sooner: the corner speed falls as 1 / length, and the torque rises as the length, so the power at the corner is the same for every stack.
    "opt"   the winding the winding study chose for each stack (OPTIMIZED_WINDINGS: optimization/winding_study.py, Storage 0acd770; 48 V, the shaft fan held to 50 W at the top speed,
            N45 at 80 C, class F, compact bundles, at most 8 strands in hand): 26 turns 4 paths 8 x 0.315 mm on 20 mm, 13 / 4 4 x 0.63 on 40 mm, 9 / 4 7 x 0.56 on 60 mm, 7 / 4 7 x 0.63 on
            80 mm. Real windings, each the smallest inverter current that still holds the continuous torque to 3,300 rpm on 48 V.
    "hold"  (an idealised family, kept for a stack the table does not have) the turns re-scaled to hold the corner speed of the anchor stack (20 mm): N' = N x anchor / stack turns of
            wire of stack / anchor times the copper area, so the copper in the slot and the fill do not change; a fraction of a turn and a wire of any diameter are allowed.

Masses and the inertia
    steel      the stator's and the rotor's visible steel areas (motorgeo.visible, one sector times the sectors) x the stack x 7,650 kg/m^3 x the lamination fill 0.98 (the maker's
               density and FEMM's lamination fill).
    magnets    the BOUGHT magnet (9.08 x 1.57 x 19.14 mm, 16 a row: 8 poles x 2) at 7,600 kg/m^3 (the N45 datasheet's value, as thermal.py); the FEMM model fills the slot, which is
               14.7 % more volume (reported). The magnets to buy: 16 x rows.
    copper     36 coils x N turns x n_s strands x S x the mean turn (2 L + 2 l_end) x 8,900 kg/m^3.
    inertia    J = integral of r^2 dm: the rotor's steel (the polygon's second moment about the axis, exact for the polygon: the formula of the sum over its edges, holes subtracted),
               the bought magnets at the radius of the modelled slots, and the shaft as a solid cylinder over the stack (15 mm, 7,850 kg/m^3). The shaft beyond the stack, the fan,
               the bearings' inner rings and any sleeve are NOT in it: the rotor-dynamics note has the whole rotor.

Run
    uv run python -m optimization.motor_map_stack       the numbers of both windings for the four stacks

Numbers to expect (this PC, 2026-11-09; copper_loss.py's full loss, the strands in a cluster, the diamond end turn)
    steel, magnets and inertia (the same for both windings)
    20 mm   16 magnets (33 g), steel 0.54 + 0.89 kg (stator + rotor), rotor inertia 0.96 g m^2
    40 mm   32 magnets (66 g), steel 1.08 + 1.78 kg, inertia 1.92 g m^2
    60 mm   48 magnets (100 g), steel 1.62 + 2.66 kg, inertia 2.87 g m^2
    80 mm   64 magnets (133 g), steel 2.15 + 3.55 kg, inertia 3.83 g m^2
    'same', 9 turns 4 paths 7 x 0.56 mm (fill 0.428), on 20 / 40 / 60 / 80 mm: R(100 C) 12.16 / 15.71 / 19.26 / 22.81 mOhm, copper 0.68 / 0.88 / 1.08 / 1.28 kg, active parts 2.14 / 3.80 / 5.46 / 7.12 kg
    'opt', the winding study's: 26 / 4 8 x 0.315 mm (fill 0.448) R 97.18 mOhm, copper 0.71 kg, active 2.17 kg on 20 mm; 13 / 4 4 x 0.63 mm (0.448) 31.38 mOhm, 0.92 kg, 3.84 kg on 40 mm;
    9 / 4 7 x 0.56 mm (0.428) 19.26 mOhm, 1.08 kg, 5.46 kg on 60 mm; 7 / 4 7 x 0.63 mm (0.422) 14.02 mOhm, 1.26 kg, 7.10 kg on 80 mm. The magnets the FEMM model holds are 14.7 % more than
    the bought ones (it fills the slot); the polygon formula of the inertia is checked against the closed form of a ring (tests/test_motor_map.py, 5e-7).

Self-check (the answers follow).
  1. Why does the same winding give the same power at the corner on every stack?
  2. In the "hold" family the turns are held by N' = N x 20 / stack. Why does the current in amperes rise with the stack, and what does the inverter see?
Answers.
  1. The torque at the current limit goes with the stack and the corner speed against it (the voltage at a speed goes with psi, which goes with the stack): their product, the
     power at the corner, does not change. A longer stack buys torque at a lower corner, not power.
  2. Fewer turns of thicker wire carry the same current density, so more amperes, at a voltage that is the same at the same speed: the inverter's current, and so its size, goes
     with the stack; its voltage does not.
"""
from __future__ import annotations

import math

import numpy as np

import motorgeo as mg
from optimization import motor_map as mm
from optimization.iron_loss import LAM_FILL, RHO_STEEL

RHO_MAGNET = 7600.0                    # kg/m^3, N45 (the datasheet value thermal.py uses)
RHO_COPPER = 8900.0                    # kg/m^3 (thermal.py's value for C11000)
RHO_SHAFT = 7850.0                     # kg/m^3, C45E
MAGNET_MM = (9.08, 1.57, 19.14)        # the bought magnet: width, thickness (the magnetised direction), length
MAGNETS_PER_ROW = 16                   # 8 poles x 2
ANCHOR_MM = 20.0
STACKS_MM = (20.0, 40.0, 60.0, 80.0)


def polar_moment_mm4(geom) -> float:
    """The integral of (x^2 + y^2) dA [mm^4] of a shapely Polygon or MultiPolygon about the origin, holes subtracted: for each ring the sum over its edges of
    (x_i^2 + x_i x_j + x_j^2 + y_i^2 + y_i y_j + y_j^2)(x_i y_j - x_j y_i) / 12 (exact for the polygon), the exterior counter-clockwise and the holes clockwise."""
    from shapely.geometry.polygon import orient
    total = 0.0
    for part in getattr(geom, "geoms", [geom]):
        poly = orient(part, 1.0)
        for ring in [poly.exterior, *poly.interiors]:
            xy = np.asarray(ring.coords)
            x, y = xy[:-1, 0], xy[:-1, 1]
            x1, y1 = np.roll(x, -1), np.roll(y, -1)
            total += float(((x * x + x * x1 + x1 * x1) + (y * y + y * y1 + y1 * y1)).dot(x * y1 - x1 * y)) / 12.0
    return total


def geometry_numbers(m: mm.Machine) -> dict:
    """Masses [kg], the rotor's inertia [kg m^2] and the magnets to buy for this machine (module docstring)."""
    p = m.p
    regions = mg.build(p)
    vis = mg.visible(regions)
    n = p.n_sectors
    L = m.L
    a_stator, a_rotor = n * vis["stator_steel"].area, n * vis["rotor_steel"].area               # mm^2
    steel_rho = RHO_STEEL * LAM_FILL
    rows = m.stack_mm / mm.ROW_MM
    n_mag = MAGNETS_PER_ROW * rows
    v_each = MAGNET_MM[0] * MAGNET_MM[1] * MAGNET_MM[2] * 1e-9
    mag_names = [r.name for r in regions if r.name.startswith("magnet_")]
    a_slot = n * sum(vis[k].area for k in mag_names)                                              # mm^2 of N45 the FEMM model has in one 20 mm row's cross-section
    slot_w, slot_t = p.slot_w, p.slot_t
    ratio = MAGNET_MM[0] * MAGNET_MM[1] / (slot_w * slot_t)                                       # the bought magnet's share of the modelled slot's cross-section
    j_steel = steel_rho * L * n * polar_moment_mm4(vis["rotor_steel"]) * 1e-12
    j_mag = RHO_MAGNET * L * ratio * n * sum(polar_moment_mm4(vis[k]) for k in mag_names) * 1e-12
    m_shaft = RHO_SHAFT * math.pi * (p.r_shaft * 1e-3) ** 2 * L
    j_shaft = 0.5 * m_shaft * (p.r_shaft * 1e-3) ** 2
    l_mt = m.copper.mean_turn_m(m) if hasattr(m.copper, "mean_turn_m") else 2.0 * L + 2.0 * 48.55e-3
    v_cu = 3 * mm.COILS_PER_PHASE * m.w.turns * m.w.copper_area_m2 * l_mt
    return {"stack_mm": m.stack_mm, "rows": rows, "m_stator_steel": steel_rho * a_stator * 1e-6 * L, "m_rotor_steel": steel_rho * a_rotor * 1e-6 * L,
            "n_magnets": n_mag, "m_magnets": n_mag * v_each * RHO_MAGNET, "m_magnets_femm": RHO_MAGNET * a_slot * 1e-6 * L, "m_copper": v_cu * RHO_COPPER, "m_shaft_active": m_shaft,
            "J_rotor": j_steel + j_mag + j_shaft, "J_rotor_steel": j_steel, "J_rotor_magnets": j_mag, "J_rotor_shaft": j_shaft,
            "m_active": steel_rho * (a_stator + a_rotor) * 1e-6 * L + n_mag * v_each * RHO_MAGNET + v_cu * RHO_COPPER, "mean_turn_m": l_mt}


# the winding study's choice at each stack (optimization/winding_study.py, Storage 0acd770, results.json "stacks"): turns a coil side, paths, strands, wire [mm], wire over the enamel [mm]
OPTIMIZED_WINDINGS = {20.0: mm.Winding(26, 4, 8, 0.315, 0.367), 40.0: mm.Winding(13, 4, 4, 0.63, 0.7042), 60.0: mm.Winding(9, 4, 7, 0.56, 0.63), 80.0: mm.Winding(7, 4, 7, 0.63, 0.7042)}


def stack_winding(base: mm.Winding, stack_mm: float, variant: str, table: dict | None = None, anchor_mm: float = ANCHOR_MM) -> mm.Winding:
    """The winding on a stack: "same" (the base winding), "opt" (the table's winding for this stack, OPTIMIZED_WINDINGS by default) or "hold" (turns x anchor / stack, the copper area per
    turn x stack / anchor)."""
    if variant == "same":
        return base
    if variant == "opt":
        t = OPTIMIZED_WINDINGS if table is None else table
        if float(stack_mm) not in {float(k) for k in t}:
            raise ValueError(f"no optimized winding for a {stack_mm:g} mm stack: the table has {sorted(t)}")
        return next(v for k, v in t.items() if float(k) == float(stack_mm))
    if variant == "hold":
        return base.scaled(anchor_mm / stack_mm)
    raise ValueError("variant: 'same', 'opt' or 'hold'")


def main() -> int:
    copper = mm.FullCopper(arrangement="cluster")
    for variant in ("same", "opt"):
        print(f"the headline winding on every stack ('same')" if variant == "same" else "the winding study's winding of each stack ('opt')")
        for stack in STACKS_MM:
            w = stack_winding(OPTIMIZED_WINDINGS[60.0], stack, variant)
            m = mm.Machine(w, stack, copper=copper)
            g = geometry_numbers(m)
            print(f"  {stack:3.0f} mm: {w.turns:6.2f} turns, {w.paths} paths, {w.strands:g} x {w.wire_mm:g} mm, fill {m.fill:.3f}, R(100 C) {1e3 * float(m.resistance(100.0)):.2f} mOhm; steel {g['m_stator_steel']:.2f} + "
                  f"{g['m_rotor_steel']:.2f} kg, magnets {g['m_magnets'] * 1e3:.0f} g ({g['n_magnets']:.0f}), copper {g['m_copper']:.2f} kg, active {g['m_active']:.2f} kg, J {g['J_rotor'] * 1e3:.3f} g m^2")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
