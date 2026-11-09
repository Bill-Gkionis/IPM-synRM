"""optimization.thermal: the lumped thermal network (LPTN) of the motor: how much loss it can shed, and so how much current it can carry.

What it is, and why
    Every watt lost in the motor ends up as heat, and the heat has to walk out through copper, enamel, varnish, steel, an air gap and a
    housing to the air or a coolant. Temperature, not torque, sets what current the motor can carry for an hour (the continuous rating)
    and for a minute (the short-term one). Two parts have a hard ceiling: the winding's insulation (class F 155 C, class H 180 C at its
    hottest spot) and the N45 magnets (80 C working temperature, a design decision: a hard limit). This module is the model that says
    which of the two is reached first, and at what loss. It takes the losses as INPUTS (the copper loss, the iron loss map and the
    magnets' eddy loss come from other modules) and returns the temperature of every part. It is fast (a fraction of a millisecond a
    call), because the winding study and the optimizer will call it tens of thousands of times.

The picture (a node is a part at one temperature; every line is a thermal resistance R [K/W]; 1/R is a conductance G [W/K])

      copper in the slot (W) --axial copper-- end turns (E) --h_end-- end-space air (A) --h_end-- end bells (B) -- flange -- housing (H)
         |   \\ liner, slot                       |                      |                          |  \\                     |
         |    \\                                  +------ rotor ends ----+                          |   ambient             ambient
       tooth (T) --- yoke (Y) --- contact --------------------------------------------------------- +                      or coolant
         |
       air gap (speed dependent)
         |
       rotor steel (R) -- magnets (M);   rotor steel -- shaft (S) -- bearings (G) -- end bells (B)

    Eleven nodes W E T Y H B A R M S G, plus the fixed ambient (and the coolant inlet of the liquid jacket). A node is the MEAN temperature
    of its part (the rotor node sits in the magnets' zone). Two outputs are not nodes: the hot spot in the slot (above W) and the nose of the
    end turn (above E). The cooling attaches at H (the cylinder over the stack) and B (the end bells and the cylinder over the end spaces).

The maths, each step built from the one before
    1. Heat flows like current. A temperature difference dT is the voltage and a heat flow Q the current:
           conduction  R = t / (k A)        a layer of thickness t, area A, conductivity k     (a pipe wall: ln(r2/r1) / (2 pi k L))
           convection  R = 1 / (h A)        h the heat transfer coefficient of the surface
           radiation   h_rad = eps sigma (Ts^2 + Tinf^2)(Ts + Tinf)    (the T^4 law, linearised, added to h)
       Parts store heat in a capacitance C = rho c V [J/K].
    2. The network. At each node, heat stored = heat made + heat arriving:   C_i dT_i/dt = P_i + sum_j (T_j - T_i) / R_ij.
       In steady state the left side is zero and the eleven equations are linear:  G T = P + G_b T_b  (G the conductance matrix, T_b the
       ambient). numpy solves the 11 x 11 system in microseconds; the outside convection depends on temperature, so it is iterated.
    3. The slot (the biggest unknown: wire, enamel, varnish and air inside a slot lined with paper). The wire bundle is replaced by one
       material, the Hashin-Shtrikman (Maxwell) lower bound for round conductors in a matrix:  k_eq = k_m (1 + f)/(1 - f), f the copper's share
       of the block inside the liner and k_m the matrix (enamel + varnish + voids, 0.22 W/(m K); an input). The heat made in the block
       conducts out through its walls: a Poisson problem, lap T = -q/k, in a rectangle. With the slot's equivalent rectangle (width a, depth h)
       and its top wall adiabatic (it faces the slot opening), mirror it about the top: a rectangle a x 2h with T = 0 on all four sides, whose
       solution is the double sine series
           T(x,y) = sum_{m,n odd} 16 q / (k pi^2 m n lam_mn) sin(m pi x/a) sin(n pi y/b),    lam_mn = (m pi/a)^2 + (n pi/b)^2 .
       The mean over the block is 64 q/(k pi^4) sum 1/(m^2 n^2 lam_mn); the hot spot, at the centre, is 16 q/(k pi^2) sum (+-)/(m n lam_mn); the share
       of the heat leaving through the long walls (the tooth flanks) and through the bottom (the yoke) follows from the wall gradients. They give
       R_mean (winding mean to the wall), R_hs (hot spot above the mean) and the two paths. The liner (Nomex, 0.25 mm) is in series on each wall.
       Checks (selfcheck): a square gives T_max = 0.07367 q a^2/k (textbook) and a long strip q a^2/(8k).
    4. The tooth is a fin: heat enters along its flank and flows radially to the yoke. Its width is w(r) = c_t r (a wedge), so the
       mean-to-root resistance is an integral of 1/(k A), done once when the model is built. Likewise the root (where the slot's round
       bottom closes) and the yoke ring, split in two halves around the yoke node: R = ln(r2/r1) / (2 pi k L) / 2 each.
    5. The copper loss falls on the temperature it makes. Copper's resistivity rises 0.393 %/K, P = P20 [1 + alpha (T - 20)], LINEAR in T, so it goes
       straight into the matrix (alpha P20 is taken off the diagonal): no iteration. A winding that would run away (alpha P20 R_th >= 1) is
       reported. Any other copper loss given as a function f(T_copper) is iterated to convergence.
    6. The air gap. Nusselt against the modified Taylor number, from Howey, Childs and Holmes (IEEE Trans. Ind. Electron. 59(3), 2012; open text
       read in this session), after Becker and Kaye (1962):
           Ta = w^2 r_m d^3 / nu^2,  F_g = pi^4 / (1697 S (1 - d/(2 r_m))),  S = 0.0571 (1 - 0.652 d/r_i) + 0.00056 (1 - 0.652 d/r_i)^-1,  Ta_m = Ta/F_g
           Ta_m < 1700          Nu = 2 (d/r_i) / ln(1 + d/r_i)         (laminar: conduction)
           1700 < Ta_m < 1e4    Nu = 0.128 Ta_m^0.367                   (Taylor vortices)
           1e4  < Ta_m < 1e7    Nu = 0.409 Ta_m^0.241                   (turbulent vortices)
       h = Nu k_air / (2 d) on the log-mean annulus area 2 pi r_lm L. The rotor's speed enters through w: the gap conducts better as it spins.
    7. The end space (end turns, rotor end faces and end bells face each other) is one air node. Every surface on it exchanges heat at
       h_end = 13.29 + 1.693 v [W/(m^2 K), v the rotor's peripheral speed in m/s]: the measured fit of Tovar-Barranco, Lopez-de-Heredia, Villar and
       Briz ("Modeling of end-space convection heat-transfer for internal and external rotor PMSMs with fractional-slot concentrated windings",
       IEEE Trans. Ind. Electron.; volume and pages: check before citing) for a smooth-rotor internal-rotor PMSM tested to about 7 m/s. Ours
       reaches 14 m/s at 3,000 rpm and 77 m/s at the top speed, so it is EXTRAPOLATED, and capped at h_end_max.
    8. The outside, three ways. (a) Natural convection plus radiation: Churchill and Chu's correlation for a horizontal cylinder (the end discs:
       a vertical plate), iterated because h depends on the surface temperature. (b) Forced air along axial fins: duct flow between the fins
       (laminar Nu below Re 2300, Gnielinski above 4000), fin efficiency tanh(mH)/(mH). (c) A liquid jacket: a rectangular water-glycol channel
       (laminar Nu, or Gnielinski), the wall's conduction in series, and the coolant's own warming (bulk mean).
    9. The rotor loses heat only through the air gap (to the stator) and through the shaft, the bearings and the end space. It is the worst-cooled
       part and the magnets sit in it: a magnet node hangs on the rotor steel through the adhesive's 0.1 mm.
   10. Hot spots. T_hot(slot) = T_W + R_hs Q_walls + the tooth's own gradient; the end turn's nose = T_E + conduction along the copper and across
       the bundle. The winding's hot spot is the higher of the two.
   11. The limits. Each limit has its own current density (every temperature rises with J), found by bisection on the full solution; the sustainable
       current is the smaller one, and the limit that sets it is the BINDING one.
   12. Transients. With heat capacities the system is C dT/dt = -G T + P. Over a step it is solved exactly for frozen coefficients,
       T = T_ss + expm(-C^-1 G dt) (T - T_ss), which cannot go unstable however short a time constant (the end-space air's is milliseconds).
   13. The shaft fan (Cooling.shaft_fan(v_ref, n_ref): a fan on the shaft over the finned housing, like a simple TEFC motor). The air speed in the 36 fin channels is
       v(n) = v_ref n / n_ref, so the cooling follows the speed and at standstill only natural convection and radiation are left. The forced conductance is the preset's
       (Gnielinski / laminar duct, fin efficiency) limited by the air stream's own capacity: the air warms as it flows, so the surfaces can give only
       G_f = m cp (1 - exp(-G_raw / (m cp))) to air at the inlet temperature (effectiveness-NTU, an isothermal wall; m cp = rho v A_flow cp, shared between the housing and the
       bell cylinders by length). The still-air convection of the finned housing is the plain cylinder's (Churchill and Chu) times a gain (2.5, range 1.5-3.5, assumed), the
       radiation the plain cylinder's times 1.15; the two kinds of convection are blended, G = (G_f^3 + G_n^3)^(1/3) (Churchill-Usagi, check before citing). The fan itself: the
       pressure the air path asks for is the channels' friction plus the rest (inlet, exit, cowl) in dynamic heads, dp = (f L_ch / d_h + K) rho v^2 / 2, with the laminar
       f = 57 / Re, the Blasius f = 0.316 Re^-0.25 above Re 4000 and K = 3.5 (assumed); the air power is dp Q with Q = v A_flow; the shaft power is that over the fan's efficiency
       (0.20, assumed). dp goes as v^2 and Q as v, so P goes as n^3 (the fan laws; the friction factor's fall with Re bends it by 3 % over the range). The fan's power is a load on
       the SHAFT and heats the AIR, not the motor: it is not a source in the network; the warming of the stream by it, dT = P / (m cp), is carried by the convective surfaces
       (external_full). size_shaft_fan(budget) finds the v_ref whose power at the top speed is the budget.

The interface
    steady_state(losses, cooling, stack_mm, *, overrides=None, speed_rpm=3000, fill=0.45, ...) -> SteadyResult
        losses   Losses(copper=..., iron_teeth=..., ...). copper may be a number [W, the total], a function f(T_C) -> W of the copper temperature
                 (the total if all the copper were at T_C), or a CopperLoss (a current density: the DC loss rho(T) J^2 V_cu).
        cooling  Cooling.natural() | .forced_air() | .liquid()         stack_mm  20, 40, 60, 80 (rows of magnets, 20 mm each)
        fill     the copper fill of the gross slot (the winding study varies it); overrides: the design (default: the optimized one).
    transient(j_rms, cooling, stack_mm, t_end_s, *, start="cold"|"hot"|array) -> TransientResult    (every node against time)
    sustainable(cooling, stack_mm, *, wind_limit_c=155, magnet_limit_c=80, ...) -> dict (the continuous J, loss, ampere-turns, binding limit)
    overload(cooling, stack_mm, factor, j_cont, ...)    how long 1.5x or 2x the continuous current can be carried, from cold and from hot
    time_constants(cooling, stack_mm, j_rms)            the modes of the network (eigenvalues of C^-1 G)
    ThermalModel(...)  one built network: build once, call .steady_state(losses) many times (the fast path)
    ThermalModel.radial_2d(...)  the network reduced to a long stack (no end turns, no rotor): what the 2D finite-element check compares with.
    Cooling.shaft_fan(v_ref_ms=6.0, n_ref_rpm=3000)   a fan on the shaft: v = v_ref n / n_ref in the fin channels (kind "shaft_fan"; the other coolings are unchanged)
    ThermalModel.fan_state(v), .fan_power()           the air path at v: Re, friction, pressure, flow, air power, shaft power [W], the stream's warming
    shaft_fan_power(rpm, cooling, stack_mm)           the fan's shaft power [W] at rpm (scalar or array): the motor map's fan load, a function rpm -> W
    size_shaft_fan(budget_w, rpm_top=16330, n_ref_rpm=3000, stack_mm=60)   the fan whose power at the top speed is the budget: {"v_top", "v_ref", "cooling", ...}
    ThermalModel.heat_flows(res)     the heat flow [W] on every branch of a steady state, the exits (housing and bells, convection and radiation apart) and the groups
                                     (slot winding to teeth and yoke, end turns through the slot or into the air, the rotor across the gap, into the air or through the shaft)
    ThermalModel.sustainable(..., iron=function of j_rms | losses_at=function of j_rms -> Losses, j_hi=...)   a current-dependent iron loss, or the whole loss model
    ThermalModel.overload(..., iron=, iron_start=)   the iron loss during the overload and before it (a loss that depends on the current)
    ThermalModel.outside_parts(t_h, t_b)             the outside's conductances to the ambient, convection and radiation apart
    The motor's own loss and torque models (flux maps, iron-loss map, copper-loss model) drive the network in thermal_machine.py; the shaft-fan study is thermal_fan_study.py.

Example (the optimized design, 3,000 rpm, a 40 C ambient)
    from optimization import thermal
    m = thermal.ThermalModel(stack_mm=60, cooling=thermal.Cooling.liquid())
    r = m.steady_state(m.default_losses(copper=m.copper(10.0)))       # 10 A/mm^2 rms in the copper, the iron placeholder, no magnet eddy loss
    print(r.hot_spot, r.magnet, r.binding, r.copper_loss)               # C, C, "winding" or "magnet", W
    print(m.sustainable())                                              # the continuous J, copper loss, ampere-turns and the binding limit
    thermal.steady_state(thermal.Losses(copper=lambda t: 20.0 * (1 + 0.00393 * (t - 20.0)), iron_teeth=5.0), thermal.Cooling.natural(), 40)

Plugging in the other loss modules (iron_loss.py's keys teeth / yoke / rotor are Losses' iron_teeth / iron_yoke / iron_rotor; the magnets' eddy loss is Losses(magnet=...))
    The copper-loss module's function of the copper temperature goes in as the copper term, and its end turn (the diamond, 48.55 mm) is given to the network so that both use one length
    (a run of this, 60 mm forced air at 118.6 ampere-turns: 69.6 W converged in 5 iterations, hot spot 69.5 C, magnets 64.2 C):
        from optimization import copper_loss
        w, op = {"N": 26, "a": 4, "n_s": 8, "d_mm": 0.315, "d_over_mm": 0.367}, {"i_peak": 118.6 * 4 / 26, "gamma_deg": 22.47, "rpm": 3000.0, "f_sw": 20000.0, "v_dc": 48.0}
        m = thermal.ThermalModel(stack_mm=60, cooling=thermal.Cooling.forced_air(), assume=thermal.Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6))
        iron = m.iron_placeholder()
        r100 = copper_loss.losses(w, op, 100.0, 60.0)                        # only to read the slot's share of the copper loss
        share = r100["p_dc_slot"] / (r100["p_dc_slot"] + r100["p_dc_end"])
        r = m.steady_state(thermal.Losses(copper=lambda t: copper_loss.losses(w, op, t, 60.0)["p_total"], copper_slot_fraction=share,
                                          iron_teeth=iron["iron_teeth"], iron_yoke=iron["iron_yoke"], iron_rotor=iron["iron_rotor"]))

Run (from the Python to FEMM folder; no solver is started)
    uv run python -m optimization.thermal         the self-checks (building blocks against textbook cases), a speed test and a one-line summary
    The study that uses it, and the finite-element checks, are optimization/thermal_study.py and optimization/thermal_fem.py.
    The tests: uv run pytest tests/test_thermal.py tests/test_thermal_fem.py -q (about 2 s).

Through-stack conduction. The laminations conduct along the shaft 13 times worse than along the sheet (NREL: 1.68 against 21.9 W/(m K)). A stack with a
    uniform source is L / (12 k A) from its end faces to its mean, and the end faces then give the heat to the end-space air: two small links (stator iron
    to the air, rotor steel to the air). They carry well under 1 % of the stator's heat, which is why the radial path through the sheet decides.

Numbers to expect (this PC, 2026-11-09)
    selfcheck: all 16 checks pass. The optimized design's slot is 72.44 mm^2 (two layers of 36.218), the yoke 4.161 mm, the tooth 3.0-3.6 mm; the reference design's DC loss at 0.45
    fill, 100 C and 100 A-turns peak is 56.6 W (the winding study's number, mean turn 160.7 mm); the optimized design's mean turn is 2 L + 2 x 61.15 mm (the end turns: half a
    circle on the 36.38 mm span, plus 2 mm legs). The copper-loss module (copper_loss.py) draws a shorter, diamond-shaped end turn (48.55 mm a coil side, about 13.6 mm of overhang):
    Assumptions(end_turn_mm=48.55, end_overhang_mm=13.6) uses it (60 mm stack: J 1.73 / 6.96 / 15.36 A/mm^2 in natural convection / forced air / liquid, the magnets binding in all three).
    Speed (60 mm, one steady_state call): natural convection 0.36 ms from the default guess (it iterates) and 0.04 ms warm-started; forced air 0.07 ms; liquid 0.08 ms; a new fill
    on every call adds 0.01 ms; sustainable() (two bisections on the full solution) takes 2.5 to 7.5 ms; building a ThermalModel takes 0.2 ms once the design's geometry is cached (the
    first build of a design costs 0.1 s: the motorgeo drawing).
    The headline (60 mm stack, class F, 3,000 rpm, 40 C): natural convection J = 1.85 A/mm^2, 10.3 W of copper; forced air 6.88 A/mm^2, 148 W; liquid jacket 14.44 A/mm^2, 714 W; the
    magnets (80 C) are the binding limit in all three. thermal_study.py has every stack and class; thermal_fem.py and thermal_ccx.py the finite-element checks of the slot.
    The shaft fan (thermal_fan_study.py; 60 mm, the optimized design with the copper-loss model's end turn): sizing (a), 6 m/s at 3,000 rpm, takes 12.4 W of shaft power at 3,000 rpm,
    449 W at 10,000 rpm and 1,941 W at 16,330 rpm (P ~ n^3), and warms its air by 0.34 / 3.7 / 9.8 K; sizing (b), 50 W at 16,330 rpm, blows 1.76 m/s at 3,000 rpm and takes 0.32 W
    there. At 3,000 rpm (a) and the fixed 6 m/s preset differ by 0.2 % in continuous ampere-turns (radiation and mixing add, the air's warming takes away). The earlier coolings are
    unchanged to the last digit (tests/test_thermal_fan.py keeps twelve of their numbers).
    Speed of the new calls: fan_state 1 us; heat_flows 15 us; a steady_state of the shaft fan 0.15 ms; one continuous limit with the real losses
    (thermal_machine.continuous_limit: about 40 steady states with the copper-loss model as the copper function) 0.1 to 0.2 s.

Assumptions
    Every coefficient that is assumed, with its range, is in the Assumptions dataclass and ASSUMPTION_TABLE. Materials come from the materials
    notes (03-thermal.md and the README), NREL's report and the datasheets (MATERIALS gives the source of each). Anything marked "check before
    citing" was written from memory and was not looked up in this session.

Self-check (the answers follow).
  1. The magnets lose heat only through the air gap and the shaft. Why does a longer stack put them CLOSER to the stator's temperature?
  2. Copper's loss rises with temperature. Why does the network not diverge, and when would it?
  3. Why are the end turns' losses a problem for the SLOT, and why most of all in a short stack?
Answers.
  1. The gap's area grows with the stack (R_gap goes as 1/L), while the shaft and end-face paths do not: the rotor is tied more tightly to the
     stator than to the outside, so it follows the stator.
  2. Each kelvin of rise adds a fixed fraction (alpha) of the cold loss: the feedback gain is alpha P R_th, below 1 for any real winding. It
     diverges only if alpha P20 R_th >= 1 (about 254 K of cold-resistance rise); the solver then reports thermal runaway.
  3. Copper conducts along the wire: the end-turn heat enters the slot through the wire itself and leaves through the slot's walls with the
     slot's own heat. With a short stack the end turns hold most of the copper (75 % at 20 mm) and a short slot has to pass all of it.
"""
from __future__ import annotations

import functools
import json
import math
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy import linalg

import motorgeo as mg
import paths

# ---------------------------------------------------------------------------------------------------- constants
SIGMA_SB = 5.670374419e-8            # W/(m^2 K^4)  Stefan-Boltzmann
G_ACC = 9.80665                      # m/s^2
P_ATM = 101325.0                     # Pa
R_AIR = 287.05                       # J/(kg K)
K_ABS = 273.15
NODES = ("W", "E", "T", "Y", "H", "B", "A", "R", "M", "S", "G")
IDX = {name: i for i, name in enumerate(NODES)}
N_NODES = len(NODES)
NODE_NAMES = {"W": "winding in the slot", "E": "end turns", "T": "stator teeth", "Y": "stator yoke", "H": "housing", "B": "end bells",
              "A": "end-space air", "R": "rotor steel", "M": "magnets", "S": "shaft", "G": "bearings"}
CLASS_F_C, CLASS_H_C = 155.0, 180.0   # IEC 60085 hot-spot limits [C]
MAGNET_LIMIT_C = 80.0                 # N45's working temperature (materials notes, Arnold datasheet, Eclipse): a hard limit
AMBIENT_C = 40.0                      # IEC 60034-1 standard ambient; also the coolant inlet of the jacket
MAX_SPEED_RPM = 16_330                # from the stress check (the maximum operating speed)
RATED_RPM = 3_000

# materials: value, source (W/(m K), kg/m^3, J/(kg K))
MATERIALS = {
    "copper": dict(k=391.0, rho=8900.0, c=394.0, src="Aurubis C11000 (03-thermal.md)"),
    "aluminium": dict(k=154.0, rho=2700.0, c=943.0, src="NIST 6061-T6 fits at 293 K; Smiths density (README)"),
    "m19": dict(k_in=21.9, k_ax=1.68, rho=7650.0, c=463.0, src="NREL 2015 table 3 (k, c; k_ax at 138 kPa); AK Steel density"),
    "n45": dict(k=6.2, rho=7600.0, c=460.0, src="Arnold 2021 datasheet (the lower of 6.2 / 6.7)"),
    "shaft": dict(k=41.7, rho=7850.0, c=461.0, src="thyssenkrupp C45E (README)"),
    "bearing_steel": dict(k=30.0, rho=7830.0, c=460.0, src="100Cr6 data sheets (README: 20-33)"),
    "liner": dict(k=0.139, rho=960.0, c=1090.0, src="Nomex 410, 0.25 mm at 150 C (README); c = Kapton's, assumed"),
    "varnish": dict(k=0.22, rho=1200.0, c=1090.0, src="Voltatex 4200 cured (README: 0.21-0.23, 1,150-1,250); c assumed"),
}
RHO_CU_20 = 1.7241e-8                 # ohm m (IEC 60028), as the winding study
ALPHA_CU = 0.00393                    # 1/K (+0.393 %/K), as the winding study


# ---------------------------------------------------------------------------------------------------- assumptions
@dataclass(frozen=True)
class Assumptions:
    """Every number that is ASSUMED (the housing, the shaft and the bearings are not designed yet), with its range.
    Change one with `Assumptions(contact_factor=5)` and pass it as ThermalModel(assume=...). ASSUMPTION_TABLE lists them with ranges."""
    liner_mm: float = 0.25               # slot liner (Nomex 410); range 0.15-0.4
    k_matrix: float = 0.22               # W/(m K): varnish + enamel + voids around the wire; range 0.15 (air voids) - 1.0 (filled resin)
    end_leg_mm: float = 2.0              # straight leg of the end turn at each end (the winding study's e)
    end_turn_mm: float | None = None     # length of ONE end turn [mm] (per coil side, both ends together); None = the winding study's half circle on the span plus two
                                         # legs (61.15 mm for the optimized design). The copper-loss module's diamond end turn is 48.55 mm
    end_overhang_mm: float | None = None # how far the end winding sticks out beyond the stack [mm]; None = the half circle's radius plus a leg (20.2 mm); a diamond: about 13.6
    end_clearance_mm: float = 5.0        # end winding to end cap; range 3-15
    contact_nrel_mm2kw: float = 193.0    # NREL 2015: lamination-to-lamination contact resistance, 138 kPa, 29 gauge M19 [mm^2 K/W]
    contact_factor: float = 2.0          # x NREL: the stack's punched OD against a machined bore; range 1 (NREL as measured) - 8
    housing_wall_mm: float = 5.0         # aluminium housing wall (plain and finned); range 3-10
    flange_mm2kw: float = 500.0          # bolted joint between end bell and housing [mm^2 K/W]; range 100-3000
    bell_mm: float = 8.0                 # end-bell thickness; range 5-12
    magnet_gap_mm: float = 0.1           # magnet to slot clearance (filled with adhesive); range 0.05-0.3
    k_adhesive: float = 0.5              # W/(m K) magnet adhesive; range 0.2 (unfilled) - 1.5 (filled); air would be 0.03
    rotor_k_fraction: float = 0.8        # rotor's radial conduction vs solid steel (the magnet pockets); range 0.5-1.0
    rotor_shaft_mm2kw: float = 500.0     # rotor core to shaft fit contact [mm^2 K/W]; range 100-2,000
    bearing_r_kw: float = 2.0            # one bearing and its seat, inner ring to end bell [K/W]; range 0.5-5
    bearing_mass_kg: float = 0.04        # one bearing (6002 size), for the capacity
    shaft_stub_mm: float = 8.0           # shaft from the end space to the bearing's centre (beyond the end bell's inside face); range 4-20
    h_end_k1: float = 13.29              # W/(m^2 K): Tovar-Barranco et al., internal-rotor PMSM (smooth rotor, no fan)
    h_end_k2: float = 1.693              # W/(m^2 K) per (m/s)
    h_end_max: float = 100.0             # W/(m^2 K): cap on the extrapolated correlation (their data stop near 7 m/s)
    end_area_factor: float = 1.0         # x the end winding's exposed-area estimate; range 0.6-1.4
    t_gap_film_c: float = 70.0           # C: air temperature for the gap and end-space properties; range 50-100
    wind_gap_factor: float = 1.0         # x smooth-cylinder windage; slotted stators run 1.2-2
    bearing_mu: float = 0.0015           # friction coefficient of a deep-groove ball bearing (SKF's simple estimate; check before citing)
    bearing_load_n: float = 60.0         # N on each bearing: half the rotor's weight plus an axial spring preload (assumed)
    # the iron-loss PLACEHOLDER (provisional until the loss map arrives)
    b_tooth_t: float = 1.5               # T peak in the teeth (rough; the tooth is 3.0 mm at its narrowest); range 1.3-1.8
    b_yoke_t: float = 1.2                # T peak in the yoke; range 1.0-1.6
    rotor_iron_fraction: float = 0.10    # rotor iron loss as a share of the stator's (slot harmonics only); range 0-0.5
    # outside, natural convection
    emissivity: float = 0.85             # painted or anodised aluminium; bare aluminium is about 0.1
    natural_fin_gain: float = 1.0        # area gain of a LIGHTLY finned housing in still air (1.0 = plain); range 1.0-1.8
    # outside, forced air over axial fins
    n_fins: int = 36
    fin_height_mm: float = 15.0
    fin_thickness_mm: float = 3.0
    air_speed_ms: float = 6.0            # in the fin channels (an external fan); range 2-12
    t_air_film_c: float = 55.0           # C: film temperature for the forced-air properties
    # outside, the shaft fan (Cooling.shaft_fan): the same fins and channels, the air driven by a fan on the shaft, v proportional to the speed
    duct_f_re: float = 57.0              # Darcy friction factor of a laminar fin channel is duct_f_re / Re (rectangle of aspect ratio 1.4: 57; range 54-64; check before citing)
    fan_k_system: float = 3.5            # pressure loss of the whole air path in dynamic heads rho v^2 / 2 of the channel, besides the channel's own friction: inlet contraction
                                         # 0.5 + exit 1.0 + cowl, fan hub and turning 2.0; range 2-6 (assumed; check before citing)
    fan_efficiency: float = 0.20         # air power / shaft power of a simple radial-blade shaft fan in its cowl; range 0.10-0.35 (assumed; check before citing)
    fan_fin_gain_conv: float = 2.5       # still-air natural convection of the finned housing as a multiple of the plain cylinder's (Churchill-Chu on the housing's diameter);
                                         # the fins add 3.5 x the wetted area but shade each other; range 1.5-3.5 (assumed: no still-air data for this housing)
    fan_fin_gain_rad: float = 1.15       # radiation of the finned housing as a multiple of the plain cylinder's (the fins see each other: about the envelope's area); range 1.0-1.3
    mixed_exponent: float = 3.0          # Churchill-Usagi blend of forced and natural convection, G = (G_forced^n + G_natural^n)^(1/n), n = 3 (check before citing); range 2-4
    air_bulk_warming: bool = True        # the air warms as it flows along the fins: the forced conductance cannot exceed the stream's m cp (shaft fan only)
    fan_air_preheat: bool = True         # the fan's whole shaft power ends as heat in the air before the fins (about 95 % does: the exit's dynamic-head loss comes after them)
    # outside, liquid jacket (50/50 water-glycol)
    jacket_wall_mm: float = 3.0          # aluminium between the stator's OD and the coolant
    jacket_channel_w_mm: float = 8.0
    jacket_channel_h_mm: float = 6.0
    jacket_rib_mm: float = 6.0           # land between channel turns
    jacket_outer_wall_mm: float = 3.0
    jacket_rib_efficiency: float = 0.9   # the lands as fins
    flow_lpm: float = 4.0                # coolant flow [L/min]; range 1-10
    coolant_in_c: float | None = None    # coolant inlet [C]; None = the cooling's ambient (40 C: the brief)
    coolant_rho: float = 1060.0          # kg/m^3  50/50 ethylene glycol-water near 45 C (check before citing)
    coolant_cp: float = 3400.0           # J/(kg K) (check before citing)
    coolant_k: float = 0.40              # W/(m K) (check before citing)
    coolant_mu: float = 2.7e-3           # Pa s (check before citing; about 3.4e-3 at 40 C and 2.4e-3 at 50 C)


ASSUMPTION_TABLE = [
    ("liner_mm", "slot liner thickness", "0.25 mm", "0.15-0.4", "Nomex 410 (README)"),
    ("k_matrix", "slot matrix conductivity (varnish, enamel, voids)", "0.22 W/(m K)", "0.15-1.0", "Voltatex 4200; measured random windings: 0.4-0.7 overall"),
    ("end_turn_mm / end_overhang_mm", "end turn per coil side, both ends / how far it sticks out", "61.15 / 20.2 mm", "48.55 / 13.6 (the copper-loss module's diamond) - 61.15 / 20.2", "the winding study's half circle on the span plus two 2 mm legs"),
    ("contact_factor", "stator-housing contact against NREL's lamination contact", "2 x 193 mm2 K/W", "1-8", "NREL 2015 (lamination faces; the punched OD is rougher)"),
    ("housing_wall_mm", "aluminium housing wall", "5 mm", "3-10", "assumed: no housing is designed"),
    ("flange_mm2kw", "end bell to housing joint", "500 mm2 K/W", "100-3000", "assumed"),
    ("rotor_shaft_mm2kw", "rotor to shaft fit", "500 mm2 K/W", "100-2000", "assumed"),
    ("bearing_r_kw", "one bearing and its seat", "2.0 K/W", "0.5-5", "assumed"),
    ("h_end_k1/k2", "end-space h = k1 + k2 v", "13.29 + 1.693 v", "x0.5-x2", "Tovar-Barranco et al. (PMSM, 100-1200 rpm; extrapolated)"),
    ("b_tooth_t / b_yoke_t", "iron-loss placeholder flux densities", "1.5 / 1.2 T", "1.3-1.8 / 1.0-1.6", "rough: PROVISIONAL"),
    ("rotor_iron_fraction", "rotor iron loss / stator iron loss", "0.10", "0-0.5", "assumed: PROVISIONAL"),
    ("emissivity", "housing surface", "0.85", "0.1-0.95", "painted or anodised"),
    ("air_speed_ms", "air speed in the fin channels", "6 m/s", "2-12", "an assumed external fan"),
    ("flow_lpm", "coolant flow", "4 L/min", "1-10", "assumed"),
    ("fan_efficiency", "shaft fan: air power / shaft power", "0.20", "0.10-0.35", "assumed (a simple radial-blade fan in a cowl); check before citing"),
    ("fan_k_system", "shaft fan: pressure loss of the air path besides the channels' friction, in dynamic heads", "3.5", "2-6", "assumed (inlet 0.5, exit 1.0, cowl and fan 2.0); check before citing"),
    ("duct_f_re", "fin channel: laminar Darcy friction factor x Re", "57", "54-64", "rectangular duct of aspect ratio 1.4 (from memory; check before citing)"),
    ("fan_fin_gain_conv / _rad", "still-air convection / radiation of the finned housing, as a multiple of a plain cylinder's", "2.5 / 1.15", "1.5-3.5 / 1.0-1.3", "assumed: no still-air data for this housing"),
    ("mixed_exponent", "forced and natural convection blended as (G_f^n + G_n^n)^(1/n)", "3", "2-4", "Churchill-Usagi form (from memory; check before citing)"),
    ("air_bulk_warming / fan_air_preheat", "the air warms along the fins (effectiveness-NTU) / the fan's loss warms it before the fins", "on / on", "off / off", "physics, switchable to show their size"),
]


# ---------------------------------------------------------------------------------------------------- air and heat transfer correlations
def air_props(t_c: float) -> dict:
    """Dry air at 1 atm and t_c [C]: density, dynamic and kinematic viscosity, conductivity, cp, Pr, beta (ideal gas). Sutherland's law
    for the viscosity and k = 0.02624 (T/300)^0.8646 (fits 300-400 K to 1 %): textbook values (check before citing)."""
    T = t_c + K_ABS
    rho = P_ATM / (R_AIR * T)
    mu = 1.458e-6 * T ** 1.5 / (T + 110.4)
    k = 0.02624 * (T / 300.0) ** 0.8646
    cp = 1007.0
    return {"rho": rho, "mu": mu, "nu": mu / rho, "k": k, "cp": cp, "Pr": mu * cp / k, "beta": 1.0 / T}


def h_radiation(eps: float, ts_c: float, tinf_c: float) -> float:
    """Linearised radiation coefficient eps sigma (Ts^2 + Tinf^2)(Ts + Tinf) [W/(m^2 K)] (Kelvin inside)."""
    ts, ti = ts_c + K_ABS, tinf_c + K_ABS
    return eps * SIGMA_SB * (ts * ts + ti * ti) * (ts + ti)


def h_natural_cylinder(d: float, ts_c: float, tinf_c: float) -> float:
    """Churchill and Chu (1975): a horizontal cylinder of diameter d [m] in still air, Ra <= 1e12,
    Nu = {0.60 + 0.387 Ra^(1/6) / [1 + (0.559/Pr)^(9/16)]^(8/27)}^2, film properties at (Ts + Tinf)/2. Returns h [W/(m^2 K)]."""
    dt = max(ts_c - tinf_c, 0.01)
    a = air_props(0.5 * (ts_c + tinf_c))
    ra = G_ACC * a["beta"] * dt * d ** 3 / (a["nu"] * (a["k"] / (a["rho"] * a["cp"])))
    nu = (0.60 + 0.387 * ra ** (1 / 6) / (1 + (0.559 / a["Pr"]) ** (9 / 16)) ** (8 / 27)) ** 2
    return nu * a["k"] / d


def h_natural_vplate(length: float, ts_c: float, tinf_c: float) -> float:
    """Churchill and Chu: a vertical plate of height `length` [m] (the end discs of a horizontal-shaft motor),
    Nu = {0.825 + 0.387 Ra^(1/6) / [1 + (0.492/Pr)^(9/16)]^(8/27)}^2. Returns h [W/(m^2 K)]."""
    dt = max(ts_c - tinf_c, 0.01)
    a = air_props(0.5 * (ts_c + tinf_c))
    ra = G_ACC * a["beta"] * dt * length ** 3 / (a["nu"] * (a["k"] / (a["rho"] * a["cp"])))
    nu = (0.825 + 0.387 * ra ** (1 / 6) / (1 + (0.492 / a["Pr"]) ** (9 / 16)) ** (8 / 27)) ** 2
    return nu * a["k"] / length


def nusselt_duct(re: float, pr: float, laminar_nu: float = 4.0) -> float:
    """Nusselt number in a duct: `laminar_nu` below Re = 2300, Gnielinski's turbulent correlation above 4000,
    Nu = (f/8)(Re - 1000) Pr / (1 + 12.7 sqrt(f/8)(Pr^(2/3) - 1)),  f = (0.79 ln Re - 1.64)^-2, and a straight blend between."""
    def turbulent(r):
        f = (0.79 * math.log(r) - 1.64) ** -2
        return (f / 8) * (r - 1000.0) * pr / (1 + 12.7 * math.sqrt(f / 8) * (pr ** (2 / 3) - 1))
    if re <= 2300:
        return laminar_nu
    if re >= 4000:
        return max(turbulent(re), laminar_nu)
    w = (re - 2300.0) / (4000.0 - 2300.0)
    return (1 - w) * laminar_nu + w * max(turbulent(4000.0), laminar_nu)


def duct_friction(re: float, f_re: float = 57.0) -> float:
    """Darcy friction factor of a fin channel: laminar f = f_re / Re up to Re = 2300 (f_re = 57 for a rectangle of aspect ratio 1.4), the smooth-pipe Blasius
    value f = 0.316 Re^(-1/4) from Re = 4000, and a straight blend between (the same thresholds as nusselt_duct). Check before citing. Zero at Re = 0."""
    if re <= 0.0:
        return 0.0
    if re <= 2300.0:
        return f_re / re
    turb_hi = 0.316 * 4000.0 ** -0.25
    if re >= 4000.0:
        return 0.316 * re ** -0.25
    w = (re - 2300.0) / (4000.0 - 2300.0)
    return (1.0 - w) * f_re / 2300.0 + w * turb_hi


def fin_efficiency(h: float, k: float, t: float, height: float) -> float:
    """Straight rectangular fin, tip-corrected: eta = tanh(m Hc)/(m Hc), m = sqrt(2 h / (k t)), Hc = H + t/2."""
    hc = height + t / 2
    m = math.sqrt(2.0 * h / (k * t))
    return math.tanh(m * hc) / (m * hc)


def air_gap_h(omega: float, r_rotor: float, gap: float, t_film_c: float = 70.0) -> dict:
    """Heat transfer coefficient across the air gap of a rotating smooth cylinder, no axial throughflow (Howey, Childs and Holmes 2012, after
    Becker and Kaye 1962). Returns {"h", "Nu", "Ta_m", "regime", "F_g"}. h = Nu k / (2 gap); the area to use is 2 pi r_lm L (r_lm the log-mean
    radius). Ta_m = Ta / F_g with Ta = w^2 r_m d^3 / nu^2 (the module docstring has F_g and S). A stationary rotor gets the laminar value."""
    a = air_props(t_film_c)
    r_o = r_rotor + gap
    r_m = 0.5 * (r_rotor + r_o)
    x = gap / r_rotor
    s = 0.0571 * (1 - 0.652 * x) + 0.00056 / (1 - 0.652 * x)
    fg = math.pi ** 4 / (1697.0 * s * (1 - gap / (2 * r_m)))
    ta = omega * omega * r_m * gap ** 3 / a["nu"] ** 2
    tam = ta / fg
    nu_lam = 2.0 * x / math.log(1 + x)
    if tam < 1700.0:
        nu, regime = nu_lam, "laminar (conduction)"
    elif tam < 1.0e4:
        nu, regime = 0.128 * tam ** 0.367, "Taylor vortices"
    else:
        nu, regime = 0.409 * tam ** 0.241, "turbulent vortices"
    nu = max(nu, nu_lam)
    return {"h": nu * a["k"] / (2.0 * gap), "Nu": nu, "Ta_m": tam, "regime": regime, "F_g": fg}


def windage_loss(omega: float, r_rotor: float, gap: float, length: float, t_film_c: float = 70.0, factor: float = 1.0) -> float:
    """Friction in the air gap, P = C_f pi rho w^3 r^4 L (Vrancik's form). C_f is the larger of the laminar Couette value 2/Re_d and the vortex
    values of Bilgen and Boulos, 0.515 (d/r)^0.3 / Re_d^0.5 (500 < Re_d < 1e4) and 0.0325 (d/r)^0.3 / Re_d^0.2 (Re_d > 1e4), with
    Re_d = w r d / nu (check before citing the constants). `factor` scales for a slotted stator (1.2-2)."""
    a = air_props(t_film_c)
    re = omega * r_rotor * gap / a["nu"]
    if re <= 0:
        return 0.0
    cf_lam = 2.0 / re
    if re < 500.0:
        cf = cf_lam
    elif re < 1.0e4:
        cf = max(cf_lam, 0.515 * (gap / r_rotor) ** 0.3 / math.sqrt(re))
    else:
        cf = max(cf_lam, 0.0325 * (gap / r_rotor) ** 0.3 / re ** 0.2)
    return factor * cf * math.pi * a["rho"] * omega ** 3 * r_rotor ** 4 * length


def bearing_loss(omega: float, assume: Assumptions) -> float:
    """Two ball bearings: M = 0.5 mu F d (the simple friction-torque estimate; check before citing), P = 2 M w, with the 15 mm bore."""
    return 2.0 * 0.5 * assume.bearing_mu * assume.bearing_load_n * 0.015 * omega


def end_space_h(v_tip: float, assume: Assumptions) -> float:
    """Equivalent heat transfer coefficient of the end space, h = k1 + k2 v (Tovar-Barranco et al.; v the rotor's peripheral speed), capped."""
    return min(assume.h_end_k1 + assume.h_end_k2 * v_tip, assume.h_end_max)


# ---------------------------------------------------------------------------------------------------- the slot: Poisson in a rectangle
def rect_poisson(a: float, b: float, terms: int = 400) -> tuple[float, float, float]:
    """Solve  lap T = -1  on the rectangle [0,a] x [0,b] with T = 0 on all four sides (k = q = 1): the double sine series
    T = sum_{m,n odd} c_mn sin sin, c_mn = 16 / (pi^2 m n lam_mn). Returns (mean T, T at the centre = the maximum, the share of the heat that
    leaves through the two walls of length b, the x = 0 and x = a ones). Lengths in any unit u: T comes out in u^2. The wall shares are
    normalised so that the two shares sum to 1 (the series for each converges slowly, and the sum is known)."""
    m = np.arange(1, 2 * terms, 2, dtype=float)[:, None]
    n = np.arange(1, 2 * terms, 2, dtype=float)[None, :]
    lam = (m * math.pi / a) ** 2 + (n * math.pi / b) ** 2
    c = 16.0 / (math.pi ** 2) / (m * n * lam)
    sgn = np.cos(0.5 * math.pi * (m - 1)) * np.cos(0.5 * math.pi * (n - 1))       # (-1)^((m-1)/2) (-1)^((n-1)/2)
    t_max = float(np.sum(c * sgn))
    t_mean = float(np.sum(c * (2.0 / (m * math.pi)) * (2.0 / (n * math.pi))))
    qx = float(np.sum(c * (m * math.pi / a) * (2.0 * b / (n * math.pi))))            # flux through the walls of length b (one wall)
    qy = float(np.sum(c * (n * math.pi / b) * (2.0 * a / (m * math.pi))))
    return t_mean, t_max, qx / (qx + qy)


# ---------------------------------------------------------------------------------------------------- the design's geometry
@dataclass(frozen=True)
class Geom:
    """Numbers taken from motorgeo.Params (all SI). Built once per design (cached by build_geom)."""
    slots: int
    pole_pairs: int
    r_so: float
    r_si: float
    r_ro: float
    r_shaft: float
    airgap: float
    r_shoe: float
    r_wall: float
    a_s: float                # slot half angle [rad]
    r_b: float                # radius of the slot's round bottom
    r_yb: float               # outer radius of the slot bottoms = inner radius of the yoke ring
    a_slot: float             # gross slot (winding) area, both layers
    a_layer: float            # one coil layer
    perim: float              # slot's wetted perimeter (two flanks and the bottom arc)
    flank: float              # length of one flank
    w_slot: float             # mean slot width (equivalent rectangle)
    h_slot: float             # equivalent depth = a_slot / w_slot
    a_steel_stator: float
    a_yoke: float
    a_teeth: float
    a_steel_rotor: float
    c_tooth: float            # tooth width = c_tooth r
    f_tooth_mean: float       # mean tooth temperature above its root = Q f / (k L), see build_geom
    f_tooth_tip: float
    f_root: float
    f_ring_half: float
    r_centroid: float         # radius of the slot's centroid
    span: float               # end-turn span (chord of the coil pitch at the centroid)
    r_magnet_zone: float      # radius of the magnet slots' centres
    magnets_per_row: int
    mag_w: float
    mag_t: float
    mag_len: float            # magnet length in a 20 mm row
    row_mm: float


@functools.lru_cache(maxsize=16)
def build_geom(p: mg.Params) -> Geom:
    """Everything the network needs from the cross-section: the closed forms of motorgeo._slot for the slot, shapely (motorgeo.visible) for the
    steel areas. The tooth, root and ring factors are the integrals of 1/(k A) along the heat path per unit of k L (so that R = f / (k L))."""
    a_s = math.radians(p.slot_half_deg)
    r_b = p.r_wall * math.tan(a_s)
    c_b = p.r_wall / math.cos(a_s)
    r_yb = c_b + r_b
    a_slot = 0.5 * (p.r_wall ** 2 - p.r_shoe ** 2) * math.sin(2 * a_s) + 0.5 * r_b ** 2 * (math.pi + 2 * a_s + math.sin(2 * a_s))
    flank = p.r_wall - p.r_shoe
    perim = 2 * flank + r_b * (math.pi + 2 * a_s)
    w_slot = (p.r_shoe + p.r_wall) * math.tan(a_s)
    regions = mg.build(p)
    vis = mg.visible(regions)
    a_steel_s = p.n_sectors * vis["stator_steel"].area * 1e-6
    a_steel_r = p.n_sectors * vis["rotor_steel"].area * 1e-6
    a_yoke = math.pi * (p.r_so ** 2 - r_yb ** 2) * 1e-6
    metrics, _ = mg.check(p, regions)
    a_layer = metrics["coil_layer_area_mm2"] * 1e-6
    from shapely.ops import unary_union
    centroid = unary_union([r.polygon() for r in regions if r.name in ("coil_0_0", "coil_0_1")]).centroid
    r_cen = math.hypot(centroid.x, centroid.y) * 1e-3
    span = 2 * r_cen * math.sin(math.pi * p.coil_pitch / p.slots)
    tau = 2 * math.pi / p.slots
    c_t = tau - 2 * math.tan(a_s)
    q = p.slots
    # tooth: heat injected uniformly between r_shoe and r_wall, width c_t r:  T(r) - T_root = 1/(k L) int_r^root Q(s)/(q c_t s) ds,  Q(s) = (s - ra)/len
    ra, rr = p.r_shoe * 1e-3, p.r_wall * 1e-3
    s = np.linspace(ra, rr, 4001)
    integrand = ((s - ra) / (rr - ra)) / (q * c_t * s)
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(s)
    t_rel = np.concatenate([np.cumsum(seg[::-1])[::-1], [0.0]])                       # T(s) - T_root on the grid
    f_tip = float(t_rel[0])
    f_mean = float(np.trapezoid(t_rel, s) / (rr - ra))
    # root: from r_wall to r_yb the steel's width between the slot's round bottoms is w(r) = r (tau - 2 phi(r))
    rc, rb = c_b * 1e-3, r_b * 1e-3
    rs = np.linspace(rr, r_yb * 1e-3, 2001)
    phi = np.arccos(np.clip((rs ** 2 + rc ** 2 - rb ** 2) / (2 * rs * rc), -1.0, 1.0))
    f_root = float(np.trapezoid(1.0 / (q * rs * (tau - 2 * phi)), rs))
    f_ring_half = math.log(p.r_so / r_yb) / (2 * math.pi) / 2
    return Geom(
        slots=p.slots, pole_pairs=p.pole_pairs, r_so=p.r_so * 1e-3, r_si=p.r_si * 1e-3, r_ro=p.r_ro * 1e-3, r_shaft=p.r_shaft * 1e-3,
        airgap=p.airgap * 1e-3, r_shoe=p.r_shoe * 1e-3, r_wall=p.r_wall * 1e-3, a_s=a_s, r_b=r_b * 1e-3, r_yb=r_yb * 1e-3,
        a_slot=a_slot * 1e-6, a_layer=a_layer, perim=perim * 1e-3, flank=flank * 1e-3, w_slot=w_slot * 1e-3, h_slot=a_slot / w_slot * 1e-3,
        a_steel_stator=a_steel_s, a_yoke=a_yoke, a_teeth=a_steel_s - a_yoke, a_steel_rotor=a_steel_r, c_tooth=c_t,
        f_tooth_mean=f_mean, f_tooth_tip=f_tip, f_root=f_root, f_ring_half=f_ring_half,
        r_centroid=r_cen, span=span, r_magnet_zone=math.hypot(p.mag_x, p.mag_y) * 1e-3, magnets_per_row=4 * p.pole_pairs,
        mag_w=p.mag_w * 1e-3, mag_t=p.mag_t * 1e-3, mag_len=19.14e-3, row_mm=20.0)


@functools.lru_cache(maxsize=1)
def _final_overrides_cached() -> str:
    with open(paths.path("out") / "optimization" / "final" / "final.json", encoding="utf-8") as f:
        return json.dumps(json.load(f)["overrides"], sort_keys=True)


def final_overrides() -> dict:
    """The optimized design: the "overrides" of out/optimization/final/final.json (motorgeo.Params fields on top of the band-free rotor)."""
    return json.loads(_final_overrides_cached())


def design_params(overrides: dict | None = None) -> mg.Params:
    """design.params(overrides): motorgeo.Params(barrier="none", **overrides). None means the optimized design."""
    from optimization import design
    return design.params(final_overrides() if overrides is None else overrides)


# ---------------------------------------------------------------------------------------------------- inputs and results
@dataclass(frozen=True)
class WindingSpec:
    """The winding as the network needs it: the copper fill of the GROSS slot (the winding study's k: copper area / the layer's area)."""
    fill: float = 0.45


@dataclass(frozen=True)
class Cooling:
    """How the housing is cooled. Build with Cooling.natural(), .forced_air(), .shaft_fan() or .liquid(); the geometry behind each is in Assumptions.
    forced_air blows a fixed air speed (Assumptions.air_speed_ms, an external fan) whatever the shaft does; shaft_fan drives the same fin channels with a fan
    on the shaft, so the air speed is v_ref_ms x speed / n_ref_rpm (and natural convection alone at standstill)."""
    kind: str                         # "natural" | "forced_air" | "liquid" | "shaft_fan"
    ambient_c: float = AMBIENT_C
    label: str = ""
    v_ref_ms: float = 6.0             # shaft_fan: the air speed in the fin channels at n_ref_rpm [m/s]
    n_ref_rpm: float = 3000.0         # shaft_fan: the speed v_ref_ms belongs to [rpm]

    @classmethod
    def natural(cls, ambient_c: float = AMBIENT_C):
        return cls("natural", ambient_c, "natural convection + radiation (plain housing, no fan)")

    @classmethod
    def forced_air(cls, ambient_c: float = AMBIENT_C):
        return cls("forced_air", ambient_c, "forced air over a finned housing")

    @classmethod
    def shaft_fan(cls, v_ref_ms: float = 6.0, n_ref_rpm: float = 3000.0, ambient_c: float = AMBIENT_C, label: str | None = None):
        """A fan on the shaft over the finned housing (like a simple TEFC motor): the air speed in the fin channels is v(n) = v_ref_ms x n / n_ref_rpm. At
        standstill the fan stops and only natural convection and radiation are left. The fan's own shaft power (ThermalModel.fan_power, shaft_fan_power) follows
        the fan laws, P proportional to n^3; it is a mechanical loss on the shaft and heats the air, not the motor. Defaults: 6 m/s at 3,000 rpm (sized for the rated
        speed). size_shaft_fan(budget_w) returns the v_ref_ms whose power at the top speed stays within a budget."""
        lab = label or f"shaft fan, {v_ref_ms:.2f} m/s in the fin channels at {n_ref_rpm:.0f} rpm"
        return cls("shaft_fan", ambient_c, lab, float(v_ref_ms), float(n_ref_rpm))

    @classmethod
    def liquid(cls, ambient_c: float = AMBIENT_C):
        return cls("liquid", ambient_c, "liquid jacket (water-glycol)")


def copper_resistivity(t_c: float) -> float:
    """rho(T) = rho20 (1 + alpha (T - 20)) [ohm m]."""
    return RHO_CU_20 * (1.0 + ALPHA_CU * (t_c - 20.0))


@dataclass(frozen=True)
class CopperLoss:
    """The DC copper loss of the whole winding at the current density j_rms [A/mm^2 rms in the copper], at the copper's temperature:
    P = rho(T) J^2 V_cu, split between the slot part (length L) and the end turns (the winding study's mean turn 2 L + pi span + 4 e).
    Build it with ThermalModel.copper(j_rms)."""
    j_rms: float
    v_slot: float                     # copper volume in the slots [m^3]
    v_end: float                      # copper volume in the end turns [m^3]

    def slot(self, t_c: float) -> float:
        return copper_resistivity(t_c) * (self.j_rms * 1e6) ** 2 * self.v_slot

    def end(self, t_c: float) -> float:
        return copper_resistivity(t_c) * (self.j_rms * 1e6) ** 2 * self.v_end

    def total(self, t_c: float) -> float:
        return self.slot(t_c) + self.end(t_c)

    def p20(self) -> tuple[float, float]:
        return self.slot(20.0), self.end(20.0)


@dataclass
class Losses:
    """The heat sources [W] (whole motor, at the stack length of the model). copper: a number (the total), a function f(T_C) -> total W, or a
    CopperLoss. bearing and windage are estimated by the model when None (windage splits between the rotor and the stator bore; the bearing
    loss goes to the bearings). copper_slot_fraction: the share of a number's or a function's copper loss in the slots (None: by length)."""
    copper: float | Callable[[float], float] | CopperLoss = 0.0
    iron_teeth: float = 0.0
    iron_yoke: float = 0.0
    iron_rotor: float = 0.0
    magnet: float = 0.0
    bearing: float | None = None
    windage: float | None = None
    copper_slot_fraction: float | None = None


@dataclass
class SteadyResult:
    """The answer of steady_state: temperatures [C], the winding's hot spot, the sources actually used, the margins and the binding limit."""
    T: dict
    hot_spot: float
    hot_spot_slot: float
    hot_spot_end: float
    magnet: float
    copper_loss: float
    copper_slot: float
    copper_end: float
    sources: dict
    heat_out: dict
    margin_winding: float
    margin_magnet: float
    binding: str
    iterations: int
    runaway: bool = False
    converged: bool = True


@dataclass
class TransientResult:
    t: np.ndarray                     # s
    T: np.ndarray                     # (n_t, 11) C
    hot_spot: np.ndarray
    magnet: np.ndarray
    copper_loss: np.ndarray
    j_rms: np.ndarray


# ---------------------------------------------------------------------------------------------------- iron-loss placeholder
# Bertotti three-term fit of 29-gauge M-19 (materials/current.py, F2; the maker's typical data, 50-400 Hz): energy per kg per cycle
# P/f = a + b f + c sqrt(f), at 1.0 T and 1.5 T. Interpolated by power laws in B (a ~ B^2.03, b ~ B^2, c ~ B^2.37).
_BERT = {1.0: (0.0151886, 5.17774e-5, 5.71845e-4), 1.5: (0.034612, 1.16499e-4, 1.49175e-3)}


def iron_w_per_kg(b_peak: float, f_hz: float) -> float:
    """Datasheet-based core loss of 29-gauge M-19 [W/kg] at b_peak [T] and f_hz [Hz] (placeholder-grade: sinusoidal B, no building factor)."""
    a1, b1, c1 = _BERT[1.0]
    a2, _, c2 = _BERT[1.5]
    ea, ec = math.log(a2 / a1) / math.log(1.5), math.log(c2 / c1) / math.log(1.5)
    return f_hz * (a1 * b_peak ** ea + b1 * b_peak ** 2 * f_hz + c1 * b_peak ** ec * math.sqrt(f_hz))


# ---------------------------------------------------------------------------------------------------- the model
class ThermalModel:
    """One built network: a design, a stack length, a cooling and a speed. Build once; call steady_state() as often as needed.

        ThermalModel(overrides=None, stack_mm=20, cooling=Cooling.natural(), speed_rpm=3000, winding=WindingSpec(), assume=Assumptions())
    """

    def __init__(self, overrides: dict | None = None, stack_mm: float = 20.0, cooling: Cooling | None = None, speed_rpm: float = RATED_RPM,
                 winding: WindingSpec | None = None, assume: Assumptions | None = None):
        self.p = design_params(overrides)
        self.g = build_geom(self.p)
        self.stack_mm = float(stack_mm)
        self.L = stack_mm * 1e-3
        self.rows = stack_mm / self.g.row_mm
        self.cooling = cooling or Cooling.natural()
        self.speed_rpm = float(speed_rpm)
        self.omega = 2 * math.pi * speed_rpm / 60.0
        self.winding = winding or WindingSpec()
        self.a = assume or Assumptions()
        self.amb = self.cooling.ambient_c
        self.t_in = self.a.coolant_in_c if self.a.coolant_in_c is not None else self.amb     # the jacket's coolant inlet
        self._slot_cache: dict = {}
        self._rect_cache: dict = {}
        self._build_constant()

    # ---------------------------------------------------------------- the constant part of the network
    def _build_constant(self) -> None:
        g, a, L = self.g, self.a, self.L
        al, fe, mg_, sh = MATERIALS["aluminium"], MATERIALS["m19"], MATERIALS["n45"], MATERIALS["shaft"]
        # --- the winding's lengths (the winding study's mean turn 2 L + pi span + 4 e)
        e = a.end_leg_mm * 1e-3
        self.l_end_side = a.end_turn_mm * 1e-3 if a.end_turn_mm is not None else math.pi * g.span / 2 + 2 * e     # end-turn copper per coil side, both ends
        self.l_overhang = a.end_overhang_mm * 1e-3 if a.end_overhang_mm is not None else g.span / 2 + e         # how far the end winding sticks out beyond the stack
        self.l_es = self.l_overhang + a.end_clearance_mm * 1e-3     # axial length of the end space, stack face to end cap
        self.f_slot = L / (L + self.l_end_side)                     # share of the copper loss that falls in the slots
        # --- the stator iron: tooth mean to yoke mean, and the yoke's two halves
        k_fe = fe["k_in"]
        self.g_ty = 1.0 / ((g.f_tooth_mean + g.f_root + g.f_ring_half) / (k_fe * L))
        self.r_tip_minus_mean = (g.f_tooth_tip - g.f_tooth_mean) / (k_fe * L)             # K/W: the tooth's tip above its mean per watt it passes
        self.r_yo = g.f_ring_half / (k_fe * L)
        self.r_contact = a.contact_factor * a.contact_nrel_mm2kw * 1e-6 / (2 * math.pi * g.r_so * L)       # stator-housing contact [K/W]
        # --- the housing
        kind = self.cooling.kind
        t_w = a.housing_wall_mm * 1e-3
        if kind == "liquid":
            self.r_ho = g.r_so + (a.jacket_wall_mm + a.jacket_channel_h_mm + a.jacket_outer_wall_mm) * 1e-3
            self.r_hmid = g.r_so + a.jacket_wall_mm * 1e-3 / 2        # the H node sits in the middle of the wall under the channel
        else:
            self.r_ho = g.r_so + t_w
            self.r_hmid = g.r_so + t_w / 2
        self.r_h1 = math.log(self.r_hmid / g.r_so) / (2 * math.pi * al["k"] * L)           # stator OD to the housing's mid-wall
        self.g_yh = 1.0 / (self.r_yo + self.r_contact + self.r_h1)                          # yoke mean to housing node
        a_wall_ax = math.pi * (self.r_ho ** 2 - g.r_so ** 2)
        self.r_bh = (self.l_es / 2) / (al["k"] * a_wall_ax) / 2 + a.flange_mm2kw * 1e-6 / a_wall_ax / 2   # bells to housing, two joints in parallel
        # --- the end space
        self.A_ew_air = 2 * a.end_area_factor * (math.pi * (g.r_wall ** 2 - g.r_shoe ** 2) + 2 * math.pi * g.r_shoe * self.l_overhang)
        self.A_ew_out = 2 * a.end_area_factor * 2 * math.pi * (g.r_wall + 0.5 * (g.r_yb - g.r_wall)) * self.l_overhang
        self.A_endwall = 2 * (math.pi * (g.r_so ** 2 - g.r_shaft ** 2) + 2 * math.pi * g.r_so * self.l_es)
        self.A_rot_end = 2 * math.pi * (g.r_ro ** 2 - g.r_shaft ** 2)
        self.A_shaft_end = 2 * math.pi * g.r_shaft * self.l_es
        self.h_end = end_space_h(self.omega * g.r_ro, a)
        self.g_ea = self.h_end * self.A_ew_air
        self.g_ab = self.h_end * self.A_endwall
        # through-stack conduction (the laminations conduct along the shaft 13 times worse than along the sheet: NREL's 1.68 against 21.9 W/(m K)):
        # the mean of a stack with a uniform source to its two end faces is L / (12 k A) away, and then the end space's convection
        k_ax = fe["k_ax"]
        r_ax_stator = L / (12.0 * k_ax * g.a_steel_stator)
        r_ax_rotor = L / (12.0 * k_ax * g.a_steel_rotor)
        self.g_ta = 1.0 / (r_ax_stator + 1.0 / (self.h_end * 2.0 * g.a_steel_stator))          # stator iron to the end-space air, through the end faces
        self.g_ar = 1.0 / (1.0 / (self.h_end * self.A_rot_end) + r_ax_rotor)                    # rotor steel to the end-space air, through its end faces
        self.g_as = self.h_end * self.A_shaft_end
        gap_to_wall = max(g.r_so - g.r_yb, 0.5e-3)
        self.g_eb_direct = (air_props(a.t_gap_film_c)["k"] / gap_to_wall + h_radiation(0.5, 100.0, 60.0)) * self.A_ew_out * 0.5
        # --- the axial copper path (the wires run through slot and end turn): in metres, divided by (k_cu A_cu n) in _slot_branch
        self.r_ax_geom = (self.l_end_side / 2) / 12.0 + (L / 2) / 3.0
        # --- the air gap, the rotor chain (surface -> magnet zone -> shaft), the magnets, the shaft and the bearings
        self.gap = air_gap_h(self.omega, g.r_ro, g.airgap, a.t_gap_film_c)
        r_lm = g.airgap / math.log(g.r_si / g.r_ro)
        self.a_gap = 2 * math.pi * r_lm * L
        k_rot = a.rotor_k_fraction * k_fe
        r_surface_to_zone = math.log(g.r_ro / g.r_magnet_zone) / (2 * math.pi * k_rot * L)
        self.r_tr = 1.0 / (self.gap["h"] * self.a_gap) + r_surface_to_zone                 # tooth to the rotor node (gap + the 4 mm of steel)
        self.r_rs = (math.log(g.r_magnet_zone / g.r_shaft) / (2 * math.pi * k_rot * L)
                     + a.rotor_shaft_mm2kw * 1e-6 / (2 * math.pi * g.r_shaft * L))
        self.n_mag = self.rows * g.magnets_per_row
        self.v_mag = self.n_mag * g.mag_w * g.mag_t * g.mag_len
        mag_face = self.n_mag * 2 * (g.mag_w + g.mag_t) * g.mag_len                         # the four long faces
        mag_large = self.n_mag * 2 * g.mag_w * g.mag_len
        self.r_mr = (a.magnet_gap_mm * 1e-3 / a.k_adhesive) / mag_face + g.mag_t / (12.0 * mg_["k"] * mag_large)   # adhesive + inside the magnet
        a_sh = math.pi * g.r_shaft ** 2
        stub = self.l_es + a.shaft_stub_mm * 1e-3
        self.r_sg = 0.5 * (stub / (sh["k"] * a_sh) + (L / 2) / (3.0 * sh["k"] * a_sh))      # two ends in parallel
        self.r_gb = a.bearing_r_kw / 2.0                                                    # two bearings in parallel
        self._setup_external()
        self._setup_capacities()
        self.p_bearing = bearing_loss(self.omega, a)
        self.p_windage = windage_loss(self.omega, g.r_ro, g.airgap, L, a.t_gap_film_c, a.wind_gap_factor)
        G = np.zeros((N_NODES, N_NODES))
        for i, j, gc in (("T", "Y", self.g_ty), ("Y", "H", self.g_yh), ("E", "A", self.g_ea), ("E", "B", self.g_eb_direct), ("A", "B", self.g_ab),
                         ("A", "R", self.g_ar), ("A", "S", self.g_as), ("T", "A", self.g_ta), ("T", "R", 1.0 / self.r_tr), ("R", "M", 1.0 / self.r_mr),
                         ("R", "S", 1.0 / self.r_rs), ("S", "G", 1.0 / self.r_sg), ("G", "B", 1.0 / self.r_gb), ("B", "H", 1.0 / self.r_bh)):
            self._link_into(G, i, j, gc)
        self.G0 = G

    # ---------------------------------------------------------------- the outside
    def _setup_external(self) -> None:
        a, L = self.a, self.L
        al = MATERIALS["aluminium"]
        self.d_ho = 2 * self.r_ho
        self.A_cylH = 2 * math.pi * self.r_ho * L
        self.A_cylB = 2 * math.pi * self.r_ho * 2 * self.l_es
        self.A_disc = 2 * math.pi * self.r_ho ** 2
        self.mass_fins_per_len = 0.0
        self.v_air = 0.0
        if self.cooling.kind in ("forced_air", "shaft_fan"):
            af = air_props(a.t_air_film_c)
            hf, tf, n = a.fin_height_mm * 1e-3, a.fin_thickness_mm * 1e-3, a.n_fins
            gap = 2 * math.pi * (self.r_ho + hf / 2) / n - tf
            dh = 2 * gap * hf / (gap + hf)
            # the air speed in the channels: fixed for forced_air, proportional to the shaft speed for the shaft fan
            self.v_air = a.air_speed_ms if self.cooling.kind == "forced_air" else self.cooling.v_ref_ms * self.speed_rpm / self.cooling.n_ref_rpm
            self.fin_gap, self.fin_dh, self.fin_a_flow = gap, dh, n * gap * hf         # channel width at mid-height, hydraulic diameter, free flow area
            self.l_channel = L + 2 * self.l_es                                           # the air runs along the housing and both end-bell cylinders
            self.fin_re = self.v_air * dh / af["nu"]
            self.h_fin = nusselt_duct(self.fin_re, af["Pr"], laminar_nu=4.0) * af["k"] / dh
            self.fin_eta = fin_efficiency(self.h_fin, al["k"], tf, hf)
            per_len_base = 2 * math.pi * self.r_ho - n * tf
            per_len_fin = n * (2 * hf + tf)
            self.g_fin_per_len = self.h_fin * (per_len_base + self.fin_eta * per_len_fin)    # W/(m K) per metre of finned length
            self.mass_fins_per_len = n * hf * tf * al["rho"]
            if self.cooling.kind == "shaft_fan":
                # the finned surfaces (housing and bell cylinders) would pass g_raw to air at the inlet temperature if the air did not warm; the stream can carry only
                # m cp (T_wall - T_in) away, so the effective conductance is the effectiveness times m cp, eps = 1 - exp(-g_raw / (m cp)) (an isothermal wall).
                # Shared between the housing and the bells by length. No flow, no forced convection.
                g_raw = self.g_fin_per_len * self.l_channel
                self.mdot_cp_air = af["rho"] * self.v_air * self.fin_a_flow * af["cp"]
                if self.v_air <= 0.0:
                    g_eff = 0.0
                elif a.air_bulk_warming:
                    g_eff = self.mdot_cp_air * (1.0 - math.exp(-g_raw / self.mdot_cp_air))
                else:
                    g_eff = g_raw
                self.g_forced_h = g_eff * L / self.l_channel
                self.g_forced_b = g_eff * 2 * self.l_es / self.l_channel
                self.dT_fan = self.fan_state()["dT_air_k"] if a.fan_air_preheat else 0.0  # the fan's loss warms the stream (external_full)
        if self.cooling.kind == "liquid":
            self._setup_jacket()

    def _setup_jacket(self) -> None:
        a, g, L = self.a, self.g, self.L
        al = MATERIALS["aluminium"]
        w, h = a.jacket_channel_w_mm * 1e-3, a.jacket_channel_h_mm * 1e-3
        dh = 2 * w * h / (w + h)
        vel = (a.flow_lpm / 60000.0) / (w * h)
        pr = a.coolant_mu * a.coolant_cp / a.coolant_k
        self.jacket_re = a.coolant_rho * vel * dh / a.coolant_mu
        self.jacket_nu = nusselt_duct(self.jacket_re, pr, laminar_nu=3.7)
        self.h_jacket = self.jacket_nu * a.coolant_k / dh
        self.jacket_phi = (w + 2 * h * a.jacket_rib_efficiency) / (w + a.jacket_rib_mm * 1e-3)    # wetted area / the stator OD's area
        a_in = 2 * math.pi * g.r_so * L
        self.r_jacket = a.jacket_wall_mm * 1e-3 / 2 / (al["k"] * a_in) + 1.0 / (self.h_jacket * self.jacket_phi * a_in)   # H node to the coolant
        self.mdot_cp = a.coolant_rho * a.flow_lpm / 60000.0 * a.coolant_cp

    def external(self, t_h: float, t_b: float, q_cool: float) -> tuple[float, float, float, float]:
        """(G_H_amb, G_B_amb, G_H_cool, T_coolant) for the housing at t_h and the end bells at t_b [C]; q_cool [W] is the heat the jacket takes."""
        a, kind, amb, L = self.a, self.cooling.kind, self.amb, self.L
        eps = a.emissivity
        th, tb = max(t_h, amb + 0.01), max(t_b, amb + 0.01)
        r_h2 = math.log(self.r_ho / self.r_hmid) / (2 * math.pi * MATERIALS["aluminium"]["k"] * L)        # mid-wall to the outer surface
        g_hc, tc = 0.0, amb
        if kind == "liquid":
            tc = self.t_in + q_cool / (2.0 * self.mdot_cp)
            g_hc = 1.0 / self.r_jacket
        hb_d = h_natural_vplate(self.d_ho, tb, amb) + h_radiation(eps, tb, amb)
        if kind == "natural":
            fin = a.natural_fin_gain
            g_h_ext = (h_natural_cylinder(self.d_ho, th, amb) + h_radiation(eps, th, amb)) * self.A_cylH * fin
            g_b_ext = (h_natural_cylinder(self.d_ho, tb, amb) + h_radiation(eps, tb, amb)) * self.A_cylB * fin + hb_d * self.A_disc
        elif kind == "forced_air":
            g_h_ext = self.g_fin_per_len * L
            g_b_ext = self.g_fin_per_len * 2 * self.l_es + hb_d * self.A_disc
        elif kind == "shaft_fan":
            p = self._fan_parts(th, tb)
            g_h_ext = p["h_conv"] + p["h_rad"]
            g_b_ext = p["b_conv"] + p["b_rad"]
        else:                                                                        # liquid: the jacket takes the stack part; the bells stand in still air
            g_h_ext = 0.0
            g_b_ext = (h_natural_cylinder(self.d_ho, tb, amb) + h_radiation(eps, tb, amb)) * self.A_cylB + hb_d * self.A_disc
        g_ha = 1.0 / (r_h2 + 1.0 / g_h_ext) if g_h_ext > 0 else 0.0
        return g_ha, g_b_ext, g_hc, tc

    def external_full(self, t_h: float, t_b: float, q_cool: float) -> tuple[float, float, float, float, float, float]:
        """external() plus the two constant heat flows [W] into the housing and the end bells that the shaft fan's air preheat adds: (G_H_amb, G_B_amb, G_H_cool,
        T_coolant, x_h, x_b). The fan's loss warms the air before it reaches the fins by dT = P_shaft / (m cp) (fan_state), so the convective part of the
        outside gives its heat to air at ambient + dT, not at the ambient: Q = G (T - T_amb) - G_conv dT, and the second term is the source x here (the heat the
        fan's loss keeps the surfaces from giving away). Zero for every other cooling and at standstill; external() is the same without it."""
        g_ha, g_ba, g_hc, tc = self.external(t_h, t_b, q_cool)
        x_h = x_b = 0.0
        if self.cooling.kind == "shaft_fan" and self.v_air > 0.0 and self.dT_fan > 0.0:
            amb = self.amb
            p = self._fan_parts(max(t_h, amb + 0.01), max(t_b, amb + 0.01))
            g_tot = p["h_conv"] + p["h_rad"]
            x_h = g_ha * self.dT_fan * p["h_conv"] / g_tot if g_tot > 0 else 0.0
            x_b = p["b_conv_fins"] * self.dT_fan
        return g_ha, g_ba, g_hc, tc, x_h, x_b

    def _fan_parts(self, th: float, tb: float) -> dict:
        """The shaft fan's outside at the housing surface th and the bell surface tb [C], as conductances to the ambient [W/K], convection and radiation apart:
        h_conv, h_rad (the housing's cylinder), b_conv, b_rad (the bells: the finned cylinders and the end discs). The finned surfaces' convection is the forced
        conductance (set up by _setup_external, with the air's warming) blended with the still-air natural convection of the finned housing,
        G = (G_forced^n + G_natural^n)^(1/n); the discs stand in still air as in the other coolings."""
        a, amb = self.a, self.amb
        n, eps = a.mixed_exponent, a.emissivity
        nat_h = h_natural_cylinder(self.d_ho, th, amb) * self.A_cylH * a.fan_fin_gain_conv
        nat_b = h_natural_cylinder(self.d_ho, tb, amb) * self.A_cylB * a.fan_fin_gain_conv
        conv_h = (self.g_forced_h ** n + nat_h ** n) ** (1.0 / n)
        conv_b = (self.g_forced_b ** n + nat_b ** n) ** (1.0 / n)
        return {"h_conv": conv_h, "h_rad": h_radiation(eps, th, amb) * self.A_cylH * a.fan_fin_gain_rad,
                "b_conv": conv_b + h_natural_vplate(self.d_ho, tb, amb) * self.A_disc,
                "b_rad": h_radiation(eps, tb, amb) * (self.A_cylB * a.fan_fin_gain_rad + self.A_disc),
                "b_conv_fins": conv_b, "b_conv_discs": h_natural_vplate(self.d_ho, tb, amb) * self.A_disc,
                "h_forced": self.g_forced_h, "h_natural": nat_h}

    def outside_parts(self, t_h: float, t_b: float) -> dict:
        """The outside's conductances to the ambient [W/K] at the housing surface t_h and the end bells t_b [C], convection and radiation apart (the sum is what
        external() uses): {"h_conv", "h_rad"} for the housing's cylinder (zero for the liquid jacket, which has "h_cool" instead: the jacket's conductance to the
        coolant) and {"b_conv", "b_rad"} for the end bells (cylinders and discs). The forced_air preset has no radiation term on its fins (an old simplification:
        a few per cent of its conductance), so its "h_rad" is zero."""
        a, kind, amb, L = self.a, self.cooling.kind, self.amb, self.L
        eps = a.emissivity
        th, tb = max(t_h, amb + 0.01), max(t_b, amb + 0.01)
        hb_conv = h_natural_vplate(self.d_ho, tb, amb) * self.A_disc
        hb_rad = h_radiation(eps, tb, amb) * self.A_disc
        out = {"h_conv": 0.0, "h_rad": 0.0, "b_conv": 0.0, "b_rad": 0.0, "h_cool": 0.0}
        if kind == "natural":
            fin = a.natural_fin_gain
            out["h_conv"] = h_natural_cylinder(self.d_ho, th, amb) * self.A_cylH * fin
            out["h_rad"] = h_radiation(eps, th, amb) * self.A_cylH * fin
            out["b_conv"] = h_natural_cylinder(self.d_ho, tb, amb) * self.A_cylB * fin + hb_conv
            out["b_rad"] = h_radiation(eps, tb, amb) * self.A_cylB * fin + hb_rad
        elif kind == "forced_air":
            out["h_conv"] = self.g_fin_per_len * L
            out["b_conv"] = self.g_fin_per_len * 2 * self.l_es + hb_conv
            out["b_rad"] = hb_rad
        elif kind == "shaft_fan":
            p = self._fan_parts(th, tb)
            out.update({k: p[k] for k in ("h_conv", "h_rad", "b_conv", "b_rad")})
        else:
            out["h_cool"] = 1.0 / self.r_jacket
            out["b_conv"] = h_natural_cylinder(self.d_ho, tb, amb) * self.A_cylB + hb_conv
            out["b_rad"] = h_radiation(eps, tb, amb) * self.A_cylB + hb_rad
        return out

    # ---------------------------------------------------------------- the fan
    def fan_state(self, v: float | None = None) -> dict:
        """The air path of the fin channels at the air speed v [m/s] (default: the cooling's own at the model's speed). The pressure the fan must build is the
        channel's friction plus the rest of the path, dp = (f L_ch / d_h + K) rho v^2 / 2 (f the Darcy friction factor of the channel, d_h its hydraulic diameter,
        L_ch the finned length, K = Assumptions.fan_k_system); the air power is dp times the flow Q = v A_flow; the shaft power is that over the fan's efficiency.
        With K and f nearly constant, dp goes as v^2, Q as v and so the power as v^3: the fan laws. Returns {"v", "re", "f_darcy", "dp_pa", "q_m3s", "p_air_w",
        "p_shaft_w", "mdot_cp_w_per_k", "dT_air_k", "a_flow_m2", "d_h_m", "l_channel_m"}; dT_air_k is the warming of the stream if all the shaft power ended as heat in the
        air (the fan's loss heats the air, not the motor; it is not in the network's sources)."""
        if self.cooling.kind not in ("forced_air", "shaft_fan"):
            raise ValueError(f"fan_state: a {self.cooling.kind} cooling has no fin channels")
        a = self.a
        af = air_props(a.t_air_film_c)
        v = self.v_air if v is None else max(float(v), 0.0)
        re = v * self.fin_dh / af["nu"]
        f = duct_friction(re, a.duct_f_re)
        dp = (f * self.l_channel / self.fin_dh + a.fan_k_system) * 0.5 * af["rho"] * v * v
        q = v * self.fin_a_flow
        p_air = dp * q
        mdot_cp = af["rho"] * q * af["cp"]
        p_shaft = p_air / a.fan_efficiency
        return {"v": v, "re": re, "f_darcy": f, "dp_pa": dp, "q_m3s": q, "p_air_w": p_air, "p_shaft_w": p_shaft, "mdot_cp_w_per_k": mdot_cp,
                "dT_air_k": (p_shaft / mdot_cp) if mdot_cp > 0 else 0.0, "a_flow_m2": self.fin_a_flow, "d_h_m": self.fin_dh, "l_channel_m": self.l_channel}

    def fan_power(self) -> float:
        """The power [W] the air path needs at the model's speed: the shaft power of a shaft fan, or (for the forced_air preset, whose fan is external) the power
        that fixed 6 m/s would take. Zero for the natural and liquid coolings (the jacket's pump is not modelled)."""
        if self.cooling.kind in ("forced_air", "shaft_fan"):
            return self.fan_state()["p_shaft_w"]
        return 0.0

    # ---------------------------------------------------------------- where the heat goes
    def heat_flows(self, res: SteadyResult, fill: float | None = None) -> dict:
        """The heat flow [W] on every branch of the network at the steady state `res` (from steady_state), and the groups that were asked about.
        Every link is q = G (T_from - T_to), written "from>to"; a negative value is the heat going the other way. Returns
          made    the sources: copper_slot, copper_end, iron_teeth, iron_yoke, iron_rotor, magnet, bearing, windage; total = their sum (the fan's power is not
                  here: it is a shaft load that heats the air, not a source in the motor)
          links   W>T, W>Y (the slot winding to the teeth and to the yoke), E>W (the end turns' heat through the slot), E>A (end turns to the end-space air),
                  E>B, T>Y, Y>H, T>A, T>R (stator teeth to the rotor across the air gap), R>M, R>S, A>R, A>S, A>B, S>G, G>B, B>H
          exits   housing_convection, housing_radiation, bells_convection, bells_radiation (the surfaces that give the heat to the air) and housing_to_coolant
                  (the liquid jacket); their sum is `total`, which is checked against made (balance_w, ~1e-9)
          groups  slot_to_teeth, slot_to_yoke (W>T, W>Y; they add up to the slot copper plus the end turns' heat that came through the slot);
                  end_turns_via_slot (E>W), end_turns_to_air (E>A), end_turns_to_bells (E>B) (they add up to the end-turn copper loss);
                  rotor_made (the rotor's own heat: iron, half the windage, magnets), rotor_across_gap (positive: from the rotor to the stator teeth),
                  rotor_to_end_space_air, rotor_through_shaft (shaft to bearings) (they add up to rotor_made); bearings_to_bells (the shaft's heat plus the bearing
                  loss, entering the end bells through the bearings)
        The shaft's far end (the stub beyond the bearings) and the coupling carry no heat in this network: the shaft and the bearings give their heat to the
        end bells, which give it to the air."""
        fill = self.winding.fill if fill is None else fill
        s = self._slot_branch(fill)
        T = res.T
        links = {
            "W>T": s["g_wt"] * (T["W"] - T["T"]), "W>Y": s["g_wy"] * (T["W"] - T["Y"]), "E>W": s["g_we"] * (T["E"] - T["W"]),
            "E>A": self.g_ea * (T["E"] - T["A"]), "E>B": self.g_eb_direct * (T["E"] - T["B"]), "T>Y": self.g_ty * (T["T"] - T["Y"]),
            "Y>H": self.g_yh * (T["Y"] - T["H"]), "T>A": self.g_ta * (T["T"] - T["A"]), "T>R": (T["T"] - T["R"]) / self.r_tr,
            "R>M": (T["R"] - T["M"]) / self.r_mr, "R>S": (T["R"] - T["S"]) / self.r_rs, "A>R": self.g_ar * (T["A"] - T["R"]),
            "A>S": self.g_as * (T["A"] - T["S"]), "A>B": self.g_ab * (T["A"] - T["B"]), "S>G": (T["S"] - T["G"]) / self.r_sg,
            "G>B": (T["G"] - T["B"]) / self.r_gb, "B>H": (T["B"] - T["H"]) / self.r_bh,
        }
        parts = self.outside_parts(T["H"], T["B"])
        q_h, q_b, q_c = res.heat_out["housing_to_ambient"], res.heat_out["bells_to_ambient"], res.heat_out["housing_to_coolant"]
        h_tot, b_tot = parts["h_conv"] + parts["h_rad"], parts["b_conv"] + parts["b_rad"]
        g_ha, g_ba, _, _ = self.external(T["H"], T["B"], q_c)
        # radiation gives its heat to the ambient surfaces (the share of the conductance times the rise); the convection is what is left of the surface's heat
        rad_h = g_ha * parts["h_rad"] / h_tot * (T["H"] - self.amb) if h_tot > 0 else 0.0
        rad_b = g_ba * parts["b_rad"] / b_tot * (T["B"] - self.amb) if b_tot > 0 else 0.0
        exits = {"housing_convection": q_h - rad_h, "housing_radiation": rad_h, "bells_convection": q_b - rad_b, "bells_radiation": rad_b, "housing_to_coolant": q_c}
        made = dict(res.sources)
        total = float(sum(made.values()))
        p_rotor = made["iron_rotor"] + 0.5 * made["windage"] + made["magnet"]
        groups = {"slot_to_teeth": links["W>T"], "slot_to_yoke": links["W>Y"], "end_turns_via_slot": links["E>W"], "end_turns_to_air": links["E>A"],
                  "end_turns_to_bells": links["E>B"], "rotor_made": p_rotor, "rotor_across_gap": -links["T>R"],
                  "rotor_to_end_space_air": -links["A>R"] - links["A>S"], "rotor_through_shaft": links["S>G"], "bearings_to_bells": links["G>B"]}
        return {"made": made, "total_made": total, "links": links, "exits": exits, "total_out": float(sum(exits.values())), "groups": groups,
                "balance_w": total - float(sum(exits.values())), "outside": parts}

    # ---------------------------------------------------------------- heat capacities
    def _setup_capacities(self) -> None:
        g, a, L = self.g, self.a, self.L
        al, fe, mg_, sh = MATERIALS["aluminium"], MATERIALS["m19"], MATERIALS["n45"], MATERIALS["shaft"]
        c = np.zeros(N_NODES)
        c[IDX["T"]] = g.a_teeth * L * fe["rho"] * fe["c"]
        c[IDX["Y"]] = g.a_yoke * L * fe["rho"] * fe["c"]
        wall_vol = math.pi * (self.r_ho ** 2 - g.r_so ** 2) * L
        if self.cooling.kind == "liquid":                              # the housing around the channels is aluminium except the channels' volume
            w = a.jacket_channel_w_mm * 1e-3
            r1 = g.r_so + a.jacket_wall_mm * 1e-3
            r2 = r1 + a.jacket_channel_h_mm * 1e-3
            wall_vol -= w / (w + a.jacket_rib_mm * 1e-3) * math.pi * (r2 ** 2 - r1 ** 2) * L
        fins_h = self.mass_fins_per_len * L
        c[IDX["H"]] = (wall_vol * al["rho"] + fins_h) * al["c"]
        bell_vol = 2 * math.pi * self.r_ho ** 2 * a.bell_mm * 1e-3 + math.pi * (self.r_ho ** 2 - g.r_so ** 2) * 2 * self.l_es
        c[IDX["B"]] = (bell_vol * al["rho"] + self.mass_fins_per_len * 2 * self.l_es) * al["c"]
        c[IDX["A"]] = 2 * math.pi * g.r_so ** 2 * self.l_es * 1.1 * 1007.0
        c[IDX["R"]] = g.a_steel_rotor * L * fe["rho"] * fe["c"]
        c[IDX["M"]] = self.v_mag * mg_["rho"] * mg_["c"]
        c[IDX["S"]] = math.pi * g.r_shaft ** 2 * (L + 2 * (self.l_es + a.shaft_stub_mm * 1e-3)) * sh["rho"] * sh["c"]
        c[IDX["G"]] = 2 * a.bearing_mass_kg * MATERIALS["bearing_steel"]["c"]
        self.C0 = c

    # ---------------------------------------------------------------- the slot, as a function of the fill
    def _slot_branch(self, fill: float) -> dict:
        """The fill-dependent part of the network: the winding block's equivalent conductivity, its resistances to the tooth and the yoke, the
        hot-spot resistance, the axial copper path and the heat capacities of W and E. Cached per fill."""
        key = round(fill, 9)
        if key in self._slot_cache:
            return self._slot_cache[key]
        if len(self._slot_cache) > 512:
            self._slot_cache.clear()
        g, a, L = self.g, self.a, self.L
        q = g.slots
        cu, ln, vr = MATERIALS["copper"], MATERIALS["liner"], MATERIALS["varnish"]
        t_l = a.liner_mm * 1e-3
        a_liner = t_l * g.perim
        a_block = g.a_slot - a_liner
        a_cu = fill * g.a_slot                                            # copper of one slot (both layers)
        f_loc = a_cu / a_block
        if f_loc >= 0.9:
            raise ValueError(f"fill {fill} leaves no room: copper is {f_loc:.2f} of the block inside the liner")
        k_eq = a.k_matrix * (1 + f_loc) / (1 - f_loc)                      # Hashin-Shtrikman / Maxwell lower bound, copper >> matrix
        # the equivalent rectangle of the block: depth h_b (no liner on the top), width w_b; the top wall is adiabatic, so mirror it
        h_b = g.h_slot - t_l
        w_b = a_block / h_b
        rk = (round(w_b, 9), round(h_b, 9))
        if rk not in self._rect_cache:
            self._rect_cache[rk] = rect_poisson(w_b, 2 * h_b)
        t_mean, t_max, f_x = self._rect_cache[rk]
        r_mean = t_mean / (q * k_eq * a_block * L)                          # per slot R = theta / (k A L), q slots in parallel
        r_hs = (t_max - t_mean) / (q * k_eq * a_block * L)
        share_fl = 2 * h_b / (2 * h_b + w_b)                                # the liner's area splits like the rectangle's walls
        s_tot = q * g.perim * L
        r_wt = r_mean / f_x + t_l / (ln["k"] * share_fl * s_tot)            # mean winding to the tooth: the flank path (block + liner)
        r_wy = r_mean / (1 - f_x) + t_l / (ln["k"] * (1 - share_fl) * s_tot)     # ... to the yoke: the bottom path
        a_cu_side = a_cu / 2
        r_ew = self.r_ax_geom / (cu["k"] * a_cu_side * 2 * q)               # K/W: end-turn mean to slot mean along the copper
        v_cu_slot = a_cu * q * L
        v_cu_end = a_cu * q * self.l_end_side
        ins_ratio = (1 - f_loc) / f_loc
        c_w = v_cu_slot * cu["rho"] * cu["c"] + v_cu_slot * ins_ratio * vr["rho"] * vr["c"] + a_liner * q * L * ln["rho"] * ln["c"]
        c_e = v_cu_end * cu["rho"] * cu["c"] + v_cu_end * ins_ratio * vr["rho"] * vr["c"]
        out = {"k_eq": k_eq, "f_loc": f_loc, "a_block": a_block, "g_wt": 1.0 / r_wt, "g_wy": 1.0 / r_wy, "r_hs": r_hs, "r_mean": r_mean,
               "f_x": f_x, "g_we": 1.0 / r_ew, "c_w": c_w, "c_e": c_e, "v_slot": v_cu_slot, "v_end": v_cu_end, "a_cu_side": a_cu_side,
               "r_wt": r_wt, "r_wy": r_wy, "w_b": w_b, "h_b": h_b, "t_mean": t_mean, "t_max": t_max}
        self._slot_cache[key] = out
        return out

    # ---------------------------------------------------------------- losses
    def copper(self, j_rms: float, fill: float | None = None) -> CopperLoss:
        """The DC copper loss at the current density j_rms [A/mm^2 rms in the copper] for the model's winding (or for `fill`)."""
        s = self._slot_branch(self.winding.fill if fill is None else fill)
        return CopperLoss(j_rms, s["v_slot"], s["v_end"])

    def ampere_turns(self, j_rms: float, fill: float | None = None) -> dict:
        """The ampere-turns per coil side at j_rms: rms and peak (the project's 100 A-turns at 10 A and 10 turns are peak values)."""
        s = self._slot_branch(self.winding.fill if fill is None else fill)
        at_rms = j_rms * 1e6 * s["a_cu_side"]
        return {"rms": at_rms, "peak": math.sqrt(2) * at_rms, "a_cu_side_mm2": s["a_cu_side"] * 1e6}

    def iron_placeholder(self, rpm: float | None = None) -> dict:
        """The PROVISIONAL iron loss [W] at `rpm` (default: the model's): datasheet W/kg at the tooth and yoke flux densities (Assumptions), the
        electrical frequency f = p n / 60, M-19 at 7,650 kg/m^3; the rotor's as a share of the stator's. Replace it with the loss map."""
        rpm = self.speed_rpm if rpm is None else rpm
        f = self.g.pole_pairs * rpm / 60.0
        rho = MATERIALS["m19"]["rho"]
        p_t = iron_w_per_kg(self.a.b_tooth_t, f) * self.g.a_teeth * self.L * rho
        p_y = iron_w_per_kg(self.a.b_yoke_t, f) * self.g.a_yoke * self.L * rho
        return {"iron_teeth": p_t, "iron_yoke": p_y, "iron_rotor": self.a.rotor_iron_fraction * (p_t + p_y), "f_hz": f}

    def default_losses(self, copper=0.0, rpm: float | None = None, magnet: float = 0.0) -> Losses:
        """Losses with the iron placeholder, the model's bearing and windage estimates, and the copper given (a number, a function or a CopperLoss)."""
        ir = self.iron_placeholder(rpm)
        return Losses(copper=copper, iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"], magnet=magnet)

    @staticmethod
    def _link_into(G: np.ndarray, i: str, j: str, gc: float) -> None:
        ii, jj = IDX[i], IDX[j]
        G[ii, ii] += gc
        G[jj, jj] += gc
        G[ii, jj] -= gc
        G[jj, ii] -= gc

    # ---------------------------------------------------------------- the steady state
    def steady_state(self, losses: Losses, fill: float | None = None, wind_limit_c: float = CLASS_F_C, magnet_limit_c: float = MAGNET_LIMIT_C,
                     t_guess: np.ndarray | None = None, max_iter: int = 60, tol: float = 1e-4) -> SteadyResult:
        """Solve the network at the model's cooling, stack length and speed. `fill` overrides the winding's copper fill (the winding study varies it).
        The outside's temperature dependence (and a copper loss given as a function) is iterated: Picard iteration, which converges in a few
        steps because h goes as dT^0.25-0.33; a CopperLoss is linear in T and enters the matrix instead."""
        fill = self.winding.fill if fill is None else fill
        s = self._slot_branch(fill)
        amb = self.amb
        cu = losses.copper
        lin, func, fs = None, None, self.f_slot if losses.copper_slot_fraction is None else losses.copper_slot_fraction
        fixed_s = fixed_e = 0.0
        if isinstance(cu, CopperLoss):
            lin = CopperLoss(cu.j_rms, s["v_slot"], s["v_end"]).p20()              # the volumes of THIS fill
        elif callable(cu):
            func = cu
        else:
            fixed_s, fixed_e = float(cu) * fs, float(cu) * (1 - fs)
        p_bear = self.p_bearing if losses.bearing is None else losses.bearing
        p_wind = self.p_windage if losses.windage is None else losses.windage
        P = np.zeros(N_NODES)
        P[IDX["T"]] = losses.iron_teeth + 0.5 * p_wind
        P[IDX["Y"]] = losses.iron_yoke
        P[IDX["R"]] = losses.iron_rotor + 0.5 * p_wind
        P[IDX["M"]] = losses.magnet
        P[IDX["G"]] = p_bear
        base = self.G0.copy()
        self._link_into(base, "W", "T", s["g_wt"])
        self._link_into(base, "W", "Y", s["g_wy"])
        self._link_into(base, "W", "E", s["g_we"])
        iw, ie, ih, ib = IDX["W"], IDX["E"], IDX["H"], IDX["B"]
        T = np.full(N_NODES, amb + 30.0) if t_guess is None else np.array(t_guess, float)
        q_cool = float(P.sum() + (sum(lin) if lin else fixed_s + fixed_e))
        runaway, converged, iters = False, False, 0
        for iters in range(1, max_iter + 1):
            g_ha, g_ba, g_hc, tc, x_h, x_b = self.external_full(T[ih], T[ib], q_cool)
            G = base.copy()
            b = P.copy()
            G[ih, ih] += g_ha + g_hc
            b[ih] += g_ha * amb + g_hc * tc + x_h
            G[ib, ib] += g_ba
            b[ib] += g_ba * amb + x_b
            if lin is not None:                                              # P = P20 (1 + alpha (T - 20)): linear in T
                G[iw, iw] -= ALPHA_CU * lin[0]
                G[ie, ie] -= ALPHA_CU * lin[1]
                b[iw] += lin[0] * (1 - ALPHA_CU * 20.0)
                b[ie] += lin[1] * (1 - ALPHA_CU * 20.0)
            elif func is not None:
                b[iw] += func(T[iw]) * fs
                b[ie] += func(T[ie]) * (1 - fs)
            else:
                b[iw] += fixed_s
                b[ie] += fixed_e
            try:
                Tn = np.linalg.solve(G, b)
            except np.linalg.LinAlgError:
                runaway = True
                break
            if lin is not None and (Tn.max() > 800.0 or Tn.min() < min(amb, self.t_in) - 50.0):
                runaway = True
                break
            d = float(np.max(np.abs(Tn - T)))
            T = Tn
            if self.cooling.kind == "liquid":
                q_cool = max(0.0, g_hc * (T[ih] - tc))
            if d < tol:
                converged = True
                break
        if runaway:
            T = np.full(N_NODES, 800.0)
        if lin is not None:
            p_s = lin[0] * (1 + ALPHA_CU * (T[iw] - 20.0))
            p_e = lin[1] * (1 + ALPHA_CU * (T[ie] - 20.0))
        elif func is not None:
            p_s, p_e = func(T[iw]) * fs, func(T[ie]) * (1 - fs)
        else:
            p_s, p_e = fixed_s, fixed_e
        return self._result(T, s, p_s, p_e, p_bear, p_wind, losses, iters, runaway, converged, wind_limit_c, magnet_limit_c, q_cool)

    def _hot_spot(self, T: np.ndarray, s: dict, p_end: float) -> tuple[float, float, float]:
        """(winding hot spot, slot hot spot, end-turn nose) in C. Slot: T_W + R_hs Q_walls + the tooth's gradient at the slot's top;
        end turn: T_E + conduction along the copper (a parabola over half the end turn) and across the bundle (3 mm half thickness)."""
        iw, ie, it, iy = IDX["W"], IDX["E"], IDX["T"], IDX["Y"]
        q_walls = s["g_wt"] * (T[iw] - T[it]) + s["g_wy"] * (T[iw] - T[iy])
        q_ty = self.g_ty * (T[it] - T[iy])
        t_slot = T[iw] + s["r_hs"] * q_walls + self.r_tip_minus_mean * q_ty
        qv = p_end / s["v_end"] if s["v_end"] > 0 else 0.0                    # W per m^3 of end-turn copper
        d_nose = qv * (self.l_end_side / 2) ** 2 / (6.0 * MATERIALS["copper"]["k"]) + qv * s["f_loc"] * (3e-3) ** 2 / (6.0 * s["k_eq"])
        t_end = T[ie] + d_nose
        return float(max(t_slot, t_end)), float(t_slot), float(t_end)

    def _result(self, T, s, p_s, p_e, p_bear, p_wind, losses, iters, runaway, converged, wl, ml, q_cool) -> SteadyResult:
        hot, t_slot, t_end = self._hot_spot(T, s, p_e)
        t_m = float(T[IDX["M"]])
        g_ha, g_ba, g_hc, tc, x_h, x_b = self.external_full(T[IDX["H"]], T[IDX["B"]], q_cool)
        heat_out = {"housing_to_ambient": float(g_ha * (T[IDX["H"]] - self.amb) - x_h), "bells_to_ambient": float(g_ba * (T[IDX["B"]] - self.amb) - x_b),
                    "housing_to_coolant": float(g_hc * (T[IDX["H"]] - tc))}
        sources = {"copper_slot": float(p_s), "copper_end": float(p_e), "iron_teeth": float(losses.iron_teeth), "iron_yoke": float(losses.iron_yoke),
                   "iron_rotor": float(losses.iron_rotor), "magnet": float(losses.magnet), "bearing": float(p_bear), "windage": float(p_wind)}
        m_w, m_m = wl - hot, ml - t_m
        return SteadyResult(T={n: float(T[i]) for n, i in IDX.items()}, hot_spot=hot, hot_spot_slot=t_slot, hot_spot_end=t_end, magnet=t_m,
                            copper_loss=float(p_s + p_e), copper_slot=float(p_s), copper_end=float(p_e), sources=sources, heat_out=heat_out,
                            margin_winding=float(m_w), margin_magnet=float(m_m), binding="winding" if m_w < m_m else "magnet",
                            iterations=iters, runaway=runaway, converged=converged)

    # ---------------------------------------------------------------- limits
    def sustainable(self, wind_limit_c: float = CLASS_F_C, magnet_limit_c: float = MAGNET_LIMIT_C, fill: float | None = None,
                    rpm: float | None = None, iron: dict | Callable[[float], dict] | None = None, magnet_loss: float = 0.0, j_hi: float = 40.0, tol: float = 1e-4,
                    losses_at: Callable[[float], Losses] | None = None) -> dict:
        """The continuous current density at which the first limit is reached, by bisection on the full solution (every temperature rises with J).
        The iron losses are the placeholder at `rpm` unless `iron` (a dict like iron_placeholder's) is given; `iron` may also be a FUNCTION of the current density
        j_rms [A/mm^2] returning such a dict (an iron loss that grows with the current: the loss map at the MTPA point). `losses_at(j_rms) -> Losses` replaces the
        whole loss model (a copper loss from another module as a function of the copper temperature, iron loss, magnet loss...); then rpm, iron and magnet_loss are
        not used. j_hi is the largest current density searched (a model of the current that exists only up to some value, such as the flux maps, sets it).
        Returns J_rms [A/mm^2], the copper loss at the operating temperature, the ampere-turns per coil side, the binding limit and the temperatures. If a limit is
        passed with no copper loss at all (the iron, windage and bearing heat alone), J_rms is 0, feasible is False and the margins in "result" are negative."""
        fill = self.winding.fill if fill is None else fill
        iron_fn = iron if callable(iron) else None
        ir = self.iron_placeholder(rpm) if iron is None else (iron_fn(0.0) if iron_fn is not None else iron)
        base = dict(iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"], magnet=magnet_loss)
        guess = [None]

        def losses_for(j: float) -> Losses:
            if losses_at is not None:
                return losses_at(j)
            if iron_fn is None:
                return Losses(copper=self.copper(j, fill), **base)
            w = iron_fn(j)
            return Losses(copper=self.copper(j, fill), iron_teeth=w["iron_teeth"], iron_yoke=w["iron_yoke"], iron_rotor=w["iron_rotor"], magnet=magnet_loss)

        def solve(j: float) -> SteadyResult:
            r = self.steady_state(losses_for(j), fill=fill, wind_limit_c=wind_limit_c, magnet_limit_c=magnet_limit_c, t_guess=guess[0])
            if not r.runaway:
                guess[0] = np.array([r.T[n] for n in NODES])
            return r

        def margin(r: SteadyResult, which: str) -> float:
            if r.runaway:
                return -1.0e3
            return r.margin_winding if which == "winding" else r.margin_magnet

        found = {}
        for which in ("winding", "magnet"):
            r0 = solve(0.0)
            if margin(r0, which) <= 0:                                        # over the limit with no copper loss at all
                found[which] = (0.0, r0)
                continue
            r_hi = solve(j_hi)
            if margin(r_hi, which) > 0:
                found[which] = (j_hi, r_hi)
                continue
            lo, hi = 0.0, j_hi
            while hi - lo > tol * max(1.0, hi):
                mid = 0.5 * (lo + hi)
                if margin(solve(mid), which) > 0:
                    lo = mid
                else:
                    hi = mid
            found[which] = (lo, solve(lo))
        j_w, j_m = found["winding"][0], found["magnet"][0]
        binding = "winding" if j_w < j_m else "magnet"
        j = min(j_w, j_m)
        res = found[binding][1]
        at = self.ampere_turns(j, fill)
        return {"J_rms": j, "feasible": j > 0.0, "J_winding_limit": j_w, "J_magnet_limit": j_m, "binding": binding, "P_cu": res.copper_loss, "AT_pk": at["peak"],
                "AT_rms": at["rms"], "T": res.T, "hot_spot": res.hot_spot, "hot_spot_slot": res.hot_spot_slot, "hot_spot_end": res.hot_spot_end,
                "magnet": res.magnet, "sources": res.sources, "P_total": float(sum(res.sources.values())), "result": res, "fill": fill,
                "wind_limit_c": wind_limit_c, "magnet_limit_c": magnet_limit_c}

    # ---------------------------------------------------------------- the dynamic system
    def _dynamic(self, T: np.ndarray, j_rms: float, fill: float, base_losses: Losses, q_cool: float):
        """(G, b, C, slot branch) of  C dT/dt = -G T + b, linearised at T: the housing's convection at its present temperature, copper's
        linear temperature feedback at the current density j_rms, the constant sources and the heat capacities."""
        s = self._slot_branch(fill)
        amb = self.amb
        G = self.G0.copy()
        self._link_into(G, "W", "T", s["g_wt"])
        self._link_into(G, "W", "Y", s["g_wy"])
        self._link_into(G, "W", "E", s["g_we"])
        g_ha, g_ba, g_hc, tc, x_h, x_b = self.external_full(T[IDX["H"]], T[IDX["B"]], q_cool)
        b = np.zeros(N_NODES)
        ih, ib, iw, ie = IDX["H"], IDX["B"], IDX["W"], IDX["E"]
        G[ih, ih] += g_ha + g_hc
        b[ih] += g_ha * amb + g_hc * tc + x_h
        G[ib, ib] += g_ba
        b[ib] += g_ba * amb + x_b
        p_wind = self.p_windage if base_losses.windage is None else base_losses.windage
        p_bear = self.p_bearing if base_losses.bearing is None else base_losses.bearing
        b[IDX["T"]] += base_losses.iron_teeth + 0.5 * p_wind
        b[IDX["Y"]] += base_losses.iron_yoke
        b[IDX["R"]] += base_losses.iron_rotor + 0.5 * p_wind
        b[IDX["M"]] += base_losses.magnet
        b[IDX["G"]] += p_bear
        p20s, p20e = CopperLoss(j_rms, s["v_slot"], s["v_end"]).p20()
        G[iw, iw] -= ALPHA_CU * p20s
        G[ie, ie] -= ALPHA_CU * p20e
        b[iw] += p20s * (1 - ALPHA_CU * 20.0)
        b[ie] += p20e * (1 - ALPHA_CU * 20.0)
        C = self.C0.copy()
        C[iw], C[ie] = s["c_w"], s["c_e"]
        return G, b, C, s

    def _iron_losses(self, rpm: float | None, iron: dict | None) -> Losses:
        ir = self.iron_placeholder(rpm) if iron is None else iron
        return Losses(copper=0.0, iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"], magnet=0.0)

    def time_constants(self, j_rms: float, fill: float | None = None, rpm: float | None = None, t_min_s: float = 2.0) -> dict:
        """The modes of the network at the operating point of current density j_rms [A/mm^2]: eigenvalues of C^-1 G with G linearised at the steady
        state, copper's temperature feedback included. Returns the time constants tau_i [s] above t_min_s (slowest first) with each mode's weight in
        the response of the winding node to a step in its own loss (v_W^2 tau_i, normalised: how much of the winding's rise each mode carries), and the
        63 % rise times of the winding's hot spot and of the magnets after a step of the current from cold."""
        fill = self.winding.fill if fill is None else fill
        base = self._iron_losses(rpm, None)
        r = self.steady_state(Losses(copper=self.copper(j_rms, fill), iron_teeth=base.iron_teeth, iron_yoke=base.iron_yoke,
                                     iron_rotor=base.iron_rotor), fill=fill)
        T = np.array([r.T[n] for n in NODES])
        G, _, C, _ = self._dynamic(T, j_rms, fill, base, float(sum(r.sources.values())))
        lam, vec = linalg.eigh(G, np.diag(C))
        tau = 1.0 / lam
        order = [i for i in np.argsort(-tau) if tau[i] > t_min_s]
        iw = IDX["W"]
        w = np.array([vec[iw, i] ** 2 * tau[i] for i in order])
        taus = [float(tau[i]) for i in order]
        # 63 % rise times from the step response of the real (nonlinear) transient, from cold
        tr = self.transient(8 * 3600.0, j_rms, start="cold", fill=fill, rpm=rpm, dt_s=1.0, dt_max_s=30.0)
        out = {"tau_s": taus, "winding_weight": [float(x) for x in (w / w.sum())], "operating_point": r}
        for name, series in (("hot_spot", tr.hot_spot), ("magnet", tr.magnet)):
            t0, tf = series[0], series[-1]
            target = t0 + 0.632 * (tf - t0)
            i = int(np.argmax(series >= target))
            out[f"t63_{name}_s"] = float(tr.t[i])
            out[f"rise_{name}_K"] = float(tf - t0)
        return out

    def transient(self, t_end_s: float, j_rms: float | Callable[[float], float], start: str | np.ndarray = "cold", fill: float | None = None,
                  rpm: float | None = None, iron: dict | None = None, dt_s: float = 1.0, dt_max_s: float = 30.0, stop_when=None) -> TransientResult:
        """Heating curves. The copper's current density j_rms [A/mm^2] may be a function of time (an overload profile). start: "cold" (everything at
        the ambient, or the coolant's inlet), "hot" (the steady state at j_rms(0)) or an array of 11 temperatures. The step starts at dt_s and grows
        to dt_max_s once the curves are smooth; each step is exact for frozen coefficients (module docstring, step 12).
        stop_when(t, T, hot_spot) -> True ends the run early."""
        fill = self.winding.fill if fill is None else fill
        base = self._iron_losses(rpm, iron)
        jf = j_rms if callable(j_rms) else (lambda t, v=float(j_rms): v)
        if isinstance(start, str):
            if start == "cold":
                T = np.full(N_NODES, self.t_in if self.cooling.kind == "liquid" else self.amb)
            elif start == "hot":
                r = self.steady_state(Losses(copper=self.copper(jf(0.0), fill), iron_teeth=base.iron_teeth, iron_yoke=base.iron_yoke,
                                             iron_rotor=base.iron_rotor), fill=fill)
                T = np.array([r.T[n] for n in NODES])
            else:
                raise ValueError("start: 'cold', 'hot' or an array of 11 temperatures")
        else:
            T = np.array(start, float)
        s = self._slot_branch(fill)
        ts, Ts, hots, mags, pcus, js = [0.0], [T.copy()], [], [], [], []
        j0 = jf(0.0)
        cl0 = self.copper(j0, fill)
        hots.append(self._hot_spot(T, s, cl0.end(T[IDX["E"]]))[0])
        mags.append(T[IDX["M"]])
        pcus.append(cl0.slot(T[IDX["W"]]) + cl0.end(T[IDX["E"]]))
        js.append(j0)
        t, step = 0.0, dt_s
        q_cool = 100.0
        while t < t_end_s - 1e-9:
            dt = min(step, t_end_s - t)
            j = jf(t + 0.5 * dt)
            G, b, C, _ = self._dynamic(T, j, fill, base, q_cool)
            T_ss = np.linalg.solve(G, b)
            T = T_ss + linalg.expm(-(G / C[:, None]) * dt) @ (T - T_ss)
            t += dt
            cl = CopperLoss(j, s["v_slot"], s["v_end"])
            ts.append(t)
            Ts.append(T.copy())
            hots.append(self._hot_spot(T, s, cl.end(T[IDX["E"]]))[0])
            mags.append(T[IDX["M"]])
            pcus.append(cl.slot(T[IDX["W"]]) + cl.end(T[IDX["E"]]))
            js.append(j)
            if self.cooling.kind == "liquid":
                q_cool = max(1.0, pcus[-1] + sum(base.iron_teeth + base.iron_yoke + base.iron_rotor for _ in (0,)))
            if t > 120.0:
                step = min(dt_max_s, step * 1.15)
            if stop_when is not None and stop_when(t, T, hots[-1]):
                break
        return TransientResult(np.array(ts), np.array(Ts), np.array(hots), np.array(mags), np.array(pcus), np.array(js))

    def overload(self, factor: float, j_cont: float, start: str = "hot", wind_limits=(CLASS_F_C, CLASS_H_C), magnet_limit_c: float = MAGNET_LIMIT_C,
                 fill: float | None = None, rpm: float | None = None, t_max_s: float = 4 * 3600.0, iron: dict | None = None,
                 iron_start: dict | None = None) -> dict:
        """How long the winding can carry `factor` x the continuous current density j_cont, from "cold" (everything at the ambient) or "hot" (the
        steady state at j_cont). For each limit, the first time its temperature is exceeded (None if it never is within t_max_s), by linear
        interpolation between steps (1 s while the curve moves fast, up to 10 s later). The magnets are reported on their own: they heat slowly.
        `iron` (a dict like iron_placeholder's) is the iron loss DURING the overload (default: the placeholder); `iron_start` the iron loss of the hot state before the
        step (default: `iron`): an iron loss that depends on the current, such as the loss map's at the MTPA point, is different at the two currents."""
        fill = self.winding.fill if fill is None else fill
        base = self._iron_losses(rpm, iron)
        if start == "hot":
            base0 = base if iron_start is None else self._iron_losses(rpm, iron_start)
            r0 = self.steady_state(Losses(copper=self.copper(j_cont, fill), iron_teeth=base0.iron_teeth, iron_yoke=base0.iron_yoke,
                                          iron_rotor=base0.iron_rotor), fill=fill)
            start_state = np.array([r0.T[n] for n in NODES])
        elif start == "cold":
            start_state = "cold"
        else:
            raise ValueError("start: 'cold' or 'hot'")

        def done(t, T, hot):
            return hot > wind_limits[-1] + 1.0 and (T[IDX["M"]] > magnet_limit_c + 0.5 or t > 3600.0)
        r = self.transient(t_max_s, factor * j_cont, start=start_state, fill=fill, rpm=rpm, iron=iron, dt_s=1.0, dt_max_s=10.0, stop_when=done)
        if start == "hot":
            r.hot_spot[0] = r0.hot_spot            # the state BEFORE the step: the hot spot's offsets above the node means follow the current at once

        def first(x, limit):
            idx = np.nonzero(x > limit)[0]
            if len(idx) == 0:
                return None
            i = int(idx[0])
            if i == 0:
                return 0.0
            f = (limit - x[i - 1]) / (x[i] - x[i - 1])
            return float(r.t[i - 1] + f * (r.t[i] - r.t[i - 1]))
        out = {"factor": factor, "start": start, "j_rms": factor * j_cont, "hot_spot_start": float(r.hot_spot[0]), "magnet_start": float(r.magnet[0])}
        for lim in wind_limits:
            out[f"winding_{lim:.0f}C_s"] = first(r.hot_spot, lim)
        out[f"magnet_{magnet_limit_c:.0f}C_s"] = first(r.magnet, magnet_limit_c)
        out["magnet_rise_60s"] = float(np.interp(60.0, r.t, r.magnet) - r.magnet[0])
        out["transient"] = r
        return out

    # ---------------------------------------------------------------- the network reduced to a long stack (for the 2D finite-element check)
    def radial_2d(self, p_slot: float, p_teeth: float, p_yoke: float, h_outer: float, fill: float | None = None, t_amb: float | None = None) -> dict:
        """The radial chain only (winding W, teeth T, yoke Y, housing H, the ambient): what a 2D finite-element model of the slot pitch contains.
        No end turns, no air gap, no rotor, nothing leaves axially. p_slot, p_teeth, p_yoke are the heat made [W per metre of stack] in ALL the
        slots, teeth and yoke; h_outer [W/(m^2 K)] the heat transfer coefficient on the housing's outer cylinder (radius r_ho). Returns temperatures
        [C] of W, the slot's hot spot, T, Y, H and the heat split. Uses the same slot, tooth, yoke, contact and housing resistances as the full network."""
        fill = self.winding.fill if fill is None else fill
        s = self._slot_branch(fill)
        amb = self.amb if t_amb is None else t_amb
        L = self.L
        al = MATERIALS["aluminium"]
        r_h2 = math.log(self.r_ho / self.r_hmid) / (2 * math.pi * al["k"] * L)
        g_ext = 1.0 / (r_h2 + 1.0 / (h_outer * 2 * math.pi * self.r_ho * L))
        n = 4
        iw, it, iy, ih = 0, 1, 2, 3
        G = np.zeros((n, n))

        def link(i, j, gc):
            G[i, i] += gc
            G[j, j] += gc
            G[i, j] -= gc
            G[j, i] -= gc
        link(iw, it, s["g_wt"])
        link(iw, iy, s["g_wy"])
        link(it, iy, self.g_ty)
        link(iy, ih, self.g_yh)
        G[ih, ih] += g_ext
        b = np.array([p_slot * L, p_teeth * L, p_yoke * L, g_ext * amb])
        T = np.linalg.solve(G, b)
        q_walls = s["g_wt"] * (T[iw] - T[it]) + s["g_wy"] * (T[iw] - T[iy])
        q_ty = self.g_ty * (T[it] - T[iy])
        hot = T[iw] + s["r_hs"] * q_walls + self.r_tip_minus_mean * q_ty
        t_tip = T[it] + self.r_tip_minus_mean * q_ty                           # the tooth's tip (the bore's end of the fin) above its mean
        return {"T_W": float(T[iw]), "T_hot": float(hot), "T_T": float(T[it]), "T_tip": float(t_tip), "T_Y": float(T[iy]), "T_H": float(T[ih]),
                "q_flank": float(s["g_wt"] * (T[iw] - T[it])), "q_bottom": float(s["g_wy"] * (T[iw] - T[iy])), "k_eq": s["k_eq"], "f_loc": s["f_loc"],
                "r_hs_per_w": s["r_hs"], "G_ext": g_ext}

    # ---------------------------------------------------------------- the picture of the network
    def edges(self, fill: float | None = None, dt_outside: float = 40.0) -> list[tuple[str, str, float]]:
        """Every link as (node, node, R [K/W]) at the model's present coefficients (the outside evaluated at a `dt_outside` K rise)."""
        fill = self.winding.fill if fill is None else fill
        s = self._slot_branch(fill)
        e = [("W", "T", 1.0 / s["g_wt"]), ("W", "Y", 1.0 / s["g_wy"]), ("W", "E", 1.0 / s["g_we"]), ("T", "Y", 1.0 / self.g_ty),
             ("Y", "H", 1.0 / self.g_yh), ("E", "A", 1.0 / self.g_ea), ("E", "B", 1.0 / self.g_eb_direct), ("A", "B", 1.0 / self.g_ab),
             ("A", "R", 1.0 / self.g_ar), ("A", "S", 1.0 / self.g_as), ("T", "A", 1.0 / self.g_ta), ("T", "R", self.r_tr), ("R", "M", self.r_mr),
             ("R", "S", self.r_rs), ("S", "G", self.r_sg), ("G", "B", self.r_gb), ("B", "H", self.r_bh)]
        g_ha, g_ba, g_hc, _ = self.external(self.amb + dt_outside, self.amb + dt_outside, 100.0)
        if g_ha > 0:
            e.append(("H", "ambient", 1.0 / g_ha))
        if g_hc > 0:
            e.append(("H", "coolant", 1.0 / g_hc))
        e.append(("B", "ambient", 1.0 / g_ba))
        return e


# ---------------------------------------------------------------------------------------------------- module-level interface
@functools.lru_cache(maxsize=64)
def _model(overrides_key: str, stack_mm: float, cooling: Cooling, speed_rpm: float, fill: float, assume: Assumptions) -> ThermalModel:
    ov = json.loads(overrides_key) if overrides_key else None
    return ThermalModel(ov, stack_mm, cooling, speed_rpm, WindingSpec(fill), assume)


def get_model(overrides: dict | None = None, stack_mm: float = 20.0, cooling: Cooling | None = None, speed_rpm: float = RATED_RPM,
              fill: float = 0.45, assume: Assumptions | None = None) -> ThermalModel:
    """The cached ThermalModel for these inputs (the module's functions use it, so repeated calls cost only the solve)."""
    key = "" if overrides is None else json.dumps(overrides, sort_keys=True)
    return _model(key, float(stack_mm), cooling or Cooling.natural(), float(speed_rpm), float(fill), assume or Assumptions())


def steady_state(losses: Losses, cooling: Cooling, stack_mm: float, *, overrides: dict | None = None, speed_rpm: float = RATED_RPM,
                 fill: float = 0.45, assume: Assumptions | None = None, **kw) -> SteadyResult:
    """steady_state(losses, cooling, stack_mm, ...): see the module docstring. The copper loss in `losses` may be a number, a function of the
    copper temperature, or a CopperLoss (get_model(...).copper(j_rms)). Extra keywords go to ThermalModel.steady_state."""
    return get_model(overrides, stack_mm, cooling, speed_rpm, fill, assume).steady_state(losses, **kw)


def transient(j_rms, cooling: Cooling, stack_mm: float, t_end_s: float, *, overrides: dict | None = None, speed_rpm: float = RATED_RPM,
              fill: float = 0.45, assume: Assumptions | None = None, start="cold", **kw) -> TransientResult:
    """transient(j_rms, cooling, stack_mm, t_end_s, ...): heating curves at a current density j_rms [A/mm^2] (a number or a function of time),
    from `start` ('cold', 'hot' or an array of 11 temperatures)."""
    return get_model(overrides, stack_mm, cooling, speed_rpm, fill, assume).transient(t_end_s, j_rms, start=start, **kw)


def sustainable(cooling: Cooling, stack_mm: float, *, overrides: dict | None = None, speed_rpm: float = RATED_RPM, fill: float = 0.45,
                assume: Assumptions | None = None, **kw) -> dict:
    """The continuous current density, copper loss, ampere-turns and binding limit (ThermalModel.sustainable)."""
    return get_model(overrides, stack_mm, cooling, speed_rpm, fill, assume).sustainable(**kw)


def overload(cooling: Cooling, stack_mm: float, factor: float, j_cont: float, *, overrides: dict | None = None, speed_rpm: float = RATED_RPM,
             fill: float = 0.45, assume: Assumptions | None = None, **kw) -> dict:
    """How long `factor` x j_cont can be carried (ThermalModel.overload)."""
    return get_model(overrides, stack_mm, cooling, speed_rpm, fill, assume).overload(factor, j_cont, **kw)


def time_constants(cooling: Cooling, stack_mm: float, j_rms: float, *, overrides: dict | None = None, speed_rpm: float = RATED_RPM,
                   fill: float = 0.45, assume: Assumptions | None = None, **kw) -> dict:
    """The network's time constants at j_rms (ThermalModel.time_constants)."""
    return get_model(overrides, stack_mm, cooling, speed_rpm, fill, assume).time_constants(j_rms, **kw)


def shaft_fan_power(rpm, cooling: Cooling | None = None, stack_mm: float = 60.0, *, overrides: dict | None = None, assume: Assumptions | None = None):
    """The shaft power [W] of a shaft fan at rpm (a number or an array): the air speed in the fin channels is v = v_ref rpm / n_ref, the pressure the channels and the
    cowl ask for times the flow, over the fan's efficiency (ThermalModel.fan_state); P goes as rpm^3 (the fan laws) up to the small change of the friction factor with
    Re. cooling: Cooling.shaft_fan(...) (default 6 m/s at 3,000 rpm). A load on the shaft like the windage: pass it to the motor map as the fan's power
    (`fan=lambda r: thermal.shaft_fan_power(r, cooling, 60.0)`)."""
    c = cooling or Cooling.shaft_fan()
    if c.kind != "shaft_fan":
        raise ValueError(f"shaft_fan_power: a Cooling.shaft_fan(...) is needed, not {c.kind!r}")
    m = get_model(overrides, stack_mm, c, RATED_RPM, 0.45, assume)
    r = np.asarray(rpm, float)
    out = np.array([m.fan_state(c.v_ref_ms * float(x) / c.n_ref_rpm)["p_shaft_w"] for x in np.ravel(r)]).reshape(r.shape)
    return float(out) if out.ndim == 0 else out


def size_shaft_fan(budget_w: float, rpm_top: float = MAX_SPEED_RPM, n_ref_rpm: float = RATED_RPM, stack_mm: float = 60.0, *, overrides: dict | None = None,
                   assume: Assumptions | None = None, ambient_c: float = AMBIENT_C) -> dict:
    """The shaft fan whose power at the top speed rpm_top is `budget_w`: the air speed v_top at which fan_state's shaft power equals the budget (bisection), and
    v_ref = v_top n_ref / rpm_top, the speed it gives at n_ref_rpm (a fan big enough for the rated speed costs a lot at the top speed, P ~ n^3; a fan held to a
    budget at the top speed is weak at the rated speed). Returns {"budget_w", "v_top", "v_ref", "p_ref_w" (its power at n_ref_rpm), "cooling": Cooling.shaft_fan(v_ref, n_ref)}."""
    probe = Cooling.shaft_fan(1.0, n_ref_rpm, ambient_c)
    m = get_model(overrides, stack_mm, probe, RATED_RPM, 0.45, assume)
    lo, hi = 0.0, 400.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if m.fan_state(mid)["p_shaft_w"] < budget_w:
            lo = mid
        else:
            hi = mid
    v_top = 0.5 * (lo + hi)
    v_ref = v_top * n_ref_rpm / rpm_top
    cooling = Cooling.shaft_fan(v_ref, n_ref_rpm, ambient_c)
    return {"budget_w": budget_w, "v_top": v_top, "v_ref": v_ref, "p_ref_w": m.fan_state(v_ref)["p_shaft_w"], "cooling": cooling}


# ---------------------------------------------------------------------------------------------------- self-checks
def selfcheck(verbose: bool = True) -> list[tuple[str, bool, str]]:
    """The building blocks against textbook cases and against the project's earlier numbers. Returns (name, passed, detail) for each."""
    out: list[tuple[str, bool, str]] = []

    def check(name: str, ok: bool, detail: str) -> None:
        out.append((name, bool(ok), detail))
    tm, tx, fx = rect_poisson(1.0, 1.0)
    check("rectangle series, square: T_max = 0.07367 q a^2/k", abs(tx - 0.073671) < 2e-5, f"T_max = {tx:.6f}, T_mean = {tm:.6f}, wall share {fx:.4f}")
    tm, tx, fx = rect_poisson(1.0, 100.0)
    check("rectangle series, long strip: T_max -> 1/8, T_mean -> 1/12", abs(tx - 0.125) < 5e-4 and abs(tm - 1 / 12) < 3e-3, f"T_max = {tx:.5f}, T_mean = {tm:.5f}")
    nu = 0.128 * 6858.0 ** 0.367
    check("Becker-Kaye vortex branch at Ta_m = 6858 gives Nu = 3.27 (Howey et al.'s worked example)", abs(nu - 3.27) < 0.01, f"Nu = {nu:.3f}")
    n1a, n1b = 2 * 0.0178 / math.log(1.0178), 0.128 * 1700.0 ** 0.367
    n2a, n2b = 0.128 * 1.0e4 ** 0.367, 0.409 * 1.0e4 ** 0.241
    check("air-gap Nusselt branches join (laminar / vortex / turbulent)", abs(n1a / n1b - 1) < 0.05 and abs(n2a / n2b - 1) < 0.01,
          f"{n1a:.3f} vs {n1b:.3f} at 1,700; {n2a:.3f} vs {n2b:.3f} at 1e4")
    hnc = h_natural_cylinder(0.138, 80.0, 40.0)
    check("natural convection of a 138 mm cylinder at dT = 40 K: h = 4-8 W/(m2 K)", 4.0 < hnc < 8.0, f"h = {hnc:.2f}")
    g = build_geom(design_params())
    check("slot area 72.44 mm2 (two layers of 36.218)", abs(g.a_slot * 1e6 - 72.436) < 0.01, f"{g.a_slot * 1e6:.3f} mm2, layer {g.a_layer * 1e6:.3f}")
    check("yoke 4.161 mm at its thinnest", abs((g.r_so - g.r_yb) * 1e3 - 4.1608) < 1e-3, f"{(g.r_so - g.r_yb) * 1e3:.4f} mm")
    check("tooth 3.0 mm at the shoe, 3.6 mm at the root", abs(g.c_tooth * g.r_shoe * 1e3 - 2.995) < 0.01 and abs(g.c_tooth * g.r_wall * 1e3 - 3.624) < 0.01,
          f"{g.c_tooth * g.r_shoe * 1e3:.3f} / {g.c_tooth * g.r_wall * 1e3:.3f} mm")
    gr = build_geom(design_params({}))
    l_mt = 2 * 0.020 + math.pi * gr.span + 4 * 2e-3
    p_cu = copper_resistivity(100.0) * (100.0 / (0.45 * gr.a_layer * 1e6) / math.sqrt(2) * 1e6) ** 2 * 72 * 0.45 * gr.a_layer * (l_mt / 2)
    check("reference design: DC loss at 0.45 fill, 100 C, 100 A-turns = the winding study's 56.6 W", abs(p_cu - 56.604) < 0.1,
          f"{p_cu:.3f} W (mean turn {l_mt * 1e3:.1f} mm)")
    m = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air())
    r = m.steady_state(m.default_losses(copper=50.0))
    q_in, q_out = sum(r.sources.values()), sum(r.heat_out.values())
    check("energy balance (forced air, 40 mm): heat out = heat made", abs(q_out - q_in) < 1e-6 * q_in, f"{q_out:.4f} W out, {q_in:.4f} W in")
    cl = m.copper(4.0)
    r1 = m.steady_state(Losses(copper=cl, iron_teeth=2.0, iron_yoke=2.0))
    r2 = m.steady_state(Losses(copper=lambda t: cl.total(t), iron_teeth=2.0, iron_yoke=2.0))
    check("copper loss as a matrix term = the same loss iterated as a function of T", abs(r1.T["W"] - r2.T["W"]) < 1e-3 and abs(r1.copper_loss - r2.copper_loss) < 1e-3,
          f"T_W {r1.T['W']:.4f} vs {r2.T['W']:.4f}; P {r1.copper_loss:.4f} vs {r2.copper_loss:.4f}")
    ir = m.iron_placeholder()
    tr = m.transient(8 * 3600.0, 4.0, start="cold", dt_max_s=60.0)
    ss = m.steady_state(Losses(copper=m.copper(4.0), iron_teeth=ir["iron_teeth"], iron_yoke=ir["iron_yoke"], iron_rotor=ir["iron_rotor"]))
    check("transient (8 h, forced air) ends at the steady state", abs(tr.T[-1, IDX["W"]] - ss.T["W"]) < 0.2, f"{tr.T[-1, IDX['W']]:.3f} vs {ss.T['W']:.3f} C")
    # a single-resistor limit: huge internal conductances, everything at one temperature, T = T_amb + P / G_outside
    big = ThermalModel(stack_mm=40.0, cooling=Cooling.forced_air())
    big.G0 = big.G0 * 1.0e6
    for key in list(big._slot_cache):
        big._slot_cache.pop(key)
    s = big._slot_branch(0.45)
    s["g_wt"] *= 1e6
    s["g_wy"] *= 1e6
    s["g_we"] *= 1e6
    rr = big.steady_state(Losses(copper=30.0, copper_slot_fraction=0.5, bearing=0.0, windage=0.0))
    g_ha, g_ba, _, _ = big.external(rr.T["H"], rr.T["B"], 0.0)
    spread = max(rr.T.values()) - min(rr.T.values())
    check("single-resistor limit: with huge internal conductances every node sits at T_amb + P / G_outside", spread < 0.01 and abs((rr.T["W"] - big.amb) - 30.0 / (g_ha + g_ba)) < 0.05,
          f"spread {spread:.4f} K, rise {rr.T['W'] - big.amb:.3f} K against {30.0 / (g_ha + g_ba):.3f} K")
    # the shaft fan: the fan laws, the standstill, the energy balance with the air's warming, and the heat flows
    sa = ThermalModel(stack_mm=60.0, cooling=Cooling.shaft_fan(6.0, 3000.0), speed_rpm=3000.0)
    p1, p2 = sa.fan_state(8.0)["p_shaft_w"], sa.fan_state(16.0)["p_shaft_w"]
    check("shaft fan: doubling the air speed (and so the shaft speed) multiplies the power by 7.4 to 8.0 (n^3, the friction factor falling a little)", 7.4 < p2 / p1 < 8.0,
          f"P(8 m/s) = {p1:.2f} W, P(16 m/s) = {p2:.2f} W, ratio {p2 / p1:.3f}")
    st = ThermalModel(stack_mm=60.0, cooling=Cooling.shaft_fan(6.0, 3000.0), speed_rpm=0.0)
    pp = st.outside_parts(80.0, 80.0)
    want = h_natural_cylinder(st.d_ho, 80.0, st.amb) * st.A_cylH * st.a.fan_fin_gain_conv
    check("shaft fan at standstill: no forced convection; the housing's convection is the plain cylinder's times the fin gain", st.v_air == 0.0 and abs(pp["h_conv"] / want - 1) < 1e-12,
          f"v {st.v_air} m/s, {pp['h_conv']:.4f} against {want:.4f} W/K")
    rr = sa.steady_state(Losses(copper=60.0, iron_teeth=10.0, iron_yoke=6.0, iron_rotor=2.0))
    hf = sa.heat_flows(rr)
    check("shaft fan: heat out = heat made, the air's warming by the fan included; the exits split into convection and radiation", abs(hf["balance_w"]) < 1e-3 and hf["exits"]["housing_radiation"] > 0,
          f"imbalance {hf['balance_w']:.2e} W; housing radiation {hf['exits']['housing_radiation']:.2f} W, convection {hf['exits']['housing_convection']:.2f} W; stream warmed {sa.dT_fan:.2f} K")
    if verbose:
        for name, ok, detail in out:
            print(("PASS  " if ok else "FAIL  ") + name + "  [" + detail + "]")
    return out


def main() -> int:
    import time
    results = selfcheck()
    failed = [n for n, ok, _ in results if not ok]
    print()
    m = ThermalModel(stack_mm=60.0, cooling=Cooling.natural())
    t0 = time.perf_counter()
    for _ in range(200):
        m.steady_state(m.default_losses(copper=m.copper(3.0)))
    print(f"steady_state call (natural convection, 60 mm): {(time.perf_counter() - t0) / 200 * 1e3:.3f} ms")
    for c in (Cooling.natural(), Cooling.forced_air(), Cooling.liquid()):
        s = ThermalModel(stack_mm=60.0, cooling=c).sustainable()
        print(f"{c.label:60s} J = {s['J_rms']:.2f} A/mm2  Pcu = {s['P_cu']:.1f} W  binding: {s['binding']}  hot spot {s['hot_spot']:.1f} C  magnets {s['magnet']:.1f} C")
    print("selfcheck:", "all passed" if not failed else "FAILED: " + ", ".join(failed))
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
