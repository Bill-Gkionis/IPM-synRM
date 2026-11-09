"""optimization.thermal_fan_figures: the four pictures of the shaft-fan study (matplotlib PNGs drawn from the rows thermal_fan_study.py computes; no solver).

What it is, and why
    The shaft fan's story is told by four pictures, each with its numbers in a CSV beside it (the table view the data-visualization method asks for): what the fan is, what the
    cooling gives at each speed, where the heat goes, and what each cooling lets the motor do. thermal_fan_study.main() builds the rows and calls make_all.

The pictures (each has a title, a note under it saying what is drawn, and a legend or direct labels)
    fig_fan_model.png          left: the air speed in the fin channels against the shaft speed for the two shaft fans and the fixed-speed fan; middle: the fans' shaft power
                               (log scale: the cube law) with the 50 W budget; right: the sizing trade, the continuous ampere-turns at four speeds against the budget at the top speed.
    fig_fan_limit_vs_speed.png the continuous limit (ampere-turns a coil side, peak) against speed from standstill to 16,330 rpm for the five coolings, N45 and the 150 C grade, and the
                               MTPA torque at the limit. A cross on the axis: a limit is passed with no current at all.
    fig_fan_heat_paths.png     where the heat goes at each cooling's limit (60 mm): where it is made, through which surfaces it leaves, and the flows inside (slot winding to teeth and yoke,
                               end turns through the slot or into the end-space air, the rotor across the air gap, into the end-space air or through the shaft), all in per cent of the heat made.
    fig_fan_capability.png     per cooling against stack length at 3,000 rpm: the continuous MTPA torque, the mechanical power (an estimate) and the time from cold at 2x the
                               current to the first limit; N45 solid, the 150 C grade dashed.

The look (the palette of optimization/winding_figures.py, as thermal_figures.py)
    Surface #fcfcfb, ink #0b0b0b and #52514e, grid #e1e0d9, axes #c3c2b7; series blue #2a78d6, orange #eb6834, aqua #1baf7a; context gray #898781. Natural convection is blue,
    the shaft fan (a) orange, the shaft fan (b) a lighter orange, the fixed-speed forced air dark gray (dotted) and the liquid jacket aqua; every series also has its own marker shape
    (aqua is 2.7:1 on the surface), and the lines are named in the legend and the CSVs.

Run
    Not on its own: thermal_fan_study.py calls make_all (uv run python -m optimization.thermal_fan_study).
Numbers to expect: four PNGs, 1430 px wide, drawn in a few seconds.

Self-check (the answers follow).
  1. Why is the fan's power a straight line of slope 3 on the log-log plot?
  2. Why does the 50 W fan do nothing for the motor at 3,000 rpm that the finned housing without a fan does not?
Answers.
  1. P goes as n^3 (the fan laws: the flow as n, the pressure as n^2), so log P = 3 log n + a constant; the small bend is the channel's friction factor falling with Re.
  2. At 3,000 rpm it blows 1.8 m/s, laminar flow in the channels: its forced convection is the same order as the natural convection of the fins, so the gain is a fraction;
     the heat it can carry away is limited by what the air stream can take at that flow, m cp (T_wall - T_air).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paths
from optimization import thermal
from optimization import winding_figures as wf
from optimization.thermal import MAX_SPEED_RPM, RATED_RPM, Cooling

SURFACE, INK, INK2, MUTED, GRID, AXIS = wf.SURFACE, wf.INK, wf.INK2, wf.MUTED, wf.GRID, wf.AXIS
BLUE, ORANGE, AQUA = wf.BLUE, wf.ORANGE, wf.AQUA


def mix(hex_a: str, hex_b: str, w: float) -> str:
    """w of colour a over colour b (both '#rrggbb')."""
    a = np.array([int(hex_a[i:i + 2], 16) for i in (1, 3, 5)], float)
    b = np.array([int(hex_b[i:i + 2], 16) for i in (1, 3, 5)], float)
    c = np.round(w * a + (1 - w) * b).astype(int)
    return "#" + "".join(f"{v:02x}" for v in c)


ORANGE_LIGHT = mix(ORANGE, SURFACE, 0.55)
STYLE = {"natural": (BLUE, "o", "-", "natural convection (plain housing)"),
         "shaft_a": (ORANGE, "s", "-", "shaft fan (a): 6 m/s at 3,000 rpm"),
         "shaft_b": (ORANGE_LIGHT, "D", "-", "shaft fan (b): 50 W at 16,330 rpm"),
         "forced": (INK2, "v", ":", "fixed-speed forced air (6 m/s)"),
         "liquid": (AQUA, "^", "--", "liquid jacket (reference)")}
KEYS = tuple(STYLE)
OUT = lambda: paths.path("out") / "optimization" / "thermal"
MAIN = 60.0


def _legend(ax, handles, loc="upper left", ncol=1, fs=8.2) -> None:
    leg = ax.legend(handles=handles, loc=loc, frameon=False, fontsize=fs, ncol=ncol)
    for t in leg.get_texts():
        t.set_color(INK2)


def _fig_legend(fig, handles, y: float, ncol: int, fs: float = 8.4) -> None:
    leg = fig.legend(handles=handles, loc="center", bbox_to_anchor=(0.5, y), frameon=False, fontsize=fs, ncol=ncol)
    for t in leg.get_texts():
        t.set_color(INK2)


def _spread(ys: list[float], gap: float) -> list[float]:
    """Label heights pushed apart to at least `gap` (in data units), keeping their order."""
    order = sorted(range(len(ys)), key=lambda i: ys[i])
    out = list(ys)
    for a, b in zip(order, order[1:]):
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    return out


def _cool_handles(extra=()) -> list:
    h = [Line2D([0], [0], color=c, lw=2.2, marker=m, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.2, linestyle=ls, label=lab) for c, m, ls, lab in STYLE.values()]
    return h + list(extra)


# ---------------------------------------------------------------------------------------------------- 1. the fan and its sizing trade
def _fan_p(fan_rows: list[dict], sizing: str, rpm: int) -> float:
    return next(r["p_shaft_w"] for r in fan_rows if r["sizing"] == sizing and r["rpm"] == rpm)


def fig_fan_model(path: Path, fan_rows: list[dict], sizing_rows: list[dict], budget_w: float) -> Path:
    pa3, pa16 = _fan_p(fan_rows, "shaft_a", RATED_RPM), _fan_p(fan_rows, "shaft_a", MAX_SPEED_RPM)
    pb3, pb16 = _fan_p(fan_rows, "shaft_b", RATED_RPM), _fan_p(fan_rows, "shaft_b", MAX_SPEED_RPM)
    fig = wf._figure(14.3, 5.2, "The shaft fan: its air speed and power against shaft speed, and what its sizing trades (60 mm stack, the forced-air preset's fins)",
                     "Left: the air speed in the 36 fin channels is proportional to the shaft speed (the fixed-speed fan is the preset's 6 m/s). Middle: the fan's shaft power, the pressure the "
                     "channels and the cowl ask for times the flow over a fan efficiency of 0.20 (assumed): it goes as the cube of the speed. Right: the continuous ampere-turns a coil side "
                     "(peak; N45, class F, MTPA, real losses) against the fan's power budget at 16,330 rpm; the point at 5 W is the finned housing with no fan. A fan big enough for 3,000 rpm "
                     f"(a) takes {pa3:.1f} W there and {pa16:,.0f} W at the top speed; one held to {pb16:.0f} W at the top speed (b) takes {pb3:.2f} W at 3,000 rpm and is a weak fan there.")
    gs = fig.add_gridspec(1, 3, left=0.05, right=0.985, top=0.84, bottom=0.225, wspace=0.24)
    ax1, ax2, ax3 = (fig.add_subplot(gs[0, i]) for i in range(3))
    sa = Cooling.shaft_fan(6.0, float(RATED_RPM))
    sb = next(r for r in fan_rows if r["sizing"] == "shaft_b")["v_ref_ms"]
    cb = Cooling.shaft_fan(sb, float(RATED_RPM))
    rpm = np.linspace(0.0, MAX_SPEED_RPM, 120)
    wf._style(ax1, "shaft speed [rpm]", "air speed in the fin channels [m/s]")
    ax1.plot(rpm, sa.v_ref_ms * rpm / sa.n_ref_rpm, color=STYLE["shaft_a"][0], lw=2.2)
    ax1.plot(rpm, cb.v_ref_ms * rpm / cb.n_ref_rpm, color=STYLE["shaft_b"][0], lw=2.2)
    ax1.plot(rpm, np.full_like(rpm, 6.0), color=STYLE["forced"][0], lw=1.8, linestyle=":")
    for c, v_ref, col, mk in ((sa, 6.0, STYLE["shaft_a"][0], "s"), (cb, sb, STYLE["shaft_b"][0], "D")):
        for n in (RATED_RPM, MAX_SPEED_RPM):
            ax1.plot([n], [v_ref * n / RATED_RPM], marker=mk, color=col, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5)
            ax1.text(n - 250, v_ref * n / RATED_RPM + (1.2 if n == MAX_SPEED_RPM else 0.9), f"{v_ref * n / RATED_RPM:.1f}", color=INK2, fontsize=8.5, ha="right" if n == MAX_SPEED_RPM else "center")
    ax1.set_xlim(0, MAX_SPEED_RPM * 1.03)
    ax1.set_ylim(0, 36)
    ax1.text(MAX_SPEED_RPM * 0.52, 3.0, "fixed 6 m/s (external fan)", color=INK2, fontsize=8.5)
    ax1.text(MAX_SPEED_RPM * 0.30, 15.0, "(a)", color=STYLE["shaft_a"][0], fontsize=10, fontweight="bold")
    ax1.text(MAX_SPEED_RPM * 0.80, 9.2, "(b)", color=mix(ORANGE, INK, 0.75), fontsize=10, fontweight="bold")
    # power
    wf._style(ax2, "shaft speed [rpm]", "fan shaft power [W] (log)")
    n_arr = np.linspace(300.0, MAX_SPEED_RPM, 100)
    pa = thermal.shaft_fan_power(n_arr, sa, MAIN, assume=_assume())
    pb = thermal.shaft_fan_power(n_arr, cb, MAIN, assume=_assume())
    ax2.plot(n_arr, pa, color=STYLE["shaft_a"][0], lw=2.2)
    ax2.plot(n_arr, pb, color=STYLE["shaft_b"][0], lw=2.2)
    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.axhline(budget_w, color=MUTED, lw=1.2, linestyle=(0, (4, 3)))
    ax2.text(330, budget_w * 1.18, f"budget {budget_w:.0f} W at 16,330 rpm", color=INK2, fontsize=8.5)
    for n in (RATED_RPM, MAX_SPEED_RPM):
        for col, mk, c in ((STYLE["shaft_a"][0], "s", sa), (STYLE["shaft_b"][0], "D", cb)):
            p = float(thermal.shaft_fan_power(n, c, MAIN, assume=_assume()))
            ax2.plot([n], [p], marker=mk, color=col, markersize=8, markeredgecolor=SURFACE, markeredgewidth=1.5)
            lab = f"{p:,.0f} W" if p >= 100 else f"{p:.3g} W"
            if n == RATED_RPM:
                ax2.text(n * 0.93, p * 1.0, lab, color=INK2, fontsize=8.5, ha="right", va="center")
            else:
                ax2.text(n * 0.93, p * (1.0 if col == STYLE["shaft_a"][0] else 0.7), lab, color=INK2, fontsize=8.5, ha="right", va="center")
    ax2.set_xlim(300, MAX_SPEED_RPM * 1.15)
    ax2.set_xticks([500, 1000, 3000, 10000])
    ax2.set_xticklabels(["500", "1,000", "3,000", "10,000"])
    ax2.set_ylim(0.05, 6000)
    ax2.text(5000, 1.6, "slope 3: P ~ n^3", color=INK2, fontsize=8.5, rotation=0)
    # the trade
    wf._style(ax3, "fan power budget at 16,330 rpm [W] (log)", "continuous ampere-turns a coil side, peak")
    spd = {1000: (MUTED, "o"), 3000: (ORANGE, "s"), 6000: (BLUE, "^"), 10000: (AQUA, "v")}
    rows = [r for r in sizing_rows if r["label"] != "no fan"]
    nofan = next(r for r in sizing_rows if r["label"] == "no fan")
    x = [r["budget_w"] for r in rows]
    ends = {}
    for n, (col, mk) in spd.items():
        y = [r[f"AT_{n}"] for r in rows]
        ax3.plot(x, y, color=col, lw=2.0, marker=mk, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.2)
        ax3.plot([5.0], [nofan[f"AT_{n}"]], marker=mk, color=col, markersize=7, markerfacecolor=SURFACE, markeredgewidth=1.6)
        ends[n] = y[-1]
    lab_y = _spread([ends[n] for n in spd], 9.0)
    for (n, (col, mk)), ly in zip(spd.items(), lab_y):
        ax3.text(x[-1] * 1.12, ly, f"{n:,} rpm", color=INK2, fontsize=8.5, va="center")
    ax3.set_xscale("log")
    ax3.axvline(budget_w, color=MUTED, lw=1.2, linestyle=(0, (4, 3)))
    ax3.text(budget_w * 1.08, 4, "(b) 50 W", color=INK2, fontsize=8.5, rotation=90, va="bottom")
    ax3.text(x[-1] * 0.93, 4, "(a)", color=INK2, fontsize=8.5, rotation=90, va="bottom", ha="right")
    ax3.set_xlim(4.0, x[-1] * 3.2)
    ax3.set_ylim(bottom=0)
    ax3.set_xticks([5, 10, 50, 200, 1000])
    ax3.set_xticklabels(["no fan", "10", "50", "200", "1,000"])
    return wf._save(fig, path)


def _assume():
    from optimization import thermal_machine as tm
    return tm.assumptions()


# ---------------------------------------------------------------------------------------------------- 2. the limit against speed
def _series(rows: list[dict], key: str, magnets: str, stack: float, field: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rr = sorted([r for r in rows if r["cooling"] == key and r["magnets"] == magnets and r["stack_mm"] == stack], key=lambda r: r["rpm"])
    x = np.array([r["rpm"] for r in rr], float)
    y = np.array([0.0 if r[field] is None else r[field] for r in rr], float)
    feas = np.array([r["feasible"] for r in rr], bool)
    return x, y, feas


def fig_limit_vs_speed(path: Path, rows: list[dict], at_max: float, fan_rows: list[dict]) -> Path:
    pa3, pa16 = _fan_p(fan_rows, "shaft_a", RATED_RPM), _fan_p(fan_rows, "shaft_a", MAX_SPEED_RPM)
    pb3, pb16 = _fan_p(fan_rows, "shaft_b", RATED_RPM), _fan_p(fan_rows, "shaft_b", MAX_SPEED_RPM)
    fig = wf._figure(14.3, 5.4, "The continuous limit against speed, standstill to 16,330 rpm (60 mm stack, class F, 40 C ambient, MTPA currents, real losses)",
                     "Left and middle: the ampere-turns a coil side (peak) at which the first limit is reached, for the N45 magnets (80 C) and for the 150 C grade; a cross on the axis: a limit is "
                     "passed with no current at all (the iron, windage and bearing heat alone); the dotted line is the edge of the flux maps (3x rated). Right: the MTPA torque at the N45 limit. "
                     f"The fan of (a) takes {pa3:.1f} W at 3,000 rpm and {pa16:,.0f} W at 16,330 rpm, (b) {pb3:.2f} W and {pb16:.0f} W. Operating point: MTPA at every speed; above the corner speed "
                     "the real motor weakens the field and has less iron loss, so these are thermal limits at MTPA currents there (exact at low speed). Copper loss from the copper-loss model, "
                     "iron loss from the iron-loss map.")
    gs = fig.add_gridspec(1, 3, left=0.052, right=0.985, top=0.855, bottom=0.305, wspace=0.2)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, mag, field, ylab in ((axes[0], "N45", "AT_pk", "continuous ampere-turns a coil side, peak (N45, 80 C)"),
                                 (axes[1], "150C", "AT_pk", "continuous ampere-turns a coil side, peak (150 C grade)"),
                                 (axes[2], "N45", "torque_Nm", "continuous MTPA torque [N.m] (N45)")):
        wf._style(ax, "shaft speed [rpm]", ylab)
        for key in KEYS:
            col, mk, ls, _ = STYLE[key]
            x, y, feas = _series(rows, key, mag, MAIN, field)
            ax.plot(x, y, color=col, lw=2.0 if key != "forced" else 1.8, linestyle=ls, marker=mk, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1.2)
            bad = ~feas
            if bad.any():
                ax.plot(x[bad], np.zeros(bad.sum()), linestyle="none", marker="x", color=INK, markersize=8, markeredgewidth=1.8)
        ax.set_xlim(-300, MAX_SPEED_RPM + 500)
        ax.set_xticks([0, 4000, 8000, 12000, 16330])
        ax.set_xticklabels(["0", "4,000", "8,000", "12,000", "16,330"])
        if field == "AT_pk":
            ax.axhline(at_max, color=MUTED, lw=1.1, linestyle=(0, (1.5, 3)))
            ax.set_ylim(-12, at_max * 1.06)
        else:
            ax.set_ylim(-0.5, None)
    _fig_legend(fig, _cool_handles([Line2D([0], [0], color=INK, lw=0, marker="x", markersize=8, markeredgewidth=1.8, label="a limit is passed with no current")]), 0.2, 3)
    return wf._save(fig, path)


# ---------------------------------------------------------------------------------------------------- 3. where the heat goes
MADE_GROUPS = [("copper in the slots", ("made_copper_slot_w",), ORANGE), ("copper in the end turns", ("made_copper_end_w",), mix(ORANGE, SURFACE, 0.55)),
               ("iron: teeth and yoke", ("made_iron_teeth_w", "made_iron_yoke_w"), BLUE), ("iron: rotor", ("made_iron_rotor_w",), mix(BLUE, SURFACE, 0.5)),
               ("windage and bearings", ("made_windage_w", "made_bearing_w"), MUTED), ("magnets (eddy loss, 0)", ("made_magnet_w",), AQUA)]
OUT_GROUPS = [("housing: convection", ("out_housing_convection_w",), INK2), ("housing: radiation", ("out_housing_radiation_w",), mix(INK2, SURFACE, 0.45)),
              ("end bells: convection", ("out_bells_convection_w",), AQUA), ("end bells: radiation", ("out_bells_radiation_w",), mix(AQUA, SURFACE, 0.45)),
              ("housing to coolant (jacket)", ("out_housing_to_coolant_w",), BLUE)]
INSIDE_SLOT = [("slot winding to the teeth", ("grp_slot_to_teeth_w",), BLUE), ("slot winding to the yoke", ("grp_slot_to_yoke_w",), mix(BLUE, SURFACE, 0.5))]
INSIDE_END = [("end turns through the slot", ("grp_end_turns_via_slot_w",), ORANGE), ("end turns into the end-space air", ("grp_end_turns_to_air_w",), mix(ORANGE, SURFACE, 0.55)),
              ("end turns to the bells", ("grp_end_turns_to_bells_w",), MUTED)]
INSIDE_ROT = [("rotor across the air gap to the stator", ("grp_rotor_across_gap_w",), BLUE), ("rotor into the end-space air", ("grp_rotor_to_end_space_air_w",), mix(ORANGE, SURFACE, 0.55)),
              ("rotor through the shaft to the bearings", ("grp_rotor_through_shaft_w",), AQUA)]
CASE_LABEL = {"natural": "natural convection", "shaft_a": "shaft fan (a)", "shaft_b": "shaft fan (b)", "forced": "fixed forced air", "liquid": "liquid jacket"}


def _stack(ax, rows, groups, y0, height, total_key="total_made_w", diverging=False, label_min=7.0) -> None:
    """Draw one stacked bar per row; segments in per cent of the heat made; negative segments go to the left of zero when diverging."""
    for i, r in enumerate(rows):
        tot = r[total_key]
        pos, neg = 0.0, 0.0
        for name, keys, col in groups:
            v = 100.0 * sum(r.get(k, 0.0) for k in keys) / tot if tot > 0 else 0.0
            if abs(v) < 1e-9:
                continue
            if diverging and v < 0:
                ax.barh(y0[i], v, left=neg, height=height, color=col, edgecolor=SURFACE, linewidth=0.8)
                if abs(v) >= label_min:
                    ax.text(neg + v / 2, y0[i], f"{v:.0f}", ha="center", va="center", fontsize=7, color=INK)
                neg += v
            else:
                ax.barh(y0[i], v, left=pos, height=height, color=col, edgecolor=SURFACE, linewidth=0.8)
                if v >= label_min:
                    ax.text(pos + v / 2, y0[i], f"{v:.0f}", ha="center", va="center", fontsize=7, color=INK if col not in (INK2, ORANGE, BLUE, AQUA) else SURFACE)
                pos += v


def fig_heat_paths(path: Path, heat_rows: list[dict]) -> Path:
    fig = wf._figure(14.3, 9.8, "Where the heat goes at each cooling's continuous limit (60 mm stack, N45, class F), in per cent of the heat made in the motor",
                     "One row per case (cooling, shaft speed, and the ampere-turns a coil side at the limit). Left: where the heat is made. Second: through which surfaces it leaves the motor "
                     "(all of it; the end bells carry the heat that reaches them through the end-space air, the bearings and the housing joint). Third: inside the stator, the slot "
                     "winding's heat to the teeth and the yoke, and the end turns' heat through the slot, into the end-space air, or to the bells. Right: the rotor's heat (its own and what "
                     "the stator gives it across the gap, a negative bar: the stator heats the rotor) into the end-space air and through the shaft to the bearings. The fan's power is a "
                     "shaft load that heats the air and is not in these sources. Numbers in the segments are per cent; the table is fan_heat_paths.csv.")
    n = len(heat_rows)
    gs = fig.add_gridspec(1, 4, left=0.155, right=0.985, top=0.92, bottom=0.335, wspace=0.07, width_ratios=[1.0, 1.0, 1.0, 0.75])
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    ys = np.arange(n)[::-1].astype(float)
    titles = ["heat made", "heat leaving the motor", "inside the stator", "the rotor's heat"]
    for ax, t in zip(axes, titles):
        wf._style(ax, "per cent of the heat made", "" if ax is not axes[0] else "")
        ax.set_title(t, color=INK, fontsize=10, fontweight="bold", loc="left")
        ax.set_ylim(-0.6, n - 0.4)
        ax.grid(axis="y", visible=False)
    _stack(axes[0], heat_rows, MADE_GROUPS, ys, 0.62)
    _stack(axes[1], heat_rows, OUT_GROUPS, ys, 0.62)
    # the stator: two bars a row (slot to teeth/yoke above, end turns below), as per cent of the total heat made
    _stack(axes[2], heat_rows, INSIDE_SLOT, ys + 0.17, 0.31, diverging=True, label_min=9.0)
    _stack(axes[2], heat_rows, INSIDE_END, ys - 0.17, 0.31, diverging=True, label_min=9.0)
    _stack(axes[3], heat_rows, INSIDE_ROT, ys, 0.62, diverging=True, label_min=1.5)
    axes[0].set_xlim(0, 100)
    axes[1].set_xlim(0, 100)
    axes[2].set_xlim(-12, 100)
    axes[3].set_xlim(-6, 8)
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels([f"{CASE_LABEL[r['cooling']]}, {r['rpm']:,} rpm\n{r['AT_pk']:.0f} A-turns, {r['total_made_w']:.0f} W" for r in heat_rows], fontsize=8, color=INK2)
    for ax in axes[1:]:
        ax.set_yticks([])
    axes[2].axvline(0, color=AXIS, lw=0.8)
    axes[3].axvline(0, color=AXIS, lw=0.8)
    for ax, groups in zip(axes, (MADE_GROUPS, OUT_GROUPS, INSIDE_SLOT + INSIDE_END, INSIDE_ROT)):
        h = [Patch(facecolor=c, edgecolor=SURFACE, label=name) for name, _, c in groups]
        leg = ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, -0.085), frameon=False, fontsize=7.6, ncol=1)
        for t in leg.get_texts():
            t.set_color(INK2)
    return wf._save(fig, path)


# ---------------------------------------------------------------------------------------------------- 4. the capability
def fig_capability(path: Path, cap_rows: list[dict]) -> Path:
    fig = wf._figure(14.3, 5.3, "What each cooling lets the motor do at 3,000 rpm (40 C ambient, class F): the continuous point and 2x the current from cold",
                     "Solid: N45 magnets (80 C); dashed: the 150 C grade (an open decision; the flux maps are N45's, so its torque would be a little lower than shown). Left: the continuous MTPA "
                     "torque from the flux maps, scaled with the stack. Middle: the mechanical power T x omega, an ESTIMATE (MTPA, sinusoidal currents, no voltage limit, the maps' d-q torque). "
                     "Right: the time from cold at 2x the continuous current until the first limit (class F winding or magnets); a ring marks a case where the winding limit comes first, "
                     "and a point at the top is a case that never reaches a limit in four hours. The torque at 2x the current (from the maps; beyond them for the liquid jacket) is in fan_capability.csv.")
    gs = fig.add_gridspec(1, 3, left=0.055, right=0.985, top=0.84, bottom=0.315, wspace=0.2)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    wf._style(axes[0], "stack length [mm] (rows of magnets x 20 mm)", "continuous MTPA torque [N.m]")
    wf._style(axes[1], "stack length [mm]", "mechanical power T x omega [kW] (estimate)")
    wf._style(axes[2], "stack length [mm]", "time to the first limit at 2x the current, from cold [s] (log)")
    # the fixed-speed forced air is drawn first and wide: the shaft fan (a) lies on top of it at 3,000 rpm (the same 6 m/s), and both must be seen
    for key in ("forced", "natural", "shaft_b", "shaft_a", "liquid"):
        col, mk, ls, _ = STYLE[key]
        wide = key == "forced"
        lw = 5.0 if wide else 2.0
        for mag, dash in (("N45", "-"), ("150C", (0, (4, 2.5)))):
            rr = sorted([r for r in cap_rows if r["cooling"] == key and r["magnets"] == mag], key=lambda r: r["stack_mm"])
            x = [r["stack_mm"] for r in rr]
            ms = 9 if wide else 7
            fc = col if mag == "N45" else SURFACE
            axes[0].plot(x, [r["torque_Nm"] or 0.0 for r in rr], color=col, lw=lw, linestyle=dash, marker=mk, markersize=ms, markeredgecolor=SURFACE if mag == "N45" else col,
                         markeredgewidth=1.2 if mag == "N45" else 1.6, markerfacecolor=fc, alpha=0.55 if wide else 1.0)
            axes[1].plot(x, [r["p_em_w"] / 1000.0 for r in rr], color=col, lw=lw, linestyle=dash, marker=mk, markersize=ms, markeredgecolor=SURFACE if mag == "N45" else col,
                         markeredgewidth=1.2 if mag == "N45" else 1.6, markerfacecolor=fc, alpha=0.55 if wide else 1.0)
            t = [r["two_x_t_first_s"] if r["two_x_t_first_s"] is not None else 4 * 3600.0 for r in rr]
            axes[2].plot(x, t, color=col, lw=lw, linestyle=dash, marker=mk, markersize=ms, markeredgecolor=SURFACE if mag == "N45" else col,
                         markeredgewidth=1.2 if mag == "N45" else 1.6, markerfacecolor=fc, alpha=0.55 if wide else 1.0)
            for xi, r, ti in zip(x, rr, t):
                if r["two_x_first"] == "winding":
                    axes[2].plot([xi], [ti], marker="o", markersize=13, markerfacecolor="none", markeredgecolor=col, markeredgewidth=1.0)
    axes[2].set_yscale("log")
    for ax in axes:
        ax.set_xticks([20, 40, 60, 80])
        ax.set_xlim(15, 85)
    axes[0].set_ylim(bottom=0)
    axes[1].set_ylim(bottom=0)
    extra = [Line2D([0], [0], color=INK2, lw=2.0, linestyle="-", label="N45 (80 C)"), Line2D([0], [0], color=INK2, lw=2.0, linestyle=(0, (4, 2.5)), label="150 C grade")]
    _fig_legend(fig, _cool_handles(extra), 0.2, 4)
    return wf._save(fig, path)


def make_all(fan_rows, sizing_rows, speed_rows, heat_rows, cap_rows, budget_w: float = 50.0, at_max: float = 355.88) -> list[Path]:
    OUT().mkdir(parents=True, exist_ok=True)
    return [fig_fan_model(OUT() / "fig_fan_model.png", fan_rows, sizing_rows, budget_w),
            fig_limit_vs_speed(OUT() / "fig_fan_limit_vs_speed.png", speed_rows, at_max, fan_rows),
            fig_heat_paths(OUT() / "fig_fan_heat_paths.png", heat_rows),
            fig_capability(OUT() / "fig_fan_capability.png", cap_rows)]
