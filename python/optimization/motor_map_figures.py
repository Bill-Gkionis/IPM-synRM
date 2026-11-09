"""optimization.motor_map_figures: the figures of the motor map (matplotlib PNGs, no FEMM): the torque-speed envelopes, the highest speed each fan allows, the continuous limit against the thermal
study's table, the losses along the envelope, the efficiency maps, the current plane, the stack table, the winding study's front and the 150 C magnet grade.

What it is, and why
    motor_map_run.py solves and checks; this file draws what it found. Each picture has its numbers in a CSV beside it (fig_<name>.csv, the "table view" a picture alone cannot give), a title
    that says what it shows and a note under it that says what is drawn, what to look at, the inputs and what is not modelled. The inputs in the note come from the run's own record
    (results.pkl), so a run on other windings says so without a change here; a run's inputs are put back into motor_map_run (restore_inputs) before anything is drawn.

The pictures (written to out\\optimization\\motor_map\\<label>\\figures\\; 130 dpi), for the winding study's winding of the headline stack (60 mm) unless said
    fig_envelopes_fans.png        torque and power against speed for the three sizings of the shaft fan, both buses, continuous (the thermal network's limit at every speed) and short-term
                                  (twice the winding's own continuous current at 3,000 rpm, from a cold start), the corner marked, the rotor's safe speed (16,330 rpm) as the mechanical limit.
    fig_envelopes_comparisons.png the same for natural convection, the fixed 6 m/s external fan and the liquid jacket (parked).
    fig_speed_limit.png           the highest speed each fan allows: the motor's air-gap torque less its iron drag against the fan + windage + bearings load, the crossing marked; a cross where a
                                  continuous envelope ends because no current both holds the voltage and keeps the temperatures.
    fig_limits.png                the continuous ampere-turns of each fan against speed: this module's with each bus's voltage limit, and the thermal study's table (MTPA, no voltage limit)
                                  as squares: where and why they differ.
    fig_losses_<fan>.png          the losses stacked against speed along the envelope: copper DC, copper AC (skin, proximity, PWM ripple, circulating currents), iron by region, mechanical
                                  (bearings, windage, the shaft fan).
    fig_efficiency_<fan>.png      the efficiency over (speed, torque) at minimum-loss control, the continuous envelope solid, the short-term one dashed.
    fig_saving_<fan>.png          what minimum-loss control saves over the MTPA control, in watts (for the headline fan).
    fig_plane_<fan>_<stack>.png   the (i_d, i_q) plane of the 20 mm stack's own winding: the current circle, constant-torque lines, the voltage curve at several speeds, the MTPA path, the
                                  envelope's operating points by region.
    fig_stack_table_<fan>.png     the stack table of each fan: magnets to buy, masses, rotor inertia, and for the winding study's winding of each stack and for the headline winding on every
                                  stack the corner speed, peak torque and power, the continuous torque and where it ends, and the overload time from a cold start.
    fig_front_48V.png, fig_front_60V.png  the winding study's front at 60 mm (the headline winding and the 8 / 4, 5 / 4, 4 / 4 and 3 / 4 windings of the 100, 150, 200 and 300 A inverter
                                  classes) with the headline fan: continuous and short-term torque and power on that bus.
    fig_grade.png                 the 8 / 4, 4 x 0.8 mm winding for a magnet grade that may run at 150 C against the headline winding on the N45 grade, both buses.
    (<fan> is fanA, fanB or fanC: the three sizings; the headline fan is fanB, held to 50 W at the top speed.)

The look (the palette of optimization/maps_figures.py: surface #fcfcfb, ink #0b0b0b and #52514e; series blue #2a78d6, orange #eb6834, aqua #1baf7a (with a second channel: its own marker), context gray
    #898781; one sequential blue ramp for the maps and for the front's windings (darker = more current); lines 2 px, markers at least 8 px with a 2 px surface ring; text is never in a series
    colour; one axis a panel)

Run (from the Python to FEMM folder; after motor_map_run has written results.pkl; seconds)
    uv run python -m optimization.motor_map_figures [label]       (the label is the run's folder under out\\optimization\\motor_map\\, default "final")

Self-check (the answers follow).
  1. In fig_plane the voltage curves shrink as the speed rises. Around which point, and why there?
  2. The continuous torque curves of the three shaft fans cross. What decides on which side of a crossing each fan is better?
  3. In fig_speed_limit the load of the fan, the windage and the bearings is drawn as a torque and it grows as the square of the speed, while the fan's power grows as the cube. Why?
  4. The continuous limit of fan (b) falls to a valley near 5,900 rpm and rises again (fig_limits). Why there?
Answers.
  1. Around the point where the flux linkage is zero (i_d of about -psi_PM / L_d): the speed voltage is w_e |psi|, so |v| = V_max is the curve |psi| = V_max / w_e, a smaller and smaller
     contour around the zero of |psi|.
  2. The fan on the shaft blows in proportion to the speed and costs power as its cube: the one sized for 3,000 rpm cools best at low and middle speeds and eats the shaft power at the top
     speed; the ones held to 50 W and 200 W at the top speed are cheap there and cool weakly below it. Which is better at a speed is decided by whether the heat that the extra air takes
     away there is worth more than the shaft power the fan takes.
  3. Torque is power over angular speed: a power that goes as n^3 asks for a torque that goes as n^3 / n = n^2.
  4. The air in the fin channels moves in proportion to the speed (1.76 m/s at 3,000 rpm for fan (b)), and the flow is laminar below Re 2,300 (3.45 m/s: 5,880 rpm): there the heat-transfer
     coefficient does not rise with the speed (Nu is constant) while the iron loss does, so the limit falls; above it the flow turns turbulent (fully at Re 4,000, 6 m/s, 10,230 rpm), the
     coefficient rises several times, and the limit rises again. Fan (a) has its valley at 1,725 rpm and fan (c) at 3,690 rpm: the same 3.45 m/s in the channels.
"""



from __future__ import annotations

import csv
import math
import pickle
import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.tri as mtri
import numpy as np

from optimization import motor_map as mm
from optimization.maps_figures import AQUA, AXIS, BLUE, CONTEXT, FONT, GRID, INK, INK2, MUTED, ORANGE, PT, SEQ, SURFACE, TOP, _cbar, _dot, _legend, _line, _style
from optimization.motor_map_run import OUT_ROOT

COOL_KEYS: list = []
LABELS: dict = {}                      # cooling key -> label, filled from the run's own record by load()
SHORT: dict = {}
TOP_RPM = mm.TOP_RPM
_INPUTS: dict = {}


def provisional_text() -> str:
    """The statement of the inputs, from the run's own record (results['inputs']), printed under every figure."""
    i = _INPUTS
    ws = i["windings"]["opt"]
    wtxt = "; ".join(f"{k} mm {_wl(w)}" for k, w in ws.items())
    return (f"Inputs: the winding study's windings ({wtxt}), stack {i['headline_stack']:g} mm unless said, bus {' and '.join(f'{v:g}' for v in i['vdc'])} V with {100 * i['k_margin']:.0f} % of the "
            f"modulation limit (an assumption), copper: {i['copper']}, iron loss of the maps at the CAD's {i['lamination_mm']:g} mm lamination, mechanical loss estimated, the maps' magnets at N45 "
            f"({i['magnet_limit_c']:g} C). The power figures are estimates.")


def not_modelled() -> str:
    ac = "the winding's AC loss (this run prices the copper as DC loss)" if _INPUTS.get("copper_kind", "dc") == "dc" else "the circulating currents of strands that are not in a compact bundle"
    return f"Not modelled: the magnets' eddy loss, the inverter's loss, the end-winding leakage, demagnetisation, {ac}."


def bus_of(vdc: float) -> mm.Bus:
    return mm.Bus(vdc, _INPUTS["k_margin"])


def load(out: Path) -> dict:
    with open(out / "results.pkl", "rb") as f:
        res = pickle.load(f)
    _INPUTS.clear()
    _INPUTS.update(res["inputs"])
    from optimization import motor_map_run as run
    run.restore_inputs(res["inputs"])
    COOL_KEYS[:] = _INPUTS["coolings"]
    LABELS.clear()
    LABELS.update({k: v["label"] for k, v in _INPUTS["cool_labels"].items()})
    SHORT.clear()
    SHORT.update({k: v["short"] for k, v in _INPUTS["cool_labels"].items()})
    return res


def _figure(width: float, height: float, title: str, note: str):
    """A figure with its title (top left) and its note (bottom left, wrapped by hand to the width); fig._bottom is the bottom edge the panels must keep for the note and the axis labels."""
    plt.rcParams["font.family"] = FONT
    fig = plt.figure(figsize=(width, height), facecolor=SURFACE)
    chars = width * 10.0                                                         # the characters of a bold 11 pt title that fit the width; a longer title shrinks to fit
    fig.text(0.012, 0.985, title, ha="left", va="top", fontsize=11.0 if len(title) <= chars else max(8.0, 11.0 * chars / len(title)), color=INK, fontweight="bold")
    lines = textwrap.wrap(note, width=int(width * 16.0))
    fig.text(0.012, 0.012, "\n".join(lines), ha="left", va="bottom", fontsize=7.6, color=INK2, linespacing=1.35)
    fig._bottom = min(0.42, (0.012 * height + len(lines) * 0.135 + 0.62) / height)
    return fig


def _save(fig, out: Path, name: str, rows: list[dict] | None = None) -> Path:
    folder = out / "figures"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / f"fig_{name}.png"
    fig.savefig(path, dpi=130, facecolor=SURFACE)
    plt.close(fig)
    if rows:
        with open(folder / f"fig_{name}.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(dict.fromkeys(k for r in rows for k in r)), restval="")
            w.writeheader()
            w.writerows(rows)
    return path


def _grid(fig, nrow: int, ncol: int, left: float = 0.055, right: float = 0.985, gap_x: float = 0.05, gap_y: float = 0.105, cbar: float = 0.0) -> list[list]:
    """A grid of axes above the figure's reserved bottom margin; cbar: room at each panel's right for a colour bar and its labels (the gap grows with it)."""
    gap_x = gap_x + (0.045 if cbar else 0.0)
    if cbar:
        right = min(right, 0.935)                                  # the last colour bar and its labels sit to the right of its panel
    w = (right - left - gap_x * (ncol - 1)) / ncol
    h = (TOP - fig._bottom - gap_y * (nrow - 1)) / nrow
    return [[fig.add_axes([left + c * (w + gap_x), TOP - (r + 1) * h - r * gap_y, w - cbar, h]) for c in range(ncol)] for r in range(nrow)]


def _vline(ax, x: float, text: str | None = None) -> None:
    ax.axvline(x, color=CONTEXT, linewidth=1.2 * PT, linestyle=":")
    if text:
        ax.text(x, 0.985, text, transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=7.2, color=INK2, rotation=90)


def _thousands(ax) -> None:
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:,.0f}"))


def _case(res: dict, cool: str, stack: float, variant: str, vdc: float) -> dict:
    return res["cases"][f"{cool}|{stack:g}|{variant}|{vdc:g}"]


def _buses(res: dict) -> tuple:
    """The run's bus voltages and their colours: the lowest in blue, the highest in orange (a run with one bus draws it once)."""
    v = list(res["inputs"]["vdc"])
    return tuple((x, c) for x, c in zip((v[0], v[-1]) if len(v) > 1 else (v[0],), (BLUE, ORANGE)))


def _head(res: dict) -> float:
    return float(res["inputs"]["headline_stack"])


# ====================================================================================================== figure 1: the envelopes
def _hv(res: dict) -> str:
    return res["inputs"]["headline_variant"]


def _wl(w: dict) -> str:
    """A winding of the run's record in the winding study's notation: N / a, strands x wire."""
    return f"{w['turns']:g} / {w['paths']}, {w['strands']:g} x {w['wire_mm']:g} mm"


def _head_winding(res: dict) -> str:
    return _wl(res["inputs"]["windings"][_hv(res)][f"{_head(res):g}"])


def fig_envelopes(res: dict, out: Path, group: str = "fans") -> Path:
    H, buses, variant = _head(res), _buses(res), _hv(res)
    wanted = _INPUTS["fan_keys"] if group == "fans" else [k for k in COOL_KEYS if k not in _INPUTS["fan_keys"]]
    keys = [k for k in wanted if f"{k}|{H:g}|{variant}|{buses[0][0]:g}" in res["cases"]]
    what = "the three sizings of the shaft fan" if group == "fans" else "the comparisons"
    title = f"Torque and power against speed, {what}: the winding study's winding of the {H:g} mm stack ({_head_winding(res)})"
    note = (f"Top row: shaft torque; bottom row: shaft power. Solid: continuous (the thermal network's limit at every speed: magnets at or under {_INPUTS['magnet_limit_c']:g} C, winding under 155 C). "
            f"Dashed: short-term, twice the winding's own continuous current at 3,000 rpm with no voltage limit, from a COLD start (how long it lasts: stack table). Blue {buses[0][0]:g} V, orange "
            f"{buses[-1][0]:g} V; the dot is the corner speed of the short-term curve. Dotted vertical line: the rotor's safe speed, 16,330 rpm, the end of the range. A continuous curve that drops "
            f"to zero and comes back is a band of speeds where no current both holds the voltage limit and keeps the temperatures. "
            + ("" if group == "fans" else "The liquid jacket's current is held at the edge of the flux maps (3 times the placeholder's rated current), not by its temperature. ")
            + f"{provisional_text()} {not_modelled()}")
    fig = _figure(max(10.0, 4.4 * len(keys)), 8.4, title, note)
    axes = _grid(fig, 2, len(keys))
    rows = []
    t_max = max(float(np.nanmax(_case(res, k, H, variant, v)[kind]["d"]["T"])) for k in keys for v, _ in buses for kind in ("cont", "short"))
    p_max = max(float(np.nanmax(np.nan_to_num(_case(res, k, H, variant, v)[kind]["d"]["P_out"]))) for k in keys for v, _ in buses for kind in ("cont", "short")) / 1e3
    for c, k in enumerate(keys):
        st = res["static"][f"{k}|{H:g}|{variant}"]
        for r in range(2):
            ax = axes[r][c]
            _style(ax, "speed [rpm]", "shaft torque [N.m]" if r == 0 else "shaft power [kW]")
            for vdc, col in buses:
                case = _case(res, k, H, variant, vdc)
                for kind, ls in (("cont", "-"), ("short", "--")):
                    d = case[kind]["d"]
                    y = d["T"] if r == 0 else np.nan_to_num(d["P_out"]) / 1e3
                    ax.plot(d["rpm"], y, color=col, linestyle=ls, linewidth=2 * PT)
                    if r == 0:
                        for j in range(len(y)):
                            rows.append({"cooling": LABELS[k], "vdc": vdc, "kind": kind, "rpm": float(d["rpm"][j]), "T_Nm": float(d["T"][j]), "P_kW": float(np.nan_to_num(d["P_out"][j]) / 1e3)})
                sm = case["short"]["summary"]
                if sm["corner_rpm"] is not None:
                    yc = sm["corner_T"] if r == 0 else sm["corner_T"] * sm["corner_rpm"] * 2 * math.pi / 60 / 1e3
                    ax.plot([sm["corner_rpm"]], [yc], marker="o", markersize=6, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT, linestyle="")
            _vline(ax, TOP_RPM, "rotor's safe speed" if (r == 0 and c == 0) else None)
            ax.set_xlim(0, 17200)
            ax.set_ylim(0, (t_max if r == 0 else p_max) * 1.06)
            if r == 0:
                ax.set_title(f"{SHORT[k]}\nthermal rating {st['rating']['I']:.1f} A, peak {st['i_peak']:.1f} A{' (maps edge)' if st['peak_clipped'] else ''}", fontsize=8.6, color=INK, loc="left")
            _thousands(ax)
    _legend(axes[0][-1], [_line(c, f"{v:g} V") for v, c in buses] + [_line(INK2, "continuous"), _line(INK2, "short-term (cold start)", "--")], loc="upper right")
    return _save(fig, out, f"envelopes_{group}", rows)


# ====================================================================================================== figure 1b: the highest speed each fan allows
def fig_speed_limit(res: dict, out: Path) -> Path:
    H, buses, variant = _head(res), _buses(res), _hv(res)
    keys = _INPUTS["fan_keys"]
    stat = res["static"][f"{keys[0]}|{H:g}|{variant}"]
    w = mm.Winding(stat["winding"]["turns"], stat["winding"]["paths"], stat["winding"]["strands"], stat["winding"]["wire_mm"], stat["winding"]["wire_over_mm"])
    title = f"The highest speed each fan allows: the motor's torque against the loads on its shaft ({_head_winding(res)}, {H:g} mm)"
    note = ("Colored: the air-gap torque of the best operating point of the envelope (field weakening above the corner) less the drag of its iron loss, solid for the continuous envelope and dashed for "
            "the short-term one. Black: the torque the fan, the windage and the bearings ask for (the fan's grows as the square of the speed). Where the two meet (a dot) the shaft torque is zero: that "
            "speed is the highest the motor can reach with that fan on its shaft. A curve that reaches the right-hand edge is not stopped by the fan. A cross (a t after the speed in the "
            "title) is a continuous curve that ends for another reason: above that speed no current both holds the voltage limit and keeps the magnets and the winding inside their "
            f"temperatures. {provisional_text()} {not_modelled()}")
    fig = _figure(max(10.0, 4.4 * len(keys)), 5.8, title, note)
    axes = _grid(fig, 1, len(keys))
    rows = []
    speeds = np.linspace(0.0, TOP_RPM, 400)
    for c, k in enumerate(keys):
        ax = axes[0][c]
        _style(ax, "speed [rpm]", "torque [N.m]")
        cool = next(x for x in _fan_functions(res) if x[0] == k)
        m = mm.Machine(w, H, mech=mm.Mechanical(fan=cool[1]))
        load = m.mech_torque(2.0 * math.pi * speeds / 60.0)
        ax.plot(speeds, load, color=INK, linewidth=2 * PT)
        ymax = 0.0
        for vdc, col in buses:
            case = _case(res, k, H, variant, vdc)
            for kind, ls in (("cont", "-"), ("short", "--")):
                d = case[kind]["d"]
                on = d["region"] != "off"
                net = (d["T_em"] - d["P_fe"] / (2.0 * math.pi * np.maximum(d["rpm"], 1e-9) / 60.0))
                net = np.where(on, net, np.nan)
                ax.plot(d["rpm"], net, color=col, linestyle=ls, linewidth=2 * PT)
                ymax = max(ymax, float(np.nanmax(net[on])) if on.any() else 0.0)
                lp, stop = case["limit_" + kind], case[kind]["summary"]["stop"]
                if lp is not None:
                    ax.plot([lp["rpm"]], [lp["T_em"] - lp["T_iron"]], marker="o", markersize=6.5, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT, linestyle="")
                elif stop == "limit" and on.any():
                    j_end = int(np.nonzero(on)[0].max())
                    ax.plot([d["rpm"][j_end]], [net[j_end]], marker="X", markersize=7.5, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT, linestyle="")
                for j in range(0, len(d["rpm"]), 4):
                    if on[j]:
                        rows.append({"fan": LABELS[k], "vdc": vdc, "kind": kind, "rpm": float(d["rpm"][j]), "T_air_gap_less_iron_Nm": float(net[j]), "T_load_Nm": float(np.interp(d["rpm"][j], speeds, load))})
        ax.set_xlim(0, TOP_RPM * 1.03)
        ax.set_ylim(0, max(ymax, float(load.max()) * 1.05) * 1.05)
        _thousands(ax)
        _vline(ax, TOP_RPM)
        lims = []
        for vdc, _ in buses:
            case = _case(res, k, H, variant, vdc)
            ends = []
            for kd in ("cont", "short"):
                sm = case[kd]["summary"]
                ends.append("top" if sm["stop"] == "top" else (f"{sm['zero_rpm']:,.0f} t" if sm["stop"] == "limit" else f"{sm['zero_rpm']:,.0f}"))
            lims.append(f"{vdc:g} V: " + " / ".join(ends))
        ax.set_title(f"{SHORT[k]}\nhighest speed, continuous / short-term [rpm]\n" + "   ".join(lims), fontsize=8.2, color=INK, loc="left")
    _legend(axes[0][-1], [_line(INK, "fan + windage + bearings")] + [_line(c, f"{v:g} V") for v, c in buses] + [_line(INK2, "continuous"), _line(INK2, "short-term", "--")], loc="upper right")
    return _save(fig, out, "speed_limit", rows)


def _fan_functions(res: dict) -> list:
    """(key, the fan's power function of the speed on the headline stack) for the shaft fans of the run."""
    from optimization import motor_map_run as run
    H = _head(res)
    return [(c["key"], run.fan_of(c, H)) for c in run.coolings() if c["key"] in _INPUTS["fan_keys"]]


# ====================================================================================================== figure 1c: the continuous limit against the thermal study's table
def fig_limits(res: dict, out: Path) -> Path | None:
    """The continuous ampere-turns of each shaft fan against speed: this module's with the voltage limit of each bus, and the thermal study's (MTPA, no voltage limit) at its six speeds."""
    ref = res.get("thermal_reference")
    H, buses, variant = _head(res), _buses(res), _hv(res)
    keys = _INPUTS["fan_keys"]
    if ref is None:
        return None
    title = f"The continuous limit of each shaft fan against speed: this module (voltage limit) and the thermal study (MTPA, no voltage limit), {H:g} mm, {_head_winding(res)}"
    note = ("Peak ampere-turns a coil side that the cooling carries for ever (magnets at or under " + f"{_INPUTS['magnet_limit_c']:g} C, winding under 155 C). Lines: this module, at the operating point the voltage limit leaves (MTPA below the "
            "corner speed, field weakening above it, so less flux and less iron loss), blue and orange for the two buses; a line that falls to zero is a speed where no current both holds the voltage "
            "and keeps the temperatures. Squares: the thermal study's table, thermal_machine.continuous_limit at MTPA currents with no voltage limit (it keeps the full flux and so the full iron loss "
            "above the corner: a lower bound there), at 0, 1,000, 3,000, 6,000, 10,000 and 16,330 rpm. Where the lines are above the squares the field weakening has made room; where they "
            f"fall to zero the voltage cannot be held within the temperatures. {provisional_text()}")
    fig = _figure(max(10.0, 4.4 * len(keys)), 5.4, title, note)
    axes = _grid(fig, 1, len(keys))
    rows = []
    for c, k in enumerate(keys):
        ax = axes[0][c]
        _style(ax, "speed [rpm]", "continuous ampere-turns, peak [A]")
        stat = res["static"][f"{k}|{H:g}|{variant}"]
        ymax = 1.0
        for vdc, col in buses:
            case = _case(res, k, H, variant, vdc)
            at = np.where(case["limits"]["limited"] == "none", 0.0, case["limits"]["I"] * stat["winding"]["turns"] / stat["winding"]["paths"])
            rpm = case["cont"]["d"]["rpm"]
            ax.plot(rpm, at, color=col, linewidth=2 * PT)
            ymax = max(ymax, float(at.max()))
            for j in range(0, len(rpm), 4):
                rows.append({"fan": LABELS[k], "vdc": vdc, "rpm": float(rpm[j]), "AT_pk": float(at[j])})
        pts = [r for r in ref["rows"] if r["cooling"] == SHORT[k]]
        ax.plot([r["rpm"] for r in pts], [r["tm_headline"] for r in pts], linestyle="", marker="s", markersize=7, markerfacecolor=INK, markeredgecolor=SURFACE, markeredgewidth=2 * PT)
        ymax = max([ymax] + [r["tm_headline"] for r in pts])
        for r in pts:
            rows.append({"fan": LABELS[k], "vdc": "thermal study (no voltage limit)", "rpm": r["rpm"], "AT_pk": r["tm_headline"]})
        ax.set_xlim(0, TOP_RPM * 1.03)
        ax.set_ylim(0, ymax * 1.08)
        _thousands(ax)
        _vline(ax, TOP_RPM)
        ax.set_title(SHORT[k], fontsize=8.6, color=INK, loc="left")
    _legend(axes[0][0], [_line(c, f"{v:g} V, with the voltage limit") for v, c in buses] + [_dot(INK, "s", 7, "thermal study: MTPA, no voltage limit")], loc="lower left")
    return _save(fig, out, "limits", rows)


# ====================================================================================================== figure 2: the losses along the envelope
LAYERS = [("P_cu_dc", "copper, DC", ORANGE), ("P_cu_ac", "copper, AC, PWM and circulating currents", "#f6cdb8"), ("P_fe_teeth", "iron, teeth", BLUE),
          ("P_fe_yoke", "iron, yoke", "#9ec5f2"), ("P_fe_rotor", "iron, rotor", "#14478a"), ("P_bearing", "mechanical, bearings", "#c3c2b7"), ("P_windage", "mechanical, windage", MUTED),
          ("P_fan", "mechanical, shaft fan", INK2)]


def fig_losses(res: dict, out: Path, cool: str) -> Path:
    H, buses, variant = _head(res), _buses(res), _hv(res)
    title = f"Where the loss goes along the envelope: {LABELS[cool]}, {_head_winding(res)}, {H:g} mm"
    ac = ("The copper AC band is empty on purpose: this run prices the copper as DC loss." if _INPUTS.get("copper_kind", "dc") == "dc" else
          "The copper AC band (skin, proximity, PWM ripple, circulating currents) is the thin light-orange line: with a compact bundle of thin strands it is a fraction of a percent of the copper loss.")
    note = ("Each band is a loss at the operating point the envelope uses, stacked: copper (orange), iron by region (blues), mechanical (grays; the shaft fan is the dark one). " + ac +
            f" Top row continuous, bottom row short-term (cold start); each row has its own scale. {provisional_text()} {not_modelled()}")
    fig = _figure(12.5, 8.2, title, note)
    axes = _grid(fig, 2, 2)
    rows = []
    ymax = {}
    for kind in ("cont", "short"):
        top = 0.0
        for vdc, _ in buses:
            d = _case(res, cool, H, variant, vdc)[kind]["d"]
            on = d["region"] != "off"
            if on.any():
                top = max(top, float(np.sum([np.nan_to_num(d[k]) for k, _, _ in LAYERS], axis=0)[on].max()))
        ymax[kind] = max(top, 1.0)
    for r, kind in enumerate(("cont", "short")):
        for c, (vdc, _) in enumerate(buses):
            ax = axes[r][c]
            _style(ax, "speed [rpm]", "loss [W]")
            d = _case(res, cool, H, variant, vdc)[kind]["d"]
            on = d["region"] != "off"
            x = d["rpm"][on]
            base = np.zeros(on.sum())
            for key, label, col in LAYERS:
                v = np.nan_to_num(d[key][on])
                ax.fill_between(x, base, base + v, color=col, linewidth=0, label=label)
                if key == "P_cu_ac":
                    ax.plot(x, base, color=ORANGE, linewidth=1 * PT, linestyle=(0, (3, 2)))
                base = base + v
                for j in range(len(x)):
                    rows.append({"cooling": LABELS[cool], "vdc": vdc, "kind": kind, "rpm": float(x[j]), "layer": key, "W": float(v[j])})
            ax.set_xlim(0, TOP_RPM * 1.05)
            ax.set_ylim(0, ymax[kind] * 1.08)
            ax.set_title(f"{vdc:.0f} V, {'continuous' if kind == 'cont' else 'short-term (cold start)'}", fontsize=9, color=INK, loc="left")
            _thousands(ax)
            if r == 0 and c == 0:
                ax.legend(loc="upper left", frameon=False, fontsize=7.3, handlelength=1.3)
    return _save(fig, out, f"losses_{cool}", rows)


# ====================================================================================================== figure 3: the efficiency maps
def _quads(valid: np.ndarray) -> np.ndarray:
    """Triangles over a grid of points (the index grid of valid[i, j]): two for every cell whose four corners are valid."""
    nS, nF = valid.shape
    idx = np.arange(nS * nF).reshape(nS, nF)
    tri = []
    for i in range(nS - 1):
        for j in range(nF - 1):
            if valid[i, j] and valid[i + 1, j] and valid[i, j + 1] and valid[i + 1, j + 1]:
                tri.append([idx[i, j], idx[i + 1, j], idx[i, j + 1]])
                tri.append([idx[i + 1, j], idx[i + 1, j + 1], idx[i, j + 1]])
    return np.array(tri)


def _map_axes(ax, em, z, levels, cmap, lines, fmt_lines):
    X = np.broadcast_to(em["rpm"][:, None], em["T"].shape)
    Y = em["T"]
    valid = em["minloss"]["found"] & np.isfinite(z)
    tri = mtri.Triangulation(X.ravel(), Y.ravel(), _quads(valid))
    zz = np.where(valid, z, np.nanmin(z[valid])).ravel()
    cf = ax.tricontourf(tri, zz, levels=levels, cmap=cmap, extend="both")
    if lines is not None:
        cs = ax.tricontour(tri, zz, levels=lines, colors=INK, linewidths=0.6)
        ax.clabel(cs, fmt=fmt_lines, fontsize=7)
    return cf


def fig_efficiency(res: dict, out: Path, cool: str) -> Path | None:
    H, buses, variant = _head(res), _buses(res), _hv(res)
    keys = [f"{cool}|{H:g}|{variant}|{v:g}" for v, _ in buses]
    if not all(k in res.get("maps", {}) for k in keys):
        return None
    title = f"Efficiency map, minimum-loss control: {LABELS[cool]}, {_head_winding(res)}, {H:g} mm"
    gain = max(float(np.nanmax(res["maps"][k]["minloss"]["eta"] - res["maps"][k]["mtpa"]["eta"])) for k in keys)
    note = (f"Efficiency (shaft power over the machine's terminal power) at the current vector that minimises the total loss for each speed and torque; the map ends at the short-term "
            f"envelope (twice the thermal rating). Solid line: the continuous envelope; dashed: the short-term one. Minimum-loss control beats the plain MTPA control (the smallest current that holds "
            f"the voltage) by at most {100 * gain:.3f} percentage points anywhere on the map. Motor only: the inverter's loss is not in it. The copper is at {_INPUTS['t_cu_short']:g} C on this map "
            f"(the continuous operating points run at the network's temperature: a copper loss changes 0.4 % a kelvin). {provisional_text()} {not_modelled()}")
    fig = _figure(12.5, 5.8, title, note)
    axes = _grid(fig, 1, len(buses), cbar=0.012)
    levels = [0, 40, 60, 70, 75, 80, 85, 88, 90, 92, 94, 96, 100]
    rows = []
    for c, (vdc, _) in enumerate(buses):
        ax = axes[0][c]
        _style(ax, "speed [rpm]", "shaft torque [N.m]")
        em = res["maps"][keys[c]]
        cf = _map_axes(ax, em, 100 * em["minloss"]["eta"], levels, SEQ, [80, 90, 93], "%d %%")
        case = _case(res, cool, H, variant, vdc)
        d, ds = case["cont"]["d"], case["short"]["d"]
        ax.plot(d["rpm"], d["T"], color=INK, linewidth=2 * PT)
        ax.plot(ds["rpm"], ds["T"], color=INK, linewidth=1.4 * PT, linestyle="--")
        ax.set_xlim(0, TOP_RPM * 1.03)
        ax.set_ylim(0, float(np.nanmax(ds["T"])) * 1.05)
        ax.set_title(f"{vdc:.0f} V", fontsize=9, color=INK, loc="left")
        _thousands(ax)
        _vline(ax, TOP_RPM)
        _cbar(fig, ax, cf, "efficiency [%]", fmt="%d")
        for i in range(0, em["rpm"].size, 4):
            for j in range(0, em["T"].shape[1], 4):
                if em["minloss"]["found"][i, j]:
                    rows.append({"vdc": vdc, "rpm": float(em["rpm"][i]), "T_Nm": float(em["T"][i, j]), "eta_minloss_pct": float(100 * em["minloss"]["eta"][i, j]),
                                 "eta_mtpa_pct": float(100 * em["mtpa"]["eta"][i, j]), "loss_W": float(em["minloss"]["P_loss"][i, j])})
    return _save(fig, out, f"efficiency_{cool}", rows)


def fig_saving(res: dict, out: Path, cool: str) -> Path | None:
    H, buses, variant = _head(res), _buses(res), _hv(res)
    keys = [f"{cool}|{H:g}|{variant}|{v:g}" for v, _ in buses]
    if not all(k in res.get("maps", {}) for k in keys):
        return None
    title = f"What minimum-loss control saves over the MTPA control: {LABELS[cool]}, {_head_winding(res)}, {H:g} mm"
    mx = max(float(np.nanmax(np.maximum(res["maps"][k]["saving"], 0.0))) for k in keys)
    note = (f"The total loss of the MTPA control (the smallest current that makes the torque within the voltage limit) less that of the minimum-loss control, in watts, at every speed and torque of the map. "
            f"It is at most {mx:.2f} W here: well above the corner the voltage limit leaves one current vector for a torque, and below it the iron loss is small against the copper loss, so the best "
            f"vector is the MTPA one; only just above the corner is there a little to choose. {provisional_text()}")
    fig = _figure(12.5, 5.6, title, note)
    axes = _grid(fig, 1, len(buses), cbar=0.012)
    hi = max(mx, 1e-6)
    rows = []
    for c, (vdc, _) in enumerate(buses):
        ax = axes[0][c]
        _style(ax, "speed [rpm]", "shaft torque [N.m]")
        em = res["maps"][keys[c]]
        cf = _map_axes(ax, em, np.maximum(em["saving"], 0.0), np.linspace(0, hi, 13), SEQ, None, None)
        d = _case(res, cool, H, variant, vdc)["short"]["d"]
        ax.plot(d["rpm"], d["T"], color=INK, linewidth=1.4 * PT, linestyle="--")
        ax.set_xlim(0, TOP_RPM * 1.03)
        ax.set_ylim(0, float(np.nanmax(d["T"])) * 1.05)
        ax.set_title(f"{vdc:.0f} V", fontsize=9, color=INK, loc="left")
        _thousands(ax)
        _cbar(fig, ax, cf, "loss saved [W]", fmt="%.2f")
        for i in range(0, em["rpm"].size, 4):
            for j in range(0, em["T"].shape[1], 4):
                if np.isfinite(em["saving"][i, j]):
                    rows.append({"vdc": vdc, "rpm": float(em["rpm"][i]), "T_Nm": float(em["T"][i, j]), "saving_W": float(em["saving"][i, j])})
    return _save(fig, out, f"saving_{cool}", rows)


# ====================================================================================================== figure 4: the current plane
def fig_plane(res: dict, out: Path, cool: str, stack: float = 20.0) -> Path:
    variant = _hv(res)
    stat = res["static"][f"{cool}|{stack:g}|{variant}"]
    w = mm.Winding(stat["winding"]["turns"], stat["winding"]["paths"], stat["winding"]["strands"], stat["winding"]["wire_mm"], stat["winding"]["wire_over_mm"])
    kind = _INPUTS.get("copper_kind", "dc")
    m = mm.Machine(w, stack, copper=mm.FullCopper(arrangement=_INPUTS.get("arrangement", "cluster")) if kind == "full" else mm.DcCopper())
    v0 = _buses(res)[0][0]
    case = _case(res, cool, stack, variant, v0)
    bus = bus_of(v0)
    d = case["short"]["d"]
    I_pk = stat["i_peak"]
    sc = mm.Scan(m, i_max=I_pk, n_gamma=361, n_radius=241)
    sm = case["short"]["summary"]
    corner = sm["corner_rpm"] or 3000.0
    mtpv = d["rpm"][d["region"] == "MTPV"]
    title = (f"The current plane: what limits the torque, and the path the controller takes ({stack:g} mm, {_wl(stat['winding'])}, {v0:g} V, short-term current {I_pk:.1f} A)")
    note = ("The motoring half of the (i_d, i_q) plane; d is the magnet axis, so the magnets' flux is cancelled to the left. Thin gray lines: constant air-gap torque. Gray dashed: the MTPA "
            "path (the most torque per ampere at each current). The circle is the current limit. Each colored curve is the voltage limit |v| = V_max at one speed: the machine can work only "
            "inside it and the circle, and as the speed rises it shrinks around the point where the flux linkage is zero. Dots, one for every 500 rpm: blue at the MTPA point (below the corner "
            f"{corner:,.0f} rpm the circle is the only limit, and the point does not move), orange in field weakening (circle and curve both limit), aqua on MTPV"
            + (f" (from {mtpv.min():,.0f} rpm)" if len(mtpv) else "") + " (only the curve limits: more current would give less torque). " + provisional_text())
    fig = _figure(10.5, 8.8, title, note)
    ax = fig.add_axes([0.07, fig._bottom, 0.72, TOP - fig._bottom])
    _style(ax, "i_d [A]", "i_q [A]")
    ax.set_aspect("equal")
    T_em = 1.5 * mm.PP * (sc.psi_d * sc.i_q - sc.psi_q * sc.i_d)
    cs = ax.contour(sc.i_d, sc.i_q, T_em, levels=np.linspace(0.0, float(T_em.max()), 9)[1:-1], colors="#c3c2b7", linewidths=0.7)
    ax.clabel(cs, fmt="%.1f N.m", fontsize=6.5, colors=INK2)
    th = np.linspace(0.0, math.pi / 2, 300)
    ax.plot(-I_pk * np.sin(th), I_pk * np.cos(th), color=INK, linewidth=1.6 * PT)
    ramp = ["#9ec5f2", "#5d9be6", "#2a78d6", "#14478a", "#0e2f5c", "#0b0b0b"]
    speeds = sorted({round(s, 0) for s in (corner, 1.6 * corner, 2.5 * corner, 4.0 * corner, 6.0 * corner, TOP_RPM) if s <= TOP_RPM})
    handles = []
    for rpm, col in zip(speeds, ramp):
        c = sc.ctx(rpm, bus, 100.0, "air_gap")
        _, vv = sc.grid(c)
        ax.contour(sc.i_d, sc.i_q, np.sqrt(vv) - bus.v_max, levels=[0.0], colors=col, linewidths=1.4 * PT)
        handles.append(_line(col, f"|v| = V_max at {rpm:,.0f} rpm"))
    g, a, b, _ = m.maps.mtpa(np.linspace(0.02, 1.0, 100) * I_pk, turns=m.w.turns, paths=m.w.paths, stack_mm=stack)
    ax.plot(a, b, color=INK2, linewidth=1.2 * PT, linestyle="--")
    on = d["region"] != "off"
    reg_col = {"MTPA": (BLUE, "o", "constant torque: the MTPA point"), "FW": (ORANGE, "s", "field weakening"), "MTPV": (AQUA, "^", "MTPV")}
    rows = []
    for reg, (col, mk, lab) in reg_col.items():
        sel = (d["region"] == reg) & (np.round(d["rpm"] / 500.0) * 500.0 == d["rpm"])
        sel_all = d["region"] == reg
        if not sel_all.any():
            continue
        ax.plot(d["i_d"][sel_all], d["i_q"][sel_all], color=col, linewidth=1.0 * PT, alpha=0.8)
        ax.plot(d["i_d"][sel], d["i_q"][sel], linestyle="", marker=mk, markersize=6.5, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=1.5 * PT)
        handles.append(_dot(col, mk, 6.5, lab))
    for j in np.nonzero(on)[0]:
        rows.append({"rpm": float(d["rpm"][j]), "i_d_A": float(d["i_d"][j]), "i_q_A": float(d["i_q"][j]), "region": str(d["region"][j]), "T_Nm": float(d["T"][j])})
    handles.append(_line(INK2, "MTPA path", "--"))
    ax.set_xlim(-I_pk * 1.04, I_pk * 0.12)
    ax.set_ylim(-I_pk * 0.03, I_pk * 1.04)
    ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=7.8)
    return _save(fig, out, f"plane_{cool}_{stack:g}", rows)


# ====================================================================================================== figure 5: the stack table
def _table(fig, y_top: float, height: float, title: str, cols: list[str], cells: list[list[str]], fontsize: float = 8.6) -> float:
    """One table in the band [y_top - height, y_top] of the figure (fractions); returns the next band's top."""
    ax = fig.add_axes([0.012, y_top - height, 0.976, height - 0.034])
    ax.axis("off")
    ax.set_title(title, fontsize=9.5, color=INK, loc="left", pad=3)
    tb = ax.table(cellText=cells, colLabels=cols, bbox=[0.0, 0.0, 1.0, 1.0], cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(fontsize)
    for (r, c), cell in tb.get_celld().items():
        cell.set_edgecolor(GRID)
        cell.set_facecolor("#f1f6fd" if r == 0 else (SURFACE if r % 2 else "#f6f6f3"))
        cell.get_text().set_color(INK)
    return y_top - height


def fig_stack_table(res: dict, out: Path, cool: str) -> Path:
    stacks = res["inputs"]["stacks"]
    v0, v1 = res["inputs"]["vdc"][0], res["inputs"]["vdc"][-1]
    title = f"The stack table: 1 to 4 rows of magnets ({LABELS[cool]})"
    note = ("The winding study's winding: the winding the winding study chose for each stack (the smallest inverter current that still holds the continuous torque to 3,300 rpm on 48 V). Same "
            "winding: the headline winding on every stack; a longer stack links more flux, so its corner speed falls as 1 / length. "
            "Peak = twice the winding's own thermal rating at 3,000 rpm; torque and power are at the shaft; the overload time is from a COLD start at that current at 3,000 rpm (magnets to "
            f"{res['inputs']['magnet_limit_c']:g} C / winding to 155 C); the continuous figures come from the thermal network at the lower bus voltage (a t after a speed: the continuous envelope ends "
            "there because no current both holds the voltage and keeps the temperatures, not because the loads take all the torque). The fan is sized on each stack. Masses are active parts only. "
            + provisional_text())
    fig = _figure(15.5, 9.0, title, note)
    rows_csv = []
    cols_g = ["stack\n[mm]", "magnets\nto buy", "steel: stator\n+ rotor [kg]", "magnets\n[g]", "active\ntotal [kg]", "rotor inertia\n[g m^2]"]
    cells_g = []
    hv = _hv(res)
    for s in stacks:
        g = res["static"][f"{cool}|{s:g}|{hv}"]["geometry"]
        cells_g.append([f"{s:.0f}", f"{g['n_magnets']:.0f}", f"{g['m_stator_steel']:.2f} + {g['m_rotor_steel']:.2f}", f"{1e3 * g['m_magnets']:.0f}", f"{g['m_active']:.2f}", f"{1e3 * g['J_rotor']:.2f}"])
    fmt_s = lambda x: "-" if x is None else f"{x:,.0f}"
    blocks = [("Active parts (steel, magnets and inertia are the same for both windings; the copper and so the active total are of the winding study's windings)", cols_g, cells_g)]
    clipped = False
    for variant, name in (("opt", "The winding study's winding of each stack"), ("same", "The headline winding on every stack")):
        cols = ["stack\n[mm]", "winding\nN / a, strands x wire", "peak\ncurrent [A]", f"corner {v0:g} V\n[rpm]", f"corner {v1:g} V\n[rpm]", "peak torque\n[N.m]", f"peak power\n{v0:g} V [kW]", f"peak power\n{v1:g} V [kW]",
                "continuous T at\n3,000 rpm [N.m]", "continuous: torque\nends [rpm]", "overload from cold:\nmagnets / winding [s]", "copper\n[kg]"]
        cells = []
        for s in stacks:
            stt = res["static"][f"{cool}|{s:g}|{variant}"]
            c48, c60 = res["cases"][f"{cool}|{s:g}|{variant}|{v0:g}"], res["cases"][f"{cool}|{s:g}|{variant}|{v1:g}"]
            s48, s60, k48 = c48["short"]["summary"], c60["short"]["summary"], c48["cont"]["summary"]
            ov = c48["overload"]
            clipped |= bool(stt["peak_clipped"])
            when = lambda x: ">3600" if x is None else f"{x:,.0f}"
            ovs = f"{when(ov.get('magnet_s'))} / {when(ov.get('winding_s'))}" if ov.get("feasible") else "-"
            cells.append([f"{s:.0f}", _wl(stt["winding"]), f"{stt['i_peak']:.1f}{'*' if stt['peak_clipped'] else ''}", fmt_s(s48["corner_rpm"]), fmt_s(s60["corner_rpm"]), f"{s48['T_peak']:.2f}",
                          f"{s48['P_peak'] / 1e3:.2f}", f"{s60['P_peak'] / 1e3:.2f}", f"{k48['T_rated']:.2f}", fmt_s(k48["zero_rpm"]) + (" t" if k48.get("stop") == "limit" else ""), ovs,
                          f"{stt['geometry']['m_copper']:.2f}"])
            rows_csv.append({"cooling": LABELS[cool], "variant": variant, "stack_mm": s, "winding": _wl(stt["winding"]), "peak_I_A": stt["i_peak"], "corner_48_rpm": s48["corner_rpm"], "corner_60_rpm": s60["corner_rpm"],
                             "T_peak_Nm": s48["T_peak"], "P_peak_48_W": s48["P_peak"], "P_peak_60_W": s60["P_peak"], "short_zero_48_rpm": s48["zero_rpm"], "cont_T_3000_Nm": k48["T_rated"],
                             "cont_zero_48_rpm": k48["zero_rpm"], "cont_stop_48": k48.get("stop"), "overload_magnet_s": ov.get("magnet_s"), "overload_winding_s": ov.get("winding_s")})
        blocks.append((name, cols, cells))
    avail = TOP - fig._bottom - 0.03
    unit = avail / (sum(len(b[2]) + 1.6 for b in blocks) + 0.45 * len(blocks))
    y = TOP
    for name, cols, cells in blocks:
        y = _table(fig, y, unit * (len(cells) + 1.6 + 0.45), name, cols, cells) - 0.0
    if clipped:
        fig.text(0.012, fig._bottom + 0.004, "* the peak current is clipped at the edge of the flux maps (3 times the rated current of the 10-turn placeholder).", fontsize=7.5, color=INK2)
    return _save(fig, out, f"stack_table_{cool}", rows_csv)


# ====================================================================================================== figure 6: the front's windings, and the magnet grade
def _extra_key(res: dict, variant: str) -> str:
    return f"{res['inputs']['extra_fan']}|{_head(res):g}|{variant}"


def _front_members(res: dict) -> list[tuple[str, str, dict]]:
    """(key of the case, label, static record) of the headline winding and of the extras on the N45 grade, in order of the winding's peak current."""
    members = [(_extra_key(res, _hv(res)), _head_winding(res), res["static"][_extra_key(res, _hv(res))])]
    for e in res["inputs"]["extras"]:
        if e["magnet_limit_c"] == res["inputs"]["magnet_limit_c"]:
            key = _extra_key(res, e["key"])
            members.append((key, _wl(e["winding"]), res["static"][key]))
    return sorted(members, key=lambda t: t[2]["i_peak"])


def fig_front(res: dict, out: Path, vdc: float) -> Path | None:
    members = _front_members(res)
    if len(members) < 2:
        return None
    ramp = ["#9ec5f2", "#5d9be6", "#2a78d6", "#14478a", "#0b0b0b"][-len(members):] if len(members) <= 5 else SEQ(np.linspace(0.3, 1.0, len(members)))
    H = _head(res)
    title = f"The winding study's front on {vdc:g} V, {H:g} mm, {SHORT[res['inputs']['extra_fan']]}: what the inverter's current class buys"
    note = ("Each line is one winding of the front (the smallest inverter current that holds the continuous torque to 3,300 rpm, and the larger-current windings that trade it for a wider speed range); "
            "the darker the line, the larger the current. Left: continuous (the thermal network's limit at every speed); right: short-term, twice the winding's own continuous current at 3,000 rpm, "
            "from a COLD start, which is the current the inverter must supply. Top row shaft torque, bottom row shaft power. The torque at low speed is the same for all (the continuous torque is held "
            "the same); what the larger current buys is the corner speed and the torque and power above it. "
            f"{provisional_text()} {not_modelled()}")
    fig = _figure(12.0, 8.4, title, note)
    axes = _grid(fig, 2, 2)
    rows = []
    t_max = max(float(np.nanmax(res["cases"][f"{k}|{vdc:g}"]["short"]["d"]["T"])) for k, _, _ in members)
    p_max = max(float(np.nanmax(np.nan_to_num(res["cases"][f"{k}|{vdc:g}"]["short"]["d"]["P_out"]))) for k, _, _ in members) / 1e3
    for c, (kind, name) in enumerate((("cont", "continuous"), ("short", "short-term (cold start)"))):
        for r in range(2):
            ax = axes[r][c]
            _style(ax, "speed [rpm]", "shaft torque [N.m]" if r == 0 else "shaft power [kW]")
            for (key, label, stt), col in zip(members, ramp):
                case = res["cases"][f"{key}|{vdc:g}"]
                d = case[kind]["d"]
                y = d["T"] if r == 0 else np.nan_to_num(d["P_out"]) / 1e3
                ax.plot(d["rpm"], y, color=col, linewidth=2 * PT)
                if r == 0:
                    for j in range(0, len(y), 2):
                        rows.append({"winding": label, "vdc": vdc, "kind": kind, "rpm": float(d["rpm"][j]), "T_Nm": float(d["T"][j]), "P_kW": float(np.nan_to_num(d["P_out"][j]) / 1e3)})
                sm = case["short"]["summary"]
                if kind == "short" and sm["corner_rpm"] is not None:
                    yc = sm["corner_T"] if r == 0 else sm["corner_T"] * sm["corner_rpm"] * 2 * math.pi / 60 / 1e3
                    ax.plot([sm["corner_rpm"]], [yc], marker="o", markersize=6, markerfacecolor=col, markeredgecolor=SURFACE, markeredgewidth=2 * PT, linestyle="")
            _vline(ax, TOP_RPM, "rotor's safe speed" if (r == 0 and c == 0) else None)
            ax.set_xlim(0, 17200)
            ax.set_ylim(0, (t_max if r == 0 else p_max) * 1.06)
            if r == 0:
                ax.set_title(name, fontsize=9, color=INK, loc="left")
            _thousands(ax)
    _legend(axes[0][0], [_line(col, f"{label} ({stt['i_peak']:.0f} A peak)") for (key, label, stt), col in zip(members, ramp)], loc="upper center")
    return _save(fig, out, f"front_{vdc:g}V", rows)


def fig_grade(res: dict, out: Path) -> Path | None:
    """The 150 C magnet grade's winding against the headline winding on the N45 grade (the same fan, stack and buses)."""
    grade = [e for e in res["inputs"]["extras"] if e["magnet_limit_c"] != res["inputs"]["magnet_limit_c"]]
    if not grade:
        return None
    g = grade[0]
    buses = _buses(res)
    H = _head(res)
    base_key, g_key = _extra_key(res, _hv(res)), _extra_key(res, g["key"])
    title = f"The 150 C magnet grade ({_wl(g['winding'])}) against the headline winding ({_head_winding(res)}), {H:g} mm, {SHORT[res['inputs']['extra_fan']]}"
    note = (f"Blue: the headline winding with N45 magnets held at {res['inputs']['magnet_limit_c']:g} C; orange: the winding study's pick for a magnet grade that may run at {g['magnet_limit_c']:g} C "
            f"(the thermal limit is the only change: the flux maps stay N45's, and a grade that runs hotter has a lower remanence, so its real torque would be a few percent lower). Left: "
            f"{buses[0][0]:g} V; right: {buses[-1][0]:g} V. Solid: continuous; dashed: short-term (twice the winding's own continuous current, from a COLD start). Top row shaft torque, bottom row shaft power. "
            f"{provisional_text()} {not_modelled()}")
    fig = _figure(max(10.0, 5.4 * len(buses)), 8.4, title, note)
    axes = _grid(fig, 2, len(buses))
    rows = []
    t_max = max(float(np.nanmax(res["cases"][f"{k}|{v:g}"][kind]["d"]["T"])) for k in (base_key, g_key) for v, _ in buses for kind in ("cont", "short"))
    p_max = max(float(np.nanmax(np.nan_to_num(res["cases"][f"{k}|{v:g}"][kind]["d"]["P_out"]))) for k in (base_key, g_key) for v, _ in buses for kind in ("cont", "short")) / 1e3
    for c, (vdc, _) in enumerate(buses):
        for r in range(2):
            ax = axes[r][c]
            _style(ax, "speed [rpm]", "shaft torque [N.m]" if r == 0 else "shaft power [kW]")
            for key, col in ((base_key, BLUE), (g_key, ORANGE)):
                case = res["cases"][f"{key}|{vdc:g}"]
                for kind, ls in (("cont", "-"), ("short", "--")):
                    d = case[kind]["d"]
                    y = d["T"] if r == 0 else np.nan_to_num(d["P_out"]) / 1e3
                    ax.plot(d["rpm"], y, color=col, linestyle=ls, linewidth=2 * PT)
                    if r == 0:
                        for j in range(0, len(y), 2):
                            rows.append({"case": "N45 headline" if key == base_key else "150 C grade", "vdc": vdc, "kind": kind, "rpm": float(d["rpm"][j]), "T_Nm": float(d["T"][j]),
                                         "P_kW": float(np.nan_to_num(d["P_out"][j]) / 1e3)})
            _vline(ax, TOP_RPM, "rotor's safe speed" if (r == 0 and c == 0) else None)
            ax.set_xlim(0, 17200)
            ax.set_ylim(0, (t_max if r == 0 else p_max) * 1.06)
            if r == 0:
                ax.set_title(f"{vdc:g} V", fontsize=9, color=INK, loc="left")
            _thousands(ax)
    st_b, st_g = res["static"][base_key], res["static"][g_key]
    _legend(axes[0][-1], [_line(BLUE, f"{_head_winding(res)}, magnets {res['inputs']['magnet_limit_c']:g} C ({st_b['i_peak']:.0f} A peak)"),
                          _line(ORANGE, f"{_wl(g['winding'])}, magnets {g['magnet_limit_c']:g} C ({st_g['i_peak']:.0f} A peak)"), _line(INK2, "continuous"), _line(INK2, "short-term", "--")], loc="upper right")
    return _save(fig, out, "grade", rows)


def make_all(out: Path) -> list[Path]:
    res = load(out)
    fans, main_fan = _INPUTS["fan_keys"], _INPUTS["main_fan"]
    paths_ = [fig_envelopes(res, out, "fans"), fig_envelopes(res, out, "comparisons"), fig_speed_limit(res, out), fig_limits(res, out)]
    for cool in fans:
        paths_.append(fig_losses(res, out, cool))
        paths_.append(fig_efficiency(res, out, cool))
        paths_.append(fig_stack_table(res, out, cool))
    paths_.append(fig_saving(res, out, main_fan))
    paths_.append(fig_plane(res, out, main_fan, min(_INPUTS["stacks"])))
    for vdc, _ in _buses(res):
        paths_.append(fig_front(res, out, vdc))
    paths_.append(fig_grade(res, out))
    return [p for p in paths_ if p]


def main(argv=None) -> int:
    label = (argv or sys.argv[1:] or ["final"])[0]
    for p in make_all(OUT_ROOT / label):
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
