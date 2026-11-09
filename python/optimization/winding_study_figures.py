"""optimization.winding_study_figures: the pictures of the winding study (matplotlib PNGs drawn from what optimization/winding_study.py computed; no FEMM).

What it is, and why
    The winding study enumerates about twenty thousand windings, keeps the few that satisfy every constraint and reads off the exact front of three objectives. These figures are
    what that says: where the front lies, what the losses of the shortlisted windings are made of, what the strand arrangement costs, how the answer moves with the stack, the bus,
    the magnet grade and the cap on strands, and what the speed range of each winding is. Every figure has a title, a note under it that says what is drawn, units on its axes,
    and its numbers in a table beside it (the report.txt, the CSV files and results.json of the study's folder): the 'table view'.

The pictures (written by winding_study.py into out/optimization/winding_study/)
    fig_ws_front.png          the windings that pass the thermal limits and can hold the high-speed demand, as points: copper loss at the high-speed point against the short-term current
                              (left) and the corner speed on 48 V against it (right); the ones that satisfy the voltage rule, the exact front and the recommended winding marked.
    fig_ws_loss.png           what the copper loss of the shortlisted windings is made of at 3,000 rpm and at the high-speed point: the proximity, skin, PWM and circulating terms
                              (logarithmic axis) against the DC loss.
    fig_ws_arrangements.png   the recommended winding in each strand arrangement: the copper loss at the high-speed point split into the circulating part and the rest, and what the
                              arrangement does to the continuous rating.
    fig_ws_stacks.png         the exact fronts of the four stacks in the plane of the corner speed against the short-term current, and, at 60 mm, what the bus, the magnet grade and
                              the fan do to the continuous torque and the recommended winding.
    fig_ws_caps.png           the cap on strands in hand: the number of allowed windings and the mean cost in copper loss at the high-speed point, for each stack.
    fig_ws_speed.png          the torque-speed envelope on 48 V of the shortlisted windings (motor_map.py's engine): the short-term torque of each and the continuous torque of the
                              recommended winding.
    The look follows optimization/winding_figures.py and ac_slot_figures.py: surface #fcfcfb, ink #0b0b0b and #52514e, grid #e1e0d9, axes #c3c2b7; series blue #2a78d6, orange #eb6834,
    aqua #1baf7a; context gray #898781.

Run
    Not on its own: winding_study.py calls draw_all() after its stages (uv run python -m optimization.winding_study). A few seconds in all.

Self-check (the answers follow).
  1. Why is the copper loss of the front drawn on a logarithmic axis, and why does it fall so steeply to the right of the voltage rule?
  2. The AC terms are drawn on a logarithmic axis next to the DC loss. What would a linear axis show?
Answers.
  1. The loss at the high-speed point spans a factor of fifteen along the front (the windings of few turns need little field-weakening current, the ones of many need much): on a
     linear axis the cheap end would be a line on the floor. It falls to the right because fewer turns a path mean less voltage for the same speed, hence less field-weakening current, and
     the loss is the current squared.
  2. Nothing but the DC loss: the proximity, skin, PWM and circulating terms are a few tenths of a watt or less against some tens of watts. That is the finding (the first study's: the loss does
     not move along the front), and it is why the figure needs the logarithm to show it.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

from optimization.ac_slot_figures import AXIS, AQUA, BLUE, GRID, INK, INK2, MUTED, ORANGE, PT, SURFACE, _figure, _legend, _save, _style

STACK_COLOUR = {20.0: MUTED, 40.0: AQUA, 60.0: BLUE, 80.0: ORANGE}
TERMS = (("p_prox", "proximity", BLUE), ("p_skin", "skin", MUTED), ("p_pwm", "PWM ripple", ORANGE), ("p_circ", "circulating (compact bundle)", AQUA))


def _plain_ticks(ax, axis: str, ticks) -> None:
    """Ticks of a logarithmic axis written as plain numbers (1,000, not 1 x 10^3), and no minor ticks."""
    if axis == "x":
        ax.set_xticks(list(ticks))
        ax.set_xticklabels([f"{t:,.0f}" for t in ticks])
    else:
        ax.set_yticks(list(ticks))
        ax.set_yticklabels([f"{t:,.0f}" for t in ticks])
    ax.minorticks_off()


def _tag(d: dict) -> str:
    """Short label of a winding: N, a, strands x wire."""
    return f"N {d['N']}  a {d['a']}\n{d['n_s']} x {d['d_mm']:g}"


def _marks(ax, x, y, run_masks: dict, rec_xy, front_xy, front_labels, label_dy: float) -> None:
    """The points of the front figure: the grey cloud, the windings that satisfy every constraint, the exact front (joined) and the recommended winding."""
    ax.scatter(x[run_masks["thermal"]], y[run_masks["thermal"]], s=7, color=MUTED, alpha=0.4, linewidths=0, label="thermal limits kept, high-speed demand held", zorder=2)
    ax.scatter(x[run_masks["allowed"]], y[run_masks["allowed"]], s=44, color=BLUE, edgecolors=SURFACE, linewidths=0.6, label="and corner speed on 48 V >= 3,300 rpm", zorder=3)
    fx, fy = front_xy
    ax.plot(fx, fy, color=ORANGE, lw=1.6 * PT, marker="D", ms=4.6, mec=SURFACE, mew=0.8, label="exact front (3 objectives)", zorder=4)
    for px, py, text in front_labels:
        ax.annotate(text, (px, py), textcoords="offset points", xytext=(6, label_dy), fontsize=6.6, color=INK2, ha="left", va="bottom")
    if rec_xy is not None:
        ax.scatter([rec_xy[0]], [rec_xy[1]], s=170, facecolors="none", edgecolors=INK, linewidths=1.3, label="recommended winding", zorder=5)


def fig_front(h: dict, path: Path) -> None:
    run, case = h["run"], h["case"]
    c, o = run.cand, run.obj
    ok_thermal = o["thermal_ok"] & (c["n_s"] <= case.n_s_cap) & o["hs_feasible"]
    ok_all = run.allowed(True)
    front = run.front(True)
    front = front[np.argsort(o["I_st"][front], kind="stable")]
    masks = {"thermal": ok_thermal & ~ok_all, "allowed": ok_all}
    groups: dict = {}
    for j in front:                                                          # windings of one current (one N / a) share a label: they lie on top of each other
        groups.setdefault(round(float(o["I_st"][j]), 3), []).append(j)
    labels = [(float(o["I_st"][g[0]]), float(o["P_hs"][g[0]]), f"{c['N'][g[0]]}/{c['a'][g[0]]}: " + " or ".join(f"{c['n_s'][j]} x {c['d_mm'][j]:g}" for j in g)) for g in groups.values()]
    j0 = h["rec_index"]
    fig = _figure(11.4, 5.6, f"The windings that can be chosen, and the exact front: {case.stack_mm:g} mm stack, {case.v_dc:g} V",
                  f"Every point is one winding (turns N a coil side, parallel paths a, strands in hand n_s of d mm) that keeps the hot spot under {case.wind_limit_c:g} C and the magnets under "
                  f"{case.magnet_limit_c:g} C at the common continuous rating and can hold the high-speed demand within its short-term current (at most {case.n_s_cap} strands, {case.arrangement} bundle). "
                  "Short-term current = twice the continuous one; copper at 100 C. Labels: N / a: strands x wire [mm]. Table view: report.txt, front_60mm.csv.")
    for k, (ykey, ylabel, ylim) in enumerate((("P_hs", "copper loss at 10,000 rpm, 1.24 N.m, 48 V [W]", None), ("corner_bus", f"corner speed at the continuous torque on {case.v_dc:g} V [rpm]", None))):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, "short-term phase current I_st [A peak]  (sizes the inverter)", ylabel)
        ax.set_position([0.07 + 0.50 * k, 0.24, 0.41, 0.64])
        x, y = o["I_st"], o[ykey]
        _marks(ax, x, y, masks, (x[j0], y[j0]) if j0 is not None else None, (x[front], y[front]), labels if k == 0 else [], 3.0)
        ax.set_xscale("log")
        ax.set_yscale("log")
        _plain_ticks(ax, "x", (40, 60, 80, 100, 150, 200, 300, 500))
        ax.set_xlim(35, 600)
        if k == 0:
            _plain_ticks(ax, "y", (5, 10, 20, 50, 100, 200))
            ax.set_ylim(4.5, 260)
            ax.axvline(x[ok_all].min(), color=ORANGE, lw=1 * PT, ls="--", zorder=1)
            ax.text(x[ok_all].min() * 1.03, 4.9, "voltage rule: no winding to\nthe left of the line satisfies it", fontsize=7, color=INK2, va="bottom")
            _legend(ax, "upper right")
        else:
            _plain_ticks(ax, "y", (1000, 2000, 3000, 5000, 10000, 20000))
            ax.set_ylim(1200, 15000)
            ax.axhline(1.1 * 3000.0, color=ORANGE, lw=1 * PT, ls="--", zorder=1)
            ax.text(570, 1.1 * 3000.0 * 0.96, "3,300 rpm: 3,000 rpm with 10 %\nof the voltage to spare", fontsize=7, color=INK2, va="top", ha="right")
    _save(fig, path)


def fig_loss(h: dict, results: dict, path: Path) -> None:
    rows = list(results["headline"]["shortlist"])
    if results["headline"]["recommended"] is not None:
        rows = [results["headline"]["recommended"]] + [r for r in rows if (r["N"], r["a"], r["n_s"], r["d_mm"]) != (results["headline"]["recommended"]["N"], results["headline"]["recommended"]["a"],
                                                                                      results["headline"]["recommended"]["n_s"], results["headline"]["recommended"]["d_mm"])]
    rows.sort(key=lambda d: -d["N"] / d["a"])
    fig = _figure(11.4, 5.6, "What the copper loss of the shortlisted windings is made of",
                  "Bars: the AC terms (watts, logarithmic axis); the black diamond and its number: the DC loss of the slots and the end turns, which the AC terms are small against (copper at 100 C; 3,000 rpm "
                  "at each winding's own continuous rating; 10,000 rpm at 1.24 N.m on 48 V, the smallest current that holds the voltage). The recommended winding is the first. Table view: report.txt, shortlist_60mm.csv.")
    for k, (pre, title) in enumerate((("cont", "3,000 rpm, continuous rating"), ("hs", "10,000 rpm, 1.24 N.m (a short-term demand)"))):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, "", "loss [W]")
        ax.set_position([0.07 + 0.49 * k, 0.30, 0.42, 0.58])
        use = [d for d in rows if pre == "cont" or d["hs_feasible"]]
        x = np.arange(len(use))
        for m, (key, name, colour) in enumerate(TERMS):
            vals = np.array([max(d[f"{pre}_{key}"], 1e-7) for d in use])
            ax.bar(x + (m - 1.5) * 0.19, vals, width=0.18, color=colour, label=name, zorder=3)
        dc = np.array([d[f"{pre}_p_dc_slot"] + d[f"{pre}_p_dc_end"] for d in use])
        ax.scatter(x, dc, marker="D", s=34, color=INK, zorder=4, label="DC loss (slots + end turns)")
        for xi, v, d in zip(x, dc, use):
            ac = d[f"{pre}_p_total"] - v
            ax.annotate(f"{v:.1f} W\nAC {100 * ac / d[f'{pre}_p_total']:.1f} %", (xi, v), textcoords="offset points", xytext=(0, 7), fontsize=6.8, color=INK2, ha="center", va="bottom")
        ax.set_yscale("log")
        ax.set_ylim(1e-5, 4e2)
        ax.set_xticks(x)
        ax.set_xticklabels([_tag(d) for d in use], fontsize=7.2, color=INK2)
        ax.set_title(title, fontsize=9.2, color=INK2, loc="left")
        if k == 0:
            handles, names = ax.get_legend_handles_labels()
    leg = fig.legend(handles, names, loc="lower center", bbox_to_anchor=(0.5, 0.085), ncol=5, frameon=False, fontsize=8.2)
    for text in leg.get_texts():
        text.set_color(INK2)
    _save(fig, path)


def fig_arrangements(h: dict, results: dict, path: Path) -> None:
    rows = results["rules"]["arrangements"]
    rec = results["headline"]["recommended"]
    names = [r["arrangement"] for r in rows]
    fig = _figure(10.6, 5.2, f"What the strand arrangement costs the recommended winding (N {rec['N']}, a {rec['a']}, {rec['n_s']} x {rec['d_mm']:g} mm)",
                  "Left: the copper loss at 10,000 rpm, 1.24 N.m, 48 V (copper 100 C), the circulating-current part on top. Right: the change of the winding's own continuous rating (3,000 rpm, the full copper "
                  "loss) against the compact bundle's, which is the cost in continuous torque. Transposed and row are options; random and line are the layouts to avoid. Table view: report.txt.")
    ax = fig.add_subplot(1, 2, 1)
    _style(ax, "", "copper loss at 10,000 rpm [W]")
    ax.set_position([0.08, 0.24, 0.40, 0.62])
    x = np.arange(len(rows))
    rest = np.array([r["hs_p_total"] - r["hs_p_circ"] for r in rows])
    circ = np.array([r["hs_p_circ"] for r in rows])
    ax.bar(x, rest, width=0.6, color=BLUE, label="DC, proximity, skin and PWM", zorder=3)
    ax.bar(x, circ, width=0.6, bottom=rest, color=ORANGE, label="circulating currents between strands", zorder=3)
    for xi, r in zip(x, rows):
        ax.annotate(f"{r['hs_p_total']:.1f} W\n({r['d_hs_loss_pct']:+.1f} %)", (xi, r["hs_p_total"]), textcoords="offset points", xytext=(0, 4), fontsize=7.2, color=INK2, ha="center", va="bottom")
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8, color=INK2)
    ax.set_ylim(0, 1.18 * max(r["hs_p_total"] for r in rows))
    _legend(ax, "upper left")
    ax = fig.add_subplot(1, 2, 2)
    _style(ax, "", "change of the continuous rating against the compact bundle [%]")
    ax.set_position([0.57, 0.24, 0.40, 0.62])
    dat = np.array([r["d_AT_pct"] for r in rows])
    cols = [BLUE if n in ("cluster", "row", "transposed") else ORANGE for n in names]
    ax.bar(x, dat, width=0.6, color=cols, zorder=3)
    for xi, r in zip(x, rows):
        v = r["d_AT_pct"]
        ax.annotate(f"{r['AT_cont']:.2f} A-turns\n({v:+.2f} %)", (xi, v), textcoords="offset points", xytext=(0, 5 if v >= 0 else -5), fontsize=7.2, color=INK2, ha="center",
                    va="bottom" if v >= 0 else "top")
    ax.axhline(0, color=AXIS, lw=1 * PT, zorder=2)
    ax.set_ylim(min(dat.min() * 1.9, -0.5), max(dat.max() + 0.5, 0.4))
    ax.set_xticks(x)
    ax.set_xticklabels(names, fontsize=8, color=INK2)
    _save(fig, path)


def fig_stacks(h: dict, results: dict, path: Path) -> None:
    var = results.get("variants")
    if not var:
        return
    fig = _figure(11.6, 5.6, "The answer against the stack, the bus, the magnet grade and the fan",
                  "Left: the exact front of each stack (thermal limits at that stack's common rating, voltage rule, at most 8 strands) in the plane of the corner speed at the continuous torque on 48 V "
                  "against the short-term current; the star is the recommended winding of the stack. Right (60 mm): the continuous torque each case holds its windings to (the common rating) and the recommended winding with its continuous and short-term current at it. "
                  "Table view: report.txt, shortlist_<stack>mm.csv.")
    ax = fig.add_subplot(1, 2, 1)
    _style(ax, "short-term phase current I_st [A peak]", "corner speed at the continuous torque, 48 V [rpm]")
    ax.set_position([0.07, 0.24, 0.40, 0.62])
    for stack, pts in sorted(var["fronts"].items()):
        pts = sorted(pts, key=lambda p: p["I_st"])
        ax.plot([p["I_st"] for p in pts], [p["corner_bus"] for p in pts], color=STACK_COLOUR[stack], lw=1.6 * PT, marker="o", ms=4.2, mec=SURFACE, mew=0.7, label=f"{stack:g} mm", zorder=3)
        row = var["stacks"][stack]
        if "I_st_common" in row:
            ax.scatter([row["I_st_common"]], [row["corner48_common"]], marker="*", s=190, color=STACK_COLOUR[stack], edgecolors=INK, linewidths=0.8, zorder=5)
    ax.axhline(3300.0, color=ORANGE, lw=1 * PT, ls="--", zorder=1)
    ax.text(440, 3300 * 0.97, "3,300 rpm: the voltage rule", fontsize=7.2, color=INK2, va="top", ha="right")
    ax.set_xscale("log")
    ax.set_yscale("log")
    _plain_ticks(ax, "x", (20, 30, 40, 60, 100, 200, 400))
    _plain_ticks(ax, "y", (3000, 5000, 10000, 20000, 40000))
    ax.set_xlim(20, 450)
    ax.set_ylim(2500, 4e4)
    _legend(ax, "upper left", title="stack")
    ax = fig.add_subplot(1, 2, 2)
    _style(ax, "continuous torque every winding of the case is held to [N.m]", "")
    ax.set_position([0.62, 0.24, 0.35, 0.62])
    names = list(var["variants"].keys())
    y = np.arange(len(names))[::-1]
    for yi, name in zip(y, names):
        r = var["variants"][name]
        ax.barh(yi, r["T_req"], height=0.56, color=BLUE if name.startswith("48") else (AQUA if "V bus" in name else ORANGE), zorder=3)
        ax.annotate(f"N {r['N']}, a {r['a']}, {r['n_s']} x {r['d_mm']:g} mm\n{r['I_cont_common']:.0f} A continuous, {r['I_st_common']:.0f} A short-term", (r["T_req"], yi),
                    textcoords="offset points", xytext=(6, 0), fontsize=7, color=INK2, ha="left", va="center")
    ax.set_yticks(y)
    ax.set_yticklabels(names, fontsize=7.6, color=INK2)
    ax.set_xlim(0, 17)
    ax.set_xticks([0, 2, 4, 6, 8, 10, 12])
    _save(fig, path)


def fig_caps(h: dict, results: dict, path: Path) -> None:
    var = results.get("variants")
    if not var:
        return
    fig = _figure(10.6, 5.2, "The cap on strands in hand",
                  "Left: the number of windings that satisfy every constraint when at most `cap` strands are wound in hand. Right: how much more copper loss at the high-speed point the best winding "
                  "at each N / a has under the cap than without one (mean over the N / a values that exist at the cap). The cap of 8 is marked. Table view: report.txt.")
    for k, (key, ylabel) in enumerate((("allowed", "windings that satisfy every constraint"), ("mean_excess_pct", "mean extra copper loss at 10,000 rpm over no cap [%]"))):
        ax = fig.add_subplot(1, 2, k + 1)
        _style(ax, "cap on strands in hand", ylabel)
        ax.set_position([0.08 + 0.49 * k, 0.24, 0.40, 0.62])
        for stack, rows in sorted(var["cap"].items()):
            ax.plot([r["cap"] for r in rows], [r[key] if key in r else np.nan for r in rows], color=STACK_COLOUR[stack], lw=1.6 * PT, marker="o", ms=4.5, mec=SURFACE, mew=0.7, label=f"{stack:g} mm")
        ax.axvline(8, color=ORANGE, lw=1 * PT, ls="--", zorder=1)
        ax.set_xscale("log")
        ax.set_xticks([4, 6, 8, 12, 30])
        ax.set_xticklabels(["4", "6", "8", "12", "30"])
        ax.minorticks_off()
        if k == 0:
            ax.set_yscale("log")
            _legend(ax, "upper left", title="stack")
        else:
            ax.set_ylim(bottom=0)
    _save(fig, path)


def fig_speed(h: dict, results: dict, path: Path) -> None:
    env = results["headline"].get("envelopes")
    if not env:
        return
    palette = [INK, BLUE, ORANGE, AQUA, MUTED, INK2]
    fig = _figure(10.8, 5.4, "Torque against speed on 48 V: the shortlisted windings",
                  "d-q torque (air-gap basis) at the largest torque inside the voltage limit V_dc / sqrt 3 and the winding's current limit (motor_map.py's solver on the same flux maps, copper at 100 C): the "
                  "short-term current (solid) and, for the recommended winding, the continuous current (dashed). Each winding is held to the same continuous torque at 3,000 rpm; the turns decide how far it reaches. "
                  "Table view: report.txt.")
    ax = fig.add_subplot(1, 1, 1)
    _style(ax, "speed [rpm]", "torque [N.m]")
    ax.set_position([0.07, 0.22, 0.90, 0.66])
    for e, colour in zip(env, palette):
        ax.plot(e["rpm"], e["T_st"], color=colour, lw=1.8 * PT, marker="o", ms=3.6, mec=SURFACE, mew=0.6, label=f"{e['name']}: N {e['N']}, a {e['a']}, {e['n_s']} x {e['d_mm']:g} mm, {e['I_st']:.0f} A short-term")
        if e["name"] == "recommended":
            ax.plot(e["rpm"], e["T_cont"], color=colour, lw=1.6 * PT, ls="--", label=f"recommended, continuous ({e['I_cont']:.0f} A)")
    ax.axvline(3000, color=AXIS, lw=1 * PT, zorder=0)
    ax.text(3000 + 120, 0.2, "3,000 rpm", fontsize=7.2, color=INK2, va="bottom")
    ax.set_xlim(0, 16500)
    ax.set_ylim(0, 12.4)
    _legend(ax, "upper center", ncol=2)
    _save(fig, path)


def draw_all(h: dict, results: dict, out: Path) -> list[Path]:
    """Every figure of the study into `out`; returns the paths written."""
    written = []
    for name, fn, args in (("fig_ws_front.png", fig_front, (h,)), ("fig_ws_loss.png", fig_loss, (h, results)), ("fig_ws_arrangements.png", fig_arrangements, (h, results)),
                           ("fig_ws_stacks.png", fig_stacks, (h, results)), ("fig_ws_caps.png", fig_caps, (h, results)), ("fig_ws_speed.png", fig_speed, (h, results))):
        path = Path(out) / name
        try:
            fn(*args, path)
        except KeyError as error:
            print(f"({name}: not drawn, {error} is missing from the results)")
            continue
        if path.exists():
            written.append(path)
            print(f"  figure {path}")
    return written
