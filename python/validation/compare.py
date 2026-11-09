"""compare.py: FEMM 4.2 against FEMM 5.0, number by number, with a verdict for every number.

What this proves. The motor-design scripts ran for weeks on FEMM 4.2. They now run on FEMM 5.0, which is a different
program (its mesher, `tangle`, replaces `Triangle`, and its solver replaces FEMM 4.2's). Having the same model (see
same_model.py) is not enough: the same model can give different numbers. This script reads every number the old run
left behind and every number the new run made, puts them side by side, and says for each one whether it moved by more
than it should.

"The report" below is the project note "4.2 vs 5.0": it holds the
numbers this script prints, explains every flag (section 2.12) and says what is still open (section 7).

Run (from the Python to FEMM folder, after `validation/run_all.py` has run every script on FEMM 5.0):
    uv run python validation/compare.py femm42-out [out] [--csv FILE] [--detail FILE] [--audit]
femm42-out (OLD_OUT) is FEMM 4.2's run: a copy of the `out` folder that FEMM 4.2 left in D:/dev/scripts/automation/proc/out,
checked file by file (femm42-out.sha256 lists the SHA-256 of all 634 files). [out] (NEW_OUT) is FEMM 5.0's (default: out).
Nothing in femm42-out is written to. It prints one table per study and writes
    validation/femm5_vs_femm42.csv          one row per compared number: study, quantity, femm42, femm50,
                                            delta_pct, verdict (and kind, limit_pct, plain_pct, n, ... after them)
    NEW_OUT/validation/femm5_vs_femm42_detail.csv   one row per state, solve or sample (thousands; for a spreadsheet)
With --audit it also runs the comparison five more times (twice with the safeguard below switched off, three times with it set a
little lower or higher) and lists the rows whose verdict depends on it (see audit(); about 20 seconds).

How to read a row.
  femm42, femm50   the two values. FEMM 4.2's come from the old run's files (csv, logs, .ans), or, for what the old run
                   only printed, from the numbers the README recorded on 2026-10-28.
  delta_pct        (femm50 - femm42) / femm42 * 100, with one safeguard (below). For a group of numbers that belong
                   together (the 2,040 states of a sweep, the 120 no-load samples) the row is the median, the 95th
                   percentile (p95) and the maximum of |delta|, and the group's verdict sits on the max row.
  verdict          identical  equal to the last digit that is stored.
                   ok         within the tolerance for its kind (the table below).
                   FLAG       beyond it. A flag is not a failure. It is a number that has to be explained, with a
                              check, in the report (section 2.12). Nothing is accepted silently.
                   info       there is no tolerance for this number (a node count, a time); it is only shown.
                   DIFFERS / MISSING / FAIL   something that should have been equal is not, a file is absent, or
                              a check that must hold on FEMM 5.0 itself does not.

The tolerances (from the plan; the same ones for every study).
  kind                          fine mesh   coarse mesh
  a mean torque, a flux linkage    1 %          2 %          (a state's torque and every quantity worked out from them:
                                                              EMF, inductance, saliency, air-gap B, the MTPA angle)
  a ripple, a cogging torque,     10 %         10 %          (judged on its own value, however small it is)
  a harmonic, a tail statistic
The mesh is "fine" where the study used FEMM's smart mesh (the barrier study, the slice of validate_full) and "coarse"
where it did not (the flux explorer, the optimizer, check_femm without --fine).

The one safeguard. A state's flux linkage or torque can be nearly zero (a flux linkage that the current has almost
cancelled, a torque at the angle where it should vanish). A percentage of almost nothing says nothing, so a FEMM 4.2
value smaller than a tenth of the largest of its kind in the same study (FLOOR) is judged against that tenth instead. The
plain percentage, which has no safeguard, is kept in the column plain_pct for a single number. For a group (a sweep's
states, a waveform's samples) that column is empty and the column n_plain_flagged counts the members that are beyond the
tolerance by the plain rule (the members themselves are in the detail table), so nothing hides behind the floor.
A second rule works the same way: the samples of a waveform that crosses zero are judged as a share of the waveform's
own peak (or of the cogging curve's peak-to-peak), because the percentage of one sample near a zero crossing says nothing.
--audit lists the rows that pass only through the floor or only through that rule, with how close each one is to its limit.

Exit status. 0 if every check that must hold on FEMM 5.0 itself holds and nothing to compare is missing. 1 if one of those fails
(FAIL) or an input is missing (MISSING). 2 if a folder is not there. FLAG and DIFFERS rows do not make the run fail: they are
printed at the end, and the report (section 2.12) explains every one of them.

What to expect (FEMM 5.0 of 2026-11-04 against the FEMM 4.2 run of 2026-10-28; 350 verdicts, which are 462 data lines in the
CSV, because a group of numbers writes three lines, its median, its p95 and its maximum): 51 identical, 206 ok, 73 info, 19
FLAG and one DIFFERS (a picture, check_vs_fem.png: 850 of 1.9 million pixels, along the stator's slot-bottom arcs). The flags
are the ones the report explains: flux linkage where the d current has cancelled most of the magnet flux, a few single torques
and flux linkages of the coarse flux-explorer sweep, three small torque harmonics (the 24th and 30th), the air-gap flux density
at the edges of its waveform, the full model's torque, one optimizer angle on a flat maximum, and the tail statistics of the materials.
What --audit prints on the same data: 13 rows that pass only through the floor (three of them within 4 % of their limit: the with-band
mean torque of the 52 (I, gamma) points, 0.997 % of 1 %; the no-band torque at 312 single positions, 0.964 %; the 120-degree psi_d,
1.992 % of 2 %), 3 FLAG rows whose members the floor hides (2 of 312 become 18, 2 of 2,040 become 384, 11 of 2,040 become 217),
17 rows that pass only through the waveform rule, and 22, 19 and 16 FLAG rows with the floor at 0.09, 0.12 and 0.15 (19 at 0.10).

Self-check (the answers are in the report, section "Self-check", questions 1, 2 and 4).
  1. A flux linkage changes from 0.023 mWb to 0.031 mWb, and the machine's flux linkage is 12 mWb. What are the plain
     percentage and the judged one, and which one is the fair description of FEMM 5.0?
  2. Why is a harmonic that grows from 0.00066 to 0.00104 N.m flagged, although it is under a tenth of a percent of
     the mean torque?
  3. Why does the optimizer's generation 0 compare row by row, and generation 3 not?
"""
from __future__ import annotations

import argparse
import contextlib
import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
P = HERE.parent
for _folder in (str(P), str(HERE)):                       # femfile.py is in P, compare_materials.py is next to this file
    if _folder not in sys.path:
        sys.path.insert(0, _folder)
import compare_materials  # noqa: E402
import femfile  # noqa: E402

# ----------------------------------------------------------------------------------------------- 1. the rules
FLOOR = 0.10                  # a FEMM 4.2 value below this fraction of the biggest of its kind is judged against that fraction
WAVE_RULE = True              # the samples of a waveform are judged as a share of the waveform's own peak (Report.group, wave=True)
LIMITS = {                    # percent: (fine mesh, coarse mesh)
    "mean_torque": (1.0, 2.0), "flux": (1.0, 2.0), "angle": (1.0, 2.0),
    "ripple": (10.0, 10.0), "cogging": (10.0, 10.0), "harmonic": (10.0, 10.0), "tail": (10.0, 10.0),
}
FLOORED = {"mean_torque", "flux", "angle"}                # kinds whose small values get the safeguard
KINDS = set(LIMITS) | {"exact", "info"}                   # exact: must be equal. info: no tolerance, only shown


def limit_for(kind: str, mesh: str) -> float | None:
    """The tolerance in percent for a kind of number on a mesh ("fine" or "coarse"); None where there is none."""
    if kind not in LIMITS:
        return None
    return LIMITS[kind][0 if mesh == "fine" else 1]


def plain_percent(old: float, new: float) -> float | None:
    """(new - old) / |old| * 100, or None where old is zero."""
    return None if old == 0 else 100.0 * (new - old) / abs(old)


def judged_percent(old: float, new: float, kind: str, scale: float | None = None) -> float | None:
    """The percentage the verdict is based on: the plain one, except that a small old value is judged against FLOOR * scale."""
    if kind in FLOORED and scale:
        return 100.0 * (new - old) / max(abs(old), FLOOR * scale)
    return plain_percent(old, new)


def verdict_for(kind: str, mesh: str, delta: float | None, old=None, new=None) -> str:
    """identical / ok / FLAG for a number that has a tolerance, info for one that has none, DIFFERS for exact ones."""
    if kind == "exact":
        return "identical" if old == new else "DIFFERS"
    limit = limit_for(kind, mesh)
    if limit is None or delta is None:
        return "info"
    if old is not None and old == new:
        return "identical"
    return "ok" if abs(delta) <= limit + 1e-12 else "FLAG"


def harmonic_amplitudes(y) -> np.ndarray:
    """Amplitudes of harmonics 0..N/2 of one period sampled evenly (the last point is not repeated)."""
    c = np.fft.rfft(np.asarray(y, float)) / len(y)
    c[1:] *= 2
    return np.abs(c)


# ----------------------------------------------------------------------------------------------- 2. rows
@dataclass
class Row:
    """One compared number (or one group of them). Written to femm5_vs_femm42.csv."""
    study: str
    quantity: str
    femm42: float | str | None
    femm50: float | str | None
    delta_pct: float | None
    verdict: str
    kind: str = ""
    limit_pct: float | None = None
    plain_pct: float | None = None
    n: int = 1                                  # members of a group
    n_flagged: int = 0                          # members beyond the tolerance (judged percentage)
    n_plain_flagged: int = 0                    # members beyond the tolerance by the plain percentage (no floor)
    note: str = ""
    stats: tuple | None = None                  # (median, p95, max) of |delta| for a group


class Report:
    """Collects rows and detail lines. Every comparison in this file goes through value(), group() or check()."""

    def __init__(self):
        self.rows: list[Row] = []
        self.detail: list[dict] = []

    # one number against one number
    def value(self, study, quantity, old, new, kind, mesh="fine", scale=None, note=""):
        if kind == "exact" or kind == "info" and (isinstance(old, str) or isinstance(new, str)):
            self.rows.append(Row(study, quantity, old, new, None, verdict_for("exact", mesh, None, old, new) if kind == "exact"
                                 else "info", kind, None, None, note=note))
            return self.rows[-1]
        delta = judged_percent(old, new, kind, scale)
        plain = plain_percent(old, new)
        verdict = verdict_for(kind, mesh, delta, old, new)
        n_flag = int(verdict == "FLAG")
        n_plain = int(plain is not None and limit_for(kind, mesh) is not None and abs(plain) > limit_for(kind, mesh))
        self.rows.append(Row(study, quantity, old, new, delta, verdict, kind, limit_for(kind, mesh), plain,
                             n_flagged=n_flag, n_plain_flagged=n_plain, note=note))
        return self.rows[-1]

    # many numbers against as many (a sweep's states, a waveform's samples)
    def group(self, study, quantity, old, new, kind, mesh="coarse", scale=None, wave=False, keys=None, note=""):
        """Compare arrays element by element. The row carries median, p95 and max of |delta|.

        scale: the biggest value of this kind in the study (the safeguard's base). wave=True judges every sample against
        the waveform's own scale (its peak, or its peak-to-peak for a cogging curve): for a signal that crosses zero
        the percentage of one sample's value says nothing, the percentage of the waveform's size does.
        """
        old, new = np.asarray(old, float), np.asarray(new, float)
        if wave and WAVE_RULE:
            base = np.full_like(old, float(scale))
        elif kind in FLOORED and scale:
            base = np.maximum(np.abs(old), FLOOR * scale)
        else:
            base = np.abs(old)
        with np.errstate(divide="ignore", invalid="ignore"):
            d = 100.0 * (new - old) / base
            plain = 100.0 * (new - old) / np.abs(old)
        limit = limit_for(kind, mesh)
        flagged = np.abs(d) > limit + 1e-12 if limit is not None else np.zeros(len(d), bool)
        plain_flagged = (np.abs(plain) > limit + 1e-12) & np.isfinite(plain) if limit is not None else np.zeros(len(d), bool)
        a = np.abs(d[np.isfinite(d)])
        stats = (float(np.median(a)), float(np.percentile(a, 95)), float(a.max())) if len(a) else (math.nan,) * 3
        verdict = "info" if limit is None else ("identical" if not np.any(new != old) else
                                                ("FLAG" if flagged.any() else "ok"))
        self.rows.append(Row(study, quantity, None, None, stats[2], verdict, kind, limit, None, len(d), int(flagged.sum()),
                             int(plain_flagged.sum()), note, stats))
        for i in range(len(d)):
            self.detail.append({"study": study, "quantity": quantity, "member": keys[i] if keys is not None else i,
                                "femm42": old[i], "femm50": new[i], "delta_pct": d[i],
                                "plain_pct": plain[i] if np.isfinite(plain[i]) else "", "beyond": bool(flagged[i])})
        return self.rows[-1]

    # a yes/no statement about one side only (FEMM 5.0 alone, or the old numbers alone)
    def check(self, study, quantity, value, limit, note="", old=None, higher_is_worse=True, verdict=None):
        """A check that must hold: ok when value <= limit (or >= limit), FAIL when not."""
        holds = value <= limit if higher_is_worse else value >= limit
        self.rows.append(Row(study, quantity, old, value, None, verdict or ("ok" if holds else "FAIL"), "must hold", limit,
                             note=note))
        return holds

    def info(self, study, quantity, old, new, note=""):
        """Two numbers with no tolerance (a node count, a time): shown with their plain percentage."""
        plain = plain_percent(old, new) if isinstance(old, (int, float)) and isinstance(new, (int, float)) else None
        self.rows.append(Row(study, quantity, old, new, plain, "info", "info", None, plain, note=note))

    def missing(self, study, what):
        self.rows.append(Row(study, what, None, None, None, "MISSING", "exact", note="file not found or unreadable"))

    # ---- output
    def csv_rows(self):
        for r in self.rows:
            if r.stats is None:
                yield [r.study, r.quantity, _num(r.femm42), _num(r.femm50), _num(r.delta_pct), r.verdict, r.kind,
                       _num(r.limit_pct), _num(r.plain_pct), r.n, r.n_flagged, r.n_plain_flagged, r.note]
            else:
                for label, v in zip(("median |delta|", "p95 |delta|", "max |delta|"), r.stats):
                    last = label.startswith("max")
                    yield [r.study, f"{r.quantity}: {label}", "", "", _num(v), r.verdict if last else "info", r.kind,
                           _num(r.limit_pct), "", r.n if last else "", r.n_flagged if last else "",
                           r.n_plain_flagged if last else "", r.note if last else ""]


CSV_HEADER = ["study", "quantity", "femm42", "femm50", "delta_pct", "verdict", "kind", "limit_pct", "plain_pct", "n",
              "n_flagged", "n_plain_flagged", "note"]
DETAIL_HEADER = ["study", "quantity", "member", "femm42", "femm50", "delta_pct", "plain_pct", "beyond"]


def _num(x):
    """A value for a csv cell: numbers with 10 significant digits, text as it is, None as an empty cell."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    if isinstance(x, (float, np.floating)):
        return f"{float(x):.10g}"
    return x


def _fmt(x, width=11):
    if x is None or x == "":
        return " " * width
    if isinstance(x, str):
        return x[:width].ljust(width)
    if isinstance(x, (int, np.integer)):
        return f"{x:>{width},}"
    ax = abs(x)
    if 1000 <= ax < 1e7:
        return f"{x:,.0f}".rjust(width)
    s = f"{x:.4g}" if (ax >= 1e-3 or ax == 0) and ax < 1e7 else f"{x:.3e}"
    return s.rjust(width)


def print_report(rep: Report, out=print):
    """One block per study: the verdict counts, then every row that is not simply ok, then the ok ones as a short list."""
    studies = list(dict.fromkeys(r.study for r in rep.rows))
    for study in studies:
        rows = [r for r in rep.rows if r.study == study]
        counts = {}
        for r in rows:
            counts[r.verdict] = counts.get(r.verdict, 0) + 1
        out(f"\n=== {study}   ({', '.join(f'{v} {k}' for k, v in counts.items())})")
        for r in rows:
            if r.stats is not None:
                med, p95, mx = r.stats
                lim = f"of {r.limit_pct:g} %" if r.limit_pct is not None else ""
                extra = f"  {r.n_flagged} of {r.n} beyond {r.limit_pct:g} %" if r.n_flagged else ""
                if r.n_plain_flagged and r.n_plain_flagged != r.n_flagged:
                    extra += f" ({r.n_plain_flagged} by the plain rule)"
                out(f"  {r.quantity[:58]:58s} n={r.n:<5d} median {med:7.3f}  p95 {p95:7.3f}  max {mx:7.3f} %   {r.verdict}{extra}"
                    f"{('  | ' + r.note) if r.note else ''}")
            else:
                dp = "" if r.delta_pct is None else f"{r.delta_pct:+9.3f} %"
                out(f"  {r.quantity[:58]:58s} {_fmt(r.femm42)} {_fmt(r.femm50)} {dp:>12s}   {r.verdict}"
                    f"{('  | ' + r.note) if r.note else ''}")


# ----------------------------------------------------------------------------------------------- 3. reading things
def read_text(path: Path) -> str:
    """A text file whatever its encoding: the old Windows logs are UTF-16 or UTF-8 with a byte-order mark."""
    raw = Path(path).read_bytes()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        return raw.decode("utf-16")
    return raw.decode("utf-8-sig", errors="replace")


def find_log(folder: Path, name: str) -> Path | None:
    """The log run_all.py wrote for a step: NN-name.log in NEW_OUT/validation/logs."""
    hits = sorted((Path(folder) / "validation" / "logs").glob(f"*-{name}.log"))
    return hits[0] if hits else None


def log_lines(path: Path | None) -> list[str]:
    """A step's output, without run_all's own header and footer lines (those start with '#')."""
    if path is None or not Path(path).is_file():
        return []
    return [ln for ln in read_text(path).splitlines() if not ln.startswith("#")]


def read_csv(rep: Report, study: str, folder: Path, relative: str, **kw) -> pd.DataFrame | None:
    path = Path(folder) / relative
    if not path.is_file():
        rep.missing(study, f"{relative} in {folder}")
        return None
    return pd.read_csv(path, **kw)


def align(rep: Report, study: str, old: pd.DataFrame, new: pd.DataFrame, keys: list[str], what: str) -> pd.DataFrame:
    """Both tables on their key columns (a missing key is filled with -1). Rows that only one side has are reported."""
    old, new = old.copy(), new.copy()
    for k in keys:
        old[k] = old[k].fillna(-1)
        new[k] = new[k].fillna(-1)
    m = old.merge(new, on=keys, suffixes=("_o", "_n"), how="outer", indicator=True)
    odd = int((m["_merge"] != "both").sum())
    if odd:
        rep.rows.append(Row(study, f"{what}: rows that only one run has", None, odd, None, "MISSING", "exact"))
    return m[m["_merge"] == "both"].reset_index(drop=True)


# ----------------------------------------------------------------------------------------------- 4. FEMM 4.2's printed numbers
# What FEMM 4.2 printed on 2026-10-28 and the README kept (git show "cd2c977:Python to FEMM/README.md", "Expected numbers".
# The path is the one at that commit, where the code sat in a folder called Python to FEMM, and it is quoted because it has
# spaces; the README has moved since, so "cd2c977:README.md" fails). They are the only
# record of these runs: the old out folder has no log for check_geometry, check_vs_fem, check_femm, example, evaluate or mesh_study.
README_GEOMETRY = [          # (quantity, the README's value, regex into the new log, decimals the README gave)
    ("regions", 40, r"^(\d+) regions", 0),
    ("centre rib [mm]", 1.500, r"centre_rib\s+([\d.]+)", 3),
    ("q-rib [mm]", 1.794, r"q_rib\s+([\d.]+)", 3),
    ("outer bridge [mm]", 1.400, r"outer_bridge\s+([\d.]+)", 3),
    ("web, thinnest point [mm]", 1.460, r"^\s*web\s+([\d.]+)", 3),
    ("tooth [mm]", 3.969, r"tooth_min\s+([\d.]+)", 3),
    ("yoke [mm]", 5.698, r"yoke\s+([\d.]+)", 3),
    ("magnet area [mm2]", 261.57, r"magnet_area_mm2\s+([\d.]+)", 2),
    ("winding factor kw1", 0.9452, r"kw1 = ([\d.]+)", 4),
    ("DXF round trip, worst [%]", 0.030, r"worst \S+: kernel .*\(\+([\d.]+)%\)", 3),
    ("phase axis A [deg mech]", 65.0, r"'A': ([\d.]+)", 1),
    ("phase axis B [deg mech]", 5.0, r"'B': ([\d.]+)", 1),
    ("phase axis C [deg mech]", 35.0, r"'C': ([\d.]+)", 1),
]
README_VS_FEM = [
    ("stator steel area, motorgeo against .FEM [%]", -0.009, r"stator_steel .*\(([-+][\d.]+)%\)", 3),
    ("rotor steel area [%]", 0.000, r"rotor_steel .*\(([-+][\d.]+)%\)", 3),
    ("magnet area [%]", 0.000, r"magnet_0_up .*\(([-+][\d.]+)%\)", 3),
    ("slot area, motorgeo [mm2]", 54.127, r"slot 0: motorgeo ([\d.]+) mm2", 3),
    ("slot area, .FEM [mm2]", 54.152, r"\.FEM ([\d.]+) mm2 = faces", 3),
    ("magnet direction 1, .FEM [deg]", 170.5, r"magdir ([\d.]+), the pattern says 165\.5", 1),
    ("magnet direction 2, .FEM [deg]", 305.5, r"magdir ([\d.]+), the pattern says 300\.5", 1),
]
CHECK_FEMM_UNITS = {"psi_d": "mWb", "psi_q": "mWb", "T_age": "N.m", "T_wst": "N.m", "T_dq": "N.m"}     # the order of a check_femm log line
README_CHECK_FEMM = {        # coarse mesh, FEMM 4.2: (id, iq) -> {quantity: value, in CHECK_FEMM_UNITS' units}; 11,644 nodes
    (0, 0): {"psi_d": 11.941, "psi_q": 0.023},
    (0, 10): {"T_age": 0.8061, "T_wst": 0.8059},
    (-6, 8): {"T_age": 0.9110},
    (6, 8): {"T_age": 0.3886},
}
README_NODES = {"coarse": 11644, "fine": 27939}
README_BARRIER = {"with band psi_pm [mWb]": 12.09, "with band T at 10 A MTPA [N.m]": 0.857,
                  "no band psi_pm [mWb]": 16.57, "no band T at 10 A MTPA [N.m]": 1.125}
README_REFERENCE_OPEN = {"no-load psi_d [mWb]": 19.18, "best mean torque at 10 A [N.m]": 1.266}
CHECK_FEMM_LINE = re.compile(r"th_e=\s*0 id=\s*(-?\d+) iq=\s*(-?\d+) \| psi_d=\s*(-?[\d.]+) psi_q=\s*(-?[\d.]+) mWb \| "
                             r"T_age=\s*(-?[\d.]+) T_wst=\s*(-?[\d.]+) T_dq=\s*(-?[\d.]+)")


def parse_check_femm(lines: list[str]) -> dict:
    """The lines of a check_femm log as {(id, iq): {"psi_d": ..., "psi_q": ..., "T_age": ..., "T_wst": ..., "T_dq": ...}}.

    A line looks like  `th_e=   0 id= -6 iq=  8 | psi_d=  5.188 psi_q= 15.883 mWb | T_age= 0.9104 T_wst= 0.9100 T_dq= 0.8208 N.m | 0.3 s`
    (units in CHECK_FEMM_UNITS). Lines that do not look like that are skipped, so a complete log gives four keys. This is the
    only reader of those lines: the coarse and the fine log and the self-consistency check all use it.
    """
    got = {}
    for ln in lines:
        m = CHECK_FEMM_LINE.search(ln)
        if m:
            got[(int(m.group(1)), int(m.group(2)))] = dict(zip(CHECK_FEMM_UNITS, map(float, m.groups()[2:])))
    return got


def match_log(lines: list[str], pattern: str) -> float | None:
    """The first number the regex captures in the log's lines, or None."""
    rx = re.compile(pattern)
    for ln in lines:
        m = rx.search(ln)
        if m:
            return float(m.group(1))
    return None


# ----------------------------------------------------------------------------------------------- 5. the studies
def study_not_femm(rep: Report, old: Path, new: Path):
    """Everything that does not use FEMM must come out exactly as before: the geometry checks and the files they write."""
    s = "not FEMM: identical"
    geo = log_lines(find_log(new, "check_geometry"))
    vs = log_lines(find_log(new, "check_vs_fem"))
    if not geo or not vs:
        rep.missing(s, "the logs of check_geometry and check_vs_fem in NEW_OUT/validation/logs")
    for label, expected, rx, dec in README_GEOMETRY:
        got = match_log(geo, rx)
        if got is None:
            rep.missing(s, f"check_geometry: {label}")
        else:
            rep.value(s, f"check_geometry: {label}", round(expected, dec), round(got, dec), "exact", note="README, 10-28")
    for label, expected, rx, dec in README_VS_FEM:
        got = match_log(vs, rx)
        if got is None:
            rep.missing(s, f"check_vs_fem: {label}")
        else:
            rep.value(s, f"check_vs_fem: {label}", round(expected, dec), round(got, dec), "exact", note="README, 10-28")
    text = "\n".join(geo)
    for label, needle in (("v_half_deg 55 is rejected", "v_half_deg 55 rejected: magnet slot reaches the bridge"),
                          ("mag_x 38.4 is rejected", "mag_x 38.4 violations: ['q_rib = 0.00 < 1.0', 'web = 0.45 < 1.0']")):
        rep.value(s, f"check_geometry: {label}", "yes", "yes" if needle in text else "no", "exact", note="README, 10-28")
    # the files the scripts write
    for name in ("equations.txt", "ipm_sector.msh"):
        a, b = old / name, new / name
        if a.is_file() and b.is_file():
            rep.value(s, f"file {name}", "bytes", "bytes" if a.read_bytes() == b.read_bytes() else "other bytes", "exact")
        else:
            rep.missing(s, name)
    if (old / "sector.png").is_file() and (new / "sector.png").is_file():
        n, total, _ = pixel_difference(old / "sector.png", new / "sector.png")
        rep.value(s, "file sector.png (pixels that differ)", 0, n, "exact", note=f"of {total:,}")
    else:
        rep.missing(s, "sector.png")
    if (old / "ipm_full.dxf").is_file() and (new / "ipm_full.dxf").is_file():
        a, b = dxf_entities(old / "ipm_full.dxf"), dxf_entities(new / "ipm_full.dxf")
        rep.value(s, "file ipm_full.dxf (lines and arcs that differ)", 0, sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b)),
                  "exact", note=f"{len(b)} entities; only dates, ids and the order of two housekeeping objects differ in the text")
    else:
        rep.missing(s, "ipm_full.dxf")
    if (old / "check_vs_fem.png").is_file() and (new / "check_vs_fem.png").is_file():
        n, total, box = pixel_difference(old / "check_vs_fem.png", new / "check_vs_fem.png")
        rep.value(s, "file check_vs_fem.png (pixels that differ)", 0, n, "exact",
                  note=f"of {total:,}; all inside x {box[0]}-{box[1]}, y {box[2]}-{box[3]}" if n else "")
    else:
        rep.missing(s, "check_vs_fem.png")


def pixel_difference(a: Path, b: Path):
    """(pixels that differ, pixels in all, bounding box x0, x1, y0, y1 of the differences)."""
    from PIL import Image
    x = np.asarray(Image.open(a).convert("RGBA")).astype(int)
    y = np.asarray(Image.open(b).convert("RGBA")).astype(int)
    if x.shape != y.shape:
        return x.shape[0] * x.shape[1], x.shape[0] * x.shape[1], (0, 0, 0, 0)
    d = (x != y).any(axis=2)
    if not d.any():
        return 0, int(d.size), (0, 0, 0, 0)
    ys, xs = np.nonzero(d)
    return int(d.sum()), int(d.size), (int(xs.min()), int(xs.max()), int(ys.min()), int(ys.max()))


def dxf_entities(path: Path) -> list:
    """Every line and arc of a DXF with its layer and its numbers (nothing about dates, ids or handles)."""
    import ezdxf
    out = []
    for e in ezdxf.readfile(str(path)).modelspace():
        attrs = e.dxf.all_existing_dxf_attribs()
        out.append((e.dxftype(), tuple(sorted((k, tuple(round(c, 9) for c in v) if hasattr(v, "__iter__") and not isinstance(v, str)
                                               else (round(v, 9) if isinstance(v, float) else v))
                                              for k, v in attrs.items() if k not in ("handle", "owner")))))
    return out


def count_nodes(path: Path):
    """(nodes, triangles, the parsed solution) of a .ans file."""
    sol = femfile.read_ans(path)
    return len(sol["xy"]), len(sol["tri"]), sol


def mean_triangle_area_by_material(sol) -> dict:
    """The mean area [mm2] of the triangles in each material of a solved model."""
    xy, tri = sol["xy"], sol["tri"]
    x, y = xy[tri, 0], xy[tri, 1]
    area = 0.5 * np.abs((x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (x[:, 2] - x[:, 0]) * (y[:, 1] - y[:, 0]))
    material = np.array([sol["labels"][k]["material"] for k in sol["label"]])
    return {m: float(area[material == m].mean()) for m in sorted(set(material)) if m != "<No Mesh>"}


def study_mesh(rep: Report, old: Path, new: Path):
    """How many nodes each FEMM made from the same model. This is the one real difference between the two (see the report, section 5)."""
    s = "mesh: nodes and triangles"
    cases = [("coarse mesh, the default model", "scaling/n1_w0.ans", "coarse"),
             ("fine mesh (smart mesh on), the default model", "barrier/fields/noload_with_band.ans", "fine"),
             ("the hand-drawn full model (360 deg), one position", "validate/hand_000.ans", "fine")]
    areas = {}
    for label, rel, _ in cases:
        a, b = old / rel, new / rel
        if not (a.is_file() and b.is_file()):
            rep.missing(s, rel)
            continue
        (na, ta, sa), (nb, tb, sb) = count_nodes(a), count_nodes(b)
        rep.info(s, f"nodes, {label}", na, nb)
        rep.info(s, f"triangles, {label}", ta, tb)
        if rel.startswith("scaling"):
            areas = (mean_triangle_area_by_material(sa), mean_triangle_area_by_material(sb))
    for material in (areas[0] if areas else {}):
        if material in areas[1]:
            rep.info(s, f"mean triangle area in {material}, coarse [mm2]", areas[0][material], areas[1][material])
    for name, expected in README_NODES.items():                         # what the README recorded on 10-28, against the old files
        rel = {"coarse": "scaling/n1_w0.ans", "fine": "barrier/fields/noload_with_band.ans"}[name]
        if (old / rel).is_file():
            rep.value(s, f"nodes of the {name} mesh in the old files, against the README's {expected:,}", expected, count_nodes(old / rel)[0],
                      "exact", note="README, 10-28")
    # the flux explorer's 17 rotors, one mesh each, as the sweep stored them
    ratios = []
    for g in range(17):
        a, b = old / f"flux_explorer/g{g:02d}_air.npz", new / f"flux_explorer/g{g:02d}_air.npz"
        if a.is_file() and b.is_file():
            ratios.append((len(np.load(a)["xy"]), len(np.load(b)["xy"])))
    if ratios:
        rep.info(s, "nodes, flux explorer, the reference design (rotor 0)", ratios[0][0], ratios[0][1])
        rep.info(s, "nodes, flux explorer, smallest of the 17 rotors", min(r[0] for r in ratios), min(r[1] for r in ratios))
        rep.info(s, "nodes, flux explorer, largest of the 17 rotors", max(r[0] for r in ratios), max(r[1] for r in ratios))
        rep.info(s, "nodes, flux explorer, mean of the 17 rotors", sum(r[0] for r in ratios) / len(ratios),
                 sum(r[1] for r in ratios) / len(ratios))


def study_check_femm(rep: Report, old: Path, new: Path):
    """check_femm.py: four solves at one rotor position. FEMM 4.2's coarse numbers are the README's; the fine ones are in the study."""
    s = "check_femm, coarse mesh"
    got = parse_check_femm(log_lines(find_log(new, "check_femm")))
    if len(got) != 4:
        rep.missing(s, "the four lines of check_femm's log")
        return
    for point, expected in README_CHECK_FEMM.items():
        for q, v in expected.items():
            scale = README_CHECK_FEMM[(0, 0)]["psi_d"] if q.startswith("psi") else None
            rep.value(s, f"id={point[0]:+d} iq={point[1]}: {q} [{CHECK_FEMM_UNITS[q]}]", v, got[point][q],
                      "flux" if q.startswith("psi") else "mean_torque", "coarse", scale=scale,
                      note="psi_q at no load should be 0: judged against psi_d" if q == "psi_q" else "")
    # fine: the same model, the same four points, are in the barrier study's rows (theta_e = 0)
    s = "check_femm --fine (against the barrier study's rows)"
    po = read_csv(rep, s, old, "barrier/phase1.csv")
    pn = read_csv(rep, s, new, "barrier/phase1.csv")
    if po is None or pn is None:
        return
    fine = parse_check_femm(log_lines(find_log(new, "check_femm_fine")))
    for label, pick in (("no load", lambda d: d[(d.study == "noload") & (d.variant == "with band") & (d.theta_e == 0)]),
                        ("iq = 10 A", lambda d: d[(d.study == "gamma") & (d.variant == "with band") & (d.theta_e == 0)
                                                   & (d.I == 10) & (d.gamma == 0)])):
        a, b = pick(po), pick(pn)
        if len(a) != 1 or len(b) != 1:
            rep.missing(s, f"the {label} row")
            continue
        a, b = a.iloc[0], b.iloc[0]
        point = (0, 0) if label == "no load" else (0, 10)
        for q, kind, scale in (("psi_d", "flux", a.psi_d), ("psi_q", "flux", a.psi_d), ("T_age", "mean_torque", None),
                               ("T_wst", "mean_torque", None)):
            if label == "no load" and q.startswith("T"):
                continue
            f = 1e3 if q.startswith("psi") else 1.0                  # flux linkages in mWb, torques in N.m
            rep.value(s, f"{label}: {q} [{'N.m' if q.startswith('T') else 'mWb'}]", a[q] * f, b[q] * f, kind, "fine",
                      scale=None if scale is None else scale * f)
            shown = fine.get(point, {}).get(q)
            if shown is not None and abs(shown - b[q] * (1e3 if q.startswith("psi") else 1)) > 6e-4:
                rep.rows.append(Row(s, f"{label}: {q}, the log of check_femm --fine", b[q], shown, None, "DIFFERS", "exact",
                                    note="the printed number is not the study's number"))
    # the second, third and fourth lines need no old counterpart: they are the README's coarse numbers above


def study_printed_only(rep: Report, old: Path, new: Path):
    """example.py, evaluate.py and mesh_study.py print their numbers and keep nothing, and FEMM 4.2's run of them was not recorded
    anywhere (the old out folder has no log, the README gives no number). So they have a FEMM 5.0 column and no FEMM 4.2 one.
    The nearest FEMM 4.2 numbers are in the barrier study and are named in the note; they are on the fine mesh, so they are a
    neighbour, not a counterpart."""
    s = "no FEMM 4.2 record (FEMM 5.0 numbers only)"
    po = read_csv(rep, s, old, "barrier/phase1.csv")
    near = None
    if po is not None:
        g = po[(po.study == "gamma") & (po.variant == "with band") & (po.I == 10) & (po.gamma == 0)]
        near = float(g.T_age.mean()) if len(g) else None
    example = log_lines(find_log(new, "example"))
    m = re.search(r"mean ([\d.]+) N\.m, ripple ([\d.]+) N\.m", "\n".join(example))
    if m:
        for label, v in (("example.py: mean torque, 30 positions, 10 A on the q axis [N.m]", float(m.group(1))),
                         ("example.py: torque ripple, peak to peak [N.m]", float(m.group(2)))):
            rep.rows.append(Row(s, label, None, v, None, "info", "info",
                                note=f"FEMM 4.2's fine-mesh mean torque for these currents (barrier study, 6 positions): {near:.4f}" if near and "mean" in label else ""))
    for ln in log_lines(find_log(new, "evaluate")):
        d = re.search(r"'T_mean': ([\d.]+), 'T_ripple': ([\d.]+)", ln)
        if d:
            label = "baseline" if "mag_y" not in ln and "web" not in ln else ("mag_y 5.0" if "mag_y" in ln else "web 1.2")
            rep.rows.append(Row(s, f"evaluate.py: {label}: mean torque [N.m]", None, float(d.group(1)), None, "info", "info"))
            rep.rows.append(Row(s, f"evaluate.py: {label}: torque ripple [N.m]", None, float(d.group(2)), None, "info", "info"))
    for ln in log_lines(find_log(new, "mesh_study")):
        d = re.match(r"(coarse|fine): .*mean ([\d.]+), p-p ([\d.]+) N\.m", ln)
        if d:
            rep.rows.append(Row(s, f"mesh_study.py, {d.group(1)} mesh: mean torque [N.m]", None, float(d.group(2)), None, "info", "info"))
            rep.rows.append(Row(s, f"mesh_study.py, {d.group(1)} mesh: torque ripple, peak to peak [N.m]", None, float(d.group(3)), None, "info", "info"))
        h = re.search(r"harmonics per electrical cycle: 6th ([\d.]+), 12th ([\d.]+), 18th ([\d.]+), 24th ([\d.]+)", ln)
        if h:
            rep.rows.append(Row(s, "mesh_study.py: harmonics 6, 12, 18, 24 of the same sweep [N.m]", None, "/".join(h.groups()), None, "info", "info"))


def study_barrier(rep: Report, old: Path, new: Path):
    """barrier_study.py: with and without the air band, fine mesh, 1,070 solves. summary.csv holds the headline numbers."""
    mesh = "fine"
    s = "barrier: summary.csv"
    so = read_csv(rep, s, old, "barrier/summary.csv", index_col=0)
    sn = read_csv(rep, s, new, "barrier/summary.csv", index_col=0)
    if so is not None and sn is not None:
        for variant in ("with band", "no band"):
            for q in so.index:
                if q not in sn.index:
                    rep.missing(s, f"{q} in NEW")
                    continue
                rep.value(s, f"{variant}: {q}", float(so.loc[q, variant]), float(sn.loc[q, variant]), summary_kind(q), mesh)
        readme = {"with band psi_pm [mWb]": ("with band", "psi_pm_mWb"), "with band T at 10 A MTPA [N.m]": ("with band", "T_mtpa_10A"),
                  "no band psi_pm [mWb]": ("no band", "psi_pm_mWb"), "no band T at 10 A MTPA [N.m]": ("no band", "T_mtpa_10A")}
        for label, (variant, q) in readme.items():        # the README's headline, as FEMM 4.2 wrote it, against the old csv
            rep.value("barrier: the README's headline", label, README_BARRIER[label], round(float(so.loc[q, variant]), 3 if "T" in q else 2),
                      "exact", note="README, 10-28, against the old summary.csv")
    po = read_csv(rep, "barrier: phase 1", old, "barrier/phase1.csv")
    pn = read_csv(rep, "barrier: phase 1", new, "barrier/phase1.csv")
    if po is not None and pn is not None:
        barrier_phase1(rep, po, pn, mesh)
    qo = read_csv(rep, "barrier: phase 2", old, "barrier/phase2.csv")
    qn = read_csv(rep, "barrier: phase 2", new, "barrier/phase2.csv")
    if qo is not None and qn is not None:
        barrier_phase2(rep, qo, qn, mesh)
    gaps(rep, old, new, mesh)


def summary_kind(q: str) -> str:
    if q.startswith(("mtpa_gamma",)):
        return "angle"
    if q == "cogging_pp_mNm":
        return "cogging"
    if q in ("ripple_pp_Nm", "ripple_pct"):
        return "ripple"
    if q in ("emf_thd_pct", "B_gap_thd_pct"):
        return "harmonic"
    if q.startswith(("T_", "T_mean")):
        return "mean_torque"
    return "flux"                                       # psi, EMF, ke, L, saliency, air-gap B1


def barrier_phase1(rep: Report, po: pd.DataFrame, pn: pd.DataFrame, mesh: str):
    s = "barrier: phase 1 (no load, cogging, torque vs current angle)"
    m = align(rep, s, po, pn, ["study", "variant", "theta_e", "I", "gamma"], "phase1.csv")
    for variant in ("with band", "no band"):
        v = m[m.variant == variant]
        nl = v[v.study == "noload"].sort_values("theta_e")
        flux_scale = float(np.abs(nl.psi_d_o).max())
        for q in ("psi_a", "psi_d", "psi_q"):
            peak = float(np.abs(nl[q + "_o"]).max()) if q == "psi_a" else flux_scale
            rep.group(s, f"{variant}, no load: {q} at 120 rotor positions", nl[q + "_o"], nl[q + "_n"], "flux", mesh, scale=peak,
                      wave=True, keys=[f"theta_e={t:g}" for t in nl.theta_e],
                      note="% of the " + ("waveform's peak" if q == "psi_a" else "no-load flux linkage"))
        pp = float(np.ptp(nl.T_age_o))
        rep.group(s, f"{variant}, no load: torque at 120 rotor positions (cogging)", nl.T_age_o, nl.T_age_n, "cogging", mesh, scale=pp,
                  wave=True, keys=[f"theta_e={t:g}" for t in nl.theta_e], note="% of the curve's peak-to-peak")
        cg = v[v.study == "cogging"].sort_values("theta_e")
        pp = float(np.ptp(cg.T_age_o))
        rep.group(s, f"{variant}, cogging: torque at 41 rotor positions", cg.T_age_o, cg.T_age_n, "cogging", mesh, scale=pp, wave=True,
                  keys=[f"theta_e={t:g}" for t in cg.theta_e], note="% of the curve's peak-to-peak")
        g = v[v.study == "gamma"]
        grid = g.groupby(["I", "gamma"])[["T_age_o", "T_age_n", "psi_d_o", "psi_d_n", "psi_q_o", "psi_q_n"]].mean().reset_index()
        keys = [f"I={i:g} gamma={a:g}" for i, a in zip(grid.I, grid.gamma)]
        rep.group(s, f"{variant}, torque vs current angle: mean torque, {len(grid)} (I, gamma) points", grid.T_age_o, grid.T_age_n,
                  "mean_torque", mesh, scale=float(np.abs(grid.T_age_o).max()), keys=keys, note="mean over 6 rotor positions")
        flux = float(max(np.abs(grid.psi_d_o).max(), np.abs(grid.psi_q_o).max()))
        for q in ("psi_d", "psi_q"):
            rep.group(s, f"{variant}, torque vs current angle: {q}, {len(grid)} points", grid[q + "_o"], grid[q + "_n"], "flux", mesh,
                      scale=flux, keys=keys, note="mean over 6 rotor positions")
        rep.group(s, f"{variant}, torque vs current angle: torque at {len(g)} single rotor positions", g.T_age_o, g.T_age_n, "mean_torque",
                  mesh, scale=float(np.abs(g.T_age_o).max()),
                  keys=[f"I={i:g} gamma={a:g} theta_e={t:g}" for i, a, t in zip(g.I, g.gamma, g.theta_e)],
                  note="one position: includes the ripple")


def barrier_phase2(rep: Report, qo: pd.DataFrame, qn: pd.DataFrame, mesh: str):
    s = "barrier: phase 2 (ripple at the best current angle, snapshots)"
    for variant in ("with band", "no band"):
        ro = qo[(qo.study == "ripple") & (qo.variant == variant)].sort_values("theta_e")
        rn = qn[(qn.study == "ripple") & (qn.variant == variant)].sort_values("theta_e")
        if len(ro) != len(rn) or len(ro) == 0:
            rep.missing(s, f"{variant}: ripple rows")
            continue
        rep.value(s, f"{variant}: current angle of the best torque at 10 A [deg]", float(ro.gamma.iloc[0]), float(rn.gamma.iloc[0]), "angle",
                  mesh, note="each run found its own; the torque is flat there")
        To, Tn = ro.T_age.values[:-1], rn.T_age.values[:-1]                 # 0..59 deg: one 60-degree period
        rep.group(s, f"{variant}: torque at 61 rotor positions", ro.T_age, rn.T_age, "ripple", mesh, scale=float(np.ptp(ro.T_age)), wave=True,
                  keys=[f"theta_e={t:g}" for t in ro.theta_e], note="% of the curve's peak-to-peak")
        rep.value(s, f"{variant}: mean torque over the period [N.m]", float(To.mean()), float(Tn.mean()), "mean_torque", mesh)
        rep.value(s, f"{variant}: torque ripple, peak to peak [N.m]", float(np.ptp(To)), float(np.ptp(Tn)), "ripple", mesh)
        ho, hn = harmonic_amplitudes(To), harmonic_amplitudes(Tn)
        for k in range(1, 7):
            rep.value(s, f"{variant}: torque harmonic {6 * k} per electrical cycle [N.m]", float(ho[k]), float(hn[k]), "harmonic", mesh,
                      note=f"{100 * ho[k] / ho[1:].max():.1f} % of the largest harmonic")
        for q in ("psi_d", "psi_q"):
            rep.value(s, f"{variant}: mean {q} over the period [mWb]", float(ro[q].mean() * 1e3), float(rn[q].mean() * 1e3), "flux", mesh,
                      scale=float(abs(ro.psi_d.mean() * 1e3)))
    for variant in ("with band", "no band"):
        a = qo[(qo.study == "snapshot") & (qo.variant == variant)]
        b = qn[(qn.study == "snapshot") & (qn.variant == variant)]
        if len(a) == 1 and len(b) == 1:
            a, b = a.iloc[0], b.iloc[0]
            rep.value(s, f"{variant}: all-q current, 10 A: torque [N.m]", float(a.T_age), float(b.T_age), "mean_torque", mesh)
            rep.value(s, f"{variant}: all-q current, 10 A: psi_d [mWb]", float(a.psi_d * 1e3), float(b.psi_d * 1e3), "flux", mesh,
                      scale=float(max(abs(a.psi_d), abs(a.psi_q)) * 1e3))
            rep.value(s, f"{variant}: all-q current, 10 A: psi_q [mWb]", float(a.psi_q * 1e3), float(b.psi_q * 1e3), "flux", mesh,
                      scale=float(max(abs(a.psi_d), abs(a.psi_q)) * 1e3))


def gaps(rep: Report, old: Path, new: Path, mesh: str):
    """The air-gap flux density samples (360 points per 90-degree sector), six files."""
    s = "barrier: air-gap flux density"
    for case in ("noload", "qaxis", "mtpa"):
        for variant in ("with_band", "no_band"):
            name = f"barrier/gap_{case}_{variant}.csv"
            a, b = read_csv(rep, s, old, name), read_csv(rep, s, new, name)
            if a is None or b is None:
                continue
            peak = float(np.abs(a.br).max())
            rep.group(s, f"{case}, {variant.replace('_', ' ')}: radial flux density at {len(a)} points", a.br, b.br, "flux", mesh,
                      scale=peak, wave=True, keys=[f"deg={d:g}" for d in a.deg], note="% of the waveform's peak")
            if case == "noload":
                ha, hb = harmonic_amplitudes(a.br.values), harmonic_amplitudes(b.br.values)
                top = np.argsort(ha[2:60])[::-1][:5] + 2
                for k in sorted(top):
                    rep.value(s, f"noload, {variant.replace('_', ' ')}: harmonic {k} of the radial flux density [T]", float(ha[k]), float(hb[k]),
                              "harmonic", mesh, note=f"{100 * ha[k] / ha[1]:.1f} % of the fundamental")


def study_reference_open(rep: Report, old: Path, new: Path):
    """reference_open.py: the magnet inside the band (the 10-22 prototype), same tests at 10 A, fine mesh."""
    s = "reference_open (magnet inside the band)"
    ro = read_csv(rep, s, old, "barrier/reference_open.csv")
    rn = read_csv(rep, s, new, "barrier/reference_open.csv")
    if ro is None or rn is None:
        return
    mesh = "fine"
    def headline(df):
        nl = df[df.study == "noload"]
        g = df[df.study == "gamma"].groupby("gamma").T_age.mean()
        return float(nl.psi_d.mean() * 1e3), float(g.max()), float(g.idxmax())
    (pa, ta, ga), (pb, tb, gb) = headline(ro), headline(rn)
    rep.value(s, "no-load psi_d [mWb] (the old log: 19.18)", pa, pb, "flux", mesh, note="the README says 19.18")
    rep.value(s, "best mean torque at 10 A [N.m]", ta, tb, "mean_torque", mesh,
              note=f"the old log says {ta:.4f}; the README says {README_REFERENCE_OPEN['best mean torque at 10 A [N.m]']}")
    rep.value(s, "current angle of the best torque [deg]", ga, gb, "angle", mesh, note="from the 7.5-degree grid")
    m = align(rep, s, ro, rn, ["study", "theta_e", "I", "gamma"], "reference_open.csv")
    nl = m[m.study == "noload"].sort_values("theta_e")
    rep.group(s, "no-load psi_a at 24 rotor positions", nl.psi_a_o, nl.psi_a_n, "flux", mesh, scale=float(np.abs(nl.psi_a_o).max()), wave=True,
              keys=[f"theta_e={t:g}" for t in nl.theta_e], note="% of the waveform's peak")
    g = m[m.study == "gamma"]
    grid = g.groupby(["I", "gamma"])[["T_age_o", "T_age_n", "psi_d_o", "psi_d_n", "psi_q_o", "psi_q_n"]].mean().reset_index()
    keys = [f"gamma={a:g}" for a in grid.gamma]
    rep.group(s, "mean torque at 13 current angles", grid.T_age_o, grid.T_age_n, "mean_torque", mesh, scale=float(np.abs(grid.T_age_o).max()), keys=keys)
    flux = float(max(np.abs(grid.psi_d_o).max(), np.abs(grid.psi_q_o).max()))
    for q in ("psi_d", "psi_q"):
        rep.group(s, f"{q} at 13 current angles", grid[q + "_o"], grid[q + "_n"], "flux", mesh, scale=flux, keys=keys)


def study_validate(rep: Report, old: Path, new: Path):
    """validate_full.py: the hand-drawn 360-degree model against the 90-degree slice, no load and 10 A, fine mesh."""
    s = "validate_full (slice against the full model)"
    vo = read_csv(rep, s, old, "validate/validate.csv")
    vn = read_csv(rep, s, new, "validate/validate.csv")
    if vo is None or vn is None:
        return
    mesh = "fine"
    m = align(rep, s, vo, vn, ["model", "theta_e", "I"], "validate.csv")
    for model in ("hand-drawn 360", "motorgeo 90"):
        racy = " (FEMM 4.2's slice rows came from the script whose jobs shared file names: see README)" if model == "motorgeo 90" else ""
        d0, d10 = m[(m.model == model) & (m.I == 0)].sort_values("theta_e"), m[(m.model == model) & (m.I == 10)].sort_values("theta_e")
        flux_scale = float(max(np.abs(d0.psi_a_o).max(), np.abs(d10.psi_a_o).max()))
        for label, d in (("no load", d0), ("10 A", d10)):
            keys = [f"theta_e={t:g}" for t in d.theta_e]
            rep.group(s, f"{model}, {label}: phase-A flux linkage at {len(d)} rotor positions", d.psi_a_o, d.psi_a_n, "flux", mesh,
                      scale=float(np.abs(d.psi_a_o).max()), wave=True, keys=keys, note="% of the peak" + racy)
        rep.group(s, f"{model}, 10 A: torque at 24 rotor positions", d10.T_o, d10.T_n, "mean_torque", mesh, scale=float(np.abs(d10.T_o).max()), wave=True,
                  keys=[f"theta_e={t:g}" for t in d10.theta_e], note="% of the peak torque" + racy)
        rep.group(s, f"{model}, no load: torque at 24 rotor positions (cogging)", d0.T_o, d0.T_n, "cogging", mesh, scale=float(np.ptp(d0.T_o)),
                  wave=True, keys=[f"theta_e={t:g}" for t in d0.theta_e], note="% of the curve's peak-to-peak" + racy)
        rep.value(s, f"{model}: peak phase-A flux linkage, 10 A [mWb]", float(np.abs(d10.psi_a_o).max() * 1e3), float(np.abs(d10.psi_a_n).max() * 1e3),
                  "flux", mesh)
        rep.value(s, f"{model}: peak torque, 10 A [N.m]", float(np.abs(d10.T_o).max()), float(np.abs(d10.T_n).max()), "mean_torque", mesh)
        rep.value(s, f"{model}: torque peak to peak, no load [N.m] (cogging)", float(np.ptp(d0.T_o)), float(np.ptp(d0.T_n)), "cogging", mesh)
    # FEMM 4.2's slice rows came from the racing script; validate_fix_zeros re-solved the 11 no-load rows whose flux linkage was 0.0
    # (it left a fix_NNN.fem for each). If the race had corrupted others, they would disagree more with FEMM 5.0's clean rows.
    fixed = sorted(int(x.group(1)) for x in (re.match(r"fix_(\d+)\.fem$", f.name) for f in (Path(old) / "validate").glob("fix_*.fem")) if x)
    s0 = m[(m.model == "motorgeo 90") & (m.I == 0)]
    if fixed and len(s0):
        gap = 100 * (s0.psi_a_n - s0.psi_a_o).abs() / float(s0.psi_a_o.abs().max())
        mask = s0.theta_e.isin(fixed)
        rep.info(s, f"slice, no load: largest |psi_a difference| in the {int(mask.sum())} rows FEMM 4.2 had to re-solve, and in the other "
                    f"{int((~mask).sum())} [% of peak]", float(gap[~mask].max()), float(gap[mask].max()),
                 note="femm42 column: the other rows; femm50 column: the re-solved ones. Neither set stands out")
    for tag, col in (("nodes, hand-drawn model, mean over the 48 solves", "nodes"),):
        h = vo[vo.model == "hand-drawn 360"][col].mean(), vn[vn.model == "hand-drawn 360"][col].mean()
        rep.info(s, tag, float(h[0]), float(h[1]))
    for model in ("hand-drawn 360", "motorgeo 90"):
        rep.info(s, f"seconds per solve, {model} (8 workers)", float(vo[vo.model == model].seconds.mean()), float(vn[vn.model == model].seconds.mean()))


def slice_vs_full(vdf: pd.DataFrame) -> dict:
    """How far the 90-degree slice is from the full model, in one table (old or new): the plan's acceptance numbers.

    flux: the largest |slice - full| of the phase-A flux linkage at any rotor position, as a percentage of the full model's peak,
          for no load and for 10 A (the README's 1.1 %).
    torque: peak torque of the slice against the full model's at 10 A (the README's 1.5 %).
    """
    out = {}
    for label, current in (("no load", 0.0), ("10 A", 10.0)):
        full = vdf[(vdf.model == "hand-drawn 360") & (vdf.I == current)].sort_values("theta_e")
        part = vdf[(vdf.model == "motorgeo 90") & (vdf.I == current)].sort_values("theta_e")
        out[f"flux, {label}"] = 100 * float(np.abs(part.psi_a.values - full.psi_a.values).max() / np.abs(full.psi_a).max())
        out[f"peak flux, {label}"] = 100 * float(np.abs(part.psi_a).max() / np.abs(full.psi_a).max() - 1)
        if current:
            out["peak torque, 10 A"] = 100 * float(np.abs(part["T"]).max() / np.abs(full["T"]).max() - 1)
    return out


def study_scaling(rep: Report, old: Path, new: Path):
    """worker_scaling.py: 48 coarse solves on 1, 2, 4, 6, 8 and 10 workers."""
    s = "worker_scaling"
    so, sn = read_csv(rep, s, old, "scaling/scaling.csv"), read_csv(rep, s, new, "scaling/scaling.csv")
    if so is None or sn is None:
        return
    for _, a in so.iterrows():
        b = sn[sn.workers == a.workers]
        if b.empty:
            rep.missing(s, f"{int(a.workers)} workers")
            continue
        b = b.iloc[0]
        rep.info(s, f"seconds for 48 solves on {int(a.workers)} workers", float(a.wall_s), float(b.wall_s))
        if a.workers > 1:
            rep.info(s, f"speed-up on {int(a.workers)} workers", float(so.wall_s.iloc[0] / a.wall_s), float(sn.wall_s.iloc[0] / b.wall_s))
    rep.check(s, "largest |T - T(1 worker)| over every worker count [N.m]", float(sn.max_diff_Nm.max()), 0.0, old=float(so.max_diff_Nm.max()),
              note="the plan: bit-identical, so exactly 0")


def study_optimize(rep: Report, old: Path, new: Path):
    """optimize_demo.py: a small NSGA-II run. Generations 0 and 1 are the same candidates (same seed), so they compare row by row."""
    s = "optimize_demo"
    eo, en = read_csv(rep, s, old, "opt/evaluations.csv"), read_csv(rep, s, new, "opt/evaluations.csv")
    if eo is None or en is None:
        return
    mesh = "coarse"                                    # the optimizer uses the coarse mesh
    cols = ["mag_y", "mag_x", "v_half_deg", "web"]
    for gen in range(0, 7):
        a, b = eo[eo.gen == gen].reset_index(drop=True), en[en.gen == gen].reset_index(drop=True)
        same = len(a) == len(b) and np.allclose(a[cols].values, b[cols].values, rtol=0, atol=1e-12)
        label = f"generation {gen}"
        if gen <= 1:
            rep.value(f"optimize_demo: {label}, row by row", "the 8 candidates are the same designs", "yes", "yes" if same else "no", "exact")
            if not same:
                continue
            ok = a.feasible.values & b.feasible.values
            rep.value(f"optimize_demo: {label}, row by row", "feasible in both runs", int(a.feasible.sum()), int(b.feasible.sum()), "exact")
            if ok.any():
                tb = float(np.abs(a.T_mean[ok]).max())
                rep.group(f"optimize_demo: {label}, row by row", "mean torque of the feasible candidates", a.T_mean[ok], b.T_mean[ok], "mean_torque",
                          mesh, scale=tb, keys=[f"candidate {i}" for i in np.nonzero(ok)[0]])
                rep.group(f"optimize_demo: {label}, row by row", "torque ripple of the feasible candidates", a.ripple[ok], b.ripple[ok], "ripple",
                          mesh, scale=float(a.ripple[ok].max()), keys=[f"candidate {i}" for i in np.nonzero(ok)[0]])
                rep.group(f"optimize_demo: {label}, row by row", "best current angle of the feasible candidates", a.gamma[ok], b.gamma[ok], "angle",
                          mesh, scale=float(a.gamma[ok].max()), keys=[f"candidate {i}" for i in np.nonzero(ok)[0]])
        else:
            rep.rows.append(Row("optimize_demo: later generations (not row by row)", f"{label}: the same 8 designs in both runs?", None, None, None,
                                "info", "info", note="yes" if same else "no: selection depends on the torques, so the two searches part ways"))
    for gen in range(0, 7):
        a, b = eo[(eo.gen == gen) & eo.feasible], en[(en.gen == gen) & en.feasible]
        if len(a) and len(b):
            rep.info("optimize_demo: later generations (not row by row)", f"best torque in generation {gen} [N.m]", float(a.T_mean.max()),
                     float(b.T_mean.max()))
    rep.info("optimize_demo: later generations (not row by row)", "best torque of all 56 evaluations [N.m]", float(eo.T_mean.max()), float(en.T_mean.max()))
    # the geometry columns of every row are computed without FEMM: they must be equal for equal designs
    for gen in (0, 1):
        a, b = eo[eo.gen == gen].reset_index(drop=True), en[en.gen == gen].reset_index(drop=True)
        if len(a) == len(b) and np.allclose(a[cols].values, b[cols].values, rtol=0, atol=1e-12):
            geo = ["centre_rib", "q_rib", "outer_bridge", "web_mm"]
            same = bool(np.allclose(a[geo].fillna(-1).values, b[geo].fillna(-1).values, rtol=0, atol=1e-12))
            rep.value(f"optimize_demo: generation {gen}, row by row", "geometry numbers (ribs, bridge, web: no FEMM)", "equal",
                      "equal" if same else "different", "exact")


def study_flux_explorer(rep: Report, old: Path, new: Path):
    """flux_explorer_sweep.py and sweep120: 2,040 + 306 coarse-mesh solves, one row per state."""
    mesh = "coarse"
    for name, label in (("flux_explorer/states.csv", "flux explorer: 2,040 states"), ("flux_explorer/states_p120.csv", "flux explorer: 306 states at 120 deg")):
        a, b = read_csv(rep, label, old, name), read_csv(rep, label, new, name)
        if a is None or b is None:
            continue
        m = align(rep, label, a, b, ["gid", "band", "I", "gamma", "theta_e"], name)
        keys = [f"g{int(g):02d} {bd} I={i:g} gamma={ga:g} theta_e={t:g}" for g, bd, i, ga, t in zip(m.gid, m.band, m.I, m.gamma, m.theta_e)]
        t_scale = float(np.abs(m.T_age_o).max())
        phase = ["psi_a", "psi_b", "psi_c"]
        dq_scale = float(max(np.abs(m.psi_d_o).max(), np.abs(m.psi_q_o).max()))
        # worst rotors first, so that the report can say where the large differences are
        d = 100 * (m.T_age_n - m.T_age_o) / np.maximum(np.abs(m.T_age_o), FLOOR * t_scale)
        worst = m.assign(d=d).groupby(["gid", "band"]).d.apply(lambda x: x.abs().max()).sort_values(ascending=False).head(3)
        rep.group(label, "torque T_age", m.T_age_o, m.T_age_n, "mean_torque", mesh, scale=t_scale, keys=keys,
                  note=f"largest torque {t_scale:.3f} N.m; worst rotors " + ", ".join(f"g{int(g):02d} {bd} {v:.2f} %" for (g, bd), v in worst.items()))
        for c in ("psi_d", "psi_q"):
            rep.group(label, f"flux linkage {c}", m[c + "_o"], m[c + "_n"], "flux", mesh, scale=dq_scale, keys=keys,
                      note=f"largest d or q flux linkage {dq_scale * 1e3:.2f} mWb")
        # the three phase flux linkages are sines of the rotor angle: every value is judged against the waveform's peak
        peak = float(max(np.abs(m[c + "_o"]).max() for c in phase))
        rep.group(label, "phase flux linkages psi_a, psi_b, psi_c (3 values per state)", np.concatenate([m[c + "_o"] for c in phase]),
                  np.concatenate([m[c + "_n"] for c in phase]), "flux", mesh, scale=peak, wave=True,
                  keys=[f"{c} {k}" for c in phase for k in keys], note=f"% of the largest, {peak * 1e3:.2f} mWb")
        flux_cols = ["psi_d", "psi_q", *phase]
        rep.value(label, "flux-linkage retries in the run", int(a.psi_retries.sum()) if "psi_retries" in a else 0, int(b.psi_retries.sum()), "exact")
        rep.value(label, "states with a flux linkage of exactly 0", int(((a[flux_cols] == 0).any(axis=1)).sum()), int(((b[flux_cols] == 0).any(axis=1)).sum()), "exact")


def study_materials(rep: Report, old: Path, new: Path):
    """out/materials: the numbers the notes quote. Those that need no FEMM must be identical; those read from a solved field may move."""
    s = "materials: numbers that need no FEMM"
    mo, mn = Path(old) / "materials", Path(new) / "materials"
    if not (mo.is_dir() and mn.is_dir()):
        rep.missing(s, "the materials folder")
        return
    for name in ("current.json", "thermal.json", "structural.json", "magnetic.json", "fields.json"):
        a, b = mo / name, mn / name
        if not (a.is_file() and b.is_file()):
            rep.missing(s, name)
            continue
        found, count = compare_materials.compare_json(a, b, 0.0)
        if name in ("magnetic.json", "fields.json"):
            fem_stats(rep, name, a, b, found, count)
        else:
            rep.rows.append(Row(s, name, count, count - len(found), None, "identical" if not found else "DIFFERS", "exact", n=count,
                                n_flagged=len(found), note="values equal / values compared"))
    a, b = mo / "magnetic_m19_curve.csv", mn / "magnetic_m19_curve.csv"
    if a.is_file() and b.is_file():
        found, count = compare_materials.compare_csv(a, b, 0.0)
        rep.rows.append(Row(s, "magnetic_m19_curve.csv", count, count - len(found), None, "identical" if not found else "DIFFERS", "exact", n=count,
                            n_flagged=len(found), note="values equal / values compared"))


def flatten(obj, prefix=""):
    """{'a/b[2]': 1.5, ...} for every number in a parsed JSON value."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from flatten(v, f"{prefix}/{k}")
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            yield from flatten(v, f"{prefix}[{i}]")
    elif isinstance(obj, (int, float)) and not isinstance(obj, bool):
        yield prefix, float(obj)


def materials_kind(path: str) -> str:
    """The kind of a number in fields.json or magnetic.json, from its name. Statistics over many triangles are 'flux' (1 %);
    anything decided by the weakest few triangles, by a threshold on the tail of a distribution, or by one triangle's position
    is a 'tail' (10 %). The bridge zone is the triangles within 1 mm of the single most loaded triangle of the steel, so which
    triangles are in it depends on where the mesh happens to put that peak: a tail statistic too."""
    if "area_pct_below_B" in path or "weakest_element" in path or "B_parallel_p1_T" in path or "B_parallel_min_T" in path:
        return "tail"
    if path.startswith("/F4") or path.startswith("/bridge_zone_B_T"):
        return "tail"
    return "flux"


def fem_stats(rep: Report, name: str, a: Path, b: Path, found: list, count: int):
    """The numbers of a materials file that come from a solved field: which moved, and by how much."""
    s = f"materials: {name} (numbers read from solved fields)"
    old = dict(flatten(json.loads(a.read_text(encoding="utf-8"))))
    new = dict(flatten(json.loads(b.read_text(encoding="utf-8"))))
    paths = [p for p in old if p in new and old[p] != new[p]]
    rep.info(s, "values compared / values that moved", float(len(old)), float(len(paths)))
    for kind in ("flux", "tail"):
        sel = [p for p in paths if materials_kind(p) == kind]
        if not sel:
            continue
        pairs = [(old[p], new[p]) for p in sel]
        rep.group(s, f"{'statistics over many triangles (medians, percentiles)' if kind == 'flux' else 'tail statistics (weakest triangles, areas below a level)'}",
                  [x for x, _ in pairs], [y for _, y in pairs], kind, "fine", keys=sel,
                  note="only the values that moved are listed here; the rest are equal")
    # the headline numbers the notes quote
    for case in ("noload_with_band", "noload_no_band", "mtpa_with_band", "mtpa_no_band"):
        for q, kind in (("B_parallel_median_T", "flux"), ("B_parallel_p5_T", "flux"), ("B_parallel_p1_T", "tail"), ("B_parallel_min_T", "tail"),
                        ("steel_B_p99_T", "flux")):
            p = f"/{case}/{q}"
            if p in old and p in new:
                rep.value(s, f"{case}: {q}", old[p], new[p], kind, "fine")
        p = f"/{case}/weakest_element/n_magnet_elements"
        if p in old and p in new:
            rep.info(s, f"{case}: triangles in the magnets", old[p], new[p])


# ----------------------------------------------------------------------------------------------- 6. the page
def study_page(rep: Report, new: Path, page: Path | None):
    """pages/Rotor Flux Explorer.html: it must load nothing from the network, and its data must be FEMM 5.0's sweep."""
    s = "flux page"
    if page is None or not page.is_file():
        rep.missing(s, "pages/Rotor Flux Explorer.html")
        return
    text = page.read_text(encoding="utf-8")
    blob = re.compile(r'(<script type="application/octet-stream" id="g\d+">)[^<]*(</script>)')
    shell = blob.sub(r"\1\2", text)                     # the page without its 17 MB of data
    meta_match = re.search(r'<script type="application/json" id="meta">(.*?)</script>', shell, re.S)
    code = re.sub(r'<script type="application/json" id="meta">.*?</script>', "", shell, flags=re.S)
    rep.check(s, "src=\"http... (scripts, images, frames loaded from the network)", len(re.findall(r'\bsrc\s*=\s*["\']?(?:https?:)?//', code)), 0)
    rep.check(s, "<link ... href=\"http...\" (style sheets, fonts, icons from the network)",
              len([t for t in re.findall(r"<link\b[^>]*>", code, re.I) if re.search(r'href\s*=\s*["\']?(?:https?:)?//', t, re.I)]), 0)
    rep.check(s, "fetch( or XMLHttpRequest or import( or @import url(http", len(re.findall(r"\bfetch\s*\(|XMLHttpRequest|\bimport\s*\(|@import|url\(\s*[\"']?https?:", code)), 0)
    rep.check(s, "any \"http://\" or \"https://\" at all in the page's code", len(re.findall(r"https?://", code)), 0,
              note="(a link in the visible text would be listed here)" if re.search(r"https?://", code) else "")
    if not meta_match:
        rep.missing(s, "the page's data block (meta)")
        return
    meta = json.loads(meta_match.group(1))
    footer = " ".join(meta["footer"])
    rep.value(s, "the footer names FEMM 4.2", "no", "yes" if "FEMM 4.2" in footer else "no", "exact")
    rep.value(s, "the footer names FEMM 5.0", "yes", "yes" if "FEMM 5.0" in footer else "no", "exact", note=meta["footer"][0][:150])
    rep.rows.append(Row(s, "the footer's first sentence", None, meta["footer"][0], None, "info", "info"))
    page_data_vs_sweep(rep, s, new, text, meta)


def page_data_vs_sweep(rep: Report, s: str, new: Path, text: str, meta: dict):
    """The page's torque and flux-linkage tables are the sweep's states.csv, and its meshes are the sweep's meshes."""
    import base64
    import gzip
    import flux_explorer_sweep as sweep
    states = pd.read_csv(Path(new) / "flux_explorer/states.csv")
    p120 = pd.read_csv(Path(new) / "flux_explorer/states_p120.csv")
    n_solves = len(states) + len(p120)
    rep.value(s, "the footer's number of solves", f"{n_solves:,}", f"{n_solves:,}" if f"{n_solves:,} solves" in " ".join(meta["footer"]) else "other",
              "exact")
    gammas = [0, 15, 30, 45, 60, 75, 90]
    worst, n = 0.0, 0
    for gid in range(17):
        for bi, band in enumerate(sweep.BANDS):
            s_ = states[(states.gid == gid) & (states.band == band)]
            for name, col in (("T", "T_age"), ("psi_d", "psi_d"), ("psi_q", "psi_q")):
                table = np.array(meta["scal"][str(gid)][bi][name])
                for li, current in enumerate((0.0, *sweep.FIELD_I)):
                    for gi, g in enumerate(gammas):
                        sel = s_[(s_.I == current) & ((s_.gamma == g) if current > 0 else (s_.gamma == 0))].set_index("theta_e")[col]
                        want = np.array([sel[th] for th in sweep.POS])
                        worst = max(worst, float(np.abs(table[li, gi] - want).max()))
                        n += 1
    rep.check(s, f"largest |page table - states.csv| over {n * 4:,} torque and flux-linkage entries", worst, 5e-7,
              note="the page stores 6 decimals, so up to 5e-7 is rounding")
    blobs = dict(re.findall(r'<script type="application/octet-stream" id="g(\d+)">([^<]*)</script>', text))
    nodes_page, nodes_npz = [], []
    for gid in range(17):
        raw = gzip.decompress(base64.b64decode(blobs[str(gid)]))
        nodes_page.append(int(np.frombuffer(raw[:12], "<u4")[0]))
        nodes_npz.append(len(np.load(Path(new) / f"flux_explorer/g{gid:02d}_air.npz")["xy"]))
    rep.value(s, "mesh nodes of the 17 rotors: the page's, then the sweep's files", sum(nodes_npz), sum(nodes_page), "exact",
              note=f"the reference design: {nodes_page[0]:,} nodes")


# ----------------------------------------------------------------------------------------------- 7. what must hold on FEMM 5.0 itself
def study_self(rep: Report, old: Path, new: Path):
    """Checks that need no FEMM 4.2: the new numbers must agree with themselves."""
    s = "FEMM 5.0 on its own terms"
    for mesh_name, log_name in (("coarse", "check_femm"), ("fine", "check_femm_fine")):
        got = parse_check_femm(log_lines(find_log(new, log_name)))
        if len(got) == 4:
            worst = max(abs(v["T_age"] - v["T_wst"]) / abs(v["T_age"]) * 100 for k, v in got.items() if k != (0, 0))
            rep.check(s, f"|T_age - T_wst| at the three loaded points of check_femm, {mesh_name} mesh [% of the torque]", worst, 0.1)
            rep.check(s, f"L_d < L_q: torque at (-6, 8) A above torque at (+6, 8) A, {mesh_name} mesh [N.m difference]",
                      got[(-6, 8)]["T_age"] - got[(6, 8)]["T_age"], 0.0, higher_is_worse=False)
    results = {}
    for label, base in (("femm42", Path(old)), ("femm50", Path(new))):
        frames = []
        for name in ("barrier/phase1.csv", "barrier/phase2.csv", "barrier/reference_open.csv"):
            p = base / name
            if p.is_file():
                frames.append(pd.read_csv(p))
        if frames:
            allrows = pd.concat(frames, ignore_index=True)
            allrows = allrows[allrows.T_age.notna() & allrows.T_wst.notna()]
            scale = float(np.abs(allrows.T_age).max())
            diff = np.abs(allrows.T_age - allrows.T_wst)
            big = np.abs(allrows.T_age) >= 0.3 * scale
            mid = (np.abs(allrows.T_age) >= 0.1 * scale) & ~big
            results[label] = (len(allrows), 100 * float(diff.max() / scale), 100 * float((diff / np.abs(allrows.T_age))[big].max()),
                              100 * float((diff / np.abs(allrows.T_age))[mid].max()), float(diff.max()) * 1e3)
    if "femm50" in results:
        n, of_scale, big, mid, abs_mnm = results["femm50"]
        old_vals = results.get("femm42", (0, None, None, None, None))
        rep.rows.append(Row(s, f"T_age against T_wst over {n:,} solves: largest difference [% of the largest torque]", old_vals[1], of_scale, None,
                            "ok" if of_scale <= 0.1 else "FAIL", "must hold", 0.1, note=f"{abs_mnm:.2f} mN.m at most"))
        rep.rows.append(Row(s, "... and as a % of the torque itself, where it is at least 30 % of the largest", old_vals[2], big, None,
                            "ok" if big <= 0.1 else "FAIL", "must hold", 0.1))
        rep.rows.append(Row(s, "... and where it is between 10 % and 30 % of the largest (no limit: the two methods agree to about 1 mN.m, not to a share)",
                            old_vals[3], mid, None, "info", "info"))
    # Ld < Lq in the study
    sn = pd.read_csv(Path(new) / "barrier/summary.csv", index_col=0) if (Path(new) / "barrier/summary.csv").is_file() else None
    if sn is not None:
        ok = all(sn.loc[f"Ld_{i}A_mH", v] < sn.loc[f"Lq_{i}A_mH", v] for i in (5, 10, 15, 20) for v in ("with band", "no band"))
        smallest = min(sn.loc[f"saliency_{i}A", v] for i in (5, 10, 15, 20) for v in ("with band", "no band"))
        rep.check(s, "L_d < L_q at 5, 10, 15 and 20 A, with and without the band: the smallest saliency Lq/Ld", float(smallest), 1.0, higher_is_worse=False,
                  verdict="ok" if ok and smallest > 1 else "FAIL")
    # the slice against the full model
    vo = pd.read_csv(Path(old) / "validate/validate.csv") if (Path(old) / "validate/validate.csv").is_file() else None
    vn = pd.read_csv(Path(new) / "validate/validate.csv") if (Path(new) / "validate/validate.csv").is_file() else None
    if vn is not None:
        a, b = (slice_vs_full(vo) if vo is not None else {}), slice_vs_full(vn)
        worst_flux, worst_flux_old = max(b["flux, no load"], b["flux, 10 A"]), max(a.get("flux, no load", math.nan), a.get("flux, 10 A", math.nan))
        peak, peak_old = abs(b["peak torque, 10 A"]), (abs(a["peak torque, 10 A"]) if vo is not None else None)
        # "within 1.1 % / 1.5 %, or better": the plan's numbers are what FEMM 4.2's README rounded to, so when FEMM 4.2's own
        # figure is there (1.146 % and 1.479 %) FEMM 5.0 must be no worse than that; without it the plan's number is the limit.
        # The limit_pct column shows the limit that the verdict used.
        limit_flux = worst_flux_old if vo is not None else 1.1
        limit_peak = peak_old if vo is not None else 1.5
        rep.rows.append(Row(s, "the slice against the full model: worst |flux linkage difference| at any rotor position [% of the full peak]",
                            worst_flux_old if vo is not None else None, worst_flux, None, "ok" if worst_flux <= limit_flux + 1e-9 else "FAIL",
                            "must hold", limit_flux, note="the plan: within 1.1 % (FEMM 4.2: the README's 1.1 %), or better"))
        rep.rows.append(Row(s, "the slice against the full model: peak torque at 10 A [% difference]", peak_old, peak, None,
                            "ok" if peak <= limit_peak + 1e-9 else "FAIL", "must hold", limit_peak, note="the plan: within 1.5 %, or better"))
    # no flux linkage came back as exactly zero, anywhere. One row of these tables is one solve, so the number of solves in
    # the label is counted from the tables that were there (1,070 + 102 + 2,040 + 306 = 3,518 when all five are).
    zeros, solves = 0, 0
    for name in ("barrier/phase1.csv", "barrier/phase2.csv", "barrier/reference_open.csv", "flux_explorer/states.csv", "flux_explorer/states_p120.csv"):
        p = Path(new) / name
        if p.is_file():
            d = pd.read_csv(p)
            solves += len(d)
            zeros += int((d[["psi_a", "psi_b", "psi_c"]] == 0).any(axis=1).sum()) + (int(d.psi_retries.sum()) if "psi_retries" in d else 0)
    rep.check(s, f"flux linkages of exactly 0, plus read-back retries, in {solves:,} solves", float(zeros), 0.0)
    # the same scripts run twice print the same numbers
    pairs = same_logs(Path(new))
    if pairs:
        rep.check(s, f"logs of two complete runs that differ in more than their times ({pairs[0]} pairs)", float(pairs[1]), 0.0,
                  note="run_all.py keeps the last log of each step as NN-name.log.prev")


TIMES = re.compile(r"\d+(?:\.\d+)?\s?(?:ms|s/solve|s\b|x\b)|PIDs: \[[^\]]*\]|PID \d+")    # what two runs of one script may differ in


def same_logs(new: Path):
    """(pairs of logs from two complete runs, how many differ once the times are blanked out), or None if there are no pairs."""
    folder = Path(new) / "validation" / "logs"
    pairs = sorted(folder.glob("*.log.prev"))
    if not pairs:
        return None
    different = 0
    for prev in pairs:
        cur = prev.with_suffix("")
        if not cur.is_file():
            continue
        a = [TIMES.sub("<t>", ln) for ln in log_lines(prev)]
        b = [TIMES.sub("<t>", ln) for ln in log_lines(cur)]
        different += a != b
    return len(pairs), different


# ----------------------------------------------------------------------------------------------- 8. speed
def parse_old_times(old: Path) -> dict:
    """Seconds that FEMM 4.2's runs wrote into their own logs (the old out folder holds a few of them)."""
    t = {}
    def grab(name, pattern):
        p = Path(old) / name
        if p.is_file():
            m = re.search(pattern, read_text(p))
            return m
    m = grab("barrier_study.log", r"phase 1: (\d+) solves in (\d+) s")
    n = grab("barrier_study.log", r"phase 2: (\d+) solves in (\d+) s")
    if m and n:
        t["barrier_study"] = (float(m.group(2)) + float(n.group(2)), int(m.group(1)) + int(n.group(1)))
    for key, name, pat in (("reference_open", "reference_open.log", r"(\d+) solves in (\d+) s"), ("validate_full", "validate_full.log", r"(\d+) solves in (\d+) s"),
                           ("flux_explorer_sweep", "flux_explorer_sweep.log", r"(\d+) solves in (\d+) s"),
                           ("flux_explorer_sweep120", "flux_explorer_sweep120.log", r"(\d+) solves in (\d+) s")):
        m = grab(name, pat)
        if m:
            t[key] = (float(m.group(2)), int(m.group(1)))
    m = grab("optimize_demo.log", r"(\d+) generations in (\d+) s")
    if m:
        t["optimize_demo"] = (float(m.group(2)), None)
    p = Path(old) / "worker_scaling.log"
    if p.is_file():
        seconds = [float(x) for x in re.findall(r"workers:\s+([\d.]+) s", read_text(p))]
        if seconds:
            t["worker_scaling"] = (sum(seconds), 6 * 48)
    return t


def study_speed(rep: Report, old: Path, new: Path):
    """How long each study took: FEMM 4.2's own logs against run_all.csv."""
    s = "speed (seconds)"
    table = Path(new) / "validation" / "run_all.csv"
    if not table.is_file():
        rep.missing(s, "validation/run_all.csv")
        return
    runs = pd.read_csv(table).set_index("name")
    old_t = parse_old_times(old)
    total_old = total_new = 0.0
    for name, (seconds, solves) in old_t.items():
        if name in runs.index:
            now = float(runs.loc[name, "seconds"])
            rep.info(s, f"{name}" + (f" ({solves:,} solves)" if solves else ""), seconds, now, note=f"{seconds / now:.1f} times faster")
            total_old += seconds
            total_new += now
    rep.info(s, "the studies above, together", total_old, total_new, note=f"{total_old / total_new:.1f} times faster")
    rep.info(s, "the whole sweep of 33 steps (run_all.csv, FEMM 5.0; FEMM 4.2: about 2 h 13 min by the README's step times)", 8000.0,
             float(runs.seconds.sum()), note="FEMM 4.2: the README's step times added up, approximate")
    if "femm_left_running" in runs:
        rep.check(s, "FEMM programs left running after any step (run_all.csv)", float(runs.femm_left_running.sum()), 0.0)
    rep.check(s, "steps that did not finish ok (run_all.csv)", float((runs.status != "ok").sum()), 0.0)


# ----------------------------------------------------------------------------------------------- 9. the command
def compare_all(old: Path, new: Path, page: Path | None = None) -> Report:
    rep = Report()
    study_not_femm(rep, old, new)
    study_mesh(rep, old, new)
    study_check_femm(rep, old, new)
    study_printed_only(rep, old, new)
    study_barrier(rep, old, new)
    study_reference_open(rep, old, new)
    study_validate(rep, old, new)
    study_scaling(rep, old, new)
    study_optimize(rep, old, new)
    study_flux_explorer(rep, old, new)
    study_materials(rep, old, new)
    study_page(rep, new, page if page is not None else P / "pages" / "Rotor Flux Explorer.html")
    study_self(rep, old, new)
    study_speed(rep, old, new)
    return rep


# ----------------------------------------------------------------------------------------------- 10. the audit of the two extra rules
@contextlib.contextmanager
def rules(floor: float | None = None, wave: bool | None = None):
    """Run a block with the floor and/or the waveform rule changed, and put them back afterwards (even after an error)."""
    global FLOOR, WAVE_RULE
    saved = FLOOR, WAVE_RULE
    FLOOR = saved[0] if floor is None else floor
    WAVE_RULE = saved[1] if wave is None else wave
    try:
        yield
    finally:
        FLOOR, WAVE_RULE = saved


def judged_size(r: Row) -> float | None:
    """The percentage a row's verdict rests on: |delta| of a single number, the largest |delta| of a group."""
    if r.stats is not None:
        return r.stats[2]
    return None if r.delta_pct is None else abs(r.delta_pct)


NEAR = 4.0      # audit(): a row whose judged percentage is within this share (percent) of its limit is called close to it
FLOORS_TRIED = (0.09, 0.12, 0.15)      # audit(): other values of the floor, to see how much the FLAGs depend on 0.10


def audit(old: Path, new: Path, page: Path | None = None, out=print) -> dict:
    """Which verdicts rest on the two rules of this script that the plan does not have: the floor, and the waveform's own peak?

    The comparison runs three times and the rows are matched by their place in it (it is the same code and the same files, so
    the rows come in the same order): as it is; with FLOOR = 0, so a small FEMM 4.2 value is judged by its own plain percentage;
    and with FLOOR = 0 and the waveform rule off, so every sample of a waveform is judged by its own plain percentage too.
    Returns, and prints, three lists (and, below, the floor's sensitivity):
      floor      rows that are ok as they are and FLAG without the floor: they pass ONLY through it. Each carries its judged
                 percentage, its limit, how far below the limit it is (margin, in percent of the limit) and whether that is
                 within NEAR percent (close). A reader who doubts the floor should look at these rows first.
      more       rows that are FLAG either way, where the floor hides some of the members that are beyond the limit.
      wave       rows that are ok without the floor and FLAG once the waveform rule is off as well: they pass only because
                 a sample near a zero crossing is judged as a share of the waveform's peak, not of its own tiny value.
    It also runs the comparison with the floor a little lower and a little higher (FLOORS_TRIED) and names the rows that then
    become FLAG or stop being FLAG ("floors" in the result): how much the 19 flags depend on the value 0.10.
    """
    base = compare_all(old, new, page)
    with rules(floor=0.0):
        plain_floor = compare_all(old, new, page)
    with rules(floor=0.0, wave=False):
        plain_wave = compare_all(old, new, page)
    assert [(r.study, r.quantity) for r in base.rows] == [(r.study, r.quantity) for r in plain_floor.rows] == \
        [(r.study, r.quantity) for r in plain_wave.rows], "the three runs must produce the same rows in the same order"
    result = {"floor": [], "more": [], "wave": []}
    for r0, r1, r2 in zip(base.rows, plain_floor.rows, plain_wave.rows):
        size, plain = judged_size(r0), judged_size(r1)
        if r0.verdict == "ok" and r1.verdict == "FLAG":
            margin = 100 * (r0.limit_pct - size) / r0.limit_pct
            result["floor"].append({"study": r0.study, "quantity": r0.quantity, "judged": size, "limit": r0.limit_pct, "margin": margin,
                                    "close": margin < NEAR, "plain": plain, "n_beyond": r0.n_flagged, "n_beyond_plain": r1.n_flagged})
        elif r0.verdict == "FLAG" and r1.verdict == "FLAG" and r1.n_flagged > r0.n_flagged:
            result["more"].append({"study": r0.study, "quantity": r0.quantity, "n": r0.n, "n_beyond": r0.n_flagged, "n_beyond_plain": r1.n_flagged})
        if r1.verdict == "ok" and r2.verdict == "FLAG":
            result["wave"].append({"study": r0.study, "quantity": r0.quantity, "judged": judged_size(r1), "limit": r1.limit_pct,
                                   "plain": judged_size(r2), "n_beyond": r1.n_flagged, "n_beyond_plain": r2.n_flagged})
    out(f"\naudit: the verdicts that rest on the two rules that are not the plan's (FLOOR = {FLOOR}, waveform rule on)")
    out(f"\n{len(result['floor'])} rows pass only through the floor (ok as they are, FLAG with FLOOR = 0); "
        f"{sum(x['close'] for x in result['floor'])} of them are within {NEAR:g} % of their limit (marked <<):")
    for x in result["floor"]:
        out(f"  {'<<' if x['close'] else '  '} {x['study']} | {x['quantity'][:70]}: judged {x['judged']:.3f} % of a limit of {x['limit']:g} %"
            f" (margin {x['margin']:.1f} % of the limit); the plain percentage is {x['plain']:.1f} %"
            + (f"; members beyond the limit by the plain rule: {x['n_beyond_plain']}" if x["n_beyond_plain"] else ""))
    out(f"\n{len(result['more'])} rows are FLAG either way, and the floor hides some of their members:")
    for x in result["more"]:
        out(f"     {x['study']} | {x['quantity'][:70]}: {x['n_beyond']} of {x['n']} members beyond the limit with the floor, {x['n_beyond_plain']} without")
    out(f"\n{len(result['wave'])} further rows pass only through the waveform rule (ok with FLOOR = 0, FLAG once each sample is judged by its own value):")
    for x in result["wave"]:
        out(f"     {x['study']} | {x['quantity'][:70]}: judged {x['judged']:.3f} % of a limit of {x['limit']:g} %; "
            f"by its own value the largest is {x['plain']:.1f} %, {x['n_beyond_plain']} members beyond")
    result["floors"] = {}
    n_flag = sum(r.verdict == "FLAG" for r in base.rows)
    out(f"\nhow much the FLAGs depend on the value of the floor (now {FLOOR}, {n_flag} FLAG rows):")
    for f in FLOORS_TRIED:
        with rules(floor=f):
            other = compare_all(old, new, page)
        became = [f"{a.study} | {a.quantity}" for a, b in zip(base.rows, other.rows) if a.verdict != "FLAG" and b.verdict == "FLAG"]
        stopped = [f"{a.study} | {a.quantity}" for a, b in zip(base.rows, other.rows) if a.verdict == "FLAG" and b.verdict != "FLAG"]
        result["floors"][f] = {"flag": sum(r.verdict == "FLAG" for r in other.rows), "became": became, "stopped": stopped}
        out(f"  FLOOR = {f}: {result['floors'][f]['flag']} FLAG rows; {len(became)} become FLAG, {len(stopped)} stop being FLAG")
        for name in became:
            out(f"       becomes FLAG: {name[:150]}")
        for name in stopped:
            out(f"       stops being FLAG: {name[:150]}")
    return result


def write_csv(path: Path, header: list, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="FEMM 4.2 against FEMM 5.0: every number, with a verdict.")
    ap.add_argument("old", help="FEMM 4.2's out folder (read only)")
    ap.add_argument("new", nargs="?", default="out", help="FEMM 5.0's out folder (default: out)")
    ap.add_argument("--csv", default=str(HERE / "femm5_vs_femm42.csv"), help="where the table goes")
    ap.add_argument("--detail", default=None, help="where the per-state table goes (default: NEW_OUT/validation/femm5_vs_femm42_detail.csv)")
    ap.add_argument("--page", default=None, help="the flux page to check (default: pages/Rotor Flux Explorer.html in P)")
    ap.add_argument("--quiet", action="store_true", help="print the totals only")
    ap.add_argument("--audit", action="store_true", help="after the comparison, list the verdicts that rest on the floor and on the waveform rule")
    args = ap.parse_args(argv)
    old, new = Path(args.old), Path(args.new)
    for label, folder in (("OLD_OUT", old), ("NEW_OUT", new)):
        if not folder.is_dir():
            print(f"{label} is not a folder: {folder}", file=sys.stderr)
            return 2
    rep = compare_all(old, new, Path(args.page) if args.page else None)
    if not args.quiet:
        print_report(rep)
    write_csv(Path(args.csv), CSV_HEADER, rep.csv_rows())
    detail = Path(args.detail) if args.detail else new / "validation" / "femm5_vs_femm42_detail.csv"
    write_csv(detail, DETAIL_HEADER, ([d[h] if not isinstance(d[h], (float, np.floating)) else _num(d[h]) for h in DETAIL_HEADER] for d in rep.detail))
    counts = {}
    for r in rep.rows:
        counts[r.verdict] = counts.get(r.verdict, 0) + 1
    for verdict, header in (("FLAG", "FLAGGED: beyond the tolerance, each one needs an explanation (the report, section 2.12)"),
                            ("FAIL", "FAILED: a check that must hold on FEMM 5.0 itself"), ("MISSING", "MISSING: nothing to compare"),
                            ("DIFFERS", "DIFFERS: should have been equal")):
        rows = [r for r in rep.rows if r.verdict == verdict]
        if rows and not args.quiet:
            print(f"\n{header}")
            for r in rows:
                print(f"  {r.study} | {r.quantity}" + (f": max {r.stats[2]:.3g} %, {r.n_flagged} of {r.n}" if r.stats else
                                                       (f": {r.delta_pct:+.3g} %" if r.delta_pct is not None else "")))
    print(f"\n{len(rep.rows)} compared quantities or groups: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
    print(f"wrote {args.csv}\nwrote {detail}  ({len(rep.detail):,} members)")
    if args.audit:
        audit(old, new, Path(args.page) if args.page else None)
    # Exit status: 1 if a check that must hold on FEMM 5.0 failed or something to compare is missing. FLAG and DIFFERS rows do not
    # fail the run: they are listed above for explanation, and the report (section 2.12) explains each of them.
    return 1 if any(r.verdict in ("FAIL", "MISSING") for r in rep.rows) else 0


if __name__ == "__main__":
    sys.exit(main())
