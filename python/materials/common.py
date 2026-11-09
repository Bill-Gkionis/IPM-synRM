"""common.py: what every materials script shares.

The materials notes (IPM-synRM/docs/1-materials/) quote numbers. Every one of them is either
computed by a script in this folder or read from a source named next to it here. This module
holds the three things the scripts share:

1. the physical constants (CODATA 2018, the values SI fixed in 2019);
2. readers for the material data that FEMM actually uses: the block properties inside our
   reference model `IPM synRM analysis.FEM`, and FEMM's own library `matlib.dat`;
3. one plot style, so every figure in the notes looks like part of one set, and a `save()` that
   writes each figure next to the scripts' outputs, in out/materials/, and copies it anywhere else only when asked, through the
   environment variable MOTOR_FIGURES_DIR (see save()).

Run any script from the `Python to FEMM` folder (the one above this one), for example:
    uv run python materials/magnetic.py                   # the figures land in out/materials/ and nowhere else
    $env:MOTOR_FIGURES_DIR = "out/materials/figures"      # (PowerShell) also copy every figure into that folder
    $env:MOTOR_FIGURES_DIR = "docs"                       # also copy every figure into the repository's docs figures folder (config.toml, docs_figures)
Unset or empty, nothing is copied into the docs figures folder (the rule: no run writes there; figures are copied by hand). The word `docs` is
the one way to have a run write there, on purpose, when the repository's figures are being regenerated. The numbers (JSON, CSV) always land in out/materials/.
"""
from __future__ import annotations

import json
import math
import os
import re
import shutil
import sys
from pathlib import Path

import numpy as np

# ----------------------------------------------------------------------------- 1. where things are
PROC = Path(__file__).resolve().parents[1]                 # Storage\Python to FEMM
if str(PROC) not in sys.path:                               # `python materials/x.py` starts with materials\ first on sys.path, not P,
    sys.path.insert(0, str(PROC))                           # and paths.py (like motorgeo.py) is in P
import paths  # noqa: E402  (config.toml: the one place that names this PC's folders and programs)

OUT = PROC / "out" / "materials"                            # numbers (JSON, CSV) and figures land here
OUT.mkdir(parents=True, exist_ok=True)
FIGURES = paths.path("docs_figures")                        # the repository's docs figures folder (save() copies figures here only when MOTOR_FIGURES_DIR says `docs`)
DOCS = "docs"                                               # the one word MOTOR_FIGURES_DIR may hold instead of a folder: "copy into FIGURES, the repository's docs figures folder"
REFERENCE_FEM = paths.path("reference_fem")                 # our reference model, the copy in P\reference (read_fem_blocks() reads it)
MATLIB = paths.path("femm_matlib")                          # FEMM 5.0's material library (Air, N45 and the ASM "M-19" as in 4.2)
M19_STEEL_FOLDER = "Soft Magnetic Materials / Silicon Iron"   # where FEMM 4.2's library kept the block "M-19 Steel"

# ----------------------------------------------------------------------------- 2. constants (CODATA 2018)
# Source: NIST, "The NIST Reference on Constants, Units, and Uncertainty", https://physics.nist.gov/cuu/Constants/
# Since 2019 e, h, k_B and N_A are exact by definition; mu0 is measured (it was exactly 4*pi*1e-7 before).
E = 1.602176634e-19          # elementary charge [C], exact
ME = 9.1093837015e-31        # electron mass [kg]
HBAR = 1.054571817e-34       # reduced Planck constant [J s]  (h / 2 pi, h exact)
KB = 1.380649e-23            # Boltzmann constant [J/K], exact
NA = 6.02214076e23           # Avogadro constant [1/mol], exact
MU0 = 1.25663706212e-6       # vacuum permeability [N/A^2]; 4*pi*1e-7 = 1.2566370614e-6 differs by 5e-10 relative
MU_B = 9.2740100783e-24      # Bohr magneton e*hbar/(2 me) [A m^2] = [J/T]
A0 = 5.29177210903e-11       # Bohr radius [m]
ALPHA = 7.2973525693e-3      # fine-structure constant [-]
C = 299792458.0              # speed of light [m/s], exact
LORENZ0 = math.pi ** 2 / 3 * (KB / E) ** 2   # Sommerfeld's Lorenz number = 2.443e-8 W Ohm/K^2

# ----------------------------------------------------------------------------- 3. what FEMM uses
def read_fem_blocks(path: Path = REFERENCE_FEM) -> dict[str, dict]:
    """Return every material block of a FEMM .fem file: {name: {property: value, ..., 'bh': array}}.

    A .fem file is plain text. Its [BlockProps] section lists one <BeginBlock> ... <EndBlock> per
    material, with lines like `<Mu_x> = 1.05`, and for nonlinear materials `<BHPoints> = 47`
    followed by 47 lines of `B  H` (tesla, A/m).
    """
    text = Path(path).read_text(encoding="latin-1")
    blocks = {}
    for body in re.findall(r"<BeginBlock>(.*?)<EndBlock>", text, flags=re.S | re.I):
        props, bh = {}, []
        for line in body.strip().splitlines():
            line = line.strip()
            m = re.match(r"<(\w+)>\s*=\s*(.*)", line)
            if m:
                key, val = m.group(1), m.group(2).strip()
                props[key] = val.strip('"') if val.startswith('"') else float(val)
            elif line:                                     # a "B  H" pair after <BHPoints>
                bh.append([float(v) for v in line.split()])
        props["bh"] = np.array(bh) if bh else np.empty((0, 2))
        blocks[props["BlockName"]] = props
    return blocks


def read_matlib(path: Path = MATLIB) -> list[dict]:
    """Every block of FEMM's material library, each tagged with the folder it sits in."""
    text = Path(path).read_text(encoding="latin-1")
    out, folders = [], []
    for token in re.finditer(r"<BeginFolder>|<EndFolder>|<FolderName>\s*=\s*\"([^\"]*)\"|<BeginBlock>(.*?)<EndBlock>",
                             text, flags=re.S | re.I):
        t = token.group(0)
        if t.lower() == "<beginfolder>":
            folders.append("?")
        elif t.lower() == "<endfolder>":
            folders.pop()
        elif token.group(1) is not None:
            folders[-1] = token.group(1)
        else:
            block = read_fem_blocks_text(t)
            block["folder"] = " / ".join(folders)
            out.append(block)
    return out


def m19_steel_block() -> dict:
    """FEMM 4.2's library block "M-19 Steel", in the shape read_matlib() returns for a block.

    FEMM 5.0's library has no "M-19 Steel" any more (its "M-19 14 mil" is a different curve), but the motor model
    still uses FEMM 4.2's, which motorgeo carries as M19_STEEL_FEMM42. This turns that dictionary into what
    read_matlib() would have returned from FEMM 4.2's matlib.dat: the same keys (BlockName, Mu_x, ..., BHPoints) as
    numbers, "bh" as a (47, 2) array of (B [T], H [A/m]) and "folder" as the library folder FEMM 4.2 kept it in.
    """
    if str(PROC) not in sys.path:
        sys.path.insert(0, str(PROC))
    import motorgeo
    constant = motorgeo.M19_STEEL_FEMM42
    block = {k: v if isinstance(v, str) else float(v) for k, v in constant.items() if k != "BH"}
    block["bh"] = np.array(constant["BH"], dtype=float)
    block["folder"] = M19_STEEL_FOLDER
    return block


def read_fem_blocks_text(text: str) -> dict:
    """Same as read_fem_blocks, for a single <BeginBlock>...<EndBlock> string."""
    tmp = {}
    for body in re.findall(r"<BeginBlock>(.*?)<EndBlock>", text, flags=re.S | re.I):
        props, bh = {}, []
        for line in body.strip().splitlines():
            line = line.strip()
            m = re.match(r"<(\w+)>\s*=\s*(.*)", line)
            if m:
                key, val = m.group(1), m.group(2).strip()
                props[key] = val.strip('"') if val.startswith('"') else float(val)
            elif line:
                bh.append([float(v) for v in line.split()])
        props["bh"] = np.array(bh) if bh else np.empty((0, 2))
        tmp = props
    return tmp


# ----------------------------------------------------------------------------- 4. the plot style
# The validated reference palette of a data-visualization method (light surface). Categorical slots are used
# in this fixed order and never cycled; ordered series (temperatures) use one blue ramp, light to dark.
SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE_RAMP = {250: "#86b6ef", 300: "#6da7ec", 350: "#5598e7", 400: "#3987e5", 450: "#2a78d6",
             500: "#256abf", 550: "#1c5cab", 600: "#184f95", 650: "#104281", 700: "#0d366b"}
CRITICAL = "#d03b3b"         # status colour: reserved for limits that break something (always with a label)


def style():
    """Apply the shared matplotlib style and return pyplot."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": ["Segoe UI", "DejaVu Sans"], "font.size": 10,
        "text.color": INK, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.titlesize": 11, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "axes.axisbelow": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "lines.linewidth": 1.8, "lines.markersize": 5.5, "legend.frameon": False, "legend.fontsize": 9,
        "axes.prop_cycle": matplotlib.cycler(color=CAT),
    })
    return plt


def figures_copy_destination() -> Path | None:
    """Where save() copies a figure besides out/materials/, or None: nowhere. The environment variable MOTOR_FIGURES_DIR decides, read each time:
    unset or empty, None (a run writes nothing into the docs figures folder); a folder, that folder; the word `docs` (in any case), FIGURES, the repository's docs figures folder.
    A folder that is really called docs is written ./docs."""
    asked = os.environ.get("MOTOR_FIGURES_DIR")
    if not asked:
        return None
    return FIGURES if asked.strip().lower() == DOCS else Path(asked)


def save(fig, name: str):
    """Write a figure as PNG to out/materials/ and, only if asked, copy it to a second folder.

    The copy is opt-in (the rule: no run writes into the repository's docs figures folder; figures are copied by hand). MOTOR_FIGURES_DIR unset or empty:
    no copy. MOTOR_FIGURES_DIR a folder: the figure is copied there too, and the folder is made if it is missing. MOTOR_FIGURES_DIR=docs:
    the copy goes to FIGURES, the repository's docs figures folder (config.toml, docs_figures). The variable is read here, each time, so it counts however
    late it was set. Every figure of every materials script goes through this function: it is the only place that copies into the docs figures folder.

    MOTOR_FIGURES_DIR may be the numbers' folder itself (out/materials): the figure is already there, so it is left
    alone. (Copying a file onto itself fails: on Windows with PermissionError, "being used by another process",
    because Python 3.13's shutil.copy2 hands the two names straight to the operating system; elsewhere with
    shutil.SameFileError.) Both paths are resolved before they are compared, so it does not matter how the folder is
    spelled, or whether out is a junction.
    """
    path = OUT / name
    fig.savefig(path, dpi=170, bbox_inches="tight")
    destination = figures_copy_destination()
    if destination is not None:
        destination.mkdir(parents=True, exist_ok=True)
        if (destination / name).resolve() != path.resolve():
            shutil.copy2(path, destination / name)
    print(f"  figure: {name}")
    return path


def write_numbers(name: str, numbers: dict):
    """Save the numbers a note quotes as JSON (rounded for reading, not for computing)."""
    def clean(v):
        if isinstance(v, (np.floating, float)):
            return float(f"{float(v):.6g}")
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, dict):
            return {k: clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple, np.ndarray)):
            return [clean(x) for x in v]
        return v
    path = OUT / f"{name}.json"
    path.write_text(json.dumps(clean(numbers), indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"  numbers: {path.name}")
    return path


def label_end(ax, x, y, text, color=INK2, dx=4, dy=0, **kw):
    """A direct label at the end of a line (text in ink, never in the series colour)."""
    ax.annotate(text, (x, y), xytext=(dx, dy), textcoords="offset points", va="center",
                color=color, fontsize=9, **kw)
