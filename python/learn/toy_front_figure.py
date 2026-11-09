"""Figure for the project note on the optimizer: what NSGA-II returns on the two-cooks toy problem.

Run, from the learn folder: uv run python toy_front_figure.py <output.png>   (the output folder must exist)
"""
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.callback import Callback
from pymoo.core.problem import Problem
from pymoo.optimize import minimize

INK, INK2, MUTED, GRID, AXIS, SURF = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"
BLUE, ORANGE = "#2a78d6", "#eb6834"


class TwoCooks(Problem):
    def __init__(self):
        super().__init__(n_var=2, n_obj=2, n_ieq_constr=1, xl=[-2, -2], xu=[2, 2])

    def _evaluate(self, X, out, *args, **kwargs):
        x, y = X[:, 0], X[:, 1]
        out["F"] = np.column_stack([(x - 1) ** 2 + y ** 2, (x + 1) ** 2 + y ** 2])
        out["G"] = np.column_stack([x - 0.5])


class Keep(Callback):                   # remember the first generation, to show where the search started
    def notify(self, algorithm):
        if algorithm.n_gen == 1:
            self.first = algorithm.pop.get("X").copy(), algorithm.pop.get("F").copy()


cb = Keep()
res = minimize(TwoCooks(), NSGA2(pop_size=40), ("n_gen", 150), seed=1, callback=cb, verbose=False)
X0, F0 = cb.first

fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=SURF)
for ax in axs:
    ax.set_facecolor(SURF); ax.grid(color=GRID, lw=0.6); ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(AXIS)
    ax.tick_params(colors=MUTED)
ax = axs[0]
ax.axvspan(0.5, 2, color="#f0efec", zorder=0)
ax.text(1.25, 1.75, "not allowed\n(x > 0.5)", ha="center", va="top", fontsize=8.5, color=INK2)
ax.plot(X0[:, 0], X0[:, 1], "o", ms=4.5, color=AXIS, label="generation 1 (random)")
ax.plot(res.X[:, 0], res.X[:, 1], "o", ms=5, color=BLUE, label="generation 150: the front")
for x, name in ((1, "cook 1's favourite"), (-1, "cook 2's favourite")):
    ax.plot(x, 0, "*", ms=13, color=ORANGE, mec=SURF)
    ax.annotate(name, (x, 0), textcoords="offset points", xytext=(0, -18), ha="center", fontsize=8.5, color=INK)
ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_aspect("equal")
ax.set_xlabel("recipe setting x", color=INK2); ax.set_ylabel("recipe setting y", color=INK2)
ax.set_title("The recipes (what the optimizer changes)", loc="left", fontsize=10, color=INK)
ax.legend(frameon=False, loc="upper left", fontsize=8)
ax = axs[1]
ax.plot(F0[:, 0], F0[:, 1], "o", ms=4.5, color=AXIS, label="generation 1")
ax.plot(res.F[:, 0], res.F[:, 1], "o", ms=5, color=BLUE, label="generation 150")
ax.set_xlim(0, 12); ax.set_ylim(0, 12)
ax.set_xlabel("cook 1's unhappiness (lower is better)", color=INK2)
ax.set_ylabel("cook 2's unhappiness (lower is better)", color=INK2)
ax.set_title("The scores: no point beats a blue one on both", loc="left", fontsize=10, color=INK)
ax.legend(frameon=False, loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(sys.argv[1], dpi=150, facecolor=SURF)
print("wrote", sys.argv[1])
