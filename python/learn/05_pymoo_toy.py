"""NSGA-II on a toy problem with two goals that disagree. No FEMM.
Run, from the learn folder: uv run python 05_pymoo_toy.py

Cook 1 likes recipes near x = +1, cook 2 near x = -1, and both dislike large y. No recipe pleases both,
so the answer is not one point but a front of best compromises: the line from (-1, 0) to (+1, 0).
"""
import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize


class TwoCooks(Problem):
    def __init__(self):
        super().__init__(n_var=2, n_obj=2, n_ieq_constr=1, xl=[-2, -2], xu=[2, 2])

    def _evaluate(self, X, out, *args, **kwargs):       # X: one row per candidate in this generation
        x, y = X[:, 0], X[:, 1]
        out["F"] = np.column_stack([(x - 1) ** 2 + y ** 2, (x + 1) ** 2 + y ** 2])   # both to be minimised
        out["G"] = np.column_stack([x - 0.5])            # a rule: x <= 0.5 (G <= 0 means allowed)


res = minimize(TwoCooks(), NSGA2(pop_size=40), ("n_gen", 150), seed=1, verbose=False)
order = np.argsort(res.X[:, 0])
print(f"{len(res.X)} best compromises found. x from {res.X[order[0], 0]:.2f} to {res.X[order[-1], 0]:.2f}, "
      f"largest |y| {np.abs(res.X[:, 1]).max():.3f}")
for i in order[::8]:
    print(f"  x = {res.X[i, 0]:+.2f}, y = {res.X[i, 1]:+.3f}  ->  cook 1 unhappiness {res.F[i, 0]:.2f}, cook 2 {res.F[i, 1]:.2f}")
