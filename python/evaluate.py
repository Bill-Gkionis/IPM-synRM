"""evaluate(x): one optimizer candidate in, one result row out (the unit that optimize_demo.py repeats).

A candidate is a dict of Params overrides. It is drawn and checked without FEMM first; if it can be built and no rib
is too thin, one FEMM of its own (own port, below-normal priority, see femm_session.py) solves 15 rotor positions at
i_q = 10 A, and the row holds mean torque, ripple, magnet mass and the q-rib width. Every candidate gets private
files in runs/, so workers never share a .fem or .ans.
Run: uv run python evaluate.py   (from this folder; 5 candidates on 4 workers, 6 s with FEMM 5.0)
Expect with FEMM 5.0 (coarse mesh, 2026-11-04): the baseline {} gives T_mean 0.7289 and T_ripple 0.1579 N.m (the same
as example.py: 15 positions or 30 over the period give the same mean), mag_y 5.0 gives T_mean 1.0013, web 1.2 gives 0.6635,
and the last two are rejected without FEMM ("magnet slot reaches the bridge", "q_rib = 0.00 < 1.0; web = 0.45 < 1.0").
"""
import os, uuid
from concurrent.futures import ProcessPoolExecutor
import femm_session
import motorgeo as mg

def evaluate(x: dict) -> dict:
    """One optimizer candidate -> one CSV row. Runs inside a worker process."""
    p = mg.Params(**x)
    try:
        regions = mg.build(p)
    except mg.GeometryError as e:
        return {**x, "feasible": False, "reason": str(e)}
    metrics, problems = mg.check(p, regions)
    if problems:
        return {**x, "feasible": False, "reason": "; ".join(problems)}

    import femm
    fem = os.path.abspath(f"runs/{uuid.uuid4().hex}.fem")      # private files per candidate
    femm_session.open_femm()                                   # this worker's own FEMM
    try:
        mg.to_femm(p, regions, fem)
        femm.mi_close()
        model = mg.FemmModel(p, fem)
        T = [model.solve(th, 0, 10)["T_age"] for th in range(0, 60, 4)]  # 15 positions, one ripple period
    finally:
        femm_session.close_femm()
    return {**x, "feasible": True, "T_mean": sum(T) / len(T), "T_ripple": max(T) - min(T),
            "magnet_g": metrics["magnet_mass_g"], "q_rib": metrics["q_rib"]}

if __name__ == "__main__":
    os.makedirs("runs", exist_ok=True)
    candidates = [{}, {"mag_y": 5.0}, {"web": 1.2}, {"v_half_deg": 55}, {"mag_x": 38.4}]
    with ProcessPoolExecutor(4) as pool:
        for row in pool.map(evaluate, candidates):
            print(row)
