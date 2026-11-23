"""Rioptimizimi i fundit te tures.

Per nje ture te mire mbahen fikse hyrja dhe m planetet e para, dhe kerkohen me dendur
(me shume fillime, me shume zgjidhje per dritare) kembet e mbetura me nje beam te vogel.
"""
import numpy as np
import beam_search as BS
from heuristic import states_of


def suffix_task(args):
    """Nje detyre: tura x, prerja pas planetit m, beam Bm. Kthen turat e reja (N x 34)."""
    x, m, Bm, nstart, keep, seed = args
    x = np.asarray(x, float)
    p, t, vin, dv = states_of(x)[m]
    s0 = dict(seq=list(x[27:27 + m + 1].astype(int)), t=t, dv=dv, vin=vin, x=list(x[:3 + 4 * m]))
    cur = [s0]
    for lev in range(6 - m):
        kids = []
        for i, s in enumerate(cur):
            kids += BS.expand((s, nstart, seed + 31 * i + lev, keep))
        if lev == 5 - m:
            cur = kids
            break
        cur = BS.select(kids, Bm)
    return np.array([list(s['x']) + list(map(float, s['seq'])) for s in cur]) if cur else np.zeros((0, 34))


def run_suffix(pool, X, cuts=(1, 2, 3), Bm=40, nstart=6, keep=3, seed=0):
    """Rioptimizon te gjitha turat X per cdo prerje ne `cuts`, paralelisht."""
    tasks = [(x, m, Bm, nstart, keep, seed + 7 * j + m) for j, x in enumerate(X) for m in cuts]
    res = [r for r in pool.map(suffix_task, tasks, chunksize=1) if len(r)]
    return np.vstack(res) if res else np.zeros((0, 34))
