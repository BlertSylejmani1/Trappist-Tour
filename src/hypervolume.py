"""Hypervolume, verifikimi zyrtar dhe zgjedhja e <= 100 zgjidhjeve per dorezim.

Hypervolume = siperfaqja e "mbuluar" nga fronti Pareto deri te pika e referimit (2500 m/s, 4000 dite).
Cdo zgjidhje (DV, T) mbulon drejtkendeshin [T, 4000] x [DV, 2500]; score = siperfaqja e bashkimit.
"""
import io
import contextlib
import numpy as np
import pygmo as pg
import evaluator as F

REF = np.array([2500., 4000.])


def hv(points):
    """Hypervolume i nje bashkesie pikash (DV, T); pikat jashte referimit injorohen."""
    P = [p for p in points if p[0] < 2500 and p[1] < 4000]
    if not P:
        return 0.0
    return pg.hypervolume(np.array(P)).compute(REF)


def eval_all(X):
    """Vleresim i shpejte (numba) i shume turave: rreshtat [DV, T, ineq]."""
    return np.array([F.fitness(np.asarray(x, float))[:3] for x in X])


def official_eval(x):
    """Vleresim me kodin ZYRTAR te ESA-s. Kthen (DV, T, e_vlefshme)."""
    import official as O
    with contextlib.redirect_stdout(io.StringIO()):
        try:
            f = O.udp.fitness(list(map(float, x)))
        except Exception:
            return np.inf, np.inf, False
    if len(f) != 4 or not np.isscalar(f[0]):
        return np.inf, np.inf, False
    dv, T, eq, iq = map(float, f)
    return dv, T, (eq == 0 and iq <= 0 and dv < 2500 and T < 4000)


def official_score(points):
    """Score-i me funksionin zyrtar combine_scores (kthehet pozitiv)."""
    import official as O
    return -O.combine_scores([list(p) for p in points])


def valid_front(X, R=None):
    """Fronti Pareto i zgjidhjeve qe jane te VLEFSHME sipas kodit zyrtar.

    Disa tura jane ne kufi me distancen minimale nga ylli: kalojne kontrollin e shpejte, por jo
    ate zyrtar (dallime numerike ~1e-13). Nese nuk hiqen, ato "fshehin" (dominojne) tura te
    vlefshme. Prandaj: llogaritet fronti -> kontrollohet zyrtarisht -> hiqen te pavlefshmet ->
    perseritet derisa i gjithe fronti te jete i vlefshem.
    Kthen indekset e frontit (pa dublikata) dhe vlerat zyrtare [DV, T]."""
    if R is None:
        R = eval_all(X)
    pool = np.where((R[:, 2] <= 0) & (R[:, 0] < 2500) & (R[:, 1] < 4000))[0]
    checked = {}
    while True:
        nd = pool[pg.non_dominated_front_2d(R[pool, :2])]
        for i in nd:
            if i not in checked:
                checked[i] = official_eval(X[i])
        bad = set(i for i in nd if not checked[i][2])
        if not bad:
            break
        pool = np.array([i for i in pool if i not in bad])
    idx = np.array(list(nd))
    Ro = np.array([[checked[i][0], checked[i][1]] for i in idx])
    _, u = np.unique(np.round(Ro, 6), axis=0, return_index=True)
    return idx[u], Ro[u]


def select_greedy(P, n=100):
    """Zgjedhje grykese: heq me radhe piken me kontributin me te vogel ne hypervolume. O(n^2 log n)."""
    idx = list(range(len(P)))
    while len(idx) > n:
        c = pg.hypervolume(P[idx]).contributions(REF)
        idx.pop(int(np.argmin(c)))
    return np.array(idx)


def select_exact(P, n=100):
    """Zgjedhja OPTIMALE e n pikave nga nje front 2D, me programim dinamik
    (ideja: Bringmann, Friedrich & Klitzke, GECCO 2014).

    Pikat renditen sipas T (ne front DV zbret). Hypervolume i nje zinxhiri p_1..p_k eshte
        sum_m (4000 - T_m) * (DV_{m-1} - DV_m),   me DV_0 = 2500,
    pra varet vetem nga ciftet e njepasnjeshme. Rekurrenca:
        f[c][j] = max_{i<j} f[c-1][i] + (4000 - T_j) * (DV_i - DV_j)
    jep zgjidhjen optimale ne O(n * m^2), ku m = numri i pikave te frontit."""
    P = np.asarray(P, float)
    o = np.argsort(P[:, 1])
    m = len(P)
    if m <= n:
        return o
    y = np.r_[2500.0, P[o, 0]]               # nyja 0 = fillimi virtual
    a = np.r_[0.0, 4000.0 - P[o, 1]]
    NEG = -1e30
    f = np.full(m + 1, NEG); f[0] = 0.0
    back = np.zeros((n + 1, m + 1), dtype=np.int32)
    best_val, best_c, best_j = NEG, 0, 0
    tri = np.tril(np.ones((m + 1, m + 1), bool), -1)         # vetem i < j
    for c in range(1, n + 1):
        cand = f[None, :] + a[:, None] * (y[None, :] - y[:, None])   # cand[j, i]
        cand = np.where(tri, cand, NEG)
        back[c] = np.argmax(cand, axis=1)
        f = cand[np.arange(m + 1), back[c]]
        f[0] = NEG
        j = int(np.argmax(f))
        if f[j] > best_val:
            best_val, best_c, best_j = f[j], c, j
    sel = []; j = best_j
    for c in range(best_c, 0, -1):                               # rindertimi i zgjidhjes
        sel.append(j - 1); j = back[c][j]
    return o[np.array(sel[::-1])]


def best_subset(X, R=None, n=100):
    """Nga cdo bashkesi turash: fronti i vlefshem zyrtarisht -> n pikat optimale.
    Kthen (vektoret e zgjedhur, score zyrtar)."""
    idx, Ro = valid_front(X, R)
    sel = select_exact(Ro, n)
    return X[idx[sel]], official_score(Ro[sel])


def to_submission(X):
    """Formati JSON i dorezimit te ESA-s."""
    return [{"decisionVector": [list(map(float, x)) for x in X],
             "problem": "trappist-tour", "challenge": "spoc-trappist-tour"}]
