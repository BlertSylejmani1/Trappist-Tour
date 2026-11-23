"""Beam search kembe pas kembe, me heuristike per te ardhmen.

Gjendja = (sekuenca e planeteve deri tani, koha t, shpejtesia e mberritjes, DV i grumbulluar,
vektori i pjesshem i vendimeve). Nje ture ka 7 nivele: hyrja + 6 kembe.

Cikli:
  1. Gjendjet fillestare: mijera drejtime hyrjeje (u, v), secila me T0 qe jep DV_init ~ 0.
  2. Zgjerimi: per cdo gjendje dhe cdo planet te pavizituar optimizohet kemba tjeter
     (leg_optimizer.py) ne shume dritare kohore; ruhen 1-2 zgjidhjet me te mira per dritare,
     qe te mbetet kompromisi DV / kohe.
  3. Perzgjedhja: mbahen B gjendjet me te mira sipas DV + lambda*(t + h) per 10 vlera lambda
     (heuristic.py), plus nje pjese sipas DV + lambda*t per diversitet.
Kompleksiteti: O(6 * B * (7-k) * S * I), ku S = numri i fillimeve, I = iteracionet e Nelder-Mead.
"""
import time
import numpy as np
from multiprocessing import Pool
import evaluator as F
import entry as E
import leg_optimizer as L
from heuristic import min_period, h_of

# Dritaret kohore te kembes [dite]: te trasha per fillimet e rastit, te imta per perzgjedhjen.
WINDOWS = [(5, 8), (8, 12), (12, 18), (18, 27), (27, 40), (40, 60), (60, 90), (90, 130), (130, 190), (190, 280)]
FINE = [5, 6.5, 8, 10, 12, 14.5, 17, 20, 23.5, 27, 31, 36, 42, 48, 55, 63, 72, 82, 93, 105, 118, 132, 146,
        161, 177, 195, 215, 240, 280]
FINE_W = list(zip(FINE[:-1], FINE[1:]))
# Peshat lambda [m/s per dite]: secila perfaqeson nje pike tjeter te frontit Pareto.
LAMS = [0.3, 0.7, 1.2, 2.0, 3.0, 4.5, 7.0, 10.0, 15.0, 25.0]


def make_entry(args):
    """Gjendjet fillestare: n drejtime te rastit (u, v) drejt planetit p, me DV_init < 2 m/s."""
    p, n, seed = args
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        u, v = rng.random(), rng.uniform(.25, .75)
        T0, dv0, vin, ok = E.solve_T0(u, v, p)
        if ok and dv0 < 2.0:
            out.append(dict(seq=[p], t=T0, dv=dv0, vin=vin.copy(), x=[u, v, T0]))
    return out


def expand(args):
    """Zgjeron nje gjendje: per cdo planet q te pavizituar gjen kembet me te mira p -> q.
    Fillimet vijne nga dy burime: fizike (seeded_candidates) dhe te rastit (multistart_leg).
    Per cdo dritare te imte kohore ruhen `keep` femijet me DV me te vogel (pa dublikata)."""
    st, nstart, seed, keep = args
    p = st['seq'][-1]
    Tmax = WINDOWS[-1][1]
    kids = []
    for q in range(7):
        if q in st['seq']:
            continue
        res, vas = [], []
        C = L.seeded_candidates(p, q, st['t'], st['vin'], 5.0, Tmax, 24, 8, 4000)
        if len(C):
            starts = []
            for (a, b) in FINE_W:                         # 3 kandidatet me te mire per dritare
                m = np.where((C[:, 3] >= a) & (C[:, 3] < b))[0]
                if len(m):
                    starts += list(m[np.argsort(C[m, 4])[:3]])
            if starts:
                out, va = L.refine_list(p, q, st['t'], st['vin'], C[starts], 5.0, Tmax, 200)
                res.append(out); vas.append(va)
        for wi, (a, b) in enumerate(WINDOWS):
            out, va = L.multistart_leg(p, q, st['t'], st['vin'], 0.0, a, b, nstart, seed + 97 * q + wi, 250)
            res.append(out); vas.append(va)
        out = np.vstack(res); va = np.vstack(vas)
        for (a, b) in FINE_W:
            m = np.where((out[:, 3] >= a) & (out[:, 3] < b) & (out[:, 6] > 0) & (out[:, 4] < 3000))[0]
            if len(m) == 0:
                continue
            m = m[np.argsort(out[m, 4])]
            taken = []
            for i in m:
                if len(taken) >= keep:
                    break
                if any(abs(out[i, 3] - out[j, 3]) < 0.5 and abs(out[i, 4] - out[j, 4]) < 1.0 for j in taken):
                    continue                              # pothuajse e njejta zgjidhje
                taken.append(i)
                beta, rp, eta, T, dv = out[i, :5]
                kids.append(dict(seq=st['seq'] + [q], t=st['t'] + T, dv=st['dv'] + dv, vin=va[i].copy(),
                                 x=st['x'] + [beta, rp, eta, T]))
    return kids


def add_heur(states):
    """Shton P_min dhe kohen e mbetur te parashikuar h ne cdo gjendje."""
    for s in states:
        if 'h' not in s:
            k = len(s['seq']) - 1
            s['P'] = min_period(s['seq'][-1], s['t'], s['vin'])
            s['h'] = h_of(k, s['P']) if k <= 5 else 0.0
    return states


def _lam_pick(states, n, tkey, seen):
    """Merr n gjendje duke alternuar renditjet DV + lambda*tkey per te gjitha lambda-t."""
    dv = np.array([s['dv'] for s in states]); t = np.array([tkey(s) for s in states])
    orders = [np.argsort(dv + l * t) for l in LAMS]
    chosen = []; k = 0
    while len(chosen) < n and k < len(states):
        for o in orders:
            i = o[k]
            key = (tuple(states[i]['seq']), round(states[i]['dv'], 1), round(states[i]['t'], 2))
            if key not in seen:
                seen.add(key); chosen.append(i)
                if len(chosen) >= n:
                    break
        k += 1
    return chosen


def select(states, B, frac_heur=0.6):
    """Perzgjedhja e beam-it: frac_heur e vendeve sipas DV + lambda*(t + h) (heuristika),
    pjesa tjeter sipas DV + lambda*t (si pa heuristike), per te ruajtur diversitetin."""
    if len(states) <= B:
        return states
    add_heur(states)
    seen = set()
    a = _lam_pick(states, int(B * frac_heur), lambda s: s['t'] + s['h'], seen)
    b = _lam_pick(states, B - len(a), lambda s: s['t'], seen)
    return [states[i] for i in a + b]


def select_entries(states, n):
    """Hyrjet: gjysma me kapjen me te mire (P_min me i vogel), gjysma te rastit."""
    add_heur(states)
    o = np.argsort([s['P'] for s in states])
    best = list(o[:n // 2])
    rest = np.setdiff1d(np.arange(len(states)), best)
    rng = np.random.default_rng(len(states))
    rnd = list(rng.choice(rest, size=min(len(rest), n - len(best)), replace=False))
    return [states[i] for i in best + rnd]


def run(B=400, n_entry=1500, nstart=4, seed=0, procs=2, firsts=(0,), keep=2, frac_heur=0.6, verbose=True):
    """Nje kerkim i plote. B = gjeresia e beam-it, n_entry = drejtime hyrjeje per proces,
    firsts = planetet e mundshme te para (0 = b). Kthen matricen e turave (N x 34)."""
    t0 = time.time()
    with Pool(procs) as pool:
        ent = pool.map(make_entry, [(p, n_entry, seed * 100 + p * 7 + k) for p in firsts for k in range(procs)])
        st = select_entries([s for e in ent for s in e], 2 * B)
        if verbose:
            print('entry states', len(st), 'best P_min %.1f d' % min(s['P'] for s in st), flush=True)
        for lev in range(6):
            kids = pool.map(expand, [(s, nstart, seed * 1000 + 13 * i + lev, keep) for i, s in enumerate(st)],
                            chunksize=1)
            kids = [k for kk in kids for k in kk]
            if lev == 5:
                st = kids
                break
            st = select(kids, B, frac_heur)
            if verbose:
                print(f'level {lev} kids {len(kids)} beam {len(st)} minDV {min(s["dv"] for s in st):.1f} '
                      f'minT {min(s["t"] for s in st):.1f}  {time.time() - t0:.0f}s', flush=True)
    return np.array([list(s['x']) + list(map(float, s['seq'])) for s in st])
