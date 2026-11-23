"""Optimizimi lokal i NJE kembe (planet p -> planet q).

Nje kembe ka 4 variabla: beta (orientimi i fly-by-it), rp (distanca e fly-by-it),
eta (kur ndizet motori, si pjese e kohes) dhe T (kohezgjatja). Kostoja = DV i manovres.
Problemi 34-dimensional zberthehet keshtu ne shume probleme te vogla 4-dimensionale.

Variablat optimizohen ne trajten z = [beta, log(rp), eta, T]; logaritmi e ben rp me te
"lehte" per optimizuesin, sepse rp ndryshon nga 1.1 ne 100 rreze planetare.
"""
import numpy as np
from numba import njit
import evaluator as F

LRMIN = np.log(1.1)      # kufijte e log(rp), sipas rregullave te ESA-s
LRMAX = np.log(100.0)


@njit(cache=True, error_model='numpy')
def _clip(z, Tlo, Thi):
    """I con variablat brenda kufijve (beta eshte periodik)."""
    y = z.copy()
    y[0] = y[0] % (2 * np.pi)
    y[1] = min(max(y[1], LRMIN), LRMAX)
    y[2] = min(max(y[2], 1e-3), 1 - 1e-3)
    y[3] = min(max(y[3], Tlo), Thi)
    return y


@njit(cache=True, error_model='numpy')
def leg_obj(z, p, q, t0, vin, lam, Tlo, Thi):
    """Funksioni objektiv i kembes: DV + lam*T (1e9 nese trajektorja eshte e pavlefshme)."""
    y = _clip(z, Tlo, Thi)
    dv, va, ok = F.leg(p, q, t0, vin, y[0], np.exp(y[1]), y[2], y[3])
    if not ok:
        return 1e9
    return dv + lam * y[3]


@njit(cache=True, error_model='numpy')
def nm_leg(z0, step, p, q, t0, vin, lam, Tlo, Thi, maxit):
    """Nelder-Mead ne 4 dimensione (simpleks me 5 kulme).
    Ne cdo iteracion kulmi me i keq zevendesohet me reflektim / zgjerim / tkurrje;
    nese asgje nuk ndihmon, i gjithe simpleksi tkurret drejt kulmit me te mire."""
    n = 4
    S = np.empty((n + 1, n))
    fS = np.empty(n + 1)
    S[0] = z0
    for i in range(n):
        S[i + 1] = z0
        S[i + 1, i] += step[i]
    for i in range(n + 1):
        fS[i] = leg_obj(S[i], p, q, t0, vin, lam, Tlo, Thi)
    for it in range(maxit):
        order = np.argsort(fS)
        S = S[order]; fS = fS[order]
        if fS[n] - fS[0] < 1e-9 * (1 + abs(fS[0])) and it > 20:
            break                                         # konvergjenca
        c = np.zeros(n)                                   # qendra e kulmeve pervec me te keqit
        for i in range(n):
            c += S[i]
        c /= n
        xr = c + (c - S[n]); fr = leg_obj(xr, p, q, t0, vin, lam, Tlo, Thi)       # reflektim
        if fr < fS[0]:
            xe = c + 2 * (c - S[n]); fe = leg_obj(xe, p, q, t0, vin, lam, Tlo, Thi)  # zgjerim
            if fe < fr:
                S[n] = xe; fS[n] = fe
            else:
                S[n] = xr; fS[n] = fr
        elif fr < fS[n - 1]:
            S[n] = xr; fS[n] = fr
        else:                                                                    # tkurrje
            if fr < fS[n]:
                xc = c + 0.5 * (xr - c)
            else:
                xc = c + 0.5 * (S[n] - c)
            fc = leg_obj(xc, p, q, t0, vin, lam, Tlo, Thi)
            if fc < min(fr, fS[n]):
                S[n] = xc; fS[n] = fc
            else:                                                                # tkurrje e plote
                for i in range(1, n + 1):
                    S[i] = S[0] + 0.5 * (S[i] - S[0])
                    fS[i] = leg_obj(S[i], p, q, t0, vin, lam, Tlo, Thi)
    b = np.argmin(fS)
    return _clip(S[b], Tlo, Thi), fS[b]


@njit(cache=True, error_model='numpy')
def multistart_leg(p, q, t0, vin, lam, Tlo, Thi, nstart, seed, maxit):
    """Fillime te rastit: per secilin nga `nstart` fillime merret me i miri nga 30 pika te rastit,
    pastaj dy kalime Nelder-Mead (hap i madh, pastaj i vogel) brenda dritares kohore [Tlo, Thi].
    Kthen tabelen (nstart, 8) = [beta, rp, eta, T, dv, f, ok, -] dhe shpejtesite e mberritjes."""
    np.random.seed(seed)
    out = np.full((nstart, 8), 1e9)
    varr = np.zeros((nstart, 3))
    for k in range(nstart):
        best = 1e18
        zb = np.zeros(4)
        for s in range(30):
            z = np.array([np.random.uniform(0, 2 * np.pi), np.random.uniform(LRMIN, LRMAX * 0.6),
                          np.random.uniform(0.02, 0.98), np.random.uniform(Tlo, Thi)])
            f = leg_obj(z, p, q, t0, vin, lam, Tlo, Thi)
            if f < best:
                best = f; zb = z
        step = np.array([0.5, 0.5, 0.15, 0.1 * (Thi - Tlo) + 0.2])
        z, f = nm_leg(zb, step, p, q, t0, vin, lam, Tlo, Thi, maxit)
        step = np.array([0.05, 0.05, 0.02, 0.01 * (Thi - Tlo) + 0.02])
        z, f = nm_leg(z, step, p, q, t0, vin, lam, Tlo, Thi, maxit)
        dv, va, ok = F.leg(p, q, t0, vin, z[0], np.exp(z[1]), z[2], z[3])
        out[k, 0] = z[0]; out[k, 1] = np.exp(z[1]); out[k, 2] = z[2]; out[k, 3] = z[3]
        out[k, 4] = dv; out[k, 5] = f; out[k, 6] = 1.0 if ok else 0.0
        varr[k] = va
    return out, varr


@njit(cache=True, error_model='numpy')
def _lon(r):
    """Gjatesia (kendi ne planin xy) e nje pozicioni."""
    return np.arctan2(r[1], r[0])


@njit(cache=True, error_model='numpy')
def seeded_candidates(p, q, t0, vin, Tmin, Tmax, nb, nr, maxc):
    """Fillime "fizike" per kemben p -> q.

    Per nje rrjete (beta, rp) llogaritet orbita pas fly-by-it. Nese ajo pret rrezen e orbites se q,
    gjenden momentet kur anija kalon ate rreze (ne cdo revolucion, hyrje dhe dalje) dhe koha kur
    planeti q arrin ne te njejten gjatesi. Keto kohe T jane pikerisht vendet ku DV(T) ka lugina
    te ngushta, qe fillimet e rastit shpesh i humbasin. Cdo kandidat merr nje vleresim te shpejte
    (DV me i mire per 4 vlera te eta). Kthen tabelen [beta, log rp, eta, T, DV_i_perafert]."""
    rP, vP = F.eph(p, t0)
    rq = F.ELEM[q, 0]
    Pq = 2 * np.pi / F.NMOT[q] / F.DAY2SEC            # periudha e planetit q [dite]
    cands = np.zeros((maxc, 5))
    nc = 0
    for i in range(nb):
        beta = 2 * np.pi * (i + 0.5) / nb
        for j in range(nr):
            lrp = LRMIN + (LRMAX * 0.75 - LRMIN) * j / max(1, nr - 1)
            vo = F.fb_prop(vin, vP, np.exp(lrp) * F.EARTH_RADIUS, beta, F.PMU[p])
            R = np.sqrt(rP @ rP)
            eps = 0.5 * (vo @ vo) - F.MU / R          # energjia orbitale pas fly-by-it
            if eps >= 0:
                continue                               # orbite e hapur: s'kthehet me
            a = -F.MU / (2 * eps)
            ev = ((vo @ vo - F.MU / R) * rP - (rP @ vo) * vo) / F.MU
            e = np.sqrt(ev @ ev)
            if a * (1 - e) > rq or a * (1 + e) < rq:
                continue                               # orbita nuk e pret rrezen e q
            n = np.sqrt(F.MU / a ** 3)
            P = 2 * np.pi / n / F.DAY2SEC
            cE0 = min(1.0, max(-1.0, (1 - R / a) / e))  # anomalia ekscentrike aktuale
            E0 = np.arccos(cE0)
            if rP @ vo < 0:
                E0 = -E0
            M0 = E0 - e * np.sin(E0)
            Ec = np.arccos(min(1.0, max(-1.0, (1 - rq / a) / e)))   # anomalia kur r = rq
            for sgn in (-1.0, 1.0):                    # kalimi hyres dhe dales
                Ex = sgn * Ec
                Mx = Ex - e * np.sin(Ex)
                dt0 = ((Mx - M0) % (2 * np.pi)) / n / F.DAY2SEC
                k = 0
                while True:                            # nje kalim per cdo revolucion
                    tc = dt0 + k * P
                    k += 1
                    if tc > Tmax + Pq:
                        break
                    if tc < Tmin - Pq:
                        continue
                    rc, vc, ok = F.propagate(rP, vo, tc * F.DAY2SEC, F.MU)
                    if not ok:
                        break
                    rqc, vqc = F.eph(q, t0 + tc)
                    dth = (_lon(rc) - _lon(rqc)) % (2 * np.pi)
                    tq = tc + dth / (2 * np.pi) * Pq   # kur planeti arrin ne piken e kalimit
                    for T in (tq, tq - Pq):
                        if T < Tmin or T > Tmax:
                            continue
                        if nc < maxc:
                            cands[nc, 0] = beta; cands[nc, 1] = lrp
                            cands[nc, 2] = 0.5; cands[nc, 3] = T
                            best = 1e18
                            for et in (0.05, 0.3, 0.6, 0.9):
                                dv, va, ok2 = F.leg(p, q, t0, vin, beta, np.exp(lrp), et, T)
                                if ok2 and dv < best:
                                    best = dv; cands[nc, 2] = et
                            cands[nc, 4] = best
                            nc += 1
    return cands[:nc]


@njit(cache=True, error_model='numpy')
def refine_list(p, q, t0, vin, starts, Tlo, Thi, maxit):
    """Nelder-Mead nga nje liste fillimesh (p.sh. nga seeded_candidates). Format si multistart_leg."""
    n = starts.shape[0]
    out = np.full((n, 8), 1e9)
    varr = np.zeros((n, 3))
    for k in range(n):
        z0 = starts[k, :4].copy()
        step = np.array([0.1, 0.1, 0.1, 0.5])
        z, f = nm_leg(z0, step, p, q, t0, vin, 0.0, Tlo, Thi, maxit)
        step = np.array([0.02, 0.02, 0.02, 0.05])
        z, f = nm_leg(z, step, p, q, t0, vin, 0.0, Tlo, Thi, maxit)
        dv, va, ok = F.leg(p, q, t0, vin, z[0], np.exp(z[1]), z[2], z[3])
        out[k, 0] = z[0]; out[k, 1] = np.exp(z[1]); out[k, 2] = z[2]; out[k, 3] = z[3]
        out[k, 4] = dv; out[k, 5] = f; out[k, 6] = 1.0 if ok else 0.0
        varr[k] = va
    return out, varr
