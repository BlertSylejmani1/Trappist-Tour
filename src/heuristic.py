"""Heuristika per te ardhmen (ideja kryesore e algoritmit, si funksioni h tek A*).

Pas mberritjes te nje planet, llogaritet P_min: periudha me e shkurter e orbites qe mund te
jape fly-by-i ne kete planet (pa motor). Sa me e vogel P_min, aq me e "kapur" eshte anija dhe
aq me te shkurtra do te jene kembet e mbetura. Mbi turat tona me te mira, koha e mbetur e
tures parashikohet shume mire nga regresioni linear

    h(k, P_min) = a_k + b_k * P_min

ku k = numri i kembeve te pershkuara pas hyrjes. Korrelacioni ishte 0.88 - 0.97 (k = 1..5) mbi
arkivin e turave kur u ndertua heuristika. Per k = 0 heuristika ka pak informacion, sepse P_min
pas hyrjes eshte pothuajse i njejte per te gjitha drejtimet (~134 dite, kufiri fizik i kapjes). Beam search-i i rendit gjendjet sipas
DV + lambda * (t + h) ne vend te DV + lambda * t, ashtu si A* perdor f = g + h.
"""
import numpy as np
from numba import njit
import evaluator as F

# Koeficientet (a_k, b_k), te nxjerra me fit_heuristic() nga turat tona me te mira.
HEUR = {0: (-7.8, 2.43), 1: (28.2, 2.93), 2: (18.8, 2.86), 3: (11.6, 2.47), 4: (4.8, 1.97), 5: (1.9, 1.04)}


def states_of(x):
    """Rindjek nje ture te plote dhe kthen gjendjen ne mberritje te cdo planeti:
    lista me (planeti, koha [dite], shpejtesia e mberritjes, DV i grumbulluar)."""
    x = np.asarray(x, float); seq = x[27:].astype(int)
    dv0, vin, ok = F.entry_leg(x[0], x[1], x[2], seq[0])
    out = [(seq[0], x[2], vin.copy(), dv0)]
    t = x[2]; dv = dv0
    for i in range(6):
        b, rp, eta, T = x[3 + 4 * i:7 + 4 * i]
        d, va, ok = F.leg(seq[i], seq[i + 1], t, vin, b, rp, eta, T)
        t += T; dv += d; vin = va
        out.append((seq[i + 1], t, va.copy(), dv))
    return out


@njit(cache=True, error_model='numpy')
def min_period(p, t, vin, nb=36, nr=16):
    """P_min [dite]: periudha me e shkurter e arritshme me fly-by te planeti p.
    Kerkohet ne nje rrjete nb x nr te (beta, rp) energjia orbitale me e ulet pas fly-by-it;
    nga ajo, gjysmaksi a = -mu / (2 eps) dhe periudha P = 2 pi sqrt(a^3 / mu).
    Kthen 1e6 nese anija mbetet ne orbite te hapur (nuk kapet)."""
    rP, vP = F.eph(p, t); R = np.sqrt(rP @ rP)
    best = 1e18
    for i in range(nb):
        beta = 2 * np.pi * i / nb
        for j in range(nr):
            rp = 1.1 * np.exp(np.log(100 / 1.1) * j / (nr - 1))
            vo = F.fb_prop(vin, vP, rp * F.EARTH_RADIUS, beta, F.PMU[p])
            e = 0.5 * (vo @ vo) - F.MU / R
            if e < best:
                best = e
    if best >= 0:
        return 1e6
    a = -F.MU / (2 * best)
    return 2 * np.pi * np.sqrt(a ** 3 / F.MU) / 86400.0


def h_of(k, P):
    """Koha e mbetur e parashikuar [dite] pas k kembeve, kur periudha minimale eshte P."""
    a, b = HEUR[k]
    return a + b * min(P, 400.0)


def fit_heuristic(X, R):
    """Rillogarit koeficientet HEUR nga turat X (me rezultatet R = [DV, T, ...]).
    Per cdo ture dhe cdo planet te ndermjetem: P_min aty dhe koha qe mbeti realisht deri ne fund.
    Perdoren vetem turat qe harxhojne pak karburant pas asaj pike (< 100 m/s)."""
    rows = []
    for x, r in zip(X, R):
        for k, (p, t, v, dv) in enumerate(states_of(x)[:-1]):
            rows.append((k, min_period(p, t, v), r[1] - t, r[0] - dv))
    rows = np.array(rows)
    coef, corr = {}, {}
    for k in range(6):
        m = (rows[:, 0] == k) & (rows[:, 3] < 100) & ((rows[:, 1] < 300) | (k == 0))
        P, y = rows[m, 1], rows[m, 2]
        A = np.c_[np.ones_like(P), P]
        c, *_ = np.linalg.lstsq(A, y, rcond=None)
        coef[k] = (round(float(c[0]), 1), round(float(c[1]), 2))
        corr[k] = float(np.corrcoef(P, y)[0, 1]) if len(P) > 2 else float('nan')
    return coef, corr
