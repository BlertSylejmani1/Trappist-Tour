"""Kemba e hyrjes: nga sfera 10 AU deri te planeti i pare.

Anija nis me shpejtesi 10 km/s; cdo devijim nga 10 km/s kushton DV_init = | |v1| - 10 km/s |.
Per nje drejtim nisjeje (u, v) ekziston nje kohezgjatje "natyrale" T0 (~1500 dite) per te cilen
DV_init ~ 0. Ky modul e gjen ate T0.
"""
import numpy as np
from numba import njit
import evaluator as F


@njit(cache=True, error_model='numpy')
def solve_T0(u, v, p):
    """Gjen T0 qe minimizon DV_init per drejtimin (u, v) dhe planetin e pare p.
    1) skanim i trashe me hap 1 dite ne [1440, 1560];
    2) kerkim me prerjen e arte (golden section) brenda +-1 dite rreth me te mirit.
    Kthen (T0, DV_init, shpejtesia e mberritjes te p, ok)."""
    best = 1e18; Tb = 1500.0
    for i in range(121):
        T0 = 1440.0 + i
        dv, va, ok = F.entry_leg(u, v, T0, p)
        if ok and dv < best:
            best = dv; Tb = T0
    gr = 0.6180339887                      # 1/phi, raporti i prerjes se arte
    a = Tb - 1.0; b = Tb + 1.0
    c = b - gr * (b - a); d = a + gr * (b - a)
    fc, va, okc = F.entry_leg(u, v, c, p)
    fd, va, okd = F.entry_leg(u, v, d, p)
    for _ in range(60):
        if fc < fd:                        # minimumi eshte ne [a, d]
            b = d; d = c; fd = fc
            c = b - gr * (b - a); fc, va, okc = F.entry_leg(u, v, c, p)
        else:                              # minimumi eshte ne [c, b]
            a = c; c = d; fc = fd
            d = a + gr * (b - a); fd, va, okd = F.entry_leg(u, v, d, p)
    T0 = 0.5 * (a + b)
    dv, va, ok = F.entry_leg(u, v, T0, p)
    return T0, dv, va, ok
