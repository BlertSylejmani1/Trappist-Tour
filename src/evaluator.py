"""Vleresuesi i shpejte (numba) i problemit Trappist Tour.

Eshte rishkrim i funksionit zyrtar `fitness` (official.py), por rreth 4 here me i shpejte,
sepse kerkimi e therret miliona here. Riprodhon saktesisht pykep 2.6:
  - efemeridat e planeteve nga elementet Kepleriane,
  - problemin e Lambert-it (algoritmi i Izzo-s, pa revolucione, drejtim prograd),
  - propagimin Lagrangian te orbites,
  - fly-by-in (fb_prop) dhe kontrollin e distances minimale nga ylli.
Validimi: ne 3000 vektore te rastit, gabimi relativ mesatar ndaj pykep eshte ~1e-13.

Njesite: distanca [m], shpejtesia [m/s], koha [dite] ne nderfaqe dhe [s] brenda llogaritjeve.
Planetet indeksohen 0..6 = b, c, d, e, f, g, h.
"""
import numpy as np
from numba import njit

# Konstantet fizike, te njejtat si ne official.py
DAY2SEC = 86400.0
AU = 149597870691.0
DEG2RAD = np.pi / 180.0
EARTH_RADIUS = 6378137.0
G = 6.67430E-11
SM = 1.989E30
EM = 5.9722E24
MU = G * 8.98266512E-2 * SM
SAFE_R = 10 * 83484000.0        # distanca minimale e lejuar nga ylli [m]
R_START = 10 * AU                # anija nis nga sfera me rreze 10 AU
V_START = 1e4                    # me shpejtesi 10 km/s (devijimi kushton DV)

# Masat dhe elementet orbitale [a, e, i, Omega, omega, M] te planeteve b..h
_masses = np.array([1.36004499e+00, 1.29688971e+00, 3.85884170e-01, 6.88627613e-01,
                    1.03510927e+00, 1.31679654e+00, 3.20988718e-01]) * EM
_el = np.array([
    [0.01154035, 0.0158073, 0.02226718, 0.02926548, 0.0385054, 0.04683828, 0.06191385],
    [0.00455247, 0.00107698, 0.00624362, 0.00577129, 0.00861154, 0.00400202, 0.00378361],
    [1.1, 1.3, 0.5, -0.3, 0.01, 1.2, 2.3],
    [238.35940145, 221.02196847, 90.26523049, 164.16001523, 171.42668515, 166.63349296, 282.808809],
    [126.22225948, 80.78036351, 147.44525146, 315.09376191, 179.82855007, 25.31242216, 176.19682451],
    [33.40998384, -89.79506054, -28.38482427, 45.00185452, -145.16404128, 12.59883066, -12.12072197]]).T.copy()
_el[:, 0] *= AU
_el[:, 2:] *= DEG2RAD
ELEM = _el
PMU = G * _masses                # parametri gravitacional i cdo planeti
NMOT = np.sqrt(MU / ELEM[:, 0] ** 3)   # levizja mesatare [rad/s]


@njit(cache=True, error_model='numpy')
def eph(p, t_days):
    """Pozicioni dhe shpejtesia e planetit p ne kohen t_days (dite pas fillimit te misionit).
    Anomalia mesatare rritet linearisht; ekuacioni i Kepler-it zgjidhet me Newton."""
    a = ELEM[p, 0]; e = ELEM[p, 1]; inc = ELEM[p, 2]; W = ELEM[p, 3]; w = ELEM[p, 4]
    M = ELEM[p, 5] + NMOT[p] * t_days * DAY2SEC
    # pykep m2e: Newton on E - e sinE - M
    E = M + e * np.sin(M)
    for _ in range(50):
        dE = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E -= dE
        if abs(dE) < 1e-15:
            break
    b = a * np.sqrt(1.0 - e * e)
    n = np.sqrt(MU / a ** 3)
    xper = a * (np.cos(E) - e)
    yper = b * np.sin(E)
    xdotper = -(a * n * np.sin(E)) / (1.0 - e * np.cos(E))
    ydotper = (b * n * np.cos(E)) / (1.0 - e * np.cos(E))
    cW = np.cos(W); sW = np.sin(W); cw = np.cos(w); sw = np.sin(w); ci = np.cos(inc); si = np.sin(inc)
    R11 = cW * cw - sW * sw * ci; R12 = -cW * sw - sW * cw * ci
    R21 = sW * cw + cW * sw * ci; R22 = -sW * sw + cW * cw * ci
    R31 = sw * si; R32 = cw * si
    r = np.empty(3); v = np.empty(3)
    r[0] = R11 * xper + R12 * yper; r[1] = R21 * xper + R22 * yper; r[2] = R31 * xper + R32 * yper
    v[0] = R11 * xdotper + R12 * ydotper; v[1] = R21 * xdotper + R22 * ydotper; v[2] = R31 * xdotper + R32 * ydotper
    return r, v


@njit(cache=True, error_model='numpy')
def _cross(a, b):
    c = np.empty(3)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c


@njit(cache=True, error_model='numpy')
def _norm(a):
    return np.sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2])


@njit(cache=True, error_model='numpy')
def _hypergeo(z, tol):
    Sj = 1.0; Cj = 1.0; err = 1.0; j = 0
    while err > tol:
        Cj1 = Cj * (3.0 + j) * (1.0 + j) / (2.5 + j) * z / (j + 1)
        Sj = Sj + Cj1
        err = abs(Cj1)
        Cj = Cj1
        j += 1
        if j > 1000:
            break
    return Sj


@njit(cache=True, error_model='numpy')
def _x2tof2(x, lam):
    a = 1.0 / (1.0 - x * x)
    if a > 0:
        alfa = 2.0 * np.arccos(x)
        beta = 2.0 * np.arcsin(np.sqrt(lam * lam / a))
        if lam < 0.0:
            beta = -beta
        return a * np.sqrt(a) * ((alfa - np.sin(alfa)) - (beta - np.sin(beta))) / 2.0
    else:
        alfa = 2.0 * np.arccosh(x)
        beta = 2.0 * np.arcsinh(np.sqrt(-lam * lam / a))
        if lam < 0.0:
            beta = -beta
        return -a * np.sqrt(-a) * ((beta - np.sinh(beta)) - (alfa - np.sinh(alfa))) / 2.0


@njit(cache=True, error_model='numpy')
def _x2tof(x, lam):
    battin = 0.01; lagrange = 0.2
    dist = abs(x - 1.0)
    if dist < lagrange and dist > battin:
        return _x2tof2(x, lam)
    K = lam * lam
    E = x * x - 1.0
    rho = abs(E)
    z = np.sqrt(1 + K * E)
    if dist < battin:
        eta = z - lam * x
        S1 = 0.5 * (1.0 - lam - x * eta)
        Q = _hypergeo(S1, 1e-11)
        Q = 4.0 / 3.0 * Q
        return (eta * eta * eta * Q + 4.0 * lam * eta) / 2.0
    else:
        y = np.sqrt(rho)
        g = x * z - lam * E
        if E < 0:
            l = np.arccos(g)
            d = l
        else:
            f = y * (z - lam * x)
            d = np.log(f + g)
        return (x - lam * z - d / y) / E


@njit(cache=True, error_model='numpy')
def lambert0(r1, r2, tof, mu):
    """Problemi i Lambert-it: cila orbite e con anijen nga r1 ne r2 per kohen tof [s]?
    Kthen shpejtesine ne nisje (v1), ne mberritje (v2) dhe `ok` (False nese nuk konvergjon).
    Algoritmi i Izzo-s (2015), kopje e pykep 2.6: iteracione Householder te rendit te 3-te,
    me formulat Battin/Lagrange/Lancaster sipas zones se variablit x."""
    v1 = np.zeros(3); v2 = np.zeros(3)
    if tof <= 0:
        return v1, v2, False
    c = r2 - r1
    cn = _norm(c); R1 = _norm(r1); R2 = _norm(r2)
    s = (cn + R1 + R2) / 2.0
    ir1 = r1 / R1; ir2 = r2 / R2
    ih = _cross(ir1, ir2)
    nh = _norm(ih)
    if nh == 0.0 or ih[2] == 0.0:
        return v1, v2, False
    ih = ih / nh
    lambda2 = 1.0 - cn / s
    lam = np.sqrt(lambda2)
    if ih[2] < 0.0:
        lam = -lam
        it1 = _cross(ir1, ih); it2 = _cross(ir2, ih)
    else:
        it1 = _cross(ih, ir1); it2 = _cross(ih, ir2)
    it1 = it1 / _norm(it1); it2 = it2 / _norm(it2)
    lambda3 = lam * lambda2
    T = np.sqrt(2.0 * mu / s / s / s) * tof
    T00 = np.arccos(lam) + lam * np.sqrt(1.0 - lambda2)
    T1 = 2.0 / 3.0 * (1.0 - lambda3)
    if T >= T00:
        x0 = -(T - T00) / (T - T00 + 4.0)
    elif T <= T1:
        x0 = T1 * (T1 - T) / (2.0 / 5.0 * (1.0 - lambda2 * lambda3) * T) + 1.0
    else:
        x0 = (T / T00) ** (0.69314718055994529 / np.log(T1 / T00)) - 1.0
    # householder
    it = 0; err = 1.0
    l2 = lambda2; l3 = lambda3
    while err > 1e-5 and it < 15:
        tofx = _x2tof(x0, lam)
        umx2 = 1.0 - x0 * x0
        y = np.sqrt(1.0 - l2 * umx2)
        y2 = y * y; y3 = y2 * y
        DT = 1.0 / umx2 * (3.0 * tofx * x0 - 2.0 + 2.0 * l3 * x0 / y)
        DDT = 1.0 / umx2 * (3.0 * tofx + 5.0 * x0 * DT + 2.0 * (1.0 - l2) * l3 / y3)
        DDDT = 1.0 / umx2 * (7.0 * x0 * DDT + 8.0 * DT - 6.0 * (1.0 - l2) * l2 * l3 * x0 / y3 / y2)
        delta = tofx - T
        DT2 = DT * DT
        xnew = x0 - delta * (DT2 - delta * DDT / 2.0) / (DT * (DT2 - delta * DDT) + DDDT * delta * delta / 6.0)
        err = abs(x0 - xnew)
        x0 = xnew
        it += 1
    x = x0
    if not np.isfinite(x):
        return v1, v2, False
    gamma = np.sqrt(mu * s / 2.0)
    rho = (R1 - R2) / cn
    sigma = np.sqrt(1.0 - rho * rho)
    y = np.sqrt(1.0 - lambda2 + lambda2 * x * x)
    vr1 = gamma * ((lam * y - x) - rho * (lam * y + x)) / R1
    vr2 = -gamma * ((lam * y - x) + rho * (lam * y + x)) / R2
    vt = gamma * sigma * (y + lam * x)
    vt1 = vt / R1; vt2 = vt / R2
    for k in range(3):
        v1[k] = vr1 * ir1[k] + vt1 * it1[k]
        v2[k] = vr2 * ir2[k] + vt2 * it2[k]
    ok = np.isfinite(v1[0]) and np.isfinite(v2[0]) and err < 1e-3
    return v1, v2, ok


@njit(cache=True, error_model='numpy')
def propagate(r0, v0, t, mu):
    """Propagimi Kepler: ku ndodhet anija pas kohes t [s], nese nis nga (r0, v0)?
    Orbita eliptike: zgjidhet ekuacioni i Kepler-it per anomaline ekscentrike (Newton).
    Orbita hiperbolike: Newton i mbrojtur me bisection (gjithmone konvergjon)."""
    R = _norm(r0); V = _norm(v0)
    energy = V * V / 2.0 - mu / R
    a = -mu / 2.0 / energy
    sigma0 = (r0[0] * v0[0] + r0[1] * v0[1] + r0[2] * v0[2]) / np.sqrt(mu)
    if a > 0:
        sqrta = np.sqrt(a)
        DM = np.sqrt(mu / a ** 3) * t
        k = np.floor((DM + np.pi) / (2 * np.pi))
        DMr = DM - 2 * np.pi * k
        DE = DMr
        for _ in range(100):
            f = DE + sigma0 / sqrta * (1 - np.cos(DE)) - (1 - R / a) * np.sin(DE) - DMr
            df = 1 + sigma0 / sqrta * np.sin(DE) - (1 - R / a) * np.cos(DE)
            d = f / df
            DE -= d
            if abs(d) < 1e-14:
                break
        r = a + (R - a) * np.cos(DE) + sigma0 * sqrta * np.sin(DE)
        F = 1 - a / R * (1 - np.cos(DE))
        Gg = a * sigma0 / np.sqrt(mu) * (1 - np.cos(DE)) + R * np.sqrt(a / mu) * np.sin(DE)
        Ft = -np.sqrt(mu * a) / (r * R) * np.sin(DE)
        Gt = 1 - a / r * (1 - np.cos(DE))
    else:
        sqrtma = np.sqrt(-a)
        DN = np.sqrt(-mu / a ** 3) * t
        # f(DH) = -DH + sigma0/sqrtma*(cosh-1) + (1-R/a) sinh - DN, increasing in DH
        lo = -1.0; hi = 1.0
        # bracket
        for _ in range(200):
            fl = -lo + sigma0 / sqrtma * (np.cosh(lo) - 1) + (1 - R / a) * np.sinh(lo) - DN
            if fl < 0:
                break
            lo *= 2
        for _ in range(200):
            fh = -hi + sigma0 / sqrtma * (np.cosh(hi) - 1) + (1 - R / a) * np.sinh(hi) - DN
            if fh > 0:
                break
            hi *= 2
        DH = 0.5 * (lo + hi)
        for _ in range(200):
            f = -DH + sigma0 / sqrtma * (np.cosh(DH) - 1) + (1 - R / a) * np.sinh(DH) - DN
            if f > 0:
                hi = DH
            else:
                lo = DH
            df = -1 + sigma0 / sqrtma * np.sinh(DH) + (1 - R / a) * np.cosh(DH)
            DHn = DH - f / df
            if not (DHn > lo and DHn < hi):
                DHn = 0.5 * (lo + hi)
            if abs(DHn - DH) < 1e-14 * max(1.0, abs(DH)):
                DH = DHn
                break
            DH = DHn
        r = a + (R - a) * np.cosh(DH) + sigma0 * sqrtma * np.sinh(DH)
        F = 1 - a / R * (1 - np.cosh(DH))
        Gg = a * sigma0 / np.sqrt(mu) * (1 - np.cosh(DH)) + R * np.sqrt(-a / mu) * np.sinh(DH)
        Ft = -np.sqrt(-mu * a) / (r * R) * np.sinh(DH)
        Gt = 1 - a / r * (1 - np.cosh(DH))
    rr = F * r0 + Gg * v0
    vv = Ft * r0 + Gt * v0
    return rr, vv, np.isfinite(rr[0]) and np.isfinite(vv[0])


@njit(cache=True, error_model='numpy')
def fb_prop(v_in, v_pla, rp, beta, mu):
    """Fly-by: shpejtesia relative ndaj planetit rrotullohet me kendin delta (madhesia mbetet e njejte).
    rp = distanca minimale nga qendra e planetit, beta = orientimi i planit te fly-by-it."""
    vrel = v_in - v_pla
    v2 = vrel[0] ** 2 + vrel[1] ** 2 + vrel[2] ** 2
    vn = np.sqrt(v2)
    ecc = 1 + rp / mu * v2
    delta = 2 * np.arcsin(1.0 / ecc)
    i_hat = vrel / vn
    j_hat = _cross(i_hat, v_pla)
    j_hat = j_hat / _norm(j_hat)
    k_hat = _cross(i_hat, j_hat)
    return v_pla + vn * np.cos(delta) * i_hat + vn * np.cos(beta) * np.sin(delta) * j_hat \
        + vn * np.sin(beta) * np.sin(delta) * k_hat


@njit(cache=True, error_model='numpy')
def _anom(r, v, mu):
    """Elementet orbitale si te pykep.ic2par: |a|, ekscentriciteti e, anomalia ekscentrike E
    (hiperbolike per e > 1) dhe distanca e periapsis rp."""
    R = _norm(r)
    h = _cross(r, v)
    rv = r[0] * v[0] + r[1] * v[1] + r[2] * v[2]
    V2 = v[0] ** 2 + v[1] ** 2 + v[2] ** 2
    ev = ((V2 - mu / R) * r - rv * v) / mu
    e = _norm(ev)
    p = (h[0] ** 2 + h[1] ** 2 + h[2] ** 2) / mu
    a = abs(p / (1 - e * e))
    c = (ev[0] * r[0] + ev[1] * r[1] + ev[2] * r[2]) / e / R
    if c > 1: c = 1.0
    if c < -1: c = -1.0
    ni = np.arccos(c)
    if rv < 0:
        ni = 2 * np.pi - ni
    if e < 1:
        E = 2.0 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(ni / 2.0))
        rp = a * (1 - e)
    else:
        E = 2.0 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(ni / 2.0))
        rp = -a * (1 - e)
    return a, e, E, rp


@njit(cache=True, error_model='numpy')
def check_distance(r0, v0, dt_days, mu):
    """Kontrolli i distances nga ylli, me te njejten logjike si official._check_distance:
    nese gjate harkut anija kalon periapsis-in, kthen (True, SAFE_R - rp). Kufizimi plotesohet kur d <= 0."""
    a, e, E0, rp = _anom(r0, v0, mu)
    r, v, ok = propagate(r0, v0, dt_days * DAY2SEC, mu)
    if not ok:
        return True, 1e17, False
    period = 2 * np.pi * np.sqrt(a ** 3 / mu) / DAY2SEC
    _, _, E, _ = _anom(r, v, mu)
    if dt_days > period:
        passed = True
    elif E0 > 0:
        passed = (0 < E) and (E < E0)
    else:
        passed = (E > 0) or (dt_days > period / 2)
    return passed, SAFE_R - rp, True


@njit(cache=True, error_model='numpy')
def fitness(x):
    """Vlereson nje vektor te plote 34-dimensional. Kthen (DV total [m/s], T total [dite], ineq, ok).
    ineq <= 0 do te thote se distanca nga ylli respektohet; ok = False kur Lambert/propagimi deshton."""
    seq = np.empty(7, np.int64)
    for i in range(7):
        seq[i] = int(x[27 + i])
    u = x[0]; vv = x[1]
    theta = 2 * np.pi * u
    phi = np.arccos(2 * vv - 1) - np.pi / 2
    rinit = np.array([R_START * np.cos(phi) * np.cos(theta), R_START * np.cos(phi) * np.sin(theta),
                      R_START * np.sin(phi)])
    T = np.empty(7)
    T[0] = x[2]
    for i in range(6):
        T[i + 1] = x[6 + 4 * i]
    iq = -10e16
    t = 0.0
    tP = np.empty(7)
    for i in range(7):
        t += T[i]
        tP[i] = t
    rP0, vP0 = eph(seq[0], tP[0])
    v1, v2, ok = lambert0(rinit, rP0, T[0] * DAY2SEC, MU)
    if not ok:
        return 1e10, 1e10, 1e17, False
    DV = abs(_norm(v1) - V_START)
    c, d, ok = check_distance(rinit, v1, T[0], MU)
    if not ok:
        return 1e10, 1e10, 1e17, False
    if c:
        iq = max(iq, d)
    vend = v2
    rP = rP0; vP = vP0
    for i in range(6):
        beta = x[3 + 4 * i]; rpr = x[4 + 4 * i]; eta = x[5 + 4 * i]
        p = seq[i]
        vout = fb_prop(vend, vP, rpr * EARTH_RADIUS, beta, PMU[p])
        c, d, ok = check_distance(rP, vout, eta * T[i + 1], MU)
        if not ok:
            return 1e10, 1e10, 1e17, False
        if c:
            iq = max(iq, d)
        r, v, ok = propagate(rP, vout, eta * T[i + 1] * DAY2SEC, MU)
        rN, vN = eph(seq[i + 1], tP[i + 1])
        l1, l2, ok = lambert0(r, rN, (1 - eta) * T[i + 1] * DAY2SEC, MU)
        if not ok:
            return 1e10, 1e10, 1e17, False
        c, d, ok = check_distance(r, l1, (1 - eta) * T[i + 1], MU)
        if not ok:
            return 1e10, 1e10, 1e17, False
        if c:
            iq = max(iq, d)
        DV += _norm(l1 - v)
        vend = l2
        rP = rN; vP = vN
    return DV, T.sum(), iq, True


# ---------- Funksione per nje kembe te vetme (perdoren nga kerkimi) ----------

@njit(cache=True, error_model='numpy')
def entry_leg(u, vv, T0, p):
    """Kemba e hyrjes: nga pika (u, vv) ne sferen 10 AU te planeti p pas T0 ditesh.
    Kthen (DV_init = | |v1| - 10 km/s |, shpejtesia e mberritjes, ok)."""
    theta = 2 * np.pi * u
    phi = np.arccos(2 * vv - 1) - np.pi / 2
    rinit = np.array([R_START * np.cos(phi) * np.cos(theta), R_START * np.cos(phi) * np.sin(theta),
                      R_START * np.sin(phi)])
    rP, vP = eph(p, T0)
    v1, v2, ok = lambert0(rinit, rP, T0 * DAY2SEC, MU)
    if not ok:
        return 1e10, v2, False
    c, d, ok2 = check_distance(rinit, v1, T0, MU)
    if (not ok2) or (c and d > 0):
        return 1e10, v2, False
    return abs(_norm(v1) - V_START), v2, True


@njit(cache=True, error_model='numpy')
def leg(p, q, t0, vin, beta, rpr, eta, T):
    """Nje kembe planet-planet: fly-by te p (beta, rp), udhetim pa motor per eta*T dite,
    pastaj manovra DSM dhe hark Lambert deri te q per (1-eta)*T dite.
    Kthen (DV i manovres, shpejtesia e mberritjes te q, ok)."""
    rP, vP = eph(p, t0)
    vout = fb_prop(vin, vP, rpr * EARTH_RADIUS, beta, PMU[p])
    c, d, ok = check_distance(rP, vout, eta * T, MU)
    if (not ok) or (c and d > 0):
        return 1e10, vout, False
    r, v, ok = propagate(rP, vout, eta * T * DAY2SEC, MU)
    if not ok:
        return 1e10, vout, False
    rN, vN = eph(q, t0 + T)
    l1, l2, ok = lambert0(r, rN, (1 - eta) * T * DAY2SEC, MU)
    if not ok:
        return 1e10, l2, False
    c, d, ok = check_distance(r, l1, (1 - eta) * T, MU)
    if (not ok) or (c and d > 0):
        return 1e10, l2, False
    return _norm(l1 - v), l2, True
