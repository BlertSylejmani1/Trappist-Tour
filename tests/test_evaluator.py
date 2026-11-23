"""Teste: (1) vleresuesi i shpejte perputhet me kodin zyrtar; (2) dorezimi eshte i vlefshem.

    cd src
    python ../tests/test_evaluator.py
"""
import os, sys, io, json, contextlib, warnings
import numpy as np
warnings.filterwarnings('ignore')
SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'src')
sys.path.insert(0, SRC)
import official as O
import evaluator as F


def test_evaluator_matches_official(n=600, seed=0):
    """DV dhe T te vleresuesit te shpejte duhet te perputhen me kodin zyrtar (pykep).
    Testohet ne zonen qe ka rendesi: tura te mira nga dorezimi, me zhurme te rastit."""
    rng = np.random.default_rng(seed)
    lb, ub = map(np.array, O.udp.get_bounds())
    base = np.array(json.load(open(os.path.join(SRC, '..', 'submission', 'best_submission.json')))[0]['decisionVector'])
    worst = 0.0; compared = 0
    for _ in range(n):
        x = base[rng.integers(len(base))].copy()
        x[:27] += rng.normal(0, 1, 27) * np.r_[1e-6, 1e-6, 0.01, np.tile([0.01, 0.01, 0.002, 0.01], 6)]
        x[:27] = np.clip(x[:27], lb[:27], ub[:27])
        with contextlib.redirect_stdout(io.StringIO()):
            try:
                f = O.udp.fitness(list(x))
            except Exception:
                continue
        if len(f) != 4 or f[0] > 1e4:          # zonat me DV shume te madh nuk kane rendesi
            continue
        g = F.fitness(x)
        assert g[3], 'vleresuesi i shpejte deshtoi aty ku pykep punon'
        worst = max(worst, abs(f[0] - g[0]) / max(1.0, f[0]), abs(f[1] - g[1]))
        compared += 1
    assert compared > 100 and worst < 1e-6, (compared, worst)
    print(f'OK: {compared} vektore, gabimi maksimal relativ {worst:.1e}')


def test_submission_is_valid():
    path = os.path.join(SRC, '..', 'submission', 'best_submission.json')
    X = json.load(open(path))[0]['decisionVector']
    pts = []
    for x in X:
        with contextlib.redirect_stdout(io.StringIO()):
            dv, T, eq, iq = O.udp.fitness(x)
        assert eq == 0 and iq <= 0 and dv < 2500 and T < 4000
        pts.append([dv, T])
    score = -O.combine_scores(pts)
    assert score > 5.8e6
    print(f'OK: {len(X)} tura te vlefshme, hypervolume {score:,.1f}')


if __name__ == '__main__':
    test_evaluator_matches_official()
    test_submission_is_valid()
