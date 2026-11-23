"""Kontrollon nje skedar dorezimi me kodin ZYRTAR te ESA-s (official.py, pykep 2.6):
cdo ture duhet te kete eq == 0, ineq <= 0, DV < 2500 m/s dhe T < 4000 dite.
Ne fund printon score-in me funksionin zyrtar combine_scores.

    python verify_submission.py ../submission/best_submission.json
"""
import sys, json, io, contextlib, warnings
warnings.filterwarnings('ignore')
import official as O

path = sys.argv[1] if len(sys.argv) > 1 else '../submission/best_submission.json'
sub = json.load(open(path))[0]
assert sub['problem'] == 'trappist-tour' and sub['challenge'] == 'spoc-trappist-tour'
X = sub['decisionVector']
assert len(X) <= 100, 'maksimumi 100 tura'
lb, ub = O.udp.get_bounds()
pts, bad = [], 0
for k, x in enumerate(X):
    in_bounds = len(x) == 34 and all(l - 1e-12 <= v <= u + 1e-12 for v, l, u in zip(x, lb, ub)) \
        and sorted(int(v) for v in x[27:]) == list(range(7))
    with contextlib.redirect_stdout(io.StringIO()):
        dv, T, eq, iq = O.udp.fitness(x)
    ok = in_bounds and eq == 0 and iq <= 0 and dv < 2500 and T < 4000
    bad += not ok
    pts.append([dv, T])
    print(f'{k:3d}  DV = {dv:8.2f} m/s   T = {T:8.2f} dite   {"OK" if ok else "GABIM"}')
print(f'\nTura: {len(X)}, te pavlefshme: {bad}')
print(f'Score zyrtar (combine_scores): {O.combine_scores(pts):,.1f}  ->  hypervolume {-O.combine_scores(pts):,.1f}')
