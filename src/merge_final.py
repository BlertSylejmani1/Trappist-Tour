"""Krijon dorezimin final nga TE GJITHA zgjidhjet e ruajtura:
  1) bashkon arkivin fillestar, run-et (runs/*.npy) dhe dorezimet ekzistuese,
  2) mban frontin Pareto qe eshte i vlefshem sipas kodit zyrtar,
  3) zgjedh <= 100 pikat me hypervolume maksimal (programim dinamik i sakte),
  4) e ruan ne ../submission/best_submission.json vetem nese eshte me i mire se ai ekzistues.

    python merge_final.py
"""
import os, sys, glob, json, warnings
import numpy as np

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE); sys.path.insert(0, HERE)
import hypervolume as H

OUT = os.path.join('..', 'submission', 'best_submission.json')

parts = [np.load('seed_archive.npy')]
parts += [np.load(f) for f in sorted(glob.glob(os.path.join('runs', 'run_*.npy')))]
for f in [OUT, 'best_found.json']:
    if os.path.exists(f):
        parts.append(np.array(json.load(open(f))[0]['decisionVector']))
X = np.vstack(parts)
print('tours in total:', len(X), flush=True)

Xb, score = H.best_subset(X)
if os.path.exists(OUT):
    old = H.official_score([H.official_eval(x)[:2] for x in json.load(open(OUT))[0]['decisionVector']])
    if old >= score - 1e-6:
        print(f'existing submission is already as good ({old:,.1f}); not overwritten')
        sys.exit(0)
json.dump(H.to_submission(Xb), open(OUT, 'w'))
print(f'OFFICIAL score: {score:,.1f}  ({len(Xb)} tours)  ->  {OUT}')
