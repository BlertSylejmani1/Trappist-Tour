"""Kerkimi i gjate (p.sh. 12-24 ore): perserit beam search me konfigurime te ndryshme dhe,
here pas here, rioptimizimin e fundit te turave me te vlefshme. Pas cdo run-i:
  - rezultati ruhet ne runs/run_XXX.npy,
  - arkivi bashkohet dhe pastrohet (vetem fronti i vlefshem zyrtarisht),
  - nese score-i rritet, ruhet best_found.json dhe shenohet ne log_search.txt.
Nese nderpritet, nuk humbet asgje; me te njejten komande vazhdon nga run-i i fundit.

    python run_search.py --hours 20
"""
import os, sys, time, json, glob, argparse, warnings
import numpy as np

warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE); sys.path.insert(0, HERE)

# Konfigurimet qe perseriten ne cikel:
#   ('beam', B, hyrje per proces, planetet e para, keep, nstart, pjesa e heuristikes)
#   ('suffix', numri i turave, prerjet, B i vogel, keep, nstart, -)
CONFIGS = [
    ('beam', 400, 1500, [0], 2, 4, 0.6),
    ('beam', 600, 2500, [0], 2, 4, 0.8),
    ('beam', 400, 1500, [0, 1], 2, 4, 0.6),
    ('suffix', 10, (1, 2, 3), 40, 3, 6, None),
    ('beam', 500, 2000, [0], 3, 3, 0.7),
    ('beam', 300, 1200, [1, 5], 2, 4, 0.6),
    ('beam', 800, 3000, [0], 2, 3, 0.6),
    ('suffix', 10, (1, 2, 3), 40, 3, 6, None),
]


def log(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S ') + msg
    print(line, flush=True)
    with open('log_search.txt', 'a', encoding='utf-8') as f:
        f.write(line + '\n')


def shrink(archive, n=600):
    """Mban vetem frontin e vlefshem zyrtarisht (maksimumi n pika)."""
    import hypervolume as H
    idx, Ro = H.valid_front(archive)
    if len(idx) > n:
        idx = idx[H.select_exact(Ro, n)]
    return archive[idx]


def most_valuable(archive, n):
    """Turat me kontributin me te madh ne hypervolume (per rioptimizimin e fundit)."""
    import hypervolume as H, pygmo as pg
    R = H.eval_all(archive)
    ok = np.where((R[:, 2] <= 0) & (R[:, 0] < 2500) & (R[:, 1] < 4000))[0]
    c = pg.hypervolume(R[ok, :2]).contributions(H.REF)
    return archive[ok[np.argsort(-c)[:n]]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--hours', type=float, default=1e9)
    ap.add_argument('--procs', type=int, default=max(1, (os.cpu_count() or 2) // 2))
    a = ap.parse_args()
    t_end = time.time() + a.hours * 3600
    import beam_search as BS, suffix_reopt as SR, hypervolume as H
    from multiprocessing import Pool
    os.makedirs('runs', exist_ok=True)

    # arkivi fillestar: zgjidhjet e meparshme + te gjitha run-et e ruajtura
    parts = [np.load('seed_archive.npy')]
    sub = os.path.join('..', 'submission', 'best_submission.json')
    if os.path.exists(sub):
        parts.append(np.array(json.load(open(sub))[0]['decisionVector']))
    parts += [np.load(f) for f in sorted(glob.glob(os.path.join('runs', 'run_*.npy')))]
    archive = shrink(np.vstack(parts))
    Xb, best = H.best_subset(archive)
    json.dump(H.to_submission(Xb), open('best_found.json', 'w'))
    log(f'start: archive {len(archive)} tours, best score {best:,.0f}, processes {a.procs}')

    k = len(glob.glob(os.path.join('runs', 'run_*.npy')))
    while time.time() < t_end:
        cfg = CONFIGS[k % len(CONFIGS)]
        t0 = time.time()
        if cfg[0] == 'beam':
            _, B, ne, firsts, keep, nstart, fh = cfg
            log(f'run {k}: heuristic beam B={B} entries={ne} first planets={firsts} keep={keep} heur={fh} ...')
            X = BS.run(B=B, n_entry=ne, nstart=nstart, seed=2000 + k, procs=a.procs, firsts=firsts,
                       keep=keep, frac_heur=fh)
        else:
            _, ntours, cuts, Bm, keep, nstart, _ = cfg
            log(f'run {k}: suffix re-optimisation of the {ntours} most valuable tours, cuts {cuts} ...')
            with Pool(a.procs) as pool:
                X = SR.run_suffix(pool, most_valuable(archive, ntours), cuts=cuts, Bm=Bm,
                                  nstart=nstart, keep=keep, seed=3000 + k)
        np.save(os.path.join('runs', f'run_{k:03d}.npy'), X)
        archive = shrink(np.vstack([archive, X]))
        Xb, score = H.best_subset(archive)
        dt = (time.time() - t0) / 60
        if score > best + 1e-6:
            best = score
            json.dump(H.to_submission(Xb), open('best_found.json', 'w'))
            log(f'run {k} done in {dt:.0f} min  ->  NEW BEST {best:,.0f} saved')
        else:
            log(f'run {k} done in {dt:.0f} min  ->  score {score:,.0f} (best stays {best:,.0f})')
        k += 1
    log(f'finished. best score {best:,.0f}, file best_found.json (then run merge_final.py)')


if __name__ == '__main__':
    main()
