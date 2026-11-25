# exponent/n at the top (mythr13.top_cost, r=4.26, III=1, S2q) for each solver table, at given log2 m; ratio to balanced.
import os, sys, json, math, numpy as np
os.environ['III'] = '1'
sys.path.insert(0, '/home/user/pnp-recherche/code/gen26/p13')
src = open('/home/user/pnp-recherche/code/gen26/p13/mythr13.py').read().split('if __name__')[0]
exec(src)
def table(t):
    d = json.load(open(f'sol_{t}.json')); X = np.array(d['x']); E = np.array(d['e']); Q = np.array(d['q']); K = np.array(d['kind'])
    o = np.argsort(X); X, E, Q, K = X[o], E[o], Q[o], K[o]; keep = []; em = np.inf
    for i in range(len(X)):
        if E[i] < em: keep.append(i); em = E[i]
    return X[keep], E[keep], Q[keep], K[keep]
def f(sol, lm, r=4.26):
    best = None
    for x, e, q, kd in zip(*sol):
        if 2*x + 2 > lm: break
        c = top_cost(x, e, lm, r, 'S2q', qk=(q, kd))
        if c and (best is None or c[0] < best): best = c[0]
    return best
tabs = {t: table(t) for t in sys.argv[1:]}
LMS = [100, 200, 400, 800, 1600, 3200, 6400, 12800, 25600]
base = {lm: f(tabs['cinf'], lm) for lm in LMS}
print('log2m: ' + ' '.join(f'{lm:>9d}' for lm in LMS))
print('bal exp/n: ' + ' '.join(f'{base[lm]:9.4f}' for lm in LMS))
for t in sys.argv[1:]:
    if t == 'cinf': continue
    v = [f(tabs[t], lm) for lm in LMS]
    print(f'{t:8s} exp/n: ' + ' '.join(f'{a:9.4f}' for a in v) + ' | bal/corner: ' + ' '.join(f'{base[lm]/a:.2f}' for lm, a in zip(LMS, v)), flush=True)
