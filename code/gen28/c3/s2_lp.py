# Non-symmetric search in the site-pattern model (n=2): degree-d g on {0,1}^k, measure mu_{1/4}; maximise agreement with parity.
# Compare: symmetric window optimum, Theorem S'' bound.
import numpy as np, itertools, sys
from math import comb
from scipy.optimize import linprog
def bound(k, d, b=0.25):
    s = 0.0
    for m in range(k+1):
        pm = comb(k, m) * (2*b)**m * (1-2*b)**(k-m)
        s += pm * min(1.0, sum(comb(m, i) for i in range((m+d)//2 + 1)) / 2**m)
    return s
def sym_best(k, d, b=0.25):
    pm = [comb(k, m) * b**m * (1-b)**(k-m) for m in range(k+1)]
    best = max(sum(pm[a:a+d+1]) for a in range(k+1))
    return max(best, sum(pm[0::2]))
def lp_search(k, d, iters=8, seed=0):
    rng = np.random.default_rng(seed)
    pats = np.array(list(itertools.product([0,1], repeat=k)))
    mu = np.prod(np.where(pats == 1, 0.25, 0.75), axis=1)
    par = (-1.0)**pats.sum(1)
    mons = [U for r in range(d+1) for U in itertools.combinations(range(k), r)]
    A = np.array([[np.prod(T[list(U)]) if U else 1.0 for U in mons] for T in pats])
    P, M = A.shape
    best = 0
    # init from best symmetric window interpolant
    pmw = [sum(mu[pats.sum(1)==m]) for m in range(k+1)]
    a0 = max(range(k+1), key=lambda a: sum(pmw[a:a+d+1]))
    win = list(range(a0, min(k, a0+d)+1))
    ms = pats.sum(1)
    from numpy.polynomial import polynomial as Pn
    coef = np.polyfit(win, [(-1.0)**m for m in win], len(win)-1) if len(win)>1 else np.array([(-1.0)**win[0]])
    gs = np.polyval(coef, ms)
    res0 = np.abs(gs - par); best = mu[res0 < 1e-7].sum()
    wt = mu / (res0 + 1e-3)
    for it in range(iters):
        # min sum wt*(s+ + s-), A a + s+ - s- = par
        c = np.concatenate([np.zeros(M), wt, wt])
        Aeq = np.hstack([A, np.eye(P), -np.eye(P)])
        bnds = [(None, None)]*M + [(0, None)]*(2*P)
        r = linprog(c, A_eq=Aeq, b_eq=par, bounds=bnds, method='highs')
        a = r.x[:M]; res = np.abs(A @ a - par)
        agree = mu[res < 1e-7].sum()
        best = max(best, agree)
        wt = mu / (res + 1e-3) * (1 + 0.05*rng.random(P))
    return best
if __name__ == '__main__':
    for k in map(int, sys.argv[1:]):
        for d in range(0, k+1):
            sb = sym_best(k, d); bd = bound(k, d)
            lp = lp_search(k, d) if (k <= 12 and comb(k, d) < 2000) else float('nan')
            print(f"k={k} d={d} sym={sb:.4f} lp={lp:.4f} bound={bd:.4f}  q_max={sum(comb(k,j)*3**j for j in range(d+1))}", flush=True)
