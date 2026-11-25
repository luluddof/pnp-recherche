# hill-climb the network c (entries in [0,1]) to maximise Var(Y12-Y34)/(Var Y12+Var Y34) exactly (K=5,t=2 and K=6,t=2)
import numpy as np, sys
from enum14 import stats
rng = np.random.default_rng(3)
for K, t, it in ((5, 2, 1500), (6, 2, 120), (8, 1, 3000)):
    g = K*t; eps = (4*K-2)/((K-2)*(K-3))
    if t == 1 and K == 8:
        pass
    bestall = 0
    for restart in range(3):
        c = rng.random((g, g)); c = np.triu(c, 1); c = c + c.T
        cur = stats(c, K, t)[0]
        for i in range(it):
            a, b = rng.integers(g, size=2)
            if a == b: continue
            c2 = c.copy(); v = rng.random() if rng.random() < 0.5 else float(rng.integers(2)); c2[a, b] = c2[b, a] = v
            r = stats(c2, K, t)[0]
            if r > cur: c, cur = c2, r
        bestall = max(bestall, cur)
    print(f"K={K} t={t}: hill-climb max ratio {bestall:.4f}; 1+2/(K-2) = {1+2/(K-2):.4f}; 1+eps_K = {1+eps:.4f}", flush=True)
