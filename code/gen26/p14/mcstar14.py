# Monte Carlo (stratified grouping) for the star network at larger K: ratio vs 1+2/(K-2) and 1+eps_K
import numpy as np
rng = np.random.default_rng(11)
for K, t, reps in ((8, 3, 40000), (16, 3, 40000), (32, 2, 40000)):
    g = K*t; c = np.zeros((g, g)); v = 0; c[v, :] = 1; c[:, v] = 1; c[v, v] = 0
    Y = np.zeros((reps, 2))
    for r in range(reps):
        G = np.array([rng.permutation(K) for _ in range(t)]).T + np.arange(t)*K   # G[group] = node indices
        Y[r] = (c[np.ix_(G[0], G[1])].sum(), c[np.ix_(G[2], G[3])].sum())
    rv = (Y[:, 0]-Y[:, 1]).var()/(Y[:, 0].var()+Y[:, 1].var()); eps = (4*K-2)/((K-2)*(K-3))
    print(f"star K={K} t={t}: ratio {rv:.4f}  1+2/(K-2)={1+2/(K-2):.4f}  1+eps={1+eps:.4f}", flush=True)
