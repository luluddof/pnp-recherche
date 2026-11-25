# G26-P13 checks of Step E: (1) |E F(a,b)H(c,d)| <= eps_K rms F rms H for a uniform injection (exact, random F,H, worst of many);
# (2) Monte Carlo: Var(Y_e - Y_e') vs (1+eps_K)(Var Y_e + Var Y_e') for stratified grouping on random networks.
import numpy as np, itertools
rng = np.random.default_rng(1)
for K in (5, 6, 8, 12):
    eps = (4*K-2)/((K-2)*(K-3)); worst = 0
    for trial in range(300):
        F = rng.standard_normal((K, K)); H = rng.standard_normal((K, K))
        if trial % 3 == 0: H = F.copy()          # adversarial-ish: correlated
        if trial % 3 == 1: H = F.T.copy()
        np.fill_diagonal(F, 0); np.fill_diagonal(H, 0)
        off = ~np.eye(K, dtype=bool); F[off] -= F[off].mean(); H[off] -= H[off].mean()
        rF = np.sqrt((F[off]**2).mean()); rH = np.sqrt((H[off]**2).mean())
        # exact E over injections
        S = 0.0; n = 0
        for a, b, c, d in itertools.permutations(range(K), 4):
            S += F[a, b]*H[c, d]; n += 1
        worst = max(worst, abs(S/n)/(rF*rH))
    print(f"(1) K={K}: max |E FH|/(rmsF rmsH) = {worst:.4f}  eps_K = {eps:.4f}  ratio {worst/eps:.3f}")
def mc(K, t, kind, reps=20000):
    g = K*t; c = rng.random((g, g)) if kind == 'unif' else (rng.random((g, g)) < 0.3).astype(float)
    if kind == 'rank2':
        e = np.where(np.arange(g) % 2 == 0, 1.0, -1.0); c = 0.5*(1 + np.outer(e, e))
    c = np.triu(c, 1); c = c + c.T
    order = np.argsort(c.sum(1)); strata = order.reshape(t, K)
    Y = np.zeros((reps, 3))
    for r in range(reps):
        groups = np.array([rng.permutation(s) for s in strata])   # groups[k, gamma] = node of stratum k in group gamma
        G = groups.T
        def y(i, j): return c[np.ix_(G[i], G[j])].sum()
        Y[r] = (y(0, 1), y(2, 3), y(0, 2))
    v1, v2 = Y[:, 0].var(), Y[:, 1].var(); vd = (Y[:, 0] - Y[:, 1]).var()
    eps = (4*K-2)/((K-2)*(K-3))
    print(f"(2) K={K} t={t} {kind}: Var(Ye-Ye')/(VarYe+VarYe') = {vd/(v1+v2):.3f} (bound 1+eps = {1+eps:.3f}); "
          f"E|Ye-Ye'| / (E|Ye-EY|+E|Ye'-EY|) = {np.abs(Y[:,0]-Y[:,1]).mean()/(np.abs(Y[:,0]-Y[:,0].mean()).mean()+np.abs(Y[:,1]-Y[:,1].mean()).mean()):.3f}")
for K, t, kind in ((8, 8, 'unif'), (8, 8, 'bern'), (8, 8, 'rank2'), (16, 8, 'bern')):
    mc(K, t, kind, 6000)
