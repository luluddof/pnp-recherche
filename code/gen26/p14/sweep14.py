# Plain clique sweep with per-edge widths: cost(lam) = min over orderings max_j sum_{a<=j<=b} lam_{o_a o_b}.
# Claim (P13 Step D): cost <= e(T + q D), q = 1, e = s/K^2, T = K^2 c', D = min_c (2K/(K-1)) sum|lam - c|.
# q_needed = (cost - s c') / (beta * min_c sum|lam-c|), beta = 2s/(K(K-1)). Hill-climb lam to maximise it.
import numpy as np, itertools, math
rng = np.random.default_rng(7)
def run(K, iters):
    s = (K+1)**2//4 - 1; beta = 2*s/(K*(K-1))
    E = [(a, b) for a in range(K) for b in range(a+1, K)]; M = len(E)
    orders = np.array(list(itertools.permutations(range(K))))
    pos = np.argsort(orders, axis=1)                         # pos[o, node]
    # incidence: for each ordering and j, which edges are in cut_j
    pa = pos[:, [a for a, b in E]]; pb = pos[:, [b for a, b in E]]
    lo = np.minimum(pa, pb); hi = np.maximum(pa, pb)
    inc = np.stack([(lo <= j) & (j <= hi) for j in range(K)], axis=1).astype(float)   # (n_ord, K, M)
    def q(lam):
        cost = (inc @ lam).max(axis=1).min()
        cbar = lam.mean(); med = np.median(lam); dev = np.abs(lam - med).sum()
        return (cost - s*cbar)/(beta*dev) if dev > 1e-12 else -np.inf
    best = -np.inf
    for restart in range(4):
        lam = rng.random(M) if restart else np.zeros(M)
        if restart == 0: lam[:K-1] = 1.0
        cur = q(lam)
        for it in range(iters):
            l2 = lam.copy(); i = rng.integers(M)
            l2[i] = max(0.0, l2[i] + rng.normal()*0.5) if rng.random() < 0.7 else float(rng.integers(2))
            v = q(l2)
            if v >= cur: lam, cur = l2, v
        best = max(best, cur)
    # explicit star
    lam = np.array([1.0 if 0 in e else 0.0 for e in E]); qs = q(lam)
    print(f"K={K}: s={s}, hill-climb max q_needed = {best:.4f}; star q = {qs:.4f} (formula K(K-1)/(2s)-1 = {K*(K-1)/(2*s)-1:.4f})", flush=True)
for K, it in ((4, 3000), (5, 3000), (6, 2000), (7, 600)):
    run(K, it)
