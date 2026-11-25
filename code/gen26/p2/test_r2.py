# Own check of R2 (stratified grouping): exact Doob martingale on a tiny instance, decomposition, linear-part bound.
import numpy as np, itertools, math
rng = np.random.default_rng(7)
def mknet(g, P, kind):
    if kind == 'unif': c = rng.random((g, g))
    elif kind == 'hetero': u = rng.random(g)**4; c = np.outer(u, u)*rng.random((g, g))
    elif kind == 'star': c = np.zeros((g, g)); hub = rng.choice(g, max(1, g//8), replace=False); c[hub, :] = 1; c[:, hub] = 1
    c = np.triu(c, 1); c = c + c.T
    c = c/ c.max() * P/g
    # enforce node sums <= P
    s = c.sum(1).max()
    if s > P: c *= P/s
    return c
def Yval(c, groups_a, groups_b):  # groups: list of node ids
    return c[np.ix_(groups_a, groups_b)].sum()
# tiny exact: K=3, t=3
K, t = 3, 3; g = K*t; P = 1.0
for kind in ('unif', 'hetero', 'star'):
    c = mknet(g, P, kind); ca = c.sum(1); order = np.argsort(ca)
    strata = [order[i*K:(i+1)*K] for i in range(t)]
    perms = list(itertools.permutations(range(K)))
    gam, gam2 = 0, 1
    # all bijections: pi_i maps stratum position -> group; group gamma gets strata[i][p.index(gamma)]
    def grp(pis, gm): return [strata[i][pis[i].index(gm)] for i in range(len(pis))]
    def condE(prefix):
        tot = 0; cnt = 0
        for rest in itertools.product(perms, repeat=t-len(prefix)):
            pis = list(prefix) + list(rest)
            tot += Yval(c, grp(pis, gam), grp(pis, gam2)); cnt += 1
        return tot/cnt
    maxerr = 0; maxlinsum = 0; maxrest = 0
    for trial in range(20):
        pis = [perms[rng.integers(len(perms))] for _ in range(t)]
        M = [condE(pis[:k]) for k in range(t+1)]
        linsum = 0
        for k in range(1, t+1):   # stratum index k-1 (0-based)
            kk = k-1; S = strata[kk]
            a = strata[kk][pis[kk].index(gam)]; b = strata[kk][pis[kk].index(gam2)]
            Gp = [strata[j][pis[j].index(gam2)] for j in range(kk)]; G = [strata[j][pis[j].index(gam)] for j in range(kk)]
            later = [v for j in range(kk+1, t) for v in strata[j]]
            h = lambda v, Gx: c[v, Gx].sum() + c[v, later].sum()/K
            Dk = (h(a, Gp) - np.mean([h(v, Gp) for v in S])) + (h(b, G) - np.mean([h(v, G) for v in S])) \
                 + (c[a, b] - np.mean([c[u, v] for u in S for v in S if u != v]))
            maxerr = max(maxerr, abs(Dk - (M[k]-M[k-1])))
            lin = lambda v: (ca[v] - c[v, [w for w in S if w != v]].sum())/K
            la = lin(a) - np.mean([lin(v) for v in S]); lb = lin(b) - np.mean([lin(v) for v in S])
            linsum += abs(la) + abs(lb)
        maxlinsum = max(maxlinsum, linsum)
    print(f"{kind}: max |D_k - (M_k - M_(k-1))| = {maxerr:.2e};  max sum|lin| = {maxlinsum:.4f} vs bound 4P/K = {4*P/K:.4f}; "
          f"E Y={M[0]:.4f} <= U/(K(K-1)) = {c.sum()/(K*(K-1)):.4f}")
