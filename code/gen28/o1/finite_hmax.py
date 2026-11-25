# Task 3: can any known finite (p,q)@N point set, used as a tensor-powered piece chain (entropy formula
#   h(mu) = log_{N^2} sum_s r_s^mu c_s^(1-mu), B2/finite_chain), beat App-B's circle h_B(mu) = ||(mu,1-mu)||_2 at ANY mu?
# If h_fin >= h_B for all mu, no level (balanced or unbalanced, any e_a/e_b) can gain: e' = min_mu h(mu)/(mu/e_a+(1-mu)/e_b).
import itertools, math, numpy as np
pts = {  # realizable (p,q)@N (exact/upper points from G26-P11, G27-T2, BRIEF); symmetric (q,p) added below
 2: [(1,1),(1,2)], 3: [(3,3),(3,1)], 4: [(4,4),(5,3),(6,1)], 5: [(6,6),(7,4),(10,1)],
 6: [(7,7),(8,6),(11,5),(12,3),(15,1)], 7: [(8,8),(11,7),(12,6),(14,5),(17,4),(20,3),(21,1)],
 8: [(8,8),(19,6),(20,6),(22,5),(28,1)], 9: [(17,9),(22,6),(36,1)], 10: [(25,5),(45,1)], 12: [(24,18),(66,1)],
 16: [(120,1)], 32: [(496,1)], 64: [(2016,1)]}
mus = np.linspace(0, 1, 2001); hB = np.hypot(mus, 1 - mus)
def hchain(N, chain, mode):
    N2 = N*N; full = [(0, N2)] + chain + [(N2, 0)]
    rs, cs = [], []
    for (p0, q0), (p1, q1) in zip(full[:-1], full[1:]):
        if mode == 'safe': r, c = min(N2, p0 + p1), min(N2, q0 + q1)
        else: r, c = max(p1, 1), max(q0, 1)          # optimistic: nested chain, piece = (p_s, q_{s-1})
        rs.append(r); cs.append(c)
    rs = np.array(rs, float); cs = np.array(cs, float)
    S = (rs[None, :]**mus[:, None] * cs[None, :]**(1 - mus[:, None])).max(1)
    return np.log(S)/math.log(N2)
best = {}
for N, P in pts.items():
    allp = sorted(set(P + [(q, p) for p, q in P]))
    # Pareto-filter: keep points not dominated
    allp = [a for a in allp if not any(b != a and b[0] <= a[0] and b[1] <= a[1] for b in allp)]
    allp.sort()
    for mode in ('opt', 'safe'):
        gmin = None
        for k in range(1, len(allp) + 1):
            for ch in itertools.combinations(allp, k):
                h = hchain(N, list(ch), mode); gap = (h - hB).min()
                if gmin is None or gap < gmin[0]: gmin = (gap, ch, mus[(h - hB).argmin()], h[1000])
        print(f"N={N:2d} {mode:4s}: min_mu [h_fin - h_B] = {gmin[0]:+.4f} at mu={gmin[2]:.3f} chain={gmin[1]}  (h(1/2)={gmin[3]:.4f} vs 0.7071)")
