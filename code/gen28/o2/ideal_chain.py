# "Ideal" merged chain formula: pieces of a chain (0,N^2)=(p_0,q_0) < (p_1,q_1) < ... < (N^2,0), merged type classes of the
# k-th power if the row (column) operators M_s were commuting projections of rank d_s = p_s - p_{s-1} (e_s = q_{s-1} - q_s):
#   h_id(mu) = log_{N^2} sum_s d_s^mu e_s^(1-mu).   For ONE point this is exactly h_LB (proved lower bound, attained).
# Compare with circle; also O1 optimistic-entropy (p_s, q_{s-1}) and multiplicity-free (max).
import itertools, math, numpy as np
pts = {2: [(1,1),(1,2)], 3: [(3,3),(3,1)], 4: [(4,4),(5,3),(6,1)], 5: [(6,6),(7,4),(10,1)],
 6: [(7,7),(8,6),(11,5),(12,3),(15,1)], 7: [(8,8),(11,7),(12,6),(14,5),(17,4),(20,3),(21,1)],
 8: [(8,8),(19,6),(20,6),(22,5),(28,1)], 9: [(17,9),(22,6),(36,1)], 10: [(25,5),(45,1)], 12: [(24,18),(66,1)],
 16: [(120,1)], 32: [(496,1)], 64: [(2016,1)]}
mus = np.linspace(0, 1, 2001); hB = np.hypot(mus, 1-mus); inn = slice(1, -1); W = mus[inn]*(1-mus[inn])
def h(N, chain, kind):
    N2 = N*N; full = [(0, N2)] + list(chain) + [(N2, 0)]
    if kind == 'ideal': rs = [b[0]-a[0] for a, b in zip(full, full[1:])]; cs = [a[1]-b[1] for a, b in zip(full, full[1:])]
    else: rs = [max(b[0], 1) for a, b in zip(full, full[1:])]; cs = [max(a[1], 1) for a, b in zip(full, full[1:])]
    rs = np.array(rs, float); cs = np.array(cs, float)
    if (rs <= 0).any() or (cs <= 0).any(): return None
    T = rs[None]**mus[:, None]*cs[None]**(1-mus[:, None])
    return np.log(T.sum(1) if kind != 'max' else T.max(1))/math.log(N2)
for N, P in pts.items():
    allp = sorted(set(P + [(q, p) for p, q in P]))
    allp = [a for a in allp if not any(b != a and b[0] <= a[0] and b[1] <= a[1] for b in allp)]; allp.sort()
    res = {}
    for kind in ('ideal', 'O1', 'max'):
        best = None
        for k in range(1, len(allp)+1):
            for ch in itertools.combinations(allp, k):
                if any(b[1] >= a[1] for a, b in zip(ch, ch[1:])): continue
                hh = h(N, ch, kind)
                if hh is None: continue
                g = ((hh - hB)[inn]/W).min(); 
                if best is None or g < best[0]: best = (g, ch, hh[1000], mus[inn][((hh-hB)[inn]/W).argmin()])
        res[kind] = best
    print(f"N={N:2d} " + " | ".join(f"{k}: min (h-circ)/(mu(1-mu)) {v[0]:+.4f} at mu={v[3]:.3f}, h(1/2) {v[2]:.4f}, chain {v[1]}" for k, v in res.items()))
