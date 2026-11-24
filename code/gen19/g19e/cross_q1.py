# Exhaustive: min rank(F - G (x) G') over all G, G' in M_N(F_p)  (crossing, q = 1), vs Theorem 2 bound.
# Reduction: (g(x)h) F (h^-1 (x) g^-1) = F, so G may be taken in rank normal form diag(1_rho,0).
import sys, itertools, numpy as np
sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from ff import brank
N, p = int(sys.argv[1]), int(sys.argv[2])
n2 = N*N
F = np.zeros((n2, n2), dtype=np.int64)   # rows (x,y), cols (u,v): [x=v][y=u]
for x in range(N):
    for y in range(N):
        F[x*N+y, y*N+x] = 1
def allmats(chunk=200000):
    tot = p**n2
    for s in range(0, tot, chunk):
        ids = np.arange(s, min(tot, s+chunk))
        digs = np.stack([(ids // p**i) % p for i in range(n2)], axis=1)
        yield digs.reshape(-1, N, N)
bound = max(k*(N-k*1) for k in range(N+1))
best = {}
for rho in range(1, N+1):
    G = np.zeros((N, N), dtype=np.int64); G[:rho,:rho] = np.eye(rho, dtype=np.int64)
    mn = n2; hist = {}
    for Gp in allmats():
        # (G (x) G')[(x,y),(u,v)] = G[x,u] G'[y,v]
        K = np.einsum('xu,byv->bxyuv', G, Gp).reshape(-1, n2, n2)
        r = brank(F[None] - K, p)
        u, c = np.unique(r, return_counts=True)
        for a, b in zip(u, c): hist[int(a)] = hist.get(int(a), 0) + int(b)
    best[rho] = min(hist)
    print(f"N={N} p={p} rank(G)={rho}: min rank(F-G(x)G') = {min(hist)}  (count {hist[min(hist)]})", flush=True)
print(f"N={N} p={p}: exact min p at q=1 is {min(best.values())}; Theorem-2 bound {bound}; N(N-1)/2 = {N*(N-1)//2}")
assert min(best.values()) >= bound
