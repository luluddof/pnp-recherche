"""Exact (integer/float) checks of the two identities used in Theorem F, on the P of hexlib (N=2,3)."""
import numpy as np, itertools
from hexlib import perm_P
rng = np.random.default_rng(0)
for N in (2, 3):
    d = N*N; P = perm_P(N); P6 = P.reshape((d,)*6)            # [o0,o1,o2,i0,i1,i2]
    assert np.allclose(P @ P.T, np.eye(d**3))
    # (i) Schmidt rank across site 0 | sites 12 of P v  ==  Schmidt rank of v across (strands moving INTO site 0 | rest)
    #     input strands: site0=(a,a2) site1=(b,b2) site2=(c,c2); out site0 = (b,c)
    for s in (1, 2, 3):
        Phi = rng.normal(size=(d, s))                         # s-dim subspace of site 0
        v = np.einsum('ps,sjk->pjk', Phi, rng.normal(size=(s, d, d)))     # v in Phi_s (x) H_12
        Pv = (P @ v.reshape(-1)).reshape(d, d, d)
        sr_out = np.linalg.matrix_rank(Pv.reshape(d, d*d))
        v6 = v.reshape(N, N, N, N, N, N)                      # a,a2,b,b2,c,c2
        sr_in = np.linalg.matrix_rank(v6.transpose(2, 4, 0, 1, 3, 5).reshape(N*N, -1))   # (b,c) | rest
        assert sr_out == sr_in, (sr_out, sr_in)
        print(f"N={N} s={s}: SR_0|12(Pv)={sr_out} = SR_(b,c)|rest(v)={sr_in}  (generic value min(N^2, s N^2)={min(d, s*d)})")
    # (ii) Lambda_P(w) = Tr_0[(w^dagger (x) 1) P] has rank = rank(w) * N^2
    for rho in (1, 2, 3):
        w = rng.normal(size=(d, rho)) @ rng.normal(size=(rho, d))
        Lam = np.einsum('pq,pjkqmn->jkmn', w, P6).reshape(d*d, d*d)
        print(f"N={N} rank(w)={rho}: rank Lambda_P(w) = {np.linalg.matrix_rank(Lam)} (= rho*N^2 = {rho*d})")
        assert np.linalg.matrix_rank(Lam) == rho*d
# (iii) the proven lower bound r_min(b) = max_s s Q (M - b s) for the hexagon, vs known points
def rmin(N, b):
    M = N*N; Q = N*N
    cand = [s*Q*(M - b*s) for s in range(1, M+1)] + [Q*(rho*M - b)/rho**2 for rho in range(1, M+1)]
    return max(0, int(np.ceil(max(cand) - 1e-9)))
print("hexagon N=2: b -> proven min r:", {b: rmin(2, b) for b in range(1, 17)})
print("hexagon N=4: b -> proven min r:", {b: rmin(4, b) for b in (1, 2, 4, 8, 16, 32, 64, 128, 200, 255, 256)})
print("known points N=2 (r;b): (28;1) ok" , 28 >= rmin(2, 1), " (10;7)", 10 >= rmin(2, 7), " (1;16)", 1 >= rmin(2, 16), "  [B=identity: r=(N^6-N^3)/2]")
# crossing: r >= max_rho (rho N - b)/rho^2
def rmin_cross(N, b): return max(0, int(np.ceil(max((rho*N - b)/rho**2 for rho in range(1, N+1)) - 1e-9)))
print("crossing: N=2 b=1 ->", rmin_cross(2, 1), "(point (1,1) exists: tight);  N=3:", [rmin_cross(3, b) for b in (1, 2, 3)],
      " N=4:", [rmin_cross(4, b) for b in range(1, 9)], " N=8 b=8:", rmin_cross(8, 8))
