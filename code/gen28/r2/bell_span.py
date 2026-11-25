# Bell model, Pauli support in a subspace S of F_2^{2k} with dim r: optimum p = sum over cosets of S^perp of min(#+,#-) of beta.
# Compare with my bound (N^2 - N 2^r)/2 and C2's N^2 Pr[L1-L2 > 2r] (d=2).
import numpy as np, itertools
from math import comb
rng = np.random.default_rng(4)
def c2bound(k, u):
    pr = {m: (comb(k,m)-(comb(k,m-1) if m else 0))*(k-2*m+1)/2**k for m in range(k//2+1)}
    return 4**k*sum(pr[a]*pr[b] for a in pr for b in pr if b-a > 2*u)
for k in [3, 4, 5]:
    N = 2**k; V = np.array(list(itertools.product([0,1], repeat=2*k)))
    beta = (-1)**((V[:, :k]*V[:, k:]).sum(1) % 2)
    for r in range(1, 2*k):
        worst = None
        for t in range(30):
            B = rng.integers(0, 2, size=(r, 2*k))
            if np.linalg.matrix_rank(B) < r: continue  # (real rank proxy, fine for filtering)
            # phi depends on v only through the r bits <s_j, v>_symp ? Bell: chi_s(v) = (-1)^{<s,v>} (any fixed bilinear pairing);
            key = (V @ B.T) % 2
            keys = key @ (2**np.arange(r))
            p = 0
            for c in np.unique(keys):
                bb = beta[keys == c]; p += min((bb == 1).sum(), (bb == -1).sum())
            worst = p if worst is None else min(worst, p)
        print(f"k={k} r={r}: min p over random S = {worst}, my bound (N^2-N2^r)/2 = {(N*N-N*2**r)/2:.0f}, C2 bound = {c2bound(k,r):.1f}")
