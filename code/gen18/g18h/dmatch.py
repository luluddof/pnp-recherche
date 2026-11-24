"""G18h: D-direction box at lambda=1.  T_eps = prod over the 2^(D-1) antipodal pairs of eps(w_p, w_pbar)
(eps = antisymmetric 2x2; T_eps is the D-box crossing up to an invertible change of basis at half of the
corners, which preserves every flattening rank).  Question: is T_eps a combination of products of eps over
PERFECT MATCHINGS OF THE CUBE GRAPH Q_D by edges?  Such a term has rank 2^(#edges in direction i) in
direction i, total cost = 2^(D-1) bits = floor 1/D exactly.   Exact arithmetic (sympy Rational / mod-p + check)."""
import itertools, sys
import numpy as np
from fractions import Fraction
import sympy

def matchings(pts, allowed):
    if not pts:
        yield (); return
    a = pts[0]
    for b in pts[1:]:
        if allowed(a, b):
            rest = [p for p in pts if p not in (a, b)]
            for m in matchings(rest, allowed):
                yield ((a, b),) + m

def tensor(m, n):
    """prod eps(w_a,w_b) as a vector of length 2^n (integer entries)."""
    T = np.zeros(1 << n, dtype=np.int64)
    k = len(m)
    for bits in range(1 << k):
        idx = 0; sign = 1
        for j, (a, b) in enumerate(m):
            wa = (bits >> j) & 1; wb = 1 - wa
            idx |= (wa << a) | (wb << b)
            if wa: sign = -sign
        T[idx] = sign
    return T

def run(D, maxdist=1):
    n = 1 << D
    full = n - 1
    anti = tuple((p, full ^ p) for p in range(n // 2))
    Teps = tensor(anti, n)
    ham = lambda a, b: bin(a ^ b).count("1")
    ms = list(matchings(list(range(n)), lambda a, b: ham(a, b) <= maxdist))
    costs = sorted(set(sum(ham(a, b) for a, b in m) for m in ms))
    print(f"D={D}: corners {n}, matchings with all pair distances <= {maxdist}: {len(ms)}, total costs {costs}, floor {n//2}, full {D*n//2}")
    A = np.array([tensor(m, n) for m in ms], dtype=np.int64).T   # 2^n x k
    # exact solve over Q on the support rows
    rows = np.nonzero(np.any(A != 0, axis=1) | (Teps != 0))[0]
    M = sympy.Matrix(A[rows].tolist()); b = sympy.Matrix(Teps[rows].tolist())
    aug = M.row_join(b)
    r1 = M.rank(); r2 = aug.rank()
    print(f"   rank of matching tensors = {r1}; T_eps in span: {r1 == r2}")
    if r1 == r2:
        sol, params = M.gauss_jordan_solve(b)
        sol0 = sol.subs({p: 0 for p in params})
        assert M * sol0 == b
        nz = [(ms[i], sol0[i]) for i in range(len(ms)) if sol0[i] != 0]
        print(f"   one exact solution uses {len(nz)} matchings; coefficients {sorted(set(c for _, c in nz))}")
        for m, c in nz[:20]: print("     ", c, m, "cost", sum(ham(a, b) for a, b in m))
    return r1 == r2

if __name__ == "__main__":
    run(2, 1)
    run(3, 1)
    run(3, 2)
