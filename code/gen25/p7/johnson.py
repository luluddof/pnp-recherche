# Johnson-scheme regrading of Hrubes's construction:  in the Clifford-tensor framework
#   T1 = [F1(A,B) e_A e_B^t]  (rr = n rank F1),  T2 = [F2(A,B) e_B e_A^t]^{swap}  (rc = n rank F2),
#   2X = T1 + T2  <=>  F1(A,B) + (-1)^{k-|A cap B|} F2(A,B) = [A = B]  (since e_B e_A^t = (-1)^{k-|AnB|} e_A e_B^t).
# Johnson-invariant choice F_i(A,B) = phi_i(|A cap B|), phi_i polynomial of degree <= j_i  => rank <= C(t, j_i).
import sys, itertools, flint, json
from fractions import Fraction
from math import comb, log2, log
import numpy as np
def solve(k, j1, j2):
    """phi1 deg<=j1, phi2 deg<=j2 (in binomial basis C(c,m)), phi1(c) + (-1)^(k-c) phi2(c) = [c==k], c=0..k. returns (phi1, phi2) values or None"""
    nun = (j1 + 1) + (j2 + 1)
    rows = []; rhs = []
    for c in range(k + 1):
        rows.append([comb(c, m) for m in range(j1 + 1)] + [(-1) ** (k - c) * comb(c, m) for m in range(j2 + 1)])
        rhs.append(1 if c == k else 0)
    A = flint.fmpq_mat(k + 1, nun, [x for r in rows for x in r]); b = flint.fmpq_mat(k + 1, 1, rhs)
    Ab = flint.fmpq_mat(k + 1, nun + 1, [x for i, r in enumerate(rows) for x in r + [rhs[i]]])
    if A.rank() != Ab.rank(): return None
    # least-norm-ish particular solution: use rref
    R = Ab.rref()[0]
    sol = [Fraction(0)] * nun
    piv = []
    for i in range(k + 1):
        for j in range(nun + 1):
            if R[i, j] != 0: piv.append((i, j)); break
    for i, j in piv:
        if j == nun: return None
        sol[j] = Fraction(int(R[i, nun].p), int(R[i, nun].q))
    phi1 = [sum(sol[m] * comb(c, m) for m in range(j1 + 1)) for c in range(k + 1)]
    phi2 = [sum(sol[j1 + 1 + m] * comb(c, m) for m in range(j2 + 1)) for c in range(k + 1)]
    for c in range(k + 1): assert phi1[c] + (-1) ** (k - c) * phi2[c] == (c == k)
    return phi1, phi2
def frank(phi, fam):
    N = len(fam); S = [set(A) for A in fam]
    den = 1
    for x in phi: den = den * x.denominator // __import__('math').gcd(den, x.denominator)
    M = flint.fmpz_mat(N, N, [int(phi[len(S[a] & S[b])] * den) for a in range(N) for b in range(N)])
    return M.rank()
if __name__ == "__main__":
    t, k, N = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); n = N
    fam = list(itertools.combinations(range(t), k))[:N]
    print(f"t={t} k={k} N=n={N}:  j1 j2  rank F1|A  rank F2|A  (C(t,j1),C(t,j2))   p = n r1, q = n r2   log_N p, log_N q")
    out = []
    for j1 in range(-1, k + 1):
        for j2 in range(-1, k + 1):
            if j1 + j2 != k - 1 and not (j1 + j2 == k - 2): continue
            s = solve(k, j1, j2) if j1 >= 0 and j2 >= 0 else None
            if j1 < 0: s = ([Fraction(0)] * (k + 1), [Fraction((-1) ** (k - c) * (c == k)) for c in range(k + 1)])
            if j2 < 0: s = ([Fraction(int(c == k)) for c in range(k + 1)], [Fraction(0)] * (k + 1))
            if s is None: print(f"   {j1:2} {j2:2}  no solution"); continue
            r1, r2 = frank(s[0], fam), frank(s[1], fam)
            p, q = n * r1, n * r2
            out.append((j1, j2, r1, r2, p, q))
            print(f"   {j1:2} {j2:2}  {r1:6} {r2:6}   ({comb(t,j1) if j1>=0 else 0},{comb(t,j2) if j2>=0 else 0})   ({p}, {q})   {log(max(p,1))/log(N):.3f} {log(max(q,1))/log(N):.3f}")
    json.dump(out, open(f"johnson_t{t}_k{k}_N{N}.json", "w"))
