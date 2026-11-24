# "Diagonal (Pauli) family" for N = 2^k: in the Bell basis vec(P_b) the swap F is diag(s_b), s_b = +-1 (P_b^T = s_b P_b);
# C = sum_a c_a P_a (x) conj(P_a) is diag(f(b)), f(b) = sum_a c_a w(a,b), w = commutation sign.  rc(C) = #supp(c),
# rr(F - C) = #{b : f(b) != s_b}.  Exact frontier of this family by brute force (N = 4): for every support A (|A| = q) and
# every set T0 of q rows with M[T0,A] invertible, solve and count disagreements.  Integer arithmetic via exact rational solve
# at the end for the reported optimum.
import numpy as np, itertools, sys
from fractions import Fraction
k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
n = 4**k
def digs(a): return [(a >> (2*i)) & 3 for i in range(k)]   # 0=I,1=X,2=Y,3=Z
def sgn_T(a): return (-1)**sum(1 for d in digs(a) if d == 2)
def comm(a, b): return (-1)**sum(1 for x, y in zip(digs(a), digs(b)) if x != 0 and y != 0 and x != y)
s = np.array([sgn_T(b) for b in range(n)], dtype=float)
M = np.array([[comm(a, b) for a in range(n)] for b in range(n)], dtype=float)
assert np.allclose(M @ (s/ (2**k)), s)       # F = (1/N) sum_a s_a P_a (x) conj(P_a)
best = {}
for q in range(1, int(sys.argv[2]) + 1 if len(sys.argv) > 2 else 5):
    bp = n; arg = None
    rows = np.array(list(itertools.combinations(range(n), q)))
    for A in itertools.combinations(range(n), q):
        if 0 not in A and q < n: pass
        sub = M[:, A]                                   # n x q
        B = sub[rows]                                   # (#T0, q, q)
        det = np.linalg.det(B); ok = np.abs(det) > 0.5
        if not ok.any(): continue
        c = np.linalg.solve(B[ok], s[rows[ok]][..., None])[..., 0]   # (#ok, q)
        f = c @ sub.T                                   # (#ok, n)
        dis = (np.abs(f - s) > 1e-9).sum(1)
        i = dis.argmin()
        if dis[i] < bp: bp = int(dis[i]); arg = (A, c[i].copy())
    A, c = arg
    cf = [Fraction(float(x)).limit_denominator(64) for x in c]
    # exact verification of the count
    f = [sum(cf[j]*int(M[b, A[j]]) for j in range(q)) for b in range(n)]
    p = sum(1 for b in range(n) if f[b] != int(s[b]))
    assert p == bp
    best[q] = p
    print(f"N={2**k} diagonal family: q={q}: least p = {p};  support A={A} coefficients={[str(x) for x in cf]}", flush=True)
