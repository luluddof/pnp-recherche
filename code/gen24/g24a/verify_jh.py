# exact check of an UNBALANCED Johnson-regraded point by the full evaluator: 2X = T1 + T2 with T2 built from F2 = U V^t
import numpy as np, scipy.sparse as sp, itertools, sys, flint, math
from fractions import Fraction
import ev, hr, johnson
r, k, j1, j2 = map(int, sys.argv[1:5])
S = hr.setup(r, k); N, n = S['N'], S['n']; fam = S['fam']; eAs = np.array(S['eAs'])
phi1, phi2 = johnson.solve(k, j1, j2)
den = 1
for x in phi1 + phi2: den = math.lcm(den, x.denominator)
F2 = np.array([[int(phi2[len(set(a) & set(b))] * den) for b in fam] for a in fam], dtype=np.int64)
# factor F2 = U V^t over Q via rref (integer-scaled): use column basis
M = flint.fmpz_mat(F2.tolist()); r2 = M.rank()
# T2[x=A,y=i,u=j,v=B] = F2(A,B) sum_r e_B[i,r] e_A[j,r]; as C-type terms we need its rc; build tensor directly
T2 = np.einsum('ab,bir,ajr->aijb', F2, eAs, eAs)
X = np.zeros((N, n, n, N), dtype=np.int64)
for x in range(N):
    for y in range(n): X[x, y, y, x] = den
R = X - T2
rr = ev.rank(sp.csr_matrix(R.reshape(N*N, N*N)))[0]
rc = ev.rank(sp.csr_matrix(T2.transpose(0, 2, 1, 3).reshape(N*N, N*N)))[0]
print(f"N={N} k={k} (j1,j2)=({j1},{j2}): rank F2 = {r2};  exact p = rr(X - T2) = {rr}, q = rc(T2) = {rc}")
