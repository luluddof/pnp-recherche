# Exact Jordan data at eigenvalues 0 and 1 of M_row = R_row S^-1 and M_col = C_col S^-1 (ranks of powers mod a prime).
# p_s = N^2 - algmult_0(M_row) (rank of the semisimple part), q_s likewise; LB h >= log_{N^2}(p_s^mu (N^2-q_s)^(1-mu) + (N^2-p_s)^mu q_s^(1-mu))
# by lower semicontinuity of rank (the semisimple part lies in the closure of the conjugation orbit) + the 'tops' argument.
import sys, glob, numpy as np, flint, math
from pieces import *
PR = 1000003
def nm(A): return flint.nmod_mat(A.shape[0], A.shape[1], [int(v) % PR for v in A.ravel()], PR)
def weyr(M, n):
    rs = [n]; P = flint.nmod_mat(n, n, [1 if i == j else 0 for i in range(n) for j in range(n)], PR)
    while True:
        P = P*M; r = P.rank(); rs.append(r)
        if r == rs[-2]: break
    return rs   # rank of M^j, j=0..; alg mult of 0 = n - final rank
mus = np.linspace(0, 1, 401); hB = np.hypot(mus, 1-mus); inn = slice(1, -1); W = mus[inn]*(1-mus[inn])
def doit(name, N, R, C, s):
    n = N*N; S = swap(N)
    Mr = rowflat(R) @ rowflat(S).T; Mc = colflat(C) @ colflat(S).T     # = s * M (integer)
    I = np.eye(n, dtype=np.int64)
    w0r = weyr(nm(Mr), n); w1r = weyr(nm(Mr - s*I), n); w0c = weyr(nm(Mc), n); w1c = weyr(nm(Mc - s*I), n)
    p, q = w0r[1], w0c[1]; ps, qs = w0r[-1], w0c[-1]
    L = math.log(n)
    hLB = np.log(ps**mus*(n-qs)**(1-mus) + (n-ps)**mus*qs**(1-mus))/L
    print(f"{name:20s} N={N} p={p} q={q} | rank M_row^j {w0r} (M_row-1)^j {w1r[:4]}.. | rank M_col^j {w0c} (M_col-1)^j {w1c[:4]}.."
          f" | p_s={ps} q_s={qs} | h_LB(1/2)={hLB[200]:.4f} min (h_LB-circle)/(mu(1-mu)) = {((hLB-hB)[inn]/W).min():+.4f}", flush=True)
for f in sys.argv[1:]:
    if f == 'oct': N, (R, C, s) = octonion()
    elif f.startswith('F1_'): N, (R, C, s) = F1(int(f[3:]))
    else: N, (R, C, s) = load_json(f)
    doit(f.split('/')[-1], N, R, C, s)
