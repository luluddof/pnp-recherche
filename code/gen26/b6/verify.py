# Exact (over Q) construction and verification of a relaxation point: X = R + C, rr(R) <= dim U_W, rc(C) <= dim U'_W'.
import sys, math
sys.dont_write_bytecode = True
import flint
from relax import setup, g_rows, gp_rows
from appb_exact import rr, rc
h, W, Wp = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); Nreq = int(sys.argv[4]) if len(sys.argv) > 4 else None
N, mus, bs = setup(h, Nreq); n = N * N
G = g_rows(h, N, mus, bs, max(W, 0)); Gp = gp_rows(h, N, mus, bs, max(Wp, 0))
Uv = [G[s] for s in G if sum(s) <= W]; Upv = [Gp[s] for s in Gp if sum(s) <= Wp]
# basis of U' over Q
Mp = flint.fmpz_mat(len(Upv), n, [x for v in Upv for x in v])
Bp_r, rk = flint.fmpq_mat(Mp).rref(); Bp = [[Bp_r[i, j] for j in range(n)] for i in range(rk)]
# annihilator of U over Q
if Uv:
    MU = flint.fmpz_mat(len(Uv), n, [x for v in Uv for x in v]); Xn, nul = MU.nullspace()
    K = [[Xn[i, j] for i in range(n)] for j in range(nul)]
else: K = [[int(i == j) for j in range(n)] for i in range(n)]
dU = n - len(K)
# unknowns F_m(x,u): for each x solve M f = r_x
rows = len(K) * N; cols = len(Bp) * N
M = flint.fmpq_mat(rows, cols + N)
for j, k in enumerate(K):
    for y in range(N):
        for m, g in enumerate(Bp):
            for u in range(N):
                M[j*N+y, m*N+u] = sum(g[y*N+v] * k[u*N+v] for v in range(N))
        for x in range(N): M[j*N+y, cols+x] = k[y*N+x]
Rr, r = M.rref()
# read off a particular solution: pivot columns
piv = []; i = 0
for c in range(cols + N):
    if i < r and Rr[i, c] != 0: piv.append(c); i += 1
assert all(c < cols for c in piv), "infeasible"
F = [[0]*cols for _ in range(N)]   # F[x][(m,u)]
for i, c in enumerate(piv):
    for x in range(N): F[x][c] = Rr[i, cols + x]
C = [[[[flint.fmpq(0)]*N for _ in range(N)] for _ in range(N)] for _ in range(N)]
for x in range(N):
    for m, g in enumerate(Bp):
        for u in range(N):
            f = F[x][m*N+u]
            if f == 0: continue
            for y in range(N):
                for v in range(N): C[x][y][u][v] = flint.fmpq(C[x][y][u][v]) + f * g[y*N+v]
den = 1
for x in range(N):
    for y in range(N):
        for u in range(N):
            for v in range(N):
                den = den * C[x][y][u][v].q // math.gcd(den, C[x][y][u][v].q)
Ci = [[[[int(C[x][y][u][v]*den) for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
Ri = [[[[den*(x == v)*(y == u) - Ci[x][y][u][v] for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
print(f"h={h} N={N} W={W} W'={Wp}: dimU={dU} dimU'={len(Bp)};  exact rr(R)={rr(Ri, N)}  rc(C)={rc(Ci, N)}  (den={den})")
