# Exact verification of a decomposition X = R + C given integer (or rational) column factors G_l, G'_l:
#   C(x,y,u,v) = sum_l G_l[x,u] G'_l[y,v];  R := X - C;  exact ranks (python-flint, over Q);  explicit rank factorisation
#   R = sum_j f_j(x,y) f'_j(u,v) (rref over Q) re-expanded and compared entry by entry with X - C.
import numpy as np, sys, flint, json, itertools
from fractions import Fraction
def verify(N, G, Gp, name=""):
    q = len(G)
    X = [[[[Fraction(int(x == v and y == u)) for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
    C = [[[[sum(Fraction(G[l][x][u]) * Fraction(Gp[l][y][v]) for l in range(q)) for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
    Rrow = [[X[x][y][u][v] - C[x][y][u][v] for u in range(N) for v in range(N)] for x in range(N) for y in range(N)]
    Ccol = [[C[x][y][u][v] for y in range(N) for v in range(N)] for x in range(N) for u in range(N)]
    den = 1
    for row in Rrow + Ccol:
        for e in row: den = np.lcm(den, e.denominator)
    den = int(den)
    Rm = flint.fmpq_mat(N*N, N*N, [flint.fmpq(e.numerator, e.denominator) for row in Rrow for e in row])
    Cm = flint.fmpq_mat(N*N, N*N, [flint.fmpq(e.numerator, e.denominator) for row in Ccol for e in row])
    p, qq = Rm.rank(), Cm.rank()
    # rank factorisation of R: R = B * E where E = nonzero rows of rref, B = pivot columns of R
    E, rk = Rm.rref(); piv = []
    for i in range(rk):
        j = next(j for j in range(N*N) if E[i, j] != 0); piv.append(j)
    ok = True
    for a in range(N*N):
        for b in range(N*N):
            val = sum((Rm[a, piv[i]] * E[i, b] for i in range(rk)), flint.fmpq(0))
            xx, yy, uu, vv = a // N, a % N, b // N, b % N
            tot = val + flint.fmpq(C[xx][yy][uu][vv].numerator, C[xx][yy][uu][vv].denominator)
            if tot != (1 if (xx == vv and yy == uu) else 0): ok = False
    print(f"{name}: N={N}: rr(R) = {p}, rc(C) = {qq} (number of column terms {q}); explicit rank-{rk} factorisation of R re-expanded, R + C == X entry by entry: {ok}")
    assert ok
    return p, qq
if __name__ == "__main__":
    z = np.load(sys.argv[2]); N, q = int(sys.argv[1]), int(sys.argv[3])
    zr = np.rint(z).astype(int); print("max rounding error:", np.abs(z - zr).max())
    G = zr[:q*N*N].reshape(q,N,N).tolist(); Gp = zr[q*N*N:].reshape(q,N,N).tolist()
    verify(N, G, Gp, sys.argv[2]); print("G =", G); print("G' =", Gp)
