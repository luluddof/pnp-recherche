# exact check: given rational G[l][x][u], H[l][y][v] (Fractions/ints), return (rank_Q(X - C) in (xy|uv), rank_Q(C) in (xu|yv))
import flint, numpy as np
from fractions import Fraction
def ranks(N, G, H):
    q = len(G)
    def fq(a): a = Fraction(a); return flint.fmpq(a.numerator, a.denominator)
    R = [[0] * (N * N) for _ in range(N * N)]; Cr = [[0] * (N * N) for _ in range(N * N)]
    for x in range(N):
        for y in range(N):
            for u in range(N):
                for v in range(N):
                    c = sum(Fraction(G[l][x][u]) * Fraction(H[l][y][v]) for l in range(q))
                    R[x * N + y][u * N + v] = Fraction(int(x == v and y == u)) - c
                    Cr[x * N + u][y * N + v] = c
    Rm = flint.fmpq_mat(N * N, N * N, [fq(e) for row in R for e in row])
    Cm = flint.fmpq_mat(N * N, N * N, [fq(e) for row in Cr for e in row])
    return Rm.rank(), Cm.rank()
