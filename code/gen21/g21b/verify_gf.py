# Independent exact verification of a crossing decomposition over GF(P):  X = R + C, R = X - C,
# rr(R) = rank of R as (x,y)|(u,v), rc(C) = rank of C as (x,u)|(y,v), computed by own Gaussian elimination mod P.
# Also checks R == U V entry by entry when U, V are present.
import sys, json, itertools
def rank_mod(M, P):
    M = [[e % P for e in row] for row in M]; r = 0; rows = len(M); cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = next((i for i in range(r, rows) if M[i][c]), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; inv = pow(M[r][c], P-2, P)
        M[r] = [(e*inv) % P for e in M[r]]
        for i in range(rows):
            if i != r and M[i][c]:
                f = M[i][c]; M[i] = [(a - f*b) % P for a, b in zip(M[i], M[r])]
        r += 1
    return r
def verify(dec):
    N, P = dec["N"], dec["field"]; G, H = dec["G"], dec["H"]; q = len(G)
    C = {}
    for x, y, u, v in itertools.product(range(N), repeat=4):
        C[(x, y, u, v)] = sum(G[l][x][u]*H[l][y][v] for l in range(q)) % P
    R = {k: ((k[0] == k[3] and k[1] == k[2]) - C[k]) % P for k in C}
    if "U" in dec:
        U, V = dec["U"], dec["V"]; p = len(V)
        for x, y, u, v in itertools.product(range(N), repeat=4):
            assert sum(U[x*N+y][k]*V[k][u*N+v] for k in range(p)) % P == R[(x, y, u, v)], "R != UV"
    for k in C:  # re-expansion R + C == X
        assert (R[k] + C[k]) % P == (k[0] == k[3] and k[1] == k[2])
    rr = rank_mod([[R[(x, y, u, v)] for u in range(N) for v in range(N)] for x in range(N) for y in range(N)], P)
    rc = rank_mod([[C[(x, y, u, v)] for y in range(N) for v in range(N)] for x in range(N) for u in range(N)], P)
    return rr, rc
if __name__ == "__main__":
    for f in sys.argv[1:]:
        d = json.load(open(f)); print(f, "field GF(%d) N=%d: rr(R)=%d rc(C)=%d  (re-expansion OK)" % ((d["field"], d["N"]) + verify(d)))
