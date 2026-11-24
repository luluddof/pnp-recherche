# "Laser" test: S_{2^k} = sum_T R^T (x) I^{T^c} (R = S_2 - I, pieces (rank,osr) = (4^{k-|T|}, 4^{|T|}) ).
# Restrict to a sub-swap on a code X' in F_2^k and threshold by |T| (and by weighted |T|); exact ranks over Q (flint).
import itertools, numpy as np, flint, sys
def rank(M): return flint.fmpz_mat(M.astype(int).tolist()).rank() if M.size else 0
def code_from_gen(G):
    k = G.shape[1]; G = np.array(G)
    return sorted({tuple(np.array(c) @ G % 2) for c in itertools.product([0,1], repeat=G.shape[0])})
def pieces_by_j(X, k):
    X = np.array(X); n = len(X)
    # For pair (x,y)->(u,v): piece T nonzero iff for i not in T: u_i=x_i, v_i=y_i ; for i in T: x_i!=y_i and value is
    # R(x_i,y_i,u_i,v_i) = [x=v][y=u]-[x=u][y=v] (per coordinate). Sum over |T|=j.
    P = np.zeros((k+1, n, n, n, n), dtype=np.int64)
    for a in range(n):
        for b in range(n):
            x, y = X[a], X[b]
            for c in range(n):
                for d in range(n):
                    u, v = X[c], X[d]
                    # per coordinate weights: I-factor value [u=x][v=y]; R-factor value
                    Iv = ((u == x) & (v == y)).astype(int)
                    Rv = ((x == v) & (y == u)).astype(int) - Iv
                    # generating polynomial prod_i (Iv_i + t Rv_i): coefficient of t^j
                    poly = np.array([1], dtype=np.int64)
                    for i in range(k):
                        poly = np.convolve(poly, np.array([Iv[i], Rv[i]]))
                        if not poly.any(): break
                    if poly.any(): P[:len(poly), a, b, c, d] += poly
    return P
def osr(M4):  # M4[a,b,c,d] rows (a,b) cols (c,d); OSR = rank of realigned (a,c)x(b,d)
    n = M4.shape[0]; return rank(M4.transpose(0, 2, 1, 3).reshape(n*n, n*n))
def rk(M4):
    n = M4.shape[0]; return rank(M4.reshape(n*n, n*n))
def run(name, X, k):
    n = len(X); P = pieces_by_j(X, k)
    S = P.sum(0); n2 = n*n
    chk = rk(S), osr(S)
    best = []
    for j0 in range(k+2):
        Rm = P[j0:].sum(0); Cm = P[:j0].sum(0)
        best.append((max(rk(Rm), osr(Cm)), rk(Rm), osr(Cm), j0))
    b = min(best)
    print(name, 'N\'=', n, 'check swap rank/osr', chk, 'per-j0 (p,q):', [(t[1], t[2]) for t in best], 'balanced', b[0],
          'exp %.3f' % (np.log(b[0]) / np.log(n)), flush=True)
k = int(sys.argv[1]) if len(sys.argv) > 1 else 6
I = np.eye
codes = {
  'full3': (np.eye(3, dtype=int), 3),
  'even4': (np.array([[1,1,0,0],[0,1,1,0],[0,0,1,1]]), 4),
  'rep2x3': (np.array([[1,1,0,0,0,0],[0,0,1,1,0,0],[0,0,0,0,1,1]]), 6),
  'ham743': (np.array([[1,0,0,0,0,1,1],[0,1,0,0,1,0,1],[0,0,1,0,1,1,0],[0,0,0,1,1,1,1]]), 7),
  'rm13': (np.array([[1,1,1,1,1,1,1,1],[0,0,0,0,1,1,1,1],[0,0,1,1,0,0,1,1],[0,1,0,1,0,1,0,1]]), 8),
  'simplex3': (np.array([[1,0,0,1,1,0,1],[0,1,0,1,0,1,1],[0,0,1,0,1,1,1]]), 7),
  'even5_4': (np.array([[1,1,0,0,0],[0,1,1,0,0],[0,0,1,1,0],[0,0,0,1,1]]), 5),
}
for nm, (G, kk) in codes.items():
    run(nm, code_from_gen(G), kk)
