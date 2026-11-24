# T7: strengthening found during review:  osr(B) < m*n  ==>  rank(A) >= Q   (Theorem C'), and tightness at b = mn-1.
import sys, numpy as np
sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from ff import rank1
rng = np.random.default_rng(11)
def buildP(m, n, Q):
    P = np.zeros((n, m, Q, m, n, Q), dtype=np.int64)
    for e in range(m):
        for g in range(n):
            for w in range(Q): P[g, e, w, e, g, w] = 1
    return P
def osr(T, m, n, Q, p): return rank1(T.transpose(0, 3, 1, 2, 4, 5).reshape(n*m, -1), p)
def rk(T, p): return rank1(T.reshape(T.shape[0]*T.shape[1]*T.shape[2], -1), p)
cnt = 0
for (m, n, Q) in [(2,2,2),(2,2,4),(3,2,2),(2,3,3),(3,3,2),(2,2,8),(4,4,4)]:
    P = buildP(m, n, Q); Dm = m*n*Q
    for p in (2, 3, 5):
        # tight example: A = E_{g0 e0} (x) P_{g0 e0}: rank Q, osr(P-A) = mn-1
        A = np.zeros_like(P); A[0, :, :, 0, :, :] = P[0, :, :, 0, :, :]
        assert rk(A, p) == Q and osr((P - A) % p, m, n, Q, p) == m*n - 1
        # more generally c of rank rho: A = L (x) (c^T (x) 1_W) with c(L)=1: rank(A) = rank(L) rho Q
        # every A with rank < Q must leave osr(P-A) = mn:
        for r in range(0, Q):
            for t in range(40 if Dm <= 48 else 6):
                mode = t % 4
                if mode == 0: A = (rng.integers(0, p, (Dm, r)) @ rng.integers(0, p, (r, Dm))) % p
                elif mode == 1:   # rank-r piece of the tight block
                    blk = P[0, :, :, 0, :, :].reshape(m*Q, n*Q); sel = np.zeros((n*Q, n*Q), dtype=np.int64)
                    idx = rng.choice(n*Q, r, replace=False); sel[idx, idx] = 1
                    A = np.zeros_like(P); A[0, :, :, 0, :, :] = (blk @ sel).reshape(m, Q, n, Q); A = A.reshape(Dm, Dm)
                elif mode == 2:   # P on r basis inputs
                    A = np.zeros((Dm, Dm), dtype=np.int64); c = rng.choice(Dm, r, replace=False); A[:, c] = P.reshape(Dm, Dm)[:, c]
                else:             # product-structured low rank: x = x0 (x) x1, y = y0 (x) y1
                    A = np.zeros((Dm, Dm), dtype=np.int64)
                    for _ in range(r):
                        x = np.kron(rng.integers(0, p, n), rng.integers(0, p, m*Q)); y = np.kron(rng.integers(0, p, m), rng.integers(0, p, n*Q))
                        A += np.outer(x, y)
                A = (A % p).reshape(P.shape)
                assert rk(A, p) < Q and osr((P - A) % p, m, n, Q, p) == m*n, (m, n, Q, p, r, mode)
                cnt += 1
print(f"T7 ok: tight block example (rank Q, osr mn-1) verified; {cnt} matrices A of rank < Q all leave osr(P-A) = mn (full)")
