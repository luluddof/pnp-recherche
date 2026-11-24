# T4: the two mechanism identities of Theorem C with m != n and a spectator W, exact mod p; T5: random decompositions, general (m,n,Q);
# T6: integer optimisation of the inequality system (constant 2/9).
import sys, numpy as np
sys.path.insert(0, sys.argv[0].rsplit('/',1)[0])
from ff import rank1, nullspace
rng = np.random.default_rng(7)
def buildP(m, n, Q):
    # P[(g,e,w),(e',g',w')] : E(x)G(x)W -> G(x)E(x)W
    P = np.zeros((n, m, Q, m, n, Q), dtype=np.int64)
    for e in range(m):
        for g in range(n):
            for w in range(Q): P[g, e, w, e, g, w] = 1
    return P
def osr(T, m, n, Q, p):   # T[g,e,w,e',g',w'] -> rows (g,e') cols (e,w,g',w')
    return rank1(T.transpose(0, 3, 1, 2, 4, 5).reshape(n*m, -1), p)
def rk(T, p): return rank1(T.reshape(T.shape[0]*T.shape[1]*T.shape[2], -1), p)
cnt = 0
for (m, n, Q) in [(2,2,1),(3,3,1),(2,3,1),(3,2,2),(2,2,2),(4,3,1),(3,4,2),(4,4,1)]:
    P = buildP(m, n, Q)
    for p in (2, 3, 7):
        for t in range(15):
            # (a): pi: E->K^k surjective, S random subspace of G(x)W of dim d: compressed P has rank k*d
            k = int(rng.integers(1, m+1)); d = int(rng.integers(0, n*Q+1))
            while True:
                pi = rng.integers(0, p, (k, m))
                if rank1(pi, p) == k: break
            while True:
                S = rng.integers(0, p, (n*Q, d)) if d else np.zeros((n*Q, 0), dtype=np.int64)
                if d == 0 or rank1(S, p) == d: break
            Y = np.einsum('ke,gewfhv->gkwfhv', pi, P).reshape(n, k, Q, m, n*Q) @ S   # (g,k,w,e',d)
            assert rank1(Y.reshape(n*k*Q, m*d) % p, p) == k*d if d else True
            # (b): Phi in E of dim k, Psi in G(x)W of dim d: osr(P|Phi(x)Psi) = k * dim(G-support of Psi) >= k*d/Q
            Phi = pi.T                                                 # m x k
            R = np.einsum('gewfhv,fk->gewkhv', P, Phi).reshape(n, m, Q, k, n*Q) @ S     # (g,e,w,k,d)
            o = rank1(R.transpose(0, 3, 1, 2, 4).reshape(n*k, m*Q*max(d,1))[:, :m*Q*d] % p, p) if d else 0
            supp = rank1(S.reshape(n, Q*d) % p, p) if d else 0     # rows g: dimension of the G-support
            assert o == k*supp and o*Q >= k*d, (m, n, Q, p, k, d, o, supp)
            cnt += 1
print(f"T4 ok: {cnt} random (pi,S) / (Phi,Psi): rank Y(P) = k dim S ; osr(P|Phi(x)Psi) = k dim supp_G(Psi) >= k dim Psi / Q")
def check(m, n, Q, r, b):
    bad = []
    for k in range(m+1):
        if k*Q*n > r + k*k*Q*b: bad.append(('a', k))
        if k*n*Q > b*Q + k*k*r: bad.append(('b', k))
    for k in range(n+1):
        if k*Q*m > r + k*k*Q*b: bad.append(("a'", k))
        if k*m*Q > b*Q + k*k*r: bad.append(("b'", k))
    return bad
cnt = 0; tight = {}
for (m, n, Q) in [(2,2,1),(2,3,1),(3,2,2),(2,2,2),(4,3,1),(3,3,2),(2,2,4)]:
    P = buildP(m, n, Q); Dm = m*n*Q
    for p in (2, 3, 5):
        for b in range(0, m*n+1):
            for t in range(12):
                B = np.zeros_like(P)
                for _ in range(b):
                    mode = rng.integers(0, 3)
                    if mode == 0: L = rng.integers(0, p, (n, m)); Mk = rng.integers(0, p, (m, Q, n, Q))
                    elif mode == 1:
                        L = np.zeros((n, m), dtype=np.int64); g, e = rng.integers(0, n), rng.integers(0, m); L[g, e] = 1
                        Mk = P[g, :, :, e, :, :].copy()          # the matching realignment row of P
                    else:
                        L = (rng.random((n, m)) < .5).astype(np.int64); Mk = np.einsum('gewfhv,gf->ewhv', P, L)
                    B += np.einsum('gf,ewhv->gewfhv', L, Mk)
                B %= p
                bb = osr(B, m, n, Q, p); r = rk((P - B) % p, p); v = check(m, n, Q, r, bb); cnt += 1
                assert not v, (m, n, Q, p, bb, r, v)
                tight[(m,n,Q,bb)] = min(tight.get((m,n,Q,bb), 10**9), r)
        for r in range(0, Dm+1):
            for t in range(6):
                A = (rng.integers(0, p, (Dm, r)) @ rng.integers(0, p, (r, Dm))) % p if r else np.zeros((Dm, Dm), dtype=np.int64)
                if t % 2: A = (P.reshape(Dm, Dm) @ np.diag((np.arange(Dm) < r).astype(np.int64))[rng.permutation(Dm)][:, rng.permutation(Dm)]) % p
                A = A.reshape(P.shape); rr = rk(A, p); bb = osr((P - A) % p, m, n, Q, p); v = check(m, n, Q, rr, bb); cnt += 1
                assert not v, (m, n, Q, p, rr, bb, v)
print(f"T5 ok: {cnt} exact decompositions with general (m,n,Q), no violation of (a),(b),(a'),(b')")
# T6: the inequality system of Theorem 2 alone: min of p*q over integer points with p,q >= 1
worst = 1
for N in range(2, 241):
    best = None
    for q in range(1, N*N+1):
        pmin = max(max(k*(N - k*q) for k in range(N+1)), 1)
        # smallest p >= pmin with q >= k(N-kp) for all k
        pp = pmin
        while any(q < k*(N - k*pp) for k in range(N+1)): pp += 1
        if best is None or pp*q < best[0]: best = (pp*q, pp, q)
        if q > N: break
    worst = min(worst, best[0]/N**2)
    if N in (2,3,4,5,6,9,12,30,240): print(f"   N={N}: min p*q allowed by Theorem 2 = {best[0]} at (p,q)=({best[1]},{best[2]}); /N^2 = {best[0]/N**2:.4f}")
print(f"T6: inf over N<=240 of min(pq)/N^2 = {worst:.5f}  (2/9 = {2/9:.5f}, 3/16 = {3/16})")
