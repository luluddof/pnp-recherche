# Exact checks (integer arithmetic, python-flint ranks over Q):
#  (1) base identities R + C = d X for N = 2 (one-bit), 2,4,8 (Hurwitz), general N (I - F);
#  (2) tensor powers: A_i = coeff of T^i in prod_blocks (R + T C); thresholds R_th = sum_{i<=i0} A_i,
#      C_th = sum_{i>i0} A_i give exactly feasible pairs (p,q) = (rr(R_th), rc(C_th));
#  (3) every feasible pair is tested against the hyperbola bound p >= k(N-kq), q >= k(N-kp);
#  (4) telescoping pieces A_i: max_i rr*rc and the implied rho.
import numpy as np, itertools, math, sys, flint
def rank(M): return flint.fmpz_mat(M.astype(object).tolist()).rank() if M.size else 0
def rowm(T): N = T.shape[0]; return T.reshape(N*N, N*N)                       # (x,y)|(u,v)
def colm(T): N = T.shape[0]; return T.transpose(0,2,1,3).reshape(N*N, N*N)    # (x,u)|(y,v)
def X(N):
    T = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): T[x,y,y,x] = 1
    return T
def onebit():      # X = [x=u][y=v] - [x!=y](-1)^x [u!=v](-1)^u   (paper, App. H); R = second term, C = first
    N = 2; R = np.zeros((N,)*4, dtype=np.int64); C = np.zeros((N,)*4, dtype=np.int64)
    for x,y,u,v in itertools.product(range(N), repeat=4):
        C[x,y,u,v] = (x==u)*(y==v); R[x,y,u,v] = -(x!=y)*(-1)**x*(u!=v)*(-1)**u
    return R, C, 1
def ImF(N):        # X = C + R, C = [x=u][y=v] (rc 1), R = X - C (rr N(N-1)/2)
    C = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): C[x,y,x,y] = 1
    return X(N) - C, C, 1
def hurwitz(n):    # N = 2^n, n <= 3: sigma from eq. (27) of the paper restricted to the first n bits
    N = 2**n
    def bits(p): b = [(p >> (n-1-i)) & 1 for i in range(n)]; return b + [0]*(3-n)
    def sig(p, q):
        P, Q = bits(p), bits(q)
        phi = sum(P[i]*Q[j] for i in range(3) for j in range(i,3)) + Q[0]*P[1]*P[2] + P[0]*Q[1]*P[2] + P[0]*P[1]*Q[2]
        return (-1)**phi
    R = np.zeros((N,)*4, dtype=np.int64); C = np.zeros((N,)*4, dtype=np.int64)
    for x,y,u,v in itertools.product(range(N), repeat=4):
        if x ^ y == v ^ u: R[x,y,u,v] = sig(x,y)*sig(v,u)
        if x ^ u == v ^ y: C[x,y,u,v] = sig(x,u)*sig(v,y)
    return R, C, 2
def hyper_ok(N, p, q):
    return p >= max(k*(N-k*q) for k in range(N+1)) and q >= max(k*(N-k*p) for k in range(N+1))
def power(name, base, t):
    R, C, d = base; N0 = R.shape[0]; N = N0**t
    assert np.array_equal(R + C, d*X(N0)), "base identity fails"
    out = []
    for arr, tag in ((rowm, 'row'), (colm, 'col')):
        Rm, Cm = arr(R), arr(C)
        polys = {0: np.ones((1,1), dtype=np.int64)}
        for b in range(t):
            new = {}
            for i, M in polys.items():
                for add, B in ((0, Rm), (1, Cm)):
                    K = np.kron(M, B); new[i+add] = new.get(i+add, 0) + K
            polys = new
        out.append(polys)
    prow, pcol = out
    tot = sum(prow.values()); Xr = rowm(X(N0)); Xt = np.ones((1,1), dtype=np.int64)
    for b in range(t): Xt = np.kron(Xt, Xr)
    assert np.array_equal(tot, d**t * Xt), "power identity fails"
    print(f"{name}^(x){t}: N={N}, identity sum_i A_i = {d**t} X verified exactly")
    worst = 0
    for i in range(t+1):
        rr, rc = rank(prow[i]), rank(pcol[i]); worst = max(worst, rr*rc)
        print(f"   piece i={i}: rr={rr} rc={rc} rr*rc={rr*rc}")
    if N > 1: print(f"   max rr*rc = {worst} vs N^2 = {N*N}, N^4 = {N**4};  rho = log(max)/(4 log N) = {math.log(worst)/(4*math.log(N)):.4f}  (t+1 = {t+1} pieces, not charged)")
    for i0 in range(-1, t+1):
        Rth = sum((prow[i] for i in range(0, i0+1)), np.zeros_like(tot)); Cth = sum((pcol[i] for i in range(i0+1, t+1)), np.zeros_like(tot))
        p, q = rank(Rth), rank(Cth)
        print(f"   threshold i0={i0}: feasible (p,q)=({p},{q})  p*q={p*q}  p*q/N^2={p*q/N**2:.3f}  p+q={p+q}  hyperbola bound respected: {hyper_ok(N,p,q)}")
        assert hyper_ok(N, p, q)
if __name__ == "__main__":
    power("onebit", onebit(), 1); power("onebit", onebit(), 2); power("onebit", onebit(), 3); power("onebit", onebit(), 4)
    for N in (3,4,5): power(f"I-F(N={N})", ImF(N), 1)
    power("I-F(N=3)", ImF(3), 2)
    for n in (1,2,3): power(f"hurwitz{2**n}", hurwitz(n), 1)
    power("hurwitz4", hurwitz(2), 2)
