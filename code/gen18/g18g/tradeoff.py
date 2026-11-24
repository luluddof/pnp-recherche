# Numerical trade-off (p,q): X = sum_{j<=p} f_j(x,y) f'_j(u,v) + sum_{k<=q} g_k(x,u) g'_k(y,v)
# ALS with restarts; reports the least p with residual < tol for each q.  [N] numerical only.
import numpy as np, sys, time
def crossing(N):
    T = np.zeros((N,N,N,N))
    for x in range(N):
        for y in range(N): T[x,y,y,x] = 1.0   # [x=v][y=u], index order (x,y,u,v)
    return T
def als(N, p, q, iters, rng, cplx=False, tol=1e-11):
    T = crossing(N); n2 = N*N
    Mr = T.reshape(n2, n2)                              # (x,y)|(u,v)
    Mc = T.transpose(0,2,1,3).reshape(n2, n2)           # (x,u)|(y,v)
    def rnd(a): 
        z = rng.standard_normal((n2, a))
        return z + 1j*rng.standard_normal((n2,a)) if cplx else z
    F, Fp, G, Gp = rnd(p), rnd(p), rnd(q), rnd(q)
    def to_r(Mc_like): return Mc_like.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n2,n2)
    res = None
    for it in range(iters):
        Cc = G @ Gp.T if q else np.zeros((n2,n2))
        Er = Mr - to_r(Cc)                              # target for row part
        if p:
            F = np.linalg.lstsq(Fp, Er.T, rcond=None)[0].T
            Fp = np.linalg.lstsq(F, Er, rcond=None)[0].T
        Rr = F @ Fp.T if p else np.zeros((n2,n2))
        Ec = Mc - to_r(Rr)                              # to_r is an involution (swap y<->u)
        if q:
            G = np.linalg.lstsq(Gp, Ec.T, rcond=None)[0].T
            Gp = np.linalg.lstsq(G, Ec, rcond=None)[0].T
        if it % 25 == 0 or it == iters-1:
            Cc = G @ Gp.T if q else np.zeros((n2,n2))
            res = np.linalg.norm(Mr - (F @ Fp.T if p else 0) - to_r(Cc))
            if res < tol: break
    return res
def lower(N, q):
    return max(k*(N-k*q) for k in range(0, N+1))
if __name__ == "__main__":
    N = int(sys.argv[1]); cplx = len(sys.argv) > 2 and sys.argv[2] == "c"
    restarts = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    iters = int(sys.argv[4]) if len(sys.argv) > 4 else 3000
    rng = np.random.default_rng(1)
    t0 = time.time()
    print(f"N={N} field={'C' if cplx else 'R'}  q : hyperbola bound max_k k(N-kq) : least p found (best residual at p-1)")
    pprev = N*N
    for q in range(0, N*N+1):
        lb = lower(N, q); found = None; below = None
        p = pprev
        while p >= max(lb,0):
            best = min(als(N, p, q, iters, rng, cplx) for _ in range(restarts))
            if best < 1e-9: found = p; p -= 1
            else: below = best; break
        pprev = found if found is not None else pprev
        print(f"  q={q:2d}  bound={lb:2d}  p_found={found}  p+q={None if found is None else found+q}  resid(p-1)={below}  t={time.time()-t0:.0f}s", flush=True)
        if found == 0: break
        if time.time()-t0 > 420: print("  time limit"); break
