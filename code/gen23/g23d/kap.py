# [N] kappa(P) = sup_{G,G'} |tr(F (G(x)G') (1-P))|^2 / ||(G(x)G')(1-P)||^2, P = rank-p projection on V(x)V.
# Inner sup over G closed form (generalized Rayleigh); outer over G' (and optionally P) by BFGS with restarts.
import numpy as np, sys
from scipy.optimize import minimize

def swapmat(N):
    F = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N): F[x*N+y, y*N+x] = 1
    return F

def inner(Gp, Pc, F, N):
    # Pc = 1 - P (N^2 x N^2).  basis E_ab (x) G'
    Bs = np.einsum('ax,by,uv->abxuyv', np.eye(N), np.eye(N), Gp)  # E_ab[x,u]=d(a,x)d(b,u)
    Bs = Bs.reshape(N*N, N*N, N*N)  # index ab, then operator rows (x,u)?? fix below
    return Bs

def ops(Gp, N):
    # O[k] = E_ab (x) G' as N^2 x N^2 operator with rows (x,y), cols (u,v): O[x,y,u,v] = E_ab[x,u] G'[y,v]
    O = np.zeros((N*N, N, N, N, N), complex)
    for a in range(N):
        for b in range(N):
            O[a*N+b, a, :, b, :] = Gp
    return O.reshape(N*N, N*N, N*N)

def kap_fixed(Gp, Pc, F, N):
    O = ops(Gp, N)
    c = np.einsum('ij,kji->k', Pc @ F, O)           # tr(Pc F O_k)
    OP = O @ Pc                                       # O_k Pc
    Q = np.einsum('kij,lij->kl', OP.conj(), OP)     # <O_k Pc, O_l Pc>
    w = np.linalg.eigvalsh(Q)
    Q = Q + 1e-13*max(w)*np.eye(len(Q))
    # sup_g |c.g|^2 / g^* Q g = c^T Q^{-1} conj(c)... with g complex: = conj(c)^T Q^{-T}... use c^H? compute directly
    x = np.linalg.solve(Q.T, c)   # maximiser g = Q^{-1} conj(c) ; value = c^T Q^{-1} conj(c)
    return float(np.real(c @ np.linalg.solve(Q, c.conj())))

def projP(Z, N, p):
    Qm, _ = np.linalg.qr(Z.reshape(N*N, p))
    return Qm @ Qm.conj().T

def run(N, p, restarts, seed=0, fixP=None):
    rng = np.random.default_rng(seed); F = swapmat(N); best = (-1, None)
    nG = 2*N*N; nP = 2*N*N*p
    def unpackz(z):
        Gp = (z[:N*N] + 1j*z[N*N:nG]).reshape(N, N)
        if fixP is not None: return Gp, fixP
        Z = z[nG:nG+N*N*p] + 1j*z[nG+N*N*p:]
        return Gp, projP(Z, N, p)
    for r in range(restarts):
        z0 = rng.standard_normal(nG + (0 if fixP is not None else nP))
        f = lambda z: -kap_fixed(unpackz(z)[0], np.eye(N*N) - unpackz(z)[1], F, N)
        s = minimize(f, z0, method='BFGS', options=dict(maxiter=3000, gtol=1e-10))
        v = -s.fun
        if v > best[0]: best = (v, s.x)
        print(f" r{r} kappa={v:.8f}", flush=True)
    return best, unpackz

if __name__ == "__main__":
    N, p, R = map(int, sys.argv[1:4])
    best, up = run(N, p, R)
    m = max(k for k in range(1, 50) if k*(k-1)//2 <= p)
    print(f"RESULT N={N} p={p} kappa_max={best[0]:.8f}  law p+m={p+m}")
