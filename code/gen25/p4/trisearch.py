# search for a true-N solution whose P_b E-family (side=E) or C-family (side=C) is upper triangular.
# By the symmetry (X (x) Y) . (Y^-1 (x) X^-1) of the target, "triangularisable" = "upper triangular after a gauge", so WLOG triangular.
import sys, numpy as np, time
sys.path.insert(0, '/home/user/pnp-recherche/code/gen24/g24f')
import numsearch
from modes import Setup, npf
from scipy.optimize import least_squares
n, E0, N, side, trials, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), int(sys.argv[6])
W = float(sys.argv[7]) if len(sys.argv) > 7 else 1.0
St = Setup(n, E0)
Bn, On, dA, dB, dM = numsearch.prep(n, E0, 0, 0)
gs = [tuple(int(x) for x in p) for (p, m) in St.A]
Dn = np.concatenate([np.array([npf(St.mode(g, E0, E)) for g in gs]) for E in range(E0)], axis=1)  # g, m, w
Up = np.concatenate([np.array([npf(St.mode(g, E, E0)) for g in gs]) for E in range(E0)], axis=2)  # g, w, m
T0 = np.einsum('xv,uy->xvuy', np.eye(N), np.eye(N))
low = np.tril(np.ones((N, N)), -1).astype(bool)
sizes = [N*dA, N*dB, N*dM, dM*N]; shapes = [(N, dA), (N, dB), (N, dM), (dM, N)]
def unpack(z):
    o = 0; out = []
    for s, sh in zip(sizes, shapes): out.append(z[o:o+s].reshape(sh)); o += s
    return out
def fun(z):
    A, Bm, U, Y = unpack(z)
    beta = np.einsum('xi,gij,vj->gxv', A, Bn, Bm); Om = np.einsum('uk,gkl,ly->guy', U, On, Y)
    r1 = (np.einsum('gxv,guy->xvuy', beta, Om) - T0).ravel()
    if side == 'E': F = np.einsum('xg,gmw,wy->mxy', A, Dn, Y)
    else: F = np.einsum('vg,uw,gwm->mvu', Bm, U, Up)
    return np.concatenate([r1, W*F[:, low].ravel(), 1e-7*z])
def jac(z):
    A, Bm, U, Y = unpack(z); I = np.eye(N)
    BB = np.einsum('gij,vj->giv', Bn, Bm); AB = np.einsum('xi,gij->gxj', A, Bn)
    beta = np.einsum('xi,giv->gxv', A, BB); Om = np.einsum('uk,gkl,ly->guy', U, On, Y)
    OY = np.einsum('gkl,ly->gky', On, Y); UO = np.einsum('uk,gkl->gul', U, On)
    JA = np.einsum('xa,giv,guy->xvuyai', I, BB, Om).reshape(N**4, -1)
    JB = np.einsum('gxj,vb,guy->xvuybj', AB, I, Om).reshape(N**4, -1)
    JU = np.einsum('gxv,uc,gky->xvuyck', beta, I, OY).reshape(N**4, -1)
    JY = np.einsum('gxv,guk,yd->xvuykd', beta, UO, I).reshape(N**4, -1)
    J1 = np.hstack([JA, JB, JU, JY]); k = Dn.shape[1] if side == 'E' else Up.shape[2]
    Z = lambda sh: np.zeros((k, N, N) + sh)
    if side == 'E':
        fA = np.einsum('xa,gmw,wy->mxyag', I, Dn, Y); fY = np.einsum('xg,gmw,yc->mxywc', A, Dn, I)
        blocks = [fA.reshape(k, N, N, -1), Z((N*dB,)), Z((N*dM,)), fY.reshape(k, N, N, -1)]
    else:
        fB = np.einsum('vb,uw,gwm->mvubg', I, U, Up); fU = np.einsum('vg,uc,gwm->mvucw', Bm, I, Up)
        blocks = [Z((N*dA,)), fB.reshape(k, N, N, -1), fU.reshape(k, N, N, -1), Z((dM*N,))]
    J2 = np.concatenate(blocks, axis=3)[:, low].reshape(-1, len(z))
    return np.vstack([J1, W*J2, 1e-7*np.eye(len(z))])
rng = np.random.default_rng(seed); best = 1e9
for t in range(trials):
    t0 = time.time()
    r = least_squares(fun, rng.standard_normal(sum(sizes)), jac=jac, method='lm', xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=20000)
    A, Bm, U, Y = unpack(r.x); f = fun(r.x)
    r1 = np.linalg.norm(f[:N**4]); r2 = np.linalg.norm(f[N**4:-len(r.x)])
    print("trial %d: identity residual %.2e  triangularity residual %.2e  maxabs %.2f (%.0fs)" % (t, r1, r2, np.abs(r.x).max(), time.time()-t0), flush=True)
