# Spread measurements for decompositions F = R + C (F = swap, X[x,y,u,v] = [x=v][y=u]).
# Fl1: (x,y)|(u,v); Fl2: (x,u)|(y,v).  Group G = GL x GL: T -> (a (x) b) T (b^-1 (x) a^-1) fixes F, preserves rank Fl1 / rank Fl2.
import numpy as np, scipy.sparse.linalg as sla
def swapT(N):
    F = np.zeros((N, N, N, N))
    for x in range(N):
        for y in range(N): F[x, y, y, x] = 1
    return F
Fl1 = lambda T: T.reshape(T.shape[0] ** 2, -1)
Fl2 = lambda T: T.transpose(0, 2, 1, 3).reshape(T.shape[0] ** 2, -1)
def unFl2(M, N): return M.reshape(N, N, N, N).transpose(0, 2, 1, 3)
def act(T, a, b):
    ai, bi = np.linalg.inv(a), np.linalg.inv(b)
    return np.einsum('xX,yY,XYUV,Uu,Vv->xyuv', a, b, T, bi, ai, optimize=True)
def opn(M):
    if M.shape[0] <= 1200: return np.linalg.norm(M, 2)
    return sla.svds(M, k=1, return_singular_vectors=False, tol=1e-8)[0]
def rk(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False); return int((s > max(tol * s[0], 1e-9)).sum())
def specrad(M):
    if M.shape[0] <= 1500: return np.abs(np.linalg.eigvals(M)).max()
    return np.abs(sla.eigs(M, k=1, which='LM', return_eigenvectors=False, tol=1e-8)).max()
def analyze(R, C, p=None, q=None, inv=True, check=True):
    N = R.shape[0]; S1 = Fl1(swapT(N)); S2 = Fl2(swapT(N))
    if check:
        err = np.abs(R + C - swapT(N)).max(); assert err < 1e-8, err
    M1, M2 = Fl1(R), Fl2(C)
    if p is None: p = rk(M1)
    if q is None: q = rk(M2)
    oR, oC = opn(M1), opn(M2)
    d = dict(N=N, p=p, q=q, opR=oR, opC=oC, frR=np.linalg.norm(M1), frC=np.linalg.norm(M2),
             FR=float((R * swapT(N)).sum()), FC=float((C * swapT(N)).sum()))
    d['L'] = max(oR / np.sqrt(max(q, 1)), oC / np.sqrt(max(p, 1)))
    d['ineq_op'] = N * N / (p * oR + q * oC)          # <= 1 (proved); 1 = tight
    d['ineq_fr'] = N * N / (np.sqrt(p) * d['frR'] + np.sqrt(q) * d['frC'])
    if inv:   # G-invariant lower bounds: rho(Fl1 R . S1), rho(Fl2 C . S2^T)  (similarity invariants)
        d['rhoR'] = specrad(M1 @ S1); d['rhoC'] = specrad(M2 @ S2.T)
        d['Linv'] = max(d['rhoR'] / np.sqrt(max(q, 1)), d['rhoC'] / np.sqrt(max(p, 1)))
    return d
def fmt(d):
    s = f"N={d['N']:4d} p={d['p']:6d} q={d['q']:6d} |R|op={d['opR']:10.4g} |C|op={d['opC']:10.4g} L={d['L']:9.4g}"
    s += f" N2/(p|R|+q|C|)={d['ineq_op']:.4f} frob={d['ineq_fr']:.4f}"
    if 'Linv' in d: s += f" rhoR={d['rhoR']:.4g} rhoC={d['rhoC']:.4g} Linv={d['Linv']:.4g}"
    return s
# ---- minimisation over G (and optionally over Z = sum_k rank-one product operators moved R -> C) ----
def schatten(M, m):
    U, s, Vt = np.linalg.svd(M, full_matrices=False)
    sm = s / s[0]; val = s[0] * (sm ** m).sum() ** (1 / m)
    w = sm ** (m - 1) * (sm ** m).sum() ** (1 / m - 1)
    return val, (U * w) @ Vt, s[0]
def grads_XY(G4, T):
    gX = np.einsum('xyuv,zyuv->xz', G4, T, optimize=True) - np.einsum('xyuv,xyuw->wv', G4, T, optimize=True)
    gY = np.einsum('xyuv,xzuv->yz', G4, T, optimize=True) - np.einsum('xyuv,xywv->wu', G4, T, optimize=True)
    return gX, gY
def objective(R, C, p, q, m, beta):
    N = R.shape[0]
    vR, GR, sR = schatten(Fl1(R), m); vC, GC, sC = schatten(Fl2(C), m)
    a1, a2 = np.log(vR / np.sqrt(q)), np.log(vC / np.sqrt(p))
    mx = max(a1, a2); w1, w2 = np.exp(beta * (a1 - mx)), np.exp(beta * (a2 - mx))
    f = mx + np.log(w1 + w2) / beta; w1, w2 = w1 / (w1 + w2), w2 / (w1 + w2)
    G4R = GR.reshape(N, N, N, N) * (w1 / vR); G4C = unFl2(GC, N) * (w2 / vC)
    Lt = max(sR / np.sqrt(q), sC / np.sqrt(p))
    return f, G4R, G4C, Lt
def minimise(R, C, p, q, iters=300, m=24, beta=30, lr=0.05, nz=0, seed=0, verbose=False, init=0.0):
    """minimise L over (a,b) in GL x GL (+ nz product rank-one moves). returns (best true L, R', C', history)"""
    N = R.shape[0]; rng = np.random.default_rng(seed)
    a, b = np.eye(N) + init * rng.standard_normal((N, N)), np.eye(N) + init * rng.standard_normal((N, N))
    Zf = [rng.standard_normal((4, N)) * 1e-3 for _ in range(nz)]
    def build(a, b, Zf):
        R1, C1 = act(R, a, b), act(C, a, b)
        for al, be, ga, de in Zf:
            Z = np.einsum('x,y,u,v->xyuv', al, be, ga, de); R1 = R1 + Z; C1 = C1 - Z
        return R1, C1
    R1, C1 = build(a, b, Zf); f, G4R, G4C, Lt = objective(R1, C1, p + nz, q + nz, m, beta)
    best = (Lt, R1, C1); hist = [Lt]
    for it in range(iters):
        gXR, gYR = grads_XY(G4R, R1); gXC, gYC = grads_XY(G4C, C1)
        gX, gY = gXR + gXC, gYR + gYC
        gZ = []
        for al, be, ga, de in Zf:
            Gd = G4R - G4C
            gZ.append(np.array([np.einsum('xyuv,y,u,v->x', Gd, be, ga, de), np.einsum('xyuv,x,u,v->y', Gd, al, ga, de),
                                np.einsum('xyuv,x,y,v->u', Gd, al, be, de), np.einsum('xyuv,x,y,u->v', Gd, al, be, ga)]))
        nrm = np.sqrt((gX ** 2).sum() + (gY ** 2).sum() + sum((g ** 2).sum() for g in gZ)) + 1e-300
        step = lr
        while step > 1e-6:
            from scipy.linalg import expm
            a2 = expm(-step * gX / nrm) @ a; b2 = expm(-step * gY / nrm) @ b
            Zf2 = [z - step * g / nrm for z, g in zip(Zf, gZ)]
            R2, C2 = build(a2, b2, Zf2); f2, G4R2, G4C2, Lt2 = objective(R2, C2, p + nz, q + nz, m, beta)
            if f2 < f: break
            step /= 2
        if step <= 1e-6: break
        a, b, Zf, R1, C1, f, G4R, G4C, Lt = a2, b2, Zf2, R2, C2, f2, G4R2, G4C2, Lt2
        lr = min(step * 1.5, 0.5); hist.append(Lt)
        if Lt < best[0]: best = (Lt, R1, C1)
        if verbose and it % 20 == 0: print(it, f, Lt, flush=True)
    return best[0], best[1], best[2], hist
