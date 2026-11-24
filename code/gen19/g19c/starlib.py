#!/usr/bin/env python3
"""G19c: 'three stars + row-rank r' for the hexagon permutation operator, general strand dims n=(n01,n02,n12).
Sites j=0,1,2 (columns), d_j = product of the two strand dims at j, K_j = d_j^2 = dim Op(H_j).
P = S_01 (x) S_02 (x) S_12 as a D x D matrix (D = d0 d1 d2), rows = (o0,o1,o2), cols = (i0,i1,i2).
Operator-space view: T3[p0,p1,p2], p_j = (o_j,i_j).
Star centred at c with leaves l,m:  T3 = sum_{a,b} C[p_c,a,b] LA[p_l,a] LB[p_m,b],  a < bond(c,l) = n_{cl}^2, b < bond(c,m) = n_{cm}^2
(= a tensor of multilinear rank (*, n_cl^2, n_cm^2); the bond opposite to c is 1).
Model: P = A + star_0 + star_1 + star_2, rank(A) <= r as a D x D matrix (row cut).
Solver: variable projection (A eliminated by Eckart-Young), f(B) = sum_{k>r} sigma_k(P-B)^2,
Levenberg-Marquardt with the projected Gauss-Newton matrix J^T Pi_N J (Pi_N = normal space of the rank-r manifold), matrix-free CG."""
import numpy as np, itertools, time
from scipy.linalg import svd as _svd

class Inst:
    def __init__(self, n01, n02, n12):
        self.n = (n01, n02, n12)
        d = (n01*n02, n01*n12, n02*n12); self.d = d
        self.K = tuple(x*x for x in d); self.D = d[0]*d[1]*d[2]
        P = np.zeros(d + d)
        for a, a2, b, b2, c, c2 in itertools.product(range(n01), range(n02), range(n01), range(n12), range(n02), range(n12)):
            # in: site0=(a,a2) site1=(b,b2) site2=(c,c2); swaps a<->b (01), a2<->c (02), b2<->c2 (12)
            P[b*n02+c, a*n12+c2, a2*n12+b2, a*n02+a2, b*n12+b2, c*n12+c2] = 1.0
        self.P = P.reshape(self.D, self.D)
        b01, b02, b12 = n01*n01, n02*n02, n12*n12
        # (centre, leafA, leafB, bondA, bondB)
        self.stars = [(0, 1, 2, b01, b02), (1, 0, 2, b01, b12), (2, 0, 1, b02, b12)]
        self.shapes = []
        for (c, l, m, ba, bb) in self.stars:
            self.shapes += [(self.K[c], ba, bb), (self.K[l], ba), (self.K[m], bb)]
        self.sizes = [int(np.prod(s)) for s in self.shapes]; self.npar = sum(self.sizes)
        self.prod_r = (n01*(n01-1)//2)*(n02*(n02-1)//2)*(n12*(n12-1)//2)
    def to3(self, M):
        d = self.d
        return np.ascontiguousarray(M.reshape(d + d).transpose(0, 3, 1, 4, 2, 5)).reshape(self.K)
    def tomat(self, T3):
        d = self.d
        return np.ascontiguousarray(T3.reshape(d[0], d[0], d[1], d[1], d[2], d[2]).transpose(0, 2, 4, 1, 3, 5)).reshape(self.D, self.D)
    def unpack(self, x):
        o = 0; out = []
        for s, k in zip(self.shapes, self.sizes):
            out.append(x[o:o+k].reshape(s)); o += k
        return out
    # ---- star maps (in operator-space ordering (p0,p1,p2)) ----
    @staticmethod
    def _arr(c, t):   # t has axes (centre, leafA, leafB) -> (p0,p1,p2)
        if c == 0: return t
        if c == 1: return t.transpose(1, 0, 2)
        return t.transpose(1, 2, 0)
    @staticmethod
    def _unarr(c, G):  # (p0,p1,p2) -> (centre, leafA, leafB)
        if c == 0: return G
        if c == 1: return G.transpose(1, 0, 2)
        return G.transpose(2, 0, 1)
    def build3(self, T):
        out = np.zeros(self.K)
        for s, (c, l, m, ba, bb) in enumerate(self.stars):
            C, LA, LB = T[3*s:3*s+3]
            t = np.tensordot(C, LB, axes=([2], [1]))          # (c,a,m)
            t = np.tensordot(t, LA, axes=([1], [1]))          # (c,m,l)
            out += self._arr(c, t.transpose(0, 2, 1))
        return out
    def star3(self, T, s):
        c, l, m, ba, bb = self.stars[s]; C, LA, LB = T[3*s:3*s+3]
        t = np.tensordot(np.tensordot(C, LB, axes=([2], [1])), LA, axes=([1], [1]))
        return self._arr(c, t.transpose(0, 2, 1))
    def jvec3(self, T, dT):
        out = np.zeros(self.K)
        for s, (c, l, m, ba, bb) in enumerate(self.stars):
            C, LA, LB = T[3*s:3*s+3]; dC, dLA, dLB = dT[3*s:3*s+3]
            t1 = np.tensordot(dC, LB, axes=([2], [1])) + np.tensordot(C, dLB, axes=([2], [1]))   # (c,a,m)
            t = np.tensordot(t1, LA, axes=([1], [1]))
            t += np.tensordot(np.tensordot(C, LB, axes=([2], [1])), dLA, axes=([1], [1]))
            out += self._arr(c, t.transpose(0, 2, 1))
        return out
    def jtvec3(self, T, G3):
        g = []
        for s, (c, l, m, ba, bb) in enumerate(self.stars):
            C, LA, LB = T[3*s:3*s+3]
            G = np.ascontiguousarray(self._unarr(c, G3))       # (c,l,m)
            GB = np.tensordot(G, LB, axes=([2], [0]))          # (c,l,b)
            gC = np.tensordot(GB, LA, axes=([1], [0])).transpose(0, 2, 1)   # (c,b,a)->(c,a,b)
            gLA = np.tensordot(GB, C, axes=([0, 2], [0, 2]))   # (l,a)
            GA = np.tensordot(G, LA, axes=([1], [0]))          # (c,m,a)
            gLB = np.tensordot(GA, C, axes=([0, 2], [0, 1]))   # (m,b)
            g += [gC.ravel(), gLA.ravel(), gLB.ravel()]
        return np.concatenate(g)
    def regauge(self, x):
        """orthonormalise the leaves (QR), absorb into the cores.  Does not change B."""
        T = [t.copy() for t in self.unpack(x)]
        for s in range(3):
            C, LA, LB = T[3*s:3*s+3]
            QA, RA = np.linalg.qr(LA); QB, RB = np.linalg.qr(LB)
            C = np.einsum('cab,xa,yb->cxy', C, RA, RB)
            T[3*s:3*s+3] = [C, QA, QB]
        return np.concatenate([t.ravel() for t in T])
    def random_x(self, rng, core_scale=0.1):
        T = []
        for s, (c, l, m, ba, bb) in enumerate(self.stars):
            T += [core_scale * rng.normal(size=(self.K[c], ba, bb)),
                  np.linalg.qr(rng.normal(size=(self.K[l], ba)))[0], np.linalg.qr(rng.normal(size=(self.K[m], bb)))[0]]
        return np.concatenate([t.ravel() for t in T])
    def product_point(self):
        """exact product point: P = prod_e (1 - Y_e); stars absorb every term with >=1 identity factor. returns x (stars) and A."""
        n01, n02, n12 = self.n; d = self.d
        def swap(n):
            S = np.zeros((n, n, n, n))
            for a in range(n):
                for b in range(n): S[b, a, a, b] = 1
            return S  # [o_left, o_right, i_left, i_right]
        def ident(n): return np.einsum('ac,bd->abcd', np.eye(n), np.eye(n))
        S01, S02, S12 = swap(n01), swap(n02), swap(n12)
        # NOTE P's convention: out site0=(b,c), out site1=(a,c2), out site2=(a2,b2): the swap is 'in the strands'; with the
        # relabelling of outputs used in P, a rectangle acts as IDENTITY wiring in[site j strand] -> out[site k strand].
        def assemble(F01, F02, F12):
            # F_e[x_out_at_j, x_out_at_k, x_in_at_j, x_in_at_k] for edge e=(j,k), j<k
            # site0 in (a,a2) out (u,u2); site1 in (b,b2) out (v,v2); site2 in (c,c2) out (w,w2)
            T = np.einsum('uvab,UwAc,VWBC->uUvVwWaAbBcC', F01, F02, F12)
            return T.reshape(self.D, self.D)
        return swap, ident, assemble

def svd_top(R, r):
    if R.shape[0] < 200:
        U, s, Vt = _svd(R, full_matrices=False, lapack_driver='gesdd', check_finite=False)
        return U[:, :r], s, Vt[:r].T
    from scipy.linalg import eigh
    w, V = eigh(R.T @ R, check_finite=False)
    V = V[:, ::-1][:, :r]; s = np.sqrt(np.maximum(w[::-1], 0))
    U = (R @ V) / np.maximum(s[:r], 1e-300)
    U, _ = np.linalg.qr(U) if r > 0 else (U, None)   # orthonormal basis of the left space (only the projector is used)
    return U, s, V

def lm_solve(inst, r, x0, maxit=200, tol=1e-13, cg_max=40, cg_tol=1e-2, mu0=1e-2, tlim=1e9, verbose=False, ridge=0.0, log=None):
    """minimise f(x) = sum_{k>r} sigma_k(P - B(x))^2 (+ ridge*|x|^2).  returns x, f, history"""
    t0 = time.time(); P = inst.P
    def state(x):
        T = inst.unpack(x)
        R = P - inst.tomat(inst.build3(T))
        U, s, V = svd_top(R, r)
        E = R - (R @ V) @ V.T
        f = float((E**2).sum()) + ridge * float(x @ x)
        return T, R, U, s, V, E, f
    x = inst.regauge(x0); T, R, U, s, V, E, f = state(x); mu = mu0; hist = [f]
    for it in range(maxit):
        if f < tol or time.time() - t0 > tlim: break
        g = inst.jtvec3(T, inst.to3(E)) - ridge * x      # = -grad/2
        def H(v):
            X = inst.tomat(inst.jvec3(T, inst.unpack(v)))
            if r > 0:
                UX = U.T @ X; XV = X @ V
                X = X - U @ UX - XV @ V.T + U @ ((UX @ V) @ V.T)
            return inst.jtvec3(T, inst.to3(X)) + (mu + ridge) * v
        # CG
        dx = np.zeros_like(x); res = g.copy(); p = res.copy(); rs = res @ res; rs0 = rs
        for k in range(cg_max):
            Hp = H(p); al = rs / (p @ Hp); dx += al * p; res -= al * Hp
            rsn = res @ res
            if rsn < cg_tol**2 * rs0: break
            p = res + (rsn / rs) * p; rs = rsn
        xn = x + dx
        Tn, Rn, Un, sn, Vn, En, fn = state(xn)
        if fn < f:
            ratio = (f - fn) / max(dx @ (g + mu * dx + res), 1e-300) if False else 1.0
            x = inst.regauge(xn); T, R, U, s, V, E, f = state(x) if False else (inst.unpack(x), Rn, Un, sn, Vn, En, fn)
            mu = max(mu * 0.4, 1e-12)
        else:
            mu = min(mu * 6.0, 1e8)
        hist.append(f)
        if verbose or log:
            msg = f"   it {it} f={f:.6e} mu={mu:.1e} cg={k+1} |dx|={np.linalg.norm(dx):.2e} t={time.time()-t0:.0f}s"
            if verbose: print(msg, flush=True)
        if mu >= 1e8: break
        if len(hist) > 30 and hist[-30] - f < 1e-9 * max(f, 1e-30) and f > 1e-6: break
    return x, f, hist

def norms(inst, x, r):
    T = inst.unpack(x); out = []
    B = np.zeros(inst.K)
    for s in range(3):
        t = inst.star3(T, s); out.append(float(np.linalg.norm(t))); B += t
    R = inst.P - inst.tomat(B); sv = np.linalg.svd(R, compute_uv=False)
    return dict(star_norms=out, B_norm=float(np.linalg.norm(B)), A_norm=float(np.sqrt((sv[:r]**2).sum())), sig_r=float(sv[r-1]) if r > 0 else None,
                sig_r1=float(sv[r]), tail=float((sv[r:]**2).sum()))
