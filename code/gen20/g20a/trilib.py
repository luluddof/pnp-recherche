#!/usr/bin/env python3
"""G20a: feasible tuples (p_0, p_1, ..., p_{D-1}) for a strand tensor X (product of identities on strands) with D cuts.
X = A_0 + A_1 + ... + A_{D-1},  rank of A_j across cut j <= p_j.
Solver (design of G19c/starlib): A_0 eliminated by Eckart-Young in flattening 0 (variable projection),
A_j = U_j V_j^T (j >= 1) in flattening j; Levenberg-Marquardt with the projected Gauss-Newton matrix
J^T Pi_N J (Pi_N = normal space of the rank-p_0 manifold), matrix-free CG; QR gauge balancing at each accepted step.

Box tensors:
  D = 2 crossing (G18g): strands s0=(a0,b0), s1=(a1,b1); legs order a0 a1 b0 b1;
      cut 0 (rows) = (a0 a1 | b0 b1), cut 1 (cols) = (a0 b1 | b0 a1).
  D = 3 box (G18g 1bis): 4 body diagonals = strands delta=0..3 ; legs a0 a1 a2 a3 b0 b1 b2 b3;
      orientation of strand delta=(d1 d2) in cut 0: 0, cut 1: d1, cut 2: d2.
      cut 0 = (a0 a1 a2 a3 | b0 b1 b2 b3), cut 1 = (a0 a1 b2 b3 | b0 b1 a2 a3), cut 2 = (a0 b1 a2 b3 | b0 a1 b2 a3).
"""
import numpy as np, time, itertools

class Box:
    def __init__(self, dims, cuts):
        """dims: list of strand dims; legs = a_0..a_{S-1}, b_0..b_{S-1}; cuts: list of row-leg index lists."""
        self.dims = list(dims); S = len(dims); self.S = S
        self.legdims = self.dims + self.dims
        L = 2 * S
        X = np.ones([1] * 0)
        X = np.einsum(','.join(f'{chr(97+i)}{chr(97+S+i)}' for i in range(S)) + '->' + ''.join(chr(97+i) for i in range(L)),
                      *[np.eye(d) for d in dims])
        self.X = X
        self.cuts = []
        for rows in cuts:
            cols = [l for l in range(L) if l not in rows]
            perm = list(rows) + cols
            inv = np.argsort(perm)
            R = int(np.prod([self.legdims[l] for l in rows])); C = int(np.prod([self.legdims[l] for l in cols]))
            self.cuts.append((perm, inv, R, C))
        self.D = len(cuts)
    def flat(self, T, j):
        perm, inv, R, C = self.cuts[j]
        return np.ascontiguousarray(T.transpose(perm)).reshape(R, C)
    def unflat(self, Mx, j):
        perm, inv, R, C = self.cuts[j]
        shp = [self.legdims[l] for l in perm]
        return np.ascontiguousarray(Mx.reshape(shp).transpose(inv))

def crossing(N):
    return Box([N, N], [[0, 1], [0, 3]])

def box3(dims):
    # legs: a0 a1 a2 a3 = 0..3, b0..b3 = 4..7
    return Box(dims, [[0, 1, 2, 3], [0, 1, 6, 7], [0, 5, 2, 7]])

class Prob:
    """X = A_0 + sum_{j>=1} U_j V_j^T (flattening j);  minimise tail_{p0}(flat_0(X - B))."""
    def __init__(self, box, ps):
        self.box = box; self.ps = list(ps); self.X = box.X
        self.shapes = []
        for j in range(1, box.D):
            perm, inv, R, C = box.cuts[j]
            self.shapes += [(R, ps[j]), (C, ps[j])]
        self.sizes = [int(np.prod(s)) for s in self.shapes]; self.npar = sum(self.sizes)
    def unpack(self, x):
        o = 0; out = []
        for s, k in zip(self.shapes, self.sizes):
            out.append(x[o:o+k].reshape(s)); o += k
        return out
    def pack(self, T): return np.concatenate([t.ravel() for t in T]) if T else np.zeros(0)
    def build(self, T):
        B = np.zeros_like(self.X)
        for j in range(1, self.box.D):
            U, V = T[2*(j-1)], T[2*(j-1)+1]
            B += self.box.unflat(U @ V.T, j)
        return B
    def jvec(self, T, dT):
        B = np.zeros_like(self.X)
        for j in range(1, self.box.D):
            U, V = T[2*(j-1)], T[2*(j-1)+1]; dU, dV = dT[2*(j-1)], dT[2*(j-1)+1]
            B += self.box.unflat(dU @ V.T + U @ dV.T, j)
        return B
    def jtvec(self, T, G):
        out = []
        for j in range(1, self.box.D):
            U, V = T[2*(j-1)], T[2*(j-1)+1]
            Gj = self.box.flat(G, j)
            out += [(Gj @ V).ravel(), (Gj.T @ U).ravel()]
        return np.concatenate(out)
    def regauge(self, x):
        T = [t.copy() for t in self.unpack(x)]
        for j in range(len(T) // 2):
            U, V = T[2*j], T[2*j+1]
            if U.shape[1] == 0: continue
            # balanced SVD gauge
            Qu, Ru = np.linalg.qr(U); Qv, Rv = np.linalg.qr(V)
            W, s, Zt = np.linalg.svd(Ru @ Rv.T)
            sq = np.sqrt(s)
            T[2*j] = Qu @ W * sq; T[2*j+1] = Qv @ Zt.T * sq
        return self.pack(T)
    def random_x(self, rng, scale=1.0):
        return scale * rng.normal(size=self.npar) / np.sqrt(max(1, max(s[0] for s in self.shapes) if self.shapes else 1))

def tail_state(prob, x):
    T = prob.unpack(x); p0 = prob.ps[0]
    R = prob.box.flat(prob.X - prob.build(T), 0)
    U, s, Vt = np.linalg.svd(R, full_matrices=False)
    U = U[:, :p0]; V = Vt[:p0].T
    E = R - (R @ V) @ V.T
    return T, R, U, s, V, E, float((E**2).sum())

def lm_solve(prob, x0, maxit=500, tol=1e-24, cg_max=60, cg_tol=1e-3, mu0=1e-2, tlim=1e9, verbose=False):
    t0 = time.time(); p0 = prob.ps[0]; box = prob.box
    x = prob.regauge(x0); T, R, U, s, V, E, f = tail_state(prob, x); mu = mu0; hist = [f]
    for it in range(maxit):
        if f < tol or time.time() - t0 > tlim: break
        g = prob.jtvec(T, box.unflat(E, 0))
        def H(v):
            Xm = box.flat(prob.jvec(T, prob.unpack(v)), 0)
            if p0 > 0:
                UX = U.T @ Xm; XV = Xm @ V
                Xm = Xm - U @ UX - XV @ V.T + U @ ((UX @ V) @ V.T)
            return prob.jtvec(T, box.unflat(Xm, 0)) + mu * v
        dx = np.zeros_like(x); res = g.copy(); p = res.copy(); rs = res @ res; rs0 = rs
        for k in range(cg_max):
            Hp = H(p); al = rs / (p @ Hp); dx += al * p; res -= al * Hp
            rsn = res @ res
            if rsn < cg_tol**2 * rs0: break
            p = res + (rsn / rs) * p; rs = rsn
        xn = x + dx
        st = tail_state(prob, xn)
        if st[-1] < f:
            x = prob.regauge(xn); T, R, U, s, V, E, f = tail_state(prob, x)
            mu = max(mu * 0.3, 1e-14)
        else:
            mu = min(mu * 8.0, 1e10)
        hist.append(f)
        if verbose and it % 20 == 0:
            print(f"   it {it} f={f:.4e} mu={mu:.1e} cg={k+1} |x|={np.linalg.norm(x):.2e} t={time.time()-t0:.0f}s", flush=True)
        if mu >= 1e10: break
        if len(hist) > 60 and hist[-60] - f < 1e-7 * max(f, 1e-30) and f > 1e-8: break
    return x, f, hist

def solve_restarts(prob, nstart, seed=0, scale=1.0, **kw):
    rng = np.random.default_rng(seed); out = []
    for i in range(nstart):
        x0 = prob.random_x(rng, scale)
        x, f, h = lm_solve(prob, x0, **kw)
        out.append((f, len(h), float(np.linalg.norm(x)), x))
    return out
