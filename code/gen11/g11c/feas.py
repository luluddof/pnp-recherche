# G11c: numerical feasibility of rank profiles in the GENERAL (non-commutative) per-point problem.
# J(t) 2^l x 2^l, deg<=d, J(0)=X^{(x)l}; one point tau_k per direction k.
# Reformulation: unknown off-nice blocks B01^k, B10^k (2^{l-1} x 2^{l-1}, rank<=r_k) of J(tau_k).
# J exists  <=>  for every entry (x,y) with D=D(x,y)={k: x_k!=y_k}, |D|>=d+1, the data
#   (0, [y==~x]) , (tau_k, E_k[x,y])_{k in D}  lie on a polynomial of degree <= d.
import numpy as np, sys, itertools
from scipy.optimize import least_squares
from scipy.linalg import null_space

class Prob:
    def __init__(s, l, d, tau, ranks):
        s.l, s.d, s.tau, s.ranks = l, d, np.array(tau, float), list(ranks)
        n = 1 << l; h = n >> 1; s.n = n; s.h = h
        # parameter layout
        s.off = {}; p = 0
        for k in range(l):
            for side in (0, 1):
                r = ranks[k]
                s.off[(k, side, 'U')] = p; p += h * r
                s.off[(k, side, 'V')] = p; p += h * r
        s.np_ = p
        # entries
        s.eidx = {}; s.elist = []
        def loc(x, k):  # index of x among words with bit k removed
            lo = x & ((1 << k) - 1); hi = x >> (k + 1)
            return lo | (hi << k)
        s.loc = loc
        for k in range(l):
            for x in range(n):
                for y in range(n):
                    if ((x ^ y) >> k) & 1:
                        s.eidx[(k, x, y)] = len(s.elist)
                        s.elist.append((k, (x >> k) & 1, loc(x, k), loc(y, k)))
        ne = len(s.elist)
        rows = []; c = []
        for x in range(n):
            for y in range(n):
                D = [k for k in range(l) if ((x ^ y) >> k) & 1]
                if len(D) < d + 1: continue
                nodes = np.array([0.0] + [s.tau[k] for k in D])
                V = np.vander(nodes, d + 1, increasing=True)
                W = null_space(V.T).T  # (|D|-d) x (|D|+1)
                for w in W:
                    row = np.zeros(ne)
                    for i, k in enumerate(D):
                        row[s.eidx[(k, x, y)]] = w[i + 1]
                    rows.append(row); c.append(w[0] * (1.0 if y == (n - 1) ^ x else 0.0))
        s.W = np.array(rows); s.c = np.array(c); s.nres = len(c)
        s.ek = np.array([e[0] for e in s.elist]); s.es = np.array([e[1] for e in s.elist])
        s.ei = np.array([e[2] for e in s.elist]); s.ej = np.array([e[3] for e in s.elist])
    def blocks(s, p):
        B = {}
        for k in range(s.l):
            r = s.ranks[k]
            for side in (0, 1):
                U = p[s.off[(k, side, 'U')]:s.off[(k, side, 'U')] + s.h * r].reshape(s.h, r)
                V = p[s.off[(k, side, 'V')]:s.off[(k, side, 'V')] + s.h * r].reshape(s.h, r)
                B[(k, side)] = (U, V)
        return B
    def evals(s, p):
        B = s.blocks(p); e = np.zeros(len(s.elist))
        for k in range(s.l):
            for side in (0, 1):
                U, V = B[(k, side)]
                m = (s.ek == k) & (s.es == side)
                M = U @ V.T
                e[m] = M[s.ei[m], s.ej[m]]
        return e
    def res(s, p): return s.W @ s.evals(p) + s.c
    def jac(s, p):
        B = s.blocks(p); ne = len(s.elist); Je = np.zeros((ne, s.np_))
        for k in range(s.l):
            r = s.ranks[k]
            for side in (0, 1):
                U, V = B[(k, side)]
                m = np.where((s.ek == k) & (s.es == side))[0]
                oU = s.off[(k, side, 'U')]; oV = s.off[(k, side, 'V')]
                for j in range(r):
                    Je[m, oU + s.ei[m] * r + j] = V[s.ej[m], j]
                    Je[m, oV + s.ej[m] * r + j] = U[s.ei[m], j]
        return s.W @ Je
    def solve(s, starts=20, seed=0, scale=1.0, verbose=False):
        rng = np.random.default_rng(seed); best = None
        if s.np_ == 0:
            r = s.res(np.zeros(0)); return 0.5 * float(r @ r), np.zeros(0)
        for it in range(starts):
            p0 = rng.normal(size=s.np_) * scale
            try:
                o = s.lm(p0)
            except Exception as ex:
                continue
            if best is None or o.cost < best[0]: best = (o.cost, o.x)
            if verbose: print("  start", it, "cost %.3e" % o.cost, "maxabs %.2e" % np.abs(o.x).max(), flush=True)
            if o.cost < 1e-24: break
        return best
    def lm(s, p, iters=600, tol=1e-26):
        # plain Levenberg-Marquardt (handles the GL(r) gauge degeneracy)
        class O: pass
        lam = 1e-2; r = s.res(p); c = 0.5 * r @ r
        for it in range(iters):
            J = s.jac(p); g = J.T @ r; H = J.T @ J
            while True:
                try:
                    dp = np.linalg.solve(H + lam * np.eye(len(p)), -g)
                except np.linalg.LinAlgError:
                    lam *= 10; continue
                r2 = s.res(p + dp); c2 = 0.5 * r2 @ r2
                if c2 < c:
                    p = p + dp; r = r2; dec = c - c2; c = c2; lam = max(lam / 5, 1e-14); break
                lam *= 4
                if lam > 1e12: dec = 0; break
            if c < tol or lam > 1e12: break
        o = O(); o.cost = c; o.x = p; o.nit = it
        return o
    def J_of(s, p):
        """reconstruct J(t) coefficient matrices (minimal: free entries set by lowest-degree interpolation) and verify"""
        e = s.evals(p); n = s.n; l = s.l; d = s.d
        A = np.zeros((d + 1, n, n))
        for x in range(n):
            for y in range(n):
                D = [k for k in range(l) if ((x ^ y) >> k) & 1]
                nodes = [0.0] + [s.tau[k] for k in D]
                vals = [1.0 if y == (n - 1) ^ x else 0.0] + [e[s.eidx[(k, x, y)]] for k in D]
                deg = min(d, len(nodes) - 1)
                co = np.linalg.lstsq(np.vander(np.array(nodes), deg + 1, increasing=True), np.array(vals), rcond=None)[0]
                A[:deg + 1, x, y] = co
        return A
    def check(s, A, tol=1e-8):
        n = s.n; out = []
        for k in range(s.l):
            Jt = sum(A[i] * s.tau[k] ** i for i in range(s.d + 1))
            r0 = [x for x in range(n) if not (x >> k) & 1]; r1 = [x for x in range(n) if (x >> k) & 1]
            sv0 = np.linalg.svd(Jt[np.ix_(r0, r1)], compute_uv=False); sv1 = np.linalg.svd(Jt[np.ix_(r1, r0)], compute_uv=False)
            out.append((int((sv0 > tol).sum()), int((sv1 > tol).sum())))
        return out

if __name__ == "__main__":
    l = int(sys.argv[1]); d = int(sys.argv[2]); ranks = [int(c) for c in sys.argv[3].split(',')]
    starts = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    tau = [1.0, 1.7, 2.9, 4.3, 5.9, 7.1][:l] if len(sys.argv) <= 5 else [float(v) for v in sys.argv[5].split(',')]
    P = Prob(l, d, tau, ranks)
    print("l=%d d=%d ranks=%s params=%d residuals=%d" % (l, d, ranks, P.np_, P.nres))
    cost, p = P.solve(starts=starts, verbose=True)
    print("best cost %.3e" % cost, "max|param| %.2e" % (np.abs(p).max() if len(p) else 0))
    if cost < 1e-20:
        A = P.J_of(p); print("verified (rank01,rank10) per direction:", P.check(A))
