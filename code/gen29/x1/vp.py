# Variable-projection numerics: min over C (column format) of sum_{i>p} sigma_i(T - C)^2, row flattening.
# Operators: 'rect' X_{N1,N2}[x,y,u,v]=[x=v][y=u] (rows (x,y)|(u,v), cols (x,u)|(y,v), C = sum_l G_l(x,u) H_l(y,v));
#            'cyc3' P_3[a1,a2,a3,b1,b2,b3]=[a1=b2][a2=b3][a3=b1] (rows a|b, cols = ring L,M,R, node i holds (a_i,b_i)),
#             C = sum f1[a1,b1,t31,t12] f2[a2,b2,t12,t23] f3[a3,b3,t23,t31]  (bonds b12,b23,b31).
#            'copy' T[x,y,u,v,w]=[x=v][y=u][x=w], w in bottom-left cell: rows (x,y)|(u,v,w), cols (x,u,w)|(y,v), OSR q.
import numpy as np, sys
from scipy.optimize import minimize

def rect_T(N1, N2):
    T = np.zeros((N1, N2, N2, N1))
    for x in range(N1):
        for y in range(N2): T[x, y, y, x] = 1
    return T

class Rect:
    def __init__(s, N1, N2, q):
        s.N1, s.N2, s.q = N1, N2, q; s.T = rect_T(N1, N2)
        s.nG = q * N1 * N2; s.n = 2 * s.nG
    def C(s, z):
        G = z[:s.nG].reshape(s.q, s.N1, s.N2); H = z[s.nG:].reshape(s.q, s.N2, s.N1)
        return np.einsum('lxu,lyv->xyuv', G, H), G, H
    def rows(s, A): return A.reshape(s.N1 * s.N2, s.N2 * s.N1)
    def grad(s, z, D):  # D = dF/dC (shape of T)
        _, G, H = s.C(z)
        gG = np.einsum('xyuv,lyv->lxu', D, H); gH = np.einsum('xyuv,lxu->lyv', D, G)
        return np.concatenate([gG.ravel(), gH.ravel()])

class Cyc3:
    def __init__(s, N, b):
        s.N, s.b = N, b; b12, b23, b31 = b
        T = np.zeros((N,) * 6)
        for a1 in range(N):
            for a2 in range(N):
                for a3 in range(N): T[a1, a2, a3, a3, a1, a2] = 1  # b1=a3, b2=a1, b3=a2
        s.T = T
        s.sh = [(N, N, b31, b12), (N, N, b12, b23), (N, N, b23, b31)]
        s.sz = [int(np.prod(x)) for x in s.sh]; s.n = sum(s.sz)
    def parts(s, z):
        o = 0; F = []
        for sh, sz in zip(s.sh, s.sz): F.append(z[o:o + sz].reshape(sh)); o += sz
        return F
    def C(s, z):
        f1, f2, f3 = s.parts(z)
        return np.einsum('ipcm,jqmn,krnc->ijkpqr', f1, f2, f3), f1, f2, f3
    def rows(s, A): return A.reshape(s.N ** 3, s.N ** 3)
    def grad(s, z, D):
        _, f1, f2, f3 = s.C(z)
        g1 = np.einsum('ijkpqr,jqmn,krnc->ipcm', D, f2, f3)
        g2 = np.einsum('ijkpqr,ipcm,krnc->jqmn', D, f1, f3)
        g3 = np.einsum('ijkpqr,ipcm,jqmn->krnc', D, f1, f2)
        return np.concatenate([g1.ravel(), g2.ravel(), g3.ravel()])

class Copy:  # T[x,y,u,v,w] = [x=v][y=u][x=w]; cols (x,u,w)|(y,v): C = sum_l G_l[x,u,w] H_l[y,v]
    def __init__(s, N, q):
        s.N, s.q = N, q; T = np.zeros((N,) * 5)
        for x in range(N):
            for y in range(N): T[x, y, y, x, x] = 1
        s.T = T; s.nG = q * N ** 3; s.n = s.nG + q * N * N
    def C(s, z):
        G = z[:s.nG].reshape(s.q, s.N, s.N, s.N); H = z[s.nG:].reshape(s.q, s.N, s.N)
        return np.einsum('lxuw,lyv->xyuvw', G, H), G, H
    def rows(s, A): return A.reshape(s.N ** 2, s.N ** 3)
    def grad(s, z, D):
        _, G, H = s.C(z)
        return np.concatenate([np.einsum('xyuvw,lyv->lxuw', D, H).ravel(), np.einsum('xyuvw,lxuw->lyv', D, G).ravel()])

def fun(z, prob, p):
    Cz = prob.C(z)[0]; M = prob.rows(prob.T - Cz)
    U, S, Vt = np.linalg.svd(M, full_matrices=False)
    Mp = (U[:, :p] * S[:p]) @ Vt[:p]
    R = M - Mp; f = (R ** 2).sum()
    D = -2 * R.reshape(prob.T.shape)
    return f, prob.grad(z, D)

def solve(prob, p, starts=20, seed=0, scale=0.5, maxiter=3000):
    rng = np.random.default_rng(seed); best = (np.inf, None)
    for st in range(starts):
        z0 = rng.normal(size=prob.n) * scale
        r = minimize(fun, z0, args=(prob, p), jac=True, method='L-BFGS-B',
                     options={'maxiter': maxiter, 'gtol': 1e-14, 'ftol': 1e-16})
        if r.fun < best[0]: best = (r.fun, r.x)
        if r.fun < 1e-20: break
    return best

if __name__ == '__main__':
    kind = sys.argv[1]
    if kind == 'rect':
        N1, N2, q = map(int, sys.argv[2:5]); prob = Rect(N1, N2, q); pmax = N1 * N2
    elif kind == 'cyc3':
        N = int(sys.argv[2]); b = tuple(map(int, sys.argv[3].split(','))); prob = Cyc3(N, b); pmax = N ** 3
    elif kind == 'copy':
        N, q = map(int, sys.argv[2:4]); prob = Copy(N, q); pmax = N * N
    import os; starts = int(os.environ.get('STARTS', '20'))
    for p in range(pmax, -1, -1):
        f, z = solve(prob, p, starts=starts)
        print(kind, sys.argv[2:], 'p =', p, 'tail = %.3e' % f, flush=True)
        if f > 1e-8: break
