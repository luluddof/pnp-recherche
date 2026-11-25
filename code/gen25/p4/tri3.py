# numerical test: is the family {F_m} of 3x3 matrices simultaneously triangularisable by equivalence over C?
# N=3 criterion: exists c, r != 0 with rank[F_m c]_m <= 1, rank[r^T F_m]_m <= 1, r^T F_m c = 0 for all m.
import sys, pickle, numpy as np
from scipy.optimize import minimize
def cvec(z, k): return z[:k] + 1j*z[k:2*k]
def tail(M):  # sigma_2^2 + sigma_3^2 relative
    s = np.linalg.svd(M, compute_uv=False); return (s[1:]**2).sum()/max((s**2).sum(), 1e-300)
def span_sv(F):
    return np.linalg.svd(F.reshape(len(F), -1), compute_uv=False)
def test(F, starts=60, seed=0):
    rng = np.random.default_rng(seed)
    F = F/np.sqrt((np.abs(F)**2).sum())
    def parts(z):
        c = cvec(z[:6], 3); r = cvec(z[6:], 3); c = c/np.linalg.norm(c); r = r/np.linalg.norm(r)
        Mc = np.einsum('mab,b->am', F, c); Mr = np.einsum('a,mab->bm', r, F)
        cross = np.einsum('a,mab,b->m', r, F, c)
        nn = (np.abs(Mc)**2).sum() + (np.abs(Mr)**2).sum()
        return tail(Mc), tail(Mr), (np.abs(cross)**2).sum()/max(nn, 1e-300)
    best = {}
    for name, fn in [('col', lambda z: parts(z)[0]), ('row', lambda z: parts(z)[1]), ('all', lambda z: sum(parts(z)))]:
        b = 1e9
        for s in range(starts):
            r = minimize(fn, rng.standard_normal(12), method='BFGS', options=dict(gtol=1e-12, maxiter=3000))
            b = min(b, r.fun)
        best[name] = b
    return best
if __name__ == "__main__":
    res = pickle.load(open(sys.argv[1], 'rb'))
    for i, R in enumerate(res):
        for side in ['Es', 'Cs']:
            F = R[side]; sv = span_sv(F)
            b = test(F, seed=i)
            print("sol %d %s: span dim %d (sv %s) | min col-cond %.2e  row-cond %.2e  full %.2e" %
                  (i, side, (sv > 1e-8*sv[0]).sum(), np.round(sv/sv[0], 3), b['col'], b['row'], b['all']), flush=True)
