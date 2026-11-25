# ALS search: min rank(M - C'), C' = sum_{P in A, Q in B} kappa(P,Q) g(P,Q)  (general coefficient matrix on a product
# support, ring <= |A||B|).  Agreement-subspace formulation: minimise ||(M - C')V||, V orthonormal n x (n-r).
import numpy as np, sys, json
from core import Mperm, pauli_M
from bellwin import majorana_windows
rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 0)
def paulis_from(k, A, B):
    msk = 2**k - 1; ops = []; labs = []
    for t1 in A:
        a1, c = t1 & msk, t1 >> k
        for t2 in B:
            d, a2 = t2 & msk, t2 >> k
            ops.append(pauli_M(k, a1, a2, c, d).astype(complex)); labs.append((t1, t2))
    return ops, labs
def als(Mm, ops, r, iters=400, kappa0=None, tol=1e-12):
    n = Mm.shape[0]; m = len(ops)
    kap = kappa0.copy() if kappa0 is not None else (rng.standard_normal(m) + 1j*rng.standard_normal(m))/np.sqrt(n)
    Opsarr = np.array(ops)
    hist = []
    for it in range(iters):
        C = np.tensordot(kap, Opsarr, 1); E = Mm - C
        U, s, Vh = np.linalg.svd(E); V = Vh.conj().T[:, r:]   # smallest n-r right singular vectors
        res = np.sum(s[r:]**2); hist.append(res)
        if res < tol: break
        Dm = np.stack([(O @ V).ravel() for O in ops], 1); b = (Mm @ V).ravel()
        kap, *_ = np.linalg.lstsq(Dm, b, rcond=None)
    C = np.tensordot(kap, Opsarr, 1); s = np.linalg.svd(Mm - C, compute_uv=False)
    return kap, s, hist
if __name__ == '__main__':
    k = int(sys.argv[1]); i1, i2 = map(int, sys.argv[2].split(',')); rs = list(map(int, sys.argv[3].split(',')))
    fr, beta = majorana_windows(k); d1, d2 = fr[i1], fr[i2]
    Mm = Mperm(k).astype(complex)
    ops, labs = paulis_from(k, d1['supp'], d2['supp'])
    # product start: kappa = kappa1 (x) kappa2 (normalised Pauli coefficients)
    n = 4**k
    k1 = {t: d1['wh'][t]/(d1['L']*n) for t in d1['supp']}; k2 = {t: d2['wh'][t]/(d2['L']*n) for t in d2['supp']}
    kap0 = np.array([float(k1[a])*float(k2[b])*(-1)**bin((a>>k)&(b>>k)).count('1') for a, b in labs], complex)
    C0 = np.tensordot(kap0, np.array(ops), 1); s0 = np.linalg.svd(Mm - C0, compute_uv=False)
    print(f"k={k} A=q{d1['q']} B=q{d2['q']} |AxB|={len(ops)} composite rank(num)={np.sum(s0>1e-8)} (p1+p2={d1['p']+d2['p']})", flush=True)
    for r in rs:
        best = None
        for trial in range(int(sys.argv[5]) if len(sys.argv) > 5 else 6):
            start = kap0 + (0.3*(rng.standard_normal(len(ops)) + 1j*rng.standard_normal(len(ops)))/n if trial % 2 else 0) if trial < 2 else None
            kap, s, hist = als(Mm, ops, r, kappa0=start)
            tail = np.sum(s[r:]**2)
            if best is None or tail < best: best = tail
            print(f'  r={r} trial={trial} tail={tail:.3e} s[r-1..r+2]={np.round(s[max(r-1,0):r+3],5)}', flush=True)
        print(f' r={r} BEST tail {best:.3e}', flush=True)
