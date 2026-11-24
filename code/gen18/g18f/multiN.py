"""N=2 hexagon: P = A + sum_i B_i, B_i triangle operators with bond triples Dlist[i] = (D01, D12, D02); minimise the tail
sum_{k>r} sigma_k(P - sum B_i)^2 (A eliminated).  usage: multi.py "r1,r2,.." "D01.D12.D02;D01.D12.D02;..." starts maxiter"""
import sys, numpy as np, time
from hexlib import *
N = int(sys.argv[5]); d = N*N; P = perm_P(N); TL = float(sys.argv[6])
rs = [int(a) for a in sys.argv[1].split(',')]
Dl = [tuple(int(x) for x in t.split('.')) for t in sys.argv[2].split(';')]
starts = int(sys.argv[3]); maxiter = int(sys.argv[4])
shapes = []
for (D01, D12, D02) in Dl: shapes += [(d, d, D01, D02), (d, d, D01, D12), (d, d, D12, D02)]
sizes = [int(np.prod(s)) for s in shapes]; n = sum(sizes)
def unpack(x):
    o = 0; out = []
    for s, k in zip(shapes, sizes): out.append(x[o:o+k].reshape(s)); o += k
    return out
def make(r):
    def fg(x):
        T = unpack(x); R = P.copy()
        for i in range(len(Dl)): R -= build(*T[3*i:3*i+3])
        f, G, s = tail(R, r); G6 = (-2*G).reshape((d,)*6); g = []
        for i in range(len(Dl)):
            B0, B1, B2 = T[3*i:3*i+3]
            g += [np.einsum('abcijk,bjAB,ckBG->aiAG', G6, B1, B2, optimize=True).ravel(),
                  np.einsum('abcijk,aiAG,ckBG->bjAB', G6, B0, B2, optimize=True).ravel(),
                  np.einsum('abcijk,aiAG,bjAB->ckBG', G6, B0, B1, optimize=True).ravel()]
        return f, np.concatenate(g)
    return fg
t0 = time.time()
cost = max(np.log2(np.prod(D)) for D in Dl)
for r in rs:
    fg = make(r); rng = np.random.default_rng(7 + r); vals = []; best = (9e9, None)
    for s in range(starts):
        x0 = rng.normal(size=n) * [0.3, 0.6, 1.0][s % 3]
        res = minimize(fg, x0, jac=True, method='L-BFGS-B', options=dict(maxiter=maxiter, maxfun=2*maxiter, ftol=1e-16, gtol=1e-11))
        vals.append(res.fun)
        if res.fun < best[0]: best = (res.fun, res.x.copy())
        print(f"   start {s}: f={res.fun:.4e} nit={res.nit} [{time.time()-t0:.0f}s]", flush=True)
        if res.fun < 1e-10 or time.time()-t0 > TL: break
    vals = np.sort(vals)
    print(f"r={r} terms={Dl} (max term cost {cost:.2f} bits; log2 r + cost = {np.log2(max(r,1))+cost:.2f} of {6*np.log2(N):.2f}): best residual^2={best[0]:.3e}; "
          f"lowest minima {np.round(vals[:5], 4)} starts={len(vals)} [{time.time()-t0:.0f}s]", flush=True)
    if best[0] < 1e-8: np.save(f"multi_r{r}_{sys.argv[2].replace(';','_')}.npy", best[1])
