import sys, numpy as np, itertools, time
from hexlib import *
N = 2; d = 4; P = perm_P(N)
# gradient check
fg, unpack, n = make_obj(P, d, (2, 2, 2), 5)
x = np.random.default_rng(0).normal(size=n); f, g = fg(x); e = np.random.default_rng(1).normal(size=n)*1e-6
print("gradcheck", (fg(x+e)[0]-fg(x-e)[0])/2, g@e)
Ds = [(1,1,1),(2,1,1),(2,2,1),(2,2,2),(3,1,1),(4,1,1),(3,2,1),(4,2,1),(3,3,1),(4,4,1),(3,2,2),(4,2,2),(3,3,2),(3,3,3),(4,3,3),(4,4,2),(4,4,3)]
t0 = time.time()
for D in Ds:
    lim = int(np.ceil(64/np.prod(D)))        # gain needs r*prod(D) < 64
    out = []
    rmin = None
    for r in range(lim+2, 0, -1) if lim <= 12 else list(range(min(lim+2, 40), 0, -3)):
        f, x, vals = search(P, d, D, r, starts=8, maxiter=400, seed=r)
        out.append((r, f))
        if f < 1e-7: rmin = r
        elif f > 0.5 and rmin is not None: break
    print(f"D={D} prodD={np.prod(D)} gain needs r<{64/np.prod(D):.1f}; min feasible r found={rmin}; (r,res)=" +
          " ".join(f"({r},{f:.1e})" for r, f in out), f"[{time.time()-t0:.0f}s]", flush=True)
