# continue L-BFGS from a saved point (vp objective), report tail and max|coeff| every chunk.  [N]
import numpy as np, sys, time
from scipy.optimize import minimize
from vp import swap, fg
N, p, q = map(int, sys.argv[1].split(',')); z = np.load(sys.argv[2]); tmax = float(sys.argv[3]); F = swap(N); t0 = time.time()
while time.time() - t0 < tmax:
    r = minimize(fg, z, args=(N, p, q, F), jac=True, method='L-BFGS-B', options=dict(maxiter=10000, maxfun=20000, ftol=1e-30, gtol=1e-16))
    z = r.x; print(f"{time.time()-t0:.0f}s tail {r.fun:.3e} max|coeff| {np.abs(z).max():.2f}", flush=True)
    if r.fun < 1e-22: break
np.save(f"pol_{N}_{p}_{q}.npy", z)
print(f"N={N} (p,q)=({p},{q}) polished: tail {r.fun:.3e} max|coeff| {np.abs(z).max():.2f} -> {'FEASIBLE(num)' if r.fun < 1e-20 else 'unresolved'}")
