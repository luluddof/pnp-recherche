# continue real L-BFGS from a saved point; print d2 and max coefficient every chunk (border vs exact diagnostic)
import numpy as np, sys
from scipy.optimize import minimize
exec(open('vpreal.py').read().split('N,p,q,R,seed=')[0])
fn, N, p, q, chunks, MI = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
z = np.load(fn); F = swap(N)
for c in range(chunks):
    res = minimize(lambda z: fg(z, F, N, p, q), z, jac=True, method='L-BFGS-B', options=dict(maxiter=MI, ftol=1e-18, gtol=1e-14, maxcor=50))
    z = res.x; print(f"chunk {c}: d2={res.fun:.3e} coef={np.abs(z).max():.2f} nit={res.nit}", flush=True)
    np.save(fn.replace('.npy', '_c.npy'), z)
    if res.fun < 1e-24: break
