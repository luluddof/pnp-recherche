# min-norm Gauss-Newton (lstsq steps, backtracking), mu = 0.  Quadratic convergence <=> regular exact solution nearby.
import numpy as np, sys
from lm import *
N,p,q,f,its = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4],int(sys.argv[5])
rc = float(sys.argv[6]) if len(sys.argv)>6 else 1e-8
z = np.load(f); X = swapmat(N).astype(complex)
d2 = lambda z: float(np.linalg.norm(resid(z,N,p,q,X))**2)
f0 = d2(z)
for it in range(its):
    r = resid(z,N,p,q,X).ravel(); J = jac(z,N,p,q)
    dz = np.linalg.lstsq(J, -r, rcond=rc)[0]; t = 1.0
    while d2(z+t*dz) >= f0 and t > 1e-6: t /= 2
    if t <= 1e-6: print("no descent"); break
    z = z + t*dz; f0 = d2(z)
    print(f"it{it} d2={f0:.3e} t={t:.3g} |z|^2={np.vdot(z,z).real:.1f} max|z|={np.abs(z).max():.2f}", flush=True)
    if f0 < 1e-26: break
np.save("gn_"+f, z)
