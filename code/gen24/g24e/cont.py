# Tikhonov continuation: min d2 + mu ||z||^2, mu decreasing. exact => ||z|| bounded; border => ||z|| -> inf, d2 ~ mu^alpha
import numpy as np, sys
from lm import *
N,p,q,f = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
z = np.load(f)
if not f.startswith('lm') and not f.startswith('ct'): z = from_vp(z, N, p, q)
X = swapmat(N)
for e in np.arange(float(sys.argv[5]) if len(sys.argv)>5 else 1, 12.1, 0.5):
    mu = 10**(-e)
    z, h = lm(z, N, p, q, mu=mu, iters=int(sys.argv[6]) if len(sys.argv)>6 else 60)
    d2 = h[-1][1]
    U,V,G,H = unpack(z,N,p,q)
    print(f"mu=1e-{e:.1f} d2={d2:.3e} |z|^2={np.vdot(z,z).real:.2f} max|z|={np.abs(z).max():.2f} |G||H|={np.linalg.norm(G)*np.linalg.norm(H):.2f} |U||V|={np.linalg.norm(U)*np.linalg.norm(V):.2f} its={len(h)}", flush=True)
    np.save(f"ct_{N}_{p}_{q}.npy", z)
