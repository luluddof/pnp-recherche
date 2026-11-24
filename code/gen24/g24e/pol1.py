import numpy as np, sys
from lm import *
N,p,q,f = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
z = from_vp(np.load(f), N, p, q)
print("start d2", np.linalg.norm(resid(z,N,p,q,swapmat(N)))**2, "max", np.abs(z).max())
z, h = lm(z, N, p, q, iters=int(sys.argv[5]) if len(sys.argv)>5 else 150, verbose=True)
np.save(f"lm_{N}_{p}_{q}.npy", z)
