# warm VP: embed exact json (N0,q0) into N>=N0 with q>=q0 terms; optimise tail at target p (all entries free, or border only).
# usage: extvp.py json N q p seeds it [border]
import sys, json, numpy as np, time
from scipy.optimize import minimize
from vp11 import fg, swap
from lib import load
f,N,q,p,ns,it=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]); border=len(sys.argv)>7 and sys.argv[7]=='border'
N0,A,B=load(f); F=swap(N); s=N*N; best=None; t0=time.time()
for seed in range(ns):
    rng=np.random.default_rng(seed)
    G=np.zeros((q,N,N)); H=np.zeros((q,N,N))
    for l in range(len(A)): G[l,:N0,:N0]=A[l]; H[l,:N0,:N0]=B[l]
    sc=rng.choice([0.3,0.6,1.0])
    mask=np.ones((q,N,N),bool); mask[:len(A),:N0,:N0]=False
    G[mask]=rng.standard_normal(mask.sum())*sc; H[mask]=rng.standard_normal(mask.sum())*sc
    z0=np.concatenate([G.ravel(),H.ravel()]); mz=np.concatenate([mask.ravel(),mask.ravel()])
    if border:
        def obj(w):
            z=z0.copy(); z[mz]=w; f_,g_=fg(z,F,N,p,q,'f2',0,0); return f_,g_[mz]
        w=minimize(obj,z0[mz],jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x; z=z0.copy(); z[mz]=w
    else:
        z0=z0+0.02*rng.standard_normal(z0.size)
        z=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x
    K=np.einsum('kxu,kyv->xyuv',z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)).reshape(s,s); sv=np.linalg.svd(F.reshape(s,s)-K,compute_uv=False)
    tail=float(np.sum(sv[p:]**2)); print(f"seed {seed}: tail {tail:.2e} sv[p-1]={sv[p-1]:.2e} sv[p]={sv[p]:.2e} coef {abs(z).max():.1f} t={time.time()-t0:.0f}s",flush=True)
    if best is None or tail<best[0]: best=(tail,z)
np.save(f"ext_{N}_{p}_{q}.npy",best[1]); print("best",best[0])
