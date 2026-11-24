# real-valued variant of vp.py (g24g): B = sum_l G_l (x) H_l real; saves real z = [G.ravel(), H.ravel()]
import numpy as np, sys, time
from scipy.optimize import minimize
def swap(N):
    X = np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): X[x,y,y,x] = 1
    return X
def fg(z, F, N, p, q):
    s=N*N; G=z[:q*s].reshape(q,N,N); H=z[q*s:].reshape(q,N,N)
    B=np.einsum('kxu,kyv->xyuv',G,H); M=(F-B).reshape(s,s)
    U,S,Vh=np.linalg.svd(M); E=((U[:,p:]*S[p:])@Vh[p:,:]).reshape(N,N,N,N)
    gG=-2*np.einsum('xyuv,kyv->kxu',E,H); gH=-2*np.einsum('xyuv,kxu->kyv',E,G)
    return float(np.sum(S[p:]**2)), np.concatenate([gG.ravel(),gH.ravel()])
N,p,q,R,seed=map(int,sys.argv[1:6]); MI=int(sys.argv[6]) if len(sys.argv)>6 else 4000; F=swap(N); rng=np.random.default_rng(seed); best=1e9
for r in range(R):
    z0=rng.standard_normal(2*q*N*N)*rng.choice([.5,1.,1.5])
    res=minimize(lambda z: fg(z,F,N,p,q),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=MI,ftol=1e-16,gtol=1e-12,maxcor=30))
    print(f"r{r} d2={res.fun:.3e} coef={np.abs(res.x).max():.2f} nit={res.nit}",flush=True)
    if res.fun<best: best=res.fun; np.save(f"real_{N}_{p}_{q}.npy",res.x)
print("BEST",best)
