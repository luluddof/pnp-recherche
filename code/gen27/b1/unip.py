# KT5: unipotent Theta (1-Theta nilpotent of index 2): D = S - K, rank D <= p and D S D = 0.  VP: tail_p(D) + lam ||D S D||_F^2
import numpy as np, sys, time
from scipy.optimize import minimize
def swapm(N):
    S=np.zeros((N,)*4)
    for x in range(N):
        for y in range(N): S[x,y,y,x]=1
    return S.reshape(N*N,N*N)
def fg(z,Sm,N,p,q,lam):
    s=N*N; G=z[:q*s].reshape(q,N,N); H=z[q*s:].reshape(q,N,N)
    K=np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s); D=Sm-K
    U,Sv,Vh=np.linalg.svd(D); f=float(np.sum(Sv[p:]**2)); E=(U[:,p:]*(2*Sv[p:]))@Vh[p:,:]
    Z=D@Sm@D; f+=lam*float(np.sum(Z**2)); E+=lam*2*(Z@D.T@Sm+Sm@D.T@Z)
    E=E.reshape(N,N,N,N)   # df/dD ; dK = -dD
    gG=-np.einsum('xyuv,kyv->kxu',E,H); gH=-np.einsum('xyuv,kxu->kyv',E,G)
    return f,np.concatenate([gG.ravel(),gH.ravel()])
N,p,q=map(int,sys.argv[1:4]); lam=float(sys.argv[4]) if len(sys.argv)>4 else 1.0; seeds=int(sys.argv[5]) if len(sys.argv)>5 else 6
Sm=swapm(N); s=N*N
for sd in range(seeds):
    z=np.random.default_rng(sd).standard_normal(2*q*s)/np.sqrt(N); t=time.time()
    res=minimize(fg,z,args=(Sm,N,p,q,lam),jac=True,method='L-BFGS-B',options=dict(maxiter=6000,maxcor=50,gtol=1e-14,ftol=1e-18))
    G=res.x[:q*s].reshape(q,N,N); H=res.x[q*s:].reshape(q,N,N); D=Sm-np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s)
    sv=np.linalg.svd(D,compute_uv=False); nil=np.linalg.norm(D@Sm@D); ev=np.linalg.eigvals(np.eye(s)-D@Sm)
    print(f"N={N} p={p} q={q} lam={lam} seed={sd}: f={res.fun:.2e} sv[p-1],sv[p]={sv[p-1]:.2e},{sv[p]:.2e} ||DSD||={nil:.2e} max|ev-1|={np.max(np.abs(ev-1)):.2e} t={time.time()-t:.0f}s",flush=True)
