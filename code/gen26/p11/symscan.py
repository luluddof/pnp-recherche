# descending-p scan with the symmetric-Theta ansatz (see symR.py).  usage: symscan.py N q p0 R seed tmax it
import numpy as np, sys, time
from scipy.optimize import minimize
N,q,p0,Rr,seed,tmax,it=map(int,sys.argv[1:8]); n=N*N
S=np.zeros((n,n))
for x in range(N):
    for y in range(N): S[x*N+y,y*N+x]=1
re=lambda M: M.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n,n)
def mk(p):
    def fg(w):
        U=w[:n*p].reshape(n,p); lam=w[n*p:]
        R=(U*lam)@U.T@S; Uu,sg,Vt=np.linalg.svd(re(S-R))
        f=float(np.sum(sg[q:]**2)); Z=re((Uu[:,q:]*sg[q:])@Vt[q:])@S
        return f,np.concatenate([(-2*(Z+Z.T)@U*lam).ravel(),-2*np.einsum('ip,ij,jp->p',U,Z,U)])
    return fg
rng=np.random.default_rng(seed); t0=time.time(); p=p0
while p>0 and time.time()-t0<tmax:
    fg=mk(p); found=False; bt=1e9
    for r in range(Rr):
        if time.time()-t0>tmax: break
        w0=np.concatenate([rng.standard_normal(n*p)/np.sqrt(n)*rng.choice([.5,1,2]),rng.choice([1.,2.,3.],p)])
        res=minimize(fg,w0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30))
        U=res.x[:n*p].reshape(n,p); lam=res.x[n*p:]; W=(U*lam)@U.T; ev=np.linalg.eigvalsh(W); ev=ev[abs(ev)>1e-6]
        sg=np.linalg.svd(re(S-W@S),compute_uv=False); gap=sg[q-1]/max(sg[q],1e-300); bt=min(bt,res.fun)
        if res.fun<1e-10 and gap>1e5:
            print(f"N={N} q={q} p={p} r{r} FOUND tail={res.fun:.1e} osr-gap={gap:.1e} eig(1-Theta)={np.round(np.sort(ev),3).tolist()} t={time.time()-t0:.0f}",flush=True)
            np.save(f"sym_{N}_{p}_{q}_s{seed}.npy",res.x); found=True; break
    if not found: print(f"N={N} q={q} p={p}: not found ({r+1} restarts) best tail {bt:.2e} t={time.time()-t0:.0f}",flush=True); break
    p-=1
