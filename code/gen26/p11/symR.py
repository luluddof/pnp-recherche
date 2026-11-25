# Structured search suggested by the spectra of all exact optima (1-Theta = R S symmetric, PSD, spectrum {0,2,3}):
# ansatz R = U diag(lam) U^T S (so Theta = K S = 1 - U diag(lam) U^T is SYMMETRIC), U: n x p.  Minimise the operator-Schmidt
# tail sum_{i>q} sigma_i(realign(S - R))^2 over U, lam.  usage: symR.py N p q R seed it
import numpy as np, sys, time
from scipy.optimize import minimize
N,p,q,Rr,seed,it=map(int,sys.argv[1:7]); n=N*N
S=np.zeros((n,n))
for x in range(N):
    for y in range(N): S[x*N+y,y*N+x]=1
re=lambda M: M.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n,n)   # (x,y),(u,v) -> (x,u),(y,v); self-inverse
def fg(w):
    U=w[:n*p].reshape(n,p); lam=w[n*p:]
    R=(U*lam)@U.T@S; M=re(S-R); Uu,sg,Vt=np.linalg.svd(M)
    f=float(np.sum(sg[q:]**2)); Et=re((Uu[:,q:]*sg[q:])@Vt[q:]); Z=Et@S
    gU=-2*(Z+Z.T)@U*lam; gl=-2*np.einsum('ip,ij,jp->p',U,Z,U)
    return f,np.concatenate([gU.ravel(),gl])
rng=np.random.default_rng(seed); best=None; t0=time.time()
for r in range(Rr):
    w0=np.concatenate([rng.standard_normal(n*p)/np.sqrt(n)*rng.choice([.5,1,2]),rng.choice([1.,2.,3.],p)])
    res=minimize(fg,w0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30))
    U=res.x[:n*p].reshape(n,p); lam=res.x[n*p:]; R=(U*lam)@U.T@S
    sg=np.linalg.svd(re(S-R),compute_uv=False); W=(U*lam)@U.T; ev=np.linalg.eigvalsh(W)
    print(f"r{r} tail={res.fun:.2e} osr-gap={sg[q-1]/max(sg[q],1e-300):.1e} rank R={np.linalg.matrix_rank(R,1e-8)} eig(1-Theta) nonzero={np.round(ev[abs(ev)>1e-6],3)} t={time.time()-t0:.0f}",flush=True)
    if best is None or res.fun<best[0]: best=(res.fun,res.x)
np.save(f"symR_{N}_{p}_{q}_s{seed}.npy",best[1]); print("BEST",best[0])
