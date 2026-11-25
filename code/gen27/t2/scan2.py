# R restarts (it iters), keep best K by tail, polish each in balanced gauge with long L-BFGS; report sv gap and Schmidt rank.
import numpy as np, sys, time
from scipy.optimize import minimize
from vp11 import fg, diag, swap
import gauge
N,p,q,R,seed,it,K,tmax=map(int,sys.argv[1:9]); F=swap(N); rng=np.random.default_rng(seed); t0=time.time(); s=N*N; res=[]
for r in range(R):
    if time.time()-t0>tmax*0.6: break
    z0=rng.standard_normal(2*q*s)*rng.choice([.4,.7,1.])
    z=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x
    d=diag(z,N,p,q); res.append((d['tail'],r,z)); print(f"r{r} tail={d['tail']:.2e} coef={d['coef']:.1f}",flush=True)
res.sort(key=lambda t:t[0])
for tail,r,z in res[:K]:
    if time.time()-t0>tmax: break
    G,H=z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)
    sc=[np.sqrt(np.linalg.norm(h)/max(np.linalg.norm(g),1e-300)) for g,h in zip(G,H)]; G=np.array([g*c for g,c in zip(G,sc)]); H=np.array([h/c for h,c in zip(H,sc)])
    G2,H2=gauge.run(G,H,N,1); zb=np.concatenate([G2.ravel(),H2.ravel()])
    zb=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),zb,jac=True,method='L-BFGS-B',options=dict(maxiter=4*it,ftol=1e-30,gtol=1e-18,maxcor=50)).x
    Km=np.einsum('kxu,kyv->xyuv',zb[:q*s].reshape(q,N,N),zb[q*s:].reshape(q,N,N)).reshape(s,s); sv=np.linalg.svd(F.reshape(s,s)-Km,compute_uv=False)
    sch=np.linalg.svd(np.einsum('kxu,kyv->xuyv',zb[:q*s].reshape(q,N,N),zb[q*s:].reshape(q,N,N)).reshape(s,s),compute_uv=False)
    print(f"POLISH r{r}: tail {tail:.2e} -> sv[p-1]={sv[p-1]:.3e} sv[p]={sv[p]:.3e} gap={sv[p-1]/sv[p]:.2e} coef={abs(zb).max():.2f} Schmidt[q-1]={sch[q-1]:.3e} Schmidt[q]={sch[q]:.1e} t={time.time()-t0:.0f}",flush=True)
    if sv[p]<1e-10: np.save(f"ex_{N}_{p}_{q}_P{seed}r{r}.npy",zb)
