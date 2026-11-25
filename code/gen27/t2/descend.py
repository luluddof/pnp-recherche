# descending-p VP from a warm start z (npy) at (N,q): for p = pstart-1, ... try R noisy restarts (parallel), accept if
# gauge-balanced polish gives sv[p] < 1e-11 and gap > 1e10.  Saves dz_{N}_{p}_{q}.npy for each accepted p.
# usage: descend.py z.npy N q pstart pmin R it tmax [noise]
import os
os.environ.setdefault('OMP_NUM_THREADS','1'); os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys, numpy as np, time
from multiprocessing import Pool
from scipy.optimize import minimize
from vp11 import fg, swap
import gauge
f,N,q,ps,pmin,R,it,tmax=sys.argv[1],*map(int,sys.argv[2:9]); noise=float(sys.argv[9]) if len(sys.argv)>9 else 0.05
F=swap(N); s=N*N; t0=time.time()
def polish(z,p):
    G,H=z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)
    sc=[np.sqrt(np.linalg.norm(h)/max(np.linalg.norm(g),1e-300)) for g,h in zip(G,H)]; G=np.array([g*c for g,c in zip(G,sc)]); H=np.array([h/c for h,c in zip(H,sc)])
    try: G2,H2=gauge.run(G,H,N,1)
    except Exception: G2,H2=G,H
    zb=np.concatenate([G2.ravel(),H2.ravel()])
    zb=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),zb,jac=True,method='L-BFGS-B',options=dict(maxiter=3*it,ftol=1e-30,gtol=1e-18,maxcor=50)).x
    Km=np.einsum('kxu,kyv->xyuv',zb[:q*s].reshape(q,N,N),zb[q*s:].reshape(q,N,N)).reshape(s,s); sv=np.linalg.svd(F.reshape(s,s)-Km,compute_uv=False)
    return zb,sv
def job(args):
    z0,p,seed=args; rng=np.random.default_rng(seed)
    z0=z0+noise*rng.standard_normal(z0.size)*(seed>0)
    z=minimize(lambda z: fg(z,F,N,p,q,'f2',0,0),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x
    Km=np.einsum('kxu,kyv->xyuv',z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)).reshape(s,s); sv=np.linalg.svd(F.reshape(s,s)-Km,compute_uv=False)
    tail=float(np.sum(sv[p:]**2))
    if tail<1e-8:
        r=min(int(np.sum(sv>1e-6)),p)
        zb,sv2=polish(z,r)
        if sv2[r]<1e-11 and sv2[r-1]/sv2[r]>1e9: return (tail,seed,zb,sv2[r-1],sv2[r],True,r)
        return (tail,seed,z,sv2[r-1],sv2[r],False,r)
    return (tail,seed,z,sv[p-1],sv[p],False,p)
z=np.load(f); p=ps
with Pool(int(os.environ.get('NPROC','6'))) as pool:
    while p>pmin and time.time()-t0<tmax:
        p-=1
        res=pool.map(job,[(z,p,sd) for sd in range(R)])
        res.sort(key=lambda r:(not r[5],r[6],r[0]))
        for r in res: print(f"  p={p} r={r[6]} seed{r[1]} tail={r[0]:.2e} sv[p-1]={r[3]:.2e} sv[p]={r[4]:.2e} exact-like={r[5]}",flush=True)
        if res[0][5]:
            p=res[0][6]; z=res[0][2]; np.save(f"dz_{N}_{p}_{q}.npy",z); print(f"ACCEPT ({p},{q})@{N} t={time.time()-t0:.0f}s",flush=True)
        else:
            print(f"FAIL at p={p} t={time.time()-t0:.0f}s",flush=True); break
