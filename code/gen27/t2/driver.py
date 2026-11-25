# G27-T2 driver: for each (N,q): warm starts (border extension of best exact at N-1; best exact at (N, q-1) + one term),
# descending VP, then exactification; exact points saved in exact/.  usage: driver.py "N:q1-q2,..." R it
import os
os.environ['OMP_NUM_THREADS']='1'; os.environ['OPENBLAS_NUM_THREADS']='1'
import sys, glob, json, time, numpy as np, subprocess
from multiprocessing import Pool
from scipy.optimize import minimize
from vp11 import fg, swap
import gauge
from lib import load, SRC, exrank, modrank, swapm, Kmat, osr, save, frank
R=int(sys.argv[2]); IT=int(sys.argv[3]); NP=int(os.environ.get('NPROC','8'))
def db():
    pts=[]
    for f in SRC+glob.glob('cons/c_*.json')+glob.glob('exact/exact_*.json'):
        try: N,E,C=load(f)
        except Exception: continue
        b=os.path.basename(f).replace('.json','').split('_'); p,q=int(b[-2]),int(b[-1])
        if f in SRC or 'Q8' in f: pass
        pts.append((N,p,q,f))
    return pts
def best(N,q):
    c=[t for t in db() if t[0]==N and t[2]<=q]
    return min(c,key=lambda t:(t[1],t[2])) if c else None
F={}
def Fm(N):
    if N not in F: F[N]=swap(N)
    return F[N]
def svals(z,N,q):
    s=N*N; K=np.einsum('kxu,kyv->xyuv',z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)).reshape(s,s); return np.linalg.svd(Fm(N).reshape(s,s)-K,compute_uv=False)
def opt(z0,N,p,q,it,mask=None):
    if mask is None:
        return minimize(lambda z: fg(z,Fm(N),N,p,q,'f2',0,0),z0,jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x
    def obj(w):
        z=z0.copy(); z[mask]=w; f_,g_=fg(z,Fm(N),N,p,q,'f2',0,0); return f_,g_[mask]
    w=minimize(obj,z0[mask],jac=True,method='L-BFGS-B',options=dict(maxiter=it,ftol=1e-17,gtol=1e-14,maxcor=30)).x; z=z0.copy(); z[mask]=w; return z
def polish(z,N,p,q,it):
    s=N*N; G,H=z[:q*s].reshape(q,N,N),z[q*s:].reshape(q,N,N)
    sc=[np.sqrt(np.linalg.norm(h)/max(np.linalg.norm(g),1e-300)) if np.linalg.norm(g)>1e-12 else 1 for g,h in zip(G,H)]
    G=np.array([g*c for g,c in zip(G,sc)]); H=np.array([h/c for h,c in zip(H,sc)])
    try: G,H=gauge.run(G,H,N,1)
    except Exception: pass
    zb=np.concatenate([G.ravel(),H.ravel()])
    zb=minimize(lambda z: fg(z,Fm(N),N,p,q,'f2',0,0),zb,jac=True,method='L-BFGS-B',options=dict(maxiter=3*it,ftol=1e-30,gtol=1e-18,maxcor=50)).x
    return zb,svals(zb,N,q)
def embed(f,N,q,rng,sc):
    N0,A,B=load(f); s=N*N; G=np.zeros((q,N,N)); H=np.zeros((q,N,N))
    if len(A)>q:   # file has more terms than its osr: re-decompose by SVD of the realignment (float warm start)
        Rm=sum(np.outer(a.ravel(),b.ravel()) for a,b in zip(A,B)).astype(float); U,S_,Vh=np.linalg.svd(Rm)
        A=[np.sqrt(S_[l])*U[:,l].reshape(N0,N0) for l in range(q)]; B=[np.sqrt(S_[l])*Vh[l].reshape(N0,N0) for l in range(q)]
    k=min(len(A),q)
    for l in range(k): G[l,:N0,:N0]=A[l]; H[l,:N0,:N0]=B[l]
    mask=np.ones((q,N,N),bool); mask[:k,:N0,:N0]=False
    G[mask]=rng.standard_normal(mask.sum())*sc; H[mask]=rng.standard_normal(mask.sum())*sc
    return np.concatenate([G.ravel(),H.ravel()]), np.concatenate([mask.ravel(),mask.ravel()])
def job_ext(a):
    f,N,q,p,seed,it=a; rng=np.random.default_rng(seed); z0,m=embed(f,N,q,rng,rng.choice([0.3,0.6,1.0]))
    z=opt(z0,N,p,q,it,m); sv=svals(z,N,q); return float(np.sum(sv[p:]**2)),z
def job_desc(a):
    z0,N,q,p,seed,it,noise=a; rng=np.random.default_rng(seed)
    z=opt(z0+noise*rng.standard_normal(z0.size)*(seed>0),N,p,q,it); sv=svals(z,N,q); tail=float(np.sum(sv[p:]**2))
    if tail<1e-8:
        r=min(int(np.sum(sv>1e-6)),p); zb,sv2=polish(z,N,r,q,it)
        if sv2[r]<1e-11 and sv2[r-1]/sv2[r]>1e9: return (tail,seed,zb,True,r)
    return (tail,seed,z,False,p)
def exactify(zf,N,p,q,R=4):
    out=subprocess.run(['python3','exactify.py',zf,str(N),str(p),str(q),str(R),'0'],capture_output=True,text=True,timeout=900).stdout
    for line in out.splitlines():
        if line.startswith('SAVED'):
            fn=line.split()[1]; os.replace(fn,'exact/'+fn); return 'exact/'+fn, out
    return None,out
def descend(pool,z,N,q,p,tag,tmax=int(os.environ.get('TMAX','300'))):
    t0=time.time(); acc=[]
    while time.time()-t0<tmax:
        res=pool.map(job_desc,[(z,N,q,p-1,sd,IT,0.05) for sd in range(R)])
        res.sort(key=lambda r:(not r[3],r[4],r[0]))
        if not res[0][3]:
            print(f"   [{tag}] fail at p={p-1}: best tail {res[0][0]:.2e} t={time.time()-t0:.0f}",flush=True); break
        p=res[0][4]; z=res[0][2]; zf=f"work/d_{N}_{p}_{q}_{tag}.npy"; np.save(zf,z); acc.append((p,zf)); print(f"   [{tag}] accept ({p},{q})@{N} t={time.time()-t0:.0f}",flush=True)
    return acc
def task(pool,N,q):
    t0=time.time(); acc=[]
    b=best(N-1,q)
    if b:
        N0,p0,q0,f=b
        for pt in [p0+N0,p0+N0+3,p0+N0+6,p0+2*N0+1]:
            res=pool.map(job_ext,[(f,N,q,pt,sd,IT) for sd in range(min(R,4))]); res.sort(key=lambda r:r[0])
            print(f" ({N},{q}) ext from {os.path.basename(f)} ({p0},{q0})@{N0} at p={pt}: tail {res[0][0]:.2e}",flush=True)
            if res[0][0]<1e-10: acc+=descend(pool,res[0][1],N,q,pt+1,'ext'); break
    b2=best(N,q-1)
    if b2 and b2[0]==N and os.environ.get('ADDQ','1')=='1' and b2[1]<N*(N-1)//2-1:
        N0,p0,q0,f=b2; rng=np.random.default_rng(7); z0,m=embed(f,N,q,rng,0.05)
        print(f" ({N},{q}) from {os.path.basename(f)} ({p0},{q0}) + term",flush=True)
        acc+=descend(pool,z0,N,q,p0,'addq')
    cur=best(N,q); curp=cur[1] if cur else 10**9
    acc.sort()
    for p,zf in acc:
        if p>=curp: break
        fn,out=exactify(zf,N,p,q)
        print(f" exactify ({p},{q})@{N}: {fn}  {out.strip().splitlines()[-1] if out.strip() else ''}",flush=True)
        if fn: break
    print(f"DONE ({N},{q}) best now {best(N,q)[1:3]} t={time.time()-t0:.0f}",flush=True)
if __name__=='__main__':
    todo=[]
    for item in sys.argv[1].split(','):
        N,qs=item.split(':'); a,b=map(int,qs.split('-')); todo+=[(int(N),q) for q in range(a,b+1)]
    with Pool(NP) as pool:
        for N,q in todo: task(pool,N,q)
