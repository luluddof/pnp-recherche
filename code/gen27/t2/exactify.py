# exactify a numerical (N,p,q) solution: sparsify (L1 over GL_N x GL_N x GL_q gauge) -> diagonal log-fit -> round -> exact rank (flint)
# usage: exactify.py z.npy N p q restarts seed
import os
os.environ.setdefault('OMP_NUM_THREADS','2')
import sys, numpy as np, subprocess, json
from lib import exrank, modrank, swapm, Kmat, osr, save
f,N,p,q,R,seed=sys.argv[1],*map(int,sys.argv[2:7]); s=N*N
from scipy.optimize import minimize
z=np.load(f); G=z[:q*s].reshape(q,N,N); Gp=z[q*s:].reshape(q,N,N)
def tr(w):
    g=w[:N*N].reshape(N,N); h=w[N*N:2*N*N].reshape(N,N); m=w[2*N*N:].reshape(q,q)
    gi,hi,mit=np.linalg.inv(g),np.linalg.inv(h),np.linalg.inv(m).T
    A=np.einsum('lk,kab->lab',m,np.einsum('xa,kab,bu->kxu',g,G,hi)); B=np.einsum('lk,kab->lab',mit,np.einsum('ya,kab,bv->kyv',h,Gp,gi)); return A,B
def obj(w,d):
    try: A,B=tr(w)
    except np.linalg.LinAlgError: return 1e9
    return np.sqrt(A**2+d).sum()+np.sqrt(B**2+d).sum()
def diagfit(A,B,thr):
    G=A.copy(); H=B.copy(); G[abs(G)<thr]=0; H[abs(H)<thr]=0
    rows=[];rhs=[]; nv=2*N+q
    for l in range(q):
        for a,b in zip(*np.nonzero(G[l])):
            r=np.zeros(nv); r[2*N+l]=1; r[a]+=1; r[N+b]-=1; rows.append(r); rhs.append(-np.log(abs(G[l,a,b])))
        for a,b in zip(*np.nonzero(H[l])):
            r=np.zeros(nv); r[2*N+l]=-1; r[N+a]+=1; r[b]-=1; rows.append(r); rhs.append(-np.log(abs(H[l,a,b])))
    if not rows: return None
    Am=np.array(rows); w=np.linalg.lstsq(Am,np.array(rhs),rcond=None)[0]; xi,eta,tau=w[:N],w[N:2*N],w[2*N:]
    G2=np.array([np.exp(tau[l])*np.diag(np.exp(xi))@G[l]@np.diag(np.exp(-eta)) for l in range(q)])
    H2=np.array([np.exp(-tau[l])*np.diag(np.exp(eta))@H[l]@np.diag(np.exp(-xi)) for l in range(q)])
    return G2,H2,float(np.abs(Am@w-np.array(rhs)).max())
def tryround(A,B):
    best=None
    for thr in [0.05,0.1,0.2,0.3,0.4,0.5]:
        r=diagfit(A,B,thr)
        if r is None: continue
        G2,H2,res=r
        Gi=np.rint(G2).astype(int); Hi=np.rint(H2).astype(int)
        E=list(Gi); C=list(Hi); M=swapm(N)-Kmat(E,C)
        pr=int(np.linalg.matrix_rank(M.astype(float)))
        if best is None or pr<best[0]: best=(pr,thr,res,E,C)
    return best
rng=np.random.default_rng(seed); bestL=(1e18,None); found=None
for r in range(R):
    w=np.concatenate([np.eye(N).ravel(),np.eye(N).ravel(),np.eye(q).ravel()])+(0.5 if r else 0.0)*rng.standard_normal(2*N*N+q*q)
    for d in (1e-2,1e-4,1e-6,1e-8):
        w=minimize(obj,w,args=(d,),method='L-BFGS-B',options=dict(maxiter=3000)).x
    v=obj(w,1e-10); A,B=tr(w)
    tb=tryround(A,B); print(f"restart {r}: L1={v:.2f} round -> p={tb[0] if tb else None} (thr {tb[1] if tb else None}, logfit res {tb[2] if tb else 0:.1e})",flush=True)
    if tb and tb[0]<=p:
        E,C=tb[3],tb[4]; M=swapm(N)-Kmat(E,C); rk=exrank(M); o=exrank(osr(E,C))
        mods=[(modrank(M,pr),modrank(osr(E,C),pr)) for pr in (1000003,998244353)]
        print(f"EXACT CHECK: rank_Q={rk} osr={o} mod-primes={mods} entries={sorted(set(np.concatenate([np.ravel(x) for x in E+C]).tolist()))}",flush=True)
        if rk<=p and o<=q:
            fn=f"exact_{N}_{rk}_{o}.json"; save(fn,N,E,C,{"rank":rk,"osr":o}); print("SAVED",fn); found=fn; break
