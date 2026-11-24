"""vectorised version of ncsearch (P-mode only = the relaxation the theorem's proof actually bounds; plus optional PQ)."""
import numpy as np, sys, itertools, time, scipy.sparse as sp
from scipy.optimize import least_squares
from math import comb, ceil
def build(l,d,tau,mode):
    n=1<<l; h=n>>1
    pos=[({x:i for i,x in enumerate([x for x in range(n) if not (x>>k)&1])},{y:i for i,y in enumerate([y for y in range(n) if (y>>k)&1])}) for k in range(l)]
    rows=[];cols=[];vals=[];b=[]; nr=0
    proj={}
    for x in range(n):
        for y in range(n):
            D=[k for k in range(l) if ((x^y)>>k)&1]
            if mode=='P': D=[k for k in D if not (x>>k)&1]
            if len(D)<=d: continue
            key=tuple(D)
            if key not in proj:
                t=np.array([tau[k] for k in D],float)
                V=np.vander(t,d,increasing=True) if d>0 else np.zeros((len(D),0))
                Qm=np.linalg.qr(V)[0] if d>0 else np.zeros((len(D),0))
                w,vec=np.linalg.eigh(np.eye(len(D))-Qm@Qm.T); C=vec[:,w>0.5].T
                proj[key]=(C/t[None,:],C@(1/t))
            C,c1=proj[key]; xv=1.0 if y==(n-1)^x else 0.0
            cid=[]
            for k in D:
                if not (x>>k)&1: cid.append((2*k)*h*h+pos[k][0][x]*h+pos[k][1][y])
                else: cid.append((2*k+1)*h*h+pos[k][1][x]*h+pos[k][0][y])
            for i in range(C.shape[0]):
                for j,cc in enumerate(cid): rows.append(nr+i);cols.append(cc);vals.append(C[i,j])
                b.append(xv*c1[i])
            nr+=C.shape[0]
    A=sp.csr_matrix((vals,(rows,cols)),shape=(nr,2*l*h*h)); return A,np.array(b),h
def solve(l,d,tau,ranks,mode,restarts=5,seed=0,maxnfev=300,built=None):
    A,b,h=built if built else build(l,d,tau,mode); rng=np.random.default_rng(seed)
    sizes=[(ranks[bb//2] if (mode=='PQ' or bb%2==0) else 0) for bb in range(2*l)]
    offs=np.cumsum([0]+[2*h*r for r in sizes]); Ac=A.tocsc(); Abs={bb:Ac[:,bb*h*h:(bb+1)*h*h].toarray().reshape(-1,h,h) for bb in range(2*l) if sizes[bb]}
    def UV(z,bb):
        r=sizes[bb]; return z[offs[bb]:offs[bb]+h*r].reshape(h,r), z[offs[bb]+h*r:offs[bb+1]].reshape(h,r)
    def vecB(z):
        out=np.zeros(2*l*h*h)
        for bb in range(2*l):
            if sizes[bb]:
                U,V=UV(z,bb); out[bb*h*h:(bb+1)*h*h]=(U@V.T).ravel()
        return out
    def res(z): return A@vecB(z)-b
    def jac(z):
        Jm=np.zeros((A.shape[0],offs[-1]))
        for bb in range(2*l):
            r=sizes[bb]
            if not r: continue
            U,V=UV(z,bb); Ab=Abs[bb]
            Jm[:,offs[bb]:offs[bb]+h*r]=np.einsum('rij,js->ris',Ab,V).reshape(-1,h*r)
            Jm[:,offs[bb]+h*r:offs[bb+1]]=np.einsum('rij,is->rjs',Ab,U).reshape(-1,h*r)
        return Jm
    best=(np.inf,None)
    for it in range(restarts):
        z0=rng.normal(size=offs[-1])*rng.choice([0.3,1.0])
        s=least_squares(res,z0,jac=jac,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=maxnfev)
        f=np.linalg.norm(s.fun)
        if f<best[0]: best=(f,vecB(s.x).reshape(2*l,h,h))
        if f<1e-9: break
    return best
def bound(l,d,ma,mc): return sum(comb(l,w)*(ma+1-w) for w in range(mc+1))
def bestb(l,d): return max(bound(l,d,l-1-d-mc,mc) for mc in range(l) if l-1-d-mc>=mc)
def lemma(l,d): return max(ceil(l*(s-d)*2**(l-s)/s) for s in range(d+1,l+1))
if __name__=="__main__":
    l=int(sys.argv[1]); d=int(sys.argv[2]); mode=sys.argv[3]; tlim=float(sys.argv[5]); restarts=int(sys.argv[4])
    profiles=[tuple(int(c) for c in a.split(',')) for a in sys.argv[6:]]
    tau=[float(k+1) for k in range(l)]; t0=time.time(); built=build(l,d,tau,mode)
    m=(l-d-1)//2
    print(f"(l,d)=({l},{d}) mode={mode}: N_m={bound(l,d,m,m)} asym={bestb(l,d)} subcube-lemma={lemma(l,d)}; #residuals={built[0].shape[0]}",flush=True)
    for pr in profiles:
        if time.time()-t0>tlim: print(" time limit"); break
        f,B=solve(l,d,tau,pr,mode,restarts=restarts,built=built,seed=sum(pr))
        mx=np.abs(B).max() if B is not None else 0
        rk=[int((np.linalg.svd(B[2*k],compute_uv=False)>1e-8*max(1,mx)).sum()) for k in range(l)] if B is not None else None
        print(f"  ranks={pr} sum={sum(pr)}: residual {f:.2e} max|entry|={mx:.1e} numeric rank P_k={rk} -> {'FEASIBLE' if f<1e-9 and mx<1e5 else 'not found'}  [{time.time()-t0:.0f}s]",flush=True)
