"""G13a: numerical search (reals) for NON-commutative instances with small rank profile.
mode 'P'  : only the blocks P_k=J(tau_k)[x_k=0,y_k=1] are forced to rank<=r_k  (relaxation; the theorem's proof bounds sum rank P_k)
mode 'PQ' : also Q_k=J(tau_k)[x_k=1,y_k=0] rank<=r_k   (then min rank E_k = max(rank P_k, rank Q_k) = r_k is attainable)
Unknowns: factors U_k V_k^T (and U'_k V'_k^T). Entry (x,y) lying in blocks k in D must be interpolable by a polynomial
of degree<=d with value X[x,y] at 0: residual = (I-Pi_D)((v_k - X[x,y])/tau_k)_k, Pi_D = projector on deg<=d-1 Vandermonde at tau_D."""
import numpy as np, sys, itertools
from scipy.optimize import least_squares
def setup(l,d,tau,mode):
    n=1<<l; half=n>>1
    idx=[]  # per direction: lists xs (x_k=0), ys (y_k=1)
    for k in range(l):
        xs=[x for x in range(n) if not (x>>k)&1]; ys=[y for y in range(n) if (y>>k)&1]
        idx.append((xs,ys))
    cons=[]  # (list of (block id, i, j)), Xval, projector N
    for x in range(n):
        for y in range(n):
            D=[k for k in range(l) if ((x^y)>>k)&1]
            if mode=='P': D=[k for k in D if not (x>>k)&1]
            if len(D)<=d: continue
            ent=[]
            for k in D:
                if not (x>>k)&1: ent.append((2*k, idx[k][0].index(x), idx[k][1].index(y)))
                else: ent.append((2*k+1, idx[k][1].index(x), idx[k][0].index(y)))
            t=np.array([tau[k] for k in D],float)
            V=np.vander(t,d,increasing=True) if d>0 else np.zeros((len(D),0))
            Qm=np.linalg.qr(V)[0] if d>0 else np.zeros((len(D),0))
            Nn=np.eye(len(D))-Qm@Qm.T
            # reduce to an orthonormal basis of the complement
            w,vec=np.linalg.eigh(Nn); C=vec[:,w>0.5].T
            cons.append((ent,1.0 if y==(n-1)^x else 0.0,C/t[None,:],C@(1/t)))
    return idx,cons,half
def solve(l,d,tau,ranks,mode,restarts=20,seed=0,verbose=False):
    idx,cons,half=setup(l,d,tau,mode); rng=np.random.default_rng(seed)
    nb=2*l if mode=='PQ' else 2*l
    sizes=[]
    for b in range(2*l):
        r=ranks[b//2] if (mode=='PQ' or b%2==0) else 0
        sizes.append(r)
    offs=np.cumsum([0]+[2*half*r for r in sizes])
    def blocks(z):
        Bs=[]
        for b in range(2*l):
            r=sizes[b]
            if r==0: Bs.append(np.zeros((half,half))); continue
            U=z[offs[b]:offs[b]+half*r].reshape(half,r); V=z[offs[b]+half*r:offs[b+1]].reshape(half,r)
            Bs.append(U@V.T)
        return Bs
    def res(z):
        Bs=blocks(z); out=[]
        for ent,xv,C,c1 in cons:
            v=np.array([Bs[b][i,j] for b,i,j in ent])
            out.append(C@v-xv*c1)
        return np.concatenate(out) if out else np.zeros(1)
    best=(np.inf,None)
    for it in range(restarts):
        z0=rng.normal(size=offs[-1])*rng.choice([0.3,1.0,2.0])
        if offs[-1]==0:
            f=np.linalg.norm(res(z0)); return f,None
        s=least_squares(res,z0,method='trf',xtol=1e-15,ftol=1e-15,gtol=1e-15,max_nfev=400)
        f=np.linalg.norm(s.fun)
        if f<best[0]: best=(f,s.x)
        if f<1e-9: break
    return best[0],(blocks(best[1]) if best[1] is not None else None)
if __name__=="__main__":
    l=int(sys.argv[1]); d=int(sys.argv[2]); mode=sys.argv[3]; total=int(sys.argv[4]); restarts=int(sys.argv[5]) if len(sys.argv)>5 else 15
    tau=[float(k+1) for k in range(l)]; half=1<<(l-1)
    profs=sorted(set(tuple(sorted(p)) for p in itertools.product(range(half+1),repeat=l) if sum(p)==total))
    for pr in profs:
        # all distinct permutations matter in principle (tau differ) but by symmetry of the problem class under relabeling directions+points, sorted suffices
        f,_=solve(l,d,tau,pr,mode,restarts=restarts)
        print(f"l={l} d={d} mode={mode} ranks={pr} total={total}: best residual {f:.3e} -> {'FEASIBLE (numerically)' if f<1e-8 else 'no solution found'}",flush=True)
