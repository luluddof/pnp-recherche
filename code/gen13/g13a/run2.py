"""min feasible totals (numerical, reals) + verification of the theorem's injectivity on the found low-rank instances."""
import ncsearch, numpy as np, sys, itertools, time
from math import comb
def bound(l,d,ma,mc): return sum(comb(l,w)*(ma+1-w) for w in range(mc+1))
def Nm(l,d):
    m=(l-d-1)//2; return bound(l,d,m,m)
def best(l,d):
    return max((bound(l,d,l-1-d-mc,mc),(l-1-d-mc,mc)) for mc in range(l) if l-1-d-mc>=mc)
def lemma(l,d):
    import math
    return max(math.ceil(l*(s-d)*2**(l-s)/s) for s in range(d+1,l+1))
def reconstruct(l,d,tau,Bs,mode):
    n=1<<l; idx=[([x for x in range(n) if not (x>>k)&1],[y for y in range(n) if (y>>k)&1]) for k in range(l)]
    co=np.zeros((d+1,n,n))
    for x in range(n):
        for y in range(n):
            xv=1.0 if y==(n-1)^x else 0.0; co[0,x,y]=xv
            D=[k for k in range(l) if ((x^y)>>k)&1]
            if mode=='P': D=[k for k in D if not (x>>k)&1]
            if not D or d==0: continue
            vals=[]
            for k in D:
                if not (x>>k)&1: vals.append(Bs[2*k][idx[k][0].index(x),idx[k][1].index(y)])
                else: vals.append(Bs[2*k+1][idx[k][1].index(x),idx[k][0].index(y)])
            t=np.array([tau[k] for k in D]); V=np.vander(t,d+1,increasing=True)[:,1:]
            c=np.linalg.lstsq(V,np.array(vals)-xv,rcond=None)[0]; co[1:,x,y]=c
    return co
def verify(l,d,tau,co,tol=1e-7):
    n=1<<l; ranksP=[];ranksQ=[]; Ps=[]
    for k in range(l):
        Jt=sum(co[i]*tau[k]**i for i in range(d+1))
        xs=[x for x in range(n) if not (x>>k)&1]; ys=[y for y in range(n) if (y>>k)&1]
        P=Jt[np.ix_(xs,ys)]; Q=Jt[np.ix_(ys,xs)]; Ps.append((xs,ys,P))
        sc=max(1,np.abs(Jt).max())
        ranksP.append(int((np.linalg.svd(P,compute_uv=False)>tol*sc).sum())); ranksQ.append(int((np.linalg.svd(Q,compute_uv=False)>tol*sc).sum()))
    out=[]
    for mc in range(l):
        ma=l-1-d-mc
        if ma<mc: break
        rows=[]
        for x in range(n):
            w=bin(x).count('1')
            if w>mc: continue
            for i in range(ma-w+1):
                row=[]
                for k in range(l):
                    xs,ys,P=Ps[k]
                    if (x>>k)&1: row+= [0.0]*len(ys)
                    else:
                        v=tau[k]**i*np.prod([1-tau[k]/tau[j] for j in range(l) if (x>>j)&1])
                        row+=list(v*P[xs.index(x)])
                rows.append(row)
        M=np.array(rows); sv=np.linalg.svd(M,compute_uv=False)
        out.append((ma,mc,M.shape[0],int((sv>1e-9*sv[0]).sum()),sv[-1]/sv[0]))
    return ranksP,ranksQ,out
if __name__=="__main__":
    l=int(sys.argv[1]); d=int(sys.argv[2]); mode=sys.argv[3]; lo=int(sys.argv[4]); hi=int(sys.argv[5]); restarts=int(sys.argv[6]); tlim=float(sys.argv[7])
    tau=[float(k+1) for k in range(l)]; half=1<<(l-1); t0=time.time()
    print(f"(l,d)=({l},{d}) mode={mode}: N_m={Nm(l,d)}  asym bound={best(l,d)}  subcube-lemma LB={lemma(l,d)}",flush=True)
    for total in range(lo,hi+1):
        profs=sorted(set(tuple(sorted(p)) for p in itertools.product(range(half+1),repeat=l) if sum(p)==total),key=lambda p:max(p)-min(p))
        found=False
        for pr in profs:
            if time.time()-t0>tlim: print("  time limit"); sys.exit()
            f,Bs=ncsearch.solve(l,d,tau,pr,mode,restarts=restarts,seed=total)
            if f<1e-9 and Bs is not None and max(np.abs(b).max() for b in Bs)<1e4:
                co=reconstruct(l,d,tau,Bs,mode); rP,rQ,inj=verify(l,d,tau,co)
                print(f"  total={total} ranks={pr}: FEASIBLE res={f:.1e}; J rebuilt: rank P_k={rP} (sum {sum(rP)}), rank Q_k={rQ}; injectivity (ma,mc,dimA,rank,relmin sv)={[(a,b,c,e,float(f'{g:.2e}')) for a,b,c,e,g in inj]}",flush=True)
                found=True; break
            else: print(f"  total={total} ranks={pr}: not found (res {f:.1e})",flush=True)
        if found: break
