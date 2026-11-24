# Matrix-level check of the per-point design: P(0)=I, degree, niceness (support of nice part), kernel condition by ranks.
import numpy as np, itertools, sys, time
from gf import GF
from design import *

def rank_mod3(A):
    A=np.array(A,dtype=np.int64)%3
    if A.size==0: return 0
    rk=0; rows,cols=A.shape
    for c in range(cols):
        piv=np.nonzero(A[rk:,c])[0]
        if len(piv)==0: continue
        p=piv[0]+rk
        if p!=rk: A[[rk,p]]=A[[p,rk]]
        if A[rk,c]==2: A[rk]=(A[rk]*2)%3
        nz=np.nonzero(A[:,c])[0]; nz=nz[nz!=rk]
        A[nz]=(A[nz]-np.outer(A[nz,c],A[rk]))%3
        rk+=1
        if rk==rows: break
    return rk

def run(r,Q,m,shift,do_rank=True):
    F=GF(2); t0=time.time()
    taus=list(range(1,Q+1))  # elements 1..Q of GF(9) (ints 1..8, distinct nonzero)
    des=Design(F,Q,m,shift,taus)
    maxdeg,nbp,acc=des.scalar_checks()
    D=r**(2*Q); r2=r*r
    Sw=np.zeros((r2,r2),dtype=np.int64)
    for a in range(r):
        for d in range(r): Sw[a*r+d,d*r+a]=1
    I=np.eye(r2,dtype=np.int64)
    pi=[((I+Sw)*2)%3, ((I-Sw)*2)%3]     # (I +- Sw)/2, 1/2 = 2 mod 3
    def kron_z(z):
        K=np.ones((1,1),dtype=np.int8)
        for s in range(Q): K=(np.kron(K,pi[z[s]])%3).astype(np.int8)
        return K
    # F_9 element e = re + 3*im
    def accumulate(vals):   # vals: dict z -> element ; returns (Re,Im) int8 matrices mod 3
        Re=np.zeros((D,D),dtype=np.int8); Im=np.zeros((D,D),dtype=np.int8)
        for z,e in vals.items():
            if e==0: continue
            K=kron_z(z)
            if e%3: Re=((Re+(e%3)*K)%3).astype(np.int8)
            if e//3: Im=((Im+(e//3)*K)%3).astype(np.int8)
        return Re,Im
    zs=list(itertools.product((0,1),repeat=Q))
    # P(0) = I
    Re,Im=accumulate({z:F.peval(des.Pz(z),0) for z in zs})
    assert (Re==np.eye(D,dtype=np.int8)).all() and not Im.any(), "P(0)!=I"
    # index helpers: u -> a_u(t), d_u(t)
    idx=np.arange(D)
    def digit(t): return (idx//(r2**(Q-1-t)))%r2
    worst=0; out=[]
    for t in range(Q):
        good={}; bad={}
        for z in zs:
            w=sum(z)-z[t]
            (bad if w in des.tops else good)[z]=F.peval(des.Pz(z),taus[t])
        NRe,NIm=accumulate(good)
        pt=digit(t); a_t=pt//r; d_t=pt%r
        allowed=(d_t[:,None]==a_t[None,:])       # row u, column y: need a_y(t)=d_u(t)
        assert not (NRe[~allowed].any() or NIm[~allowed].any()), ("nice part not nice", t)
        nice_nonzero = int((NRe!=0).sum()+(NIm!=0).sum())
        # bad eigenvectors: rows of C_t^B
        ev=eigvecs(r)
        frac=float(des.bad_fraction(r)); ct=round(frac*D)
        res=dict(t=t, ct=ct, frac=frac, nice_nnz=nice_nonzero)
        if do_rank:
            ERe,EIm=accumulate(bad)
            rows=[]
            for choice in itertools.product(range(len(ev)),repeat=Q):
                w=sum(ev[c][1] for s,c in enumerate(choice) if s!=t)
                if w in des.tops:
                    v=np.ones(1,dtype=np.int64)
                    for c in choice: v=np.kron(v,ev[c][0])
                    rows.append(v)
            C=np.array(rows,dtype=np.int64).reshape(len(rows),D)
            assert len(rows)==ct, (len(rows),ct)
            assert rank_mod3(C)==ct
            sumrk=0; rks=[]
            for g in range(r):
                cols=np.nonzero(a_t==g)[0]
                Cg=C[:,cols]; rk=rank_mod3(Cg); rks.append(rk); sumrk+=rk
                urows=np.nonzero(d_t!=g)[0]
                # full rows of P(tau_t) = N + E restricted to cols (N vanishes there), check in rowspace of Cg
                PR=((NRe.astype(np.int64)+ERe)[np.ix_(urows,cols)])%3
                PI=((NIm.astype(np.int64)+EIm)[np.ix_(urows,cols)])%3
                assert not NRe[np.ix_(urows,cols)].any()
                stack=np.vstack([Cg,PR,PI])
                assert rank_mod3(stack)==rk, ("kernel condition fails", t, g)
                # non-vacuity: is the tested block nonzero, and is kernel nontrivial?
                res.setdefault("nonvac",[]).append((int(PR.any() or PI.any()), len(cols)-rk))
            res["rk_per_group"]=rks; res["sum_rk"]=sumrk
            rankE=rank_mod3(np.vstack([(ERe.astype(np.int64))%3,(EIm.astype(np.int64))%3]))
            res["rankE_F3span"]=rankE
        out.append(res); worst=max(worst,frac)
    print(f"r={r} Q={Q} m={m} shift={shift} blocks={des.blocks} D={D} deg={maxdeg} "
          f"badpairs={nbp} accidentally_good={acc} max c_t/D={worst:.4f} (1/m={1/m:.4f}) time={time.time()-t0:.1f}s")
    for res in out: print("   ",res)
    return worst

if __name__=="__main__":
    cfgs=[tuple(int(x) for x in a.split(",")) for a in sys.argv[1:]]
    for (r,Q,m,dr) in cfgs:
        ws=[run(r,Q,m,sh,bool(dr)) for sh in range(1,m+1)]
        print(f"== r={r} Q={Q} m={m}: best shift max c_t/D = {min(ws):.4f}  <= 1/m ? {min(ws)<=1/m+1e-12}")
