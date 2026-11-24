#!/usr/bin/env python3
"""G18b exp_aux: per-point piecewise construction over F_3 (coefficients in F_9), auxiliary level.
Checks, for r=2 and small Q:  P(0)=I, deg<=m-1, Def-3.1 kernel condition at every node t with a per-node C_t,
ranks of C_t restricted to groups, and counts reachable measurement vectors (over Z and mod 3)."""
import itertools, sys, math
import numpy as np
rng = np.random.default_rng(1)
# ---------- F_9 = F_3[i]/(i^2+1), element = a+3b -> pair (a,b)
def f9(a,b): return (a%3, b%3)
def add(x,y): return ((x[0]+y[0])%3,(x[1]+y[1])%3)
def mul(x,y): return ((x[0]*y[0]-x[1]*y[1])%3,(x[0]*y[1]+x[1]*y[0])%3)
def inv(x):
    n=(x[0]*x[0]+x[1]*x[1])%3  # norm, nonzero
    ni=1 if n==1 else 2
    return ((x[0]*ni)%3,(-x[1]*ni)%3)
ZERO=(0,0); ONE=(1,0)
def padd(p,q):
    n=max(len(p),len(q)); p=p+[ZERO]*(n-len(p)); q=q+[ZERO]*(n-len(q))
    return [add(a,b) for a,b in zip(p,q)]
def pmul(p,q):
    r=[ZERO]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): r[i+j]=add(r[i+j],mul(a,b))
    return r
def pscal(c,p): return [mul(c,a) for a in p]
def peval(p,t):
    v=ZERO
    for a in reversed(p): v=add(mul(v,t),a)
    return v
def pdeg(p):
    d=-1
    for i,a in enumerate(p):
        if a!=ZERO: d=i
    return d
def rank_mod3(A):
    A=np.array(A,dtype=np.int64)%3
    if A.size==0: return 0, A
    A=A.copy(); r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        if A[r,c]==2: A[r]=(A[r]*2)%3
        nz=np.nonzero(A[:,c])[0]
        for i in nz:
            if i!=r: A[i]=(A[i]-A[i,c]*A[r])%3
        r+=1
        if r==rows: break
    return r, A[:r]
def rank_Q(A):
    A=np.array(A,dtype=float)
    return 0 if A.size==0 else int(np.linalg.matrix_rank(A))

def build(Q, m, nodes, r=2, verbose=True):
    """r=2. coordinates u in ([2]x[2])^Q ; index = sum 4^s*(2a+d). Returns dict of objects."""
    # position eigenvectors on pairs (a,d)->index 2a+d : e00,e11,e01+e10 (sign +), e01-e10 (sign -)
    pos=[(np.array([1,0,0,0]),+1),(np.array([0,0,0,1]),+1),(np.array([0,1,1,0]),+1),(np.array([0,1,-1,0]),-1)]
    D=4**Q
    Psi=np.zeros((D,D),dtype=np.int64); signs=[]
    for k,choice in enumerate(itertools.product(range(4),repeat=Q)):
        v=np.array([1])
        for s in range(Q):           # position s is digit s (least significant first)
            v=np.kron(pos[choice[s]][0],v) if s>0 else pos[choice[s]][0].copy()
        Psi[:,k]=v; signs.append(tuple(pos[c][1] for c in choice))
    norms=(Psi*Psi).sum(axis=0)       # powers of 2
    # blocks of consecutive weights of length m (weight = number of -1), offset chosen = 0
    blocks=[list(range(b,min(b+m,Q+1))) for b in range(0,Q+1,m)]
    blk_of={w:bi for bi,B in enumerate(blocks) for w in B}
    # coefficients a_i (i<m) with sum_i a_i C(w,i) = (-1)^w mod 3 on the block
    coef=[]
    for B in blocks:
        sol=None
        for a in itertools.product(range(3),repeat=m):
            if all(sum(a[i]*math.comb(w,i) for i in range(m))%3==((-1)**w)%3 for w in B): sol=a;break
        assert sol is not None, "no F_3 polynomial for block"
        coef.append(sol)
    # P_z(T) = 1 + chi(z) [F_B(y(T)) - F_B(y(0))],  y_j = x_j (1 - T/tau_j)
    lin=[[ONE, mul(f9(2,0),inv(nodes[j]))] for j in range(Q)]   # 1 - T/tau_j
    Pz={}
    for z in set(signs):
        supp=[j for j in range(Q) if z[j]==-1]; w=len(supp); a=coef[blk_of[w]]
        # elementary symmetric polys of the linear factors on supp
        e=[[ONE]]+[[ZERO] for _ in range(m)]
        for j in supp:
            for i in range(min(m-1,len(supp)),0,-1):
                e[i]=padd(e[i],pmul(e[i-1],lin[j]))
        F=[ZERO]
        for i in range(m): F=padd(F,pscal(f9(a[i],0),e[i]))
        F0=peval(F,ZERO)
        chi=f9((-1)**w,0)
        G=padd(F,[mul(f9(2,0),F0)])
        Pz[z]=padd([ONE],pscal(chi,G))
        assert peval(Pz[z],ZERO)==ONE
        assert F0==chi, (z,F0,chi)
    maxdeg=max(pdeg(p) for p in Pz.values())
    ninv=np.array([pow(int(x),-1,3) for x in norms])
    def Pmat(t):
        vals=[peval(Pz[z],t) for z in signs]
        out=[]
        for comp in (0,1):
            d=np.array([v[comp] for v in vals])*ninv%3
            out.append(((Psi%3)*d[None,:]%3) @ (Psi.T%3) %3)
        return out   # (real, imag) D x D mod 3
    return dict(D=D,Psi=Psi,signs=signs,Pz=Pz,Pmat=Pmat,maxdeg=maxdeg,blocks=blocks,blk_of=blk_of,norms=norms)

def bad_rows(obj,Q,t):
    """indices of basis vectors psi whose sign pattern z_{-t} is a bad t-edge (computed from P_z, not assumed)."""
    Pz=obj['Pz']; bad=[]
    return bad

def run(Q,m,nodes):
    print(f"=== r=2 Q={Q} block length m={m} (degree <= {m-1}), field F_9, nodes={nodes[:Q]}")
    obj=build(Q,m,nodes); D=obj['D']; signs=obj['signs']; Psi=obj['Psi']
    print("  D =",D,"; max deg_T P_z =",obj['maxdeg'],"; blocks",obj['blocks'])
    P0=obj['Pmat'](ZERO); assert (P0[0]==np.eye(D,dtype=np.int64)).all() and not P0[1].any(); print("  P(0)=I  OK")
    idx=np.arange(D)
    a_of=lambda s:(idx>>(2*s+1))&1; d_of=lambda s:(idx>>(2*s))&1
    tot=[]
    for t in range(Q):
        tau=nodes[t]
        # bad edges: (P_z+P_{z^t})(tau) != 0
        badset=set()
        for z in set(signs):
            z2=tuple(-z[j] if j==t else z[j] for j in range(Q))
            if add(peval(obj['Pz'][z],tau),peval(obj['Pz'][z2],tau))!=ZERO: badset.add(z)
        # predicted: crossing block boundary
        pred={z for z in set(signs) if obj['blk_of'][sum(1 for j in range(Q) if z[j]==-1)]!=
              obj['blk_of'][sum(1 for j in range(Q) if (z[j]==-1)!=(j==t))]}
        assert badset<=pred
        rows=[k for k,z in enumerate(signs) if z in badset]
        C=Psi[:,rows].T                       # c_t x D, entries 0,+-1
        ct=len(rows)
        Pt=obj['Pmat'](tau)
        A=a_of(t); Dd=d_of(t)
        ok=True; rk3=[]; rkQ=[]; nzrows=[]
        for g in (0,1):
            cols=np.nonzero(A==g)[0]
            Cg=C[:,cols]
            r3,basis=rank_mod3(Cg); rk3.append(r3); rkQ.append(rank_Q(Cg)); nzrows.append(int((np.abs(Cg).sum(axis=1)>0).sum()))
            urows=np.nonzero(Dd!=g)[0]
            for comp in (0,1):
                R=Pt[comp][np.ix_(urows,cols)]
                R=np.unique(R,axis=0)
                r2,_=rank_mod3(np.vstack([basis,R]) if r3>0 else R)
                if r2!=r3: ok=False
        tot.append((ct,rk3,rkQ,nzrows))
        print(f"  t={t}: |bad|={len(badset)} (pred {len(pred)}) c_t={ct} rho_t={ct/D:.4f}  nonzero rows per group {nzrows}"
              f"  rank_Q per group {rkQ}  rank_3 per group {rk3}   kernel condition (Def 3.1(4)) over F_9: {'OK' if ok else 'FAIL'}")
        assert ok
    return obj

def reach(Q,m,nodes,nfill,seed=0):
    """reachable measurement vectors for group g=0 at t=0, single copy, nfill random coordinates of S_{t,g} filled
    with independent Boolean variables (others 0)."""
    obj=build(Q,m,nodes,verbose=False); D=obj['D']; signs=obj['signs']; Psi=obj['Psi']
    t=0; idx=np.arange(D); A=(idx>>(2*t+1))&1
    badrows=[k for k,z in enumerate(signs) if obj['blk_of'][sum(1 for j in range(Q) if z[j]==-1)]!=
              obj['blk_of'][sum(1 for j in range(Q) if (z[j]==-1)!=(j==t))]]
    cols=np.nonzero(A==0)[0]
    r=np.random.default_rng(seed); cols=r.choice(cols,size=nfill,replace=False)
    C=Psi[:,badrows].T[:,cols]
    C=C[np.abs(C).sum(axis=1)>0]
    rkQ=rank_Q(C); rk3,_=rank_mod3(C)
    X=((np.arange(2**nfill)[:,None]>>np.arange(nfill))&1).astype(np.int64)
    Y=X@C.T
    nZ=len(np.unique(Y,axis=0)); n3=len(np.unique(Y%3,axis=0))
    # greedy basis with smallest l1 norms
    order=np.argsort(np.abs(C).sum(axis=1),kind='stable'); basis=[]; cur=0
    for i in order:
        if rank_Q(C[basis+[i]])>cur: basis.append(i); cur+=1
        if cur==rkQ: break
    prodnorm=sum(math.log2(int(np.abs(C[i]).sum())+1) for i in basis)
    Dp=4**Q
    print(f"  Q={Q} m={m} filled={nfill}: nonzero rows={len(C)} rank_Q={rkQ} rank_3={rk3} | log2 #reachable over Z = {math.log2(nZ):.2f}"
          f" (2^rank: {rkQ}, basis-norm bound: {prodnorm:.2f}, cube: {nfill}, P's box (2D+1)^rows: {len(C)*math.log2(2*Dp+1):.0f})"
          f" | log3 #reachable mod 3 = {math.log(n3,3):.2f} (<= rank_3 = {rk3})")
    assert nZ>=2**rkQ-0 and math.log2(nZ)<=prodnorm+1e-9 and n3<=3**rk3

if __name__=="__main__":
    nodes=[f9(1,0),f9(2,0),f9(0,1),f9(0,2),f9(1,1),f9(1,2),f9(2,1),f9(2,2)]
    run(3,3,nodes)
    run(4,3,nodes)
    run(5,3,nodes)
    print("=== reachable measurement vectors (group 0 at t=0)")
    for Q,nf in ((3,16),(4,18),(4,20),(5,20)):
        reach(Q,3,nodes,nf)
