#!/usr/bin/env python3
"""G18b exp_reach: number of distinct measurement vectors C x, x Boolean, for (1) the tightness examples and
(2) dense sub-boxes of the auxiliary coordinates (rows = bad basis vectors psi restricted to a group)."""
import math, itertools, numpy as np
from exp_aux import build, rank_Q, rank_mod3, f9
def count(C,label):
    C=C[np.abs(C).sum(axis=1)>0]; nf=C.shape[1]
    rkQ=rank_Q(C); rk3,_=rank_mod3(C)
    C=np.unique(C,axis=0)
    X=((np.arange(2**nf)[:,None]>>np.arange(nf))&1).astype(np.int16)
    _,piv=rank_mod3(C.T)  # not used
    # keep a rational row basis for the Z-count (other rows are determined by it)
    bas=[]; cur=0
    for i in range(len(C)):
        if rank_Q(C[bas+[i]])>cur: bas.append(i); cur+=1
    YZ=X@C[bas].T.astype(np.int16); nZ=len(np.unique(YZ,axis=0)); del YZ
    r3,B3=rank_mod3(C)
    Y=(X@B3.T.astype(np.int16)) if r3 else np.zeros((len(X),1),dtype=np.int16); n3=len(np.unique(Y%3,axis=0))
    order=np.argsort(np.abs(C).sum(axis=1),kind='stable'); basis=[]; cur=0
    for i in order:
        if rank_Q(C[basis+[i]])>cur: basis.append(i); cur+=1
        if cur==rkQ: break
    pn=sum(math.log2(int(np.abs(C[i]).sum())+1) for i in basis)
    Dmax=int(np.abs(C).sum(axis=1).max())
    print(f"  {label}: cols={nf} rows={len(C)} rank={rkQ} rank_3={rk3} maxnorm={Dmax} | log2#Z={math.log2(nZ):.2f}"
          f"  [rank={rkQ} <= . <= basis-norm bound {pn:.2f} <= rank*log2(D+1)={rkQ*math.log2(Dmax+1):.2f};"
          f" box rows*log2(2D+1)={len(C)*math.log2(2*Dmax+1):.1f}] | log2 #mod3 = {math.log2(n3):.2f} <= rank_3*log2(3) = {rk3*math.log2(3):.2f}")
    assert 2**rkQ<=nZ and math.log2(nZ)<=pn+1e-9 and n3<=3**rk3
print("(1) tightness: c blocks of D ones (image exactly (D+1)^c); identity (image 2^c)")
for c,Dn in ((2,8),(3,6),(4,4)):
    C=np.kron(np.eye(c,dtype=np.int64),np.ones((1,Dn),dtype=np.int64)); count(C,f"blocks c={c} D={Dn}")
    assert len(np.unique(((np.arange(2**(c*Dn))[:,None]>>np.arange(c*Dn))&1)@C.T,axis=0))==(Dn+1)**c
count(np.eye(6,dtype=np.int64),"identity c=6")
C=np.kron(np.eye(2,dtype=np.int64),np.array([[1]*5+[-1]*5])); count(C,"signed blocks c=2 (+5,-5)")
print("(2) auxiliary rows, r=2, group 0 at t=0, dense sub-box of coordinates (all positions s>=1 in {01,10,00,11})")
nodes=[f9(1,0),f9(2,0),f9(0,1),f9(0,2),f9(1,1),f9(1,2),f9(2,1),f9(2,2)]
for Q,m,box in ((4,3,[(0,1),(1,2),(1,2),(1,2)]),(4,3,[(1,),(1,2),(1,2,0),(1,2,0)]),(5,3,[(1,),(1,2),(1,2),(1,2),(1,2)]),
                (5,3,[(0,),(1,2),(1,2),(1,2),(1,2)]),(5,3,[(1,),(1,2),(1,2),(1,2,0),(0,)])):
    obj=build(Q,m,nodes)
    signs=obj['signs']; Psi=obj['Psi']; t=0
    bad=[k for k,z in enumerate(signs) if obj['blk_of'][sum(1 for j in range(Q) if z[j]==-1)]!=
              obj['blk_of'][sum(1 for j in range(Q) if (z[j]==-1)!=(j==t))]]
    cols=[sum(4**s*c[s] for s in range(Q)) for c in itertools.product(*box)]
    C=Psi[np.ix_(cols,bad)].T
    count(C,f"Q={Q} m={m} per-point bad rows, box {len(cols)} coords")
