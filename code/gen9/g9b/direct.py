"""Direct check, lambda=1 one-bit identity X = [x=u][y=v] - g(x,y)g(u,v), g=[x!=y](-1)^x  (c=1, valid over Z).
Pieces A_i = coeff of T^i in (R + T C)^{(x) t}.  Ranks of row/col flattenings over Q, GF(2), GF(3)."""
import numpy as np, itertools, sys
from math import comb, log2
def crossing():
    X=np.zeros((2,2,2,2),dtype=np.int64)
    for x in range(2):
        for y in range(2): X[x,y,y,x]=1
    return X
X=crossing()
C=np.zeros((2,2,2,2),dtype=np.int64); R=np.zeros((2,2,2,2),dtype=np.int64)
g=lambda x,y:(x!=y)*(-1)**x
for x,y,u,v in itertools.product(range(2),repeat=4):
    C[x,y,u,v]=(x==u)*(y==v); R[x,y,u,v]=-g(x,y)*g(u,v)
assert (R+C==X).all()
def rowm(A): return A.transpose(0,1,3,2).reshape(4,4)   # (x,y),(v,u)
def colm(A): return A.transpose(0,2,3,1).reshape(4,4)   # (x,u),(v,y)
def rank_p(M,p):
    M=(M%p).astype(np.int64).copy(); r=0; rows,cols=M.shape
    for c in range(cols):
        piv=None
        nz=np.nonzero(M[r:,c])[0]
        if len(nz)==0: continue
        piv=r+nz[0]
        if piv!=r: M[[r,piv]]=M[[piv,r]]
        inv=pow(int(M[r,c]),p-2,p); M[r]=(M[r]*inv)%p
        idx=np.nonzero(M[:,c])[0]; idx=idx[idx!=r]
        M[idx]=(M[idx]-np.outer(M[idx,c],M[r]))%p
        r+=1
        if r==rows: break
    return r
def pieces(Rm,Cm,t):
    # polynomial in T with matrix coefficients: kron product
    polys=[np.array([[1]],dtype=np.int64)]  # list index = power of T
    for s in range(t):
        new=[None]*(len(polys)+1)
        for i,P in enumerate(polys):
            a=np.kron(P,Rm); b=np.kron(P,Cm)
            new[i]=a if new[i] is None else new[i]+a
            new[i+1]=b if new[i+1] is None else new[i+1]+b
        polys=new
    return polys
for t in range(1,6):
    Pr=pieces(rowm(R),rowm(C),t); Pc=pieces(colm(R),colm(C),t)
    assert (sum(Pr)==np.eye(4**t,dtype=np.int64)).all() and (sum(Pc)==np.eye(4**t,dtype=np.int64)).all()
    out=[]
    for i in range(t+1):
        rq=(np.linalg.matrix_rank(Pr[i].astype(float)),np.linalg.matrix_rank(Pc[i].astype(float)))
        r2=(rank_p(Pr[i],2),rank_p(Pc[i],2)); r3=(rank_p(Pr[i],3),rank_p(Pc[i],3))
        out.append((i,rq,r2,r3))
    print("t=%d"%t)
    for i,rq,r2,r3 in out: print("   i=%d  Q %s  GF2 %s  GF3 %s"%(i,rq,r2,r3))
    for name,k in (("Q",1),("GF2",2),("GF3",3)):
        m=max(log2(o[k][0])+log2(o[k][1]) for o in out); print("   %s: max_i (crow+ccol)/(4t) = %.4f"%(name,m/(4*t)))
