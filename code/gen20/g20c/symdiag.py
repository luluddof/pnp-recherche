# Symmetric Pauli-diagonal family.  N = 2^n.  Points a in F_2^{2n} = (F_2^2)^n, classes t=(n0,ns,nns) (counts of
# coordinate types 00 / {01,10} / 11).  chi(a) = (-1)^{a_x.a_z} = (-1)^{nns}.  Characters chi_c(a)=(-1)^{omega(c,a)},
# classes of c likewise.  K[t_a][t_c] = sum_{c in class t_c} chi_c(a) = coeff of v^{ms} w^{mns} in
#  (1+2v+w)^{n0} (1-w)^{ns} (1-2v+w)^{nns}.
# Problem: choose a set Qc of character classes (q = total size) and a set A of point classes such that
# chi|_A is in span{K[.,C] : C in Qc}|_A ; p = 4^n - |A|.  (Symmetrisation: an invariant solution exists iff any does,
# for invariant A, Qc.)  Greedy over classes by decreasing size, consistency checked mod a large prime.
import sys, math, itertools
from math import comb, factorial
P=(1<<61)-1
def classes(n):
    return [(n0,ns,n-n0-ns) for n0 in range(n+1) for ns in range(n+1-n0)]
def size(t):
    n0,ns,nn=t; n=n0+ns+nn
    return factorial(n)//(factorial(n0)*factorial(ns)*factorial(nn))*2**ns
def polypow(base,e,n):  # base: dict (i,j)->coef ; truncated deg<=n
    res={(0,0):1}
    for _ in range(e):
        new={}
        for (i,j),c in res.items():
            for (k,l),d in base.items():
                if i+k<=n and j+l<=n:
                    new[(i+k,j+l)]=new.get((i+k,j+l),0)+c*d
        res=new
    return res
def mul(a,b,n):
    new={}
    for (i,j),c in a.items():
        for (k,l),d in b.items():
            if i+k<=n and j+l<=n: new[(i+k,j+l)]=new.get((i+k,j+l),0)+c*d
    return new
def Kmatrix(n):
    cl=classes(n)
    B0={(0,0):1,(1,0):2,(0,1):1}; Bs={(0,0):1,(0,1):-1}; Bn={(0,0):1,(1,0):-2,(0,1):1}
    P0=[polypow(B0,e,n) for e in range(n+1)]
    Ps=[polypow(Bs,e,n) for e in range(n+1)]
    Pn=[polypow(Bn,e,n) for e in range(n+1)]
    K={}
    for t in cl:
        pol=mul(mul(P0[t[0]],Ps[t[1]],n),Pn[t[2]],n)
        K[t]={c: pol.get((c[1],c[2]),0) for c in cl}
    return cl,K
def greedy(n,cl,K,Qc,order=None):
    # incremental elimination mod P on rows [K[t][Qc] | chi]
    cols=list(Qc); m=len(cols)
    basis={}  # pivot col -> row (list length m+1)
    A=[]
    if order is None: order=sorted(cl,key=lambda t:-size(t))
    for t in order:
        row=[K[t][c]%P for c in cols]+[(1 if t[2]%2==0 else P-1)]
        for pc,br in basis.items():
            if row[pc]:
                f=row[pc]
                row=[(x-f*y)%P for x,y in zip(row,br)]
        piv=next((i for i in range(m) if row[i]),None)
        if piv is None:
            if row[m]: continue   # inconsistent: reject class
            A.append(t); continue
        inv=pow(row[piv],P-2,P); row=[x*inv%P for x in row]
        for pc in list(basis):
            br=basis[pc]
            if br[piv]:
                f=br[piv]; basis[pc]=[(x-f*y)%P for x,y in zip(br,row)]
        basis[piv]=row; A.append(t)
    return A
if __name__=="__main__":
    n=int(sys.argv[1])
    cl,K=Kmatrix(n)
    N2=4**n; lN=n
    print("n",n,"classes",len(cl))
    out=[]
    # Q shapes: support <= D (ms+mns<=D), and boxes (ms<=D1, mns<=D2) with m0 free
    shapes=[]
    for D in range(0,n+1): shapes.append(("supp<=%d"%D,[c for c in cl if c[1]+c[2]<=D]))
    for D1 in range(0,n+1,max(1,n//6)):
        for D2 in range(0,n+1,max(1,n//6)):
            shapes.append(("ms<=%d,mns<=%d"%(D1,D2),[c for c in cl if c[1]<=D1 and c[2]<=D2]))
    for name,Qc in shapes:
        q=sum(size(c) for c in Qc)
        A=greedy(n,cl,K,Qc)
        p=N2-sum(size(t) for t in A)
        th=math.log2(max(p,1))/lN; T=math.log2(max(q,1))/lN
        circ=math.sqrt(max(0,4-th*th))
        print(f"{name:18s} p=2^{math.log2(max(p,1)):.2f} q=2^{math.log2(max(q,1)):.2f} theta={th:.3f} T={T:.3f} th+T={th+T:.3f} circleT={circ:.3f} {'BELOW' if T<circ-1e-9 and p>0 and q>0 else ''}",flush=True)
