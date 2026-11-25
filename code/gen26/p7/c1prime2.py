import sys, numpy as np; sys.path.insert(0,'.')
import flat7
from smith7 import predicted
flat7.NB=90; flat7.MOD=1<<90
NB=90; MOD=1<<90
def smith_obj(A):
    A=A.astype(object)%MOD; out=[]
    while A.size:
        best=None
        n,m=A.shape
        for r in range(n):
            for c in range(m):
                x=int(A[r,c])
                if x:
                    v=(x&-x).bit_length()-1
                    if best is None or v<best[0]: best=(v,r,c)
                    if v==0: break
            if best and best[0]==0: break
        if best is None: break
        v,r,c=best; out.append(v); u=int(A[r,c])>>v; ui=pow(u,-1,MOD)
        f=np.array([((int(A[k,c])>>v)*ui)%MOD for k in range(n)],dtype=object); f[r]=0
        A=(A-np.outer(f,A[r]))%MOD
        A=np.delete(np.delete(A,r,0),c,1)
    return sorted(out)
def flatten_obj(a,i,t):
    # rebuild with python ints mod 2^90
    import itertools
    from fractions import Fraction
    from smith7 import red
    L=1<<a; T=t+Fraction(2**(a+1))
    C=Fraction(1)
    for q in range(-L,L): C*=(T-q)
    Cp0=Fraction(1)
    for q in range(-L,L):
        if q: Cp0*=-q
    W=[[red(C/(Cp0*(T-(y-x))),NB) for y in range(L)] for x in range(L)]
    rg=[x for x in range(L) if (x>>i)&1==0]; cg=[y for y in range(L) if (y>>i)&1==1]
    rh=[x for x in range(L) if (x>>i)&1==1]; ch=[y for y in range(L) if (y>>i)&1==0]
    def asg(rows,cols):
        return [(dict(zip(rows,b[:len(rows)])),dict(zip(cols,b[len(rows):]))) for b in itertools.product((0,1),repeat=len(rows)+len(cols))]
    Ag=asg(rg,cg); Ah=asg(rh,ch)
    F=np.zeros((len(Ag),len(Ah)),dtype=object)
    for I,(Yg,Xg) in enumerate(Ag):
        for J,(Yh,Xh) in enumerate(Ah):
            Y={**Yg,**Yh}; X={**Xg,**Xh}; pr=1
            for x in range(L):
                z=sum(W[x][y] for y in range(L) if X[y])%MOD
                pr=pr*(z if Y[x] else 1-z)%MOD
            F[I,J]=pr
    return F
def sqfree(modes,p):
    cnt=[0]*p; cnt[0]=1
    for s in modes:
        for w in range(p-1,s-1,-1): cnt[w]+=cnt[w-s]
    return sum(cnt)
a,i,t=[int(x) for x in sys.argv[1:4]]
F=flatten_obj(a,i,t); sv=smith_obj(F); modes=predicted(a,i,1)*2
print("max inv",max(sv),"count",len(sv))
bad=[]
for p in range(1,88):
    R=sum(1 for s in sv if s<p); lo=sqfree(modes,p); hi=flat7.monomial_count(modes,p)
    if not (lo<=R<=hi): bad.append((p,lo,R,hi))
print(a,i,t,"violations of sqfree<=R_p<=monomials:",bad[:10])
print("Smith(F):",sv)
