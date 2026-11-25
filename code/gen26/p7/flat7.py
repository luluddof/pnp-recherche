# Own check of C1 counterexample: flattening of one packet tensor across the two groups
import sys, itertools, numpy as np
from fractions import Fraction
sys.path.insert(0,'.')
from smith7 import v2, red
NB=30; MOD=1<<NB
def smith_np(A):
    A=A.copy()%MOD; n,m=A.shape; out=[]
    while A.size:
        nz=A!=0
        if not nz.any(): break
        low=A & (-A); low[~nz]=1<<62
        r,c=np.unravel_index(np.argmin(low),low.shape); v=int(low[r,c]).bit_length()-1
        out.append(v)
        u=int(A[r,c])>>v; ui=pow(u,-1,MOD)
        col=(A[:,c]>>v)%MOD; f=(col*ui)%MOD; f[r]=0
        # A -= outer(f, A[r]) mod MOD ; avoid overflow: values <2^30, products <2^60
        A=(A-(np.outer(f,A[r])%MOD))%MOD
        A=np.delete(np.delete(A,r,0),c,1)
    return sorted(out)
def flatten(a,i,t,zeta):
    L=1<<a; T=t+Fraction(2**(a+1))*zeta
    C=Fraction(1)
    for q in range(-L,L): C*=(T-q)
    Cp0=Fraction(1)
    for q in range(-L,L):
        if q: Cp0*=-q
    W=np.array([[red(C/(Cp0*(T-(y-x))),NB) for y in range(L)] for x in range(L)],dtype=np.int64)
    rows_g=[x for x in range(L) if (x>>i)&1==0]; cols_g=[y for y in range(L) if (y>>i)&1==1]
    rows_h=[x for x in range(L) if (x>>i)&1==1]; cols_h=[y for y in range(L) if (y>>i)&1==0]
    def assignments(rows,cols):
        res=[]
        for bits in itertools.product((0,1),repeat=len(rows)+len(cols)):
            res.append((dict(zip(rows,bits[:len(rows)])),dict(zip(cols,bits[len(rows):]))))
        return res
    Ag=assignments(rows_g,cols_g); Ah=assignments(rows_h,cols_h)
    F=np.zeros((len(Ag),len(Ah)),dtype=np.int64)
    for I,(Yg,Xg) in enumerate(Ag):
        for J,(Yh,Xh) in enumerate(Ah):
            Y={**Yg,**Yh}; X={**Xg,**Xh}
            Xv=np.array([X[y] for y in range(L)],dtype=np.int64)
            WX=(W@Xv)%MOD  # entries <2^30*8 fine
            prod=1
            for x in range(L):
                z=int(WX[x]); fac=z if Y[x] else (1-z)
                prod=(prod*fac)%MOD
            F[I,J]=prod
    return F
def monomial_count(svals,p):
    # number of beta in N^modes with sum s beta < p
    cnt=[0]*p; cnt[0]=1
    for s in svals:
        for w in range(s,p): cnt[w]+=cnt[w-s]
    return sum(cnt)
if __name__=="__main__":
    a,i=int(sys.argv[1]),int(sys.argv[2])
    for t in [int(x) for x in sys.argv[3].split(',')]:
        F=flatten(a,i,t,1); sv=smith_np(F)
        from smith7 import predicted
        modes=predicted(a,i,1)*2
        print("a",a,"i",i,"t",t,"modes",modes)
        print(" p : R_p(F)  monomials(C1)")
        for p in range(1,17):
            print(" ",p, sum(1 for s in sv if s<p), monomial_count(modes,p))
