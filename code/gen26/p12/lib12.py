# Independent library (P12): exact Cauchy blocks, Smith valuations over Z_(2) via mod 2^N elimination.
import itertools, math
from fractions import Fraction
import numpy as np
def v2(n):
    n=abs(int(n)); return 10**9 if n==0 else (n&-n).bit_length()-1
def v2q(q):
    q=Fraction(q); return 10**9 if q==0 else v2(q.numerator)-v2(q.denominator)
def toZ(q,MOD):
    q=Fraction(q); assert q.denominator%2==1, q; return q.numerator*pow(q.denominator,-1,MOD)%MOD
def smith_vals(A,N):
    """A: list of lists of ints (2-adic integers mod 2^N). returns sorted valuations of invariants, capped at N (rank-deficient => N)."""
    MOD=1<<N
    M=np.array([[int(x)%MOD for x in row] for row in A],dtype=object)
    r,c=M.shape; out=[]
    for k in range(min(r,c)):
        sub=M[k:,k:]
        best=None
        # find min valuation entry
        nz=np.nonzero(sub!=0)
        if len(nz[0])==0:
            out+= [N]*(min(r,c)-k); break
        bv=N+1;bp=None
        for (ii,jj) in zip(*nz):
            vv=v2(sub[ii,jj])
            if vv<bv:
                bv=vv;bp=(ii+k,jj+k)
                if vv==0: break
        pi,pj=bp
        M[[k,pi],:]=M[[pi,k],:]; M[:,[k,pj]]=M[:,[pj,k]]
        p=int(M[k,k]); u=p>>bv; uinv=pow(u,-1,MOD)
        M[k,:]=(M[k,:]*uinv)%MOD   # pivot now 2^bv
        col=M[k+1:,k]
        f=np.array([int(x)>>bv for x in col],dtype=object)
        if len(f):
            M[k+1:,:]=(M[k+1:,:]-np.outer(f,M[k,:]))%MOD
        out.append(bv)
    return sorted(out)
def cauchy_block(a,t,zeta,b):
    L=1<<a; i=v2(t); T=t+(1<<(a+1))*Fraction(zeta)
    Cp0=1
    for q in range(-L,L):
        if q: Cp0*=-q
    CT=Fraction(1)
    for q in range(-L,L): CT*=(T-q)
    B=[x for x in range(L) if (x>>i)&1==b]
    c=CT/Cp0
    E=[[c/(T+x-y) for y in B] for x in B]
    return E,c,T,B
def per(Mx):
    n=len(Mx)
    if n==0: return Fraction(1)
    # Ryser-free simple expansion via DP over subsets of columns
    dp={0:Fraction(1)}
    for r in range(n):
        nd={}
        for mask,val in dp.items():
            for j in range(n):
                if not mask>>j&1:
                    nm=mask|1<<j; nd[nm]=nd.get(nm,0)+val*Mx[r][j]
        dp=nd
    return dp[(1<<n)-1]
def det(Mx):
    n=len(Mx); A=[list(map(Fraction,r)) for r in Mx]; d=Fraction(1)
    for k in range(n):
        p=next((i for i in range(k,n) if A[i][k]!=0),None)
        if p is None: return Fraction(0)
        if p!=k: A[k],A[p]=A[p],A[k]; d=-d
        d*=A[k][k]
        for i in range(k+1,n):
            f=A[i][k]/A[k][k]
            if f: A[i]=[x-f*y for x,y in zip(A[i],A[k])]
    return d
def Mmat(E):
    """M[S,X'] = prod_{x in S} (E X')_x, rows/cols indexed by bitmasks."""
    n=len(E); N=1<<n
    lin=[[Fraction(0)]*N for _ in range(n)]
    for x in range(n):
        for X in range(1,N):
            low=X&-X; k=low.bit_length()-1; lin[x][X]=lin[x][X^low]+E[x][k]
    M=[[None]*N for _ in range(N)]
    for S in range(N):
        xs=[x for x in range(n) if S>>x&1]
        for X in range(N):
            p=Fraction(1)
            for x in xs: p*=lin[x][X]
            M[S][X]=p
    return M
def compound(E,k,f):
    n=len(E); subs=list(itertools.combinations(range(n),k))
    return [[f([[E[r][cc] for cc in C] for r in R]) for C in subs] for R in subs]
