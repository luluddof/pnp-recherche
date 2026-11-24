"""G13a: independent exact checks (GF(p) and Q) of the dual lower bound and of a FIELD-FREE bilinear variant.
Setting: J(t)=X^{(x)l}+sum_{i=1..d} t^i B_i (2^l x 2^l), points tau_k (direction k), P_k := J(tau_k)[x_k=0, y_k=1].
For any decomposition J(tau_k)=K_k+E_k (K_k zero on x_k!=y_k):  rank E_k >= rank P_k.
Space A(ma,mc): a_x(t)=prod_{k:x_k=1}(1-t/tau_k) r_x(t), deg a_x<=ma, r_x=0 for |x|>mc.   ma+mc+d<=l-1, ma>=mc.
Claim (bilinear variant): the map a -> (a(tau_k)[x_k=0] P_k)_k is injective on A(ma,mc)  (any field).
Hence sum_k rank P_k >= dim A(ma,mc) = sum_{w<=mc} C(l,w)(ma+1-w).   (ma=mc=m gives N_m.)
"""
import random, sys, itertools
from math import comb
def rank_mod(rows,p):
    M=[r[:] for r in rows]; r=0
    if not M: return 0
    nc=len(M[0])
    for c in range(nc):
        piv=None
        for i in range(r,len(M)):
            if M[i][c]%p: piv=i;break
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]; inv=pow(M[r][c],p-2,p); M[r]=[x*inv%p for x in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]%p:
                f=M[i][c]; M[i]=[(a-f*b)%p for a,b in zip(M[i],M[r])]
        r+=1
        if r==len(M): break
    return r
def polymul(a,b,p):
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]=(r[i+j]+x*y)%p
    return r
def peval(a,t,p):
    v=0
    for c in reversed(a): v=(v*t+c)%p
    return v
def basisA(l,ma,mc,tau,p):
    """list of (x, poly coeffs of a_x) basis elements"""
    B=[]
    for x in range(1<<l):
        w=bin(x).count('1')
        if w>mc or w>ma: continue
        pi=[1]
        for k in range(l):
            if (x>>k)&1: pi=polymul(pi,[1,(-pow(tau[k],p-2,p))%p],p)
        for i in range(ma-w+1):
            B.append((x,[0]*i+pi))
    return B
def Jeval(Bc,t,p,n):
    # Bc[i][x][y], i=0..d
    return [[sum(Bc[i][x][y]*pow(t,i,p) for i in range(len(Bc)))%p for y in range(n)] for x in range(n)]
def bound(l,d,ma,mc): return sum(comb(l,w)*(ma+1-w) for w in range(mc+1))
def best_bound(l,d):
    best=(0,None)
    for mc in range(0,l):
        ma=l-1-d-mc
        if ma<mc: break
        b=bound(l,d,ma,mc)
        if b>best[0]: best=(b,(ma,mc))
    return best
def Nm(l,d):
    m=(l-d-1)//2
    return bound(l,d,m,m) if m>=0 else 0
def randomJ(l,d,p,lowrank=None):
    n=1<<l
    Bc=[[[0]*n for _ in range(n)]]+[[[random.randrange(p) for _ in range(n)] for _ in range(n)] for _ in range(d)]
    for x in range(n): Bc[0][x][(n-1)^x]=1
    return Bc
def check_instance(l,d,p,tau,Bc,ma,mc):
    """returns (dim A, rank of map a->(a(tau_k)P_k)_k, sum_k rank P_k, max violation flag of Phi-identity)"""
    n=1<<l; B=basisA(l,ma,mc,tau,p); rows=[[] for _ in B]; srank=0
    for k in range(l):
        Jt=Jeval(Bc,tau[k],p,n)
        xs=[x for x in range(n) if not (x>>k)&1]; ys=[y for y in range(n) if (y>>k)&1]
        Pk=[[Jt[x][y] for y in ys] for x in xs]; srank+=rank_mod(Pk,p)
        for bi,(x,pol) in enumerate(B):
            if (x>>k)&1: v=0; assert peval(pol,tau[k],p)==0
            else: v=peval(pol,tau[k],p)
            # row vector a(tau_k) restricted to x_k=0 has single entry v at x ; times Pk
            if (x>>k)&1: rows[bi]+= [0]*len(ys)
            else:
                ix=xs.index(x); rows[bi]+=[v*Pk[ix][j]%p for j in range(len(ys))]
    rk=rank_mod(rows,p) if B else 0
    return len(B),rk,srank
def phi_identity(l,d,p,tau,Bc,ma,mc,trials=5):
    """Phi(t)=a(t)J(t)c(t)^T with a in A(ma,mc) random and c=pi_{x0} e_{bar x0}, |x0|<=mc.  Check:
       Phi(tau_k) == a(tau_k) E_k c(tau_k) for E_k := J(tau_k) restricted to x_k!=y_k  (K-part kills), Phi(0)=a_{x0}(0)."""
    n=1<<l; B=basisA(l,ma,mc,tau,p); ok=True
    for _ in range(trials):
        co=[random.randrange(p) for _ in B]
        def a(x,t): return sum(c*peval(pol,t,p) for c,(xx,pol) in zip(co,B) if xx==x)%p
        x0=random.choice([x for x in range(n) if bin(x).count('1')<=mc])
        def pi0(t):
            v=1
            for k in range(l):
                if (x0>>k)&1: v=v*(1-t*pow(tau[k],p-2,p))%p
            return v
        y0=(n-1)^x0
        # Phi(0)
        J0=Jeval(Bc,0,p,n)
        if sum(a(x,0)*J0[x][y0] for x in range(n))*pi0(0)%p!=a(x0,0): ok=False
        for k in range(l):
            Jt=Jeval(Bc,tau[k],p,n)
            full=sum(a(x,tau[k])*Jt[x][y0] for x in range(n))*pi0(tau[k])%p
            Epart=sum(a(x,tau[k])*Jt[x][y0] for x in range(n) if ((x^y0)>>k)&1)*pi0(tau[k])%p
            Kpart=(full-Epart)%p
            if Kpart!=0: ok=False
    return ok
if __name__=="__main__":
    random.seed(13)
    print("bounds: l d | N_m (claimed) | best asymmetric bound (ma,mc)")
    for l in range(2,9):
        for d in range(0,l):
            print(f"  {l} {d} | {Nm(l,d)} | {best_bound(l,d)}")
    for p in (5,7,11,13,101,(1<<31)-1):
        allok=True; cnt=0
        for l in (2,3,4):
            if p-1<l: continue
            for d in range(0,l):
                for mc in range(0,l):
                    ma=l-1-d-mc
                    if ma<mc: break
                    for tr in range(6 if l<4 else 2):
                        tau=random.sample(range(1,p),l)
                        Bc=randomJ(l,d,p)
                        dim,rk,sr=check_instance(l,d,p,tau,Bc,ma,mc)
                        idok=phi_identity(l,d,p,tau,Bc,ma,mc,2)
                        cnt+=1
                        if not(rk==dim and sr>=dim and idok):
                            allok=False; print("  VIOLATION",p,l,d,ma,mc,dim,rk,sr,idok)
        print(f"GF({p}): random J, all (l<=4,d,ma,mc): injective & sum rank P_k >= dim A & Phi-identity: {allok} ({cnt} instances)")
