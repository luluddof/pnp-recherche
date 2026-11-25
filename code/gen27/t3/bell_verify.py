# explicit verification: build D*(SWAP-K) as integer matrix for the Bell polynomial family, rank mod two primes
import numpy as np, itertools, sys
from fractions import Fraction
from math import lcm
from bell_poly import interp
from bell_poly_big import best_window
def rank_mod(M,p):
    M=M%p; n,m=M.shape; r=0
    for c in range(m):
        piv=np.nonzero(M[r:,c])[0]
        if len(piv)==0: continue
        i=r+piv[0]; M[[r,i]]=M[[i,r]]
        inv=pow(int(M[r,c]),p-2,p); M[r]=(M[r]*inv)%p
        col=M[:,c].copy(); col[r]=0
        nz=np.nonzero(col)[0]
        if len(nz): M[nz]=(M[nz]-np.outer(col[nz],M[r]))%p
        r+=1
        if r==n: break
    return r
k=int(sys.argv[1]); N=2**k; n=N*N
xs=np.array(list(itertools.product([0,1],repeat=k)))
ints=xs@(1<<np.arange(k)[::-1])
dot=lambda u,v:(u@v.T)
for d in range(0,k+1):
    A=best_window(k,d); f=interp(A)
    # phi on points v=(c,d') : phi = f(c.d')
    M=dot(xs,xs); fv={m:f(m) for m in range(k+1)}
    D=1
    for v in fv.values(): D=lcm(D,v.denominator)
    phiD=np.vectorize(lambda m:int(fv[m]*D))(M).astype(np.int64)   # index [c,d']
    # kappa(a,b) = (1/n) sum_{c,d'} phi(c,d') (-1)^{a.d'+b.c}  -> n*D*kappa integer
    Hk=(-1)**(dot(xs,xs)%2)   # Hk[a,d']=(-1)^{a.d'}
    kapnD=(Hk@phiD.T@Hk.T).T   # [a? ] compute carefully below
    # kap[a,b] = sum_{c,d'} phi[c,d'] H[a,d'] H[b,c] = (H_b . phi . H_a^T) -> matrix [b,a] = H phi H^T
    KB=Hk@phiD@Hk.T   # KB[b,a]
    q=int((KB!=0).sum())
    # build n*D*K : (P_ab (x) P_ab)|y1,y2> = (-1)^{b.(y1^y2)} |y1^a,y2^a>
    big=np.zeros((n,n),dtype=object)
    Kint=np.zeros((n,n),dtype=np.int64)
    y1,y2=np.meshgrid(np.arange(N),np.arange(N),indexing='ij'); y1=y1.ravel(); y2=y2.ravel()
    par=lambda z: np.array([bin(t).count('1')&1 for t in range(N)])[z]
    for ai in range(N):
        for bi in range(N):
            kv=KB[bi,ai]
            if kv==0: continue
            sgn=1-2*par(np.bitwise_and(ints[bi],np.bitwise_xor(y1,y2)))
            rows=(np.bitwise_xor(y1,ints[ai]))*N+np.bitwise_xor(y2,ints[ai]); cols=y1*N+y2
            Kint[rows,cols]+=kv*sgn
    SW=np.zeros((n,n),dtype=np.int64); SW[y1*N+y2,y2*N+y1]=1
    Mt=n*D*SW-Kint
    r1=rank_mod(Mt.copy(),2147483629); r2=rank_mod(Mt.copy(),2147483587)
    beta=(-1)**(M%2)
    pBell=int((phiD!=beta*D).sum())
    print("k=%d d=%d q=%d  pBell=%d  rank mod p1=%d p2=%d %s"%(k,d,q,pBell,r1,r2,"OK" if r1==r2==pBell else "CHECK"),flush=True)
