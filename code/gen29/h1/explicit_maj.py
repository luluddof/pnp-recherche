# Explicit operator check of Majorana points: phi from majorana window -> Fourier coefficients (rationalised, exact
# re-check of phi on all points) -> K, rank(SWAP-K) and osr mod 2 primes.
import numpy as np, sys
from fractions import Fraction
from bell import popc, fwht, beta
from explicit_check import check
k=int(sys.argv[1]); lo=int(sys.argv[2]); hi=int(sys.argv[3]); n=4**k; msk=(1<<k)-1
A=[]
for j in range(k):
    for P in 'XY':
        a=1<<j; b=((1<<j)-1) | ((1<<j) if P=='Y' else 0); A.append(a|(b<<k))
t=0
for v in A: t^=v
A.append(t)
w=np.arange(n)
def Qv(x): return popc((x&msk)&(x>>k))&1
def Bv(x,s): return (popc((x&msk)&(s>>k)) + popc((s&msk)&(x>>k)))&1
J=np.zeros(n,np.int64)
for s in A: J+=(Qv(np.array([s]))[0] + Bv(w,s))&1
bt=1-2*Qv(w); par=int(J[0]%2); u=(J-par)//2
F={x: int(bt[u==x][0]) for x in set(u.tolist())}
W=list(range(lo,hi+1))
def f(x):
    s=Fraction(0)
    for i in W:
        tt=Fraction(F[i])
        for jj in W:
            if jj!=i: tt*=Fraction(x-jj,i-jj)
        s+=tt
    return s
fv={x:f(x) for x in range(k+1)}
phi=np.array([float(fv[x]) for x in u]); coef=fwht(phi)/n
S=[int(s) for s in np.where(np.abs(coef)>1e-9)[0]]
cf=[Fraction(float(coef[s])).limit_denominator(1<<20) for s in S]
den=1
for c in cf: den=den*c.denominator//np.gcd(den,c.denominator)
ci=[int(c*den) for c in cf]
# exact recheck of phi*den on all points
val=np.zeros(n,dtype=object)
for s,c in zip(S,ci): val=val+c*(1-2*(popc(w&s)&1)).astype(object)
p=int(sum(1 for i in range(n) if val[i]!=den*int(bt[i])))
print(f'k={k} window {W}: q={len(S)} p(exact recount)={p} den={den}',flush=True)
# check() expects integer coef2 with K = coef2/2 ...: use generic: scale SWAP by den
import explicit_check as ec
N=2**k; D=N*N
SW=np.zeros((D,D),dtype=np.int64)
for i in range(N):
    for j in range(N): SW[i*N+j, j*N+i]=1
for conv in (0,1):
    K=np.zeros((D,D),dtype=np.int64)
    for s,c in zip(S,ci):
        a=s&(N-1); b=s>>k
        P=ec.pauli(a,b,k) if conv==0 else ec.pauli(b,a,k)
        K+=c*np.kron(P,P)
    rk=ec.rank_mod(den*SW-K,1000003); rk2=ec.rank_mod(den*SW-K,998244353)
    R=K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(D,D); osr=ec.rank_mod(R,1000003)
    print(f'  conv {conv}: rank(SWAP-K) = {rk},{rk2}  osr(K) = {osr}',flush=True)
