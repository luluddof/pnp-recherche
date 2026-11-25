# Explore the core (a,0): alternant candidates A_T = sum_S det(Binv[T,S]) prod_{x in S} l_x ; report min v2 on cube vs w(T)
import sys, itertools
from fractions import Fraction
sys.path.insert(0,'.')
from verify_am2 import v2, v2q, toZ, MOD, K
def Ecore(a,t,zeta,b=0):
    L=1<<a; i=v2(t)
    T=t+(1<<(a+1))*Fraction(zeta); Cp0=1
    for q in range(-L,L):
        if q: Cp0*=-q
    CT=Fraction(1)
    for q in range(-L,L): CT*=(T-q)
    B=[x for x in range(L) if (x>>i)&1==b]
    return [[toZ(CT/(Cp0*(T-(y-x)))) for y in B] for x in B], v2q(CT/Cp0)-i
def lins(E):
    n=len(E); lin=[[0]*(1<<n) for _ in range(n)]
    for x in range(n):
        for Xm in range(1,1<<n):
            low=Xm&-Xm; k=low.bit_length()-1
            lin[x][Xm]=(lin[x][Xm^low]+E[x][k])%MOD
    return lin
from math import comb
def det(M):
    M=[[Fraction(v) for v in r] for r in M]; n=len(M); d=Fraction(1)
    for c in range(n):
        p=next((r for r in range(c,n) if M[r][c]),None)
        if p is None: return 0
        if p!=c: M[c],M[p]=M[p],M[c]; d=-d
        d*=M[c][c]
        for r in range(c+1,n):
            f=M[r][c]/M[c][c]
            for k in range(c,n): M[r][k]-=f*M[c][k]
    return int(d)
if __name__=='__main__':
    a=int(sys.argv[1]); E,w0=Ecore(a,1,1); n=len(E); lin=lins(E)
    h=lambda j: j+sum(v2(m) for m in range(1,j+1))
    Binv=[[(-1)**(j-x)*comb(j,x) for x in range(n)] for j in range(n)]
    for k in range(1,n+1):
        for T in itertools.combinations(range(n),k):
            f=[0]*(1<<n)
            for S in itertools.combinations(range(n),k):
                d=det([[Binv[j][x] for x in S] for j in T])
                if d==0: continue
                pr=[1]*(1<<n)
                for x in S: pr=[p*q%MOD for p,q in zip(pr,lin[x])]
                f=[(F+d*p)%MOD for F,p in zip(f,pr)]
            mv=min(v2(z) for z in f)
            print(f"T={T}: min v2={mv}  target w(T)={sum(w0+2*h(j) for j in T)}  rowdiff-only={sum(w0+h(j) for j in T)}")
