# (U) and (A) for varied t, zeta, b: compute Smith(M(E)) directly, Per_k Smiths, C_k (k-subset sums of Thm-1 profile).
import sys, itertools
from fractions import Fraction
sys.path.insert(0,'.')
from percomp import perm, smith_vals, MOD2
from verify_am2 import v2, v2q
def toZ(fr): fr=Fraction(fr); return (fr.numerator*pow(fr.denominator,-1,MOD2))%MOD2
def Eblk(a,t,zeta,b):
    L=1<<a; i=v2(t); T=t+(1<<(a+1))*Fraction(zeta); Cp0=1
    for q in range(-L,L):
        if q: Cp0*=-q
    CT=Fraction(1)
    for q in range(-L,L): CT*=(T-q)
    B=[x for x in range(L) if (x>>i)&1==b]
    return [[toZ(CT/(Cp0*(T-(y-x)))) for y in B] for x in B]
def Mmat(E):
    n=len(E); lin=[[0]*(1<<n) for _ in range(n)]
    for x in range(n):
        for Xm in range(1,1<<n):
            low=Xm&-Xm; k=low.bit_length()-1
            lin[x][Xm]=(lin[x][Xm^low]+E[x][k])%MOD2
    out=[]
    for Sm in range(1<<n):
        row=[1]*(1<<n)
        for x in range(n):
            if Sm>>x&1: row=[r*l%MOD2 for r,l in zip(row,lin[x])]
        out.append(row)
    return out
def prof(a,i,vz): return [(a-i+1)+vz+2*k+2*sum(v2(m) for m in range(1,k+1)) for k in range(1<<(a-i-1))]*(1<<i)
if __name__=='__main__':
    a=int(sys.argv[1]); doM=sys.argv[2]=='1'; kmax=int(sys.argv[3])
    L=1<<a
    cases=[]
    for i in range(a):
        ts=[t for t in range(-L,L) if t and v2(t)==i]
        for t in (ts[0],ts[-1]):
            for zeta in (1,-3,Fraction(1,3),2,Fraction(-5,7)):
                cases.append((i,t,zeta,(t>0)*1))
    for (i,t,zeta,b) in cases:
        E=Eblk(a,t,zeta,b); n=len(E); s=prof(a,i,v2q(zeta))
        U=[0]; Aok=True
        for k in range(1,min(n,kmax)+1):
            subs=list(itertools.combinations(range(n),k))
            P=[[perm([[E[x][y] for y in Us] for x in S]) for Us in subs] for S in subs]
            sp=smith_vals(P,900); U+=sp
            sq=sorted(sum(c) for c in itertools.combinations(s,k))
            Aok&=all(p<=q for p,q in zip(sp,sq)) and len(sp)==len(sq)
        msg=f"a={a} i={i} t={t} zeta={zeta} b={b}: (A) k<={min(n,kmax)}: {Aok}"
        if doM:
            sM=smith_vals(Mmat(E),900); msg+=f"  (U): {sorted(U)==sM}"
        print(msg,flush=True)
