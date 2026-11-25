# Smith of M restricted to rows of one residue class (all columns), in multilinear coefficient basis (degree <= #rows).
import sys, itertools
from fractions import Fraction as Fr
from lib12 import *
NB=300; MOD=1<<NB
def mult_coeffs(lins, n, maxdeg):
    """lins: list of coefficient vectors (mod MOD). returns dict monomial(frozenset)->coef of multilinearized product"""
    poly={frozenset():1}
    for l in lins:
        new={}
        for mon,cf in poly.items():
            for y in range(n):
                if l[y]==0: continue
                m2=mon|{y}
                new[m2]=(new.get(m2,0)+cf*l[y])%MOD
        poly=new
    return poly
def perclass(a,t,zeta,b,rho):
    E,c,T,B=cauchy_block(a,t,zeta,b); n=len(E); i=v2(t)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    rows=[k for k,x in enumerate(B) if x%(1<<i)==rho] if i>0 else list(range(n))
    nr=len(rows)
    mons=[frozenset(m) for d in range(nr+1) for m in itertools.combinations(range(n),d)]
    mi={m:j for j,m in enumerate(mons)}
    Mx=[]
    for d in range(nr+1):
        for S in itertools.combinations(rows,d):
            p=mult_coeffs([Ez[r] for r in S],n,d)
            v=[0]*len(mons)
            for m,cf in p.items(): v[mi[m]]=cf
            Mx.append(v)
    sm=smith_vals(Mx,NB)[:1<<nr]
    sE=smith_vals(Ez,NB)
    return sm, nr
if __name__=="__main__":
    a=int(sys.argv[1]);t=int(sys.argv[2]);zeta=Fr(sys.argv[3]);rho=int(sys.argv[4])
    sm,nr=perclass(a,t,zeta,0,rho)
    i=v2(t); s0=a-i+1+v2q(zeta)
    sk=[s0+2*k+2*v2(math.factorial(k)) for k in range(nr)]
    tgt=sorted(sum(x) for d in range(nr+1) for x in itertools.combinations(sk,d))
    print(a,t,zeta,rho,"perclass",sm,"\n target",tgt, "OK" if all(p<=q for p,q in zip(sm,tgt)) else "FAIL")
