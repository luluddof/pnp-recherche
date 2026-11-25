# Independent check of P8 Theorem A2 (level i=a-2): (1) factor claims symbolically for a<=8; (2) direct Smith(M) vs bound for a<=4.
import sys, itertools
from fractions import Fraction as Fr
from lib12 import *
K=400; MOD=1<<K
def vv(z): z%=MOD; return K if z==0 else v2(z)
def factors(a,t,zeta,b):
    E,c,T,B=cauchy_block(a,t,zeta,b); i=a-2; n=len(E)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    w0=3+v2q(zeta); u=toZ(T/2**i,MOD)
    idx={x:k for k,x in enumerate(B)}
    row=lambda r,j: idx[r+b*(1<<i)+(j<<(i+1))]
    ok=True
    for r in range(1<<i):
        l0=Ez[row(r,0)]; l1=Ez[row(r,1)]
        nu={row(r,0),row(r,1)}  # columns of class r (cols indexed same as rows)
        rp=(r^(1<<(i-1))) if i>=1 else None
        nup={row(rp,0),row(rp,1)} if i>=1 else set()
        # F1
        if min(map(vv,l0))!=w0 or {y for y in range(n) if vv(l0[y])==w0}!=nu: ok=False; print(" F1 bad",r)
        F2=[(u*p-(u+2)*q)%MOD for p,q in zip(l0,l1)]
        pred={row(r,1)}|nup
        if min(map(vv,F2))!=w0+2 or {y for y in range(n) if vv(F2[y])==w0+2}!=pred: ok=False; print(" F2 bad",r)
        # F3 multilinear coefficients
        w3=2*w0+1; quad_odd=set(); 
        for y in range(n):
            for yp in range(y+1,n):
                co=(l0[y]*l1[yp]+l0[yp]*l1[y])%MOD
                if vv(co)<w3: ok=False; print(" F3 quad not divisible",r,y,yp)
                elif vv(co)==w3: quad_odd.add((y,yp))
            co=(l0[y]*l1[y]+(l0[y]<<w0))%MOD
            if vv(co)<w3: ok=False; print(" F3 lin not divisible",r,y)
        if quad_odd!={tuple(sorted((row(r,0),row(r,1))))}: ok=False; print(" F3 quad reduction",quad_odd)
    return ok
def tensor_sorted(lst,reps):
    out=[0]
    for _ in range(reps): out=sorted(p+q for p in out for q in lst)
    return out
def direct(a,t,zeta,b):
    E,c,T,B=cauchy_block(a,t,zeta,b); i=a-2; w0=3+v2q(zeta)
    NB=200; M=Mmat(E); sm=smith_vals([[toZ(z,1<<NB) for z in r] for r in M],NB)
    bound=tensor_sorted([0,w0,w0+2,2*w0+1],1<<i); sq=tensor_sorted([0,w0,w0+2,2*w0+2],1<<i)
    return all(p<=q for p,q in zip(sm,bound)), sm==bound, all(p<=q for p,q in zip(sm,sq))
mode=sys.argv[1]
tot=0;good=0
for a in ([3,4,5,6,7,8] if mode=='sym' else [2,3]):
    L=1<<a; i=a-2
    ts=[t for t in range(-L,L) if t!=0 and v2(t)==i]
    for t in ts[:8]:
        for zeta in (1,-3,Fr(1,3),2,-12,Fr(5,7),4):
            for b in (0,1):
                if mode=='sym': r=factors(a,t,zeta,b); tot+=1; good+=r
                else:
                    r=direct(a,t,zeta,b); tot+=1; good+=r[0]
                    print(a,t,zeta,b,"bound<=",r[0],"equal",r[1],"sqfree",r[2],flush=True)
print(mode,"ok",good,"/",tot)
