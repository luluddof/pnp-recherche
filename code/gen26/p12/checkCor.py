# P8 Corollary: sub-block E_sub (x',y' in {0,1} in each class) at ANY level i<=a-2: same factor claims as Thm A2 with w0=s0=a-i+1+v2 zeta.
import sys, itertools
from fractions import Fraction as Fr
from lib12 import *
K=400; MOD=1<<K
def vv(z): z%=MOD; return K if z==0 else v2(z)
def factors_sub(a,t,zeta,b):
    E,c,T,B=cauchy_block(a,t,zeta,b); i=v2(t)
    sub=[k for k,x in enumerate(B) if (x>>(i+1))<2]
    Bs=[B[k] for k in sub]; n=len(sub)
    Ez=[[toZ(E[r][u],MOD) for u in sub] for r in sub]
    w0=a-i+1+v2q(zeta); u=toZ(T/2**i,MOD)
    idx={x:k for k,x in enumerate(Bs)}
    row=lambda r,j: idx[r+b*(1<<i)+(j<<(i+1))]
    ok=True
    for r in range(1<<i):
        l0=Ez[row(r,0)]; l1=Ez[row(r,1)]
        nu={row(r,0),row(r,1)}
        nup={row(r^(1<<(i-1)),0),row(r^(1<<(i-1)),1)} if i>=1 else set()
        if min(map(vv,l0))!=w0 or {y for y in range(n) if vv(l0[y])==w0}!=nu: ok=False
        F2=[(u*p-(u+2)*q)%MOD for p,q in zip(l0,l1)]
        if min(map(vv,F2))!=w0+2 or {y for y in range(n) if vv(F2[y])==w0+2}!=({row(r,1)}|nup): ok=False
        w3=2*w0+1; qo=set()
        for y in range(n):
            for yp in range(y+1,n):
                co=(l0[y]*l1[yp]+l0[yp]*l1[y])%MOD
                if vv(co)<w3: ok=False
                elif vv(co)==w3: qo.add((y,yp))
            if vv((l0[y]*l1[y]+(l0[y]<<w0))%MOD)<w3: ok=False
        if qo!={tuple(sorted((row(r,0),row(r,1))))}: ok=False
    return ok
tot=good=0
for a in (4,5,6,7):
    L=1<<a
    for i in range(0,a-2):
        ts=[t for t in range(-L,L) if t and v2(t)==i][:3]
        for t in ts:
            for z in (1,Fr(-1,3),2):
                tot+=1; good+=factors_sub(a,t,z,0)
print("corollary factor claims ok",good,"/",tot)
