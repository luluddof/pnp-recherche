# Check (A) at level a-2 with exact profile: Smith(Per_k) == {sum_rho w(T_rho): sum |T_rho| = k}, w in {0:0,{0}:w0,{1}:w0+2,{01}:2w0+1}
import sys,itertools
from fractions import Fraction as Fr
from lib12 import *
NB=200;MOD=1<<NB
def run(a,t,zeta,b):
    E,c,T,B=cauchy_block(a,t,zeta,b); n=len(E); i=a-2; w0=3+v2q(zeta)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    ok=True
    for k in range(1,n+1):
        sp=smith_vals([[toZ(z,MOD) for z in r] for r in compound(E,k,per)],NB)
        prof=sorted(sum(y[1] for y in x) for x in itertools.product([(0,0),(1,w0),(1,w0+2),(2,2*w0+1)],repeat=1<<i) if sum(y[0] for y in x)==k)
        sc=sorted(sum(x) for x in itertools.combinations(smith_vals(Ez,NB),k))
        eq=(sp==prof); A=all(p<=q for p,q in zip(sp,sc))
        if not (eq and A): ok=False; print(" k",k,sp,prof,sc)
    print(a,t,zeta,b,"Per_k==graded profile & (A):",ok,flush=True); return ok
r=[]
for a in (2,3,4):
    L=1<<a; i=a-2
    for t in [t for t in range(-L,L) if t and v2(t)==i][:3]:
        for zeta in ([1,Fr(-1,3),2] if a<4 else [1]):
            r.append(run(a,t,zeta,0))
print(sum(r),"/",len(r))
