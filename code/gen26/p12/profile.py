import sys,itertools
from fractions import Fraction as Fr
from lib12 import *
NB=200;MOD=1<<NB
a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3]) if len(sys.argv)>3 else 1
E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E)
vc=v2q(c)
sE=smith_vals([[toZ(z,MOD) for z in r] for r in E],NB)
print("a,t,i",a,t,v2(t),"v2c",vc,"Smith E",sE)
for k in range(2,n+1):
    sp=smith_vals([[toZ(z,MOD) for z in r] for r in compound(E,k,per)],NB)
    sc=sorted(sum(x) for x in itertools.combinations(sE,k))
    full=sorted(sum(x) for x in itertools.combinations_with_replacement(sE,k))[:len(sc)]
    print("k",k,"\n per ",[x-k*vc for x in sp],"\n det ",[x-k*vc for x in sc],"\n full",[x-k*vc for x in full])
