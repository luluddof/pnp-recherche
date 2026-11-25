import sys,itertools
from fractions import Fraction as Fr
from lib12 import *
from levelm3 import perZ
NB=300;MOD=1<<NB
def perprofile(E):
    Ez=[[toZ(z,MOD) for z in r] for r in E]; n=len(E); out={0:[0]}
    for k in range(1,n+1):
        subs=list(itertools.combinations(range(n),k))
        out[k]=smith_vals([[perZ([[Ez[r][u] for u in U] for r in R]) for U in subs] for R in subs],NB)
    return out
a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3])
E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E)
E0,c0,_,_=cauchy_block(a-1,t//2,zeta,0); sh=v2q(c)-1-v2q(c0)
big=perprofile(E); small=perprofile(E0)
Mb=sorted(x for k in big for x in big[k])
Ms=sorted(x+y+(k0+k1)*sh for k0 in small for k1 in small for x in small[k0] for y in small[k1])
sq=sorted(sum(x) for k in range(n+1) for x in itertools.combinations(smith_vals([[toZ(z,MOD) for z in r] for r in E],NB),k))
print(a,t,zeta,"M-level TD <=:",all(p<=q for p,q in zip(Mb,Ms)),">=:",all(p>=q for p,q in zip(Mb,Ms)),"C1'right:",all(p<=q for p,q in zip(Mb,sq)))
d=[(j,p,q) for j,(p,q) in enumerate(zip(Mb,Ms)) if p!=q]; print("  diffs (pos,big,tens):",d[:12])
