import sys,itertools
from fractions import Fraction as Fr
from lib12 import *
NB=300;MOD=1<<NB
a=int(sys.argv[1]); t=int(sys.argv[2]); k=int(sys.argv[3]); zeta=Fr(sys.argv[4]) if len(sys.argv)>4 else 1
E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E); vc=v2q(c); i=v2(t)
Ez=[[toZ(z,MOD) for z in r] for r in E]
sE=[x-vc+i for x in smith_vals(Ez,NB)]
# permanent mod MOD via DP over integers
def perZ(Mx):
    m=len(Mx); dp={0:1}
    for r in range(m):
        nd={}
        for mask,val in dp.items():
            for j in range(m):
                if not mask>>j&1:
                    nm=mask|1<<j; nd[nm]=(nd.get(nm,0)+val*Mx[r][j])%MOD
        dp=nd
    return dp[(1<<m)-1]
subs=list(itertools.combinations(range(n),k))
P=[[perZ([[Ez[r][u] for u in U] for r in R]) for U in subs] for R in subs]
sp=[x-k*(vc-i) for x in smith_vals(P,NB)]
sc=sorted(sum(x) for x in itertools.combinations(sE,k))
print("a t i k",a,t,i,k,"modes(normalized, -v2c+i)",sE)
print(" per",sp); print(" det",sc); print(" slack",[q-p for p,q in zip(sp,sc)])
