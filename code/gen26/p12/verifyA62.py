# Independent re-check of the (A) violation at a=6, t=2, k=2.
import sys,itertools,time
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,v2q,smith_vals
from fastA import smith_fast
a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3]); NB=int(sys.argv[4]); b=int(sys.argv[5]) if len(sys.argv)>5 else 0
MOD=1<<NB
E,c,T,B=cauchy_block(a,t,zeta,b); n=len(E)
# entries exactly as Fractions -> mod
Ez=[[toZ(x,MOD) for x in r] for r in E]
pairs=list(itertools.combinations(range(n),2))
P=[[ (Ez[r][u]*Ez[rp][up]+Ez[r][up]*Ez[rp][u])%MOD for (u,up) in pairs] for (r,rp) in pairs]
t0=time.time()
sE=smith_vals(Ez,NB)
sp=smith_vals(P,NB)
tgt=sorted(x+y for x,y in itertools.combinations(sE,2))
viol=[(j,p,q) for j,(p,q) in enumerate(zip(sp,tgt)) if p>q]
print("Smith E",sE)
print("violations (pos, Per2, C2):",viol[:20], "count",len(viol), f"{time.time()-t0:.0f}s")
def R(s,p): return sum(1 for x in s if x<p)
bad=[(p,R(sp,p),R(tgt,p)) for p in range(max(tgt)+2) if R(sp,p)<R(tgt,p)]
print("p with R_p(Per2) < R_p(C2):",bad[:20])
