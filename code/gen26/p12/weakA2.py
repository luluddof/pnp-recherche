# Check components of weak (A)_2: Smith(Sym2)<=pair sums; Smith(2K) >= Smith(E)+3+v2(zeta); and the final count inequality.
import sys,itertools
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,v2q
from fastA import smith_fast, NB, MOD
def R(sm,p): return sum(1 for x in sm if x<p)
for (a,t,z) in [(3,1,1),(3,2,1),(3,3,Fr(1,3)),(4,1,1),(4,2,1),(4,4,2),(4,6,Fr(-5,3)),(5,2,1),(5,4,1)]:
    z=Fr(z); E,c,T,B=cauchy_block(a,t,z,0); n=len(E)
    Ez=[[toZ(x,MOD) for x in r] for r in E]
    pairs=list(itertools.combinations(range(n),2))
    P=[[ (Ez[r][u]*Ez[rp][up]+Ez[r][up]*Ez[rp][u])%MOD for (u,up) in pairs] for (r,rp) in pairs]
    K2=[[ (2*Ez[r][u]*Ez[r][up])%MOD for (u,up) in pairs] for r in range(n)]
    sE=smith_fast(Ez); sP=smith_fast(P); sS=smith_fast(P+K2)[:len(pairs)]; sK=smith_fast(K2)
    tgt=sorted(x+y for x,y in itertools.combinations(sE,2)); vz=v2q(z)
    c1=all(p<=q for p,q in zip(sS,tgt)); c2=all(p>=q+3+vz for p,q in zip(sK,sE))
    c3=all(R(sP,p)>=R(tgt,p)-R(sE,p-3-vz) for p in range(0,max(tgt)+2))
    print(a,t,z,"Sym2<=tgt",c1,"Smith(2K)>=Smith(E)+3+v(zeta)",c2,"weakA2",c3, "minSmith(2K)-s0",sK[0]-sE[0])
