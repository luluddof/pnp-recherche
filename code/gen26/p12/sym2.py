# Is 2K (rows (2 E_ru E_ru')_{u<u'}) inside rowspace(Per_2 E)?  Compare Smith(Per_2) with Smith([Per_2;2K]).
import sys,itertools
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,v2q
from fastA import smith_fast, NB, MOD
from levelm3 import perZ
def run(a,t,zeta):
    E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E)
    Ez=[[toZ(z,MOD) for z in r] for r in E]
    pairs=list(itertools.combinations(range(n),2))
    P=[[ (Ez[r][u]*Ez[rp][up]+Ez[r][up]*Ez[rp][u])%MOD for (u,up) in pairs] for (r,rp) in pairs]
    K2=[[ (2*Ez[r][u]*Ez[r][up])%MOD for (u,up) in pairs] for r in range(n)]
    s1=smith_fast(P); s2=smith_fast(P+K2)[:len(pairs)]
    sE=smith_fast(Ez); tgt=sorted(x+y for x,y in itertools.combinations(sE,2))
    print(a,t,zeta,"Per2==Sym2:",s1==s2,"| Sym2<=tgt:",all(p<=q for p,q in zip(s2,tgt)),"| Per2<=tgt:",all(p<=q for p,q in zip(s1,tgt)))
    if s1!=s2: print("   Per2",s1[:16],"\n   Sym2",s2[:16])
import sys

for (a,t,z) in [(2,1,1),(2,2,1),(3,1,1),(3,2,1),(3,4,1),(4,1,1),(4,2,1),(4,4,1),(4,8,1),(3,3,Fr(1,3))]:
    run(a,t,Fr(z))
