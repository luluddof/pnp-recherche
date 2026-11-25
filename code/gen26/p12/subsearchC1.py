# Search Cauchy sub-blocks of crossing blocks for C1'-right (union over k, = Smith M by Thm U) violations; also count (A) violations.
import sys,itertools,random,time
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2
from fastA import smith_fast
from levelm3 import perZ
NB=260; MOD=1<<NB
a=int(sys.argv[1]); t=int(sys.argv[2]); nsub=int(sys.argv[3]); random.seed(int(sys.argv[4])); tl=float(sys.argv[5])
E,c,T,B=cauchy_block(a,t,1,0); n=len(E); Ez=[[toZ(x,MOD) for x in r] for r in E]
t0=time.time(); tried=Aviol=C1viol=0
def R(s,p): return sum(1 for x in s if x<p)
while time.time()-t0<tl:
    R0=sorted(random.sample(range(n),nsub)); C0=sorted(random.sample(range(n),nsub)) if random.random()<0.5 else R0
    Es=[[Ez[r][u] for u in C0] for r in R0]; sE=smith_fast(Es); tried+=1
    un=[0]; sq=[0]; av=False
    for k in range(1,nsub+1):
        subs=list(itertools.combinations(range(nsub),k))
        sp=smith_fast([[perZ([[Es[r][u] for u in U] for r in Rr]) for U in subs] for Rr in subs])
        sc=sorted(sum(x) for x in itertools.combinations(sE,k))
        if any(p>q for p,q in zip(sp,sc)): av=True
        un+=sp; sq+=sc
    un.sort(); sq.sort(); Aviol+=av
    if any(p>q for p,q in zip(un,sq)):
        C1viol+=1; print("C1' VIOLATION rows",[B[r] for r in R0],"cols",[B[u] for u in C0],"SmithE'",sE,flush=True)
print(f"a={a} t={t} nsub={nsub}: tried {tried}, (A) violated in {Aviol}, C1'-right violated in {C1viol}")
