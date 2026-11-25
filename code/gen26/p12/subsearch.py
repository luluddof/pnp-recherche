# Search for small Cauchy sub-blocks E'=E[R0,C0] of crossing blocks violating (A) for some k; then test C1'-right (union) on them.
import sys,itertools,random,time
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,v2q
from fastA import smith_fast
from levelm3 import perZ
NB=260; MOD=1<<NB
random.seed(int(sys.argv[4]) if len(sys.argv)>4 else 1)
a=int(sys.argv[1]); t=int(sys.argv[2]); nsub=int(sys.argv[3])
E,c,T,B=cauchy_block(a,t,1,0); n=len(E)
Ez=[[toZ(x,MOD) for x in r] for r in E]
t0=time.time(); found=0; tried=0
while time.time()-t0<float(sys.argv[5]) if len(sys.argv)>5 else 200:
    R0=sorted(random.sample(range(n),nsub)); C0=sorted(random.sample(range(n),nsub)) if random.random()<0.5 else R0
    Es=[[Ez[r][u] for u in C0] for r in R0]
    sE=smith_fast(Es); tried+=1
    viol=[]
    for k in range(2,nsub+1):
        subs=list(itertools.combinations(range(nsub),k))
        sp=smith_fast([[perZ([[Es[r][u] for u in U] for r in Rr]) for U in subs] for Rr in subs])
        sc=sorted(sum(x) for x in itertools.combinations(sE,k))
        bad=[(j,p,q) for j,(p,q) in enumerate(zip(sp,sc)) if p>q]
        if bad: viol.append((k,bad[:3]))
    if viol:
        found+=1
        print("VIOLATION rows",[B[r] for r in R0],"cols",[B[u] for u in C0],"SmithE'",sE,"viol",viol,flush=True)
        if found>=3: break
print("tried",tried,"found",found)
