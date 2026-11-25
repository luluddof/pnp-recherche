# C1'-right (union level) for given sub-blocks; Smith(M) computed directly (M 2^n' x 2^n') and via union Per_k (Thm U), both.
import sys,itertools
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,Mmat
from fastA import smith_fast
from levelm3 import perZ
NB=260; MOD=1<<NB
E,c,T,B=cauchy_block(6,2,1,0); idx={x:k for k,x in enumerate(B)}
cases=[([1,24,25,37,48,49,52,56],[5,13,16,20,36,45,49,56]),([0,1,24,33,40,41,57,60],[0,13,28,29,32,45,48,53]),
       ([5,20,25,32,36,44,45,57],[5,20,25,32,36,44,45,57])]
for R0,C0 in cases:
    Es=[[E[idx[r]][idx[u]] for u in C0] for r in R0]; ns=len(R0)
    Ez=[[toZ(x,MOD) for x in r] for r in Es]
    sE=smith_fast(Ez)
    un=[0]
    for k in range(1,ns+1):
        subs=list(itertools.combinations(range(ns),k))
        un+=smith_fast([[perZ([[Ez[r][u] for u in U] for r in Rr]) for U in subs] for Rr in subs])
    un.sort()
    M=Mmat(Es); sm=smith_fast([[toZ(x,MOD) for x in r] for r in M])
    sq=sorted(sum(x) for k in range(ns+1) for x in itertools.combinations(sE,k))
    def R(s,p): return sum(1 for x in s if x<p)
    bad=[(p,R(sm,p),R(sq,p)) for p in range(max(sq)+2) if R(sm,p)<R(sq,p)]
    print("U check (Smith M == union Per_k):",sm==un,"| C1'-right holds:",not bad,"| deficits (p,R_p(M),sqfree):",bad[:8],flush=True)
