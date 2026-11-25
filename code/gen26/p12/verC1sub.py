# Independent verification of a C1'-right violation for a Cauchy sub-block (direct Smith of M, lib12.smith_vals).
import sys,itertools
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,v2q,Mmat,smith_vals
a,t=6,4
R0=[0,2,3,9,18,33,35,40]; C0=[10,24,26,33,34,35,51,56]
E,c,T,B=cauchy_block(a,t,1,0); idx={x:k for k,x in enumerate(B)}
assert all(r in idx for r in R0+C0)
Es=[[E[idx[r]][idx[u]] for u in C0] for r in R0]
# hypothesis H for Thm U
print("H:",all(v2q(c)>=v2(r-rp) for r in R0 for rp in R0 if r!=rp))
NB=200; MOD=1<<NB
sE=smith_vals([[toZ(x,MOD) for x in r] for r in Es],NB)
sm=smith_vals([[toZ(x,MOD) for x in r] for r in Mmat(Es)],NB)
sq=sorted(sum(x) for k in range(9) for x in itertools.combinations(sE,k))
def R(s,p): return sum(1 for x in s if x<p)
print("Smith E'",sE)
print("deficits (p, R_p(M), sqfree count):",[(p,R(sm,p),R(sq,p)) for p in range(max(sq)+2) if R(sm,p)<R(sq,p)])
print("sorted Smith(M) vs sqfree where M>sq:",[(j,x,y) for j,(x,y) in enumerate(zip(sm,sq)) if x>y][:10])
# also left side: full monomials <= Smith(M)?
full=sorted(sum(x) for k in range(9) for x in itertools.combinations_with_replacement(sE,k))[:256]
print("C1'-left (full<=M):",all(x<=y for x,y in zip(full,sm)))
