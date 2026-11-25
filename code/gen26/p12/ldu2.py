# LDU of E in order (x',rho); Per_2(E) = Per_2(L) [Delta_2 + 2 sum_j d_j^2 alpha_j beta_j^T] Per_2(U)  (exact over Q)
import sys,itertools
from fractions import Fraction as Fr
from lib12 import *
a=int(sys.argv[1]); t=int(sys.argv[2]); zeta=Fr(sys.argv[3]) if len(sys.argv)>3 else Fr(1)
E,c,T,B=cauchy_block(a,t,zeta,0); n=len(E); i=v2(t)
order=sorted(range(n),key=lambda k:(B[k]>>(i+1), B[k]%(1<<i)))
A=[[E[r][u]/c for u in order] for r in order]   # normalized by c
# LDU no pivoting
Lm=[[Fr(int(p==q)) for q in range(n)] for p in range(n)]; Um=[r[:] for r in A]; d=[None]*n
for k in range(n):
    piv=Um[k][k]; assert piv!=0
    for r in range(k+1,n):
        f=Um[r][k]/piv; Lm[r][k]=f
        Um[r]=[x-f*y for x,y in zip(Um[r],Um[k])]
for k in range(n):
    d[k]=Um[k][k]; Um[k]=[x/d[k] for x in Um[k]]
print("L integral:",all(v2q(x)>=0 for r in Lm for x in r if x),"U integral:",all(v2q(x)>=0 for r in Um for x in r if x))
print("v2 d (normalized by c):",[v2q(x) for x in d])
pairs=list(itertools.combinations(range(n),2))
P2L=[[per([[Lm[r][m] for m in M] for r in R]) for M in pairs] for R in pairs]
# solve P2L alpha = a_j
def solve(Mx,b):
    m=len(Mx); A_=[row[:]+[bb] for row,bb in zip(Mx,b)]
    for k in range(m):
        p=next(r for r in range(k,m) if A_[r][k]!=0); A_[k],A_[p]=A_[p],A_[k]
        for r in range(m):
            if r!=k and A_[r][k]!=0:
                f=A_[r][k]/A_[k][k]; A_[r]=[x-f*y for x,y in zip(A_[r],A_[k])]
    return [A_[k][m]/A_[k][k] for k in range(m)]
vd=[v2q(x) for x in d]
for j in range(min(n,4)):
    aj=[Lm[r][j]*Lm[rp][j] for r,rp in pairs]
    al=solve(P2L,aj)
    # relative excess: v(2 d_j^2 alpha_j[mm']) - (v d_m + v d_m')/2 ... print entries where correction (with beta similarly) could matter
    ex=[(pairs[q], v2q(al[q]) if al[q] else None) for q in range(len(pairs))]
    print("j",j,"v(alpha_j[mm'])  (need vs  vd_m+vd_m' - 2vd_j - 1):",[(M,v,vd[M[0]]+vd[M[1]]-2*vd[j]-1) for M,v in ex if v is not None][:12])
