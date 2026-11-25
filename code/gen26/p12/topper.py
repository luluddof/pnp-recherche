# v2(per E) - n*(v2 c - i)  vs  v2(det E) - n*(v2 c - i)  (normalized), Ryser formula mod 2^K
import sys,itertools
from fractions import Fraction as Fr
from lib12 import cauchy_block,toZ,v2,v2q
K=200;MOD=1<<K
def ryser(A):
    n=len(A); tot=0
    # Gray-code Ryser
    rowsum=[0]*n; sgn=1 if n%2==0 else -1
    prev=0
    for g in range(1,1<<n):
        gc=g^(g>>1); diff=gc^prev; j=diff.bit_length()-1; prev=gc
        if gc>>j&1: 
            for r in range(n): rowsum[r]+=A[r][j]
        else:
            for r in range(n): rowsum[r]-=A[r][j]
        p=1
        for r in range(n): p=p*rowsum[r]%MOD
        k=bin(gc).count('1')
        tot+= p if (n-k)%2==0 else -p
    return tot%MOD
for a in (2,3,4,5):
    L=1<<a
    for i in range(0,a):
        t=1<<i
        E,c,T,B=cauchy_block(a,t,1,0); n=len(E)
        if n>16: continue
        Ez=[[toZ(z,MOD) for z in r] for r in E]
        vp=v2(ryser(Ez)); norm=n*(v2q(c)-i)
        D=1<<(a-i); import math
        vdet=(1<<i)*sum(2*j+2*v2(math.factorial(j)) for j in range(D//2))
        print(f"a={a} i={i} n={n} D={D}: v(per)-norm={vp-norm}  v(det)-norm={vdet}",flush=True)
