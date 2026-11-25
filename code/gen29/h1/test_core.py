import numpy as np, itertools, sys
from bell import *
# check: singlet window family k=3 (N=8): full ball weight<=1 -> q=1+3*3=10
k=3; b=beta(k)
def sitew(s,k): return sum(1 for i in range(k) if ((s>>i)&1) or ((s>>(i+k))&1))
for w in range(k+1):
    S=[s for s in range(4**k) if sitew(s,k)<=w]
    M=charmat(k,S); e,x=rw_lad(M,b)
    print(k,w,len(S),'rwLAD errors',e,'exact',exact_errors(M,b,x))
