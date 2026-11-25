# Singlet window family: phi = f(m), m = #{i: c_i=d_i=1}, f interpolates (-1)^m on window W (all windows). Exact p; q by WHT
# (rational values -> float WHT with tolerance 1e-6; values have denominators | prod).  Prints Pareto frontier (q,p).
import numpy as np, sys
from fractions import Fraction
from bell import popc, fwht
k=int(sys.argv[1]); n=4**k
z=np.arange(n); m=popc((z&((1<<k)-1)) & (z>>k))
pts=[]
for lo in range(k+1):
    for hi in range(lo,k+1):
        W=list(range(lo,hi+1))
        def f(x):
            s=Fraction(0)
            for i in W:
                t=Fraction((-1)**i)
                for j in W:
                    if j!=i: t*=Fraction(x-j,i-j)
                s+=t
            return s
        vals=np.array([float(f(x)) for x in range(k+1)])
        phi=vals[m]; p=int(sum(1 for x in range(k+1) if f(x)!=(-1)**x and True for _ in [0]) )
        p=int(np.sum([np.sum(m==x) for x in range(k+1) if f(x)!=(-1)**x]))
        q=int(np.sum(np.abs(fwht(phi))/n>1e-9))
        pts.append((q,p,lo,hi))
pts.sort()
front=[];bp=10**9
for q,p,lo,hi in pts:
    if p<bp: front.append((q,p,lo,hi)); bp=p
print('k=%d singlet frontier (q,p,window):'%k, front)
