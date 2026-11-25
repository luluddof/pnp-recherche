# Theorem S' curve for d = 3 (qutrit sites), N = 3^k. u(q) from two-row minimiser f^{(k-j,j)} (valid j <= k/4).
import numpy as np, sys
from math import lgamma, log, comb
def logf(nu):
    r=len(nu); l=[nu[i]+(r-1-i) for i in range(r)]; k=sum(nu); s=lgamma(k+1)
    for i in range(r):
        for j in range(i+1,r): s+=log(l[i]-l[j])
    for i in range(r): s-=lgamma(l[i]+1)
    return s
def first_row_dist(k):
    pr=np.zeros(k+1)
    for a in range(k, -1, -1):
        for b in range(min(a, k-a), -1, -1):
            c=k-a-b
            if c<0 or c>b: continue
            dimW=(a-b+1)*(b-c+1)*(a-c+2)/2
            pr[a]+=np.exp(logf([a,b,c])+log(dimW)-k*log(3))
    return pr
for k in [int(x) for x in sys.argv[1:]]:
    lnN=k*log(3); pr=first_row_dist(k); tail=np.concatenate([np.cumsum(pr[::-1])[::-1],[0]])
    print(f"== d=3 k={k} (sum pr = {pr.sum():.6f})")
    for eta in [0.05,0.1,0.2,0.3]:
        lq=eta*lnN; u=0
        while u+1<=k//4 and log(comb(k,u+1)*(k-2*u-1)/(k-u)) <= lq: u+=1
        b=sum(pr[a]*tail[a+2*u+1] for a in range(k+1) if a+2*u+1<=k)
        delta=-log(b)/lnN if b>0 else float('nan')
        print(f"  eta={eta}: u={u} delta<= {delta:.5f}  delta/eta^2={delta/eta**2:.3f}")
