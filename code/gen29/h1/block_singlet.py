# Block generalisation of the singlet family: n = k/b blocks, y_B = IP_b(block) in {0,1}, Pr = rho = (1-2^-b)/2,
# centred-variable cost 4^b - 1.  phi = f(m), m = sum y_B, f interpolates (-1)^m on a window W of d+1 integers.
# q <= sum_{j<=d} C(n,j)(4^b-1)^j ; error fraction = Pr[m notin W], m ~ Bin(n,rho). eta = log2 q / k, delta = -log2(err)/k.
import math, sys
from math import lgamma, log
def lbin(n,j): return lgamma(n+1)-lgamma(j+1)-lgamma(n-j+1)
def lse(xs):
    m=max(xs); return m+math.log(sum(math.exp(x-m) for x in xs))
def curve(k,b):
    n=k//b; kk=n*b; rho=(1-2.0**-b)/2; c=4**b-1
    lp=[lbin(n,j)+j*log(rho)+(n-j)*log(1-rho) for j in range(n+1)]
    lq=[lbin(n,j)+j*log(c) for j in range(n+1)]
    out=[]
    for d in range(0,n):
        q=lse(lq[:d+1])
        # best window of d+1 consecutive integers maximising mass
        best=None
        for s in range(0,n-d+1):
            inside=lse(lp[s:s+d+1])
            outside=[x for i,x in enumerate(lp) if i<s or i>s+d]
            if not outside: continue
            e=lse(outside)
            if best is None or e<best: best=e
        if best is None: break
        out.append((q/log(2)/kk, -best/log(2)/kk))
    return out
k=int(sys.argv[1]) if len(sys.argv)>1 else 600
for target in [0.05,0.1,0.2,0.3,0.5]:
    row=[]
    for b in [1,2,3,4,5,6]:
        cv=curve(k,b)
        # interpolate delta at eta=target (take best delta with eta<=target)
        ds=[dl for (et,dl) in cv if et<=target]
        row.append(max(ds) if ds else float('nan'))
    print(f"k={k} eta<={target}: delta by b=1..6:", " ".join(f"{x:.5f}" for x in row), " delta/eta^2:", " ".join(f"{x/target**2:.3f}" for x in row))
