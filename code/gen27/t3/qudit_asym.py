# qudit-site singlet(antisym)-number family: SWAP_{n^k} = prod_i (1-2A_i) = (-1)^M, K = f(M), deg f = d = tau k
# q <= sum_{u<=d} C(k,u)(n^2-1)^u ; p = dim of M-eigenspaces outside window; pi = Pr[A] = (n-1)/(2n)
from math import log2
from scipy.optimize import brentq
def H(t): return 0 if t<=0 or t>=1 else -t*log2(t)-(1-t)*log2(1-t)
def D(a,p):
    if a<=0: return log2(1/(1-p))
    return a*log2(a/p)+(1-a)*log2((1-a)/(1-p))
def delta_bits(tau,p):
    f=lambda a: D(a,p)-D(a+tau,p)
    lo=max(0.0,p-tau); hi=p
    a=lo if f(lo)<0 else brentq(f,lo,hi)
    return min(D(a,p),D(a+tau,p))
for n in [2,3,5,8]:
    pi=(n-1)/(2*n); L=log2(n)
    out=[]
    for tau in [0.01,0.03,0.1]:
        eta=(H(tau)+tau*log2(n*n-1))/L; d=delta_bits(tau,pi)/L
        out.append("eta=%.3f delta=%.5f d/eta^2=%.4f"%(eta,d,d/eta**2))
    print("n=%d:"%n," | ".join(out))
