from math import log2
from scipy.optimize import brentq
def H(t): return 0 if t<=0 or t>=1 else -t*log2(t)-(1-t)*log2(1-t)
def D(a,p=0.25):
    if a<=0: return log2(1/(1-p))
    if a>=1: return log2(1/p)
    return a*log2(a/p)+(1-a)*log2((1-a)/(1-p))
def delta(tau):
    # best window [a,a+tau]: maximize min tail exponent
    if tau>=0.75: return D(tau)  # window from 0
    f=lambda a: D(a)-D(a+tau)
    lo=max(0.0,0.25-tau); hi=0.25
    if f(lo)<0: a=lo
    else: a=brentq(f,lo,hi)
    return min(D(a),D(a+tau))
print(" tau   eta      delta    eta^2/4   delta/eta  delta/eta^2")
for tau in [0.001,0.003,0.01,0.03,0.05,0.1,0.15,0.2,0.3,0.4,0.5]:
    eta=H(tau)+tau*log2(3); d=delta(tau)
    print("%.3f  %.4f  %.6f  %.6f  %.4f  %.4f"%(tau,eta,d,eta*eta/4,d/eta,d/eta**2))
