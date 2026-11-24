import math, numpy as np
K=np.arange(1,20000); S=np.arange(-30,32)
bs=np.exp(np.linspace(math.log(60),math.log(0.004),800))
LP=np.array([-np.log1p(-np.exp(-b*K)).sum() for b in bs])
lZM=np.array([math.log(np.exp(-b*(S**2)).sum()) for b in bs])+LP
def s(x): return float(np.min(lZM+bs*x))
for u in [0.003,0.01,0.03,0.06,0.1,0.116,0.2]:
    for v in [0.0,0.03]:
        cost=2*s(u+v); lo=min(s(u),s(2*u+2*v)/2)
        print(u,v,"sM(u)=%.4f approx %.4f  sM(2u+2v)/2=%.4f  ratio=%.4f"%(s(u),u*math.log(3*math.e/u),s(2*u+2*v)/2,cost/lo))
