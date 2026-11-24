import math, numpy as np
from chars import *
from scaling import run, mk_phi
from modscale import sl2, entropy
N=300
for k,j2 in [(1,1),(2,2),(4,4),(8,8),(8,4)]:
    eps,r=run("sl2 k=%d M=L(%g), own profile"%(k,j2/2),sl2(N,k,j2))
    i=int(np.argmin(r)); print("     -> min r=%.4f at eps=%.3g"%(r[i],eps[i]))
# psi scan for L(1/2) level 1
phi,mean=mk_phi(sl2(N,1,1),1)
print([ (b, round(phi(b)-b*mean(b)[0],4)) for b in [0.3,0.5,0.7,0.83,1,1.5,2,4,8]])
# does V dominate? compare s_V, s_M
sV,_=entropy(sl2(N,1,0)); sM,_=entropy(sl2(N,1,1))
print("eps, s_V, s_M:",[(e,round(float(sV(e)),3),round(float(sM(e)),3)) for e in [0.1,0.5,1,2,5,10]])
