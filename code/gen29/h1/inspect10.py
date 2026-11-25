import numpy as np
from bell import *
k=4; b=beta(k); n=4**k
S=[0]; rng=np.random.default_rng(0)
while len(S)<=19:
    M=charmat(k,S); e,x=rw_lad(M,b,iters=4)
    r=M@x-b; E=np.abs(r)>1e-7
    sc=np.abs(fwht(np.where(E,np.sign(r),0.0))); sc[S]=-1; sc=sc+1e-9*rng.random(n)
    S=S+[int(t) for t in np.argsort(-sc)[:2]]
M=charmat(k,S[:-2]); e,x=rw_lad(M,b,iters=4)
print('errors',e)
def pauli(s):
    out=''
    for i in range(k):
        a=(s>>i)&1; bb=(s>>(i+k))&1
        out+='IXZY'[a+2*bb]
    return out
for s,v in zip(S[:-2],x):
    if abs(v)>1e-9: print(pauli(s), round(v,6))
r=M@x-b; E=np.where(np.abs(r)>1e-7)[0]
z=E; c=z&15; d=z>>4
print('error points (c,d):',[(int(a),int(bb)) for a,bb in zip(c,d)][:60])
