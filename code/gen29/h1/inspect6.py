import numpy as np
from bell import *
k=3; b=beta(k); n=64; rng=np.random.default_rng(0)
S=[0]
while len(S)<=6:
    M=charmat(k,S); e,x=rw_lad(M,b,iters=4)
    p,nz=exact_errors(M,b,x)
    if len(S)==6: break
    r=M@x-b; E=np.abs(r)>1e-7
    sc=np.abs(fwht(np.where(E,np.sign(r),0.0))); sc[S]=-1; sc=sc+1e-9*rng.random(n)
    S=S+[int(np.argmax(sc))]
print('S',S,'p',p,'coeffs',np.round(x,4))
msk=7
Q=lambda v: bin((v&msk)&(v>>3)).count('1')&1
Bf=lambda u,v:(bin((u&msk)&(v>>3)).count('1')+bin((v&msk)&(u>>3)).count('1'))&1
print('Q',[Q(s) for s in S]); print(np.array([[Bf(u,v) for v in S] for u in S]))
