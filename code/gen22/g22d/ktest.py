import numpy as np, sys
from d2 import *
N,p,q=3,2,2; n=9; F=swap(N); Fm=F.reshape(n,n)
e=lambda i,j: np.eye(n)[i*N+j]
cands={'01,10':[e(0,1),e(1,0)], '00,11':[e(0,0),e(1,1)], '01-10,00':[(e(0,1)-e(1,0))/2**.5,e(0,0)],
 '00,01':[e(0,0),e(0,1)], '00,10':[e(0,0),e(1,0)], '01-10,22':[(e(0,1)-e(1,0))/2**.5,e(2,2)], '01,02':[e(0,1),e(0,2)],
 '00+11, 01':[ (e(0,0)+e(1,1))/2**.5, e(0,1)], '01,12':[e(0,1),e(1,2)]}
rng=np.random.default_rng(0)
for name,vs in cands.items():
    K=np.array(vs).T; PK=K@K.conj().T
    T=(Fm@(np.eye(n)-PK)).reshape(N,N,N,N).astype(complex)
    best=min(solve(T,N,p,q,rng.standard_normal(4*(p+q)*N*N)*0.7,2500)[0] for r in range(4))
    print(name, 'dist^2(F(1-P_K), J(2,2)) =', round(best,6), flush=True)
