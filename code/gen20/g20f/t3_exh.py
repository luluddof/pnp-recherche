# Theorem 3 exhaustive: min over G (rank normal form) and all G' of rank(F - G(x)G') over GF(p).
import numpy as np, itertools, sys
from brank import batch_rank
def run(N,p):
    F=np.zeros((N*N,N*N),dtype=np.int64)
    for x in range(N):
        for y in range(N): F[x*N+y, y*N+x]=1    # rows (x,y), cols (u,v): [x=v][y=u]
    best=None
    allGp=np.array(list(itertools.product(range(p),repeat=N*N)),dtype=np.int64).reshape(-1,N,N)
    for rk in range(N+1):
        G=np.zeros((N,N),dtype=np.int64); G[range(rk),range(rk)]=1
        res=[]
        for s in range(0,len(allGp),20000):
            Gp=allGp[s:s+20000]
            K=np.einsum('xu,byv->bxyuv',G,Gp).reshape(-1,N*N,N*N)
            res.append(batch_rank(F[None]-K,p))
        res=np.concatenate(res); mn=res.min()
        print(f'N={N} p={p} rankG={rk}: min rank = {mn}', flush=True)
        best=mn if best is None else min(best,mn)
    print(f'== N={N} GF({p}): min p at q=1 = {best}; N(N-1)/2 = {N*(N-1)//2}',flush=True)
for N,p in [(2,2),(2,3),(2,5),(3,2),(3,3),(4,2)]: run(N,p)
