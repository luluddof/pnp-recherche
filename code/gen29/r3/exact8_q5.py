# q=5 at N=8: supports {0} u T (translation), T up to O(Q) (reflections); batched vertex enumeration.
import itertools, time, numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
k=3; n2=64; msk=7; q=5
POP=np.array([bin(i).count('1') for i in range(256)])
Q=lambda x: POP[(x&msk)&(x>>k)]&1
B=lambda x,v: (POP[(x&msk)&(v>>k)]+POP[(x>>k)&(v&msk)])&1
W=np.arange(n2); beta=1-2*Q(W)
combs=np.array(list(itertools.combinations(range(1,n2),q-1)),dtype=np.int64)
def mk(c):
    m=np.zeros(len(c),dtype=np.uint64)
    for j in range(c.shape[1]): m|=(np.uint64(1)<<c[:,j].astype(np.uint64))
    return m
masks=mk(combs); order=np.argsort(masks); sm=masks[order]
rows=[];cols=[]
for v in range(1,n2):
    if Q(v)!=1: continue
    g=W^(B(W,v)*v); m=mk(g[combs]); idx=order[np.searchsorted(sm,m)]; assert np.all(masks[idx]==m)
    rows.append(np.arange(len(combs))); cols.append(idx)
r=np.concatenate(rows); c=np.concatenate(cols)
nc,lab=connected_components(coo_matrix((np.ones(len(r)),(r,c)),shape=(len(combs),)*2),directed=True,connection='weak')
reps=[np.where(lab==i)[0][0] for i in range(nc)]
print(f'q=5: {len(combs)} sets {{0}}uT, {nc} O(Q)-orbits',flush=True)
signs=np.array(list(itertools.product((1,-1),repeat=q)),float).T  # q x 32
def best(S):
    X=np.stack([B(W,s) for s in S],1); keys=X@(1<<np.arange(q)); im,inv=np.unique(keys,return_inverse=True)
    a=np.bincount(inv,weights=(beta==1)); b=np.bincount(inv,weights=(beta==-1))
    Chi=(1-2*((im[:,None]>>np.arange(q))&1)).astype(float); m=len(im)
    idx=np.array(list(itertools.combinations(range(m),q))); bestag=-1
    for st in range(0,len(idx),20000):
        A=Chi[idx[st:st+20000]]; d=np.linalg.det(A); ok=np.abs(d)>1e-9; A=A[ok]
        if len(A)==0: continue
        K=np.linalg.inv(A)@signs; Gv=np.einsum('mq,Bqs->Bms',Chi,K)
        ag=((np.abs(Gv-1)<1e-9)*a[None,:,None]+(np.abs(Gv+1)<1e-9)*b[None,:,None]).sum(1)
        bestag=max(bestag,int(ag.max()))
    return n2-bestag
t0=time.time(); res=[]
for ri in reps:
    S=[0]+[int(x) for x in combs[ri]]; res.append((best(S),S))
res.sort(); print(f'q=5: p_min={res[0][0]} support {res[0][1]}; p-values {sorted(set(x[0] for x in res))}; {time.time()-t0:.0f}s',flush=True)
