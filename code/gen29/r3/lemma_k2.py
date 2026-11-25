# exhaustive N=4 (k=2) Bell optimum for q=1..5 (all supports, vertex enumeration) -> test "q<=2k => p>=N^2/4"
import itertools, numpy as np
k=2; n2=16; msk=3
POP=np.array([bin(i).count('1') for i in range(64)])
Q=lambda x: POP[(x&msk)&(x>>k)]&1
B=lambda x,v: (POP[(x&msk)&(v>>k)]+POP[(x>>k)&(v&msk)])&1
W=np.arange(n2); beta=1-2*Q(W)
for q in range(1,6):
    signs=np.array(list(itertools.product((1,-1),repeat=q)),float).T; pm=99
    for S in itertools.combinations(range(n2),q):
        X=np.stack([B(W,s) for s in S],1); keys=X@(1<<np.arange(q)); im,inv=np.unique(keys,return_inverse=True)
        a=np.bincount(inv,weights=(beta==1)); b=np.bincount(inv,weights=(beta==-1))
        Chi=(1-2*((im[:,None]>>np.arange(q))&1)).astype(float)
        idx=np.array(list(itertools.combinations(range(len(im)),q)))
        A=Chi[idx]; ok=np.abs(np.linalg.det(A))>1e-9; A=A[ok]
        if not len(A): continue
        Gv=np.einsum('mq,Bqs->Bms',Chi,np.linalg.inv(A)@signs)
        ag=((np.abs(Gv-1)<1e-9)*a[None,:,None]+(np.abs(Gv+1)<1e-9)*b[None,:,None]).sum(1).max()
        pm=min(pm,n2-int(ag))
    print(f'N=4 q={q}: p_min={pm}  (N^2/4={n2//4})',flush=True)
