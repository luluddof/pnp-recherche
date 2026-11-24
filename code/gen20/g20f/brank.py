import numpy as np
def batch_rank(A,p):
    A=np.array(A,dtype=np.int64)%p; B,n,m=A.shape; r=np.zeros(B,dtype=np.int64); ar=np.arange(B); rows=np.arange(n)
    inv=np.zeros(p,dtype=np.int64)
    for a in range(1,p): inv[a]=pow(a,p-2,p)
    for c in range(m):
        mask=(A[:,:,c]!=0)&(rows[None,:]>=r[:,None])
        has=mask.any(1)
        if not has.any(): continue
        b=ar[has]; piv=mask[has].argmax(1); rb=r[has]
        tmp=A[b,piv].copy(); A[b,piv]=A[b,rb]; A[b,rb]=tmp
        prow=(A[b,rb]*inv[A[b,rb,c]][:,None])%p; A[b,rb]=prow
        f=A[b,:,c].copy(); f[np.arange(len(b)),rb]=0
        A[b]=(A[b]-f[:,:,None]*prow[:,None,:])%p
        r[has]+=1
    return r
