# Theorem 3 at N=3 over GF(4) (char 2, outside the proof): min rank(F - G(x)G').
import numpy as np, itertools
# GF(4) = {0,1,w,w^2} encoded 0,1,2,3 ; add = xor
MUL=np.zeros((4,4),dtype=np.int64); lg={1:0,2:1,3:2}; ex=[1,2,3]
for a in range(1,4):
    for b in range(1,4): MUL[a,b]=ex[(lg[a]+lg[b])%3]
INV=np.array([0,1,3,2])
def batch_rank4(A):
    A=A.copy(); B,n,m=A.shape; r=np.zeros(B,dtype=np.int64); ar=np.arange(B); rows=np.arange(n)
    for c in range(m):
        mask=(A[:,:,c]!=0)&(rows[None,:]>=r[:,None]); has=mask.any(1)
        if not has.any(): continue
        b=ar[has]; piv=mask[has].argmax(1); rb=r[has]
        tmp=A[b,piv].copy(); A[b,piv]=A[b,rb]; A[b,rb]=tmp
        prow=MUL[INV[A[b,rb,c]][:,None],A[b,rb]]; A[b,rb]=prow
        f=A[b,:,c].copy(); f[np.arange(len(b)),rb]=0
        A[b]=A[b]^MUL[f[:,:,None],prow[:,None,:]]
        r[has]+=1
    return r
N=3
F=np.zeros((9,9),dtype=np.int64)
for x in range(3):
    for y in range(3): F[x*3+y,y*3+x]=1
allGp=np.array(list(itertools.product(range(4),repeat=9)),dtype=np.int64).reshape(-1,3,3)
best=99
for rk in range(4):
    G=np.zeros((3,3),dtype=np.int64); G[range(rk),range(rk)]=1
    mn=99
    for s in range(0,len(allGp),32768):
        Gp=allGp[s:s+32768]
        K=MUL[G[None,:,None,:,None],Gp[:,None,:,None,:]].reshape(-1,9,9)   # [x,y,u,v] = G[x,u]Gp[y,v]
        mn=min(mn,batch_rank4(F[None]^K).min())
    print('GF(4) N=3 rankG',rk,'min rank',mn,flush=True); best=min(best,mn)
print('GF(4) N=3 q=1 min p =',best)
