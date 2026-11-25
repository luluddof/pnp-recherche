# Explicit operator check: K = sum_s kappa_s P_s (x) P_s (P_s = X^a Z^b, real), rank(SWAP - K) and osr(K) mod 2 primes.
import numpy as np, sys
from bell import *
def pauli(a,b,k):
    X=np.array([[0,1],[1,0]],dtype=np.int64); Z=np.array([[1,0],[0,-1]],dtype=np.int64); I=np.eye(2,dtype=np.int64)
    M=np.array([[1]],dtype=np.int64)
    for i in range(k):
        P=I
        if (a>>i)&1: P=P@X
        if (b>>i)&1: P=P@Z
        M=np.kron(M,P)
    return M
def rank_mod(A,p):
    A=(A%p).astype(np.int64).copy(); r=0; R,C=A.shape
    for c in range(C):
        piv=np.where(A[r:,c]!=0)[0]
        if len(piv)==0: continue
        i=r+piv[0]; A[[r,i]]=A[[i,r]]; inv=pow(int(A[r,c]),p-2,p); A[r]=(A[r]*inv)%p
        nz=np.where(A[:,c]!=0)[0]
        for j in nz:
            if j!=r: A[j]=(A[j]-A[j,c]*A[r])%p
        r+=1
        if r==R: break
    return r
def check(k,S,coef2):  # coef2 = 2*kappa (integers)
    N=2**k; D=N*N
    SW=np.zeros((D,D),dtype=np.int64)
    for i in range(N):
        for j in range(N): SW[i*N+j, j*N+i]=1
    res=[]
    for conv in (0,1):
        K2=np.zeros((D,D),dtype=np.int64); terms=[]
        for s,c in zip(S,coef2):
            a=s&(N-1); b=s>>k
            P=pauli(a,b,k) if conv==0 else pauli(b,a,k)
            K2+=c*np.kron(P,P); terms.append(P)
        A=2*SW-K2
        rk=[rank_mod(A,p) for p in (1000003,998244353)]
        # osr: realignment rank of K2
        R=K2.reshape(N,N,N,N).transpose(0,2,1,3).reshape(N*N,N*N)
        osr=[rank_mod(R,p) for p in (1000003,998244353)]
        res.append((rk,osr))
    return res
if __name__=='__main__':
    # (6,19)@8 and Majorana (10,46)@16 (greedy solution earlier), coefficients +-1/2 -> coef2 = +-1
    print('(q,p)=(6,19)@8:', check(3,[0,26,13,52,48,19],[1,-1,-1,-1,1,-1]))
    P=['YIYX','IYZY','ZIII','IZYX','YYZI','XZXY','YYXX','XZYI','IXXZ','XIZY']; sg=[1,1,1,-1,1,-1,1,-1,1,-1]
    def vec(p):
        a=b=0
        for i,ch in enumerate(p):
            if ch in 'XY': a|=1<<i
            if ch in 'ZY': b|=1<<i
        return a|(b<<4)
    print('(q,p)=(10,46)@16:', check(4,[vec(p) for p in P],sg))
