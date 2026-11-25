# Theorem G, non-regular action: Gamma = Z_k acting freely on [N] by a fixed-point-free permutation iota of order k
# (orbits of size k); irreps d=1, multiplicities m = N/k each (k of them): bound (N^2 - k (N/k)^2)/2 = (N^2 - N^2/k)/2.
import numpy as np, random
random.seed(3)
def commuting_signed_perm(N,k):
    # [N] = Z_k x [N/k], iota(c,i) = (c+1,i); commuting perms: (c,i) -> (c + f(i), pi(i)), sign s(i)
    M=N//k; pi=list(range(M)); random.shuffle(pi); f=[random.randrange(k) for _ in range(M)]; s=[random.choice([1,-1]) for _ in range(M)]
    A=np.zeros((N,N),dtype=int)
    for c in range(k):
        for i in range(M):
            if random.random()<0.85: A[((c+f[i])%k)*M+pi[i], c*M+i]=s[i]
    return A
def partial_fix(A,k,N):
    # make partial deletions orbit-consistent (delete whole orbits) to keep commuting
    M=N//k; B=A.copy()
    for i in range(M):
        cols=[c*M+i for c in range(k)]
        if any(np.all(B[:,col]==0) for col in cols): B[:,cols]=0
    return B
for N,k in [(8,2),(12,3),(12,2),(16,4)]:
    S=np.zeros((N*N,N*N))
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1
    bound=(N*N-N*N/k)/2; best=10**9
    for t in range(150):
        q=random.choice([1,2,3,4,6])
        K=sum(np.kron(partial_fix(commuting_signed_perm(N,k),k,N),partial_fix(commuting_signed_perm(N,k),k,N)) for _ in range(q))
        p=np.linalg.matrix_rank(S-K); best=min(best,p); assert p>=bound-1e-9,(N,k,p,bound)
    print('N',N,'Z_%d free action'%k,'bound',bound,'min p over 150 random commuting instances',best,' N(N-1)/2=',N*(N-1)//2)
