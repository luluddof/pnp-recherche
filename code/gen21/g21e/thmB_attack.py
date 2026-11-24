# Attack Theorem B: ratio R^2/(N^2 c1 c2) for (i) single class-B blocks, all 2^4 bipartitions, coordinate/random A,B, X = F and random X;
# (ii) maps OUTSIDE the stated class: Phi(X) = P Fl_S(C X D) Q with C, D general (non-local) -- notes claim compositions are in B.
import numpy as np, itertools, sys
from ff import rank_mod, flat
p=10007; rng=np.random.default_rng(3)
def swap(N):
    F=np.zeros((N*N,N*N),dtype=np.int64)
    for x in range(N):
        for y in range(N): F[x*N+y,y*N+x]=1
    return F
def rnd(*s): return rng.integers(0,p,s)
def prof(phi,N,X,t=2):
    c1=max(rank_mod(phi(np.outer(rnd(N*N),rnd(N*N))),p) for _ in range(t))
    c2=max(rank_mod(phi(np.kron(rnd(N,N),rnd(N,N))),p) for _ in range(t))
    R=rank_mod(phi(X),p); return c1,c2,R
splits=[S for r in range(0,5) for S in itertools.combinations(range(4),r)]
for N in [2,3]:
    D=N*N; F=swap(N); best=(0,None); best2=(0,None)
    for it in range(400 if N==2 else 150):
        S=splits[rng.integers(len(splits))]
        m,n=N**len(S),N**(4-len(S))
        mode=rng.integers(3)
        if mode==0: A=np.eye(m,dtype=np.int64)[rng.permutation(m)[:rng.integers(1,m+1)]]; B=np.eye(n,dtype=np.int64)[:,rng.permutation(n)[:rng.integers(1,n+1)]]
        elif mode==1: a=rng.integers(1,m+1); b=rng.integers(1,n+1); A=rnd(a,a)@rnd(a,m)%p; B=rnd(n,b)@rnd(b,b)%p
        else:
            ka=[rng.integers(1,N+1) for _ in S]; kb=[rng.integers(1,N+1) for _ in range(4-len(S))]
            A=np.ones((1,1),dtype=np.int64)
            for k in ka: A=np.kron(A,rnd(k,N))
            B=np.ones((1,1),dtype=np.int64)
            for k in kb: B=np.kron(B,rnd(N,k))
        for X,name in [(F,'F'),(rnd(D,D),'rand')]:
            phi=lambda Y: (A@flat(Y,N,S)%p)@B%p
            c1,c2,R=prof(phi,N,X)
            if c1*c2==0: continue
            r=R*R/(N*N*c1*c2)
            if r>best[0]: best=(r,(S,mode,name,c1,c2,R))
        # outside-class composition
        C=rnd(D,D) if rng.random()<.5 else np.kron(rnd(N,N),rnd(N,N)); Dm=rnd(D,D) if rng.random()<.5 else np.eye(D,dtype=np.int64)
        S2=splits[rng.integers(len(splits))]; m,n=N**len(S2),N**(4-len(S2))
        P=np.eye(m,dtype=np.int64)[rng.permutation(m)[:rng.integers(1,m+1)]]; Q=np.eye(n,dtype=np.int64)[:,rng.permutation(n)[:rng.integers(1,n+1)]]
        phi2=lambda Y: P@flat(C@Y%p@Dm%p,N,S2)%p@Q%p
        c1,c2,R=prof(phi2,N,F)
        if c1*c2:
            r=R*R/(N*N*c1*c2)
            if r>best2[0]: best2=(r,(S2,c1,c2,R))
    print('N',N,'class-B single block max R^2/(N^2c1c2)=%.3f'%best[0],best[1])
    print('N',N,'composition P Fl_S(CXD) Q, X=F: max ratio=%.3f'%best2[0],best2[1],flush=True)
