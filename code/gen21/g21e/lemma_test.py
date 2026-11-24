# exact test of key lemma: max_{G,G'} rank A(G(x)G')B >= rank A rank B / N^2, over F_p
import numpy as np, sys, itertools
from ff import rank_mod
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 1)
p = int(sys.argv[1])
def rand_rank(m, n, r):
    return (rng.integers(0,p,(m,r)) @ rng.integers(0,p,(r,n))) % p
def gens(N):
    D=N*N; I=np.eye(N,dtype=np.int64)
    # random prescribed ranks
    for _ in range(150):
        a=rng.integers(1,D+1); b=rng.integers(1,D+1)
        yield 'rand', rand_rank(a,D,a), rand_rank(D,b,b)
    # Kronecker coordinate subspaces
    for k1 in range(1,N+1):
     for k2 in range(1,N+1):
      for j1 in range(1,N+1):
       for j2 in range(1,N+1):
        A=np.kron(I[:k1],I[:k2]); B=np.kron(I[:,:j1],I[:,:j2])
        yield 'kron',A,B
    # symmetric / antisymmetric (projection-free: spanning vectors)
    sym=[];alt=[]
    for i in range(N):
        for j in range(i,N):
            v=np.zeros(D,dtype=np.int64); v[i*N+j]+=1; v[j*N+i]+=1; sym.append(v%p)
            if i<j:
                w=np.zeros(D,dtype=np.int64); w[i*N+j]=1; w[j*N+i]=p-1; alt.append(w)
    S=np.array(sym); Al=np.array(alt)
    for X in [S,Al]:
        for Y in [S,Al]:
            yield 'symalt', X, Y.T
    # graph of a matrix: subspace {(x, Mx)} type / diagonal span
    diag=np.array([np.eye(D,dtype=np.int64)[i*N+i] for i in range(N)])
    yield 'diag', diag, diag.T
    yield 'diag-full', diag, np.eye(D,dtype=np.int64)
    # random coordinate subsets
    for _ in range(100):
        R=rng.permutation(D)[:rng.integers(1,D+1)]; C=rng.permutation(D)[:rng.integers(1,D+1)]
        yield 'coord', np.eye(D,dtype=np.int64)[R], np.eye(D,dtype=np.int64)[:,C]
worst=(9,None)
for N in [2,3,4]:
    wN=9
    for name,A,B in gens(N):
        a=rank_mod(A,p); b=rank_mod(B,p)
        if a*b==0: continue
        trials = 3 if p>1000 else 30
        s=0
        for _ in range(trials):
            G=rng.integers(0,p,(N,N)); H=rng.integers(0,p,(N,N))
            s=max(s, rank_mod(A@np.kron(G,H)%p@B, p))
            if s*N*N>=a*b and s==min(a,b): break
        ratio=s*N*N/(a*b)
        if ratio<wN: wN=ratio; wc=(name,a,b,s)
    print('p',p,'N',N,'min s N^2/(ab) =',round(wN,3),'at',wc, flush=True)
