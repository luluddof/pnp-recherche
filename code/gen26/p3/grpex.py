# exhaustive min rank over group-based MC instances (A_l=+-L(g), B_l=L(h)), small q
import numpy as np, itertools, sys
sys.argv+=[]
from grpthm import G, irrep_dims
name=sys.argv[1]; qmax=int(sys.argv[2])
tab=G[name]; ds,L=irrep_dims(tab); N=len(tab)
S=np.zeros((N*N,N*N))
for x in range(N):
    for y in range(N): S[x*N+y,y*N+x]=1
terms=[(s,a,b) for s in (1,-1) for a in range(N) for b in range(N)]
KR=[s*np.kron(L[a],L[b]) for s,a,b in terms]
bound=(N*N-sum(d**3 for d in ds))/2
for q in range(1,qmax+1):
    best=10**9;arg=None
    for comb in itertools.combinations_with_replacement(range(len(terms)),q):
        K=sum(KR[i] for i in comb); p=np.linalg.matrix_rank(S-K)
        if p<best: best=p;arg=[terms[i] for i in comb]
    print(name,'N',N,'q',q,'exhaustive min p',best,'bound',bound,'N(N-1)/2',N*(N-1)//2,'arg',arg,flush=True)
