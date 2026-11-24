import numpy as np, json, itertools
d=json.load(open('exact_6_12_3.json')); E=np.array(d['E']); C=np.array(d['C']); N=6
K=sum(np.kron(e,c) for e,c in zip(E,C))
S=np.zeros((36,36),int)
for x in range(N):
    for y in range(N): S[x*N+y,y*N+x]=1
M=S-K
print("entries of K off supp(Swap):", int(((K!=0)&(S==0)).sum()), " on supp:", int(((K!=0)&(S!=0)).sum()))
print("K values on supp(Swap):", sorted(set(K[S!=0].tolist())))
cnt=0; grp=[]
for perm in itertools.permutations(range(N)):
    for signs in itertools.product([1,-1],repeat=N):
        g=np.zeros((N,N),int)
        for i,pi in enumerate(perm): g[pi,i]=signs[i]
        G=np.kron(g,g)
        if np.array_equal(G@K, K@G): cnt+=1; grp.append((perm,signs))
print("signed-permutation symmetries g (g(x)g commuting with K):", cnt)
print(grp[:12])
Th=K@S; ev=np.linalg.eigvals(Th.astype(float)); 
from collections import Counter; print("Theta spectrum:", Counter(np.round(ev.real,6)+0))
