# verify Bell-model formula: for qubit N=2^k, K = sum_S kappa P(x)P, rank(SWAP-K) == #{phi != beta}
import numpy as np, itertools, sys
from fractions import Fraction
X=np.array([[0,1],[1,0]]); Z=np.array([[1,0],[0,-1]]); I2=np.eye(2,dtype=int)
def pauli(a,b,k):
    M=np.array([[1]])
    for i in range(k):
        P=I2
        if a[i]: P=X
        if b[i]: P=P@Z if a[i] else Z
        M=np.kron(M,P)
    return M
k=int(sys.argv[1]) if len(sys.argv)>1 else 2
N=2**k
vecs=list(itertools.product([0,1],repeat=k))
SW=np.zeros((N*N,N*N),dtype=int)
for x in range(N):
    for y in range(N): SW[x*N+y,y*N+x]=1
# check SWAP = (1/N) sum eps P(x)P
tot=np.zeros((N*N,N*N))
paulis={}
for a in vecs:
    for b in vecs:
        P=pauli(a,b,k); paulis[(a,b)]=P
        tot+= (-1)**(np.dot(a,b)%2)*np.kron(P,P)
print("SWAP identity ok:", np.allclose(tot/N,SW))
rng=np.random.default_rng(1)
keys=list(paulis)
def chi(ab,cd): a,b=ab; c,d=cd; return (-1)**((np.dot(a,d)+np.dot(b,c))%2)
pts=[(c,d) for c in vecs for d in vecs]
for trial in range(20):
    q=rng.integers(1,len(keys)); S=[keys[i] for i in rng.choice(len(keys),q,replace=False)]
    kap={s:int(rng.integers(-2,3)) for s in S}
    K=sum(kap[s]*np.kron(paulis[s],paulis[s]) for s in S)
    r=np.linalg.matrix_rank(SW-K)
    cnt=sum(1 for v in pts if (-1)**(np.dot(v[0],v[1])%2)!=sum(kap[s]*chi(s,v) for s in S))
    print(q,r,cnt, r==cnt)
