# check (q,p)=(6,19)@8: S={0} u 5-vector Majorana frame of a 4-dim nondegenerate subspace (qubits 1,2), kappa=(1/2)(-1)^{Q(s)}
import numpy as np, itertools, flint
k=3; N=8; n2=64; msk=7
POP=np.array([bin(i).count('1') for i in range(256)])
Q=lambda x: POP[(x&msk)&(x>>k)]&1
B=lambda x,v: (POP[(x&msk)&(v>>k)]+POP[(x>>k)&(v&msk)])&1
W=np.arange(n2); beta=1-2*Q(W)
sub=[w for w in range(n2) if (w&4)==0 and (w&32)==0]   # qubits 0,1 only
best=None
for fr in itertools.combinations([w for w in sub if w],4):
    if all(B(a,b)==1 for a,b in itertools.combinations(fr,2)):
        t=fr[0]^fr[1]^fr[2]^fr[3]; S=[0,*fr,t]
        phi=sum(0.5*(1-2*Q(s))*(1-2*B(W,s)) for s in S)
        p=int(np.sum(np.abs(phi-beta)>1e-9)); best=(p,S) if best is None or p<best[0] else best
print('Bell count p =',best)
p,S=best
def pauli(s):
    c,d=s&msk,s>>k; M=np.zeros((N,N),dtype=np.int64)
    for x in range(N): M[x^c,x]=-1 if bin(d&x).count('1')%2 else 1
    return M
SW=np.zeros((n2,n2),dtype=np.int64)
for i in range(N):
    for j in range(N): SW[i*N+j,j*N+i]=1
for pr in (1000003,998244353):
    inv2=pow(2,-1,pr); K=np.zeros((n2,n2),dtype=np.int64)
    for s in S: K=(K+(int(1-2*Q(s))*inv2%pr)*np.kron(pauli(s),pauli(s)))%pr
    rk=lambda M: flint.nmod_mat(M.shape[0],M.shape[1],[int(x) for x in (M%pr).ravel()],pr).rank()
    print(pr,'osr',rk(K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(n2,n2)),'rank(SWAP-K)',rk(SW-K))
