# SA over signed partial-perm MC instances (q fixed) minimizing rank(Swap-K); compare with Theorem-A bounds (own certs).
import numpy as np, sys, time
from common import swap
N=int(sys.argv[1]); q=int(sys.argv[2]); T=float(sys.argv[3]); seed=int(sys.argv[4])
rng=np.random.default_rng(seed); S=swap(N)
def rand(): 
    M=np.zeros((N,N)); pr=rng.permutation(N)
    for i in range(N):
        if rng.random()<0.8: M[i,pr[i]]=rng.choice([-1,1])
    return M
def cost(A,B):
    s=np.linalg.svd(S-sum(np.kron(a,b) for a,b in zip(A,B)),compute_uv=False)
    r=int(np.sum(s>1e-8)); return r+ 0.5*np.sum(s[s>1e-8][-max(1,r//4):])/ (1+np.sum(s)) ,r
def mutate(M):
    M=M.copy(); i,j=rng.integers(N,size=2); t=rng.random()
    if t<0.4: M[[i,j]]=M[[j,i]]
    elif t<0.7: M[i]*=-1
    else:
        if M[i].any(): M[i]=0
        else:
            free=[c for c in range(N) if not M[:,c].any()]
            if free: M[i,rng.choice(free)]=rng.choice([-1,1])
    return M
A=[rand() for _ in range(q)]; B=[rand() for _ in range(q)]; c,r=cost(A,B); best=r; t0=time.time(); temp=1.0
while time.time()-t0<T:
    A2=[a for a in A]; B2=[b for b in B]; l=rng.integers(q)
    if rng.random()<0.5: A2[l]=mutate(A2[l])
    else: B2[l]=mutate(B2[l])
    c2,r2=cost(A2,B2)
    if c2<=c or rng.random()<np.exp((c-c2)/temp): A,B,c,r=A2,B2,c2,r2
    best=min(best,r); temp=max(0.05,temp*0.9995)
print('N',N,'q',q,'best p',best,'seed',seed)
