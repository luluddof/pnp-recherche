from load import *
rng=np.random.default_rng(1)
def bounds(N,A,Bm):
    S,K=mats(N,A,Bm); Th=(K@S).astype(float); n=N*N; X=np.eye(n)-Th
    s=np.linalg.svd(X,compute_uv=False); p=int(np.sum(s>1e-9)); trX=np.trace(X)
    S4=(s**4).sum(); KtK=(K.T@K).astype(float)
    energy=np.sum(KtK**2)/n
    return p, trX**(4/3)/S4**(1/3), energy
print('exact optima: p, Thm C bound (trX)^{4/3}/||X||_S4^{4/3}, Energy=||K^T K||_F^2/n')
for f in FILES:
    N,A,Bm=load(f); p,b,e=bounds(N,A,Bm); print(' ',f.split('/')[-1],p,round(b,2),'energy/q^2=%.2f'%(e/len(A)**2))
def randmono(N,dens=1.0):
    P=np.zeros((N,N),dtype=int); perm=rng.permutation(N)
    for i in range(N):
        if rng.random()<dens: P[perm[i],i]=rng.choice([1,-1])
    return P
print('random full signed perms: N, q, p, bound, energy/q^2')
for N in (8,12):
    for q in (2,4,8):
        A=[randmono(N) for _ in range(q)];Bm=[randmono(N) for _ in range(q)]
        p,b,e=bounds(N,A,Bm); print(' ',N,q,p,round(b,1),round(e/q**2,2))
