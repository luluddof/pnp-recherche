# S_k-invariant K = sum_r c_r (X_r (x) Y_r)^{(x)k} (pair-site product), osr <= R. For R <= 2 at k=4, u(q)=0 => p >= 73.
import numpy as np
rng=np.random.default_rng(5); k=4; N=16
def swap(N):
    S=np.zeros((N*N,N*N))
    for x in range(N):
        for y in range(N): S[y*N+x,x*N+y]=1
    return S
S=swap(N)
def kp(M,k):
    o=np.array([[1.0+0j]])
    for _ in range(k): o=np.kron(o,M)
    return o
def rk(M,tol=1e-8):
    s=np.linalg.svd(M,compute_uv=False); return int((s>tol*max(1,s[0])).sum())
P=[np.eye(2),np.array([[0,1],[1,0]]),np.array([[0,-1j],[1j,0]]),np.diag([1.,-1])]
mins=99999
for t in range(400):
    R=rng.integers(1,3); K=0
    for r in range(R):
        if rng.random()<0.5:
            X=P[rng.integers(4)]; Y=X.conj() if rng.random()<.5 else P[rng.integers(4)]
        else:
            X=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2)); Y=rng.normal(size=(2,2))
        c=rng.normal()+1j*rng.normal() if rng.random()<.5 else rng.choice([1,-1,0.5])
        K=K+c*np.kron(kp(X,k),kp(Y,k))
    mins=min(mins,rk(S-K))
print("min p over 400 random S_4-invariant K with osr<=2:",mins," (S' bound 73, Theorem-G-type/q=1 optimum 120)")
