import json, numpy as np
B='/home/user/pnp-recherche/code/'
OPT=[('gen19/g19b/exact_4_5_3.json',5),('gen19/g19b/exact_5_7_4.json',7),('gen19/g19b/exact_5_6_6.json',6),
     ('gen19/g19b/exact_6_7_7.json',7),('gen24/g24g/exact_6_12_3.json',12)]
def load(f):
    d=json.load(open(B+f))
    a,b=('G','Gp') if 'G' in d else ('E','C')
    A=[np.array([[int(x) for x in r] for r in m]) for m in d[a]]
    Bm=[np.array([[int(x) for x in r] for r in m]) for m in d[b]]
    return int(d['N']),A,Bm
def swap(N):
    S=np.zeros((N*N,N*N))
    for x in range(N):
        for y in range(N): S[x*N+y,y*N+x]=1   # row (x,y), col (u,v)=(y,x): [x=v][y=u]
    return S
def Kmat(A,Bm): return sum(np.kron(a,b) for a,b in zip(A,Bm)).astype(float)
def rand_mono(N,rng,dens=1.0,signed=True,weights=None):
    P=np.zeros((N,N)); perm=rng.permutation(N)
    for i in range(N):
        if rng.random()<dens:
            P[i,perm[i]]=(rng.choice([-1,1]) if signed else 1)*(1 if weights is None else weights(rng))
    return P
