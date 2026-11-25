# Union of the degree<=d Majorana supports of r frames (frame 1 = Jordan-Wigner, others = images under random O(Q)
# elements = products of random orthogonal transvections).  Coefficients by reweighted LAD, EXACT re-solve & count.
import numpy as np, sys, itertools
from bell import *
k=int(sys.argv[1]); d=int(sys.argv[2]); r=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
rng=np.random.default_rng(seed); n=4**k; b=beta(k); msk=(1<<k)-1
def Qs(v): return bin((v&msk)&(v>>k)).count('1')&1
def Bs(u,v): return (bin((u&msk)&(v>>k)).count('1')+bin((v&msk)&(u>>k)).count('1'))&1
A=[]
for j in range(k):
    for P in 'XY':
        a=1<<j; bb=((1<<j)-1) | ((1<<j) if P=='Y' else 0); A.append(a|(bb<<k))
t=0
for v in A: t^=v
A.append(t)
nons=[v for v in range(1,n) if Qs(v)==1]
def randg(x, steps=60, rs=rng):
    vs=[int(rs.choice(nons)) for _ in range(steps)]
    out=[]
    for y in x:
        for v in vs:
            if Bs(y,v): y^=v
        out.append(y)
    return out
def Jm(v): return ((v&msk)<<k) | (v>>k)
def supp(frame):
    S=set([0])
    for rr in range(1,d+1):
        for T in itertools.combinations(range(2*k+1),rr):
            v=0
            for i in T: v^=Jm(frame[i])
            S.add(v)
    return S
frames=[A]+[randg(A) for _ in range(r-1)]
for f in frames: assert all(Bs(u,v)==1 for u,v in itertools.combinations(f,2))
S1=sorted(supp(A)); M1=charmat(k,S1); e1,x1=rw_lad(M1,b,iters=3); p1,_=exact_errors(M1,b,x1)
S=set(S1)
for f in frames[1:]: S|=supp(f)
S=S1+sorted(S-set(S1)); M=charmat(k,S)
x=np.concatenate([x1,np.zeros(len(S)-len(S1))]); best=(p1,len(S1))
for tau in [1e-3,1e-2,1e-1,0.3]:
    for it in range(4):
        res=np.abs(M@x-b); w=1.0/(res+tau); x=lad(M,b,w)
        p,nz=exact_errors(M,b,x)
        if p<best[0]: best=(p,nz)
print('  single frame p',p1,' best union',best)
p,nz=best
print(f'k={k} d={d} frames={r} seed={seed}: |S|={len(S)} -> exact p={p}, nnz={nz}',flush=True)
