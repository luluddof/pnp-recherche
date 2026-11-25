# Does a non-symmetric support beat the Majorana family?  S = Majorana support of degree <= d (all products of <= d of the
# 2k+1 Majorana characters) [+ extra characters chosen greedily by residual correlation]; coefficients by reweighted LAD;
# EXACT rational re-solve and count (rigorous points).
import numpy as np, sys, itertools
from bell import *
k=int(sys.argv[1]); d=int(sys.argv[2]); extra=int(sys.argv[3]); seed=int(sys.argv[4]) if len(sys.argv)>4 else 0
rng=np.random.default_rng(seed); n=4**k; b=beta(k)
A=[]
for j in range(k):
    for P in 'XY':
        a=1<<j; bb=((1<<j)-1) | ((1<<j) if P=='Y' else 0); A.append(a|(bb<<k))
t=0
for v in A: t^=v
A.append(t)
# characters chi_s with s = J(sum of subset)?  chi_{s}(w) = (-1)^{s.w};  Majorana bit y_a(w) = Q(a)+B(a,w) = Q(a) + (Ja).w
def Jm(v): return ((v&((1<<k)-1))<<k) | (v>>k)
S=set([0])
for r in range(1,d+1):
    for T in itertools.combinations(range(2*k+1),r):
        v=0
        for i in T: v^=Jm(A[i])
        S.add(v)
S=sorted(S); q0=len(S)
M=charmat(k,S); e,x=rw_lad(M,b,iters=4); p,nz=exact_errors(M,b,x)
print(f'k={k} d={d}: Majorana-support |S|={q0}: rwLAD p={p} nnz={nz}',flush=True)
for it in range(extra):
    r=M@x-b; E=np.abs(r)>1e-7
    sc=np.abs(fwht(np.where(E,np.sign(r),0.0))); sc[S]=-1; sc+=1e-9*rng.random(n)
    S=S+[int(np.argmax(sc))]
    M=charmat(k,S); e,x=rw_lad(M,b,iters=4); p,nz=exact_errors(M,b,x)
    print(f'   +{it+1} extra: |S|={len(S)} p={p} nnz={nz}',flush=True)
