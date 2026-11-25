# structured larger-N MC tests: A(x)I_m, A(x)P (P random signed perm), direct sums; compare p with best verified q=3 bound
import json, glob, numpy as np
from fractions import Fraction as F
from common import *
certs=[json.load(open(f)) for f in glob.glob('cert_q3_*.json')]
def best_bound(N): return max((float(F(c['a']))*N*N-float(F(c['L']))*N)/(1+float(F(c['a']))) for c in certs)
rng=np.random.default_rng(3)
for f in ['gen19/g19b/exact_4_5_3.json','gen24/g24g/exact_6_12_3.json']:
    N0,A,Bm=load(f)
    for m in (2,3,4):
        for kind in ('I','P'):
            if kind=='I': X=[np.eye(m)]*3; Y=[np.eye(m)]*3
            else: X=[rand_mono(m,rng) for _ in range(3)]; Y=[rand_mono(m,rng) for _ in range(3)]
            AA=[np.kron(a,x) for a,x in zip(A,X)]; BB=[np.kron(b,y) for b,y in zip(Bm,Y)]; N=N0*m
            p=np.linalg.matrix_rank(swap(N)-Kmat(AA,BB),tol=1e-8)
            print(f'{f.split("/")[-1]} (x){kind}_{m}: N={N} q=3 p={p}  N^2/2={N*N/2}  best ThmA bound={best_bound(N):.1f}  ok={p>=best_bound(N)}')
