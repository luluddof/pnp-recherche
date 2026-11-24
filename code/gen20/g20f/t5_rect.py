# Theorem 5 rectangular step: S: K^k (x) K^m -> K^m (x) K^k swap, g~ : K^k->K^m (m x k), gam: K^m->K^k (k x m).
# claim: rank(S - g~ (x) gam) >= k m - k'(k'+1)/2, k'=min(k,m), char != 2; charpoly of Theta = S^{-1}(g~(x)gam).
import numpy as np, itertools, sympy as sp, random
from brank import batch_rank
def Smat(k,m):
    S=np.zeros((m*k,k*m),dtype=np.int64)
    for i in range(k):
        for j in range(m): S[j*k+i, i*m+j]=1   # e_i(x)f_j -> f_j(x)e_i
    return S
# 1) charpoly identity, exact over Q with random integers
t=sp.symbols('t')
for (k,m) in [(1,3),(2,3),(3,2),(2,4),(3,5),(4,2),(3,3)]:
    for trial in range(3):
        gt=sp.Matrix(m,k,lambda a,b:random.randint(-3,3)); ga=sp.Matrix(k,m,lambda a,b:random.randint(-3,3))
        S=sp.Matrix(Smat(k,m)); Th=S.inv()*sp.kronecker_product(gt,ga)  # g~(x)gam maps K^k(x)K^m -> K^m(x)K^k
        cp=sp.expand(Th.charpoly(t).as_expr())
        kp=min(k,m); small=(ga*gt) if k<=m else (gt*ga)
        cs=small.charpoly(t).as_expr()     # roots mu_i, degree kp
        # build prod (t-mu_i) prod_{i<j}(t^2-mu_i mu_j) symbolically via resultant-free route: use power sums -> compare numerically
        mus=[complex(r) for r in sp.Poly(cs,t).nroots(n=30)]
        val=lambda z: complex(cp.subs(t,z))
        pred=lambda z: z**(k*m-kp*kp)*np.prod([z-a for a in mus])*np.prod([z*z-mus[i]*mus[j] for i in range(kp) for j in range(i+1,kp)])
        zs=[0.3+0.7j,1.1-0.2j,-0.9+0.4j,2.0]
        err=max(abs(val(z)-pred(z))/(1+abs(pred(z))) for z in zs)
        assert sp.degree(cp,t)==k*m
        print(f'k={k} m={m}: charpoly relative err {err:.1e}')
# 2) exhaustive rank bound over GF(p)
def exh(k,m,p):
    S=Smat(k,m); best=10**9
    allg=np.array(list(itertools.product(range(p),repeat=m*k)),dtype=np.int64).reshape(-1,m,k)
    allga=np.array(list(itertools.product(range(p),repeat=k*m)),dtype=np.int64).reshape(-1,k,m)
    for ga in allga:
        K=np.einsum('bab,cd->bacbd'.replace('bab','zab').replace('bacbd','zacbd'),allg,ga).reshape(-1,m*k,k*m)
        r=batch_rank(S[None]-K,p); best=min(best,r.min())
    kp=min(k,m); print(f'GF({p}) k={k} m={m}: min rank={best}, bound km-k\'(k\'+1)/2={k*m-kp*(kp+1)//2}',flush=True)
for (k,m,p) in [(1,2,2),(1,2,3),(2,2,2),(2,2,3),(2,3,2),(2,3,3),(3,2,3),(3,3,2)]:
    exh(k,m,p)
