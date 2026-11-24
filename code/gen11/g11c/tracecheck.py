# Sanity check of the Trace Lemma on numerically found general (non-commutative) solutions.
# alpha_x = [t^1] J(t)[x, ~x];  for a saturated block Q(T,sigma):  sum_{x in Q} alpha_x = - sum_{k in T} n_k(Q)/tau_k.
import numpy as np, itertools
from feas import Prob
def run(l,d,tau,prof,seed):
    P=Prob(l,d,tau,prof); rng=np.random.default_rng(seed)
    for it in range(50):
        o=P.lm(rng.normal(size=P.np_),iters=600,tol=1e-28)
        if o.cost<1e-24: break
    assert o.cost<1e-24
    A=P.J_of(o.x); n=1<<l; alpha=np.array([A[1][x,(n-1)^x] for x in range(n)])
    print("(l,d)=(%d,%d) profile %s, verified ranks %s"%(l,d,prof,P.check(A)))
    for s in range(d+1,l+1):
        for T in itertools.combinations(range(l),s):
            mask=sum(1<<k for k in T)
            for sig in range(n):
                if sig & ~mask: continue
                Q=[x for x in range(n) if x&mask==sig]; nQ=[(n-1)^x for x in Q]
                nul=[]
                for k in T:
                    M=sum(A[i]*tau[k]**i for i in range(d+1))[np.ix_(Q,nQ)]
                    nul.append(len(Q)-int((np.linalg.svd(M,compute_uv=False)>1e-8).sum()))
                assert sum(nul)<=d*len(Q)
                if sum(nul)==d*len(Q):
                    lhs=alpha[Q].sum(); rhs=-sum(m/tau[k] for m,k in zip(nul,T))
                    print("  saturated block T=%s sigma=%s nullities=%s: sum alpha=%.10f  -sum n_k/tau_k=%.10f"%(T,bin(sig),nul,lhs,rhs))
                    assert abs(lhs-rhs)<1e-6
run(3,1,[1.0,1.7,2.9],[1,1,2],0)
run(4,2,[1.0,1.7,2.9,4.3],[1,1,1,1],1)
run(4,1,[1.0,1.7,2.9,4.3],[3,3,3,3],2)
print("trace lemma checks passed")
