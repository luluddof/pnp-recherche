# plain Pauli truncation K = (1/N) sum_{P in S} eps(P) P(x)P (qubit, Bell-diagonal): p(S) = #{v: sum_{P notin S} eps chi_P(v) != 0}
import numpy as np, itertools, sys
k=int(sys.argv[1]); N=2**k; n=N*N
pts=np.array(list(itertools.product([0,1],repeat=2*k)))
a=pts[:,:k]; b=pts[:,k:]
eps=(-1)**((a*b).sum(1)%2)
# chi_P(v) = (-1)^{a.d + b.c}, v=(c,d)
c=pts[:,:k]; d=pts[:,k:]
CH=(-1)**(((a@d.T)+(b@c.T))%2)   # [P,v]
best={}
if k==2:
    for mask in range(1<<n):
        inS=np.array([(mask>>i)&1 for i in range(n)],bool)
        lam=(eps*(~inS))@CH
        p=int((lam!=0).sum()); q=int(inS.sum())
        if q not in best or p<best[q]: best[q]=p
else:
    rng=np.random.default_rng(0)
    # greedy: grow complement T (Paulis kept out of K) is S^c ; start S=all, remove Paulis greedily to keep p minimal
    for trial in range(20):
        inS=np.ones(n,bool)
        order=rng.permutation(n)
        while inS.sum()>0:
            cand=np.nonzero(inS)[0]; rng.shuffle(cand); cand=cand[:64]
            bestp=None
            for j in cand:
                inS[j]=False; lam=(eps*(~inS))@CH; p=int((lam!=0).sum()); inS[j]=True
                if bestp is None or p<bestp[0]: bestp=(p,j)
            inS[bestp[1]]=False; q=int(inS.sum())
            if q not in best or bestp[0]<best[q]: best[q]=bestp[0]
for q in sorted(best): 
    if q%max(1,n//32)==0 or q<N+2 or q>n-N-2: print("N=%d q=%d min p=%d"%(N,q,best[q]))
