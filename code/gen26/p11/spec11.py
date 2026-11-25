import sys; sys.path.insert(0,'/home/user/pnp-recherche/code/gen26/p3')
from load import *
import numpy as np
from collections import Counter
np.set_printoptions(linewidth=200,precision=3,suppress=True)
def report(name,N,K,S,Qn=None):
    n=N*N; Th=K@S; p=np.linalg.matrix_rank(np.eye(n)-Th,tol=1e-8); m=n-p
    Ts=(Th+Th.T)/2; ev=np.linalg.eigvalsh(Ts)
    r=lambda x: round(float(x),3)
    c=Counter(np.round(ev,3))
    print(f"== {name} N={N} p={p} m={m} trTheta={r(np.trace(Th))} ||Ts||={r(max(abs(ev)))} ||K||={r(np.linalg.norm(K,2))} Q={Qn}")
    print("  #ev>1+e",int(np.sum(ev>1+1e-8)),"#ev=1",int(np.sum(abs(ev-1)<1e-8)),"#ev in(0,1)",int(np.sum((ev>1e-8)&(ev<1-1e-8))),"#ev=0",int(np.sum(abs(ev)<1e-8)),"#ev<0",int(np.sum(ev<-1e-8)))
    print("  Ts eigen multiset (non-1):",sorted([(k,v) for k,v in c.items() if abs(k-1)>1e-6]))
    rest=ev[abs(ev-1)>1e-8]
    print("  sum(rest)=",r(rest.sum()),"  -(n-p-#extra1) check; odd moments rest j=1,3,5:",[r((rest**j).sum()) for j in (1,3,5)],"tr Ts^j:",[r((ev**j).sum()) for j in (1,3,5)])
    eT=np.linalg.eigvals(Th); cT=Counter(np.round(eT,3))
    print("  Theta eigen multiset:",sorted(cT.items(),key=lambda t:(-t[1],t[0].real)))
    sv=np.linalg.svd(np.eye(n)-Th,compute_uv=False); print("  sv(1-Theta) nonzero:",sv[:p+1])
    svK=np.linalg.svd(K,compute_uv=False); print("  sv(K) multiset:",sorted(Counter(np.round(svK,3)).items(),reverse=True))
for f in FILES:
    N,A,Bm=load(f); S,K=mats(N,A,Bm); S=S.astype(float); K=K.astype(float)
    Q=sum(np.linalg.norm(a,2)*np.linalg.norm(b,2) for a,b in zip(A,Bm))
    report(f.split('/')[-1]+f" q={len(A)}",N,K,S,round(Q,3))
