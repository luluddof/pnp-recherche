import sys, glob, numpy as np
sys.path.insert(0,'/home/user/pnp-recherche/code/gen26/p3')
from load import load, mats, FILES
from collections import Counter
fs=FILES+sorted(glob.glob('exact_*.json'))+['q8_8_8_8.json']
print("file N p q | m #ev(Ts)>=1 #ev>1 #ev<0 #ev=0 | ||Ts|| ||K|| Q trTheta symTheta | Theta-eigs(non-1)")
for f in fs:
    N,A,B=load(f); S,K=mats(N,A,B); T=(K@S).astype(float); n=N*N; p=np.linalg.matrix_rank(np.eye(n)-T); Ts=(T+T.T)/2; ev=np.linalg.eigvalsh(Ts)
    Q=sum(np.linalg.norm(a,2)*np.linalg.norm(b,2) for a,b in zip(A,B)); e=np.linalg.eigvals(T)
    c=Counter([round(x.real,2) for x in e if abs(x-1)>1e-6])
    print(f"{f.split('/')[-1]:18s} {N} {p:2d} {len(A)} | {n-p:2d} {int(np.sum(ev>=1-1e-9)):2d} {int(np.sum(ev>1+1e-9))} {int(np.sum(ev<-1e-9)):2d} {int(np.sum(abs(ev)<1e-9))} | {abs(ev).max():.2f} {np.linalg.norm(K,2):.2f} {Q:.1f} {T.trace():.0f} {bool((T==T.T).all())} | {dict(sorted(c.items()))}")
