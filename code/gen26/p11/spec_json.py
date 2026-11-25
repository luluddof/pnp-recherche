import sys, json, numpy as np
sys.path.insert(0,'/home/user/pnp-recherche/code/gen26/p3')
from load import load, mats
from spec11 import report
for f in sys.argv[1:]:
    N,A,B=load(f); S,K=mats(N,A,B); Q=sum(np.linalg.norm(a,2)*np.linalg.norm(b,2) for a,b in zip(A,B))
    report(f+f" q={len(A)}",N,K.astype(float),S.astype(float),round(Q,3))
