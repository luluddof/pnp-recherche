# warm start: Q8 ST array (8,8)@8 (g23a/sol_Q8_8_8_8.json): S-classes A_k give C = sum_k A_k (x) A_k^T (ST convention, X*2 = R + C with weight 2
# on S colours => G_k = A_k, G'_k = A_k^T after dividing by 2 on both sides).  Embed in N >= 8 with q >= 8 terms (+ small noise terms).
import json, ast, numpy as np, sys
from vp11 import swap
Nt,q,seed=int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]); rng=np.random.default_rng(seed)
d=json.load(open('/home/user/pnp-recherche/code/gen23/g23a/sol_Q8_8_8_8.json')); cells=ast.literal_eval(d['cells']) if isinstance(d['cells'],str) else d['cells']
N=8; cls={}
for a in range(N):
    for b in range(N):
        P,k,s=cells[a][b]; cls.setdefault((P,k),[]).append((a,b,s))
S_cls=[v for k,v in cls.items() if k[0]=='S']
G=[];H=[]
for v in S_cls:
    A=np.zeros((N,N))
    for a,b,s in v: A[a,b]=s
    G.append(A); H.append(A.T.copy())
s=N*N; K=sum(np.kron(g,h) for g,h in zip(G,H)); Sm=swap(N).reshape(s,s)
print("Q8 check: #S classes",len(G),"rank(Swap-K) =",np.linalg.matrix_rank(Sm-K))
Gt=np.zeros((q,Nt,Nt)); Ht=np.zeros((q,Nt,Nt))
for l in range(min(q,len(G))): Gt[l,:N,:N]=G[l]; Ht[l,:N,:N]=H[l]
z=np.concatenate([Gt.ravel(),Ht.ravel()])+0.05*rng.standard_normal(2*q*Nt*Nt)
np.save(f"warmQ8_{Nt}_{q}_s{seed}.npy",z)
