# embed an exact json (N0) into N >= N0 with q >= q0 terms (extra terms random small) + noise. usage: warmEmb.py json N q seed out
import sys, json, numpy as np
sys.path.insert(0,'/home/user/pnp-recherche/code/gen26/p3')
from load import load
f,N,q,seed,out=sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),sys.argv[5]; rng=np.random.default_rng(seed)
N0,A,B=load(f); G=np.zeros((q,N,N)); H=np.zeros((q,N,N))
for l in range(len(A)): G[l,:N0,:N0]=A[l]; H[l,:N0,:N0]=B[l]
for l in range(len(A),q): G[l]=rng.standard_normal((N,N))*0.3; H[l]=rng.standard_normal((N,N))*0.3
z=np.concatenate([G.ravel(),H.ravel()])+0.03*rng.standard_normal(2*q*N*N); np.save(out,z)
