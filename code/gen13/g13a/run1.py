import ncsearch, numpy as np, sys
l,d=3,1; tau=[1.,2.,3.]
for seed in range(3):
    f,B=ncsearch.solve(l,d,tau,(1,1,1),'PQ',restarts=60,seed=seed)
    print("PQ (1,1,1) seed",seed,"residual",f)
    if B is not None:
        print("  singular values per block:",[np.round(np.linalg.svd(b,compute_uv=False),3).tolist() for b in B])
