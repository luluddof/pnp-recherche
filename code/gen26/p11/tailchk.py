import numpy as np, sys
from vp11 import swap
f,N,p,q=sys.argv[1],*map(int,sys.argv[2:5]); z=np.load(f); s=N*N; G=z[:q*s].reshape(q,N,N); H=z[q*s:].reshape(q,N,N)
K=np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s); sv=np.linalg.svd(swap(N).reshape(s,s)-K,compute_uv=False)
print(f, 'sv around p:',np.round(sv[max(0,p-3):p+4],5),'tail',np.sum(sv[p:]**2))
