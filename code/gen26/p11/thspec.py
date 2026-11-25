import numpy as np, sys, re
from collections import Counter
from vp11 import swap
for f in sys.argv[1:]:
    m=re.search(r'_(\d+)_(\d+)_(\d+)_',f); N,p,q=map(int,m.groups()); z=np.load(f); s=N*N
    G=z[:q*s].reshape(q,N,N); H=z[q*s:].reshape(q,N,N); K=np.einsum('kxu,kyv->xyuv',G,H).reshape(s,s); Th=K@swap(N).reshape(s,s)
    ev=np.linalg.eigvals(Th); c=Counter([complex(round(e.real,2),round(e.imag,2)) for e in ev])
    Ts=(Th+Th.T)/2; print(f, f"N={N} p={p} q={q} trTheta={np.trace(Th):.3f} ||Theta-Theta^T||={np.linalg.norm(Th-Th.T):.2e}", sorted(c.items(),key=lambda t:-t[1])[:8])
