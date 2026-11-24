# dim J(p,q) by Terracini (rank of the differential of (a,b,G,H) -> A+B at a random point, mod PR), and tangent checks.
import numpy as np, sys
from common import *
def jac(N, p, q, rng):
    a = rnd((p,N,N), rng); b = rnd((p,N,N), rng); G = rnd((q,N,N), rng); H = rnd((q,N,N), rng)
    rows = []
    E = np.eye(N, dtype=np.int64)
    for k in range(p):
        for i in range(N):
            for j in range(N):
                e = np.zeros((N,N), dtype=np.int64); e[i,j] = 1
                rows.append(np.einsum('xy,uv->xyuv', e, b[k]).ravel()); rows.append(np.einsum('xy,uv->xyuv', a[k], e).ravel())
    for k in range(q):
        for i in range(N):
            for j in range(N):
                e = np.zeros((N,N), dtype=np.int64); e[i,j] = 1
                rows.append(np.einsum('xu,yv->xyuv', e, H[k]).ravel()); rows.append(np.einsum('xu,yv->xyuv', G[k], e).ravel())
    return np.array(rows) % PR
if __name__ == "__main__":
    rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
    for N,p,q in [(3,2,2),(3,3,1),(3,1,3),(3,2,1),(3,3,3),(4,3,3),(4,4,3),(4,3,4),(4,2,2),(4,4,2),(4,5,2)]:
        d = max(rank(jac(N,p,q,rng)) for _ in range(2))
        exp = min(N**4, (2*p*N*N-p*p)+(2*q*N*N-q*q))
        print(f"N={N} (p,q)=({p},{q}) dim J = {d}  expected {exp}  codim {N**4-d}", flush=True)
