# Check the trace-duality inequality  N^2 = <F,R> + <F,C> <= ||R||_tr + ||Fl2 C||_tr <= p||R||_op + q||Fl2 C||_op
# on a monomial ST array (STU convention of g22b, U empty): C(W)=sum_S K W K, R(W)=sum_T tr(K W) K.
import json, sys, numpy as np
d = json.load(open(sys.argv[1])); N = d['N']; cells = d['cells']
cls = {}
for a in range(N):
    for b in range(N):
        pal, k, s = cells[a][b]; cls.setdefault((pal,k), np.zeros((N,N)))[a,b] = s
X4 = lambda f: f  # placeholder
# build 4-tensors X[x,y,u,v]: map W -> sum over (u,v) X[x,y,u,v] W[u,v]
C = np.zeros((N,N,N,N)); R = np.zeros((N,N,N,N))
for (pal,k),K in cls.items():
    if pal == 'S': C += np.einsum('xu,vy->xyuv', K, K)      # (K W K)[x,y] = sum K[x,u] W[u,v] K[v,y]
    else:          R += np.einsum('vu,xy->xyuv', K, K)      # tr(K W) K[x,y] = sum K[v,u] W[u,v] K[x,y]
F = np.zeros((N,N,N,N))
for x in range(N):
    for y in range(N): F[x,y,y,x] = 1                      # W -> W^T : (W^T)[x,y] = W[y,x] => u=y, v=x
print("identity error:", np.abs(F - R - C).max())
Fl1 = lambda X: X.reshape(N*N, N*N)                         # (x,y)|(u,v)
Fl2 = lambda X: X.transpose(0,2,1,3).reshape(N*N, N*N)      # (x,u)|(y,v)
sR = np.linalg.svd(Fl1(R), compute_uv=False); sC = np.linalg.svd(Fl2(C), compute_uv=False)
p = (sR > 1e-9).sum(); q = (sC > 1e-9).sum()
print("p,q =", p, q, " <F,R>, <F,C> =", (F*R).sum(), (F*C).sum())
print("trace norms:", sR.sum(), sC.sum(), " p*op + q*op =", p*sR[0] + q*sC[0], " N^2 =", N*N)
print("Frobenius form sqrt(p)|R|+sqrt(q)|C| =", np.sqrt(p)*np.linalg.norm(R) + np.sqrt(q)*np.linalg.norm(C))
