# [E] tensor family: K = A_k (x) (I_m (x) I_m) on N = k m, A_k = Swap_k - Phi Phi^T (osr k(k-1)/2, rank(Swap_k - A_k) = 1).
# Predicted rank(Swap_N - K) = m(m+1)/2 * 1 + m(m-1)/2 * rank(2 Swap_k - Phi Phi^T)  (= k^2 unless k = 2 in the field).
# Also K = Swap_k (x) (I_m (x) I_m): osr k^2, rank = k^2 m(m-1)/2 = N(N-k)/2.
import numpy as np, flint, itertools
fz = lambda A: flint.fmpz_mat([[int(v) for v in r] for r in A])
def swap4(n):
    X = np.zeros((n,)*4, int)
    for x in range(n):
        for y in range(n): X[x,y,y,x] = 1
    return X
def build(A4, B4, k, m):
    # A4[x,y,u,v] on K^k, B4 on K^m; combined index (x,x') with x in K^k, x' in K^m
    T = np.einsum('abcd,efgh->aebfcgdh', A4, B4)   # [x,x',y,y',u,u',v,v']
    N = k*m; return T.reshape(N, N, N, N)
for k, m in [(3,2), (3,3), (2,3), (4,2), (2,4)]:
    N = k*m; Sk = swap4(k); Phi = np.einsum('xy,uv->xyuv', np.eye(k, dtype=int), np.eye(k, dtype=int))
    Ak = Sk - Phi; Im = np.einsum('xu,yv->xyuv', np.eye(m, dtype=int), np.eye(m, dtype=int))
    for name, A in [("A_k", Ak), ("Swap_k", Sk)]:
        K = build(A, Im, k, m); SN = swap4(N)
        rk = fz((SN - K).reshape(N*N, N*N)).rank(); osr = fz(K.transpose(0,2,1,3).reshape(N*N, N*N)).rank()
        print(f"k={k} m={m} N={N} K={name}(x)I: rank(Swap-K)={rk}, osr={osr}, N(N-1)/2={N*(N-1)//2}, gain={N*(N-1)//2-rk}")
