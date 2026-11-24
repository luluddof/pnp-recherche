# Lemma (Schur averaging): max_{G,G'} rank A (G(x)G') B >= rank(A) rank(B) / N^2.  Random sparse integer A, B (float ranks).
import numpy as np, sys
rng = np.random.default_rng(5)
rk = lambda M: np.linalg.matrix_rank(M, tol=1e-8)
for N in [2,3,4]:
    D=N*N; worst = 9
    for t in range(3000 if N<4 else 800):
        a = rng.integers(1,D+1); b = rng.integers(1,D+1); d = rng.uniform(0.05,0.6)
        # structured: sometimes Kronecker-structured supports
        if rng.random()<0.5:
            A = np.kron(rng.integers(-1,2,size=(rng.integers(1,N+1),N))*(rng.random((1,N))<0.9), rng.integers(-1,2,size=(rng.integers(1,N+1),N)))
            B = np.kron(rng.integers(-1,2,size=(N,rng.integers(1,N+1))), rng.integers(-1,2,size=(N,rng.integers(1,N+1))))
            if rng.random()<0.5: A = A @ np.kron(np.eye(N), np.eye(N))[:, rng.permutation(D)]
        else:
            A = rng.integers(-1,2,size=(a,D))*(rng.random((a,D))<d); B = rng.integers(-1,2,size=(D,b))*(rng.random((D,b))<d)
        ra, rb = rk(A), rk(B)
        if ra*rb == 0: continue
        s = max(rk(A@np.kron(rng.standard_normal((N,N)),rng.standard_normal((N,N)))@B) for _ in range(2))
        worst = min(worst, s*D/(ra*rb))
    print('N',N,' min over samples of s N^2/(rank A rank B) = %.3f (lemma says >= 1)'%worst, flush=True)
