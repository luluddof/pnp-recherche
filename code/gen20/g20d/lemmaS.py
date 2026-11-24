# Lemma S test: single compression Phi(X) = A X B (c1 = 1). Is rank(AFB)^2 <= N^2 * s, s = generic rank A(G(x)G')B ?
import numpy as np, sys
N = int(sys.argv[1]); trials = int(sys.argv[2]); rng = np.random.default_rng(int(sys.argv[3]))
D = N*N
F = np.zeros((D,D))
for x in range(N):
    for y in range(N): F[x*N+y, y*N+x] = 1
rk = lambda M: np.linalg.matrix_rank(M, tol=1e-8)
best = {}
for t in range(trials):
    r = rng.integers(1, D+1)
    dens = rng.uniform(0.1, 0.9)
    A = rng.integers(-1,2,size=(r,D)) * (rng.random((r,D))<dens)
    B = rng.integers(-1,2,size=(D,r)) * (rng.random((D,r))<dens)
    R = rk(A@F@B)
    if R == 0: continue
    s = max(rk(A@np.kron(rng.standard_normal((N,N)), rng.standard_normal((N,N)))@B) for _ in range(2))
    ratio = R*R/(D*s)
    if ratio > best.get(R, (0,))[0]: best[R] = (ratio, s)
for R in sorted(best): print('N',N,'R',R,'max ratio R^2/(N^2 s) = %.3f'%best[R][0],'s =',best[R][1])
