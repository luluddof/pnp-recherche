import numpy as np, itertools
exec(open('c3_check.py').read().split('for n, k in')[0])
rng = np.random.default_rng(12)
for n, k in [(2,3), (2,4), (3,3)]:
    N = n**k; Sw = swap(n); I = np.eye(n*n); Pa = (I - Sw)/2; Ps = (I + Sw)/2
    pts = list(itertools.product([0,1], repeat=k)); S = site_order([Sw]*k, n, k)
    for trial in range(3):
        g = {T: ((-1)**sum(T) if rng.random() < 0.6 else rng.normal()) for T in pts}
        K = sum(g[T]*site_order([Pa if T[i] else Ps for i in range(k)], n, k) for T in pts)
        p = np.linalg.matrix_rank(S - K, tol=1e-8)
        mult = {T: np.prod([n*(n-1)/2 if t else n*(n+1)/2 for t in T]) for T in pts}
        pform = sum(mult[T] for T in pts if abs(g[T] - (-1)**sum(T)) > 1e-9)
        # coefficients c_U in tau basis: tau_i = a - 2T_i with a = 1-1/n ; solve linear system
        Us = [U for r in range(k+1) for U in itertools.combinations(range(k), r)]
        Amat = np.array([[np.prod([1-1/n-2*T[i] for i in U]) for U in Us] for T in pts])
        c = np.linalg.solve(Amat, np.array([g[T] for T in pts]))
        qform = sum((n*n-1)**len(U) for U, cc in zip(Us, c) if abs(cc) > 1e-9)
        C = K.reshape(N,N,N,N).transpose(0,2,1,3).reshape(N*N,N*N)
        q = np.linalg.matrix_rank(C, tol=1e-8)
        print(f"n={n} k={k}: p={p} formula {int(pform)}; osr={q} formula {qform}")
