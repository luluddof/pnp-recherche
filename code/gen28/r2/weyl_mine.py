import numpy as np, itertools
rng = np.random.default_rng(5)
def W(d, c, e):
    w = np.exp(2j*np.pi/d); X = np.roll(np.eye(d), 1, axis=0); Z = np.diag(w**np.arange(d))
    return np.linalg.matrix_power(X, c % d) @ np.linalg.matrix_power(Z, e % d)
def kl(ms):
    o = np.eye(1, dtype=complex)
    for m in ms: o = np.kron(o, m)
    return o
def symp(d, k, u, v):  # omega(u,v) = c.e' - e.c'
    return (np.dot(u[:k], v[k:]) - np.dot(u[k:], v[:k])) % d
def test(d, k, M):
    N = d**k; labs = list(itertools.product(range(d), repeat=2*k))
    # check anti-symplectic involution
    for u in labs[:20]:
        for v in labs[:20]:
            assert (symp(d,k,M@np.array(u)%d, M@np.array(v)%d) + symp(d,k,np.array(u),np.array(v))) % d == 0
    ops = []
    for u in labs:
        u = np.array(u); mu = M @ u % d
        ops.append(np.kron(kl([W(d,u[i],u[k+i]) for i in range(k)]), kl([W(d,mu[i],mu[k+i]) for i in range(k)])))
    G = sum(rng.normal()*o for o in ops)
    w, V = np.linalg.eig(G)
    S = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N): S[y*N+x, x*N+y] = 1
    Sb = np.linalg.solve(V, S @ V)
    A = np.abs(Sb) > 1e-6
    rowcnt = A.sum(1); fixed = sum(A[i,i] for i in range(N*N))
    # optimum over full span: diagonal K; 2-cycles each cost 1, fixed points 0
    print(f"d={d} k={k} N={N}: SWAP monomial in eigenbasis: {np.all(rowcnt==1)}, fixed points {fixed}, p_opt=(N^2-fixed)/2={(N*N-fixed)//2}, (N^2-N)/2={(N*N-N)//2}")
    # brute random check
    best = N*N
    for t in range(200):
        dvec = np.where(rng.random(N*N) < 0.5, np.diag(Sb), rng.normal(size=N*N))
        # choose rank-1 completions on 2-cycles
        K = np.diag(dvec).astype(complex)
        for i in range(N*N):
            j = np.argmax(np.abs(Sb[i]))
            if j > i:
                a = rng.normal()+1j*rng.normal(); K[i,i] = a; K[j,j] = Sb[i,j]*Sb[j,i]/a
        best = min(best, np.linalg.matrix_rank(Sb - K, tol=1e-8))
    print("   heuristic best rank:", best)
test(3, 1, np.diag([1, -1]))
test(3, 1, np.diag([-1, 1]))
test(5, 1, np.diag([1, -1]))
test(3, 2, np.diag([1, 1, -1, -1]))
M2 = np.array([[0,0,1,0],[0,0,0,1],[1,0,0,0],[0,1,0,0]])  # (c,e)->(e,c): omega(Mu,Mv) = -omega(u,v)
test(3, 2, M2)
