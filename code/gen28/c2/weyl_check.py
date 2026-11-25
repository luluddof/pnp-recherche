# Qudit Weyl "Bell-diagonal" models for odd d: K in span{ W_u (x) W_{Mu} } with M an anti-symplectic involution
# (M(c,e) = (c,-e) <-> W_u (x) conj(W_u) up to phases).  Claim: p >= (N^2 - N)/2 for every such K.
import numpy as np, itertools, sys
rng = np.random.default_rng(3)
def weyl(d, c, e):
    w = np.exp(2j*np.pi/d); Z = np.diag([w**j for j in range(d)]); X = np.roll(np.eye(d), 1, axis=0)
    return np.linalg.matrix_power(Z, c) @ np.linalg.matrix_power(X, e)
def kl(ms):
    o = np.array([[1.0+0j]])
    for m in ms: o = np.kron(o, m)
    return o
def swap(N):
    S = np.zeros((N*N, N*N))
    for x in range(N):
        for y in range(N): S[y*N+x, x*N+y] = 1
    return S
def rank(M, tol=1e-8):
    s = np.linalg.svd(M, compute_uv=False); return int((s > tol*max(1, s[0])).sum())
for d, k in [(3,1), (3,2), (5,1)]:
    N = d**k; S = swap(N); labels = list(itertools.product(range(d), repeat=2*k))
    ops = []
    for lab in labels:
        cs, es = lab[:k], lab[k:]
        A = kl([weyl(d, c, e) for c, e in zip(cs, es)])
        Bm = kl([weyl(d, c, (-e) % d) for c, e in zip(cs, es)])  # M(c,e) = (c,-e)
        ops.append(np.kron(A, Bm))
    # commutation check
    comm = max(np.abs(ops[i] @ ops[j] - ops[j] @ ops[i]).max() for i in range(0, len(ops), 3) for j in range(0, len(ops), 5))
    best = N*N
    for t in range(300):
        q = rng.integers(1, len(ops)+1)
        sel = rng.choice(len(ops), size=q, replace=False)
        K = sum((rng.normal()+1j*rng.normal())*ops[i] for i in sel)
        best = min(best, rank(S - K))
    # also least-squares best K in the full span (q = d^{2k}): project SWAP onto span
    Mmat = np.array([o.ravel() for o in ops]).T
    coef, *_ = np.linalg.lstsq(Mmat, S.ravel().astype(complex), rcond=None)
    Kls = (Mmat @ coef).reshape(S.shape)
    # exact optimum over the span: eigenbasis is common; compute SWAP in that eigenbasis
    w, V = np.linalg.eig(sum(rng.normal()*o for o in ops))
    Sb = np.linalg.solve(V, S @ V)
    offdiag = np.abs(Sb - np.diag(np.diag(Sb))).max()
    print(f"d={d} k={k} N={N}: max commutator {comm:.1e}; random K min p = {best}; LS-projection p = {rank(S-Kls)};"
          f" (N^2-N)/2 = {(N*N-N)//2}; SWAP off-diagonal in common eigenbasis (max) = {offdiag:.2f}")
