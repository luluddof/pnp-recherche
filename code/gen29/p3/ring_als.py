# General ring (periodic MPO) ALS for P_3: minimise ||(sigma - C) V||, V = n x (n-r) agreement subspace, C ring with
# bonds (b12,b23,b31); node tensors A1[l31,i12] (N x N op), A2[i12,j23], A3[j23,l31].
import numpy as np, sys, itertools
N = int(sys.argv[1]); b12, b23, b31 = map(int, sys.argv[2].split(',')); rs = list(map(int, sys.argv[3].split(',')))
starts = int(sys.argv[4]) if len(sys.argv) > 4 else 5; iters = int(sys.argv[5]) if len(sys.argv) > 5 else 300
rng = np.random.default_rng(int(sys.argv[6]) if len(sys.argv) > 6 else 0)
n = N**3
S = np.zeros((n, n))
for t in itertools.product(range(N), repeat=3):
    S[(t[2]*N + t[0])*N + t[1], (t[0]*N + t[1])*N + t[2]] = 1
def build(T1, T2, T3):
    # C = sum_{l,i,j} T1[l,i] (x) T2[i,j] (x) T3[j,l]
    C = np.zeros((n, n), complex)
    for l in range(b31):
        for i in range(b12):
            for j in range(b23):
                C += np.kron(np.kron(T1[l, i], T2[i, j]), T3[j, l])
    return C
def design(node, T1, T2, T3, V):
    # columns: derivative of C V w.r.t. each entry of the node tensor
    cols = []
    if node == 0:
        for l in range(b31):
            for i in range(b12):
                R = sum(np.kron(T2[i, j], T3[j, l]) for j in range(b23))
                for a in range(N):
                    for c in range(N):
                        E = np.zeros((N, N)); E[a, c] = 1
                        cols.append((np.kron(E, R) @ V).ravel())
    elif node == 1:
        for i in range(b12):
            for j in range(b23):
                for a in range(N):
                    for c in range(N):
                        E = np.zeros((N, N)); E[a, c] = 1
                        Cc = sum(np.kron(np.kron(T1[l, i], E), T3[j, l]) for l in range(b31))
                        cols.append((Cc @ V).ravel())
    else:
        for j in range(b23):
            for l in range(b31):
                for a in range(N):
                    for c in range(N):
                        E = np.zeros((N, N)); E[a, c] = 1
                        Cc = sum(np.kron(np.kron(T1[l, i], T2[i, j]), E) for i in range(b12))
                        cols.append((Cc @ V).ravel())
    return np.stack(cols, 1)
for r in rs:
    best = np.inf
    for st in range(starts):
        T1 = rng.standard_normal((b31, b12, N, N)) + 1j*rng.standard_normal((b31, b12, N, N))
        T2 = rng.standard_normal((b12, b23, N, N)) + 1j*rng.standard_normal((b12, b23, N, N))
        T3 = rng.standard_normal((b23, b31, N, N)) + 1j*rng.standard_normal((b23, b31, N, N))
        for it in range(iters):
            C = build(T1, T2, T3); s, Vh = np.linalg.svd(S - C)[1:]
            V = Vh.conj().T[:, r:]; tail = np.sum(s[r:]**2)
            if tail < 1e-20: break
            for node in range(3):
                Dm = design(node, T1, T2, T3, V); bvec = (S @ V).ravel()
                x, *_ = np.linalg.lstsq(Dm, bvec, rcond=None)
                if node == 0: T1 = x.reshape(b31, b12, N, N)
                elif node == 1: T2 = x.reshape(b12, b23, N, N)
                else: T3 = x.reshape(b23, b31, N, N)
        s = np.linalg.svd(S - build(T1, T2, T3), compute_uv=False); tail = np.sum(s[r:]**2)
        best = min(best, tail)
        print(f'N={N} bonds=({b12},{b23},{b31}) r={r} start={st} tail={tail:.2e} s[r-1:r+2]={np.round(s[max(r-1,0):r+2],4)}', flush=True)
    print(f'== N={N} bonds=({b12},{b23},{b31}) r={r} BEST {best:.2e}', flush=True)
