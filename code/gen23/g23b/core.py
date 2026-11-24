# G23b independent core: array = (N, col[a][b] int, pal: dict colour->'S'/'T'/'U', sg[a][b] in {+1,-1})
import numpy as np, itertools, math

def classes(N, col):
    M = {}
    for a in range(N):
        for b in range(N): M.setdefault(col[a][b], []).append((a, b))
    return M

def identity_tensor_ok(N, col, pal, sg):
    """Independent check of (*) via matrix products on every basis W=E_uv (doubled to stay integral)."""
    M = classes(N, col)
    K = {}
    for k, cells in M.items():
        A = np.zeros((N, N), dtype=np.int64)
        for (a, b) in cells: A[a, b] = sg[a][b]
        K[k] = A
    for u in range(N):
        for v in range(N):
            W = np.zeros((N, N), dtype=np.int64); W[u, v] = 1
            rhs = np.zeros((N, N), dtype=np.int64)
            for k, A in K.items():
                wc = {'S': 2, 'T': 0, 'U': 1}[pal[k]]; wr = 2 - wc
                rhs += wc * (A @ W @ A) + wr * int(np.trace(A @ W)) * A
            if not (rhs == 2 * W.T).all(): return False
    return True

def ranks(N, col, pal, sg):
    """exact ranks of R (flattening (x,y)|(u,v)) and C (flattening (x,u)|(y,v)), via sympy-free fraction rank (flint)."""
    import flint
    M = classes(N, col)
    R = np.zeros((N,)*4, dtype=np.int64); C = np.zeros((N,)*4, dtype=np.int64)
    for k, cells in M.items():
        wc = {'S': 2, 'T': 0, 'U': 1}[pal[k]]; wr = 2 - wc
        for (a, b) in cells:
            for (c, d) in cells:
                s = sg[a][b] * sg[c][d]
                # C[x,y,u,v] += wc K[x,u] K[v,y] : (x,u)=(a,b), (v,y)=(c,d)
                C[a, d, b, c] += wc * s
                # R[x,y,u,v] += wr K[x,y] K[v,u] : (x,y)=(a,b), (v,u)=(c,d)
                R[a, b, d, c] += wr * s
    X = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): X[x, y, y, x] = 2
    ok = (R + C == X).all()
    Rm = flint.fmpz_mat(N*N, N*N, [int(e) for e in R.reshape(-1)])
    Cm = flint.fmpz_mat(N*N, N*N, [int(e) for e in C.transpose(0, 2, 1, 3).reshape(-1)])
    return bool(ok), Rm.rank(), Cm.rank()

def predicate(N, col, pal, sg):
    """(i) every class a partial permutation; (ii) each 2x2 with a mono diagonal: other diagonal mono, palettes {S,T} or {U,U},
    product of 4 signs = -1."""
    for a in range(N):
        if len(set(col[a])) < N: return False
        if len(set(col[x][a] for x in range(N))) < N: return False
    for x, v in itertools.combinations(range(N), 2):
        for y, u in itertools.combinations(range(N), 2):
            d1 = col[x][y] == col[v][u]; d2 = col[x][u] == col[v][y]
            if d1 != d2: return False
            if d1:
                P = {pal[col[x][y]], pal[col[x][u]]}
                if not (P == {'S', 'T'} or (P == {'U'})): return False
                if sg[x][y] * sg[v][u] * sg[x][u] * sg[v][y] != -1: return False
    return True

def gstar(t, N):
    g = 1
    while g + 1 <= N and (g + 1) * g // 2 <= t: g += 1
    return g

def key_lemma_check(N, col, pal):
    """For every S/T colour of size h: off-diagonal cells of its h x h submatrix are of the opposite palette, each opposite colour
    occupies exactly one symmetric pair there; returns (ok, details)."""
    M = classes(N, col)
    nS = sum(1 for k in M if pal[k] == 'S'); nT = sum(1 for k in M if pal[k] == 'T')
    for k, cells in M.items():
        if pal[k] == 'U': continue
        opp = 'S' if pal[k] == 'T' else 'T'
        h = len(cells); seen = {}
        for i in range(h):
            for j in range(h):
                if i == j: continue
                c = col[cells[i][0]][cells[j][1]]
                if pal[c] != opp: return False, ('offdiag palette', k)
                seen.setdefault(c, set()).add(frozenset((i, j)))
        if any(len(s) != 1 for s in seen.values()): return False, ('colour in two pairs', k)
        if len(seen) != h * (h - 1) // 2: return False, ('count', k)
        if (nS if pal[k] == 'T' else nT) < h * (h - 1) // 2: return False, ('palette count', k)
    return True, None

def load_json(path):
    import json
    d = json.load(open(path)); N = d['N']
    col = [[None]*N for _ in range(N)]; sg = [[0]*N for _ in range(N)]; pal = {}
    for a in range(N):
        for b in range(N):
            P, k, s = d['cells'][a][b]; key = (P, k); col[a][b] = key; pal[key] = P; sg[a][b] = s
    return N, col, pal, sg
