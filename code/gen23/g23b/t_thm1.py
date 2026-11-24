# Exhaustive test of Theorem 1 (identity <=> predicate) using the GENERAL definition C = sum G (x) G', R = sum H (x) H'
# (tensor built from scratch), vectorised over all sign vectors.
import numpy as np, itertools, sys, random, time
from core import predicate

def partitions(n, pp_only, N):
    """restricted growth strings over n = N*N cells (reading order); if pp_only, classes must be partial permutations."""
    rg = [0]*n
    def rec(i, m):
        if i == n: yield list(rg); return
        a, b = divmod(i, N)
        for c in range(m + 1):
            if pp_only and c < m:
                bad = False
                for j in range(i):
                    if rg[j] == c:
                        a2, b2 = divmod(j, N)
                        if a2 == a or b2 == b: bad = True; break
                if bad: continue
            rg[i] = c
            yield from rec(i + 1, max(m, c + 1))
    yield from rec(0, 0)

def tensors_all_signs(N, rg, pal, signs):
    """signs: (B, N*N) array of +-1. Returns (B,) bool: identity 2X == 2C + 2R."""
    n = N * N; B = signs.shape[0]
    T = np.zeros((B, N, N, N, N), dtype=np.int64)
    for c1 in range(n):
        for c2 in range(n):
            if rg[c1] != rg[c2]: continue
            k = rg[c1]; wc = {'S': 2, 'T': 0, 'U': 1}[pal[k]]; wr = 2 - wc
            a, b = divmod(c1, N); c, d = divmod(c2, N)
            s = signs[:, c1] * signs[:, c2]
            # C(x,y,u,v) = sum_l G_l[x,u] G'_l[y,v], G_l = w K, G'_l = K^T:  K[x,u] K[v,y]  -> x=a,u=b,v=c,y=d
            if wc: T[:, a, d, b, c] += wc * s
            # R(x,y,u,v) = sum_j H_j[x,y] H'_j[u,v], H' = K^T: K[x,y] K[v,u] -> x=a,y=b,v=c,u=d
            if wr: T[:, a, b, d, c] += wr * s
    X = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): X[x, y, y, x] = 2
    return (T == X).reshape(B, -1).all(axis=1)

def pred_nosign(N, col, pal):
    for a in range(N):
        if len(set(col[a])) < N or len(set(col[x][a] for x in range(N))) < N: return False, None
    quads = []
    for x, v in itertools.combinations(range(N), 2):
        for y, u in itertools.combinations(range(N), 2):
            d1 = col[x][y] == col[v][u]; d2 = col[x][u] == col[v][y]
            if d1 != d2: return False, None
            if d1:
                P = {pal[col[x][y]], pal[col[x][u]]}
                if not (P == {'S', 'T'} or P == {'U'}): return False, None
                quads.append((x*N+y, v*N+u, x*N+u, v*N+y))
    return True, quads

def run(N, pp_only, sample=None, seed=1):
    n = N * N
    allsigns = np.array(list(itertools.product([1, -1], repeat=n)), dtype=np.int64) if n <= 9 else None
    rnd = random.Random(seed)
    stats = dict(structs=0, id_true=0, mismatch=0)
    for rg in partitions(n, pp_only, N):
        m = max(rg) + 1
        for palt in itertools.product('STU', repeat=m):
            if sample is not None and rnd.random() > sample: continue
            pal = dict(enumerate(palt)); stats['structs'] += 1
            sg = allsigns if allsigns is not None else np.array([[rnd.choice([1, -1]) for _ in range(n)] for _ in range(64)])
            ids = tensors_all_signs(N, rg, pal, sg)
            col = [[rg[a*N+b] for b in range(N)] for a in range(N)]
            ok0, quads = pred_nosign(N, col, pal)
            if not ok0:
                pr = np.zeros(sg.shape[0], dtype=bool)
            else:
                pr = np.ones(sg.shape[0], dtype=bool)
                for q in quads: pr &= (sg[:, q[0]] * sg[:, q[1]] * sg[:, q[2]] * sg[:, q[3]] == -1)
            # spot-check the vectorised predicate against core.predicate on a few sign vectors
            for t in rnd.sample(range(sg.shape[0]), 2):
                s2 = [[int(sg[t, a*N+b]) for b in range(N)] for a in range(N)]
                assert predicate(N, col, pal, s2) == bool(pr[t])
            mm = int((ids != pr).sum())
            if mm: stats['mismatch'] += mm; print('MISMATCH', N, rg, palt)
            stats['id_true'] += int(ids.sum())
    return stats

if __name__ == "__main__":
    t = time.time()
    print('N=2 all partitions (incl. non-pp), all palettes, all 16 signs:', run(2, False), f'{time.time()-t:.1f}s', flush=True)
    t = time.time()
    print('N=3 all partial-permutation partitions, all palettes, all 512 signs:', run(3, True), f'{time.time()-t:.1f}s', flush=True)
    t = time.time()
    print('N=3 all partitions (incl. non-pp), 2% sample of (partition,palette), all 512 signs:', run(3, False, sample=0.02), f'{time.time()-t:.1f}s', flush=True)
