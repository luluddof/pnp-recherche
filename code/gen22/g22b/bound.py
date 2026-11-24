# Theorem 3 (ST case): N^2 <= q*g*(p) + p*g*(q), g*(t) = max{g <= N : g(g-1)/2 <= t}; and the lemma checked on arrays.
import json, glob, itertools, math
def gstar(t, N):
    g = 1
    while g + 1 <= N and (g + 1) * g // 2 <= t: g += 1
    return g
def pmin(N, q):
    p = 0
    while N*N > q*gstar(p, N) + p*gstar(q, N): p += 1
    return p
def check_lemma(N, cells):
    cls = {}
    for a in range(N):
        for b in range(N): cls.setdefault((cells[a][b][0], cells[a][b][1]), []).append((a, b))
    nS = sum(1 for k in cls if k[0] == 'S'); nT = sum(1 for k in cls if k[0] == 'T')
    ok = True
    for k, M in cls.items():
        if k[0] == 'U': continue
        h = len(M); rows = [a for a, b in M]; cols = [b for a, b in M]
        opp = 'S' if k[0] == 'T' else 'T'
        seen = {}
        for i, a in enumerate(rows):
            for j, b in enumerate(cols):
                if i == j: continue
                c = (cells[a][b][0], cells[a][b][1])
                if c[0] != opp: ok = False
                seen.setdefault(c, set()).add((min(i, j), max(i, j)))
        if any(len(v) != 1 for v in seen.values()) or len(seen) != h*(h-1)//2: ok = False
        if (nS if k[0] == 'T' else nT) < h*(h-1)//2: ok = False
    return ok, nS, nT, sorted(len(v) for v in cls.values())
if __name__ == "__main__":
    for N in range(2, 13):
        print(f"N={N}: bound p_min(q) for q=1..: ", [pmin(N, q) for q in range(1, 13)], " balanced >=", next(p for p in range(1, N*N) if pmin(N, p) <= p))
    for N in (16, 32, 64, 128):
        pb = next(p for p in range(1, N*N) if pmin(N, p) <= p)
        print(f"N={N}: balanced ST bound p >= {pb};  N^(4/3)/2 = {N**(4/3)/2:.1f}")
    for f in sorted(glob.glob('shape_*.json') + glob.glob('st_*.json')):
        if not f.endswith('.json'): continue
        d = json.load(open(f)); N = d['N']; cells = [[tuple(c) for c in r] for r in d['cells']]
        if any(c[0] == 'U' for r in cells for c in r): print(f, 'has U colours, skip'); continue
        ok, nS, nT, sizes = check_lemma(N, cells)
        print(f, f"N={N} p={nT} q={nS}: lemma (off-diagonal of each class submatrix = opposite palette, symmetric pairs, distinct) {ok};",
              "bound slack", nS*gstar(nT, N) + nT*gstar(nS, N) - N*N, "sizes", sizes)
