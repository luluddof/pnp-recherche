"""Exhaustive, exact over GF(2), N = 4:  m* = min over F, G in GF(2)^{4x4} of rank( X - f(x,u) g(y,v) ) in the ROW flattening
(= rank(I - (F(x)G) Swap)).  WLOG F = diag(I_r, 0) (W -> A W B^T maps (F,G) -> (A F B^-1, B G A^-1)); G runs over all 2^16.
Rows of the 16x16 matrix are bit masks; rank by elimination on ints."""
import itertools, sys
N = 4
def rank_masks(rows):
    rows = list(rows); r = 0
    for bit in range(16):
        piv = None
        for i in range(r, 16):
            if rows[i] >> bit & 1: piv = i; break
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        for i in range(16):
            if i != r and rows[i] >> bit & 1: rows[i] ^= rows[r]
        r += 1
    return r
best = {}
for r in range(N + 1):
    F = [[1 if (i == j and i < r) else 0 for j in range(N)] for i in range(N)]
    mn = 99; cnt = 0
    for gbits in range(1 << 16):
        G = [[(gbits >> (4 * i + j)) & 1 for j in range(N)] for i in range(N)]
        rows = []
        for x in range(N):
            for y in range(N):
                m = 0
                for u in range(N):
                    if not F[x][u] and not (x == 0 and False): pass
                    for v in range(N):
                        val = (1 if (x == v and y == u) else 0) ^ (F[x][u] & G[y][v])
                        if val: m |= 1 << (4 * u + v)
                rows.append(m)
        k = rank_masks(rows)
        if k < mn: mn = k; cnt = 1
        elif k == mn: cnt += 1
    best[r] = (mn, cnt)
    print(f"rank F = {r}: min row rank of X - f g = {mn} (attained by {cnt} matrices G)"); sys.stdout.flush()
print("m* =", min(v[0] for v in best.values()))
