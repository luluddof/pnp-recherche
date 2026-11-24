# Exhaustive enumeration of ST arrays (no U) for small N, up to renaming of colours (canonical first-occurrence labels),
# all palette choices. Sign-free intercalate condition checked incrementally; signs solved over GF(2) afterwards.
# For every sign-free-valid structure: key lemma + counting inequality checked; frontier computed for signable ones.
import sys, itertools, time
from core import gstar, key_lemma_check, predicate, identity_tensor_ok

def gf2_solvable(nvars, eqs):
    """eqs: list of (list of var indices, rhs bit). Gaussian elimination over GF(2) with bitmasks."""
    rows = []
    for vs, r in eqs:
        m = 0
        for v in vs: m ^= 1 << v
        rows.append((m, r))
    piv = {}
    sol_rows = []
    for m, r in rows:
        for p, (pm, pr) in piv.items():
            if m >> p & 1: m ^= pm; r ^= pr
        if m == 0:
            if r: return None
            continue
        p = m.bit_length() - 1
        # eliminate p from existing pivots
        for q in list(piv):
            qm, qr = piv[q]
            if qm >> p & 1: piv[q] = (qm ^ m, qr ^ r)
        piv[p] = (m, r)
    # build a solution: free vars 0, pivot var p = r xor (other vars in row, all free -> 0)
    x = [0] * nvars
    for p, (m, r) in piv.items(): x[p] = r
    return x

def enumerate_st(N, report):
    n = N * N
    col = [[-1] * N for _ in range(N)]
    pal = []
    out = []
    def ok_cell(a, b):
        c = col[a][b]
        for x in range(a):
            for y in range(N):
                if y == b: continue
                # 2x2 rows x<a, cols {y,b}; all four assigned iff (a, y) assigned or y < b ... (a,y) assigned iff y < b
                if y > b: continue
                d1 = col[x][y] == c          # diagonal (x,y),(a,b)
                d2 = col[x][b] == col[a][y]  # antidiagonal
                if d1 != d2: return False
                if d1 and pal[col[x][y]] == pal[col[x][b]]: return False
        return True
    def rec(i):
        if i == n:
            out.append(([r[:] for r in col], pal[:])); return
        a, b = divmod(i, N)
        used_row = set(col[a][:b]); used_col = set(col[x][b] for x in range(a))
        for c in range(len(pal) + 1):
            if c < len(pal):
                if c in used_row or c in used_col: continue
                col[a][b] = c
                if ok_cell(a, b): rec(i + 1)
            else:
                for P in 'ST':
                    pal.append(P); col[a][b] = c
                    if ok_cell(a, b): rec(i + 1)
                    pal.pop()
        col[a][b] = -1
    rec(0)
    return out

if __name__ == "__main__":
    N = int(sys.argv[1])
    t = time.time(); structs = enumerate_st(N, None)
    print(f"N={N}: sign-free valid ST structures (canonical labels): {len(structs)}  [{time.time()-t:.1f}s]", flush=True)
    front = {}; nsign = 0; lemma_bad = 0; ineq_bad = 0; tight = {}
    for col, palL in structs:
        pal = dict(enumerate(palL))
        q = palL.count('S'); p = palL.count('T')
        ok, info = key_lemma_check(N, col, pal)
        if not ok: lemma_bad += 1; print('LEMMA FAIL', col, palL, info)
        slack = q * gstar(p, N) + p * gstar(q, N) - N * N
        if slack < 0: ineq_bad += 1; print('INEQ FAIL', col, palL)
        # signs
        eqs = []
        for x, v in itertools.combinations(range(N), 2):
            for y, u in itertools.combinations(range(N), 2):
                if col[x][y] == col[v][u]:
                    eqs.append(([x*N+y, v*N+u, x*N+u, v*N+y], 1))   # product -1 <=> xor of 'minus' bits = 1
        sol = gf2_solvable(N*N, eqs)
        if sol is None: continue
        nsign += 1
        sg = [[-1 if sol[a*N+b] else 1 for b in range(N)] for a in range(N)]
        assert predicate(N, col, pal, sg)
        if nsign % 50 == 1: assert identity_tensor_ok(N, col, pal, sg)
        front[q] = min(front.get(q, 10**9), p)
        if slack == 0: tight[(p, q)] = tight.get((p, q), 0) + 1
    print(f"  signable: {nsign}; key-lemma failures: {lemma_bad}; counting-inequality failures: {ineq_bad}")
    print("  exact ST frontier q -> min p:", dict(sorted(front.items())))
    bnd = {}
    for q in sorted(front):
        p = 0
        while N*N > q*gstar(p, N) + p*gstar(q, N): p += 1
        bnd[q] = p
    print("  counting bound   q -> p_min:", bnd)
    print("  tight (slack 0) signable arrays by (p,q):", tight)
