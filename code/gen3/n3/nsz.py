"""NS degree over Z, Q (exact integer HNF via flint) and GF(p) (nmod_mat rank).
Matrix: rows = generators g*p_C (g multilinear monomial on vars outside C),
cols = multilinear monomials of degree <= d, ordered degree-descending, constant LAST.
Row-HNF: e_const in Q-rowspace iff HNF has a pivot in the last column; then that
pivot value D is the order of e_const in Z^N/L, i.e. minimal D with D*1 = sum g_C p_C over Z.
Z-refutation iff D == 1."""
import sys, itertools, time, json
sys.path.insert(0, '/home/user/pnp-recherche/code/step3d')
from ns import rand3cnf, unsat, clause_poly, tseitin
import flint

def build(n, polys, d):
    mons = [sum(1 << i for i in S) for k in range(d, -1, -1) for S in itertools.combinations(range(n), k)]
    idx = {m: i for i, m in enumerate(mons)}
    rows = []
    for p in polys:
        vs = 0
        for m in p: vs |= m
        dp = max(bin(m).count('1') for m in p)
        free = [i for i in range(n) if not vs >> i & 1]
        for k in range(0, d - dp + 1):
            for S in itertools.combinations(free, k):
                g = sum(1 << i for i in S)
                rows.append({idx[m | g]: a for m, a in p.items()})
    return rows, len(mons)

def hnf_D(rows, nc, block=None):
    """Blocked row HNF over Z. Returns D (pivot in last col) or None (not in Q-span)."""
    block = block or max(nc, 200)
    H = None
    for s in range(0, len(rows), block):
        chunk = rows[s:s + block]
        M = [[0] * nc for _ in chunk]
        for r, row in zip(M, chunk):
            for c, a in row.items(): r[c] = a
        if H is not None: M = H + M
        Hm = flint.fmpz_mat(M).hnf()
        H = [r for r in Hm.tolist() if any(r)]
    last = H[-1]
    if all(x == 0 for x in last[:-1]) and last[-1] != 0:
        return int(abs(last[-1])), max(abs(int(x)) for r in H for x in r)
    return None, max(abs(int(x)) for r in H for x in r)

def gfp_ok(rows, nc, p):
    M = flint.nmod_mat(len(rows) + 1, nc, p)
    for i, row in enumerate(rows):
        for c, a in row.items(): M[i, c] = a % p
    r1 = M.rank()
    M[len(rows), nc - 1] = 1
    return M.rank() == r1

def degs(n, polys, dmax, primes=(2, 3, 5, 1000003)):
    dp = max(max(bin(m).count('1') for m in p) for p in polys)
    out = {str(p): None for p in primes}; out.update(Q=None, Z=None, D={}, maxentry={})
    for d in range(dp, dmax + 1):
        rows, nc = build(n, polys, d)
        r2, nc2 = unit_reduce(rows, nc)   # unimodular; pivots are units mod every p too
        for p in primes:
            if out[str(p)] is None and r2 and gfp_ok(r2, nc2, p): out[str(p)] = d
        D, mx = hnf_D(r2, nc2) if r2 else (None, 0)
        out.setdefault('resid', {})[d] = (len(r2), nc2)
        out['D'][d] = D; out['maxentry'][d] = mx
        if D is not None and out['Q'] is None: out['Q'] = d
        if D == 1: out['Z'] = d; break
    return out

def unit_reduce(rows, nc):
    """Unimodular sparse elimination on unit (+-1) pivots in non-constant columns.
    Pivot rows are discarded: since the target e has zeros on pivot columns and every
    later row is zero there, e in L  <=>  e in span(residual rows). Exact over Z."""
    rows = [dict(r) for r in rows if r]
    colrows = {}
    for i, r in enumerate(rows):
        for c in r: colrows.setdefault(c, set()).add(i)
    alive = set(range(len(rows)))
    for c in range(nc - 1):
        cand = [i for i in colrows.get(c, ()) if i in alive and c in rows[i]]
        units = [i for i in cand if abs(rows[i][c]) == 1]
        if not units: continue
        piv = min(units, key=lambda i: len(rows[i])); P = rows[piv]; s = P[c]
        for i in cand:
            if i == piv: continue
            R = rows[i]; f = R[c] * s  # R - f*P kills column c
            for cc, a in P.items():
                v = R.get(cc, 0) - f * a
                if v: 
                    if cc not in R: colrows.setdefault(cc, set()).add(i)
                    R[cc] = v
                else:
                    R.pop(cc, None)
            if not R: alive.discard(i)
        alive.discard(piv)
    res = [rows[i] for i in alive if rows[i]]
    cols = sorted({c for r in res for c in r} | {nc - 1})
    remap = {c: j for j, c in enumerate(cols)}
    return [{remap[c]: a for c, a in r.items()} for r in res], len(cols)

def hnf_D_fast(rows, nc):
    r2, nc2 = unit_reduce(rows, nc)
    if not r2: return None, 0
    return hnf_D(r2, nc2)

if __name__ == '__main__':
    n = int(sys.argv[1]); dens = float(sys.argv[2]); K = int(sys.argv[3]); dmax = int(sys.argv[4]) if len(sys.argv) > 4 else n
    m = round(dens * n); res = []; seed = 0
    while len(res) < K:
        seed += 1; s = 1000 * n + seed; cl = rand3cnf(n, m, s)
        if not unsat(n, cl): continue
        t = time.time(); r = degs(n, [clause_poly(c) for c in cl], dmax)
        r['seed'] = s; r['t'] = round(time.time() - t, 1); res.append(r)
        print(n, dens, json.dumps(r), flush=True)
    json.dump(res, open(f'/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/n3/res_n{n}_d{dens}.json', 'w'))

