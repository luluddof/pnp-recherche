# Signed two-palette intercalate matrices ("the shape").
# A solution = N x N array: cell (a,b) -> (palette, colour, sign), palette in {'S','T','U'}:
#   S-colours l: G_l = sum sign E_ab  (column terms, weight 1),   T-colours j: H_j (row terms, weight 1),
#   U-colours k: weight 1/2 on both sides (2X = ... ; counted in p and in q).
# Identity:  X(x,y,u,v) = sum_l lam_l G_l[x,u]G_l[v,y] + sum_j mu_j H_j[x,y]H_j[v,u].
import numpy as np, itertools, flint, sys, time
from fractions import Fraction
def verify_shape(N, cells, quiet=False):
    """cells[a][b] = (pal, colour, sign). EXACT check of the identity by full expansion over [N]^4 (doubled, integers),
    and exact ranks (flint over Q) of R (rows (x,y), cols (u,v)) and of C (rows (x,u), cols (y,v))."""
    cls = {}
    for a in range(N):
        for b in range(N):
            pal, k, s = cells[a][b]; cls.setdefault((pal, k), []).append((a, b, s))
    R2 = np.zeros((N,)*4, dtype=np.int64); C2 = np.zeros((N,)*4, dtype=np.int64)   # 2R, 2C
    for (pal, k), lst in cls.items():
        wC = {'S': 2, 'T': 0, 'U': 1}[pal]; wR = 2 - wC
        for (a, b, s) in lst:
            for (c, d, t) in lst:
                if wC: C2[a, d, b, c] += wC * s * t      # G[x,u]G[v,y]: x=a,u=b,v=c,y=d
                if wR: R2[a, b, d, c] += wR * s * t      # H[x,y]H[v,u]: x=a,y=b,v=c,u=d
    X2 = np.zeros((N,)*4, dtype=np.int64)
    for x in range(N):
        for y in range(N): X2[x, y, y, x] = 2
    ok = bool((R2 + C2 == X2).all())
    Rm = flint.fmpz_mat(N*N, N*N, [int(e) for e in R2.reshape(-1)])
    Cm = flint.fmpz_mat(N*N, N*N, [int(e) for e in C2.transpose(0, 2, 1, 3).reshape(-1)])
    p, q = Rm.rank(), Cm.rank()
    nS = sum(1 for (pal, k) in cls if pal == 'S'); nT = sum(1 for (pal, k) in cls if pal == 'T'); nU = sum(1 for (pal, k) in cls if pal == 'U')
    if not quiet: print(f"  verify N={N}: identity R+C==X exact: {ok}; rr(R)={p}, rc(C)={q}; colours S,T,U = {nS},{nT},{nU} (predicted p={nT+nU}, q={nS+nU})")
    return ok, p, q
def show(N, cells):
    for a in range(N):
        print('   ' + ' '.join(f"{'+' if cells[a][b][2] > 0 else '-'}{cells[a][b][0]}{cells[a][b][1]:<2}" for b in range(N)))
def sat_shape(N, q, p, c=0, timeout=None, fix=None, extra=None, symbreak=True, seed=None, maxsize=None):
    """Find an N x N two-palette signed intercalate matrix with q S-colours, p T-colours, c U-colours. Returns cells or None."""
    from pysat.solvers import Cadical153, Glucose4
    from pysat.card import CardEnc
    n = q + p + c; pal = ['S']*q + ['T']*p + ['U']*c
    vid = itertools.count(1)
    col = [[[next(vid) for k in range(n)] for b in range(N)] for a in range(N)]
    sg = [[next(vid) for b in range(N)] for a in range(N)]
    cl = []
    for a in range(N):
        for b in range(N):
            cl.append([col[a][b][k] for k in range(n)])
            for k in range(n):
                for k2 in range(k+1, n): cl.append([-col[a][b][k], -col[a][b][k2]])
    for k in range(n):
        for a in range(N):
            for b in range(N):
                for b2 in range(b+1, N): cl.append([-col[a][b][k], -col[a][b2][k]]); cl.append([-col[b][a][k], -col[b2][a][k]])
    E = {}
    for a in range(N):
        for c2 in range(a+1, N):
            for b in range(N):
                for d in range(N):
                    if b != d: E[(a, b, c2, d)] = {X: next(vid) for X in 'STU'}
    for (a, b, c2, d), e in E.items():
        for k in range(n):
            cl.append([-col[a][b][k], -col[c2][d][k], e[pal[k]]])
            for X in 'STU':
                if pal[k] == X:
                    cl.append([-e[X], -col[a][b][k], col[c2][d][k]]); cl.append([-e[X], col[a][b][k], -col[c2][d][k]])
                else:
                    cl.append([-e[X], -col[a][b][k]])
        e2 = E[(a, d, c2, b)]
        cl.append([-e['S'], e2['T']]); cl.append([-e['T'], e2['S']]); cl.append([-e['U'], e2['U']])
        if b < d:
            s4 = [sg[a][b], sg[c2][d], sg[a][d], sg[c2][b]]
            for X in 'STU':
                for bits in itertools.product([0, 1], repeat=4):
                    if sum(bits) % 2 == 0:   # forbid even parity (xor must be 1)
                        cl.append([-e[X]] + [(-v if bt else v) for v, bt in zip(s4, bits)])
    if symbreak:
        for a in range(N): cl.append([sg[a][0]]); cl.append([sg[0][a]])
        # colour precedence inside each palette along the cells in reading order (weak form: colour k of a palette may be used
        # at cell i only if colour k-1 of the same palette is used at an earlier cell)
        cellsl = [(a, b) for a in range(N) for b in range(N)]
        for X, ks in (('S', range(0, q)), ('T', range(q, q+p)), ('U', range(q+p, n))):
            ks = list(ks)
            if len(ks) < 2: continue
            used = {}
            for i, (a, b) in enumerate(cellsl):
                for k in ks:
                    v = next(vid); used[(i, k)] = v
                    # v <-> col at cell i or used before
                    if i == 0: cl.append([-v, col[a][b][k]]); cl.append([v, -col[a][b][k]])
                    else:
                        pv = used[(i-1, k)]
                        cl.append([-v, col[a][b][k], pv]); cl.append([v, -col[a][b][k]]); cl.append([v, -pv])
                for j in range(1, len(ks)):
                    if i == 0: cl.append([-col[a][b][ks[j]]])
                    else: cl.append([-col[a][b][ks[j]], used[(i-1, ks[j-1])]])
    if maxsize is not None:
        top = next(vid) - 1
        for k in range(n):
            enc = CardEnc.atmost(lits=[col[a][b][k] for a in range(N) for b in range(N)], bound=maxsize, top_id=top)
            top = max(top, enc.nv); cl += enc.clauses
    if fix:
        for (a, b, k) in fix: cl.append([col[a][b][k]])
    if extra: cl += extra(col, sg)
    S = Cadical153(bootstrap_with=cl)
    t0 = time.time()
    r = S.solve()
    if not r: return None
    m = set(l for l in S.get_model() if l > 0)
    cells = [[None]*N for _ in range(N)]
    cnt = {}
    for a in range(N):
        for b in range(N):
            k = next(k for k in range(n) if col[a][b][k] in m)
            cells[a][b] = (pal[k], k, 1 if sg[a][b] in m else -1)
    return cells
if __name__ == "__main__":
    N, q, p, c = map(int, sys.argv[1:5])
    t0 = time.time(); cells = sat_shape(N, q, p, c)
    print(f"N={N} (p,q)=({p+c},{q+c}) [S={q},T={p},U={c}]:", "SAT" if cells else "UNSAT", f"{time.time()-t0:.1f}s", flush=True)
    if cells: show(N, cells); verify_shape(N, cells)
