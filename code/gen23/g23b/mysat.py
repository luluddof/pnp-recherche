# Independent SAT encoding of STU arrays with AT MOST sS S-colours, sT T-colours, sU U-colours (colours may be unused).
# Usage: python3 mysat.py N sS sT sU [cap]   (cap = optional max class size)
import sys, itertools, time
from pysat.solvers import Cadical153
from pysat.card import CardEnc, EncType
from core import predicate, identity_tensor_ok, ranks

def solve(N, sS, sT, sU, cap=None, timeout=None):
    pal = ['S'] * sS + ['T'] * sT + ['U'] * sU; n = len(pal)
    cnt = itertools.count(1)
    x = [[[next(cnt) for k in range(n)] for b in range(N)] for a in range(N)]
    sg = [[next(cnt) for b in range(N)] for a in range(N)]   # true = sign -1
    cl = []
    for a in range(N):
        for b in range(N):
            cl.append(x[a][b][:])
            for k1, k2 in itertools.combinations(range(n), 2): cl.append([-x[a][b][k1], -x[a][b][k2]])
    for k in range(n):
        for a in range(N):
            for b1, b2 in itertools.combinations(range(N), 2):
                cl.append([-x[a][b1][k], -x[a][b2][k]]); cl.append([-x[b1][a][k], -x[b2][a][k]])
    mono = {}
    def M(c1, c2, X):
        key = (min(c1, c2), max(c1, c2), X)
        if key not in mono:
            v = next(cnt); mono[key] = v
            (a1, b1), (a2, b2) = key[0], key[1]
            ks = [k for k in range(n) if pal[k] == X]
            for k in range(n):
                if pal[k] == X:
                    cl.append([-x[a1][b1][k], -x[a2][b2][k], v])       # both k -> mono_X
                    cl.append([-v, -x[a1][b1][k], x[a2][b2][k]])       # mono_X & c1=k -> c2=k
                else:
                    cl.append([-v, -x[a1][b1][k]])                     # mono_X -> c1 not outside X
        return mono[key]
    partner = {'S': 'T', 'T': 'S', 'U': 'U'}
    for r1, r2 in itertools.combinations(range(N), 2):
        for c1, c2 in itertools.combinations(range(N), 2):
            D = ((r1, c1), (r2, c2)); A = ((r1, c2), (r2, c1))
            for P, Q in ((D, A), (A, D)):
                for X in 'STU':
                    if X == 'S' and sS == 0 or X == 'T' and sT == 0 or X == 'U' and sU == 0: continue
                    if partner[X] == 'S' and sS == 0 or partner[X] == 'T' and sT == 0:
                        cl.append([-M(*P, X)]); continue
                    m = M(*P, X); cl.append([-m, M(*Q, partner[X])])
                    # sign: odd number of minus signs among the 4 cells
                    s4 = [sg[r1][c1], sg[r2][c2], sg[r1][c2], sg[r2][c1]]
                    for bits in itertools.product([0, 1], repeat=4):
                        if sum(bits) % 2 == 0:
                            cl.append([-m] + [(-v if bt else v) for v, bt in zip(s4, bits)])
    # symmetry: first row and column signs +
    for a in range(N): cl.append([-sg[a][0]]); cl.append([-sg[0][a]])
    # colour first-occurrence order inside each palette (reading order)
    cells = [(a, b) for a in range(N) for b in range(N)]
    for X in 'STU':
        ks = [k for k in range(n) if pal[k] == X]
        if len(ks) < 2: continue
        prev = {}
        for i, (a, b) in enumerate(cells):
            cur = {}
            for k in ks:
                v = next(cnt); cur[k] = v      # v <-> colour k occurs at a cell <= i
                if i == 0: cl += [[-v, x[a][b][k]], [v, -x[a][b][k]]]
                else: cl += [[-v, x[a][b][k], prev[k]], [v, -x[a][b][k]], [v, -prev[k]]]
            for j in range(1, len(ks)):
                if i == 0: cl.append([-x[a][b][ks[j]]])
                else: cl.append([-x[a][b][ks[j]], prev[ks[j-1]]])
            prev = cur
    if cap is not None:
        top = next(cnt)
        for k in range(n):
            e = CardEnc.atmost([x[a][b][k] for a in range(N) for b in range(N)], bound=cap, top_id=top, encoding=EncType.seqcounter)
            top = max(top, e.nv) + 1; cl += e.clauses
    S = Cadical153(bootstrap_with=cl)
    r = S.solve()
    if not r: return None
    m = set(l for l in S.get_model() if l > 0)
    col = [[next(k for k in range(n) if x[a][b][k] in m) for b in range(N)] for a in range(N)]
    sgn = [[-1 if sg[a][b] in m else 1 for b in range(N)] for a in range(N)]
    return col, dict(enumerate(pal)), sgn

if __name__ == "__main__":
    N, sS, sT, sU = map(int, sys.argv[1:5]); cap = int(sys.argv[5]) if len(sys.argv) > 5 else None
    t = time.time(); r = solve(N, sS, sT, sU, cap)
    msg = f"N={N} S<={sS} T<={sT} U<={sU} cap={cap}: " + ("SAT" if r else "UNSAT") + f" {time.time()-t:.1f}s"
    if r:
        col, pal, sgn = r
        used = set(c for row in col for c in row)
        nS = sum(1 for k in used if pal[k] == 'S'); nT = sum(1 for k in used if pal[k] == 'T'); nU = sum(1 for k in used if pal[k] == 'U')
        ok, rr, rc = ranks(N, col, pal, sgn)
        msg += f" used S,T,U={nS},{nT},{nU} pred={predicate(N, col, pal, sgn)} R+C=X {ok} rk(R)={rr} rk(C)={rc}"
    print(msg, flush=True)
