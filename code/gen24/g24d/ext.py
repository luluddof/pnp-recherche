# G24d: SAT-HUNTING ONLY (not used for UNSAT claims): ST arrays N x N with a FIXED corner sub-array (e.g. the 8x8 (8,8) array),
# no symmetry breaking except first-occurrence order of the non-corner colours.  usage: ext.py N p q corner.json timeout [solver]
import sys, itertools, time, json, os
from pysat.solvers import Solver
import st

def build(N, p, q, fixed, caps=True, extra=None):
    E = st.Enc(); n = p + q
    x = [[[E.new() for k in range(n)] for b in range(N)] for a in range(N)]
    P = [[E.new() for b in range(N)] for a in range(N)]
    for (a, b), k in fixed.items(): E.add([x[a][b][k]])
    for a in range(N):
        for b in range(N):
            E.add(x[a][b][:]); E.amo(x[a][b])
            for k in range(n): E.add([-x[a][b][k], P[a][b] if k < q else -P[a][b]])
    for k in range(n):
        for a in range(N):
            E.amo([x[a][b][k] for b in range(N)]); E.amo([x[b][a][k] for b in range(N)])
    eqv = {}
    def eq(c1, c2):
        key = (min(c1, c2), max(c1, c2))
        if key in eqv: return eqv[key]
        v = E.new(); eqv[key] = v
        (a1, b1), (a2, b2) = key
        for k in range(n):
            E.add([-v, -x[a1][b1][k], x[a2][b2][k]]); E.add([v, -x[a1][b1][k], -x[a2][b2][k]])
        return v
    inter = {}
    for r1, r2 in itertools.combinations(range(N), 2):
        for c1, c2 in itertools.combinations(range(N), 2):
            eD = eq((r1, c1), (r2, c2)); eA = eq((r1, c2), (r2, c1))
            E.add([-eD, eA]); E.add([eD, -eA])
            E.add([-eD, P[r1][c1], P[r1][c2]]); E.add([-eD, -P[r1][c1], -P[r1][c2]])
            inter[(r1, r2, c1, c2)] = eD
    if int(os.environ.get('ZSYM', '0')): usedS = list(range(q)); usedT = list(range(q, n))   # no first-occurrence with ZSYM
    else: usedS = None
    if usedS is None: usedS = sorted(set(k for k in fixed.values() if k < q)); usedT = sorted(set(k for k in fixed.values() if k >= q))
    scan = [(a, b) for a in range(N) for b in range(N) if (a, b) not in fixed]
    for ks in ([k for k in range(q) if k not in usedS], [k for k in range(q, n) if k not in usedT]):
        prev = None
        for i, (a, b) in enumerate(scan):
            cur = {}
            for k in ks:
                v = E.new(); cur[k] = v
                if prev is None: E.add([-v, x[a][b][k]]); E.add([v, -x[a][b][k]])
                else: E.add([-v, x[a][b][k], prev[k]]); E.add([v, -x[a][b][k]]); E.add([v, -prev[k]])
            for j in range(1, len(ks)):
                if prev is None: E.add([-x[a][b][ks[j]]])
                else: E.add([-x[a][b][ks[j]], prev[ks[j - 1]]])
            prev = cur
    m = int(os.environ.get('ZSYM', '0'))   # SAT-HUNTING ansatz: invariance under (a,b,k) -> (a+N/m, b+N/m, pi(k)), pi of order m
    if m:
        sh = N // m
        pi = lambda k: (k // m) * m + (k % m + 1) % m if k < q else q + ((k - q) // m) * m + ((k - q) % m + 1) % m
        for a in range(N):
            for b in range(N):
                for k in range(n):
                    E.add([-x[a][b][k], x[(a + sh) % N][(b + sh) % N][pi(k)]])
    if caps:
        cS, cT = st.gstar(N, p), st.gstar(N, q)
        for k in range(n): E.atmost([x[a][b][k] for a in range(N) for b in range(N)], cS if k < q else cT)
    return E, x, inter

if __name__ == "__main__":
    N, p, q = map(int, sys.argv[1:4]); cf = sys.argv[4]; tl = float(sys.argv[5]); solver = sys.argv[6] if len(sys.argv) > 6 else 'cadical195'
    c = json.load(open(cf)) if cf != 'none' else dict(N=0, q=0, col=[]); m = c['N']; cq = c['q']
    # recolour corner: its S colours -> 0.., T colours -> q..
    mapS, mapT = {}, {}
    fixed = {}
    for a in range(m):
        for b in range(m):
            k = c['col'][a][b]
            if k < cq: fixed[(a, b)] = mapS.setdefault(k, len(mapS))
            else: fixed[(a, b)] = q + mapT.setdefault(k, len(mapT))
    E, x, inter = build(N, p, q, fixed)
    print(f"  enc vars {E.nv} clauses {len(E.cl)}", flush=True)
    S = Solver(name=solver, bootstrap_with=E.cl)
    import threading
    timer = threading.Timer(tl, lambda: S.interrupt()); timer.start()
    t0 = time.time(); it = 0
    while True:
        it += 1
        r = S.solve_limited(expect_interrupt=True)
        if r is None: print(f"ext N={N} ({p},{q}) corner {cf}: TIMEOUT {tl}s", flush=True); break
        if not r: print(f"ext N={N} ({p},{q}) corner {cf}: UNSAT (with this corner only) [{time.time()-t0:.1f}s]", flush=True); break
        ms = set(l for l in S.get_model() if l > 0)
        col = [[next(k for k in range(p + q) if x[a][b][k] in ms) for b in range(N)] for a in range(N)]
        sg, conf = st.gf2_signs(N, col)
        if sg is None: S.add_clause([-inter[t] for t in conf]); print("  sign conflict blocked", flush=True); continue
        v = st.verify(N, p, q, col, sg)
        print(f"ext N={N} ({p},{q}): SAT it={it} [{time.time()-t0:.1f}s] " + json.dumps(v), flush=True)
        json.dump(dict(N=N, p=p, q=q, col=col, sg=sg, pal='k<q:S else T', verify=v, corner=cf), open(f"solext_{N}_{p}_{q}_z{os.environ.get('ZSYM', '0')}.json", 'w'))
        break
    timer.cancel()
