# cube-and-conquer for MC(N,p,q) branch r: one GF(2) SAT instance per orbit representative of C_1.
# usage: python3 cube.py N p q r [workers] [tl_per_cube] ; env SOLVER=kissat|cadical, ORD=1, DRAT=1, START,STOP (cube idx)
import sys, os, json, time, subprocess, itertools
from multiprocessing import Pool
import mc, orbits
EXT = '/tmp/claude-0/-home-user-pnp-recherche/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/gen24/g24d/ext'
BIN = {'kissat': EXT + '/kissat/build/kissat', 'cadical': EXT + '/cadical/build/cadical'}
N, p, q, r = map(int, sys.argv[1:5])
W = int(sys.argv[5]) if len(sys.argv) > 5 else 3
TL = int(sys.argv[6]) if len(sys.argv) > 6 else 300
SOLVER = os.environ.get('SOLVER', 'kissat'); ORD = bool(int(os.environ.get('ORD', '0'))); DRAT = bool(int(os.environ.get('DRAT', '0')))
D = f"cubes/N{N}_{p}_{q}_r{r}"; os.makedirs(D, exist_ok=True)

def work(arg):
    idx, C1 = arg
    Cm = [[0] * N for _ in range(N)]
    for y, v in C1: Cm[y][v] = 1
    E, e, c, U, V = mc.build(N, p, q, r, order=ORD, extra=Cm, c22=False)
    ev = [v for l in range(q) for X in (e[l], c[l]) for row in X for v in row]
    t0 = time.time(); nblock = 0
    if SOLVER == 'cms':
        import pycryptosat
        S = pycryptosat.Solver(threads=1, time_limit=TL)
        for x in E.cl: S.add_clause(x)
        for t, rr in E.xors: S.add_xor_clause(t, bool(rr))
        while True:
            sat, sol = S.solve()
            if sat is None or time.time() - t0 > TL: return idx, 'TIMEOUT', nblock, None, time.time() - t0
            if not sat: return idx, 'UNSAT', nblock, None, time.time() - t0
            Es = [[[int(sol[e[l][x][u]]) for u in range(N)] for x in range(N)] for l in range(q)]
            Cs = [[[int(sol[c[l][y][v]]) for v in range(N)] for y in range(N)] for l in range(q)]
            res = mc.sign_search(N, p, Es, Cs)
            if res is not None and res != 'MAXITER':
                E2, C2, rq = res
                return idx, 'SAT', nblock, {"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}, time.time() - t0
            if res == 'MAXITER': return idx, 'SIGN-MAXITER', nblock, None, time.time() - t0
            nblock += 1; S.add_clause([(-v if sol[v] else v) for v in ev])
    nv, cl = mc.xor_cnf(E)
    while True:
        fn = f"{D}/c{idx}.cnf"
        with open(fn, 'w') as f:
            f.write(f"p cnf {nv} {len(cl)}\n")
            for x in cl: f.write(' '.join(map(str, x)) + ' 0\n')
        rem = TL - (time.time() - t0)
        if rem <= 1: return idx, 'TIMEOUT', nblock, None, time.time() - t0
        cmd = [BIN[SOLVER], '-q', fn]
        if DRAT and SOLVER == 'kissat': cmd += [f"{D}/c{idx}.drat", '--no-binary']
        if DRAT and SOLVER == 'cadical': cmd += [f"{D}/c{idx}.drat"]
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=rem).stdout
        except subprocess.TimeoutExpired:
            return idx, 'TIMEOUT', nblock, None, time.time() - t0
        if 's UNSATISFIABLE' in out:
            if not DRAT: os.remove(fn)
            return idx, 'UNSAT', nblock, None, time.time() - t0
        if 's SATISFIABLE' not in out: return idx, 'ERR', nblock, out[-300:], time.time() - t0
        ms = set()
        for line in out.splitlines():
            if line.startswith('v'): ms.update(int(t) for t in line.split()[1:] if int(t) > 0)
        Es = [[[int(e[l][x][u] in ms) for u in range(N)] for x in range(N)] for l in range(q)]
        Cs = [[[int(c[l][y][v] in ms) for v in range(N)] for y in range(N)] for l in range(q)]
        res = mc.sign_search(N, p, Es, Cs)
        if res is not None and res != 'MAXITER':
            E2, C2, rq = res
            return idx, 'SAT', nblock, {"N": N, "p": p, "q": q, "rankQ": rq, "E": E2, "C": C2}, time.time() - t0
        if res == 'MAXITER': return idx, 'SIGN-MAXITER', nblock, None, time.time() - t0
        nblock += 1
        cl.append([(-v if v in ms else v) for v in ev])   # support provably unsignable (GF(3) rank > p)

if __name__ == "__main__":
    R = orbits.reps(N, r)
    a = int(os.environ.get('START', 0)); b = int(os.environ.get('STOP', len(R)))
    jobs = [(i, R[i][1]) for i in range(a, b)]
    if os.environ.get('LIST'): jobs = [(int(i), R[int(i)][1]) for i in os.environ['LIST'].split(',')]
    print(f"MC N={N} (p,q)=({p},{q}) r={r}: {len(R)} cubes, running {a}..{b-1}, {W} workers, TL {TL}s, {SOLVER}", flush=True)
    T0 = time.time(); stat = {}
    with Pool(W) as P:
        for idx, st, nb, dec, dt in P.imap_unordered(work, jobs):
            stat[st] = stat.get(st, 0) + 1
            print(f"  cube {idx} {R[idx][1]}: {st} (blocked {nb}) [{dt:.1f}s]", flush=True)
            if st == 'SAT':
                fn = f"mc_N{N}_{p}_{q}.json"; json.dump(dec, open(fn, 'w'))
                print(f"  SAT -> rankQ {dec['rankQ']} saved {fn}", flush=True)
                P.terminate(); break
    print(f"SUMMARY N={N} (p,q)=({p},{q}) r={r} cubes {a}..{b-1}: {stat} [{time.time()-T0:.1f}s]", flush=True)
