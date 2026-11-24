# G24d: run the st.py encoding with an external binary (kissat / cadical) on DIMACS; optional DRAT proof.
# usage: python3 stk.py N p q s timeout [kissat|cadical] [proof]
import sys, os, subprocess, time, json
import st
N, p, q, s, tl = map(int, sys.argv[1:6]); binname = sys.argv[6] if len(sys.argv) > 6 else 'kissat'
proof = len(sys.argv) > 7
H = os.path.dirname(os.path.abspath(__file__))
B = {'kissat': H + '/ext/kissat/build/kissat', 'cadical': H + '/ext/cadical/build/cadical'}[binname]
bal = (p == q)
E, x, inter = st.build(N, p, q, s, bal)
base = f"k_{N}_{p}_{q}_s{s}_{binname}"
extra = []
it = 0
while True:
    it += 1
    with open(base + '.cnf', 'w') as f:
        f.write(f"p cnf {E.nv} {len(E.cl) + len(extra)}\n")
        for c in E.cl + extra: f.write(' '.join(map(str, c)) + ' 0\n')
    cmd = [B, base + '.cnf'] + ([base + '.drat'] if proof else [])
    if binname == 'kissat' and proof: cmd += ['--no-binary']
    t0 = time.time()
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=tl).stdout
    except subprocess.TimeoutExpired:
        print(f"N={N} (p,q)=({p},{q}) s={s} [{binname}]: TIMEOUT {tl}s", flush=True); sys.exit(0)
    if '\ns UNSATISFIABLE' in '\n' + out:
        print(f"N={N} (p,q)=({p},{q}) s={s} [{binname}]: UNSAT it={it} [{time.time()-t0:.1f}s]" + (f" proof {base}.drat" if proof else ''), flush=True)
        sys.exit(0)
    lits = set()
    for line in out.splitlines():
        if line.startswith('v '): lits.update(int(t) for t in line[2:].split() if int(t) > 0)
    n = p + q
    col = [[next(k for k in range(n) if x[a][b][k] in lits) for b in range(N)] for a in range(N)]
    sg, conf = st.gf2_signs(N, col)
    if sg is None:
        extra.append([-inter[t] for t in conf]); print("  sign conflict blocked", flush=True); continue
    v = st.verify(N, p, q, col, sg)
    print(f"N={N} (p,q)=({p},{q}) s={s} [{binname}]: SAT it={it} [{time.time()-t0:.1f}s] " + json.dumps(v), flush=True)
    json.dump(dict(N=N, p=p, q=q, s=s, col=col, sg=sg, pal='k<q:S else T', verify=v), open(f"sol_{N}_{p}_{q}_s{s}.json", 'w'))
    os.remove(base + '.cnf')
    break
