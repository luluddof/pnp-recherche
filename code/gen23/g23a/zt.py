# "Z-transversal" ansatz (generalises the Q8 point):  Gamma with a central subgroup Z (|Z| = m = number of blocks),
# Z0 < Z of index 2 (palette).  Blocks are labelled by Z; each Z-coset meets every block once:
#   element c = rep[cb] * w  (cb = coset, w in Z) lies in block  lam(cb)^-1 w ; palette S iff block in Z0.
# t is constant on Z-cosets: t(c) = tb(cb).  Cell (a, a c) gets colour (block(c), a t(c)).  (free colour action)
# Unknowns: tb (n cosets x N), lam (n x m).  Conditions (Theorem 1, reduced by Gamma- and Z-invariance), see notes.
import sys, itertools, time, json, os
from pysat.solvers import Cadical195
from orbit import sign_solve, check_colouring, ev_check

def central_subgroups(mul, inv, order):
    N = len(mul); cen = [z for z in range(N) if all(mul[z][x] == mul[x][z] for x in range(N))]
    out = []
    for S in itertools.combinations(cen[1:], order - 1):
        Z = frozenset((0,) + S)
        if all(mul[a][b] in Z for a in Z for b in Z) and Z not in out: out.append(Z)
    return out

class ZT:
    def __init__(self, mul, inv, Z, Z0, allowed=None, nopal=False):
        self.nopal = nopal
        self.mul, self.inv = mul, inv; N = self.N = len(mul)
        Zl = self.Zl = sorted(Z); m = len(Zl); self.Z0 = Z0
        # cosets
        cos = {}; reps = []
        for x in range(N):
            k = min(mul[x][z] for z in Zl)
            if k not in cos: cos[k] = len(reps); reps.append(k)
        self.reps = reps; n = self.n = len(reps)
        cb = self.cb = [cos[min(mul[x][z] for z in Zl)] for x in range(N)]
        wz = self.wz = [next(z for z in Zl if mul[reps[cb[x]]][z] == x) for x in range(N)]   # x = rep[cb(x)] * wz(x)
        self.top = 0
        def nv(): self.top += 1; return self.top
        tb = self.tb = [[nv() for x in range(N)] for c in range(n)]
        lam = self.lam = [[nv() for w in Zl] for c in range(n)]   # lam[c][j]: lam(c) = Zl[j]
        cl = self.cl = []
        for c in range(n):
            cl.append(tb[c][:]); cl.append(lam[c][:])
            for x, y in itertools.combinations(range(N), 2): cl.append([-tb[c][x], -tb[c][y]])
            for i, j in itertools.combinations(range(m), 2): cl.append([-lam[c][i], -lam[c][j]])
        c0 = cb[0]; cl.append([lam[c0][Zl.index(0)]])
        if allowed is None:
            allowed = [list(range(N)) for c in range(n)]; cl.append([tb[c0][0]])
        self.allowed = allowed
        for c in range(n):
            for x in range(N):
                if x not in allowed[c]: cl.append([-tb[c][x]])
        V = self.V = {}
        for c in range(n):
            for d in range(n):
                if c == d: continue
                for v in range(N): V[(c, d, v)] = nv()
                for x in allowed[c]:
                    for y in allowed[d]:
                        cl.append([-tb[c][x], -tb[d][y], V[(c, d, mul[x][inv[y]])]])
        zi = {z: j for j, z in enumerate(Zl)}
        for c in range(n):
            for d in range(n):
                if c == d: continue
                posv = set(mul[x][inv[y]] for x in allowed[c] for y in allowed[d])
                for v in range(N):
                    q = V[(c, d, v)]
                    if v not in posv: cl.append([-q]); continue
                    if v == 0: cl.append([-q]); continue
                    vd = mul[v][reps[d]]; e = cb[vd]; w1 = wz[vd]           # v rep[d] = rep[e] w1
                    vc = mul[inv[v]][reps[c]]; ep = cb[vc]; w2 = wz[vc]     # v^-1 rep[c] = rep[ep] w2
                    # t(e) = v t(ep)
                    if e == ep:
                        for x in range(N):
                            if mul[v][x] != x: cl.append([-q, -tb[e][x]])
                        # v x = x impossible for v != 1 -> forced false
                    else:
                        for x in range(N): cl.append([-q, -tb[ep][x], tb[e][mul[v][x]]])
                    # lambda conditions: over assignments of lam(c), lam(d), lam(e), lam(ep)
                    vars_ = sorted(set([c, d, e, ep]))
                    for ass in itertools.product(range(m), repeat=len(vars_)):
                        L = {u: Zl[a] for u, a in zip(vars_, ass)}
                        # c_el = rep[c] L[c] ; d_el = rep[d] L[d]  (both in block 1 = identity of Z, palette S)
                        # column distinct: v d_el != c_el
                        bad = mul[vd][L[d]] == mul[reps[c]][L[c]]
                        be = mul[mul[inv[L[e]]][w1]][L[d]]          # block(e_el)
                        bep = mul[mul[inv[L[ep]]][w2]][L[c]]        # block(ep_el)
                        bad = bad or be != bep or (be in Z0 if not getattr(self, 'nopal', False) else False)
                        if bad: cl.append([-q] + [-lam[u][a] for u, a in zip(vars_, ass)])
    def decode(self, model):
        s = set(l for l in model if l > 0)
        tb = [next(x for x in range(self.N) if self.tb[c][x] in s) for c in range(self.n)]
        lam = [self.Zl[next(j for j in range(len(self.Zl)) if self.lam[c][j] in s)] for c in range(self.n)]
        return tb, lam
    def array(self, tb, lam):
        N, mul, inv = self.N, self.mul, self.inv; col = [[None] * N for _ in range(N)]
        for c in range(N):
            blk = mul[inv[lam[self.cb[c]]]][self.wz[c]]
            pal = 'S' if blk in self.Z0 else 'T'
            for a in range(N): col[a][mul[a][c]] = (pal, blk * N + mul[a][tb[self.cb[c]]])
        return col
    def block_clause(self, tb, lam):
        return [-self.tb[c][tb[c]] for c in range(self.n)] + [-self.lam[c][self.Zl.index(lam[c])] for c in range(self.n)]

def runzt(name, G, Z, Z0, tries=200, outdir='.', verify=True, allowed=None):
    mul, inv = G[0], G[1]; N = len(mul); t0 = time.time()
    E = ZT(mul, inv, Z, Z0, allowed); t1 = time.time()
    m = len(Z); nS = N * len(Z0); nT = N * (m - len(Z0))
    print(f"{name} N={N} Z={sorted(Z)} Z0={sorted(Z0)} -> (p,q)=({nT},{nS}) clauses={len(E.cl)} build={t1-t0:.1f}s", flush=True)
    S = Cadical195(bootstrap_with=E.cl); res = None
    for it in range(tries):
        if not S.solve():
            print("  UNSAT" if it == 0 else f"  SIGNFAIL after {it} colourings", f"{time.time()-t1:.1f}s", flush=True); res = 'UNSAT' if it == 0 else 'SIGNFAIL'; break
        tb, lam = E.decode(S.get_model()); col = E.array(tb, lam); sg = sign_solve(col)
        if sg is None: S.add_clause(E.block_clause(tb, lam)); continue
        bad = check_colouring(col)
        nSc = len(set(col[a][b] for a in range(N) for b in range(N) if col[a][b][0] == 'S'))
        nTc = len(set(col[a][b] for a in range(N) for b in range(N) if col[a][b][0] == 'T'))
        print(f"  SAT ({it} sign-retries) {time.time()-t1:.1f}s: colour violations {bad}, T={nTc}, S={nSc}; tb={tb} lam={lam}", flush=True)
        cells = [[(col[a][b][0], col[a][b][1], sg[a][b]) for b in range(N)] for a in range(N)]
        info = {'name': name, 'N': N, 'Z': sorted(Z), 'Z0': sorted(Z0), 'tb': tb, 'lam': lam, 'reps': E.reps, 'cells': cells, 'p': nTc, 'q': nSc, 'violations': bad}
        if verify:
            if N <= 32:
                sys.path.insert(0, '/home/user/pnp-recherche/code/gen22/g22b'); from shape import verify_shape
                ok, p, q = verify_shape(N, cells, quiet=True); info['verify_shape'] = (ok, p, q); print(f"  verify_shape: identity {ok} rr(R)={p} rc(C)={q}", flush=True)
            e = ev_check(col, sg, dense=(N <= 16)); info['ev'] = e; print("  ev.py:", e, flush=True)
        json.dump(info, open(os.path.join(outdir, f"zt_{name}_{N}_{nTc}_{nSc}.json".replace(' ', '').replace(':', '-').replace('^', '')), 'w'))
        res = 'SAT'; break
    S.delete(); return res

if __name__ == "__main__":
    import multiprocessing as mp
    from groups import catalogue
    N, m, to = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]); gl = sys.argv[4].split(',')
    def child(g, Z, Z0):
        G = catalogue(N)[g] if N != 64 else __import__('g64').G64[g]()
        runzt(g, G, Z, Z0)
    for g in gl:
        G = catalogue(N)[g] if N != 64 else __import__('g64').G64[g]()
        mul, inv = G[0], G[1]
        for Z in central_subgroups(mul, inv, m)[:4]:
            for Z0 in central_subgroups(mul, inv, m // 2):
                if not Z0 <= Z: continue
                pr = mp.Process(target=child, args=(g, Z, Z0)); pr.start(); pr.join(to)
                if pr.is_alive(): pr.kill(); pr.join(); print(f"  TIMEOUT {to}s ({g} Z={sorted(Z)} Z0={sorted(Z0)})", flush=True)
