# Semiregular generalisation: Gamma (order M) acts freely on rows Gamma x [k] and cols Gamma x [k] (N = M k), by left mult.
# "difference" delta = (c, al, be): cell ((a, al), (a c, be)).  Colour of that cell = (i(delta), a t_delta)  (free colour orbits).
# Q(d1, d2, v) <=> cells ((1,al1),(c1,be1)) and ((v,al2),(v c2,be2)) have the same colour.
import sys, itertools, time, json, os
from pysat.solvers import Cadical195
from pysat.card import CardEnc, EncType
from orbit import gstar, sign_solve, check_colouring, ev_check

class Enc2:
    def __init__(self, mul, inv, k, blocks, caps=None, sizes=None):
        self.mul, self.inv, self.k, self.blocks = mul, inv, k, blocks
        M = self.M = len(mul); N = self.N = M * k; nb = len(blocks); self.top = 0
        D = self.D = [(c, al, be) for c in range(M) for al in range(k) for be in range(k)]
        di = self.di = {d: n for n, d in enumerate(D)}; nD = len(D)
        def nv(): self.top += 1; return self.top
        T = self.T = [[nv() for x in range(M)] for d in range(nD)]
        B = self.B = [[nv() for i in range(nb)] for d in range(nD)]
        P = self.P = [nv() for d in range(nD)]
        cl = self.cl = []
        for d in range(nD):
            cl.append(T[d][:]); cl.append(B[d][:])
            for x, y in itertools.combinations(range(M), 2): cl.append([-T[d][x], -T[d][y]])
            for i, j in itertools.combinations(range(nb), 2): cl.append([-B[d][i], -B[d][j]])
            for i, pal in enumerate(blocks): cl.append([-B[d][i], P[d] if pal == 'S' else -P[d]])
        cl.append([T[0][0]])
        SB = {}
        for d1 in range(nD):
            for d2 in range(d1 + 1, nD):
                s = nv(); SB[(d1, d2)] = SB[(d2, d1)] = s
                for i in range(nb):
                    cl.append([-B[d1][i], -B[d2][i], s]); cl.append([-s, -B[d1][i], B[d2][i]]); cl.append([-s, -B[d2][i], B[d1][i]])
        Q = {}
        for d1 in range(nD):
            for d2 in range(nD):
                for v in range(M):
                    if d1 == d2 and v == 0: continue
                    alt = (d2, d1, inv[v])
                    Q[(d1, d2, v)] = Q[alt] if alt in Q else nv()
        self.Q = Q
        for d1 in range(nD):
            c1, a1, b1 = D[d1]
            for d2 in range(nD):
                c2, a2, b2 = D[d2]
                for v in range(M):
                    if d1 == d2:
                        if v: cl.append([-Q[(d1, d2, v)]])     # free colours: (1,c),(v,vc) never same colour
                        continue
                    q = Q[(d1, d2, v)]
                    if (v == 0 and a1 == a2) or (mul[v][c2] == c1 and b1 == b2):
                        cl.append([-q])
                        for x in range(M): cl.append([-SB[(d1, d2)], -T[d1][x], -T[d2][mul[inv[v]][x]]])
                        continue
                    cl.append([-q, SB[(d1, d2)]])
                    e = di[(mul[v][c2], a1, b2)]; ep = di[(mul[inv[v]][c1], a2, b1)]
                    cl.append([-q, Q[(e, ep, v)]] if e != ep or v else [-q])
                    cl.append([-q, P[d1], P[e]]); cl.append([-q, -P[d1], -P[e]])
                    for x in range(M):
                        cl.append([-q, -T[d1][x], T[d2][mul[inv[v]][x]]])
                        # forward: same block, t1 = x, t2 = y  ->  Q(d1,d2, x y^-1)
                    for x in range(M):
                        y = mul[inv[v]][x]
                        cl.append([-SB[(d1, d2)], -T[d1][x], -T[d2][y], q])
        for i in range(nb):
            lits = [B[d][i] for d in range(nD)]
            if sizes and sizes[i] is not None: enc = CardEnc.equals(lits=lits, bound=sizes[i], top_id=self.top, encoding=EncType.seqcounter)
            elif caps and caps[i] is not None: enc = CardEnc.atmost(lits=lits, bound=caps[i], top_id=self.top, encoding=EncType.seqcounter)
            else: continue
            self.top = max(self.top, enc.nv); cl += enc.clauses
        for i in range(nb - 1):
            if blocks[i] == blocks[i + 1]:
                for d in range(nD): cl.append([-B[d][i + 1]] + [B[dd][i] for dd in range(d)])
    def decode(self, model):
        m = set(l for l in model if l > 0)
        t = [next(x for x in range(self.M) if self.T[d][x] in m) for d in range(len(self.D))]
        bl = [next(i for i in range(len(self.blocks)) if self.B[d][i] in m) for d in range(len(self.D))]
        return t, bl
    def array(self, t, bl):
        M, k, mul = self.M, self.k, self.mul; N = self.N
        col = [[None] * N for _ in range(N)]
        for n, (c, al, be) in enumerate(self.D):
            for a in range(M):
                col[a * k + al][mul[a][c] * k + be] = (self.blocks[bl[n]], bl[n] * M + mul[a][t[n]])
        return col
    def block_clause(self, t, bl):
        return [-self.T[d][t[d]] for d in range(len(self.D))] + [-self.B[d][bl[d]] for d in range(len(self.D))]

def run2(name, G, k, blocks, caps=None, sizes=None, tries=300, outdir='.'):
    mul, inv = G[0], G[1]; M = len(mul); N = M * k
    t0 = time.time(); E = Enc2(mul, inv, k, blocks, caps, sizes); t1 = time.time()
    nT = M * blocks.count('T'); nS = M * blocks.count('S')
    print(f"{name} M={M} k={k} N={N} blocks={''.join(blocks)} caps={caps} sizes={sizes} -> (p,q)<=({nT},{nS}) clauses={len(E.cl)} build={t1-t0:.1f}s", flush=True)
    S = Cadical195(bootstrap_with=E.cl); res = None
    for it in range(tries):
        if not S.solve():
            print("  UNSAT" if it == 0 else f"  no colouring with consistent signs ({it} tried)", f"{time.time()-t1:.1f}s", flush=True)
            res = 'UNSAT' if it == 0 else 'SIGNFAIL'; break
        t, bl = E.decode(S.get_model()); col = E.array(t, bl); sg = sign_solve(col)
        if sg is None: S.add_clause(E.block_clause(t, bl)); continue
        bad = check_colouring(col)
        nSc = len(set(col[a][b] for a in range(N) for b in range(N) if col[a][b][0] == 'S'))
        nTc = len(set(col[a][b] for a in range(N) for b in range(N) if col[a][b][0] == 'T'))
        print(f"  SAT ({it} sign-retries) {time.time()-t1:.1f}s: violations {bad}, T={nTc} S={nSc}", flush=True)
        cells = [[(col[a][b][0], col[a][b][1], sg[a][b]) for b in range(N)] for a in range(N)]
        info = {'name': name, 'M': M, 'k': k, 'N': N, 'blocks': blocks, 't': t, 'bl': bl, 'cells': cells, 'p': nTc, 'q': nSc}
        if N <= 32:
            sys.path.insert(0, '/home/user/pnp-recherche/code/gen22/g22b'); from shape import verify_shape
            ok, p, q = verify_shape(N, cells, quiet=True); info['verify_shape'] = (ok, p, q); print(f"  verify_shape: {ok} rr(R)={p} rc(C)={q}", flush=True)
        e = ev_check(col, sg, dense=(N <= 16)); info['ev'] = e; print("  ev.py:", e, flush=True)
        fn = os.path.join(outdir, f"sol2_{name}_k{k}_{N}_{nTc}_{nSc}.json".replace(':', '-').replace('^', ''))
        json.dump(info, open(fn, 'w')); res = 'SAT'; break
    S.delete(); return res

if __name__ == "__main__":
    from groups import catalogue
    import multiprocessing as mp
    M, k, blocks, to = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4])
    gl = sys.argv[5].split(',')
    nT = M * blocks.count('T'); nS = M * blocks.count('S')
    caps = [gstar(nT if b == 'S' else nS) for b in blocks]
    def child(g): run2(g, catalogue(M)[g], k, list(blocks), caps=caps)
    for g in gl:
        pr = mp.Process(target=child, args=(g,)); pr.start(); pr.join(to)
        if pr.is_alive(): pr.kill(); pr.join(); print(f"  TIMEOUT {to}s ({g})", flush=True)
