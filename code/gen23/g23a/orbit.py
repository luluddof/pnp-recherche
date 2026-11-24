# Orbit (group-invariant) SAT encoding of ST intercalate arrays.
# Rows = cols = group Gam (order N), left-multiplication invariant: cell (a, a c) has colour (i(c), a t_c H_i)
# (block i(c) with palette S/T and normal subgroup H_i acting as colour stabiliser), t_c in Gam.
# Q(c,d,v)  <=>  cells (1,c) and (v, v d) have the same colour.  Theorem-1 conditions become:
#   Q(c,d,1) false (c != d); Q(c,d,v) false if v d = c; Q(c,d,v) -> Q(vd, v^-1 c, v) and pal(c) != pal(vd).
# Signs: solved afterwards on the full array (F_2 linear system), blocking clause on failure.
import sys, itertools, time, json, os
sys.path.insert(0, '/home/user/pnp-recherche/code/gen22/g22b'); sys.path.insert(0, '/home/user/pnp-recherche/code/gen22/g22c')
from pysat.solvers import Cadical195
from pysat.card import CardEnc, EncType

def gstar(t):
    g = 1
    while (g + 1) * g // 2 <= t: g += 1
    return g

def normal_subgroups(mul, inv, order):
    N = len(mul); out = []
    for h in range(N):   # cyclic subgroups <h> of given order + products if needed (enough for |H| <= 4 here)
        S = {0}; x = h
        while x not in S: S.add(x); x = mul[x][h]
        if len(S) == order and all(mul[mul[g][s]][inv[g]] in S for g in range(N) for s in S):
            fs = frozenset(S)
            if fs not in out: out.append(fs)
    if order == 4:
        for a in range(N):
            for b in range(N):
                S = {0, a, b, mul[a][b]}
                if len(S) == 4 and all(mul[x][y] in S for x in S for y in S) and all(mul[mul[g][s]][inv[g]] in S for g in range(N) for s in S):
                    fs = frozenset(S)
                    if fs not in out: out.append(fs)
    return out

class Enc:
    def __init__(self, mul, inv, blocks, caps=None, exact=None, fixt0=True, zpair=None, autos=()):
        """blocks: list of (pal, H) with H a frozenset (normal subgroup). caps[i]: max #differences in block i."""
        self.mul, self.inv, self.blocks = mul, inv, blocks
        N = self.N = len(mul); nb = len(blocks); self.top = 0
        def nv(): self.top += 1; return self.top
        T = self.T = [[nv() for x in range(N)] for c in range(N)]
        B = self.B = [[nv() for i in range(nb)] for c in range(N)]
        P = self.P = [nv() for c in range(N)]
        cl = self.cl = []
        Hs = []; cls_of = []
        for (pal, H) in blocks:
            if H not in Hs: Hs.append(H)
            cls_of.append(Hs.index(H))
        self.Hs = Hs
        canon = [[min(mul[x][h] for h in H) == x for x in range(N)] for H in Hs]
        for c in range(N):
            cl.append(T[c][:]); cl.append(B[c][:])
            for x, y in itertools.combinations(range(N), 2): cl.append([-T[c][x], -T[c][y]])
            for i, j in itertools.combinations(range(nb), 2): cl.append([-B[c][i], -B[c][j]])
            for i, (pal, H) in enumerate(blocks):
                cl.append([-B[c][i], P[c] if pal == 'S' else -P[c]])
                for x in range(N):
                    if not canon[cls_of[i]][x]: cl.append([-B[c][i], -T[c][x]])
        if fixt0: cl.append([T[0][0]])
        # SBH_k(c,d)
        SBH = self.SBH = {}
        for k, H in enumerate(Hs):
            ik = [i for i in range(nb) if cls_of[i] == k]
            for c in range(N):
                for d in range(c, N):
                    s = nv(); SBH[(k, c, d)] = SBH[(k, d, c)] = s
                    if c == d:
                        cl.append([-s] + [B[c][i] for i in ik])
                        for i in ik: cl.append([-B[c][i], s])
                    else:
                        for i in ik:
                            cl.append([-B[c][i], -B[d][i], s]); cl.append([-s, -B[c][i], B[d][i]]); cl.append([-s, -B[d][i], B[c][i]])
                        cl.append([-s] + [B[c][i] for i in ik])
        # Q variables, symmetric Q(c,d,v) = Q(d,c,v^-1)
        Q = self.Q = {}
        for c in range(N):
            for d in range(N):
                for v in range(N):
                    key = (c, d, v); alt = (d, c, inv[v])
                    if alt in Q: Q[key] = Q[alt]
                    else: Q[key] = nv()
        for c in range(N):
            for d in range(N):
                for v in range(N):
                    q = Q[(c, d, v)]
                    if c == d and v == 0: continue
                    if (v == 0) or mul[v][d] == c: cl.append([-q]); continue
                    if c <= d or True:
                        cl.append([-q] + [SBH[(k, c, d)] for k in range(len(Hs))])
                        e, ep = mul[v][d], mul[inv[v]][c]
                        cl.append([-q, Q[(e, ep, v)]])
                        cl.append([-q, P[c], P[e]]); cl.append([-q, -P[c], -P[e]])
                        for k, H in enumerate(Hs):
                            s = SBH[(k, c, d)]
                            for x in range(N):
                                if not canon[k][x]: continue
                                vx = mul[inv[v]][x]
                                cl.append([-q, -s, -T[c][x]] + [T[d][mul[vx][h]] for h in H])
        # forward: same colour -> Q
        for k, H in enumerate(Hs):
            for c in range(N):
                for d in range(N):
                    s = SBH[(k, c, d)]
                    if c == d:
                        for h in H:
                            if h: cl.append([-s, Q[(c, c, h)]])
                        continue
                    for x in range(N):
                        if not canon[k][x]: continue
                        for y in range(N):
                            if not canon[k][y]: continue
                            for h in H:
                                v = mul[mul[x][h]][inv[y]]
                                cl.append([-s, -T[c][x], -T[d][y], Q[(c, d, v)]])
        # block caps, symmetry (identical block types ordered by min element)
        for i in range(nb):
            lits = [B[c][i] for c in range(N)]
            if exact and exact[i] is not None:
                enc = CardEnc.equals(lits=lits, bound=exact[i], top_id=self.top, encoding=EncType.seqcounter)
            elif caps and caps[i] is not None:
                enc = CardEnc.atmost(lits=lits, bound=caps[i], top_id=self.top, encoding=EncType.seqcounter)
            else: continue
            self.top = max(self.top, enc.nv); cl += enc.clauses
        for phi in autos:   # equivariance: t(phi c) = phi(t(c)), block(phi c) = block(c)
            for c in range(N):
                for x in range(N): cl.append([-T[c][x], T[phi[c]][phi[x]]])
                for i in range(nb): cl.append([-B[c][i], B[phi[c]][i]])
        if zpair is not None:   # z central involution; blocks paired 2j <-> 2j+1 (S <-> T); t(zc) = t(c), block(zc) = partner
            z = zpair
            for c in range(N):
                zc = mul[z][c]
                for x in range(N): cl.append([-T[c][x], T[zc][x]])
                for i in range(nb): cl.append([-B[c][i], B[zc][i ^ 1]])
            for i in range(0, nb - 2, 2):
                if blocks[i] == blocks[i + 2]:
                    for c in range(N): cl.append([-B[c][i + 2]] + [B[cc][i] for cc in range(c)])
        for i in range(nb - 1):
            if zpair is not None: break
            if blocks[i] == blocks[i + 1]:
                for c in range(N): cl.append([-B[c][i + 1]] + [B[cc][i] for cc in range(c)])
                cl.append([B[c][i] for c in range(N)])   # nonempty
    def decode(self, model):
        m = set(l for l in model if l > 0); N = self.N; mul = self.mul
        t = [next(x for x in range(N) if self.T[c][x] in m) for c in range(N)]
        bl = [next(i for i in range(len(self.blocks)) if self.B[c][i] in m) for c in range(N)]
        return t, bl
    def array(self, t, bl):
        N = self.N; mul = self.mul; col = [[None] * N for _ in range(N)]
        for a in range(N):
            for c in range(N):
                pal, H = self.blocks[bl[c]]
                g = min(mul[mul[a][t[c]]][h] for h in H)
                col[a][mul[a][c]] = (pal, bl[c] * N + g)
        return col
    def block_clause(self, t, bl):
        return [-self.T[c][t[c]] for c in range(self.N)] + [-self.B[c][bl[c]] for c in range(self.N)]

def sign_solve(col):
    """F_2 system: each intercalate (both diagonals mono) has sign product -1. Returns signs or None."""
    N = len(col); cls = {}
    for a in range(N):
        for b in range(N): cls.setdefault(col[a][b], []).append((a, b))
    rows = []
    for k, lst in cls.items():
        for (a, b), (a2, b2) in itertools.combinations(lst, 2):
            v = 0
            for (x, y) in ((a, b), (a2, b2), (a, b2), (a2, b)): v ^= 1 << (x * N + y)
            rows.append((v, 1))
    piv = {}
    for v, r in rows:
        while v:
            h = v.bit_length() - 1
            if h in piv: pv, pr = piv[h]; v ^= pv; r ^= pr
            else: piv[h] = (v, r); break
        if v == 0 and r == 1: return None
    x = 0
    for h in sorted(piv):
        pv, pr = piv[h]; val = pr ^ (bin(pv & x & ~(1 << h)).count('1') & 1)
        if val: x |= 1 << h
    return [[-1 if (x >> (a * N + b)) & 1 else 1 for b in range(N)] for a in range(N)]

def check_colouring(col):
    """direct check of Theorem 1 (i),(ii) colour part; returns #violations."""
    N = len(col); bad = 0; cls = {}
    for a in range(N):
        for b in range(N): cls.setdefault(col[a][b], []).append((a, b))
    for k, lst in cls.items():
        if len(set(a for a, b in lst)) < len(lst) or len(set(b for a, b in lst)) < len(lst): bad += 1
        for (a, b), (a2, b2) in itertools.combinations(lst, 2):
            k1, k2 = col[a][b2], col[a2][b]
            if k1 != k2 or k1[0] == k[0]: bad += 1
    return bad

def ev_check(col, sg, dense=False):
    import numpy as np, scipy.sparse as sp, ev
    N = len(col); cls = {}
    for a in range(N):
        for b in range(N): cls.setdefault(col[a][b], []).append((a, b, sg[a][b]))
    terms = []
    for k, lst in cls.items():
        if k[0] != 'S': continue
        G = np.zeros((N, N), dtype=np.int64)
        for a, b, s in lst: G[a, b] = s
        terms.append((G, G.T.copy()))
    p, mth = ev.p_of(N, terms)
    out = {'q': len(terms), 'p_ev': p, 'method': mth}
    if dense: out['p_dense'] = ev.dense_check(N, terms)
    return out

def run(name, G, blocks, caps=None, exact=None, timeout=500, tries=50, verify=True, outdir='.', zpair=None, autos=()):
    mul, inv = G[0], G[1]; N = len(mul)
    t0 = time.time(); E = Enc(mul, inv, blocks, caps, exact, zpair=zpair, autos=autos)
    t1 = time.time()
    from pysat.solvers import Solver
    S = Solver(name=os.environ.get('SOLVER', 'cadical195'), bootstrap_with=E.cl)
    nT = sum(N // len(H) for pal, H in blocks if pal == 'T'); nS = sum(N // len(H) for pal, H in blocks if pal == 'S')
    tag = f"{name} N={N} blocks={[ (p, len(H)) for p, H in blocks]} caps={caps} exact={exact} -> target (p,q)<=({nT},{nS})"
    print(tag, f"clauses={len(E.cl)} build={t1-t0:.1f}s", flush=True)
    res = None
    for k in range(tries):
        r = S.solve()
        if not r: print("  UNSAT" if k == 0 else f"  no more colourings after {k} sign failures", f"{time.time()-t1:.1f}s", flush=True); res = 'UNSAT' if k == 0 else 'SIGNFAIL'; break
        t, bl = E.decode(S.get_model()); col = E.array(t, bl)
        bad = check_colouring(col); sg = sign_solve(col)
        if sg is None:
            S.add_clause(E.block_clause(t, bl)); continue
        cells = [[(col[a][b][0], col[a][b][1], sg[a][b]) for b in range(N)] for a in range(N)]
        nSc = len(set(col[a][b] for a in range(N) for b in range(N) if col[a][b][0] == 'S'))
        nTc = len(set(col[a][b] for a in range(N) for b in range(N) if col[a][b][0] == 'T'))
        print(f"  SAT after {k} sign-retries, {time.time()-t1:.1f}s: colour violations {bad}; colours T={nTc} S={nSc}; t={t} blocks={bl}", flush=True)
        info = {'name': name, 'N': N, 'blocks': [(p, sorted(H)) for p, H in blocks], 't': t, 'bl': bl, 'cells': cells, 'p': nTc, 'q': nSc}
        if verify:
            if N <= 32:
                from shape import verify_shape
                ok, p, q = verify_shape(N, cells, quiet=True); info['verify_shape'] = (ok, p, q)
                print(f"  verify_shape: identity {ok}, rr(R)={p}, rc(C)={q}", flush=True)
            e = ev_check(col, sg, dense=(N <= 16)); info['ev'] = e; print("  ev.py:", e, flush=True)
        fn = os.path.join(outdir, f"sol_{name}_{N}_{nTc}_{nSc}.json".replace(' ', '').replace(':', '-').replace('^', ''))
        json.dump(info, open(fn, 'w')); res = 'SAT'; break
    S.delete()
    return res
