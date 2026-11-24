"""G18S: small DF instances.  Instance = (u, wires), wires = [(row party, col party)], parties in {0..2^u-1}.
End bit 2e = row end of wire e (at party p), 2e+1 = column end (at party q).  Direction d: group = bit d of the party."""
import itertools

def bounds(u, wires):
    """returns dict: thmA, A' (vertex cut), A'' (oriented wire subset) -- all as log2 N lower bounds."""
    M = len(wires); P = 2 ** u
    def cs(sel):  # sel: list of (party of row-index end, party of column-index end)
        out = []
        for d in range(u):
            c = 0
            for g in (0, 1):
                e1 = sum(1 for (a, b) in sel if (a >> d) & 1 == g); e2 = sum(1 for (a, b) in sel if (b >> d) & 1 == g)
                c += min(e1, e2)
            out.append(c)
        return out
    bestA = bestAp = 0; argAp = None
    for asg in itertools.product((0, 1, 2), repeat=P):
        sel = []
        for (p, q) in wires:
            if asg[p] == 1 and asg[q] == 2: sel.append((p, q))
            elif asg[p] == 2 and asg[q] == 1: sel.append((q, p))
        if not sel: continue
        c = cs(sel); w = len(sel)
        if w - max(c) > bestAp: bestAp = w - max(c); argAp = asg
        # Thm A hypothesis: no group of any direction meets both C1 and C2
        okA = all(not (asg[x] == 1 and asg[y] == 2 and (x >> d) & 1 == (y >> d) & 1) for x in range(P) for y in range(P) for d in range(u))
        if okA: bestA = max(bestA, w)
    bestApp = 0; argApp = None
    for o in itertools.product((0, 1, 2), repeat=M):
        sel = [(p, q) if oo == 1 else (q, p) for (p, q), oo in zip(wires, o) if oo]
        if not sel: continue
        v = len(sel) - max(cs(sel))
        if v > bestApp: bestApp = v; argApp = o
    return dict(thmA=bestA, Ap=bestAp, App=bestApp, argAp=argAp, argApp=argApp)

def groups(u, wires, d):
    party = [x for w in wires for x in w]
    g0 = [k for k in range(2 * len(wires)) if (party[k] >> d) & 1 == 0]
    g1 = [k for k in range(2 * len(wires)) if (party[k] >> d) & 1 == 1]
    return g0, g1

def sat_search(u, wires, dirs, timeout=None, sym=True):
    """GF(2) search: I == sum_k tau0^(k)(bits of group 0 of dirs[k]) * tau1^(k)(bits of group 1)  (mod 2).
    returns ('SAT', tables) / ('UNSAT', None) / ('TIMEOUT', None)."""
    from pysat.solvers import Cadical153
    from threading import Timer
    M = len(wires); nb = 2 * M; N = len(dirs)
    cnt = itertools.count(1); var = {}
    def V(*key):
        if key not in var: var[key] = next(cnt)
        return var[key]
    s = Cadical153()
    for b in itertools.product((0, 1), repeat=nb):
        tgt = int(all(b[2 * e] == b[2 * e + 1] for e in range(M)))
        ands = []
        for k, d in enumerate(dirs):
            g0, g1 = groups(u, wires, d)
            a = V('t', k, 0, tuple(b[i] for i in g0)); c = V('t', k, 1, tuple(b[i] for i in g1)); z = V('z', k, b)
            s.add_clause([-z, a]); s.add_clause([-z, c]); s.add_clause([z, -a, -c]); ands.append(z)
        # xor(ands) == tgt
        for signs in itertools.product((0, 1), repeat=N):
            if sum(signs) % 2 != tgt:   # forbid assignments with wrong parity
                s.add_clause([(-z if sg else z) for z, sg in zip(ands, signs)])
    res = s.solve()
    if res is None: return 'TIMEOUT', None
    if not res: return 'UNSAT', None
    model = set(x for x in s.get_model() if x > 0)
    tabs = {key: int(v in model) for key, v in var.items() if key[0] == 't'}
    return 'SAT', tabs

def verify_gf2(u, wires, dirs, tabs):
    M = len(wires); nb = 2 * M
    for b in itertools.product((0, 1), repeat=nb):
        tot = 0
        for k, d in enumerate(dirs):
            g0, g1 = groups(u, wires, d)
            tot ^= tabs[('t', k, 0, tuple(b[i] for i in g0))] & tabs[('t', k, 1, tuple(b[i] for i in g1))]
        if tot != int(all(b[2 * e] == b[2 * e + 1] for e in range(M))): return False
    return True
