"""G14a: GENERAL (non-commutative) class over GF(2^e): SAT search for J(t) with a given rank profile
r_k = max(rank B01_k, rank B10_k) (G11c (1a),(1c)), followed by an independent exact verification.
usage: python3 gen_char2.py l d e r_1,..,r_l [tau indices as field elements, comma sep] [timeout s]"""
import sys, itertools, time, threading
from pysat.solvers import Cadical153
POLY = {1: 0b11, 2: 0b111, 3: 0b1011, 4: 0b10011}

class GF:
    def __init__(s, e):
        s.e, s.q, s.poly = e, 1 << e, POLY[e]
        s.inv = {a: next(b for b in range(1, s.q) if s.mul(a, b) == 1) for a in range(1, s.q)}
    def mul(s, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> s.e & 1: a ^= s.poly
        return r
    def rank(s, M):
        M = [row[:] for row in M]; r = 0
        for c in range(len(M[0]) if M else 0):
            piv = next((i for i in range(r, len(M)) if M[i][c]), None)
            if piv is None: continue
            M[r], M[piv] = M[piv], M[r]; iv = s.inv[M[r][c]]
            M[r] = [s.mul(iv, v) for v in M[r]]
            for i in range(len(M)):
                if i != r and M[i][c]:
                    f = M[i][c]; M[i] = [a ^ s.mul(f, b) for a, b in zip(M[i], M[r])]
            r += 1
        return r
    def polyeval(s, co, t):
        r = 0
        for c in reversed(co): r = s.mul(r, t) ^ c
        return r
    def interp(s, pts):  # Lagrange -> coefficient list
        n = len(pts); co = [0] * n
        for j, (xj, yj) in enumerate(pts):
            num = [1]; den = 1
            for i, (xi, _) in enumerate(pts):
                if i == j: continue
                num = [0] + num                      # multiply by (t + xi)
                for a in range(len(num) - 1): num[a] ^= s.mul(xi, num[a + 1])
                den = s.mul(den, xj ^ xi)
            f = s.mul(yj, s.inv[den])
            for a in range(len(num)): co[a] ^= s.mul(f, num[a])
        return co

def search(l, d, e, prof, tau, tmo=300):
    F = GF(e); n = 1 << l; full = n - 1
    S = Cadical153(); cnt = [0]
    def new(k=1):
        v = list(range(cnt[0] + 1, cnt[0] + k + 1)); cnt[0] += k; return v
    def mulgate(u, v):
        w = new(e)
        for a in range(F.q):
            for b in range(F.q):
                c = F.mul(a, b)
                pre = [(-u[i] if a >> i & 1 else u[i]) for i in range(e)] + [(-v[i] if b >> i & 1 else v[i]) for i in range(e)]
                for i in range(e): S.add_clause(pre + [w[i] if c >> i & 1 else -w[i]])
        return w
    def xor(lits, rhs):
        lits = list(lits)
        while len(lits) > 4:
            a = new()[0]; chunk = lits[:3]; lits = [a] + lits[3:]
            for sg in itertools.product([1, -1], repeat=4):      # a = xor(chunk)
                if (sum(1 for x in sg if x < 0) % 2) == 1: S.add_clause([sg[0] * a] + [sg[i + 1] * chunk[i] for i in range(3)])
        m = len(lits)
        if m == 0:
            if rhs: S.add_clause([])
            return
        for sg in itertools.product([1, -1], repeat=m):
            # forbid assignments with wrong parity: assignment x_i = (sg_i<0) ; clause = OR(sg_i * lit_i) is violated exactly there
            if (sum(1 for x in sg if x < 0) % 2) != rhs: S.add_clause([sg[i] * lits[i] for i in range(m)])
    # entries: ent[(k,x,y)] = list of products (each e bit-vars) for x_k != y_k
    ent = {}
    for k in range(l):
        for side in (0, 1):
            rows = [x for x in range(n) if (x >> k & 1) == side]; cols = [y for y in range(n) if (y >> k & 1) != side]
            terms = []
            for i in range(prof[k]):
                U = {x: new(e) for x in rows}; V = {y: new(e) for y in cols}; terms.append((U, V))
            for x in rows:
                for y in cols:
                    ent[(k, x, y)] = [mulgate(U[x], V[y]) for U, V in terms]
    def constmul_bits(c):  # matrix: out bit b = xor of in bits i with (c*2^i)>>b &1
        return [[(F.mul(c, 1 << i) >> b) & 1 for i in range(e)] for b in range(e)]
    for x in range(n):
        for y in range(n):
            D = [k for k in range(l) if (x ^ y) >> k & 1]
            if len(D) <= d: continue
            J0 = 1 if y == x ^ full else 0
            base = D[:d]
            for kk in D[d:]:
                nodes = [0] + [tau[k] for k in base] + [tau[kk]]; ks = [None] + base + [kk]
                lam = []
                for j in range(len(nodes)):
                    den = 1
                    for i in range(len(nodes)):
                        if i != j: den = F.mul(den, nodes[j] ^ nodes[i])
                    lam.append(F.inv[den])
                const = F.mul(lam[0], J0)
                for b in range(e):
                    lits = []
                    for j in range(1, len(nodes)):
                        Mb = constmul_bits(lam[j])[b]
                        for prod in ent[(ks[j], x, y)]:
                            lits += [prod[i] for i in range(e) if Mb[i]]
                    xor(lits, (const >> b) & 1)
    timer = threading.Timer(tmo, lambda: S.interrupt()); timer.start()
    res = S.solve_limited(expect_interrupt=True); timer.cancel()
    if res is None: return "TIMEOUT", None
    if not res: return "UNSAT", None
    model = set(v for v in S.get_model() if v > 0)
    val = lambda bits: sum(1 << i for i in range(e) if bits[i] in model)
    E = {key: 0 for key in ent}
    for key, prods in ent.items():
        for p_ in prods: E[key] ^= val(p_)
    # reconstruct J(t) entrywise and verify independently
    J = [[None] * n for _ in range(n)]
    for x in range(n):
        for y in range(n):
            D = [k for k in range(l) if (x ^ y) >> k & 1]; J0 = 1 if y == x ^ full else 0
            pts = [(0, J0)] + [(tau[k], E[(k, x, y)]) for k in D[:d]]
            co = F.interp(pts); co += [0] * (d + 1 - len(co)); assert len(co) <= d + 1
            for k in D: assert F.polyeval(co, tau[k]) == E[(k, x, y)], "interpolation check failed"
            assert F.polyeval(co, 0) == J0
            J[x][y] = co
    ranks = []
    for k in range(l):
        Jk = [[F.polyeval(J[x][y], tau[k]) for y in range(n)] for x in range(n)]
        r01 = F.rank([[Jk[x][y] for y in range(n) if y >> k & 1] for x in range(n) if not x >> k & 1])
        r10 = F.rank([[Jk[x][y] for y in range(n) if not y >> k & 1] for x in range(n) if x >> k & 1])
        ranks.append((r01, r10))
    return "SAT", (J, ranks)

if __name__ == "__main__":
    l, d, e = map(int, sys.argv[1:4]); prof = list(map(int, sys.argv[4].split(',')))
    tau = list(map(int, sys.argv[5].split(','))) if len(sys.argv) > 5 else list(range(1, l + 1))
    tmo = float(sys.argv[6]) if len(sys.argv) > 6 else 300
    assert len(set(tau)) == l and 0 not in tau and max(tau) < (1 << e)
    t0 = time.time(); st, out = search(l, d, e, prof, tau, tmo)
    print(f"(l,d)=({l},{d}) GF(2^{e}) tau={tau} profile={prof} sum={sum(prof)}: {st} {time.time()-t0:.1f}s", flush=True)
    if st == "SAT":
        J, ranks = out
        print("   verified (r01,r10) per k:", ranks, " sum max =", sum(max(a) for a in ranks))
        if '--dump' in sys.argv:
            import json; json.dump(J, open(f"sol_{l}_{d}_gf{1<<e}_{'-'.join(map(str,prof))}.json", "w"))
