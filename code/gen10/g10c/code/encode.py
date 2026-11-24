# G10c: independent implementation of the 3-CNF Psi(Phi,p,w) of Lemma 1 for the three skeletons
# (Steiner trees, l-cycles, l-vertex paths), following the text of consequences.md Section 3.
class CNF:
    def __init__(self): self.nv = 0; self.cl = []; self.gates = 0
    def new(self): self.nv += 1; return self.nv
    def add(self, *lits): assert 1 <= len(lits) <= 3; self.cl.append(tuple(lits))
    def XOR(self, a, b): z = self.new(); self.gates += 1; self.add(-z, a, b); self.add(-z, -a, -b); self.add(z, -a, b); self.add(z, a, -b); return z
    def AND(self, a, b): z = self.new(); self.gates += 1; self.add(-z, a); self.add(-z, b); self.add(z, -a, -b); return z
    def OR(self, a, b):  z = self.new(); self.gates += 1; self.add(z, -a); self.add(z, -b); self.add(-z, a, b); return z
    def counter(self, bits):
        """outputs list of literals o_0..o_k (binary value of the number of ones among bits); O(len(bits)) gates."""
        buckets = [list(bits)]; out = []; i = 0
        while i < len(buckets):
            bk = buckets[i]
            while len(bk) >= 2:
                if len(buckets) == i + 1: buckets.append([])
                if len(bk) >= 3:
                    a, b, c = bk.pop(), bk.pop(), bk.pop()
                    s1 = self.XOR(a, b); s = self.XOR(s1, c); c1 = self.AND(a, b); c2 = self.AND(s1, c); co = self.OR(c1, c2)
                else:
                    a, b = bk.pop(), bk.pop(); s = self.XOR(a, b); co = self.AND(a, b)
                bk.insert(0, s); buckets[i + 1].append(co)
            out.append(bk[0] if bk else None); i += 1
        return out
    def equals_const(self, out, val):
        for i, o in enumerate(out):
            bit = val >> i & 1
            if o is None:
                if bit: self.add_false()
            else: self.add(o if bit else -o)
        if val >> len(out): self.add_false()
    def add_false(self):
        z = self.new(); self.add(z); self.add(-z)
    def leq_const(self, bits, q):
        """returns literal r with r <-> value(bits) <= q  (bits little-endian), O(len(bits)) gates; None means constant True"""
        if q >= (1 << len(bits)) - 1: return None
        r = None  # None = True
        for i, c in enumerate(bits):
            if q >> i & 1:  # r_i = not c_i or r_{i-1}
                r = None if r is None else self.OR(-c, r)
            else:           # r_i = not c_i and r_{i-1}
                r = -c if r is None else self.AND(-c, r)
        return r

def build(n, E, p, w, kind, l, K=()):
    """w[e] in 0..p ; w[e]==0 means edge deleted (unit clause not x_e). returns CNF, primary var info"""
    F = CNF(); m = len(E)
    x = [F.new() for _ in E]; y = [F.new() for _ in range(n)]
    for j, (u, v) in enumerate(E): F.add(-x[j], y[u]); F.add(-x[j], y[v])          # (i)
    if kind == 'ST':
        for v in K: F.add(y[v])
        F.equals_const(F.counter(y), l + 1); F.equals_const(F.counter(x), l)
    elif kind in ('CYC', 'PATH'):
        incl = [[] for _ in range(n)]
        for j, (a_, b_) in enumerate(E): incl[a_].append(x[j]); incl[b_].append(x[j])
        for v in range(n):
            inc = incl[v]
            o = F.counter(inc) if inc else []
            if kind == 'CYC':   # y_v -> exactly two incident selected edges
                if len(o) < 2: F.add(-y[v])
                else:
                    for i, oi in enumerate(o):
                        if oi is None:
                            if i == 1: F.add(-y[v])
                        else: F.add(-y[v], oi if i == 1 else -oi)
            else:               # degree <= 2
                if len(o) >= 2 and o[0] is not None and o[1] is not None: F.add(-o[0], -o[1])
                for oi in o[2:]:
                    if oi is not None: F.add(-oi)
        if kind == 'CYC': F.equals_const(F.counter(y), l)
        else: F.equals_const(F.counter(y), l); F.equals_const(F.counter(x), l - 1)
    M0 = len(F.cl)
    for j in range(m):
        if w[j] == 0: F.add(-x[j])
    b = max(1, (p - 1).bit_length())            # ceil(log2 p)
    assert (1 << b) >= p
    c = [[F.new() for _ in range(b)] for _ in range(n)]
    for v in range(n):
        r = F.leq_const(c[v], p - 1)                                  # (a)
        if r is not None: F.add(r)
        for i in range(b): F.add(y[v], -c[v][i])                      # (b)
    for j, (u, v) in enumerate(E):
        for i in range(b):                                            # (c)
            F.add(-x[j], -c[u][i], c[v][i]); F.add(-x[j], c[u][i], -c[v][i])
    g = [[F.new() for _ in range(b)] for _ in E]
    for j in range(m):
        for i in range(b): F.add(x[j], -g[j][i])                      # (d)
        if w[j] >= 1:
            r = F.leq_const(g[j], w[j] - 1)                           # (e)
            if r is not None: F.add(-x[j], r)
    return F, M0, b

def dimacs(F):
    return "p cnf %d %d\n" % (F.nv, len(F.cl)) + "".join(" ".join(map(str, c)) + " 0\n" for c in F.cl)
