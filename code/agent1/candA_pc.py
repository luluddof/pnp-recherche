"""Candidate A: degree-d Polynomial Calculus over GF(2) (Gaussian elimination on a
truncated Macaulay matrix, closed under multiplication by variables).
Subsumes Gaussian elimination on XOR substructure (degree 3 suffices for 3-XOR /
Tseitin) and width-d resolution reasoning restricted to degree d.
Runtime: #monomials N = sum_{i<=d} C(n,i) = O(n^d). Basis has <= N rows; each
closure step multiplies a row by n variables and reduces (O(N) row ops, each O(N)
bits). Total O(n * N^3) = O(n^{3d+1}) bit ops: polynomial for fixed d.
Sound: returns 'UNSAT' only if 1 is in the derived ideal (a PC refutation).
"""
import itertools


class PC:
    def __init__(self, n, d):
        self.n, self.d = n, d
        self.idx = {}
        self.mons = []
        for k in range(d + 1):
            for S in itertools.combinations(range(1, n + 1), k):
                self.idx[frozenset(S)] = len(self.mons)
                self.mons.append(frozenset(S))
        self.basis = {}   # leading bit -> row (python int bitmask over monomials)
        self.lowdeg = []  # rows of degree < d to multiply by variables
        self.refuted = False

    def deg(self, row):
        return len(self.mons[row.bit_length() - 1]) if row else -1

    def poly_from(self, monsets):
        r = 0
        for S in monsets:
            S = frozenset(S)
            if len(S) > self.d:
                return None
            r ^= 1 << self.idx[S]
        return r

    def reduce(self, r):
        b = self.basis
        while r:
            lb = r.bit_length() - 1
            p = b.get(lb)
            if p is None:
                return r
            r ^= p
        return 0

    def add(self, r):
        r = self.reduce(r)
        if not r:
            return
        self.basis[r.bit_length() - 1] = r
        if r == 1:
            self.refuted = True
        self.lowdeg.append(r)

    def clause_poly(self, c):
        # product over literals of "literal false": x false -> (1+x); -x false -> x
        polys = [frozenset()]
        for l in c:
            v = abs(l)
            if l > 0:   # (1+x)
                new = {}
                for S in polys:
                    for T in (S, S | {v}):
                        new[T] = new.get(T, 0) ^ 1
                polys = [T for T, c_ in new.items() if c_]
            else:
                polys = list({S | {v} for S in polys})
        return polys

    def mul_var(self, r, v):
        out = 0
        x = r
        while x:
            lb = x.bit_length() - 1
            x ^= 1 << lb
            T = self.mons[lb] | {v}
            if len(T) > self.d:
                return None
            out ^= 1 << self.idx[T]
        return out

    def close(self, budget=None):
        i = 0
        while i < len(self.lowdeg) and not self.refuted:
            r = self.lowdeg[i]; i += 1
            if self.deg(r) >= self.d:
                continue
            for v in range(1, self.n + 1):
                q = self.mul_var(r, v)
                if q is not None:
                    self.add(q)
                    if self.refuted:
                        break
        self.lowdeg = []
        return self.refuted


def load(pc, cls):
    for c in cls:
        p = pc.clause_poly(c)
        r = pc.poly_from(p)
        if r is None:
            continue  # clause wider than d: skip (weaker but sound)
        pc.add(r)
        # also multiply clause by all monomials so total degree <= d (Macaulay rows)
        dc = max(len(S) for S in p)
        for k in range(1, pc.d - dc + 1):
            for M in itertools.combinations(range(1, pc.n + 1), k):
                q = pc.poly_from([S | set(M) for S in p])
                if q is not None:
                    pc.add(q)


def refute(n, cls, d):
    pc = PC(n, d)
    load(pc, cls)
    return pc.close()


def copy(pc):
    q = PC.__new__(PC)
    q.__dict__.update(pc.__dict__)
    q.basis = dict(pc.basis); q.lowdeg = list(pc.lowdeg)
    return q


def decide(n, cls, d):
    """Returns ('UNSAT', None) | ('SAT', assign) | ('FAIL', None).
    Greedy decimation using the degree-d PC oracle as lookahead."""
    base = PC(n, d)
    load(base, cls)
    if base.close():
        return 'UNSAT', None
    cur = base
    assign = {}
    forced_only = True
    for v in range(1, n + 1):
        ok = None
        for val in (0, 1):
            t = copy(cur)
            t.add(t.poly_from([{v}] if val == 0 else [set(), {v}]))  # x=0 : x ; x=1 : x+1
            if not t.close():
                ok = (val, t); break
        if ok is None:
            return ('UNSAT', None) if forced_only else ('FAIL', None)
        # if the other value was refutable, this choice was forced
        val, t = ok
        if val == 0:
            t2 = copy(cur); t2.add(t2.poly_from([set(), {v}]))
            if not t2.close():
                forced_only = False
        assign[v] = bool(val)
        cur = t
    return 'SAT', assign
