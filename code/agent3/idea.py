"""Test the key step of the 'p-evaluation' attack on AC0[p]-Frege.

Plan of attack (Razborov-Smolensky meets PC pseudo-expectations):
  1. A hard formula F has no F_p degree-D PC refutation  <=>  there is a linear functional
     L : P_D -> F_p with L(1)=1 and L(V_D)=0  (V_D = degree-D PC closure of F's axioms).
  2. Replace each AC0[p] gate g in a Frege proof by its Smolensky approximator A_g (degree ~ polylog).
  3. WANT: L(error_g) 'small' for every gate, so a union bound shows the approximated proof is
     still a valid-enough low-degree derivation, contradicting step 1.
We compute (a) PC degree over F_3 vs F_2 for Tseitin (step 1 is real), and (b) the set of values
L(error) can take over all valid L -- this tests whether 'small' in step 3 has any meaning.
Then (c) the real-field horn: the approximate degree of MOD_2/MOD_3 over R (needed for positivity).
"""
import itertools, random, sys
import numpy as np
import networkx as nx
from scipy.optimize import linprog

P = 3

def mul_var(poly, x):
    out = {}
    for m, c in poly.items():
        mm = m | {x}
        out[mm] = (out.get(mm, 0) + c) % P
    return {m: c for m, c in out.items() if c}

def poly_add(a, b, s=1):
    out = dict(a)
    for m, c in b.items():
        out[m] = (out.get(m, 0) + s * c) % P
    return {m: c for m, c in out.items() if c}

def indicator_wrong_parity(vs, charge):
    """Multilinear 0/1 polynomial that is 1 on assignments violating sum(vs)=charge mod 2."""
    poly = {}
    for bits in itertools.product([0, 1], repeat=len(vs)):
        if sum(bits) % 2 == charge:
            continue
        # prod x_i^{b_i} (1-x_i)^{1-b_i}
        term = {frozenset(): 1}
        for v, b in zip(vs, bits):
            if b:
                term = mul_var(term, v)
            else:
                term = poly_add(term, mul_var(term, v), -1)
        poly = poly_add(poly, term)
    return poly

class Echelon:
    """Fully-reduced echelon basis over F_P with a graded monomial order (leading = highest degree)."""
    def __init__(self):
        self.rows = {}  # lead monomial -> poly (lead coef 1)
    @staticmethod
    def key(m):
        return (len(m), tuple(sorted(m)))
    def lead(self, poly):
        return max(poly, key=self.key)
    def reduce(self, poly):
        poly = dict(poly)
        changed = True
        while poly:
            # reduce the highest monomial that is a pivot
            cand = [m for m in poly if m in self.rows]
            if not cand:
                break
            m = max(cand, key=self.key)
            poly = poly_add(poly, self.rows[m], -poly[m])
        return poly
    def add(self, poly):
        r = self.reduce(poly)
        if not r:
            return None
        l = self.lead(r)
        inv = pow(r[l], P - 2, P)
        r = {m: (c * inv) % P for m, c in r.items()}
        self.rows[l] = r
        return r

def pc_closure(axioms, nvars, D):
    E = Echelon()
    queue = []
    for a in axioms:
        r = E.add(a)
        if r is not None:
            queue.append(r)
    while queue:
        f = queue.pop()
        # we multiply every element ever added whose degree < D; since the span of such elements
        # equals V cap P_{D-1} only up to cancellation, we also re-add multiples of reduced forms.
        if max(len(m) for m in f) >= D:
            continue
        for x in range(1, nvars + 1):
            g = mul_var(f, x)
            r = E.add(g)
            if r is not None:
                queue.append(r)
    # second pass: guarantee closure w.r.t. V cap P_{D-1} using fully reduced basis
    changed = True
    while changed:
        changed = False
        for l, f in list(E.rows.items()):
            if len(l) < D:
                for x in range(1, nvars + 1):
                    if E.add(mul_var(f, x)) is not None:
                        changed = True
    return E

def tseitin_axioms(nv, seed):
    G = nx.random_regular_graph(3, nv, seed=seed)
    eid = {frozenset(e): i + 1 for i, e in enumerate(G.edges())}
    ax = []
    for u in G:
        vs = sorted(eid[frozenset((u, w))] for w in G[u])
        ax.append(indicator_wrong_parity(vs, 1 if u == 0 else 0))
    return ax, len(eid)

def main():
    global P
    out = []
    print("(a) PC degree of Tseitin (3-regular, random) over F_2 and F_3")
    for nv in [6, 8, 10]:
        for p in [2, 3]:
            P = p
            ax, ne = tseitin_axioms(nv, 1)
            for D in range(1, ne + 1):
                E = pc_closure(ax, ne, D)
                if frozenset() in E.rows:  # constant 1 derivable
                    print(f"  |V|={nv} edges={ne} F_{p}: PC degree = {D}", flush=True)
                    break
    print("(b) Pseudo-expectation values on a Smolensky error polynomial (F_3)")
    P = 3
    for nv in [8, 10]:
        ax, ne = tseitin_axioms(nv, 1)
        # largest D with no refutation
        D = 1
        while frozenset() not in pc_closure(ax, ne, D + 1).rows:
            D += 1
        E = pc_closure(ax, ne, D)
        nmon = sum(1 for k in range(D + 1) for _ in itertools.combinations(range(ne), k))
        print(f"  |V|={nv}: no F_3 PC refutation in degree {D}; dim P_D={nmon}, dim V_D={len(E.rows)},"
              f" pseudo-expectations form an affine space of dim {nmon - len(E.rows) - 1}")
        rng = random.Random(0)
        vals_seen = set()
        for trial in range(20):
            k = 3
            xs = rng.sample(range(1, ne + 1), k)
            S = [x for x in xs if rng.random() < 0.5] or xs[:1]
            # OR(x) = 1 - prod(1-x_i)
            OR = {frozenset(): 1}
            prod = {frozenset(): 1}
            for x in xs:
                prod = poly_add(prod, mul_var(prod, x), -1)
            OR = poly_add(OR, prod, -1)
            lin = {frozenset([x]): 1 for x in S}
            A = {}
            for m1, c1 in lin.items():
                for m2, c2 in lin.items():
                    A[m1 | m2] = (A.get(m1 | m2, 0) + c1 * c2) % P
            A = {m: c for m, c in A.items() if c}
            err = poly_add(OR, A, -1)
            if not err or max(len(m) for m in err) > D:
                continue
            r = E.reduce(err)
            # L(err) is forced iff r is a constant (multiple of 1); otherwise L(err) ranges over all of F_3
            forced = (not r) or set(r) == {frozenset()}
            vals_seen.add('forced=%s' % (r.get(frozenset(), 0) if forced else None) if forced else 'free')
        print("  values of L(error) over valid pseudo-expectations:", vals_seen)
    print("(c) Real-field horn: best uniform error of degree-d real polynomial for MOD_2 / MOD_3 on m bits")
    for modq in [2, 3]:
        m = 24
        row = []
        for d in [2, 4, 8, 12, 16, 20, 23]:
            # symmetrized: q univariate of degree d on k=0..m ; target f(k)=1 iff k%modq==0
            ks = np.arange(m + 1)
            t = (ks % modq == 0).astype(float)
            V = np.vstack([np.polynomial.chebyshev.chebvander(2 * ks / m - 1, d)]).astype(float)
            nc = d + 1
            # variables: coeffs (nc), eps ; minimize eps s.t. -eps <= Vc - t <= eps
            c = np.zeros(nc + 1); c[-1] = 1
            A = np.vstack([np.hstack([V, -np.ones((m + 1, 1))]), np.hstack([-V, -np.ones((m + 1, 1))])])
            b = np.concatenate([t, -t])
            res = linprog(c, A_ub=A, b_ub=b, bounds=[(None, None)] * nc + [(0, None)])
            row.append((d, round(res.x[-1], 4)))
        print(f"  MOD_{modq}, m={m}:", row)

if __name__ == '__main__':
    main()
