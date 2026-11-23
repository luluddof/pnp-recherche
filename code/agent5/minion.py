"""Check the polymorphism-minion conditions characterising each relaxation (Boolean templates):
 BLP     : symmetric polymorphism of arity k            (f = g(weight))
 AIP     : alternating polymorphism of arity 2L+1       (f = g(sum_odd - sum_even))
 BLP+AIP : 2-block-symmetric polymorphism, blocks L+1,L (f = g(w1, w2))
Existence of such g is itself a CSP(Gamma) instance on the 'free structure'; solved by SAT."""
import itertools
from pysat.solvers import Glucose3
from relax import TEMPLATES

def multisets(R, k):
    return itertools.combinations_with_replacement(sorted(R), k)

def vectors(R, kind, k):
    r = len(next(iter(R))); out = set()
    if kind == 'BLP':
        for ms in multisets(R, k):
            out.add(tuple(sum(t[i] for t in ms) for i in range(r)))
    else:
        L = k
        for A in multisets(R, L + 1):
            for B in multisets(R, L):
                a = tuple(sum(t[i] for t in A) for i in range(r))
                b = tuple(sum(t[i] for t in B) for i in range(r))
                out.add(tuple(x - y for x, y in zip(a, b)) if kind == 'AIP' else tuple(zip(a, b)))
    return out

def exists(T, kind, k):
    var = {}
    def v(key):
        if key not in var: var[key] = len(var) + 1
        return var[key]
    cls = []
    for R in T:
        r = len(next(iter(R)))
        for vec in vectors(R, kind, k):
            for t in itertools.product((0, 1), repeat=r):
                if t not in R:  # forbid (g(vec_1..vec_r)) == t
                    cls.append([(-v(vec[i]) if t[i] else v(vec[i])) for i in range(r)])
    s = Glucose3(bootstrap_with=cls); ok = s.solve(); s.delete(); return ok

if __name__ == '__main__':
    for name, T in TEMPLATES.items():
        blp = [k for k in range(1, 8) if exists(T, 'BLP', k)]
        aip = [L for L in range(0, 5) if exists(T, 'AIP', L)]
        ba = [L for L in range(1, 5) if exists(T, 'BA', L)]
        print(f'{name:7s} BLP-sym arities ok {blp} | AIP alt L ok {aip} (arity 2L+1) | BLP+AIP block L ok {ba}')
