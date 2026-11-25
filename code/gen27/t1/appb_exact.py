# Exact implementation of the pieces A_s of Appendix B (eqs. (47),(48),(53)) for small h, N = p(h) (or a prefix of
# the partitions of h, N < p(h)).  Checks sum_s A_s = h! X exactly, then computes EXACT ranks over Q (python-flint)
# of the threshold decompositions  R_{s0} = sum_{s<s0} A_s,  C_{s0} = h! X - R_{s0}:  rr(R_{s0}), rc(C_{s0}),
# and of each single piece (rr(A_s), rc(A_s)).  Compares with the counting bounds P(2h+s0-1), P(2h-s0).
import sys, math, functools, flint, json, time
from math import comb
def gbinom(k, j):  # generalized binomial k choose j, k any integer, j >= 0
    if j < 0: return 0
    num = 1
    for t in range(j): num *= (k - t)
    return num // math.factorial(j)
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k): yield (k,) + rest
def mono_mul_t(sig, l):  # t_l * t^sig
    return tuple(sorted(sig + (l,), reverse=True))
def nabla(sig, k):  # nabla_k t^sig = k mult_k(sig) t^sig / t_k ; returns (coef, sig') or None
    m = sig.count(k)
    if m == 0: return None
    lst = list(sig); lst.remove(k)
    return k*m, tuple(lst)
def addto(d, key, c):
    if c:
        v = d.get(key, 0) + c
        if v: d[key] = v
        else: d.pop(key, None)
@functools.lru_cache(maxsize=None)
def apply_mono(mu, r, sig):
    """ (t^mu)[r] applied to the monomial t^sig; returns tuple of (sig', coef).  Recursion (O2) on the first part. """
    E = sum(sig)
    if E + r < 0: return ()
    if len(mu) == 0: return ((sig, 1),) if r == 0 else ()
    n, rest = mu[0], mu[1:]
    out = {}
    # multiplications: sum_{l >= n} C(l-1, n-1) t_l a'[r-l]
    for l in range(n, E + r + 1):
        c = comb(l-1, n-1)
        for s2, c2 in apply_mono(rest, r - l, sig):
            addto(out, mono_mul_t(s2, l), c*c2)
    # differentiations: sum_{k >= 1} C(-k-1, n-1) a'[r+k] nabla_k
    for k in set(sig):
        c = gbinom(-k-1, n-1); cc, s1 = nabla(sig, k)
        for s2, c2 in apply_mono(rest, r + k, s1):
            addto(out, s2, c*cc*c2)
    return tuple(out.items())
def apply_poly(a, r, f):  # a, f dicts
    out = {}
    for mu, ca in a.items():
        for sig, cf in f.items():
            for s2, c2 in apply_mono(mu, r, sig): addto(out, s2, ca*cf*c2)
    return out
def poly_mul(f, g):
    out = {}
    for a, ca in f.items():
        for b, cb in g.items(): addto(out, tuple(sorted(a+b, reverse=True)), ca*cb)
    return out
def zmu(mu):
    z = 1
    for k in set(mu): m = mu.count(k); z *= k**m * math.factorial(m)
    return z
def pieces(h, N=None):
    parts = list(partitions(h))   # lexicographic, nonincreasing sequences, decreasing lex order; order irrelevant
    parts = sorted(parts)          # increasing lex order as in the paper (numbering from 0)
    if N is None: N = len(parts)
    mus = parts[:N]
    def y(n): return {(j,): comb(n, j) for j in range(1, n+1)}
    A = {}
    for mu in mus:
        a = {(): 1}
        for m in mu:
            a = poly_mul(a, {k: (-1)**(m-1)*c for k, c in y(m).items()})
        A[mu] = a
    def a_bracket(mu, i):  # a<i> = sum_h' C(h'+h, i) abar_h'
        return {k: comb(sum(k)+h, i)*c for k, c in A[mu].items() if comb(sum(k)+h, i)}
    B = {nu: {nu: math.factorial(h)//zmu(nu)} for nu in mus}
    I = list(range(-2*h, 2*h+1))
    # T[s][x][y][u][v]
    T = {s: [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)] for s in I}
    for vi, nu in enumerate(mus):
        for yi, muy in enumerate(mus):
            for s in I:
                for i in range(0, 2*h+1):
                    w = apply_poly(B[nu], h - s - i, {muy: 1})
                    if not w: continue
                    for xi, mux in enumerate(mus):
                        res = apply_poly(a_bracket(mux, i), s - h + i, w)
                        for ui, muu in enumerate(mus):
                            c = res.get(muu, 0)
                            if c: T[s][xi][yi][ui][vi] += (-1)**i * c
    return N, mus, I, T
def rr(T4, N):
    return flint.fmpz_mat(N*N, N*N, [T4[x][y][u][v] for x in range(N) for y in range(N) for u in range(N) for v in range(N)]).rank()
def rc(T4, N):
    return flint.fmpz_mat(N*N, N*N, [T4[x][y][u][v] for x in range(N) for u in range(N) for y in range(N) for v in range(N)]).rank()
if __name__ == "__main__":
    h = int(sys.argv[1]); Nreq = int(sys.argv[2]) if len(sys.argv) > 2 else None
    t0 = time.time(); N, mus, I, T = pieces(h, Nreq); d = math.factorial(h)
    ok = all(sum(T[s][x][y][u][v] for s in I) == d*(x == v)*(y == u) for x in range(N) for y in range(N) for u in range(N) for v in range(N))
    print(f"h={h} N={N} d=h!={d}: sum_s A_s == d X exactly: {ok}   ({time.time()-t0:.1f}s)", flush=True)
    assert ok
    p = [len(list(partitions(w))) for w in range(4*h+2)]; P = [sum(p[:w+1]) for w in range(4*h+2)]
    Pc = lambda w: 0 if w < 0 else P[w]
    print(" s   rr(A_s) rc(A_s) | s0  rr(R_s0) rc(C_s0)  [count bounds P(2h+s0-1), P(2h-s0)]")
    res = []
    Z = lambda: [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    R = Z()
    for s0 in range(-2*h, 2*h+2):
        C = [[[[d*(x == v)*(y == u) - R[x][y][u][v] for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
        a, b = rr(R, N), rc(C, N)
        sa = (rr(T[s0], N), rc(T[s0], N)) if s0 in T else ("-", "-")
        print(f"{s0:>3} {str(sa[0]):>7} {str(sa[1]):>7} | {s0:>3} {a:>8} {b:>8}   [{min(Pc(2*h+s0-1), N*N)}, {min(Pc(2*h-s0), N*N)}]", flush=True)
        res.append((s0, a, b))
        if s0 in T:
            for x in range(N):
                for y in range(N):
                    for u in range(N):
                        for v in range(N): R[x][y][u][v] += T[s0][x][y][u][v]
    bal = min(max(a, b) for _, a, b in res)
    print(f"EXACT threshold frontier of App B at h={h}, N={N}: {sorted(set((a,b) for _,a,b in res))}; balanced value {bal}; N^2 = {N*N}")
    json.dump({"h": h, "N": N, "frontier": res, "balanced": bal}, open(f"appb_exact_h{h}_N{N}.json", "w"))
