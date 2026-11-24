# Exact test of the GENERALISED Appendix-B pieces on a MODULE of the Heisenberg vertex algebra:
#   Fock module pi_gamma (zero mode nabla_0 = gamma; gamma = 0 is the vacuum module), labels y,u = monomials of weight E0 (any E0),
#   states a_mu(x), b_mu(v) of weight h as in Lemma 10.  Pieces A_s = [t^u] Q_s(a_x, b_v) t^y, |s| <= E0 + h.
# Checks sum_s A_s = h! [x=v][y=u] EXACTLY and compares exact ranks of the threshold decompositions with the window
# dimensions P(E0+h+s0-1), P(E0+h-s0)  (budget B = 2 E0 + H, H = 2h).   Adapted from gen20/g20b/appb_exact.py.
import sys, math, functools, flint, time
from math import comb
GAMMA = 0
def gbinom(k, j):
    if j < 0: return 0
    num = 1
    for t in range(j): num *= (k - t)
    return num // math.factorial(j)
def partitions(n, maxpart=None):
    if maxpart is None: maxpart = n
    if n == 0: yield (); return
    for k in range(min(n, maxpart), 0, -1):
        for rest in partitions(n-k, k): yield (k,) + rest
def mono_mul_t(sig, l): return tuple(sorted(sig + (l,), reverse=True))
def nabla(sig, k):
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
    E = sum(sig)
    if E + r < 0: return ()
    if len(mu) == 0: return ((sig, 1),) if r == 0 else ()
    n, rest = mu[0], mu[1:]
    out = {}
    for l in range(n, E + r + 1):
        c = comb(l-1, n-1)
        for s2, c2 in apply_mono(rest, r - l, sig): addto(out, mono_mul_t(s2, l), c*c2)
    for k in set(sig):
        c = gbinom(-k-1, n-1); cc, s1 = nabla(sig, k)
        for s2, c2 in apply_mono(rest, r + k, s1): addto(out, s2, c*cc*c2)
    if GAMMA:   # zero mode: coefficient C(-1, n-1) = (-1)^(n-1), nabla_0 = GAMMA (scalar)
        c = gbinom(-1, n-1)*GAMMA
        for s2, c2 in apply_mono(rest, r, sig): addto(out, s2, c*c2)
    return tuple(out.items())
def apply_poly(a, r, f):
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
def pieces(h, E0):
    mus = sorted(partitions(h)); labs = sorted(partitions(E0)); N, Np = len(mus), len(labs)
    def y(n): return {(j,): comb(n, j) for j in range(1, n+1)}
    A = {}
    for mu in mus:
        a = {(): 1}
        for m in mu: a = poly_mul(a, {k: (-1)**(m-1)*c for k, c in y(m).items()})
        A[mu] = a
    def a_bracket(mu, i): return {k: comb(sum(k)+h, i)*c for k, c in A[mu].items() if comb(sum(k)+h, i)}
    B = {nu: {nu: math.factorial(h)//zmu(nu)} for nu in mus}
    I = list(range(-(E0+h), E0+h+1))
    T = {s: {} for s in I}
    for vi, nu in enumerate(mus):
        for yi, muy in enumerate(labs):
            for s in I:
                for i in range(0, 2*h+1):
                    w = apply_poly(B[nu], h - s - i, {muy: 1})
                    if not w: continue
                    for xi, mux in enumerate(mus):
                        res = apply_poly(a_bracket(mux, i), s - h + i, w)
                        for ui, muu in enumerate(labs):
                            c = res.get(muu, 0)
                            if c: addto(T[s], (xi, yi, ui, vi), (-1)**i * c)
    return N, Np, I, T
def rr(T, N, Np):
    M = [[0]*(Np*N) for _ in range(N*Np)]
    for (x, y, u, v), c in T.items(): M[x*Np+y][u*N+v] = c
    return flint.fmpz_mat(M).rank()
def rc(T, N, Np):
    M = [[0]*(Np*N) for _ in range(N*Np)]
    for (x, y, u, v), c in T.items(): M[x*Np+u][y*N+v] = c
    return flint.fmpz_mat(M).rank()
if __name__ == "__main__":
    h = int(sys.argv[1]); E0 = int(sys.argv[2]); GAMMA = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    t0 = time.time(); N, Np, I, T = pieces(h, E0); d = math.factorial(h)
    tot = {}
    for s in I:
        for k, c in T[s].items(): addto(tot, k, c)
    X = {(x, y, y, x): d for x in range(N) for y in range(Np)}
    ok = (tot == X)
    print(f"h={h} E0={E0} gamma={GAMMA}: N={N} states, N'={Np} labels, budget B=2E0+2h={2*E0+2*h}; sum_s A_s == h! [x=v][y=u] exactly: {ok}  ({time.time()-t0:.1f}s)", flush=True)
    assert ok
    p = [len(list(partitions(w))) for w in range(2*E0+2*h+2)]; P = [sum(p[:w+1]) for w in range(len(p))]
    Pc = lambda w: 0 if w < 0 else P[w]
    R = {}; rows = []
    for s0 in range(-(E0+h), E0+h+2):
        C = dict(X)
        for k, c in R.items(): addto(C, k, -c)
        a, b = rr(R, N, Np), rc(C, N, Np)
        rows.append((s0, a, b, min(Pc(E0+h+s0-1), N*Np), min(Pc(E0+h-s0), N*Np)))
        if s0 in T:
            for k, c in T[s0].items(): addto(R, k, c)
    print(" s0: rr(R_s0) rc(C_s0) [window bounds P(E0+h+s0-1), P(E0+h-s0), capped at N N'] ")
    print("   " + "  ".join(f"{s0}:{a},{b}[{pa},{pb}]" for s0, a, b, pa, pb in rows))
    viol = [r for r in rows if r[1] > r[3] or r[2] > r[4]]
    print(f"   all exact ranks <= window dimensions: {not viol};  min over s0 of rr*rc = {min(a*b for _,a,b,_,_ in rows)} ; N*N' = {N*Np};  balanced max(rr,rc) = {min(max(a,b) for _,a,b,_,_ in rows)}")
