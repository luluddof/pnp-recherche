# App B pieces for the rank-m Heisenberg algebra (m colours). Parts are (n,c). Exact over Z; ranks over Q (flint)
# Checks sum_s A_s = L*X exactly, then exact threshold frontier (rank R_s0 in (xy)|(uv), osr C_s0 in (xu)|(yv)).
import sys, math, functools, flint, json, time
from math import comb
def gbinom(k, j):
    if j < 0: return 0
    num = 1
    for t in range(j): num *= (k - t)
    return num // math.factorial(j)
def cparts(n, m, maxpart=None):
    # coloured partitions: nonincreasing sequences of (size,colour)
    if maxpart is None: maxpart = (n, m-1)
    if n == 0: yield (); return
    for k in range(min(n, maxpart[0]), 0, -1):
        for c in range(m-1, -1, -1):
            if (k, c) > maxpart: continue
            for rest in cparts(n-k, m, (k, c)): yield ((k, c),) + rest
def srt(t): return tuple(sorted(t, reverse=True))
def wt(sig): return sum(k for k, c in sig)
def addto(d, key, c):
    if c:
        v = d.get(key, 0) + c
        if v: d[key] = v
        else: d.pop(key, None)
@functools.lru_cache(maxsize=None)
def apply_mono(mu, r, sig):
    E = wt(sig)
    if E + r < 0: return ()
    if len(mu) == 0: return ((sig, 1),) if r == 0 else ()
    (n, col), rest = mu[0], mu[1:]
    out = {}
    for l in range(n, E + r + 1):
        cf = comb(l-1, n-1)
        for s2, c2 in apply_mono(rest, r - l, sig):
            addto(out, srt(s2 + ((l, col),)), cf*c2)
    for part in set(sig):
        k, c = part
        if c != col: continue
        cf = gbinom(-k-1, n-1); mult = sig.count(part)
        lst = list(sig); lst.remove(part); s1 = tuple(lst)
        for s2, c2 in apply_mono(rest, r + k, s1):
            addto(out, s2, cf*k*mult*c2)
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
        for b, cb in g.items(): addto(out, srt(a+b), ca*cb)
    return out
def zmu(mu):
    z = 1
    for p in set(mu): mlt = mu.count(p); z *= p[0]**mlt * math.factorial(mlt)
    return z
def pieces(h, m):
    mus = sorted(cparts(h, m)); N = len(mus)
    L = 1
    for mu in mus: L = L*zmu(mu)//math.gcd(L, zmu(mu))
    A = {}
    for mu in mus:
        a = {(): 1}
        for (n, c) in mu:
            a = poly_mul(a, {((j, c),): (-1)**(n-1)*comb(n, j) for j in range(1, n+1)})
        A[mu] = a
    def a_br(mu, i): return {k: comb(wt(k)+h, i)*c for k, c in A[mu].items() if comb(wt(k)+h, i)}
    B = {nu: {nu: L//zmu(nu)} for nu in mus}
    I = list(range(-2*h, 2*h+1))
    T = {s: [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)] for s in I}
    for vi, nu in enumerate(mus):
        for yi, muy in enumerate(mus):
            for s in I:
                for i in range(0, 2*h+1):
                    w = apply_poly(B[nu], h - s - i, {muy: 1})
                    if not w: continue
                    for xi, mux in enumerate(mus):
                        res = apply_poly(a_br(mux, i), s - h + i, w)
                        for ui, muu in enumerate(mus):
                            c = res.get(muu, 0)
                            if c: T[s][xi][yi][ui][vi] += (-1)**i * c
    return N, mus, I, T, L
def rr(T4, N): return flint.fmpz_mat(N*N, N*N, [T4[x][y][u][v] for x in range(N) for y in range(N) for u in range(N) for v in range(N)]).rank()
def rc(T4, N): return flint.fmpz_mat(N*N, N*N, [T4[x][y][u][v] for x in range(N) for u in range(N) for y in range(N) for v in range(N)]).rank()
if __name__ == "__main__":
    h = int(sys.argv[1]); m = int(sys.argv[2])
    t0 = time.time(); N, mus, I, T, d = pieces(h, m)
    ok = all(sum(T[s][x][y][u][v] for s in I) == d*(x == v)*(y == u) for x in range(N) for y in range(N) for u in range(N) for v in range(N))
    print(f"h={h} m={m} N={N} L={d}: sum_s A_s == L X exactly: {ok} ({time.time()-t0:.1f}s)", flush=True)
    assert ok
    pm = [len(list(cparts(w, m))) for w in range(4*h+2)]; P = [sum(pm[:w+1]) for w in range(4*h+2)]
    Pc = lambda w: 0 if w < 0 else P[min(w, 4*h)]
    R = [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)]; res = []
    for s0 in range(-2*h, 2*h+2):
        C = [[[[d*(x == v)*(y == u) - R[x][y][u][v] for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
        a, b = rr(R, N), rc(C, N)
        print(f"s0={s0:>3} rank R={a:>6} osr C={b:>6}  counts [{min(Pc(2*h+s0-1), N*N)}, {min(Pc(2*h-s0), N*N)}]", flush=True)
        res.append((s0, a, b))
        if s0 in T:
            for x in range(N):
                for y in range(N):
                    for u in range(N):
                        for v in range(N): R[x][y][u][v] += T[s0][x][y][u][v]
    print(f"frontier h={h} m={m} N={N}: {sorted(set((a,b) for _,a,b in res))}  N^2={N*N}")
    json.dump({"h": h, "m": m, "N": N, "frontier": res}, open(f"appb_mc_h{h}_m{m}.json", "w"))
