# Appendix B over GF(p) with h >= p (h! = 0 mod p).  Replacement of Lemma 10: Lemma 9(c) is an integer identity,
#   sum_s Q_s(a,b) = (a[-h] b) Id,
# so it holds mod p; choose a_x (combinations of monomials of weight <= h) and b_v = t^{nu_v} (nu_v in S) with
# a_x[-h] b_v = [x = v] mod p.  This is possible exactly for N = rank_p of the pairing = d_p(h) (appb_gram.py).
# y-states t^{nu_y} (nu_y in S), u-functionals f_u(w) = a_u[-h] w.  Then sum_s A_s = X (mod p) with d = 1, where
#   A_s(x,y,u,v) = f_u( Q_s(a_x, b_v) t^{nu_y} ).
# Computes exact ranks mod p of all thresholds R_{s0} = sum_{s<s0} A_s, C = X - R, and checks sum_s A_s = X mod p.
import sys, math, itertools, importlib.util, json, time
spec = importlib.util.spec_from_file_location("ab", "/home/user/pnp-recherche/code/gen20/g20b/appb_exact.py")
ab = importlib.util.module_from_spec(spec); spec.loader.exec_module(ab)
from verify_gf import rank_mod
def inv_mod_matrix(A, p):
    n = len(A); M = [[a % p for a in row] + [int(i == j) for j in range(n)] for i, row in enumerate(A)]
    for c in range(n):
        piv = next(i for i in range(c, n) if M[i][c]); M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], p-2, p); M[c] = [(e*iv) % p for e in M[c]]
        for i in range(n):
            if i != c and M[i][c]:
                f = M[i][c]; M[i] = [(a - f*b) % p for a, b in zip(M[i], M[c])]
    return [row[n:] for row in M]
def run(h, p):
    t0 = time.time()
    rows = [mu for w in range(h+1) for mu in ab.partitions(w)]
    cols = sorted(ab.partitions(h))
    Mg = [[dict(ab.apply_mono(mu, -h, nu)).get((), 0) % p for nu in cols] for mu in rows]
    d = rank_mod(Mg, p)
    # greedy choice of independent columns S and rows R0 (mod p)
    S = []
    for j in range(len(cols)):
        if rank_mod([[Mg[r][c] for c in S + [j]] for r in range(len(rows))], p) > len(S): S.append(j)
    R0 = []
    for r in range(len(rows)):
        if rank_mod([[Mg[rr][c] for c in S] for rr in R0 + [r]], p) > len(R0): R0.append(r)
    assert len(S) == d == len(R0)
    Inv = inv_mod_matrix([[Mg[r][c] for c in S] for r in R0], p)   # Inv * M[R0][S] = I  => rows a_x = Inv[x] . t^{R0}
    # a_x as dict over monomials rows[R0[k]]: coefficient Inv[x][k] ... check: sum_k Inv[x][k] M[R0[k]][S[v]] = [x=v]
    # NB: Inv computed for M[R0][S] as a matrix (rows R0, cols S): we need left inverse L with L M = I.
    Mrs = [[Mg[r][c] for c in S] for r in R0]
    L = inv_mod_matrix(Mrs, p)  # inverse of square matrix: L*Mrs = I
    N = d
    a = [{rows[R0[k]]: L[x][k] for k in range(N) if L[x][k]} for x in range(N)]
    nu = [cols[j] for j in S]
    f = [[sum(L[u][k]*Mg[R0[k]][c] for k in range(N)) % p for c in range(len(cols))] for u in range(N)]
    cidx = {mu: c for c, mu in enumerate(cols)}
    monos = sorted(set(m for ax in a for m in ax))
    I = list(range(-2*h, 2*h+1))
    A = {s: [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)] for s in I}
    for mono in monos:
        hp = sum(mono)
        for s in I:
            for yi in range(N):
                for vi in range(N):
                    vec = {}
                    for i in range(2*h+1):
                        cb = math.comb(hp+h, i) % p
                        if not cb: continue
                        w = ab.apply_mono(nu[vi], h - s - i, nu[yi])
                        for sig, c1 in w:
                            for sig2, c2 in ab.apply_mono(mono, s - h + i, sig):
                                vec[sig2] = (vec.get(sig2, 0) + (-1)**i * cb * c1 * c2) % p
                    if not any(vec.values()): continue
                    for ui in range(N):
                        val = sum(f[ui][cidx[sg]]*cv for sg, cv in vec.items() if sg in cidx) % p
                        if not val: continue
                        for xi in range(N):
                            co = a[xi].get(mono, 0)
                            if co: A[s][xi][yi][ui][vi] = (A[s][xi][yi][ui][vi] + co*val) % p
    ok = all(sum(A[s][x][y][u][v] for s in I) % p == int(x == v and y == u) for x, y, u, v in itertools.product(range(N), repeat=4))
    print(f"h={h} p={p}: N = d_p(h) = {N} (p(h) = {len(cols)}); sum_s A_s == X mod p: {ok}  [{time.time()-t0:.0f}s]", flush=True)
    assert ok
    def rr(T): return rank_mod([[T[x][y][u][v] for u in range(N) for v in range(N)] for x in range(N) for y in range(N)], p)
    def rc(T): return rank_mod([[T[x][y][u][v] for y in range(N) for v in range(N)] for x in range(N) for u in range(N)], p)
    # restricted counts d_p(w), cumulative D_p
    dp = [sum(1 for mu in ab.partitions(w) if ab.zmu(mu) % p) for w in range(4*h+2)]
    Dp = [sum(dp[:w+1]) for w in range(4*h+2)]
    R = [[[[0]*N for _ in range(N)] for _ in range(N)] for _ in range(N)]
    res = []
    for s0 in range(-2*h, 2*h+2):
        C = [[[[(int(x == v and y == u) - R[x][y][u][v]) % p for v in range(N)] for u in range(N)] for y in range(N)] for x in range(N)]
        r1, c1 = rr(R), rc(C); res.append((s0, r1, c1))
        Dr = Dp[2*h+s0-1] if 2*h+s0-1 >= 0 else 0; Dc = Dp[2*h-s0] if 2*h-s0 >= 0 else 0
        print(f"  s0={s0:>3}: rr(R)={r1:>4} rc(C)={c1:>4}   [restricted counts D_p(2h+s0-1)={Dr}, D_p(2h-s0)={Dc}]", flush=True)
        if s0 in A:
            for x, y, u, v in itertools.product(range(N), repeat=4): R[x][y][u][v] = (R[x][y][u][v] + A[s0][x][y][u][v]) % p
    bal = min(max(r, c) for _, r, c in res)
    print(f"h={h} p={p} N={N}: balanced threshold value {bal} (N^2 = {N*N})", flush=True)
    json.dump({"h": h, "p": p, "N": N, "frontier": res, "balanced": bal}, open(f"appb_modp_h{h}_p{p}.json", "w"))
if __name__ == "__main__":
    run(int(sys.argv[1]), int(sys.argv[2]))
