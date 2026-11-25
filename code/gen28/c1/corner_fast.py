# Fast exact corner thresholds of (multi-colour) App B: C_{s0} with s0 = 2h - w0 only needs C-channel weights <= w0:
#   C(x,y,u,v) = sum_{w'<=w0} sum_{|sig|=w'} <u| at_x^{(w')}[h-w'] sig> <sig| b_v[w'-h] y>,  at^{(w')} = sum_{i<=w0-w'} (-1)^i a<i>
# R = L X - C.  Ranks mod a prime P (flint nmod_mat): rank R in (xy)|(uv), osr C in (xu)|(yv).
import sys, math, time, flint, json
from appb_mc import cparts, apply_poly, poly_mul, zmu, wt
from math import comb
def run(h, m, w0max, P=2147483629, check_full=False):
    mus = sorted(cparts(h, m)); N = len(mus); idx = {mu: k for k, mu in enumerate(mus)}
    L = 1
    for mu in mus: L = L*zmu(mu)//math.gcd(L, zmu(mu))
    A = {}
    for mu in mus:
        a = {(): 1}
        for (n, c) in mu: a = poly_mul(a, {((j, c),): (-1)**(n-1)*comb(n, j) for j in range(1, n+1)})
        A[mu] = a
    def at(mu, w0, wp):  # sum_{i=0}^{w0-wp} (-1)^i a<i>,  a<i> = C(h'+h, i) abar_h'
        out = {}
        for k, c in A[mu].items():
            cc = sum((-1)**i*comb(wt(k)+h, i) for i in range(0, min(w0-wp, 2*h)+1))
            if cc*c: out[k] = cc*c
        return out
    res = []
    # G[w'][sig][(y,v)] = [t^sig] b_v[w'-h] t^y   (independent of w0)
    G = {}
    for wp in range(0, w0max+1):
        sigs = list(cparts(wp, m)); sidx = {s: k for k, s in enumerate(sigs)}
        Gw = [[0]*(N*N) for _ in sigs]
        for vi, nu in enumerate(mus):
            bv = {nu: L//zmu(nu)}
            for yi, muy in enumerate(mus):
                for s, c in apply_poly(bv, wp-h, {muy: 1}).items(): Gw[sidx[s]][yi*N+vi] += c
        G[wp] = (sigs, Gw)
    for w0 in range(1, w0max+1):
        t0 = time.time()
        Cm = [[0]*(N*N) for _ in range(N*N)]   # rows (x,u), cols (y,v)
        for wp in range(0, w0+1):
            sigs, Gw = G[wp]
            for xi, mux in enumerate(mus):
                ax = at(mux, w0, wp)
                for si, sg in enumerate(sigs):
                    img = apply_poly(ax, h-wp, {sg: 1})
                    if not img: continue
                    g = Gw[si]
                    for uu, c in img.items():
                        row = Cm[xi*N+idx[uu]]
                        for col in range(N*N):
                            if g[col]: row[col] += c*g[col]
        osr = flint.nmod_mat(N*N, N*N, [Cm[r][c] % P for r in range(N*N) for c in range(N*N)], P).rank()
        # R in (xy)|(uv): R[(x,y),(u,v)] = L[x=v][y=u] - C[(x,u),(y,v)]
        Rl = []
        for x in range(N):
            for y in range(N):
                for u in range(N):
                    for v in range(N):
                        Rl.append((L*(x == v)*(y == u) - Cm[x*N+u][y*N+v]) % P)
        rk = flint.nmod_mat(N*N, N*N, Rl, P).rank()
        eta = math.log(osr)/math.log(N) if osr > 0 else 0; dl = 2 - math.log(rk)/math.log(N)
        print(f"h={h} m={m} N={N} w0={w0}: rank R={rk} osr C={osr}  (eta={eta:.3f}, delta={dl:.4f}, N^2-p={N*N-rk})  {time.time()-t0:.1f}s", flush=True)
        res.append((w0, rk, osr))
    return N, res
if __name__ == "__main__":
    h, m, w0max = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    N, res = run(h, m, w0max)
    json.dump({"h": h, "m": m, "N": N, "corner": res}, open(f"corner_h{h}_m{m}.json", "w"))
