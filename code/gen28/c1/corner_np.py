# numpy version of corner_fast.py (mod prime P ~ 2^20, int64 matmul), w0 = 0..w0max.
import sys, math, time, flint, json
import numpy as np
from appb_mc import cparts, apply_poly, poly_mul, zmu, wt
from math import comb
def run(h, m, w0max, P=1048573):
    mus = sorted(cparts(h, m)); N = len(mus); idx = {mu: k for k, mu in enumerate(mus)}
    L = 1
    for mu in mus: L = L*zmu(mu)//math.gcd(L, zmu(mu))
    A = {}
    for mu in mus:
        a = {(): 1}
        for (n, c) in mu: a = poly_mul(a, {((j, c),): (-1)**(n-1)*comb(n, j) for j in range(1, n+1)})
        A[mu] = a
    G = {}; F = {}
    for wp in range(0, w0max+1):
        sigs = list(cparts(wp, m)); sidx = {s: k for k, s in enumerate(sigs)}
        Gw = np.zeros((len(sigs), N*N), dtype=np.int64)
        for vi, nu in enumerate(mus):
            bv = {nu: L//zmu(nu)}
            for yi, muy in enumerate(mus):
                for s, c in apply_poly(bv, wp-h, {muy: 1}).items(): Gw[sidx[s], yi*N+vi] = (Gw[sidx[s], yi*N+vi] + c) % P
        G[wp] = Gw
        # F[wp][i] : (x,u) x sig  for abar parts with coefficient C(h'+h,i) -- store per h' to combine later
        Fw = {}
        for xi, mux in enumerate(mus):
            for k, c in A[mux].items():
                hp = wt(k)
                for si, sg in enumerate(sigs):
                    for uu, c2 in apply_poly({k: c}, h-wp, {sg: 1}).items():
                        Fw.setdefault(hp, np.zeros((N*N, len(sigs)), dtype=np.int64))
                        Fw[hp][xi*N+idx[uu], si] = (Fw[hp][xi*N+idx[uu], si] + c2) % P
        F[wp] = Fw
    Xm = np.zeros((N*N, N*N), dtype=np.int64)  # in (xy)|(uv)
    for x in range(N):
        for y in range(N): Xm[x*N+y, y*N+x] = L % P
    res = []
    for w0 in range(0, w0max+1):
        t0 = time.time()
        Cm = np.zeros((N*N, N*N), dtype=np.int64)  # (x,u) | (y,v)
        for wp in range(0, w0+1):
            for hp, Fh in F[wp].items():
                cc = sum((-1)**i*comb(hp+h, i) for i in range(0, min(w0-wp, 2*h)+1)) % P
                if cc: Cm = (Cm + cc * ((Fh @ G[wp]) % P)) % P
        osr = flint.nmod_mat(Cm.tolist(), P).rank()
        # reorder C to (x,y)|(u,v): Cm[x*N+u, y*N+v] -> C4[x,u,y,v] -> [x,y,u,v]
        C4 = Cm.reshape(N, N, N, N).transpose(0, 2, 1, 3).reshape(N*N, N*N)
        Rm = (Xm - C4) % P
        rk = flint.nmod_mat(Rm.tolist(), P).rank()
        eta = math.log(osr)/math.log(N) if osr > 0 else 0; dl = 2 - math.log(rk)/math.log(N)
        print(f"h={h} m={m} N={N} w0={w0}: rank R={rk} osr C={osr} (eta={eta:.3f}, delta={dl:.4f}, N^2-p={N*N-rk}) {time.time()-t0:.1f}s", flush=True)
        res.append((w0, rk, osr))
    return N, res
if __name__ == "__main__":
    h, m, w0max = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    P = int(sys.argv[4]) if len(sys.argv) > 4 else 1048573
    N, res = run(h, m, w0max, P)
    json.dump({"h": h, "m": m, "N": N, "P": P, "corner": res}, open(f"cornernp_h{h}_m{m}_P{P}.json", "w"))
