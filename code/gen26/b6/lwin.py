# Idea 2: second grading (number of parts l(sigma)) windows inside the envelope: S = {|sigma|<=W, l(sigma)<=L}.
# Looks for any feasible (dimU, dimU') with dimU + dimU' < N^2 (below the trivial antidiagonal).
import sys, time
sys.dont_write_bytecode = True
from relax import setup, g_rows, gp_rows, feasible, basis_mod, PR
h = int(sys.argv[1]); N, mus, bs = setup(h); n = N*N; p = PR[0]
wmax = h + 6
G = g_rows(h, N, mus, bs, wmax); Gp = gp_rows(h, N, mus, bs, wmax)
def sub(D, W, L): return [D[s] for s in D if sum(s) <= W and len(s) <= L]
best = None; pts = []
t0 = time.time()
for W in range(0, wmax+1):
    for L in range(1, W+1):
        Uv = sub(G, W, L); dU = len(basis_mod(Uv, n, p))
        if dU >= n: continue
        for Lp in range(1, wmax+1):
            lo, hi = 0, wmax
            if not feasible(N, Uv, sub(Gp, wmax, Lp), p)[0]: continue
            while lo < hi:
                mid = (lo+hi)//2
                if feasible(N, Uv, sub(Gp, mid, Lp), p)[0]: hi = mid
                else: lo = mid+1
            dUp = len(basis_mod(sub(Gp, lo, Lp), n, p))
            pts.append((dU+dUp-n, dU, dUp, W, L, lo, Lp))
    if pts:
        m = min(pts); print(f"W<={W} ({time.time()-t0:.0f}s): best excess dU+dU'-N^2 = {m[0]}  at dims ({m[1]},{m[2]}) (W,L,W',L')={m[3:]}", flush=True)
pts.sort(); print("10 best:", pts[:10])
