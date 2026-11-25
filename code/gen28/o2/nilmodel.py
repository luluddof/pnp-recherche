# Idealized Jordan-heavy seed (hypothetical): M_row nilpotent with p Jordan 2-blocks (+ N^2-2p zero 1-blocks), M_col semisimple
# with q nonzero eigenvalues (none equal 1).  Exact merged ranks of type class (k positions, j = k-i R-positions):
#   row = sum_{m>=j} C(k,m) p^m (N^2-2p)^(k-m) * sum_w min(C(m,w), C(m,w-j))     [inclusion matrices have full rank]
#   col = sum_{a'<=j} C(k,a') (N^2-q)^a' q^(k-a')
# Asymptotic exponents; compare h(mu) with the circle, scanning p <= N^2/2, q with Theorem C p q >= 2N^2/9.
import numpy as np, math
G = 300
def H(t):
    t = np.clip(t, 1e-300, 1); u = np.clip(1-t, 1e-300, 1); return -(t*np.log(t) + u*np.log(u))
def phis(N, p, q, xs):
    N2 = N*N; L0 = math.log(N2-2*p) if N2 > 2*p else -np.inf
    nu = np.linspace(0, 1, G+1)
    pr = []; pc = []
    for x in xs:
        y = 1-x; v = nu[nu >= y - 1e-12]
        with np.errstate(divide='ignore', invalid='ignore'):
            inner = np.where(v > 0, v*H(0.5 + np.where(v > 0, y/np.maximum(v, 1e-300), 0)/2), 0)
            val = H(v) + v*math.log(p) + np.where(v < 1, (1-v)*L0, 0) + inner
        pr.append(np.nanmax(val))
        a = nu[nu <= y + 1e-12]
        pc.append(np.max(H(a) + a*math.log(N2-q) + (1-a)*math.log(q)))
    return np.array(pr), np.array(pc)
xs = np.linspace(0, 1, 201); mus = np.linspace(0.005, 0.995, 199); hB = np.hypot(mus, 1-mus)
worst = None
for N in [2, 3, 4, 5, 6, 8, 10, 16, 32]:
    N2 = N*N; bestN = None
    for p in sorted(set(np.unique(np.round(np.geomspace(1, N2//2, 25)).astype(int)))):
        for q in sorted(set(np.unique(np.round(np.geomspace(1, N2-1, 25)).astype(int)))):
            if p*q < 2*N2/9: continue
            pr, pc = phis(N, p, q, xs)
            h = np.array([np.max(m*pr + (1-m)*pc) for m in mus])/math.log(N2)
            g = ((h-hB)/(mus*(1-mus))).min()
            if bestN is None or g < bestN[0]: bestN = (g, p, q, h[99])
    print(f"N={N}: nilpotent-row model, min over (p,q) [pq>=2N^2/9] of min_mu (h-circle)/(mu(1-mu)) = {bestN[0]:+.4f} at (p,q)=({bestN[1]},{bestN[2]}), h(1/2)={bestN[3]:.4f}", flush=True)
