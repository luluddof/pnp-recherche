# For which hypothetical single seeds (p,q)@N (semisimple) does h_LB(mu) = log_{N^2}(p^mu (N^2-q)^(1-mu) + (N^2-p)^mu q^(1-mu))
# drop below the circle for some mu?  Scan all integer p,q (or a log grid for large N); report the region and Theorem C status.
import numpy as np, math
mus = np.linspace(0.0025, 0.9975, 399); hB = np.hypot(mus, 1-mus)
def lb(N, p, q):
    N2 = N*N; return np.log(p**mus*(N2-q)**(1-mus) + (N2-p)**mus*q**(1-mus))/math.log(N2)
for N in [2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 32, 64, 128, 1024, 2**16, 2**32]:
    N2 = N*N
    if N <= 128: grid = np.arange(1, N2)
    else: grid = np.unique(np.round(np.exp(np.linspace(0, math.log(N2-1), 3000))).astype(np.int64))
    best = None; beat = []
    P, Q = np.meshgrid(grid.astype(float), grid.astype(float), indexing='ij')
    # vectorised over mu in chunks
    gmin = np.full(P.shape, np.inf); argmu = np.zeros(P.shape)
    for m, hb in zip(mus, hB):
        h = np.log(P**m*(N2-Q)**(1-m) + (N2-P)**m*Q**(1-m))/math.log(N2); g = h - hb
        upd = g < gmin; gmin[upd] = g[upd]; argmu[upd] = m
    viol = gmin < -1e-12
    if viol.any():
        prod = (P*Q)[viol]/N2
        i = np.argmax(prod); pv, qv = P[viol][i], Q[viol][i]
        print(f"N={N}: {viol.sum()} (p,q) pairs beat the circle; max pq/N^2 among them = {prod.max():.4f} at (p,q)=({pv:.0f},{qv:.0f})"
              f" [Theorem C needs pq/N^2 >= 2/9 = 0.2222]; max (log_N p + log_N q) among them = {np.max(np.log(P[viol]*Q[viol]))/math.log(N):.4f}")
    else:
        print(f"N={N}: no (p,q) in [1,N^2-1]^2 beats the circle")
