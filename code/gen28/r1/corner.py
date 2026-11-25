# Corner of the singlet family from exact integer counts at large k (q taken as the closed-form upper bound; p exact
# outside-window mass; for a window of d+1 integers the interpolant cannot agree off the window except by accident, which
# only lowers p).  delta = 2 - log2(p)/k, eta = log2(q)/k.
from math import comb, log2, log
import sys
for k in [300, 1000, 3000]:
    mass = [comb(k, m)*3**(k-m) for m in range(k+1)]
    pre = [0]
    for x in mass: pre.append(pre[-1] + x)
    tot = pre[-1]
    qcum = 0; rows = []
    for d in range(0, k//2):
        qcum += comb(k, d)*3**d
        best = max(range(0, k-d+1), key=lambda a: pre[a+d+1] - pre[a])
        p = tot - (pre[best+d+1] - pre[best])
        eta = log2(qcum)/k; delta = 2 - log2(p)/k
        rows.append((d/k, eta, delta))
    print(f"k={k}")
    for tau in [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]:
        d = max(1, int(tau*k)); t, eta, delta = rows[d]
        print(f"  tau={t:.4f} eta={eta:.4f} delta={delta:.6f} delta/eta^2={delta/eta**2:.4f} "
              f"delta*log2(1/eta)^2/eta^2={delta*log2(1/eta)**2/eta**2:.3f}  eta^2/4={eta**2/4:.5f}", flush=True)
