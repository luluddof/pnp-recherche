import math, sys
from profiles import series, cum, analyse
W = 1300; U = 2
for n in (8, 16):
    NS = cum(series([('f', 1, 2, n)], W)); R = cum(series([('f', 2, 2, n)], W, pref=2**(n//2)))
    for nm, M in (("NS", NS), ("Ramond", R)):
        print(f"n={n} V=NS M={nm}:  h, lam, exact ratio, leading prediction sqrt(2 + n/(48u)), u = h - n/48")
        for h2 in (40, 80, 160, 240, 320):
            N = NS[h2]; lam = math.log2(N); E0 = next(E for E in range(len(M)) if M[E] >= N); B = 2*E0 + 2*h2
            if B >= len(M): break
            r = max(math.log2(M[x]) + math.log2(M[B-x]) for x in range(B+1))/(2*lam); u = h2/2 - n/48
            print(f"   h={h2/2:6.1f} lam={lam:7.1f} E0={E0/2:6.1f}  ratio={r:.4f}  pred={math.sqrt(2 + n/(48*u)):.4f}")
