# Effective per-level gain of the App B pieces at finite N = p(h): rho_eff = max_s log(rr(A_s) rc(A_s)) / (4 log N),
# with the counting bounds rr(A_s) <= min(N^2, sum_{w=s}^{s+2h} p(w)), rc(A_s) <= min(N^2, sum_{w=-s}^{2h-s} p(w)), and |I| = 4h+1.
import math
M = 2600; p = [1] + [0]*M
for k in range(1, M+1):
    for n in range(k, M+1): p[n] += p[n-k]
pc = lambda w: p[w] if 0 <= w <= M else 0
print(" h        N=p(h)    rho_eff   (1/sqrt2 = 0.7071)")
for h in [2,4,6,8,10,12,15,20,30,50,100,200,400,600]:
    N = p[h]; best = 0
    for s in range(-2*h, 2*h+1):
        rr = min(N*N, sum(pc(w) for w in range(s, s+2*h+1))); rc = min(N*N, sum(pc(w) for w in range(-s, 2*h-s+1)))
        best = max(best, math.log(rr*rc))
    print(f"{h:>3} {N:>14.4g} {best/(4*math.log(N)):10.4f}")
