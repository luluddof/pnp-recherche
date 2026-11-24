# Appendix B at finite N: the feasible pairs obtained by thresholding the pieces A_s (s in [-2h,2h]) at s0.
# Counting formula (exact upper bound on the ranks, from (54)/(55)):
#  in (54) the (u,v)-factor [t^mu(u)] b[h-|sigma|] t^sigma depends on sigma only (not on i,s separately), so
#  R_{s0} = sum_{s<s0} A_s has rr <= #{sigma : |sigma| <= 2h+s0-1} = P(2h+s0-1), P(w) = sum_{j<=w} p(j);
#  in (55) the (x,u)-factor [t^mu(u)] a<i>[..] t^sigma ... sorry: (y,v)-factor [t^sigma] b[h-s-i] t^mu(y) depends on sigma only,
#  |sigma| = 2h-s-i, so C_{s0} = sum_{s>=s0} A_s has rc <= P(2h-s0).  Both also <= N^2 trivially.
# N = p(h) (all partitions of h used), or N = 2^lam (h least with p(h) >= N; restriction to a subset of indices).
import math
M = 2400
p = [1] + [0]*M
for k in range(1, M+1):
    for n in range(k, M+1): p[n] += p[n-k]
P = [0]*(M+1); acc = 0
for w in range(M+1): acc += p[w]; P[w] = acc
def Pc(w): return 0 if w < 0 else P[w]
def frontier(h, N):
    pts = []
    for s0 in range(-2*h, 2*h+2):
        pp, qq = min(Pc(2*h+s0-1), N*N), min(Pc(2*h-s0), N*N)
        pts.append((pp, qq))
    return pts
def balanced(h, N):
    return min(max(a, b) for a, b in frontier(h, N))
def minpq(h, N):  # min of p*q/N^2 over the frontier excluding endpoints
    return min(a*b for a, b in frontier(h, N) if a > 0 and b > 0) / N**2
print("N = p(h): balanced value b_AppB = min_s0 max(P(2h+s0-1), P(2h-s0)) (capped at N^2)")
print(f"{'h':>3} {'N':>12} {'b_AppB':>14} {'cap?':>5} {'log b/log N':>12} {'b/N':>12}")
rows = []
for h in list(range(2, 13)) + [15, 20, 30, 40, 60, 80, 100, 150, 200, 300, 450, 600]:
    N = p[h]; b = balanced(h, N); raw = min(max(Pc(2*h+s0-1), Pc(2*h-s0)) for s0 in range(-2*h, 2*h+2))
    rows.append((h, N, b)); print(f"{h:>3} {N:>12} {b:>14} {'yes' if raw > N*N else '':>5} {math.log(b)/math.log(N):12.4f} {b/N:12.4g}")
print()
print("N = 2^lam: h = least with p(h) >= N; restriction to N of the p(h) partitions")
print(f"{'lam':>3} {'N':>8} {'h':>3} {'b_AppB':>12} {'cap?':>5} {'log b/log N':>12}  uncapped P(2h)  exponent")
for lam in list(range(1, 13)) + [16, 20, 24, 32, 40, 50, 64, 100, 200, 400, 800]:
    N = 2**lam; h = next(h for h in range(M) if p[h] >= N)
    if 4*h > M: break
    raw = min(max(Pc(2*h+s0-1), Pc(2*h-s0)) for s0 in range(-2*h, 2*h+2)); b = min(raw, N*N)
    print(f"{lam:>3} {N:>8} {h:>3} {b:>12} {'yes' if raw > N*N else '':>5} {math.log(b)/math.log(N):12.4f}  {raw:>14} {math.log(raw)/math.log(N):8.4f}")
print()
print("first N = p(h) at which the uncapped App B balanced value is < N^2 (i.e. App B beats the trivial point at all):")
for h in range(2, 100):
    N = p[h]; raw = min(max(Pc(2*h+s0-1), Pc(2*h-s0)) for s0 in range(-2*h, 2*h+2))
    if raw < N*N: print("h =", h, "N =", N, "b =", raw, "exp =", round(math.log(raw)/math.log(N), 4)); break
for target in (1.9, 1.8, 1.7, 1.6, 1.5, 1.45, 1.42):
    for h in range(2, M//4):
        N = p[h]; raw = min(max(Pc(2*h+s0-1), Pc(2*h-s0)) for s0 in range(-2*h, 2*h+2))
        if math.log(raw)/math.log(N) < target: print(f"exponent < {target}: first at h = {h}, N = p(h) = {N:.4g}"); break
    else: print(f"exponent < {target}: not reached for h <= {M//4} (N up to {p[M//4]:.3g})")
