# asymptotic and finite-N (counting) curves of the Johnson-regraded Hrubes construction (JH) vs the circle and App. B
import math
from math import comb, log, sqrt
from cliff import rho
def H(x): return 0 if x <= 0 or x >= 1 else -x*math.log2(x) - (1-x)*math.log2(1-x)
def Hinv(y):
    lo, hi = 0, 0.5
    for _ in range(100):
        m = (lo+hi)/2
        if H(m) < y: lo = m
        else: hi = m
    return lo
al = Hinv(0.5)
print(f"alpha = H^-1(1/2) = {al:.5f};  JH curve: theta = 1 + 2H(b), T = 1 + 2H(alpha - b), b in [0, alpha]")
print("   b      theta     T_JH    T_circle   theta+T_JH")
for i in range(0, 11):
    b = al * i / 10; th, T = 1 + 2*H(b), 1 + 2*H(al - b)
    print(f" {b:.4f}  {th:.4f}  {T:.4f}  {sqrt(max(0,4-th*th)):.4f}  {th+T:.4f}")
# T_JH at theta = sqrt2
lo, hi = 0, al/2
for _ in range(100):
    m = (lo+hi)/2
    if 1 + 2*H(m) < sqrt(2): lo = m
    else: hi = m
print(f"at theta = sqrt2: b = {lo:.5f}, T_JH = {1+2*H(al-lo):.4f} (circle {sqrt(2):.4f}); balanced JH = {1+2*H(al/2):.4f}")
print(f"sup(theta+T) on JH = {max(2+2*H(al*i/1000)+2*H(al-al*i/1000) for i in range(1001)):.4f}  (circle 2 sqrt2 = {2*sqrt(2):.4f})")
# finite N counting table (upper bounds n*min(C(t,j),N), exact for the full slice; exact ranks computed for N <= 256 in johnson.py)
import sys
sys.setrecursionlimit(10000)
M = 2400
p = [1] + [0]*M
for k in range(1, M+1):
    for n in range(k, M+1): p[n] += p[n-k]
P = [0]*(M+1); acc = 0
for w in range(M+1): acc += p[w]; P[w] = acc
Pc = lambda w: 0 if w < 0 else P[w]
def appb(N):
    h = next((h for h in range(M//4) if p[h] >= N), None)
    if h is None: return None
    if 4*h + 2 > M: return None
    return [(min(Pc(2*h+s0-1), N*N), min(Pc(2*h-s0), N*N)) for s0 in range(-2*h, 2*h+2)]
def appb_T(N, theta):  # min q over App B frontier with p <= N^theta
    F = appb(N)
    if F is None: return None
    qs = [b for a, b in F if a <= N**theta]
    return min(qs) if qs else None
print("\nfinite N (N = n = 2^r): JH frontier (count) vs App. B count frontier;  log_N values")
print("   r   t   k  bal_JH   bal_AppB |  T_JH(theta<=sqrt2)  T_AppB(theta<=sqrt2)")
for r in [3,4,5,6,7,8,10,12,16,20,24,32,36,40,44,50,64,80,100,128,160]:
    n = 2**r; N = n; t = rho(n) - 1
    k = next(k for k in range(t+1) if comb(t, k) >= N)
    pts = []
    for j1 in range(-1, k+1):
        j2 = k - 1 - j1
        r1 = 0 if j1 < 0 else min(comb(t, j1), N); r2 = 0 if j2 < 0 else min(comb(t, j2), N)
        pts.append((n*r1, n*r2))
    pts.append((N*(N-1)//2, 1))
    lg = lambda x: log(x)/log(N) if x > 0 else float('-inf')
    bal = min(max(a, b) for a, b in pts)
    TJ = min((b for a, b in pts if a <= N**sqrt(2)), default=None)
    F = appb(N); balB = min(max(a, b) for a, b in F) if F else None
    TB = appb_T(N, sqrt(2))
    print(f" {r:3} {t:3} {k:3}  {lg(bal):.4f}   {lg(balB) if balB else float('nan'):.4f}   |  {lg(TJ) if TJ else float('nan'):.4f}   {lg(TB) if TB else float('nan'):.4f}")
