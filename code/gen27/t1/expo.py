import math
from counting import p, P
pts = [  # (label, h, W, W')
 ("D=h  h=11 bal (k=3)", 11, 21, 19),
 ("D=inf h=4 (4,5) k=6", 4, 4, 5), ("D=inf h=5 (6,6) k=7", 5, 6, 6), ("D=inf h=6 (8,8) k=7", 6, 8, 8),
 ("D=inf h=7 (9,10) k=8", 7, 9, 10),
 ("AppB thr h=7 (13,14)", 7, 13, 14),
]
for lab, h, W, Wp in pts:
    N = p(h); a, b = min(P(W), N*N), min(P(Wp), N*N)
    th, T = math.log(a)/math.log(N), math.log(b)/math.log(N)
    ta, Ta = math.sqrt(W/h), math.sqrt(Wp/h)
    print(f"{lab:24s} N={N:>3}: counting dims ({a},{b}) -> finite (theta,T)=({th:.2f},{T:.2f}) th^2+T^2={th*th+T*T:.2f} th+T={th+T:.2f} | "
          f"asymptotic sqrt(W/h): ({ta:.2f},{Ta:.2f}) sum sq={ta*ta+Ta*Ta:.2f} sum={ta+Ta:.2f}")
print("hypothetical k = c*h:  theta^2+T^2 = 4-c; balanced 4rho = sqrt(2(4-c)); alpha' = log2(1/rho)")
for c in [0, 0.5, 1, 8/7, 1.5, 2]:
    s = math.sqrt(2*(4-c)); rho = s/4
    print(f"  c={c:.3f}: 4rho={s:.3f} rho={rho:.3f} alpha'={math.log2(1/rho):.3f}  inside circle: {c>0}")
