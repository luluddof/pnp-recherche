"""G4b d-box DP (schemes.dbox_floor_best) with the corrected floor rho_d >= 1/d instead of 2^{1-d}.
E(x) = min_{y: d = x//y >= 2} rho(d) E(y) + add 2^{-y}, x = log2 K."""
import math
def best(X, rho, x0=6, add=4.0):
    E = {}
    for x in range(1, X + 1):
        if x <= x0: E[x] = 1.0; continue
        b = E[x - 1]
        for y in range(x0, x):
            d = x // y
            if d >= 2: b = min(b, rho(d) * E[y] + add * 2.0 ** (-y))
        E[x] = min(b, 1.0)
    return E[X]
for X in [64, 256, 1024, 4096]:
    a = best(X, lambda d: 2.0 ** (1 - d)); b = best(X, lambda d: 1.0 / d)
    c = best(X, lambda d: 0.5 if d == 2 else 1.0)
    print(f"X={X:5d}  log2(1/e): floor 2^(1-d) {-math.log2(a):6.2f} | floor 1/d {-math.log2(b):6.2f}"
          f" | d=2 only rho=1/2 {-math.log2(c):6.2f} | log2(X/6)={math.log2(X/6):5.2f} sqrt(X)={math.sqrt(X):5.1f}")
