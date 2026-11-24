# asymptotic exponent of the Kronecker-graded hybrid (App-B factor x Clifford/JH factor), certified upper bounds:
# piece (rr, rc) exponents are gamma*(factor 1) + (1-gamma)*(factor 2); cross ranks of App-B thresholds and of JH parts = 2 (full).
import math, itertools
from math import sqrt
from curves import H, al
def circle(th): return sqrt(max(0, 4 - th*th))
JH = [(1 + 2*H(b), 1 + 2*H(al - b)) for b in [al*i/200 for i in range(201)]] + [(2, 0), (0, 2)]
AB = [(th, circle(th)) for th in [2*i/400 for i in range(401)]]
best = {}
for g in [i/100 for i in range(0, 101)]:
    for (t1, T1) in AB:
        for (t2, T2) in JH:
            RR = g*t1 + (1-g)*t2                      # R1 x T1'  -> R
            CC = g*T1 + (1-g)*T2                      # C1 x T2'  -> C
            m1 = (g*t1 + 2*(1-g), 2*g + (1-g)*T2)     # R1 x T2'
            m2 = (2*g + (1-g)*t2, g*T1 + 2*(1-g))     # C1 x T1'
            for o1, o2 in itertools.product((0, 1), repeat=2):
                th = max([RR] + ([m1[0]] if o1 == 0 else []) + ([m2[0]] if o2 == 0 else []))
                T = max([CC] + ([m1[1]] if o1 == 1 else []) + ([m2[1]] if o2 == 1 else []))
                key = round(th, 2)
                if T < best.get(key, (9,))[0]: best[key] = (T, g, (t1, T1), (t2, T2), (o1, o2))
print(" theta   T_hybrid  T_circle   (gamma, AppB pt, JH pt)")
for th in [1.0, 1.2, 1.3, 1.4, 1.41, 1.42, 1.5, 1.6, 1.7, 1.8, 1.9]:
    cands = [v for k, v in best.items() if k <= th + 1e-9]
    T, g, a, b, o = min(cands)
    print(f" {th:.2f}   {T:.4f}   {circle(th):.4f}   g={g:.2f} AB=({a[0]:.2f},{a[1]:.2f}) JH=({b[0]:.2f},{b[1]:.2f}) opt={o}")
