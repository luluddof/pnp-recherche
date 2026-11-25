# Numerical sanity check of Theorem 2: sum_{k<D/2} f(th*s_k) <= zeta(2)/(4 th), s_k = m+1+2k+2v2(k!), D=2^m;
# popcount bound P(K) <= (K/2) log2 K; and C1 counterexample at odd p from Smith(M) (F ~ M (x) M').
import math
Z2 = math.pi**2/6
worst = 0
for m in range(1, 19):
    D = 2**m; s = []; v = 0
    for k in range(D//2):
        if k: v += (k & -k).bit_length() - 1 + 0  # placeholder replaced below
    s = []; vf = 0
    for k in range(D//2):
        if k:
            kk = k
            while kk % 2 == 0: kk //= 2; vf += 1
        s.append(m + 1 + 2*k + 2*vf)
    for e in range(-60, 30):
        th = 2.0**(e/3)
        G = sum(-math.log1p(-math.exp(-th*x)) for x in s if th*x < 50)
        worst = max(worst, G/(Z2/(4*th)))
print("max ratio G/(zeta2/(4th)) over m<=18, th grid:", worst)
P = 0; ok = True
for K in range(1, 1 << 20):
    P += bin(K-1).count("1")
    if P > (K/2)*math.log2(K) + 1e-9: ok = False; print("popcount bound fails", K); break
print("P(K) <= (K/2)log2 K for K < 2^20:", ok)
M = [0, 4, 6, 9, 10, 10, 12, 13, 13, 14, 15, 17, 19, 20, 20, 22]    # Smith(M_{3,0}) from m_profile_out.txt
F = sorted(u + v for u in M for v in M)
modes = [4, 6, 10, 12] * 2
def monos(ws, cap):
    res = [0]
    for w in ws:
        res = res + [r + e*w for r in res for e in range(1, cap) if r + e*w < cap]
    return res
for p in range(2, 14):
    print(f"p={p}: R_p(F)={sum(1 for x in F if x < p)}  Smith-monomials={len(monos(modes, p))}")
