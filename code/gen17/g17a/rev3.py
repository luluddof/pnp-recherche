#!/usr/bin/env python3
"""G17a part 3: product witnesses for link (i), and the mode expansion (D.12)."""
from fractions import Fraction as Fr
from rev1_lib import *
def witness(a, u, M, levels=None):
    """F_m(T) = prod_{q in S\\0} ((T-q)/(-q))^{m_q}, m_q = h'(u - v2(q))_+ with h' maximal s.t. deg <= M(2L-1)."""
    L, S = grid(a); al = 2 ** (a + 1); N = M * (2 * L - 1)
    w = sum(Fr(L, 2 ** i) * (u - i) for i in range(u)); h = int(N // w)
    m = {q: h * max(0, u - v2(q)) for q in S if q != 0}
    F = [Fr(1)]
    for q, mq in m.items(): F = pmul(F, ppow([Fr(1), Fr(-1, q)], mq))
    assert deg(F) <= N and peval(F, 0) == 1
    out = {}
    for t in S:
        g = gauss(pcompose_affine(F, t, al)); g2 = gauss(pcompose_affine(F, t, 2 ** a))
        if t == 0: assert g == 0; continue
        i = v2(t); D = Fr(L, 2 ** i)
        prof = [h * max(0, u - k) for k in range(a + 2)]
        pred = D * sum(Fr(j, 2 ** j) * (prof[i] - prof[i + j]) for j in range(1, a - i + 1)) + m[t] * (a + 1 - i)
        assert g == pred, (t, g, pred)
        assert g2 == pred - m[t]          # still integral on the larger disc
        out[i] = (int(g), float(g / D), float(Fr(g) / D / Fr(M, u)))
    return h, out
print("W1 product witness F_m (deg<=M(2L-1), F_m(0)=1, integral on all discs, ALL Taylor coefficients at node t divisible by 2^g_t):")
print("   level i -> (g_t, g_t/D_t, g_t/(D_t*M/u));  any recovery rule from Taylor data at nodes of levels <u needs p_t > g_t at some used node")
for a, u, M in [(2, 1, 3), (2, 2, 4), (3, 2, 4), (3, 3, 6), (4, 2, 6), (4, 3, 6), (4, 4, 8)]:
    h, out = witness(a, u, M)
    print(f"   a={a} u={u} M={M} h'={h}:", {i: out[i] for i in sorted(out) if i < u})
# asymptotic value of the witness bound: min_i  sum_j j/2^j min(j,u-i)  (times h' ~ M/u)
for a, u in [(8, 4), (16, 8), (40, 20), (60, 30), (30, 30)]:
    c = [sum(j / 2.0 ** j * min(j, u - i) for j in range(1, a - i + 1)) for i in range(u)]
    print(f"   a={a} u={u}: min_i c_i={min(c):.4f} (level {c.index(min(c))}), c_0={c[0]:.4f};  bound p_t/D_t > ~{min(c):.3f} M/u  (Theorem B's formula bound: 0.5 M/u, LP ~5.9 M/u, paper uses 24 M/u)")

# (D.12): geometric/mode expansion, exact valuations of the mode coefficients
def d12(a, i, P=40):
    L, S = grid(a); al = 2 ** (a + 1); ok = True; minextra = None
    Cp0 = Fr(1)
    for q in S:
        if q != 0: Cp0 *= -q
    def C(T):
        p = Fr(1)
        for q in S: p *= (T - q)
        return p
    for t in [s for s in S if s != 0 and v2(s) == i]:
        for zeta in (0, 1, 3, -5):
            T = t + al * zeta
            for x in range(L):
                for y in range(L):
                    if ((x >> i) & 1) != ((y >> i) & 1): continue      # crossing pairs only
                    y0 = y % 2 ** (i + 1); j = (y - y0) // 2 ** (i + 1)
                    v = v2(t - y0 + x); ok &= (v <= i)
                    exact = C(T) / (Cp0 * (T - y + x))
                    ser = sum(C(T) / (Cp0 * (T - y0 + x)) * (Fr(2 ** (i + 1), T - y0 + x)) ** b * j ** b for b in range(P))
                    ok &= v2(exact - ser) >= P if exact != ser else True
                    if zeta != 0:
                        for b in range(4):
                            val = v2(C(T) / (Cp0 * (T - y0 + x)) * (Fr(2 ** (i + 1), T - y0 + x)) ** b)
                            ok &= (val >= (a + 1 - v) + b * (i + 1 - v))
                            ex = val - (b + 1); minextra = ex if minextra is None else min(minextra, ex)
    return ok, minextra
for a, i in [(2, 0), (2, 1), (3, 0), (3, 1), (3, 2), (4, 1)]:
    ok, me = d12(a, i)
    print(f"W2 (D.12) a={a} level i={i}: expansion correct mod 2^40 and coefficient valuation >= (a+1-v)+b(i+1-v): {ok}; "
          f"min over rows/tiles/b of (valuation - w_Z) = {me}  (paper only uses >= 0; a-i = {a-i})")
