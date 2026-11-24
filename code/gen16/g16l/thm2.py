import itertools
from fractions import Fraction
def f(l, a, n):  # sum_x max(0, n - 2 a.x)
    return sum(max(0, n - 2*sum(a[c] for c in range(l) if x[c])) for x in itertools.product([0,1], repeat=l))
def comps(N, l):
    if l == 1: yield (N,); return
    for a in range(N+1):
        for r in comps(N-a, l-1): yield (a,)+r
bad = 0; tot = 0
for l in range(1, 6):
    for N in range(0, 13):
        bal = tuple(N//l + (1 if i < N % l else 0) for i in range(l))
        for d in range(0, N+1):
            n = N-d
            fb = f(l, bal, n); fm = f(l, [Fraction(N, l)]*l, n)
            # cyclic-average integer form:  l*f(m) >= sum_x max(0, l*n - 2N|x|)
            cyc = sum(max(0, l*n - 2*N*sum(x)) for x in itertools.product([0,1], repeat=l))
            assert cyc == l*fm
            for m in comps(N, l):
                if list(m) != sorted(m, reverse=True): continue
                v = f(l, m, n); tot += 1
                if v < fb or v < fm or l*v < cyc: bad += 1; print("VIOLATION", l, N, d, m, v, fb, fm)
print("Theorem 2 enumeration: cases", tot, "violations", bad)
