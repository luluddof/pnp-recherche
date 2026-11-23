# E1: verify identity (exact + truncated) vs brute force and vs closed form (*)
import random
from math import comb
from ie import *
rng = random.Random(1)
bad = 0; checks = 0
for trial in range(40):
    n = rng.randint(4, 9); m = rng.randint(3, 14)
    F = rand3sat(n, m, rng)
    a = hist(n, F)
    for t in range(0, m + 1):
        Td = T_direct(n, F, t); Th = T_from_hist(a, t)
        closed = int(a[0]) + (-1) ** t * sum(int(a[r]) * comb(r - 1, t) for r in range(t + 1, m + 1))
        checks += 1
        if not (Td == Th == closed): bad += 1
    assert T_direct(n, F, m) == int(a[0])
print(f"E1: {checks} (formula,t) checks, mismatches = {bad}; exact IE == brute-force #SAT on all 40 formulas")
# Bonferroni error profile example
rng = random.Random(7); n = 12; F = rand3sat(n, int(4.26 * n), rng); a = hist(n, F)
print("example n=12 m=51 #SAT=", a[0], " rmax=", max(i for i in range(len(a)) if a[i]))
for t in range(0, 26, 2):
    print(f"  t={t:2d} T_t={T_from_hist(a,t):>22d}  T_t+1={T_from_hist(a,t+1):>22d}")
