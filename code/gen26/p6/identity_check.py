# Check the decomposition in step (4) of Theorem 1 on random subsets:
# V(R)+V(C)-sum_{RxC} min(i,v2(x-y)) == Q(e) - i k + sum_rho [g(R_rho)+g(C_rho)]
import random, itertools
def v2(n): n = abs(n); return (n & -n).bit_length() - 1
random.seed(1); bad = 0; tot = 0
for a in range(2, 8):
    L = 1 << a
    for i in range(a):
        B = [x for x in range(L) if (x >> i) & 1 == 0]
        for _ in range(200):
            k = random.randint(1, len(B)); R = random.sample(B, k); C = random.sample(B, k)
            mn = lambda x, y: i if x == y else min(i, v2(x - y))
            lhs = sum(v2(x - y) for x, y in itertools.combinations(R, 2)) + sum(v2(x - y) for x, y in itertools.combinations(C, 2)) \
                  - sum(mn(x, y) for x in R for y in C)
            res = range(1 << i); m = lambda r, s: i if r == s else v2(r - s)
            e = {r: sum(1 for x in R if x % (1 << i) == r) - sum(1 for y in C if y % (1 << i) == r) for r in res}
            Q = sum(m(r, s) * e[r] * e[s] for r in res for s in res) / 2
            def g(X):
                xs = [x >> (i + 1) for x in X]
                return len(xs) * (len(xs) - 1) // 2 + sum(v2(p - q) for p, q in itertools.combinations(xs, 2))
            rhs = Q - i * k + sum(g([x for x in R if x % (1 << i) == r]) + g([y for y in C if y % (1 << i) == r]) for r in res)
            tot += 1; bad += (lhs != rhs)
print("identity checks:", tot, "failures:", bad)
