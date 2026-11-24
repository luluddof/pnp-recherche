# Weighted identity behind Theorem 2 (definitions only):
#   W = sum over (Y,X) [K<=Y, X<=E(G[Y]), |Y|=l+1, |X|=l, no deleted edge in X] of prod_{e in X} w_e * p^{cc(Y,X)}
#   claim: (W/p) mod p == f(u) = sum_{T in ST(l)} prod_{e in T} u_e mod p   (u_e in 0..p-1; u_e=0 <=> deleted; w_e=u_e)
# f(u) computed independently by the weighted matrix-tree theorem.
import itertools, random, sys
sys.path.insert(0, '.')
from identity_all import cc_of, tau
rnd = random.Random(7)
tests = 0; nz = 0
for trial in range(int(sys.argv[1])):
    n = rnd.randint(4, 8); m = rnd.randint(n - 1, min(n * (n - 1) // 2, 12))
    allE = list(itertools.combinations(range(n), 2)); rnd.shuffle(allE); E = sorted(allE[:m])
    p = rnd.choice([2, 3, 5, 7, 11])
    u = [rnd.randrange(p) for _ in E]
    K = rnd.randrange(1, 1 << n)
    for l in range(1, n):
        W = 0; f = 0
        for Yt in itertools.combinations(range(n), l + 1):
            Y = sum(1 << v for v in Yt)
            if Y & K != K: continue
            idx = [j for j, e in enumerate(E) if Y >> e[0] & 1 and Y >> e[1] & 1 and u[j] != 0]
            for Xi in itertools.combinations(idx, l):
                c = cc_of(Y, [E[j] for j in Xi], n)
                pr = 1
                for j in Xi: pr *= u[j]
                W += pr * p ** c
            f += tau(Yt, E, wt=u)
        assert W % p == 0 and (W // p) % p == f % p, (n, E, K, l, p, u)
        tests += 1; nz += (f % p != 0)
print(f"weighted identity: {tests} checks passed ({nz} with nonzero f(u) mod p)")
