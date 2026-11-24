# Brute-force check of the combinatorial identity behind Lemmas 1 and 3 (not of the CNF encoding):
#   S(p) = sum over (Y,X): K<=Y, X<=E(G[Y]), |Y|=l+1, |X|=l of p^{cc(Y,X)}  satisfies  S/p = #SteinerTrees(l) mod p
# and the analogous statement for l-cycles (2-regular (Y,X) with |Y| = l).
import itertools, random
def cc(Y, X):
    par = {v: v for v in Y}
    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for u, v in X: par[f(u)] = f(v)
    return len({f(v) for v in Y})
def run(n, m, seed):
    rnd = random.Random(seed)
    allE = list(itertools.combinations(range(n), 2)); rnd.shuffle(allE); E = allE[:m]
    K = set(rnd.sample(range(n), 2))
    for p in (2, 3, 5):
        for l in range(1, n):
            S = 0; trees = 0
            for Y in itertools.combinations(range(n), l + 1):
                Ys = set(Y)
                if not K <= Ys: continue
                EY = [e for e in E if e[0] in Ys and e[1] in Ys]
                for X in itertools.combinations(EY, l):
                    c = cc(Y, X); S += p ** c; trees += (c == 1)
            assert S % p == 0 and (S // p) % p == trees % p, (n, m, p, l)
        for l in range(3, n + 1):
            S = 0; cyc = 0
            for Y in itertools.combinations(range(n), l):
                Ys = set(Y); EY = [e for e in E if e[0] in Ys and e[1] in Ys]
                for X in itertools.combinations(EY, l):
                    deg = {v: 0 for v in Y}
                    for u, v in X: deg[u] += 1; deg[v] += 1
                    if all(d == 2 for d in deg.values()):
                        c = cc(Y, X); S += p ** c; cyc += (c == 1)
            assert S % p == 0 and (S // p) % p == cyc % p, ('cyc', n, m, p, l)
for s in range(6): run(7, 11, s)
print("identity verified on 6 random graphs (n=7, m=11), p in {2,3,5}, all l")
