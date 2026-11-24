# (i) naive contraction: Y* 2/3-separator of F*, |Y*|<=n*/9, but pi(Y*) NOT a 2/3-separator of F (cardinality), while d-weight balance holds.
import random, itertools
src = open('transfer_test.py').read().replace("\nmain()\n", "\n"); exec(src)
exec(open('transfer_test2.py').read().split("rnd = random.Random(7)")[0].split("exec(src)")[1])
def compsF(n, clauses, Y):
    adj = [set() for _ in range(n)]
    for e in clauses:
        for a in e:
            for b in e:
                if a != b: adj[a].add(b)
    return components(n, adj, Y)
rnd = random.Random(3); found = 0; tested = 0; wviol = 0
for it in range(300):
    # dense core on few variables + many low-degree variables in one sparse chain
    core = rnd.randint(4, 5); tail = rnd.randint(9, 14); n = core + tail
    clauses = [tuple(rnd.sample(range(core), 3)) for _ in range(rnd.randint(8, 12))]
    t = list(range(core, n))
    for i in range(0, tail-2, 2): clauses.append((t[i], t[i+1], t[i+2]))   # chain: tail connected, degree <= 2
    d, pi, fib, adj, ns, n0 = expand(n, clauses, rnd, keep=0)
    for k in range(0, min(ns//9, 3)+1):
        for Ys in itertools.combinations(range(3*8), k) if k else [()]:
            cs = components(ns, adj, Ys)
            if any(3*len(c) > 2*ns for c in cs): continue
            tested += 1
            Y = set(pi[x] for x in Ys)
            cF = compsF(n, clauses, Y)
            for K in cF:
                if sum(d[v] for v in K) * 3 > 2*ns: wviol += 1
            if any(3*len(K) > 2*n for K in cF):
                found += 1
                if found == 1:
                    print("example: n=%d m=%d n*=%d clauses=%s" % (n, len(clauses), ns, clauses))
                    print("  Y* =", Ys, "-> Y =", sorted(Y), "; component sizes of Gr(F*)-Y*:", sorted(map(len, cs)), "(2n*/3=%.1f)" % (2*ns/3))
                    print("  component sizes of Gr(F)-Y:", sorted(map(len, cF)), "(2n/3=%.1f)" % (2*n/3))
print("tested separators:", tested, " pi(Y*) not a 2/3-separator of F:", found, " occurrence-weight balance violated:", wviol)
