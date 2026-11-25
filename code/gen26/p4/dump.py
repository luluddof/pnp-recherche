# write GF(2) MC instance (branch r) as pure CNF (Tseitin xor chains) for kissat/cadical + DRAT
import sys, os, mc
N, p, q, r = map(int, sys.argv[1:5]); fn = sys.argv[5]
E, e, c, U, V = mc.build(N, p, q, r, order=bool(int(os.environ.get('ORD', '0'))))
nv, cl = mc.xor_cnf(E)
with open(fn, 'w') as f:
    f.write(f"p cnf {nv} {len(cl)}\n")
    for x in cl: f.write(' '.join(map(str, x)) + ' 0\n')
print(fn, nv, len(cl))
