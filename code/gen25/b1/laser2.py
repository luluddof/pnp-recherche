import itertools, numpy as np, sys
sys.argv = ['x', '0']
exec(open('laser.py').read().split("k = int(sys.argv[1])")[0])
def run2(name, X, k):
    n = len(X); P = pieces_by_j(X, k)
    js = [j for j in range(k+1) if P[j].any()]
    best = None
    for r in range(len(js)+1):
        for J in itertools.combinations(js, r):
            Rm = P[list(J)].sum(0) if J else np.zeros_like(P[0]); Cm = P[[j for j in js if j not in J]].sum(0)
            a, b = rk(Rm), osr(Cm)
            if best is None or max(a, b) < best[0]: best = (max(a, b), a, b, J)
    print(name, "N'=", n, 'best over all J-subsets:', best, 'exp %.3f' % (np.log(best[0]) / np.log(n)), flush=True)
def gen_even(m):  # even weight code length m
    return np.array([[1 if (j == i or j == i+1) else 0 for j in range(m)] for i in range(m-1)])
todo = [('full3', np.eye(3, dtype=int), 3), ('full4', np.eye(4, dtype=int), 4), ('even5', gen_even(5), 5), ('even6', gen_even(6), 6),
        ('ham743', np.array([[1,0,0,0,0,1,1],[0,1,0,0,1,0,1],[0,0,1,0,1,1,0],[0,0,0,1,1,1,1]]), 7),
        ('punct_even', np.array([[1,0,0,0,0,1],[0,1,0,0,0,1],[0,0,1,0,0,1],[0,0,0,1,0,1],[0,0,0,0,1,1]]), 6)]
for nm, G, kk in todo:
    run2(nm, code_from_gen(G), kk)
