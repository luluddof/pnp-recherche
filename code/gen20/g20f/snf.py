import json, itertools, sys, sympy as sp
from sympy.matrices.normalforms import invariant_factors
for fn in sys.argv[1:]:
    d=json.load(open(fn)); N=d['N']; G=[[[int(v) for v in r] for r in g] for g in d['G']]; Gp=[[[int(v) for v in r] for r in g] for g in d['Gp']]
    Rr=range(N)
    R=sp.Matrix([[int(x==v and y==u)-sum(G[l][x][u]*Gp[l][y][v] for l in range(len(G))) for u in Rr for v in Rr] for x in Rr for y in Rr])
    f=[a for a in invariant_factors(R,domain=sp.ZZ) if a!=0]
    print(fn.split('/')[-1],'nonzero invariant factors of R:',f)
