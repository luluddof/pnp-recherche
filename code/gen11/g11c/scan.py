import sys, numpy as np, itertools
from feas import Prob
l=int(sys.argv[1]); d=int(sys.argv[2]); starts=int(sys.argv[3])
profs=[tuple(int(c) for c in a.split(',')) for a in sys.argv[4:]]
tau=[1.0,1.7,2.9,4.3,5.9,7.1][:l]
for pr in profs:
    P=Prob(l,d,tau,pr)
    cost,p=P.solve(starts=starts)
    ok = cost<1e-20
    extra=""
    if ok: extra=str(P.check(P.J_of(p)))
    print(pr,"total",sum(pr),"params",P.np_,"res",P.nres,"best cost %.2e"%cost,"max|p| %.1e"%np.abs(p).max(),"FEASIBLE" if ok else "not found",extra,flush=True)
