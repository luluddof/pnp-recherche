# Bell-model optimum at N=8 via MILP (big-M, |kappa|<=M): min #errors s.t. #support<=q. Lower bounds only valid for |kappa|<=M.
import numpy as np, itertools, pulp, sys
k=3; n=64; M=float(sys.argv[2]) if len(sys.argv)>2 else 4.0
pts=np.array(list(itertools.product([0,1],repeat=2*k)))
beta=(-1)**((pts[:,:k]*pts[:,k:]).sum(1)%2)
CH=(-1)**((pts@pts.T)%2)
for q in [int(a) for a in sys.argv[1].split(',')]:
    prob=pulp.LpProblem("bell",pulp.LpMinimize)
    kap=[prob.add_variable("k%d"%s,lowBound=-M,upBound=M) for s in range(n)]
    z=[prob.add_variable("z%d"%s,cat="Binary") for s in range(n)]
    e=[prob.add_variable("e%d"%v,cat="Binary") for v in range(n)]
    prob+=pulp.lpSum(e)
    prob+=pulp.lpSum(z)<=q
    # symmetry breaking: translation invariance lets us assume 0 in S (char twist) -- valid: S -> S+l
    prob+=z[0]==1
    for s in range(n):
        prob+=kap[s]<=M*z[s]; prob+=kap[s]>=-M*z[s]
    BM=q*M+1
    for v in range(n):
        expr=pulp.lpSum(int(CH[s,v])*kap[s] for s in range(n))-int(beta[v])
        prob+=expr<=BM*e[v]; prob+=expr>=-BM*e[v]
    prob.solve(pulp.HiGHS(msg=True,timeLimit=int(sys.argv[3]) if len(sys.argv)>3 else 100,threads=1))
    val=pulp.value(prob.objective); bound=prob.solverModel if False else None
    S=[s for s in range(n) if z[s].value()>0.5]
    kv=np.array([kap[s].value() for s in range(n)])
    phi=kv@CH; perr=int((np.abs(phi-beta)>1e-6).sum())
    print("N=8 q=%d M=%g obj=%s verified_errors=%d S=%s kappa=%s"%(q,M,val,perr,S,np.round(kv[S],3)),flush=True)
