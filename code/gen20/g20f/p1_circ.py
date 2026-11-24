import numpy as np, random, itertools, networkx as nx, time
from p1_lib import Aexact
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
random.seed(5); rng=np.random.default_rng(5)
# (1d) brute force on small graphs, arbitrary (possibly empty / unequal) classes
viol=0; tests=0; eq=0
for trial in range(300):
    u=random.randint(1,4); n=random.randint(2,9); M=random.randint(1,9)
    wires=[(random.randrange(n),random.randrange(n)) for _ in range(M)]
    cls=np.array([random.randrange(2*u) for _ in range(n)])
    parts=[np.isin(cls,[(i+j)%(2*u) for j in range(u)]).astype(int) for i in range(u)]
    A,_,_=Aexact(wires,n,parts)
    cutK=[sum((cls[h]==c)!=(cls[t]==c) for h,t in wires) for c in range(2*u)]
    bnd=0.5*max(cutK[j]+cutK[(j+u)%(2*u)] for j in range(u))
    tests+=1; viol+= A>bnd+1e-9; eq+= A==bnd
print('(1d) Prop C brute force: tests',tests,'violations',viol,'equalities',eq)
# random 4-regular n=48, equal classes, exact MILP A''
def circ_exact(G,u,seed):
    n=G.number_of_nodes(); wires=list(G.edges()); M=len(wires)
    perm=np.random.default_rng(seed).permutation(n); cls=np.empty(n,int); cls[perm]=np.arange(n)%(2*u)
    D=np.zeros((u,M))
    for i in range(u):
        H=set((i+j)%(2*u) for j in range(u))
        for e,(h,t) in enumerate(wires): D[i,e]=(cls[h] in H)-(cls[t] in H)
    Mi=(D!=0).sum(1)
    best_lp=-1; best_int=-1
    for s in itertools.product((-1,1),repeat=u-1):
        s=np.array((1,)+s)
        # vars: x+ (M), x- (M), t ; maximize t ; s_i * D_i (x+ - x-) >= t
        A=np.hstack([-(s[:,None]*D), (s[:,None]*D), np.ones((u,1))])
        c=np.zeros(2*M+1); c[-1]=-1
        lp=linprog(c,A_ub=A,b_ub=np.zeros(u),A_eq=None,bounds=[(0,1)]*(2*M)+[(None,None)],method='highs')
        best_lp=max(best_lp,-lp.fun)
        if -lp.fun<=best_int+1e-9: continue
        cons=[LinearConstraint(A,-np.inf,0), LinearConstraint(np.hstack([np.eye(M),np.eye(M),np.zeros((M,1))]),0,1)]
        res=milp(c,constraints=cons,integrality=np.r_[np.ones(2*M),0],bounds=Bounds(np.r_[np.zeros(2*M),-1e3],np.r_[np.ones(2*M),1e3]),options=dict(time_limit=20))
        if res.x is not None: best_int=max(best_int,round(-res.fun))
    return M,Mi,best_lp,best_int
G=nx.random_regular_graph(4,48,seed=11)
for u in (3,4,6,8):
    t0=time.time(); M,Mi,lp,ip=circ_exact(G,u,u)
    print(f'u={u}: M={M} M/u={M/u:.2f} min M_i={Mi.min()} max M_i={Mi.max()}  A_LP={lp:.3f} (A_LP/M={lp/M:.4f})  A_int={ip}  ({time.time()-t0:.0f}s)',flush=True)
