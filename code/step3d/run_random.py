import sys, json, time, statistics
sys.path.insert(0,'/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/step3d')
from ns import *
n=int(sys.argv[1]); dens=float(sys.argv[2]); K=int(sys.argv[3]); dmax=int(sys.argv[4]) if len(sys.argv)>4 else n
m=round(dens*n); res=[]; seed=0
while len(res)<K:
    seed+=1; cl=rand3cnf(n,m,1000*n+seed)
    if not unsat(n,cl): continue
    P=[clause_poly(c) for c in cl]; t=time.time()
    r={q:min_deg(n,P,q,dmax) for q in (2,3,4)}
    r[6]=None if r[2] is None or r[3] is None else max(r[2],r[3])
    res.append(dict(seed=1000*n+seed,deg=r,t=round(time.time()-t,1)))
    print(n,res[-1],flush=True)
def med(q):
    v=[x['deg'][q] for x in res]; return statistics.median([x if x is not None else 99 for x in v])
summary={'n':n,'m':m,'density':dens,'median':{q:med(q) for q in (2,3,4,6)},'runs':res}
print('SUMMARY',json.dumps(summary['median']))
json.dump(summary,open(f'/tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/step3d/res_n{n}_d{dens}.json','w'))
