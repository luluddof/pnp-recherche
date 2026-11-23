"""Threshold-density random instances; counts correct decisions and false accepts per relaxation."""
import random, sys, json
from relax import *
dens = {'2-SAT':1.0,'XORSAT':0.9,'Horn':1.2,'3-SAT':4.3,'1-in-3':0.65,'NAE-3':2.1}
algs = {'BLP': blp, 'AIP_Z': aip, 'AIP_GF2': aip_gf2, 'BLP+AIP': blp_aip}
rng = random.Random(int(sys.argv[2]) if len(sys.argv)>2 else 7)
trials = int(sys.argv[1])
res = {}
for name,T in TEMPLATES.items():
    st = {a:[0,0] for a in algs}; nsat=0; ex={}
    for _ in range(trials):
        n = rng.randint(4,9)
        m = max(2,int(round(dens[name]*n*rng.uniform(0.7,1.4))))
        cons = rand_instance(T,n,m,rng)
        tr = brute(n,cons); nsat+=tr
        for a,f in algs.items():
            r=f(n,cons)
            if r==tr: st[a][0]+=1
            else:
                st[a][1]+=1
                if a not in ex: ex[a]=(n,[(sc,sorted(R)) for sc,R in cons])
    res[name]=dict(trials=trials,sat=nsat,stats=st)
    print(name,f'sat={nsat}/{trials}',' '.join(f'{a}:{c}/{trials}(FA={fa})' for a,(c,fa) in st.items()),flush=True)
    for a,e in ex.items():
        if a=='BLP+AIP' or (name=='1-in-3' and a=='AIP_Z'): print('  example false-accept',a,e)
json.dump(res,open('search_results.json','w'),indent=1)
