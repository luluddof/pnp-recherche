# check levH/hsup of dp_corner against brute force over boundary samples (no DP run)
import sys, os, numpy as np, math
for two in ('0','1'):
  for cs in ('2','4','8'):
    os.environ['TWO']=two; sys.argv=['x',cs]
    src=open('dp_corner.py').read().split('cache = {}')[0]
    g={}; exec(src,g)
    rng=np.random.default_rng(1); n=4000
    ea=rng.uniform(0.01,1,n); eb=ea*np.exp(rng.uniform(-4,4,n)); al=rng.uniform(0,0.3,n); be=rng.uniform(0,0.3,n)
    v,w=g['levH'](ea,eb,al,be)
    tau=np.linspace(0,g['LT'],200001); u,vv=g['bpt'](tau)
    bf=np.array([np.max(np.minimum(a*(u+x),b*(vv+y))) for a,b,x,y in zip(ea[:300],eb[:300],al[:300],be[:300])])
    Ew=g['hsup'](w*ea,(1-w)*eb); Ow=w*ea*al+(1-w)*eb*be
    print(f"TWO={two} c={cs}: max|v-bruteforce|/v={np.max(np.abs(v[:300]-bf)/bf):.2e}  max (Ew+Ow-v)/v={np.max((Ew+Ow-v)/v):.2e} min={np.min((Ew+Ow-v)/v):.2e}")
