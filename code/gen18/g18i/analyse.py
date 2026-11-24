import json, glob, math, statistics as st
from collections import defaultdict
rows=[json.loads(l) for f in sorted(glob.glob('res_*.jsonl')) for l in open(f)]
g=defaultdict(list)
for r in rows:
    if r.get('full_unsat') and 'mus_m' in r: g[(r['n'],r['D'])].append(r)
med=lambda xs: st.median(xs)
print("n    D    #  | pre_n/n  sqrt(4.27/D) | peel_n/n peel_m/peel_n | mus_n/n mus_m/n mus_m/mus_n  LB=1/(e^2 D) | conflicts full/pre/peel/mus | mus_done")
T={}
for (n,D) in sorted(g):
    R=g[(n,D)]; f=lambda k: med([r[k] for r in R])
    T[(n,D)]=dict(pre=f('pre_n')/n, peel=f('peel_n')/n, musn=f('mus_n')/n, musm=f('mus_m')/n, prem=f('pre_m')/n, peelm=f('peel_m')/n)
    print(f"{n:4d} {D:5.0f} {len(R):2d} | {f('pre_n')/n:.3f}   {math.sqrt(4.27/D):.3f}      | {f('peel_n')/n:.3f}   {med([r['peel_m']/r['peel_n'] for r in R]):.2f}     | {f('mus_n')/n:.3f}  {f('mus_m')/n:.3f}  {med([r['mus_m']/r['mus_n'] for r in R]):.2f}    {1/(math.e**2*D):.4f} | {f('full_conf'):.0f}/{f('pre_conf'):.0f}/{f('peel_conf'):.0f}/{f('mus_conf'):.0f} | {sum(r['mus_done'] for r in R)}/{len(R)}")
def fit(xs,ys):
    lx=[math.log(x) for x in xs]; ly=[math.log(y) for y in ys]; mx=sum(lx)/len(lx); my=sum(ly)/len(ly)
    b=sum((a-mx)*(c-my) for a,c in zip(lx,ly))/sum((a-mx)**2 for a in lx); return b, math.exp(my-b*mx)
print("\nlog-log slopes in D (quantity ~ C * D^slope), per n:")
for n in sorted({k[0] for k in T}):
    Ds=sorted(D for (nn,D) in T if nn==n)
    if len(Ds)<3: continue
    print(n, 'D in', Ds[0], '..', Ds[-1], {k: tuple(round(v,3) for v in fit(Ds,[T[(n,D)][k] for D in Ds])) for k in ['pre','peel','musn','musm','prem','peelm']})
