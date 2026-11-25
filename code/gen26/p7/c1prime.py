import sys; sys.path.insert(0,'.')
from flat7 import *
from smith7 import predicted
def sqfree(modes,p):
    cnt=[0]*p; cnt[0]=1
    for s in modes:
        for w in range(p-1,s-1,-1): cnt[w]+=cnt[w-s]
    return sum(cnt)
for a,i,t in [(3,0,1),(3,1,2),(3,2,4),(2,0,1),(2,1,2)]:
    F=flatten(a,i,t,1); sv=smith_np(F); modes=predicted(a,i,1)*2
    ok=True; rows=[]
    for p in range(1,40):
        R=sum(1 for s in sv if s<p); lo=sqfree(modes,p); hi=monomial_count(modes,p)
        if not (lo<=R<=hi): ok=False; rows.append((p,lo,R,hi))
    print(a,i,t,"sqfree<=R_p<=monomials for p<40:",ok,rows[:5])
