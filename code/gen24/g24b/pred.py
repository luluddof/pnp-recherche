# predicted flattening rank = sum over label-weight sequences and charge blocks q of min(rows_q, cols_q)
import itertools
from functools import lru_cache
@lru_cache(None)
def part(n):
    if n<0: return 0
    p=[1]+[0]*n
    for k in range(1,n+1):
        for m in range(k,n+1): p[m]+=p[m-k]
    return p[n]
def dV(h,s): return part(h-s*s+s)          # A1, xi=alpha/2: sector s weight s^2-s
dM=dV                                       # M = V'
CH=range(-6,8)
def r1(ha,e,q): return sum(dV(ha,s)*dM(e,q+s) for s in CH)     # rows (a in s, m' in q+s)
def c1(hb,e,q): return sum(dV(hb,t)*dM(e,q-t) for t in CH)     # cols (b in t, m in q-t)
def comps(tot,n):
    if n==1: yield (tot,); return
    for k in range(tot+1):
        for rest in comps(tot-k,n-1): yield (k,)+rest
def pred(n,E0,ha,hb):
    tot=0
    for es in comps(E0,n):
        for q in itertools.product(CH,repeat=n):
            R=sum(__import__('math').prod(r1(a,e,qq) for a,e,qq in zip(has,es,q)) for has in comps(ha,n))
            C=sum(__import__('math').prod(c1(b,e,qq) for b,e,qq in zip(hbs,es,q)) for hbs in comps(hb,n))
            tot+=min(R,C)
    return tot
if __name__=="__main__":
    for (e,h) in [(0,0),(2,2),(4,4),(6,6),(4,0),(6,2)]: print("1copy",e,h,pred(1,e,h,h))
    for (e,ha,hb) in [(0,0,0),(1,1,1),(2,0,0),(0,2,2),(2,1,1),(2,2,2),(1,2,2),(2,1,3),(3,1,1),(3,2,2),(2,3,3),(3,3,3),(4,2,2),(4,1,1),(4,3,3),(4,4,4)]:
        print("2copy",e,ha,hb,pred(2,e,ha,hb))
