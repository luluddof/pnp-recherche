# formula-only evaluation of the Bell polynomial family for large k; optimal window by mass; exact integers
import sys
from fractions import Fraction
from math import comb, log2
from bell_poly import formula
def best_window(k,d):
    best=None
    for a0 in range(0,k-d+1):
        mass=sum(comb(k,m)*3**(k-m) for m in range(a0,a0+d+1))
        if best is None or mass>best[0]: best=(mass,a0)
    return list(range(best[1],best[1]+d+1))
for k in ([int(a) for a in sys.argv[1:]] if __name__=="__main__" else []):
    rows=[]
    for d in range(0,k+1):
        A=best_window(k,d); p,q=formula(k,A)
        if p==0: break
        th=log2(p)/k; T=log2(q)/k
        rows.append((d,A[0],q,p,T,th))
        print("k=%d d=%d a0=%d  eta=T=%.4f theta=%.4f delta=%.4f  eta^2/4=%.4f"%(k,d,A[0],T,th,2-th,T*T/4),flush=True)
