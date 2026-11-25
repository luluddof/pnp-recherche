from math import comb, log2
from bell_poly import formula
from bell_poly_big import best_window
import sys
# balanced point: min over d of max(p,q); compare N^2/8
for k in range(3,25):
    best=None
    for d in range(0,k+1):
        A=best_window(k,d); p,q=formula(k,A)
        m=max(p,q)
        if best is None or m<best[0]: best=(m,d,p,q)
    m,d,p,q=best; N=2**k
    print("k=%2d N=2^%d bal max(p,q)=%d (p=%d,q=%d,d=%d) exponent %.4f ; N^2/8 exponent %.4f"%(k,k,m,p,q,d,log2(m)/k,2-3/k),flush=True)
