from math import log2
from bell_poly import formula
for k in [25,31,41,51,61]:
    d=(k-1)//2; A=list(range(0,d+1)); p,q=formula(k,A)
    print("k=%d d=%d p=%d q=%d p==q:%s exponent %.4f  (Hrubes asympt 1.614, App B 1.414, N^2/8 %.3f)"%(k,d,p,q,p==q,log2(max(p,q))/k,2-3/k),flush=True)
