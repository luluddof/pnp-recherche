# kappa = u * max_i Lambda_i / M for (a) paper taper h=ceil(3M/u), budget 8hD -> 24 (stated), loss 6h -> 18
# (b) linear taper with minimal h = M/(u-1): exact max_i Lambda_i ; (c) LP-optimal taper.  Units M=1, L -> inf.
import numpy as np
from scipy.optimize import linprog
def lam_profile(q,a):
    u=len(q); Q=lambda k: q[k] if k<u else 0.0
    return [sum(j/2**j*(Q(i)-Q(i+j)) for j in range(1,a-i+1)) for i in range(u)]
for u in (10,30,100,300):
    a=u+20
    h=1.0/(u-1); q=[h*(u-i) for i in range(u)]
    # check sum_i 2^-i q_i >= 2 (deg bound M(2L-1)/L -> 2)
    s=sum(2**-i*q[i] for i in range(u))
    lin=max(lam_profile(q,a))*u
    # LP: min Lam s.t. Lambda_i <= Lam, sum 2^-i q_i >= 2, q>=0
    A=[];b=[]
    for i in range(u):
        row=np.zeros(u+1)
        for j in range(1,a-i+1):
            w=j/2**j; row[i]+=w
            if i+j<u: row[i+j]-=w
        row[u]=-1; A.append(row); b.append(0)
    row=np.zeros(u+1); row[:u]=[-2.0**-i for i in range(u)]; A.append(row); b.append(-2)
    c=np.zeros(u+1); c[u]=1
    r=linprog(c,A_ub=A,b_ub=b,bounds=[(0,None)]*(u+1),method='highs')
    print(f"u={u:4d}  sum2^-i q={s:.4f}  linear-min-h kappa={lin:.4f}   LP kappa={r.x[u]*u:.4f}   ratio sqrt={np.sqrt(r.x[u]*u/lin):.4f}")
