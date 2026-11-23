from math import log, sqrt
H=20000
# partition numbers via Euler pentagonal recurrence
p=[0]*(4*H+1); p[0]=1
for n in range(1,4*H+1):
    s=0;k=1
    while True:
        g1=k*(3*k-1)//2
        if g1>n: break
        sg = 1 if k%2 else -1
        s+=sg*p[n-g1]
        g2=k*(3*k+1)//2
        if g2<=n: s+=sg*p[n-g2]
        k+=1
    p[n]=s
pre=[0]
for v in p: pre.append(pre[-1]+v)
def S(lo,hi):
    lo=max(lo,0)
    return 0 if hi<lo else pre[hi+1]-pre[lo]
for h in [5,10,20,50,100,200,500,1000,2000,5000,10000,20000]:
    if 4*h>4*H: break
    B=sum(min(S(s,s+2*h),S(-s,2*h-s)) for s in range(-2*h,2*h+1))
    N=p[h]
    print(f"h={h:6d} log2 p(h)={log(N,2):10.1f}  log_p(h) B = {log(B)/log(N):.4f}   B<N^2: {B<N*N}  (log B - sqrt2 log N)/log2(h) = {(log(B)-sqrt(2)*log(N))/log(h):.3f}")
