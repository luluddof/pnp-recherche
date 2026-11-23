from math import log, sqrt
H=6000
p=[1]+[0]*(4*H)
for k in range(1,4*H+1):          # coin-change partition count, independent of pentagonal recurrence
    for n in range(k,4*H+1): p[n]+=p[n-k]
pre=[0]
for v in p: pre.append(pre[-1]+v)
S=lambda lo,hi: 0 if hi<max(lo,0) else pre[hi+1]-pre[max(lo,0)]
first2=first162=None
for h in list(range(2,60))+list(range(60,1500,10))+[2000,3000,4000,6000]:
    B=sum(min(S(s,s+2*h),S(-s,2*h-s)) for s in range(-2*h,2*h+1)); N=p[h]
    e=log(B)/log(N)
    if first2 is None and B<N*N: first2=(h,log(N,2))
    if first162 is None and e<1.62: first162=(h,log(N,2))
    if h in (10,50,100,500,1000,2000,6000): print(h, round(log(N,2),1), round(e,4), round(log(p[2*h])/log(N),4))
print("first h with B<N^2:",first2,"first with exponent<1.62:",first162)
