import mpmath as mp
mp.mp.dps=60
H=lambda x: -x*mp.log(x,2)-(1-x)*mp.log(1-x,2)
print('2H(1/4)=',2*H(mp.mpf(1)/4))
for e in [2,3,5,10,20,40,80]:
    t=mp.mpf(10)**(-e); eta=2*H(t); de=2-2*H(mp.mpf(1)/2-t)
    print(f'tau=1e-{e}: eta={mp.nstr(eta,5)} delta={mp.nstr(de,5)}  delta*ln^2(1/eta)/eta^2={mp.nstr(de*mp.log(1/eta)**2/eta**2,6)}  (ln2={mp.nstr(mp.log(2),6)})')
# finite k exact balanced exponent of the family (log-binomials)
from math import lgamma, log, exp
def lc(n,r): return (lgamma(n+1)-lgamma(r+1)-lgamma(n-r+1))/log(2)
def lsum(ls):
    m=max(ls); return m+log(sum(2**(x-m) for x in ls),2)
for k in [100,400,1600,6400]:
    n=2*k+1; par=k%2; um=(n-par)//2; best=1e99
    cq=[lc(n,j) for j in range(n+1)]; pre=[]; acc=[]
    for d in range(0,k+1):
        acc.append(cq[d]); lq=lsum(acc[-60:]) if d>60 else lsum(acc)
        # centered window of d+1 u-values
        lo=(um-d)//2; hi=lo+d
        ls=[lc(n,2*u+par) for u in range(um+1) if (u<lo or u>hi) and 2*u+par<=n]
        lp=lsum(ls) if ls else -1e9
        best=min(best,max(lp,lq))
    print(f'k={k}: balanced exponent log_N max(p,q) = {best/k:.4f}')
