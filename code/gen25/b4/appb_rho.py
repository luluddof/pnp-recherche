# Effective per-level rho of Appendix B pieces at crossing size lambda bits (N = 2^lambda).
# 4 rho_eff(lambda) ~ max_s [log2 P(2h+s-1) + log2 P(2h-s)] / lambda,  lambda = log2 p(h)  (g20b counting formula)
import math
def partitions(nmax):
    p=[0]*(nmax+1); p[0]=1
    for n in range(1,nmax+1):
        s=0;k=1
        while True:
            g1=k*(3*k-1)//2
            if g1>n: break
            sg = 1 if k%2 else -1
            s += sg*p[n-g1]
            g2=k*(3*k+1)//2
            if g2<=n: s+= sg*p[n-g2]
            k+=1
        p[n]=s
    return p
NMAX=6000
p=partitions(NMAX)
P=[0]*(NMAX+1); acc=0
for i in range(NMAX+1): acc+=p[i]; P[i]=acc
l2=lambda x: math.log2(x)
def lp_asym(n): return (math.pi*math.sqrt(2*n/3)-math.log(4*n*math.sqrt(3)))/math.log(2)
def lP_asym(n): return lp_asym(n)+l2(math.sqrt(6*n)/math.pi)
def lP(n): return l2(P[n]) if n<=NMAX else lP_asym(n)
def lpp(n): return l2(p[n]) if n<=NMAX else lp_asym(n)
print("check asym at n=6000: exact log2 P", l2(P[6000]), "asym", lP_asym(6000))
print(" lambda      h        4rho_eff   rho_eff   rho_eff/(1/sqrt2)")
for lam in [8,16,32,64,128,256,512,1024,4096,2**14,2**16,2**20,2**31,2**40,2**64]:
    # least h with log2 p(h) >= lam
    lo,hi=1,10
    while (lpp(hi) if hi<=NMAX else lp_asym(hi))<lam: hi*=2
    while lo<hi:
        mid=(lo+hi)//2
        if lpp(mid)>=lam: hi=mid
        else: lo=mid+1
    h=lo
    best=max(lP(2*h+s-1)+lP(2*h-s) for s in range(-3,5))
    r4=min(best,4*lam)/lam
    print(f"{lam:>10} {h:>12} {r4:10.5f} {r4/4:9.5f} {r4/4*math.sqrt(2):9.5f}")
