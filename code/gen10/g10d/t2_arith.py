"""G10d test 2: Cor B3 arithmetic (own first-moment + mixing), fine grid in d; B1 threshold; H1 growth rate."""
import math
def H2(a): return 0.0 if a<=0 or a>=1 else -(a*math.log2(a)+(1-a)*math.log2(1-a))
def fm(a,d):  # log2 E[#IS of size a n]/n in pairing model (derived independently: see report)
    return H2(a)+d*((1-a)*math.log2(1-a)-0.5*(1-2*a)*math.log2(1-2*a))
def alpha_fm(d):
    # largest a in (0,1/2) with fm(a)>=0 : scan from 1/2 downward
    lo=None; N=200000
    for i in range(N-1,0,-1):
        a=0.5*i/N
        if fm(a,d)>0: lo=a; break
    hi=lo+0.5/N
    for _ in range(60):
        mid=(lo+hi)/2
        if fm(mid,d)>0: lo=mid
        else: hi=mid
    return hi
best=(0,None); rows=[]
ds=sorted(set([3,4,5,6,7,8,10,12,16,20,24,32,40,48,56,64,72,80,96,128,192,256,512,1024,4096,65536,2**20]+list(range(40,121,4))))
for d in ds:
    a=alpha_fm(d); lam=2*math.sqrt(d-1); mm=d/2
    nb_full=a*(d*a+lam*(1-a))/2          # one-sided mixing, exact form
    nb_loose=a*(d*a+lam)/2               # prover's form (drops 1-a)
    U=lambda nb: 2*nb*max(1.0,math.log2(3*mm/nb))
    rows.append((d,a,nb_full,U(nb_full),U(nb_loose)))
    if U(nb_loose)>best[0]: best=(U(nb_loose),d)
for r in rows:
    if r[0] in (3,4,8,16,32,48,56,64,72,80,128,256,1024,4096,65536,2**20): print("d=%7d alpha/n=%.5f nbar/n=%.5f U_A/n(exact mixing)=%.3f U_A/n(prover form)=%.3f  asym 2ln d/d=%.5f"%(r+(2*math.log(r[0])/r[0],)))
print("max over grid of prover-form U_A/n: %.3f at d=%d"%best)
# B1 threshold: alpha = n/(log 2d loglog 4d): value of h* log(4m/h*) / n
print("B1/B4 threshold: (h*/n) log2(4m/h*) at alpha/n = 1/(log 2d loglog 4d)")
for d in [2,4,16,256,2**16,2**32,2**64,2.0**256,2.0**1000]:
    a=1/(math.log2(2*d)*math.log2(math.log2(4*d)))
    a=min(a,0.5); h=H2(a)
    print("  d=2^%-5g alpha/n=%.4f h/n=%.4f  exponent/(e n)=%.3f"%(math.log2(d),a,h,h*math.log2(2*d/h)))
# B2: Hoffman
print("B2: (h*/n) log(4m/h*) at alpha/n = lam/(d+lam), lam = d/log^2 d")
for d in [16,256,2**16,2**32,2.0**64,2.0**256]:
    lam=d/math.log2(d)**2; a=lam/(d+lam); h=H2(a)
    print("  d=2^%-4g alpha/n=%.5f  exponent/(e n)=%.3f"%(math.log2(d),a,h*math.log2(2*d/h)))
print("H1 growth: log2(x)*D/log2(D), x^(D+1)=x^D+1")
for D in [2,3,4,8,16,64,256,1024,65536]:
    lo,hi=1.0,2.0
    for _ in range(200):
        x=(lo+hi)/2
        if D*math.log(x)+math.log(x-1)>0: hi=x
        else: lo=x
    print("  D=%6d  %.4f"%(D, math.log2(x)*D/math.log2(D)))
