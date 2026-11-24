# Own re-derivation of the cost exponent (natural logs, per n) for 3-occurrence 3-CNF (b=w=3, M<=3n, hD<=2M<=6n).
#   E/n = 2ln2/r + lnC_q/n + ERR,   q = 3 n deg/Q,  lnC_q/n = theta ln2 + H(theta) (theta=q/n<2/3; exact exponent of sum_{s<=q}2^s C(n,s))
#   constraint (Lemma 4.5): 12 r^{4Q} <= n   =>  Q <= (L - ln 12)/(4 ln r)
#   V3 (this work): deg = m-1, ERR = 12 ln3 * rho, rho <= 1/m (sum_g rank <= 2 h c_t <= 2hD/m <= 12 n/m)
#   V0 (P as written): deg = R-1, ERR = 2 tau ln(2D+1), 2 tau n = r h c^B <= 6 r eps n, eps = Pr[Bin(Q,(r-1)/2r) not in I], |I|=R, D=r^{2Q}
import math, numpy as np
from scipy.stats import binom
def taylor(theta):
    if theta>=2/3: return math.log(3)
    if theta<=0: return 0.0
    return theta*math.log(2)-theta*math.log(theta)-(1-theta)*math.log(1-theta)
def best_V3(L, pow3=False):
    best=(9e9,)
    r=2
    while True:
        Q=math.floor((L-math.log(12))/(4*math.log(r)))
        if Q<2: break
        ms=[float(3**a) for a in range(1,400) if 3**a<=Q] if pow3 else sorted(set(float(round(x)) for x in np.geomspace(2.0,float(max(3,Q)),600)))
        for m in ms:
            if m<2 or m>Q: continue
            e=2*math.log(2)/r+taylor(3*(m-1)/Q)+12*math.log(3)/m
            if e<best[0]: best=(e,r,Q,m)
        r=max(r+1,int(r*1.03))
        if r>1e60: break
    return best
def eps_tail(Q,r,R):
    p=(r-1)/(2*r); mu=Q*p; lo=int(round(mu-R/2)); lo=max(0,min(lo,Q-R+1)) if R<=Q else 0
    best=1.0
    for l in (lo-1,lo,lo+1):
        if l<0: continue
        best=min(best, binom.cdf(l-1,Q,p)+binom.sf(l+R-1,Q,p))
    return best
def best_V0(L):
    best=(9e9,)
    r=2
    while True:
        Q=int((L-math.log(12))/(4*math.log(r)))
        if Q<2 or r>1e9: break
        lnD=2*Q*math.log(r)
        for R in sorted(set(int(round(x)) for x in np.geomspace(2,Q,200))):
            if R>=Q: continue
            e=2*math.log(2)/r+taylor(3*(R-1)/Q)+6*r*eps_tail(Q,r,R)*(lnD+math.log(2))
            if e<best[0]: best=(e,r,Q,R)
        r=max(r+1,int(r*1.05))
    return best
if __name__=="__main__":
  print("log2n        L2/sqrtL   V3 E/n     E/n /(L2/sqrtL)  r Q m  | m=3^a: E/n ratio m | V0(P-style) E/n  /(L2^2/sqrtL)  r Q R | V0/V3 | E/n<ln2?")
  for lg in [1e3,1e4,1e5,1e6,1e8,1e10,1e12,1e15,1e20,1e30,1e50,1e100]:
      L=lg*math.log(2); L2=math.log(L); u=L2/math.sqrt(L)
      b=best_V3(L); b3=best_V3(L,True)
      v0=best_V0(L) if lg<=1e15 else None
      s=f"{lg:8.0e}  {u:10.3e} {b[0]:10.3e} {b[0]/u:8.3f}   {b[1:]}  | {b3[0]:10.3e} {b3[0]/u:7.3f} {b3[3] if len(b3)>1 else None}"
      if v0: s+=f" | {v0[0]:10.3e} {v0[0]/(L2*u):7.3f} {v0[1:]} | {v0[0]/b[0]:6.2f}"
      s+=f" | {b[0]<math.log(2)}"
      print(s,flush=True)
  # asymptotic constant check: min_c 1.5c + 12ln3/c times sqrt2
  c=math.sqrt(12*math.log(3)/1.5); print("predicted constant:", math.sqrt(2)*(1.5*c+12*math.log(3)/c), " c=",c)
