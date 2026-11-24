# Crossover model for Claude's Theorem A recursion (Appendix A), paper constants vs re-optimised constants.
# exponent (bits) of the main term = 2 e_l * W, W = #wires ~ (3r-1) n (paper bounds W by k m = 3 r n).
# x_l = log2 K_l ; x_{l+1} = 2 (x_l + log2(g_l/K_l)) ; e_{l+1} = (rho + 3 zeta_l) e_l + 3/K_l.
# Top level: log2 m >= a_top * x_l + b_top.  Need 2 e_l W/n < tau.
import math
rho=1/math.sqrt(2)
def run(e0,x0,zeta,ratio,a_top,b_top,wn,tau,lmax=40,verbose=False):
    e,x=e0,x0
    for l in range(lmax):
        if 2*e*wn<tau:
            return l,x,e,a_top*x+b_top
        z=zeta(l); rr=ratio(l,x,z)
        e=(rho+3*z)*e+(3/2**x if x<1000 else 0.0)
        x=2*(x+math.log2(rr))
    return None
def paper_ratio(l,x,z): return 2**(3*l+18)
def hoeff_ratio(l,x,z):
    # minimal g/K with 2 zeta^2 g >= K ln(2 g K): g/K = t, need t >= ln(2 t K^2)/(2 z^2)
    t=1.0
    for _ in range(60): t=max(1.0,math.log(2*t)+2*x*math.log(2))/(2*z*z) if True else t
    return t
def gran_ratio(C):
    return lambda l,x,z: C/z
def sched(z0,q): return lambda l: z0*q**l
for dens,wn_paper,wn_chain in [(1,6,2),(4.26,6*4.26,3*4.26-1)]:
  print(f"\n=== density r={dens}: W/n paper bound {wn_paper:.1f}, chain wires {wn_chain:.2f}")
  for tau,name in [(1,'brute 2^n'),(0.716,'#3SAT 1.642^n'),(0.386,'PPSZ 2^0.386n')]:
    res=run(0.5,18,sched(2**-5,0.5),paper_ratio,9,0,wn_paper,tau)
    l,x,e,lm=res
    out=[f"paper: l={l} log2K={x:.0f} log2m*={lm:.3g}"]
    # (b) same recursion, base K0=2 (e0=1/4), Hoeffding-minimal ratio, best geometric zeta schedule, top a=2,b=15
    best=None
    for z0 in [2**-k for k in range(2,9)]:
      for q in [0.5,0.6,0.7,0.8]:
        for e0,x0 in [(0.25,1),(0.5,18)]:
          r=run(e0,x0,sched(z0,q),hoeff_ratio,2,15,wn_chain,tau)
          if r and (best is None or r[3]<best[0]): best=(r[3],r,z0,q,e0)
    out.append(f"opt-Hoeffding: l={best[1][0]} log2K={best[1][1]:.0f} log2m*={best[0]:.3g} (z0={best[2]},q={best[3]},e0={best[4]})")
    best=None
    for z0 in [2**-k for k in range(2,9)]:
      for q in [0.5,0.6,0.7,0.8]:
        r=run(0.25,1,sched(z0,q),gran_ratio(4),2,15,wn_chain,tau)
        if r and (best is None or r[3]<best[0]): best=(r[3],r,z0,q)
    out.append(f"opt-granularity(g/K=4/zeta): l={best[1][0]} log2K={best[1][1]:.0f} log2m*={best[0]:.3g}")
    # idealised: no balancing overhead at all (ratio 1, zeta 0), e0=1/4, top a=2
    r=run(0.25,1,lambda l:0.0,lambda l,x,z:1.0,2,15,wn_chain,tau)
    out.append(f"ideal(no overhead): l={r[0]} log2K={r[1]:.0f} log2m*={r[3]:.3g}")
    print(f" target {name}:\n   "+"\n   ".join(out))
