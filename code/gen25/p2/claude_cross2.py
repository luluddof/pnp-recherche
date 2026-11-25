# Ablation of the crossover of Claude's Theorem A recursion (Appendix A).
# Main exponent (bits) = 2 e_l W, W = #wires (paper: W <= k m = 3 r n; chain: 3m - n).
# x_{l+1} = 2 (x_l + log2 t_l), t_l = g_l/K_l ; e_{l+1} = (rho + 3 zeta_l) e_l + 3/K_l ; top: log2 m >= a x_l + b.
import math, itertools
rho=1/math.sqrt(2)
def run(e0,x0,zeta,ratio,a,b,wn,tau,lmax=60):
    e,x=e0,x0
    for l in range(lmax):
        if 2*e*wn<tau: return dict(l=l,x=x,e=e,lm=a*x+b)
        z=zeta(l); t=ratio(l,x,z)
        e=(rho+3*z)*e+(3/2**x if x<1000 else 0.0)
        x=2*(x+math.log2(t))
    return None
paper_ratio=lambda l,x,z: 2**(3*l+18)
def hoeff(l,x,z):            # least t=g/K with 2 z^2 t >= ln(2 t K^2), K=2^x
    t=1.0
    for _ in range(100): t=max(1.0,(math.log(2*t)+2*x*math.log(2))/(2*z*z))
    return t
def best(e0f,x0s,ratio,a,b,wn,tau):
    B=None
    for x0 in x0s:
      for k in range(2,10):
        for q in [0.5,0.6,0.7,0.8,0.9]:
          r=run(e0f(x0),x0,lambda l,z0=2.0**-k,q=q: z0*q**l,ratio,a,b,wn,tau)
          if r and (B is None or r['lm']<B['lm']): B=dict(r,x0=x0,z0=2.0**-k,q=q)
    return B
brute_e0=lambda x0: 0.5
cut_e0=lambda x0: 0.25+0.5/2**x0
for dens in [1,4.26]:
  wpap=3*dens; wch=3*dens-1
  print(f"\n=== density r={dens} (W/n: paper bound {wpap:.2f}, chain {wch:.2f})")
  for tau,name in [(1,'brute 2^n'),(0.716,'#3SAT 1.642^n'),(0.386,'PPSZ 2^.386n')]:
    rows=[]
    r=run(0.5,18,lambda l:2.0**(-l-5),paper_ratio,9,0,wpap,tau); rows.append(("A paper as written",r))
    r=run(0.5,18,lambda l:2.0**(-l-5),paper_ratio,9,0,wch,tau); rows.append(("B +chain wire count",r))
    r=best(brute_e0,range(4,21),hoeff,9,0,wch,tau); rows.append(("C +exact Hoeffding g/K, x0 & zeta optimised",r))
    r=best(cut_e0,range(4,21),hoeff,9,0,wch,tau); rows.append(("D +base by cut-width sweep (e0=1/4)",r))
    r=best(cut_e0,range(4,21),hoeff,2,15,wch,tau); rows.append(("E +top level m~K^2 (random clustering)",r))
    r=best(cut_e0,range(4,21),lambda l,x,z:4/z,2,15,wch,tau); rows.append(("F +granularity balancing g/K=4/zeta [speculative]",r))
    print(f" target {name}")
    for nm,r in rows: print(f"   {nm:52s} levels={r['l']:2d} log2K={r['x']:9.0f} log2 m* = {r['lm']:.3g}"+(f"  (x0={r['x0']},z0={r['z0']},q={r['q']})" if 'x0' in r else ''))
