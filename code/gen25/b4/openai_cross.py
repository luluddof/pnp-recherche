# Crossover model for OpenAI App. D (tensor contraction), root-dominated:
# exponent (bits) ~ gamma * c * M / sqrt(u0),  u0 = log2(s)/beta,  M ~ s ~ w n  ->  need gamma c w / sqrt(u0) < tau
# => log2 s* = beta * (gamma c w / tau)^2.   (leaves 2^{s/sqrt(lambda)} are negligible at these sizes)
import math
def ls(beta,gamma,c,w,tau): return beta*(gamma*c*w/tau)**2
rows=[
 ("A paper: eps=2^-17, kappa/taper as written (26.7 nat=38.5 bit), nu=6, depth sum 1.215", 2**17, 1.215, 26.7/math.log(2)),
 ("B +u0 at its legal max under (D.6) with r=64, Delta=s^0.9 (beta=960)", 960, 1.215, 26.7/math.log(2)),
 ("C +LP-optimal taper & full packets nu=2 (6.2*sqrt2 nat = 12.6 bit)", 960, 1.215, 6.2*math.sqrt(2)/math.log(2)),
 ("D +Smith/Mahler modes (n->n/4: c/2) [unverified as algorithm]", 960, 1.215, 3.1*math.sqrt(2)/math.log(2)),
 ("E +r=32, Delta=s^0.5 (beta=(16*5)/0.5=160, depth sum 1/(1-6sqrt2/32)=1.36)", 160, 1.36, 3.1*math.sqrt(2)/math.log(2)),
 ("F idealised (K1): u0 <= log2 M/(1+2log2 r), r=32 (beta=11)", 11, 1.36, 3.1*math.sqrt(2)/math.log(2)),
]
for dens,w in [(1,2.0),(4.26,11.78)]:
  print(f"\n=== density {dens}, wires/variable w={w}")
  for tau,name in [(1,'brute'),(0.716,'#3SAT 1.642^n'),(0.386,'PPSZ')]:
    print(f" target {name}")
    for nm,beta,g,c in rows:
      print(f"   {nm:85s} log2 s* = {ls(beta,g,c,w,tau):.3g}")
# Prop B' design lower bound check: c >= 0.3 bit (crude), so even F cannot go below:
print("\nDesign floor with Prop.B' crude constant c=0.3 bit, gamma=1, beta=11, w=11.78:", ls(11,1,0.3,11.78,1))
