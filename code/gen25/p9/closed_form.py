# G25-P9 task 1: closed form of Theorem A' exponent, and the "figure of merit" Phi/x.
# Phi_l := 1/e_l^2, x_l := log2 K_l, mu_l := Phi_l / x_l.  If the top uses K_l <= m^{1/a}, the exponent is
#   2 e_l W = 2 W / sqrt(mu_l x_l)  >=  2 W sqrt(a/(mu_l log m)),
# and c1 = limsup 2 e_l sqrt(a x_{l+1}) = 2 sqrt(a * (x_{l+1}/x_l) / mu_l)  ->  2 sqrt(2a/mu_inf)  (sawtooth sqrt2).
import math, thmA_explicit as T
x0, z0, q = 9, 2.0**-6, 0.5
L = 14
ts = T.sched_ts(x0, L, z0, q, 'hoeff')
lv = T.levels(x0, ts, 'sweep', True, True, T.RHO)
print("l  x_l      e_l      Phi_l     mu_l=Phi/x  logt_l/x_l  eff.rho_l   c1(a=2) at l  c1(a=4)")
for l in range(L):
    x, e, t, zs = lv[l]; xn, en = lv[l+1][0], lv[l+1][1]
    Phi = 1/e**2; mu = Phi/x
    lt = math.log2(ts[l])
    rho_eff = (en - 2.0**-x)/e
    print(f"{l:2d} {x:8.1f} {e:.5f} {Phi:10.1f} {mu:8.4f}   {lt/x:8.4f}   {rho_eff:.5f}   {2*e*math.sqrt(2*xn):7.3f}  {2*e*math.sqrt(4*xn):7.3f}")
# loss decomposition of mu: mu_{l+1}/mu_l = (Phi_{l+1}/Phi_l) / (x_{l+1}/x_l) = (2/(1+d_l)^2) / (2(1+logt/x))
print("\nideal (no overhead): base K0=2 (sweep e0=1/4, x0=1): mu=16 -> c1 = 2 sqrt(2a/16) = sqrt(a/2):",
      "a=2:", math.sqrt(1), " a=4:", math.sqrt(2), "  (and without the sawtooth: sqrt(a)/2... = ", math.sqrt(2)/2, math.sqrt(4)/2, ")")
mu_inf = 1/(lv[L][1]**2*lv[L][0])
print(f"mu at l={L}: {mu_inf:.4f}; 2 sqrt(2a/mu): a=2 {2*math.sqrt(4/mu_inf):.3f}  a=4 {2*math.sqrt(8/mu_inf):.3f};"
      f"  liminf (no sawtooth) 2 sqrt(a/mu): a=2 {2*math.sqrt(2/mu_inf):.3f} a=4 {2*math.sqrt(4/mu_inf):.3f}")
