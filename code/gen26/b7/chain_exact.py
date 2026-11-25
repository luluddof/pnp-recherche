# Exact (non-binned) long chains, asymptotic mu_inf = lim Phi_l/x_l, for
#  (i) p2 chain (clique-sweep base K0 = 2^x0, Lemma 8 grouping at every level),
#  (ii) row-sweep level 1 (g x g grid, no grouping) then the p2 recursion,
# and the resulting constants: single schedule (sawtooth) c1 = limsup 2 e_l sqrt(a x_{l+1});
# no sawtooth (unbalanced top levels, Lemma U) c1 -> 2 sqrt(a/mu_inf).
import math, thmA_explicit as T
RHO = T.RHO; LN2 = math.log(2)

def zs_sweep(lg):
    g = 2.0**lg
    return math.sqrt(math.log(2*g*g + 2*g + 1)/(4*g))*(1 + 1/g) + 1/(4*(g-1)) + 2/g

def chain(x_start, e_start, z0, q, L=26):
    out = [(x_start, e_start)]; x, e = x_start, e_start
    for l in range(L):
        z = z0*q**l; t = T.t_for(z, x)
        zs = T.zsum(t, x, True); g = t*2.0**x if x < 1000 else float('inf')
        inv_g = 1/g if x < 1000 else 0.0; invK = 2.0**-x if x < 1000 else 0.0
        e = e*(RHO + inv_g + zs) + invK
        x = 2*(x + math.log2(t)); out.append((x, e))
    return out

def summarize(ch, lo=18):
    mu = [1/(e*e*x) for x, e in ch]
    saw2 = max(2*ch[l][1]*math.sqrt(2*ch[l+1][0]) for l in range(lo, len(ch)-1))
    saw4 = max(2*ch[l][1]*math.sqrt(4*ch[l+1][0]) for l in range(lo, len(ch)-1))
    return mu[-1], saw2, saw4

best = {}
for kind in ('clique', 'row'):
    starts = []
    if kind == 'clique':
        for x0 in range(4, 16):
            K = 2**x0; starts.append((f"x0={x0}", float(x0), (((K+1)**2)//4 - 1)/K**2))
    else:
        for lg in [x/2 for x in range(24, 81)]:
            e1 = (RHO + 2.0**-lg)/4 + zs_sweep(lg)
            starts.append((f"g=2^{lg}", 2*lg, e1))
    for name, xs, es in starts:
        for k in range(2, 12):
            for q in (0.4, 0.5, 0.6, 0.7, 0.8):
                ch = chain(xs, es, 2.0**-k, q)
                mu, s2, s4 = summarize(ch)
                if kind not in best or mu > best[kind][0]:
                    best[kind] = (mu, s2, s4, name, 2.0**-k, q, ch)
for kind, (mu, s2, s4, name, z0, q, ch) in best.items():
    print(f"{kind:7s} best mu_inf = {mu:.4f} ({name}, z0={z0}, q={q})")
    print(f"        single schedule (sawtooth) c1: a=2 {s2:.3f}  a=4 {s4:.3f}")
    print(f"        no sawtooth 2 sqrt(a/mu):      a=2 {2*math.sqrt(2/mu):.3f}  a=4 {2*math.sqrt(4/mu):.3f}")
    print("        levels (x, e, mu):", [(round(x, 1), round(e, 5), round(1/(e*e*x), 4)) for x, e in ch[:9]])
