"""Growth profiles p_S(E) and the two conditions of the window framework:
 (cost)  2*L(2h)/(4*L(h))  = 'rho' the Appendix-B calculation would give,
 (rank)  sum_s rank_row(A_s) >= N^2 needs (4h+1)(2h+1)*max_{w<=4h} p_S(w) >= p_S(h)^2.
"""
import math
def counts(S, Emax):
    p = [0]*(Emax+1); p[0] = 1
    for k in S:
        if k > Emax: continue
        for e in range(k, Emax+1): p[e] += p[e-k]
    return p
def lg(x): return math.log2(x) if x > 0 else float('-inf')
def table(name, Sgen, hs):
    Emax = 4*max(hs)
    p = counts(Sgen(Emax), Emax)
    first_viol = None
    pm = []; mx = 0
    for v in p:
        mx = max(mx, v); pm.append(mx)
    for h in range(2, max(hs)+1):
        if p[h] < 2: continue
        if (4*h+1)*(2*h+1)*pm[4*h] < p[h]**2:
            first_viol = h; break
    print(f"== S = {name}; first h with rank-count violation: {first_viol}")
    print("      h   lam=log2 N   L(2h)/L(h)  rho_calc=L(2h)/2L(h)   L(4h)/L(h)  slack=log2[(4h+1)(2h+1)D(4h)]-2lam")
    for h in hs:
        L1, L2, L4 = lg(p[h]), lg(p[2*h]), lg(max(p[:4*h+1]))
        if L1 <= 0: continue
        print(f" {h:6d} {L1:10.2f} {L2/L1:10.3f} {L2/(2*L1):12.3f} {L4/L1:14.3f} {lg((4*h+1)*(2*h+1))+L4-2*L1:12.1f}")
hs = [10, 30, 100, 300, 1000, 3000, 10000]
table("cubes", lambda E: [k**3 for k in range(1, int(E**(1/3))+3)], hs + [30000])
table("powers of 2", lambda E: [2**k for k in range(0, 40)], hs + [30000])
table("parts <= 8", lambda E: range(1, 9), hs)
table("parts <= 30", lambda E: range(1, 31), hs)
table("odd parts", lambda E: range(1, E+1, 2), hs)
