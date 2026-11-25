# Explicit Lemma-U chain with a BRUTE-FORCE small child (e_b = 1/2, K_b = 2^y) and corner slope c (H inside v <= c(1-u)).
# Overheads per Lemma 8(ii): t = group factor, zeta = sqrt(ln(2 g K)/(2t)); straight edges 1/h; entry 2^{1-y}.
import math
for c in (1.0,2.0,6.6,50.0):
    x=20.0; e=0.5; out=[]; nxt=1e3
    while x<1e9:
        y=math.ceil(2*math.log2(x)+math.log2(c)+4)
        t=x**3; ln=math.log(2)*(2*x+3*math.log2(x))       # ln(2 g K) <= ln 2^(2x')
        z=math.sqrt(ln/(2*t))
        xa=x; x=xa+y+2*math.log2(t)                         # x' = x_a + x_b + log t_a + log t_b
        eb=0.5
        Fm=e*c/(c+e/eb)                                     # ray exit on v=c(1-u)
        e=Fm*(1+2*z+2.0**-y+(2*z+2.0**-xa)/c)+2.0**(1-y)
        if x>=nxt: out.append(f"x=1e{math.log10(x):.0f}: e*x={e*x:.3g} e*x/log2x={e*x/math.log2(x):.3g}"); nxt*=10
    print(f"c={c}: "+' | '.join(out))
