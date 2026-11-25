import math, numpy as np, b7lib as B, chain_b7 as C
mu, s2, name, ch = C.best_chain('row-R3', 'new')
print(name, mu)
G=1.0; A=ch[0][0]
for l in range(len(ch)-1):
    x,e = ch[l]; x2,e2 = ch[l+1]; lt = x2/2 - x
    t=2**lt; z=float(B.zeta_new(lt,x)); lin=4/t
    G *= (C.RHO+z+2.0**-(x+lt))/C.RHO
    A += 2.0**-(l+1)*0  # placeholder
    if l<10: print(f"l={l} x={x:.1f} lt={lt:.2f} zeta={z:.2e} (lin 4/t={lin:.2e}, 1/K={2.0**-x:.1e}) entry/ (e rho)={2.0**-x/(e*C.RHO):.2e} Gamma_cum={G:.4f} mu={1/(e2*e2*x2):.3f}")
