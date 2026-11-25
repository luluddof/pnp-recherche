# G25-P10: coset window (G24f) read as a slice of the (p,q) trade-off.
# Dilute regime: p ~ n^a, q ~ n^b, N ~ n^kappa, a = E0+hb, b = E0+ha-1, cost 4rho = (2E0+H)/kappa = theta+T+1/kappa.
import math
R2=math.sqrt(2)
print("E0 ha hb | ThmC: 4rho>= | labels: 4rho>= | best proved | open=>subcircle? | spreading closes?")
for E0 in range(2,11):
    for H in range(0,10):
        if H >= (2*R2-2)*E0: continue
        for ha in range(0,H+1):
            hb=H-ha; a=E0+hb; b=E0+ha-1
            thmC=(2*E0+H)/((a+b)/2)          # kappa <= (a+b)/2 from theta+T>=2
            lab=(2*E0+H)/E0                   # kappa <= E0
            best=max(thmC,lab)
            if best>=2*R2: continue           # already closed
            sub = (a-b)**2 <= 2*(a+b)+1       # theta* <= theta_circle(slope)
            # spreading: 2*max+min >= 4  => kappa <= (2*max(a,b)+min(a,b))/4
            ksp=(2*max(a,b)+min(a,b))/4
            closes = (2*E0+H)/ksp >= 2*R2-1e-12
            print(f"{E0:2d} {ha:2d} {hb:2d} | {thmC:.3f} | {lab:.3f} | {best:.3f} | {sub} | {closes} ({(2*E0+H)/ksp:.3f})")
