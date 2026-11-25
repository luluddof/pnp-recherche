# alpha achievable by self-similar ansatz e(x)=x^-alpha iff min_w F(w^-a,(1-w)^-a) < 1 (F geometric, no duality)
# prints slack = 1 - min ratio (positive => alpha achievable with margin)
import numpy as np, sys
sys.path.insert(0,'.'); from r1dp import Hull, fam
names=sys.argv[1:] or ['circle','line','lp1.5','p2q','circ+lin0.3','circ+lin0.1']
for name in names:
    H=Hull(fam(name))
    out=[]
    for a in [0.45,0.5,0.55,0.6,0.667,0.7,0.8,0.9,0.95,0.99]:
        lk=np.linspace(0,-300,60001); k=10.0**lk
        ea=(1+k)**a; eb=((1+k)/k)**a          # x_a=1/(1+k), x_b=k/(1+k), e=x^-a
        r=np.minimum(H.s_star(ea,eb),H.s_star(eb,ea)); i=np.argmin(r)
        out.append(f"a={a}:{1-r[i]:+.2e}@k=1e{lk[i]:.1f}")
    print(name, ' '.join(out),flush=True)
