# Independent check of P13's level q (fixed minimax w, closed form) against the EXACT model chord:
# cost(d) = max_th min(ea((1-d')cos th + d' + (1+d) al), eb((1-d')sin th + d' + (1+d) be)) + (1+d) ent,  d' = min(d,1)
# q_exact = sup_d (cost(d) - cost(0))/(d cost(0)).  Must have q_exact <= q_P13 (P13's is a linear upper bound).
import numpy as np, math
rng = np.random.default_rng(5)
TH = np.linspace(0, math.pi/2, 40001)
D = np.concatenate([np.geomspace(1e-6, 1, 300), np.geomspace(1, 1e3, 100)[1:]])
def cost(ea, eb, al, be, ent, d):
    # exact max-min: f1 = ea(A cos th + B) decreasing, f2 = eb(A sin th + C) increasing; crossing R A cos(th+phi) = eb C - ea B
    dp = min(d, 1.0); A = 1 - dp; B = dp + (1+d)*al; C = dp + (1+d)*be
    if ea*(A + B) <= eb*C: return ea*(A + B) + (1+d)*ent
    if ea*B >= eb*(A + C): return eb*(A + C) + (1+d)*ent
    R = math.hypot(ea, eb); phi = math.atan2(eb, ea)
    th = math.acos(max(-1.0, min(1.0, (eb*C - ea*B)/(R*A)))) - phi
    th = min(max(th, 0.0), math.pi/2)
    return min(ea*(A*math.cos(th) + B), eb*(A*math.sin(th) + C)) + (1+d)*ent
def p13(ea, eb, al, be, ent):
    R = math.hypot(ea, eb); phi = math.atan2(eb, ea)
    th = min(max(math.acos(max(-1, min(1, (eb*be - ea*al)/R))) - phi, 0), math.pi/2)
    v = min(ea*(math.cos(th) + al), eb*(math.sin(th) + be))
    w = eb*math.cos(th)/(ea*math.sin(th) + eb*math.cos(th))
    if ea*(1+al) <= eb*be: v, w = ea*(1+al), 1.0
    if ea*al >= eb*(1+be): v, w = eb*(1+be), 0.0
    Ew = math.hypot(w*ea, (1-w)*eb); Sw = w*ea + (1-w)*eb; Ow = w*ea*al + (1-w)*eb*be
    e2 = max(v, Ew + Ow) + ent
    return e2, (Sw - Ew + Ow + ent)/e2
worst = -1; rows = []
for it in range(400):
    ea, eb = rng.uniform(0.01, 0.3, 2)
    if it % 4 == 0: eb = ea*rng.uniform(0.9, 1.1)
    al, be = rng.uniform(0, 0.08, 2); ent = 2.0**-rng.uniform(6, 20)
    e2, q = p13(ea, eb, al, be, ent)
    c0 = cost(ea, eb, al, be, ent, 0.0)
    qe = max((cost(ea, eb, al, be, ent, d) - c0)/(d*c0) for d in D)
    rows.append((q, qe, e2, c0)); worst = max(worst, qe - q)
r = np.array(rows)
print(f"level: max(q_exact - q_P13) = {worst:.2e} (must be <= ~0);  |e_P13 - cost(0)| max {np.abs(r[:,2]-r[:,3]).max():.2e}")
print(f"  q_P13 range {r[:,0].min():.3f}..{r[:,0].max():.3f}; q_exact range {r[:,1].min():.3f}..{r[:,1].max():.3f}; mean gap {np.mean(r[:,0]-r[:,1]):.3f}")
