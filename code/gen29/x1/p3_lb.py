# Compression (Theorem-C type) lower bound for P_3 vs ring bonds, in exponents base N.
# strands: s1 = a1->b2 (L->M), s2 = a2->b3 (M->R), s3 = a3->b1 (R->L); compress s_i to N^{x_i}.
# node L: top x1, bottom x3; M: top x2, bottom x1; R: top x3, bottom x2. Bonds b12 (L-M), b23 (M-R), b31 (R-L).
# rank_rows(C) <= N^{min_sigma cost}; bound theta >= sum x when sum x > min cost. Returns max sum x (normalized /3).
import numpy as np, itertools, sys
g=np.linspace(0,1,51)
X1,X2,X3=np.meshgrid(g,g,g,indexing='ij')
def lb(b12,b23,b31):
    S=X1+X2+X3; best=np.full(S.shape,np.inf)
    legs={'L':(X1,X3),'M':(X2,X1),'R':(X3,X2)}
    for sig in itertools.product((0,1),repeat=3):  # 0 = top, 1 = bottom
        sL,sM,sR=sig; c=np.zeros(S.shape)
        for s,(top,bot) in zip(sig,(legs['L'],legs['M'],legs['R'])): c=c+(bot if s==0 else top)
        c=c+b12*(sL!=sM)+b23*(sM!=sR)+b31*(sR!=sL)
        best=np.minimum(best,c)
    ok=S>best+1e-12
    return (S[ok].max() if ok.any() else 0)/3
for v in [0.05,0.1,0.2,1/3,0.4,0.5]:
    b=v  # symmetric: q = N^{3b}, v = b
    print(f"symmetric bonds b=N^{b:.3f} (v={v:.3f}): u >= {lb(b,b,b):.3f}   line 1-v={1-v:.3f}  1-4v/3={1-4*v/3:.3f}")
for v in [0.1,0.2,1/3]:
    print(f"bonds (3v,0,0) v={v:.3f}: u >= {lb(3*v,0,0):.3f};  (1.5v,1.5v,0): {lb(1.5*v,1.5*v,0):.3f}; line {1-v:.3f}")
print("--- min over bond patterns with b12+b23+b31 = 3v (grid step 0.05 in pattern) ---")
for v in [0.05,0.1,0.2,0.3,0.4,0.45,0.5,0.6]:
    worst=(9,None)
    for w1 in np.linspace(0,1,21):
        for w2 in np.linspace(0,1-w1+1e-9,int(round((1-w1)*20))+1):
            w3=max(0,1-w1-w2); r=lb(3*v*w1,3*v*w2,3*v*w3)
            if r<worst[0]: worst=(r,(round(w1,2),round(w2,2),round(w3,2)))
    print(f"v={v:.2f}: min_patterns u_lower={worst[0]:.3f} at pattern {worst[1]}  (line {1-v:.3f})")
