import numpy as np
rng = np.random.default_rng(5)
g=30; c=np.zeros((g,g)); hub=[3]; c[hub,:]=1; c[:,hub]=1; np.fill_diagonal(c,0)
U=c.sum(); mu=U/(g*(g-1)); ca=c.sum(1); al=(ca-(g-1)*mu)/(g-2); ga=c-mu-al[:,None]-al[None,:]; np.fill_diagonal(ga,0)
print("row sums", np.abs(ga.sum(1)).max())
print("gamma values: hub row", np.unique(np.round(ga[3],6)), " other", np.unique(np.round(ga[5],6)))
p=rng.permutation(g)
w=lambda j:(g-j)*(g-j-1)/(g*(g-1))
for j in range(0,6):
    S=list(p[:j]); T=list(p[j:]); G0=ga[np.ix_(S,T)].sum()
    Gm=max(abs(ga[v,S].sum()) for v in T) if j else 0
    print(j, "hub in S" if 3 in S else "", "Gamma",round(G0,6),"G_j",round(Gm,6))
