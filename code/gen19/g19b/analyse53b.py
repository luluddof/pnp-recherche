import numpy as np
np.set_printoptions(precision=3, suppress=True, linewidth=200)
z = np.load('vp_4_5_3.npy'); q=3
G = z[:q*16].reshape(q,4,4); Gp = z[q*16:].reshape(q,4,4)
for name, S in (("G", G), ("G'", Gp)):
    print(name, "rank vstack (common right kernel if <4):", np.linalg.matrix_rank(np.vstack(list(S)), tol=1e-8), " rank hstack (common left kernel if <4):", np.linalg.matrix_rank(np.hstack(list(S)), tol=1e-8))
# gauge: x->g, u->h^{-T}:  G_l -> g G_l h^{-1}, G'_l -> h G'_l g^{-1}
# find left kernels / right kernels
def lk(S): U,s,Vt = np.linalg.svd(np.hstack(list(S))); return U[:,-1], s
def rk(S): U,s,Vt = np.linalg.svd(np.vstack(list(S))); return Vt[-1], s
print("G  hstack sv", lk(G)[1], " vstack sv", rk(G)[1])
print("G' hstack sv", lk(Gp)[1], " vstack sv", rk(Gp)[1])
