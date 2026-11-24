import numpy as np
np.set_printoptions(precision=3, suppress=True, linewidth=200)
N=4; z = np.load('vp_4_5_3.npy'); q=3
G = z[:q*16].reshape(q,4,4); Gp = z[q*16:].reshape(q,4,4)
F = np.zeros((16,16))
for x in range(4):
    for y in range(4): F[x*4+y, y*4+x] = 1
C = np.einsum('kxu,kyv->xyuv', G, Gp).reshape(16,16)
print("sing(F-C):", np.linalg.svd(F-C, compute_uv=False))
print("eig(CF):", np.sort_complex(np.linalg.eigvals(C@F)))
# kernel K of F - C: decompose into symmetric / antisymmetric parts
U,s,Vt = np.linalg.svd(F-C); K = Vt[5:].T   # 16 x 11
P = (np.eye(16)+F)/2
print("dim K∩Sym:", 11 - np.linalg.matrix_rank(K - F@K, tol=1e-8), " rank of sym-projection of K:", np.linalg.matrix_rank(P@K, tol=1e-8), " antisym-proj:", np.linalg.matrix_rank((np.eye(16)-P)@K, tol=1e-8))
# pencil invariants: min rank in span G, span G'
for name, S in (("G", G), ("G'", Gp)):
    rng = np.random.default_rng(0); 
    print(name, "generic rank in span:", np.linalg.matrix_rank(np.tensordot(rng.standard_normal(3), S, 1), tol=1e-8), " ranks of the three:", [np.linalg.matrix_rank(S[i], tol=1e-8) for i in range(3)])
    # cubic det(sum t_i G_i) quartic form: look for singular members: sample
print("products G_l G'_m traces:"); print(np.array([[np.trace(G[l]@Gp[m]) for m in range(3)] for l in range(3)]))
# span of G_l G'_m and commutation
A = np.array([ (G[l]@Gp[m]).ravel() for l in range(3) for m in range(3)]); print("dim span{G_l G'_m} =", np.linalg.matrix_rank(A, tol=1e-8))
A = np.array([ (Gp[l]@G[m]).ravel() for l in range(3) for m in range(3)]); print("dim span{G'_l G_m} =", np.linalg.matrix_rank(A, tol=1e-8))
