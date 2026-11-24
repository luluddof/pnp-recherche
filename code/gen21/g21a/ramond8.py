# Smallest Ramond instance: n = 8 fermions, fields psi^x of weight 1/2 (H = 1), labels in the chiral ground space S+ (dim 8).
# {psi^x_0, psi^v_0} = delta: R = <u|g_v g_x|y>, C = <u|g_x g_v|y>, R + C = 2 [x=v][y=u].  Exact ranks (flint).
import flint, itertools
I2 = [[1,0],[0,1]]; X = [[0,1],[1,0]]; Z = [[1,0],[0,-1]]; E = [[0,1],[-1,0]]
def kron(*ms):
    out = [[1]]
    for m in ms:
        out = [[a*b for a in ra for b in rb] for ra in out for rb in m]
    return out
def mul(A, B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
L = [kron(E,X,I2), kron(E,Z,I2), kron(I2,E,X), kron(I2,E,Z), kron(X,I2,E), kron(Z,I2,E), kron(E,E,E)]
I8 = kron(I2,I2,I2)
G = [kron(E, l) for l in L] + [kron(X, I8)]
n = 16
for i, j in itertools.product(range(8), repeat=2):
    A = mul(G[i], G[j]); B = mul(G[j], G[i])
    assert all(A[a][b] + B[a][b] == (2 if (i == j and a == b) else 0) for a in range(n) for b in range(n))
print("8 real symmetric anticommuting 16x16 gamma matrices: Clifford relations verified exactly")
N = 8  # labels y,u in S+ = first 8 coordinates
R = {}; C = {}
for x in range(8):
    for v in range(8):
        A = mul(G[v], G[x]); B = mul(G[x], G[v])
        for y in range(8):
            for u in range(8):
                R[x,y,u,v] = A[u][y]; C[x,y,u,v] = B[u][y]
assert all(R[k] + C[k] == 2*(k[0] == k[3])*(k[1] == k[2]) for k in R)
def rr(T): return flint.fmpz_mat(64, 64, [T[x,y,u,v] for x in range(8) for y in range(8) for u in range(8) for v in range(8)]).rank()
def rc(T): return flint.fmpz_mat(64, 64, [T[x,y,u,v] for x in range(8) for u in range(8) for y in range(8) for v in range(8)]).rank()
print("R + C = 2 [x=v][y=u] exactly.  rr(R), rc(R) =", rr(R), rc(R), "  rr(C), rc(C) =", rr(C), rc(C))
print("=> feasible pair (p,q) = (rr(R), rc(C)) =", (rr(R), rc(C)), "at N = 8 (Hurwitz point); as PIECES (Def. 7): log2 rr + log2 rc =", 3+6, "of 12 bits (rho = 3/4)")
