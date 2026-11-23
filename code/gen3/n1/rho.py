"""N1 brick (re-targeted): crossing decompositions and the per-level factor rho.
Paper: claude_3sat.txt. Crossing X(x,y,u,v) = [x=v][y=u].
Row split (x,y)|(u,v) ; column split (x,u)|(y,v).  A function costs c bits in a direction
if its matrix there factors through 2^c (rank <= 2^c over Q, integer after scaling)."""
import itertools, math
import numpy as np

def sigma3(p, q):
    phi = sum(p[i]*q[j] for i in range(3) for j in range(i, 3)) \
        + q[0]*p[1]*p[2] + p[0]*q[1]*p[2] + p[0]*p[1]*q[2]
    return -1 if phi % 2 else 1

B3 = list(itertools.product([0, 1], repeat=3))

def verify_eight_square():
    bad = 0
    for x, y, u, v in itertools.product(B3, repeat=4):
        lhs = 2 * (x == v) * (y == u)
        xor = all((a ^ b ^ c ^ d) == 0 for a, b, c, d in zip(x, y, u, v))
        rhs = xor * (sigma3(x, y)*sigma3(v, u) + sigma3(x, u)*sigma3(v, y))
        bad += (lhs != rhs)
    return bad

def crossing(N):
    X = np.zeros((N,)*4)
    for x in range(N):
        for y in range(N):
            X[x, y, y, x] = 1
    return X

# matrices arranged so that the crossing is the IDENTITY in both directions:
def rowmat(A):   # rows (x,y), cols (v,u)
    N = A.shape[0]; return A.transpose(0, 1, 3, 2).reshape(N*N, N*N)
def colmat(A):   # rows (x,u), cols (v,y)
    N = A.shape[0]; return A.transpose(0, 2, 3, 1).reshape(N*N, N*N)

def rank(M, tol=1e-8):
    return int(np.linalg.matrix_rank(M, tol=tol))

def classes(M, kappa):
    """(z, k, m): geometric mult. of eigenvalue 0, of kappa, rest (assumes diagonalizable)."""
    n = M.shape[0]
    z = n - rank(M); k = n - rank(M - kappa*np.eye(n))
    ev = np.linalg.eigvals(M)
    az = int(np.sum(abs(ev) < 1e-7)); ak = int(np.sum(abs(ev - kappa) < 1e-7))
    diag = (az == z and ak == k)     # necessary check on the two special eigenvalues
    return (z, k, n - z - k), diag

# exponent of rank of grouped piece A_tau (fraction tau of blocks taking piece C), per block:
# eigen type 0 -> factor kappa*T (degree 1); type kappa -> degree 0; mixed -> degree 0 or 1.
_g = np.linspace(0, 1, 401)
PZ, PM = np.meshgrid(_g, _g, indexing='ij'); OK = PZ + PM <= 1 + 1e-12
PK = np.clip(1 - PZ - PM, 0, 1)
def _H(p): return np.where(p > 0, -p*np.log2(np.where(p > 0, p, 1)), 0)
def f_curve(cls, taus):
    z, k, m = cls
    val = _H(PZ) + _H(PM) + _H(PK)
    for p, w in ((PZ, z), (PK, k), (PM, m)):
        val = np.where(p > 0, np.where(w > 0, val + p*np.log2(max(w, 1)), -np.inf), val)
    val = np.where(OK, val, -np.inf)
    out = []
    for t in taus:
        msk = (PZ <= t + 1e-9) & (t <= PZ + PM + 1e-9)
        out.append(np.max(np.where(msk, val, -np.inf)))
    return np.array(out)

TAUS = np.linspace(0, 1, 201)
def rho_two_piece(rcls, ccls, lam):
    s = f_curve(rcls, TAUS) + f_curve(ccls, TAUS)
    return float(np.max(s)) / (4*lam), float(TAUS[np.argmax(s)])

def rho_of(R, kappa, lam):
    rc, d1 = classes(rowmat(R), kappa); cc, d2 = classes(colmat(R), kappa)
    r, tau = rho_two_piece(rc, cc, lam)
    return r, rc, cc, (d1 and d2), tau

# ---------- direct check at finite t: exact ranks of grouped pieces A_i^{(t)} ----------
def kron4(A, B):
    N1, N2 = A.shape[0], B.shape[0]
    return np.einsum('abcd,efgh->aebfcgdh', A, B).reshape((N1*N2,)*4)
def grouped_ranks(R, C, t):
    polys = [R, C]               # coefficient list in T
    for _ in range(t - 1):
        new = [None]*(len(polys)+1)
        for i, P in enumerate(polys):
            for j, Q in enumerate((R, C)):
                term = kron4(P, Q)
                new[i+j] = term if new[i+j] is None else new[i+j] + term
        polys = new
    return [(rank(rowmat(A)), rank(colmat(A))) for A in polys]
