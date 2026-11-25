# [E] Strassen's rank-7 decomposition of <2,2,2> read in P_3 language, and its Kronecker powers.
# P_3 parties: node i holds (a_i,b_i); a CP term = local product operator U1(x)U2(x)U3 (U_i: a_i -> b_i), ring cost 1,
# row cost rank(U1)rank(U2)rank(U3).  sigma[b,a] = [b1=a3][b2=a1][b3=a2].
import numpy as np, flint, itertools, sys
def rk(M): return flint.fmpz_mat([[int(v) for v in r] for r in M]).rank()
E=lambda i,j: np.array([[1 if (r,c)==(i,j) else 0 for c in range(2)] for r in range(2)])
# standard Strassen for C=AB (A,B,C 2x2, 0-based): M_r = (sum alpha A)(sum beta B), C += gamma M_r
A=lambda *t: sum(s*E(i,j) for s,i,j in t)
terms=[ (A((1,0,0),(1,1,1)), A((1,0,0),(1,1,1)), A((1,0,0),(1,1,1))),
        (A((1,1,0),(1,1,1)), A((1,0,0)),         A((1,1,0),(-1,1,1))),
        (A((1,0,0)),         A((1,0,1),(-1,1,1)),A((1,0,1),(1,1,1))),
        (A((1,1,1)),         A((1,1,0),(-1,0,0)),A((1,0,0),(1,1,0))),
        (A((1,0,0),(1,0,1)), A((1,1,1)),         A((-1,0,0),(1,0,1))),
        (A((1,1,0),(-1,0,0)),A((1,0,0),(1,0,1)), A((1,1,1))),
        (A((1,0,1),(-1,1,1)),A((1,1,0),(1,1,1)), A((1,0,0))) ]
# tensor T[A_ij,B_jk,C_ik] = sum_r alpha_r[ij] beta_r[jk] gamma_r[ik] must equal [matmul]
T=np.einsum('rij,rjk,rik->ijjkik'.replace('ijjkik','ijxkyz'),np.array([t[0] for t in terms]),np.array([t[1] for t in terms]),np.array([t[2] for t in terms])) if False else None
Tm=np.zeros((2,)*6,int)
for t in terms: Tm+=np.einsum('ab,cd,ef->abcdef',t[0],t[1],t[2])
ref=np.zeros((2,)*6,int)
for i,j,k in itertools.product(range(2),repeat=3): ref[i,j,j,k,i,k]=1
assert (Tm==ref).all(), "Strassen check failed"
# map to P_3: L=(a1,b1): A_ij with i=b1, j=a1; M=(a2,b2): B_jk, j=b2, k=a2; R=(a3,b3): C_ik with i=a3, k=b3
# local operator U_i[b_i,a_i]: U1 = alpha (rows i=b1, cols j=a1); U2[b2,a2] = beta[j=b2,k=a2]; U3[b3,a3]=gamma[i=a3,k=b3]^T
ops=[(t[0],t[1],t[2].T) for t in terms]
def op3(U):  # operator on (a1,a2,a3)->(b1,b2,b3), matrix [b,a]
    return np.kron(np.kron(U[0],U[1]),U[2])
sig2=np.zeros((8,8),int)
for a in itertools.product(range(2),repeat=3):
    b=(a[2],a[0],a[1]); sig2[(b[0]*2+b[1])*2+b[2],(a[0]*2+a[1])*2+a[2]]=1
S=sum(op3(U) for U in ops); assert (S==sig2).all(), "P_3 mapping failed"
print("Strassen terms as local operators: row ranks", [int(np.prod([np.linalg.matrix_rank(u) for u in U])) for U in ops])
print("term 1 == identity operator:", (op3(ops[0])==np.eye(8,dtype=int)).all(), "; rank(sigma - I) =", rk(sig2-np.eye(8,dtype=int)))
# Kronecker powers: sigma_{2^k} = sigma_2^{(x)k} (site i of the big system = (site i of each factor)) -> reorder
def kron_power_op(mats,k):
    # mats: list of k 8x8 operators on 3 qubit-sites; return operator on 3 sites of dim 2^k each (site-major order)
    M=mats[0]
    for m in mats[1:]: M=np.kron(M,m)
    # current index order: (b1_1,b2_1,b3_1,b1_2,...) ; reorder to (b1_1..b1_k, b2_..., b3_...)
    d=[2]*(3*k); Mt=M.reshape(d+d)
    perm=[3*f+s for s in range(3) for f in range(k)]
    Mt=Mt.transpose(perm+[3*k+p for p in perm]); return Mt.reshape(8**k,8**k)
for k in [1,2,3]:
    N=2**k; sig=kron_power_op([sig2]*k,k)
    # groups G_j = sum of Strassen^{(x)k} terms with exactly j non-identity factors (j = k - n1)
    T6=sig2-np.eye(8,dtype=int); I8=np.eye(8,dtype=int)
    Gs=[]
    for j in range(k+1):
        G=np.zeros((8**k,8**k),int)
        for S in itertools.combinations(range(k),j):
            G+=kron_power_op([T6 if f in S else I8 for f in range(k)],k)
        Gs.append(G)
    assert (sum(Gs)==sig).all()
    out=[]
    for m in range(k+1):
        C=sum(Gs[:m+1]); nterms=sum(len(list(itertools.combinations(range(k),j)))*6**j for j in range(m+1))
        out.append((m,nterms,rk(sig-C)))
    print(f"k={k} N={N} N^3={N**3}: C = sum_(j<=m) G_j -> (m, #CP terms of C, rank(sigma-C)):",out,
          " row ranks of G_j:",[rk(G) for G in Gs])
# flattening ranks (node i | other two) of T6 = sigma_2 - I, sigma_2, I, and of C_m at k=2,3
def flat(M,N):
    T=M.reshape((N,)*6)   # (b1,b2,b3,a1,a2,a3)
    return [rk(T.transpose(0,3,1,2,4,5).reshape(N*N,-1)), rk(T.transpose(1,4,0,2,3,5).reshape(N*N,-1)), rk(T.transpose(2,5,0,1,3,4).reshape(N*N,-1))]
print("flattening ranks node|rest: sigma_2",flat(sig2,2)," I",flat(np.eye(8,dtype=int),2)," T6",flat(sig2-np.eye(8,dtype=int),2))
for k in [2,3]:
    N=2**k; T6=sig2-np.eye(8,dtype=int); I8=np.eye(8,dtype=int); sig=kron_power_op([sig2]*k,k)
    Gs=[sum(kron_power_op([T6 if f in S else I8 for f in range(k)],k) for S in itertools.combinations(range(k),j)) for j in range(k+1)]
    print(f"k={k}: flattening ranks of C_m:",[flat(sum(Gs[:m+1]),N) for m in range(k+1)], " N^2 =",N*N)
