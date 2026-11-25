# [E] Composite construction: P_3 = S12*S23. If Swap = R + K (rank R = p, osr K = q) then
#   P_3 - K12*K23 = R12*S23 + K12*R23  has rank <= 2pN, and K12*K23 is a ring (MPO) with bonds (q,q,1).
# Uses exact SWAP decompositions from code/gen26/p11 (read only). Operators as matrices M[out,in], site order (1,2,3).
import numpy as np, flint, json, sys, itertools
def rk(M):
    Z=flint.fmpz_mat([[int(v) for v in r] for r in M]); return Z.rank()
def load(fn):
    d=json.load(open(fn)); N=d['N']; E=np.array(d['E']); C=np.array(d['C']); return N,d['rank'],d['osr'],E,C
def ops(N,E,C):
    # K[x,y,u,v]=sum_l E_l[x,u] C_l[y,v]; inputs (x,y), outputs (u,v): operator O[(u,v),(x,y)] = sum_l E_l[x,u] C_l[y,v]
    K2=np.einsum('lxu,lyv->uvxy',E,C).reshape(N*N,N*N)
    S2=np.zeros((N,N,N,N),int)
    for x in range(N):
        for y in range(N): S2[y,x,x,y]=1     # out (u,v)=(y,x)
    S2=S2.reshape(N*N,N*N); I=np.eye(N,dtype=int)
    return K2,S2,I
def run(fn):
    N,p,q,E,C=load(fn); K2,S2,I=ops(N,E,C)
    assert rk(S2-K2)==p
    K12=np.kron(K2,I); K23=np.kron(I,K2); S12=np.kron(S2,I); S23=np.kron(I,S2)
    P3=S12@S23
    # check P3 = cyclic shift b1=a3,b2=a1,b3=a2
    Q=np.zeros((N,)*6,int)
    for a in itertools.product(range(N),repeat=3): Q[a[2],a[0],a[1],a[0],a[1],a[2]]=1
    assert (Q.reshape(N**3,N**3)==P3).all()
    Cring=K12@K23
    # ring bonds: check osr across site1|rest and site3|rest
    T=Cring.reshape((N,)*6)            # (b1,b2,b3,a1,a2,a3)
    f1=T.transpose(0,3,1,2,4,5).reshape(N*N,-1); f3=T.transpose(2,5,0,1,3,4).reshape(N*N,-1)
    pP=rk(P3-Cring)
    print(f"N={N} SWAP(p={p},q={q}) -> P3: rank(P3-C)={pP}  (bound 2pN={2*p*N}); ring bonds (12,23,31)=({q},{q},1): osr site1|23={rk(f1)}, site3|12={rk(f3)}; "
          f"X1 lower bound N^3/(4*1*q)={N**3/(4*q):.1f}; C=I gives {(2*(N**3-N))//3}; trivial ring N^3={N**3}")
for fn in sys.argv[1:]: run(fn)
