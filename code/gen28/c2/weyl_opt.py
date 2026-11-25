# exact optimum of the odd-d Weyl-diagonal model: build K in the common eigenbasis matching SWAP on 2-cycles optimally
import numpy as np
exec(open('weyl_check.py').read().split('for d, k in')[0])
for d,k in [(3,1),(3,2),(5,1)]:
    N=d**k; S=swap(N)
    # Bell basis Phi_v = (W_v (x) 1)|Phi>
    phi = np.eye(N).ravel()/np.sqrt(N)
    B=[]
    for lab in itertools.product(range(d), repeat=2*k):
        cs,es=lab[:k],lab[k:]; W=kl([weyl(d,c,e) for c,e in zip(cs,es)])
        B.append(np.kron(W,np.eye(N))@phi)
    B=np.array(B).T  # columns
    Sb = B.conj().T @ S @ B
    perm = np.abs(Sb).argmax(axis=0); fixed = sum(perm[i]==i for i in range(N*N))
    # optimal diagonal: fixed v: phi=Sb[v,v]; 2-cycle (v,w): phi_v = Sb[v,w], phi_w = Sb[w,v]  -> block [[-a,s],[s',-b]] rank 1
    diag = np.zeros(N*N, complex)
    for v in range(N*N):
        w = perm[v]
        diag[v] = Sb[v,v] if w==v else Sb[w,v]
    K = B @ np.diag(diag) @ B.conj().T
    print(d,k,"fixed points",fixed,"(N =",N,") p_opt =",rank(S-K),"(N^2-N)/2 =",(N*N-N)//2, "is SWAP monomial in Bell basis:", np.allclose(np.abs(Sb).sum(axis=0),1))
