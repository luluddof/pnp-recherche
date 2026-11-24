# Code-loop (Moufang, Griess) twisted group algebras on F_2^k from doubly-even codes; composition ansatz (G22c famA) -> p exactly.
import sys, itertools, numpy as np
sys.path.insert(0, '/home/user/pnp-recherche/code/gen22/g22c')
import ev, fams

def simplex_basis(k):          # generator rows of the [2^k-1, k] simplex code (weights 2^{k-1}: doubly even for k>=3)
    cols = [c for c in itertools.product([0,1], repeat=k) if any(c)]
    return np.array(cols).T     # k x n

def rm1(m):                    # RM(1,m) without all-ones? include all-ones: [2^m, m+1], weights 0, 2^{m-1}, 2^m
    pts = np.array(list(itertools.product([0,1], repeat=m)))
    rows = [pts[:, i] for i in range(m)] + [np.ones(2**m, int)]
    return np.array(rows)

def cocycle(Bs):
    k, n = Bs.shape
    w = Bs.sum(1); inter2 = Bs @ Bs.T
    tau = np.einsum('in,jn,ln->ijl', Bs, Bs, Bs) % 2
    el = np.array(list(itertools.product([0,1], repeat=k)))
    N = len(el)
    th = np.zeros((N, N), int)
    for A in range(N):
        a = el[A]
        for B in range(N):
            b = el[B]
            t = sum(a[i]*b[i]*(w[i]//4) for i in range(k))
            t += sum(a[i]*b[j]*(inter2[i, j]//2) for i in range(k) for j in range(k) if i > j)
            t += sum(b[l]*a[i]*a[j]*tau[i, j, l] for l in range(k) for i in range(k) for j in range(i+1, k))
            th[A, B] = t % 2
    # checks: theta(a,a)=|a|/4, theta(a,b)+theta(b,a)=|a&b|/2 , delta theta = |a&b&c|
    words = (el @ Bs) % 2
    ok1 = all(th[A, A] == (words[A].sum()//4) % 2 for A in range(N))
    ok2 = all((th[A, B] + th[B, A]) % 2 == ((words[A]*words[B]).sum()//2) % 2 for A in range(N) for B in range(N))
    idx = {tuple(e): i for i, e in enumerate(el)}
    add = lambda A, B: idx[tuple((el[A]+el[B]) % 2)]
    ok3 = all((th[B, C] + th[add(A, B), C] + th[A, add(B, C)] + th[A, B]) % 2 == (words[A]*words[B]*words[C]).sum() % 2
              for A in range(N) for B in range(N) for C in range(N)) if N <= 32 else None
    T = np.array([[add(A, B) for B in range(N)] for A in range(N)])
    return T, (-1) ** th, (ok1, ok2, ok3)

for name, Bs in [('simplex3', simplex_basis(3)), ('simplex4', simplex_basis(4)), ('rm1_4', rm1(4)), ('simplex5', simplex_basis(5)), ('rm1_5', rm1(5))]:
    T, f, ok = cocycle(Bs)
    N = T.shape[0]
    terms = fams.struct_slices(T, f)
    for d in (2, 1):
        r = ev.p_of(N, terms, d=d)
        print(name, 'N=', N, 'checks', ok, 'd', d, 'p=', r, 'q=', N, flush=True)
