# Symplectic maps of the M-picture (2k qubits): vector (a1,a2,zx,zy) in F2^{4k}; omega((a,z),(a',z')) = a.z' + z.a'.
# F_M: a=(a1,a2) -> (a2, a1+a2), z=(zx,zy) -> (zx+zy, zx).   F_12 (M12|x,y>=|x+y,y>): a -> (a1+a2,a2), z -> (zx, zx+zy).
# F_23 (M23|x,y>=|x,x+y>): a -> (a1, a1+a2), z -> (zx+zy, zy).
# Question: is F_12 a product of m conjugates of F_M in Sp(4k,2)?  (=> Pauli-sparse rigidity of M controls that of M12)
import numpy as np, sys, itertools
k = int(sys.argv[1]); n = 4*k; rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
I = np.eye(k, dtype=np.int64); Z = np.zeros((k, k), dtype=np.int64)
def blk(rows): return np.block(rows) % 2
FM = blk([[Z, I, Z, Z], [I, I, Z, Z], [Z, Z, I, I], [Z, Z, I, Z]])
F12 = blk([[I, I, Z, Z], [Z, I, Z, Z], [Z, Z, I, Z], [Z, Z, I, I]])
F23 = blk([[I, Z, Z, Z], [I, I, Z, Z], [Z, Z, I, I], [Z, Z, Z, I]])
# column convention: new = F @ old
J = np.block([[np.zeros((2*k, 2*k), dtype=np.int64), np.eye(2*k, dtype=np.int64)], [np.eye(2*k, dtype=np.int64), np.zeros((2*k, 2*k), dtype=np.int64)]])
def issymp(F): return np.array_equal((F.T @ J @ F) % 2, J)
assert issymp(FM) and issymp(F12) and issymp(F23)
assert np.array_equal((F12 @ F23) % 2, FM), 'FM = F12 F23'
Id = np.eye(n, dtype=np.int64)
def inC(g):  # g^2 + g + 1 = 0  (fixed-point-free order 3)
    return not np.any((g @ g + g + Id) % 2)
def inv(g):  # symplectic inverse: J g^T J
    return (J @ g.T @ J) % 2
def transvection(v):  # x -> x + omega(x,v) v
    return (Id + np.outer(v, v @ J)) % 2
def rand_symp(steps=60):
    g = Id.copy()
    for _ in range(steps):
        v = rng.integers(0, 2, n)
        if v.any(): g = (transvection(v) @ g) % 2
    return g
assert inC(FM)
if k == 1:
    # enumerate Sp(4,2) by BFS over transvections
    gens = [transvection(np.array(v)) for v in itertools.product((0, 1), repeat=n) if any(v)]
    key = lambda g: g.tobytes(); seen = {key(Id): Id}; frontier = [Id]
    while frontier:
        nf = []
        for g in frontier:
            for t in gens:
                h = (t @ g) % 2
                if key(h) not in seen: seen[key(h)] = h; nf.append(h)
        frontier = nf
    G = list(seen.values()); print('|Sp(4,2)| =', len(G))
    C = {key((g @ FM @ inv(g)) % 2): (g @ FM @ inv(g)) % 2 for g in G}
    allC = [g for g in G if inC(g)]
    print('class of F_M:', len(C), ' #elements with g^2+g+1=0:', len(allC), '(equal => unique class)')
    for m in (1, 2, 3):
        prods = {key(Id): Id}
        for _ in range(m):
            prods = {key((p @ c) % 2): (p @ c) % 2 for p in prods.values() for c in C.values()}
        print(f'm={m}: |C^m| = {len(prods)};  F12 in C^m: {key(F12) in prods}')
else:
    # sampling: random A in C (conjugate of FM); test whether A^{-1} F12 in C (class assumed = {g^2+g+1=0})
    hits = 0; T = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
    for t in range(T):
        g = rand_symp(); A = (g @ FM @ inv(g)) % 2
        if inC((inv(A) @ F12) % 2):
            hits += 1
            if hits == 1:
                B = (inv(A) @ F12) % 2
                assert np.array_equal((A @ B) % 2, F12) and inC(A) and inC(B) and issymp(A) and issymp(B)
                print('found F12 = A B with A,B fixed-point-free order 3 (both in C if class unique)')
    print(f'k={k}: {hits}/{T} random A in C have A^-1 F12 in C')
