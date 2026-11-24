# Hurwitz-Radon block construction of an integer formula [N,N,m] as a signed intercalate matrix, then (p,q) in the S/T/U model:
# per connected component of the intercalation graph: bipartite -> (|A|,|B|) as T/S; else U (counts in both).  Exact verify.
import numpy as np, itertools, sys, networkx as nx
from shape import verify_shape
from convert import intercal_graph
X = np.array([[0,1],[1,0]]); Z = np.array([[1,0],[0,-1]]); I2 = np.eye(2, dtype=int); Y = X @ Z   # Y = XZ real, skew
P = {'I': I2, 'X': X, 'Z': Z, 'Y': Y}
def mat(s):
    M = np.array([[1]])
    for ch in s: M = np.kron(M, P[ch])
    return M
def hr_family(k):
    """max set of mutually anticommuting skew real Pauli strings on k qubits (signed perms), plus identity first."""
    strs = [''.join(t) for t in itertools.product('IXZY', repeat=k)]
    sk = [s for s in strs if s.count('Y') % 2 == 1]
    def anti(s, t): return sum(1 for a, b in zip(s, t) if a != 'I' and b != 'I' and a != b) % 2 == 1
    g = nx.Graph(); g.add_nodes_from(sk)
    for s, t in itertools.combinations(sk, 2):
        if anti(s, t): g.add_edge(s, t)
    best = max(nx.find_cliques(g), key=len)
    return [np.eye(2**k, dtype=int)] + [mat(s) for s in best]
def block(r, k, colour0):
    """r x n array (n = 2^k) from z = sum_i x_i A_i y : cell (i,j) gets colour of the nonzero row of column j of A_i."""
    A = hr_family(k)[:r]; n = 2**k
    out = [[None]*n for _ in range(r)]
    for i in range(r):
        for j in range(n):
            kk = int(np.nonzero(A[i][:, j])[0][0]); out[i][j] = (colour0 + kk, int(A[i][kk, j]))
    return out
def build(N, k):
    n = 2**k; rho = len(hr_family(k))
    rows_blocks = []; r0 = 0
    while r0 < N: rows_blocks.append((r0, min(rho, N - r0))); r0 += rho
    arr = [[None]*N for _ in range(N)]; c0 = 0
    for (r0, r) in rows_blocks:
        for cb in range(0, N, n):
            B = block(r, k, c0); c0 += n
            for i in range(r):
                for j in range(min(n, N - cb)): arr[r0 + i][cb + j] = B[i][j]
    return arr, rho
def palettes(N, arr):
    cells = [[('U', arr[a][b][0], arr[a][b][1]) for b in range(N)] for a in range(N)]
    g, bad = intercal_graph(N, cells); assert bad == 0
    used = set(arr[a][b][0] for a in range(N) for b in range(N))
    out = {}; S = T = U = 0
    for comp in nx.connected_components(g):
        sub = g.subgraph(comp)
        if nx.is_bipartite(sub):
            A, B = nx.bipartite.sets(sub) if len(comp) > 1 else (set(comp), set())
            for c in A: out[c[1]] = 'T'
            for c in B: out[c[1]] = 'S'
            T += len(A); S += len(B)
        else:
            for c in comp: out[c[1]] = 'U'
            U += len(comp)
    return out, S, T, U, g
if __name__ == "__main__":
    N, k = int(sys.argv[1]), int(sys.argv[2])
    arr, rho = build(N, k)
    m = len(set(arr[a][b][0] for a in range(N) for b in range(N)))
    pal, S, T, U, g = palettes(N, arr)
    print(f"N={N} HR blocks rho({2**k})={rho}: integer formula [N,N,{m}]; components {len(list(nx.connected_components(g)))}, S,T,U = {S},{T},{U} -> (p,q) = ({T+U},{S+U})")
    cells = [[(pal[arr[a][b][0]], arr[a][b][0], arr[a][b][1]) for b in range(N)] for a in range(N)]
    if N <= 32: verify_shape(N, cells)
    # all-U reading (always valid): (m,m)
    cellsU = [[('U', arr[a][b][0], arr[a][b][1]) for b in range(N)] for a in range(N)]
    if N <= 32: verify_shape(N, cellsU)
