# Convert the exact solutions of G19b/G20b to the two-palette array and check the structure theorem on them.  [E]
import numpy as np, networkx as nx, itertools
from anat import load, Rmat, files
from blocks import blocks
from shape import verify_shape, show
def to_cells(N, G, Gp):
    q = len(G)
    if not all((G[l] == Gp[l].T).all() for l in range(q)): return None
    cells = [[None]*N for _ in range(N)]
    for l in range(q):
        for a in range(N):
            for b in range(N):
                if G[l, a, b]:
                    assert cells[a][b] is None; cells[a][b] = ('S', l, int(G[l, a, b]))
    R, C = Rmat(N, G, Gp)
    for j, (rows, cols, sub) in enumerate(blocks(N, R)):
        if np.linalg.matrix_rank(sub) != 1 or not (sub != 0).all(): return None
        if sorted((a % N) * N + a // N for a in rows) != cols: return None
        s0 = sub[:, [cols.index((rows[0] % N) * N + rows[0] // N)]].ravel()   # column of the transposed first row cell: = h * h[0]
        for i, a in enumerate(rows):
            assert cells[a // N][a % N] is None; cells[a // N][a % N] = ('T', q + j, int(s0[i]))
    assert all(cells[a][b] is not None for a in range(N) for b in range(N))
    return cells
def intercal_graph(N, cells):
    g = nx.Graph(); bad = 0
    for a in range(N):
        for b in range(N): g.add_node(cells[a][b][:2])
    for a, c in itertools.combinations(range(N), 2):
        for b in range(N):
            for d in range(N):
                if b != d and cells[a][b][:2] == cells[c][d][:2]:
                    if cells[a][d][:2] != cells[c][b][:2]: bad += 1
                    elif cells[a][b][2]*cells[c][d][2]*cells[a][d][2]*cells[c][b][2] != -1: bad += 1
                    else: g.add_edge(cells[a][b][:2], cells[a][d][:2])
    return g, bad
if __name__ == "__main__":
    for f in files:
        N, G, Gp = load(f); cells = to_cells(N, G, Gp)
        print("=====", f.split('/')[-1])
        if cells is None: print("  NOT of the strict shape (G' != G^T or R not a sum of rank-one blocks)"); continue
        show(N, cells); verify_shape(N, cells)
        g, bad = intercal_graph(N, cells)
        sizes = {}
        for a in range(N):
            for b in range(N): sizes[cells[a][b][:2]] = sizes.get(cells[a][b][:2], 0) + 1
        print("  violations of the signed intercalate condition:", bad, "| intercalation graph bipartite S|T:", all(u[0] != v[0] for u, v in g.edges()),
              "| components:", [len(c) for c in nx.connected_components(g)], "| edges", g.number_of_edges())
        print("  class sizes S:", sorted((v for k, v in sizes.items() if k[0] == 'S'), reverse=True), " T:", sorted((v for k, v in sizes.items() if k[0] == 'T'), reverse=True))
