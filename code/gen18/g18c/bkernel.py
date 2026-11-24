"""G18c: OpenAI App. D Cauchy kernel on ONE packet (L = 2^a slots, one-bit wires), read as a family of
Claude-type pieces (substitution (ii)): piece (t, j) = [z^j] K(t + z), K(T)[Y,X] = prod_x l_{Y_x}((W(T)X)_x).
Exact arithmetic mod a 61-bit prime (ranks over GF(P): w.h.p. = ranks over Q, never larger).
Checks: K(0) = I; K(t) has rank 1 in direction v2(t); rank profile of the jets in every direction;
exact per-packet Hermite identity I = sum_{t,j} c_{t,j} piece(t,j) for a chosen multiplicity profile."""
import sys, itertools, math
import flint
P = (1 << 61) - 1
a = int(sys.argv[1]); L = 1 << a
S = list(range(-L, L))
def poly(c): return flint.nmod_poly(c, P)
T = poly([0, 1])
def Cwithout(skip):
    f = poly([1])
    for t in S:
        if t != skip: f = f * poly([(-t) % P, 1])
    return f
Cp0 = 1
for t in S:
    if t != 0: Cp0 = Cp0 * (-t)
inv = pow(Cp0 % P, P - 2, P)
W = [[Cwithout(y - x) * inv for y in range(L)] for x in range(L)]   # W_xy(T), degree 2L-1
deg = L * (2 * L - 1)
bits = list(itertools.product([0, 1], repeat=L))
# K[Y][X] as polynomials
rowform = {}
for xi, X in enumerate(bits):
    forms = []
    for x in range(L):
        f = poly([0])
        for y in range(L):
            if X[y]: f = f + W[x][y]
        forms.append(f)
    rowform[xi] = forms
K = [[None] * len(bits) for _ in bits]
one = poly([1])
for xi, X in enumerate(bits):
    forms = rowform[xi]
    for yi, Y in enumerate(bits):
        f = one
        for x in range(L):
            f = f * (forms[x] if Y[x] else (one - forms[x]))
        K[yi][xi] = f
ok0 = all(int(K[yi][xi](0)) == (1 if yi == xi else 0) for yi in range(len(bits)) for xi in range(len(bits)))
print(f"a={a} L={L} deg K <= {deg}; K(0)=I: {ok0}; max deg found {max(K[y][x].degree() for y in range(len(bits)) for x in range(len(bits)))}")
def flatten(mat, i):
    """mat[yi][xi] scalars. direction i: group of row port x = bit i of x, of column port y = 1 - bit i of y.
    rows of the flattening: bits of group 0, columns: bits of group 1."""
    g0 = [('Y', x) for x in range(L) if not (x >> i) & 1] + [('X', y) for y in range(L) if (y >> i) & 1]
    g1 = [('Y', x) for x in range(L) if (x >> i) & 1] + [('X', y) for y in range(L) if not (y >> i) & 1]
    n0, n1 = len(g0), len(g1)
    M = [[0] * (1 << n1) for _ in range(1 << n0)]
    for yi, Y in enumerate(bits):
        for xi, X in enumerate(bits):
            v = mat[yi][xi]
            if v == 0: continue
            r = 0
            for k, (s, p) in enumerate(g0): r |= ((Y[p] if s == 'Y' else X[p]) << k)
            c = 0
            for k, (s, p) in enumerate(g1): c |= ((Y[p] if s == 'Y' else X[p]) << k)
            M[r][c] = v
    return flint.nmod_mat(M, P)
def v2(t):
    t = abs(t); k = 0
    while t % 2 == 0: t //= 2; k += 1
    return k
print("identity ranks per direction:", [flatten([[1 if y == x else 0 for x in range(len(bits))] for y in range(len(bits))], i).rank() for i in range(a)])
jmax = int(sys.argv[2]) if len(sys.argv) > 2 else deg
prof = {}
for t in S:
    if t == 0: continue
    sh = poly([t % P, 1])
    Kt = [[K[y][x](sh) for x in range(len(bits))] for y in range(len(bits))]   # K(t+z)
    i0 = v2(t)
    line = []
    for j in range(0, jmax + 1):
        mat = [[int(Kt[y][x][j]) for x in range(len(bits))] for y in range(len(bits))]
        rk = [flatten(mat, i).rank() for i in range(a)]
        prof[(t, j)] = rk
        line.append(rk)
    print(f"t={t:+d} v2={i0}: jets j=0..{jmax}: ranks per direction (dir0..dir{a-1}) =", line)
# per-packet Claude accounting: for each node t (valuation < a), jets j < q; cost of piece = sum of log2 ranks over 2 best?? report dir-cost and total
print("\nsummary (log2 rank in own direction / L) by jet order, averaged over nodes of each valuation:")
for i in range(a):
    nodes = [t for t in S if t != 0 and v2(t) == i]
    row = []
    for j in range(0, jmax + 1):
        row.append(sum(math.log2(max(1, prof[(t, j)][i])) for t in nodes) / len(nodes) / L)
    print(f" v2={i} ({len(nodes)} nodes):", " ".join(f"{x:.2f}" for x in row))
if a == 2:
    print("\nClaude accounting (two directions): rho(t,j) = (log2 r0 + log2 r1)/(2L)  [crossing: 2L bits per direction]")
    for t in S:
        if t == 0 or v2(t) >= a: continue
        print(f" t={t:+d}:", " ".join(f"{(math.log2(max(1,prof[(t,j)][0]))+math.log2(max(1,prof[(t,j)][1])))/(2*L):.2f}" for j in range(jmax + 1)))
