import itertools
from grpans import groups, table, solve
N = 8; Gp = groups(8)['Q8']; mul, inv = table(Gp)
# name elements: find -1 (central involution), i, j, k
e = 0; m1 = next(x for x in range(8) if x and mul[x][x] == 0)
names = {0: '1', m1: '-1'}
i_ = next(x for x in range(8) if x not in names); names[i_] = 'i'; names[mul[m1][i_]] = '-i'
j_ = next(x for x in range(8) if x not in names); names[j_] = 'j'; names[mul[m1][j_]] = '-j'
k_ = mul[i_][j_]; names[k_] = 'k'; names[mul[m1][k_]] = '-k'
S, F, B = solve(mul, inv, 1, 1, 4); sols = 0; seen = set()
while S.solve() and sols < 400:
    m = set(l for l in S.get_model() if l > 0)
    f = [next(v for v in range(N) if F[c][v] in m) for c in range(N)]; bl = [next(i for i in range(2) if B[c][i] in m) for c in range(N)]
    key = tuple(sorted(names[c] for c in range(N) if bl[c] == 0))
    if key not in seen and sols < 400:
        seen.add(key); print("G-block", key, " f:", {names[c]: names[f[c]] for c in range(N)})
    sols += 1
    S.add_clause([-F[c][f[c]] for c in range(N)] + [-B[c][bl[c]] for c in range(N)])
print("colourings (f,B) total (capped 400):", sols, " distinct G-blocks:", len(seen))
