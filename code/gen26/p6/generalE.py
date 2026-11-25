# Is "Smith(M(E)) <= square-free subset sums of Smith(E), termwise" true for ARBITRARY E over Z_2?  Random test.
import random, sys
src = open(sys.argv[1]).read()
exec(src.split("def sform")[0].split("def Eblock")[0])   # v2, smith_vals, MOD
exec("def Mmat" + src.split("def Mmat")[1].split("def sform")[0])
random.seed(7)
def randGL(n):
    while True:
        A = [[random.randint(-3, 3) for _ in range(n)] for _ in range(n)]
        # det odd check via mod-2 rank
        Bm = [[x & 1 for x in r] for r in A]; rk = 0; rows = Bm[:]
        for c in range(n):
            piv = next((r for r in range(rk, n) if rows[r][c]), None)
            if piv is None: continue
            rows[rk], rows[piv] = rows[piv], rows[rk]
            for r in range(n):
                if r != rk and rows[r][c]: rows[r] = [(x ^ y) for x, y in zip(rows[r], rows[rk])]
            rk += 1
        if rk == n: return A
def matmul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(len(B))) % MOD for j in range(len(B[0]))] for i in range(len(A))]
viol = 0; tot = 0; lowviol = 0
for trial in range(300):
    n = random.choice([2, 3, 4, 5])
    s = sorted(random.randint(0, 5) for _ in range(n))
    P, Q = randGL(n), randGL(n)
    E = matmul(matmul(P, [[(1 << s[i]) if i == j else 0 for j in range(n)] for i in range(n)]), Q)
    E = [[x % MOD for x in r] for r in E]
    sm = smith_vals(Mmat(E), 200)
    sq = [0]
    for w in s: sq = sq + [x + w for x in sq]
    sq.sort()
    full = [0]
    cap = max(sm) + 1 if sm else 1
    for w in s: full = full + [x + e * w for x in full for e in range(1, 60) if x + e * w < cap]
    full = sorted(full)[:len(sm)]
    tot += 1
    if len(sm) < len(sq) or any(u > v for u, v in zip(sm, sq)):
        viol += 1
        if viol <= 5: print("upper-side violation (Smith > sqfree):", n, s, sm, sq)
    if any(u < v for u, v in zip(sm, full)): lowviol += 1; print("Smith < full monomial?", s, sm)
print("trials", tot, "violations of Smith<=sqfree:", viol, " violations of Smith>=full:", lowviol)
