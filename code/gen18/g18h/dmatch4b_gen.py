ham = [[bin(a ^ b).count("1") for b in range(16)] for a in range(16)]
def gen(pts, budget):
    if not pts:
        yield (); return
    a = pts[0]; k = len(pts) // 2
    for b in pts[1:]:
        d = ham[a][b]
        if d + (k - 1) <= budget:
            rest = [q for q in pts if q != a and q != b]
            for m in gen(rest, budget - d):
                yield ((a, b),) + m
