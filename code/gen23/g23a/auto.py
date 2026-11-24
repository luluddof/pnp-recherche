# find group automorphisms of a given order (small groups as tables), by generator-image backtracking
import itertools, random
def order(mul, a):
    k, x = 1, a
    while x != 0: x = mul[x][a]; k += 1
    return k
def gens_of(mul):
    N = len(mul); gens = []; S = {0}
    def close(gs):
        S = {0}; fr = [0]
        while fr:
            nf = []
            for x in fr:
                for g in gs:
                    y = mul[x][g]
                    if y not in S: S.add(y); nf.append(y)
            fr = nf
        return S
    cur = {0}
    for a in sorted(range(1, N), key=lambda a: -order(mul, a)):
        if a not in cur: gens.append(a); cur = close(gens)
        if len(cur) == N: break
    return gens
def extend(mul, gens, imgs):
    """try to extend gens -> imgs to a homomorphism; return map list or None"""
    N = len(mul); phi = {0: 0}; fr = [0]
    while fr:
        nf = []
        for x in fr:
            for g, h in zip(gens, imgs):
                y = mul[x][g]; z = mul[phi[x]][h]
                if y in phi:
                    if phi[y] != z: return None
                else: phi[y] = z; nf.append(y)
        fr = nf
    if len(set(phi.values())) != N: return None
    # full homomorphism check
    for a in range(N):
        for b in range(N):
            if phi[mul[a][b]] != mul[phi[a]][phi[b]]: return None
    return [phi[x] for x in range(N)]
def automorphisms_of_order(mul, k, limit=5, seed=0):
    N = len(mul); gens = gens_of(mul); ords = [order(mul, g) for g in gens]
    cand = [[x for x in range(N) if order(mul, x) == o] for o in ords]
    out = []; rnd = random.Random(seed)
    for c in cand: rnd.shuffle(c)
    def rec(i, imgs):
        if len(out) >= limit: return
        if i == len(gens):
            phi = extend(mul, gens, imgs)
            if phi is not None:
                # order of phi
                p = list(range(N)); o = 0
                while True:
                    p = [phi[x] for x in p]; o += 1
                    if p == list(range(N)): break
                if o == k: out.append(phi)
            return
        for x in cand[i]:
            # quick pairwise commutation check
            ok = True
            for j in range(i):
                if (mul[gens[i]][gens[j]] == mul[gens[j]][gens[i]]) != (mul[x][imgs[j]] == mul[imgs[j]][x]): ok = False; break
            if ok: rec(i + 1, imgs + [x])
            if len(out) >= limit: return
    rec(0, [])
    return out
if __name__ == "__main__":
    from groups import catalogue
    import time
    for g in ('ES-', 'ES+'):
        G = catalogue(32)[g]; t = time.time()
        for k in (5, 3):
            A = automorphisms_of_order(G[0], k, limit=1)
            print(g, 'order', k, 'found', len(A), f'{time.time()-t:.1f}s')
