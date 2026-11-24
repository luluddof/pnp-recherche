# groups of order q^3: Sylow 2-subgroups of SU(3,q) (q = 2: Q8; q = 4: order 64), and a few order-64 comparisons
from groups import from_elems, direct, dicyclic, abelian, catalogue
def gf(n, poly):
    def mul(a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            b >>= 1; a <<= 1
            if a >> n: a ^= poly
        return r
    return mul
def su3_sylow(q):
    n = {2: 2, 4: 4, 8: 6}[q]; poly = {2: 0b111, 4: 0b10011, 8: 0b1000011}[q]   # F_{q^2}
    m = gf(n, poly)
    def pw(a, k):
        r = 1
        for _ in range(k): r = m(r, a)
        return r
    Q2 = 1 << n
    elems = [(a, b) for a in range(Q2) for b in range(Q2) if b ^ pw(b, q) == pw(a, q + 1)]
    assert len(elems) == q ** 3, len(elems)
    return from_elems(elems, lambda x, y: (x[0] ^ y[0], x[1] ^ y[1] ^ m(x[0], pw(y[0], q))), (0, 0))
G64 = {
    'SU3(4)syl': lambda: su3_sylow(4),
    'Q8xQ8': lambda: direct(dicyclic(2), dicyclic(2)),
    'Q8xZ2^3': lambda: direct(dicyclic(2), abelian([2, 2, 2])),
    'Q8xZ8': lambda: direct(dicyclic(2), abelian([8])),
    'Q16xZ4': lambda: direct(dicyclic(4), abelian([4])),
    'Q8xQ8?': None,
}
del G64['Q8xQ8?']
if __name__ == "__main__":
    from groups import invariants
    print('q=2', invariants(su3_sylow(2)), invariants(dicyclic(2)))
    G = su3_sylow(4); print('q=4', invariants(G)[1:])
