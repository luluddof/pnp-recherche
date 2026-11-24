"""Piecewise-OpenAI per-point family (G7a), exact rational verification.
Weights w=|z| in {0..l} (z in {0,1}^l; s = l - 2w = sum of +-1 spins).  Partition {0..l} into blocks of
<= d+1 consecutive weights.  For block B choose p_B in Q[s], deg <= d, p_B(l-2w) = (-1)^w for w in B
(Lagrange).  For z with |z| in B:
     g_z(t) = chi(z) + p_B( sum_j (1 - t/tau_j) s_j(z) ) - p_B( sum_j s_j(z) ),   s_j = (-1)^{z_j}.
deg_t g_z <= d, g_z(0) = chi(z).  Edge (z,z^k) good iff g_z(tau_k) = g_{z^k}(tau_k).
Claim: bad edges = exactly the edges joining different blocks (checked here); hence
rank E_k = #{y in {0,1}^{l-1}: the edge from weight |y| to |y|+1 crosses a block boundary}."""
import itertools, sys
from fractions import Fraction as F
from math import comb
def lagrange_poly(xs, ys):
    def p(x):
        tot = F(0)
        for i, (a, b) in enumerate(zip(xs, ys)):
            term = F(b)
            for j, c in enumerate(xs):
                if j != i: term *= F(x - c, a - c)
            tot += term
        return tot
    return p
def blocks_for(l, d, r):
    """boundaries after weights r, r+d+1, ...: blocks [.., r], [r+1, r+d+1], ..."""
    B = []; start = 0; cut = r
    while start <= l:
        end = min(cut, l); B.append(list(range(start, end + 1))); start = end + 1; cut = end + d + 1
    return [b for b in B if b]
def check(l, d, blocks, taus):
    blk = {}
    ps = []
    for i, b in enumerate(blocks):
        assert len(b) <= d + 1
        for w in b: blk[w] = i
        ps.append(lagrange_poly([l - 2 * w for w in b], [(-1) ** w for w in b]))
    def g(z, t):
        sp = [(-1) ** ((z >> j) & 1) for j in range(l)]; w = sum((z >> j) & 1 for j in range(l)); p = ps[blk[w]]
        return F((-1) ** w) + p(sum((1 - F(t) / taus[j]) * sp[j] for j in range(l))) - p(sum(sp))
    per = [0] * l; mism = 0
    for z in range(1 << l):
        for k in range(l):
            if z >> k & 1: continue
            u = z | (1 << k); w = bin(z).count('1')
            bad = g(z, taus[k]) != g(u, taus[k])
            cross = blk[w] != blk[w + 1]
            per[k] += bad; mism += (bad and not cross)   # good-inside-block must hold; crossing edges may accidentally be good
    return per, mism
if __name__ == "__main__":
    for l in range(3, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 10):
        for d in range(1, min(l, 6)):
            best = None
            for r in range(d + 1):
                bl = blocks_for(l, d, r)
                pred = sum(comb(l - 1, w) for w in range(l) if any(w == b[-1] for b in bl[:-1]))
                if best is None or pred < best[0]: best = (pred, r, bl)
            pred, r, bl = best
            per, mism = check(l, d, bl, [F(k + 1) for k in range(l)])
            print(f"l={l} d={d} r={r}: predicted bad/dir={pred} actual per-dir={per} inside-block failures={mism} "
                  f"rho_max={max(per)/2**l:.4f}  (1/(2(d+1))={1/(2*(d+1)):.4f})", flush=True)
