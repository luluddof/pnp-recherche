"""G14a: check the refined bound mechanism on the dumped char-2 solutions:
S = {a in A_L : a(tau_k)|_{x_k=0} B01_k = 0 all k}, L = l-1-d.  Claim: dim S <= sum_{j<L/2} (ball(L-j)-ball(j)),
hence sum_k rank B01_k >= N_L - that = N_m+N_{m-1} (L=2m) or 2N_m (L=2m+1).  Same for B10 (complemented)."""
import json, sys, glob
from math import comb
from gen_char2 import GF
def ball(l, r): return sum(comb(l, i) for i in range(r + 1)) if r >= 0 else 0
def check(fn, l, d, e, tau):
    F = GF(e); J = json.load(open(fn)); n = 1 << l; L = l - 1 - d
    def pw(a, i):
        r = 1
        for _ in range(i): r = F.mul(r, a)
        return r
    res = []
    for side in (0, 1):          # side 0: a supported on x_k=0 (B01); side 1: mirrored (B10)
        wt = lambda x: bin(x).count('1') if side == 0 else l - bin(x).count('1')
        params = [(x, i) for x in range(n) for i in range(L + 1 - wt(x))]
        rows = []; R = 0
        for (x, i) in params:
            row = []
            for k in range(l):
                tk = tau[k]
                pref = 1
                for kk in range(l):
                    if (x >> kk & 1) == (1 - side): pref = F.mul(pref, 1 ^ F.mul(F.inv[tau[kk]], tk))
                ax = F.mul(pref, pw(tk, i))        # value of a_x(tau_k) for this basis parameter
                for y in range(n):
                    if (y >> k & 1) == (1 - side):
                        row.append(F.mul(ax, F.polyeval(J[x][y], tk)) if (x >> k & 1) == side else 0)
            rows.append(row)
        for k in range(l):
            tk = tau[k]
            B = [[F.polyeval(J[x][y], tk) for y in range(n) if (y >> k & 1) == 1 - side] for x in range(n) if (x >> k & 1) == side]
            R += F.rank(B)
        dimS = len(params) - F.rank(rows)
        cap = sum(ball(l, L - j) - ball(l, j) for j in range(L + 1) if 2 * j < L)
        res.append((R, len(params), dimS, cap, len(params) - cap))
        assert dimS <= cap and R >= len(params) - cap
    return res
for fn, l, d, e, tau in [("sol_3_1_gf4_1-1-1.json", 3, 1, 2, [1, 2, 3]), ("sol_3_1_gf4_1-1-2.json", 3, 1, 2, [1, 2, 3])] + \
        [(f, 4, int(f.split('_')[2]), 3, None) for f in sorted(glob.glob("sol_4_*_gf8_*.json"))]:
    if tau is None: continue
    print(fn, "[side: (sum rank B, dim A_L, dim S, cap, bound)]", check(fn, l, d, e, tau))
