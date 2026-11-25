# Route H: counted App-B (GTC) curve for the untwisted symmetric orbifold Sym^M(Heisenberg) = (H^{(x)M})^{S_M}.
# dims d_M(w) = # multisets of <= M non-vacuum Heisenberg states (p(n) states at weight n) of total weight w.
# Counting as in gen28/c1/appb_curve.py with p -> d_M.  Status: [C] (counted; identity = VOA locality + dual bases).
import math, sys
from math import comb
sys.path.insert(0, '/home/user/pnp-recherche/code/gen28/c1')
def parts(M):
    p=[1]+[0]*M
    for n in range(1,M+1):
        k=1; s=0
        while True:
            g1=k*(3*k-1)//2
            if g1>n: break
            sg = 1 if k%2 else -1
            s+=sg*p[n-g1]
            g2=k*(3*k+1)//2
            if g2<=n: s+=sg*p[n-g2]
            k+=1
        p[n]=s
    return p
def symdims(W, M):
    p = parts(W)
    dp = [[0]*(W+1) for _ in range(M+1)]; dp[0][0] = 1
    for n in range(1, W+1):
        c = p[n]
        new = [row[:] for row in dp]
        for j in range(M+1):
            for w in range(W+1):
                v = dp[j][w]
                if not v: continue
                r = 1
                while j + r <= M and w + n*r <= W:
                    new[j+r][w+n*r] += v * comb(c + r - 1, r); r += 1
        dp = new
    return [sum(dp[j][w] for j in range(M+1)) for w in range(W+1)]
def curve(d, h):
    N = d[h]; L = math.log(N)
    P = lambda w: d[w] if 0 <= w <= 4*h else 0
    rowc = {s: sum(P(s+i) for i in range(2*h+1)) for s in range(-2*h, 2*h+1)}
    colc = {s: sum(P(2*h-s-i) for i in range(2*h+1)) for s in range(-2*h, 2*h+1)}
    R = 0; Rs = {}
    for s0 in range(-2*h, 2*h+2):
        Rs[s0] = R
        if s0 <= 2*h: R += rowc[s0]
    C = 0; Cs = {}
    for s0 in range(2*h+1, -2*h-1, -1):
        if s0 <= 2*h: C += colc[s0]
        Cs[s0] = C
    Cs[2*h+1] = 0
    pts = []
    for s0 in range(-2*h, 2*h+2):
        th = math.log(Rs[s0])/L if Rs[s0] > 0 else 0.0
        T = math.log(Cs[s0])/L if Cs[s0] > 0 else 0.0
        pts.append((s0, th, T))
    th_end = math.log(Rs[2*h+1])/L
    # jump: smallest T among points with theta < 2
    jT = min((T for s0, th, T in pts if th < 2), default=None)
    # balanced
    bal = min(max(min(th, 2), T) for s0, th, T in pts)
    return N, th_end, jT, bal, pts
if __name__ == "__main__":
  for M in [1, 2, 3, 4, 6, 10, 1000]:
    for h in [6, 10, 16, 24]:
        d = symdims(4*h+2, min(M, 4*h+2))
        N, te, jT, bal, pts = curve(d, h)
        print(f"M={M:5d} h={h:3d} N=d(h)={N:.3e} log2N={math.log2(N):6.1f}  theta_end={te:.3f}  jump T*={jT:.3f}  balanced={bal:.3f}", flush=True)
