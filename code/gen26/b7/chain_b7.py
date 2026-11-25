# G26-B7: long balanced chains (p9 chain_exact model) with the new lemmas.
#  bases: clique sweep (exp. space, e0 = s(K)/K^2), clique D&C (poly space, e0 = d(K)/K^2),
#         row sweep level 1 with p9 Lemma RS (zs_p9) or new ordering lemma R3 (zs_R3),
#         poly-space row D&C level 1 (e1 = (rho+1/g)/3 + zs_DC, zs_DC := 2*zs_R3 [H constant]).
#  groupings: p2/p9 Hoeffding (zeta_old) or stratified (R2, zeta_new).
#  lt per level: greedy (maximise next mu) over a grid, and also grids of geometric schedules.
import math, sys, functools
import numpy as np
import b7lib as B
RHO = B.RHO

_DC = [0, 0]
def dC(K):  # poly-space D&C exponent (per lam) of a K-clique: min_j j(K-j) + max(d(j), d(K-j))
    while len(_DC) <= K:
        n = len(_DC); _DC.append(min(j*(n-j) + max(_DC[j], _DC[n-j]) for j in range(1, n//2+1)))
    return _DC[K]

def s_sweep(K): return ((K+1)**2)//4 - 1

def step(x, e, lt, zfun):
    z = float(zfun(lt, x)); g = 2.0**(x+lt) if x+lt < 1000 else float('inf')
    inv_g = 0 if g == float('inf') else 1/g; invK = 2.0**-x if x < 1000 else 0.0
    return 2*(x+lt), e*(RHO + inv_g + z) + invK

LTS = np.arange(1, 60.01, 0.25)
def chain_greedy(x, e, zfun, L=24, memo=False):
    out = [(x, e)]
    for l in range(L):
        best = None
        for lt in LTS:
            x2, e2 = step(x, e, lt, zfun)
            mu = 1/(e2*e2*x2)
            if best is None or mu > best[0]: best = (mu, x2, e2, lt)
        _, x, e, lt = best; out.append((x, e))
    return out

def chain_param(x, e, zfun, a0, b0, L=24):
    # lt_l = a0 + b0*log2(x_l): t ~ 2^a0 x^b0
    out = [(x, e)]
    for l in range(L):
        lt = max(1.0, a0 + b0*math.log2(x)); x, e = step(x, e, lt, zfun); out.append((x, e))
    return out

def summarize(ch, lo=16):
    mu = 1/(ch[-1][1]**2*ch[-1][0])
    s2 = max(2*ch[l][1]*math.sqrt(2*ch[l+1][0]) for l in range(lo, len(ch)-1))
    return mu, s2

def starts(kind):
    S = []
    if kind == 'clique-sweep':
        for x0 in range(3, 16): K = 2**x0; S.append((f"x0={x0}", float(x0), s_sweep(K)/K**2))
    elif kind == 'clique-DC':
        for x0 in range(3, 12): K = 2**x0; S.append((f"x0={x0}", float(x0), dC(K)/K**2))
    elif kind in ('row-p9', 'row-R3', 'row-DC'):
        for lg in np.arange(6, 40.01, 0.5):
            if kind == 'row-p9': e1 = (RHO + 2.0**-lg)/4 + B.zs_p9(lg)
            elif kind == 'row-R3': e1 = (RHO + 2.0**-lg)/4 + B.zs_R3(lg)
            else: e1 = (RHO + 2.0**-lg)/3 + 2*B.zs_R3(lg)
            S.append((f"g=2^{lg}", 2*lg, e1))
    return S

def best_chain(kind, zname):
    zfun = B.zeta_new if zname == 'new' else (lambda lt, x: B.zeta_old(lt, x))
    best = None
    for name, x, e in starts(kind):
        cands = [chain_greedy(x, e, zfun)]
        for a0 in np.arange(-4, 16.01, 1.0):
            for b0 in (0.0, 0.5, 1.0, 1.5, 2.0):
                cands.append(chain_param(x, e, zfun, a0, b0))
        for ch in cands:
            mu, s2 = summarize(ch)
            if best is None or mu > best[0]: best = (mu, s2, name, ch)
    return best

if __name__ == "__main__":
    kinds = sys.argv[1:] or ['clique-sweep', 'row-p9', 'row-R3', 'clique-DC', 'row-DC']
    print("d(K)/K^2 (poly-space D&C):", [(2**x, round(dC(2**x)/4**x, 4)) for x in range(1, 12)])
    for kind in kinds:
        for zname in ('old', 'new'):
            mu, s2, name, ch = best_chain(kind, zname)
            print(f"{kind:13s} zeta={zname}: mu_inf={mu:.4f} ({name}); single-schedule c1(a=2)={s2:.3f}; "
                  f"no-sawtooth c1 a=2 {2*math.sqrt(2/mu):.3f}, a=4 {2*math.sqrt(4/mu):.3f}")
            print("     levels (x, e, mu):", [(round(x, 1), round(e, 5), round(1/(e*e*x), 3)) for x, e in ch[:8]], flush=True)
