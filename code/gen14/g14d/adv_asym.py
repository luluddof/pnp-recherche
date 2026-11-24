# G14d: (B) adversary optimisation over placements, (C) matrix check on multi-point layered J, (D) asymptotics.
import itertools, math, sys
from fractions import Fraction as Fr
from math import comb, log, sqrt, exp
from mp_check import bound_old, bound_new

def compositions(N, l):
    if l == 1: yield (N,); return
    for a in range(N+1):
        for r in compositions(N-a, l-1): yield (a,)+r

print("== B. minimum over ALL placements (m_c) of N points on l directions: attained at the balanced placement?")
ok = True
for l in range(1, 6):
    for N in range(1, 13):
        for d in range(0, N):
            bal = tuple(N//l + (1 if c < N % l else 0) for c in range(l))
            for f in (bound_old, bound_new):
                mn = min(f(l, N, d, mc) for mc in compositions(N, l))
                if mn != f(l, N, d, bal): ok = False; print("FAIL", l, N, d, f.__name__, mn, f(l, N, d, bal))
            # Jensen relaxation: 2^l * (N/l) * E[(S_l - d l/N)_+]  <= new bound at balanced
            relax = sum(comb(l, j)*max(Fr(0), Fr(N-d) - Fr(2*N*j, l)) for j in range(l+1))
            if relax > bound_new(l, N, d, bal): ok = False; print("FAIL relax", l, N, d)
print("balanced placement minimises both bounds for all l<=5, N<=12, d<N; real relaxation <= integer value:", ok)

# ---------- C. matrix check
def rank_fr(M):
    M = [r[:] for r in M]; rk = 0; rows = len(M); cols = len(M[0]) if M else 0
    for c in range(cols):
        p = None
        for i in range(rk, rows):
            if M[i][c] != 0: p = i; break
        if p is None: continue
        M[rk], M[p] = M[p], M[rk]
        for i in range(rk+1, rows):
            if M[i][c] != 0:
                f = M[i][c]/M[rk][c]
                M[i] = [a - f*b for a, b in zip(M[i], M[rk])]
        rk += 1
    return rk

def lagrange_eval(xs, ys, x):
    tot = Fr(0)
    for i in range(len(xs)):
        term = Fr(ys[i])
        for j in range(len(xs)):
            if j != i: term *= Fr(x - xs[j], xs[i] - xs[j])
        tot += term
    return tot

def layered_J_at(l, D, blocks, pts, t):
    """g_z(t) = chi(z) + p_B(sum_j s_j pi_j(t)) - p_B(sum_j s_j), |z| in B; returns matrix J(t) = H diag(g) H / 2^l (exact)."""
    n = 1 << l
    pi = []
    for j in range(l):
        v = Fr(1)
        for (tau, c) in pts:
            if c == j: v *= (1 - Fr(t)/tau)
        pi.append(v)
    g = []
    for z in range(n):
        w = bin(z).count('1'); B = next(b for b in blocks if w in b)
        xs = [l - 2*ww for ww in B]; ys = [(-1)**ww for ww in B]
        s = [1 - 2*((z >> j) & 1) for j in range(l)]
        y_t = sum(s[j]*pi[j] for j in range(l)); y_0 = sum(s)
        g.append((-1)**w + lagrange_eval(xs, ys, y_t) - lagrange_eval(xs, ys, y_0))
    J = [[sum(g[z]*(-1)**(bin((x ^ y) & z).count('1')) for z in range(n))/n for y in range(n)] for x in range(n)]
    return J

print("== C. general-class check on genuine matrices: multi-point layered J (deg_t = D*max m_c), r_k = rank J(tau_k)[x_c=0, y_c=1] (exact over Q)")
okC = True
tests = [(2, 1, (2, 1)), (2, 1, (2, 2)), (3, 1, (2, 1, 1)), (3, 1, (2, 2, 1)), (3, 1, (2, 2, 2)), (3, 2, (1, 1, 1)), (3, 2, (2, 1, 0)), (3, 1, (3, 1, 1)),
         (4, 1, (2, 1, 1, 1)), (4, 2, (1, 1, 1, 1)), (4, 1, (1, 1, 1, 1)), (4, 3, (1, 1, 1, 1)), (4, 2, (2, 2, 1, 1)), (4, 1, (2, 2, 2, 2)), (4, 1, (3, 0, 2, 1))]
for l, D, mc in tests:
    pts = []; t = 1
    for c in range(l):
        for _ in range(mc[c]): pts.append((t, c)); t += 1
    N = len(pts); d = D*max(mc)
    best = None
    # all partitions of weights 0..l into consecutive blocks of size <= D+1 (try a few: shifts)
    for shift in range(D+1):
        cuts = [0] + list(range(shift if shift > 0 else D+1, l+1, D+1)) + [l+1]
        cuts = sorted(set(cuts))
        blocks = [list(range(cuts[i], cuts[i+1])) for i in range(len(cuts)-1)]
        J0 = layered_J_at(l, D, blocks, pts, 0)
        n = 1 << l
        assert all(J0[x][y] == (1 if x ^ y == n-1 else 0) for x in range(n) for y in range(n))
        tot01 = 0; totmax = 0
        for (tau, c) in pts:
            Jt = layered_J_at(l, D, blocks, pts, tau)
            r0 = [x for x in range(n) if not x >> c & 1]; r1 = [x for x in range(n) if x >> c & 1]
            b01 = rank_fr([[Jt[x][y] for y in r1] for x in r0]); b10 = rank_fr([[Jt[x][y] for y in r0] for x in r1])
            tot01 += b01; totmax += max(b01, b10)
        if best is None or tot01 < best[0]: best = (tot01, totmax, blocks)
    bn = bound_new(l, N, d, mc); bo = bound_old(l, N, d, mc)
    good = best[0] >= bn; okC &= good
    print(f"l={l} D={D} m_c={mc} N={N} d={d}: sum rank B01 = {best[0]} (sum max(B01,B10) = {best[1]}), new bound {bn}, old {bo}  {'OK' if good else 'VIOLATION'}")
    sys.stdout.flush()
print("C all OK:", okC)

# ---------- D. asymptotics
def Epos(l, delta):      # E[(S_l - delta)_+], S_l = l - 2 Bin(l,1/2), exact float via logs
    tot = 0.0
    for j in range(l+1):
        k = l - 2*j
        if k > delta:
            tot += (k - delta)*exp(math.lgamma(l+1) - math.lgamma(j+1) - math.lgamma(l-j+1) - l*log(2))
    return tot
def star(l, delta, w):   # explicit lower bound (*) for E[(S-delta)_+]
    if w < 2 or delta + 2*w > l/sqrt(2): return None
    k = delta + 2*w
    return (w*w/4)/sqrt(2*l)*exp(-k*k/(2*l) - k**4/(6*l**3))
print("== D1. check explicit inequality (*): E[(S-delta)_+] >= (w^2/4)(2l)^(-1/2) exp(-(delta+2w)^2/2l - (delta+2w)^4/6l^3), w>=2, delta+2w<=l/sqrt2")
okD = True; worst = 1e9
for l in list(range(1, 60)) + [100, 200, 500, 1000, 3000]:
    for i in range(0, 41):
        delta = l*i/40*0.7
        for w in [2, 2.5, 3, 4, 6, 10, 20, 50, l/max(delta, 1e-9), sqrt(l)]:
            s = star(l, delta, w)
            if s is None: continue
            e = Epos(l, delta)
            if e < s*(1-1e-12): okD = False; print("FAIL", l, delta, w, e, s)
            if s > 1e-280: worst = min(worst, e/s)
print("(*) holds on grid:", okD, " min ratio true/(*) =", round(worst, 3))

def lay_rho(l, d):       # G7a one central block of d+1 weights: rho = (1/2) Pr[Bin(l-1,1/2) in boundary levels]
    a = (l - d)//2       # block = weights a..a+d ; boundaries: edges between a-1,a and a+d,a+d+1
    tot = 0
    if a - 1 >= 0: tot += comb(l-1, a-1)
    if a + d <= l - 1: tot += comb(l-1, a+d)
    return tot/2**(l-1)/2*1.0 if False else Fr(tot, 2**l)
def lb_new(l, d): return Epos(l, d)/l
def lb_old(l, d):
    m = (l-d-1)//2
    return sum(comb(l, j)*(m+1-j) for j in range(m+1))/(l*2**l) if m >= 0 else 0
print("== D2. threshold degree d*(l, l^-C): smallest d with bound <= l^-C  (old G11b | new G14d | layered UB | formula sqrt((2C-1) l ln l - 2 l lnln l))")
for C in (1, 2, 3):
    for l in (100, 400, 1000, 4000):
        tgt = l**(-C)
        do = next(d for d in range(l) if lb_old(l, d) <= tgt)
        dn = next(d for d in range(l) if lb_new(l, d) <= tgt)
        dl = next(d for d in range(l) if float(lay_rho(l, d)) <= tgt)
        form = sqrt(max(0, (2*C-1)*l*log(l) - 2*l*log(log(l))))
        print(f"C={C} l={l}: old {do} | new {dn} | layered {dl} | main-term formula {form:.1f} | sqrt((2C-1) l ln l) = {sqrt((2*C-1)*l*log(l)):.1f}")
print("== D3. ratio layered/new bound vs 2 d^2/l")
for l, d in ((1000, 83), (1000, 118), (1000, 144), (4000, 252), (4000, 354)):
    print(l, d, round(float(lay_rho(l, d))/lb_new(l, d), 2), round(2*d*d/l, 2), " old ratio", round(float(lay_rho(l, d))/lb_old(l, d), 2))
print("== D4. explicit corollary: rho_avg <= l^-C and sqrt(l) <= delta <= sqrt(l)*sqrt(l)/5  =>  delta^2 >= (2C-1) l ln l - 2 l ln(2C ln l) - 11.5 l - 108 C^2 ln^2 l")
okE = True
for C in (0.75, 1, 2, 3, 5):
    for l in (30, 100, 400, 1000, 4000, 20000):
        tgt = l**(-C)
        dn = next((d for d in range(l) if lb_new(l, d) <= tgt), None)
        rhs = (2*C-1)*l*log(l) - 2*l*log(2*C*log(l)) - 11.5*l - 108*C*C*log(l)**2
        if dn is not None and dn*dn < rhs: okE = False; print("FAIL", C, l, dn, rhs)
print("corollary consistent with exact thresholds:", okE)
