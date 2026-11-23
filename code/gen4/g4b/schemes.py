#!/usr/bin/env python3
"""G4b: compare recursion geometries for Claude's 2^{o(n)} 3-SAT algorithm.

Cost model (units: log2).  An instance is a (K, lam)-clique with S = K^2 lam / 2 bits
on its edges; level l has exponent e_l, i.e. time 2^{e_l S} (times branching).
One level takes K_big nodes to K_small nodes:
    e_new = (1 - (1-rho) * phi) * e_old + add / K_small
  phi  = fraction of edge bits that are rewritten by gadgets (crossings, boxes),
  rho  = factor by which a rewritten bit's cost is multiplied (crossings: rho >= 1/2),
  add/K_small = cost of evaluating the W of one merged group (2^{~2S/K_small}).
Branching: (#pieces)^{#gadgets}; for the clique, #gadgets ~ K_big^2/4, #pieces = poly(lam).
On a formula with m clauses: S ~ 3m and lam >= 1, so K_top <= sqrt(6m) ALWAYS.
log2(time) ~ e_top * 3m + branching_log2(top level).
"""
import math

LOG = math.log2


def run(seq_x, gain, x0, add=4.0):
    """seq_x: list of log2 K_l (l=0..L).  gain(xp, xn) -> (rho_eff multiplier).
    Returns e_L."""
    e = 1.0
    for xp, xn in zip(seq_x, seq_x[1:]):
        e = gain(xp, xn) * e + add / 2 ** xp
    return min(e, 1.0)


def time_log2(m, seq_x, gain, branch_bits_per_gadget):
    e = run(seq_x, gain, seq_x[0])
    x = seq_x[-1]
    lam = max(1.0, 6 * m / 2 ** (2 * x))
    branching = 2 ** (2 * x - 2) * branch_bits_per_gadget(lam)
    return e * 3 * m + branching, e


def best_over_levels(m, step, gain, bb, x0=6.0, xmax=None):
    """step(x) -> next x; choose number of levels L minimising time."""
    xmax = xmax or 0.5 * LOG(6 * m)
    best = (3 * m, 1.0, 0)
    seq = [x0]
    while True:
        nx = step(seq[-1])
        if nx > xmax:
            break
        seq.append(nx)
        t, e = time_log2(m, seq, gain, bb)
        if t < best[0]:
            best = (t, e, len(seq) - 1)
    return best


def poly_branch(lam):            # Lemma 2 / App. B: log2(#pieces) ~ log2(lam)+O(1)
    return LOG(lam) + 2


def phi_two_partition(xp, xn, C=2.0):
    """Proved upper bound (two-partition model, see PROOFS at end of file): phi <= min(1, C^2 K_big/K_small^2,
    (D+1) D^2 / K_big^2) with D = K_big - K_small."""
    if xn <= xp:
        return 0.0
    a = 2 * LOG(C) + xn - 2 * xp                      # log2 of C^2 K_big/K_small^2
    lD = xn + LOG(1 - 2 ** (xp - xn))                 # log2 D
    b = lD + LOG(2 ** lD + 1) + lD - 2 * xn if lD < 60 else 3 * lD - 2 * xn
    return min(1.0, 2 ** min(a, 0.0), 2 ** min(b, 0.0))


def scheme_table(ms):
    s2 = 2 ** 0.5
    schemes = {
        # name: (step, gain)
        "grid, rho=11/12, squaring (paper main)":
            (lambda x: 2 * x, lambda a, b: 11 / 12),
        "grid, rho=1/sqrt2, squaring (App. B)":
            (lambda x: 2 * x, lambda a, b: 1 / s2),
        "grid, rho=1/2 floor, squaring":
            (lambda x: 2 * x, lambda a, b: 0.5),
        "3-box, product pieces rho_3=1/sqrt2, cube":
            (lambda x: 3 * x, lambda a, b: 1 / s2),
        "3-box, FLOOR rho_3=1/4 (spec.), cube":
            (lambda x: 3 * x, lambda a, b: 0.25),
        "6-box, FLOOR rho_6=1/32 (spec.), ^6":
            (lambda x: 6 * x, lambda a, b: 1 / 32),
        "const growth K*2, rho=1/2, phi<=bound (proved)":
            (lambda x: x + 1, lambda a, b: 1 - 0.5 * phi_two_partition(a, b)),
        "growth K^1.5, rho=1/2, phi<=bound (proved)":
            (lambda x: 1.5 * x, lambda a, b: 1 - 0.5 * phi_two_partition(a, b)),
        "HYPOTHETICAL const growth K*2, rho=1/2, phi=1":
            (lambda x: x + 1, lambda a, b: 0.5),
    }
    print(f"{'scheme':50s}" + "".join(f"  log2m={LOG(m):5.0f}" for m in ms))
    for name, (step, gain) in schemes.items():
        row = f"{name:50s}"
        for m in ms:
            t, e, L = best_over_levels(m, step, gain, poly_branch)
            # report  t/m  (fraction of brute-force-ish exponent 3m -> 3) and L
            row += f"  {t / m:7.2e}/L{L:<3d}"
        print(row)
    print("\nreference curves, value = log2(time)/m:")
    for name, f in [("m/log m", lambda m: 1 / LOG(m)),
                    ("m/(log m)^2", lambda m: 1 / LOG(m) ** 2),
                    ("m/2^sqrt(log m)", lambda m: 2 ** -math.sqrt(LOG(m))),
                    ("m^{0.9}  (target 2^{n^{1-eps}})", lambda m: m ** -0.1)]:
        print(f"{name:50s}" + "".join(f"  {f(m):7.2e}     " for m in ms))


# ---------------------------------------------------------------- Theorem 1 check
def max_total_gain_two_partition(X, x0=6.0, rho=0.5, h=0.05, C=2.0):
    """DP over all sequences x0 < x1 < ... <= X (grid step h) maximising
    sum log2 1/(1-(1-rho) phi) with phi the proved bound.  Compare with
    log2(1/rho) * log2(X/x0)."""
    n = int(round((X - x0) / h)) + 1
    xs = [x0 + i * h for i in range(n)]
    best = [-1e9] * n
    best[0] = 0.0
    for j in range(1, n):
        bj = -1e9
        for i in range(j):
            phi = phi_two_partition(xs[i], xs[j], C)
            g = -LOG(1 - (1 - rho) * phi) if phi > 0 else 0.0
            bj = max(bj, best[i] + g)
        best[j] = bj
    return best[-1]


# ---------------------------------------------------------------- d-box at the floor
def dbox_floor_best(X, x0=6, add=4.0):
    """Speculative: d-box gadgets reaching the restriction floor rho_d = 2^{1-d},
    any d per level (K_big = K_small^d).  DP over integer x = log2 K.
    E(x) = min_{y*d <= x} 2^{1-d} E(y) + add/2^y ."""
    E = {}
    for x in range(1, X + 1):
        if x <= x0:
            E[x] = 1.0
            continue
        b = E[x - 1]
        for y in range(x0, x):
            d = x // y
            if d >= 2:
                b = min(b, 2.0 ** (1 - d) * E[y] + add * 2.0 ** (-y))
        E[x] = min(b, 1.0)
    return E[X]


if __name__ == "__main__":
    ms = [2.0 ** 40, 2.0 ** 100, 2.0 ** 300, 2.0 ** 1000]
    print("value printed = log2(time)/m (brute force ~ 3), /L = best number of levels\n")
    scheme_table(ms)

    print("\nTheorem 1 (two-partition model): max total log-gain vs log2(1/rho)*log2(X/x0)")
    for X in [12, 24, 48, 96]:
        g = max_total_gain_two_partition(X, h=0.25 if X > 48 else 0.1)
        print(f"  X=log2 K_top={X:4d}: DP max gain = {g:6.2f}   "
              f"doubling gives {LOG(X / 6):6.2f}")

    print("\nSpeculative d-box at floor rho_d=2^{1-d}: log2(1/e) vs sqrt(X)")
    for X in [64, 256, 1024, 4096]:
        e = dbox_floor_best(X)
        print(f"  X={X:5d}: log2(1/e)={-LOG(e):7.2f}   sqrt(X)={math.sqrt(X):6.2f}"
              f"   log2(X)={LOG(X):5.2f}")

# ============================================================================ PROOFS
# Model M2 ("two-partition"): one level takes K = K_big nodes, two partitions P1 (rows),
# P2 (cols), each with >= K_s parts of size <= s = C K/K_s (needed: W of a part costs
# 2^{~2 s S/K}, must be O(S/K_s)).  Gadgets = crossings [x=v][y=u] on edges {a,b'},{a',b}
# with a,a' same row, a,b same col (cells may hold several nodes / be empty; every equality
# bit used in one gadget).  A rewritten bit is multiplied by >= rho (rho >= 1/2 floor).
#
# Lemma 1 (geometry).  phi := rewritten fraction of edge bits satisfies
#   phi <= sum_a (|row a|-1)(|col a|-1) / K^2 <= min( C^2 K/K_s^2 , (D+1) D^2 / K^2 ),
#   D = K - K_s.
# Proof: a rewritten bit of an edge {a,b'} is paired with a bit of an edge {a',b},
# a' in row(a)\{a}, b in col(a)\{a}; pairing is injective, so the rewritten bits at a are
# <= (|row a|-1)(|col a|-1) lam.  Sum over a (each bit counted at 2 ends), divide by
# S = K^2 lam/2.  First bound: parts <= CK/K_s.  Second: AM-GM and sum_r(|r|-1) <= D.
#
# Lemma 2 (merging never helps).  Bits that become internal to a part are paid in the
# additive term h = bits held by the largest part >= 2 * (sum of internal bits removed)/...;
# for e <= 2 one gets e_new >= e_old * (1 - (1-rho) phi).  (Pure merges: e_new >= e_old.)
#
# Theorem 1.  In M2, with x_l = log2 K_l <= X = log2 K_top <= (1/2) log2(6m):
#   sum_l phi_l <= log2(X/x_0) + O(log C),  hence  e_top >= 2^{-O(1)} (X)^{-log2(1/rho)}
#   >= Omega(1/log m) for rho >= 1/2, i.e. time >= 2^{Omega(m/log m)}.
# Proof sketch: convexity: -log2(1-(1-rho)phi) <= phi log2(1/rho).  Per level with
# u = (x_l - x_{l-1})/x_{l-1}:  phi_l <= log2(1+u) + H/x_{l-1} for u > 1/2 (H = O(log C)),
# phi_l <= C^2 2^{-x_{l-1}/2} for 1/x <= u <= 1/2 (x increases by >= 1: geometric sum),
# and levels with K_l <= 2 K_{l-1}: by the second bound of Lemma 1 a dyadic block
# K in [2^j, 2^{j+1}] contributes <= O(C^2 2^{-j/3}) (summable).  Levels with u > 1/2 are
# <= log_{1.5} X many with x growing geometrically, so sum H/x = O(1).  Telescoping
# sum log2(1+u_l) = log2(X/x_0).  The DP max_total_gain_two_partition confirms:
# max gain - log2(X/x_0) stays ~0.25-0.4 for X = 12..96.
#
# Consequences: (i) constant-factor growth of K (idea 2) gives total gain O(1);
# (ii) rewriting a fraction / K^{1+a} growth is worse than squaring;
# (iii) branching (idea 3) is irrelevant for the exponent: lam >= 1 forces K_top <= sqrt(6m)
#      whatever the branching, so poly(K) branching changes e only by a constant factor.
#
# Model Md ("d-box", SPECULATIVE pieces): nodes on [g]^d, merge along one coordinate,
# K_big = K_s^d, gadget = the D = 2^{d-1} main diagonals of a box.  Restricting all but one
# diagonal gives sum_pieces 2^{min_i c_i} >= 2^lam, so with poly(lam) pieces and the
# Definition-7 style per-piece bound sum_i c_i <= d rho_d D lam:  rho_d >= 2^{1-d} (floor).
# Pairing diagonals into crossings only gives rho_d = rho_2 (worse than iterating d=2,
# which gives rho_2^{log2 d}); beating d=2 needs a genuinely d-way identity with
# rho_d < rho_2^{log2 d} (d=4: < 1/4 if rho_2 -> 1/2; floor 1/8).
# Theorem 2.  Even at the floor, with free choice of d per level, the recursion
#   log2(1/e_l) <= min(x_{l-1}, (d_l - 1) + log2(1/e_{l-1})),  x_l = d_l x_{l-1}
# (first term: additive cost 2^{-x_{l-1}}) gives log2(1/e_top) <= O(sqrt(X)), i.e.
# time >= 2^{m / 2^{O(sqrt(log m))}}; DP dbox_floor_best: log2(1/e) ~ sqrt(X).
