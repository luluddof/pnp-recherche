#!/usr/bin/env python3
"""G9d: escapes from the two-partition / product framework.  Cost model + exact small checks.
Units: instance = K nodes, S edge bits, time 2^{e(K) S};  u(K) := 1/e(K);  brute force e=1."""
import math, itertools
import numpy as np
LOG = math.log2

# ------------------------------------------------------------------ A. exact small-tensor checks
def eight_square():
    def phi(p, q):
        s = sum(p[k] * q[k2] for k in range(3) for k2 in range(k, 3))
        return (s + q[0]*p[1]*p[2] + p[0]*q[1]*p[2] + p[0]*p[1]*q[2]) % 2
    B = list(itertools.product([0, 1], repeat=3)); n = 8
    sg = np.array([[(-1) ** phi(p, q) for q in B] for p in B])
    xor = np.array([[int(''.join(map(str, [a ^ b for a, b in zip(p, q)])), 2) for q in B] for p in B])
    R = np.zeros((n,)*4, dtype=np.int64); C = np.zeros_like(R); T = np.zeros_like(R)
    for x, y, u, v in itertools.product(range(n), repeat=4):
        par = int(xor[xor[x, y], xor[u, v]] == 0)
        R[x, y, u, v] = par * sg[x, y] * sg[v, u]
        C[x, y, u, v] = par * sg[x, u] * sg[v, y]
        T[x, y, u, v] = 2 * int(x == v and y == u)
    assert (R + C == T).all()
    row = lambda A: np.linalg.matrix_rank(A.reshape(n*n, n*n).astype(float))
    col = lambda A: np.linalg.matrix_rank(A.transpose(0, 2, 1, 3).reshape(n*n, n*n).astype(float))
    print("A1 eight-square: R+C=2*crossing OK; (row,col) ranks  R:", (row(R), col(R)), " C:", (row(C), col(C)),
          " crossing:", (row(T), col(T)))
    # restriction x=v=0 -> identity in (y|u): each piece has rank <= min(row,col) there
    print("   restriction x=v=0: rank crossing|=", np.linalg.matrix_rank(T[0, :, :, 0].astype(float)),
          " R|=", np.linalg.matrix_rank(R[0, :, :, 0].astype(float)),
          " C|=", np.linalg.matrix_rank(C[0, :, :, 0].astype(float)))

def triangle(n=2):
    # MM tensor on parts A=(i,j) B=(j',k) C=(k',i'):  [j=j'][k=k'][i=i'];  3 edges of L=log2 n bits
    T = np.zeros((n*n,)*3)
    for i, j, k in itertools.product(range(n), repeat=3):
        T[i*n+j, j*n+k, k*n+i] = 1
    fl = [np.linalg.matrix_rank(np.moveaxis(T, a, 0).reshape(n*n, -1)) for a in range(3)]
    print(f"A2 triangle <{n},{n},{n}>: flattening ranks (A|BC),(B|AC),(C|AB) = {fl} = n^2 each; "
          f"plain bits/edge L={LOG(n):.2f}; any cut crossed by 2L={2*LOG(n):.2f} bits;"
          f" rank-one hyperedge needs >= log2 R(MM)  (n=2: log2 7={LOG(7):.3f} > 2L)")

# ------------------------------------------------------------------ B1. product-piece model, any g x h grid
def pp_level(eg, eh, g, h, rho, add=4.0):
    """One level K=g*h -> (rows: g nodes) or (cols: h nodes).  Product pieces on the line
    R+C = 2 rho f_x (plain: R=C=f_x), R,C in [(2rho-1) f_x, f_x]; adversarial global term."""
    K = g * h
    fr, fc = (h - 1) / (K - 1), (g - 1) / (K - 1)
    fx = 1 - fr - fc
    lo, hi = (2 * rho - 1) * fx, fx
    best = 0.0
    # cost_row(R) = eg*(fc+R)+add/g increasing in R ; cost_col = eh*(fr+2rho fx-R)+add/h decreasing
    cr = lambda R: eg * (fc + R) + add / g
    cc = lambda R: eh * (fr + 2 * rho * fx - R) + add / h
    Rstar = (eh * (fr + 2 * rho * fx) + add / h - eg * fc - add / g) / (eg + eh)
    Rstar = min(max(Rstar, lo), hi)
    return min(1.0, min(cr(Rstar), cc(Rstar)))   # max over [lo,hi] of min(increasing, decreasing)

def nonproduct_floor_level(eg, eh, g, h, add=4.0):
    """General two-direction model (arbitrary, non-product terms; balanced instances).
    Only proved constraint (simultaneous-diagonal restriction, quadrant X1Y1-X2Y2):
    some term has R >= kr*S_x (if row term) or C >= kc*S_x (if col term)."""
    K = g * h
    fr, fc = (h - 1) / (K - 1), (g - 1) / (K - 1); fx = 1 - fr - fc
    kr, kc = h / (4 * (h - 1)), g / (4 * (g - 1))
    return min(1.0, min(eg * (fc + kr * fx) + add / g, eh * (fr + kc * fx) + add / h))

def dp(X, level, x0=6):
    """e(2^x), x<=X integer; all splits x=a+b (a<=b), or no-op monotone."""
    e = {x: 1.0 for x in range(1, x0 + 1)}
    arg = {}
    for x in range(x0 + 1, X + 1):
        b, ba = e[x - 1], None
        for a in range(1, x // 2 + 1):
            v = level(e[a], e[x - a], 2.0 ** min(a, 500), 2.0 ** min(x - a, 500))
            if v < b: b, ba = v, a
        e[x], arg[x] = b, ba
    return e, arg

def report_pp():
    print("\nB1 product-piece model (any unbalanced grid allowed): u(K)=1/e(K), DP over all splits")
    print("   harmonic law at the floor: u(gh) <= u(g)+u(h)  =>  u = O(log K)")
    for rho, name in [(0.5, "rho=1/2 (floor)"), (2 ** -0.5, "rho=1/sqrt2 (App.B)"), (11 / 12, "rho=11/12")]:
        e, arg = dp(1024, lambda a, b, g, h: pp_level(a, b, g, h, rho))
        row = "   %-20s" % name
        for X in [24, 64, 256, 1024]:
            row += f" X={X}: u={1/e[X]:8.2f} u/X={1/e[X]/X:.3f} split={arg[X]}/{X} |"
        print(row)
    print("   (X = log2 K.  u/X bounded => e >= Omega(1/log K) => time 2^{Omega(m/log m)}; best split is balanced)")

def report_np():
    print("\nB3 general two-direction model, NON-product terms: what the provable floor still allows")
    e, arg = dp(512, nonproduct_floor_level)
    for X in [32, 64, 128, 256, 512]:
        print(f"   X={X}: floor e={e[X]:.3e}  -log2(e)/X={-LOG(e[X])/X:.3f}  best split a={arg[X]}")
    print("   variant with brute-force part tables (both directions forced): log2(1/e) vs sqrt(X)")
    e2, arg2 = dp(1024, nonproduct_floor_bruteW)
    for X in [64, 256, 1024]:
        print(f"   X={X}: floor e={e2[X]:.3e}  log2(1/e)={-LOG(e2[X]):.2f}  sqrt(X)={X**.5:.1f}  2*sqrt(X)={2*X**.5:.1f}  split a={arg2[X]}")
    print("   variant with recursive part evaluation (row terms of negligible cost):")
    for g in [2, 3, 4, 5, 6, 8, 16]:
        f = 1 / g + 0.25
        print(f"   constant g={g}: per-level factor 1/g+1/4={f:.3f}, e ~ K^-{LOG(1/f)/LOG(g):.3f}"
              f" -> time 2^(m^{1-0.5*LOG(1/f)/LOG(g):.3f})")

def nonproduct_floor_bruteW(eg, eh, g, h, add=2.0):
    """Same, but BOTH directions must occur (rank subadditivity across row cuts / column cuts:
    column terms alone cannot sum to T) and a merged part's W is tabulated by brute force
    (>= 2^{add*S/g} for a row part): e >= max(add/g, add/h, gain*e_h)."""
    K = g * h
    fr, fc = (h - 1) / (K - 1), (g - 1) / (K - 1); fx = 1 - fr - fc
    kr, kc = h / (4 * (h - 1)), g / (4 * (g - 1))
    return min(1.0, max(add / g, add / h, min(eg * (fc + kr * fx), eh * (fr + kc * fx))))

# ------------------------------------------------------------------ B2. fast matrix multiplication
def report_mm():
    print("\nB2 fast MM.  Bottom level on K0 nodes: 3 groups, trace(W_A W_B W_C), dimension 2^{2S/9}:")
    print("   e0 = 2*omega/9  (AKK single-network lower bound 4/9; rank-one floor 1/2; brute force 1)")
    for om, nm in [(2.0, "omega=2 (=ARC)"), (2.371339, "omega=2.371339"), (LOG(7), "Strassen"), (3.0, "naive")]:
        e0 = 2 * om / 9
        # full recursion, App.B pieces, squaring: e_l = (rho+2/K) e_{l-1} + 5/K
        out = []
        for m_log in [40, 100, 1000]:
            m = 2.0 ** m_log; xmax = 0.5 * LOG(6 * m)
            best = None
            for x0 in [6, 8, 10, 12]:
                x, e = x0, e0
                bst = 3 * m * e
                while 2 * x <= xmax:
                    e = (2 ** -0.5 + 2 / 2 ** x) * e + 5 / 2 ** x; x *= 2
                    bst = min(bst, 3 * m * e)
                best = bst if best is None else min(best, bst)
            out.append(f"log2m={m_log}: log2T/m={best/m:.4f}")
        print(f"   {nm:16s} e0={e0:.4f} (x{e0:.3f} vs brute)   " + "  ".join(out))
    print("   triangle gadget INSIDE a level (3 edges of L bits -> one rank-one hyperedge of omega*L bits):")
    for om in [2.0, 2.371339]:
        print(f"     omega={om}: 3 parts kept separate: {om/3:.3f} per bit but instance is a 3-uniform hypergraph (leaves class);"
              f" after any merge of 2 parts: omega*L/(2L) = {om/2:.3f} >= 1  (no gain)")
    print("   slice-type triangle pieces f(a)g(b,c): slice rank(MM) >= subrank = n^{2-o(1)} => rho_triangle = 1-o(1)")

if __name__ == "__main__":
    eight_square(); triangle(2); triangle(3)
    report_pp(); report_mm(); report_np()

# ============================================================================ NOTES / PROOFS (G9d)
# Status tags: [P] proof written here, unreviewed; [E] sketch; [M] holds only in the stated cost model;
# [K] known result quoted (not re-derived); [O] open.
#
# 0. GENERAL CLASS  G ("signed sums of networks, each with its own tree").
#    T = T_{K,lam}: the K-linear form (W_1..W_K) -> Z of eq.(2)/(5) (complete-graph tensor, N=2^lam,
#    S = C(K,2) lam).  An algorithm in G:  c*T = sum_i D_i, D_i a tensor network on the W_a plus fixed
#    tensors, contracted along a tree tau_i; cost = sum_i cost(D_i,tau_i).  G is contained in
#    set-multilinear arithmetic circuits (one bilinear gate per contraction); input size M = K N^{K-1}
#    = K 2^{2S/K}, brute force 2^S = M^{K/2}; target 2^{eS} = M^{eK/2}.
#    Floors:  (F1)[K] rank-one terms (product depth 1): >= N^{floor(K^2/4)} ~ 2^{S/2} (flattening);
#             asymptotic-rank conjecture (K-party form: asymptotic rank = max flattening) gives exactly
#             this, e -> 1/2; grouped in 3 blocks: MM tensor <n,n,n>, n = 2^{2S/9}: e0 = 2 omega/9.
#             (F2)[K, via the paper's App. G] one network / one tree (AKK22): >= 2^{4S/9}; = F1(3 blocks)
#             at omega=2.  So ARC / omega=2 buy a CONSTANT factor in the exponent, never 2^{o(S)}.
#             (F3)[M] product-piece schemes: 2^{Omega(S/log K)} (section 2).
#             (F4)[P] whole class: only K 2^{2S/K} (every table entry must be read; also = subrank
#             bound N^{K-1}); with lam>=1, K<=sqrt(6m): 2^{Theta(sqrt m)}.  [O] nothing better: a
#             flattening rank across a node cut is NOT a progress measure in G, because a term
#             contracted along another partition has unconstrained rank there (the crossing identity
#             itself is the witness).  2^{S^{1-eps}} is consistent with every bound I know.
#
# 1. (a) TRIANGLES / FAST MM.
#    Lemma T [P].  Parts A,B,C pairwise joined by L bits (tensor <n,n,n>, n=2^L).
#      (i) a term f(z_A,th) g(z_B,th) h(z_C,th) summed over th: sum needs >= R(<n,n,n>) >= n^2 values,
#          i.e. a 3-hyperedge of >= omega L >= 2L bits.  After merging any two of the three parts the
#          plain edges cost 2L.  So inside a recursion on graph instances: factor omega/2 >= 1, no gain.
#          Kept unmerged: one-time factor omega/3 (0.667 / 0.790), instance becomes a GHZ hypergraph
#          (rank = flattening rank: no second application).  This IS Williams/Nesetril-Poljak.
#      (ii) slice terms f(z_A,th) g(z_B,z_C,th) (low rank across one of the 3 cuts): slice rank of MM
#          >= subrank Q(<n,n,n>) = n^{2-o(1)} (Strassen; slice rank(<r>) = r, Tao) [K]; exact value
#          ceil(3n^2/4) is quoted from memory (BCCGNSU17), not re-checked.  So with poly(L) terms some
#          term costs 2L - O(log L): rho_triangle = 1 - o(1).
#      Reason: the crossing is compressible because it has a third cut of rank 1; the triangle is a
#      tight tensor with full subrank on every cut.
#    Combination that does work: fast MM at the BOTTOM level only: e_0 = 2 omega/9 instead of 1
#    (0.444 for omega=2; 0.527 for 2.371339), e_l = e_0 * (paper's e_l): time
#    2^{(2 omega/9) c m / sqrt(log m)}; the shape and the m/log m cap are unchanged (table B2).
#
# 2. HARMONIC LAW  [M, E] (generalises G4b Thm 1 and the G6b conclusion; covers unbalanced grids,
#    d directions, gadgets of any arity, any instance graph).
#    Model PP: gadgets = disjoint sets of equality bits; c_gam T_gam = sum_{i<=poly} A_{gam,i};
#    a global term = tensor product of one piece per gadget, evaluated after merging along ONE of
#    the partitions pi_1..pi_d (g_j parts); an instance with g nodes and R bits costs 2^{R/u(g)}.
#    Claim: u(K) <= sum_j u(g_j) on the rewritten bits.
#    Proof: rank subadditivity across pi_j-aligned cuts: for each j some piece A_{gam,i_j} has full
#    pi_j-cost F - O(log #pieces).  The global term using i_j on a fraction th_j of the gadgets has
#    pi_j-cost >= th_j F (ranks multiply under tensor product), so its best direction costs
#    >= F max_th min_j th_j/u(g_j) = F / sum_j u(g_j).   QED
#    Geometry (G4b Lemma 1): rewriting a constant fraction needs prod_j g_j <= K (cells non-empty).
#    Hence u(K) <= c log K, e >= Omega(1/log K), time 2^{Omega(m/log m)}; at rho = 1/2 every split is
#    equally good, for rho > 1/2 balanced squaring is optimal (DP B1: u/log2 K -> 0.16, 0, 0).
#    (b) sparse / non-clique instances: the law does not use the instance graph; K may grow to
#        Theta(m) instead of sqrt(6m), which only removes floor F4, not the 1/log K cap.
#    (c) reused equalities ([x=y] = [x=y]^k): pieces are then multiplied entrywise, bits add
#        (rank is submultiplicative for Hadamard products); k uses cost k*rho >= 1 for k >= 2:
#        break-even exactly at the floor, a loss otherwise [M].  The PROOF of the law breaks
#        (Hadamard products can lose rank), but I found no mechanism that gains.  [O]
#
# 3. WHAT IS LEFT: NON-PRODUCT ("entangled") TERMS  [E].
#    Drop "one piece per gadget": c*T = sum_j B_j, each B_j an arbitrary tensor evaluated along one
#    of d partitions.  Only bound I can prove: choose cuts X_1..X_d (one per partition), keep the
#    EPR pairs joining one orthant class to its opposite, fix all other strings: the restriction is
#    the same identity matrix in all d flattenings, so some term has own-direction cut bits >= that:
#    per-level gain >= 4^{1-d} (d=2, g=h=2: 1/2 = the P(N) >= N bound; large g,h: 1/4), against
#    1/d for product pieces (G6b).  Also both directions must occur (column terms alone violate rank
#    subadditivity across column cuts).  Consequences (table B3):
#      - merged parts tabulated by brute force (2^{2S/g}): log(1/e) <= ~2 sqrt(log2 K), i.e. the
#        floor is 2^{m / 2^{O(sqrt(log m))}} (same shape as G4b Thm 2, now for d=2 non-product).
#      - merged parts evaluated recursively: floor e >= K^{-1/2} (g=4), i.e. 2^{Omega(m^{3/4})}.
#    So 2^{m^{1-eps}} is NOT excluded here, but needs BOTH
#      (E1) a joint decomposition of the 4 x h grid tensor into poly(lam)^{poly(h)} terms: column terms
#           with <= (1/2+o(1)) S balanced bits, row terms of row rank 2^{o(S)};  and
#      (E2) factor functions structured enough that a merged row is evaluated recursively.
#    For a single crossing, E1-type statements live between P(N) >= N and P(N) <= N^{sqrt 2}
#    (Naslund Problem 12 / bilinear Hurwitz, HWY11); P(N) >= N^{1+eps} would raise this floor.  [O]
#
# ADVERSARIAL SELF-REVIEW.
#  - Section 2 is a cost-model statement: "cost = 2^{bits/u(g)}" is postulated, the instances reached
#    may be easier than worst case; additive terms and balancing are ignored in the proof (kept in DP).
#  - The floors of section 3 are only what one restriction argument gives; they say what is not yet
#    excluded, not what is achievable.  No scheme was constructed.  Balanced terms are assumed.
#  - "K-party ARC" is the natural extension of Strassen's conjecture, not its standard statement.
#  - AKK22 and App. G numbers are taken from the (unrefereed) paper text.
