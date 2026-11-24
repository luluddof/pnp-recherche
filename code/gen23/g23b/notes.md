# G23b — adversarial review of gen22/g22b/notes.md (monomial / intercalate model)
Tags: [P] re-derived by hand here, [E] exact computation here (own code, independent of shape.py), [L] read in HWY ECCC TR10-021.

## Files (this folder)
core.py (tensor identity via matrix products, predicate, flint ranks, key-lemma checker), t_thm1.py (exhaustive Thm 1 test),
t_corrupt.py (corrupted arrays), enum_st.py (exhaustive ST enumeration N<=4), mysat.py (own SAT encoding, at-most colour counts),
front5.py, drv5.sh (+ drv5.txt, sat666.txt, enum4.txt), t_misc.py.

## Claim (1) correspondence — CONFIRMED (hypothesis redundant)
[P] Coefficient at x=v, y!=u: C-term K[x,u]K[x,y], R-term K[x,y]K[x,u] — the SAME pair of cells, total weight wC+wR = 1.
So the identity alone forces EVERY class (any palette) to be a partial permutation; the hypothesis "S,U classes are partial
permutations" is unnecessary. 2x2 case: E1 wC(D)sD + wR(A)sA = 0, E2 wC(A)sA + wR(D)sD = 0, solved case by case: OK.
Ranks: disjoint supports -> exact: OK.
[E] t_thm1.py, tensor built from the GENERAL definition C = sum G(x)G', R = sum H(x)H' (G = wK, G' = K^T):
 N=2: all 309 (partition, palette) x 16 signs (incl. non-partial-permutation classes): 0 mismatches.
 N=3: all 586,656 partial-permutation structures x STU palettes x all 512 signs: 0 mismatches (11.76M identities true).
 N=3: 2% sample of ALL partitions (incl. non-pp) x 512 signs: 0 mismatches.
 t_corrupt.py: 1500 corruptions (sign flip, palette change, merge, swap, recolour, all-U) of the 10 arrays on file: 0 mismatches.
 All arrays on file (shape_*, grp_8_Q8, transl_*) satisfy identity, R+C=X, ranks = colour counts (flint).
Genuine special case: G_l = w K_l, G'_l = K_l^T, R = X - C of rank exactly |T|+|U|. Integer identity for ST (any ring);
U needs 1/2 (char != 2). So model constructions are general constructions; model lower bounds are not.

## Claim (2) key lemma — CONFIRMED
Uses: partial-permutation property (cells of the class have distinct rows/cols; the two S-cells (a_i,b_k),(a_i',b_k') then have
i!=i', k!=k', so the 2x2 is genuine) and BICHROMATICITY in both directions (T-diag => S-antidiag; S-diag => T-antidiag).
The SIGN condition is NOT used. [E] enum_st.py: all sign-free-valid ST structures N=2,3,4 (18 / 1100 / 612,736; all happen to
be signable): 0 lemma failures, 0 counting failures. g*(t) := max{g <= N : g(g-1)/2 <= t} is exactly the max class size of a
T-colour when |S| = t (resp. S-colour when |T| = t) allowed by the lemma; g*(t) <= sqrt(2t) + 1, and <= 1/2 + sqrt(2t + 1/4).

## Claim (3) counting — CONFIRMED (ST) / CONFIRMED modulo HWY Prop 7.8 (STU; Lemma 7.7 + Prop 7.8 re-checked by hand here)
N^2 = sum |class| <= q g*(p) + p g*(q). Balanced table recomputed: N=8:8, 16:21, 32:52, 64:128, 128:325. q<=p trade-off
algebra OK (q(sqrt(2p)+1) <= p(sqrt(2q)+1)). q=1: N^2 <= g*(p) + 2p <= N + 2p, construction verified N=2..7 (ranks 1, N(N-1)/2).
STU: U-class submatrix is a full HWY matrix in U-colours; 2 mu^2 >= C(s,3) is HWY's Lemma 7.4 proof (each ordered colour pair at
<= 2 positions, Prop 7.8). Note HWY's "at least k^3/8 positions" is false for k <= 11 but the notes use C(s,3), fine. Prop 7.8
checked [E] on the U-classes of the N=9,10 arrays (max 2 positions per pair). Constant check: for q <= N^{6/5}/5,
mu s*(mu) <= 0.157 N^2 + 0.4 N^{6/5} <= N^2/2 for all N >= 2, so the STU trade-off holds for all N, not only N large.
Could not break (3): impossible by the proof; exhaustive N<=4 and SAT frontier N=5 respect it.

## Claim (4) data — CONFIRMED
[E] own SAT (mysat.py, at-most colour counts):
 N=5 ST frontier q<=1..10 -> min p 10,10,10,9,8,7,6,5,4,1; counting bound 10,9,6,5,4,3,3,3,2,1: never violated, tight only q=1.
 N=4 exhaustive frontier = SAT frontier; tight arrays (slack 0): (1,1)@2, q=1 families, (8,8)@8 file array slack 0.
 (6,6)@6 ST: UNSAT (no class-size cap, 65.6 s) — independent of the cap used in g22b.
 N=5, (p,q)=(7,4): (S,T,U)<=(4-m,7-m,m), m=0..4 all UNSAT; (6,6): (6-m,6-m,m), m=0..6 all UNSAT; (6,7) SAT.
 => both exact N=5 optima lie outside the STU model.

## Lean candidates (sign-free, finite combinatorics)
 col : Fin N -> Fin N -> κ, pal : κ -> Bool (true = S)
 hrow : col a b = col a b' -> b = b';  hcol : col a b = col a' b -> a = a'
 hint : x != v -> y != u -> col x y = col v u -> col x u = col v y /\ pal (col x u) != pal (col x y)
 KEY LEMMA: for any colour j, with Mj = {(a,b) | col a b = j}:  (Mj.card).choose 2 <=
   card {c | (exists a b, col a b = c) /\ pal c != pal j}.   (injection {(a,b),(a',b')} |-> col a b')
 COUNTING: N^2 <= q * gstar N p + p * gstar N q, p = #{used c | pal c = false}, q = #{used c | pal c = true},
   gstar N t = max {g <= N | g(g-1)/2 <= t}; corollary q <= p -> N^2 <= 2 p (sqrt(2q) + 1).
