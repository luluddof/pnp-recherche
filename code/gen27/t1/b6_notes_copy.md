# G26-B6 — local / perturbative improvement of Appendix B (working notes, saved as I go)
Tags: [P] short proof here (unreviewed), [E] exact computation (mod-p ranks unless stated), [H] heuristic.

## 0. Structure of App B used here [P, from (47)-(55) of the paper]
Sub-terms: beta_{s,i} = (-1)^i <b[h-s-i] a<i>[s-h+i]>  (row-cheap: factors through V_{s+i}),
           alpha_{s,i} = (-1)^i <a<i>[s-h+i] b[h-s-i]>  (col-cheap: factors through V_{2h-s-i}).
sum_i beta_{s,i} = sum_i alpha_{s,i} = A_s (Lemma 9a holds only for the i-sum, not termwise).
KEY: the (u,v)-factor of beta_{s,i}, g_sigma(u,v) = [t^mu(u)] b_v[h-|sigma|] t^sigma, depends on sigma only (not on s, i, x).
Likewise the (y,v)-factor of alpha, g'_sigma(y,v) = [t^sigma] b_v[|sigma|-h] t^mu(y), depends on sigma only.
So every R built from beta's of intermediate weight <= W lies in Q^{N^2}_(xy) (x) U_W, U_W := span{g_sigma : |sigma| <= W},
and every C built from alpha's of weight <= W' lies in Q^{N^2}_(xu) (x) U'_W'.  Threshold at s0: (W, W') = (2h+s0-1, 2h-s0),
W + W' = 4h - 1  <=>  theta^2 + T^2 = 4 (since log p(w)/log p(h) ~ sqrt(w/h)).
One-variable picture: F(x) = (1-x)^H C(x) = (1-x)^H D(x), C,D the col/row-cheap generating series, H = h + weight(a) <= 2h:
the extra 2h in "4h" is the pole order H of the generating function at x = 1.

## 1. Brainstorm (ideas, estimate at theta = sqrt2, barrier status)
LEMMA 1-D [P]: for ONE pair (a,b) with generic Laurent F = (1-x)^H C = (1-x)^H D, the row-cheap data d_r (r <= r_R, weight
h+r) see exactly F_k, k <= r_R, and the col-cheap data c_r (r >= r_C, weight h-r) see exactly F_k, k >= r_C + H (triangular).
F(1) = sum_k F_k needs r_R >= r_C + H - 1, i.e. W + W' >= 2h + H - 1 = 4h - 1: thresholding is OPTIMAL pair-by-pair.
Any gain must come from linear relations ACROSS the N^4 functions (x,y,u,v) -- i.e. from the envelope of Test A.
 1. Re-split a piece between R and C via its sub-terms (s,i) (partial (54)/(55) use) -- special case of Test A. Tested.
 2. Merge adjacent pieces / any linear recombination sharing App B factors -- special case of Test A. Tested.
 3. Windows in a second statistic (number of parts l(sigma)) -- Test B. Estimate: typical sigma of weight w has ~sqrt(w) log w
    parts and the pieces need all of them; cutting l only removes a sub-exponential fraction unless l << sqrt w. No gain. Tested.
 4. Lower-weight dual elements a_x (smaller pole order H = h + wt a): the functionals a[-h]|V_h from wt <= h* span
    <= P(h*) dims, need p(h): h* >= h - O(sqrt h log h). Gain o(h) in W+W' => o(1) in exponent (finite-N constant only).
 5. Other multiplier q(x)(1-x)^H, or truncated/reweighted binomials a<i>: Lemma 1-D, pole order fixed => none.
 6. Linear combination of App B at two h (restricted to the same N): ranks add, each side is a circle point; no better than
    the better one up to factor 2. Covered by G18g Thm 1(ii) / G20b sect. 2 (only uniform families matter).
 7. Hurwitz/8x8 blocks at the top weights only: a fixed-size block = constant factor (G19b sect. 6, G20b sect. 2).
 8. Change of index basis (Schur, monomial instead of power sums): GL_N x GL_N invariance (G25-P7 sect. 0) - ranks unchanged.
 9. Reduce App B mod a small prime p > h (h! invertible): ranks can only drop. Tested (Test C). Brief: char p no exponent gain.
10. Asymmetric weights (b in V_{h_b}, indices in V_E, a of weight h_a): N fixed by E, H = h_a + h_b; any asymmetry raises
    (W+W')/E above 4 [H]. No gain.

## 2. Exact tests
### 2a. Test A = envelope of ALL perturbations that keep App B's factor spaces (relax.py, verify.py)  [E]
Question: X = R + C with R in Q^{N^2}_(xy) (x) U_W, C in Q^{N^2}_(xu) (x) U'_W' (x-side free; this contains every threshold,
every re-splitting of a piece by sub-terms (s,i), every merge of pieces with shared factors, every partial move of A_{s0}).
Reduction: K = annihilator of U_W; decouples over x: one linear system M f = K_j-columns, M[(j,y),(m,u)] = (G'_m K_j^T)[y,u].
Ranks mod two primes (33554393, 16777213), agree everywhere.  Equivalent operator form [P]: for each x find
rho_x: V_h -> V_{<=W}, kappa_x: V_{<=W'} -> V_h with  B_v^- rho_x + kappa_x B_v^+ = [v=x] Id  for all v
(B_v^- sigma = b_v[h-|sigma|] sigma, B_v^+ y = sum_{w<=W'} b_v[w-h] y).
Results (min W' for each W; dims of U_W, U'_W'; App B threshold has W' = 4h-1-W):
 h=3 N=3 : (W,W') = (2,3) dims (4,6); (3,2) dims (6,4).          App B threshold balanced 8.
 h=4 N=5 : (4,5) dims (11,17); (5,5) (16,17); (6,0) (24,1).       App B 22.
 h=5 N=7 : (4,7) (12,39); (6,6) (27,30); (7,4) (40,12).           App B 46.
           EXACT over Q (verify.py, explicit R, C): (W,W')=(6,6): rr(R)=27, rc(C)=22;  (7,4): rr=40, rc=10.
           h=6 (8,8): exact rr(R)=62, rc(C)=65 (particular solution, not rank-minimised).
 h=6 N=11: (7,9) (42,93); (8,8) (62,67); (9,7) (89,45); (10,4) (116,12).   App B 111.
 h=7 N=15: (8,11) (64,188); (9,10) (92,139); (10,9) (131,97); (11,7) (182,45).  App B 215.  (all sums = 225 + 2..6)
READING: much better than App B's threshold at the same N (balanced 30 vs 46 at N=7, 67 vs 111 at N=11), BUT every point
satisfies dimU + dimU' = N^2 + O(10) (exact h=5: 27 + 22 = 49 = N^2): it is the TRIVIAL ANTIDIAGONAL p + q = N^2
(achieved directly by splitting the pairs (x,y) of X into two sets). App B's thresholds are WORSE than trivial at these N
(G20b), so "beating App B like-for-like" here only means reaching the trivial line. No point below p + q = N^2 found.
The balanced weights W = W' ~ h + 2 are a saturation artefact (dims reach N^2/2 there), not a sign of W + W' < 4h.

### 2b. Test B: number-of-parts windows inside the envelope (lwin.py; out_lwin_h*.txt)  [E mod p]
S = {|sigma| <= W, l(sigma) <= L}, S' likewise, full grid, minimal W' per (W,L,L').  Minimum of dimU + dimU' - N^2:
 h=4 (N=5): 0 (dims (24,1));  h=5 (N=7): 0 (dims (40,9)), then 1 at (2,48), (12,38), (27,23);  h=6: see out_lwin_h6.txt.
h=4,5: never below the antidiagonal. h=6 (partial, timeout): excess -3 at dims (2,116) -- but App B's own endpoint is
(1, 111) = (1, N^2 - N) (excess -9): near the ENDS the true baseline is slightly below p+q = N^2; nothing at the balanced region.
(Dims are upper bounds on ranks; actual ranks can be lower, e.g. h=5 dim 30 -> rank 22.)
### 2c. Test C: App B thresholds mod p (modp.py; out_modp_h5/h6.txt)  [E]
 h=5 N=7: balanced Q 46 | p=7: 43 | 11: 43 | 13: 43.   h=6 N=11: Q 111 | p=7: 100 | 11: 111 | 13: 107.
 Drops only in the capped (near-N^2) entries, a few % ; the uncapped ranks (e.g. rc = 45, 30, 19, 12 at h=6) are identical
 mod 11, 13 (p=7: 40, 29, 19, 12). No structural effect.

## 3. Verdict
 * Within App B's own factor spaces (all re-splittings, merges, partial moves, second-statistic windows), the exact optimum
   at h <= 7 (N <= 15) lies ON the trivial antidiagonal p + q = N^2 (excess 0..6), never below it. This beats App B's
   thresholds like-for-like (balanced 30/46, 67/111, 131/215 at N = 7, 11, 15) only because App B is worse than trivial there.
 * Pair-by-pair (Lemma 1-D) thresholding is optimal: W + W' = 4h - 1 is forced by the pole order H = 2h. A perturbation can
   only win through cross-(x,v) relations; nothing at h <= 7 shows any, and the decisive regime (App B below trivial)
   starts at N = 176 (h = 15), out of reach of exact linear algebra (N^4 ~ 1e9).
 * No perturbation tested bends the frontier at theta = sqrt 2. Ideas 4-10 are o(1) or covered by barriers.
 * Most promising remaining: Test A at h = 15..20 restricted to a symmetric sector (e.g. y,u,x,v in a fixed small
   subspace, or the per-x operator system B_v^- rho + kappa B_v^+ = [v=x] Id solved with structured rho, kappa that commute
   with the weight-raising operators), looking for W + W' < 4h - 1 below saturation. First step: write rho, kappa as
   polynomials in the operators b[r] themselves (an ansatz of poly(h) unknowns) and test W + W' = 4h - 1 - k for k = 1, 2 at
   h <= 10; any k growing linearly in h would move theta.
