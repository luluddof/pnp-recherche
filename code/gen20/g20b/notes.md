# G20b — Appendix B at finite N, and the balanced frontier b(N) beyond the Hurwitz range  (working notes, saved as I go)

Tags: [P] proved here (short proof, unreviewed), [E] exact computation, [N] numerical only, [H] heuristic, [O] open.
Setting and notation as G19b (X = [x=v][y=u]; (p,q) feasible: X = R + C, rr(R) <= p, rc(C) <= q; b(N) = least p with (p,p)).
Nothing here is a new algorithm or a new bound for 3-SAT.

## 1. Appendix B at finite N  (task 1)

### 1a. Exact counting formula  [P, checked exactly for h <= 6]
In (54) the (u,v)-factor [t^mu(u)] b[h-|sigma|] t^sigma depends on sigma only (not on (s,i) separately), and in (55) the
(y,v)-factor [t^sigma] b[h-s-i] t^mu(y) likewise.  Hence thresholding the pieces at s0 (R = sum_{s<s0} A_s, C = d X - R) gives
     rr(R) <= P(2h + s0 - 1),     rc(C) <= P(2h - s0),      P(w) := p(0) + p(1) + ... + p(w)  (cumulative partition numbers),
for N = p(h) (all partitions of h) and, by restriction, for every N <= p(h).  Also both <= N^2 trivially.
Balanced value b_AppB(N) = min_{s0} max(P(2h+s0-1), P(2h-s0)) ~ P(2h).

### 1b. Exact implementation  [E]  (appb_exact.py; out_appb_h*.txt; appb_exact_h*_N*.json)
Implemented (47) (via the recursion (O2)), (48), Lemma 10's a_mu, b_nu and (53) literally, with integer arithmetic.
Check: sum_s A_s = h! X holds EXACTLY for h = 2..7 (N = 2, 3, 5, 7, 11, 15; also prefixes N = 4 of h = 4, N = 8 of h = 6).
(This is an independent machine check of Lemma 9(c) + Lemma 10 for these h.)
Exact ranks over Q (python-flint) of every threshold decomposition: they coincide with the counting bound (min with N^2)
up to a few units (e.g. h = 6: rr(R_{s0}) = 1,2,4,7,12,19,29,42,62,89,116 vs P = 1,2,4,7,12,19,30,45,67,97,121).
Exact App B threshold frontiers (balanced value = min_s0 max(rr, rc)):
   N=2 (h=2): balanced 3 (N^2 = 4)          N=3 (h=3): 8 (9)         N=4 (prefix of h=4): 13 (16)
   N=5 (h=4): 22 (25)                       N=7 (h=5): 46 (49)       N=8 (prefix of h=6): 58 (64)       N=11 (h=6): 111 (121)       N=15 (h=7): 215 (225)
   i.e. at every size where b(N) is known, Appendix B is within a few units of the TRIVIAL point (N^2, 0)/(0, N^2).

### 1c. Table: N / b_AppB / exponent  (counting formula, which the exact ranks confirm where computed)  [E]
 N = 2^lam:  lam  N      h   b_AppB (capped at N^2)   log b / log N    uncapped P(2h): exponent
             1    2      2   4                         2.000            3.585
             2    4      4   16                        2.000            3.033
             3    8      6   64                        2.000            2.696
             4    16     8   256                       2.000            2.459
             5    32    10   1024                      2.000            2.281
             6    64    12   4096                      2.000            2.140
             7    128   14   16384                     2.000            2.025
             8    256   17   65536                     2.000            2.002
             9    512   20   215308                    1.968
            12    4096  29   4.84e6                    1.851
            16    65536 44   3.59e8                    1.776
            24    1.7e7 81   1.40e12                   1.681
            32    4.3e9 128  4.86e15                   1.628
            64    1.8e19 417 5.9e29                    1.545
 N = p(h): first N where App B beats the trivial balanced value N^2 at all: h = 15, N = 176 (exponent 1.985).
 exponent < 1.9 from N ~ 800;  < 1.8 from N ~ 1.2e4;  < 1.7 from N ~ 2.3e6;  < 1.6 from N ~ 5.7e11;  < 1.5: not reached for
 N <= 4.6e23 (h <= 600).  The approach to sqrt 2 is like sqrt 2 + O(log h / sqrt h) (the poly(h) factors in P(2h) vs p(h)).
 => The "circle regime" of Appendix B does not exist at any computable N; at N <= 12 App B gives b_AppB = N^2 - O(N),
    exponent 2.  The fair baseline for G19b's numbers is therefore the TRIVIAL point, not the circle.

## 2. Which single point would move rho  [P, short]
Claim. Suppose a uniform family (Definition 7-type: integer entries of poly(lambda) bits, computable) gives (N^a, N^b)
feasible for N = 2^lambda with a, b <= beta < sqrt 2.  Then, together with Appendix B, pieces exist with
     4 rho <= beta + sqrt(4 - beta^2) + o(1) < 2 sqrt 2,   e.g. beta = 1 (P(N) = N^{1+o(1)}) gives rho <= (1+sqrt 3)/4 = 0.683.
Proof. G18g Theorem 1(ii) (telescoping pieces A_j = R_j - R_{j-1}, rr rc <= (p_j+p_{j-1})(q_j+q_{j-1})) applied to the chain:
App B threshold points with theta < a (their T > sqrt(4-a^2) >= sqrt 2 > b, so the chain is monotone), then (N^a, N^b), then App B
points with T < b (theta > sqrt(4-b^2) > a).  A consecutive pair of circle points costs exponent theta_j + T(theta_{j-1}) =
theta + sqrt(4-theta^2) + O(delta), which increases on [0, sqrt 2], so it is <= a + sqrt(4-a^2) + O(delta) (resp. b + sqrt(4-b^2));
the two pairs adjacent to the new point cost max(a, a-delta) + max(b, T(a-delta)) = a + sqrt(4-a^2) + O(delta) and the mirror.
So 4 rho <= max(a + sqrt(4-a^2), b + sqrt(4-b^2)) + O(delta) <= beta + sqrt(4-beta^2) + O(delta).  O(1/delta) pieces.  QED
So "rho < 1/sqrt 2" FOLLOWS from (not: is equivalent to) a uniform "P(N) <= N^{sqrt2 - eps} for N = 2^lambda"
— i.e. beating Appendix B's own bound P(N) <= N^{sqrt 2 + o(1)} (already better than Hrubes' N^{1.62}) at the balanced point.
(G18g's "P(N) <= N^c gives only rho <= (2+c)/4" used the trivial curve around the point; with App B around it the
bound is (c + sqrt(4-c^2))/4.)  A single fixed N (a "small-case gain") does not help: G19b section 6 (gains do not tensor).

## 3. Balanced frontier b(N), N = 7..12  (task 2)
(sections 3a-3c filled below as results arrive)
### 3a. Calibration of vp.py (variable projection, L-BFGS, random restarts; 4 shared cores, other agents' jobs running)
 N=6 (7,7): found, tail 3e-26, 4 restarts.   N=7 (8,8): found, 2 restarts.   N=8 (8,8) (Hurwitz): found, 1 restart.
 For N >= 9 there is NO known feasible point near the balanced region to calibrate on (App B is trivial there, sect. 1;
 direct sum 8+1 gives only (16,17) at N=9).  Calibration run (16,17)@9: see 3b.  Failures at N >= 9 are therefore weak evidence.
### 3b. Results  [N] (out_vp*.txt)
 N=7: (7,7) tail 6.01 after 30 restarts (all restarts at 6.0x) -> not found;  (8,8) FEASIBLE.
 N=8: (7,7) tail 18.0 after 22 restarts -> not found;  (8,8) FEASIBLE.
 N=9: (9,9) tail 13.0 (14 restarts); (10,10) tail 9.07 (7); (11,11) tail 4.47 (6 restarts) -> not found.
 N=10: (10,10) tail 20.0 (10 restarts) -> not found.
 N=12: (13,13) tail 40.0 (2 restarts only) -> not found.
### 3c. Structured ansatz (msearch.py)  [E when found]
 Observation [E]: all of G19b's exact points have C(W) = sum_l G_l W G'_l^T with G_l a signed PARTIAL PERMUTATION matrix and
 G'_l (almost always) = G_l^T; the supports M_l are disjoint and cover a set S of pairs (x,y) (C(x,y,y,x) = 1 on S, 0 elsewhere,
 |S| ~ N^2/2: 7/16, 12/25, 18/36), and R = X - C has rank p although it has N^2 - |S| + garbage nonzero rows.
 In operator form: W^T - sum_l G_l W G_l (G_l signed partial permutations) has rank <= p.
 Simulated annealing over (partition of [N]^2 into partial permutations, signs), energy = rank + tail:
   recovers (5,3)@4 and (4,4)@4 in seconds (verified exactly), but NOT (6,6)@5 in 150 s (best rank 7), nor (8,8)@8 (rank 15):
   too weak to extend the frontier; no family extracted.
### 3d. Calibration at N = 9 (important)  [N]
 Known EXACT point (17,16)@9 (direct sum 8+1, see 3f).  Cold vp.py (maxit 4000, 8 restarts): best tail 1.7e-6 only ("border?").
 => cold vp.py at maxit 4000 does NOT reliably find feasible points at N = 9: all cold "not found" at N >= 9 are weak evidence.
 Warm start (vp_init.py: the numerical (8,8)@8 solution zero-padded + noise, extra terms random, maxit 20000) is much better:
 (14,14)@9 found in the FIRST restart (tail 9e-25, |coeff| <= 1.8).
### 3e. Results at N = 9, 10, 12  [N]
 N=9 : (9,9) 13.0 | (10,10) 9.07 | (11,11) 4.47 (cold)  | (12,12): cold best 8.4e-4 (6 restarts, maxit 20000), warm best 2.4e-4
       (11 restarts; coefficients 12-21, drifting) -> border candidate / unresolved | (14,14) FEASIBLE (warm, 1 restart).
 N=10: (10,10) 20.0 | (11,11) 14.6 (11 restarts) | (12,12) 0.056 (9) | (13,13) 9.5e-4 (11 restarts, coeff 4) -> unresolved |
       (14,14) FEASIBLE (cold, tail 5e-27, max|coeff| 1.0).
 N=12: (13,13) 40.0 (2 restarts) | (14,14) 29.8 (4) | (16,16) 12.3 (5, cold) | warm/cold (16,16),(18,18): see 3g.
### 3f. Exact upper bounds by direct sums  [P, elementary]
 Split [N] = A u B.  X = X_AA + X_BB + X_AB + X_BA (by the blocks of x=v and y=u).  X_AB (x,v in A; y,u in B) is a restriction
 of X_{max(|A|,|B|)}, so costs any pair feasible there, or (|A||B|, 0) / (0, |A||B|) (it has rr = rc = |A||B|).
 N = 8+1: (8,8) + (1,0) + (8,0) + (0,8) = (17,16): b(9) <= 17.   N = 8+2: (8,8)+(1,1)+(8,8)+(8,8) = (25,25): b(10) <= 25.
 N = 8+3: (27,27);  N = 8+4: (8,8)+(4,4)+(8,8)+(8,8) = (28,28): b(12) <= 28.  N = 8+8: b(16) <= 32.  (Exponent 2 in general.)
### 3g. Longer runs / polishing (polish.py: L-BFGS continued from the best point)  [N]
 (12,12)@9: polished 560 s: tail 2.4e-4 -> 1.39e-5 while max|coeff| grows 12 -> 18.9 (tail ~ coeff^-k): BORDER-type signature, unresolved.
 (13,13)@9: warm restarts 9.8e-6 (coeff 7.1) and 3.8e-4; polishing: 3.8e-4 -> 7.7e-5 with coeff 12 -> 14.6: border-type, unresolved.
 (16,17)@9 (KNOWN exactly feasible) cold, maxit 20000: tail 3.1e-8 (coeff 2.7): so even a feasible point can look like 1e-8 here.
 (13,13)@10: cold 9.5e-4 (coeff 4); warm 4.07 (1 restart): unresolved.   (14,14)@10 FEASIBLE (cold).
 N=12: (16,16) cold 12.3 (5 restarts), warm from (14,14)@10: 8.2 / 8.3 / 16.8 / 20.5 -> not found; (18,18) warm 0.23-0.52 (2 restarts);
       (20,20) and (16,16)@11: see 3h.
 Exactification: sparsify.py + torus.py (G19b pipeline) on (14,14)@9 and (14,14)@10 and (8,8)@7: one gauge-L1 restart each
 (L1 = 104, 121, 54) was not sparse enough; rounding gives rr(R) = 46..64 instead of 14: NOT exactified in the time available.
 N=12: (18,18) warm best 3.6e-3 (4 restarts); (20,20) warm best 1.4e-3 (3 restarts) -> not found.  N=11: (16,16) warm tails ~2.6-3.5 -> not found
 (warm start from a q=14 solution; cold runs not done).  So for N >= 11 the numerics give NOTHING better than the exact direct sums.

## 4. Fit: best known b(N) against line (exponent 1) and circle (sqrt 2)  (task 3)
 N  | b(N) best upper (status)                         | proved lower | log b/log N | b/N  | N^sqrt2 | App B (exact/count)
 2  | 1  [P, settled]                                  | 1            | 0           | 0.5  | 2.7     | 3
 3  | 3  [P, settled]                                  | 3            | 1.00        | 1.0  | 4.7     | 8
 4  | 4  [E] (border 3)                                | 3            | 1.00        | 1.0  | 7.1     | 13
 5  | 6  [E, G19b]                                     | 3            | 1.11        | 1.2  | 9.7     | 22
 6  | 7  [E, G19b]                                     | 3            | 1.09        | 1.17 | 12.6    | (not computed; <= 46)
 7  | 8  [restriction of the 8-square identity; N]; (7,7) not found, 30 restarts, tail 6.0 | 4 | 1.07 | 1.14 | 15.7 | 46
 8  | 8  [Hurwitz]; (7,7) not found, 22 restarts, tail 18.0 | 4       | 1.00        | 1.0  | 18.9    | 58
 9  | 14 [N]; 17 [P, direct sum]; (12,12),(13,13) border-type unresolved; <= 11 not found (cold, weak) | 5 | 1.20 (exact: 1.29) | 1.56 | 22.4 | 81 = N^2
 10 | 14 [N]; 25 [P, direct sum]; (13,13) unresolved; <= 12 not found (weak) | 5 | 1.15 (exact: 1.40) | 1.4 | 26.0 | 100 = N^2
 11 | <= b(12)                                         | 6            |             |      |         | 111
 12 | 28 [P, direct sum 8+4]; numerics: (20,20) not found (tail 1.4e-3) | 6 | 1.34 | 2.33 | 33.6 | 144 = N^2
 (proved lower bounds: ceil(N/2) from Theorem 2 of G18g / Theorem 5 of G19b at p = q; Theorems 3,4 of G19b for N <= 4.)
 Reading.  (i) Up to N = 8 the best b(N) is N + O(1) (b/N in [1, 1.2]); at N = 9 there is a JUMP (numerically b(9) >= 12,
 border-type at 12 and 13, clean feasibility only at 14): b/N goes from 1.0 to ~1.4-1.6 between N = 8 and N = 9-10.  This is the
 Hurwitz pattern (8 is special), as for sums of squares, but here only numerically.
 (ii) log b/log N at N <= 12 says nothing about the asymptotic exponent: e.g. the direct-sum family (exponent 2 asymptotically)
 gives b(16) <= 32 = N^{1.25}, below N^{sqrt 2}; and Appendix B (asymptotic exponent sqrt 2) gives N^2 - O(N) here.  Constants
 dominate completely: "line versus circle" is still undecidable from these sizes.
 (iii) Every N in this range is far below App B (which is trivial up to N ~ 176 and has exponent < 1.6 only from N ~ 10^12).

## 5. Families (task 3)
 * Matching ansatz (3c): the common shape of all exact small points (signed partial permutations, C(W) = sum_l G_l W G_l^T-type);
   my annealing search did not even recover (6,6)@5, so no family was extracted.  [O]
 * Direct sums / restriction (3f): exact, valid for all N, exponent 2 (b(2N) <= 4 b(N) at best): useless.
 * Kronecker powers: excluded (G19b section 6, mixed terms).
 * No family with exponent < sqrt 2 found.  None is claimed.

## 6. Verdict (task 4)
 * Appendix B at finite N is trivial: exact ranks (implemented and checked, h <= 7) equal the cumulative-partition counting
   bound, and its balanced value is N^2 for every N = 2^lam <= 256 by the count (exact ranks: N^2 - O(N), e.g. 58 at N = 8) and exponent > 1.9 up to N ~ 800, > 1.6 up to
   N ~ 5.7e11, > 1.5 beyond N = 4.6e23.  So G19b's small-case points (b ~ N) are NOT evidence against the circle: App B's circle is
   a statement about N >= 10^20 or so.  The fair finite-N baseline is the trivial point.
 * What would move rho (section 2): a uniform family with b(N) <= N^{beta}, beta < sqrt 2 (equivalently P(N) <= N^{sqrt2 - eps}),
   gives 4 rho <= beta + sqrt(4 - beta^2) < 2 sqrt 2 when combined with App B.  That is: beat App B's own composition-type bound.
   Nothing here does.
 * Data N <= 12: b(N) = N + O(1) up to N = 8, then a jump at N = 9 (numerical).  Too small, too dominated by constants and by the
   1,2,4,8 Hurwitz structure to say line / circle / between.  Honest status: undecided.
 * Next computation: (a) exactify (14,14)@9 and (14,14)@10 (more gauge-L1 restarts, or the matching ansatz seeded from them) to get
   EXACT points beyond the Hurwitz range; (b) decide (12,12)@9 / (13,13)@9 (border vs exact) with a border-aware method
   (degree-1 degeneration search: C(eps) with R(eps) rank p, Laurent in eps); (c) N = 16 with warm starts from 8 (+) 8 block
   structure: b(16) vs 2N = 32 is the first size where a multiplicative gain over direct sums could show; (d) asymptotically, the
   only meaningful test is a uniform construction, not more small N.

## 7. Adversarial self-review
 * App B code: an independent implementation of (47),(48),(53) from the text; the identity sum_s A_s = h! X came out exactly for
   h <= 7, which would be very unlikely with an implementation error; ranks computed exactly over Z (flint).  The counting formula
   is proved from (54)/(55) directly and is an UPPER bound; the exact ranks confirm it is nearly tight (differences <= 5% where
   uncapped).  Thresholding is the only way of using App B considered; arbitrary subsets of the 4h+1 pieces were not optimised
   (cannot beat ~N^2 here since every single piece with s in [-1,2] has rr = rc ~ N^2).
 * "Balanced" App B value = min over thresholds of max(rr, rc); it includes the near-trivial endpoints (1, N^2 - small).
 * Section 2's claim uses G18g Theorem 1(ii) and the App B circle for the other points; uniformity/bit-size of the hypothetical family
   is assumed, not discussed further.  The direction "only P(N) <= N^{sqrt2-eps} can help" is NOT claimed (other off-balance
   improvements of the circle also help; the balanced point is just where App B is worst).
 * All N >= 7 feasibility statements except the direct sums and Hurwitz are NUMERICAL (tails < 1e-20 with |coeff| <= 2: not exact).
   No decomposition at N >= 7 was verified exactly in this run.  "Not found" at N >= 9 is weak (calibration failure on (16,17)@9).
 * The exact points re-verified here: (5,3)@4 and (4,4)@4 (matching ansatz, verify.py: rank factorisation re-expanded over Q).
 * All background jobs were stopped before the report (17:10 UTC).
 * failed_rounding_exact_*.json: valid but useless exact decompositions from failed roundings (rr 15-64), kept only as a record.
