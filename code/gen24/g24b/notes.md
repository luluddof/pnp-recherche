# G24b — multi-channel, non-product constructions with xi on the Voronoi boundary (working notes)
Tags: [P] proved here (short, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [O] open.
Nothing here is a new algorithm or a new bound for 3-SAT.

## 0. Setting (G23c)
(A3) with K = 0:  sum_{g in S_xi} beta_g (x) omega_g = I_N (x) I_N,  beta_g = P B_g Q^T (pairs), omega_g = U O_g Y (labels).

## 1. General form of a solution [P]
 Let K = {kappa : sum kappa_g beta_g = 0} (relations among the beta_g), L = K^perp.  Then (A3) <=> there is c with
   P B_c Q^T = I, P B_k Q^T = 0 (k in K),   U O_kappa Y = <kappa,c> I (kappa in L)     (B_c := sum c_g B_g etc.)
 So N <= rk B_v for all v in c+K, and N <= rk O_kappa for all kappa in L with <kappa,c> = 1 (= G23c adversary form).

## 2. The key new bound: CROSS FLATTENING [P, elementary]
 The target delta_{xv} delta_{uy}, flattened (x,u)|(v,y), is I_{N^2}.  Hence for EVERY K=0 construction (any channels,
 any coefficients, any non-product choice):     N^2 <= rank F,   F := sum_{g in S_xi} B_g (x) O_g  (rows (a,m'), cols (b,m)).
 (same with O_g^T).  F couples pairs and labels; it is NOT bounded by sum_g rk B_g rk O_g cleverly — it is much smaller.
 Toy (all weights 0, A1^n): F factors through V_0^{(x)n} (2^n-dim)  =>  N <= 2^{n/2} = single-channel value exactly.
 Exact [E, one.py, blk.py]: one copy A1, xi = alpha/2: every charge block of F has FULL rank min(rows, cols);
   conserved charge q = charge(m') - charge(a) = charge(m) + charge(b) (per copy), cols(q) = rows(q-1)
   => rk F = sum_q min(r_q, r_{q-1}) = dimV' dimM' - max_q r_q.   e.g. (e,h,h) = (4,4,4): 138 (dims 14 x 14 = 196).
 CAUTION [E]: for n copies O_g (|g| >= 2) is NOT O_{g1} (x) O_{g2}: o(u1 (x) u2) = sum_k u1[k] (x) u2[-k] moves weight
   between copies (A1^2, E0 = 2: Lab_{(a,a)} = 11 > 3x3-product 9).  B_g IS a product.  Charge per copy is still conserved.
 n-copy bound [P]: rk F <= sum_{q in Z^n} min( [x^E0 y^Ha] prod R_{q_i}, [x^E0 y^Hb] prod C_{q_i} ),
   R_q(x,y) = sum_{e,ha,s} dimV^{(s)}_{ha} dimM^{(q+s)}_e x^e y^ha,  C_q = R_{q-1}.
   Checked against exact two-copy ranks (lattice A1+A1, xi = (1/2,1/2), two.py / pred2.py): bound = exact in 12 of 16 cases,
   e.g. (E0,Ha,Hb) = (4,4,4): 8470 = 8470; (3,3,3): 1640 = 1640; (2,2,2): 408 vs 410.
   Three copies [E, three.py -> three.out]: block bound >= exact everywhere, e.g. (E0,Ha,Hb) = (3,2,2): 8382 vs 8394.

## 3. Exact small-copy data  [E]  (two_ch.py, three.py; true per-channel ranks R_g, Lab_g, flattening rank)
 two copies, xi = (1/2,1/2), S = {00,01,10,11}:  (E0,h,h): R_g ; Lab_g ; max_g min(R_g,Lab_g) ; sqrt(rk F) ; loose
   (2,2,2): R = 11,16,16,28 ; Lab = 28,16,16,11 ; SC 16 ; 20.2 ; 28
   (3,3,3): R = 22,32,32,56 ; Lab = 56,32,32,22 ; SC 32 ; 40.5 ; 56
   (4,4,4): R = 52,74,74,124; Lab = 124,74,74,52; SC 74 ; 92.0 ; 124
 three copies: (2,2,2): SC 32, sqrt rkF 56.2, loose 96;  (3,2,2): SC 58, sqrt rkF 91.6, loose 96.
 Non-product effect (real, single channel): Lab_{(a,a)}(E0=2) = 11 > 9 = product of one-copy Sym's (weight transfer).
 At finite size sqrt(rk F) lies strictly between SC and loose: the exact data alone do not decide; the RATE does (sect. 4).
 True N is not computed (hard restriction problem); one-copy check (coset labels, E0 = 1, H = 0): flattening allows 2,
 the true maximum is 1 (case analysis: c_0 != 0 -> 2N <= 3; c_0 = 0 -> U Y = 0 & U O_a Y = I impossible; dependent -> R_0 = 1).

## 4. Rate level: flattening kills the dip  [P for the bound; N for the grid]   (rate2.py/rate3.py/rate4.py -> r3_*.out, r4.out)
 rk F^{(n)} <= sum_q min(X_q, Y_q) <= (H+1) [sum_q sqrt(R_q C_q)(x,y)]^n x^{-E0} y^{-H/2}   (Chernoff, valid for EVERY n)
 => 2 log N / n <= inf_{bx,by} [ log sum_q sqrt(R_q R_{q-1}) + bx u + by v ] + o(1)        (u = E0/n, v = H/(2n))
 "Bhattacharyya loss": the label-charge distribution and its shift by one channel unit overlap only partially at small
 weights per copy — exactly where the loose dip sits.  Combined criterion  min_{u,v} 2 s_M(u+v)/min(s_M(u), s_V(v),
 s_M(2u+2v)/2, flat(u,v)/2), grid u,v in [0.03, 40], b in [0.004, 30]:
   family                                 loose (old)          loose + flattening (new)
   A1^n, xi = alpha/2, labels V'          2.8160 (u=v=1.15)     2.8449 at grid edge u=v=40 (-> 2sqrt2 from above); at the
                                                               6 grid points where loose < 2sqrt2: >= 2.8864
   (A1+A1)^n, xi = (alpha1/2, 0)          2.7989 (u=v=0.33)     2.8535 at grid edge; at the 11 loose-dip points: >= 3.0224
   (A1+A1)^n, xi = (alpha1/2, alpha2/2)   = A1^n                 = A1^n
 So for these families NO K=0 construction — single- or multi-channel, product or not — reaches below 2sqrt2:
 the open case of G23c (gap 2.799 .. 3.030) is CLOSED for them [P bound + N grid].  The previous "loose" upper bound on N is
 not attainable: multi-channel N is at most sqrt(rk F), whose per-copy rate is strictly below the loose rate at the dip.

## 5. New caveat found: module labels with NEGATIVE lowest xi-weight (anti-dominance)  [N, O]
 M = V_{alpha/2+L}, xi = alpha/2: lowest xi-weight -1/4 per copy.  (A1) needs a grading E >= 0, so the upper anchor must use
 the SHIFTED grading (B = 2E0 + H with E0 measured from the bottom of M).  G23c's scan used the unshifted weights in
 rate(PC, 2u+2v): its anchor is too strong and its 2.817 for "A1 M=V_{a/2+L}" is NOT a valid loose value.
 Correct loose value: 2.4164 (u=0.116, v=0.03; pairs almost at weight 0, where V'_0 is 2-dim but M_0 is 1-dim:
 anti-dominance, escape route R2 of G22e realised by a regrading).  With flattening: 2.6740 — STILL below 2sqrt2.
 (Restricting to the old grid u >= 0.28: loose 2.526.)  The flattening is not tight there (one-copy check, sect. 3),
 so this is an upper bound on what such constructions could give, not a construction.  OPEN.
 Same caution applies to every scan row with a module whose lowest xi-weight is negative (A2, M = V_{w1+L}).

## 6. Task 2 (t-filtration) — status
 Not needed for the M = V' families: the cross flattening is a stronger, channel-free statement (it bounds every
 combination sum kappa_g B_g at once, because it never splits into channels).  For the coset regime I only have an
 H = 0 filtration bound by excitation count:  N <= min_k max( sum_{i<k} C(n,i), 2 dim M_E0^{(>= k/2 excited copies)} )
 (tau kills channels |g| >= k, or Id lies in their span and those O_g need >= k/2 excited copies on one side) [P sketch];
 numerically it does NOT beat the flattening there.  [O]
 Sanity (sanity.py): the Chernoff bound exceeds every exact rank (n = 1,2,3).

## 7. Verdict
 * Per-copy rate of the best K=0 construction in A1^n (labels V') and (A1+A1)^n at boundary xi is at most
   flat(u,v)/2 < loose rate at the dip; resulting 4rho >= 2.886 at the old dip (A1), >= 3.022 (A1+A1, xi = (a1/2,0)),
   global minimum on the grid at the edge, 2.845 / 2.854 at u = v = 40, decreasing to 2sqrt2 = 2.8284 from above.
   => no multi-channel / non-product gain below 2sqrt2: the G23c open case is closed for these families [P + N grid].
   The exact per-copy rate of the TRUE maximum N is not determined (it lies between the single-channel value and flat/2).
 * NOT closed: labels in a module with lower xi-ground than V' (coset V_{alpha/2+L}): correct loose 2.416,
   loose+flattening 2.674 < 2sqrt2.  Upper bound only; no construction.  Needs a stronger argument (or a construction).
Files: notes.md, objs.py (one-copy/n-copy B_g, O_g, exact flattening rank), one.py, blk.py, two.py, two_ch.py, three.py
 (+three.out), pred.py, pred2.py (block-count bound), chk.py, rate2.py/rate3.py/rate4.py (+r3_*.out, r4.out), sanity.py.
