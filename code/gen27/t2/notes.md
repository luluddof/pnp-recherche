# G27-T2 notes (constructive search N = 10-12 + trend analysis)
Tags [E] exact over Q (flint) + mod 2 large primes; [V] re-verified with p15 myverify.py (pure Python Fractions/Bareiss); [N] numerical.
Convention as p11: K[(x,y),(u,v)] = sum_l E_l[x,u] C_l[y,v]; p = rank(Swap - K), q = osr(K).

07:40 setup. lib.py (loader for all repo exact points + Q8 from g23a ST array), construct.py (tensor products E=kron, C=kron;
 index deletion = restriction to a coordinate subspace, f(N-1,q) <= f(N,q)).  cons/ = constructed points.
 Baselines: tensor products are weak (Theta_A (x) I_m scales pq by m^2); deletion gives (11,5)@6 EXACT from (14,5)@7 (p11 had it [N]).
07:45 IMPORTANT: numpy multithreading made VP 50x slower (54 s vs 1 s per run): always OMP_NUM_THREADS=1.
07:46 Border extension N -> N+1 (extvp.py; old block fixed, new row/col + extra terms free): generically f(N+1,q) <= f(N,q) + N
 numerically (cross cells cost N, 'perfect cross').  Exact block-diagonal version (extexact.py: E (+) e_l, C (+) c_l, needs
 M1 = sum c_l E_l, M2 = sum e_l C_l with M1 M2 = I and <e,c> = 1) finds nothing with {-1,0,1} coefficients for p11's points;
 for Q8 a term with E_l C_l = I exists, so (16,8)@9 follows.
07:48 descend.py from the (31,6)@10 extension: 29, 27, 26, 25 accepted (balanced-gauge polish gap > 1e9), 24 fails.
 exactify.py (sparsify L1 + diagonal log-fit + rounding + flint): (25,6)@10 -> EXACT (25,5)@10 (osr dropped to 5!), entries 0,+-1,
 rank mod 1000003 and 998244353 the same.  pq/N^2 = 1.25.
07:49 driver.py (ext + addq warm starts, descend, exactify) on N=10 q=3..8.  q=3,4 fail to descend from weak starts (45, 40).
07:55 (25,5)@10 block structure [E] (blocks2.py): MC (monomial, entries 0,+-1), Theta symmetric, tr Theta = 30:
 30 fixed cells + 20 triangle cliques (size 3, cost 1) + one 10-cell block of cost 5.  Same clique-block shape as N <= 9 optima.
 Verified [V] by p15's myverify (frac = bareiss = 25, osr 5) + mod 1000003, 2^31-1.
08:00 sa.py: simulated annealing directly in MC (state = q pairs of signed partial perms; cost = rank mod 1000003 of 1 - Theta).
 Calibration N=7 q=5 (opt 14 known): 150k steps reach only 19-24 (rank objective too flat).  Not competitive with VP; dropped.
08:05 driver N=9 q=7..12 (ext from (20,6)@8 / Q8).
 NOTE on the metrics: on the circle theta^2+T^2=4, min(theta+T) = 2 (at the axes) = the line's value, so min(theta+T) does NOT
 discriminate line vs circle; the discriminators are the balanced exponent log_N b(N) (line 1, circle sqrt2) and
 min(theta^2+T^2) (line 2 at theta=T=1, circle 4).  Endpoints (q=1: p=N(N-1)/2) are 'inside the circle' trivially (constant 1/2).
08:15 N=9 q=7..10 (driver): (22,7) [N] only (= known 22@6); (9,9): (20,9) [N] but construction A_3 (x) A_3 gives (17,9)@9 [E];
 (16,10)@9 [N] from border extension of Q8 (exactify > 400 s, killed).  (16,8)@9 border ext of Q8 fails numerically (tail 2e-5):
 Q8 classes are partial perms of size 4, so no term with E_l C_l = I.
08:20 lib.py: added A_k = Swap_k - Phi Phi^T (osr k(k-1)/2, p = 1) for k = 2,3,4.  Tensor products of A's / I(x)I give NEW exact
 points by pure construction:  (19,6)@8 = [I(x)I]_2 (x) A_4  (beats p11's (20,6)@8);  (17,9)@9 = A_3 (x) A_3;
 (24,18)@12 = A_3 (x) A_4, (24,18)@11 and (24,15)@10 by index deletion; (33,3)@9, (64,3)@12 (g24g's A_3 (x) I family).
 General: A_a (x) A_b at N = ab has p = a^2 + b^2 - 1, q = a(a-1)b(b-1)/4 [E for (3,3),(3,4)].
08:25 COORDINATOR request (G27-B2): corner theta -> 2: report delta_N(q) = 2 - log_N f(N,q) vs eta = log_N q -> corner.py.
 I also report delta' = log_N(f(N,1)/f(N,q)) which removes the trivial factor 2 of f(N,1) = N(N-1)/2.
08:30 (25,5)@10 does NOT extend with a perfect cross: border extension to N=11 needs p = 40 (cross cost 15 > N = 10) at q = 6
 (p = 35..38 tails 1e-3..5e-5).  The 9 -> 10 perfect cross of (22,6)@9 was special.  driver now tries a ladder of offsets.
08:35 balanced N=10: (26,10) from (17,9)@9-extension only; VP descents stall far above the G20b numerical 13-14.
08:37 (11,6) descent from the (40,6)@11 extension: (39,6) [N] accepted, 38 not reached (tails 2e-8 after 2500 its; rerun with
 8000 its killed by the deadline).
08:40 STOPPED (coordinator deadline). All processes killed.  verify_all.py (p15 myverify: Fraction + Bareiss, + ranks mod
 1000003 and 2^31-1) CONFIRMED: (25,5)@10 [search], (19,6)@8, (17,9)@9, (24,18)@12 [constructions], (11,5)@6 [deletion of (14,5)@7].
 Files exact/*.json.  Tables: table.out (frontier N=2..12, b(N), min theta+T, min pq/N^2, min theta^2+T^2), fits.out, corner.out.
## Verdict (honest)
 * Exact data N <= 10: min pq/N^2 over p,q >= 3 stays 0.94-1.25 (N=10: (25,5), pq/N^2 = 1.25): the small/moderate-q frontier hugs
   the LINE constant regime; nothing drifts away from the line where we have good data.
 * Balanced b(N): exact 4,6,7,8,8 (N=4..8); upper bounds 14,15,20,22 (ST, g24d) at N=9..12, my VP gives only 16 [N] at N=9.
   Fits of log_N b = alpha + c/lnN give alpha = 1.46..2.4 depending on the range: UNSTABLE, dominated by the Hurwitz
   anomaly at N=8 and by loose upper bounds at N >= 9.  No extrapolation is meaningful (uncertainty larger than line-vs-circle gap).
 * Inside the circle: every N has points with theta^2+T^2 < 4 (e.g. (8,8)@8: 2.0; (25,5)@10: 2.44), but only through constants
   (App B is trivial for N <= 256, so the circle itself is not realised at any computable N) -> no finite-size meaning.
 * Corner (G27-B2): delta' = log_N(f(N,1)/f(N,q)) vs eta = log_N q: linear (rms .105) and quadratic (rms .099) fits are
   indistinguishable; at finite N eta >= log_N 3 ~ 0.45, so the corner eta -> 0 is not sampled.  Values at q = N:
   delta' = 0.35 (6), 0.33 (7), 0.60 (8), 0.34 (9), 0.26 (10) >= App-B's eta^2/4 = 0.25; N = 11, 12 rows are construction-only (weak).
## Unfinished: N=11,12 small-q points (q=3..N) and balanced N=10..12 by VP; exactification of (16,10)@9 [N], (39,6)@11 [N].
