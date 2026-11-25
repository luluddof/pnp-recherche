# G26-P11 notes (spectral profile + targeted numerical search)
Tags [P] proof, [E] exact, [N] numeric, [H] heuristic.
03:25 spec11.py (spec11.out): spectra on the 5 exact optima [E-ish, float eig of integer matrices]:
 Theta = K S has INTEGER spectrum on all 5: {1 (mult m), -1, -2, 0, +2}; ||Ts|| = 2 (2.08 for (6,6)@5), ||K|| = 2..2.3, Q = q.
 tr Theta = sum tr(G_l G'_l) is LARGE: 7, 10, 12, 18, 12 (vs bound qN = 12,20,30,42,18).
 Ts: m eigenvalues >= 1 (m or m+few: 0,4,2,0,1 extra >1), rest (p-ish) negative in {-2,-1}, a few 0 / +-1.118.
 p/n = 0.31, 0.28, 0.24, 0.19, 0.33.  Count: m + sum(rest) = tr Theta, sum(rest) = -4,-4,-5,-11,-12.
appb_spec.py: Appendix B pieces thresholded (h=4,N=5; h=5,N=7): Q (via realigned SVD) 290..2e6, ||Ts|| 77..6e5,
 the non-1 spectrum of Ts is split into two ~equal clusters of huge + and - eigenvalues (e.g. h=5 s0=0: 22 below -10, 25 above +10):
 odd-trace cancellation by +- symmetry, as the design principle predicts. But at these N App B is trivial (p q ~ n or worse).
03:20 vp11.py (real VP, tail objective) + gauge.py (minimise sum ||XG_lY^-1||_F^2+||YH_lX^-1||_F^2 over X,Y; Theta similarity-invariant).
 Border-type solutions ((3,3)@4 tail 3e-3, (4,3)@4 tail 7e-5, (8,3)@5, (11,3)@6): Q ~ 1e3..1e5, Ts has a +-pair of huge
 eigenvalues (evmin ~ -evmax), BUT after gauge balancing Q drops to 6..13, ||Ts|| ~ 3-4 and the tail becomes O(10):
 the large norms are PURE GAUGE (cond(X)cond(Y) ~ 1e4); the gauge-invariant Theta has only 6 eigenvalues within 1e-2 of 1
 (need m = 12-13): border points get their small rank(1-Theta) from NON-NORMALITY (pseudospectrum), not from eigenvalue 1.
 => refined principle: Lemma 1-3 hold in every gauge, so an EXACT point needs large Q_min := min_gauge Q.  Border points at small
    N have small Q_min and are not exact; they say nothing about the needed intrinsic norm growth.
 Exact (5,3)@4 numerically re-found (tail 1e-16, gap 3e8); after gauge: Q = 3.23, ||Ts|| = 2 (integer optimum: Q=3).

## Necessary spectral profile [P-level arithmetic, from Lemmas 1-3 of P3; n = N^2, m = n-p, Q any gauge]
 Theta = 1 - W, W = E S (rank p).  Ts = 1 - (W+W^T)/2: exactly (n - p - r) eigenvalues = 1, r >= 0 eigenvalues > 1, where
 r = dim(V cap (im W + V^perp)) <= p (inertia of [[A,B],[B^T,0]]: >= r positive and >= r negative eigenvalues of the
 symmetric part of W).  Remaining <= p eigenvalues ("rest") satisfy, for odd j:  sum_rest lam^j = tr Ts^j - m - (extra >1)
 i.e.  sum_rest lam^j ~ -m + O(Q^j N).
 (a) Q >= N (i.e. T' := log_N Q >= 1): the j=1 constraint is vacuous (m <= n <= QN); higher j give at most Thm B, p >~ N^2 ln N/Q <= N ln N. In the whole region
     q >= N (above the diagonal theta <= 1 <= T) the spectral method says NOTHING, even with bounded norms.
 (b) Q < N: j=1 forces sum_rest ~ -n, so mean rest eigenvalue ~ -n/p = -N^{2-theta} and Q >= ||Ts|| >= n/p - 1: the LINE
     in terms of Q (log_N Q >= 2 - theta).  This is p >= (N^2 - QN)/(1+Q).
 (c) one-spike rest (all rest eigenvalues equal -c, the shape of all 5 optima: rest in {-2,-1}): j=3 gives
     p c^3 <= m + Q^3 N with p c ~ n  =>  log_N Q >= (5 - 2 theta)/3, strictly above the line for theta > 1
     (theta=1.2: 0.867 vs line 0.8; theta = 1.5: 0.667 vs 0.5), but still far inside the circle (1.6, 1.32).
     General j: log_N Q >= ((j-1)(2-theta)+1)/j; j=3 is the binding one for 1 < theta < 2.
 (d) to approach the line with theta > 1 the rest spectrum must be TWO-SIDED: mass ~2n/Q at ~ -Q/2 plus a few eigenvalues
     near +Q cancelling the cubes (e.g. p2 = p1 a^3/Q^3 at +Q).  Matching all odd moments up to degree k costs only log factors
     (P3 obstruction), so spectrally a point strictly inside the circle with Q ~ q^{1+o(1)} is NOT excluded: the needed profile is
       Ts: (n - p) eigenvalues >= 1 (almost all exactly 1), p eigenvalues spread two-sidedly over [-Q, Q] with sum ~ -n,
       odd moments ~ -n (up to Q^j N), ||Ts|| ~ Q >= N^{2-theta}; singular values of K: ||K|| >= ||Ts|| >= N^{2-theta}.
     The 5 optima instead sit in regime (a)/(b) border: Q = q, QN ~ n, ||Ts|| = 2, one-spike rest {-1,-2}: the j=1 count
     m + sum_rest = tr Theta with tr Theta = sum tr(G_l G'_l) as large as possible (7,10,12,18,12).
 App B pieces at computable N (h = 4, 5; N = 5, 7): all thresholds are trivial (p or q >= n - 3); the extreme ones have Q ~ 2e2..2e6,
 Ts with symmetric +- clusters of huge eigenvalues = the two-sided profile (d); X,Y gauge reduces Q only by ~25% (h=4, p=1 point:
 289 -> 214): there the norm growth is intrinsic, unlike the small-N border points.

## 03:30 numerical scans (scan.py = descending-p VP, real, random starts; EXACT? = tail<1e-12 & sv gap>1e5, then X,Y gauge
## balancing + polish: gap after polish reported).  [N] unless marked [E].
 N=6: q=3 p=12 (gap 5e14, Qmin=3.00), p=11 not found (20 restarts, best tail 4e-5: border-type as g24g).
      q=4 p=12, 11 not found (tail 9e-5).  q=5 p=11 (gap 2e14, trTheta=13.72 non-integer), p=10 not found (1e-4).
      q=6: p=10, 9, 8 found (gap 3.5e13); p=7 not found (20 restarts, 5e-5).
      (8,6)@6 EXACT [E]: sparsify.py (copy of g24g) + rounding + verify_exact.py: rank_Q(Swap-K)=8, osr=6, same mod 2,3,5,7.
      exact_6_8_6.json. Entries 0,+-1, monomial. Theta spectrum 1^28 (-2)^4 (-1)^4, Ts=Theta symmetric, ||K||=2, Q=6.
 N=7: q=3 p=20 (gap 3e14, Qmin 3.01, trTheta=13) [N]; q=4 p=19 (gap 9e13, trTheta 14) [N]; q=7 p=11 (gap 2.5e14, trTheta 21) [N].
 OBSERVATION [E]: Theta = K S is exactly SYMMETRIC for 4 of the 6 exact optima ((5,3)@4, (7,7)@6, (12,3)@6, (8,6)@6; not the two
 N=5 points, whose Theta is diagonalisable with eigenvalues {1,-1,-2} but non-normal).  Theta symmetric <=> sum G_l (x) G'_l =
 sum G'_l^T (x) G_l^T (family closed under (G,G') -> (G'^T,G^T)).  Spectrum of Theta in {1,-1,-2} (+ a 0 or +2): with a = #(-1),
 b = #(-2): p = a + b and n - tr Theta = 2a + 3b (+ corrections), i.e. p = (n - tr Theta - b)/2: maximise tr Theta = sum tr(G_l G'_l)
 and the number b of -2's (= number of singular values 2 of K).
 => symR.py: structured ansatz R = U diag(lam) U^T S (Theta symmetric), only n p + p unknowns; finds (8,6)@6 in ~1 s with
    eig(1-Theta) = {2^4, 3^4}.
 symscan.py (symmetric ansatz) at N=7: (19,3) best tail 0.86, (10,7) 0.21 after 30 restarts: WORSE than plain VP there
   (VP tails 4.5e-8 / 3.5e-7): the symmetric family is not where the N=7 near-solutions live.

## 03:35 Cell-block structure of the exact optima [E, blocks.py]
 Theta = sum_l P_l, each P_l a twisted signed partial permutation of the N^2 cells ((u,v) -> (sigma_l v, tau_l u)).
 Connected components of supp(Theta) on the cells: ALL 6 exact optima are block-diagonal with
   singletons Theta_ii = 1 (cost 0), blocks of size c = 2, 3 with Theta_B = -(J-I) up to signs, eigenvalues {1-c, 1^(c-1)}
   (cost 1 each: rank-1 defect), and a few size-4/6 blocks of cost 2/4.
   (5,3)@4: 7 fixed + 2 zero-singletons + 2x(c=2) + 1x(c=3);  (7,7)@6: 18 fixed + 3x(c=2) + 4x(c=3);
   (8,6)@6: 16 fixed + 4x(c=2) + 4x(c=3);  (12,3)@6: 12 fixed + 6x(c=2) + 2x(c=3) + one size-6 block of cost 4.
 So p = #blocks, n = #fixed + sum c_B.  Fixed cells of one twisted map <= N  =>  #fixed <= qN; a rank-1-defect block of size c needs
 c(c-1) off-diagonal +-1 entries and each P_l contributes <= c of them => c <= q+1.  Hence in this "block MC" model
   p >= (n - qN)/(q+1)  (the line; = Lemma 3 with j=1), and more generally a block of size c with defect rank r and +-1 entries has
   r >= (c - tr Theta_B)^2/||1-Theta_B||_F^2 >= (c - f_B)^2/(c(q+1)) (f_B = fixed diagonal cells in B).
 Max block size observed: 3 (any q <= 7).  To go beyond the line one needs blocks of size c >> q with low defect rank, which with
 +-1 entries is impossible by the rank >= tr^2/||.||_F^2 inequality -> needs LARGE entries (= the norm principle again, now per block).
 03:37 VP long runs (vp11.py, 15000 its): (17,5)@7 r2 tail 6.4e-15, gap 7.8e7, #ev(Ts)>=1 = 32 = m [N, exact-looking];
   (18,4)@7 tails 6e-8..9e-8, gap ~1e4, coef 2.5-4.4 (slow: unresolved);  (10,7)@7 tail 4.7e-8, gap 1.8e5, Q 454 (unresolved, border-like).
 scan.py N=8 q=8: rank 8 found (Q8-type), p=7 tail 10.0 (12 restarts): no.
 03:38 N=7 scans [N, all gaps >= 1e13 after polish]: q=6: p=15 (trTheta 17), 14, 13, 12 (trTheta 21) found; (20,3)@7 EXACT [E]
   (sparsify+round+flint: rank_Q 20, osr 3, same mod 2,3,5,7; exact_7_20_3.json).
 thspec.py: Theta eigenvalues of the numerical points (similarity invariant): mostly {1, -1, -2} (+ a few non-integers on
   non-sparsified points); trTheta always an integer (12,13,14,16,17,21).  (11,7)@7: algebraic mult of 1 = 39 but p = 11
   (geometric 38): a JORDAN BLOCK at eigenvalue 1 -- the would-be (10,7)@7 is only border (tail 4.7e-8, Q ~ 450).
 03:42 (17,5)@7: gauge-balanced polish gap 7.5e13, Qmin 22.8, Theta spectrum {1^32, (-1)^16, (-2)^1} [N]; sparsify running.
   (18,4)@7: after balanced polish sv_19 = 2.6e-5 (gap 7.6e4): unresolved/border.  (19,4)@7: sparsify failed to reach +-1
   (L1 = 43, rounding error 0.49): numerical only.  N=7 q=6 p=11 not found (tail 1.8e-5).  N=8 q=4: 24 (=Swap_2 (x) I family), 23 not found.

## 03:45 Clique-block model: combinatorial Theorem A with polynomial error  [P, short; unreviewed]
 Model: MC (signed partial perms) and Theta block-diagonal on cells with blocks = fixed singletons (Theta_ii = 1) or
 Theta_B = -D(J-I)D (size c >= 2, cost 1) -- the shape of 4 of the 6 exact optima' blocks (all of (7,7)@6, (8,6)@6).
 tr Theta = F (#fixed), tr Theta^3 = F - sum_B c_B(c_B-1)(c_B-2) (each block contributes its c(c-1)(c-2) signed triangles, all -1).
 MC: |tr Theta| <= qN, |tr Theta^3| <= q^3 N (odd words twisted).  Since c/2 - 1 <= c(c-1)(c-2)/12 for c >= 3:
     p = #blocks >= (n - F)/2 - sum_{c_B>=3} (c_B/2 - 1) >= n/2 - qN/2 - (2q^3+q)N/12  (with the 03:55 correction).
 (Checks: (7,7)@6 tr Theta^3 = 18 - 4*6 = -6 = tr Ts^3 in spec11.out.)  So in this model the strong-lemma shape p >= n/2 - O(q^3 N)
 holds; the general MC case (blocks with defect rank > 1, non-block Theta) remains Theorem A (exponential C(q,eps)).
 Mechanism: a block of size c >= 3 needs triangles = fixed cells of odd words; twisted odd words fix <= N cells each.
 Without the trace constraint the clique-block covering bound alone would allow p ~ n/(q+1) (c <= q+1): the LINE.
 03:47 (11,7)@7 EXACT [E]: sparsify + rounding (max rounding err 0.37, but the ROUNDED integer family is itself exact):
   rank_Q(Swap-K) = 11, osr 7, same mod 2,3,5,7; exact_7_11_7.json; monomial +-1.  Blocks: 20 fixed, 2x(c=2), 7x(c=3) and one
   4-cell block with eigenvalues {1,1,1,-2} but defect rank 2 (Jordan block at 1: costs 1 extra).  Theta spectrum 1^39 (-2)^8 (-1)^2.
 N=8 from random starts: q=6 p=18 not found (3e-5), q=3 p=27 not found (3.5e-4) -- VP from scratch weak at N=8 (as g24g).
 03:48 (12,6)@7 EXACT [E]: rank_Q 12, osr 6, mod 2,3,5,7 same; exact_7_12_6.json; Theta symmetric, trTheta 21; blocks 20 fixed,
   2x(c=2), 7x(c=3), one 4-block of cost 3.  (11,7)@7 has the same blocks except its 4-block costs 2 (Jordan) -> q 6->7 buys 1.
   (20,3)@7 blocks: 13 fixed, 4 zero-singletons, 10x(c=2), 2x(c=3), one 6-block cost 4 (the (12,3)@6 6-block again).
 Warm start Q8 (8,8)@8 embedded in N=9 (q=12, p=13): border-like (tail 1e-3, coef growing). N=9 q=12 p=15 random: not found (1e-3).
 03:52 Q8 (8,8)@8 (g23a/sol_Q8_8_8_8.json -> q8_8_8_8.json): Theta = 32 fixed cells + 8 clique blocks of size 4 (eigen -3), ||Ts|| = 3,
   Ts = Theta symmetric.  So ALL known exact optima are clique-block type with c in {2,3,4} (+ rare cost-2..4 blocks of size 4, 6).

## Clique-block model: exponent bound and design target  [P, short, unreviewed]
 Same model as above. Let B = blocks of size >= 3, k = |B|, C = sum_B c_B.  p >= (n - F - C)/2 + k, F <= qN, and since
 c(c-1)(c-2) >= (2/9) c^3 for c >= 3 and Hoelder (sum c^3 >= C^3/k^2):  C^3/k^2 <= 4.5 (2q^3 + q) N.  Hence
     p >= min_{0<=C<=n} [ (n - qN - C)/2 + C^{3/2} / sqrt(4.5 (2q^3+q) N) ].
 (Correction 03:55: singletons may also have Theta_ii = theta_s not in {0,1} (cost 1). Since sum_i |Theta_ii| <= sum_l fix(P_l) <= qN,
  F <= qN and sum_s |theta_s|^3 <= q^2 qN; then T_B := sum_B c(c-1)(c-2) <= F + sum theta_s^3 + q^3 N <= qN + 2q^3 N.)
 => if q^3 << N: p >= n/2 - O(q^3 N);  if q^3 >> N: p >~ N^{5/2} / (3 q^{3/2}),  i.e. in exponents  T >= (5 - 2 theta)/3  (theta >= 1).
 This EXCLUDES the line (T = 2 - theta) for theta > 1 inside this model, but leaves the whole region between (5-2theta)/3 and the
 circle open (theta = sqrt2: 0.724 < 1.414).  DESIGN TARGET (spectral profile made concrete): an MC family whose Theta is
 block-diagonal on the N^2 cells with ~n/c clique blocks Theta_B = -D(J-I)D of size c ~ q^{3/2} N^{-1/2} (needs c <= q+1 for covering,
 c^2 N <~ q^3 for triangles), rest fixed.  E.g. theta = 1.5: blocks of size N^{1/2}, q = N^{2/3}: would be (N^{1.5}, N^{0.67}), far
 inside the circle (T_circle(1.5) = 1.32).  Q8 at N = 8 is the instance c = 4, q = 8 (c^2 N = 128 <= q^3 = 512).
 Caveat: the ST model bound of g22b (p^2 q >~ N^4/8, i.e. T >= 4 - 2 theta) is stronger for theta < 7/4 where it applies; clique-block
 families inside ST cannot reach the target; a target family must use cancellations / non-ST structure.
 03:55 N=8: VP random (16,6) best tail 2e-3 (coef growth); warm from Q8 minus one term (rank 24 exactly when any term dropped)
   for (11,7): tail 5.4 (no); warm from (12,6)@7 embedded, (16,7)@8: tail 4e-8 but after gauge-balanced polish 7e-4: border/unresolved.
   log-tail objective (vp11 --mode log) on (11,3)@6: best 5e-5 = same as plain tail (border), no gain from the surrogate.
 04:05 (17,4)@7 EXACT [E]: the "(17,5)@7" VP solution had osr(K) = 4 (5th Schmidt value 0); re-expressed with 4 terms (SVD),
   balanced polish (gap 3.5e14), sparsify (weights 1.149/0.87/0.642/1.558 = diagonal gauge), NEW diaground.py (log-linear LSQ removal
   of diagonal + per-term scaling; residual 2e-5) -> all magnitudes 1 -> flint: rank_Q 17, osr 4, same mod 2,3,5,7. exact_7_17_4.json.
   So f(7,4) <= 17 (gain 4 below 21), better than the numerical (19,4).
 (11,5)@6: sparsify does not reach a sparse form (trTheta = 13.72 non-integer; diaground residual 2.5): numerical only [N],
   probably a genuinely non-monomial point (f(6,5) <= 11 [N] vs 12 [E]).  (17,4)@7 blocks: 14 fixed, 8x(c=2), 1x(c=3),
   4 four-blocks {1,1,-1,-1} of cost 2.
 04:00 N=7 q=5 (scan seed 32, 6000 its): p = 16, 15, 14 found (gaps >= 1.3e14 after polish; Schmidt rank exactly 5),
   p = 13 not found (25 restarts, 1.4e-3).  (14,5)@7 Theta spectrum {1^35, (-1)^10, (-2)^4} (integer): sparsify running.
   N=7 q=4 p=16 not found (25 restarts, 1.8e-3).  N=8 q=7: rank 19 found (gap 2.5e10), p=18 best tail 1.3e-7 (unresolved).

## TABLE (best p per q; E = exact over Q verified with flint + mod 2,3,5,7; N = numerical, sv gap >= 1e13 after balanced polish)
 N=6: q=3: 12 E (known) | q=4: 12 E (known) | q=5: 11 N (new, non-monomial?) | q=6: 8 E (NEW; was 9 N) | q=7: 7 E (known)
 N=7: q=3: 20 E (NEW; was 21) | q=4: 17 E (NEW) | q=5: 14 E (NEW) | q=6: 12 E (NEW) | q=7: 11 E (NEW) | q=8: 8 E (known, octonions)
 N=8: q=4: 24 (known family) | q=7: 19 N | q=8: 8 (Q8, known).  N=9: nothing below known (VP random/warm: border only).
 p q / N^2 of these: 1.0-1.6.  Circle/line at these N: everything satisfies theta + T > 2 (pq > N^2) except (12,3)@6 (= N^2);
 App B at these N is trivial (b_AppB = N^2 - O(N)), so "inside the circle" is not testable at N <= 9.
 04:05 (14,5)@7 EXACT [E] (sparsify 1 restart + diaground residual 1e-12 + flint; mod 2,3,5,7 same): exact_7_14_5.json.
 04:10 scan2.py (restarts, keep best 3, gauge-balanced long polish = exact/border discriminator):
   (20,6)@8: two restarts polish to sv_21 = 7e-15 / 2.6e-13 (gap 1.6e14 / 4e12), coef ~1, Schmidt rank 6 [N];
   (22,5)@8: two restarts sv_23 = 3e-15 / 1e-14 (gap 3e14) [N];  (18,7)@8, (11,6)@7: polish leaves sv_{p+1} ~ 2e-5..4e-3 = border/unresolved.
   Sparsify of (20,6)@8, (22,5)@8 running; scan2 (19,6)@8 and (21,5)@8 running.
 04:15 independent re-check (flint, own Swap/K construction) of all 6 new exact_*.json: ranks and osr confirmed, entries {0,+-1}.
 04:20 (22,5)@8 EXACT [E] (entries 0,+-1) and (20,6)@8 EXACT [E] (entries -1,0,1,2): sparsify (9 restarts) + diaground with threshold
   0.1..0.4 + flint, mod 2,3,5,7 same; exact_8_22_5.json, exact_8_20_6.json; independent re-check OK.  These are NOT clique-block:
   (22,5)@8 has a 30-cell block of cost 15 and a 9-cell block of cost 3; (20,6)@8 a 17-cell block of cost 8 and an 8-cell block of cost 3.
   (21,5)@8, (19,6)@8, (18,7)@8: polish leaves sv_{p+1} ~ 1e-4..3e-3 (border/unresolved).
 04:18 (19,3)@7, (16,4)@7 (scan2, 30 restarts): polish leaves sv_{p+1} ~ 3e-3..3e-2 with growing coefficients: not found / border.
 => N=7 exact frontier now (21,1),(21,2),(20,3),(17,4),(14,5),(12,6),(11,7),(8,8) (+ mirrors); next points down not found numerically.
 Data vs P3's fit p >= n/2 - qN/3: holds for all N=6,7,8 points EXCEPT Q8 (8,8)@8 (8 < 32 - 21.3); Q8 gives equality in
 p >= n/2 - 3qN/8.  So the constant in the "strong MC lemma" conjecture must be >= 3/8.
 04:22 table.py -> table.out: spectral table of all 14 exact points (5 given, 8 new, Q8).  N<=7 new points are all in MC
   (signed partial perms) with ||Ts|| = 2..2.11, Q = q, Theta eigenvalues in {1,-1,-2} (+ rare 0, 2, and the algebraic triple
   -2.11/0.25/1.86 of one 4-block); #eig(Ts) >= 1 = m exactly or m+1.  The two N=8 points as found are NOT in MC (rows with 2 entries;
   (20,6)@8 has one entry 2): Q = 8.6 / 7.3 > q, ||Ts|| = 2.34 / 2.59, non-integer Theta eigenvalues, big non-clique blocks.
   (Caveat: a sparser gauge might exist; sparsify is a local search.)  First exact optima outside MC in the repo, as far as I saw.

## SUMMARY (04:25)
 Spectral profile of optima [E]: Theta = KS ~ integer spectrum {1^(n-p), -1, -2 (-3 for Q8)}, Theta mostly symmetric, ||Ts|| = 2 (3),
 Q = q, sv(K) in {0,1,2}; #eig(Ts)>=1 = m (tight); remaining p eigenvalues negative (one-spike at -1/-2); combinatorially:
 fixed cells + clique blocks -D(J-I)D of size 2,3,4 (cost 1 each).  N = 8 optima found here are larger, non-clique, non-MC.
 Necessary profile for inside-circle: q >= N: no spectral constraint (j=1 vacuous). Q < N: rest mean -n/p, ||Ts|| >= n/p - 1;
 one-spike / clique-block shape => log_N Q (resp. T in the clique-block MC model, [P]) >= (5 - 2 theta)/3: beyond the line for
 theta > 1, still inside the circle; approaching the line needs a two-sided spread rest spectrum (as App B's +- clusters).
 Border solutions: large norms are pure gauge (balanced Q ~ q, tail O(10)); Theta non-normal / Jordan; use gauge-minimal Q.
 Numerics: VP (real, unregularised norms) + balanced-gauge polish + sparsify + diaground + flint => 8 new exact points
 (8,6)@6, (20,3)(17,4)(14,5)(12,6)(11,7)@7, (22,5)(20,6)@8; numerical (11,5)@6, (19,7)@8. Nothing beyond the line constant
 regime (pq/N^2 = 1.2-1.7); N = 9: nothing found (random and Q8-warm starts border-like).  log-tail surrogate: no gain.
 04:25 (27,6)@9 [N]: scan2 seed 92 r10 polished to rank 27 (sv_27=0.46, sv_28=1.3e-14, gap 3.6e13, Qmin 11.8, coef 1.2); improves the known (27,9)@9 (Swap_3 (x) I_3 family) and (33,3). Sparsify running; (25,6)@9 scan2 running. (16,9)@9: polish sv_17 = 1.7e-5 (border/unresolved).
 04:28 (25,6)@9 [N] (2 restarts, gap 2.8e13 / 6.1e13) and (24,6)@9 [N] (r9 re-polished at p=24: gap 2.8e13, Qmin 8.1, coef 0.85). scan2 (23,6),(24,6)@9 running.
 04:31 (22,6)@9 [N]: scan2 seed 94 found rank 22 when asked 23; balanced polish gap 1.7e14, Qmin 7.2, coef 0.78; Theta {1^59 (-2)^14 (-1)^5 2.56 -1.56 0} (trTheta 27). Gain n/2 - p = 18.5 = 0.34 qN (> qN/3). scan2 (20,6),(21,6)@9 running; sparsify running.
 04:33 (20,6)@9, (21,6)@9: polish leaves 3e-4..4e-2 with coef 20-34: border/not found.
 04:37 (25,5)@9 [N]: scan2 asked 26, three restarts polish to rank 25 (sv_25 ~1.3-1.6, sv_26 ~1e-15). (19,7)@9: border (3e-2). scan2 (23,5),(24,5)@9 running.
 04:40 (22,6)@9 EXACT [E]: sparsify (2 restarts) + diaground (th 0.4, residual 9e-12) + flint, mod 2,3,5,7 same, independent re-check OK, MC; exact_9_22_6.json. Blocks: 27 fixed, 5x(c=2), 12x(c=3), one 8-block cost 5. Previous best at N=9, q<=6: (33,3) and (27,9) families.
 04:43 (23,5)@9, (24,5)@9: polish leaves 7e-4..3e-2: border / not found.

## FINAL TABLE (04:45)  E = exact over Q (flint) + same ranks mod 2,3,5,7 + independent re-check; N = numerical (sv gap >= 1e13
## after gauge-balanced polish); "nf" = next p not found (border-type tails after balanced polish)
 N=6: q=3 12 E(known) nf 11 | q=4 12 E(known) nf 11 | q=5 11 N nf 10 | q=6 8 E NEW nf 7 | q=7 7 E(known)
 N=7: q=3 20 E NEW nf 19 | q=4 17 E NEW nf 16 | q=5 14 E NEW nf 13 | q=6 12 E NEW nf 11 | q=7 11 E NEW nf 10 | q=8 8 E(known)
 N=8: q=4 24 (known) nf 23 | q=5 22 E NEW (non-MC as found) nf 21 | q=6 20 E NEW (entries -1..2) nf 19 | q=7 19 N nf 18 | q=8 8 (Q8)
 N=9: q=5 25 N nf 24 | q=6 22 E NEW (MC) nf 21 | q=9 known 27 (Swap_3 (x) I_3), (14,14) ST; (16,9) border.
 pq/N^2: 1.0-1.7 everywhere (line-constant regime); no point is evidence for anything beyond the line or about the circle.
 04:47 (25,5)@9 sparsify (1 restart) not sparse: stays [N]. All my processes killed.
