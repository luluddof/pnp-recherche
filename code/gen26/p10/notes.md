# G26-P10 adversarial review of G26-P5 (notes as I go)
Paper text: /tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/papers/claude_3sat.txt (Prop. 3 at lines ~1000-1060).
## 03:10 analytic pass
- R2F: Doob increment re-derived (h(a) = c_a/K + Z_k(a) - c(a,S_k)/K): OK. Linear part sum (M_k - A_k) <= M_t - A_1 <= P: OK -> 2P/K.
  Covariance for (a,b) uniform distinct: E u(a)u'(b) = -(1/(K(K-1))) sum_a u(a)u'(a): OK (without-replacement handled exactly).
  Increment bound 4z_k + 3P/g: OK. Freedman stopped at predictable bad event: OK. Union over K(K-1)/2 unordered pairs (Y symmetric): OK.
  E1 families Z and Z' coincide by symmetry of c (so L = ln(4gK/d1) is conservative). McDiarmid for rms: OK. Popoviciu: OK.
- R3'': dX in interval of length 2(max-min)gamma/w <= 4G/w = P2's |dX| bound c -> Hoeffding-Azuma exp(-2s^2/sum c^2): Gam/2 OK.
  Serfling factor (1-(i-1)/(g-1)) on uniform i-subset of [g]-{v}: OK. A-term factor (1/2+1/g): OK.
- S2: Prop. 3 is rectangle-local (iota per rectangle, Lemma 2(a) per rectangle, (35) sum over rectangles; rows/cols chosen globally,
  min <= average still fine with per-rectangle lambda_q). BUT padding cost is miscounted by a factor 2:
  a rectangle has TWO diagonals, both padded to max(l1,l2): extra = 2max - (l1+l2) = |l1-l2|, not |l1-l2|/2.
  (In the paper's units each rectangle contributes 4 lambda to N: 2 strings per diagonal.) So padded total <= sum_e Y_e + sum_e |Y_e - c|
  and zeta_S2 = K^2 sd/N (+1/(K-1)) ~ 1.21/t, not 0.60/t.
- 03:08 P5's DP with corrected zeta_S2 (x2 on sd part), rowser base: c1 a=2 [8000,16000] = 3.312 (P5: 3.226).
- 03:09 S2 corrected thresholds (P5's mythr5, top S2c): r=4.26: 1514 / 2911 / 1.002e4 (P5: 1419/2747/9645).
- zlib10.py: own transcription of zeta_F, zeta_R2, sd_S2, R3/R3'' terms with exact discrete sums: matches P5 (P5 slightly
  conservative, e.g. Gam R3'' 0.0258 vs 0.0263 at g=2^8). Gam ratio R3''/R3 = 0.387 (Serfling alone 0.775): CONFIRMED.
- test_cov.py: covariance claim exact check (ratio sup -> 1.000, tight), sum Lambda <= 2P/K on exhaustive K=t=3: OK.
- Paper Lemma 1 (lines 833-905): child (K,L)-clique with L >= max Y and appended zero bits per edge: S2 = drop the appended bits,
  needs Def. 3/4 re-stated with per-edge widths and achievability measured by the PADDED total of the solver's own grid placement
  (induction is consistent: padding at level l is caused by the grouping at level l+1, which S2 controls for every placement).
  Second-order effects not in P5's model: P (max node sum) and c_ab <= P/g of the row network now use max widths -> factor 1+O(zeta)
  on zeta terms; negligible asymptotically.
- Counterexample to P5's padding inequality "padded total <= sum_e Y_e + (1/2) sum_e |Y_e - c|": take every rectangle with
  diagonal widths c+delta, c-delta. Padded (unordered) total = sum_q 2(c+delta) = sum_e Y_e + 2 delta |R|, while P5's bound gives
  sum_e Y_e + (1/2)(2 delta |R|) = sum_e Y_e + delta |R|. Correct: padded total <= sum_e Y_e + sum_e |Y_e - c| (tight here).
  -> zeta_S2 = K^2 sd/N + 1/(K-1) ~ 1.21/t. (Possible improvement, unproved: E|l1-l2| <= sqrt(2Var(1-corr)) ~ sqrt2 sd if the
  covariance of disjoint group pairs is O(Var/K) -> ~0.85/t.)
- 03:12 machine heavily loaded (load ~19 on 4 cores, other agents). dp10.py rewritten: closed-form level max (own derivation:
  R cos(th+phi) = eb be - ea al), vectorised over j; rowsolver on a 20001-point theta grid. Runs: R3pp+F (claim 3.824), R3+R2 (calib 4.050).
- 03:14 dp10 (independent DP + own zetas): R3+R2 = 4.050 (calib exact), R3''+R2F = 3.824 (P5 3.824): CONFIRMED.
- 03:18 dp10: R3+R2F 3.921 (P5 3.921), R3''+S2(P5's 0.5 factor) 3.226 (P5 3.226) -> P5's numerics reproduce exactly;
  R3''+S2 corrected (factor 1) 3.312; poly+R2F 5.875 (P5 5.875). Thresholds (P2/P5 mythr5 on MY solver table, top R2F):
  1943 / 3821 / 1.341e4 | r=1: 59.28 / 92.95 / 354.9 = P5 exactly.
- 03:22 S2 corrected (zeta = K^2 sd/N + 1/(K-1)), own DP: R3 base 3.433 (P5 3.352), R3'' 3.312 (P5 3.226), poly 4.853 (P5 4.611);
  thresholds (top S2c) r=4.26: 1514 / 2911 / 1.002e4 (P5 1419/2747/9645); r=1: <=50 (bisection floor) / 82.8 / 279.6.
  Corrected S2+R3'' (3.31) is now only 0.05 below floor B (3.36).
## Floors (4), brief
- F3 (zeta >= 1/t - 2/K): checked analytically, OK.
- F2 (first moment, Bernoulli network): union over K^g groupings x independent Bin(t^2,nu) pair sums: OK as an asymptotic statement;
  caveat: normal approximation of the binomial upper tail is NOT conservative for nu > 1/2 (upper tail lighter), so 3.36 is a
  heuristic floor, scope = uniform-width lemmas of the stated network-independent form.
- F1 (random stratified grouping >= ln K/t): true for typical groupings (constant-probability analyses); the algorithm takes the
  BEST grouping by exhaustive search, and on the rank-2 network perfectly balanced stratified groupings exist (excess 0), so F1
  is a floor of analyses, not of the algorithm.
## Verdicts
(1) R2F CONFIRMED. (2) R3'' CONFIRMED. (3) c1 3.921 / 3.824, poly 5.875, thresholds 1943/3821/1.341e4 CONFIRMED (independent DP).
(4) CONFIRMED with caveats (heuristic/asymptotic; F1 is about analyses of typical groupings). (5) S2 CONFIRMED WITH CORRECTION:
padding factor 2 -> zeta ~1.21/t, c1 3.31 (R3''), thresholds 1514/2911/1.002e4; remains a sketch (per-edge Def. 3/4, Lemma 1/8 with
max-width P, second-order).
No process left running.
