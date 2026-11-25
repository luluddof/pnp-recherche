# G25-P8: adversarial review of G25-P2 (explicit Theorem A') 

Files here: sweep_mine.py/.out (own sweep vs brute force), tw_linegraph.py / tw.out (exact treewidth of L(K_K)),
recur_mine.py (own recursion, paper (40), top levels), c1_mine.py, th_mine.py / th_426.out, th_1.out (thresholds).
None of p2's code was imported.

## (1) Sweep lemma (I2) -- CONFIRMED
Def. 3: (K,lam)-clique value = sum over one lam-bit string per edge of K_K of prod_w W_w(strings on the K-1 edges at w).
Step j of the node sweep enumerates C_j u N_j: j(K-j) + (K-j-1) edge strings, max over integer j = floor((K+1)^2/4)-1
= s(K) (the cut alone is j(K-j) <= floor(K^2/4) = memory; the extra K-j-1 are the new edges of node j+1).
- sweep_mine.out: 47 random integer instances (signed tables), K = 2..7, lam = 1..3: sweep = brute force; max edges
  enumerated = s(K) exactly; #terms <= K 2^{lam s(K)}.
- tw.out: exact treewidth of the line graph L(K_K) (primal graph of the edge variables) for K = 3..6 is s(K)-1, so no
  elimination order beats the sweep (matches Harvey-Wood's formula for tw(L(K_n))).
- Def. 4 / Lemma 8 only need the VALUE of the (K0,L)-clique with query access: sweep is a legal level-0 algorithm.
- Caveat (not an error): memory 2^{lam floor(K0^2/4)}; level 0 carries ~96% of the final exponent (E2 schedule), so
  the algorithm becomes exponential-space (~2^{0.95 x time-exponent}); the paper's on-demand recursion is not.

## (2) Accounting fixes (I3) vs paper text -- CONFIRMED
- Lemma 8 (ii): "let chi be uniformly random among the groupings in which every group has exactly g/K nodes"; Thm A(b):
  row setting = g(g^2-1)lam bits (Prop. 3). Entry = sum over settings of the g/K nodes of the group: 2^{(g/K)g(g^2-1)lam}
  < 2^{N/K}. The paper's 2^{2gP/K} is copied from Lemma 1 (groups <= 2g/K). (44): T 2^{2N/K+c lam}, c lam <= N/K -> 3/K.
  So e_{l+1} = e(rho+3zeta)+3/K can be replaced by e(rho_bar+1/g+z1+z2) + 1/K + c/g^4. Correct.
- Hoeffding: K group events (range P) + gK node-group events (range P/g); (39) gives Y < (U + z1 gP + z2 gP)/K^2.
  Split condition K e^{-2t z1^2} + tK^2 e^{-2t z2^2} < 1 is exactly the union bound. Correct (gain only x0.91).
- (41)-(42): U <= N/g + |R| S(lam); paper folds 1/g <= zeta. Correct.

## (3) Recurrence, c1, thresholds -- CONFIRMED (numbers slightly improvable)
Paper (40), x0=18: c1 = limsup 2 e_l sqrt(9 x_{l+1}) = 42.379 (x0 = 19, 20: 43.4, 44.4); row A 9*x_9 = 275616, row B 137430,
r=1: 8046 / 3780. Own c1 (geometric zeta schedules): C' 30.877, D 15.872, E1 10.581, E2 7.482 (p2: 30.879/15.854/10.569/7.474);
liminf/limsup = 0.7071 in every case.
Thresholds, r = 4.26, vs 2^n (own optimiser, fractional log t, no overhead terms -- shown negligible, < 2^-64 n):
 C 86526, C' 78467, D 20696, E1 9248, E2 4682   (p2: 87430, 80310, 21010, 9375, 4740).
With p2's own schedules my top-level formulas give 9374.5 / 4739.5 / 21009.6 exactly -> difference = schedule search only.
r=4.26 vs 1.6423^n: E1 18542, E2 9334; vs 1.30704^n: E1 52653, E2 26402. r=1 vs 2^n: D 270, E1 146.7, E2 106.1.

## (4) T4 and Lemma 8' -- CONFIRMED (leading order)
T4 covariances re-derived (same pair 1/K^2, shared first/second node 1/K^3, opposite roles <= 0), Var <= PU(1+2/K)/K^2,
Chebyshev + union over K^2/2 pairs as stated; g >= 6K^4/((3-1/r) zeta'^2) = 2.17 K^4/zeta'^2 at r=4.26 (notes say 2.1;
code uses exact). Lemma 8': Shannon 9 colours (Delta <= 6), Bernstein per class with s_i = d(mu_i + mu/9) gives
9 exp(-d^2 mu/(18(1+d/3))); cluster sizes Bernstein; formula correct; LV (verify + resample).

## (5) Unchecked dependencies on the paper
Prop. 3 construction ((33)-(36): rectangles, row/column networks, c_{a,a'} <= g^2 lam = P/g); Lemma 1 construction (26)
(only tiny end-to-end tests in p2); Hoeffding w/o replacement [Hoe63] (standard); Theorem B: identity sum_s A_s = h![x=v][y=u]
(Lemmas 9, 10, (53)), factor counts R^row, R^col (54)-(55), absolute constant c (evaluation 2^{c lam}, values 2^{(2lam)^c});
(46)'s C_HR is not explicit, p2's S(lam) uses the HR leading term (uncertified, effect ~1e-20 at lam >= 2^80);
section-6 network (W = 3m-n). The additive (I8) bookkeeping is a sketch, but its terms are < 2^-60 n at all thresholds.
