# G26-P2 adversarial review of G26-B7 (notes as I go)
## 01:55 start. Read b7 notes + scripts, p2/p8/p9 notes, p2 thmA_explicit.py.
### Analytic re-derivations (first pass)
- R2: Doob decomposition re-derived: D_k = [h(a_k)-avg h] + [h'(b_k)-avg h'] + [c(a_k,b_k)-E], h(a) = (c_a - c(a,S_k-a))/K + Z_k(a). OK.
  linear part |.| <= (r_k + P/t)/K, sum r_k <= P (telescoping over sorted consecutive strata). OK. Linear parts are
  mean-zero given F_{k-1} (a_k uniform in S_k) so the remainder is still a MDS. Z_k F_{k-1}-measurable; good event
  predictable -> stopped Azuma valid. Union for Z: 2gK events x d1/g^2 = 2 d1/t <= d1. Pairs: Y symmetric -> K(K-1)/2. OK.
- R3: ANOVA identity re-derived, drift re-derived (E Gamma_{j+1} = Gamma_j (1-2/(g-j))), |dX| <= 4G/w_{j+1} re-derived.
  Padding OK provided main term uses an upper bound U_bd >= nu N (nu = rho <= rho_bar: fine).
  A-term: my constant is 4x smaller than b7's (b7 conservative). 
  Probability budget: dA + dG(forward) + dG(reverse order needs its own G events!) + dGam/2 = 1/3+1/3+1/3+1/6 = 7/6 > 1. Minor (log-factor) slip.
- R1: Phi = exact cond. expectation (matching => disjoint label sets). Min over g labels <= average => nonincreasing. OK.
  Caveat: exp(s) real -> need finite precision (O(log m) bits, rounding slack from Phi(0) <= 1/2); cost O(m g^3) = poly.
- R4: continuous derivation re-done: M = 1/4 - 3 d^2/4, T = 1/4 + d^2/4 => T = (1-M)/3 for d in [0,1/2] (M in [1/16,1/4]);
  recursion with S = middle cut: M = M_h/4, T = 1/4 + T_h/4 -> same line on [0,1/16]. OK.
## Numerics (own code, nothing imported from b7)
- test_r2.py: exact Doob martingale (K=t=3, all 6^3 bijection tuples): D_k formula exact (3e-16); sum|lin| <= 4P/K; E Y <= U/(K(K-1)).
- test_r2b.py: K,t in {(16,64),(8,256),(32,32)}, 5 network kinds incl. star/ladder/blocks: stratified excess always
  <= R2 bound (bound loose: 0.47/0.11/1.04 per N vs observed <= 0.02); plain grouping much worse on star/blocks.
- test_r3.py: ANOVA identity 1e-14, drift 1e-15, |dX| <= 4G/w_{j+1} holds (max ratio 0.83) on 5 kinds.
- test_r1.py/test_r1b.py: toy derandomised clustering, m=6000 r=4.26 chain wiring, g=6, 9 colours (greedy),
  thresholds tuned to Phi(empty)=0.5: Phi monotone nonincreasing, final 0.028, no failed event. Cost per clause g*deg*g^2.
- test_r4.py: K=4,5,6 random signed tables: brute = sweep = D&C = outer-square+sweep (all d) exactly; step counts = s(K).
- count_r4.py: K=64..1024: square scheme and hybrid (middle cut + square in halves) satisfy T <= (1-M)/3 + 1/(2K). CONFIRMED.
## Own constants / thresholds (02:10)
- mylib.py: own zeta_R2 (matches b7 to 3 digits: x=49 t=2^15 6.88e-3), own R3 terms (own constants: A-range P/(g-2), gamma
  range P/g+P/(g-2), A-term 4x smaller than b7's, budget fixed dA=dG=dGam=0.33 with forward+reverse G events), row-sweep
  base with exact padding optimisation + exact max-min over the arc (minimax = optimal-w bound), Lemma U level via closed form.
- Note on padding in an unbalanced (g != h) row-sweep base: with fixed padding nu = rho the weighted-w argument can fail
  in a corner (U_row < nu N, U_col > nu N) by (1-w)(zs(g)-zs(h)); harmless because the unpadded bound x/4 + cA/x is much
  smaller there; my implementation uses min over padding levels (rigorous).
- mythr.py: own top-level (Lemma 8' derandomised) threshold code; calibration on p2's E2 chain: 4740 / 9530 (r=4.26) exact.
## DP results (own code, ratio 1.02, lt step 1, x < 2^15; c1 = max over log m/a in [Y,2Y])
- clique-sweep base + R2 + Lemma U: a=2 [8000,16000] 4.663 (b7 4.665)  CONFIRMED
- poly (D&C base) + R2: 6.162 (b7 6.160)  CONFIRMED
- row-sweep (R3, own constants, exact max-min) + R2: 3.996/4.022/4.038/4.050 for Y=1000..8000 (b7 with its constants 4.12)
  -> b7's 4.12 is a valid (slightly conservative) value; own tighter constants give 4.05.
## Thresholds (mythr.py on DP solver tables)
- r=4.26 vs 2^n: rowR3 2237 (NOT ~1.4e3 as b7's heuristic: the (c1)^2 scaling used p2's sawtooth limsup 7.47, but p2's 4740
  sat at a favourable sawtooth point; effective constant sqrt(4740)/11.78 = 5.84).
- finer DP (ratio 1.015): rowR3 4.031 at [8000,16000]; mu(16000)=0.5033 -> limit ~4.05-4.06.
- Threshold table (deterministic = LV by R1; top grouping p2-Hoeffding [R2 at top in brackets]):
  rowR3+R2+LemmaU: r=4.26: 2237 [2225] / 4385 [4372] / 1.494e4 [1.493e4]; r=1: 73.5 [67.5] / 117.2 [102.1] / 415 [400]
  clique+R2+LemmaU: r=4.26: 2988 [2976] / 5756 [5738] / 2.018e4 [2.017e4]; r=1: 99.0 [93.0] / 164.7 [147.6] / 560 [547.5]
  poly (D&C)+R2:    r=4.26: 5123 [5103] / 1.011e4 [1.01e4] / (running, hi=6e4); r=1: 131.3 [120.3] / 272.4 [259.5] / 981 [962]
  (p2 E2 LV: 4740/9530/2.68e4; 108/202/696.)
- calibration of my DP: rowb7c (b7's constants) 4.117 = b7's 4.117; rowp9old (p9 Lemma RS + Hoeffding) 4.817 = p9's 4.817.
- more thresholds (Hoeffding top): b7 constants rowR3: 2323/4453/1.553e4 | 77.6/121.2/430;
  p9 machinery + R1 (det): 3118/6158/2.133e4 | 81.8/147.6/557; poly PPSZ r=4.26: 3.495e4.
## Verdicts
R1 CONFIRMED (needs finite-precision remark). R2 CONFIRMED. R3 CONFIRMED up to a probability-budget slip
(dA+2dG+dGam/2 = 7/6 with b7's constants; fixed budget changes only logs; own constants give smaller zs).
R4 CONFIRMED (T <= (1-M)/3 + 1/(2K), exact small instances). Constants 4.12 / 4.67 / 6.16 reproduced (own: 4.05 with
tighter R3 constants). b7's heuristic m* ~1.4e3 REFUTED: actual 2.24e3 (2.32e3 with b7's constants) at r=4.26 vs 2^n.
No process left running.
