# G26-P7 adversarial review notes
02:42 start. Read p6 notes, paper D.4.3-D.4.5 (lines 3340-3640), P3 notes.
Hand re-derivation of Thm 1 (steps 1-6): all algebra checks (see below).
02:45 smith7.py (own Smith-over-Z/2^N elimination, exact Fraction entries of (D.7)): Thm 1 formula matches in
  1412 cases a=1..5 (ALL t in S with v2(t)<a incl. negative and odd multiples of 2^i, zeta in
  {1,-1,2,-3,1/3,-7/5,8,-12,20/7,96}, b=0,1) + run67: a=6 (76), a=7 (80) cases, 0 fails.
  negctl.py: wrong alpha=2^a -> 15/30 mismatches (test discriminates).
  minors.py: exact Cauchy-minor valuation formula k v2(c)+V(R)+V(C)-sum min(i,v2(x-y)) verified on 3816 random minors (a=3..5); v2(c)=a+1+v2(zeta) asserted in all.
02:46 thm2.py: P(K)/((K/2)log2K) <= 1 for K<=2^20; identity sum s_k-2K^2 = K(m-1)-2P(K) exact, min>=0 for m<=18;
  G(theta)*4theta/zeta(2) max 0.9985 (m<=18, theta in 2^[-15,10]).
02:46 flat7.py: C1 counterexample reproduced: a=3,i=0,t=1 and t=-3, zeta=1: R_9(F)=6 vs C1 8 (equal at all even p<=16;
  also p=13: 17 vs 23, p=15: 29 vs 33). a=3,i=1: R_10=35 vs 55; a=3,i=2: R_6=37 vs 45. Square-free counts at i=a-1 confirmed.
02:50 run89: a=8 (128x128 blocks): 56 cases incl. t=+-1,+-2,-(L/2),3L/4,-3L/4,5L/8, zeta in {1,-3,1/3,2,-12,5/7}: 0 fails.
02:50 c1prime.py with NB=30 showed a spurious 'violation' at p=39 -> artifact (invariants >=30 dropped). Rerun c1prime2.py
  with modulus 2^90: a=3,i=0,t=1: sqfree <= R_p(F) <= monomials for all p<88; Smith(F) max 44 = 22+22, consistent with
  P6's Smith(M_{3,0}) and Prop 3. a=2 (i=0,1), a=3 (i=1,2) sandwich also holds (NB=30 suffices there? max sqfree sums <30: yes
  for (3,1): 2*(3+3+5+5)=32 -> p<30 only checked; minor).
Thm 2 use: paper's count = exp(theta p_t + sum_modes f(theta w)); per index k at most (3M/L)*2*2^i = 6M/D Smith modes
  (same as (D.13)); real-port submatrix invariants termwise >= (B = XAY => s_k(A) | s_k(B)); zeta in {0} u units
  (binary lifts of residue elements; zeta=0 => no crossing). So (D.15) holds with 6M/D -> (6M/D)/4 exactly. CONFIRMED.
Caveats: P(K)<=(K/2)log2K cited, not proved in p6 (classical; checked K<=2^20). c=37.5 is P3's closed form with nu->nu/4;
  finite-lambda thresholds (2.25e5/1.5e6) are scaled, not recomputed (leaf term does not halve).
02:50 run9.py: a=9 (256x256): t=1 (i=0), t=-6 (i=1), t=-384 (i=7, zeta=1/3): all match. Total own checks a=1..9: 1627 blocks, 0 fails.
VERDICTS: Thm1 CONFIRMED; Thm2 CONFIRMED (popcount lemma cited, classical); C1 refutation CONFIRMED; C1' numerics confirmed a=2,3.
Row H: rigorous relative to App. D + P3 rows A-G + P3's [P] claim that the Smith-mode algorithm is correct; c=37.5 closed form only.
