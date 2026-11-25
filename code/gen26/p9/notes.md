# G26-P9 review notes (adversarial review of G26-P3)
03:00 lemma_check.py: Lemma 1 compression err ~1e-15, #eig>=1 >= m, ||Ts||<=Q, |tr Ts^j|<=Q^j N (j=1..9) on 5 optima
  (exact p recomputed by sympy: 5,7,6,7,12 ok), 400 random signed MC, 300 random real BN (Q=sum ||G|| ||G'||): 0 violations.
03:00 mycert.py (own LP, own exact Bernstein-subdivision verification + sympy Sturm cross-check):
  q=3 a=.9 d41 L=17.85 ; q=6 a=.5 d41 L=10.12 ; q=2 a=.99 d61 L=16.51 ; q=10 a=.5 d61 L=92.57  -- all verified,
  slightly BETTER than P3's table (18.01, 10.16, 17.68, 93.4). Table certificates CONFIRMED.
03:10 extra verified certs (cert_*.json): q=3 a=.3 d21 L=.90 (p>=.2308N^2-.693N); q=3 a=.5 d21 L=2.12 (N^2/3-1.41N);
  q=3 a=.7 L=4.32; q=4 a=.5 L=3.43; q=7 a=.5 L=17.6.  NB: P3's table certs are VACUOUS (negative) at N<=6, so the
  5 exact optima only test the inequality chain, not the bound. opt_test.py: chain a m-(n-m) <= tr g(Ts/Q) <= L N holds on
  all 5 optima for all applicable certs; bounds <= p. kron_test.py: (5,3)@4, (12,3)@6 (x) I_m / (x) signed perm, N<=24:
  p >= best q=3 bound everywhere (e.g. N=24: p=288 vs bound 158).
thmB.py: Fejer kernel r=1: odd, monotone, |g|<=1, a*Q/(K+1) in [0.51,0.71] >= 16/pi^4=0.164 (needs K+1<=pi Q/2), L far below
  claimed 2pi(1+ln(2K+2))(1+sqrt2)^(2K+1). Theorem B derivation OK (constant: (1+a) denominator ignored, harmless).
03:10 grp_check.py (Theorem G): Z3 circulant {-1,0,1} q<=2 exhaustive min p=3 (=bound 3); Z4 q=1 exhaustive 6 (=6);
  S3 permutation rep triv+std (NON-regular), K over the full 4-dim span (any q), 13^4 grid: min p=3 = bound 3 (tight);
  Z2 on R^4 m=(2,2) random: >=4; Z4/Z5 random group-algebra q<=4: min = bound. Proof re-derived by hand: correct
  (Theta beta(x,y)=beta(y,x)Theta; off-diagonal blocks rank>=dim; diagonal blocks Theta=F(x)Theta' via Schur, rank(1-T)+rank(1+T)>=dim).
Proof review Lemma 1/2/3: re-derived; correct over R (C with adjoints). Lemma 3 needs only g>=-1 on [-1,1] and g>=a on [1/Q,1]
  (g<=1 not needed). Caveats: (i) MC must be SIGNED (+-1) partial perms; general weighted monomial ("generalized permutation")
  matrices are not covered with Q=q, only with Q=sum ||G_l|| ||G'_l|| (gauge-optimised). (ii) char 0 only (rank over R/C).
  (iii) table bounds positive only for N > L/a (q=3,a=.9: N>=21).
SA (sa_mc.py) q=3: N=5 best p=10, N=6 best 15 (SA weak: known 12); far above bounds (3, 5).
03:12 S3 regular random group-algebra q<=4: min 15 >= bound 13. SA q=3 N=7: 21, N=8: 28 (bounds 7, 11). All processes killed.
Item 3 (BN): proof uses only sum ||G_l|| ||G'_l|| -> correct. No contradiction: f(N,1)=f(N,2)=N(N-1)/2 (Triang.lean, any norm)
  is stronger than BN bound at Q<=2; App. B has q=N^beta and super-polynomial norms (P7: L_inv); Pauli Q=N p=0 makes
  Thm A vacuous (L~e^{cQ}); j=1 alone gives (N^2-QN)/(1+Q) (MC1 line generalised to BN); Q<1 => p=N^2 trivially.
Item 5 novelty: ingredients standard (Courant-Fischer compression, trace method, odd sign-approximating polynomials as in
  QSVT/approximate degree). BN log-gain resembles bounded-coefficient circuit lower bounds (Morgenstern; Raz n^2 log n).
  Per g22a reading, no (p,q) trade-off beyond k=1 appears in Naslund/Karam/HWY/Hrubes/RZ; HWY's Z-model 6/5 bound is a
  different (sum-of-squares) model. Not checked against wider literature.
