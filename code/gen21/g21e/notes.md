# G21e — adversarial review of G20d (barrier Theorem B) and G20c (the 2/3 bound)

Tags: [P] re-derived here, [E] exact computation (my scripts in this folder), [N] numerical.

## CLAIM I (G20d Theorem B): verdict CONFIRMED (with two scope corrections)

(i) Key lemma [P, E].  Statement: for complex A (a' x N^2, rank a), B (N^2 x b', rank b):
    max_{G,G'} rank A(G(x)G')B >= ab/N^2.
 Re-derivation.  rank(AYB) = rank(P1 Y P2) (P1 = orth. proj. on (ker A)^perp, P2 on range B).  For Haar (g,h) in U(N)xU(N),
 T := E[(g(x)h)P2(g(x)h)^*] commutes with U(N)xU(N); C^N(x)C^N is (absolutely) irreducible, so T = (tr P2/N^2) I = (b/N^2) I.
 E||P1 (g(x)h) P2||_F^2 = tr(P1 T) = ab/N^2.  M = P1 U P2 has ||M||_op <= 1, so rank M >= sum sigma_i^2 = ||M||_F^2 POINTWISE;
 hence E[rank] >= ab/N^2 and some unitary pair attains >= ab/N^2.  The averaged quantity is ||.||_F^2, the pointwise inequality
 rank >= ||.||_F^2 is what makes it rigorous (no linearity of rank used).  Correct.  Char-0 fields: A, B have entries in a
 finitely generated subfield, which embeds in C; the generic rank of a polynomial matrix does not change under field extension.
 Tests (lemma_test.py, lemma_gf.py): exact ranks over F_10007, F_3, and with G,G' from GF(4096), GF(2187) (generic over the
 closure in char 2, 3): random prescribed ranks, all coordinate-subset pairs (N=2, exhaustive), Kronecker coordinate subspaces,
 Sym/Alt spans, diagonal span; N = 2,3,4: min s N^2/(ab) = 1.000 everywhere (equality at Kronecker compressions = Thm C maps).
 Apparent failures with G,G' drawn from the PRIME field F_2 (0.889, 0.844) are small-field sampling artefacts: f2_exh.py
 (exhaustive over all F_2-points at N=2, a=b=3) finds no failure, and GF(4096) gives 1.000.  So in char 2, 3 the lemma held in
 every test (with c's = generic ranks over the closure); it is unproved there, as G20d says.

(ii) Per split [P].  Fl1: c1 <= 1 (A uv^T B rank <= 1), c1 = 1 if block nonzero (rank-one operators span), c2 >= ab/N^2 by the
 lemma, R <= min(a,b): R^2 <= ab <= N^2 c1 c2.  Fl2: symmetric.  Fl3 (xv|yu): Fl3(UV^T) = U(x)V^T, Fl3(G(x)G') = (G(x)G'^T)Pi;
 c1, c2 >= max(1, ab/N^2), so N^2 c1 c2 >= ab in both cases ab >= N^2, ab < N^2.  1|3 splits: R <= N, c1, c2 >= 1.  0|4: R <= 1.
 Block sums: on the irreducible varieties {rank <= 1}, {osr <= 1} a generic point attains every block's max simultaneously, so
 c_i = sum_t c_i,t; R = sum R_t; Cauchy-Schwarz sum sqrt(c1t c2t) <= sqrt(sum c1t)(sum c2t).  Correct.
 Attack (thmB_attack.py, exact mod 10007): all 16 bipartitions, coordinate / random-rank / Kronecker A,B, X = F and random X,
 N = 2, 3: max R^2/(N^2 c1 c2) = 1.000 (attained: e.g. 1|3 split at N=2 (1,1,2); Fl2 at N=3 (9,1,9)).  No violation.
 SCOPE ERROR 1: G20d says "restrictions, compressions, partial transposes, realignments and their compositions are all of this
 form".  False for compositions with NON-LOCAL maps: X -> P Fl_S(C X D) Q with C, D general N^2 x N^2 is not A Fl_S(X) B.  Only
 local compressions (C = C1(x)C2, D = D1(x)D2) commute into the form.  Theorem B is not proved for those maps (they belong to the
 "general Phi" open case).  Random test: ratio <= 1 for them too (max 1.0 at N=2, 0.889 at N=3), so no counterexample, but the
 class must be stated without "compositions".

(iii) Corollary [P].  Excluding (p,q) by Phi means c1 p + c2 q < R.  AM-GM: c1 p + c2 q >= 2 sqrt(c1 c2 p q).  With
 R <= N sqrt(c1 c2): 2 sqrt(pq) < N, i.e. pq < N^2/4.  A family excludes a union of such half-planes: still inside pq < N^2/4.
 Correct.  Wording: "B-methods give EXACTLY pq < (1/4+o(1))N^2" — the barrier side is proved; the attainment side holds only for
 q = o(N) with N/(2q) near an integer (k = round(N/2q)); at e.g. q = N/3 the best known B-map gives 2N^2/9 (k=1,2), and whether
 another B-map reaches N^2/4 there is open (sums of Thm C maps only add inequalities; off-diagonal compressions pi_k, sigma_j give
 R = min(k,j)N, no gain).  In exponents (theta + T = 2) "exactly at the line" is right.

(iv) Theorem C as a fixed map [P].  Phi_k(X) = (pi_k(x)1) X (1(x)sigma_k): rank-one -> rank <= 1; G(x)G' -> (pi_k G)(x)(G' sigma_k),
 rank <= k^2 (and = k^2 generically); F -> delta_{xv}delta_{yu} restricted to x,v <= k: rank kN.  Profile (1,k^2,kN) gives
 p + k^2 q >= kN = Thm C(a) with Q=1, n=m=N.  (b): same compression of Fl2(X): (k^2,1,kN).  The adaptive original proof and the
 fixed map give the same inequality: CONFIRMED (also test_thmC.py).  Lean SwapRank.lean has no sorry (statements not re-read).
 Theorem C(c) (b < mn => r >= Q): for the crossing (Q=1) it says (0,q) infeasible for q < N^2, which IS a B-method (Fl2 itself,
 profile (N^2,1,N^2)); for Q > 1 it is false in the border sense (G19e), hence not any fixed linear rank method.  G19b Thm 3
 (excludes (p,1), p < N(N-1)/2: pq > N^2/4 for N >= 4; NOT at N = 3 where (2,1) has pq = 2 < 2.25) and Thm 4 ((2,2) at N=3,
 4 > 2.25) exclude pairs no B-method can, so they are outside B — conditional on those (unreviewed, char != 2) theorems.

(v) "char 0": proved; char p: open (tests support it with generic c's).  "exact or border": correct (any fixed linear rank
 method is border-valid and the barrier never uses exactness).

## CLAIM II (G20c (F1)-(F2)): verdict PARTLY

(i) Curve [P].  Paper (54)/(55): R_s^row = sum_{i=0}^{2h} p(s+i) <= (2h+1)p(2h+s), R_s^col <= (2h+1)p(2h-s); Def. 7 bounds
 c_row + c_col of EACH piece, so 4 rho = max_s [L(2h+s)+L(2h-s)]/lambda up to O(log h).  For a general graded space the honest
 quantity is the window dim V_{[s,2h+s]} <= D(2h+s) = dim V_{<=2h+s}.  G20b's threshold form (P(2h+s0-1), P(2h-s0)) is the same
 curve shifted by one weight unit: irrelevant when L varies on scale h, RELEVANT for step profiles (below).  lambda = log N with
 N = number of strings = size of the dual pair of bases of Lemma 10 inside V_h, so lambda <= L(h) (not a normalisation: an
 inequality, equality when all of V_h is used).

(ii) Optimisation (profile_lp.py: LP over f = L/lambda on 0..4h, f(0)=0, nondecreasing, f(h) >= 1, the line at EVERY threshold
 f(2h+s0-1)+f(2h-s0) >= 2, cost max_s f(2h+s)+f(2h-s)):
  * concave: optimum 4rho = 2.714, 2.678, 2.6695 (h = 5, 20, 80) -> 8/3: rho >= 2/3 CONFIRMED (only concavity on [h,4h] is needed
    for the lower bound: max_s >= value at s=0 = 2L(2h)); kinked profile N, N^{4/3}, N^2 is concave and line-consistent; its curve
    10/3 - 2|s|/(3h) re-derived.
  * WITHOUT concavity: optimum 4rho = 2 + 1/h -> 2, i.e. rho -> 1/2, attained by ramp-plateau-ramp: f = E/h on [0,h], 1 on [h,3h],
    1 + (E-3h)/h on [3h,4h] (LP output at h=20: 0, .5, 1, 1, 1.05, 1.55, 2).  It lies ON the line at every threshold, so it is
    compatible with pq >= (2/9)N^2 (up to constants/poly(h) factors, which the exponent framework does not see).
  * A pure STEP function does not work: f = 1 on [1,4h-1], 2 at the end costs 3 (rho = 3/4) once the line holds at every
    threshold; f(4h) = 2 alone costs 2 but violates the line at s0 = 2h (f(0)+f(4h-1) = 1).  So the task's guess "step function"
    is wrong; the right non-concave optimiser is the plateau.
(iii) Concavity is an ASSUMPTION, not a property of the objects.  log D(E) is concave for polynomial rings / Fock spaces /
 finitely strongly generated VAs asymptotically (Cardy growth c sqrt E), but: generalized free fields of weight w have D(E)=1 on
 [0,w) then jump (non-concave at small E); graded MODULES with a large lowest-weight space at weight E0 (e.g. V = Fock (+) M with
 dim M_{E0} = N^c at E0 ~ 3h) have L = max of two concave functions, which is not concave — exactly the plateau-then-rise shape.
 Whether Lemma 9 (locality) and Lemma 10 (dual bases at weight h) extend to such spaces is NOT known; that, not concavity, is the
 real constraint.  (G20c itself mentions the plateau giving the line, but its headline "cannot do better than 2/3" omits the
 hypothesis' status.)
(iv) Paper's remark (App. H, l. ~3009) [P]: L = c E^alpha, lambda = L(h): 4rho = 2L(2h)/lambda = 2^{1+alpha}, rho = 2^{alpha-1};
 alpha = 1/2 gives 1/sqrt2.  Floor: L(4h) >= 2 lambda <=> 4^alpha >= 2 <=> alpha >= 1/2.  Using N < dim V_h (lambda < L(h)) to
 satisfy the floor with alpha < 1/2 gives rho = 2^{-alpha} > 1/sqrt2.  So "vacuous for alpha < 1/2" CONFIRMED.

## Corrected statements (quotable)
I.  Let Phi = (+)_t A_t Fl_{S_t}(X) B_t (block-diagonal; S_t any bipartition of the four legs; A_t, B_t arbitrary matrices over a
 field of characteristic 0; amplifications included).  Then rank Phi(X)^2 <= N^2 c1(Phi) c2(Phi) for every X, c_i being generic
 ranks on rank-one / Schmidt-rank-one operators.  Consequently no such Phi, nor any union of the half-planes they exclude,
 excludes a pair (p,q) with pq >= N^2/4, exactly or in the border sense, for any operator.  The class does NOT include
 compositions with non-local maps (P Fl_S(CXD) Q with C, D not Kronecker); for those and for general linear Phi, and in positive
 characteristic, the statement is open.  Theorem C(a),(b) are B-maps with profiles (1,k^2,kN), (k^2,1,kN).
II. In the graded framework of App. B, with L(E) = log dim V_{<=E} and lambda = log N <= L(h), the pieces have exponent pairs
 (L(2h+s), L(2h-s))/lambda + O(log h/lambda); the line forces L(a) + L(4h-1-a) >= 2 lambda - O(1) for all a (in particular
 dim V_{<4h} >= N^2/2).  For L(E) = cE^alpha this forces alpha >= 1/2 and gives rho = 2^{alpha-1} >= 1/sqrt2.  IF L is concave on
 [h,4h], then rho >= 2/3, attained (as a profile) by dims N, N^{4/3}, N^2 at h, 2h, 4h.  Without concavity the dimension
 constraints alone allow rho -> 1/2 (profile flat on [h,3h]); concavity is a hypothesis, not a theorem about the framework.

Files: ff.py, lemma_test.py, lemma_gf.py, f2_exh.py, thmB_attack.py, profile_lp.py, notes.md (this folder).
