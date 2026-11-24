# G21d — an algebraic witness separating the swap from J(2,2) at N = 3 (working notes, saved as I go)
Tags: [P] proof written here (unreviewed), [E] exact computation (mod prime unless said), [N] numerical, [O] open.
Conventions: X(x,y,u,v), legs 1=x 2=y (outputs) 3=u 4=v (inputs).  Fl1 = (12|34) rank, Fl2 = (13|24) osr (realignment),
Fl3 = (14|23).  F = swap = delta_{xv} delta_{yu}.  J(p,q) = closure(sigma_p(Fl1) + sigma_p(Fl2)), N = 3 unless said.

## 1. Symmetric reformulation  [P]
Group: (a,b,c,d) in GL_N^4 acts by X -> (a(x)b) X (c(x)d) (all four legs independently); rank (Fl1) and osr (Fl2) are invariant.
Leg permutations: S_4 permutes the three pairings 12|34, 13|24, 14|23 through S_4 -> S_3 (kernel Klein V_4).  Those fixing
Fl1 and Fl2 each = V_4 (transpose X -> X^T, conjugation X -> F X F, and their product); those exchanging Fl1 <-> Fl2 (p <-> q) give
D_4 (order 8).  Full symmetry group of the pair (sigma_p, tau_q):  GL_N^4 x| V_4 (and x| D_4 for the unordered problem).
Stab(F) = {(a, b, b^-1, a^-1)} (x| V_4): dim 2N^2.  Orbit of F = {M(x,v) M'(y,u): M, M' invertible}, dense in
Z3 := S_(14|23) = rank-one tensors in the THIRD pairing.  Since J(p,q) is closed and G-stable:
   F in J(p,q)  <=>  Z3 subset J(p,q)  <=>  sigma_1(14|23) subset closure( sigma_p(12|34) + sigma_q(13|24) ).
Symmetric statement: for the three pairings alpha, beta, gamma of 4 legs, "rank-one along gamma, full support" vs
"rank <= p along alpha + rank <= q along beta".  General 3-parameter question: sigma_r(gamma) subset J_{alpha,beta}(p,q)?
Witness <=> f in I(J) with f|Z3 != 0.  Representation theory (Cauchy): C[Z3]_(d,d) = (+)_{mu,nu |- d} S_mu V_x (x) S_nu V_y (x)
S_nu V_u (x) S_mu V_v (multiplicity one each), so a separating equation exists in degree d iff some highest-weight vector of
weight (mu, nu, nu, mu) (x,y,u,v) in I(J)_d is nonzero on Z3.  Only weights with mu, nu having 3 nonzero rows can matter at N=3
(if mu has a zero row, the polynomial lives on the coordinate projection to the 2x3-format problem; see section 3).

## 2. dim J(2,2), N = 3  [E, Terracini mod 32003, terracini.py]
dim J(2,2) = 64 = expected (2*2*9-4)*2, codim 17 in the 81-dim space (NOT 28: the 2N^2 p - p^2 count is with N^2 = 9 rows).
Others: J(3,1) = J(1,3) = 62, J(2,1) = 49, J(3,3) = 81.  N = 4: J(3,3) 174 (codim 82), J(4,3) = J(3,4) 199 (codim 57), J(4,2) 172.

## 3. Non-membership of F in J(2,2), N = 3, over C  [N]
dist.py (variable projection, tail_2 of Fl1(F - B), complex L-BFGS): 20/20 restarts -> 2.0000.  BUT calibration FAILS for this
  solver: on a random point of J(2,2) (norm^2 9) it stalls at 0.07-0.4 (12 restarts) -> vp is not trustworthy here (gradient checked OK;
  landscape problem).
distlm.py (Levenberg-Marquardt on all factors a,b,G,H, complex): calibration random J(2,2) point: 6/8 restarts -> 1e-30;
  swap at (3,1): 5/5 -> 1e-30; swap at (2,1): 3.0 (Thm 3).  Swap at (2,2): 23/23 restarts -> dist^2 = 2.0000 (||F||^2 = 9).
  Two independent solvers agree: inf dist^2(F, J(2,2)) = 2 [N], relative distance sqrt(2/9) = 0.47.
Structure of the near-optimum (inspect_opt.py): factors DIVERGE (|A| ~ |B| ~ 500, A + B bounded): the nearest points of J(2,2) to F
  are BORDER points L with Fl1(L) of rank 7, singular values (1,..,1,0,0), and F - L of rank 2 in all three flattenings with singular
  values (1,1).  I.e. "F with a rank-2 piece removed" is (border) in J(2,2); F itself stays at distance sqrt 2.
Projection check (distfmt.py): the coordinate projection to format (2,3,3,3) (x restricted to 2 values) has codim 2 (proj.py) but
  CONTAINS the projected swap (LM dist^2 1e-31): coordinate projections lose the separation.  Reynolds projections onto
  K-invariants for 12 finite subgroups K of Stab(F) (reyn.py) all FILL W^K: no reduction there either.

## 4. Equations of J(2,2), N = 3
### 4a. Flattenings and Koszul/Young flattenings  [E: koszul.py, out_koszul_3_22.txt]
168 linear maps: Koszul flattenings Lambda^k A (x) S* -> Lambda^{k+1} A (x) S^c with A = one leg (k = 0,1,2) or A = a pair of legs
(dim 9, k = 0,1,2), S any subset of the remaining legs (k = 0 = ordinary flattenings, all bipartitions incl. 1|3).
Result: on J(2,2) EVERY one of them has generic rank = min(#rows, #cols) (maximal; mod-p rank at a random integer point of J,
a valid lower bound for the char-0 generic rank).  So none of them gives ANY equation of J(2,2), separating or not.
At F the ranks are often SMALLER (F is rank one in Fl3, so e.g. Fl3 and the Koszul maps built on it have rank 1 or 3 at F).
This is a concrete instance of G20d's barrier (pq = 4 >= N^2/4) and shows more: Koszul flattenings (not covered by Theorem B)
do not help here either -- they vanish identically on nothing.
### 4b. Interpolation by highest-weight vectors  [E mod 32003, probabilistic: hw.py, hw2.py; interp.py for full weight spaces]
Method: I(J) is GL_3^4-stable; a separating equation exists in degree d iff some h.w. vector of weight (mu,nu,nu,mu) in I(J)_d is
nonzero on Z3 (section 1).  h.w. vectors of weight lam in S^d(W*) are spanned by tableau polynomials (column antisymmetrisations
of X^{(x)d}, one filling per leg); their number is the generalised Kronecker coefficient g(lam) (kron.py, Murnaghan-Nakayama).
For each lam: random fillings until the rank on generic points reaches g(lam) (this independently re-checks g(lam)); evaluate on
g+25 random points of J(2,2) (rank r_J) and 20 points of Z3; eqs = g - r_J, sep = rank(J+Z3) - r_J.
Positive control: sigma_1(Fl1) (p,q) = (1,0): finds the 2x2 Fl1-minors in weight ((2),(1,1),(2),(1,1)) and the separating one in
((1,1),(2),(2),(1,1)) (sep = 1) -- the pipeline detects equations and separation.
Cross-check: for d <= 5 the full weight-space interpolation (interp.py: monomial bases up to 6860 monomials) agrees (no equations).
Rigour of "no equation": a full-column-rank evaluation matrix mod p at reductions of integer points of J, for a basis that is
independent mod p, excludes a char-0 equation in that weight (clear denominators of a primitive combination).  So the
"eqs = 0" results below are CERTIFICATES (char 0), given: (i) tableau polynomials are h.w. vectors (standard; checked: the span of
48-138 random tableau polys never exceeds g, 3 weights) and (ii) g(lam) (Murnaghan-Nakayama; independently re-checked by the rank
reached).  Positive results (eqs > 0) would be probabilistic; none occurred for J(2,2).
RESULTS (N = 3, J(2,2)):
 * ALL dominant weights, degree <= 6 (d = 3,4,5: 13+40+102 weights up to D_4; d = 6: 289 non-separating-type + 6 separating-type;
   out_hw2_all345.txt, out_hw2_all6.txt, out_hw2_6.txt): NO equation of J(2,2) at all.  [E]
   => J(2,2) (codim 17) has no polynomial equation of degree <= 6; in particular no separating one.
 * degree 7, separating weights (mu,nu,nu,mu): see out_hw2_7.txt (partial, listed in section 6).
 * control J(2,1) (codim 32, swap at distance^2 3.0 [N], infeasible exactly by Thm 3): no separating equation of degree <= 6 among
   the 15 separating weights with g <= 60 either (out_ctrl_21.txt).  Low-degree witnesses do not exist even where a spectral proof does.
### 4c. A structural fact about Z3 n J(2,2)  [N + P]
 Z3 = {M(x,v) M'(y,u)}.  If det M = 0 (or det M' = 0) then M(x)M' is in J(2,2): by GL it is the swap with one x-row deleted, i.e.
 the format-(2,3,3,3) swap, which LM puts in J(2,2) at distance 1e-31 (distfmt.py) [N]; M' singular by the D_4 symmetry.  With the
 non-membership of F this gives (numerically) Z3 n J(2,2) = {det M * det M' = 0}.  Consequently the restriction to Z3 of any
 separating equation is divisible by det M det M' -- consistent with the representation theory: a h.w. vector of weight
 (mu,nu,nu,mu) restricts to prod of leading principal minors Delta_1^{mu1-mu2} Delta_2^{mu2-mu3} Delta_3^{mu3}(M) x (same for M'),
 and mu3, nu3 >= 1 is forced (section 1).  The minimal conceivable degree (mu = nu = (1,1,1), the restriction det M det M') does not occur.
 Calibration of the border-following ability of LM [N]: N = 4, (3,3) (border-only point per G19b): restart 0 reaches dist^2 7.5e-5
 with coefficients 150 (still decreasing: border sequence followed), other restarts sit in local minima at 2.0.  At N = 3 (2,2)
 all 23 restarts plateau at exactly 2.0000 while coefficients grow to ~500 -> the inf is 2, not 0.

## 5. Generalisation test (N = 4)  [E]
No witness was found at N = 3, so there is nothing to transport.  What was tested at N = 4:
 * dims (Terracini): J(3,2) codim 109, J(4,2) codim 84, J(5,2) codim 61, J(4,3) codim 57, J(3,3) codim 82 (F border-in, G19b).
 * Koszul/ordinary flattenings (200 maps, koszul.py 4 4 2, out_koszul_4_42.txt): ALL of maximal rank on J(4,2): no equations.
 * h.w. interpolation, separating weights with 4-row partitions of lowest degree, d = 4 ((1111),(1111)) and d = 5 ((2111),(2111)):
   no equations of J(4,2), J(5,2), J(4,3), J(3,2) there (out_n4.txt, out_n4b.txt).  Degree 6 at N = 4 not completed (slow).
 Honest reading: the method (lowest-degree covariant) does not produce anything checkable at these sizes; the degree of any
 separating equation is >= 7 at N = 3 and unknown at N = 4.  Nothing here says whether a witness degree stays bounded in N; given
 that J(2,2) at N = 3 has no equation of degree <= 6 AT ALL, a bounded-degree family of witnesses looks unlikely.

## 6. Degree 7, separating weights (N = 3, J(2,2))  [E]
ALL 10 separating weights (mu,nu,nu,mu), mu,nu |- 7 with 3 rows (up to D_4), g = 16..301: eqs = 0 (out_hw2_7.txt).
=> no separating equation of degree <= 7.  (Non-separating weights of degree 7 not scanned.)

## 7. Verdict
 * Symmetric reformulation (section 1): F in J(p,q) <=> sigma_1(14|23) subset J_{12|34, 13|24}(p,q); witnesses = h.w. vectors of
   weight (mu,nu,nu,mu) in I(J) nonzero on Z3, restricting to prod of leading principal minors of M and M'.
 * dim J(2,2) = 64, codim 17 [E].  inf dist^2(F, J(2,2)) = 2 (||F||^2 = 9), approached only by divergent (border) sequences [N, two
   solvers, calibrated].  Z3 n J(2,2) = {det M det M' = 0} [N].
 * NO witness found.  Clean negatives [E, certificates mod 32003 lifting to char 0]:
   - J(2,2) has NO polynomial equation of degree <= 6 (all weights, all h.w. components).
   - In degree 7 no separating equation in the separating weights listed in out_hw2_7.txt.
   - Consequently NO linear map Phi whatsoever (any flattening, Koszul/Young or exotic, any size) can separate F from J(2,2) unless
     its generic rank on J(2,2) is >= 6 (an (r+1)-minor separating would be an equation of degree r+1 <= 6), and more generally no
     determinantal or other witness of degree <= 6 exists.
   - The 168 Koszul flattenings tried all have maximal rank on J(2,2) (no equations at all), confirming and extending G20d's barrier
     beyond class B at this size.
 * Structural proof (task 3c): not achieved.  Exact Thm 4 + Grassmannian limits reduce the border case to "F|_K0 in a limit of
   4-dim spaces (U(x)W)|_K", where flat limits can be larger; I did not close it.
 * What it suggests: separating F from J(p,q) is a genuinely high-degree phenomenon already at N = 3 (degree >= 7, beyond any
   flattening of rank <= 5 on J); a general proof is more likely to come from a semi-algebraic/analytic quantity (the inf distance
   2 = an exact-looking constant, the divergent nearest-point structure F - L of rank 2 with unit singular values) or from a
   border version of the spectral/Segre arguments than from a low-degree covariant.
 Self-review: (i) "no equation" claims rely on the standard fact that tableau polynomials are h.w. vectors and on g(lam) -- both
 cross-checked numerically, not proved here; (ii) the probability that a mod-p rank is SMALLER than the true rank only makes "no
 equation" claims safer (full rank is a certificate); (iii) degree 7: complete for separating weights only; (iv) all
 distance statements are [N]; (v) the claim Z3 n J = {det M det M' = 0} uses LM membership of the 2x3 rectangular swap [N].
Files: notes.md, common.py, terracini.py, dist.py, distlm.py, distfmt.py, inspect_opt.py, koszul.py, monos.py, interp.py,
 reyn.py, proj.py, kron.py, hw.py, hw2.py, hwN.py, out_*.txt, spec*.json.
