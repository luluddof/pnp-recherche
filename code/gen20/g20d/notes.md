# G20d — rank methods for the crossing trade-off: definitions, barrier, generic frontier (working notes, saved as I go)
Tags: [P] proof written here (unreviewed), [E] exact computation, [N] numerical, [L] literature, [O] open.

## 0. Setting
V = K^N (x) K^N, F = swap.  X ranges over End(V) = 4-tensors X(x,y,u,v), operator (u,v)->(x,y).
Three flattenings: Fl1 = (xy|uv) [rank], Fl2 = (xu|yv) [osr = realignment], Fl3 = (xv|yu) [partial transpose].
F: rank_1 = rank_2 = N^2, rank_3 = 1  (F is a rank-one tensor in the third grouping: F = Omega(x,v) Omega(y,u)).
(p,q) feasible: F = A + B, rank_1 A <= p, rank_2 B <= q.
Literature read: EGOW arXiv:1710.09502 (Thm 4.4 / Lemma 3.5 read in the pdf, dl1/ untrusted).

## 1. Definition of a rank method; Theorem C is one  [P, E]
RANK METHOD (linear, fixed): a linear map Phi : End(V) -> Mat_{a x b}(K) with
   c1(Phi) := max{rank Phi(Y) : rank_1 Y <= 1},  c2(Phi) := max{rank Phi(Y) : rank_2 Y <= 1}  (both = generic values on the
   two Segre-type varieties, by semicontinuity), R(Phi) := rank Phi(F).
   Conclusion: every feasible (p,q) has c1 p + c2 q >= R  (subadditivity).  Valid for border decompositions too (rank is lsc).
   A PROOF BY RANK METHODS = a family of such Phi; it excludes the union of the half-planes {c1 p + c2 q < R}.
Theorem C (a) = Phi_k(X) = (pi_k (x) 1) X (1 (x) sigma_k)  (compress x to k, restrict v to k):  (c1,c2,R) = (1, k^2, kN),
   i.e. p + k^2 q >= kN.  (b) = same map after realignment Fl2: (k^2, 1, kN).  The published proof is adaptive (S depends on B)
   but the inequality is exactly that of this FIXED map.  test_thmC.py: exact profiles for N = 3,4,5, all k: R^2/(N^2 c1 c2) = 1.
G18g Theorem 2 = Theorem C.  G19b Theorem 5 = compression (rank method) + Theorem 3 on the one-term remainder (not rank).
G19b Theorem 3 (q = 1, charpoly / eigenvalue multiplicity) is NOT a rank method: it excludes (p,1) for p < N(N-1)/2, beyond
   what any method of class B below can exclude (pq < N^2/4).  G19b Theorem 4 (Segre variety cut by a linear space, dimension
   count) is NOT a rank method either (nonlinear: uses dim of {ab^T in ker R}); it stays inside pq < N^2/4 as well (p >= N at q = 2).

## 2. Barrier theorem  [P over C / any char-0 field; tested E/N]
CLASS B ("flattening-compression methods"): block-diagonal direct sums Phi = (+)_t Phi_t (amplification X -> Phi(X) (x) I_m
included), each block Phi_t(X) = A_t Fl_{S_t}(X) B_t with S_t ANY bipartition of the legs {x,y,u,v} (3 balanced ones
Fl1, Fl2, Fl3, the unbalanced 1|3 ones, and the trivial one) and A_t, B_t ARBITRARY matrices (restrictions, compressions,
partial transposes, realignments and their compositions are all of this form).
THEOREM B.  For every Phi in B and EVERY operator X on C^N (x) C^N:   rank Phi(X)^2 <= N^2 c1(Phi) c2(Phi).
COROLLARY.  If a family of B-methods excludes (p,q) for ANY operator X (in particular for F, exact or border), then
   c1 p + c2 q < R <= N sqrt(c1 c2), and AM-GM gives 2 sqrt(c1 c2 p q) < N sqrt(c1 c2):  p q < N^2 / 4.
   So no proof by B-methods excludes a single pair with p q >= N^2/4; in exponent form, B-methods cannot move T(theta) above
   the line 2 - theta at any theta, for F or for any other operator.  Tight: Theorem C(a) with k = round(N/2q) excludes
   p < k(N - kq) ~ N^2/(4q).  So B-methods give EXACTLY the region p q < (1/4 + o(1)) N^2 (up to integrality).
KEY LEMMA (Schur averaging).  For any complex matrices A (a' x N^2), B (N^2 x b'):
      max_{G,G'} rank A (G (x) G') B  >=  rank(A) rank(B) / N^2.
 Proof.  rank(AYB) = rank(P1 Y P2), P1 = orth. projection on (ker A)^perp (dim a = rank A), P2 on range B (dim b).  For Haar
 (g,h) in U(N) x U(N), E[(g(x)h) P2 (g(x)h)^*] = (b/N^2) I because C^N (x) C^N is irreducible for U(N) x U(N) (Schur's lemma).
 So E||P1 (g(x)h) P2||_F^2 = ab/N^2, and rank M >= ||M||_F^2 / ||M||_op^2 >= ||M||_F^2 (||M||_op <= 1).  QED
Proof of Theorem B.  Cauchy-Schwarz: if R_t <= N sqrt(c1_t c2_t) then sum R_t <= N sqrt(sum c1_t)(sum c2_t), and for a
 block-diagonal map generic ranks add (c_i = sum c_i,t), R adds.  One block, a = rank A, b = rank B, R <= min(a,b) <= sqrt(ab):
  * Fl1: Fl1(G(x,u)G'(y,v)) = G (x) G' (rows (x,y), cols (u,v)); c1 <= 1 (rank-one operators stay rank <= 1; c1 >= 1 if the block
    is nonzero, since rank-one operators span), c2 >= ab/N^2 by the Lemma: R^2 <= ab <= N^2 c2 = N^2 c1 c2.
  * Fl2: symmetric (Fl2(u v^T) = U (x) V, Kronecker of the reshaped vectors): c2 = 1, c1 >= ab/N^2.
  * Fl3: Fl3(u v^T) = U (x) V^T and Fl3(G(x,u)G'(y,v)) = (G (x) G'^T) up to a fixed column permutation: c1, c2 >= ab/N^2, and
    c1, c2 >= 1; if ab >= N^2 then N^2 c1 c2 >= (ab)^2/N^2 >= ab; if ab < N^2 then N^2 c1 c2 >= N^2 > ab.  (For F itself R <= 1.)
  * unbalanced (one leg | three legs): R <= N, c1, c2 >= 1.  Trivial bipartition: R <= 1.                                QED
 Notes. (i) The proof never uses X = F: B-methods give the SAME bound for F and for a generic operator, whose true frontier
 (section 3) is p + q ~ 0.59 N^2 -- an exact EGOW-type barrier (EGOW: rank methods cannot exceed ~ n^{floor(d/2)} for tensor rank).
 (ii) Over a char-0 field K: c's are generic ranks of polynomial matrices, same over K and C.  Char p: the Lemma is not proved
 (Haar averaging), Theorem B open there.  (iii) Valid verbatim for border decompositions (the method does not see the difference).
 (iv) Theorem C's original adaptive proof yields the same inequality as a fixed B-map (section 1), so it is covered.
Tests: test_thmC.py (Theorem C maps: ratio R^2/(N^2 c1 c2) = 1 exactly); schur_test.py (Lemma: min s N^2/(ab) = 1.000 over
 random sparse and Kronecker-structured A,B, N = 2,3,4; equality = Theorem C's maps); lemmaS.py (single Fl1 block, F: max ratio
 R^2/(N^2 s) = 1.000, N = 2,3); search.py (hill-climb over NON-block-diagonal sums of Fl1/Fl2/Fl3 blocks, N = 3,4: max ratio 0.67).
NOT COVERED [O]: general linear Phi = sum_t A_t Fl_{S_t}(X) B_t with cancellation between blocks (Young/Koszul flattenings), whose
 c1, c2 may be smaller than the sum over blocks.  EGOW's own decomposition (their Lemma 3.5 applied to Phi on the fully decomposable
 4-tensors) shows every Phi is of this form with <= 16 c0 blocks (c0 = generic rank on a(x)b(y)c(u)d(v) <= min(c1,c2)), each block a
 general (non-product) compression of one flattening; but c0 can be much smaller than c1, c2 and the per-block ranks are not bounded
 by c1, c2, so Theorem B does not follow.  Partial: every Phi with c1 = 1 (or c2 = 1) is covered IF the classification "linear maps
 sending rank-<=1 to rank-<=1 are X -> AXB, AX^T B, or have rank-<=1 image" holds for singular maps [L, recalled, NOT checked].
 A counterexample to Theorem B for a general Phi would itself be a proof of rho > 1/2 (it would give c1 p + c2 q >= R >> N sqrt(c1c2)).
 Correction to section 1: G19b Theorem 4 excludes (2,2) at N = 3 (pq = 4 > N^2/4 = 2.25), so at N = 3 it ALSO beats class B;
 asymptotically (p >= N at q = 2) it is far inside.  The two known non-rank arguments (Thm 3 spectral, Thm 4 Segre) are exactly
 the ones that exceed the B-barrier, by constant factors only.

## 3. Generic frontier: join dimensions  [E: join.py, Terracini, exact rank mod 1000003]
J(p,q) = closure{A + B}.  Expected dim = min(N^4, (2pN^2 - p^2) + (2qN^2 - q^2)).  Computed for ALL (p,q), N = 2,3,4: J is
NEVER defective (dim = expected).  Least q such that J(p,q) fills End(V) (= generic operator exactly feasible, since the image of
the sum map is constructible and dense; below it generic X is not even border feasible):
  N=2: (1,2) (2,1) (3,1)            P(F) frontier: (1,1)                       [G19b]
  N=3: (1,5) (2,4) (3,3) (4,2) (5,1)        F: (1,3) (3,1)  [(2,2) infeasible, Thm 4]
  N=4: (1,11) (2,9) (3,7) (4,6) (5,5) (6,4) (7,3) (9,2) (11,1)   F: (1,6) (3,5) (4,4) (5,3) (6,1), border (3,3)
Asymptotically: generic needs 2(p+q)N^2 - p^2 - q^2 >= N^4: q = 1 -> p >= N^2 - sqrt2 N;  p = q -> p >= (1 - 1/sqrt2) N^2;
 in exponents the generic frontier is the boundary of the square, max(theta, T) = 2.  F's frontier is strictly below it already at
 N = 2 and asymptotically far below (circle T <= sqrt(4 - theta^2)).  So F is SPECIAL, not generic: "F is hard" must be proved
 from properties generic operators do not share -- which is why B-methods (blind to X by Theorem B) cannot do it.
 Where F sits: F's GL_N x GL_N-orbit {(a(x)b)F} = Z3 = rank-one tensors in the third grouping (xv|yu), dim 2N^2 - 1, all points
 of the dense orbit having F's frontier.  The question "is (p,q) feasible for F" = "is Z3 contained in J(p,q)" (border version:
 the closure).  Hardness of F = existence of a polynomial vanishing on J(p,q) but not on Z3.  Minors of B-flattenings are such
 polynomials only for pq < N^2/4 (Theorem B).

## 4. Border decompositions  [P/O]
Survive in the border sense: every linear rank method (any Phi: semicontinuity), so Theorem C(a),(b) and G18g Thm 2; the barrier
 Theorem B applies equally (it never distinguishes exact from border).  Fail: Theorem C(c) (G19e, explicit border example).
 Unknown: G19b Theorems 3, 4, 5 (eigenvalue / Segre proofs are exact-only).  Numerical evidence from G19b that Thm 3 holds in the
 border sense at small N (q = 1 tails stay at 3.0, no divergence), and that (2,2) at N = 3 is not border-feasible (tail 2.0).
 For the generic frontier exact = border (dense image).  For F: the border frontier is strictly lower at N = 4 ((3,3) border only).
 What would be needed: a border version of the spectral argument (Theorem 3) -- e.g. via the normalised limit B_eps/|B_eps| ->
 B_0 of osr <= q and rank <= p -- not attempted.

## 5. Beyond the barrier (non-rank techniques)
 (a) Equations of J(p,q) not vanishing on Z3 (algebraic geometry).  Determinantal equations of B-maps are capped at pq < N^2/4.
     Young/Koszul flattenings with cancellation (general Phi) are NOT covered by Theorem B: open whether they can pass it.
     Non-determinantal equations (higher-degree covariants) are beyond any rank method.
 (b) Dimension counts on F-specific varieties (G19b Thm 4: Segre ab^T cut by ker R).  Only proves p >= N-type statements so far.
 (c) Spectral/trace arguments (G19b Thm 3: charpoly of (G(x)G')F).  Beats the barrier by a factor 2 at q = 1; for q terms the traces
     tr((C F)^k) are "trace words" in the G_l, G'_l; no bound for q = N^c known.
 (d) Restriction from tensors with known super-linear bounds: impossible today -- no explicit tensor has a super-linear bound of this
     two-partition type (Karam arXiv:2310.11356 Q1.2, read by G19e) nor super-linear explicit tensor-rank bounds for 3-tensors
     (best ~3n, recalled), nor explicit rigidity beyond (n^2/r) log(n/r) (recalled).  F is NON-rigid relative to generic (circle vs
     square), as Walsh-Hadamard is non-rigid (Alman-Williams) -- explicit symmetric objects are easy at the generic scale.
 (e) Analogy with G19a: same phenomenon -- the flattening A'' is capped at M/u on the circular labelling, the truth unknown.

## 6. Verdict
 method                                  | best bound it proves               | barrier?
 Theorem C / G18g Thm 2 (B-maps)          | pq >= 2N^2/9 (line, constant 2/9)  | YES, Thm B: no B-method excludes pq >= N^2/4 (proved, char 0)
 any block-diag sum of compressed         | p q < N^2/4 region only            | YES (Thm B), for every operator X, exact and border
   flattenings (any bipartition, any A,B)
 general linear Phi (cancellation, Young) | nothing better known               | OPEN (counterexample = proof of rho > 1/2)
 Thm 3 spectral (q=1)                     | p >= N(N-1)/2 (exact; factor 2 > B) | none; exact only, border open
 Thm 4 Segre dimension (q=2)              | p >= N                              | none; exact only
 join dimension count                     | generic X: p+q >~ 0.59 N^2          | not applicable to F (F is special)
 restriction from known hard tensors      | nothing (no explicit super-linear) | --
Most promising: (a) -- find the lowest-degree equation of J(2,2) at N = 3 that does not vanish at F (exists: (2,2) infeasible and,
 numerically, not border-feasible), restricted to GL_3-covariants to make the interpolation feasible, and check whether it is a
 minor of a general (non-block) flattening.  If it is not, it is the first explicit non-rank witness and the template to scale.

## 7. Self-review
 * Theorem B: the Lemma's proof uses Haar measure (char 0).  The block-diagonal hypothesis is essential to the Cauchy-Schwarz step.
   The c1 <= 1 claim for Fl1 blocks and the Kronecker identities for Fl2/Fl3 were checked by test_thmC.py and schur_test.py profiles.
 * "B excludes nothing with pq >= N^2/4": combinations of methods only take unions of excluded half-planes; methods applied to a
   restricted instance (sub-block) are compressions, so inside B.  Methods that transform the instance non-linearly are not covered.
 * Join computations: one random point per (p,q) (Terracini gives a lower bound; it equals the upper bound = expected dim, so exact
   up to the probability of an unlucky random point mod 1000003).
 * Searches (search.py, lemmaS.py) are weak hill-climbs / random sampling: evidence only.
 * Literature: EGOW Thm 4.4/Lemma 3.5 read in the pdf; rank-one-preserver classification, explicit tensor/rigidity bounds recalled.
Files: notes.md, phi.py, test_thmC.py, search.py, lemmaS.py, schur_test.py, join.py; dl1/ (EGOW pdf, untrusted).
