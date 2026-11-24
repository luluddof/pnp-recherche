# G18S — adversarial review of G18e (Theorem A', Corollary A'1, small cases)
(interim save 0; AI reviewer, same model family as the author; not an independent human review)
Plan: (1) re-derive A' from definitions; (2) recompute corollary; (3) exact refutation attempts; (4) scope; (5) Lean-ready statement.

## [interim save 1] What has been checked so far
- Sources read: g18e/notes.md + all scripts/outputs; g14b/notes.md (model DF, Thm A); g17a/rapport.md; preprint
  paper.txt lines 3295-3360 (D.4.2, Lemma D.3) and 3490-3600 (D.4.5-D.4.6, (D.12)-(D.16)).
- A' re-derived from definitions: CORRECT (details in section 1 below).
- Corollary numbers recomputed by explicit wire enumeration (r1_cor.py, r1_out.txt): all of the notes' values reproduced
  (u=2: 1/4, u=3: 1/4, u=5: 3/16, u=9: 35/256, closed form C(u'-1,(u'-1)/2)/2^u', u' = u or u+1 odd).
- FOUND: the same flattening with an ORIENTED WIRE SUBSET instead of a vertex cut gives C(2u,u)/4^u (u=2: 3/8, u=3: 5/16),
  constant 1/sqrt(pi) = 0.564 instead of 0.399.  So "0.25 is this method's best at u=2,3" is only true for label-class
  vertex cuts.
- FOUND: on the single packet u=3 / 4 pairs the vertex-cut A' already gives N >= 4 (notes' table says proved LB 2 and
  "A' gives nothing more on a single packet type": wrong for u = 3).
- G18e's five explicit decompositions verified with my own verifier (r5_verify.py): all exact.
- Rank bound rank <= 2^c checked exactly over Q on 956 random restricted terms (r4_rank.py): 0 violations, 877 equalities.
- GF(2) SAT refutation attempts (r3_sat.py): no decomposition below the bound found; see section 3.

## 1. Model DF and Theorem A', re-derived
Model DF (g14b sec. 2).  M wires; wire e has two end bits (Y_e at a row vertex, X_e at a column vertex).  Target
I(Y,X) = prod_e [Y_e = X_e]  (a 0/1 tensor on 2M bits).  Directions i = 1..u, each a partition pi_i of the vertices.
A direction-i TERM is tau(Y,X) = prod_{g in pi_i} tau_g(the end bits sitting on vertices of g), tau_g arbitrary
(values in K = Frac(O)).  "Factoring over one direction's groups" = this product form, for ONE i per term.
A DF scheme: sum_{k<=N} kappa_k tau^(k) - I has all entries in 2^p O (p >= 1), kappa_k in K, term k of direction i(k).
N = number of terms (each term = r independent child contractions, so cost >= N).
Fidelity to App. D (read: D.4.5-D.4.6, rendered p.58 for (D.4),(D.6)): at T = t + alpha*zeta, i = v2(t) < u, "group vertices
on both sides by their symbol in coordinate i"; the child of group g carries "precisely the row bits Y and column bits
X of group g in the packet"; (D.16) F = sum_beta prod_g Z_{g,beta} mod 2^{p_t} holds for all vertex tensors, hence
entrywise on the wire kernel.  So App. D's root identity is a DF identity with k = (t, zeta, beta).  OK.

Flattening.  Choose a set W' of wires and, for each e in W', which of its two ends is the "row-index end" (the other is
the "column-index end").  Set every other end bit to 0.  I restricts to prod_{e in W'}[row bit = column bit] = identity
matrix of size 2^w, w = |W'| (the non-selected wires contribute [0 = 0] = 1).  This is why the flattening rank is 2^w.
A direction-i term restricts to prod_g tau_g(row-index bits located in g ; column-index bits located in g): each factor is
a 2^{e1(g)} x 2^{e2(g)} matrix, the bit sets of different groups are disjoint and exhaust the w row bits and the w column
bits, so the restricted term IS the Kronecker product of these matrices (up to the fixed identification of index
sets): rank = prod_g rank(tau_g) <= 2^{c_i}, c_i = sum_g min(e1(g), e2(g)).  Elementary proof without Kronecker rank:
for g with e1 <= e2 write tau_g(y_g, x_g) = sum_m [y_g = m] tau_g(m, x_g) (else symmetrically); the term factors as U*V
through an index set of size 2^{c_i}.
Counting.  A := restriction of sum_k kappa_k tau^(k) = Id + 2^p E, E over O.  det A = 1 mod 2O, 2O proper (O is a domain,
O/2O a field), so det A != 0 and rank_K A = 2^w.  Rank over the FIELD K is subadditive and invariant under the scalars
kappa_k:  2^w <= sum_k 2^{c_{i(k)}} <= N * 2^{max_i c_i}.
Each term chooses ONE direction; the bound 2^{c_i} depends only on the direction; hence the max over i.  CORRECT.
 * G18e's A' is the special case W' = all wires between disjoint vertex sets C1, C2, row-index end = the end in C1.
 * The vertex sets are NOT needed.  A'' (this review): for every oriented wire subset,
      log2 N >= w - max_i c_i = min_i (1/2) sum_{g in pi_i} |e1_i(g) - e2_i(g)|      (since sum_g e1 = sum_g e2 = w),
   and more precisely 2^w <= sum_k 2^{c_{i(k)}}.  A'' is strictly stronger (see 2b, 3).
 * Hypotheses actually used: p >= 1; the error is entrywise in 2^p O; each term is a product over the groups of one
   partition of the END BITS (a vertex in no group is harmless; a vertex in two groups of one direction is not allowed).
   Not used: bipartiteness, labels, uniformity, what the vertices' tensors are.
VERDICT CLAIM A': CONFIRMED (proof re-derived; 956 exact random rank tests r4_rank.py; no counterexample by exact
GF(2) search r3/r7).

## 2. Corollary A'1, re-derived and recomputed
Labelling of App. D (Lemma D.3, checked on rendered p.58-59): both sides labelled in [r]^a, r = 64, a = 2 sum_{j<=d} u_j >= 2u;
pair class gamma = (row label, column label), D_arr = r^{2a} classes, loads X_gamma with sum_gamma (X_gamma - M/D_arr)^2
<= 2 M Delta; (D.6) M >= 2 D_arr^2 Delta.  Directions = coordinates i = v2(t) in 0..u-1, r groups each.
Cauchy-Schwarz: for any coefficients |c_gamma| <= 1:  |sum_gamma c_gamma (X_gamma - M/D_arr)| <= sqrt(D_arr * 2 M Delta)
<= M / r^a.  OK (as in the notes).
Majority cut (u odd), b_i(x) = [x_i in S_i], s(x) = sum_i b_i(x), C1 = {s > u/2}, C2 = {s < u/2}:
 w = M/2 (+- M/r^a);  w - c_i >= #(cut wires, C1-end in S_i) - #(cut wires, C2-end in S_i)   [only 2 unions of classes]
 = w * ((1/2 + eps) - (1/2 - eps)) = 2 eps w = eps M,  eps = P(Bin(u-1,1/2) >= (u-1)/2) - 1/2 = C(u-1,(u-1)/2)/2^u;  error <= M/r^a.
 Recomputed by explicit enumeration of all wires of the uniform instance (r1_cor.py; r = 2 and r = 4 give identical values):
   u        : 2     3     4      5      6       7       8        9        10       11       12       41      401
   notes    : .25   .25   .1875  .1875  .15625  .15625  .13672   .13672   .12305   .12305   .11279   .06269  .01993
   mine     : 1/4   1/4   3/16   3/16   5/32    5/32    35/256   35/256   63/512   63/512   231/2048 .06269  .01993   (identical)
   A'' (new): 3/8   5/16  35/128 63/256 231/1024 429/2048 6435/32768 .18547 .17620  .16819   .16118   .08784  .02817
 Even u: the cut {s >= u/2} | {s < u/2} gives the value of u+1 (verified u = 2..12).  So "u odd" is not needed.
 M large: needed only through (D.6) (error M/r^a <= M/64^{2u}); otherwise the statement is exact.
 Constant: C(u-1,(u-1)/2)/2^u ~ 1/sqrt(2 pi u) = 0.3989/sqrt u.  CONFIRMED.
VERDICT COROLLARY: CONFIRMED for Lemma D.3's labelling (any network satisfying (D.6)); "Theta" is CONDITIONAL on App. D's
upper bound, which nobody here has verified.

### 2b. Improvement found: constant 1/sqrt(pi) = 0.564 (A'')
Select every wire whose two end labels have s(x) != s(y); row-index end = the end with the LARGER s.  Then for each i
 w - c_i >= sum over classes (x,y) of load * sign(s(x)-s(y)) * (b_i(x) - b_i(y)) = (M/u) E|s(x) - s(y)| = M * C(2u,u)/4^u,
so on App. D's labelling  log2 N >= M (C(2u,u)/4^u - r^{-a}) ~ 0.564 M/sqrt(u)  (u = 2: 3/8, u = 3: 5/16).
This is NOT a vertex cut (a vertex is a row-index end for some of its wires and a column-index end for others).
Optimality of the family on exactly uniform product loads (r6_brute.py): for 2 groups per direction
w - c_i = |sum_c delta_c (x_i - y_i)|, the maximisation is an LP after fixing signs, and symmetrising over coordinates gives
the optimum M C(2u,u)/4^u exactly (LP solved for u <= 5: 2/4, 6/16, 20/64, 70/256, 252/1024; exact integer DP u <= 3).
For general r and any orientation:  min_i <= average_i <= (M/u) E|T - T'| <= M/sqrt(2u), with equality case p_i = 1/2,
so the optimum is M C(2u,u)/4^u for every even r.  Hence remark (i) of G18e sec. 2 ("the family is capped at O(M/sqrt u)")
is now PROVED on exactly uniform loads, and the flattening method cannot give more than 0.564 M/sqrt(u) there.
VERDICT on "brute force over all cuts confirms 0.25 is this method's best at u = 2, 3": PARTLY.  True for vertex cuts that
are unions of label classes (re-run: 4/16 and 16/64, also when row and column vertices of a class are cut separately);
FALSE for the method: oriented wire subsets give 6/16 and 20/64, and these are the method's exact optima.

### 2c. Non-uniform loads: the sqrt(u) law is NOT a theorem of DF in general
The corollary uses near-uniform coarse loads (the 4^u sign classes).  It is not a max over cuts for every admissible
load, and it cannot be: put all wires on pair classes at Hamming distance 1 (x, x + e_i), M_i wires flipping coordinate i.
Then direction i is crossed by M_i wires only, I = sum over the 2^{M_i} values of those wires of a direction-i term, so
N <= 2^{min_i M_i} <= 2^{M/u}; and A'' gives N >= 2^{min_i M_i}: exact.  All groups are balanced.  So with loads chosen by the
designer DF reaches 2^{M/u} -- on networks that HAVE a balanced partition cutting M/u wires.  For App. D's hashed
labelling this is excluded by Lemma D.3; for an arbitrary labelling of an expander-like network no bound of the form
2^{Omega(M/sqrt u)} is proved (G18e's remark (ii): level E; the obstacle is exchanging "min over lambda" and "expectation
over random signs" when a direction has only r groups).

## 3. Refutation attempts and the small-case table
Exact search = SAT over GF(2): terms with 0/1 tables, kappa = 1, p = 1 (a SUB-model of DF; SAT => a genuine DF identity for
p = 1; UNSAT proves nothing about DF, where kappa in Frac(O) and tables in O).  dflib.py, r3_sat.py, r7_refute.py, r8, r10.
 * N = 2^bound - 1 on all seven instances of G18e's table: UNSAT for every direction profile.
 * 129 random instances (u = 2, 3; M <= 5) where A'' is strictly stronger than every vertex cut (e.g. the triangle
   (0,1),(0,3),(1,3) at u = 2: vertex cuts give N >= 2, A'' gives N >= 4): N = 2^{A''} - 1 UNSAT in all; N = 4 SAT on the
   triangle (so A'' is tight there).
 * No counterexample to A' or A''.
G18e's five explicit integer decompositions: verified with my own verifier (r5_verify.py: exact rank over Q across the two
groups <= 1 for every term, integer sum == I): 2, 4, 7, 4, 7 terms.  CONFIRMED.
Table (log2 N lower bounds by brute force over ALL vertex cuts / ALL oriented subsets, r2_small.py):
 instance                  M  ThmA  A'  A''  proved N>=   exact UB over Z   GF(2) sub-model (p = 1)
 u=1 both slots            2   2    2   2      4            4                min 4
 u=2 slots{00,01}          2   1    1   1      2            2                min 2
 u=2 slots{00,11}          2   2    2   2      4            4                min 4
 u=2 slots{00,11,01}       3   2    2   2      4            4                min 4
 u=2 full packet           4   2    2   2      4            7                N<=5 UNSAT for all profiles (r3, r8, r11; (3,2) took 254 s), so min in {6,7}
 u=3 three pairs           3   1    1   1      2            4                min 4 (N <= 3 UNSAT for all profiles)
 u=3 four pairs            4   1    2   2      4  (notes: 2) 7               N<=4 UNSAT; N = 5 SAT, profile (2,2,1), verified
Corrections to G18e's table and text:
 (1) "u=3, 4 pairs: proved LB 2" -> 4.  The majority cut C1 = {0,1,2,4}, C2 = {7,6,5,3} has w = 4, c_i = 2 for each i.  So
     "A' gives nothing more on a single packet type" (sec. 3, 4) is FALSE for u >= 3: on the full antipodal packet the
     majority cut gives log2 N >= M C(u-1,(u-1)/2)/2^{u-1} ~ 0.80 M/sqrt(u) against Thm A's 2M/2^u.  The sentence "the gain
     of Thm A' comes entirely from mixing packet types" is wrong; the HWY remark only concerns u = 2.
 (2) N = 5 at u = 3 / 4 pairs is not only a real-number border phenomenon: there is an EXACT identity mod 2 with 5 terms
     (n5_gf2_solution.json; it is not an identity over Z, and no sign choice makes it one: r10).  So for p = 1:
     4 <= N <= 5; for an identity over Z the verified upper bound stays 7.  Not lifted to p >= 2 (not attempted).
 (3) The numerics ("N=3 not found" etc.) are real-number LM searches: evidence only, as G18e says.
VERDICT SMALL CASES: PARTLY (all upper bounds and four of the lower bounds confirmed; one proved lower bound is not the
best the author's own theorem gives; one false general remark).

## 4. Scope
 * What N bounds.  The number of terms of ONE kernel identity at ONE call.  In App. D terms = (node t, zeta, beta); there
   are 2L nodes and poly many zeta, so at some retained node the number of mode tuples that must be enumerated is
   >= 2^{0.564 M/sqrt u}/poly, whatever the taper, the mode basis, the tiles, the recovery rule or the precision (all of
   these only change WHICH direction-factored terms are used).  This turns G18d's design-level Prop. B' (2^{0.3 M/sqrt u},
   "for the algorithm as written") into a model-level statement, and supplies the "rank-type bound in DF" that G17a
   listed as missing.  It also answers G17a's link (ii) for every scheme that ENUMERATES group-factored terms.
 * Running time.  A' bounds time only for algorithms that spend >= 1 step per term.  Children could be cheap (not
   needed: the count alone is the bound) -- but a scheme that does NOT enumerate (re-sums beta by a further contraction,
   i.e. treats modes as wires of a new network of r super-vertices) is outside; A'' then only says that this new
   network has >= 0.564 M/sqrt(u) bits of bond dimension across the oriented cut.  Multi-level gains are untouched.
 * Consistency with G18d ("root dominates; design-level bounds only"): since the recursion is root-dominated, a root
   bound of the right order is all one can ask of a one-level model; A'' gives it at the model level, with the SAME
   order as App. D's count (constants: lower 0.564, upper about 38-44 in log2 units per M/sqrt u -- a factor ~70).
 * Family for which "sqrt(log) is optimal" is now a theorem (lower-bound side):
   F = { algorithms that, at a call with M wires, (a) label the vertices by a product [r]^a (r even) and use as
   directions u <= a coordinates, with pair-class loads satisfying sum_gamma (X_gamma - M/r^{2a})^2 <= 2 M Delta and
   M >= 2 r^{4a} Delta  [Lemma D.3 + (D.6); only the 4^u coarse classes are needed];  (b) use ANY identity
   I == sum_{k<=N} kappa_k tau^(k) mod 2^p, each tau^(k) a product over the groups of one direction, valid for all vertex
   tensors;  (c) pay >= 1 per term }.   Then time >= 2^{M (C(2u,u)/4^u - r^{-a})} >= 2^{(0.564 - o(1)) M/sqrt(u)}, and since
   (D.6) gives u <= a <= log2(M)/(4 log2 r):  time >= 2^{Omega(M/sqrt(log M))}.
 * Still assumed / not covered: (1) App. D's upper bound (its correctness is unverified) for the word "optimal";
   (2) u = O(log M) is a property of the packing (D.6), not of DF: A'' gives M/sqrt(u) for any u and cannot exclude
   2^{o(M/sqrt(log M))} with many more directions and another packing (no such scheme is known);
   (3) near-uniform loads: FALSE in general DF (sec. 2c): on networks with small balanced separators a designer gets
   2^{M/u}; for arbitrary labellings of expander-like networks nothing is proved;
   (4) terms that do not factor over ONE direction's groups, non-enumerative re-summation, multi-level schemes,
   identities depending on the vertex tensors;  (5) nothing here is a statement about 3-SAT or ETH.
 * Side estimate (level E, constants unverified): the bound means the root count exceeds 2^n unless
   sqrt(u) > 0.564 M/n; with M/n in the tens this needs u in the hundreds to thousands, i.e. log2 M >= 48 u: astronomically
   large inputs.  Same conclusion as G18e's reading (i), now from the lower-bound side.

## 5. Lean-ready statement (not formalised)
Core (pure linear algebra over a field; kappa_k is absorbed into one factor of term k):

  theorem df_flattening
      {K : Type*} [Field K]
      {W : Type*} [Fintype W] [DecidableEq W]                      -- selected wires
      {ι : Type*} {G : ι → Type*} [∀ i, Fintype (G i)] [∀ i, DecidableEq (G i)]   -- directions, groups
      (ρ γ : ∀ i, W → G i)              -- group (direction i) of the row-index end / column-index end of each wire
      (σ : ∀ i, G i → Bool)             -- per group: true = "factor through its row bits", false = "through its column bits"
      {N : ℕ} (dir : Fin N → ι)
      (τ : ∀ k : Fin N, G (dir k) → (W → Bool) → (W → Bool) → K)
      (hloc : ∀ k g y y' x x', (∀ e, ρ (dir k) e = g → y e = y' e) → (∀ e, γ (dir k) e = g → x e = x' e) →
                τ k g y x = τ k g y' x')                            -- τ_g reads only the bits located in g
      (A : Matrix (W → Bool) (W → Bool) K)
      (hA : ∀ y x, A y x = ∑ k, ∏ g, τ k g y x)
      (hdet : A.det ≠ 0) :
      2 ^ Fintype.card W ≤
        ∑ k, 2 ^ (Fintype.card {e // σ (dir k) (ρ (dir k) e) = true} + Fintype.card {e // σ (dir k) (γ (dir k) e) = false})

  (choosing σ_i(g) = true iff e1_i(g) ≤ e2_i(g) gives the exponent c_i = ∑_g min(e1_i(g), e2_i(g)); then RHS ≤ N · 2^{max_i c_i}.)

2-adic lemma:  {R} [CommRing R] [IsDomain R] (I : Ideal R) (hI : I ≠ ⊤) (B : Matrix n n R) (h : ∀ i j, B i j - (1 : Matrix n n R) i j ∈ I)
               ⊢ (B.map (algebraMap R (FractionRing R))).det ≠ 0        (R = O, I = (2); only "2 is not a unit" is used).
Restriction lemma: setting the end bits of non-selected wires to false maps I to the identity matrix and a direction-i term
to a function of the above form (bookkeeping only).
Counting identity for the corollary:  ∑ x y : Fin u → Bool, |(#true x : ℤ) - #true y| = u * (2u).choose u.
Mathlib notions: Matrix.rank over a field; rank (U * V) ≤ card of the middle index (Matrix.rank_mul_le_left/right with
rank_le_card_width/height) -- this REPLACES any "rank of a Kronecker product" lemma (the term is U*V with middle index
(∀ e : {row ends in true-groups} ⊕ {column ends in false-groups}, Bool));  subadditivity of rank over a Finset sum
(from LinearMap.range_add_le / Submodule.finrank_add_le_finrank_add_finrank if Matrix.rank_add_le is not available);
full rank from det ≠ 0 (Matrix.isUnit_iff_isUnit_det, Matrix.rank_of_isUnit);  RingHom.map_det and Ideal.Quotient for the
2-adic lemma;  IsFractionRing.injective.  If one prefers the Kronecker route: Matrix.kroneckerMap plus
rank (A ⊗ₖ B) ≤ rank A * rank B (to be proved if absent from Mathlib; I did not check the current library).

## 6. Verdicts
 CLAIM A' (theorem): CONFIRMED.  Each term picks one direction; rank over the field Frac(O); Kronecker structure exact;
   subadditivity; det = 1 mod 2.  Generalises verbatim to oriented wire subsets (A''), which is strictly stronger.
 COROLLARY (0.399 M/sqrt u on App. D's labelling): CONFIRMED, all exact values reproduced; valid for even u too; needs only
   Lemma D.3's variance bound + (D.6).  Constant improvable to 0.564 (A''), and 0.564 is the exact limit of the method on
   uniform loads.  "Truth in DF is 2^{Theta(M/sqrt u)}": PARTLY -- lower bound proved for near-uniform product labellings
   only (false for designer-chosen loads on separable networks: 2^{M/u} exactly), upper bound = App. D at face value.
 SMALL CASES: PARTLY -- decompositions and four lower bounds confirmed; "u=3, 4 pairs: LB 2" should be 4 by A' itself;
   "A' gives nothing more on a single packet" false for u >= 3; "0.25 is this method's best at u=2,3" true only for
   label-class vertex cuts (method's optimum: 3/8, 5/16).  New: exact 5-term identity mod 2 at u=3 / 4 pairs.
 No error found that weakens a claimed lower bound; every correction goes in the direction of a stronger statement,
 except the scope restriction of sec. 2c.
Not checked here: the HWY/Hurwitz statements of G18e sec. 3; G18e's real-number numerics; App. D's correctness.

Files (all in this folder): notes.md; dflib.py; r1_cor.py/r1_out.txt (corollary values); r2_small.py/r2_bounds_out.txt;
r3_sat.py/r3_out_*.txt; r4_rank.py/r4_out.txt; r5_verify.py/r5_out.txt (+ g18e_two_block.py, copied from g18e/exact.py);
r6_brute.py/r6_out.txt; r7_refute.py/r7_out.txt, r7b_out.txt; r8_n5.py/r8_out.txt; r9_hamming.py/r9_out.txt;
r10_n5.py/r10_out.txt, n5_gf2_solution.json; r11_full.py/r11_out.txt; r12_packet.py/r12_out.txt; pdfimg/p-58.png, p-59.png.
