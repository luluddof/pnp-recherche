# G18f — E1 (entangled terms), first non-trivial size.  Tags: [P] proved here (unreviewed), [N] numerical, [E] sketch, [O] open


## 1. E1 restated (all definitions)
* Grid g x h of nodes (a,j).  A lam-bit string on every DIAGONAL edge {(a,j),(a',j')}, a != a', j != j' (N = 2^lam values).
  T = T_{g,h,lam} = product of the equality tests of all diagonal edges = EPR pairs on the graph K_g x K_h (tensor product of
  graphs); every node holds (g-1)(h-1) strings.  S = lam * C(g,2) * h(h-1) = total diagonal bits.  Rectangles = (2 rows, 2 columns);
  each carries 2 diagonals, i.e. one crossing [x=v][y=u].
* Expansion: c T = sum_{i<=m} A_i, every A_i typed "row" or "col".  cost_row(A) = total bond bits of the cheapest tensor network on
  the g merged rows (complete graph K_g, pairwise bonds only) representing A; cost_col likewise on the h merged columns.
  kappa_r = max over row-typed terms of cost_row/S; kappa_c likewise.  T itself has cost S in both directions.
* E1: is there, for some fixed grid with g or h >= 3 and all large lam, an expansion with m = 2^{o(S)} terms and
  kappa_r + kappa_c <= 1 - delta ?   (kappa_r + kappa_c >= 1 for all expansions <=> harmonic law u(gh) <= u(g)+u(h) <=> 2^{n/log n} wall.)
* g = 2: cost_row = log2 of the rank across the single row cut, rank is subadditive, so the row terms merge into ONE A at a cost of
  log2 m bits.  Column terms do NOT merge for h >= 3: the sum of m networks is a network whose bond (j,k) has sum_i D^i_{jk}, so its
  cost is sum over bonds of the per-bond max, not the max of the sums.  (I first wrote a "two-term reduction"; it is wrong for h>=3.)
* Operator form (g=2, h=3).  Diagonal graph = hexagon C_6; T = Tr(M_0...M_5) = iterated matrix multiplication of 6 N x N matrices.
  Sites = columns; site j has 2 input strands (node (0,j)) and 2 output strands (node (1,j)).  T = the permutation operator
  P = S_01 (x) S_02 (x) S_12, one SWAP_N per rectangle, each pair of sites exchanging one strand.  Row terms = low-rank operators,
  column terms = operators with a triangle (ring-of-3) bond structure.  As a 3-party tensor over the columns T = <N^2,N^2,N^2>
  (matrix multiplication), and the row cut is the split in|out INSIDE each party.

## 2. What a positive answer at lam = 2 would be, and what it would imply
* At fixed lam the inequality log2 r + max_i cost(B_i) < S is NOT a positive answer: it holds trivially at small size.
  lam=1 (S=6): P = (rank 1) + three "stars" of cost 4 bits each (found numerically, and = the product of the one-bit identity
  S_2 = x - Y, Y of rank 1, over the 3 rectangles): 4 of 6 bits.  Also B = identity gives rank(P - 1) = (N^6 - N^3)/2: a 1-bit gain
  at every lam.  These are O(1)-bit, finite-size effects.
* E1 <=> gain(lam) := S - min[log2 r + max_i cost(B_i)] >= delta*S for all large lam with 2^{o(S)} terms.  It is an ASYMPTOTIC statement
  about the Kronecker powers T_lam = (T_1)^{(x) lam} of one fixed tensor (T_1 = hexagon of qubit EPR pairs, 12 qubits).
* No amplification [P for a two-term identity, E in general]: if T_{lam0} = A + B, the t-th tensor power regrouped by the number i
  of A-factors, C_i = sum_{|I|=i} A^{(x)I} (x) B^{(x)I^c}, satisfies rk_row(C_i) >= rk(A)^i (D - rk A)^{t-i} (restrict the inputs
  outside I_0 to ker(M_T^{-1} M_A): only the term I = I_0 survives and B acts as the identity there), and symmetrically every column
  cut of C_i has rank >= (D_c - b)^i.  So at i = theta t the row cost is >= (1-theta)(S - 1) t and the column cost >= theta (S - 1.5) t:
  kappa_r + kappa_c >= 1 - O(1/S_0) whatever the base identity.  Reason: the cost is an infimal convolution of two gauges (row rank,
  column network); each is submultiplicative, but a mixed product A (x) B is expensive for both.  Unlike tensor rank, there is no
  Fekete/Schonhage amplification.  A finite gain at lam=2 would imply NOTHING asymptotically; it could only suggest an infinite
  algebraic family (as the Hurwitz/Clifford families do for one crossing).
* So a meaningful "positive at lam=2" = an expansion with log2 r + max_i cost_i < 12 that is not a product of crossing
  decompositions (products give >= 14.3 bits at N=4 from the known points (1,6),(3,3)border,(4,4)); by Theorem F below it must have
  r <= 15 and column terms whose sum has full cut rank 256 at every column.

## 3. Structural result: a cut-rank product inequality  [P, unreviewed; identities checked numerically in check_identities.py]
THEOREM F.  2 x h grid, any N, any field.  Let P = A + B (exact, or border: P = lim(A_e + B_e)), r = rank of A across the row cut,
b = rank of B across the cut (column j | other columns).  Put M = N^{h-1}, Q = N^{(h-1)(h-2)} (Q M^2 = N^{h(h-1)} = 2^S).  Then
      r >= s Q (M - b s)   for every integer s >= 1 (needs h >= 3)  and for every s = 1/rho, rho in {1,...,M}  (h >= 2).
Hence  r*b >= (2/9) 2^S  whenever b <= M^2/2,  and  r >= Q (1 - b/M^2)  always.   (B may be any sum of column terms.)
Proof, primal part (h=3 notation; site 0 = column j with strands a,a'; the strands b,c of sites 1,2 move into site 0 under P).
 U = ker A has codimension <= r and Pv = Bv on U.  For any operator B of operator-Schmidt rank b across 0|12,
 SR_{0|12}(Bv) <= b SR_{0|12}(v); and SR_{0|12}(Pv) = SR_{{b,c}|rest}(v) because P only permutes strands.  Take Phi_s, an
 s-dimensional subspace of H_0, and L = Phi_s (x) H_12.  Every v in U ∩ L, read as a matrix H_{bc} x (Phi_s (x) H_{b'c'}) of size
 M x sQ, has rank <= bs.  U ∩ L is a LINEAR space of dimension >= sMQ - r.  Flanders' theorem (a linear space of m x n matrices of
 rank <= k has dimension <= k max(m,n); all fields by Meshulam): sMQ - r <= bs * sQ.
Proof, dual part.  Write B = sum_{k<=b} L_k (x) M_k and K = {w in Op(H_0): <w,L_k> = 0 for all k} (codim <= b in M^2).
 For w in K: Tr_0[(w^+ (x) 1) P] = Tr_0[(w^+ (x) 1) A].  Left side = (w^+ transplanted to the strands b,c) (x) (P of the other
 columns), image = Im(w^+) (x) C^Q.  Right side: if Im(w) lies in a rho-dimensional F, its image lies in the FIXED space
 X_F = span{(<f| (x) 1) x_i : f in F, i <= r} of dimension <= rho r (A = sum |x_i><y_i|).  K_F = {w in K: Im w in F} has
 dimension >= rho M - b; G_F = sum_{w in K_F} Im(w^+) satisfies Q dim G_F <= rho r, and K_F lies in F (x) G_F^*, so
 rho M - b <= rho dim G_F <= rho^2 r / Q.   Border case: take limits of the subspaces (Grassmannians are compact, rank conditions
 are closed).   QED
Why it was missed: the known floor r + b >= N^{p(h-p)} is the case s = 1 restricted to one orthant; the gain comes from letting s
grow (resp. rho) and using Flanders instead of a single identity minor.
Tightness: h=2, N=2: r >= 2 - b gives (r,b) = (1,1), attained by delta delta = delta delta - eps eps.  Hexagon N=2: b=16 -> r>=0,
b=1 -> r >= 16 (best found 28).
COROLLARIES.
 (C1) One crossing (2 x 2; here everything merges, so this is fully general): r >= max_rho (rho N - b)/rho^2 and the mirror
      image; so r b >= (2/9) N^2, i.e. kappa_r + kappa_c >= 1 - O(1/lam) for EVERY expansion of a crossing into 2^{o(lam)} typed terms.
      This closes the item "[O] a 2x2 point with kappa_r + kappa_c < 1 is not excluded" of G14e (it is excluded), and is
      consistent with all data there (the bound gives N=3: b=1 -> r>=2; N=4: b=1 -> r>=3, b=3 -> r>=2; known points (1,3), (1,6), (3,3)).
      It says nothing on Naslund's problem (r = b ~ N/2 is allowed).
 (C2) 2 x h, column side = ONE network (or any B whose cut at some column j is not full): log2 r + cut_j(B) >= S - 2.17, hence
      kappa_r + kappa_c >= 1 - o(1); for one network even kappa_c >= (h/2)(1 - kappa_r).
 (C3) 2 x h, general: for every column j,  r >= Q (1 - b_j(B)/M^2), B = sum of all column terms.  So kappa_r + kappa_c < 1 - delta
      forces r < Q = 2^{(1-2/h)S} (kappa_r < 1 - 2/h) AND b_j(B) >= M^2 (1 - r/Q) for EVERY j: the column terms together must
      have essentially full rank at every single-column cut, so kappa_c >= 2/h.  For 2 x 3: the only region not excluded is
            kappa_r < 1/3   and   2/3 <= kappa_c < 1 - kappa_r        ("the corner").
      Minimal shape in the corner: P = (rank < N^2) + B_0 + B_1 + B_2, B_j a "star" centred at column j (bonds N^2, N^2, 1).
 (C4) g >= 3: restricting to the edges that cross a row bipartition R gives the same inequality with S replaced by S_R (the bits
      crossing R) - about S/2 for balanced R.  So on the grids the algorithm uses (g = h large) Theorem F only gives a floor of
      the type kappa_r + kappa_c >= ~1/2; E1 for g,h >= 3 is untouched except through these cut pairs.
* Product constructions never reach the corner asymptotically: with crossing pieces S = X + Y (OSR bits kappa_x, rank bits kappa_y,
  fractions of 2 lam) the 2x3 product gives stars of cost (2+kappa_x)/3 and A = Y(x)Y(x)Y with kappa_r = kappa_y; the sum is < 1
  iff 3 kappa_y + kappa_x < 1, impossible by (C1) unless kappa_y = 0.  (At N=2, Y has rank 1: that is the 4-of-6 point.)
* Schmidt-rank arguments are provably exhausted in the corner: there every single-column cut of B is full, and both halves of the
  proof become vacuous.  A proof would need a measure of the TOTAL triangle cost that is subadditive over terms.

## 4. Is E1 a border-rank / asymptotic-rank question?  No - and the reason is matrix multiplication.
* The tensor: T_{2,3,lam} = hexagon = <N^2,N^2,N^2> over the three columns (general: EPR pairs on K_g x K_h).  E1 asks whether T
  lies in a JOIN of tensor-network varieties: (one K_g-network variety on the rows) + (several K_h-network varieties on the columns),
  with total bond bits < S.  Membership in a join of two different varieties is not a rank: no submultiplicativity (section 2).
* The natural subadditive relaxation of the column cost is the 3-party tensor rank over the columns (a triangle network with bonds
  D has rank <= D01 D12 D02).  In that relaxation the answer is trivially "positive": A = 0, B = T has rank N^{2 omega} =
  2^{(omega/3) S} < 2^S.  This is exactly Lemma T of G9d (one hyperedge of omega*L bits: the instance leaves the graph class), and
  it shows that NO lower-bound method that only sees a rank-like (subadditive, hyperedge-insensitive) column measure can close
  the corner: fast matrix multiplication is the obstruction to the proof, without being a mechanism for the algorithm.
* Rigidity reading: E1(2x3) = "is the permutation operator P_N = P_2^{(x) lam} non-rigid with respect to triangle operators":
  rank(P - B) <= N^{6 alpha}, cost(B) = 6 beta lam, alpha + beta < 1.  The Alman-Williams/Dvir-Liu mechanism for Kronecker powers
  (exact on typical sectors of the +-1 eigenspaces of the 3 lam swaps, low-degree polynomial in the number k of antisymmetric
  factors) is a genuinely entangled term; it gives alpha = 1 - Theta(eps^2), beta = H(2 eps)/2 + 2 eps, i.e. alpha + beta > 1 [E].
* Classical (diagonal) terms: if A, B are supported on the support of T (invariant under the bond tori), B = P restricted to a set X
  of strand configurations and Shearer's inequality gives |X| <= D01 D12 D02, so r + prod D >= N^6: only the trivial points [P].

## 5. Computations (all scripts here; single-thread numpy; each run < 10 min)
* hexlib.py: P for any N; triangle operators; objective f_r = sum_{k>r} sigma_k(P - sum_i B_i)^2 (A eliminated by Eckart-Young),
  exact gradient (checked by finite differences), L-BFGS with random restarts.  "Found" = residual < 1e-9 and a gap sigma_r/sigma_{r+1}
  ~ 1e10; "not found" = evidence only.
* check_identities.py (exact): SR_{0|12}(Pv) = SR_{bc|rest}(v); rank Tr_0[(w^+ (x) 1)P] = rank(w) N^2; tables of the proven bound.
  Hexagon N=2: proven min r for b = 1,2,3,4,..,16: 16, 8, 5, 4, 4, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0.   N=4: b=1 -> 1024, 16 -> 64, 64 -> 16, 255 -> 1.
* lam = 1 (N=2, 64 x 64), one column network with bonds (D01,D12,D02)  [lam1.py, multi.py, prodB.py]:
  (1,1,1): r = 28 found (B = identity up to gauge; rank(P - 1) = (64-8)/2); r = 20: not found (residual 18.7, 30 starts); proven r >= 16.
  (2,1,1): r = 28 found; 25 not found.  (2,2,1),(3,1,1),(4,1,1),(3,2,1),(4,2,1),(3,3,1),(4,4,1),(3,2,2),(4,2,2),(3,3,2), and (2,2,2)
  with r = 5, 6: nothing found with r * prod(D) < 64 (residuals 5 to 50).  The optimiser is weak (it missed feasible points
  with larger r in the same sweep), so these are not infeasibility statements; none of them contradicts Theorem F.
* lam = 1, the corner: three stars (4,1,4),(4,4,1),(1,4,4) [multi.py]: r = 3, 2, 1 FOUND at the first start; r = 0 not found
  (residual 25; and r = 0 is impossible by the projection argument of section 6).  r = 1 is the product point (Y (x) Y (x) Y).
* lam = 2 (N=4): the generic problem has 2^24 tensor entries; one objective evaluation needs an SVD of a 4096 x 4096 matrix
  (measured: 1.4 s for 1024 on one core, so ~90 s; ~3 s with a top-r Lanczos solver when r <= 15).  A search of the corner
  (three stars with bonds 16: 3 x 73 728 unknowns, r <= 15) needs ~10^3 evaluations per start: about 1 hour per start, several
  days for a meaningful number of restarts.  NOT run.  What is decided at lam = 2 is decided by Theorem F (section 3).
* Proxy at N = 3 (729 x 729; the same question, between lam=1 and lam=2): three stars with bonds 9, target r < N^2 = 9; the product
  construction gives r = 27.  [multiN_fast.py; results below]
  r = 8: three starts of 1500 L-BFGS iterations end at residual^2 = 591, 682, 694 (of |P|^2 = 729): NOT converged, nowhere near
  even the truncated product point.  Calibration on the known feasible r = 27 (and r = 60): see n3_r27_cal.out / n3_r60_cal.out
  (700 iterations, one start): residual^2 = 621 at r = 27 and 571 at r = 60, although r = 27 IS feasible (product point).
  So the N = 3 search is INCONCLUSIVE: the optimiser fails its calibration; nothing can be read from r = 8.

## 6. The corner: what is known
* Minimal open instance: is P_N = (rank r) + B_0 + B_1 + B_2 with B_j a star centred at column j with bonds N^2, and r <= N^{2-delta}?
  (cost 4 lam + log2 r of 6 lam).  More generally any column terms of cost < 6 lam - log2 r.
* r = 0 is impossible for N >= 2 [P]: over the columns P = <n,n,n>, n = N^2.  A star centred at 0 has mode-1 and mode-2 images of
  dimension <= n.  Contract mode 1 with the phi that kill the mode-1 images of the stars centred at 0 and 2 (codimension <= 2n in
  n x n matrices): <n,n,n>(phi) = T_1(phi) must have mode-0 image in a fixed n-dimensional space, but the mode-0 images of
  <n,n,n>(phi) span C^n (x) colspace(Phi), of dimension >= n(n-2) > n.
* N = 2: r = 1 exists (product point), so min r = 1 there.  N >= 3: unknown between 1 and the product value (N(N-1)/2)^3.
* "Kill and compress" (same idea with A present) [E]: if every column term can be killed at some column with total killed cut
  rank beta_j <= N^{2-eps} at column j, then r >= ~N^{3 eps}/8.  For stars of bond exactly N^2 (beta_j = N^2) it gives nothing:
  the stars sit exactly on the threshold of every annihilator/cut-rank argument I tried.
* So neither a proof nor a construction in the corner.  Any proof must distinguish pairwise bonds from a GHZ hyperedge (section 4).

## 7. Verdict
* NEGATIVE at lam = 2 and at every lam, for the part of E1 that can be stated with cut ranks:
  (a) one crossing (2 x 2): kappa_r + kappa_c >= 1 - O(1/lam) for all expansions [P];
  (b) 2 x h: kappa_r + kappa_c >= 1 - o(1) unless kappa_r < 1 - 2/h and the column terms have full rank at every column cut [P];
  (c) 2 x 3 with one column network: kappa_c >= (3/2)(1 - kappa_r) [P];
  (d) at lam = 2 on the 2 x 3 grid a gain needs r <= 15 and column terms of cost < 12 - log2 r bits whose sum has cut rank
      >= 256 - 16 r at each column (Theorem F with s = 1/16).
* UNDECIDED: the corner kappa_r < 1/3 <= 2/3 <= kappa_c of the 2 x 3 grid (three or more column terms), and everything for g,h >= 3
  beyond a floor ~1/2.  No mechanism found; no finite computation can give a positive answer (no amplification).
* Exact next computation: N = 3 (and then N = 4) "three stars + rank r": decide min r.  N = 3: 729 x 729, 24 000 unknowns; needs a
  solver better than plain L-BFGS (Gauss-Newton on the bilinear system P U = (sum B_i) U with U a 729 x (729 - r) basis of ker A,
  or ALS: each star tensor enters linearly), ~1 CPU-hour.  N = 4: 4096 x 4096, 2.2e5 unknowns, r <= 15: ~1-3 CPU-days.
  A value min r >= c N^2 at N = 3, 4 would support the conjecture "min r >= N^{2-o(1)}", which (with Theorem F) would close 2 x 3.
* Next proof step: Theorem F type inequality for the pair (row rank of A, TOTAL bond bits of one column term), i.e. a
  bound on the triangle-network cost that is subadditive over terms and not fooled by <n,n,n> having rank n^omega.

## 8. Adversarial self-review
* Theorem F: proof written in one sitting, not reviewed.  Checked: the two identities numerically (N=2,3); tightness at the
  crossing N=2; consistency with every known point of G14e (28;1), (10;7), (1;16), crossing (1,6),(3,3),(4,4),(8,8); no numerical
  counterexample at N=2 (weak test: the optimiser rarely finds anything).  Flanders' theorem is quoted, not re-proved
  (statement used: a linear space of m x n matrices of rank <= k has dimension <= k max(m,n); Flanders 1962, Meshulam 1985).
  The constant 2/9 comes from rounding s to an integer or a reciprocal integer.  The border-case limit argument is a sketch.
* The first version of this note claimed a full negative answer for 2 x 3 through a "two-term reduction"; that reduction is false
  for h >= 3 (column terms do not merge) and the claim was withdrawn.  (C2) is the correct remnant.
* The no-amplification statement is proved only for a two-term base identity and only up to O(1/S_0); for m > 2 base terms the
  rank of a type-class sum is not controlled [E].
* Cost model inherited from G9d/G14e (balanced terms, side edges and W tables ignored); for g = 2 the row cost is exactly a rank.
  kappa_c >= cut bits at one column is all that is used about the column network.
* Numerics: lam = 2 was NOT computed (too large); N = 3 search inconclusive; all "not found" are evidence only.
* The Alman-Williams-type estimate and the kill-and-compress bound are sketches [E].
