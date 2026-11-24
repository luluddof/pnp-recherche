# G14e — entangled (non-product) expansions: smallest cases.  Status tags: [P] proved here, [N] numerical, [C] conjecture, [O] open

## 1. Precise formulation
Grid g x h, lam-bit string on every diagonal edge (endpoints differ in row and column); T = product of all diagonal
equality tests = tensor product over rectangles of crossings [x=v][y=u].  S = total diagonal bits.
Expansion: c*T = sum_i A_i, each A_i typed row or col.  cost_row(A) = bits of the cheapest network on the g merged rows
(g=2: log2 of the rank across the single row cut); cost_col likewise on the h merged columns.
kappa_r = max over row-typed terms of cost_row/S, kappa_c = max over col-typed terms of cost_col/S (T itself: 1 and 1),
and #terms must be 2^{o(S)}.  With recursion exponents e_g, e_h the level costs max(e_g kappa_r, e_h kappa_c) S.
* Product expansions (one piece per rectangle, pieces (R_i,C_i)): every mixture occurs as a term, so
  kappa_r + kappa_c >= 2 rho >= 1 (rho = hull-on-diagonal value: 11/12 Lemma 2, 1/sqrt2 App. B, floor 1/2).  [M, G9d]
  Even per-rectangle TWO-term decompositions with own costs (0,0) (they exist at lam=1, see 2) give exactly kappa_r+kappa_c = 1 - 1/#rect.
* Harmonic law u(gh) <= u(g)+u(h)  <=>  max(e_g kappa_r, e_h kappa_c) >= e_g e_h/(e_g+e_h) for all (e_g,e_h)
  <=>  kappa_r + kappa_c >= 1.
  E1 (precise): an entangled expansion beats the harmonic law iff  kappa_r + kappa_c < 1  with 2^{o(S)} terms.
* Known floor (restriction to one orthant class) [P]: for a column bipartition H1|H2 (p and h-p columns), 2 x h grid,
  T = A + B:  rk_row(A) + rk_{H1|H2}(B) >= 2^{p(h-p) lam}.  2x3: kappa_r >= 1/3 or kappa_c >= 1/2 (all three single-column
  cuts of B have >= 2 lam - 1 bits, a triangle network needs >= half the sum of its cut bits).  2x2: max >= 1/2.
  So the floor leaves the whole region {kappa_r + kappa_c < 1} minus a corner open.
* g = h = 2 (one crossing, N = 2^lam): terms of one type can be merged (rank is subadditive), so WLOG two terms and
  E1 becomes:  SWAP_N = X + Y,  r = operator-Schmidt rank(X) (= row rank), c = rank(Y) (= column rank).
  [P] both flattenings of the crossing are SWAP_N; realignment fixes SWAP.  Equivalently: a subspace W of N x N matrices
  of codimension c on which transposition M -> M^T equals M -> sum_{k<=r} F_k M G_k^T.
  Known points: (1, N(N-1)/2) and mirror [X = identity: delta_xy delta_uv row-rank 1, SWAP - I = -2 P_antisym];
  eight-square (8,8) at N=8; App. B: the circle kappa_r^2 + kappa_c^2 = 1 (so kappa_r+kappa_c >= 1 there too);
  floor r + c >= N.  A 2x2 point with kappa_r+kappa_c < 1 is not known and not excluded.  [O]
  NOTE: such a point alone would NOT beat the harmonic law (products of it over rectangles still give sum >= 1);
  the entanglement has to be across rectangles with different row/column pairs, i.e. g or h >= 3.

## 2. Data (lam = 1 unless stated)
### 2a. One crossing, two terms: SWAP_N = X + Y, (r,c) = (OSR(X), rank(Y)).  Real field, Levenberg-Marquardt, random starts [N]
(crossing.py, crossing_ap.py, c4_33.py; outputs cross4_summary.out and the task log quoted in the report)
| N | floor r+c>=N | exact points found | not found (best residual^2) |
| 2 | 2 | (1,1)  [P: delta_xv delta_yu = delta_xy delta_uv - eps_xu eps_yv] | - |
| 3 | 3 | (1,3), (3,1) | (1,2),(2,1): 3.0;  (2,2): 2.0  -> min r+c = 4 > floor |
| 4 | 4 | (1,6) [and mirror (6,1)], (4,4) four-square | (2,2): 8.0; (2,3): 5.0; (1,5): 3.0; (4,1): 5.0; (2,4),(4,2): ~3.0; (2,5): 1.02 |
N = 4, (3,3): residual^2 -> 4e-6 while |A|_F = |B|_F -> 2e6, singular values (e^-2, e^-1, e^-1), the divergent leading part is
rank one in BOTH flattenings (checked).  So (3,3) looks BORDER-feasible (limit of decompositions), not exactly feasible. [N, not proved]
A border decomposition is usable (interpolation in e, polynomial overhead).  It is non-product over the two bits of N=4 and
beats the grouped product frontier of the one-bit identity {(1,6),(4,4),(6,1)} and the four-square point (4,4).
Pattern (N-1,N-1) for N=2,4 [C: also N=8, i.e. (7,7) against the eight-square (8,8); untested].
"Not found" = local minima at integer residuals from 8-40 starts; evidence, not proof.

### 2b. 2 x 3 grid, lam = 1 (12 qubits, row rank of T = 64, three column cuts of rank 16)
(grid23.py -> grid23_h3.out; product23.py; grid23_generic.py -> generic23.out)
* Candidate family = all 10395 singlet-matching states (span = SL2 invariants, dim 132; contains every expansion obtained by
  applying the one-bit identity anywhere, in any order - not only rectangle by rectangle).
  r(M) = # pairs across the row cut, c(M) = # pairs between different columns (T: 6 and 6).
  T in span{r<=R} + span{c<=C}  iff  (R,C) >= (0,4) or (2,2) or (4,0) or C=6 or R=6.   [N, exact linear algebra in floats]
  This is EXACTLY the product frontier (terms (2i, 6-2i)): kappa_r + kappa_c = 2/3 = 1 - 1/#rect.  No entangled gain.
* Weighted-l1 LP (cost of a matching = min(2^r,2^c)): a 6-term non-product expansion of cost 24 against 26 for the 8-term
  product; its grouped ranks (row-rank(A)=10; column cut ranks of B = 7,7,7) equal those of the grouped product expansion,
  so the gain is a regrouping of the same two tensors (O(1) effect).
* Grouped product trade-off (row-rank(A); cut ranks of B): (0;16,16,16) (1;16,16,16) (10;7,7,7) (28;1,1,1) (64;0).
  Floor: row-rank(A) + cut_j(B) >= 4.  Observed minimum of that sum: 17.  The floor is far from met at this scale.
* Generic search (A = truncated SVD, B = one triangle network, ALS): nothing found at (r;bonds) = (4;2,2,2), (13;2,2,2),
  (8;2,2,4), (10;3,3,3), (1;3,3,3).  INCONCLUSIVE: the method converged only to 1e-4 on a point known to be feasible.

## 3. Verdict
* Inside ONE crossing (2x2 geometry) non-product structure over bit blocks helps a little at N=4 (border (3,3)); this is the
  same phenomenon as 11/12 -> 1/sqrt2 and cannot beat the harmonic law (section 1, NOTE).
* ACROSS rectangles (2x3, the first case where it could matter): no gain found; in the invariant family the frontier is exactly
  the product one.  All effects seen are O(1) bits; lam=1 cannot separate kappa_r+kappa_c<1 from =1 (terms count as much as ranks).
* Suggestion [C]: for the 2 x h grid, T in V_R + W_C (matching family) iff R + C >= S - 2: a skein-theoretic harmonic law.
  A proof would extend the harmonic law from product pieces to all SL2-equivariant expansions; it says nothing on non-invariant terms.

## 4. Adversarial self-review
* The reduction "E1 <=> kappa_r + kappa_c < 1" inherits the cost model of G9d (balanced terms, straight edges and W tables ignored).
* All infeasibility statements are numerical (local optimisation / float rank decisions).  No lower bound beyond the restriction floor
  was proved; N=3 min r+c = 4 and the r=1 bound c >= N(N-1)/2 are unproved.
* The (3,3) border point at N=4 is a numerical limit (max entry error 4e-4 at norm 1e6); no closed form extracted.
* Field = reals; complex or finite characteristic not examined.  3x3 grid not done.  2x3 generic search is weak.
* For g>2 the cost of a merged-direction term is a network cost, not a cut rank; only cut-rank lower bounds were used.
