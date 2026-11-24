# G22b-bis — continuation of G22b: the monomial (intercalate) model; upper and lower bounds   (working notes, saved as I go)
Tags: [P] proved here (short proof, unreviewed), [E] exact computation, [N] numerical, [H] heuristic, [L] literature (read in the paper).
Inherited from g22b/: shape.py (SAT for S/T/U-palette arrays + EXACT verifier: full 4-tensor identity, flint ranks), convert.py,
anat.py, blocks.py, fr3/4/5.txt (shape frontier incl. U palette), o_6_77 (ST (7,7)@6), o_9_U16, o_10_U16 (integer [9,9,16],[10,10,16]).
Convention: X(x,y,u,v) = [x=v][y=u]; C(W) = sum_l G_l W G'_l^T; R(W) = sum_j tr(H_j W) H_j; X(W) = W^T.

## 1. Anatomy of ALL exact solutions on file  [E, convert.py / blocks.py rerun here]
 files: g19b exact_4_5_3, g20b exact_m_4_5_3, exact_m_4_4_4, g19b exact_5_7_4, exact_5_6_6, exact_6_7_7.
 * (5,3)@4 (x2), (4,4)@4, (7,7)@6: STRICT ST TYPE: G'_l = G_l^T, all G_l, H_j signed partial permutations, supports disjoint
   and covering [N]^2, R = direct sum of p rank-one +-1 blocks. Signed intercalate condition: 0 violations; intercalation graph
   bipartite S|T. Class sizes: (5,3)@4 S 3,2,2 / T 3,2,2,1,1; (4,4)@4 S 3,2,2,1 / T 3,2,2,1; (7,7)@6 S 3,3,3,3,2,2,2 / T same.
 * (7,4)@5 and (6,6)@5 are NOT of the strict type. Precisely:
   - all G_l and G'_l ARE signed partial permutations (monomial), and G'^T_l agrees with G_l on most cells, but not all:
     (7,4)@5: l=1: G has (3,0), G'^T has (2,1) instead; l=3: G has (4,4), G'^T has (2,1) instead (other cells equal).
     (6,6)@5: l=3: G = {(2,4),(4,0)}, G'^T = {(4,0)};  l=5: G = {(3,1)}, G'^T = {(0,2),(3,1)}.
   - R = X - C is block diagonal (bipartite components of its support) with blocks: (7,4)@5: two 4x4 blocks of rank 2 (sparse),
     one 3x3 rank 1, two 2x2 rank 1 (2+2+1+1+1 = 7); (6,6)@5: one 4x4 rank 2 (sparse), one 3x3 rank 1, three 2x2 rank 1 (=6).
     So the deviation = a few "dangling" cells where G and G'^T disagree, which are absorbed by RANK-2 blocks of R (not +-1 rank one).
 * [E, g22b fr5.txt, SAT exhaustive over the S/T/U shape] at N = 5 the shape (even with the U palette) needs (7,7); exact b(5) = 6
   is reached only by the non-strict solution. So the strict model is already lossy at N = 5 (7 vs 6); at N = 4, 6 it is optimal
   among known solutions (b(4) = 4, b(6) = 7).

## 2. Correspondence theorem  [P, proof below; checked [E] on every array on file by shape.verify_shape / convert.py]
Model ("STU", contains the strict symmetric monomial model ST = U-free case).  Data: an N x N array, each cell (a,b) has a colour
k and a sign s(a,b) in {+1,-1}; colours are split into palettes S, T, U; class of k = M_k (cells of colour k), signed matrix
K_k = sum_{(a,b) in M_k} s(a,b) E_ab.  Identity:
   (*)  W^T = sum_{k in S} K_k W K_k + sum_{k in T} tr(K_k W) K_k + 1/2 sum_{k in U} ( K_k W K_k + tr(K_k W) K_k )   for all W.
 i.e. C = sum_S G_l (x) G_l^T + 1/2 sum_U (same), G'_l = G_l^T; R = sum_T H_j[x,y]H_j[v,u] + 1/2 sum_U (same).
THEOREM 1. Assume the S- and U-classes are partial permutations (no two cells in a row/column) [T-classes: arbitrary].
 Then (*) holds  <=>  (i) every class (incl. T) is a partial permutation, and (ii) for every 2x2 submatrix (rows x != v, columns
 u != y) with one diagonal monochromatic, the other diagonal is monochromatic too, the product of the four signs is -1, and the
 palettes of the two diagonal colours are {S,T} or {U,U}.  [(i)+(ii) without palettes = an HWY "intercalate matrix" with
 n = |S|+|T|+|U| colours, HWY sect. 7.1 = integer sum-of-squares formula [N,N,n] (HWY Prop 7.1); the palette condition says: the
 intercalation graph (colours adjacent iff they are the two diagonals of some 2x2 submatrix) has U as a union of components and
 is bipartite S|T on the rest.]  Moreover then rr(R) = |T|+|U| =: p and rc(C) = |S|+|U| =: q EXACTLY.
 In particular (U empty): monomial decompositions with G'_l = G_l^T, G_l signed partial permutations, R = sum_j tr(H_j W) H_j with
 +-1 matrices H_j, all supports disjoint and covering [N]^2   <=>   signed intercalate matrices [N,N,p+q] with a 2-colouring of the
 colours (q G-colours, p H-colours) such that every intercalate is bichromatic.  (Weights: if one allows coefficients +-1 in front
 of the terms, evaluating (*) on W = E_ba forces every coefficient to be +1 on the class of (a,b); so nothing is lost.)
PROOF. Write w_C(S)=1, w_C(T)=0, w_C(U)=1/2, w_R = 1 - w_C. Coefficient of (*) at (x,y,u,v):
   [x=v][y=u] = sum_k w_C(k) [ (x,u),(v,y) in M_k ] s(x,u)s(v,y) + sum_k w_R(k) [ (x,y),(v,u) in M_k ] s(x,y)s(v,u).
 x=v,y=u: both pairs are the single cell (x,y): RHS = w_C + w_R = 1 = LHS.  x=v, y!=u: LHS 0; the C-pair (x,u),(x,y) lies in one
 row, impossible for S/U (partial perms); so w_R(k)[(x,y),(x,u) in M_k] = 0: no T/U class has two cells in a row; same for columns
 (y=u, x!=v). Hence (i). x!=v, y!=u: the 2x2 submatrix rows {x,v}, cols {u,y} has diagonals D = {(x,u),(v,y)} and A = {(x,y),(v,u)};
 the tuple (x,y,u,v) gives  w_C(D) s_D + w_R(A) s_A = 0  and the tuple (x,u,y,v) gives  w_C(A) s_A + w_R(D) s_D = 0, where w(.) = 0
 if the pair is not monochromatic and s_D = product of the two signs. D mono S: first eq. forces A mono with w_R = 1 (T) and
 s_A = -s_D. D mono T: second eq. forces A mono with w_C = 1 (S), s_A = -s_D. D mono U: w_R(A) = 1/2, A mono U, s_A = -s_D.
 Conversely (i),(ii) make every coefficient vanish. Ranks: C = sum_k w_C(k) vec(K_k) vec(K_k^T)^T in the (x,u)|(y,v) flattening;
 the vec(K_k) (k in S u U) have disjoint nonempty supports, so they are independent, as are the vec(K_k^T): rank = |S|+|U|. Same
 for R with T u U. QED
 [E] verify_shape re-checks the identity by full expansion over [N]^4 and the flint ranks: all arrays below satisfy it.

## 3. Lower bound in the model (the trade-off)   THEOREM 3 [P; part (b) uses HWY Lemma 7.4 = [L], proof read, Lemma 7.7 not re-checked]
Notation: sigma = |S|, tau = |T|, mu = |U|; p = tau + mu, q = sigma + mu (Thm 1). g*(t) := max{g <= N : g(g-1)/2 <= t} <= sqrt(2t)+1.
 (a) KEY LEMMA. If a T-colour j has h cells (a_1,b_1),...,(a_h,b_h), then in the h x h submatrix rows {a_i}, cols {b_i}:
     every off-diagonal cell (a_i,b_k) is an S-colour, and each S-colour occupies at most the symmetric pair {(a_i,b_k),(a_k,b_i)}.
     Hence sigma >= h(h-1)/2.  Symmetrically every S-colour of size g forces tau >= g(g-1)/2.
     Proof: rows a_i,a_k, cols b_i,b_k: the diagonal (a_i,b_i),(a_k,b_k) is mono T, so (Thm 1 (ii)) the antidiagonal
     (a_i,b_k),(a_k,b_i) is mono S. If one S-colour l sits at two off-diagonal cells (a_i,b_k),(a_i',b_k') (then i != i', k != k'),
     the 2x2 submatrix rows a_i,a_i', cols b_k,b_k' has an S-mono diagonal, so its antidiagonal (a_i,b_k'),(a_i',b_k) is T-coloured;
     inside the submatrix T occurs only on the diagonal, so k' = i, k = i'. The h(h-1) off-diagonal cells thus need >= h(h-1)/2
     distinct S-colours. QED   (This is HWY's "full submatrix" step; bichromaticity upgrades their k^{3/2} colours to k^2/2.)
 (b) A U-colour of size s: same argument => the s x s submatrix is (after permutation / sign normalisation as in HWY Thm 7.5) a
     FULL HWY intercalate matrix whose off-diagonal colours are all U-colours => mu >= s - 1 trivially, and 2 mu^2 >= C(s,3)
     (HWY proof of Lemma 7.4: each ordered colour pair at <= 2 of the C(s,3) "positions"), i.e. s <= (12 mu^2)^{1/3} + 2.
 (c) Counting cells:   N^2 <= sigma g*(tau) + tau g*(sigma) + mu s*(mu).
COROLLARIES.
 ST model (mu = 0):  N^2 <= q g*(p) + p g*(q) <= sqrt2 (q sqrt p + p sqrt q) + p + q.
   * balanced: p = q >= (1 - o(1)) N^{4/3} / 2   (exponent 4/3);   table: N=8 -> 8, 16 -> 21, 32 -> 52, 64 -> 128, 128 -> 325.
   * trade-off, q <= p:  p >= N^2 / (2 (sqrt(2q) + 1)),  i.e. p^2 q >~ N^4/8. Strictly above the proved line p q >= (2/9) N^2 for
     every q >= 2; e.g. q = N^beta forces p >= N^{2 - beta/2}/(2 sqrt2) (1-o(1)), p q >= N^{2+beta/2}/(2 sqrt 2).
   * q = 1: p >= N(N-1)/2, EXACT for every N (construction W^T = W + sum_{a<b} tr(H_ab W) H_ab, H_ab = E_ab - E_ba).
 STU model: mu <= min(p,q); mu s*(mu) <= 2.3 mu^{5/3} + 2 mu. So for q <= p:
   either q >= mu >= c N^{6/5}, or p >= N^2/(4(sqrt(2q)+1)).  Balanced: p >= Omega(N^{6/5}) (this is just HWY: p+q >= number of
   colours >= sigma_Z(N)); trade-off p >= N^2/(4(sqrt(2q)+1)) for all q <= N^{6/5}/5 (N large) -- again above the line.
 [E] checks: bound vs exhaustive SAT frontier (g22b fr3/fr4/fr5, ST = STU frontier there): never violated; tight at q = 1 (N=3,4,5:
   p = 3, 6, 10), tight BALANCED at N = 2 (1,1) and N = 8 (8,8) (slack 0: 16 colours all of size 4 = g*(8)); N=4: bound 3 vs 4,
   N=5: 5 vs 7, N=6: 6 vs 7 (and (6,6)@6 ST is UNSAT [E, 27 s]).  Lemma (a) checked on the ST arrays (N=6 (7,7), N=8 (8,8)).

## 4. Upper bounds in the model
 * STU: ANY integer SOS formula [N,N,m] (= HWY intercalate matrix) read with all colours in U gives (p,q) = (m,m) [P, Thm 1];
   so b_STU(N) <= sigma_Z(N) = O(N^2/log N) [L: HWY p.3, Hurwitz-Radon].  hr.py (Hurwitz-Radon block construction: rho(16) = 9
   mutually anticommuting skew real Pauli strings found by max-clique, rows split into blocks of <= rho rows, columns into blocks
   of 16, fresh colours per block) [E, exact verify incl. ranks]:  N = 16: [16,16,32] -> (32,32);  N = 32: [32,32,128] -> (128,128).
   Per component of the intercalation graph one may choose S/T (if bipartite) instead of U: for these HR arrays EVERY component is
   non-bipartite (so nothing gained).  [P] reason: for a Pauli/Z_2^k-type block with row shifts v_1..v_r the intercalation graph is
   the Cayley graph of F_2^k with generators v_i + v_i'; any 3 rows give a triangle. So Clifford/HR structure is useless for ST.
 * ST: tiling by an s x s ST array with fresh colours per tile is valid (intercalates stay inside tiles), (p,q) scales by (N/s)^2.
   Best tile found: N = 8, (8,8) [E, SAT 0.9 s, shape_8_8_8_S8T8U0.json] -> b_ST(N) <= N^2/8 for 8 | N.  So best known ST exponent = 2.
   ST balanced exact values: N=2:1, 3:3, 4:4, 5:7 (b=6), 6:7 [(6,6) UNSAT], 7: <= 8 [(7,7) undecided in 400 s], 8: 8 (= lower bound).
   These equal b(N) for N = 2,3,4,6,8 and b(7) <= 8 is matched: the ST model is exact except N = 5 so far.
 * STRUCTURE OF THE N = 8 TIGHT ARRAY [E, grpans.py, grp_8_Q8.json]: it can be taken invariant under the QUATERNION group:
   rows = columns = Q8, cell (a, a c) gets colour (B(c), a f(c)) with a block map B: Q8 -> {S,T} (4 + 4) and some f: Q8 -> Q8;
   exactly verified. The same ansatz is UNSAT for Z8, Z4xZ2, D4 and F_2^3 (transl.py) [E]. (Abelian 2-groups fail for the
   triangle reason; for linear f over F_2^n I proved [P] that blocks = cosets fail likewise.)
   The group ansatz gives p = N (#T-blocks), q = N (#S-blocks), class size g = block size, and is COMPATIBLE with the 4/3 law:
   N = 64, 2+2 blocks of 16 differences would give (128,128) = N^{4/3}/2 (the lower bound). If such families exist (groups of
   order g^{3/2}), the ST model would reach balanced exponent 4/3 < sqrt2 -- i.e. beat Appendix B with monomial pieces. NOT tested
   (the direct SAT encoding at N = 64 has ~10^9 clauses); obvious product Q8 x Q8 fails (intercalation graph has triangles) [P].
 * N = 12 (18,18) and N = 16 (28,28) generic ST SAT (with the class-size caps of Thm 3): see section 6 (status).

## 5. Verdict
 * Correspondence (Thm 1) [P]: symmetric monomial decompositions (G' = G^T, signed partial permutations, R = sum tr(H_j W) H_j,
   disjoint supports) == signed intercalate matrices [N,N,p+q] whose intercalation graph is bipartite G|H; with the U palette
   (weight 1/2 terms) == arbitrary integer SOS formulas, U = non-bipartite part. Ranks are exactly the colour counts.
 * Model lower bounds [P]: ST: balanced >= ~N^{4/3}/2, trade-off p >= N^2/(2(sqrt(2q)+1)) (q <= p), exact at q = 1, tight at N = 2, 8.
   STU: balanced Omega(N^{6/5}) (HWY), trade-off p >= N^2/(4(sqrt(2q)+1)) for q <= N^{6/5}/5. Both strictly above the general
   line p q >= (2/9)N^2 (for q >= 2): any decomposition near the line (e.g. q = N^beta, p = N^{2-beta+o(1)}, beta > 0) must use
   pieces outside the model (non-symmetric G' != G^T, non-monomial factors, or R not a sum of +-1 rank-one blocks).
 * Model upper bounds: STU: O(N^2/log N) (HR/integer SOS; exact (32,32)@16, (128,128)@32). ST: N^2/8 (tiling by the Q8 array).
   Best known exponent of the model: 2 (both). Lower exponent 4/3 (ST) / 6/5 (STU). Gap open.
 * The 4/3 < sqrt2 lower bound does NOT exclude that monomial (ST) constructions beat Appendix B; the N = 8 Q8 array meets the
   4/3 counting bound exactly, and a group ansatz of the right size (order g^{3/2}, 2+2 blocks) would give exponent 4/3.
   Conversely: beating 4/3 (ST) or 6/5 (STU) requires non-monomial/non-symmetric pieces, and already at N = 5 the exact optimum 6
   needs G' != G^T (the model gives 7).
 * Lower-bound side: the model bound is NOT a bound on b(N) in general (Appendix B / Hrubes pieces are not monomial). It says
   where a lower-bound proof for the general problem cannot come from: counting arguments alone cannot beat what monomial
   structures allow, and the general problem's line p q >= (2/9) N^2 is far from the model's p^2 q >~ N^4/8.

## 6. Run status / files
 * Generic ST SAT with class-size caps: N=12 (18,18) and N=16 (28,28): UNDECIDED (killed at 600 s, no output). N=7 (7,7): undecided
   at 400 s; (8,8)@7 SAT. So whether ST beats tiling (p/N^2 < 1/8) at N = 12, 16 is open. All processes killed.
 * Files: notes.md; shape.py (S/T/U SAT + exact verifier, from g22b); run1.py, run2.py (capped ST SAT); stbal.sh; bound.py (Thm 3
   tables + lemma check); hr.py (Hurwitz-Radon arrays N=16,32); transl.py (F_2^n translation ansatz + F_2 sign solver);
   grpans.py (group ansatz, groups of order 8); q8show.py; arrays: shape_8_8_8_S8T8U0.json, grp_8_Q8.json (Q8-invariant (8,8)),
   shape_6_7_7_S7T7U0.json, shape_7_8_8_S8T8U0.json, shape_9/10_16_16 (U-only integer formulas).
 * Caveat on the last verdict bullet of sect. 5 (\"where a lower bound cannot come from\"): [H], not a theorem.
