# G19b — exact small cases of the crossing trade-off, and the G18g <-> G18h dictionary  (working notes, saved as I go)

Tags: [P] proved here (short proof written, unreviewed), [E] exact computation, [N] numerical only, [H] heuristic, [O] open.
Nothing here is a new algorithm or a lower bound for 3-SAT.

## 0. Notation (as G18g)
X(x,y,u,v) = [x=v][y=u], x,y,u,v in [N].  rr = rank in the flattening (x,y)|(u,v), rc = rank in (x,u)|(y,v).
(p,q) FEASIBLE: X = R + C, rr(R) <= p, rc(C) <= q.  Row-wise X is the swap F, and C = sum_{l<=q} G_l (x) G'_l, so
   (p,q) feasible  <=>  min over C of operator-Schmidt rank <= q of rank(F - C) <= p.
Symmetric under (p,q) <-> (q,p) (exchange y and u).

## 1. Dictionary (task 1)   [P]
G18h language: a profile L(E) = log dim W_{<=E}, budget M, piece at weight x charged (row, col) = (L(x), L(M-x)).
G18g language: thresholding the pieces at weight x gives X = R + C with p <= |I| 2^{L(x)}, q <= |I| 2^{L(M-x-step)}.
So a profile is the parametrised curve  x -> (theta, T) = (L(x)/lambda, L(M-x)/lambda)  in G18g's diagram, and
(V) / the shape of L only fix the PARAMETRISATION by weight, which the (theta,T) diagram does not see.
 * App. B: L = c sqrt(E), M = 4 h1, c sqrt(h1) = lambda: a^2 + b^2 = 4, the circle.
 * rho = 1/2 in the accounting <=> L(x) + L(M-x) = 2 lambda for all x <=> the curve is the line a + b = 2.
 * The staircase of G18h as literally defined (L = 0 below h1, lambda on [h1,3h1], 2 lambda above) has only three
   piece types (1,N^2), (N,N), (N^2,1).  Thresholding between the first two types gives a feasible pair
   (p,q) with p <= |I|, q <= |I| N, i.e. the point (theta,T) = (0,1), with theta + T = 1 < 2: BELOW the line.
   Theorem 2 (p q >= 3N^2/16) EXCLUDES it unless |I| >= 0.43 sqrt(N).  So the hard staircase does not exist
   with N^{o(1)} pieces.  [P]
 * Continuity lemma [P]: if X = sum of |I| pieces with rr_i rc_i <= N^2 (rho = 1/2 exactly), sorted by rr, then
   consecutive row ranks satisfy rr_{i+1}/rr_i <= 16 |I|^2 / 3 (threshold between i and i+1: p <= |I| rr_i,
   q <= |I| N^2/rr_{i+1}, Theorem 2).  So the row ranks fill the whole scale 1..N^2 with gaps N^{o(1)} and
   |I| >= 2 lambda / log2(16|I|^2/3): a rho = 1/2 family is necessarily a "continuous" one.
 * Repaired staircase compatible with (V): L rises continuously from 0 to lambda on [0,h1], is FLAT on [h1,3h1],
   rises to 2 lambda on [3h1,4h1], with L(x) + L(4h1-x) = 2 lambda.  Its image in the diagram is the WHOLE line
   (the plateau is the single point (1,1) = Naslund's point P(N) = N^{1+o(1)}).  So in G18g's language the two
   scenarios of G18h ("linear L without (V)" and "staircase with (V)") are the same statement:
        (N^theta, N^{2-theta+o(1)}) feasible for every theta in [0,2].
   Object that must exist: for every theta a subspace U_theta of M_N of codimension N^{2-theta+o(1)} on which
   transposition is a sum of N^{theta+o(1)} two-sided multiplications (G18g 1(f)); (V) would in addition ask that
   they come from one graded space with no new states between h1 and 3h1.

## 2. Theorem 3 (q = 1 exactly)   [P]
For all G, G' in M_N(K), any field:  rank(F - G (x) G') >= N(N-1)/2, with equality for G = G' = I.
Hence the least p at q = 1 is exactly N(N-1)/2 (Theorem 2 gives only N^2/4), and by symmetry the least q at p = 1.
Proof.  rank(F - G(x)G') = rank(1 - Theta), Theta = (G(x)G')F.  Theta^2 = GG' (x) G'G, and
tr Theta^{2k+1} = tr (GG')^{2k+1}, tr Theta^{2k} = (tr (GG')^k)^2.  With mu_1..mu_N the eigenvalues of GG', the power
sums of Theta are those of the multiset {mu_i} u {+-sqrt(mu_i mu_j) : i<j}, so
     charpoly(Theta)(t) = prod_i (t - mu_i) prod_{i<j} (t^2 - mu_i mu_j)
(char 0 by Newton; an identity between integer polynomials in the entries of G, G', hence valid over every field).
The root t = 1 has multiplicity <= N + N(N-1)/2, so dim ker(1 - Theta) <= N(N+1)/2.  QED
(to be checked numerically/exactly below)
Caveat found on review: in characteristic 2 the factor t^2 - mu_i mu_j has the double root 1, so the proof of Theorem 3
needs char != 2.  (Over GF(2), N = 4, the value 6 is G18h's exhaustive computation.)  Exact check of the charpoly identity:
thm3_check.py (N = 2,3,4, random integer G, G': identity and rank bound hold).  [E]

## 3. Theorem 4 (Segre argument; q = 2)   [P, char != 2; unreviewed]
If N >= 3 and (p,2) is feasible then p >= N.  In particular (2,2) is infeasible for N = 3, so P(3) = 4
(frontier of N = 3: (0,9),(1,3),(3,1),(9,0)); Theorem 2 alone gives P(3) >= 3.
Proof (over the algebraic closure).  R = F - C, C = G_1(x)G'_1 + G_2(x)G'_2, K = ker R = {W : W^T = sum_l G_l W G'_l^T},
codim K <= p.  S = {(a,b) in P^{N-1} x P^{N-1} : a b^T in K} = Segre variety cut by a linear space of codimension <= p:
every component of S has dimension >= 2N-2-p.  For (a,b) in S:  b a^T = A(a) B(b)^T with A(a) = [G_1 a, G_2 a],
B(b) = [G'_1 b, G'_2 b] (N x 2).  Sylvester: 1 >= rank A(a) + rank B(b) - 2, so S = S_1 u S_2, S_1 = {rank A(a) = 1},
S_2 = {rank B(b) = 1}, both closed.  A component Z lies in one of them, say S_1 (the other case is the mirror image under
(x,u) <-> (y,v)).  On S_1, b spans Im A(a): Z -> P^{N-1}, (a,b) -> a is injective.  If p <= N-1 then dim Z >= N-1, so
rank [G_1 a, G_2 a] <= 1 for EVERY a.  A linear space of matrices of rank <= 1 has a common image or a common kernel:
either G_1, G_2 are proportional (then rc(C) <= 1 and p >= N(N-1)/2 >= N by Theorem 3) or G_l = c d_l^T (then
rank C <= N and p >= N^2 - N).  Both contradict p <= N-1.  QED

## 4. Small cases: computations so far (saved while running)
Tools: vp.py (variable projection: minimise over C the tail sum_{i>p} sigma_i(F-C)^2, L-BFGS, restarts) [N];
sparsify.py (gauge-sparsify a numerical solution with the symmetry group of X, then round); verify.py (exact: python-flint
ranks over Q + explicit rank factorisation of R re-expanded, R + C == X entry by entry) [E]; diag.py (exact frontier of the
"diagonal Pauli family", N = 2, 4) [E]; lmfast.py (LM on all factors, slow, superseded by vp.py).
N = 3 [N, 25-60 restarts each]: (2,1) tail 3.000, (2,2) tail 2.000 (agrees with Theorems 3, 4: infeasible).
N = 4 [N]: (5,1) tail 3.0 [infeasible, Thm 3]; (4,2),(2,4) tail 3.006; (5,2),(2,5) tail 1.05-1.07 (inf looks like 1);
   (5,3),(3,5) FEASIBLE; (3,3),(4,3),(3,4): tail 1e-3..3e-3 with growing coefficients = BORDER only (as G14e saw for (3,3)).
   NEW EXACT POINT (5,3) for N = 4 [E, verify.py]: entries in {0,+-1}, E_ij = matrix unit (1-based):
     G_1 = -E33 + E42,  G_2 = -E14 - E23 + E41,  G_3 = E22 - E31;   G'_1 = E24 - E33,  G'_2 = E14 - E32 - E41,  G'_3 = -E13 + E22;
     C = sum_l G_l(x,u) G'_l(y,v) has rc = 3 and R = X - C has rr = 5 exactly (checked by expansion).
   Diagonal Pauli family (diag.py): least p = 6,6,6,4,4 for q = 1..5: the points (5,3) and border (3,3) are NOT of this type.
N = 5 [N, running]: (7,5),(7,6),(8,4) FEASIBLE; (6,5) tail 2e-4 (border?); (5,5) tail 0.029 (not found); (9,2) tail 1.07.
   N = 5 further [N]: (6,6), (7,7), (7,4) FEASIBLE; (8,3) tail 3e-7 with coefficients 63 (border?).
   NEW EXACT POINTS for N = 5 [E, torus.py + verify.py; entries in {0,+-1}; files exact_5_6_6.json, exact_5_7_4.json]:
   (6,6):  G  = [+E24+E33+E41, -E52, +E11-E25-E32, +E35-E51, +E14+E45+E53, -E42]
           G' = [+E14+E33+E42, -E25, +E11-E23-E52, -E15, +E35+E41+E54, -E24-E31]
   (7,4):  G  = [+E25+E44+E51, +E23-E31+E42, +E12-E25-E51, -E14+E35+E53]
           G' = [-E33+E44+E52, -E13+E24+E32, -E15+E21-E33, +E35-E41+E53]
   (C = sum_l G_l(x,u) G'_l(y,v); R = X - C; rr(R), rc(C) computed exactly over Q; R + C = X re-expanded entry by entry.)

## 5. Theorem 5 (Theorem 2 + Theorem 3 combined)   [P, char != 2; unreviewed]
If (p,q) is feasible, q >= 1, then for every integer k >= 1 with m := N - k(q-1) >= 1 and k' := min(k,m):
        p >= k m - k'(k'+1)/2          (and the mirror statement with p and q exchanged).
Proof.  As in Theorem 2 take independent a_1..a_k and B = {b : a_i^T G'_l b = 0 for all i and all l >= 2}, dim B >= m,
b_1..b_m independent in B.  Contracting X = R + C with a_i(y) b_j(v) leaves ONE column term:
  b_j(x) a_i(u) - g_1(x,u) gamma(i,j)  has rank <= p as a matrix (i,x)|(j,u),  gamma = (a_i^T G'_1 b_j).
Compose x with a left inverse of (b_j) and u with a left inverse of (a_i): the rectangular swap S: K^m (x) K^k -> K^k (x) K^m
minus gamma (x) g~ (gamma: K^m -> K^k, g~: K^k -> K^m) has rank <= p.  For Theta = S^{-1}(gamma (x) g~): Theta^2 = g~gamma (x) gamma g~,
tr Theta^{2n+1} = tr (g~ gamma)^{2n+1}, tr Theta^{2n} = (tr (g~gamma)^n)^2, so
charpoly(Theta) = t^{km-k'^2} prod_{i<=k'} (t - mu_i) prod_{i<j<=k'} (t^2 - mu_i mu_j) and the eigenvalue 1 has multiplicity
<= k'(k'+1)/2.  QED
Size: optimum k ~ N/(2q-1):  p >= (1 - o(1)) N^2/(4q-2)  instead of N^2/(4q):  p q >= (N^2/4) q/(q - 1/2).
It is exact at q = 1 (N(N-1)/2) and only a constant factor in general: the exponent line T(theta) >= 2 - theta is NOT moved.
At small N it is weak: N = 4: q = 2 -> p >= 2, q = 3 -> p >= 1 (Theorem 4 gives q = 2 -> p >= 4).

## 6. Do the small-case gains tensor?  No.   [E for the example, P in general by G18g section 3 / G18f section 2]
kron53.py (G18g's exact.py machinery): Kronecker square of the exact (5,3) identity, N = 16:
  pieces (rr,rc) = (25,256), (131,87), (196,9): max rr*rc = 11397, "rho" = 0.842 (base: 0.790);
  thresholds: (25,83) and (131,9): p q / N^2 = 8.1 and 4.6, against 0.94 for the base point (5,3).
So p q / N^2 < 1 at N = 4 becomes > 4 at N = 16: the mixed terms R (x) C destroy the gain, exactly as G18g / G18f predicted.

## 7. Table of feasible pairs (task 2).  Over Q (proofs: any field of char != 2).  "least q" by symmetry = same table with p,q exchanged.
Columns: q | least p with an EXACT verified decomposition | proved lower bound on p | numerics [N] / status
N = 2:  q=0: 4 | 4 |  exact.   q>=1: 1 | 1 (Thm 2) | exact.        Frontier {(4,0),(1,1),(0,4)}: SETTLED.  P(2) = 2.
N = 3:  q=0: 9 | 9.   q=1: 3 | 3 (Thm 3).   q=2: 3 | 3 (Thm 4).   q>=3: 1 | 1.
        Frontier {(9,0),(3,1),(1,3),(0,9)}: SETTLED [P].  (2,2), (2,1), (1,2) infeasible.  P(3) = 4  (Theorem 2 gives 3).
N = 4:  q=1: 6 | 6 (Thm 3) | settled.
        q=2: 6 | 4 (Thm 4) | p = 4: tail -> 3.0; p = 5: tail 1.05-1.07, inf looks like 1: conjecture 6.
        q=3: 5 (NEW, exact) | 3 (mirrors of Thm 3, 4) | p = 3, 4: BORDER only (tail 1e-3, ~|G|^-0.57, as G14e); exact (4,3),(3,3) not found.
        q=4: 4 (Hurwitz) | 2 | p = 3: border only; p = 2: tail 3.0 (no).
        q=5: 3 (mirror of (5,3)) | 2 | p = 2: tail 1.05 (no).
        q=6..15: 1 | 1 | settled (q <= 5 impossible for p = 1 by Thm 3).
        Exact frontier (proved points, completeness numerical): {(16,0),(6,1),(5,3),(4,4),(3,5),(1,6),(0,16)}; border: + (3,3).
        P(4): proved >= 6; = 7 exactly if (3,3),(2,4),(4,2) are infeasible [N]; border value 6.
N = 5:  q=1: 10 | 10 (Thm 3) | settled.
        q=2: 10 | 5 (Thm 4) | p = 9: tail 1.07 (no): conjecture 10.
        q=3: 10 | 3 | p = 9, 8: tails 3e-7 with coefficients 60-80 (border); p = 7: tail 2e-3 (border candidate, unresolved).
        q=4: 7 (NEW, exact) | 3 | p = 6: tail 5e-3 (border candidate, unresolved).
        q=5: 7 | 2 | p = 6: tail 2e-5 (border); p = 5: tail 0.03 (unresolved).
        q=6: 6 (NEW, exact) | 2.      q=7: 4 (mirror of (7,4)) | 2.      q=10: 1 | 1.
        Exact points known: (10,1),(7,4),(6,6),(4,7),(1,10).  P(5): proved >= 7, exact <= 11, [(7,4): 11, (6,6): 12].
Balanced point b(N) = least p with (p,p) feasible:  N=2: 1;  N=3: 3;  N=4: 4 (border 3);  N=5: 6 exact known (5 not excluded);
  b/N = 0.5, 1, 1, <= 1.2.   Minimum of p q / N^2 over the middle of the frontier: 0.25, 0.33, 0.94 (border 0.56), 1.12 (border <= 0.96).
Position: for N <= 5 everything is within a constant factor of p q = N^2, i.e. "on the line" T = 2 - theta as far as such small N
can say, but strictly above Theorem 2's constant 3/16 for N >= 3 (N = 3: (2,2),(2,1) allowed by Thm 2, infeasible; N = 4: (2,2),(3,2),
(2,3) allowed, infeasible).  The circle is an asymptotic statement (App. B at these N is far worse than the points above), so
"line / circle / between" cannot be decided for N <= 5; the exact frontier is a constant factor (about 5-6 in p q) above Theorem 2.

## 8. Conjectures suggested by the data   [N, O]
 C1. least p at q = 2 is N(N-1)/2, the same as at q = 1 (proved N = 3; N = 4, 5 numerical).  Theorem 4 gives only N.
 C2. border-feasible but not exactly feasible points exist: N = 4 (3,3); N = 5 (8,3), (6,5).  For the algorithm border = exact
     (interpolation, G18g 1(e)), so the relevant frontier is the BORDER one, and the Theorems 3, 4, 5 are proved for exact
     decompositions only (the eigenvalue and Segre arguments do not pass to limits as written); Theorem 2 does (G18f).
N = 6 [N, vp.py]: (7,7), (8,7), (9,6) FEASIBLE (tails 1e-23..1e-26, coefficients <= 2.2); (10,5) tail 2e-5 (border?).
   So b(6) <= 7 < 8 (the eight-square identity restricted to N = 6 gives (8,8)).  Exactification of (7,7): see below if it finished.

## 9. Verdict (task 4)
 * Dictionary: G18h's hard staircase is the point (theta,T) = (0,1), below the line; Theorem 2 kills it (section 1).  Its
   repaired (ramp-plateau-ramp) form and G18h's "linear profile without (V)" are the same object in G18g's diagram: the whole
   line.  (V) and the shape of L do not appear in the diagram at all; they only say how the line is parametrised by weight.
 * Exact small cases: N = 2, 3 settled by proof (P(2) = 2, P(3) = 4).  N = 4: frontier (6,1),(5,3),(4,4),(3,5),(1,6) with exact
   decompositions; lower bounds proved only for q = 1 (and p >= 4 at q = 2, p >= 3 at q = 3); the rest numerical; a border point (3,3).
   N = 5: exact (10,1),(7,4),(6,6); N = 6: (7,7) numerically.
 * Theorem 2 is not tight at small N (first failure: N = 3).  Improvements proved: Theorem 3 (q = 1: N(N-1)/2, exact, a factor 2),
   Theorem 4 (q = 2: p >= N), Theorem 5 (p q >= (N^2/4) q/(q-1/2)).  ALL are constant-factor or low-order: none moves the
   exponent, none is an unbalanced bound p q >= N^{2+eps}.  rho > 1/2 is NOT proved.
 * Nothing points below the circle asymptotically.  The small-case points with p q < N^2 ((5,3), border (3,3) at N = 4) do not
   tensor (section 6, exact example), and the balanced values b(N) = 1,3,4,<=6,<=7 for N = 2..6 are of the form N + O(1): at these
   sizes one cannot tell N^{1+o(1)} from N^{sqrt 2}.  Equally, nothing here supports rho > 1/2.
 * Most informative next computation: b(N) for N = 9..12 with vp.py (81- to 144-dimensional SVDs, minutes per restart), i.e. the
   first sizes beyond the Hurwitz range 1,2,4,8 where a multiplicative gap N^{1+eps} could show: is (N+1,N+1) or (N+2,N+2)
   still feasible at N = 9, 10, 12 (compare [12,12,18] of Zhang-Zhu over Q(i))?  Second: a border-aware lower bound (Theorems
   3-5 for limits), because (3,3) at N = 4 shows the border frontier is strictly below the exact one.

## 10. Adversarial self-review
 * Theorem 3: proof checked by an exact charpoly identity for N <= 4 (thm3_check.py).  Needs char != 2.  Exact decompositions only.
 * Theorem 4: uses the projective dimension theorem over the algebraic closure, Sylvester's inequality, and the classical fact
   that a linear space of rank <= 1 matrices has a common image or a common kernel.  Written in full in section 3; not
   independently reviewed.  It is consistent with all numerics (N = 3: tails 2.0, 3.0; N = 4: (3,2),(2,3) tail 5.0).
 * Theorem 5: the reduction to the rectangular one-term problem is Theorem 2's contraction; the rectangular charpoly formula was
   derived by hand (traces) and NOT separately machine-checked for k != m.
 * None of Theorems 3-5 is proved for border decompositions; for the algorithm the border frontier is what counts.
 * "Not found" / "tail -> 1 or 3" are numerical (20-60 restarts of a local method); tails between 1e-7 and 3e-2 with growing
   coefficients are read as border or unresolved, not as infeasible.  Completeness of the N = 4, 5 frontiers is NOT proved.
 * Every decomposition reported as exact was verified by verify.py: exact ranks over Q (python-flint) and an explicit rank
   factorisation of R re-expanded so that R + C = X entry by entry.  (5,3)@N=4, (7,4)@N=5, (6,6)@N=5.  (7,7)@N=6 is numerical
   unless section 11 says otherwise.
 * The dictionary uses the accounting "rank <= dimension of the weight window" and |I| = N^{o(1)}; the exclusion of the hard
   staircase needs |I| < 0.43 sqrt(N).  The repaired staircase is not constructed; nothing here says it exists.
 * Novelty of Theorems 3-5 not checked against the literature (the q = 1 statement is probably folklore; G14e listed it as unproved).
 * No finite-field (GF(2), GF(3)) searches were run here: the Q-proofs for N <= 3 made them unnecessary and N = 4 is out of
   reach exhaustively.  G18h's GF(2) statements are unchanged.

## 11. Files
 notes.md; thm3_check.py; vp.py (+ out_vp*.txt); lmfast.py (+ out_n4_*.txt); diag.py (+ out_diag4.txt); sparsify.py, torus.py,
 verify.py (+ exact_5_6_6.json, exact_5_7_4.json, sparse_*.npy, vp_*.npy); analyse53*.py; kron53.py.

## 12. Late additions (after sections 7-10 were written)
 * N = 6, (7,7): EXACT [E, sparsify.py + torus.py with threshold 0.3 + verify.py; file exact_6_7_7.json], entries in {0,+-1}:
     G  = [-E36-E41+E55, -E32+E43, +E13+E52, -E22-E63, -E16+E44-E65, +E26+E54-E61, +E11+E25+E34]
     G' = [-E14+E55-E63, -E23+E34, +E25+E31, -E22-E36, +E44-E56-E61, -E16+E45+E62, +E11+E43+E52]
   rr(R) = 7, rc(C) = 7, R + C = X re-expanded entry by entry.  So b(6) <= 7 and P(6) <= 14 (the trivial point (1,15) gives 16,
   the eight-square identity restricted to N = 6 gives (8,8)).
 * N = 6 [N]: (6,6) tail 2e-5 and (7,6) tail 3e-5 with coefficients 25-35 after 80 restarts: border candidates, not exact.
 * N = 7 [N]: (8,8) FEASIBLE (expected: eight-square restricted).
 * Balanced values, exact decompositions in hand:  b(N) <= 1, 3, 4, 6, 7, 8, 8 for N = 2..8  (= for N = 2, 3 [P]; N = 4: >= 3 [P]).
   Border candidates: 3 (N = 4), 6 (N = 6).  All of the form N + O(1).
 * All background computations were stopped before the report.
 * (5,3) at N = 4 also saved as exact_4_5_3.json (torus.py 4,5,3, verified exactly).
