# G20a — D = 3 label directions: feasible triples of the 3-box, beyond the smallest case  (work log, saved as I go)

Tags: [P] proved here (proof written), [E] exact computation, [N] numerical (calibrated solver), [H] heuristic, [O] open.
Nothing here is a lower bound for 3-SAT or a new algorithm.

## 1. Object and cost (task 1)

X_3 (G18g 1bis) = tensor on the 8 corners c in {0,1}^3 of a 3-box, X_3 = prod over the 4 body diagonals {c, cbar} of [z_c = z_cbar],
z_c in [N].  Equivalently: 4 strands delta = 0..3 (strand = body diagonal), legs a_delta, b_delta; X_3 = tensor product of 4
identities I_N(a_delta, b_delta).  Orientation of strand delta = (d1,d2) in cut j: cut 0: 0, cut 1: d1, cut 2: d2, so
   cut 0 = (a0 a1 a2 a3 | b0 b1 b2 b3),  cut 1 = (a0 a1 b2 b3 | b0 b1 a2 a3),  cut 2 = (a0 b1 a2 b3 | b0 a1 b2 a3).
(The 3 cut indices = the 3 box coordinates; the cube symmetries permute them, so everything is symmetric in (p1,p2,p3).)
In operator language: cut 0 = X_3 as the identity on (K^N)^{(x)4}; cut 1 = partial transpose on strands {2,3}; cut 2 = partial
transpose on strands {1,3}.  Every strand crosses every cut: rank M = N^4 in each cut (checked numerically N = 2,3).
Same scheme for D = 2: 2 strands, cut 0 = (a0 a1|b0 b1), cut 1 = (a0 b1|b0 a1) = G18g's crossing X(x,y,u,v) = [x=v][y=u].

G18h vs G18g: the SAME object.  G18h's T_eps (lambda = 1) = X_3 with the identity on each strand replaced by eps, i.e. X_3 after an
invertible change of basis on the b-end of every strand (cut ranks unchanged); G18h's cost of a piece, sum_j log rank_j / (D log M),
is G18g's rho_D.  The difference is only the CLASS of pieces (G18h: eps-matchings, ungrouped, N = 2) and the regime (smallest case).
So the true analogue of the D = 2 crossing X = R + C is X_3 itself (G18g 1bis), with D = 3 cut-rank costs.

Feasible triple [definition]: (p1,p2,p3) is feasible for N if X_3 = A_1 + A_2 + A_3 with rank_j(A_j) <= p_j (rank across cut j).
Exponents a_j = log p_j / log M.  "Line" -> hyperplane a1 + a2 + a3 = 1; "circle" -> sphere a1^2 + a2^2 + a3^2 = 1
(balanced points: a = 1/3, resp. a = 1/sqrt3 = 0.577; with M = N^4: p = M^{1/3}, resp. M^{0.577}).

How triples relate to the cost rho of pieces [P, D-version of G18g Thm 1(i)]: if X_3 = sum_{i in I} P_i with
prod_j rank_j(P_i) <= M^{3 rho} for all i, then for every (s1,s2,s3) with s1 s2 s3 = M^{3 rho} the triple (|I| s1, |I| s2, |I| s3)
is feasible (put P_i into the A_j with the least j such that rank_j(P_i) <= s_j: one exists since otherwise the product exceeds
s1 s2 s3).  Hence:  rho >= (1/3) * sup over frontier points of (a1+a2+a3) - O(log|I|/log M).
In particular the LEAST BALANCED p gives a NECESSARY condition: rho >= log(p_min / |I|) / log M.
CAUTION: the converse (G18g Thm 1(ii), telescoping) has no D = 3 analogue that I can prove: a single triple decomposition gives
pieces of cost (p_j, M, M).  Pieces need a 2-parameter nested family of triple decompositions.  So: a feasible balanced triple below
the sphere is necessary, not sufficient, for rho < 3^{-1/2}; a balanced least p ABOVE M^{1/sqrt3} (asymptotically) would be a no-go.
Also: the vertices (M,0,0) are on the frontier (rank count), so sup sum a_j >= 1 always: the frontier dipping below the hyperplane in
the middle would NOT lower rho below 1/3.

## 2. Lower bound: D = 3 version of Theorem C (task 2)

THEOREM G [P].  If (p1,p2,p3) is feasible (any field), then for all integers 1 <= e_0,e_1,e_2,e_3 <= N
     e0 e1 e2 e3  <=  p1 + p2 * min(e0 e1, e2 e3)^2 + p3 * min(e0 e2, e1 e3)^2        (target cut 0; same for every cut by symmetry)
Proof.  Compress strand delta to an e_delta-dimensional subspace V_delta at BOTH ends (injection at a_delta, a left inverse
compatible with it at b_delta, so the strand becomes I_{e_delta}); compression on single legs never increases any flattening rank.
X_3 becomes the identity of size prod e_delta in cut 0.  A_1 compressed has cut-0 rank <= p1.  A_2 = sum_{l<=p2} f_l (x) g_l across cut 1,
f_l on legs {a0,a1,b2,b3}, which cut 0 splits as {a0,a1} | {b2,b3}: rank_0(f_l) <= min(e0 e1, e2 e3); likewise g_l on {b0,b1,a2,a3}.
rank_0(f (x) g) <= rank_0(f) rank_0(g).  Same for A_3 with the split {a0,a2} | {b1,b3}.  Subadditivity of rank.  QED
* D = 2 specialisation (2 strands, e = (N,k)): N k <= p + q k^2, i.e. p >= k(N - kq): exactly Theorem C(a) / G18g Theorem 2.
* On the faces it is tight in exponent: p3 = 0, e = (N,N,t,t): p1 >= N^2 t^2 - p2 t^4 -> p1 p2 >= N^4/4 = M/4 (the line a1 + a2 >= 1).
* In the interior it is WEAK: with e_delta = N^{t_delta}, the balanced bound is p >= N^{max_t [T - 2 max(alpha,beta)]}, T = sum t,
  alpha = min(t0+t1,t2+t3), beta = min(t0+t2,t1+t3); the max is 1 (e.g. e = (N,N,N,1): p >= N^3/(1+2N^2)).  So the PROVED balanced
  floor for triples is  p >= N/2 - O(1) = M^{1/4}/2, exponent 1/4, BELOW the hyperplane value 1/3.
* Is the hyperplane the right floor for triples?  NOT PROVED.  What is proved: (i) for pieces, rho >= 1/3 (rank count, G18g) - this does
  not need the triple region; (ii) triples: sum a_j >= 1 on the faces (Theorem C), a >= 1/4 at the balanced point (Theorem G).
  Exact constant: on the faces the 2/9 of G19e (p1 p2 >= (2/9) M); balanced: p >= N^3/(2N^2+1).  [O] whether p1p2p3 >= M^{1-o(1)}.
  The loss is structural: killing two pieces needs restrictions on one side of cut 1 AND one side of cut 2, whose union contains both
  ends of two strands; for D = 2 one side of one cut suffices.

## 3. Computations (task 3)

Solver: trilib.py (variable projection: A_0 eliminated by Eckart-Young in cut 0; A_j = U_j V_j^T; LM + projected Gauss-Newton + CG,
gauge = balanced SVD) - the G19c design, rewritten for a general strand box.

### 3a. Calibration on D = 2 (out_calib2.txt) [N]: reproduces G19b exactly
 N=2 (1,1) 6/6 exact; N=3 (3,1),(1,3) exact, (2,1) plateau 3.000 (6/6), (2,2) plateau 2.000 (6/6); N=4 (6,1) 6/6, (5,1) plateau 3.0,
 (5,3) 1/6 exact (others 1e-5..1e-3), (4,4) 3/6 exact, (5,2) plateau 1.02-2.0, (3,3) 3e-3 with growing norms (G19b: border only).
### 3b. Calibration on D = 3 faces [N]: (p,q,0) at N is the crossing with N' = N^2
 N=2 (out_scan3_N2.txt): face points reproduce the G19b N'=4 table: (0,1,6),(0,3,5),(0,4,4) exact, (0,1,5) plateau 3.0, (0,2,5) 1.03-2.0,
 (0,3,3) border.  N=3: (36,1,0) exact 2/2 (= N'(N'-1)/2, G19b Thm 3), (35,1,0) plateau 3.000 2/2.
### 3c. Calibration on a D = 3 INTERIOR point with a known exact solution - the solver is WEAK there
 N=2, exact (1,2,3) exists (section 3d).  cal123.py: 15 starts (3 orderings x 5, 40 s each): 1/15 reaches 3e-27; the other 14 end at
 tails 1e-6..1e-1 with norms 60-130 (border-like drift).  So in the interior, "small tail + growing norm" does NOT distinguish
 border-only from feasible; only O(1) plateaus in all starts ((1,1,1): 9.83, (1,1,2): 5.0 at N=2) are evidence of infeasibility.
 MORE calibration failures (decisive for how to read the numerics):
 * N=3, (9,9,9) is EXACTLY feasible (product point, 3e) but a random start (120 s) ends at tail 0.92, norm 161.
 * N=4, (40,40,40) is feasible (it dominates the exact (36,16,16)) but a random start (240 s) ends at tail 0.04; (32,32,32): 0.96.
 * N=3 random starts: (16,16,16) and (20,20,20) exact fast, (15,15,15) 1e-11 (bounded norm), (14) 7e-4, (13) 4e-3, (12) 1e-2 (growing
   norms), (8) 2.3-2.7, (5) 18-20 - all of this is solver difficulty, NOT a threshold: (9,9,9) is exactly feasible.
 => Random-start numerics in the D = 3 interior are UNCALIBRATED (they fail on known exact points at N = 2, 3, 4).  Only continuation
    from an exact point (cont.py) behaves: from the N = 3 product point + noise 0.05, 3/3 return to tail ~1e-12 (slow, sublinear:
    the solution set is degenerate).

### 3d. Exact points [E]
 * N = 2 (match_n2b.py; exact_n2_2_1_3.json, exact_n2_2_1_4.json): X_3 (eps form) = signed sum (+-1) of 14 eps-matchings grouped into
   A_0, A_1, A_2 with exact ranks (2,1,3) [and (2,1,4)]; sum re-expanded and ranks computed in sympy over Q.  So (1,2,3) is feasible:
   balanced least p <= 3.  (1,1,1): plateau 9.83 in 6/6 starts, (1,1,2): plateau 5.0 6/6 [N, evidence of infeasibility];
   (2,2,2): 6/6 starts at tails 2e-4..6e-2 with sigma_{p0} 400-860 diverging (border-like) - UNDECIDED given 3c.
   N = 2 eps-matching class (match_n2.py, 20000 random sparse solutions): Pareto set {(0,0,16),(0,1,6),(0,4,4),(1,2,3)}.
 * PRODUCT POINT, all N (prod_exact.py, integer identity checked entry by entry, ranks over Z with flint): pair the strands (0,3),(1,2);
   both pairs are crossings for (cut 0) vs (cuts 1 AND 2) simultaneously.  With Y = 1 - F (F = swap of the pair):
        X_3 = [Y_03 (x) Y_12]  +  [I_03 (x) F_12]  +  [F_03 (x) Y_12]
   ranks (cut 0, cut 1, cut 2) = (t^2, N^2, N^2), t = N(N-1)/2:  N=2: (1,4,4); N=3: (9,9,9); N=4: (36,16,16).  (all cut ranks of the
   three terms: N=3 [[9,81,81],[81,9,9],[27,9,9]]).  More generally any two D=2 feasible pairs (p,q),(p',q') give (pp', N^2 q', N^2 q).
   Balanced only at N = 3, by the ACCIDENT t = N(N-1)/2 = N (the analogue of the D = 2 accident N(N-1)/2 = 1 at N = 2).
   Asymptotically the product point is (M/4, M^{1/2}, M^{1/2}): exponents (1, 1/2, 1/2); best product-family balanced exponent with the
   D = 2 circle: 4/5; even with the D = 2 line: 2/3.
 * N = 3 Brauer (delta-matching) class: span dim 91 < 105 (relation Lambda^4(K^3) = 0), X_3 in span; random search found only
   (12,18,18),(9,18,21) (worse than the product point, which is in this class).  N = 4: the 105 matching tensors are independent
   (rank 105): the class is EMPTY beyond X_3 itself - the matching mechanism of G18h is a small-N accident (needs N < 4 = #strands).

### 3e. N = 4 by continuation from the exact product point (cont.py, cont_more.py; one start each, 280 s + 200 s) [N, weak]
 target          tail after 280 s / +200 s     |x|     |A_1|,|A_2|      sigma_p0 / sigma_p0+1
 (32,32,32)      3.2e-8 / 1.6e-8               19.6    18.6, 21.2       0.35 / 3e-5       stable norms, slow decrease: probably feasible
 (28,28,28)      7.0e-8 / 2.2e-8               26.4    36.0, 42.7       0.032 / 4e-5      stable norms: possibly feasible
 (24,24,24)      7.8e-5 / 3.7e-5               42.8    128, 131         0.033 / 1.5e-3    |A_j| >> |X_3| = 16: border-like
 (30,20,20)      7.7e-6 (60 s)                  25.6                                       undecided
 None of these is verified exactly (no rounding done); convergence is sublinear even near exact points (3c), so these are hints only.

## 4. Table: least balanced p, D = 3 (M = N^4)
 N  M    exact upper bound          numerics                                proved lower (Thm G)  M^{1/3}  M^{1/sqrt3}  M^{0.745}
 2  16   3  [(1,2,3), eps-matchings]  (1,1,1) plateau 9.83 -> p >= 2 [N];      1                     2.52     4.96         7.9
                                      (2,2,2) border-like, undecided
 3  81   9  [product point, accident]  random starts uncalibrated (fail at 9)    2                     4.33     12.6         26.4
 4  256  36 [(36,16,16) product pt]    (32),(28) tails ~2e-8 stable norms;       2                     6.35     24.6         62
                                      (24) border-like
 exponent log p / log M of the exact upper bounds: 0.40 (N=2), 0.50 (N=3), 0.65 (N=4); numerical N=4 candidates 0.625 (32), 0.60 (28).
 Floor 1/3, sphere 0.577, 2-direction reuse 0.745 (G18g), product point asymptotically max-exponent 1.

## 5. Trend and product point (task 4)
* D = 3 product point = (t^2, N^2, N^2), t = N(N-1)/2 (exact identity above).  It explains the N = 3 value completely (t = N: balanced
  by accident) and dominates the N = 4 exact bound.  Asymptotically it is (M/4, M^{1/2}, M^{1/2}): NOT a balanced gain.
* The N = 2 value 3 comes from the SL2/eps-matching relations (span 14 < 105), the N = 3 Brauer class has relations from
  Lambda^4(K^3) = 0, and at N = 4 the matching tensors are independent: the matching mechanism dies exactly at N = #strands = 4.
* Exact exponents RISE 0.40 -> 0.50 -> 0.65 (upper bounds); the N = 4 numerics (not exact) sit at ~0.60-0.625, i.e. AT or ABOVE the
  sphere value 0.577, nowhere near the hyperplane 1/3.  Same pattern as D = 2 (small-N values below the line, then the product point
  takes over), with the caveat that at N = 4 the product point is already beaten numerically (32 < 36), so the frontier is not
  "everything is the product point".
* Interpretation limits: (i) triples are only a NECESSARY condition for pieces (section 1), so even a triple frontier on the hyperplane
  would not give rho = 1/3; (ii) M^{o(1)} factors are huge at M = 256 (|I| alone is a factor 2-4).

## 6. Verdict
No evidence that D = 3 beats D^{-1/2} beyond the smallest case.  The only sub-sphere balanced values (N = 2: 0.40, N = 3: 0.50) are
explained by accidents (matching relations at N < 4; t = N at N = 3).  At N = 4 exact: 0.65; numerical: ~0.60.  And the object that
can be computed (triples) is weaker than the one that matters (pieces).
Most informative next computation: decide the N = 4 balanced frontier exactly - many continuation starts at (24..31)^3 from the
product point and from the (32,32,32) numerical point, then rational rounding of a converged solution (gauge-fix with the GL(N)^4 strand
symmetry, as G19b's sparsify.py) - about 2-4 CPU-hours; then the same at N = 5 (M = 625, product point (100,25,25); ~1 CPU-hour per
start).  If the N = 5 balanced exponent stays >= 0.6, the D^{-1/2} law survives the triple test at D = 3.  Beyond that, the real
object is a 2-parameter piece family (D = 3 analogue of G18g Thm 1(ii)), for which no solver exists yet.

## 7. Adversarial self-review
* The random-start solver FAILED on known-feasible interior points at N = 2 (1/15), N = 3 (9,9,9) and N = 4 (40,40,40).  Every
  interior "not found" in this file is therefore without value; I report none of them as infeasibility except O(1) plateaus at N = 2.
* N = 4 values 32/28 are numerical with sublinear convergence and no exact rounding: hints, not results.
* Theorem G is a 5-line compression argument (P, unreviewed by a second agent); it reproduces Theorem C for D = 2.  The claim that its
  balanced optimum is exponent 1/4 is a small LP over t in [0,1]^4 done by hand (max of min(|t0+t1-t2-t3|,|t0+t2-t1-t3|) = 1).
* The "triples are necessary, not sufficient" point (section 1) is a proof for one direction and an ABSENCE of proof for the other;
  a cleverer D = 3 telescoping might exist.
* G18h's eps-matching result (rho_3 = 1/2 in that class) is about pieces; my N = 2 (1,2,3) is about triples: different statements.
* Exact checks: N=2 (1,2,3)/(1,2,4) sympy over Q; product point N = 2,3,4 integer identity + flint ranks.  Nothing else is exact.

## 8. Files (this folder)
 trilib.py (solver), calib2.py/out_calib2.txt (D=2 calibration), scan3.py/out_scan3_N2.txt/scan3_N2.json (N=2 scan), run1.py + runs.log
 (random starts), cal123.py (interior calibration), match_n2.py, match_n2b.py, exact_n2_2_1_3.json, exact_n2_2_1_4.json (exact N=2),
 match_gen.py (Brauer class N=3,4), prodpt.py, prod_exact.py (product point), cont.py, cont_more.py, cont_N4_*.npz, o_*.out.
