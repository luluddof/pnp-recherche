# G18g — the Hurwitz / Naslund problem and the gain rho of Claude's 3-SAT algorithm

Status marks: [P] proved here (short proof given), [E] exact computation, [N] numerical only,
[L] literature (search snippets / abstracts; full texts not re-read unless said), [H] heuristic, [O] open.
Nothing here improves a published bound on P(N) or the running time of any algorithm.

## 0. Notation
X = X_N = the crossing [x=v][y=u], x,y,u,v in [N].
rr(A) = rank of A as a matrix (x,y)|(u,v) ("row-wise"); rc(A) = rank as a matrix (x,u)|(y,v) ("column-wise").
rr(X) = rc(X) = N^2.  Row-wise X is the swap matrix F; rc is the operator-Schmidt rank of the row-wise matrix.
A pair (p,q) is FEASIBLE for N (over a field K) if X = R + C with rr(R) <= p and rc(C) <= q.
P(N) = min{p+q : (p,q) feasible}   (Naslund's question; = width of sum_ij x_i y_j x_i y_j, HWY11 Lemma 3.5).

## 1. Task 1 — the variant that matters

What the algorithm uses (Definition 7 and the proof of Proposition 3 of Claude's paper, re-read):
integer functions A_i, i in I, with sum_i A_i = d X (d >= 1), each A_i factoring row-wise through c_row
bits and column-wise through c_col bits, c_row + c_col <= 4 rho lambda + O(log lambda), N = 2^lambda.
One level of the recursion costs, per rectangle, sum_i (rr(A_i) rc(A_i))^e in place of N^{4e}
(e = current exponent; "rows or columns are chosen once per term"), so
   e' = rho e + log|I| / (4 lambda):
|I| must be N^{o(1)} (poly(lambda) in the paper); |I| = N^beta would leave a fixed point e* > 0.

So the object is:  X = sum of N^{o(1)} pieces, minimise  max_i rr(A_i) * rc(A_i) =: N^{4 rho}.

(a) Field. Rank over Q, with integer factorizations of poly(lambda) bits computable in 2^{O(lambda)}.
    A number field of degree D costs a factor <= D in every rank [P: if A = u v^T over K with Q-basis
    w_1 = 1,...,w_D, the w_1-coordinate of A is sum_j u_j (sum_l c_jl v_l)^T, rank <= D], so Q(i), and any
    field of degree N^{o(1)}, is as good as Q. R or C without a small field of definition: not usable as such.
(b) It is NOT a sum of squares. No positivity is needed; the right classical variant is "products of bilinear
    forms" (sum a^2)(sum b^2) = sum_k f_k f'_k (HWY11 Prop. 4.1: P(N) within a factor 4 of the least number of
    products). Consequences: Hurwitz-Radon (P = N only for N = 1,2,4,8) and the integer bound Omega(N^{6/5})
    (HWY, sign/integer SOS formulas) do NOT constrain the algorithm; e.g. over Q, (1,1) is feasible for N = 2
    with P(2) = 2 but through a non-SOS identity (section 4). Field dependence is real for SOS: [12,12,18]
    exists over Q(i), not over R (Zhang-Zhu, arXiv:2605.00590, May 2026 [L]) — by (a) the algorithm could use it.
(c) Square [N,N,P] only: x,y,u,v have the same length. Rectangular [r,s,n] formulas are not what is needed;
    what IS needed is unbalanced in the RANKS (p != q), see (d).
(d) The algorithm's object is neither weaker nor stronger than "P(N) small"; it is the whole trade-off curve.
    THEOREM 1 [P]. (i) If X = sum_{i in I} A_i with rr(A_i) rc(A_i) <= Q for all i, then for every real s > 0
    the pair (|I| s, |I| Q/s) is feasible. (ii) If (p_1,q_1),...,(p_m,q_m) are feasible, p_j increasing,
    q_j decreasing, and (p_0,q_0) = (0,N^2), (p_{m+1},q_{m+1}) = (N^2,0), then X is a sum of m+1 pieces with
    rr(A_j) rc(A_j) <= (p_j + p_{j-1})(q_j + q_{j-1}).
    Proof. (i) R = sum of the A_i with rr <= s, C = the rest (each has rc < Q/s). (ii) Let X = R_j + C_j;
    A_j = R_j - R_{j-1} = C_{j-1} - C_j, j = 1..m+1 (R_0 = 0, R_{m+1} = X); the sum telescopes to X.  QED
    COROLLARY. With T(theta) := least b such that (N^theta, N^b) is feasible (asymptotically, up to N^{o(1)}),
        4 rho* = sup_{0<=theta<=2} ( theta + T(theta) )
    for O(1/eps) pieces, up to eps (uniformity / bit length of the pieces NOT addressed: existence only).
    THEOREM 2 (hyperbola bound) [P]. Over every field, if (p,q) is feasible then p >= k (N - k q) and
    q >= k (N - k p) for every integer 0 <= k <= N. Hence p + q >= N, and p q >= 3 N^2 / 16 when p,q >= 1.
    Proof. Write R = sum_{j<=p} f_j(x,y) f'_j(u,v), C = sum_{l<=q} g_l(x,u) g'_l(y,v). Take independent
    a_1..a_k in K^N and B = {b : a_i^T g'_l b = 0 for all i,l}, dim B >= m := N - kq, basis b_1..b_m.
    Contract X = R + C with a_i(y) b_l(v) over y,v: the C-part vanishes, X gives b_l(x) a_i(u), R gives
    sum_j (f_j a_i)(x) (f'_j b_l)(u). As a matrix with rows (i,x) and columns (l,u) the right side has rank <= p
    and the left side is (a-matrix) (x) (b-matrix) composed with a swap of factors, of rank k m. The second
    inequality follows by exchanging y and u. For q <= N/2 take k = round(N/2q): k(N-kq) >= N^2/4q - q/4. QED
    (Probably folklore; k = 1 is the usual P(N) >= N. I did not find the product form stated. [L: not found])
    So T(theta) >= 2 - theta (the line a + b = 2), and Appendix B gives T(theta) <= sqrt(4 - theta^2)
    (the circle a^2 + b^2 = 4: thresholding the pieces A_s at s = sigma h gives p = N^{sqrt(2+sigma)},
    q = N^{sqrt(2-sigma)}), whose sup of a + b is 2 sqrt 2 at the balanced point: rho = 1/sqrt 2.
    CONSEQUENCES (this corrects the "equivalence" used in our earlier notes):
      * rho >= log P(N) / (2 log N)  [P] (known: N1). So P(N) >= N^{1+eps}  ==>  rho >= (1+eps)/2.
      * rho -> 1/2  <==>  the WHOLE hyperbola is attained: (N^theta, N^{2-theta+o(1)}) feasible for every theta.
        This is stronger than P(N) <= N^{1+o(1)} (the single point theta = 1). From P(N) <= N^c alone one
        gets only rho <= (2+c)/4, i.e. >= 3/4 [P: step curve T <= 2 on [0,c), T <= c on [c,2]]; with exactly
        two pieces (2+c)/4 is also a lower bound [P: rr(A_1) >= N^2/2, then rc(A_2) >= N^2/2 or trivial, and
        max(rc(A_1), rr(A_2)) >= P(N)/2].
      * a floor rho > 1/2  <==>  T(theta) >= 2 - theta + eps for SOME theta in (0,2). Naslund's
        P(N) >= N^{1+eps} is the case theta = 1 and implies it; the converse is not known. So proving
        rho > 1/2 is no harder, and possibly easier, than Naslund's problem ("unbalanced Naslund problem") [O].
(e) Approximate / piece-wise objects. A border (degenerating) decomposition of degree e converts to an exact
    one with e+1 times as many pieces by interpolation, so it gives the same rho when e = N^{o(1)} [H, standard].
    "Identity up to a low-rank error" is the same object (the error is one more piece). Identities modulo 2^k
    (OpenAI Appendix D, G14b) are a different model, not covered here. No loophole below rho = 1/2 in the
    exact model: rho >= 1/2 follows from rr(X) = N^2 alone [P].
(f) Reformulation [P]. (p,q) is feasible  <==>  there is a subspace U of M_N(K) of codimension <= q and
    matrices alpha_j, beta_j (j <= p) with  Z^T = sum_j alpha_j Z beta_j^T for all Z in U
    (contract X = R + C with Z(y,v), Z in the annihilator of span g'_l). "Transposition restricted to a
    subspace as a short sum of two-sided multiplications." U = symmetric matrices gives (1, N(N-1)/2).
    Theorem 2 reads: p >= dim A * dim B whenever A (x) B is contained in U.

### 1bis. Multi-direction (D >= 3) generalisation (answer to the coordinator's question)  [mostly H]
Natural object: X_D = tensor on the 2^D corners c in {0,1}^D of a D-box, X_D = prod over the 2^{D-1} body
diagonals {c, c-bar} of [z_c = z_cbar]. Cut j separates c_j = 0 from c_j = 1; every diagonal crosses every
cut, so the rank across each of the D cuts is M := N^{2^{D-1}}. (For large grids almost every pair of grid
points differs in all D coordinates, so this is the dominant gadget; lower-dimensional boxes exist too.)
Pieces: sum_i A_i = d X_D, cost of a piece = sum_j log r_j(A_i), rho_D := max_i sum_j log r_j(A_i) / (D log M).
 * Floor rho_D >= 1/D [P: some piece has r_1 >= M/|I|]. G18c's law D^{-1/2} vs this floor 1/D is, in the
   picture of 1(d), "sphere versus hyperplane": rank profiles M^{a_j} with sum a_j^2 = 1 (a D-direction
   Heisenberg-like family, balanced point sum a_j = sqrt D) against sum a_j = 1. For D = 2 these are exactly
   the circle and the line above. A real improvement = a feasible profile strictly inside the sphere. [H]
 * Analogue of P(N): P_D(N) = least number of terms, each of rank one across one of the D cuts, with sum X_D.
   P_2 = P.  rho_D >= log P_D / log M [P: min_j r_j <= geometric mean].
   Bounds [P]: P(N) <= P_D(N) <= P(N^{2^{D-2}}).  (Upper: use two cuts only; with respect to cuts 1,2 X_D is
   the 2-direction crossing of strings of 2^{D-2} letters. Lower: fix all legs outside one quadruple
   {c, c', cbar, cbar'} (c' = c with bit 1 flipped) to 0; on these four legs cut 1 is the column split and every
   other cut is the row split.) So N <= P_D(N) <= M^{1/sqrt2 + o(1)}; the natural conjecture is
   P_D(N) >= c M^{1/D}, i.e. prod_j p_j >= c M (the D-direction form of Theorem 2) [O, not proved].
 * Commutative shadow: one set of N variables per body diagonal; X_D becomes prod_{delta=1}^{2^{D-1}}
   (sum_i a^{(delta)}_i^2), a product of 2^{D-1} sums of squares (FOUR for D = 3, not three), and a term of
   rank one across a cut becomes a product of two 2^{D-1}-linear forms. For D = 2 this shadow is equivalent
   to P up to a factor 4 (HWY). For D >= 3 it is not expected to be: iterated Hurwitz identities give the
   shadow with 8 products at N = 8 for every D (|((ab)c)d|^2), far below M^{1/D}; the cut structure is lost.
   So for D >= 3 the right problem is the non-commutative one (P_D, or the profile region), not a
   sums-of-squares problem.
 * A 2-direction family used inside D directions (pair the diagonals by flipping bit j; each pair is a
   2-direction crossing for (cut j, any other cut)) gives profiles (b, a, ..., a): with Appendix B,
   rho = sqrt((D-1)^2 + 1)/D (0.745 for D = 3), worse than 1/sqrt2. Genuinely D-directional pieces are needed. [P]

## 2. Task 2 — state of the art (as of Nov 2026)
 * Real SOS [N,N,P]: P = N iff N in {1,2,4,8} (Hurwitz 1898); general upper bound O(N^2/log N) (Hurwitz-Radon
   based; Lam-Yiu tables for small cases); lower bound N (nothing better over R in general); integer/sign
   formulas: Omega(N^{6/5}) (Hrubes-Wigderson-Yehudayoff, "An asymptotic bound on the composition number of
   integer sums of squares formulas", Canad. Math. Bull.) [L].
 * Over C / products of bilinear forms: O(N^{1.62}) (P. Hrubes, "A subquadratic upper bound on sum-of-squares
   composition formulas", CCC 2024, LIPIcs 300 no. 12, ECCC TR24-026; companion ECCC TR24-099 for any field of
   char != 2 and non-commutative variants) [L, via Claude's paper and search].
 * HWY11 (Hrubes-Wigderson-Yehudayoff, "Non-commutative circuits and the sum-of-squares problem", JAMS 2011):
   a bound N^{1+eps} over C would give exponential non-commutative circuit lower bounds [L].
 * Naslund's question (Nas20, Problem 12): P(N) >= N^{1+eps}? Open. arXiv:2310.11356 ("The interplay between
   bounded ranks of tensors arising from partitions", 2023) announces power lower bounds for products of
   identity tensors, "a qualitative answer to a question of Naslund" — for the crossing a power bound N^c with
   c <= 1 is already trivial; I did not read the paper, so I do not know its exponent [L, unchecked].
 * 2025-2026: Zhang-Zhu arXiv:2605.00590 (field dependence, [12,12,18]); Brand-Kaski-Wang arXiv:2607.02241
   (partition rank and circuit lower bounds, d >= 3 tensors; not about the crossing); Lampert-Moshkovitz
   arXiv:2509.06294 (partition rank of det). Nothing found that improves the exponent 1.62 other than the
   consequence of Claude's Appendix B:  P(N) <= N^{sqrt2+o(1)} over Q (Siemaszko, Zenodo 10.5281/zenodo.24024689,
   7 Nov 2026; Debnath; our G4a) — these two references are taken from JOURNAL.md, not re-checked today; a web
   search did not surface them.
 * Lower bounds for the variant that matters: only P(N) >= N (and now p q >= 3N^2/16, Theorem 2).

## 3. Task 3 — original attempt
Theory: Theorems 1, 2, the corollary, the reformulation 1(f), and section 1bis.
Exact computations (exact.py -> out_exact.txt; python-flint ranks over Q; every identity checked entry by entry):
 * base identities: one-bit X = [x=u][y=v] - [x!=y](-1)^x [u!=v](-1)^u (N = 2); "I - F": X = [x=u][y=v] + (X - [x=u][y=v])
   (N = 3,4,5); Hurwitz N = 2,4,8 (sign function (27) of the paper), R + C = 2X.
 * Kronecker powers A_i = coefficient of T^i in prod_blocks (R + T C), all verified to sum to d^t X:
     family          N    pieces (rr,rc)                       max rr*rc   rho      feasible pairs (thresholds)
     one-bit^1       2    (1,4)(4,1)                           4 = N^2     0.5000   (1,1)
     one-bit^2       4    (1,16)(7,7)(16,1)                    49          0.7018   (1,6)(6,1)
     one-bit^3       8    (1,64)(10,37)(37,10)(64,1)           370         0.7109   (1,28)(10,10)(28,1)
     one-bit^4       16   (1,256)(13,175)(67,67)(175,13)(256,1) 4489       0.7583   (1,120)(13,67)(67,13)(120,1)
     I-F, N=3        3    (3,9)(9,1)                           27          0.7500   (3,1)
     (I-F, N=3)^2    9    (9,81)(45,17)(81,1)                  765         0.7555   (9,17)(36,1)
     Hurwitz 4       4    (4,16)(16,4)                         64          0.7500   (4,4)
     Hurwitz 8       8    (8,64)(64,8)                         512         0.7500   (8,8)
     Hurwitz 4^2     16   (16,256)(112,112)(256,16)            12544       0.8509   (16,96)(96,16)
   Every feasible pair satisfies Theorem 2 (asserted in the script); min p q / N^2 seen: 1/4 (N=2), 1/3 (N=3).
   Small N flatter the families: asymptotically one-bit^t tends to rho = (2 + log2 3)/4 = 0.896 and
   P-exponent 1 + (log2 3)/2 = 1.79; Hurwitz-8 powers to 0.901 (paper, App. H) and 11/6; Hurwitz-4 to 0.948.
 * Kronecker products do NOT make log P / log N subadditive [P, negative]: (R1+C1) (x) (R2+C2) contains the
   mixed terms R1 (x) C2, which are rank one in neither direction (rr = p1 N2^2, rc = N1^2 q2); a base identity
   of exponent 1 (N0 = 2, 4, 8) only yields exponent >= 1.79 in the limit. For a Hurwitz-type identity
   (row-wise, in the arrangement where X = Id, R is diagonalizable with p nonzero eigenvalues, none equal to d)
   the power pieces have rr(A_i) = sum_{z<=i} C(t,z) (N0^2 - p)^z p^{t-z} (paper, App. H; checked exactly:
   Hurwitz 4^2, i = 1: 16 + 2*12*4 = 112), so the middle piece gives
   rho >= (2 + log2((N0^2 - p) p)) / (4 log2 N0) >= 3/4 for t large (0.901 for N0 = p = 8). Real SOS
   identities do not tensor at all (the tensor product of two normed bilinear maps is not normed; classical).
   Best recursive (Kronecker) family: rho = 0.896, exponent 1.79. Best known overall: Appendix B, a graded
   (not Kronecker) family: exponent sqrt2, rho = 1/sqrt2 = 0.7071.
Numerical small cases (lm.py -> out_lm_N3.txt; Levenberg-Marquardt, random restarts) [N]:
   N = 2: feasible set = {p,q >= 1} plus (4,0),(0,4): P(2) = 2 [P by Theorem 2 + exact identity].
   N = 3: (3,1),(1,3) feasible [E]; (2,1),(1,2),(2,2) NOT FOUND over R (22-35 restarts) nor over C (7-11
   restarts), with all local minima at residual sqrt2 or sqrt3. If truly infeasible, P(3) = 4 while Theorem 2
   gives only >= 3: the hyperbola bound is not tight at N = 3. "Not found" is not a proof.
   N = 4: (6,1) feasible; (5,1) not found in 2 restarts; the runs were cut short (machine load ~17 on 4 cores).
   Open small questions: is (2,2) feasible for N = 4 (allowed by Theorem 2; would give P(4) = 4 by a non-SOS
   identity)? Is (5,1) feasible for N = 4, i.e. min rank(F - G (x) G') < N(N-1)/2?
   (tradeoff.py is a first ALS attempt; it gives false negatives — e.g. it missed (1,3) for N = 3 — do not use.)

## 4. Verdict
 * Variant: rank decompositions of the crossing over Q (number fields of small degree allowed), products of
   bilinear forms, not sums of squares; and not P(N) alone but the trade-off curve T(theta).
 * Best today: T(theta) <= sqrt(4 - theta^2) (Appendix B), hence P(N) <= N^{sqrt2+o(1)} and rho = 1/sqrt2;
   lower bound T(theta) >= 2 - theta (Theorem 2), hence rho >= 1/2. The gap circle/line is open at every theta.
 * Loophole? None below 1/2. The weaker-looking object does not help: rho = 1/2 needs MORE than
   P(N) = N^{1+o(1)} (whole hyperbola). On the lower-bound side it opens something: rho > 1/2 would follow from
   an unbalanced bound p q >= N^{2+eps} at any single theta, formally weaker than Naslund's statement.
   No such bound is proved here.

## 5. Adversarial self-review
 * Theorems 1, 2 are elementary and I checked Theorem 2 on 20 exact feasible pairs; but they are existence
   statements. The corollary "4 rho* = sup(theta + T)" ignores uniform computability and bit length
   (Definition 7(c)) and the cost of |I| = O(1/eps) pieces; it is about ranks only.
 * The claim "the cost per rectangle is sum_i (rr rc)^e" is my reading of the proof of Proposition 3
   (min(U_row,U_col) <= average); I did not re-derive Theorem A.
 * T(theta) is defined with a limsup/"up to N^{o(1)}" convention that I did not make fully rigorous
   (monotonicity in N holds by restriction: (p,q) feasible for N is feasible for N' <= N).
 * The Appendix-B circle uses the weight bounds rr <= sum_{w<=2h+s} p(w) of G4a for thresholds (sums of
   consecutive pieces pass through the same graded components); exact ranks at small h are far from asymptotic.
 * Section 1bis is a model of a D-direction algorithm that I defined myself by analogy; I did not check it
   against OpenAI's construction or G18c's derivation. P_D >= c M^{1/D} is a conjecture.
 * Novelty of Theorem 2 / Theorem 1 not established (literature check was short); Theorem 2 is the kind of
   statement that may be folklore.
 * Literature items marked [L] come from abstracts and search summaries.

## 6. Files (all in this folder)
 notes.md (this file); exact.py, out_exact.txt (exact identities, ranks, thresholds); lm.py, out_lm_N3.txt
 (numerical feasibility N = 2,3); lm2.py (analytic-Jacobian variant, not run to completion); out_lm_N4.txt
 (two N = 4 cases); tradeoff.py (ALS, unreliable); lit/ (two arXiv abstract pages, untrusted data).
