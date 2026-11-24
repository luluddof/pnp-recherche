# G20f — adversarial review of g19a (model DF, claims 1a-1e) and g19b (crossing, claims 2a-2e)
(working notes, saved as I go; [P] = re-derived by me, [E] = my own exact code)

## 2e. Explicit decompositions [E, e_verify.py — own code, Fraction Gaussian elimination + mod p]
Convention: C(x,y,u,v) = sum_l G_l[x][u] Gp_l[y][v], X = [x=v][y=u], R = X - C; rr = rank (x,y)|(u,v), rc = rank (x,u)|(y,v).
 exact_4_5_3: rr(R)=5, rc(C)=3;  exact_5_6_6: (6,6);  exact_5_7_4: (7,4);  exact_6_7_7: (7,7).  All CONFIRMED over Q.
 Entries in {0,+-1}; same ranks mod 2,3,5,7 => the identities are valid over EVERY field (integer identity, rank mod p <= rank over Q,
 and here equal), not only char != 2.

## 2a. Theorem 3 (q = 1)  — CONFIRMED (char != 2, exact), re-derived [P] + exhaustive [E]
Statement: K field, char K != 2, N >= 1, G, G' in M_N(K): rank(F - G (x) G') >= N(N-1)/2; equality at G = G' = I.
Re-derivation: rank(F - C) = rank(1 - CF) (F^2 = 1); Theta = (G(x)G')F, Theta^2 = GG' (x) G'G;
 tr Theta^{2k+1} = tr (GG')^{2k+1} (trace of e(x)f -> Pf (x) Qe is tr PQ), tr Theta^{2k} = (tr (GG')^k)^2  => power sums of
 {mu_i} u {+-sqrt(mu_i mu_j)}_{i<j}; Newton in char 0 on generic entries => integer polynomial identity => any field.
 In char != 2, t = 1 is a simple root of each t^2 - mu_i mu_j; alg. mult >= geom. mult => dim ker <= N(N+1)/2.  No gap.
 Exhaustive (t3_exh.py, G in rank normal form via (G,G') -> (AGB^-1, BG'A^-1), which fixes F): min p = N(N-1)/2 for
 N=2 over GF(2),GF(3),GF(5); N=3 over GF(2),GF(3); N=4 over GF(2).  gf4_q1.py: N=3 over GF(4): 3.
 So the bound also holds in char 2 for N <= 4 (GF(2)), N = 3 (GF(4)); the PROOF does not cover char 2.
 Lean: char-free finite instances are decidable but large (N=3/GF(3): 4*3^9 cases after normal form: feasible with
 `decide`-by-kernel only via a verified checker; not recommended).  The general theorem needs charpoly machinery: medium-hard.

## 2b. Theorem 4 (q = 2 => p >= N, N >= 3)  — CONFIRMED (char != 2, exact decompositions); full argument [P]
 Hypotheses: K any field with char != 2, exact X = R + C, rc(C) <= 2, rr(R) <= p.  Work over Kbar (ranks unchanged).
 1. K := ker R = {W in M_N : W^T = G1 W G1'^T + G2 W G2'^T}, codim <= p.  Segre Sigma = {[ab^T]} in P^{N^2-1}, dim 2N-2,
    irreducible.  If p <= N-1 (<= 2N-2) then S = Sigma n P(K) is non-empty and every component has dim >= 2N-2-p >= N-1.
 2. (a,b) in S: b a^T = A(a)B(b)^T, A(a) = [G1 a, G2 a], B(b) = [G1'b, G2'b] (N x 2).  rank 1 = rank(AB^T) >= rk A + rk B - 2
    => rk A <= 1 or rk B <= 1.  S1, S2 closed, S = S1 u S2; an irreducible component Z lies in one, say S1.
 3. On S1: Im(A B^T) c Im A (dim 1) and = span b, so b is determined by a; pr_1|Z injective, image closed irreducible of dim
    dim Z >= N-1 => image = P^{N-1} => rank [G1 a, G2 a] <= 1 for all a in Kbar^N.
 4. Linear space {A(a)} of N x 2 matrices of rank <= 1: common column line (G_l = c d_l^T => rank C <= N, p >= N^2 - N >= N)
    or common row line (w2 G1 = w1 G2, G's dependent => rc(C) <= 1 => Theorem 3 => p >= N(N-1)/2 >= N for N >= 3).
    S2 is the mirror (b determined by a via B A^T, conclusions on G').  Contradiction with p <= N-1.   No gap found.
 Where char != 2 enters: ONLY via Theorem 3 in step 4.  Improvement [P]: for N >= 4 replace Theorem 3 by Theorem C(a) of g19e
 with q = 1, k = floor(N/2): p >= k(N-k) >= N.  So Theorem 4 holds in EVERY characteristic for N >= 4; in char 2 and N = 3 it
 holds over GF(2) (g19e q2n3.c exhaustive) and the needed q=1 input holds over GF(2), GF(4) (exhaustive), not proved over F2bar.
 Consequences: (2,2), (2,1), (1,2) infeasible at N = 3 and P(3) = min(p+q) = 4 [(3,1) from G=G'=I; (0,q) needs q = 9].
 Not proved for border decompositions (as the notes say).  Lean: hard (projective dimension theorem, Segre); not a candidate.

## 2c. Theorem 5  — PARTLY: literal statement FALSE for q = 1; correct with the extra hypothesis k <= N
 Counterexample to the statement as written ("every integer k >= 1 with m = N - k(q-1) >= 1"): q = 1 gives m = N for all k,
 so k is unbounded; N = 2, k = 3: k' = 2, bound p >= 6 - 3 = 3, but (1,1) is feasible (rank(F - I(x)I) = 1, t5_exact.py).
 The proof silently needs a_1..a_k independent in K^N, i.e. k <= N.  For q >= 2, m >= 1 already forces k <= (N-1)/(q-1) < N.
 Corrected statement [P, char != 2, exact]: (p,q) feasible, q >= 1, 1 <= k <= N, m = N - k(q-1) >= 1, k' = min(k,m):
       p >= k m - k'(k'+1)/2      (and mirror).
 Full argument (re-derived): B = {b : a_i^T G'_l b = 0, i <= k, 2 <= l <= q}, dim >= m.  Contract y with a_i, v with b_j:
 b_j(x)a_i(u) - G1(x,u) gamma(i,j) has rank <= p as (x,i)|(u,j); apply left inverses of [b_j] on x and of [a_i] on u:
 S - g~ (x) gamma with S: K^k(x)K^m -> K^m(x)K^k the swap, g~ = L G1 L'^T (m x k), gamma (k x m), rank <= p.
 Theta = S^-1(g~ (x) gamma): e(x)f -> gamma f (x) g~ e; Theta^2 = gamma g~ (x) g~ gamma; tr Theta^{2n} = (tr (gamma g~)^n)^2,
 tr Theta^{2n+1} = tr (gamma g~)^{2n+1}; the k x k and m x m products share their nonzero eigenvalues mu_1..mu_k'.
 => charpoly = t^{km-k'^2} prod(t - mu_i) prod_{i<j}(t^2 - mu_i mu_j); eigenvalue 1 multiplicity <= k'(k'+1)/2.  Correct.
 "Derived by hand only": now machine-checked — t5_exact.py: EXACT power sums tr Theta^n = formula for all n <= km, 24 random
 integer cases (k,m) in {(1,3),(2,3),(3,2),(2,4),(3,5),(4,2),(3,3),(4,4)}; t5_rect.py: charpoly numerically; exhaustive
 min rank(S - g~(x)gamma) over GF(2)/GF(3) for (k,m) = (1,2),(2,2),(2,3),(3,2),(3,3): equals the bound in every case (tight;
 the bound held even over GF(2) in these sizes).
 Asymptotics: k = N/(2q-1) <= m, p >= kN - k^2(q - 1/2) - k/2 = N^2/(4q-2) - O(N/q): p q >= (1-o(1)) (N^2/4) q/(q-1/2). OK
 (asymptotic, o(1) as N -> infinity for fixed q; exact at q = 1).  Small-N values (N=4: q=2 -> 2, q=3 -> 1) recomputed: OK.

## 2d. Dictionary  — CONFIRMED (conditional on the accounting "rank <= dimension of the weight window")
 Literal staircase: x < h1: (2^0, 2^{2 lambda}) = (1, N^2); [h1,3h1]: (N,N); > 3h1: (N^2,1).  Threshold below the middle:
 p <= |I|, q <= |I| max(N,1) = |I| N  => (theta,T) = (0,1).  Theorem 2 (pq >= 3N^2/16, confirmed by g19e; p,q >= 1):
 |I|^2 N >= 3N^2/16 <=> |I| >= sqrt(3/16) sqrt N = 0.433 sqrt N.  OK.  (With g19e's 2/9: 0.471 sqrt N.)
 Continuity lemma: sorted pieces with rr_i rc_i <= N^2: R = pieces <= i (p <= i rr_i <= |I| rr_i), C = rest
 (q <= (|I|-i) N^2/rr_{i+1}); pq >= 3N^2/16 => rr_{i+1}/rr_i <= 16|I|^2/3.  OK (improvable to 9|I|^2/2, or 9 i(|I|-i)/2).
 Arithmetic check t2_arith.py: min pq/N^2 under Theorem 2's inequalities, N <= 200: 2/9 >= 3/16.

## Part 1 (g19a) — see below
Setting (DF, g14b/g18s): multigraph, M wires, each with a row-end bit and a column-end bit on its two end vertices (loops allowed),
I = prod_e [Y_e = X_e]; direction = partition of the vertices (empty groups harmless); a term = product over the groups of one
direction; identity mod 2^p with arbitrary kappa.  A'' = max over oriented wire subsets x in {-1,0,1}^E of
min_i (1/2) sum_{g in pi_i} |y(g)|, and log2 N >= A'' (Theorem A'', Lean-verified per the brief).

## 1a. Cut-conditioning  — CONFIRMED [P + E]
 [Y_e = X_e] = sum_{z_e} [Y_e = z_e][X_e = z_e] for each wire cut by i; each factor lives in one group; uncut wires lie in one
 group.  Exact over Z, kappa = 1, 2^{M_i} terms, all of direction i: a DF identity for every p.  Hence N_min <= 2^{min_i M_i}.
 p1_tests.py: full-tensor check on 30 random instances (loops, empty groups, r up to 4): 0 mismatching entries.
 Lean candidate (good, fits DFRank.lean): for every partition lab : V -> G,
   ∀ Y X, (∏ e, if Y e = X e then 1 else 0 : ℤ) = ∑ z : {e // cut e} → Bool, ∏ g, τ lab z g Y X
 with τ depending only on the end bits located in g.

## 1b. Separator / expander bound  — CONFIRMED [P + E]
 (i) sum_i M_i = sum_e d_e (double counting).  (ii) fraction >= 1-eps with d_e <= 1, d_e <= u otherwise =>
 min_i M_i <= M((1-eps)/u + eps) = M(1/u + eps(1-1/u)) <= M(1/u + eps).  (iii) phi = min_{0<vol S<=M} cut S/vol S,
 beta = max_g vol(g)/2M: 2M_i = sum_g cut(g); if all vol(g) <= M: M_i >= phi M; else the unique g0 with vol > M gives
 cut(g0) >= phi(2M - vol g0) and the rest (each vol < M) give >= phi(2M - vol g0) => M_i >= 2 phi (1-beta) M.
 So M_i >= phi' M, phi' = phi min(1, 2(1-beta)).  (iv) consequences u <= 1/phi' and fraction <= 1 - (phi'u-1)/(u-1): algebra OK.
 (v) balanced bipartition from merged groups, both sides >= min(1-beta,1/3) 2M: largest-first argument OK.
 p1_cond.py: exact conductance (all subsets, n <= 13, random multigraphs with loops) x 40 random skewed partitions (empty groups,
 beta up to 1): 12000 tests, 0 violations, min M_i/bound = 1.000 (the bound is attained).  Separator claim: 0/12000 violations.
 Scope remark: "costs >= 2^{phi'M}" is about the cut-conditioning identity only (as the notes say), not about N_min.
 Lean: (i)-(ii) easy (Finset.sum_comm); (iii) medium (pure finite combinatorics, good candidate):
   theorem dir_cut_ge (hφ : ∀ S, 0 < vol S → vol S ≤ M → φ * vol S ≤ cut S) (hβ : ∀ g, vol (lab⁻¹' {g}) ≤ β * (2*M)) :
     φ * M * min 1 (2*(1-β)) ≤ (wires.filter (fun e => lab e.1 ≠ lab e.2)).card

## 1c. Lemma Q / Theorem E  — CONFIRMED; full argument [P]; no gap
 Statement (r = 2, i.e. every direction has <= 2 groups): for every multigraph and directions pi_1..pi_u there is
 x in {-1,0,1}^E with b_i(x) >= M_i/u - u for EVERY i (stronger than the min form); hence log2 N >= min_i M_i/u - u.
 General r (any number of groups, empty allowed): b_i >= M_i/(2u) - u.
 Proof (re-derived).  delta_e(i) = 1[head in H_i] - 1[tail in H_i] in {-1,0,1}; b_i(x) = |sum_e x_e delta_e(i)| (y(V) = 0).
 s in {+-1}^u a global (hence local) maximiser of Q(s) = sum_e <s,delta_e>^2.  Flip of s_i: Delta Q = sum_{e cut by i}
 (4 - 4 s_i delta_e(i) <s,delta_e>) <= 0 => sum_e s_i delta_e(i) <s,delta_e> >= M_i.  x_e = <s,delta_e>/u in [-1,1]
 (|<s,delta_e>| <= d_e <= u) => s_i sum_e x_e delta_e(i) >= M_i/u.  The polytope P = {x in [-1,1]^E : these u inequalities}
 is non-empty and bounded, so has a vertex; a vertex has E independent tight constraints, at most u non-box => at most u
 coordinates not in {+-1}; zero them: each direction loses < u.  General r: signs sigma_i on groups with >= M_i/2 wires across
 (average), delta_e(i) = (sigma(head)-sigma(tail))/2; b_i = (1/2) sum_g |y(g)| >= (1/2)(|y(G+)|+|y(G-)|) = |y(G+)|; apply the
 r = 2 case with M_i' >= M_i/2.  Theorem E = this + (1b)(iii) + Theorem A''.  Every step checks; no hypothesis missing
 (beta-balance only enters through phi').
 Tests: p1_tests.py brute force A'' (3^M, 600 instances, u <= 4, r in {2,3,4,5}, skewed/empty groups): 0 violations of either
 bound.  p1_lemmaQ.py runs the CONSTRUCTION at scale (400 instances, M <= 200, u <= 8, r in {2,3,4,6}): fractional inequality
 always holds, every HiGHS vertex has <= u fractional coords, rounded integer x satisfies b_i >= M_i'/u - u per direction: 0 failures.
 Lean: the fractional half (local max of Q => s_i sum_e <s,delta_e> delta_e(i) >= M_i) is a clean finite statement (good);
 the vertex/rounding step needs LP-vertex theory (harder; could be replaced by an iterative rounding proof).
   theorem lemmaQ_frac (δ : E → Fin u → ℤ) (hδ : ∀ e i, δ e i ∈ ({-1,0,1} : Set ℤ)) (s : Fin u → ℤ) (hs : ∀ i, s i = 1 ∨ s i = -1)
     (hmax : ∀ s', (∀ i, s' i = 1 ∨ s' i = -1) → Q δ s' ≤ Q δ s) (i : Fin u) :
     (Finset.univ.filter (fun e => δ e i ≠ 0)).card ≤ ∑ e, s i * δ e i * ∑ j, s j * δ e j

## 1d. Prop C (circular labelling)  — CONFIRMED [P + E]
 Exact statement: classes K_0..K_{2u-1} (any sizes, may be empty), directions = bipartitions {H_i, V\H_i}, H_i = K_i u..u K_{i+u-1},
 0 <= i < u.  For every network and every x in [-1,1]^E (fractional included):
   min_i b_i <= (1/2) max_{0<=j<u} (cut K_j + cut K_{j+u}) <= (1/2) max_j (vol K_j + vol K_{j+u}),
 = M/u when all classes have volume M/u.  Proof: F_i = y(H_i), F_u = -F_0, consecutive opposite signs, min <= half the jump,
 |y(K)| <= cut(K).  Correct (the j with j+1 = u uses |F_u| = |F_0|).
 Tests: p1_circ.py brute force, 300 random small instances (unequal/empty classes): 0 violations, 158 equalities.
 Exact MILP (integer A'', not LP) on random 4-regular n = 48, M = 96, equal classes: u = 4, 6, 8, 12 -> A'' = 24, 16, 12, 8 =
 M/u exactly (two graphs for u = 12); u = 3 -> 31 (= M/u - 1); min_i M_i in [36, 44].  Numerical claim reproduced.
 Lean (very good candidate, elementary):
   theorem circ_ivt (F : ℕ → ℤ) (u : ℕ) (hu : 0 < u) (hF : F u = -F 0) :
     ∃ j < u, 2 * min |F j| |F (j+1)| ≤ |F (j+1) - F j|
 plus the bookkeeping |y(K)| ≤ cut K.

## 1e. Exchangeable loads  — CONFIRMED (exact, r = 2, exactly exchangeable loads) [P + E]
 Hypothesis made precise: for each d, the wires at distance d are spread EXACTLY uniformly over all unordered label pairs at
 distance d (so the count of each sign type is equal).  Then A'' = (1/u) sum_e |<1,delta_e>| = (M/u) sum_d p_d E|S_d| exactly
 (the >= construction is integral; the <= uses mean over i and invariance of the law of <s,delta> under sign flips).
 p1_exch.py: exact MILP A'' = formula in 17/17 instances (u = 2,3,4, M up to 320).  E|S_d| = d C(d-1,floor((d-1)/2))/2^{d-1}
 exact for d < 60; sqrt(d/2) <= E|S_d| <= sqrt(d) for d <= 2000 (Szarek's constant; equality at d = 2); E|S_d|/d is only
 NON-increasing (equal at d = 2k, 2k+1) — harmless, unused.  Bounds: upper Jensen sum p_d sqrt d <= sqrt(mu); lower
 sqrt(d) >= d/sqrt(u): 0.707 m M/sqrt u <= A'' <= sqrt(m) M/sqrt u.  OK (0.707 = 1/sqrt 2, rounded down).
 Scope: needs exact exchangeability; for approximately uniform loads add the deviation (as in Cor. A'1).

## 2e. Lean candidate (best of Part 2): the four explicit decompositions — integer identities with explicit rank
 factorisations, checkable by `decide`/`norm_num` on Fin 16..36 matrices, valid over every commutative ring:
   theorem cross_4_5_3 : ∃ (U : Matrix (Fin 4 × Fin 4) (Fin 5) ℤ) (V : Matrix (Fin 5) (Fin 4 × Fin 4) ℤ)
       (G G' : Fin 3 → Matrix (Fin 4) (Fin 4) ℤ),
     ∀ x y u v, (if x = v ∧ y = u then 1 else 0) = (U * V) (x,y) (u,v) + ∑ l, G l x u * G' l y v
 (witness: G, G' from exact_4_5_3.json; U, V any integer rank factorisation of R; R's Smith form should be checked to be
 unimodular-friendly, otherwise state over ℚ).
 snf.py: all nonzero invariant factors of R are 1 for the four files => R = U V with integer U (N^2 x p), V (p x N^2):
 the decompositions are identities over Z, hence over every commutative ring; the Lean statement above can be stated over ℤ.

## Summary of verdicts
 1a CONFIRMED; 1b CONFIRMED (bound attained); 1c CONFIRMED (proof complete; stronger per-direction form); 1d CONFIRMED
 (classes of volume M/u; general form with cut(K_j)+cut(K_{j+u})); 1e CONFIRMED under exact exchangeability.
 2a CONFIRMED (char != 2; also GF(2) N<=4, GF(4) N=3 exhaustively); 2b CONFIRMED (char != 2; all chars for N >= 4 via Thm C(a));
 2c PARTLY: literal statement false at q = 1 (N=2, k=3); true with k <= N; rectangular charpoly now machine-checked;
 2d CONFIRMED (given the window-dimension accounting); 2e CONFIRMED over Z (all characteristics).
Files: e_verify.py, snf.py, t3_exh.py, gf4_q1.py, t5_rect.py, t5_exact.py, t2_arith.py, brank.py, p1_lib.py, p1_tests.py,
 p1_lemmaQ.py, p1_cond.py, p1_circ.py, p1_circ12.py, p1_exch.py.
