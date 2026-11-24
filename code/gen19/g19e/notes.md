# G19e — adversarial review + novelty check of G18g Theorem 1/2 and G18f Theorem F

Status tags: [P] re-derived on paper by me (same model family as the authors: P+ at best), [E] exact computation (my code),
[L] literature (what I actually read is said), [O] open.  (work in progress: saved as I go)

## 0. Definitions (re-stated)
K any field. E, G, W finite-dimensional, dim E = m, dim G = n, dim W = Q.
P : E (x) G (x) W -> G (x) E (x) W,  e (x) g (x) w |-> g (x) e (x) w     ("swap with a spectator").
Bipartition:  site 0 = (input E, output G) | rest = (input G (x) W, output E (x) W).
 * rank(A) = ordinary rank of the linear map A (matrix inputs | outputs).            [= "row rank", "rr"]
 * osr(B) = operator-Schmidt rank across site 0 | rest = least b with B = sum_{k<=b} L_k (x) M_k,
   L_k in Hom(E,G), M_k in Hom(G(x)W, E(x)W) = rank of the REALIGNED matrix (in_0,out_0)|(in_rest,out_rest). [= "rc", "column-cut rank"]
Both are ranks of flattenings of the same 4-leg tensor, hence subadditive and lower semicontinuous.

Dictionary.
 (I)  G18g crossing X(x,y,u,v) = [x=v][y=u]: operator (u,v) -> (x,y) is the flip F on K^N (x) K^N; rr = rank (x,y)|(u,v);
      rc = rank (x,u)|(y,v) = osr(F-part) across factor 1 | factor 2.  This is m = n = N, Q = 1.
      NOTE: it is the REALIGNMENT, not the partial transpose: F^Gamma = |Omega><Omega| has rank 1 and I^Gamma = I.
      Equivalent form: identity operator I = A + B with rank(A) <= p and rank(B^Gamma) <= q (apply F on one side).
 (II) G18f 2 x h grid, column j = site 0: E = strands leaving column j (dim M = N^{h-1}), G = strands entering column j
      (dim M), W = the (h-1)(h-2) strands not touching column j (dim Q = N^{(h-1)(h-2)}); P = F_{E,G} (x) (a permutation of W).
      r = rank(A), b = osr(B).  h = 2: Q = 1, M = N: exactly (I).

## 1. Common theorem (implies both) [P]
THEOREM C. If P = A + B with rank(A) <= r and osr(B) <= b then, for all integers k,
  (a)  r   >= k Q (n - k b)        0 <= k <= m          (and also r >= k Q (m - k b), 0 <= k <= n)
  (b)  b Q >= k (n Q - k r)        0 <= k <= m          (and also b Q >= k (m Q - k r), 0 <= k <= n)
Proof of (a). Pick a surjection pi : E -> K^k. S := intersection over l of ker((pi (x) 1_W) M_l) in G (x) W, codim <= k Q b.
 Compress: Y := (1_G (x) pi (x) 1_W) o (.) o (1_E (x) iota_S).  Y(B) = sum L_l (x) ((pi(x)1) M_l iota_S) = 0.
 Y(P) : e (x) sigma |-> (pi e) (x) sigma (up to order of factors): image = K^k (x) S, rank = k dim S.  Y(A) has rank <= r.
 So r >= k dim S >= k (nQ - kQb).                                                              QED
Proof of (b). Pick Phi in E, dim k.  Psi := {psi in G (x) W : Phi (x) psi in ker A}, codim <= k r, so dim Psi >= nQ - kr.
 On Phi (x) Psi: P = B, so osr(P restricted to Phi (x) Psi) <= b.  But P|: e (x) psi -> ... has
 osr = k * dim(G-support of Psi) >= k * dim Psi / Q   (G-support = least G' with Psi in G' (x) W).  QED
The variants in parentheses: same proofs compressing on site 0 instead of the rest (resp. Phi in G* on the output side).
Only linear algebra; no Flanders; any field; no hypothesis on char.

Border version [P]: over K((eps)) with P = A_eps + B_eps + O(eps) (A_eps, B_eps may have poles): in (a) pi is fixed and S_eps
 is a subspace of K((eps))^{nQ}; its saturated lattice reduces mod eps to a subspace S_0 of the same dimension (Grassmannian is
 proper); Y_{S}(A_eps) = Y_S(P) + O(eps) has rank <= r over K((eps)), so its reduction Y_{S_0}(P) has rank <= r; but
 rank Y_{S_0}(P) = k dim S_0 for EVERY S_0. Same for (b) (osr is a matrix rank, the lower bound holds for every Psi_0).
 Over R or C: same with compactness.  So the border case is fine.


STRENGTHENING found during this review (not in G18f/G18g) [P, E]:
  (c)  EXACT case only:  b < m n  (osr(B) not full)   ==>   r >= Q.
Proof. osr(B) < dim Hom(E,G) gives a nonzero functional c on Hom(E,G) with c(L_k) = 0 for all k; view c as a map
 c^ : G -> E of rank rho >= 1.  Contract site 0 with c:  (c (x) 1)P = c^ (x) 1_W, of rank rho Q; (c (x) 1)B = 0;
 (c (x) 1)A = sum_i sum_{t<=rho} (alpha_t x_i)(beta_t y_i)^T has rank <= rho r (A = sum_{i<=r} x_i y_i^T, c = sum_t alpha_t beta_t^T).
 So rho Q <= rho r.  QED.   Tight: A = E_{g0 e0} (x) P_{g0 e0} has rank Q and osr(P - A) = mn - 1 (test_strong.py).
 G18f only had r >= Q (1 - b/M^2).  For the hexagon N=2: b_j <= 15 forces r >= 4 (their table says r >= 1 for b = 12..15).
 Consequence: r b >= (2/9) Q M^2 for EVERY 1 <= b < M^2 (not only b <= M^2/2); "r < Q" forces osr_j(B) = M^2 EXACTLY at every column.
 (c) is FALSE in the border sense: m = n = Q = 2, r = 1 < Q, b = 3 < 4 is a border decomposition (border_gap.py, symbolic):
   A_eps = x y^T, x = g_1(x)e_1(x)w_1 + eps^-1 g_2(x)e_1(x)w_2, y = e_1(x)g_1(x)w_1 + e_2(x)g_1(x)w_2, c_eps = diag(1, eps).
 So "exact or border" in Theorem F is right for (a),(b) but the two notions really differ here.

Arithmetic corollary [P, E: arith.py, test_mech.py T6]. m = n = M.  From (a),(b),(c): r b >= (2/9) Q M^2 for 1 <= b < M^2.
 For the crossing (Q = 1): p q >= (2/9) N^2 whenever p, q >= 1.   2/9 is the EXACT infimum of the inequality system
 (attained at (p,q) = (2N/3, N/3) when 3 | N); G18g's 3/16 comes from a lossy last step (N^2/4 - q^2/4 at q = N/2 instead of
 the vertex analysis).  From (a),(b) alone (border case) the product bound needs b <= M^2/2, as G18f states.

## 2. Exact tests [E]  (my code: ff.py, cross_q1.py, q2n3.c, grid.py, test_grid.py, test_mech.py, test_strong.py, border_gap.py, arith.py, attack_hex.py)
 * cross_q1.py: exhaustive min rank(F - G(x)G') (crossing, q = 1; G in rank normal form, G' all matrices):
   N=2 over F_2,F_3,F_5,F_7: min p = 1 (= bound).  N=3 over F_2, F_3, F_5: min p = 3 (bound 2).  N=4 over F_2: min p = 6 (bound 4).
   Always N(N-1)/2 (attained only with G invertible).  So (2,1) is NOT feasible for N=3 over F_2, F_3, F_5.
 * q2n3.c: exhaustive, N=3 over F_2, all C with rc(C) <= 2 (4*10^8 cases): min p = 3.  Hence over F_2: (2,2), (1,2), (2,1) infeasible
   and P(3) = 4 exactly over F_2 (G18g had this only numerically over R, C).  out_q2n3.txt.
 * Theorem 2/F are therefore tight at N=2 only; at N=3 the true minimum of pq is 3 = N^2/3 over F_2 (bound: 2).
   No known point has pq < N^2/4 (the "real s" version r >= max_s s(M - bs) = M^2/4b is neither proved nor refuted; N=3 does not
   refute it; first undecided case of the gap: (2,2) at N=4 is allowed by both, (4,2) at N=6 by the integer version only).
 * test_grid.py: P built two ways (permutation of strands / product of rectangle swaps); 5928 decompositions P = (P-B)+B with B of
   small osr (random / permutation / sparse / pieces of P; N=2,3,4 for h=2 and N=2 for h=3; mod 2,3,5) and 2502 decompositions
   P = A+(P-A): no violation.  Hexagon N=2 product point (A = Y(x)Y(x)Y, rank 1): osr_j(P-A) = 16 at every column.
   rank(P-1) = 28 (hexagon N=2), 3 (N=3 crossing), 6 (N=4 crossing), in char 2, 3 and 101.
 * test_mech.py: T4 the two mechanism identities with m != n and Q > 1 (360 random cases mod 2,3,7); T5 3366 more decompositions with
   general (m,n,Q), all four families (a),(b),(a'),(b'); T6 integer optimum of the system = 2/9 (N <= 240).
 * test_strong.py: 2592 matrices A of rank < Q all leave osr(P-A) = mn; tight example checked.
 * attack_hex.py (local search over F_2, evidence only): see out_attack.txt; nothing below the bounds.
 Small characteristic: every test above was run mod 2 and mod 3; the proofs use no division.

## 3. Verdicts on the two theorems
(I) G18g.
 * Theorem 1 (telescoping) CONFIRMED as two finite-N statements: (i) pieces with rr*rc <= Q' give feasible (|I|s, |I|Q'/s);
   (ii) m nested feasible pairs give m+1 pieces with rr*rc <= (p_j+p_{j-1})(q_j+q_{j-1}).  Corollary 4rho* = sup(theta+T(theta)):
   CONFIRMED as a statement about ranks, with a mesh loss delta for 2/delta+1 pieces; needs T defined uniformly in N
   (the author says so); nothing about bit length or uniformity.  "Determine each other" is up to factors |I|, |I|^2 and 4.
 * Theorem 2 CONFIRMED (any field, both inequalities; the second by the y<->u symmetry of X, which exchanges rr and rc).
   It is Theorem C with Q = 1.  The stated corollary pq >= 3N^2/16 is true and improvable to 2N^2/9.
(II) G18f Theorem F.
 * CONFIRMED: r >= sQ(M - bs) for integer s and for s = 1/rho, any field, exact and border.  Their dual (annihilator) proof is
   correct as written.  Their primal proof (Schmidt ranks + Flanders) is correct, but
     - Flanders is NOT needed: (a) above is a 5-line compression argument; and it removes the restriction "h >= 3" for integer s
       (the restriction is an artefact of max(M, sQ) in Flanders' bound; for h = 2 integer s is G18g's first inequality).
     - Flanders' theorem as quoted is right: a linear space of m x n matrices of rank <= k has dim <= k max(m,n); Flanders (J. London
       Math. Soc. 37 (1962) 10-16) assumed |F| >= k+1; Meshulam (Quart. J. Math. Oxford 36 (1985) 225-229) removed the field
       hypothesis [L: Meshulam's paper via search summary; de Seguins Pazzis arXiv:1004.0298 also].  So the quoted step holds over F_2.
     - non-integer s: NOT claimed and not proved; no counterexample found (it would be r b >= Q M^2/4).
 * "r b >= (2/9) 2^S when b <= M^2/2": CONFIRMED for b >= 1 (b = 0 must be excluded: r = 2^S, product 0).  With (c): true for all 1 <= b < M^2 (exact case).
 * Errata in the corollaries:
     (C1) "1 - O(1/lam)" is 1 - (2 log2 m + log2 4.5)/(2 lam): O(1/lam) only for O(1) terms, O(log lam/lam) for poly(lam) terms,
          o(1) for 2^{o(lam)} terms.  The multi-term step IS proved (rank and osr are subadditive; for the crossing both costs are
          logs of flattening ranks, so the row terms merge and the column terms merge).
     (C2) "any B whose cut at some column is not full": not justified by Theorem F as stated (needs b <= M^2/2), but TRUE by (c).
          "kappa_c >= (h/2)(1 - kappa_r)" must read min(1, (h/2)(1 - kappa_r)) - o(1) (for kappa_r < 1 - 2/h all cuts are full and kappa_c >= 1 - o(1)).
     (C3) the 2x3 corner: CONFIRMED and sharpened: kappa_r + kappa_c < 1 - delta forces r < Q and then osr_j(B) = M^2 exactly at every
          column (not only >= M^2(1 - r/Q)).  The corner itself stays open.  It relies on cost(B_i) >= log2 osr_j(B_i) (true for a K_h network).
 * The withdrawn "two-term reduction" is not used anywhere in Theorem F: only osr_j of the SUM of column terms enters.

## 4. Same theorem?  Yes.
 Crossing = 2 x 2 grid = Theorem C with m = n = N, Q = 1: (a) is G18g's "p >= k(N-kq)" = G18f primal (integer s);
 (b) is G18g's "q >= k(N-kp)" = G18f dual (s = 1/rho).  For h >= 3, Theorem F = Theorem C applied to (r, b) with spectator Q:
 the pair (r/Q, b) satisfies the crossing inequalities of size M.  Theorem C (section 1) is the common statement.
 Sharpness: 2/9 is the exact constant of the inequality system; 3/16 is a lossy corollary; neither is known to be attained by real
 decompositions (N=2: 1/4; N=3 over F_2: 1/3; q=1: (N(N-1)/2, 1), i.e. 1/2).

## 5. Novelty [L]  (about 10 searches + 2 arXiv papers read in part; downloads in lit_dl/, untrusted)
 * The object is exactly Naslund's Problem 12 (E. Naslund, "The partition rank of a tensor and k-right corners in F_q^n",
   arXiv:1701.04475, JCTA 2020; read: section 2): P-rank of delta(x1,x4)delta(x2,x3) for P = {12|34, 13|24}; asks n^{1+eps}.
   T. Karam, "The interplay between bounded ranks of tensors arising from partitions", arXiv:2310.11356 (Question 1.2, read):
   states the known bound r + s >= n, proved by fixing x1 = x4 = 1 - this is the case k = 1 of Theorem 2 - and says no n^{1+eps}
   bound is known for any explicit tensor.  Its Theorem 1.4 (power bounds n^{c(d)}) concerns general partitions, not this trade-off.
 * I did NOT find the product form p >= k(N - kq) / pq >= cN^2, nor "swap = low rank + low operator-Schmidt rank" (quantum
   information: only osr(SWAP) = N^2 and its use as a communication lower bound), nor a rigidity statement for the swap /
   commutation matrix against Kronecker rank.  Hrubes-Wigderson-Yehudayoff (JAMS 2011) and Hrubes (CCC 2024, ECCC TR24-026,
   TR24-099) were not opened (not on arXiv; only search summaries): the balanced lower bound there is n.
 * Assessment: not found in the literature searched; the proof is a 5-line extension of the standard k = 1 argument (restrict to k
   values instead of one, kill the other type on a subspace), so it is very likely folklore-level and at most a remark, not a
   result; it gives nothing on Naslund's question (p = q ~ N/2 is allowed).  A 20-minute search cannot certify novelty.

## 6. Lean 4 + Mathlib statement (NOT formalised, not type-checked)
```lean
import Mathlib
open Matrix
variable (K : Type*) [Field K]

/-- swap with a spectator, E ⊗ G ⊗ W → G ⊗ E ⊗ W; rows = outputs (g,e,w), columns = inputs (e,g,w). -/
def swapW (m n Q : ℕ) : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K :=
  fun o i => if o.1 = i.2.1 ∧ o.2.1 = i.1 ∧ o.2.2 = i.2.2 then 1 else 0

/-- realignment across site 0 | rest: rows (out₀,in₀) = (g,e'), columns (out_rest,in_rest) = ((e,w),(g',w')). -/
def realign {m n Q : ℕ} (T : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K) :
    Matrix (Fin n × Fin m) ((Fin m × Fin Q) × (Fin n × Fin Q)) K :=
  fun a b => T (a.1, b.1.1, b.1.2) (a.2, b.2.1, b.2.2)

theorem swap_rank_osr_tradeoff (m n Q r b : ℕ)
    (A B : Matrix (Fin n × Fin m × Fin Q) (Fin m × Fin n × Fin Q) K)
    (h : A + B = swapW K m n Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b) :
    (∀ k ≤ m, k * Q * n ≤ r + k * k * Q * b ∧ k * n * Q ≤ b * Q + k * k * r) ∧
    (∀ k ≤ n, k * Q * m ≤ r + k * k * Q * b ∧ k * m * Q ≤ b * Q + k * k * r) ∧
    (b < m * n → Q ≤ r)

/-- arithmetic corollary (m = n = M). Crossing = the case Q = 1 (then `hb'` only excludes B = swap). -/
theorem swap_rank_osr_product (M Q r b : ℕ)
    (A B : Matrix (Fin M × Fin M × Fin Q) (Fin M × Fin M × Fin Q) K)
    (h : A + B = swapW K M M Q) (hA : A.rank ≤ r) (hB : (realign K B).rank ≤ b)
    (hb : 1 ≤ b) (hb' : b < M * M) : 2 * Q * M ^ 2 ≤ 9 * (r * b)
```
 (For the crossing in G18g's indexing: X (x,y) (u,v) = [x = v ∧ y = u], rr = X-type `Matrix.rank`, rc = rank of
  `fun (x,u) (y,v) => C (x,y) (u,v)`; conclusion ∀ k ≤ N, k*N ≤ p + k*k*q ∧ k*N ≤ q + k*k*p.)
 Proof ingredients needed in Lean: rank of a submatrix/compression <= rank; rank subadditivity; rank of a Kronecker product
 of a surjection with an injection; dimension of an intersection of kernels.  No Flanders.

## 7. The two consequences
(a) from (I).
 * "rho > 1/2 follows from an unbalanced bound at a single theta": CONFIRMED.  Exact form: if for some theta0 every feasible (p,q)
   with p <= N^theta0 has q >= N^{2-theta0+eps}, then any decomposition into |I| pieces has max rr*rc >= N^{2+eps}/|I|^2, i.e.
   rho >= 1/2 + eps/4 - log|I|/(2 log N).  (Theorem 1(i) with s = N^theta0/|I|.)  Conversely rho* = 1/2 iff T(theta) = 2 - theta for all theta.
 * "from P(N) <= N^c alone only rho <= (2+c)/4": CONFIRMED as stated: the two-piece decomposition R, C gives 4rho <= 2 + c, and for
   two pieces max rr*rc >= N^2 P(N)/4 exactly as argued.  "Only" = no better generic derivation is known (restriction to smaller N
   and tensoring with trivial identities give nothing better); it is not a theorem that none exists.
 * Caveat that matters for the programme: Theorem 2/C can never give more than the line T = 2 - theta (it allows p = q ~ N/2), so
   it is the floor rho >= 1/2 again, in product form; it does not touch the open unbalanced problem.
(b) from (II).  "one crossing: kappa_r + kappa_c >= 1 - O(1/lam) for every expansion": PARTLY.
   Proved for any number m of typed terms with loss (2 log2 m + 2.17)/(2 lam) (sum the row terms into R and the column terms into C;
   both costs are logs of matrix ranks, hence subadditive; c != 0 in K).  So: 1 - o(1) for 2^{o(lam)} terms, 1 - O(1/lam) only
   for O(1) terms.  Edge cases (no row term or no column term) are fine.  In G18g's cost model (each piece pays rr*rc) the same
   inequality only gives rho >= 1/2 - o(1).

## 8. Adversarial self-review of this review
 * Same model family as the authors: the re-derivations are P+ at best.  Theorem C (a),(b),(c) were each tested exactly on
   thousands of decompositions, but random/structured decompositions are far from extremal; only q = 1 (N <= 4) and q <= 2
   (N = 3, F_2) are exhaustive.  attack_hex.py is a weak optimiser (for b <= 2 it did worse than b = 1): evidence only.
 * The border proof of (a),(b) (valuation ring / properness of the Grassmannian) is written, not tested beyond one symbolic example.
 * Novelty: HWY 2011 and Hrubes 2024 were not opened; arxiv.org search pages returned nothing parseable.
 * The Lean statements were not type-checked.
 * I did not re-check G18f section 2 (no amplification), sections 4-6, nor G18g sections 1bis, 2, 3 (tables), nor G18h.
