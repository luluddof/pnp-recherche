# G11c: per-point errors in the GENERAL (non-commutative) class

## 1. Exact formulation
J(t) in F[t]^{2^l x 2^l}, deg <= d, J(0) = X^{(x)l} (x -> ~x). Points tau_k (distinct, nonzero), direction c(k).
E_k = J(tau_k) - K_k, K_k weakly c(k)-nice (K[x,y]=0 if x_c != y_c).
(1a) MIN-RANK COMPLETION (proved; already used in G10b's code): with B01, B10 the two off-nice blocks
     (x_c=0,y_c=1), (x_c=1,y_c=0) of J(tau_k):   min_K rank(J(tau_k)-K) = max(rank B01, rank B10) =: r_k.
     (>=: submatrices. <=: r=max; B01 = P S, B10 = Q R with r columns; E := [P;Q][R,S].)  NOT rank B01 + rank B10.
(1b) Principal-submatrix form: G(t) := J(t) X^{(x)l}, G(0) = I, deg <= d. Then
     r_k = max_{eps in {0,1}} rank G(tau_k)[F,F],  F = facet {x_c = eps}  (principal submatrix).
     d = 1: G = I - tN, and r_k <= r  <=>  beta_k = 1/tau_k is an eigenvalue of BOTH facet compressions
     N[F_{k,0}], N[F_{k,1}] with geometric multiplicity >= 2^{l-1} - r.
     Commutative class = G is a Z_2^l-circulant (G[x,x'] depends on x xor x'); then r_k = #bad k-edges.
(1c) Entrywise form (used for numerics): given the blocks E_k (off-nice parts), J exists iff for every (x,y),
     D = {k : x_c(k) != y_c(k)}, the data (0,[y=~x]), (tau_k, E_k[x,y])_{k in D} lie on a polynomial of degree <= d
     (|D|-d linear conditions when |D| > d).  Unknowns: 2 low-rank blocks per point.
R_gen(l,d) := min over J, tau of max_k r_k / 2^l ;  T_gen(l,d) := min of sum_k r_k (one point per direction).

## 2. Trace Lemma (PROVED, any field)  -- matrix form of "Psi_Q is additive"
alpha_x := [t^1] J(t)[x,~x].  For T subset [l], sigma in {0,1}^T: Q = {x : x_T = sigma}, M_Q(t) = J(t)[Q,~Q]
(= G(t)[Q,Q] up to a permutation), D_Q(t) = det M_Q(t)/det M_Q(0), deg <= d|Q|, D_Q(0) = 1.
  D_Q(t) = 1 + (sum_{x in Q} alpha_x) t + O(t^2)        [det(P + tA) = det P (1 + t tr(P^{-1}A) + ...), P^{-1} = P^T]
and (t - tau_k)^{n_k(Q)} | D_Q for c(k) in T, n_k(Q) = nullity M_Q(tau_k) >= |Q| - r_k   (G6d, Lean PerPoint).
Hence if Q is SATURATED (sum_{c(k) in T} n_k(Q) = d|Q|):    sum_{x in Q} alpha_x = - sum_k n_k(Q)/tau_k.
The left side is ADDITIVE over any splitting of Q into parallel subcubes; only this first-order coefficient is
additive in the non-commutative class (the Schur-complement correction is O(t^2)).
Checked numerically on non-commutative solutions (tracecheck.py).

## 3. Theorems (one point per direction, char 0 -- char != 2 suffices for Thm A)
THM B' (rigidity of the determinant lemma).  q := 2^{l-d-1}. Let T, T' be (d+1)-sets with |T cap T'| = d,
T = U+b, T' = U+c, both with equality in the determinant lemma: sum_T r_k = q = sum_T' r_k.
Then r_k = 0 on U and r_b = r_c = q.
Proof. Every block Q(T,sigma) has |Q| = q, sum_T n_k >= (d+1)q - q = dq, so it is saturated with n_k = q - r_k.
Same for T'. r_b = r_c. The block Q(U,sigma_U) is the union of 2 T-blocks and of 2 T'-blocks:
  2 sum_{k in U}(q-r_k)beta_k + 2(q-r_b)beta_b = 2 sum_U (...) + 2(q-r_c)beta_c  =>  (q-r_b)(beta_b-beta_c) = 0. QED
THM A (general Q_{d+2} lemma).  l = d+2:  sum_k r_k >= 4, for EVERY J (non-commutative), all distinct nonzero tau.
Proof. Determinant lemma on the (d+1)-sets: every such sum >= 2. Total <= 3 forces r in {0,1}^l with exactly
three 1's a,b,c; T = [l]-c, T' = [l]-b are both saturated (sum = 2 = q) and r_a = 1 != 0 contradicts Thm B'.
(Direct: 2x2 blocks {x, x+e_c}: det = -(1-t/tau_a)(1-t/tau_b) prod_R (1-t/tau_k)^2, so
 alpha_x + alpha_{x+e_c} = -(2 rho_R + beta_a + beta_b); same with a<->c; going round the (a,c)-square gives
 2(beta_a - beta_c) = 0.)  QED.   Tight: T_gen(d+2,d) = 4 (commutative examples).
COR 1 (any l >= d+2, |S| = d+2):  sum_{k in S} r_k >= min( 2^{l-d}, ceil(2^{l-d-1}(d+2)/(d+1)) + 1 )
   (at most one of the d+2 subsets S-j is saturated unless the pattern is (0,..,0,q,q)).
   By restriction to a pair of subcubes also sum_S r_k >= 4, un-normalised.
COR 2 (d = 1): r_a >= 1 => 2 r_a + r_b + r_c >= 2^{l-1} + 1.  (4,1): T_gen >= 10, and (2,2,2,2), (1,3,3,3)
   are impossible, so max_k r_k >= 3.  Cor 1 is tight at (4,1): min(8, 7) = 7 is attained by the triple (2,2,3) below.

## 4. What FAILS in general (numerical counterexamples, feas.py / sol41.py / refine41.py)
tau = (1, 1.7, 2.9, 4.3), l = 4, d = 1 (Levenberg-Marquardt on formulation (1c)):
 * profile (2,3,3,3): solution found, total 11 < 12 = commutative optimum T(4,1) (G8d, generic tau).
   Independent reconstruction J = X + tA1: singular values of the 8 off-nice blocks have a gap 1e0 / 1e-15.
   refine41.py: residuals evaluated in EXACT rational arithmetic go 1e-14, 6e-29, 6e-42, 8e-55, 1e-67, 2e-80 under
   float Gauss-Newton steps.  N = -A1 X is neither symmetric nor a circulant (genuinely non-commutative).
 * profile (3,3,3,3): feasible => max_k r_k = 3, R_gen(4,1) = 3/16 (lower bound PROVED by Cor 2).
   Commutative class, same tau (comm41x.py, exhaustive with pruning): (3,3,3,3), (2,3,3,3), (2,2,3,4) infeasible;
   (0,4,4,4), (2,2,4,4), (1,3,4,4), (2,3,3,4) feasible: commutative min-max = 4 (rho = 1/4), total 12.
 * profile (2,2,3,4): feasible; the triple S = {1,2,3} has sum 7 < 8 = 2^{l-d}.
   So the normalised corollary of the commutative Q_{d+2} lemma ("sum_{k in S} bad_k >= 2^{l-d}", which uses
   additivity of bad edges over disjoint subcubes) is FALSE per (d+2)-set in the general class.
 * (2,2,3,3) (total 10): not found in 36 starts (cost ~1e-7 with growing parameters) -> T_gen(4,1) in {10, 11}.
 * (3,1): (1,1,1) not attained but costs -> 0 with parameters -> infinity (border phenomenon; irrelevant for exact J).

| (l,d) | commutative T (G8d) | general: proved LB | general: found UB | max_k r_k comm. / general |
|-------|----|----|----|----|
| (2,1) | 1  | 1 (det lemma) | 1 | 1 / 1 |
| (3,2) | 1  | 1 | 1 | 1 / 1 |
| (3,1) | 4  | 4 (Thm A) | 4: (0,2,2), (1,1,2) | 2 / 2 |
| (4,2) | 4  | 4 (Thm A) | 4: (1,1,1,1) | 1 / 1 (LB 1) |
| (4,1) | 12 (generic tau; all-tau LB 11) | 10 (Cor 2) | 11: (2,3,3,3), (2,2,3,4) [numerical] | 4 / 3 (general: LB proved, UB numerical) |
| (5,2) | 12..14 | 9 (det lemma with integrality gives 8; Thm B' kills (0,2,2,2,2), (1,1,2,2,2)) | not determined ((2,2,2,2,2), (2,2,2,2,3): 3 starts, not found) | ? |

## 5. Conclusion
* The commutative class is NOT optimal: first separation at (l,d) = (4,1) (total 11 vs 12, max 3 vs 4).
* What holds in general: determinant lemma (G6d), min-rank = max of the two blocks, Trace Lemma, Thm B', Thm A
  (Q_{d+2} lemma at l = d+2, same constant 4), Cor 1 (+1 over the determinant lemma per (d+2)-set).
* What does not transfer: the scaling 2^{l-d} per (d+2)-set, hence rho_avg >= 2^{-d}/(d+2) is UNPROVED in general
  (not refuted: total 11 >= 10.67 at (4,1)); general proved bound stays rho_max >= 2^{-d-1}/(d+1) (+1 unit).

## 6. Self-review
* Thm A, B', Trace Lemma, (1a): complete proofs; they use only G6d's nullity <= order lemma (Lean-verified).
  Thm B' needs char 0 (or char not dividing q - r_b); Thm A needs char != 2. NOT formalised in Lean.
* All upper bounds for (4,1) in the general class are NUMERICAL (floating point + exact-residual refinement to
  1e-80); no exact rational/algebraic certificate. The Jacobian is rank deficient at the solution (redundant
  constraints), so no Newton-Kantorovich certificate was obtained. Treat "T_gen(4,1) <= 11" as strongly supported,
  not proved.
* "not found" entries are heuristic (few random starts); only the proved LB column is rigorous.
* Commutative comparisons at (4,1) are for the same tau (which is generic for l = 4); for special tau the
  commutative optimum could be 11 (G8d's all-tau bound).
* Only one point per direction was treated in Thm A/B'; with several points per direction the Trace Lemma still
  holds, the saturation bookkeeping changes.
* Strong niceness (K_k of the special form the algorithm uses) can only increase the ranks: lower bounds stay valid,
  the numerical upper bounds are for WEAK niceness only.
