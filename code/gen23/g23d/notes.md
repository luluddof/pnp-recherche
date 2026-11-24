# G23d — proving the q = 1 law for the captured mass c(N;p,q) (working notes, saved as I go)
Tags: [P] proof here (unreviewed), [N] numerical, [E] exact.
Conventions: X[x,y,u,v] = <xy|X|uv>; F|uv> = |vu>; B = G(x)G' means B[x,y,u,v] = G[x,u]G'[y,v]; <X,Y> = tr(X^*Y).
tr(F(G(x)G')) = tr(G G')  [E: F(G(x)G')|uv> = G'|v>(x)... ; sum_{x,y} G[x,y]G'[y,x]].

## 0. Reformulations [P]
(R0) Closure is free for UPPER bounds on c.  If |<F,X>|^2 <= kappa ||X||^2 for all X in sigma_p + tau_q, the same holds on the
 closure (both sides continuous), and then ||F - X||^2 = N^2 - 2Re<F,X> + ||X||^2 >= N^2 - kappa for all X in the closure
 (min over t>=0 of N^2 - 2 sqrt(kappa) t + t^2).  So c(N;p,q) = sup_{X in sigma_p+tau_q, X != 0} |<F,X>|^2/||X||^2, and border
 sequences need no separate treatment for lower bounds on d^2.
(R1) One-sided split (exact formula for q = 1).  For a rank-p orthogonal projection P (on V(x)V) put
     kappa(P) := sup_{B product, B(1-P) != 0} |tr(F B (1-P))|^2 / ||B(1-P)||^2 ,   kappa_p := sup_P kappa(P).
 Then c(N;p,1) = p + kappa_p.
 Proof (<=): X = A + B, rank A <= p, P := projection onto (ker A)^perp (rank <= p; enlarge to rank p, kappa(P) only needs
 X(1-P) = B(1-P)).  tr(FX) = tr(FXP) + tr(FB(1-P)).  |tr(FXP)| = |<FP, XP>|... <= ||FP||_F ||XP|| = sqrt(p) ||XP||;
 |tr(FB(1-P))| <= sqrt(kappa_p) ||X(1-P)||.  Cauchy-Schwarz: |tr FX| <= sqrt(p + kappa_p) sqrt(||XP||^2 + ||X(1-P)||^2) = ... ||X||.
 (>=): given P, B, take X = t B + A with A := (F - tB)P scaled... direct: d^2 <= ||(F - tB)(1-P)||^2 (A = (F - tB)P has rank <= p);
 = (N^2 - p) - 2t Re tr(FB(1-P)) + t^2 ||B(1-P)||^2 -> N^2 - p - kappa(P) after optimising t (phase absorbed in B).  QED
 So the whole q = 1 problem is the one-sided quantity kappa_p, and the LAW c = m^2 + 2j (p = C(m,2) + j) is equivalent to
     kappa_p = p + m,   m = m(p) := max{m : m(m-1)/2 <= p}            (c(p,1) = 2p + m(p)).
 Extremal for j = 0: P = P_anti(W), B = P_W (x) P_W, dim W = m: tr(FB(1-P)) = m + C(m,2), ||B(1-P)||^2 = m^2 - C(m,2) => ratio C(m+1,2) = p + m.

## 1. Numerics on the one-sided quantity [N] (kap.py, kapD.py)
 kappa_p = sup_P kappa(P): N=3: p=1 3.0000000 (6/6), p=2 3.99999 (border), p=3 6.0000000; N=4: p=1 3.0000000, p=2 3.99998,
 p=4 6.998 (border).  All = p + m(p).  => c(p,1) = 2p + m(p) confirmed through the exact formula (R1).
 kappa(|b><b|) depends only on the Schmidt coefficients sigma of b (local-unitary symmetry, notes sec. 2):
   product b: 2 (sup, not attained);  Schmidt rank 2 (ANY sigma1, sigma2 > 0): 3;  uniform rank k: 2(k+1)^2/(k(2k-1))
   [k=3: 32/15, k=4: 25/14, exact: optimiser H = rank-2 antisymmetric, value (1+1/k)^2/(1-1/(2k))].  Worst case = rank 2.
 Strong "symmetric compression" variant sup |<Om|Y|Om>|^2/||P_S Y P_S||^2 over osr-2 Y: = 9 at N=3 (FALSE as a route).

## 2. Equivalent forms of the q = 1 law [P]
 (EY)  c(p,1) = N^2 - inf_{G,G'} sum_{i>p} sigma_i(F - G(x)G')^2, so  LAW <=>  for all G, G':
        sum_{i<=p} sigma_i(F - G(x)G')^2  <=  2p + m(p) - 2Re tr(GG') + ||G||^2 ||G'||^2.
       p = 1 in matrix form: for all N x N G, H and ||V||_F = 1:   ||V^T - G V H||_F^2 <= 4 - 2 Re tr(G H^T) + ||G||^2 ||H||^2.
 (Theta) with K = (G(x)G')F (a "twisted product"), Theta := K + K^* - K^*K = 1 - (1-K)^*(1-K) <= 1:
        LAW <=> the sum of the N^2 - p largest eigenvalues of Theta is <= p + m(p)  (p=1: tr Theta - lambda_min(Theta) <= 3).
        p = 0 (kappa_0 = 1): tr Theta = 2Re tr K - ||K||^2 <= 1 since |tr K| = |tr GG'| <= ||G|| ||G'|| = ||K||.
 Symmetries: G -> U G W^*, G' -> W G' U^*, b -> (W(x)U) b leave everything invariant => b in Schmidt form sum sigma_i |ii>.
 Singular-term form: with SVDs G = sum s_i x_i y_i^*, G' = sum t_k x'_k y'_k^*, B = sum_{ik} s_i t_k |x_i x'_k><y_i y'_k|
  (orthonormal product bases on both sides), Q = 1 - vv^*:  B Q = sum s_i t_k |w_ik><Q z_ik| has terms with ORTHOGONAL ranges, so
  ||BQ||^2 = sum s_i^2 t_k^2 ||Q z_ik||^2 and tr(FBQ) = sum s_i t_k <Q z_ik | F w_ik>;  each single term has |.| <= its norm.

## 3. THEOREM A (rigorous, all N): c(N;p,1) <= p + t + (sqrt(1+t) + sqrt(p theta))^2/(1-theta), t = floor(p/theta), any 0<theta<1.
 In particular c(N;1,1) <= 11 (theta -> 1/2+, t = 1), i.e. d^2(N;1,1) >= N^2 - 11 for EVERY N (truth: N^2 - 4);
 and c(N;p,1) <= 11.7 p + O(sqrt p) (theta ~ 0.42), i.e. d^2(N;p,1) >= N^2 - O(p) uniformly in N.   [P]
 Compare: every contraction+rank+Eckart-Young bound (Theorem C / class B) gives at best N^2 - Theta(N(p+q)) (G22d sec. 5).
Proof.  X = A + B, rank A <= p, B = G(x)G'.  P := projection on (ker A)^perp (rank <= p), Q = 1 - P; XQ = BQ.
 (1) |tr(FXP)| <= ||XP||_tr <= sqrt(p) ||XP||   (||F||_op = 1, rank XP <= p).
 (2) SVDs G = sum s_i x_i y_i^*, G' = sum t_k x'_k y'_k^*  =>  B = sum_j beta_j |w_j><z_j|, j = (i,k), beta_j = s_i t_k,
     w_j = x_i(x)x'_k, z_j = y_i(x)y'_k: both families ORTHONORMAL (products of orthonormal bases).
 (3) T := { j : ||P z_j||^2 >= theta }.  sum_j ||P z_j||^2 <= tr P <= p  =>  |T| <= p/theta.  B_T := sum_{j in T}, B_R := B - B_T.
     Ranges of B_T Q and B_R Q lie in span{w_j: j in T}, span{w_j: j notin T}: orthogonal, so ||BQ||^2 = ||B_T Q||^2 + ||B_R Q||^2;
     ||B_R P||^2 = sum_{j notin T} beta_j^2 ||P z_j||^2 <= theta ||B_R||^2,  hence ||B_R Q||^2 >= (1-theta) ||B_R||^2.
 (4) |tr(F B_T Q)| <= sqrt|T| ||B_T Q||   (rank <= |T|).
 (5) Realignment R (xu|yv) is an HS isometry with R(F) = F (unitary) and R(B_R) = E M_T E'^* with isometries E, E' and the N x N
     matrix M_T = (s_i t_k [ (i,k) notin T ]) = s t^T - (s t^T) o 1_T, rank <= 1 + |T|.  So |tr(F B_R)| = |<F, R(B_R)>| <=
     ||M_T||_tr <= sqrt(1+|T|) ||B_R||;  |tr(F B_R P)| <= sqrt(p) ||B_R P|| <= sqrt(p theta) ||B_R||.
     => |tr(F B_R Q)| <= gamma ||B_R Q||,  gamma = (sqrt(1+|T|) + sqrt(p theta))/sqrt(1-theta).
 (6) |tr FX| <= sqrt(p)||XP|| + sqrt|T| ||B_T Q|| + gamma ||B_R Q|| <= sqrt(p + |T| + gamma^2) ||X||  (Cauchy-Schwarz,
     ||X||^2 = ||XP||^2 + ||B_T Q||^2 + ||B_R Q||^2).  Extends to the closure by (R0).  p = 1, theta = 1/2 + 0: 1 + 1 + 9 = 11.  QED
 Check [N] checkA.py: every inequality (1)-(6) tested on 3000 random/near-extremal instances (N = 2..4, p = 1..3): 0 violations.
 Where it loses vs the truth (4): (5) bounds tr(F B_R) and tr(F B_R P) separately (triangle inequality), the true mechanism is a
 cancellation between them (the antisymmetric direction).  An attempted sharpening "B minus its top singular term has kappa <= 2"
 (which would give exactly 4) is FALSE numerically: restR.py gives 3.0 at N = 3.

## 4. Status of the sharp statements
 c(N;1,1) = 4:  ">= 4" [P, exact construction: F_W = P_W(x)P_W - 2 P_anti(W), dim W = 2];  "<= 4" NOT proved; proved <= 11.
 Equivalent sharp target (any one suffices) [P equivalences, sec. 0/2]:
   kappa_1 = 3  <=>  ||V^T - G V H||_F^2 <= 4 - 2Re tr(GH^T) + ||G||^2||H||^2 for all G, H, ||V||_F = 1
              <=>  tr Theta - lambda_min(Theta) <= 3 for Theta = K + K^* - K^*K, K = (G(x)G')F.
 Verified [N]: N = 3,4,5 (G22d) and N = 6 (log6_1_1.txt: 4/4 restarts d^2 = 32.000000 = 36 - 4); kappa_1 = 3.0000000 at N = 3, 4.
 Partial check only (NOT a proof of a special case): for G = G' normal, K is normal and F-invariant,
   sigma_1(F - G(x)G)^2 = max(max_{i<=j}|1-g_i g_j|^2, max_{i<j}|1+g_i g_j|^2); for two real nonzero eigenvalues a, b the
   antisymmetric entry satisfies the bound since 3 - 2u + uw >= (u-3)^2/3 >= 0 (u = a^2+b^2+ab, w = a^2+b^2-ab).  Other cases
   (more eigenvalues, complex, the symmetric entries, non-normal G) not written out.

## 5. THEOREM A' (p = 1 refinement): c(N;1,1) <= 8.081 for every N, i.e. d^2(N;1,1) >= N^2 - 8.09.   [P]
 Extra ingredient: with Pi_v := R(|Fv><v|) (realigned), tr(F Y Q) = <F - Pi_v, R(Y)> for every Y (Q = 1 - vv^*), and
 ||F - Pi_v||_op = 1 + sigma_1 sigma_2 =: 1 + y (Schmidt coefficients of v; computed in the Schmidt frame: F - D(x)D has
 eigenvalues 1 - sigma_i^2 on |ii>, -sigma_i sigma_j +- 1 on span{|ij>,|ji>}; the local-unitary symmetry acts in the realigned
 picture by unitaries on both sides, so op/trace norms are invariant).  Let w1 >= w2 >= ... be |<v|z_j>|^2 (sum = 1, w1 <= sigma_1^2).
   T = {}:      kappa <= (1+y)^2 ||R(B)||_tr^2/||BQ||^2 <= (1+y)^2/(1-w1)          (osr B = 1, ||BQ||^2 >= (1-w1)||B||^2)
   T = {j1}:    kappa <= 1 + 2(1+y)^2/(1-w2) <= 1 + 2(1+y)^2/w1                   (rank R(B_R) <= 2, w2 <= 1 - w1)
 w1 <= 1/2: first bound <= (3/2)^2/(1/2) = 4.5.  w1 >= 1/2: sigma_1^2 >= w1 >= 1/2 so y <= sqrt(w1(1-w1)); max over w1 of
 min(first, second) = 7.0807 (w1 ~ 0.70).  So kappa_1 <= 7.081, c(N;1,1) = 1 + kappa_1 <= 8.081.  QED
 Check [N] checkA2.py (Nelder-Mead maximisation of kappa/bound): worst 0.987 (N=2), 0.869 (N=3): bound valid, nearly tight in
 the regime it is designed for; the slack vs the truth (kappa_1 = 3) is in the coupling of the two regimes.

## 6. q = 2 and the general conjecture [N only]
 Data (G22d, N = 4, 5): c(p,2) = 2, 6, 8, 11, 13, 15, (16 = N^2 at N = 4) for p = 0..6;  c(p,1) = 1, 4, 6, 9, 11, 13, 16.
 Conjecture Q2 (tentative, 5 points p = 1..5): for N >= n0(p), c(p,2) = c(p,1) + 2 = 2p + m(p) + 2   (p >= 1).
 Symmetry c(p,q) = c(q,p) is exact [P: realignment fixes F and swaps sigma_p <-> tau_q], so c(1,q) = 2q + m(q) under the q=1 law.
 The naive symmetric extension c(p,q) = 2(p+q) + m(max(p,q)) - 2 fits every small data point EXCEPT (3,3) at N = 4, which is
 (border-)feasible, c = 16 > 13: so growth is superadditive for balanced (p,q); no general formula is proposed.
 Lower-bound constructions [P]: orthogonal sums of swaps F_{W_l} (dim W_l = m, l <= q) give c(q C(m,2), q) >= q m^2, i.e.
 c(p,q) >~ 2p + sqrt(2pq) along that family; the (3,3)@4 point shows better constructions exist when p q ~ N^2.

## 7. Scaling assessment (honest)
 * What is proved: for q = 1, c(N;p,1) <= C p + O(sqrt p) with C ~ 11.7 (Theorem A), uniformly in N.  This is linear in p, the
   right order (truth 2p + m(p)).  For feasibility (c = N^2) it only gives p >= N^2/11.7 when q = 1 -- weaker than the known
   exact q = 1 results (p >= N(N-1)/2, G19b Thm 3) and it says nothing for q >= 2.
 * The proof uses q = 1 essentially: the right singular vectors of G(x)G' form an orthonormal PRODUCT basis, so removing any
   set T of singular terms leaves osr <= 1 + |T| and keeps Pythagoras (orthogonal ranges).  For osr-q B the singular vectors are
   not products, removing terms destroys the osr bound, and splitting B into its q product pieces loses orthogonality (cross terms
   <B_l P, B_l' P>).  I found no way around this; a general-(p,q) version is OPEN.
 * Line question: the problem lives at p q ~ N^2.  A bound c(p,q) <= kappa p q (any constant kappa) shows infeasibility only for
   p q < N^2/kappa: ON the line theta + T = 2 (exponents), like rank methods (which are capped at N^2/4 by Theorem B).
   Beating the line needs c(p,q) = o(p q) in the balanced regime, e.g. c(p,q) <= C (p + q)^{1+delta} with delta < 1 or
   c(p,q) <= (pq)^{1-eps} N^{...}; equivalently the captured mass must be SUBLINEAR in p q.  For q = 1 the truth IS of this
   additive type (c = 2p + O(sqrt p): no pq term at all), and the q = 2 data (c(p,2) = c(p,1) + 2) also look additive.  If
   c(p,q) <= C(p + q) held for all p, q, feasibility would need p + q >= N^2/C, far beyond the line.  But (3,3) is feasible at
   N = 4 (c = 16 with p + q = 6, so C >= 2.67 there, above the q = 1 slope 2) and the balanced data are superadditive, so the
   additive behaviour is only seen for small parameters; nothing proved or computed decides whether c(p,p) grows like p,
   p^{1+delta} or p^2.
 * Does the technique hint at sublinearity?  Partially: in Theorem A the q-side contributes through osr(B_R) <= 1 + |T|, i.e.
   the product structure is "used once" and the bound is additive (p + |T| + O(1)), not multiplicative.  A proof for general q of
   the shape c(p,q) <= C(p + q + sqrt(pq)) would beat the line (it would force p + q >= Omega(N^2)), and it is
   too strong to expect (it fails as soon as some feasible decomposition has p + q = o(N^2); I know of none either way, the known
   ones, e.g. (C(N,2), 1), have p + q ~ N^2/2).  I have no evidence either way at p q ~ N^2.
Files: notes.md, kap.py (kappa_p), kapD.py (kappa by Schmidt profile), restR.py (refutes the top-term sharpening), symtest.py
 (refutes the symmetric-compression route), checkA.py / checkA2.py (tests of every inequality in Theorems A, A'), d2.py (copy of
 the G22d solver), log6_1_1.txt (N = 6: d^2(6;1,1) = 32.000000, 4/4), kaprun1.txt.
