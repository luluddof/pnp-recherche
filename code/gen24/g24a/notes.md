# G24a — Hrubes (CCC'24, ECCC TR24-026/099) read NON-symmetrically: pieces, trade-off curve, hybrid  (notes, saved as I go)

Tags: [R] read in paper, [E] exact computation, [P] proved here (short, unreviewed), [H] heuristic, [O] open.
Exponents: theta = log_N p, T = log_N q. Circle T = sqrt(4 - theta^2) (App. B). Line theta + T = 2 (floor).

## 0. Source (TR24-026, read in full here; untrusted data)
Thm 9 / Thm 12 (beta version, any field char != 2): t = rho(n)-1 real antisymmetric orthonormal anticommuting e_1..e_t in
Z^{n x n}; e_A = prod_{i in A} e_i, A in ([t] choose k); f multilinear deg <= d = floor(k/2), f(a) = 1 if |a| = k, 0 if |a| < k,
|a| = k mod 2 (Claim 8: f = c prod_{0<=i<k, i=k mod 2}(sum x - i), multilinearised, f(A) = sum_{C subset A} alpha_C);
xi(A)_C = alpha_C [C subset A], xi~(A)_C = [C subset A] (C in ([t] choose <= d));  H_A = e_A (x) xi(A), H~_A = e_A (x) xi~(A)
(n x ns).  HR conditions H_A H~_A^t = I, H_A H~_B^t + H_B H~_A^t = 0  =>  beta(n, C(t,k)) <= n C(t,<=d).
Thm 1: k ~ alpha t with H(alpha) = 1/2 (alpha = 0.110), t ~ 2 log2 n: n^{1 + 2H(alpha/2)} = n^{1.6147}.

## 1. Implementation and exact verification (hr.py, cliff.py, ev.py copied from gen22/g22c)   [E]
cliff.py: real antisymmetric pairwise anticommuting signed-permutation Paulis, family of size rho(n)-1 found by greedy
clique search, checked (e^T = -e, e e^T = I, anticommute) for n = 2..256 (t = 1,3,7,8,9,11,15,16).
hr.py: H_A, H~_A exactly (alpha scaled to integers), family = first N = n k-subsets in lex order; checks HR conditions exactly;
translation to our X (x = A, y = i, u = j, v = B):  2X = T1 + T2,  T1(A,i,B,j) = sum_C alpha_C [C<=A][C<=B] (e_A e_B^t)_ij
 (R-type, = D = sum_sigma H_A[i,s] H~_B[j,s]), T2(A,i,B,j) = D(B,i,A,j) (C-type: terms G_s(A,j) = H~_A[j,s], G'_s(i,B) = H_B[i,s]).
 (This is Rossman-Zhu L.3.15 / D + D^{T1} = 2I.)  Results (p = rr(T1), q = rc(T2), exact):
   N=8  t=7  k=2: HR ok; (56, 56)      N=16 t=8 k=2: (128,128)      N=32 t=9 k=2: (288,288)
   N=64 t=11 k=3: HR ok; (704, 704) [Hrubes's count n C(t,<=d) = 768; exact is n C(t,d) = 704: Johnson rank, see sect. 3]
 Exponents log_N: 1.94, 1.75, 1.63, 1.577.  (App. B is capped at N^2 for these N; the Pauli family has (27,37) at N=8.)

## 2. Natural pieces (C in ([t] choose <= d)) and their profiles (pieces.py, pieces_N*.json)   [E]
Piece C = the n terms sigma = (r, C).  P1_C(A,i,B,j) = alpha_C [C<=A,B](e_A e_B^t)_ij  (part of T1), P2_C = its partner in T2.
 Exact (rr, rc):  P1_C = (n, n M_C),  P2_C = (n M_C, n)   (M_C = #{A in family : A contains C}), for every C with M_C >= 3
 (N = 8, 16, 32; N=8 M_C = 2 gives (8,8)).  rr(T1) = rc(T2) = S;  rc(T1) = rr(T2) = N^2 (full!) at N = 8,16,32.
 => the "wrong-side" rank of a piece is FULL (n M_C): no saving at all from moving a piece across.
 Threshold (move P2_C to R, cheapest M_C first), exact:
   N=16: (128,128) (144,112) (160,96) (192,80) (224,64) (256,48) ... (256,0)
   N=32: (288,288) (416,256) (576,224) (736,192) (864,160) (960,128) (1024,96) ... (1024,0)
 Each move saves n in q and costs n M_C in p.  Asymptotically sum_C M_C over a level c = N C(k,c) >= N, so moving any
 constant fraction of the top level costs p >= N^2 * const: the piece curve is an L: corner (1.615,1.615), then nothing until
 the trivial point p = N^2.  [P, from the exact rank formula n M_C, which is verified only for N <= 32]

## 3. Johnson-eigenspace regrading = the right "non-symmetric reading" (johnson.py, verify_jh.py, solv.py)   [E + P]
Clifford-tensor framework [P, 3 lines]: for ANY F1, F2 : fam x fam -> Q put T1(A,i,B,j) = F1(A,B)(e_A e_B^t)_ij,
 T2(A,i,B,j) = F2(A,B)(e_B e_A^t)_ij.  Since diag(e_A) is invertible, rr(T1) = n rank F1, rc(T2) = n rank F2 EXACTLY, and
 e_B e_A^t = (-1)^{k-|AnB|} e_A e_B^t gives:   X = T1 + T2  <=>  F1 + E o F2 = I,  E(A,B) = (-1)^{k-|A cap B|}  (o = Hadamard).
 Hrubes = F1 = F2 = f(|A cap B|) (rank C(t,d)).  Johnson-invariant choice: F_i = phi_i(|A cap B|), phi_i polynomial of degree
 <= j_i in c = |A cap B|  => rank <= C(t, j_i) (eigenspaces V_0..V_{j_i} of J(t,k)).  Constraint phi1(c) + (-1)^{k-c} phi2(c) =
 [c = k], c = 0..k: square linear system when j1 + j2 = k - 1; checked SOLVABLE for all k <= 25 and all splits (solv.py) [E].
 Exact checks with the full evaluator (verify_jh.py, rank of X - T2 and of T2 directly):
   N=16 k=2: (j1,j2)=(0,1): (16,112); (1,0): (112,16).   N=64 k=3: (2,0): (3456,64); (0,2): (64,3456); (1,1): (704,704).  = n(r1,r2). 
 Exact frontier points (johnson.py; ranks of the N x N restricted F's, flint):
   N=8: (8,56),(56,8)   N=16: (16,112),(112,16)   N=32: (32,288),(288,32)   N=64: (64,3456),(704,704),(3456,64)
   N=128: (128,13312),(1920,1920),(13312,128)   N=256 k=3: (256,30720),(4096,4096) = N^{1.5} exactly,(30720,256);
   N=256 k=4: (3840,26880),(26880,3840).
ASYMPTOTIC CURVE (JH) [P, standard entropy estimates]: t ~ 2 log2 n, k ~ alpha t, H(alpha) = 1/2 (alpha = 0.11003), j1 = b t:
      theta = 1 + 2H(b),   T = 1 + 2H(alpha - b),   b in [0, alpha]   (plus trivial points).
   Endpoints (1, 2), (2, 1); balanced 1.6147 (= Hrubes); at theta = sqrt2: T_JH = 1.786 (circle 1.414);
   sup(theta+T) = 3.229 at the balanced point (circle 2.828).  Table in curves.out.  ENTIRELY OUTSIDE the circle for theta in
   (0,2): e.g. theta = 1.0/1.30/1.42/1.61/1.78: T_JH = 2.00/1.86/1.78/1.61/1.42 vs circle 1.73/1.52/1.41/1.18/0.91.
   Shape: concave, symmetric, "shifted by 1" (theta, T >= 1 always: the factor n of the Clifford block is paid on BOTH sides).
 Finite N (count bounds n min(C(t,j),N); exact for N <= 256): JH beats the App-B COUNT frontier at small N (App B is capped at
   N^2 up to N = 256) -- balanced up to N ~ 2^36-2^44, at theta <= sqrt2 up to N ~ 2^44; App B wins from N ~ 2^50 on (curves.out).
   This is a finite-size effect only (both exponents at N = 2^44 are ~1.57-1.74).
FRAMEWORK LIMITS [P]:
 (L1) any non-trivial point of the Clifford-tensor framework has theta >= 1 and T >= 1 (p = n rank F1, q = n rank F2).
 (L2) (= Hrubes Remark 11/15 in (p,q) form) if fam contains S with E = +1 on S x S (pairwise |AnB| = k mod 2), then restricting
      F1 + F2 = I to S gives rank F1 + rank F2 >= |S|.  For the slice ([t] choose k) (or any 1/poly-dense part of it), S = unions of
      k/2 disjoint pairs: |S| = C(t/2, k/2) = n^{H(alpha) - o(1)} = n^{1/2 - o(1)}  =>  max(theta, T) >= 3/2 - o(1).
      Hence at every theta in [1.32, 1.5] -- in particular theta = sqrt2 -- the framework has T >= 1.5 > sqrt(4 - theta^2):
      it CANNOT lower the envelope at sqrt2, so it cannot lower sup(theta+T) = 2 sqrt2 (rho).  Not excluded by (L1)-(L2):
      points like (1, 1.5) or (1.6, 1.1) [inside the circle, but away from sqrt2]; JH does not reach them (T_JH(1) = 2).
 (L3) For general Clifford families (not a slice; Paulis of mixed weight) the reduced problem "I = F1 + E o F2, E = commutation
      sign pattern" is a twisted cousin of our problem itself; E-cliques (commuting sets) can be as large as n, anticliques
      (anticommuting sets) only O(log n): (L2) then gives nothing.  OPEN, but (L1) still forces theta, T >= 1.

## 4. Hybrid (task 3)
H_K, "Kronecker-graded hybrid" (concrete, computed): X_{N1 n} = X_{N1} (x) X_n with App-B pieces A_s on the first factor and
 Hrubes's T1, T2 on the second; R = sum_{s<s1} A_s (x) T1 + sum_{s<s2} A_s (x) T2 (two thresholds = joint grading by the App-B
 weight s and the Clifford side).  hybrid_exact.py, N1 = 5 (h = 4), n = 8, N = 40: sum of pieces = 96 X exactly; all 18 x 18
 threshold pairs, exact ranks (hybrid_exact.out/.json).  Pareto frontier: (0,1600) (56,1576) (64,1408) (1296,1296) (1408,64)
 (1576,56) (1600,0)  -- essentially trivial (N^2 = 1600); worse than the plain trivial point (N(N-1)/2, 1) = (780, 1).
 Asymptotic exponent (hybrid_exp.py, certified bounds; piece exponents = gamma*factor1 + (1-gamma)*factor2, gamma = log N1/log N;
 the wrong-side ranks of App-B thresholds and of T1/T2 are full (exponent 2, verified exactly for T1, T2 in sect. 2)):
 optimising gamma, the App-B point on the circle, the JH point, and the 4 assignments of the mixed pieces R1(x)T2, C1(x)T1,
 the optimum is ALWAYS gamma = 1 (pure App B) for every theta in [1, 1.9]: the hybrid never beats the circle.  Reason: each
 mixed piece costs a full factor (exponent 2) on the side it is put on, and T1/T2 have full wrong-side rank.
E2 as stated ("App-B accounting with an extra Z_2^t grading") [H, short argument, not computed]: if the pieces respect a
 Z_2^t (or any abelian) grading, ranks split over sectors and sum back to the same window dimensions D(2h +- s): no gain unless
 pieces vanish on whole sectors; and the Clifford-tensor part, if kept as a tensor factor, brings (L1) theta, T >= 1 and, on
 slices, (L2) max(theta,T) >= 3/2.  So the only route left for the Clifford-parity mechanism is a NON-slice family solving
 I = F1 + E o F2 with sublinear ranks -- open (L3), and still confined to [1,2]^2.

## 5. Verdict
 * Hrubes's construction implemented and verified exactly (HR conditions + rank evaluation) for N = 8,16,32,64; exact balanced
   (S,S) = n C(t,d) (not n C(t,<=d)): 56, 128, 288, 704; JH regrading exact up to N = 256 ((4096,4096) = N^1.5 at N = 256).
 * Natural pieces (C-indexed): profile (n, n M_C) / (n M_C, n): wrong side FULL; piece-thresholding gives an L-shaped curve
   (corner at the balanced point, nothing below it until p = N^2).
 * Best non-symmetric reading found (Johnson eigenspace grading): theta = 1 + 2H(b), T = 1 + 2H(0.110 - b): a concave curve
   from (1,2) to (2,1) through (1.615,1.615); strictly OUTSIDE the circle for all theta in (0,2); T(sqrt2) = 1.786 vs 1.414.
 * No point strictly inside the circle in the middle range: proved for the whole Clifford-tensor framework on slice families
   (T >= 1.5 at theta in [1.32,1.5], incl. sqrt2), so this mechanism cannot lower rho.  Finite-size: JH is better than the
   App-B count frontier for N up to ~2^40 (App B trivial for N <= 256), which is irrelevant asymptotically.
 * Hybrid H_K: no gain (optimum = pure App B); exact N = 40 frontier trivial.

## 6. Files (this folder)
 cliff.py (Clifford family), hr.py (Hrubes beta construction, HR check, terms, exact p/q), ev.py (copy of g22c evaluator),
 pieces.py + pieces_N{8,16,32}_k2_lex.json (piece profiles, thresholds), johnson.py + johnson_*.json (JH exact frontiers),
 verify_jh.py (full-evaluator check of unbalanced JH points), solv.py (solvability k <= 25), curves.py/curves.out (asymptotic
 + finite-N counts vs App B), hybrid_exact.py/.out/.json (N = 40 hybrid), hybrid_exp.py (hybrid exponent optimisation).
## 7. Self-review
 * (L2) uses only cliques inside the chosen family; for the lex-first N sets I did not compute the clique size exactly (the
   asymptotic claim is for the full slice / dense subfamilies).  rc(T1) = rr(T2) = N^2 verified only N <= 32.
 * JH solvability for j1 + j2 = k-1 is verified k <= 25, not proved in general (a one-degree slack would not change exponents).
 * The finite-N comparison uses COUNT bounds for both JH (n C(t,j)) and App B (P(2h+-s)); exact ranks only for N <= 256 (JH).
