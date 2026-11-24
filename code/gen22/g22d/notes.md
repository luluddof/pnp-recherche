# G22d — the distance d^2(N;p,q) = inf ||F - A - B||_F^2 (working notes, saved as I go)
Tags: [P] proof here (unreviewed), [N] numerical, [E] exact.
Solver: d2.py = LM on all factors (complex), analytic Jacobian (copy of G21d distlm.py idea).  insp.py / insp2.py: structure.
Equivalent: d^2 = N^2 - c(N;p,q),  c := sup_{X in sigma_p + tau_q, ||X|| = 1} |<F,X>|^2  ("captured mass") [P: J is a cone, so the
 optimal residual E satisfies <E, A+B> = 0, d^2 = <F,E> = ||E||^2].

## 1. Validation [N]
 randJ (random point of J(2,2), N=3): 4/6 restarts 1e-30 (two local minima 0.076).  (3,1)@3 swap: 0 (5/5).  (2,2)@3: 2.000000 (6/6).

## 2. Data [N, 15-16 restarts, all restarts agree unless said]
N=3: (1,1) 5.000000 | (2,1) = (1,2) 3.000000 | (2,2) 2.000000 | (3,1) = (1,3) 0 | (p,0): 9-p (Eckart-Young, exact).
N=4: (1,1) 12.000000 | (2,1) 10.000000 | rest: see below.
Captured mass c = N^2 - d^2:  c(1,1) = 4 (N=3,4), c(2,1) = 6 (N=3,4), c(2,2)@3 = 7.
(continued after a usage-limit pause; no stale jobs found)
N=4: (2,1) 10 (12/12) | (3,1) 7 (12/12) | (4,1) 5.00001 (3/3, coef ~700: border) | (5,1) 3.0 (2/2) | (2,2) 8 (4/4) | (3,2) 5 (1)
     | (5,2) 2.0 (2/2, coef 6-18; G19b vp.py saw ~1.05: UNRESOLVED, may be a local plateau).
N=5: (1,1) 21 (8/8) | (2,1) 19 (5/5).
Residual structure [N]: at every optimum seen, F^* E (Fl1) or R(F)^* R(E) (Fl2) is Hermitian PSD; in one of the two pictures it
 is an orthogonal projection P_K of rank d^2 (so d^2 = dim K integer): (1,1): E = F(1 - P_{W(x)W}), dim W = 2, A+B = swap on W(x)W;
 (2,2)@3: Fl1 and Fl2 both projections of rank 2; (2,1)@3: Fl2 projection rank 3, Fl1 eigenvalues (1.5,.5,.5,.5).
 (2,2)@3: optimal K is a CONTINUUM (2-dim K inside W1(x)W2, dim W_i = 2, varying Schmidt spectra).  ktest.py: F(1 - P_K) is
 in J(2,2) (dist 1e-12 or less, [N]) for K = span{|01>,|10>}, {|00>,|11>}, {|00>,|01>}, {|01>,|02>}, {|01>,|12>}, ...; NOT for
 K = span{|01>-|10>, |22>} (0.25).  So d^2(3;2,2) <= 2 with an explicit K (numerically certified membership).

## 3. More data [N]
N=4 (3,3): restarts 1-3 -> 2e-5..4e-5 with coef 210-280 (border sequence FOLLOWED: calibration passes), restart 0 -> local min
 2.000000 at coef 1.4 with residual = rank-2 projection.  => local minima at integer values exist; a value seen in few restarts
 is only an upper bound.  N=4 (4,2): 3.0003 (2 restarts, coef 300-400, border).  N=5: (2,2) 17, (3,1) 16 (4/4).

## 4. Table of captured mass c = N^2 - d^2  [N]
          N=3   N=4   N=5
 (1,1)     4     4     4
 (2,1)     6     6     6
 (3,1)     9     9     9
 (4,1)     -    11
 (5,1)     -    13
 (6,1)     -    16 (feasible)
 (2,2)     7     8     8
 (3,2)     9f   11
 (4,2)     -    13 (border-type, 3.0003)
 (5,2)     -    14 (2 restarts; G19b saw 1.05 -> maybe 15: UNRESOLVED)
 (3,3)     9f   16 (border-feasible)
LAW q = 1 (fits ALL 12 q=1 data points, N = 3,4,5): write p = m(m-1)/2 + j, 0 <= j < m, m <= N.  Then
     d^2(N; p, 1) = N^2 - m^2 - 2j     (c = m^2 + 2j: swap on an m-dim W exactly for j = 0; each further unit of p adds 2).
 j = 0 values are attained at finite points (A + B = swap on W(x)W via I(x)I - 2 P_anti); j > 0 only by divergent (border)
 sequences.  Upper bound for j = 0 is a [P] construction; j > 0 upper bound and all lower bounds are [N].
STABILITY: for fixed (p,q), c(N;p,q) is nondecreasing in N [P: compression to V' (x) V' is a contraction, preserves rank/osr and
 maps F_N to F_{N'}; embedding gives >=] and empirically constant once N >= (dimension it needs): (1,1) N>=2, (2,1),(3,1) N>=3,
 (2,2) N >= 4 (8 = two orthogonal copies of the 2-dim swap: (1,1)+(1,1), superadditivity on orthogonal subspaces [P]).
 Conjecture S: c(p,q) := lim_N c(N;p,q) is finite, d^2(N;p,q) = N^2 - c(p,q) for N >= n0(p,q).
Candidate laws tested: max_k(kN - p - k^2 q) and the symmetric variant: FAIL everywhere (far too small, they grow ~N(N-p-q)
 while the truth is N^2 - O(1) for fixed p,q).  (p+1)(q+1): fits (1,1),(2,1) only.  (sqrt p + sqrt q)^2: fails (2,1).
 2(p+q) + small: c - 2(p+q) = 0,0,0,1,1,1,2 | 0,1,0 for the listed pairs: roughly right in this tiny range, cannot persist
 (feasible (p,p) with p ~ N needs c ~ p^2).

## 5. Quantitative Theorem C (rank-method distance bound)  [P]
Prop. Let Phi = (+)_t Phi_t be linear with sum_t ||Phi_t(X)||_F^2 <= ||X||_F^2 for all X, c1_t (c2_t) = max rank of Phi_t on
 rank-one (osr-one) operators, s_{t,1} >= s_{t,2} >= ... the singular values of Phi_t(F).  Then for all A in sigma_p, B in tau_q,
 and hence on the closure J(p,q):   ||F - A - B||^2 >= sum_t sum_{i > c1_t p + c2_t q} s_{t,i}^2.
Proof. ||F-A-B||^2 >= sum_t ||Phi_t(F) - Phi_t(A+B)||^2; rank Phi_t(A+B) <= c1_t p + c2_t q (subadditivity, each block separately);
 Eckart-Young-Mirsky in each block.  Both sides are continuous in (A,B) -> valid on the closure (border decompositions).  QED
Instance: partition [N] = I_1 u ... u I_s, |I_t| = k_t; Phi_t(X) = X restricted to x in I_t, v in I_t (disjoint entry sets, so
 sum ||Phi_t X||^2 <= ||X||^2).  Viewed as (xy|uv) matrix: profile (1, k_t^2, k_t N) (Thm C(a)); viewed realigned (xu|yv):
 (k_t^2, 1, k_t N) (Thm C(b)) -- same entries, so per block take the better; Phi_t(F) is a partial permutation (all s = 1):
     d^2(N;p,q) >= RB(N;p,q) := max_partitions sum_t max(0, k_t N - p - k_t^2 q, k_t N - k_t^2 p - q)  >= max_k (kN - p - k^2 q).
 (rb.py).  Consequence of Theorem B (G20d, class B, char 0): for any class-B Phi, an Eckart-Young bound of this form is positive
 only if c1 p + c2 q < R, hence only if pq < N^2/4.  So NO certificate of this (contraction + rank + Eckart-Young) type can give
 d^2 > 0 for pq >= N^2/4 -- in particular none for (2,2) at N = 3 (pq = 4 > 2.25) or (2,1) at N=3 for the single-block type.
Comparison (rb.py): RB vs truth:  (3;1,1) 3 vs 5; (4;1,1) 8 vs 12; (5;1,1) 15 vs 21; (4;2,1) 4 vs 10; (5;2,1) 10 vs 19;
 (4;3,1) 2 vs 7; all other infeasible data points: RB = 0 while d^2 = 2..8.  The law max_k(kN-p-k^2 q) is NOT the distance law:
 it is a strict lower bound at every infeasible data point with p, q >= 1.

## 6. (5,2)@4 resolved and q = 2 pattern [N]
logc4_5_2.txt: restarts 2.0, 2.0 (finite coef: local minima) and 1.000562 (coef 461, border) -> d^2(4;5,2) = 1 (c = 15), agreeing
 with G19b's 1.05.  (4,2)@4: 3.0002-3.0008 (5/5, border) -> c = 13.
q = 2 row (N = 4, 5 where available): c(p,2) = 2, 6, 8, 11, 13, 15, 16 for p = 0..6;  q = 1 row: 1, 4, 6, 9, 11, 13, 16.
 Difference 1,2,2,2,2,2,0: conjecture (tentative, 7 points) c(p,2) = min(c(p,1) + 2, saturation), for N large enough
 (N = 3: (2,2) gives 7, not 8, because 8 needs two orthogonal 2-dim subspaces).

## 7. Verdict
 * Data/table: sections 2-4, 6.  Residuals at optima: F times an orthogonal projection in Fl1 or Fl2 => d^2 INTEGER at every
   infeasible point computed (14 points; all values integral to >= 6 digits) [N].  Integrality conjecture: open.
 * Law: q = 1 closed form d^2 = N^2 - m^2 - 2j (p = C(m,2) + j) [N, 12 points; j=0 construction P]; stability c(N;p,q) = c(p,q)
   for N >= n0 [N]; c nondecreasing in N and superadditive over orthogonal subspaces [P].
 * Proven lower bound: quantitative Theorem C (block version) section 5 [P]; it is strictly below the truth at every infeasible
   point with p, q >= 1, and by Theorem B no bound of that type can be positive for pq >= N^2/4.
 * (2,2)@3 = 2 [N]: inf over a continuum of border optimizers; F(1 - P_K) in closure J(2,2) for explicit K = span{|00>,|11>} or
   {|01>,|10>} (LM 1e-12 with |A| ~ 1000, i.e. border only) -> d^2 <= 2 [N]; lower bound >= 2 only numerical.  No certificate
   found; any certificate must be outside the contraction+rank class (section 5) and of degree >= 7 if algebraic (G21d).
 * Scaling: for FIXED (p,q) the analytic distance is N^2 - O(1) vs rank methods N^2 - Theta(N(p+q)) -- a huge gain, but the
   feasibility question lives at pq ~ N^2 where one needs the growth of c(p,q); data c/pq = 4, 3, 2, 3, 2.75, 2.6, 1.83, 1.63,
   1.5 is too small-scale to extrapolate.  The route scales only if one proves c(p,q) <= kappa pq with kappa < 4 (kappa = 4 is
   the rank-method line); q = 1 gives c ~ 2p + sqrt(2p) (consistent with Thm 3), no proof of any upper bound on c found.
Files: d2.py (solver), insp*.py, ktest.py, rb.py, log*.txt, out_*.json, notes.md.
