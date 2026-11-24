# G24c — q >= 2 for the captured mass c(N;p,q), and data near feasibility (working notes, saved as I go)
Tags: [P] proof here (unreviewed), [N] numerical, [E] exact.  Conventions as G22d/G23d: X[x,y,u,v] = <xy|X|uv>, F = swap,
Fl1 = (xy|uv) (rank A <= p), Fl2 = (xu|yv) = realignment R (osr B <= q), c = N^2 - d^2, Q = 1 - P.

## 0. Solver (new, faster): vp.py  [N]
Variable projection: A is eliminated exactly (Eckart-Young of F - B in Fl1), B = sum_l G_l (x) H_l optimised by L-BFGS
(f = sum_{i>p} sigma_i(F-B)^2, gradient by the envelope theorem -2 Re<E, dB>; finite-difference check 202.0558 vs 202.0558).
Uses symmetry c(p,q) = c(q,p) to optimise over the smaller side.  10-100x faster than d2.py (LM), so N = 5, 6 are reachable.
Re-validation (val.txt; all agree with G22d/G23d):
 (3,1)@3 0 (4/4) | (2,2)@3 2.000000 (4/4) | (1,1)@4 12 (3/3) | (3,3)@4 0.0009, 0.0014, 0.0104 (border, followed), 2.0 x3 (local)
 (4,2)@4 3.0013 (border; LM 3.0002) | (5,2)@4 1.0101 (border; LM 1.0006) | (6,1)@4 0 | (2,2)@5 17 | (3,1)@5 16 | (3,2)@4 5.
 Caveat: on border points L-BFGS converges more slowly than LM (values 0.001-0.01 above the integer); local minima at integers
 above the optimum are common near feasibility => only the BEST of many restarts is an (upper-bound) estimate of d^2.

## 1. Theory for q >= 2  [P unless marked]
(R1_q) EXACT FORMULA FOR EVERY q.  c(N;p,q) = p + kappa_{p,q}(N),
        kappa_{p,q} := sup_{P rank p proj} sup_{osr B <= q, BQ != 0} |tr(F B Q)|^2 / ||BQ||^2      (Q = 1 - P),
  and (symmetry) = q + kappa_{q,p}.  The G23d proof of (R1) never uses q = 1: (<=) P := proj on (ker A)^perp, XQ = BQ,
  |tr FXP| <= sqrt(p)||XP||, Cauchy-Schwarz + Pythagoras;  (>=) A := (F - tB)P has rank <= p, ||(F - tB)Q||^2 optimised in t.
(GP) General splitting principle.  For ANY projection Pi >= P (P = support projection of A):
        |tr FX|^2 / ||X||^2  <=  rank Pi + |tr F B(1-Pi)|^2 / ||B(1-Pi)||^2.      (same proof with Pi in place of P)
(C) Concentration lemma (all q).  If ||BP||^2 <= theta ||B||^2 then |tr FX|^2/||X||^2 <= p + (sqrt q + sqrt(p theta))^2/(1-theta).
  Proof: |tr FBQ| <= |tr FB| + |tr FBP| <= ||R(B)||_tr + sqrt(p)||BP|| <= (sqrt q + sqrt(p theta))||B||, ||BQ||^2 >= (1-theta)||B||^2,
  then (R1_q).  theta = 1/2: ratio <= 3p + 4q.  => In ANY configuration capturing more than 3p + 4q, more than half of ||B||^2
  lives on the p-dimensional support of A (B is "hidden behind A").  This is the only regime a general bound must handle.
(L1) Realigned-orthogonal decomposition is stable under product-rectangle cuts.  Write R(B) = sum_l s_l e_l e_l'^T (SVD), i.e.
  B = sum_l B_l, B_l = s_l E_l (x) E'_l, <E_l,E_m> = <E'_l,E'_m> = delta (HS-orthogonal products).  For Pi0 = P_Y (x) P_Y'
  (any subspaces Y, Y'), Om = 1 - Pi0:      sum_l ||B_l Om||^2 <= 2 ||B Om||^2   (and same for left cuts (1 - P_X(x)P_X') B).
  Proof: right mult. by Pi0 is R -> L R M in the realigned picture (L = 1(x)P_Y^T, M = 1(x)P_Y'^T, projections), so
  ||B Om||^2 = ||R - LRM||^2 = ||(1-L)R||^2 + ||LR(1-M)||^2 >= max(||(1-L)R||^2, ||R(1-M)||^2) >= (1/2) sum s_l^2 (a_l + b_l),
  a_l = ||(1-L)e_l||^2, b_l = ||(1-M)e'_l||^2 (orthonormality of e'_l resp. e_l), while ||B_l Om||^2 = s_l^2 (1-(1-a_l)(1-b_l))
  <= s_l^2 (a_l + b_l).  QED.  Also osr(B Om) <= 2q for any rectangle, regardless of dim Y, dim Y'.
(A_q) THEOREM A FOR THE PRODUCT-SVD CLASS.  Suppose osr B <= q and B has an SVD B = sum_j beta_j w_j z_j^* whose left and right
  singular vectors are PRODUCT vectors.  Then for every rank-<=p A and 0<theta<1, t = floor(p/theta):
        |tr F(A+B)|^2 <= [ p + t + (sqrt(q + t) + sqrt(p theta))^2/(1-theta) ] ||A+B||^2        (theta = 1/2: <= 13p + 4q).
  Proof: G23d Theorem A verbatim; the only q-dependent step (5) becomes rank R(B_R) <= rank R(B) + sum_{j in T} osr(w_j z_j^*)
  <= q + |T|.  The class contains: q = 1; orthogonal direct sums B = (+)_l G_l (x) H_l with mutually orthogonal ranges and
  co-ranges (e.g. sum_l P_{W_l} (x) P_{W_l}, W_l orthogonal: exactly the G23d lower-bound family, c ~ 2p there, consistent);
  B = (U(x)U') D (V(x)V')^* with D diagonal (factors with common singular bases).  ADDITIVE, N-independent, but restricted.
(A_q') Rectangle variant of Theorem A (q = 1, new): T may be replaced by its hull I x K (|I|,|K| <= |T|): every index outside the
  hull is light, so ||B_R P||^2 <= theta ||B_R||^2 still holds; cost rank(B_T) <= |T|^2, gain osr(B_R) <= 2 (instead of 1+|T|).
OBSTACLE for general osr-q B (precise form).  Every route I tried reduces to the following missing lemma:
  (L2?) for B = sum_l B_l (realigned-orthogonal products) and a rank-p projection P there are subspaces X, X' (or Y, Y') of
  dimension <= f(p,q,theta) with   || (1 - P_X(x)P_X') B P ||^2 <= theta || (1 - P_X(x)P_X') B ||^2.
  Given (L2) with f polynomial, (L1)-free direct argument: GP-type split with the OUTPUT projection Pi_W = P_X (x) P_X' (left
  projections always give HS-orthogonal pieces: <Pi_W Y, (1-Pi_W) Y> = 0), rank(Pi_W B) <= f^2, osr((1-Pi_W)B) <= 2q, so
        c(N;p,q) <= p + f^2 + (sqrt(2q) + sqrt(p theta))^2/(1-theta)      -- polynomial in p, q, independent of N.
  For q = 1, (L2) holds with f = p/theta (hull of the heavy set in the singular grid).  For q >= 2 each B_l has its own singular
  grid; the hull of all heavy sets is a rectangle that is NOT aligned with any single B_l, and then (1 - Pi_W) B_l^{light} has
  non-orthogonal outputs and the light/heavy bookkeeping (sum_j beta_j^2 ||P z_j||^2) breaks ("cross-contamination").  With the
  singular vectors of B itself (orthonormal, Theorem-A bookkeeping intact) the vectors are entangled and osr(w_j z_j^*) can be
  N^2.  Writing B = sum_j w_j (B^* w_j)^* over an output product basis keeps osr(w_j (B^* w_j)^*) <= q and Pythagoras, but the
  heavy-set count needs a Bessel bound for the normalised B^* w_j, which is not available.  I did not find (L2) or a
  counterexample; NOT proved for q >= 2 in general.
  Lifting attempts (stack G_l into V (x) C^q, contract with the max. entangled Phi) turn B into a product but the contraction
  is entangled across the cut; no gain.

## 2. Data near feasibility  [N]  (vp.py; logs/ ; best of R restarts, R = 4-8 per run, some points with 2 runs)
c = N^2 - d^2 ("f" = feasible or border-feasible: d^2 -> 0 with growing coefficients, or ~1e-6).  ">=" = best value found is
non-integral / still decreasing at the iteration cap, so only a lower bound on c (upper bound on d^2).
 N=5 (25):  (3,2) 11 (8/8) | (3,3) 16 (9.0017 border, 3/8) | (4,2) 14 (8/8) | (4,3) >= 17.99 (7.0105 border, 2/8 -> 18)
            (4,4) >= 20.93 (4.1478 -> 4.0711 after 2e5 its, border -> 21?; 5.0 local min in 8 other restarts)
            (5,2) 16 (8/8) | (5,3) >= 21.87 (3.128; 2nd run 3.376: unresolved) | (5,4) >= 23.88 (1.1236, unresolved)
            (5,5) 25 f (0.0095, 7/8) | (6,2) 18 (6/8) | (6,3) 25 f (0.0152, 5/8) | (6,4) f | (6,5) f (8e-5) | (7,2) 20 (5.003, 6/6)
 N=6 (36):  (2,2) 8 | (3,2) 11 | (3,3) 16 (20.0021 border) | (4,2) 14 | (5,2) 16 (5/6) | (4,4) >= 24 (12.000000, 2/4)
            (5,5) >= 29.93 (6.0742, 1/4; 8.0 local) | (6,6) 36 f (r1: 1e-6 with coef 52, i.e. essentially an exact decomposition)
 N=7 (49):  (4,4) >= 24 (25.000000, first restart; run killed at the time limit).
Stability in N: c(2,2) = 8 (N=4,5,6), c(3,2) = 11 (4,5,6), c(3,3) = 16 (4f-border,5,6), c(4,2) = 13 (N=4) -> 14 (N=5,6),
 c(5,2) = 15 (4) -> 16 (5,6), c(4,4): >= 20.93 (5) -> >= 24 (6) -> >= 24 (7, one restart).  So the stabilisation index n0(p,q)
 grows (n0(4,2) = 5, n0(4,4) >= 6).  NOTE: G23d's tentative Conjecture Q2 c(p,2) = c(p,1) + 2 is FALSE: c(4,2) = 14 = c(4,1) + 3
 at N = 5, 6 (8/8 and 6/6 restarts, integer).
Integrality: every converged value is an integer (to 1e-3, border points converge from above); the non-integral entries
 (4,4)@5, (5,3)@5, (5,4)@5, (5,5)@6 are unconverged border runs, not counterexamples.
Feasibility frontier (smallest p+q, balanced side):  N=3: (3,1) [and (2,2) infeasible]; N=4: (3,3) pq = 9 = 0.56 N^2, p+q = 6;
 N=5: (6,3) pq = 18 = 0.72 N^2, (5,5) pq = 25 = N^2, (7,2) INFEASIBLE (c = 20), (4,4) infeasible (c ~ 21); N=6: (6,6) pq = N^2
 feasible, (5,5) infeasible (c ~ 30).
Shape of c at fixed p+q (N=5): p+q = 6: (5,1) 13, (4,2) 14, (3,3) 16;  p+q = 8: (7,1) 18, (6,2) 18, (5,3) >= 21.9, (4,4) ~ 21;
 p+q = 9: (8,1) 20, (7,2) 20, (6,3) 25f;  p+q = 10: (9,1) 22, (5,5) 25f.   => balance helps: c is NOT a function of p+q.
Diagonal (stabilised or largest N): c(p,p) = 2? (p=1: 4), 8, 16, >= 24, >= 30 (N=6, unconverged), i.e. increments 4, 8, 8, >= 6:
 c(p,p) ~ 8(p-1) for p = 2..4, ratio c/pq = 2, 1.78, 1.5, ratio c/(p+q) = 2, 2.67, 3.0, 3.0.
 (6,6)@6 third restart 2.6e-4 (coef 14): feasibility confirmed 2/3.  All jobs killed at the end.

## 3. Verdict
 * Proved for q >= 2 [P]: (R1_q) c(N;p,q) = p + kappa_{p,q}(N) exactly (every q); (C) any X capturing > 3p + 4q has more than
   half of ||B||^2 on the p-dim support of A; (L1) stability of the realigned-orthogonal splitting under product cuts;
   (A_q) Theorem A for osr-q operators whose SVD has product singular vectors (orthogonal direct sums of products, common
   singular bases): c-ratio <= p + t + (sqrt(q+t) + sqrt(p theta))^2/(1-theta) <= 13p + 4q, N-independent and ADDITIVE.
   NOT proved: any N-independent bound for general osr-q B, q >= 2.  The obstacle is reduced to one lemma (L2?) (product
   rectangle of dimension f(p,q) cutting off the part of B that sits on P); (L2) + the proved steps would give
   c <= p + f^2 + (sqrt(2q) + sqrt(p theta))^2/(1-theta).
 * Data trend [N]: c is NOT additive in p+q (balance helps at fixed p+q), but on the diagonal it grows roughly linearly,
   c(p,p) ~ 8(p-1) = 4(p+q) - 8 for p = 2..4 (c/pq falls 2 -> 1.5), with N-dependence still present ((4,4): 21 -> 24 from N = 5
   to 6).  Four data points, finite-N effects visible: no extrapolation is justified.
 * Feasibility: (3,3)@4 (pq = 0.56 N^2), (6,3)@5 (0.72 N^2), (5,5)@5 and (6,6)@6 (pq = N^2) are (border-)feasible; (4,4)@5,
   (7,2)@5, (5,5)@6 are not.  So at N = 4, 5, 6 the swap is feasible at pq <= N^2 -- at or BELOW the line pq ~ N^2.
 * Honest assessment: the distance d^2 is an exact infeasibility certificate, so whether the route can beat the line is
   whether feasibility actually requires pq >> N^2.  The small-N data say the opposite (feasible at pq = 0.56-1.0 N^2; (N,N)
   feasible for N = 4?, 5, 6 -- (4,4)@4 is feasible since (3,3)@4 is).  If (N,N) stays feasible for all N, NO method can give
   an exponent above the line; the only hope would be that the small-N feasibility is a small-N accident while the linear
   diagonal law c(p,p) ~ 8p persisted with bounded N-dependence (then feasibility would need p ~ N^2/8, far above the line).
   The observed N-growth of c(4,4) and the feasibility of (N,N) at N = 5, 6 point against that.  Next decisive measurement:
   is (N,N) (or (N-1,N-1)) feasible at N = 7, 8, and does c(p,p) at fixed p stabilise in N (c(4,4) at N = 7, 8)?
Files: notes.md, vp.py (variable-projection solver), d2.py (copy of the LM solver, not run), runlist.sh, val.txt (validation),
 logs/*.txt (all runs), vp_*.json, best_*.npy (best points for warm starts).
