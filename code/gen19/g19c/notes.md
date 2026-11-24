# G19c — "three stars + row-rank r" for the hexagon permutation operator: N = 3 scan
(work log; saved as I go.  Tags: [N] numerical evidence, [P] proved, [E] sketch)

## 1. Instance (exact)
* Strand dims n = (n01,n02,n12) (N = n01 = n02 = n12 in the square case).  Sites j = 0,1,2 = columns; H_j = the two strands at j,
  d_j = product of their dims (N^2); D = d0 d1 d2 = N^6.  P = S_01 (x) S_02 (x) S_12, a D x D permutation matrix,
  rows = outputs (o0,o1,o2) (row 1 of the grid), columns = inputs (i0,i1,i2) (row 0).  Same P as gen18/g18f/hexlib.perm_P (checked).
* Row rank r: A is a D x D matrix of rank <= r (cut inputs | outputs).
* Star centred at c, leaves l,m: in the operator-space view T[p0,p1,p2], p_j = (o_j,i_j) in [K_j], K_j = d_j^2 = N^4,
      star_c = sum_{a<N^2, b<N^2} C[p_c,a,b] LA[p_l,a] LB[p_m,b]        (triangle network with bonds N^2, N^2, 1)
  i.e. ANY tensor in Op(H_c) (x) V_l (x) V_m with V_l, V_m subspaces of dimension N^2 of the N^4-dimensional operator spaces
  (multilinear rank (*, N^2, N^2)).  Three stars = a tensor in  Op0 (x) V1 (x) V2 + W0 (x) Op1 (x) W2 + U0 (x) U1 (x) Op2.
* Question: min r such that P = A + star_0 + star_1 + star_2.
* Parameter count, N = 3: per star C: 81*9*9 = 6561, leaves 2*81*9 = 1458 -> 8019; three stars 24057 (gauge GL(9)^2 per star: -486
  -> 23571 effective; = 3*6561 linear core unknowns + 6 Grassmannians Gr(9,81) of dim 648).  A: r(2*729 - r).
  Equations: 729^2 = 531441.  Naive count would need r ~ 350: every solution with small r is non-generic (as at N = 2:
  1152 star parameters + 127 vs 4096 equations, and r = 1 exists).
* Trivially feasible r ("product point"): S_N = 1 - Y, Y = 2 Pi_antisym, rank N(N-1)/2.  P = prod_e (1 - Y_e); every term with
  at least one identity factor is inside a star (the identity bond has dimension 1, OSR(Y) <= N^2), the remaining term -Y(x)Y(x)Y
  has row rank (N(N-1)/2)^3:  r = 1 (N=2), 27 (N=3), 216 (N=4).  Mixed dims: r_prod = prod_e n_e(n_e-1)/2.
* Symmetries: GL(n)^6 on the six strands (local at each site on inputs, and on outputs after the permutation) maps
  solutions to solutions; plus the star gauges.

## 2. Solver (starlib.py)
* Variable projection: A eliminated by Eckart-Young, f(x) = sum_{k>r} sigma_k(P - B(x))^2.  Levenberg-Marquardt with the projected
  Gauss-Newton matrix J^T Pi_N J (Pi_N = projector on the normal space of the rank-r manifold at the truncation), solved
  matrix-free by CG; leaves re-orthonormalised (QR gauge) at every accepted step.  Checked: adjoint identity, finite-difference
  derivative, agreement of the star map with hexlib.build.
* Calibration N = 2 (64 x 64), random starts (random orthonormal leaves, cores 0.1*gaussian):
  r = 3, 2, 1: 10/10 starts reach f < 1e-9 (9-10/10 below 1e-13) in ~15 LM iterations (0.2-0.8 s each).  |A|_F = 8.10 ~ 2^3 at r=1.
  r = 0: 10/10 stall at f = 25.00 (one at 25.93), star norms growing (53) with |B| = 6.2: infimum 25 approached at the border;
  r = 0 is infeasible (g18f section 6) -- consistent.
* Calibration mixed dims n = (3,2,2) (144 x 144, product value r = 3): r = 4: 9/10, r = 3: 9/10 random starts reach f < 1e-9
  (the r = 4 optimum has sigma_4 = 2.6e-8: it is a rank-3 solution).
* Calibration N = 3, r = 27, random start seed 0 (core scale 0.1): f = 7.9e-14 after 41 LM iterations, 127 s; sigma_27 = 8.1,
  sigma_28 = 4.6e-8, star norms 31.9/33.8/30.7, |B| = 55.3, |A| = 52.5 (bounded: an honest exact solution, not a border one).
  => the solver PASSES the calibration the L-BFGS of g18f failed.  (ck_n333_r27_s0.npz)
* n = (3,2,2), r = 2 (one below the product value 3): 0/10; ALL starts plateau at f = 25.01-25.02 (one at 26.0) with growing
  norms (stars 59-95, |A| = 126): the infimum looks like 25.0, approached only at the border, never 0.  Same constant 25 as
  N = 2, r = 0.  [N]  (r = 1 run killed: implied by monotonicity.)
* exact_check.py: the product point is verified in exact integer arithmetic for N = 2 and N = 3:
  P = -Y(x)Y(x)Y + st0 + st1 + st2 with Y = 1 - SWAP, st0 = Y01 Y02 - Y01 + 1, st1 = Y01 Y12 - Y12, st2 = Y02 Y12 - Y02;
  rank(Y^{x3}) = 27, all 27 singular values = 8; multilinear ranks of the stars (81,9,9), (9,81,9), (9,9,81).  (product_n333.npz)
* N = 3, r = 27, random starts: seed 0 (scale 0.1) f = 7.9e-14 [127 s]; seed 2 (0.1) f = 9.9e-14 [200 s], sigma_27 = 8.10,
  sigma_28 = 6.5e-8, |A| = 48.3; seed 1 (scale 0.03) stuck on a slow plateau f ~ 83 after 520 s.  => 2/3.
* n = (3,3,2) (324 x 324, product value 9): r = 9: 3/4 starts < 1e-9 (4th at 1.5e-6 when stopped by the time limit).


## 3. N = 3 runs (729 x 729, 24057 star parameters; runs.log, ck_*.npz).  f = sum_{k>r} sigma_k(P-B)^2, |P|^2 = 729
r | seed | init | f | LM its | s | star norms | |B| | |A| | sigma_r | sigma_r+1
27 | 0 | random, core scale 0.1 | 7.853e-14 | 41 | 127 | [31.9, 33.8, 30.7] | 55.3 | 52.5 | 8.12 | 4.6e-08
27 | 1 | random, core scale 0.03 | 8.278e+01 | 284 | 521 | [52.7, 263.4, 345.5] | 447.4 | 447.2 | 5.69 | 2.9
27 | 2 | random, core scale 0.1 | 9.858e-14 | 143 | 200 | [25.0, 32.1, 28.1] | 51.3 | 48.3 | 8.1 | 6.5e-08
27 | 3 | random, core scale 0.1 | 6.911e-14 | 58 | 73 | [31.0, 35.4, 33.3] | 58.1 | 55.5 | 8.55 | 6.6e-08
27 | 4 | random, core scale 0.3 | 9.105e-14 | 140 | 157 | [39.1, 37.9, 32.8] | 64.5 | 62.2 | 8.17 | 4.9e-08
26 | 0 | random, core scale 0.1 | 3.154e+01 | 357 | 540 | [211.6, 174.4, 52.8] | 294.8 | 294.3 | 6.98 | 4.1
26 | 1 | random, core scale 0.1 | 8.808e+01 | 152 | 251 | [74.1, 249.3, 350.4] | 452.6 | 452.4 | 5.48 | 2.9
26 | 2 | random, core scale 0.1 | 3.200e+01 | 170 | 230 | [34.6, 147.8, 66.0] | 185.3 | 184.5 | 7.52 | 4.1
26 | 3 | random, core scale 0.3 | 3.170e+01 | 173 | 231 | [113.5, 114.6, 61.5] | 189.9 | 189.2 | 7.6 | 4.1
26 | 100 | product_n333.npz | 3.268e+01 | 81 | 231 | [121.6, 108.2, 36.7] | 133.1 | 132.1 | 6.93 | 4.2
26 | 101 | ck_n333_r27_s0.npz | 3.225e+01 | 79 | 231 | [76.3, 83.7, 60.9] | 139.2 | 138.2 | 6.82 | 4.2
20 | 0 | random, core scale 0.1 | 1.568e+02 | 155 | 201 | [85.3, 67.8, 26.3] | 116.5 | 115.3 | 6.23 | 3.7
9 | 0 | random, core scale 0.1 | 3.176e+02 | 167 | 202 | [42.0, 25.0, 41.7] | 41.9 | 38.6 | 11.2 | 1.9

Summary N = 3 [N]:
* r = 27 (product value): 4/5 random starts reach f < 1e-13 (73-200 s each), bounded norms (|B| ~ 51-55, |A| ~ 48-53), sigma_27 ~ 8.1,
  sigma_28 ~ 5e-8: exact solutions, not border.  The failing start (core scale 0.03) sits on a slow high plateau (f ~ 83, norms ~ 450).
* r = 26: 0/6 (4 random starts + 2 continuations: from the exact product point and from a numerical r = 27 solution).  Five of the
  six end on the SAME plateau f = 31.5-32.7, still creeping down (31.83 -> 31.54 over the last 250 iterations of the longest run;
  from the product point: 64.0 -> 37.1 -> 32.7) while the norms grow (|B| ~ |A| ~ 130-300): a border-type infimum that is
  clearly > 0 (not converged; the smaller instances suggest ~25).  One start in the bad basin (88).
* r = 20: f = 157 (one start, 200 s);  r = 9 = N^2: f = 318 (one start, 200 s, not converged).  Roughly 18-30 per missing rank unit.
* NOT done: the full r = 27..1 x 20 restarts scan (would cost ~15 CPU-hours at ~100-250 s per start); replaced, using the
  monotonicity of feasibility in r, by the decisive pair r = 27 / r = 26 plus r = 20, 9, and by the mixed-dimension ladder below.

## 4. Mixed-dimension ladder (same question, strand dims (n01,n02,n12), star bonds n_e^2; product value prod n_e(n_e-1)/2) [N]
instance | D | product r | at product r | one below
(2,2,2) | 64 | 1 | 10/10 exact | r = 0: 0/10, plateau 25.00 (infeasible, proved in g18f)
(3,2,2) | 144 | 3 | 9/10 exact | r = 2: 0/10, plateau 25.01-25.02 (norms -> 60-126)
(4,2,2) | 256 | 6 | 1/2 exact (other 1e-4 at the 50 s limit) | r = 5: 0/2, 25.9-27.6 (50 s limit, norms 186)
(3,3,2) | 324 | 9 | 3/4 exact (4th 1.5e-6 at the limit) | r = 8: 0/4, 26.2-26.3 (65 s limit, norms 185)
(3,3,3) | 729 | 27 | 4/5 exact | r = 26: 0/6, 31.5-32.7 (not converged, norms 130-300)
In every instance the threshold found is exactly the product value, and one unit below the residual plateaus at ~25 (of the 64 =
8^2 carried by one singular value of Y(x)Y(x)Y), with diverging parameters: not even a border solution appears.

## 5. Interpretation
* r*(2) = 1 [P, g18f], r*(3) = 27 [N, real field], and on the ladder r* = prod_e n_e(n_e-1)/2: the data fit r*(N) = (N(N-1)/2)^3
  ~ N^6/8, i.e. the product point Y(x)Y(x)Y is optimal.  This is far ABOVE N^2 = 9 (and N^{2-delta}); at N = 3 the cost is
  log2(27) + 2 log2(81) = 4.75 + 6.34 = 11.1 bits against S = 9.51: no gain, and the trend is r* ~ N^6 (row term of full order),
  neither ~N^2 nor ~N.  The N = 2 "gain" (4 of 6 bits) is the accident N(N-1)/2 = 1.
* Theorem F gives nothing here (b = M^2), so this is new information about the corner, but only numerical.
* The exact solution at the threshold is the known product identity (integer entries, verified exactly); the numerical r = 27
  solutions have sigma_27 ~ 8.1 like the product point (all singular values 8); I did not test whether they lie in its GL(3)^6 x gauge orbit.
* No entangled term found.

## 6. Verdict: NO GAIN AT N = 3 (numerical evidence, calibrated; not a proof)

## 7. Adversarial self-review
* Evidence, not proof: local optimiser, real field only (complex solutions not searched), 6 runs at r = 26 (not >= 20), r <= 25 only
  through monotonicity and two single runs (r = 20, 9).  The requested 27 x 20 scan was NOT run (budget).
* Calibration is honest but thin at N = 3: 4/5 at r = 27, and the r = 26 runs come from the same initial distributions, the same
  time limits or longer (540 s vs 73-200 s for the successes).  A measure-zero or badly conditioned family of r <= 26 solutions with
  a tiny basin cannot be excluded; continuation from the product point (the natural path) also fails.
* The r = 26 plateau is not converged (31.5 and creeping); "25" is an extrapolation from the smaller instances.  I do not have a
  derivation of the constant 25 (the naive diagonal count for N = 2, r = 0 gives 27).
* The Gram-matrix eigensolver used for D >= 200 limits f to ~1e-13 accuracy; successes were re-checked with a full SVD (norms()).
* Only the minimal shape (exactly three stars, bonds exactly N^2, one star per centre) was tested; more column terms, bonds slightly
  above N^2, or other triangle shapes of the same cost are untouched.
* Next step N = 4: 4096 x 4096, 3 x 73728 unknowns, product value 216 vs N^2 = 16.  With this solver: one eigendecomposition
  ~15-20 s, J products ~0.3 s, ~40 CG steps per LM step -> ~30 s per LM step, ~100 steps: ~1 CPU-hour per start; the pair
  r = 216 / r = 215 with 5 starts each ~ 10 CPU-hours.  The ladder (4,3,2) [D = 576, product 18], (4,4,2) [1024, 36], (4,3,3) [1296, 54]
  costs minutes to an hour and is the cheaper next test of r* = product value.
