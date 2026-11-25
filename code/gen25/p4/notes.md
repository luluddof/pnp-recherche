# G25-P4: does Lemma T (triangularisation) close the H = 0 coset window?  Answer: NO, not by this route.
Tags: [P] proof here (short, unreviewed), [E] exact computation (flint, rationals), [N] numerical.

## 1. The operators in P_b  [P + E]
Setting of G24f (A1^n, xi = alpha/2, regraded V' with V'_0 = span{e^{alpha_S} : S in {0,1}^n}, labels M = V_{alpha/2+L},
shifted weight: copy in sector m (charge (m+1/2)alpha) at oscillator level k has weight m^2 + k).
Data of a construction: a_x, b_v in V'_0, Y_y in M_E0, U_u in M_E0^*  (x, v, u, y = 1..N).
Borcherds at H = 0:  o(a(-1)b) = sum_{j>=0} a_(-1-j) b_(-1+j)  +  sum_{j>=0} b_(-2-j) a_(j)   (on M_E0).
  P_b(x,v;u,y) = sum_{j=0}^{E0-1} < U_u , (b_v)_(-2-j) (a_x)_(j) Y_y >  =  sum_{m in basis of M_<E0} E_m(x,y) C_m(v,u),
  E_m(x,y) = < m^* , (a_x)_(E0-1-wt m) Y_y >          (a-mode lowers weight by j+1:  M_E0 -> M_{wt m})
  C_m(v,u) = < U_u , (b_v)_(wt m - E0 - 1) m >         (b-mode raises weight by j+1: M_{wt m} -> M_E0)
 E0 = 2: m in M_0 (1 state): a_(1), b_(-3);   m in M_1 (3n states): a_(0), b_(-2).          r = 1 + 3n
 E0 = 3: m in M_0, M_1, M_2: (a_(2), b_(-4)), (a_(1), b_(-3)), (a_(0), b_(-2)).              r = 1 + 3n + (4n + 9n(n-1)/2)
 E0 = 4: j = 0..3, m in M_0..M_3: (a_(3-w), b_(w-5)) for wt m = w.
 Explicit modes (cocycle = 1, all <alpha_i,alpha_j> even): Y(e^g,z) e^p (x) f = z^{<g,p>} E^-(g,z) E^+(g,z) e^{p+g} (x) f.
 * vacuum drops: 1_(j) = 0 (j >= 0) and 1_(-2-j) = 0, so only S != {} contributes (one copy: E_m = c r_m^T, rank 1, trivially
   triangular -- consistent with one-copy N <= 2);
 * selection rule [P]: e^{alpha_S}_(j) maps sector p -> p + alpha_S (total charge +|S|), and kills e^p (x) (level k) unless
   sum_{i in S} (2 m_i + 1) <= k - 1 - j.
 modes.py builds every mode e^g_(n): M_Ein -> M_Eout exactly (generalises G24f's O_apply).  borch.py: the split
 o(a(-1)b) = P_a + P_b holds EXACTLY on M_E0 for all a, b in {e^g} at (n,E0) = (1,2),(2,2),(2,3),(3,2),(1,4) [E].
 Dims (n = 2): M_0..M_4 = 1, 6, 17, 38, 84;  (n = 3): 1, 9, 39, 120.

## 2. Triangularisability: NO in general  [E]
Obstruction used (exact, field-independent): a family that is upper triangular after an equivalence spans a space of
dimension <= N(N+1)/2.  So dim span{E_m} > N(N+1)/2 => Lemma T's hypothesis fails for the E-side (same for C).
EXACT COUNTEREXAMPLES (exact_cert.py, all rational, identity checked exactly):
  N = n, a_x = e^{alpha_x}, b_v = e^{alpha_[n] - alpha_v}  (a_x(-1)b_v = delta_xv e^{alpha_[n]}),  Y = random integer
  columns in M_E0, U = exact left inverse of O_[n] Y.  Then sum_g beta_g(x,v) U O_g Y = delta_xv I_N exactly, and
   (n,E0,N) = (2,2,2): dim span E = 4, dim span C = 4 > 3        -> both sides non-triangularisable
   (2,3,2): 4, 4 > 3;   (2,4,2): 4, 4 > 3;                        -> both sides, E0 = 3, 4
   (3,2,3): 9 > 6 (E), C = 6 (dimension test inconclusive);  (3,3,3): 9, 9 > 6  -> both sides
 So for E0 = 2, 3, 4 there are genuine solutions whose E-family AND C-family are not simultaneously triangularisable by
 any (X, Y), over any field containing Q.
Numerical solutions from G24f's own search (sol.py, tri2.py, tri3.py) [N]: n=2,E0=2,N=2: 5/5 solutions have dim span = 4
 on both sides (min rel. singular value 0.014-0.5); n=3,E0=2,N=2: 4/5 non-triangularisable (one special solution is);
 n=2,E0=2,N=3 (5 sols, res <= 3e-9): C-span dim 7 > 6 in all, best triangularity defect 0.08-0.74 (relative, 60 starts);
 n=2,E0=3,N=3 (5 sols, res <= 1e-8; new data point, N=3 at E0=3): spans dim 9 on both sides.
 Conversely trisearch.py (identity + penalty on the strictly lower part, WLOG by the gauge X(x)Y) finds solutions WITH
 triangular E-family and with triangular C-family at n=2, E0=2, N=2 and N=3 (residuals ~1e-9 / 1e-11) [N].
 => triangularisability is a property of the individual solution, not of the vertex-algebra data.

Why the grading does not do it [P]: the modes are bihomogeneous in (weight, charge), but this grading lives on M, i.e. on
the SLICE index m and on the target; triangularity is needed in the index pair (x, y) = (pair label in V'_0, label in M_E0),
on which a solution acts by arbitrary linear maps (A', Y).  The grading only gives charge(m) = charge(Y_y) + alpha_{S(x)}:
each slice has "anti-diagonal" support, slices with different m have different anti-diagonals, and no common flag results.
Even for charge-homogeneous a_x (as in the counterexamples), mixing of sectors in Y destroys any order.

## 3. Consequence for the window
 * Lemma T cannot close the H = 0 window uniformly: its hypothesis fails for actual solutions at E0 = 2, 3, 4.
 * The counterexamples have tiny N (2, 3) where N(N-1)/2 <= dim M_<=E0 holds trivially; they do NOT show the bound
   N(N-1)/2 <= dim M_<=E0 is false -- only that this proof route does not work.  Possible remaining routes [O]: (i) show
   that solutions in the regime N(N-1)/2 > dim M_<=E0 are forced to have a triangularisable (or flag-defect o(N^2)) family
   (G24g sect. 6 flag form); nothing in the structure suggests this; (ii) a Schmidt-rank-r swap lemma using the E/C coupling.
 * Status unchanged from G24f: E0 = 1 strip closed (4rho >= 3, >= 4 at H = 0); open: E0 >= 2 with H < 0.83 E0 (only the
   flattening bound, -> 2 in the dilute corner).  No 4rho >= 2sqrt2 statement for the H = 0 window is proved here; H > 0
   untouched.
Files: modes.py, borch.py, exact_cert.py (+ log_exact.txt), sol.py, tri2.py, tri3.py, trisearch.py, s_*.pkl, log_*.txt.
All my background jobs killed.
