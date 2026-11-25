# G27-B1 notes (brainstorm + quick kill-tests). Start 07:38 UTC, deadline 09:38.
Tags [P] short proof here (unreviewed), [E] exact, [N] numeric, [H] heuristic.
Convention as gen26: K = sum_l G_l (x) G'_l, Theta = K S, p = rank(1-Theta), q = osr(K).

## KT1 (07:58) Kronecker / substitution of exact optima  [E kron.py, mod 1000003]  -> DEAD
 Theta = Theta1 (x) Theta2, q = q1 q2: (8,6)@6^2 -> (496,36)@36; (8,6)x(11,7) -> (692,42)@42; (11,7)^2 -> (952,49)@49;
 (20,3)^2 -> (1416,9)@49; p/N^2 = 0.38..0.59, i.e. p = n1n2 - m1m2 - a1a2 (eigenvalue-1 pairs + (-1,-1) pairs) = naive
 p1n2+p2n1-p1p2 minus a1a2.  k-fold powers: fraction of eigenvalue 1 -> ((m+a)/n)^k/2... -> p/N^2 -> 1 (theta -> 2).
 [P] twisted version (level-2 fillers depending on the level-1 term, Theta = sum_l P_l (x) Theta2^(l), any fillers):
 a level-1 clique cell z has Theta1_zz = 0, so (1-Theta) restricted to {z} x [n2] is the identity: each level-1 clique block costs
 >= n2, p >= p1 n2.  So ANY substitution into an MC skeleton (gauge twists, crossed products) has theta >= (theta1 log N1 + 2 log N2)/log N.
 NB (496,36)@36 lies exactly on the finite-size circle (T=1, theta=sqrt3) -- coincidence; finite-size circle comparisons are
 meaningless anyway (App B trivial for N <= 256).

## Ideas (mechanism / escape / kill-test).  Machine shared, load ~18 on 4 cores: tests kept small.
I1 ORBIT (twirl) ansatz, non-monomial seeds: K = sum_{g in Gamma} c_g (L_g A P_g) (x) (L_g B P'_g), Gamma regular on [N], q <= |Gamma| R.
   Generalises the Q8 point (ST Q8 array = conjugation orbit of one seed) but with DENSE seeds (not monomial => escapes G23a's
   anisotropy obstruction, which is about monomial ZT arrays) and the symmetry permutes terms (escapes Theorem G).  Regime q ~ N
   where Lemmas 1-3 are vacuous (Q >= N) and the circle sits at theta = sqrt3.  Test: VP (orbit.py) at N = 8 (Q8 sanity), 9, 16.
I2 Kronecker / substitution into an exact MC skeleton with gauge twists: DEAD (KT1, block-Schur [P]).
I3 Unipotent Theta (1-Theta nilpotent of rank p; Theta_s then has a symmetric two-sided spread).  tr Theta = n forces sum tr(G_l G'_l) = N^2,
   so Q >= N.  In MC: every cell fixed exactly once => q >= N maps with graphs of a Latin square; for the group Latin square
   (sigma_g = L_g, P_g(u,v) = (gv, g^{-1}u)) the map preserves e = u+v and acts by reflections d -> 2g-d on d = u-v:
   abelian case = Theorem G (pairs {chi, chi-bar} each cost 1, p >= N(N-1)/2) [P].  Non-abelian / e-dependent signs: not tested.
I4 ST/MC UNBALANCED regime q ~ N: ST counting (G22b Thm 3) allows (p,q) = (N^1.5/2.8, N) (theta = 1.5 < sqrt3 at T = 1), and the
   clique-block cubic bound is vacuous for q >= N (sum c^3 <= q^2 n <= q^3 N).  Nobody searched unbalanced ST beyond tiling.
   Test: SAT (g22b shape.py) at N = 12, q = 12, p = 24, 30 (tiling-type references: (45,27)@12 from (5,3)@4 tiles).
I5 Frame decompositions phi = sum_i f_i f_i^T (M > N, overcomplete) of both EPR pairs, terms (i,j) split A|B: R = (sum_A a a^T) S,
   a_ij = f_i (x) g_j: PSD sums => p = dim span_A, q = dim span_B, p + q >= N^2 [P]. DEAD (positive). Signed frames
   (phi = sum eps_i f_i f_i^T) need cancellation = two-sided spectrum again: this is a reformulation, not a construction.
I6 q-deformed Heisenberg / quantum groups at roots of unity: graded dims unchanged (same character), locality multiplier
   prod_k (1 - w^k x)^{m_k} has the same degree H: B6 Lemma 1-D (thresholding optimal pair-by-pair for ANY multiplier of
   degree H) => same W + W' = 4h - 1. DEAD [P via B6].
I7 Signed / super swap X_eps = eps(x,y)[x=v][y=u]: X = CZ_eps^{-1} X_eps, CZ acts on parties x,y: preserves rank(xy|uv), multiplies
   osr by osr(eps) (= rank of eps as x|y matrix).  So a good X_eps transfers with factor rank(eps) on q only: constant for
   super-signs (rank 2); for eps of rank N^d, loses N^d.  No reason X_eps is easier.  DEAD as an exponent route [P].
I8 SU(2) recoupling / Racah uncertainty: S = sum_J (-1)^{2j-J} P_J, C = poly of degree d in J1.J2 (osr (d+1)^2), R = non-interpolated
   J's: sqrt p + sqrt q ~ N (theta, T -> 2). DEAD (= G25-B1 Casimir idea).
I9 Two-sided WEIGHTED clique blocks (blocks I - u v^T with big |u_i v_j|, some with eigenvalue 1-c < 0, some with big positive
   eigenvalue to cancel odd moments). [P, counting]: negative blocks of size c on ~n cells give tr Theta^3 ~ -n c^2; the positive
   side contributes <= ||Theta||^2 * (its trace) <= Q^2 * QN, and trace budget QN: so n c^2 <~ Q^3 N for ANY sign pattern =>
   p >~ N^{2.5}/Q^{1.5}: two-sidedness alone buys a constant (<= 1.5); the resource is Q (gauge-minimal norm).  Covering
   (c <= q+1 nonzero off-diagonal entries per map row) still gives the line.  So the target is: Q = N^{T'} with T' > T, i.e.
   genuinely non-normal factors -- no simple combinatorial source known.  No test (analysis only).
I10 Involutive Theta (Theta^2 = 1 <=> R S R = 2R): p = (n - tr Theta)/2 <= (n - QN)/2 bound; in MC blocks are size 2 only.
   Only useful for q >= N; subsumed by I1/I4 tests.

## KT2 (08:00-08:20) ORBIT ansatz, dense seeds  [N; orbit.py, orbit2.py (diag twists), orbit3.py (fast), hada.py (Z2^k Fourier form)]
 Sanity: Q8 conj orbit (R=1, q=8) finds (8,8)@8 from random starts (tail 1e-17, gap 4e9).
 NEW vs monomial: dense seeds give (8,8)@8 ALSO for D4 (3/6 seeds) and Z2^3 (3/6 seeds), where the monomial group ansatz
   (G22b transl.py / G23a ZT) is UNSAT. Z8: 0/6, Z4xZ2: 0/6 (tail 1.6 best).  So the monomial anisotropy obstruction is not
   the end of group-orbit families.
 Z2^k reduction [P, checked numerically == orbit3]: conj orbit, c = 1, R seeds: S - K ~ (+)_{kappa in F2^k} M_kappa,
   M_kappa = I - sum_r (X_r P_kappa) o (P_kappa Y_r)  (Hadamard product, X = N F A F, F Walsh-Hadamard).  p = sum rank M_kappa, q <= RN.
   N blocks of size N: the search scales to N = 32, 64 cheaply.
 Structure of the Z2^3 solution: X, Y are +-1 matrices (x 1/sqrt2) with Y = X^T up to sign gauge; each Theta_kappa = I - vv^T/2
   (eigen -3 = the Q8 clique, now in Fourier space).  SIGN ANSATZ [P]: Theta_kappa = I - v v^T/2 for all kappa  <=>  linear system over
   F2: x[i,j^k]+x[j,i^k]+x[0,i^k]+x[i,k]+x[0,j^k]+x[j,k] = 1 (0<i<j).  f2sys.py: N = 4, 8 solvable; N = 16 (1680 eqs, rank 189)
   and N = 32 INCONSISTENT [E].  => (N,N)@N in this sign ansatz exists exactly for N <= 8 (Hurwitz pattern again).
 N = 16 VP (hada.py, random starts): R=1 (q<=16): p = 16, 24, 32, 48 NOT found (tails 50-150, stationary points of rank 128, 144);
   R=2 (q<=32): p = 32 found 3/4 seeds instantly (blocks of rank 2) = the known (32,32)@16 (N^2/8);  p = 24, 16 not found.

## Coordinator's CORNER target (08:12): q = N^eta, p = N^{2-delta}, need delta >= c eta.
 Corner status of my ideas:
 - I2 Kronecker/substitution: no corner (p >= p1 n2: delta shrinks under products) [P].
 - I9 weighted (two-sided) clique blocks [P, model counting]: blocks of size c ~ N^delta are needed for delta > 0 (c = 2 blocks give
   p ~ n/2); the cubic budget n c^2 <~ Q^3 N then forces Q >= N^{(1+2 delta)/3}: in this block model the corner needs a JUMP to
   Q >= N^{1/3} (not just Q >= N^delta), and covering needs q >= c - 1 ~ N^delta (fine).  So the block model has a corner only if
   gauge-minimal norms N^{1/3} are available at q = N^{eta}, eta small: no known source.
 - I1 orbit ansatz: lives at q >= |Gamma|; corner would need small groups (|Gamma| = N^eta, free action, N^{1-eta} orbits).  Small test
   at N = 8 (osmall.out) vs the exact small-q frontier.
 - I6/I8 (graded, Casimir): App-B-type corner delta = eta^2/4 or worse.
 Further orbit results [N]: N=16 nonabelian Q8xZ2, R=1 (q=16), p=32: tails 70-89 (not found).  N=12 D12/Q12, R=1 (q=12), p=18, 24:
   not found (tails 16-50).  N=9 Z9, Z3xZ3, R=1 (q=9), p = 11, 14, 18: not found (best tail 0.8 at p=18).
   Small groups acting freely (corner direction): N=8, Z2 (k=4 orbits), q=2: p=28 only borderline (tail 1e-5); Z2xZ2 (k=2), q=4:
   p=24 found (= known (24,4)@8), p=22 not found.  => the orbit ansatz reproduces known points ((8,8)@8, (32,32)@16 via Z2^4 R=2,
   (24,4)@8) and finds nothing new.  Honest caveat: the D4 / Z2^3 dense (8,8)@8 solutions are presumably the octonion point seen
   through other symmetry groups (octonions are Z2^3-graded), so "escaping the monomial obstruction" here gains nothing new.
 VERDICT I1: no new point; Z2^k sign version proved impossible beyond N=8 (F2 system).  Remaining untested: nonabelian groups
   with R >= 2 at N = 16..32 (expensive), and complex/phase seeds.

## KT5 (08:20) unipotent Theta (I3): D = S - K, rank D <= p, D S D = 0 (Y = 1 - Theta, Y^2 = 0)  [N, unip.py, N=4]
 Feasible: (8,7),(7,7),(6,8),(5,10),(4,12)@4 (tails 1e-17); NOT: (5,4),(6,4) (tails 6-17), (8,6) (3e-2).  General frontier at N=4:
 (6,1),(5,3),(4,4).  So index-2 unipotent decompositions sit near the trivial antidiagonal p+q = N^2 (14-16 vs 16): the nilpotent
 (Jordan) structure is very costly.  DEAD (small-N evidence; not a proof).
 [P] small lemma for the Z2^k SIGN ansatz (R=1, X +-1, Y = X^T scaled 1/2): each block is D_kappa = T_kappa/2 with T_kappa a symmetric
   +-1 matrix with unit diagonal; tr T = N, ||T||_F^2 = N^2 => rank 1 (T = vv^T) or rank >= 3 (rank 2 would need l1 l2 = 0).
   So in the sign ansatz q = N gives p = N (only N <= 8, F2 system) or p >= 3N.  (48,16)@16 would need rank-3 blocks: VP below.

## 08:22-08:30 NEW numerical points in the Z2^k orbit ansatz (R = 1, q <= N)  [N, hada.py]
 N=16, target p=64: 5/10 seeds converge to rank 62 (block ranks 2^8, 6^7, 4^1), tail 1e-19, gap 5.8e9.
   Re-built in the ORIGINAL basis (A = F X F / N, B = F Y F, orbit3 build): rank(S-K) = 62 (sv 1.46 vs 4.5e-10), osr(K) = 14 numerically
   (15th Schmidt value 1.5e-8 / 71.8: probably exact 14, to be confirmed).  => (62, 14..16)@16  [N, not yet exact].
   Spectral profile: Theta non-integer two-sided spectrum (-4.97, 4.27, -2.94, ...), ||Ts|| = 20, Q = sum||G||||H|| = 125 >> q:
   exactly the "lesson (b)" profile (large norms, spread spectrum), unlike every N <= 9 optimum (integer spectrum, Q = q).
   Structure: index set splits into the index-2 subgroup U = {0..3, 8..11} and its coset V; X o Y^T = I on VxV, 1/2 on UxV, VxU
   (sign-ansatz part), continuous on UxU.
   Benchmarks at N=16 for q <= 16: Kronecker (8,8)@8 x (1,1)@2 = (88,8)@16 (spectral formula); 2+2 direct-sum split with q=16: 144;
   q = 1: 120.  So (62,<=16)@16 would be a new frontier point at N = 16 (if exact).
 N=32, R=1 (q <= 32), target 256: rank 255 found (seed 0).  Benchmark: (128,128)@32 (balanced), Kronecker (8,8)x(4,4) ~ 350.
 N=16 target 60: rank 138 only (not found).
 08:30-08:37 follow-up [N]:
  * seeds 6, 8, 9 of the same run have osr EXACTLY 12 numerically (Schmidt values 0.99..0.98 then 1e-10): (62,12)@16, sv_62 = 0.009..0.13,
    sv_63 = 2e-10 (gap >= 4e7), entries bounded (|x| <= 1.3).  Theta: mostly -3 (Q8-type cliques) + a few non-integer eigenvalues
    (-5.48, -3.16, -3.13, ...); #eig(Ts) >= 1 = 194 = n - p (tight); Q = 17 (> q = 12), ||Ts|| = 5.5.
  * Kronecker benchmarks at N = 16 [E mod 1000003, kron2.py]: (4,4)x(4,4) = (108,16), (4,4)x(5,3) = (120,12), (5,3)x(5,3) = (131,9);
    (8,8)@8 x (1,1)@2 = (88,8) (spectral count).  So (62,12)@16 beats every construction I know with q <= 12 at N = 16 (and the
    balanced (32,32)@16 is at q = 32).  Finite-size only (theta = 1.49, T = 0.90 at N = 16); says nothing about the corner.
  * (61,12): warm start + LBFGS + Gauss-Newton stall at sv_62 ~ 8e-9 with bounded entries: UNRESOLVED (degenerate/border-like).
    (60,.): tail 1e-8 with entries growing to 24: border-like.
  * EXACTNESS NOT ACHIEVED: progressive rounding (roundh.py, denominators 4, 8, 12, 16) fixes only ~190-230 of 512 entries; the rest
    are not near simple rationals (probably algebraic).  exacth.py (flint check of a rational X, Y) is ready for when a rational point
    is found.  Numerical solutions saved: num_16_62_12.json (seed 6), num_16_62_14.json (seed 1), h_N16_R1_p64_s*.npy.
  * N = 32, R = 1: "rank 255" at target 256 has NO clean gap (sv 3e-5..3e-7): unresolved; target 192: not found.
 ST SAT (KT3): (22,4)@8, (20,5)@8 TIMEOUT at 280 s each (machine load ~15); N=12 (24,12), (30,12): TIMEOUT 420 s. Undecided.
 Tied sub-ansatz Y = X^T/2 (hadaT.py): N=16 p=64 found 2/8 seeds (tail 1e-15, not clean); rounding blocked at tolerance. Not pursued.

## RANKING / VERDICT (08:40)
 1. I1 orbit ansatz, Z2^k Fourier form (hada.py): only live lead.  Gives (62,12)@16 [N] (beats all known q <= 12 constructions at N=16:
    Kronecker (120,12), (88,8)); large-norm, non-integer two-sided spectrum => first numerical point with the "lesson (b)" profile at
    moderate N.  NOT exact yet.  Sign version provably stops at N = 8 (F2 system).  Next: exactify (algebraic point, or impose extra
    symmetry U/V split + octonion part), then N = 32 with R = 1, 2 and longer runs / more seeds (cheap: 32 blocks of 32x32).
 2. I4 unbalanced ST/MC at q ~ N (ST counting allows theta = 1.5 at T = 1; cubic clique bound vacuous): SAT undecided (load);
    needs a dedicated SAT run (cube-and-conquer) at N = 8 (20..22, 4..5) and N = 12 (24..30, 12).
 3. I9 counting (weighted blocks): corner in block-type models needs Q >= N^{(1+2 delta)/3} -- a jump; suggests the corner can only
    come from non-block (App-B-like) structures. [P, model-level, unreviewed]
 DEAD: I2 Kronecker/substitution [E+P]; I3 unipotent (small-N numerics, near p+q = N^2); I5 positive frames [P]; I6 q-deformation
    [P via B6 Lemma 1-D]; I7 signed swap [P]; I8 Racah/Casimir; I10 involutive (subsumed).
 All my processes killed (other python jobs on the machine belong to other agents).
