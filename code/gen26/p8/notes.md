# G26-P8 notes (agent, unreviewed)
02:53 start. Copied P3 appd_explicit.py, model2.py, run5.py. Observation: P3 run5 already ran row H with smith=True (nufac/4, leaf NOT halved)
 -> 2.25e5/1.52e6, but via fixed-point iteration lam <- (E sqrt(lam)/tau)^2 (4 its, 0.5% tol) at density 4.26 only. Recomputing with bisection.
02:56 run_p8.py = P3's optG (model2.run, pack=slack, Dmode=adapt, amode=slack, leaf=depth, kappa=lin) + 3 NM starts, D in {2,3,4},
 geometric bisection to 0.2%. Density set by patching model2.LOGN_OFF=log2(1+alpha), w=km/n=3 alpha. smith=True divides nufac by 4
 (Thm 2: n_t -> n_t/4 exactly); leaf brute force 2^{M_leaf} NOT halved; optimizer re-balances.
 Results so far: alpha=4.26 G vs 2^n 8.979e5 (P3: 8.97e5 OK); alpha=1 G vs 2^n 4.922e4.
 const_p8.py: E sqrt(lam)/(km): G 73.76/74.16/74.36(4.26), 74.44(alpha=1,1e7); H 37.03/37.18/37.19 (1e5,1e6,1e7); ratio G/H 1.992,1.994,1.999.
03:00 TASK 2: PROVED C1'-right at level i=a-2 (all zeta != 0, all t with v2(t)=a-2, b=0,1). Construction (verify_am2.py, 144 cases a=2,3,4 OK):
 classes rho mod 2^i, 2 rows each l_{rho,0}, l_{rho,1}; u=(t+alpha zeta)/2^i odd; w0=3+v2(zeta).
 F_{rho,{}}=1; F_{rho,{0}}=l0 (wt w0, red. nu_rho=X_{rho0}+X_{rho1}); F_{rho,{1}}=u l0-(u+2) l1 (wt w0+2, red. X_{rho1}+nu_{rho xor 2^{i-1}});
 F_{rho,{0,1}}=l0 l1+2^{w0} l0 (wt 2w0+1, red. X_{rho0}X_{rho1}+linear). F_T=prod_rho F_{rho,T_rho} in Lambda (disjoint rows).
 Reductions form a basis of Fun(F_2^n,F_2): degree filtration, gr = exterior algebra (char 2: perm=det), unitriangular in #nu's.
 Sublattice lemma => Smith(M(E)) <= {0,w0,w0+2,2w0+1}^{(x)2^i} <= square-free sums {0,w0,w0+2,2w0+2}^{(x)2^i} termwise.
03:07 core_explore.py: naive alternants sum_S det(Binv[T,S]) prod_S l_x do NOT reach w(T) (lower-degree parts drag): dead end.
03:08 percomp.py: (A) Smith(Per_k(E)) <= Smith(C_k(E)) = k-subset sums of Smith(E) termwise: TRUE for all k, (a,i)=(3,0),(3,1),(4,0),(4,1),(4,2) (zeta=1).
  (A) FAILS for generic E: E=[[3,4],[4,16]]: Smith(E)={0,5}, det v=5, per=64 v=6 > 5.  => Cauchy structure needed.
03:10 union_check.py: (U) Smith(M(E)) == disjoint union_k Smith(Per_k(E)) EXACTLY in all 8 P6 cases (a<=4, all i, zeta=1).
  (U) fails for generic E (P6's n=2,s=(0,5) example: Smith M=(0,0,1,6), union contains 5). 
  => C1'-right <= (U<=) + (A): reduces to the permanental compounds of a Cauchy matrix vs its determinantal compounds.
03:14 === THEOREM U (PROVED here) === For a Cauchy-type block E_xy = c/(X_x - Y_y) with c/(X_x - X_x') in O for all rows x != x'
 (crossing blocks of (D.7): X_x=T+x, Y_y=y, v2(c/(x'-x)) >= (a+1+v2 zeta)-(a-1) >= 2):
   Smith(M(E)) = disjoint union_{k=0..n} Smith(Per_k(E)),  Per_k(E)[S,U] = per E[S,U] (|S|=|U|=k).
 Proof. Columns: Moebius on the cube (unimodular) turns M into N[S,U] = sum_{f:S->>U} prod_x E_{x f(x)} (coeff of X^U in ml(prod_S l_x)).
 Partial fractions: prod_{x in B} E_xy = sum_{x in B} gamma_{x,B} E_xy, gamma_{x,B} = prod_{x' in B\x} c/(X_x' - X_x) in O.
 Expand prod_{y in U} over blocks f^{-1}(y), choose representatives: (f, reps) <-> (R subset S, |R|=|U|, bijection g:R->U, phi:S\R->R).
 => N[S,U] = sum_{R subset S, |R|=|U|} Gamma_{S,R} per E[R,U], Gamma_{S,R} = sum_{phi:S\R->R} prod_{r in R} gamma_{r,{r} u phi^{-1}(r)},
 independent of U; Gamma_{S,S}=1, Gamma_{S,R}=0 unless R subset S  => Gamma unitriangular in inclusion order, in GL(O).
 So N = Gamma * (direct sum_k Per_k(E)).  QED.  (Explains the exact numerical identity (U); generic E has no such Gamma.)
 Consequence: C1'-right  <=  (A) Smith(Per_k E) <= Smith(C_k E) (= k-subset sums of Smith E) termwise, for each k;
 the lower-degree correction problem disappears; (A) is a statement about permanental vs determinantal compounds of one Cauchy matrix.
03:13 TASK 1 results (run_p8.py, files thr_*.txt):
  alpha=4.26: G vs 2^n 8.979e5 | H vs 2^n 2.251e5 | H vs PPSZ 1.517e6 | G vs PPSZ: rerun thrG_ppsz.py (first job hit 590 s)
  alpha=1.0 : G vs 2^n 4.922e4 | G vs PPSZ 3.290e5 | H vs 2^n 1.285e4 | H vs PPSZ 8.291e4
  Check of reviewer's point: P3's run5.py ALREADY computed row H with the full model (model2.run smith=True: nufac/4, leaf 2^{M_leaf}
  unchanged) -> its 2.25e5/1.52e6 were computed, not scaled; bisection here reproduces them (2.251e5, 1.517e6).
  The leaf does not halve, but the optimizer pushes M_leaf down (depth 2-4), so lambda*_G/lambda*_H = 3.99 (4.26), 3.83/3.97 (alpha=1).
03:14 (A) numerics (union_check2.py): a=3 all i, t in {first,last with v2(t)=i}, zeta in {1,-3,1/3,2,-5/7}: (A) all k AND (U) direct (Smith(M)
  vs union) TRUE in 30/30 (uc2_a3.txt); a=4 all i, all k<=8: (A) TRUE 40/40 (uc2_a4.txt); a=5 all i, k<=2: (A) TRUE 50/50 (uc2_a5.txt).
  incl.py: (A) is NOT a lattice inclusion rowspace(C_k) c rowspace(Per_k) (fails a=4 all k; a=3 k=2) -> needs termwise argument.
  ldu_test.py: naive explicit candidate F = Per_k(L)^{-1} Per_k(E) (E=LDU, L integral unit lower-triangular, order (x',rho)) fails
  (reductions dependent) except at i=a-1. Repeated-mode (divided power) terms drag; (A) needs the Mahler structure. OPEN.
=== TASK 2 SUMMARY ===
[P] Thm U (general, any level): Smith(M(E)) = union_k Smith(Per_k(E)) for every crossing block (partial fractions => N = Gamma (+)Per_k,
    Gamma unitriangular over O). With P6 Prop 3: Smith(F) = {sigma+sigma'}, sigma in union_k Smith Per_k(E), sigma' likewise for E'.
    So the KF flattening rank of a packet = a count over permanental compounds of two Cauchy matrices; C1'-right <= (A).
[P] Thm A2 (level i=a-2, all zeta!=0, t, b): Smith(M(E)) <= {0,w0,w0+2,2w0+1}^{(x)2^i} <= square-free sums termwise, w0=3+v2(zeta)=s_0.
    => C1'-right holds at i=a-2 (and P6 had i=a-1).  Verified numerically 144 cases (verify_am2.py).
[P] Corollary (any level, via submatrix monotonicity s_k(A) | s_k(XAY)): taking x',y' in {0,1} in each residue class, the same
    construction applies verbatim (u odd, w0=s_0=a-i+1+v2 zeta, off-class valuations unchanged) =>
    R_p(M(E)) >= #{T: sum over classes of weights in {0,s_0,s_0+2,2s_0+1} < p}: unconditional but uses only modes j=0,1 of each class.
[N] (A) Smith(Per_k E) <= Smith(C_k E): 120 cases a<=5 (a=5 only k<=2), no failure; false for generic E. OPEN.
Article 3: the algorithm's nodes have m=a-i >= c_d+1 ~ 11-18 (i<u, a=a_prev+u+c), so the proved levels (m=1,2) are NOT the dominant
    ones; the sqrt2-optimality of mode enumeration in KF (square-free count exp(pi sqrt(np/12)) vs exp(pi sqrt(np/6))) remains
    conditional, now on the single clean statement (A) about permanental vs determinantal compounds of a Cauchy matrix.
    Unconditional (corollary): only exp(O(M/D))-size lower bound -> far from the algorithm.
03:23 TASK 1 FINAL TABLE (lambda* = log2 s; full P3 finite-size model, bisection 0.2%):
  row | alpha=4.26 vs 2^n | alpha=4.26 vs PPSZ | alpha=1 vs 2^n | alpha=1 vs PPSZ
  G   | 8.98e5            | 6.06e6             | 4.92e4         | 3.29e5
  H   | 2.25e5            | 1.52e6             | 1.29e4         | 8.29e4
  Constant E sqrt(lam)/(km), alpha=4.26: H 37.03 (1e5), 37.18 (1e6), 37.19 (1e7), 37.25 (1e8, D=4; D=3 gives 37.51);
  G 73.76, 74.16, 74.36, 74.48 (P3) ; G/H = 1.992, 1.994, 1.999, 2.000 -> closed form 37.5 = 75/2 approached from below.
  alpha=1 gives the same constants (37.04, 37.18, 37.23 at 1e5-1e7). constH 1e9 job timed out (no value).
No processes left running. Nothing written to /home/user/pnp-recherche.
