# G26-P12 notes (adversarial review of P8 Thm U, Thm A2; work on (A))
03:24 start. Read P6/P7/P8 notes, verify_am2.py.
03:27 Thm U hand check: Moebius on cols -> N[S,U]=sum_{f:S->>U} prod E_{x f(x)} OK (multilinearization). Partial fractions
 1/prod_{x in B}(X_x-Y) = sum_x A_x/(X_x-Y), A_x=prod_{x'!=x} 1/(X_x'-X_x) (deg num 0 < |B|) OK => gamma_{x,B}=prod c/(X_x'-X_x) OK.
 (f,reps) <-> (R,g:R->U bij, phi:S\R->R) bijection OK; Gamma_{S,R} indep of U; Gamma_{S,S}=1; [R subset S]. So Gamma is lower
 unitriangular (any linear extension of inclusion), entries in O under H: v2(c)>=v2(X_x-X_x') all x!=x'; det=1 => GL(O).
 Smith(Gamma P)=Smith(P), P block diag => union. Needed pairs: ALL distinct row pairs of B_b; v2(x-x')<=a-1, v2(c)=a+1+v2 zeta>=a+1
 for zeta in O => v2(gamma factor)>=2 OK (needs zeta in O; zeta with v2<-1 can break H).
 checkU_exact.py: exact identity N = Gamma*(+)Per_k with Gamma from the formula, 12 random rational Cauchy matrices n<=4: 0 mismatches.
 checkU_smith.py a=2 (84 cases), a=3 (196 cases; zeta in {1,-3,1/3,2,-12,5/7,8}, all t, b): U, (A), C1'-right all TRUE.
 negU.py: random Cauchy-like over Z_2, n<=5: under H 78/78 Smith(M)=union; without H 17/32 differ (control discriminates).
03:30 Thm A2 hand check: same-class E = c2^-i/(u+2(j-k)) (u=(t+alpha zeta)/2^i odd). F{1}=u l0-(u+2)l1: class part = c2^-i*4/(u(u-2)) X_{rho1}
 (val w0+2); off-class column, m=v2(rho-sigma)<i: coeff c*2(u2^i-D)/(D(D+2^{i+1})), val v2c+1-m >= w0+2, equality iff m=i-1 =>
 reduction X_{rho1}+nu_{rho xor 2^{i-1}} CONFIRMED. F{0,1}=l0l1+2^{w0}l0 = 2^{2w0}(lam0 lam1+lam0); lam_j==nu mod 2 => nu^2+nu=0 Boolean
 => div by 2^{2w0+1}; (lam0+lam1)/2 == nu mod 2 [class part (1/u+1/(u+2)) val>=2; off-class m=i-1: v2(D0+D1)=i -> val 2] => the
 quadratic part nu*(mu0+mu1) collapses to nu (linear) => reduction X_{rho0}X_{rho1}+linear CONFIRMED. Basis: gr = exterior algebra,
 linear forms A_rho=nu_rho, X_{rho1}: factors A, X1+A_{rho'}, A^X1: unitriangular in #A's CONFIRMED. Sublattice lemma OK
 (Lambda' c Lambda => s_k(Lambda) | s_k(Lambda')). tensor-sum monotonicity OK. => Thm A2 CONFIRMED.
 checkA2.py sym: factor claims (exact valuations, reduction supports, F3 quadratic reduction = X_{rho0}X_{rho1} only) for a=3..8,
 8 t's, 7 zetas, b=0,1: 336/336 OK. checkA2.py direct a=2,3 (112 cases): Smith(M) == {0,w0,w0+2,2w0+1}^{(x)2^i} EXACTLY (stronger than P8).
 a=4 (checkU_smith, 31 cases t in {-8,-6,-3,1,2,4,12}, zeta in {1,-5/3,2}, b): U, (A), C1'-right all TRUE.
VERDICT: Thm U CONFIRMED (hypothesis H: v2(c) >= v2(X_x-X_x') all row pairs, true iff v2 zeta >= -2; holds for any Y).
         Thm A2 CONFIRMED (in fact equality at a=2,3). Corollary (sub-block lower bound) CONFIRMED (same computation, u odd).
03:58 KEY OBSERVATION (new, P12): Thm U implies rowspace(N) = rowspace(Gamma P) = rowspace(P) = (+)_k rowspace(Per_k(E)), graded by
 column degree |U|. So for ANY F in Lambda(M) (Moebius coords = multilinear coefficients), each homogeneous component F^{(k)} lies
 in Lambda(Per_k(E)), with v(F^{(k)}) >= v(F) and reduction = degree-k part of reduction of F.
 => (A) at level a-2 (all k): P8's F_T has X-degree <= |T| with top reduction tau_T (gr basis, graded); F_T^{(|T|)} in Lambda_k,
    weight w(T) <= sum_{j in T} s_j, reductions {tau_T: |T|=k} = basis of wedge^k. Sublattice lemma => Smith(Per_k) <= k-subset sums
    = Smith(C_k E). PROVED (modulo P8 construction, which I re-verified by hand). Also level a-1 trivial (Per_k(2G) = 2^k Per_k(G),
    Per_k(G) mod 2 = C_k(G mod 2) invertible).
 Lemma W'' (criterion, all levels): per class rho, E_rho = E[rho rows, all cols]; Thm U also holds for rectangular E_rho. If there are
    F_{rho,T} in Lambda(Per_{|T|}(E_rho)), weight w(T) <= sum_{j in T} s_j, reductions tau_{rho,T} in wedge^{|T|}(F_2^n) with
    {wedge_rho tau_{rho,T_rho}} a basis of wedge(F_2^n), then (A) for all k and C1'-right hold (sqfree product = wedge mod 2).
03:43 OBSTRUCTION (obstr.py, new): for D=2^{a-i}>=8 and i>=1, rows of ONE residue class (all columns) have Smith
 [s0,s0+2,s0+4,s0+7] (e.g. [4,6,8,11] vs class profile [4,6,10,12]) and their reductions at valuation <= s_j (j>=2) lie in the
 GLOBAL V(s_{j-1}): class-local row combinations give NO new Smith direction beyond mode 1 (cases a=4,5,6, i=1,2,3, zeta=1,2,3).
 => High Smith modes of E (j>=2) need rows from several classes; P8's class-product construction (Thm A2) cannot extend to
 level a-3 or below when i>=1 (it works at a-2 only because modes 0,1 are class-local). i=0 has no issue (single class).
 levelm3.py naive criterion at a=3,i=0 fails for T=(0,3),(1,2) (tie in index sum; needs block version) - abandoned.
 New idea: induction on i at fixed D via 2-adic split rows even/odd: E_{a,i} = [[E0, cZ],[cZ', E1]], E0,E1 ~ unit*E_{a-1,i-1}.
 Tensor-domination TD: Smith(M(E_{a,i})) <= Smith(M(E_{a-1,i-1}))^{(x)2} termwise ? If true, C1'-right at level m=a-i for all i
 reduces to the single case i=0 (finite check per D).
03:50 td.py/td2.py: TENSOR DOMINATION FAILS both ways for (a,i)=(4,1): per k (k=3,4) and at M level (union): e.g. sorted pos 77:
 Smith(M_{4,1})=22 > tensor 21, pos 88: 22 < 23. (3,1) and (4,2) equal. => induction on i at fixed D via even/odd split is
 NOT available as a termwise statement; off-class coupling changes Smith both ways. (C1'-right still holds there.)
04:00 obstr2.py: modes j>=2 at level a-3 need rows from ALL classes: (5,4) i=2: group {0,2} gives [2,2,1,1] new dirs (want 2 each),
 all 4 classes [4,4,4,4]; (6,8) i=3: classes {0,2,4,6} give [4,4,4,3] (want 4), all 8 classes full; (6,4) D=16: {0,2} [2,2,1,1,0,0,0,0].
 => Smith-adapted modes of E are global across the 2-adic class tree for D>=8: no class-product structure.
 topper.py (Ryser, exact mod 2^200): v2(per E) - n(v2 c - i) = 2^i (D/2 - 1) = 2^i v2((D/2)!) in all 14 cases a<=5 (n<=16),
 vs det: 2^i sum_{j<D/2}(2j+2v2 j!). (A) for k=n holds with huge slack; exact formula conjectural (not proved).
 fastA: a=6,t=1,k=2 (496x496): (A) TRUE, min slack 1.
03:50 gradetest.py: random F in Lambda(M) -> Moebius -> degree-k part IS in rowspace(Per_k) (Smith-index membership) for blocks
 a=2,3 (t=1,2,4) and for a Cauchy-like E with H (c=8); FAILS for a control without H (c=1). Grading claim confirmed.
 fastA: (A) a=5 k=3 (560x560) TRUE for t=1 (zeta=1) and t=3 (zeta=-1/3), min slack 5. a=6 k=2 t=1 TRUE (min slack 1).
 sym2.py: Sym_2 := rowspace[Per_2(E); 2K], K[r,{u,u'}]=E_ru E_ru' (= square-free columns of Sym^2(E)). Per_2 != Sym_2 in general
 (a=3,i=0: [9,10,13,14,17,20] vs [9,10,13,14,16,19]), so 2K not inside Lambda_2.

=== THEOREM WA2 (weak (A)_2) — PROVED 03:50 ===
(i) For ANY n x n E over O with Smith s_1..s_n: R_p(Per_2 E) >= #{a<b: s_a+s_b<p} - R_p(2K), K[r,{u,u'}] = E_ru E_ru'.
(ii) For a crossing block (Cauchy, H): R_p(2K) <= R_{p-3-v2 zeta}(E), so
     R_p(Per_2 E) >= #{a<b: s_a+s_b < p} - #{g: s_g < p-3-v2 zeta}.
Proof. (i) Let L_a = sum_r P^-1_{ar} l_r = 2^{s_a} q_a (Smith-adapted, q_a mod 2 a basis of F_2^n). sqfree(L_a L_b) =
 sum_{r<r'}(a_r b_r'+a_r' b_r) Per2-row_{rr'} + sum_r a_r b_r sqfree(l_r^2), sqfree(l_r^2) = 2K_r. So Omega := span{sqfree(L_aL_b), a<b}
 lies in Sym_2 = Lambda_2 + span{2K_r}. Omega = span{2^{s_a+s_b} sqfree(q_a q_b)}, reductions q_a ^ q_b (sqfree product of linear
 forms mod 2 = wedge) are a basis of wedge^2 => Smith(Omega) = pair sums, and Smith(Sym_2) <= Smith(Omega) termwise (sublattice).
 Minimal number of generators mod 2^p: R_p(Sym_2) <= R_p(Lambda_2) + R_p(span 2K). Combine.
(ii) E_ru E_ru' = [c/(Y_u-Y_u')](E_ru - E_ru'), so K = E Gamma_2 with Gamma_2 entries c/(Y_u-Y_u') of valuation
 >= (a+1+v zeta)-(a-1) = 2+v zeta; 2K = 2^{3+v zeta} E Gamma' (Gamma' integral) => s_j(2K) >= s_j(E)+3+v zeta.  QED
 Numerics (weakA2.py, 9 cases a<=5): all three components and the final inequality hold. (A)_2 itself = WA2 without the loss term.
 The loss is O(R_p(E)) vs main term ~R_p(E)^2/2: (A)_2 holds up to a relative O(1/R_p) error, unconditionally, at all levels.
 (Does NOT extend cheaply to all k: the non-square-free generators then number ~ n*C(n,k-2).)
03:52 checkCor.py: P8 corollary factor claims on sub-blocks (x',y' in {0,1}) at levels i<=a-3, a=4..7: 126/126 OK.
!!! 03:52 COUNTEREXAMPLE TO (A): a=6, t=2 (i=1), zeta=1, b=0, k=2 (Per_2 is 496x496): top Smith invariant of Per_2 = 115 >
 114 = top pair sum (58+56); R_115(Per_2)=494 < R_115(C_2)=495. Found by fastA.py, CONFIRMED by independent elimination
 (lib12.smith_vals, verifyA62.py, NB=200). Consistent with WA2 (loss <= R_112(E)=32).
03:56 (A) violation at a=6,i=1 is robust: b=1 (NB=300), t=-2 zeta=-1/3, t=6 zeta=3 -> same: top Per_2 invariant 115 > 114.
 fastA a=6 k=2: t=1 (i=0) TRUE slack>=1; t=4 (i=2) TRUE (318 ties); t=8 (i=3) TRUE (432 ties). a=5 k=3: t=1,3 (i=0) TRUE, t=2 (i=1) TRUE.
 subsearch.py: 8x8 Cauchy sub-blocks E[R0,C0] of the a=6,t=2 block violate (A) for k=2 (and k=3) in 3 of 6 random draws;
 subC1.py: for those 3, Smith(M)=union Per_k (Thm U again exact) and C1'-right STILL holds (compensation across k).
!!! subsearchC1.py: C1'-right (union level) VIOLATED by an 8x8 Cauchy sub-block of the a=6,t=4 (i=2) crossing block:
 rows [0,2,3,9,18,33,35,40], cols [10,24,26,33,34,35,51,56], Smith(E')=[5,5,5,5,8,9,9,12]: R_15(M(E'))=26 < 27 = sqfree count
 (sorted pos 26: 15 > 14). Verified by direct Smith of M (256x256, lib12.smith_vals, verC1sub.py); H holds; C1'-left holds.
 => C1'-right is NOT a consequence of Cauchy structure + H; it fails for Cauchy submatrices (real-port-type sub-blocks).
    Any proof must use the full crossing-block index set. (A) is false for crossing blocks themselves (a=6,i=1).
04:00 subsearchC1.py (330 s each, random 8x8 / 9x9 Cauchy sub-blocks): a=5,t=2: 156 tried, (A) violated 54, C1'-right violated 0;
 a=6,t=2 (n'=8): 157/31/0; a=6,t=2 (n'=9): 25/6/0; a=6,t=4 (i=2): 156 tried, (A) 32, C1'-right violated in 6 (list in sC1_6_4_8_14.txt).

=== FINAL STATUS (P12) ===
VERDICTS: Thm U CONFIRMED (hyp. H: v2(c) >= v2(X_x - X_x') for all row pairs; any Y; also rectangular E).
          Thm A2 CONFIRMED (and Smith(M) = {0,w0,w0+2,2w0+1}^{(x)2^i} exactly at a=2,3). P8 corollary CONFIRMED.
          "C1'-right <= (A)" valid as implication, but (A) REFUTED: a=6, i=1 (t=2,-2,6; zeta=1,-1/3,3; b=0,1), k=2.
PROVED here: grading rowspace(N) = (+)_k rowspace(Per_k) (from Thm U) => (A) at level a-2 for all k (via P8's F_T top parts);
          (A) at level a-1 trivial; Thm WA2: R_p(Per_2 E) >= #{a<b: s_a+s_b<p} - R_{p-3-v2 zeta}(E) (all levels).
NEGATIVE: per-class modes obstruction (D>=8, i>=1); tensor domination fails (4,1); (A) false (a=6,i=1);
          C1'-right false for Cauchy sub-blocks (a=6,i=2, 8x8). C1'-right for FULL crossing blocks with a>=5: untested/open.
CONJ (numerics, 14 cases): v2(per E) - n(v2 c - i) = 2^i (D/2 - 1).
