# G24f — the coset window (labels in V_{alpha/2+L}, A1^n, xi = alpha/2): true N vs flattening (working notes)
Tags: [P] proved here (short, unreviewed), [E] exact, [N] numerical, [H] heuristic, [O] open.
Nothing here is a new algorithm or a new bound for 3-SAT.

## 0. Setting
A1^n, xi = (alpha/2)^n, pairs in V' (g_xi weights n^2-n+k per copy, V'_0 = span{1,e^alpha} per copy), labels in
M = V_{alpha/2+L}, SHIFTED grading (bottom 0): copy in sector m (charge (m+1/2)alpha) at level k has weight m^2+k,
char = sum q^{m^2}/eta = 1+3q+4q^2+...  (M_0 1-dim, M_1 = {-,0,+} per copy).  S_xi = {0,1}^n.
At H = 0 the pair product is the algebra R = (k[e]/e^2)^{(x)n} (B_g(e^b,e^d) = +-[b,d disjoint, b+d=g]).
E0 = 1 explicitly [E, build.py]: o(e^{e_j}) = J_j : -_j -> 0_j -> -2 +_j (local);  o(e^{e_i+e_j}) = X_ij : -_i -> +_j,
-_j -> +_i (weight transfer between copies);  |g| >= 3 acts as 0.

## 1. Correction to the window itself  [N + P]
 G24b's "loose 2.416 / loose+flattening 2.674" are grid-edge values (grid v >= 0.03).  As u, v -> 0 (dilute regime,
 few excited copies) loose -> 2 (chk_loose.py: 2.19 at u = 0.003, v = 0), and the FLATTENING also -> 2: at E0 = 1, H = 0
 rk F = n^2 + 3n exactly for n = 1..5 [E, flatE1.py], i.e. N <= ~n, while the cost is 2 log(3n+1): 4rho >= 2 only.
 Dilute bookkeeping (fixed E0, H, n -> oo): dim M_E ~ n^E, cost K ~ (2E0+H) log n, N ~ n^kappa  =>  4rho -> (2E0+H)/kappa.
 Labels: kappa <= E0; so only H < (2sqrt2-2) E0 = 0.83 E0 can go below 2sqrt2; E0 = 1 needs H = 0 and kappa > 0.707.

## 2. Exact / numerical true N (small copies)   [E where stated, else N = LM least squares, many random starts]
 numsearch.py solves sum_g (A B_g Bm^T) (x) (U O_g Y) = I (x) I.   CAVEAT: "not found" is evidence, not proof; one-copy
 N = 2 at (E0,H)=(1,0) converges to residual 1e-6 with growing entries = BORDER solution only (true max 1, proved in G24b).
   copies (E0s; ha,hb)   true N            flattening floor sqrt(rkF)   single channel max min(R_g,Lab_g)
   1  (1; 0,0)           1  [P, G24b]       2                            1
   1  (2; 0,0)           2  [N, res 4e-16]  2
   1  (1; 0,1), (2; 0,1) 1  [N]             -
   1  (2; 1,1)           2  [N]             -
   2  (1; 0,0)           2  [E, verify_N2.py: a=(e^{a1},e^{a2}), b=(e^{a2},e^{a1}), y=(-1,-2), u=(+2,+1)]; N=3 not found (res 1.41)    flattening 3
   3  (1; 0,0)           2 found, N=3 not found (res 1.00)        flattening 4
   4  (1; 0,0)           >= 2 (the 2-copy construction embeds); 5 copies: >= 2   flattening 5, 6
   2  (2; 0,0)           >= 3 [N, res 3e-10, bounded entries; not rationalised]; N=4 not finished   flattening 5
   2  (1; 0,1)           2 [N]; N=3 not found (res 1.73)
 All jobs killed at the end of the session.
 So the flattening is not tight already at 2-3 copies (2 vs 3, 2 vs 4) and the gap grows with n (sect. 3).

## 3. Sharper bound: the a-first / b-first split  [P]
 Borcherds for weight-0 a, b (H = 0): o(a(-1)b) = sum_{j>=0} a_(-1-j) b_(-1+j) + sum_{j>=0} b_(-2-j) a_(j).
 The first sum passes through M_{E0-j} (j >= 0), the second through M_{E0-1-j}.  Hence for ANY V, M (any admissible grading):
   I (x) I = P_a + P_b,  rank_{(x,u)|(v,y)} P_a <= D := dim M_{<=E0},  rank_{(x,y)|(v,u)} P_b <= D' := dim M_{<E0}.
 (General H: first-type intermediate weights in [E0-ha+1, E0+hb], second-type in [E0-hb, E0+ha-1].)
 SWAP LEMMA [P]: for all E, C in Mat_N:  rank(Swap - E (x) C) >= N(N-1)/2  (tight at E = C = I; exact tests lemma_chk.py).
  Proof: R := (E (x) C^T)Swap, R(x(x)y) = Ey (x) C^T x, R^2 = X (x) Y with X = EC^T, Y = C^T E (same nonzero spectrum,
  E, C^T intertwine their generalized eigenspaces X_l, Y_l).  ker(1-R) lies in V = sum_{l m = 1} X_l (x) Y_m, R maps
  X_l(x)Y_m -> X_m(x)Y_l, so tr(R|V) = sum_{l = +-1} tr(E Q_l C^T P_l) = m_1 - m_{-1} <= N, dim V = sum m_l m_{1/l} <= N^2;
  eigenvalues of R on V are +-1, so dim ker(1-R) <= (dim V + tr R|V)/2 <= N(N+1)/2.  QED.
 THEOREM [P]: if dim M_{<E0} = 1 (and H = 0), then P_b = E(x,y)C(v,u) and rank(Swap - E (x) C^T) <= D gives
     N(N-1)/2 <= dim M_{<=E0}.
 Coset, E0s = 1, H = 0: N(N-1)/2 <= 3n+1, i.e. N <= (1 + sqrt(24n+9))/2 ~ sqrt(6n)   (flattening: sqrt(n^2+3n) ~ n).
 => dilute exponent kappa(1,0) <= 1/2 (flattening 1): 4rho >= 4 - o(1) for E0s = 1, H = 0; for E0s = 1, H >= 1 labels give
 4rho >= 3.  So the whole E0 = 1 strip of the window is CLOSED (true cost >= 3 > 2sqrt2), and there the true exponent is at
 most HALF the flattening exponent — this is the quantitative form of "M_0 is 1-dim while V'_0 is 2-dim": the b-first
 half of every zero mode factors through the 1-dim M_0.
 For E0 >= 2, M_{<E0} is not 1-dim (dim ~ n^{E0-1}); the lemma would be needed for Schmidt rank r = dim M_{<E0}:
 rank(Swap - K) for K of operator-Schmidt rank r.  NOT proved.  A bound of the form >= N^2/2 - c r N would give
 kappa(E0,0) <= max(E0/2, E0-1): closes E0 = 2, 3 at H = 0 but not E0 >= 4.  [O]
 (App. B shows no universal bound of the form D + D' >= N^{2-o(1)}: there D ~ D' ~ N^{sqrt2}.)

## 4. Verdict
 * No construction below 2sqrt2 exists in the data: best exact true N at small sizes is 2 (n = 2..5, E0s = 1), always
   strictly below the flattening value for n >= 2.
 * The window is LARGER than G24b reported on paper (loose and flattening both -> 2 in the dilute corner u, v -> 0), but
   its E0 = 1 strip is closed rigorously (4rho >= 3; at H = 0, >= 4).  Remaining open: E0 >= 2 with H < 0.83 E0
   (and the finite-u interior), where only the flattening (>= 2 in the corner, 2.674 at the G24b grid point) is proved.
 * Heuristic [H]: data and the E0 = 1 theorem suggest kappa(E0,0) ~ E0/2 (N ~ sqrt(dim M_E0)), which would put the whole
   window at >= 4; proving it needs the Schmidt-rank-r Swap lemma (sect. 3).  Not done.
Files: build.py (exact objects), numsearch.py (true-N search), verify_N2.py (exact N=2), flatE1.py (flattening),
 lemma_chk.py (Swap lemma tests), chk_loose.py (loose -> 2), s_*.out, objs.py/trank.py (copied from g24b).
