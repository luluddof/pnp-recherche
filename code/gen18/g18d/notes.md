# G18d — Limit of recursive coordinate-nice interpolation (OpenAI preprint, Appendix D)

Agent notes (AI; not an independent human review). Source read: paper.txt lines 3007-3655 (App. D;
the TeX file in the folder is a 14-byte "404", so the pdftotext version was used; superscripts as
confirmed by G17a). App. D is taken at face value; its correctness is NOT assessed.
Status legend: [P] proved here (paper proof, not Lean, not reviewed); [N] numerical; [C] conjecture;
[D] "design-level" = statement about the algorithm as written with free parameters, not a barrier.

## 1. The recursion, exactly (Task 1)   [intermediate save 1]

One call = Boolean tensor network with M wires, input precision p (ring A_p = O[V]/(V^l, 2^p)).
Parameters at depth d (D.4): lambda = ceil(log2 s); r = 64 groups per coordinate; eps = 2^-17;
u_d = floor(eps*lambda/2^d) retained valuation levels; a_d = 2*sum_{j<=d} u_j; L = 2^{a_d} slots per
packet; labels in [r]^{a_d}; D_arr = r^{2 a_d} pair classes; leaf size H = ceil(s/floor(sqrt(lambda))).
Constraints (D.6): M >= 2 D_arr^2 Delta, M >= L D_arr, 2L <= Delta = ceil(s^{9/10}).
What a level does:
 1. hash vertices to labels, pack wires in <= 3M/L packets of L slots (slot bit i = orientation in
    coordinate i);
 2. W(T) = Lagrange/Cauchy kernel on S = {-L..L-1}: W(0) = I, deg 2L-1, W(t) = unit * (shift by t),
    which is local in coordinate v2(t) (EXACT niceness, zero error at the node);
 3. F(T) (deg <= M(2L-1)), F(0) = contraction. Nodes t in S\{0} with v2(t) = i < u, multiplicity
    q_t = h(u-i), h = ceil(3M/u); data F(t + alpha z) mod (z^{q_t}, 2^{p_t}), alpha = 2^{a+1},
    p_t = p + 8 h D_t + 8M(a+1), D_t = L/2^i;
 4. each F(t + alpha zeta) (poly many zeta) = sum over mode tuples beta with sum w_Z beta_Z < p_t of
    prod_{g=1..r} Z_{g,beta}   (D.16); #tuples = exp(Theta(sqrt(n_t p_t))), n_t = 6M/D_t modes per
    weight (packets x tiles), so exp(O(M/sqrt u))   (D.15).
What is passed down: for each (t, zeta, beta, g) a child = Boolean tensor network on the ports of
group g (old vertices of group g + one new "packet tensor" with <= 2L ports per packet), with
<= 6M/r = 3M/32 < M/8 wires, at precision p_t (precisions only increase along a path).
Cost recurrence (calls):  T_d(M) <= poly(s) * r * exp(c M/sqrt(u_d)) * T_{d+1}(6M/r),  T(M <= H) = 2^H.
Unrolled: #calls = exp(c * sum_d M_d/sqrt(u_d)), M_d <= Cs 8^-d, u_d = Theta(lambda/2^d):
 sum = (s/sqrt(lambda)) * sum_d (sqrt2/8)^d: geometric, ROOT-DOMINATED. Depth used: log2(lambda)/6 + O(1)
 (until M <= H = s/sqrt(lambda)); leaves cost 2^H.
Where sqrt(log) comes from: it is ALREADY in the root call. Two balanced quantities:
   budget  p_t ~ M D_t/u   (2-adic loss of tapered Hermite recovery; taper gives the 1/u)
   supply  n_t ~ M/D_t     (modes per unit of 2-adic weight)
   #tuples = #{monomials of weight < p in n variables of each weight} = exp(Theta(sqrt(n p)))
           = exp(Theta(M/sqrt u)),  and u <= a <= log2(M)/(4 log2 r)  ((D.6)).
 The recursion does NOT compound savings (levels multiply tuple counts; exponents add); its only
 role is to evaluate each term on M/8-wire children instead of by brute force, i.e. to bring the
 leaf size down to s/sqrt(lambda) with constant r. The second reason for "sum of u_d <= O(log s)":
 a_d - u_d >= a_{d-1} + O(1) is needed so that D_t >= 2^{a_d-u_d} absorbs the inherited precision
 p ~ M_{d-1} 2^{a_{d-1}}/u_{d-1} (D.11); the paper wastes a factor 2 (a_d = a_{d-1} + 2u_d).

## [intermediate save 2] Numerics for Prop. B' (lp_taper.py, lp_taper_out.txt)   [N]
Optimal taper by LP + bisection, both counting regimes (saddle-point count, M -> infinity scaling,
evaluation loss omitted = favourable to the scheme). E* = min_taper max_i ln(#tuples_i)/M:
  a = u (no slack):  E* sqrt(u/nu) = 2.5, 3.2, 3.7, 4.0, 4.4, 4.7, 5.0 for u = 4, 8, 12, 16, 24, 32, 48
  a >= u + 4:        E* sqrt(u/nu) = 6.0 - 6.2 (stable; nu = 2 or 6)
So E* = Theta(sqrt(nu/u)) in every regime; a = u (cheap binomial regime at the top valuations) gains
only a constant < 2.5 which shrinks as u grows. The optimum corresponds to max Lambda = 5.84 M/u
(G14b's LP: 5.9), against the paper's actual loss 18 M/u: a factor sqrt 3 = 1.76 in the exponent
(paper, nu = 6: pi sqrt(2*6*18/3) = 26.7 M/sqrt u; optimum 6.2 sqrt 6 = 15.2 M/sqrt u; natural-log units).

## [intermediate save 3] 2-adic Smith profile of the Cauchy crossing block (smith_cauchy.py, smith_out.txt)   [N]
Full packet, L = 2^a (a = 4,5,6), t = 2^i, zeta = 1, crossing block (rows, columns with equal bit i).
Observed EXACTLY in all 15 cases: the elementary divisors are 2^{s_k}, each with multiplicity 2^i
(= number of tiles in the block),   s_k = (a-i+1) + 2k + 2 v2(k!),  k = 0..D/2-1  (so s_k ~ a-i+1+4k).
App. D's modes j^b have true weight (a-i+1) + b: one per unit weight and tile.  So
  #{invariants <= w} = (1/4 + o(1)) #{modes of weight <= w}:
the linear growth of the 2-adic rank profile is intrinsic to the kernel; the monomial basis j^b is
4x redundant (a Mahler-type basis on both sides would give n -> n/4, i.e. a factor 2 in the exponent).

## 2. Parametric k-level model and its optimum (Task 2)

Model VT_k ("App. D with free parameters, recursive").  Free per depth d: r_d (groups), u_d <= a_d,
taper q^(d)(i) >= 0 (i < u_d, sum_t q_t > M_d(2L_d - 1)), packing with fill phi_d (= M_d/(#packets*L_d)),
depth k, leaves by brute force.  Fixed: grid S = {-L..L-1}, Cauchy kernel (D.7), recovery (D.9),
tiles/modes (D.12), enumeration of ALL tuples (D.14), children = groups.
Named quantities at a node t of valuation i (D_t = L/2^i):
   p_t >= D_t * Lambda_i,   Lambda_i = sum_j (j/2^j)(q(i) - q(i+j))        (loss of (D.9), eq. (D.10))
   n_t  = #packets * 2^{i+1} = 2M/(phi D_t)                                 (modes per unit weight)
   ln #tuples_t = f(p_t, n_t) := min_theta [theta p + n sum_j -ln(1 - e^{-theta j})]   (leading order)
       = pi sqrt(2 n p/3) if p >> n ("partition regime");  ~ p ln(e n/p) if p << n ("binomial regime").
Constraints:
 (K1) phi*u <= log2(M) + 1                 [P]  (a packet has 2^a slots; phi*2^a*#packets = M; u <= a)
      with the paper's hashing: r^{2a} 2^a <= M, i.e. u <= a <= log2 M/(1 + 2 log2 r)      [D]
 (K2) max_t min(p_t, sqrt(n_t p_t)) >= 0.43 M sqrt(nu/u), nu = 2/phi  (Prop. B' below)      [P]
 (K3) a_d >= a_{d-1} + u_d + O(1) (so that 2^{a_d-u_d} absorbs the inherited precision)    [D]
      => sum_d u_d <= a_k <= U := log2(M_k)/(1 + 2 log2 r).
 (K4) some child has >= 2M_d/r_d wires (2M ports, r groups)                                 [P]
Cost exponent:  E = Theta( sum_{d<k} M_d/sqrt(phi_d u_d) + M_k ),  M_d = M prod_{j<d} (c/r_j).
Optimisation.
 * Depth -> infinity, r constant: min sum_d M_d u_d^{-1/2} s.t. sum u_d <= U gives u_d ∝ M_d^{2/3},
   E = (M/sqrt U)(1 - rho^{2/3})^{-3/2}, rho = M_{d+1}/M_d.  rho = 1/8: 1.54 M/sqrt U (paper's
   u_d = U 2^{-d-1}: 2.19 M/sqrt U).  Depth changes the CONSTANT only: the root term M/sqrt(u_0),
   u_0 <= U <= log2 M, is already of the final order.
 * Fixed depth k (leaves 2^{M_k}): need prod r_d >= sqrt(log M) while (K1)+(K3) force
   log r_d <= log M/(2 a_0): r = (log M)^{1/2k}, E_k = Theta(M sqrt(loglog M/(k log M))), 1 <= k <= loglog M.
   So recursion depth buys exactly the factor sqrt(loglog) and nothing else.  (If (K3) is dropped
   and r_d may grow with d, depth log* M suffices; the limit is the same.)
 * Other taper: LP optimum 6.2 M sqrt(nu/u) vs paper 10.9 M sqrt(nu/u): constant 0.57 (= 1/sqrt 3, as G14b).
 * a = u (no slack, binomial regime at the top valuations): constant < 2.5, vanishing gain (sec. save 2).
 * Finer tiles 2^{i+1+c}: 2^c times more modes of weight step 1+c: worse by sqrt(2^c/(1+c)).
 * Sparse packets (phi < 1): effective u is phi*u: worse.   Odd prime P: needs nodes with leading
   digit 1 only; same scaling with u <= log_P M: worse constant.                              [E, sketch]
 * Better mode basis (Smith/Mahler, sec. 3b): n -> n/4, constant 0.5.                         [N]
 Net: no power of log is gained by depth, taper, balance or tiling.  Achievable gain: the constant
 (roughly 0.57 * 0.5 * 0.71 ~ 1/5 in the exponent, estimate, unverified as an algorithm).
BINDING CONSTRAINT:   n_t * p_t >= c * nu * M^2/u  at some retained node (Prop. B'),  with
   phi * u <= log2 M + 1,  so  ln #tuples >= c' M/sqrt(log2 M).
 The binomial regime (p << n) would give M log(u)/u = n loglog n/log n, but it is unreachable:
 sum_i Lambda_i >= 3M/8 needs ~u levels each carrying Lambda ~ M/u, and only the top
 (1/2)log2 u levels (D_t <= sqrt u) are in that regime.

## 3. Lower bounds (Task 3)

### 3a. Prop. B' (per-level inequality, both regimes; extends G14b Thm B to a = u)   [P, design-level]
For every profile q (q(i) = 0 for i >= u, 1 <= u <= a, L >= 2, sum_t q_t > M(2L-1)) put X :=
max_{i<u} min(p_i, sqrt(n_i p_i)) with p_i = D_i Lambda_i, n_i = nu M/D_i.  Then
      X + u X^2/(nu M) > 3M/8,   hence   X > M min(3/16, sqrt(3 nu/(16 u))).
Proof. min(D Lambda, sqrt(nu M Lambda)) <= X gives Lambda_i <= lam_i := X/D_i + X^2/(nu M).
sigma_i q_i = Lambda_i + sum_{j>=1} (j/2^j) q_{i+j}, sigma_i = sum_{j<=a-i} j/2^j in [1/2, 2).
Downward induction: q_i <= 2 S_i, S_i := sum_{k>=i} lam_k (S decreasing, so the tail sum is
<= 2 S_{i+1} sigma_i, and lam_i/sigma_i <= 2 lam_i).  Then M(2L-1) < sum_i (L/2^i) q_i <= 4L sum_k lam_k,
and sum_k X/D_k <= X (D_k >= 2^{u-k}).  So X + uX^2/(nu M) > M(2-1/L)/4 >= 3M/8.  QED
Counting: #tuples >= 2^{min(p, sqrt(n p))/2 - log2 n - 1}  (p > n: 0/1 exponents on modes of weight
<= floor(sqrt(p/n)); p <= n: subsets of weight-1 modes).  With nu >= 2 (all formal modes of all
packets are enumerated; #packets >= M/L):  the root call enumerates >= 2^{0.3 M/sqrt(u) - 1} tuples
(u >= 11), for EVERY taper, slack a-u, fill, depth.  With (K1): >= 2^{c M/sqrt(log2 M)}.
Numerics (lp_taper.py): true optimum ln #tuples = (6.0-6.2) M sqrt(nu/u) for a >= u+4.
Caveats: p_i >= D_i Lambda_i is the loss of formula (D.9) (G17a: black-box precision is of the same
order but the formula is not forced); "enumerates all tuples" is the algorithm as written.

### 3b. 2-adic rank lemma + exact Smith profile of the Cauchy crossing block (target (b))
Lemma T [P] (standard, DVR O, uniformiser 2).  For a matrix A over O let R_w(A) := #{elementary
divisors of A not divisible by 2^w}.  (i) R_w(A) = least N with A == sum of N rank-one matrices
mod 2^w.  (ii) R_w(submatrix) <= R_w(A).  (iii) Smith valuations of A (x) B are the sums s_j(A)+s_k(B);
so R_p(A_1 (x) ... (x) A_m) = #{(k_1..k_m): sum_j s_{k_j}(A_j) < p}  ("budgeted" count).
Consequences.
 (1) Any mode system E == sum_Z c_Z (x) L_Z with v2(c_Z) >= w_Z for a crossing block E has
     #{Z: w_Z < w} >= R_w(E).
 (2) Model KF ("kernel-factored node"): at node (t, zeta) the deformed wire kernel
     K_T(Y,X) = prod_x l_{Y_x}((W(T)X)_x) is written mod 2^{p_t} as sum_{k<=N} prod_{groups g} tau_g^(k).
     For every bipartition of the groups, N >= R_{p_t}(flattening) = #{(k_P)_{P cut}: sum s_{k_P} < p_t},
     s_.(P) = Smith valuations of the flattening of packet P.  [P]
     One-crossing sub-block [P]: the flattening of a packet contains (rows (Y_g = e_x, X_g = 0),
     columns (Y_g' = e_{y-t}, X_g' = e_y)) the matrix [E_xy/(1 - E_xy)] * unit diagonal = the crossing
     Cauchy block at T' = T - C(T)/C'(0); so its invariants include {0} and those of the crossing
     block.  Alone this only gives N >= (L/2+1)^{#cut packets} ~ exp((M/2L) ln L): weak.
Numerics [N] (smith_cauchy.py, flat_packet.py):
 * crossing block, full packet, a <= 6, all valuations, several t, zeta: elementary divisors
   2^{s_k}, s_k = (a-i+1) + v2(zeta) + 2k + 2 v2(k!), k < D/2, multiplicity 2^i per block.  EXACT
   match in 23 cases.  [formula conjectured in general; a Cauchy-determinant proof looks feasible]
   => R_w(E) = (#tiles)(w-(a-i))/4 (1+o(1)): linear in w.  The paper's modes are 4x redundant, and
   NO linear mode basis for this kernel has fewer than n/4 modes per unit weight.
 * flattening of ONE packet tensor across the two groups (L = 4, 8; 16x16 and 256x256, exact
   2-adic Smith form): R_p(flattening) EQUALS the number of monomials of weight < p in the
   Smith-optimal modes for all small p (L=8, i=0: 1,1,3,5,8,14,33 for p = 2,4,6,8,10,12,16; the
   paper's modes give 1,1,5,9,20,39,141), and falls below it only when p approaches saturation
   (finite 2^L, Boolean X^2 = X, <= L/2 rows).
Conjecture C1 [C, supported only by L <= 8]: in the unsaturated range the flattening invariants
of a packet are the weights of the Smith-mode monomials; with Lemma T(iii) this would give
N_KF >= exp(c sqrt(n_t p_t)) = exp(c M/sqrt u): the mode enumeration is optimal up to the constant
among ALL term-by-term group-product expansions at a node (link (ii) of G17a, oblivious form).

### 3c. What is NOT proved
 * No lower bound for general kernels W (other nice matrices N_t, higher degree, non-dyadic nodes),
   for other recovery rules with smaller precision, or for non-oblivious re-summation of the beta-sum.
   Still open: DF-model direction rank (proved 2^{2M/4^u}, G14b Thm A; target 2^{cM/sqrt u}).
 * The 2-adic analogue of our Lean rank inequality (sum_k rank E_k >= sum_x max(0, N-d-2mu(x))) is NOT
   proved.  Natural statement [C]: for deg W <= d, W(0) = I, nodes with discs on which W is integral
   and nice mod 2, the rank profiles R_w of the crossing blocks satisfy an inequality forcing
   sum_w R_w ~ linear growth.  Only the scalar case (one fully cut wire: a degree-<= d polynomial with
   value 1 at 0 cannot be == 0 mod 2^w on all discs beyond the Lagrange/Hermite loss) is immediate.

## 4. Verdict (Task 4)

| family | best upper bound (exponent) | best lower bound PROVED | gap |
|---|---|---|---|
| VT_k: App. D design, free taper, a, u, r, fill, depth; all tuples enumerated | O(s/sqrt(log s)) (paper) | Omega(s/sqrt(log s)) [P, design-level: Prop. B' + (K1)] | constant (~5, incl. taper 1.76, Smith modes 2, a-increments 1.4) |
| VT_k, fixed depth k, hashing packing, (K3) | Theta(s sqrt(loglog s/(k log s))) [E] | same order [P given (K1),(K3),(K4) as design axioms] | constant |
| KF: oblivious group-product expansion of the Cauchy kernel at one node mod 2^{p_t} | exp(O(sqrt(n_t p_t))) = exp(O(M/sqrt u)) | budgeted Smith count (Lemma T) [P]; evaluated only to exp((M/2L) ln L) [P]; exp(c M/sqrt u) under C1 [C,N] | exponential until C1 is proved |
| black-box recovery precision on Lemma D.5's class | Theta(M D_t/u) | Theta(M D_t/u) (G17a) | constant ~12 |
| DF: any direction-factored 2-adic identity, one level | 2^{O(M/sqrt u)} | 2^{2M/4^u} (G14b Thm A) | exponential in u |
| general recursive coordinate-nice interpolation (any kernel/nodes, re-summation) | s/sqrt(log s) | none beyond DF | open |

Answer to the question.  Within every family analysed the limit is n/sqrt(log n) (bounded density):
depth, taper, slack a-u, tile size, fill, prime, and regime balance change constants only (and
sqrt(loglog) at bounded depth).  Nothing here supports n/(log n)^{2/3} or n/log n.  The sqrt is the
partition count exp(sqrt(n p)) of a LINEAR 2-adic rank profile against a budget, with
n_t p_t >= c M^2/u and u <= log2 M (packets are full antipodal subcubes: 2^a <= M/#packets).
Most promising place to gain / to close:
 (a) barrier side: prove C1 and the Smith formula for the Cauchy kernel (finite, concrete: Cauchy
     determinants + symmetric-power structure of the packet tensor), then the general-kernel 2-adic
     rank-profile inequality; this would upgrade "design-level" to "every oblivious enumeration".
 (b) algorithm side: the only lever of non-constant size is a kernel whose crossing blocks have a
     SUBLINEAR 2-adic rank profile R_w = o(w * #tiles) (then the count is no longer exp(sqrt(np))),
     or more than log2 M directions; both are excluded for the Cauchy kernel [N] and the second for
     any packet structure [P].  Constant-factor gain available now: Smith/Mahler mode basis.

## 5. Adversarial self-review
 * All "lower bounds" on the VT_k rows are statements about a DESIGN (formula (D.9), full
   enumeration); they are tightness results, not barriers.  G17a already showed (D.9)'s loss is not
   forced (10 bits vs 40 in a toy), although the black-box order M D_t/u is.
 * Prop. B' counting uses formal modes of all packets (nu >= 2); skipping empty tiles or dead modes
   is not modelled.  Constants (0.43, 0.13) are crude; LP gives 6.2 (natural log, per sqrt(nu/u)).
 * (K3) is derived from the paper's uniform-per-depth parameters; per-call adaptive parameters are
   not excluded (does not affect the limit, only the fixed-depth formula).
 * Smith formula and C1 are numerical (a <= 6; L <= 8), zeta integer, full packets only; G17a saw
   much smaller precision for sparse packets.  The flattening/monomial equality fails at saturation.
 * The saddle-point count is leading order only.  The odd-prime and tile remarks are sketches.
 * The estimate "constant ~1/5 gain" multiplies three unverified constants; not an algorithm.
 * App. D's correctness was not checked; the TeX source was unavailable (pdftotext read).
Files: notes.md, lp_taper.py, lp_taper_out.txt, smith_cauchy.py, smith_out.txt, flat_packet.py, flat_out.txt
