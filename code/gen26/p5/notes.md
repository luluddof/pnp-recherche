# G26-P5 notes: Freedman version of R2, the floor of the grouping design, structural options

Tags: [P] proof written here (unreviewed), [N] own numerics, [R] inherited (paper / p2 / p9 / b7 as reviewed by G26-P2),
[H] heuristic. Model = reviewer's (G26-P2) mylib/mydp/mythr, copied here and extended (mylib.py additions at the end:
zeta_F, zeta_floor, zeta_lin2, zeta_randfloor; mydp5.py = mydp.py + choice of zeta + zeta cache + 'rowzero' base).

## log
- 02:18 start; read b7, p2, p9 notes, p2 scripts, paper Lemma 8 (groups may have up to 2g/K nodes in case (i), exactly
  g/K in case (ii); L = ceil((U + 2 zeta g P)/K^2) = max over pairs Y).
- 02:27 zeta_F implemented; DP runs launched (rowR3 x {F, R2, floor}, rowzero x floor).
- 02:28 test_r2f.out: variance bound, increment bound, linear bound all hold on 6 network kinds (max ratios .57/.42/.23);
  stratified excess far below R2F bound at small K; rank-2 network: excess*t ~ 0.75-1.0 ln K (random stratified grouping
  cannot beat ~ ln K/t).

## 1. Lemma R2F (Freedman/Bernstein version of R2)  [P]
Setting as R2: network on g = tK nodes, symmetric c_ab in [0, P/g], c_aa = 0, c_a = sum_b c_ab <= P, U = sum_a c_a, N = gP.
Sort nodes by c_a (ascending), strata S_1..S_t of K consecutive nodes, independent uniform bijections pi_k: S_k -> [K],
group gamma = {pi_k^{-1}(gamma)}_k (exactly t nodes). Claim: with probability >= 1 - d1 - d2 - d3, for all gamma != gamma',
   K^2 Y_{gamma gamma'} <= U K/(K-1) + 2N/t + (N/t^2) s',                                   (R2F)
   s' = b' L2/3 + sqrt((b' L2/3)^2 + 2 V' L2),    b' = 4 sqrt((t-1)L/2) + 3,
   V' = sum_{k=1}^t ( s2 (m sqrt(k-1) + 1/2) + 1/2 )^2,  s2 = sqrt(2(1 + 1/(K-1))),  m = min( sqrt(L/2), 1/2 + sqrt(ell/2) ),
   L = ln(4gK/d1),  ell = ln(2Kt/d3),  L2 = ln(K(K-1)/(2 d2)).
So zeta_F = 2/t + s'/t^2 + 1/(K-1) ~ (2 + sqrt(2 m^2 * 2 L2))/t ~ (2 + sqrt(L L2))/t, against R2's (4 + sqrt(8 L L2))/t:
the martingale term shrinks by 2 sqrt2 (more with the McDiarmid m), the linear term by 2.

Proof. Fix gamma != gamma'; a_k = pi_k^{-1}gamma, b_k = pi_k^{-1}gamma'. The pairs (a_k, b_k), k = 1..t, are independent,
(a_k, b_k) uniform on ordered pairs of distinct nodes of S_k. G_k = sigma((a_i, b_i), i <= k), M_k = E[Y | G_k].
(a) Mean: E Y <= U/(K(K-1)) (R2 (a), [R]).
(b) Doob increments (re-derived; = b7/P2's formula): D_k = M_k - M_{k-1} = [h(a_k) - avg_{S_k} h] + [h'(b_k) - avg h']
    + [c(a_k,b_k) - c(S_k,S_k)/(K(K-1))], h(a) = c(a, B_{<k}) + c(a, S_{>k})/K, h' symmetric with A_{<k}.
    Write h(a) = c_a/K + f(a), f(a) = Z_k(a) - c(a,S_k)/K, Z_k(a) = sum_{j<k} [c(a,b_j) - c(a,S_j)/K]
    (because c(a,S_{<k}) + c(a,S_{>k}) = c_a - c(a,S_k)). Hence D_k = Lambda_k + Delta_k with
      Lambda_k = (c_{a_k} - avg_{S_k} c)/K + (c_{b_k} - avg_{S_k} c)/K,
      Delta_k  = u(a_k) + u'(b_k) + w(a_k,b_k),  u = f - avg_{S_k} f,  u' likewise,  w(a,b) = c(a,b) - E c(a_k,b_k).
    Each of the three pieces of Delta_k and Lambda_k has conditional mean 0 given G_{k-1} (a_k, b_k marginally uniform on
    S_k), so (Delta_k) is a martingale difference sequence. NOTE the difference with R2: the term -c(a,S_k)/K (range < P/g per
    stratum, total P/K) now sits in the martingale, where it costs O(1/t^{3/2}) instead of 2P/K deterministically.
(c) Linear part, deterministic: with M_k = max_{S_k} c, m_k = min, A_k = avg_{S_k} c and M_k <= m_{k+1} <= A_{k+1} (sorted,
    consecutive strata), sum_k (M_k - A_k) <= M_t - A_1 <= P (telescoping). So sum_k Lambda_k <= 2P/K, i.e. <= 2N/t after
    multiplying by K^2 (K P = N/t).
(d) Good events (predictable). E1: |Z_k(a)| <= z_k := (P/g) sqrt((k-1)L/2) for all k, a in S_k, gamma' (and the same for
    Z'_k(b) = sum_{i<k}[c(a_i,b) - c(S_i,b)/K], all gamma): Hoeffding (k-1 independent terms, each in an interval of length
    P/g), 2gK events x two tails: Pr(not E1) <= d1. E3: rms_{a in S_k} Z_k(a) <= y_k := (P/g) sqrt(k-1) (1/2 + sqrt(ell/2))
    for all k, gamma' (and Z' for all gamma), 2Kt events: the vector Z_k = (Z_k(a))_{a in S_k} is a function of the
    independent b_1..b_{k-1}; changing one b_j moves every coordinate by <= P/g, hence the rms norm by <= P/g; McDiarmid
    gives rms <= E rms + (P/g) sqrt((k-1) ell/2), and E rms <= sqrt(E mean_a Z_k(a)^2) = sqrt(mean_a sum_j Var c(a,b_j))
    <= (P/g) sqrt(k-1)/2 (Popoviciu). Pr(not E3) <= d3. Both events at step k are G_{k-1}-measurable.
(e) Conditional variance on E1 and E3. Let s_u^2 = mean_{S_k} u^2, s_u'^2, s_w^2 = E w^2 <= (P/2g)^2 (c in [0,P/g]).
    For (a,b) uniform on distinct pairs, E[u(a) u'(b)] = -(1/(K(K-1))) sum_a u(a)u'(a) (since sum u = 0), so
    |E u u'| <= s_u s_u'/(K-1) and E(u+u')^2 <= (1 + 1/(K-1))(s_u^2 + s_u'^2); Cauchy-Schwarz on the w cross term:
      Var(Delta_k | G_{k-1}) <= ( sqrt((1+1/(K-1))(s_u^2 + s_u'^2)) + s_w )^2.
    s_u <= sd(Z_k) + sd(c(.,S_k)/K) <= min(z_k, y_k) + P/(2g) (std <= rms <= max; range of c(a,S_k)/K < P/g).
    Hence Var(Delta_k | G_{k-1}) <= v_k := (P/g)^2 ( s2 (m sqrt(k-1) + 1/2) + 1/2 )^2, sum_k v_k = (P/g)^2 V'.
    Upper bound of the increment: Delta_k <= (max f - min f) + (max f' - min f') + P/g <= 2(2 z_k + P/g) + P/g
    <= (P/g) b'.
(f) Freedman's inequality for the martingale stopped at the first k where E1 or E3 fails at step k (predictable stopping,
    so increments, variance bound and b' remain valid): Pr(sum Delta_k >= s) <= exp(-s^2/(2(V + b s/3))), V = (P/g)^2 V',
    b = (P/g) b'. With exp(-L2) = 2 d2/(K(K-1)) and a union over unordered pairs (Y symmetric) we get, on E1 and E3, all
    pairs <= (P/g) s' except with probability d2. Multiply by K^2: K^2 (P/g) = N/t^2.  Sum: (R2F). QED
Exhaustive search over groupings (as in p2/R2) finds one meeting (R2F). Budget used in the numbers: d1 = d2 = d3 = 0.33.
Checks [N] (test_r2f.out): exact conditional variance <= v_k bound (max ratio 0.57), Delta <= b (0.42), sum Lambda <= 2P/K
(0.23) on unif/hetero/rank2/star/blocks/Bernoulli networks along random histories; observed stratified excess always far
below the bound.

### Effect of R2F (DP of P2's model, rowR3 base with P2's own R3 constants, Lemma U, bins 1.02, lt step 1)  [N]
dp_rowR3_F.out vs dp_rowR3_R2.out (same code, same run; R2 reproduces P2's dp_rowR3.out exactly):
| | [1000,2000] | [2000,4000] | [4000,8000] | [8000,16000] | mu(x<=16000) |
|---|---|---|---|---|---|
| R2  a=2 | 3.996 | 4.022 | 4.038 | 4.050 | 0.5030 |
| R2F a=2 | 3.871 | 3.899 | 3.901 | **3.921** | 0.5361 |
| R2F a=4 | 5.474 | 5.515 | 5.517 | 5.545 | |
Limit (slow drift, as P2's 4.05 -> ~4.06): R2F c1 ~ 3.93 (deterministic = Las Vegas by R1).
(b7's heuristic "Freedman ~3.84" assumed the McDiarmid-free optimal variance; the proved bound keeps one log, see §2.)

## 2. The remaining overheads and the floor of the design

### 2.1 Inventory (per level of Lemma U; model of P2/p9)
e' = max_th min(e_a(cos th + 1/h + z_a), e_b(sin th + 1/g + z_b)) + 2^{-min(x_a,x_b)},  x' = x_a + lt_a + x_b + lt_b.
Balanced form: Phi' = 2 Phi/(1 + sqrt2 (z + 1/h))^2 (entries aside), x' = 2(x + lt), so per level
    ln(mu'/mu) ~ -2 sqrt2 z(t,K) - lt/x - (entry, 1/h terms),
and mu_inf = Phi_1 / A_inf * Gamma^{-2}, A_inf = x1 + sum_l 2^{-l} lt_l (effective base size), Gamma = prod (1 + sqrt2 z_l + ...).
The optimal t balances the two: d/dlt [2 sqrt2 z + lt/x] = 0. Writing z = C(x)/t:
    2^{lt} ~ 2 sqrt2 ln2 C(x) x,   loss per level ~ (lt + 1/ln2)/x.
So the overheads are (i) the grouping excess z (the only reason t > 1), (ii) the growth lt ~ log2(C(x) x) it forces, (iii) the
straight edges 1/h, (iv) entries 2^{-x}, (v) the base (row sweep: R3 terms zs, 2/g, discretisation s(K)/K^2), (vi) top a = 2 and
the sawtooth/granularity (removed up to ~1% by unbalanced levels). (iii), (iv) are exponentially/inversely small at the optimum;
(i)+(ii) dominate: C(x) = growth exponent of the excess constant in log K:
   R2:  C ~ 4 + sqrt(8 L L2) ~ 4 ln2 sqrt(2 x(2x + lt))       -> lt ~ 2 log2 x + 3.4
   R2F: C ~ 2 + sqrt(L L2)  ~ 2 ln2 sqrt(x(2x + lt)) (x 0.85)  -> lt ~ 2 log2 x + 1.5
   random stratified grouping, ANY analysis: C >= ~ ln K = x ln2  (rank-2 example below)  -> lt ~ 2 log2 x + 0.4
   ANY grouping, network-independent bound (counting floor): C >= kappa y*, y* ~ sqrt(2 ln(K/(2t lnK)))  -> lt ~ 1.5 log2 x
   C = O(1) (no log K)                                            -> lt ~ log2 x + O(1)
Since x_l ~ 2^l x1, sum_l 2^{-l} lt_l ~ alpha(2 log2 x1 + 2) + 2 beta for lt = alpha log2 x + beta: A_inf ~ x1 + 2 alpha log2 x1 + O(1).
The power alpha of log x in the growth is what the concentration lemma controls; the floor is about alpha.

### 2.2 Floor statements
(F1) Random stratified grouping cannot beat C ~ ln K  [P-sketch + N]. Network c_ab = (P/2g)(1 + e_a e_b), e = +-1 balanced in every
  stratum (all node sums equal, so the sorted strata are arbitrary). Then Y_{gg'} - EY = (P/2g) X_g X_g', X_g = sum_k e(a_k) a sum of
  t independent signs; max over pairs of X_g X_g' ~ 2 t ln K, so K^2 Y - U ~ N ln K / t. test_r2f.out: excess*t = 0.75-1.0 x ln K
  (K = 8..64). So R2F (C ~ 2 ln2 sqrt(x(2x+lt)) x 0.85 ~ 1.7-2.3 ln K) is within a factor ~2 of the best possible analysis of the
  random stratified grouping; the leading growth 2 log2 x of lt cannot be improved by sharper martingale inequalities.
(F2) Counting floor for every grouping lemma of the form "for all admissible networks, some grouping has K^2 Y <= U K/(K-1) + zeta N"
  [P, first moment; normal approximation of the binomial tail = Cramer moderate deviations, asymptotic]. Random network
  c_ab = (P/g) Bern(nu) i.i.d.; for a fixed equal-size grouping the K(K-1)/2 pair sums are independent Bin(t^2, nu)(P/g);
  Pr(all <= nu t^2 + y kappa t) <= exp(-q_y K(K-1)/2), kappa = sqrt(nu(1-nu)); at most K^g groupings; so if
  q_y > 2 t ln K/(K-1) no grouping is good, i.e. zeta >= kappa y*/t with Phibar(y*) = 2 t ln K/(K-1) (U = nu N(1-1/g)).
  Used with nu = 1/sqrt2 (the crossing point of the arc), kappa = 0.455 (mylib.zeta_floor).
(F3) Linear part: any equal-size grouping of the network "K/2 isolated nodes, all other pairs P/g" has a pair of all-heavy groups:
  zeta >= 1/t - 2/K. So zeta = Omega(1/t) for every grouping design: lt >= log2 x + O(1) in any design that needs groups.

## 3. Structural change S2: edge-weighted cliques with per-rectangle padding  [P-sketch]
Observation: the ONLY place where the max over group pairs enters is the uniform width L = max_{g != g'} Y_{gg'} of the child
(K, L)-clique (paper Lemma 8: L = ceil((U + 2 zeta g P)/K^2)). Every later step is local:
 * Lemma 1/8: the group pair (g, g') can carry its own string length Y_{gg'} (Definition 3 with per-edge widths lambda_e; the
   factor of the pair factors through Y_{gg'} bits by observation (a)/(b) after Def. 1 — nothing needs equal lengths).
 * Prop. 3 (paper lines 994-1052): every diagonal belongs to exactly one rectangle q; the identity of Lemma 2(a) is applied
   rectangle by rectangle with its own iota(q) in {0..t}; the row/column networks are sums over rectangles (35). With unequal
   diagonal widths lambda_1, lambda_2 in q, pad the shorter string with zeros to lambda_q := max(lambda_1, lambda_2) ([x = v] is
   unchanged) and apply Lemma 2 with lambda_q, t_q = ceil(lambda_q/3) (or App. B pieces of size lambda_q). Straight edges keep their
   widths. Lemma U's Cauchy-Schwarz is per rectangle (p9), so the level bound holds with N replaced by sum over rectangles of
   2 lambda_q (+ straight), i.e. by the PADDED total width.
 * Padding cost, deterministic for ANY placement on the grid: sum_q (max(l1,l2) - (l1+l2)/2) = sum_q |l1 - l2|/2
   <= (1/2) sum_e |lambda_e - c| for any constant c (each diagonal in one rectangle).
 * Hence the grouping only has to control the AVERAGE absolute deviation of the pair sums, not their max:
   Lemma S2 [P]. For the stratified random grouping of R2F, every pair has E Y = EY (same for all pairs) and
     sd(Y_{gg'}) <= (P/K) [ 1/sqrt2 + sqrt(sum_k (s2 (sqrt(k-1)/2 + 1/2) + 1/2)^2)/t ]  (~ (P/K)(0.707 + 0.5)).
   Proof: Var Y = sum_k E D_k^2 (orthogonal Doob increments), Minkowski over k with D_k = Lambda_k + Delta_k:
   E Lambda_k^2 = (2/K^2) sd_{S_k}(c)^2 (1 - 1/(K-1)) <= r_k^2/(2K^2), sum r_k^2 <= (sum r_k)^2 <= P^2;
   E Delta_k^2 <= E v_k with E rms(Z_k)^2 <= (k-1)(P/g)^2/4 (Popoviciu) in place of the high-probability bounds of R2F.
   So E sum_{g<g'} |Y - EY| <= (K(K-1)/2) max sd, and some grouping (found by exhaustive search, the quantity is computable)
   has padded total <= sum_e Y_e + (K(K-1)/4) max sd. In the units of the level formula (N = K^2 x width):
     zeta_S2 = (1/2) K^2 sd/N + 1/(K-1) ~ 0.60/t + 1/(K-1)     (mylib.zeta_S2; no log K at all).
   (The placement for the straight-edge totals is then chosen by averaging, as in the uniform case; the padding bound holds for
   every placement, so both are met together.)
 * Top level: the top clique can likewise be edge-weighted (cluster-pair counts as widths); in the threshold runs I keep
   P2's top (Lemma 8' Bernstein deviation d unchanged) and only replace the top grouping zeta by zeta_S2 (conservative).
Status: the grouping lemma S2 is proved; the per-edge-width versions of Def. 3 / Lemma 1 / Lemma 8 / Prop. 3 / Lemma U are a
sketch (each step is rectangle- or pair-local in the paper's text, checked by reading lines 994-1052; not re-proved line by line).
Caveat: rho_bar(lambda_q) needs lambda_q -> inf for all rectangles (asymptotic model; finite-lambda rho_bar is per rectangle).
Check [N] (test_s2.out): empirical sd(Y) <= 0.44 x bound on bern/unif/rank2/jump/hetero/star (K = 16..32, t = 16..64);
empirical S2 constant (K^2 mean|Y-EY|/(2N)) * t = 0.11-0.23 vs the proved 0.60.

## log (cont.)
- 02:35 DP results: rowR3 x {floor 3.361, randfloor 3.839, S2 3.352}; rowzero x floor 3.159 (all a=2, [8000,16000]).
- 02:38 thresholds (mythr5 = P2's mythr with selectable top grouping): calibration rowR3+R2 top Hoeffding = 2237/4385/1.494e4 |
  73.47/117.2/415.1 (= P2 exactly).
- 02:39 chain5.out (balanced chain, greedy lt): ln-loss split growth / excess / 1/h / entries, e.g. R2F rowR3:
  -0.742 / -0.071 / 0 / 0; lt ~ 11.5, 14, 16.25, 18.25, ... at x = 25, 73, 174, 380 (lt ~ 2 log2 x + 2).
- 02:40 R3 base: the Gamma term dominates the base overhead (Gam ~ 10 cA). Observation: P2's Azuma uses |dX| <= c; the increment
  lies in a predictable interval of length c (= 2 (max - min) of gamma(v,S_j) over v / w), so Hoeffding-Azuma gives sqrt(S Lam/2)
  instead of sqrt(2 S Lam): **Gam -> Gam/2** ("R3'", [P], trivial). DP mode 'rowhalf'.
- 02:40 DP rowR3 x zero (no grouping excess at all, lt >= 1 forced): 2.309 = what the non-grouping overheads alone allow.
- 02:41 DP rowzero x S2: 3.005; rowR3 x lin2 (2/t only): 3.480.

## 1b. Lemma R3'' (base ordering lemma, two free sharpenings of P2's constants)  [P]
In R3 (P2's own constants) the Gamma term dominates the base overhead (Gam ~ 10 cA at g = 2^10..2^16).
(i) dX_j = X_{j+1} - X_j = -2(gamma(v,S_j) - avg_{u in S_j^c} gamma(u,S_j))/w_{j+1}, v uniform on S_j^c. Given S_j (and the
    predictable good event G_j <= bound) dX_j lies in an interval of length 2(max_v - min_v gamma(v,S_j))/w_{j+1} <= 4G_j/w_{j+1}
    = c_j. Hoeffding's lemma conditionally => Pr(X_n - X_0 >= s) <= exp(-2 s^2/sum c_j^2): deviation sqrt(S Lam/2), half of
    P2's sqrt(2 S Lam) (which used |dX| <= c, i.e. Azuma's weaker form). (A Freedman version gives the same leading term here:
    one random node per step, so no sqrt2 from independence as in R2F.)
(ii) gamma(v,S_i) for v not in S_i is the sum over a uniform i-subset of [g]-{v} of a zero-sum list with range R: Serfling's
    inequality adds the factor (1 - (i-1)/(g-1)) to the variance proxy in G_i; same for A_j (factor 1 - (j-1)/g, max at j = g/2
    -> cA x (1/2 + 1/g)). Same events and budget as P2.
Net: Gam x 0.3875, cA x 0.5 (mylib.r3_terms_v2; DP mode 'rowser'; 'rowhalf' = (i) only).

## 4. Results (DP of P2's model; c1 = limsup over log m/a in [8000,16000], a = 2; det = LV by R1)
| variant | status | c1 | r=4.26 vs 2^n / 1.6423^n / 1.30704^n | r=1 vs same |
|---|---|---|---|---|
| R2 + R3 (P2, reproduced) | reviewed | 4.050 | 2237 / 4385 / 1.494e4 | 73.5 / 117.2 / 415 |
| R2F + R3 (top Hoeffding / top R2F) | [P] | 3.921 | 2075 / 4062 / 1.408e4  (2057 / 4042 / 1.407e4) | 73.5/117.2/395 (65.4/98.0/377) |
| R2F + R3' (Gam/2) | [P] | 3.851 | | |
| **R2F + R3''** (top R2F) | [P] | **3.824** | **1943 / 3821 / 1.341e4** (top Hoeff: 1963/3842/1.343e4) | 59.3 / 93.0 / 355 |
| floor A: random stratified grouping, any analysis (zeta = ln K/t), R3'' / ideal base | [P-sketch] | 3.726 / 3.519 | | |
| floor B: ANY grouping, network-independent zeta (counting), R3 or R3'' / ideal base | [P-asympt.] | 3.361 / 3.159 | | |
| S2 (edge-weighted cliques) + R3 / R3' / R3'' / ideal base | [P-sketch] | 3.352 (3.342 at lt step 0.5) / 3.269 / **3.226** / 3.005 | R3'': **1419 / 2747 / 9645**; ideal base 1210 / 2396 / 8231 | R3'': <=50 (bisection floor) / 80.7 / 265 |
| no grouping excess at all (lt >= 1 forced), R3 | [H, unattainable: F3] | 2.309 | | |
| poly space (D&C base): R2 / R2F / S2 | [P] / [P] / [P-sketch] | 6.162 / 5.875 / 4.611 | | |
(Floor B with R3/R3'' coincide: there the DP prefers clique-sweep leaves at small x, where the counting floor is 0.)
Balanced-chain split (chain5.out): growth by lt is 85-90% of the log-loss of mu in every variant; excess 10-15%; 1/h and
entries < 0.1%. Effective growth lt at x ~ 25..1600: R2/R2F ~ 2 log2 x + (3.4 / 1.5); floor B ~ 1.5 log2 x; S2 ~ log2 x + 0.2.

## 5. Verdict on the floor
* The design (uniform-width child cliques, max over group pairs) with ANY grouping lemma of the network-independent form
  "K^2 Y <= U K/(K-1) + zeta(t,K) N" cannot go below c1 ~ 3.16 (ideal base) / 3.36 (best proved base) in this model (floor B,
  first-moment argument F2; asymptotic normal tail). With the random stratified grouping, no analysis beats ~3.73 / 3.52
  (rank-2 example F1); R2F + R3'' = 3.82 is within 0.1 of that.
* Reason: lt must grow like alpha log2 x (alpha = 2 random grouping, 1.5 any grouping) because the max over K^2/2 pairs of
  the group sums carries a sqrt(ln K) (any grouping) or ln K (random grouping) factor; the growth sum 2^{-l} lt_l inflates the
  effective base size A_inf = x1 + 2 alpha log2 x1 + O(1) (x1 ~ 16-25).
* Any grouping design has zeta = Omega(1/t) (F3), so alpha >= 1: lt >= log2 x + O(1) is the hard wall of grouping designs.
* Structural change S2 removes the max: only the mean absolute deviation of the pair sums is paid (per-rectangle padding),
  zeta ~ 0.6/t, alpha = 1: c1 3.23 with R3'' (3.00 ideal base), below floor B (3.36) — but only by 0.13: after S2 the base
  (R3 terms; x1 ~ 24) and the alpha = 1 growth are what remains. Below that: only a no-grouping recursion (t = 1) could help,
  and naive per-rectangle padding of a raw network costs O(1) (no averaging over t^2 pairs) [H].
* [H] direction not pursued: the worst iota is uniform (exact grouping, zero excess); any fluctuation of the row network
  needs angular dispersion of the pieces, which costs the adversary quadratically on the circle. A network-adaptive grouping
  bound (variance of the actual pair weights) would make the degenerate excess second order; the linear (node-sum) part stays
  first order (lin2 model: 3.48 with R3, worse than S2).
## Files
notes.md, mylib.py (P2 + zeta_F, zeta_S2, zeta_floor, zeta_randfloor, zeta_lin2, zeta_zero, r3_terms_v2), mydp5.py, mythr5.py,
chain5.py, test_r2f.py/.out, test_s2.py/.out, dp_*.out, thr_*.out, solvers_*.json. No process left running.
