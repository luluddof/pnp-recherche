# G26-B7 notes (lower rung, round 2: constants of the Oct-2026 algorithms)

Status tags: [P] proof here (unreviewed), [N] numerics here, [R] inherited (paper / p2 / p9), [H] heuristic, [C] conjecture.

## 0. Idea list (brainstorm)
1. (R1) Deterministic Lemma 8' by pessimistic estimators (Raghavan): the Bernstein/Chernoff union bound of
   Lemma 8' is a sum of products of per-edge MGFs -> method of conditional expectations, label clauses one by one.
   Deterministic a = 2 => deterministic c1 = Las Vegas c1 (10.6 -> 7.47 p2 chain; 6.81 -> 4.82 p9).
2. (R2) Stratified grouping lemma replacing Lemma 8 (ii)'s Hoeffding: sort nodes by node sum, strata of K
   consecutive nodes, independent random bijection stratum -> groups. Linear (Hajek) part of Y deterministic
   <= 4N/t; degenerate part via Doob martingale over strata + Hoeffding for node-to-partial-group sums.
   zeta = O(sqrt(log g log K)/t) instead of O(sqrt(log(tK^2)/t)).
3. (R3) Ordering lemma (p9's conjecture up to a log): random order, decomposition F(S_j) = j(g-j)mu + (g-2j)A_j + Gamma_j.
   Parabola trick: the linear term has coefficient (g-2j), which vanishes at the peak: contributes <= A^2/mu = O(P log g).
   Degenerate Gamma: martingale X_j = Gamma_j / w_j, w_j = (g-j)(g-j-1)/(g(g-1)). => max cut <= U/4 + O(P log g).
4. (R4) Time-space trade-off at level 0: outer-enumerate a set S of edge strings, sweep the rest. Square S around the
   middle of the cut gives time T = (1 - M)/3 (per lam K^2) at memory M, for every M in [0, 1/4]; M = 0 is the
   divide-and-conquer on the middle cut (T = 1/3, poly space).
5. Poly-space Lemma RS (D&C on the row network, cut U/4 (1 + 1/4 + ...) = U/3).
6. Same stratification trick at the top level (Lemma 8'/T4) -> smaller deviations (a stays 2; minor).
7. Decision vs counting: isolation (Valiant-Vazirani) costs nothing in exponent; parity / mod-p counting lets
   tables live in GF(2)/Z_p: entries of O(1) bits instead of lam K^2 bits -- only lower order. [H: no constant gain]
8. OpenAI side: derandomised Freedman labelling (P3's "not done", ~1.9x in lambda*) -- same pessimistic-estimator
   technique as R1 can make it deterministic; the gain itself is P3's heuristic.
9. Memoisation + R2: entries term 1/K_l becomes a max (p9 chain_variants) -- combine.
10. Smith-mode C1 for Cauchy blocks (OpenAI 37.5): not attempted here (needs the G18d machinery).

## log
- read brief + p2/p8/p9/p3/p10 notes.
- wrote b7lib.py (zeta_new = R2, zs_R3 = R3 explicit), chain_b7.py (chains). b7lib output (per N):
  R2 vs old zeta at x=49: t=2^12 5.4e-2 vs 1.6e-1; t=2^15 6.9e-3 vs 5.8e-2; t=2^19 4.4e-4 vs 1.5e-2.
  R3 vs p9 zs at g=2^16.5: 7.4e-4 vs 8.0e-3 (gamma part 5.3e-4).
- container restarted 01:42 UTC; deadline 06:15 UTC.
- chain1.out (balanced chains, single schedule; no-sawtooth = 2 sqrt(a/mu)):
  clique-sweep old zeta mu 0.284 (repro p9 0.2865) -> new zeta (R2) 0.374: c1(a=2) 7.50 -> 6.54
  row-p9 old 0.351 (repro 0.352) ; row-p9+R2 0.424 ; row-R3 old 0.395 ; row-R3 + R2: mu 0.484 (> p9's conjectural 0.4605)
  -> single-schedule c1 a=2 5.75 (p9: 6.74), no-sawtooth a=2 4.07 / a=4 5.75 (p9: 4.77/6.74).
  diag.py: with R2+R3 Gamma=prod(rho_l/rho)=1.034 only; remaining loss = growth by log t: A_inf=58 vs x1=27.
- space.out: poly space (no sweep): clique D&C base + R2: mu 0.215, c1 single a=2 8.62 (no-saw 6.10);
  brute.out: old poly-space (paper level 0, brute) best: 14.6 (old zeta) / 12.8 (R2). Row D&C level 1 [H const]: 7.75.
  memory M (per lam K0^2) at base: e0=(1-M)/3: c1 single a=2 = 8.58 (M=1/64), 8.18 (1/16), 7.86 (.1), 7.43 (.15), 7.00 (.2), 6.57 (.25 = sweep).
- checks_b7.out: Hajek identity exact (1e-15), drift exact; random orders max cut <= 1.011 U/4 (g=512, 3 kinds);
  stratified grouping hetero network: excess 0.455 -> 0.049 (plain vs stratified), uniform 0.024 vs 0.023;
  R4 discrete K=40: time/mem points on (1-M)/3 + ~0.006.
- dp_b7.py (p9 DP with R2+R3) launched.

## 1. Proofs (written here, unreviewed) -- all on top of p2/p9 [R] dependencies (Prop. 3, Lemma 8 (ii) usage, Thm B)

### R1. Deterministic Lemma 8' (top level, a = 2)  [P, standard technique]
Lemma 8' bounds failure by a union over events E = (cluster pair (a,b), colour class i) and cluster-size events,
each via an MGF (Bernstein) bound  Pr(X_E > thr_E) <= exp(-s_E thr_E) prod_{e in M_i} E exp(s_E 1[e crosses (a,b)]).
Edges of one colour class are vertex-disjoint (a matching of clauses), so the factors depend on disjoint label sets and
Phi(partial labelling) := sum_E exp(-s_E thr_E) prod_e E[exp(s_E 1_e) | partial] is EXACTLY the conditional expectation
of Psi := sum_E exp(s_E (X_E - thr_E)) >= #(failed events) (same for size events, factor per clause).
Phi(empty) = the p_fail bound <= 1/2. Label clauses one at a time, choosing the label (g choices) minimising Phi
(average over the g choices = current Phi, so Phi never increases). At the end Psi = Phi < 1: no event fails.
Each update touches O(g^2 * 9) factors: poly(m) time. Hence the Las Vegas top level becomes deterministic with the same
parameters (a = 2): **deterministic c1 = Las Vegas c1** (p2 chain 10.6 -> 7.47; p9 6.81 -> 4.82; all rows below).

### R2. Stratified grouping lemma (replaces the Hoeffding argument of Lemma 8 (ii))  [P]
Setting: network on g = tK nodes, symmetric c_ab in [0, P/g], node sums c_a <= P, U = sum_a c_a, N = gP.
Construction: sort nodes by c_a, strata S_1..S_t of K consecutive nodes, independent uniform bijections pi_i: S_i -> [K];
group gamma = {pi_i^{-1}(gamma)}: every group has exactly t nodes (as Lemma 8 (ii) requires).
Claim: with prob. >= 1 - d1 - d2, for all gamma != gamma':
   K^2 Y_{gamma gamma'} <= U K/(K-1) + 4N/t + N sqrt(2 Sig L2)/t^2,
   Sig = sum_{k<t} (4 sqrt(k L/2) + 1)^2,  L = ln(2g^2/d1),  L2 = ln(K(K-1)/(2 d2)).
So zeta_new = 4/t + sqrt(2 Sig L2)/t^2 + 1/(K-1) ~ (4 + sqrt(8 L L2))/t, versus Hoeffding's ~ 2 sqrt(ln(tK^2)/(2t)).
Proof. (a) Mean: for i != j, (pi_i^{-1}gamma, pi_j^{-1}gamma') is uniform on S_i x S_j; for i = j uniform on distinct pairs
of S_i. So E Y = sum_{i!=j} c(S_i,S_j)/K^2 + sum_i c(S_i,S_i)/(K(K-1)) <= U/(K(K-1)).
(b) Doob martingale over strata k = 1..t, a = pi_k^{-1}gamma, b = pi_k^{-1}gamma'. With G'_{<k} the part of gamma' built from
strata < k, h(a) := c(a, G'_{<k}) + c(a, S_{>k})/K, the increment is D_k = [h(a) - avg_{S_k} h] + [same for b, gamma]
+ [c(a,b) - E c(a,b)]. Split h(a) = (c_a - c(a, S_k - a))/K + Z_k(a), Z_k(a) = sum_{j<k} [c(a, pi_j^{-1}gamma') - c(a,S_j)/K].
 * linear part: |(c_a - c(a,S_k-a))/K - avg| <= (r_k + P/t)/K, r_k = range of c_a on S_k; sum_k r_k <= max c - min c <= P
   (strata consecutive in sorted order). Total over k, both sides: <= 4P/K, deterministically.  (stratification used here)
 * Z_k(a): sum of k-1 independent terms (over pi_j), each in an interval of length P/g, mean 0: Hoeffding + union over
   (a, gamma', k) (g^2 events): |Z_k(a)| <= z_k := (P/g) sqrt((k-1) L/2) for all, w.p. >= 1 - d1.
 * remaining martingale differences bounded by c_k = 4 z_k + P/g on that event; Azuma for the stopped martingale + union over
   K(K-1)/2 pairs: deviation <= sqrt(2 sum c_k^2 L2) w.p. >= 1 - d2.  Multiply by K^2 and use K^2 P/g = N/t^2, KP = N/t. QED
Exhaustive search over groupings (cost g log K, counted in p2) finds a grouping meeting the bound.
Check (checks_b7.out): heterogeneous network K=16, t=64: max excess 0.455 (plain random grouping) vs 0.049 (stratified).

### R3. Ordering lemma for the row sweep (p9's conjecture, up to a log)  [P]
Same setting, uniform random order, S_j = first j nodes, F_j = c(S_j, S_j^c). Off-diagonal ANOVA:
mu = U/(g(g-1)), alpha_a = (c_a - (g-1)mu)/(g-2) (sum 0), gamma_ab = c_ab - mu - alpha_a - alpha_b (zero row sums). Then
   F_j = j(g-j) mu + (g-2j) A_j + Gamma_j,  A_j = sum_{S_j} alpha,  Gamma_j = gamma(S_j, S_j^c)     (exact; checked 1e-15).
Parabola trick: with k = g/2 - j, j(g-j)mu + (g-2j)A_j = g^2 mu/4 - k^2 mu + 2kA_j <= g^2 mu/4 + A_j^2/mu.
(Pad c_ab by a constant so that U = nu N exactly: changes neither alpha nor gamma, only increases cuts.)
 * |alpha| <= P/(g-2); Hoeffding w/o replacement + union over j (using A_j = -A_{g-j}): A^2/mu <= ln(g/dA)(g-1)N/(nu (g-2)^2).
 * Gamma: E[Gamma_{j+1} | S_j] = Gamma_j (1 - 2/(g-j)) (checked), so X_j = Gamma_j/w_j, w_j = (g-j)(g-j-1)/(g(g-1)) is a
   martingale; |dX_i| <= 4 G_i/w_{i+1}, G_i = max_v |gamma(v,S_i)| <= R sqrt(i LG/2) (Hoeffding, R = P/g + 2P/(g-2),
   union over g^2/2 events); Azuma (stopped) for j <= g/2, reverse order for j >= g/2 (Gamma(S) = Gamma(S^c)).
 => max_j F_j <= U g/(4(g-1)) + O(P log g)   (explicit: b7lib.zs_R3; per N: 7.4e-4 at g=2^16.5 vs p9's 8.0e-3).
So zs(g) = O(log g/g) instead of O(sqrt(log g/g)). Empirically (checks_b7.out) random orders give <= 1.011 U/4 at g=512.

### R4. Time-space trade-off for the level-0 sweep  [P]
Outer-enumerate the edge strings of a set S, then run p2's sweep on the rest: time 2^{lam(|S| + max_j |(C_j u N_j) - S|)},
memory 2^{lam max_j |C_j - S|}. With S = {u,v}: K/2 - d <= u < K/2 <= v < K/2 + d (continuous limit, delta = d/K):
cut at x = 1/2 - eps loses (delta - eps) delta, max over eps at eps = delta/2: memory 1/4 - 3 delta^2/4, time 1/4 + delta^2/4
(per lam K^2). delta = 1/2 = the whole middle cut; recurse in both halves (they are independent once S is fixed):
   **time T(M) = (1 - M)/3 + O(1/K) at memory 2^{M lam K^2}, every M in [0, 1/4]**; M = 0: divide and conquer on the middle
cut, d(K)/K^2 -> 1/3 (d(2..2048)/K^2 = .25,.3125,.328,.332,.333...), polynomial space.
(Discrete check K = 40: points (mem, time) on (1-M)/3 + 0.006.)  Optimality of the line: OPEN (no lower bound attempted).
Whole algorithm: the paper's recursion above level 0 is on-demand (poly space) [R, p8], and level 0 carries f ~ 0.96-0.99
of the exponent, so memory/time exponent ratio ~ 3M f/(1 - M): poly space at M = 0, 0.95 at M = 1/4 (= p8's caveat).
- dp_new.out (Lemma U unbalanced levels + R2 + R3, DP of p9): limsup c1 a=2: 4.10/4.11/4.12 (ranges 2k-4k,4k-8k,8k-16k);
  a=4: 5.82. Clique-sweep base + R2 unbalanced: 4.67 (a=2). => with R1: deterministic AND Las Vegas c1 = 4.12
  (was 6.81 det / 4.82 LV in p9).
- freedman_H.out [H]: a variance-based (Freedman) R2 would give mu 0.543 (no-saw c1 a=2 3.84); floor of this route
  (zeta = 4/t only) mu 0.618 (3.60). So better grouping concentration is nearly exhausted; the remaining loss is the
  x-growth by log t (A_inf ~ 58 vs x1 = 27) and the leaf mu1 = 32/x1.
- heuristic threshold scaling: m* ~ (c1)^2: p2 E2 4740 (r=4.26 vs 2^n) * (4.12/7.47)^2 ~ 1.4e3 [H, not recomputed].
- dp_poly.out: POLY SPACE (clique D&C base, Lemma U unbalanced, R2), limsup c1 a=2: 6.16 (a=4 8.71).
  Old poly-space best (brute base, p2 accounting): 14.6 LV / 20.7 det (single schedule). With R1: deterministic 6.16.

## 2. Summary table (asymptotic model of p2/p9: lam -> inf, kappa/c_ev dropped; all [R] deps of p2/p9)
| variant | det. c1 | LV c1 | memory |
|---|---|---|---|
| p9 best (row sweep + unbalanced) | 6.81 | 4.82 | 2^{~0.95 E} |
| + R1 (derandomised Lemma 8') | 4.82 | 4.82 | same |
| + R2 + R3 (DP, Lemma U) | **4.12** | **4.12** | 2^{~0.98 E} (row sweep) |
| clique-sweep base + R2 + R1 (DP) | 4.67 | 4.67 | 2^{~0.95 E} |
| poly space: D&C base + R2 + R1 (DP) | **6.16** | 6.16 | poly |
| poly space, old (brute base, p2) | 20.7 (T4) | 14.6 | poly |
| intermediate memory M (single schedule, clique base, a=2) | 8.58 (M=1/64) .. 8.18 (1/16) .. 7.00 (.2) .. 6.57 (1/4) | | 2^{~3M/(1-M) E} |
| [H] Freedman version of R2 / floor of grouping route | 3.84 / 3.60 (balanced no-sawtooth) | | |

## 3. Ranked remainder
1. Second-agent review of R2/R3 proofs (Azuma with stopping, the linear-part bound sum r_k <= P, padding in R3).
2. Recompute thresholds m* with R1-R3 in p2's full machinery (heuristic ~1.4e3 vs 4740).
3. Freedman/Hanson-Wright version of R2 (3.84 [H]); beyond that the x-growth log t is the wall (floor ~3.6).
4. Lower bound for the time-space line T + M/3 >= 1/3 (pebbling-type) -- open.
5. Poly-space Lemma RS (nested-bisection version of R3) -- [H] 7.75 single schedule vs 8.62.
6. OpenAI side: pessimistic-estimator derandomised Freedman labelling (P3's ~1.9x in lambda*), Smith C1 (not touched).
No process left running; nothing written to /home/user/pnp-recherche.
