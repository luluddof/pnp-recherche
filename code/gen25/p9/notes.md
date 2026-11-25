# G25-P9: can the sweep idea go further? Power of log vs constant in Theorem A'

Tags: [P] proof written here (not yet reviewed by a second agent), [N] own numerics (scripts in this folder),
[R] taken from the paper or from p2 as written (not re-checked), [H] heuristic, [C] conjecture.
All statements are about the design of Theorem A (Claude's Oct-2026 preprint, App. A/B) and its variants. None of
them is a lower bound for 3-SAT. Notation as in p2: x = log2 K, W = number of tests, e = exponent per K^2 lam.

Files: closed_form.py/.out (task 1), dp_unbalanced.py/.out (levers b, d: DP over solver trees),
chain_exact.py/.out (exact long chains, mu_inf), chain_variants.py/.out (memoisation, optimistic ordering),
checks.py/.out (tw(L(K_n)), general network sweep vs brute force, ordering bound, App. B frontier).
thmA_explicit.py, asym_c1.py, sweep_check.py: unmodified copies from p2.

---------------------------------------------------------------------------------------------------------
## 1. Where sqrt(log m) comes from (task 1)  [P, given p2's recursion]

p2's recursion: e_{l+1} = rho_l e_l + 1/K_l, rho_l = rho_bar(lam_l) + 1/g_l + zeta1_l + zeta2_l;
x_{l+1} = 2 x_l + 2 log t_l (g_l = t_l K_l, K_{l+1} = g_l^2). Unrolled:
   x_L = 2^L A_L,  A_L = x0 + sum_{l<L} 2^{-l} log t_l ;
   e_L = e0 Gamma_L rho^L + sum_j K_j^{-1} prod_{j<l<L} rho_l,   Gamma_L = prod_{l<L} (rho_l/rho), rho = 1/sqrt2.
The top uses the largest L with a x_L + b <= log m (a = 9 paper, 4 Lemma T4, 2 Lemma 8'). So rho^L = (A_L/x_L)^{1/2} and
   exponent = 2 W e_L (1+o(1)) = 2 e0 Gamma sqrt(a A) * s * W / sqrt(log m),   s in [1, sqrt2] (sawtooth: x_L can sit
   anywhere in ((log m)/(2a), (log m)/a]).
With a general level (gain rho, growth x_{l+1} = beta x_l + o(x_l)) the same unrolling gives
   exponent = 2 e0 Gamma (aA)^alpha W/(log m)^alpha,  alpha = log(1/rho)/log(beta)   (= 1/2 here: rho = 2^{-1/2}, beta = 2).
**POWER is controlled only by (rho, beta)** and by the requirement log t_l = o(x_l) (the Chebyshev grouping of the
main text has log t = 7x, beta = 16: alpha drops to 1/8, as in the paper's Thm 2). **CONSTANT** = e0, Gamma (zeta and
1/g overheads), A (effective base size), a, and the sawtooth factor s.

Cleaner invariant. Put Phi := 1/e^2 and mu := Phi/x ("figure of merit"). Then exponent = 2W/sqrt(mu_top x_top) and
   c1 = limsup 2 e_l sqrt(a x_{l+1}) = 2 sqrt(2a/mu_inf)   (single schedule),  2 sqrt(a/mu_inf) without sawtooth.
closed_form.out for p2's E schedule (x0 = 9): mu = 1.76 (l=0), 0.59 (l=1), 0.41, 0.34, ... -> mu_inf = 0.2865;
reproduces c1 = 7.473 (a=2) / 10.568 (a=4). The sweep base has Phi0 = 16 for every K0; with no overhead at all
(K0 = 2, x0 = 1, mu = 16) c1 would be 2 sqrt(2a/16) = 1.0 (a=2). **So the factor ~7.5 between the actual constant
and this ideal is all overhead, mostly the grouping factor t_l: A_inf = 47 against x0 = 9 (log t_0 = 15.6 > x0).**

---------------------------------------------------------------------------------------------------------
## 2. The Pythagorean law (key structural fact; answers tasks 2, 3a, 3b and the verdict)

**Lemma U (unbalanced grid level)** [P, on top of the Prop. 3 / Lemma 8 (ii) constructions, R]. Let solvers A, B handle
(K_a, .)- and (K_b, .)-cliques with exponents e_a, e_b. Put the K' = g h nodes of a (K', lam)-clique on a g x h grid,
g = t_a K_a, h = t_b K_b, expand every rectangle with App. B pieces, and for each iota merge rows (row network on g nodes,
grouped into K_a groups by Lemma 8 (ii), solved by A) or columns (h nodes, K_b groups, B), whichever has the smaller
WEIGHTED cost e_a(U_row + z_a N) versus e_b(U_col + z_b N), N = K'^2 lam. Then for every w in [0,1]
   e' <= r(lam) sqrt((w e_a)^2 + ((1-w) e_b)^2) + w e_a (1/h + z_a) + (1-w) e_b (1/g + z_b) + 1/min(K_a, K_b) (+ lower order),
   r(lam) = max_s sqrt(c_s^row^2 + c_s^col^2)/(2 lam);  with w = e_b^2/(e_a^2+e_b^2) the main term is r e_a e_b/sqrt(e_a^2+e_b^2).
*Proof.* Nothing in Prop. 3 up to (36) uses g = h (rectangle = 2 rows x 2 columns; the identity is rectangle-local).
Row network: P = g h^2 lam = N/g, c_{a,a'} <= h lam + C(h,2) 2 lam <= P/g, so Lemma 8 (ii) applies; the straight edges give
U_row <= N/h + 2 sum_q c^row_q, similarly U_col <= N/g + 2 sum_q c^col_q, and 2|R| 2lam <= N. min(X, Y) <= wX + (1-w)Y, and per
rectangle w e_a u + (1-w) e_b v <= sqrt((w e_a)^2 + ((1-w)e_b)^2) sqrt(u^2+v^2) (Cauchy-Schwarz). Entry costs 2^{N/K_a} or
2^{N/K_b}. The algorithm needs computable upper bounds on e_a, e_b (p2's schedule numbers are). QED.
For w = 1/2, e_a = e_b this is p2's level (rho_bar <= r/sqrt2; checks.out (4): r(lam) = sqrt2 rho_bar to 5 digits,
1.0177 at lam = 2^10, 1.00004 at 2^20, i.e. the App. B costs lie on the quarter circle u^2 + v^2 = (2 lam)^2 [N, HR leading term]).

**Proposition P (no power gain inside the architecture)** [P, design-level]. Model M: solver trees whose internal nodes are
grid levels of Lemma U (any g, h, balanced or not; children of any kind: recursive solvers, sweeps of the row network,
memoised or not) and whose leaves are base solvers. In M:
 (i) every level satisfies Phi' <= (Phi_a + Phi_b)(1 + o(1)) while x' >= x_a + x_b;
 (ii) hence Phi(root) <= mu0 x(root) (1+o(1)), mu0 := sup over leaves of 1/(e_leaf^2 x_leaf);
 (iii) leaves that are variable-elimination (tree-decomposition) algorithms have e_leaf >= ~1/4 and x_leaf >= 1, so mu0 <= 16;
 (iv) the top (K, L)-clique has K^2 L <= 2W(1+o(1)) with L >= 1, so x_top <= (1/2) log(2W) and a >= 2;
 => time >= 2^{2W e_top} with e_top >= 1/sqrt(mu0 x_top): exponent >= 2 sqrt(2) W / sqrt(mu0 log(2W)) (>= W/sqrt(2 log 2W) for
    elimination leaves). The power of log is 1/2 for every tree in M.
*Proof of (i).* Lower bound for the design: the worst iota dominates (#iota = 2^{O(g^2h^2 log lam)} = 2^{o(N)}); take iota = s
on every rectangle. Then the row network is uniform with c = h lam + C(h,2) c_s^row, the column network uniform with c_s^col,
and the design widths satisfy c_s^row >= log2 p(2h_lam + s) (they are ceil-logs of partition sums), so by Hardy-Ramanujan
(c_s^row, c_s^col)/(2 lam) -> (cos phi_s, sin phi_s) with phi_s dense in [0, pi/2]. Children cost at least e times their input,
so e' >= max_s min(e_a cos phi_s, e_b sin phi_s)(1-o(1)) -> e_a e_b/sqrt(e_a^2 + e_b^2). (Same for a sweep child: with s fixed the
row network is uniform, every order's middle cut is (g^2/4) c.) (ii) by induction over the tree. (iii) the variables of a
(K,lam)-clique are its K(K-1)/2 edge strings, two interact iff they share a node: the primal graph is the line graph L(K_K),
so an elimination order costs >= 2^{lam (tw(L(K_K))+1)}; checks.out (1): tw(L(K_n)) + 1 = s(n) = floor((n+1)^2/4) - 1 exactly
for n = 3..6 (exact subset DP) - **the sweep is an optimal elimination order**; for all n, tw(L(K_n)) = n^2/4 + O(n) is, as
far as I recall, Harvey-Wood (J. Graph Theory, 2015) [not re-checked]. Since s(K) >= K^2/4 for K >= 2 (K = 2: e = 1/4,
x = 1), elimination leaves have e_leaf >= 1/4 - o(1), Phi_leaf <= 16(1+o(1)), x_leaf >= 1, so mu0 <= 16(1+o(1)).
(Algebraic leaves are not covered: e.g. K0 = 3 with fast matrix multiplication, e0 = omega/9, gives mu = 9.1 < 16.)
(iv) L >= 1 and Lemma 8's K^2 L <= U + 2 zeta N + K^2 with U = 2W. QED.

**Precise reason (the analogue of G18d's inequality):**
      1/e'^2 <= 1/e_a^2 + 1/e_b^2   and   log K' >= log K_a + log K_b,
i.e. Phi = 1/e^2 is at most additive in log K, because the App. B row/column costs lie on the circle u^2 + v^2 = (2 lam)^2
(log p(n) ~ pi sqrt(2n/3), the "gain D^{-1/2}" law of G18c in its weighted form). Depth, number of levels, growth K_l -> K_{l+1},
unbalanced splits, sweeps inside levels and memoisation only move along this law; they change mu (a constant), never the
power. Lowering e0 is a leaf improvement (mu0), also a constant. A power change needs a level whose cost frontier lies
strictly inside the quarter circle for some weighting (rho < 1/sqrt2: the crossing trade-off again), or leaves with
e0^2 log K0 -> 0 (themselves better than sqrt(log)). The same holds for D-direction levels with frontier on the unit sphere
(Phi' = sum Phi_i).

---------------------------------------------------------------------------------------------------------
## 3. The levers one by one (task 2, 3)

**Task 2, sweep inside a level with pieces.** (a) Sweep of the row network instead of grouping + recursion: e = (1/4) x
(row-network cost), i.e. Phi = 16 + 16 at x = 2 log g - it is a sweep leaf one level up, no new rho. (b) Sweep of the grid
POINTS in any order, with pieces on any subset of rectangles: for a bisection S of the g^2 points the plain cut is
|S||S^c| lam = N/4 for every S, and a rectangle is row- or column-split by S for at most |R|/2 rectangles (from
lam #(1-3 splits) + 2 lam #(row + col splits) + 0 #(diag) + straight = |S||S^c| lam); pieces only save on row/col-split
rectangles, so the middle cut is >= rho N/4 for the iota with s = 0 everywhere [P, assuming a piece costs >= the original
crossing on 1-3 and diagonal splits, H for App. B]. No lower effective rho, and no extra depth.
**Lemma RS (row sweep at level 1, no grouping)** [P]. For a g x g grid (or g x h), per iota sweep the row network on g nodes
in the best of all g! orders (search cost g log g, lower order). Random order + Hoeffding without replacement (events: prefix
and suffix node sums, 2g; node-to-prefix sums with range P/g, 2g^2) give some order with
   max_j cut_j <= U/4 g/(g-1) + N zs0(g),  zs0 = sqrt(ln(2g^2+2g+1)/(4g)) (1 + 1/g),
and a step costs 2^{cut_j + c_{j+1} + P} <= 2^{cut_j + 2N/g}. So e1 <= (rho_bar + 1/g)/4 + zs(g), zs = zs0 + 1/(4(g-1)) + 2/g.
checks.out (2): the general-network sweep equals brute force on 12 random Definition-2 networks; (3): on uniform, random,
two-block and heavy-half weights (g = 64..1024) the best of 20 random orders has max cut = (1.000-1.016) U/4, while the proved
bound is 1.24-2.97 x U/4: the bound is loose [N].
chain_exact.out: row sweep at g = 2^16.5 (x1 = 33, e1 = 0.1848, mu1 = 0.89) then p2's recursion: **mu_inf 0.2865 -> 0.3521**,
c1 (single schedule) 7.474 -> 6.741 (a=2), 10.569 -> 9.534 (a=4). The gain is the removal of the first grouping (log t_0 = 15.6).
**(a) growth K_{l+1}(K_l).** Not forced: Lemma U allows any K' = g h with g >= K_a, h >= K_b. Neutral for the main term
(Prop. P). The only real constraint is log t_l = o(x_l) (else alpha drops: Chebyshev log t = 7x gives alpha = 1/8). Slower
growth = more levels = same Phi/x law.
**(b) unbalanced splits.** Exactly neutral in the main term (Pythagoras is an equality for App. B), but they let the solver
sizes fill every x, which removes the sawtooth sqrt2 of the limsup. dp_unbalanced.out (DP over trees, bins of ratio 1.02,
log m/a in [8000, 16000]): clique-sweep base: limsup c1 = 5.35 (a=2) / 7.57 (a=4), vs 7.47/10.57 single schedule (a family of
balanced schedules with different x0 already gives 6.08/8.60 in the same run).
**(c) memoisation / DP over the recursion tree.** The main term is the number of DISTINCT cut states of one bottom sweep along
one worst iota-chain (#iota = 2^{o(N)}), which memoisation cannot reduce; it only turns the additive entry terms 1/K_l into a
max (precompute all 2^{(K-1)L} entries of a level's tables, memory 2^{rho N/K}). chain_variants.out: clique base mu_inf
0.2865 -> 0.3009 (x0 can drop to 7); with the row sweep: no change (K1 = 2^33 already). Constant only, no power.
**(d) sweep + piece expansion where best.** By Prop. P a sweep is best used as a leaf of smallest size that the balancing
overhead allows; concretely at level 1 directly on the row network (Lemma RS) - that is the lever above. Using it at a higher
level (larger leaf) lowers mu (Phi = 32 at larger x).
**Re-embedding between levels (b3 G1-2)** [P-sketch / H]. A unary node "split each of the K vertices into g pieces, solve the
(Kg)-clique" gives Phi(x) = Phi(x'')/c^2 <= mu x''/c^2 (c = bit blow-up, x'' = x + log g), so it gains iff x''/x > c^2. But a
vertex table is arbitrary on B = (K-1)L bits; a network on g pieces whose piece domains have b bits has <= g 2^b free entries,
so some piece needs b >= B - log g, i.e. its incident pairs carry ~B bits, width ~B/g per pair instead of L/g^2; padding to a
clique gives c >~ g/K, and x''/(x c^2) <~ (1 + log g/log K) K^2/g^2 << 1. (Rigorous floor: any (K,lam)-clique algorithm with
arbitrary tables reads K 2^{(K-1)lam} entries, so e >= (K-1)/K^2.) Dead, consistent with b3's guess.

---------------------------------------------------------------------------------------------------------
## 4. Best constant (all [R] dependencies of p2; asymptotic model: lam -> inf, kappa/c_ev terms dropped, as asym_c1.py)

| variant | mu_inf | c1, a=2 (Las Vegas) | c1, a=4 (deterministic) |
|---|---|---|---|
| p2 E2/E1 (single schedule) | 0.2865 | 7.474 | 10.569 |
| + memoisation (clique base) | 0.3009 | 7.29 | 10.31 |
| + row sweep at level 1 (Lemma RS) | 0.3521 | 6.741 | 9.534 |
| + unbalanced levels (Lemma U), clique base | 0.288 | 5.35 (DP range) | 7.57 |
| **row sweep + unbalanced levels** | 0.352 | **4.82** (DP: 4.817; limit ~4.83) | **6.81** (DP: 6.812; limit ~6.83) |
| (conj.) optimal-order lemma, zs = O(1/g) | 0.4605 | 4.17 | 5.89 |
| ideal: no balancing overhead (mu = 16) | 16 | 1.0 (sawtooth) / 0.71 | 1.41 / 1.0 |

So the constant improves again by a factor 0.645 (10.6 -> 6.8 deterministic, 7.47 -> 4.8 Las Vegas); "DP range" = the
limsup over log m/a in [8000, 16000]; the limit is 2 sqrt(a/mu_inf) times the granularity loss (~1.3% at bin ratio 1.02).
Thresholds m* were NOT recomputed with these levers.

---------------------------------------------------------------------------------------------------------
## 5. Verdict
No lever changes the power of log m. Precise reason: in the whole family M (any split shapes, any number of levels, sweeps
at any level, memoisation), Phi = 1/e^2 is at most additive in log K (Lemma U + Prop. P), because the App. B costs sit on the
quarter circle u^2 + v^2 = (2 lam)^2; and log K <= (1/2) log(2W) at the top. Exponent >= 2 sqrt2 W/sqrt(mu0 log 2W), with
mu0 <= 16 for elimination-type leaves (sweep = optimal elimination order, exact for K <= 6). Changing the power needs either a
piece frontier inside the circle (rho < 1/sqrt2 for some weighting: the crossing trade-off) or a leaf with e0 = o(1/sqrt log K0).
What the sweep idea DOES buy: constants. The sweep turns e0 = 1/2 into the elimination optimum 1/4 at the leaves, and it also
removes the first grouping (Lemma RS); unbalanced levels remove the sawtooth. The remaining gap to the model's ideal (c1 ~ 0.7
vs 4.8) is entirely the balancing overhead (grouping factor t_l and zeta_l, and the sweep-order concentration): a sharper
balancing/ordering lemma is the only sizable constant lever left (e.g. [C] optimal-order lemma: 4.8 -> 4.2).

## 6. Self-review / caveats
* Everything inherits p2's [R] dependencies (Prop. 3, Lemma 8, Thm B pieces, rho_bar from the HR leading term, not certified).
* Prop. P is a statement about the design (worst iota, children charged at their exponent), not a lower bound for algorithms;
  its step "c_s^row >= log2 p(2h+s)" uses the paper's definition of the widths, and the circle limit uses HR asymptotics.
* Prop. P (iii) for general n relies on a recalled theorem (Harvey-Wood); verified exactly only for n <= 6.
* Lemma U's proof is a sketch reusing Prop. 3/Lemma 8 verbatim with g != h; no end-to-end small-instance run of an
  unbalanced level was done (the identity is rectangle-local; only costs change). Lemma RS's sweep was tested end to end on
  tiny general networks; the ordering bound by sampling.
* DP numbers use bins (ratio 1.02) that keep the max-mu solver per bin (pruning: results are achievable, maybe not optimal),
  a t-grid with step 1 in log2 t, w = e_b^2/(e_a^2+e_b^2) (valid, maybe not optimal), lam -> inf, and finite x <= 2^15.
* No process of mine left running; nothing written to /home/user/pnp-recherche.
