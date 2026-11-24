# G9c (track C): Q1* — hybrid branching, other encodings, star identity. Statements and proofs

Labels: [P] proved here (relative to Appendix D / Lemma 1 of the paper, which is taken as given and was not
re-verified), [L] literature, [H] heuristic, [Q] open. log = log2. e = e_l = Theta(1/sqrt(log n)) is the exponent of the
(K, lambda)-clique recursion at K ~ n^{1/17}. "U" = total edge bits of a network; one network costs
exp2(e(U + gP/K + K^2) + 2gP/K + ...) = exp2(e U + O(n^{1-c})) ((70) of the paper), c = c(C0) > 0.
G graph, n vertices, m = nd/2 edges, max degree Delta <= n^{1/(40 C0)} (vertices of larger degree: there are
<= n^{1-delta} polylog of them, branch on them first). i_t(G) = number of independent sets of size t, alpha = alpha(G).

## 0. Summary

| item | result |
|---|---|
| NEW Thm A [P] | cut-off without parity checks: i_t(G) with U <= 2 nb max{1, log(3m/nb)} + 2g^2 for any nb >= max_{|S|=t} e(S) |
| NEW Thm B [P] | class-restricted Appendix D: i_t(G) in time exp2(O(e h_t log(4m/h_t)) + O(n^{1-c})), h_t = log C(n,t) |
| Cor B1-B4 [P] | #IS in 2^{O(en)} (no log d) for: alpha <= n/(log d loglog d); (n,d,lambda)-graphs with lambda <= d/log^2 d; random d-regular graphs (U <= 13.5 n by Thm A, U = O(n log^3 d/d) by Thm B). **The G7d/G8b benchmark family is not hard.** |
| Idea 1 hybrid | [P] obstruction: decision tree + Prop. R leaves costs Omega(min(e,1/d) n log d) on expanders (Prop H1/H2); on the new benchmark (bipartite) hybrid gains a constant 0.80 only [H] |
| Idea 2 zeta/Moebius | [P] obstruction: bits = log rank, basis independent; profile-restricted disjointness has rank C(N, min(s,s')) (Gottlieb), so a direct network pays log(1/density) per selected incidence: U >= d alpha log(n/alpha) |
| Idea 3 star identity | [P] obstruction: bits = (vertex, foreign cluster) incidences >= n d^2/(d^2 smax/n + lambda^2) - n = Theta(nd) for every clustering with >= d parts; = n(K-1)(1-(1-1/K)^d) for a random K-clustering (measured) |
| open core [Q] | independent sets of size t in (n/(log d loglog d), n/2]; canonical instance: random bipartite d-regular graph, profile |T| ~ n/2d on one side; every method found costs Theta(min(e,1/d) n log d). Equivalent wish: a bit shared by r clusters at cost O(1) instead of r-1 (direct) or log r (Appendix D) ("Q-hyper"). |

Best bound proved for general graphs: unchanged, exp2(O(e n log d)) (Prop. R); log log n is NOT removed.

## 1. Two positive results (my own idea: restrict the class of assignments)

Appendix D pays (kind 1) about k M log(3 sbar) with M = rbar + beta - 1, beta = n + log C(m, nb+1) + 1,
rbar <= 2.41 beta + 1, sbar = 16 m/nb. The "n" in beta is log of the number of assignments (2^n): it bounds
X <= 2^n C(m, nb+1) (number of checks needed) and the modulus 2^{n+1} for F. Independent sets live in a much
smaller set of assignments, {|S| <= alpha}; restricting the sums to a class shrinks beta, hence M.

**Class networks.** For t >= 0 let C_t = {S subset V : |S| = t}. In Step 3 of Appendix D add to a profile the
numbers s_a = |S cap A_a| (a in [g], sum s_a = t) and restrict S_a to settings with exactly s_a true variables.
This multiplies the number of profiles by <= (n+1)^g = 2^{O(n^{8/17} log n)} and changes neither U nor P, B, tau.
So every sum "sum over xi in {0,1}^n" of Appendix D is available as "sum over xi in C_t" at the same cost.

**Theorem A (no checks) [P].** Let nb >= nu_t := max_{|S|=t} e(S). Then
i_t(G) = sum_{S in C_t} sum_{Pi subset E(S), |Pi| <= nb} (-1)^{|Pi|}, and this is a sum of 2^{o(en)} class networks
with only kind-0 marks, each with U <= 2 nb max{1, log(3m/nb)} + 2 g^2.
*Proof.* For |S| = t, e(S) <= nb, so the truncation is void and the inner sum is [e(S) = 0]. Marks: one test per
edge (k = 2), at most nb marked tests among T0 <= m; (59) and the convexity step of the paper give
sum_{a<b} ceil(log C(t_ab, mu_ab)) <= nb max{1, log(3m/nb)} + g^2. ∎
*Certified cut-offs.* (i) d-regular G with second largest adjacency eigenvalue lambda_2:
e(S) <= (t/2)(d t/n + lambda_2 (1 - t/n)) (one-sided mixing: 1_S^T A 1_S <= d t^2/n + lambda_2 |P 1_S|^2). (ii) any G:
X_t(nb) = sum_{S in C_t} C(e(S), nb+1) is a sum of class networks with exactly nb+1 kind-0 marks, and X_t(nb) = 0
iff nb >= nu_t; doubling nb finds a valid cut-off <= 2 nu_t. (iii) alpha need not be known: compute i_0, i_1, ...
and stop at the first t with i_t = 0 (independent sets are closed under subsets, so t = alpha + 1).
Hence #IS in time exp2(O(e nu log(3m/nu)) + O(n^{1-c})), nu = max_{t <= alpha+1} nu_t.

**Theorem B (entropy form of Appendix D / Prop. R) [P].** For every t with 1 <= h_t := log C(n,t) <= m,
i_t(G) is computed in time exp2(O(e h_t log(4m/h_t)) + O(n^{1-c})).
*Proof (changes to Steps 1-4 of Appendix D, k = 2).* Put h = ceil(h_t), l = log(4m/h), nb = max{2, floor(h/l)},
v = the power of two in (nb/2, nb], beta = h + ceil(log C(m, nb+1)) + 1, rbar = min{j : 3^j 2^beta < 4^j},
M = rbar + beta - 1, sbar = 16m/v, p = 16/(15v). Define F, X as in (62) with xi ranging over C_t.
(63) holds verbatim. 0 <= X <= |C_t| C(m, nb+1) <= 2^{beta-1}; if X = 0 then 0 <= F = i_t <= 2^h, so F is determined
modulo 2^{h+1}. Lemma 12 is applied to phi 1_{C_t} with (phi, nu) = (C(|V|, nb+1), beta) and (eta, h+1); M >= nu-1+r.
Lemma 13 with v and the potential argument of Step 2 are unchanged (sets with > nb clauses have > v clauses),
so <= rbar <= 2.41 beta + 1 checks of size <= sbar suffice, load <= tbar + 2 with
tbar <= log m + 1.5 p rbar + 1 = log m + O(beta/nb) = O(log m). Costs: (nb+1) log(3m/(nb+1)) <= (h/l + 3)(l + log l + 1) = O(h + log m);
so beta, M = O(h + log m); kind 1 costs <= 2 M log(3 sbar) <= 2M log(96 m l/h) = O((h + log m) l). So
U = O(h_t log(4m/h_t) + log^2 m) + 2g^2. gP = O(n Delta (tbar+3) + rbar sbar + g) = O(n Delta log n + m log m), as in (69).
Step 4 unchanged. ∎
(Same statement for k-CNFs and any class that is a union of 2^{o(en)} products of node-local restrictions, e.g.
Hamming balls: time exp2(O(e k h log(4km/h))), h = log |class|.)

**Corollaries [P].** Let h* = log C(n, alpha+1) <= (alpha+1) log(en/(alpha+1)).
* B1. #IS(G) in time exp2(O(e h* log(4m/h*)) + o(en)). This is O(en), with no log d, as soon as
  alpha <= n/(log 2d · log log 4d); it is Prop. R when alpha = Theta(n).
* B2. (n,d,lambda)-graphs: alpha <= n lambda/(d+lambda) (Hoffman), so the exponent is O(e n (lambda/d) log^2 d);
  = O(en) whenever lambda <= d/log^2 d; for lambda = O(sqrt d): O(e n log^2 d/sqrt d) = o(en).
* B3. Random d-regular graphs (alpha <= (2+o(1)) n ln d/d whp, first moment; lambda_2 = O(sqrt d) whp [L]):
  Thm B gives U = O(n log^3 d/d); Thm A gives U = O(n log^2 d/sqrt d), and with the first-moment alpha and
  lambda = 2 sqrt(d-1) the explicit value 2 nb log(3m/nb) is <= 13.5 n for every d (maximum near d = 64; exp2 table),
  against 195 n ... 480 n for the same accounting of Prop. R. **So the benchmark family of G7d §4 / G8b
  (random d-regular, d ~ sqrt(log n)) is counted in time 2^{O(n/sqrt(log n))}; indeed 2^{o(n/sqrt(log n))} as d grows.**
  Q1' of G7d (2^{O(n/d)} at d ~ sqrt(log n)) holds on it: exponent O(n log^3 d/d^2).
* B4. Every graph: i_t(G) in 2^{O(en)} for all t <= n/(log 2d loglog 4d).
* Combined measure: exponent O(max_{t<=alpha+1} min{1, e log(4m/h_t)} h_t) (enumeration costs h_t).
**Lemma BC (broadcast) [P].** Let R0 be a vertex cover of G. If every selected vertex of R0 tells the nodes of its
neighbours that it is selected (marks on (vertex, foreign node) tests, one test per object, exact, no cut-off), a
class in which the selected vertices of R0 have total degree <= D costs U <= 2 D max{1, log(6m/D)} + 2g^2, whatever
the vertices outside R0 do (at most 2m tests, at most D marked). For bipartite G with
sides L, R, |R| = N: classes with |S cap R| <= N/(d log d) cost U = O(N) although their entropy is ~ N.

Limit of the method [P]: a family of classes covering all independent sets with 2^{o(en)} classes has
max h >= log i(G) - o(en) >= alpha - o(en). So Thm B cannot give O(en) when alpha >> n/log d.

## 2. Idea 1: hybrid (branch on B, enumerate IS cap B, network on the residual graph)

**Prop H1 [P].** Let G have max degree Delta. A decision tree on the vertex variables (leaf = partial assignment,
residual = undecided vertices not adjacent to a chosen one) all of whose consistent leaves have residuals of
<= k vertices has >= F(n-k) leaves, F(j) = 1 for j <= 0, F(j) = F(j-1) + F(j-Delta-1); F(j) >= x^{j-Delta},
x^{Delta+1} = x^Delta + 1, log x = (1-o(1)) W(Delta)/(Delta ln 2) ~ log(Delta)/Delta · (0.73..0.80 for Delta <= 2^16).
*Proof.* L_k(S) := least number of leaves for G[S]. Restricting a tree to S' subset S (answer 0 outside S') shows
L_k monotone. If |S| > k, the all-zero path queries v_1..v_q with q >= |S| - k; the 1-child of v_i is a tree for
S_i = S \ ({v_1..v_i} ∪ N(v_i)), |S_i| >= |S| - i - Delta, and it is consistent. Induction:
L_k(S) >= 1 + sum_{i=1}^{|S|-k} F(|S|-k-i-Delta) = F(|S|-k). ∎ (Exact DP on n <= 18: exp5, always >= F, equality at k = n/2.)
**Cor H2 [P, for the cost model "Prop. R bound at the leaves"].** G an (n,d,lambda)-graph, lambda <= d/4. Every residual S
has e(S) >= (|S|/2)(d|S|/n - lambda), so its Prop. R bound is >= e |S| log((d|S|/n - lambda)/2), increasing in |S|.
Hence total exponent >= min_k max{ log F(n-k), e k log((dk/n - lambda)/2) } = Omega(min(e, 1/d) n log d)
(k = n/2). At e d = 1 the minimum is ~ 0.31..0.42 (n/d) log d (exp5). **G8b's heuristic claim (iii) is correct and now
proved in this model:** the residual degree does drop to d k/n, but getting k <= n/log d costs 2^{Omega(n log d/d)} leaves.
(Lower bound on the *bound*, not on the running time; trees branching on other objects, e.g. edges
i(G) = i(G-e) - i(G - N[u] ∪ N[v]), obey the same recurrence with (m, 2d^2) in place of (n, d): n log d/(2d) [H].)
With Theorem B at the leaves nothing is needed on this family (Cor B3). On the new benchmark (Section 5), profile tau = 1:
enumerate T cap R1 (R1 a fraction rho of R), broadcast network on the rest:
exponent (N/d) log(ed) [rho + e d (1-rho) e^{-rho}]; at e d = 1 the minimum is 0.801 at rho = 0.443 [H]: a constant.
**Verdict: no setting of |B| gives O(n/sqrt(log n)) at d ~ sqrt(log n) via hybrid alone.**

## 3. Idea 2: other encodings of the pair functions (zeta/Moebius, subset convolution)

The cost of a pair (a,b) is ceil(log) of the rank over Q of Psi_ab (Definition 1), which no change of basis alters.
The Moebius factorisation DISJ[X,Y] = sum_Z (-1)^{|Z|} [Z subset X][Z subset Y] *is* inclusion-exclusion over common
subsets, i.e. Appendix D's marks. Fast zeta transforms save ring operations inside one table (irrelevant: tables
are evaluated by brute force at cost 2^{O(gP/K)} = 2^{o(en)} anyway), not bits.
**Prop Z [P given Gottlieb 1966].** Let H_ab (edges between A_a and A_b) be a perfect matching on N pairs, and let a
profile fix s = |S_a cap N_a(b)|, s' = |S_b cap N_b(a)|, s + s' <= N. Then Psi_ab restricted to the profile contains the
disjointness matrix of s-sets against s'-sets of [N], of rank C(N, min(s,s')) (checked numerically, exp4(3)).
So a network whose nodes hold vertex bits and whose pair functions are the conflict tests pays
c_ab >= log C(N_ab, min(s_ab, s_ba)) under any such profile; for selected sets spread evenly at density rho this is
U >= 2 m H(rho)(1-o(1)) ~ d t log(en/t) for |S| = t (realisability of the even profile by an independent set: [H]).
Without profiles the rank is 2^N (G8b: U = Theta(nd)). "O(1) bits per vertex instead of log d" is impossible
in this class: one pays 1 bit per potential conflict (dense) or log(1/density) per selected incidence (sparse).
Down-closedness of the family of independent sets is exactly what Theorems A/B use (classes of small sets).

## 4. Idea 3: star identity prod_{v in N+(u)} (1 - x_u x_v) = 1 - x_u + x_u prod_{v in N+(u)} (1 - x_v)

The star factor is a CP-rank-2 function of (x_u, y_{u,b_1}, ..., y_{u,b_s}), y_{u,b} = [N+(u) cap S cap A_b = empty]; the rank-2
index is x_u itself. In a pairwise network it costs one bit per cluster b != c(u) that meets N+(u) (star or AND
chain: the same count). So U_star = 2 inc+, inc(+) = sum_u #{b != c(u) : N(+)(u) meets A_b}, and the recursion applies
with lambda = L ~ U_star/K^2 (Lemma 1), cost e U_star.
**Prop S [P].** d-regular G, lambda = second largest absolute eigenvalue, clusters of size <= smax:
inc >= n d^2/(d^2 smax/n + lambda^2) - n, and inc+ >= n d^2/(4(d^2 smax/n + lambda^2)) - n for every orientation.
*Proof.* For fixed u, #{b : e(u,A_b) > 0} >= d^2/q_u, q_u = sum_b e(u,A_b)^2 (Cauchy-Schwarz). sum_u q_u =
sum_b 1_b^T A^2 1_b <= sum_b (d^2|A_b|^2/n + lambda^2 |A_b|) <= d^2 smax + lambda^2 n. Jensen for 1/x. Oriented:
sum_u (d_u^+)^2/q_u >= (nd/2)^2/sum q_u. ∎ So inc = Theta(nd) for every clustering into K >= d parts when
lambda = O(sqrt d). For a uniformly random K-clustering, E inc = n (K-1)(1 - (1-1/K)^d) (measured: exp3, agreement to
1%; local search lowers it by at most 22%). inc ~ nK for K << d, but then e = e(K) = Omega(1/sqrt(log K)) and the
cost e(K) n K >= n: useless; the recursion needs K = n^{Omega(1)} >> d. Star + inclusion-exclusion = edge marks
(sum over nonempty Q subset N+(u) cap S of (-1)^{|Q|+1} = 1), so nothing new there either.
**Reformulation [Q-hyper].** "x_u is one bit shared by the r <= d+1 clusters meeting N[u]". Pairwise networks pay r-1
bits for it (direct) or O(log r) (Appendix D's marks); Q1* asks for O(1). A (K, lambda)-recursion for strings shared by
r nodes at cost e · (shared bits) · O(1), r <= polylog, would give #IS in 2^{O(en)} (n shared bits) and, with the
clause-node network of §6 of the paper, #k-SAT in 2^{O(e n)} after sparsification [H, not checked]. It is false for
r = K (query lower bound 2^lambda for sum_theta prod_w W_w(theta)); no obstruction known for r << eK.

## 5. Refined open core and new benchmark

By B4 + Lemma BC only independent sets with n/(log d loglog d) < t <= n/2 remain, i.e. graphs with large
independent sets. Benchmark: random bipartite d-regular G, sides L, R of size N = n/2;
i(G) = sum_{T subset R} 2^{|L \ N(T)|}. Class |T| = tau N/d, Y = S cap L subset L \ N(T), |L \ N(T)| ~ N e^{-tau}:

| method | exponent for the class (leading order) |
|---|---|
| enumerate T | (N tau/d) log(ed/tau) |
| Lemma BC (T broadcasts) | e · 2 N tau log(3d/tau) |
| Theorem B on (T, Y) | h ~ N H(e^{-tau}) => e · Theta(N log d) |
| Y broadcasts / d-clause marks | e · Theta(N d) |
| Thm A with bipartite mixing | nb = Theta(N) => e · Theta(N log d) |

For tau in [1/log d, O(log d)] all are Theta(min(1/d, e) N log d): at d ~ 1/e again Theta(n log log n/sqrt(log n)).
(tau <= 1/log d: O(eN) by Lemma BC [P]; tau <= 1/sqrt(log d) via 2^{-|N(T)|} = 2^{-d|T|} prod_y #even subsets of the
marked edges at y [H]; large tau: symmetric in Y.) Reason: each y in L needs the OR of d sparse bits (density tau/d)
held by d different clusters; pairwise this is d H(tau/d) = tau log(ed/tau) bits per y, while the cut rank of a
single cluster is at most 1 bit per vertex [P] and a balanced cut has log rank O(N log^2 d/d) [H]. Single-node cut
ranks therefore cannot prove more than U = Omega(n), and no rank bound is available for signed sums of
2^{o(n)} networks (G7d §0). **[Q1**]** count independent sets of random bipartite d-regular graphs, d ~ sqrt(log n), in
time 2^{O(n/d)}.

## 6. Experiments (all in this directory)
* exp1_identities.py / exp1_out.txt: n = 12..16 random regular graphs, brute force over all subsets.
  Thm A identity exact for every t <= alpha+1 with nb = max e(S) and false with nb - 1 (all cases);
  Thm B pipeline (class C_t, cut-off nb in {1,2}, greedy random checks until X = 0, marks truncated at
  M' = nu-1+r, modulus 2^nu, nu = ceil(h)+1) returns i_t in all 42 runs; r ~ log2 X_0.
* exp2_costs.py / exp2_out.txt: U/n as a function of d (asymptotic accounting, same constants for all routes):
  Prop. R 195..484, Thm B 192 (d=4) .. 33 (d=1024) .. 0.25 (d=2^20), Thm A 5.4..13.4..1.2, direct = d.
  n = 60 graphs: alpha, log2 i, lambda_2, spectral and local-search cut-offs (the asymptotic regime nb << n
  is not visible at n = 60, d <= 16: nb_spec ~ 44..58).
* exp3_star.py / exp3_out.txt: incidences vs formula vs local search vs spectral bound (n = 2048, d = 8,16,32).
* exp4_bip.py / exp4_out.txt: bipartite, N = 20: F_t exactly, |L \ N(T)| ~ N e^{-tau}; hybrid constant 0.801;
  Gottlieb rank check.
* exp5_tree.py / exp5_out.txt: Prop H1 against exact minimal trees; growth rates; Cor H2 minimum.

## 7. Self-review (adversarial)
* Thms A/B are transcriptions of Appendix D with the sum over {0,1}^n replaced by a sum over a class. I re-read
  Steps 1-4 and Lemmas 12, 13 for every use of "n" (only beta, the modulus n+1, and (69)); the paper itself was not
  verified, and it states Appendix D for 3 <= k (2-clauses are clauses of <= 3 literals). The identities are checked
  by brute force (exp1), the network contraction (Lemma 1/8, (70)) is not simulated. Constants in "O" not optimised;
  the O(n^{1-c}) term dominates when U << n^{1-c}.
* Cor B3 uses whp facts about random regular graphs (alpha, lambda_2) [L]; the algorithm itself needs no certificate
  (stop rule + computed lambda_2, or the X_t = 0 test). The "13.5 n" is the value of the explicit formula, not a
  theorem with all rounding terms.
* Thm B does not improve the #k-SAT bound: for general k-CNFs the class is {0,1}^n. It removes log d only when
  log i(G) << n/log d is witnessed by local classes.
* Cor H2 and the table in §5 are statements about upper-bound formulas (cost models), not lower bounds on time.
  The bipartite "barrier" is the failure of the methods listed, not an impossibility; the hybrid constant and the
  even-profile realisability in Prop Z are [H]. Prop S bounds one specific encoding (plus G8b's rank bound for
  direct networks); networks with relay nodes are only constrained by cut ranks (Omega(n)).
* Q-hyper and its consequence for #k-SAT are unproved suggestions.
* Nothing here removes log log n in general; Q1* stays open, restricted now to graphs with alpha >= n/(log d loglog d).
