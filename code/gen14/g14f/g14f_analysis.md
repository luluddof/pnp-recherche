# G14f: Q1* on the bipartite benchmark — clause/#CSP encodings, Fourier, cut ranks, where the cost sits

Labels: [P] proved here relative to the paper's Lemma 1 / Appendix D (untrusted, not re-verified), [L] literature,
[H] heuristic, [Q] open. log = log2. e = Theta(1/sqrt(log n)) (exponent of the clique recursion), one network of
total U costs exp2(eU + O(n^{1-c})). Benchmark (G9c §5): random bipartite d-regular G on L ∪ R, |L| = |R| = N = n/2,
d ≈ 1/e ≈ sqrt(log n). Target: exponent O(eN) = O(N/d). Known: Theta(min(e,1/d) N log d).
**Result: no improvement. log log n is NOT removed.** What is new is listed in §0.

## 0. Summary
| item | status |
|---|---|
| T1: one-sided #CSP form (N variables, N weighted width-d "clauses", density 1) | [P] covered by Lemma 1 + the §6 network (arbitrary integer W_a), NOT by Theorem D verbatim (its cut-off needs the signed sum). Cost e·2N(d−1): exponent ≈ 2N at ed = 1, worse than the trivial 2^N. |
| T1: 3-CNF chain / tree gadgets | [P] U ≥ 2N(d−1) for every gadget (the original variables keep their d occurrences); chain: e·4N(d−1). Useless. |
| T2: Fourier Z_{d+1} / Z_2^d characters | [P] dominated: the factor 1+∏(1−x_u) has CP rank exactly 2, index = [v ∈ T]; Fourier has index of log(d+1) (resp. d) bits shared by the same d clusters: e·2N(d−1)⌈log(d+1)⌉ (resp. e·2Nd). Also complex weights are not allowed (integer tables). |
| T2: even-subgraph (Holant) form 2^{−|N(T)|} = 2^{−d|T|}·#even subgraphs | [P] identity (checked); every edge variable is shared by exactly 2 clusters, but there are τ²N/2·Θ(1) marks of cost log(d/τ): gain only for τ ≤ 1/sqrt(log d) (G9c). |
| localisation of the hard core | [P/H] with classes (|T|,|A|) = (τN/d, αN/d): everything except {τ ∈ [1/log d, O(1)], α/d = Θ(1)} (and its mirror image) costs o(e N log d). Hard core = Σ_{|T|=τN/d} 2^{−|N(T)|}, τ = Θ(1): a T-only sum with entropy h_T = (N/d)τ log(ed/τ). |
| cut ranks | [P] rank of F(A)=2^{|R∖N(A)|} across a cut (L1,L2) = rank of C[B1,B2] = #{T : L∖N(T) = B1∪B2} ≤ #traces; balanced cut: ≤ 0.04N at d = 1024 (Theta(N log² d/d) [H]); single cluster: full. Cut-rank lower bounds on U cannot exceed N (A-holding) resp. h_T (T-only). |
| obstructions | [P] (O1) zero-pattern rank: any pair function with the support of the conflict test costs ≥ log C(s+s′, s) — weights/signs/basis do not help in the hard class; (O2) raw-bit relays ("gathering"): U ≥ 2(c·d−1)·N·H(τ/d) under vertex expansion; [H] the same for OR-aggregation trees. So every network whose messages are restrictions or ORs of T pays d × entropy(T). |
| loophole left [Q] | messages with cancellation (as Appendix D's parity checks) for a T-only sum whose every term is non-zero; equivalently Q-hyper = weighted #CSP with k variables of degree r ≈ 1/e in exponent O(e·k) instead of e·2(r−1)·k. |

Best exponent on the benchmark: unchanged, Theta(min(e, 1/d)·N log d) = Theta(n log log n / sqrt(log n)).

## 1. Task 1: the one-sided form as a weighted d-CSP
i(G) = Σ_{A⊆L} ∏_{v∈R} (1 + [N(v)∩A = ∅]) = Σ_ξ 2^{|V(ξ)|}, V(ξ) = clauses of the monotone d-CNF
Φ_R = ∧_v (∨_{u∈N(v)} x_u) violated by ξ (N variables, N clauses, width d, density 1). (exp1: checked.)
* **Is the weighted form covered?** Definition 2 allows arbitrary integer W_a with |W_a| < 2^B. In the network of §6
  (nodes = clusters of clauses, one bit per occurrence, equality tests between consecutive occurrences) replace
  "[every clause of cluster a has a true occurrence]" by ∏_{E in a} w_E(bits of E), any integer weights:
  U, P unchanged, B ≤ m̄·log max w + 1. So weighted #CSP of arity k with m constraints costs exp2(2ekm + l.o.t.) as in
  (45) [P]. Theorem 1/A(c) is stated for 0/1 clauses only, but nothing else is used. Theorem D (bound in n) is NOT
  available: its cut-off (62)-(63) uses η(ξ) = Σ_{|Π|≤n̄}(−1)^{|Π|} = [V(ξ)=∅] on |V(ξ)| ≤ n̄ and kills the other ξ; here
  the sum is unsigned, Σ_{Π⊆V(ξ)} 1 = 2^{|V(ξ)|}, and the ξ with many violated clauses carry weight. Only the branch
  "no checks" survives: marks Π = T, (d−1) tests per mark.
* **Cost.** U = 2·(occurrences − variables) = 2N(d−1); exponent 2e(d−1)N ≈ 2N at ed = 1: worse than enumerating A ⊆ L.
  With k = d, m = N the m-bound is O(d·N/sqrt(log N)) = O(N): useless. Class-restricted (|T| = t): 2(d−1)t·max{1, log(3N/t)}
  = Lemma BC of G9c.
* **Gadgets.** Replace each R-factor by any gadget with auxiliary variables (prefix-OR chain y_i ≡ y_{i−1} ∨ x_i with a
  free weight-2 bit when y_d = 0; trees). In the clause network the tests are Σ_var (occ − 1); every original variable x_u
  still occurs once in each of its d gadgets, so tests ≥ N(d−1) and U ≥ 2N(d−1) for every gadget [P]. Chain: n′ = Nd
  variables, ≈ 3N(d−1) 3-clauses: m-bound exponent 2e·3·m′ ≈ 18 e N d; n-bound (Prop. R, density 3): O(e N d);
  network count U = 4N(d−1). The dual form (variables = T, factors at u ∈ L) and the 2-CNF form (Nd edge tests) give the
  same Nd. A single factor needs d−1 pairs with c ≥ 1 in ANY pairwise network (the function is not a product across any
  bipartition of its d clusters, so they lie in one component of {c_ab ≥ 1}) [P]; this does not extend to N factors.

## 2. Task 2: Fourier / character encodings
[Σ_{u∈N(v)} x_u = 0] = (1/q) Σ_{j∈Z_q} ∏_u ω^{j x_u}, q = d+1 (exp1: checked numerically). This is a rank-(d+1)
edge-factorised expansion whose index j_v (log(d+1) bits) must be known to the d clusters of N(v): pairwise cost
2(d−1)⌈log(d+1)⌉ per v. But [Σ x_u = 0] = ∏_u (1−x_u) is already rank ONE, and f_v = 1 + ∏(1−x_u) has CP rank exactly 2
(it is not a product), with index s_v = [v ∈ T]: the rank-2 expansion IS the original two-sided IS sum. So [P]:
every edge-factorised expansion of f_v has ≥ 2 terms, 2 is attained, and its index is the bit s_v shared by d clusters.
Fourier is strictly dominated (log(d+1) shared bits instead of 1), as is the Z_2^d expansion 2^{−d} Σ_Γ (−1)^{|Γ∩A|}
(one bit Γ_u per edge: Nd bits; small-bias families: index of O(log d) bits again shared by d clusters).
Truncation to low frequencies is pointless (exact rank 2 exists), and in the hard class |A ∩ N(v)| is Bin(d, Θ(1)), not small.
"log r per shared bit" does not arise here: the price is (r−1)·H(bit), r = d, H = H(τ/d) under profiles.
The one genuinely pairwise form is the Fourier transform of the whole sum (Holant): i(G) = 2^N Σ_T 2^{−d|T|} ∏_u #{even
subsets of N(u)∩T} — every edge bit is shared by 2 clusters — but for τ = Θ(1) it has Θ(Nτ²) marks of cost log(d/τ).

## 3. Where the cost sits (refined classes)
Classes (|T|,|A|) = (τN/d, αN/d), feasible iff α ≲ d e^{−τ} or τ ≲ d e^{−α}. Costs: T broadcasts e·2Nτ log(3d/τ);
A broadcasts e·2Nα log(3d/α); Theorem B e·O(h log(4Nd/h)), h = N(H(τ/d)+H(α/d)); enumeration h.
* both sparse (τ, α ≤ d^{1−δ}): Theorem B gives O(e N d^{−δ} log² d) = o(eN) [P given Thm B];
* τ ≤ 1/log d or α ≤ 1/log d: O(eN) by Lemma BC [P];
* remaining: τ ∈ [1/log d, O(1)] and α/d = Θ(1) (A a constant fraction of L∖N(T)), and the mirror image.
  exp2 §C (Thm B constant set to 1): the max–min is 7.2, 8.5, 10.2, 12.1, 18.0 × N/d at d = 16, 64, 256, 1024, 2^16,
  attained at τ ≈ 0.6–0.85, α/d ≈ 0.42–0.5: Theta((N/d) log d) [H: leading-order formulas].
In the hard classes A is free, so do not refine by α: **hard core = S_t = Σ_{|T|=t} 2^{−|N(T)|}, t = τN/d, τ = Θ(1)**,
= 2^{−dt} Σ_T 2^{c(T)}, c(T) = Σ_u (deg_T(u) − 1)^+ (a ferromagnetic lattice gas on the d-uniform hypergraph {N(u)}).
Enumeration costs h_T = (N/d) τ log(ed/τ); every network found costs e·d·h_T; target would be e·d·h_T/log d.

## 4. Cut ranks (what a rank argument can and cannot show)
**Lemma C [P].** For a cut L = L1 ∪ L2, the matrix F[A1,A2] = 2^{|R∖N(A1∪A2)|} equals ζ1 · C · ζ2^T with ζ_i the (invertible)
subset matrices and C[B1,B2] = #{T ⊆ R : L∖N(T) = B1 ∪ B2}. So rank F = rank C ≤ min_i #{N(T) ∩ L_i}.
(Proof: F = Σ_T [A1 ⊆ L1∖N(T)][A2 ⊆ L2∖N(T)].) exp1: rank F = rank C in all cases; N = 12, d = 4, balanced: 35 < 2^6.
* Balanced cuts: log rank ≤ max_τ min{N H(τ/d), (N/2) H(e^{−τ})}: 0.41N, 0.23N, 0.10N, 0.041N at d = 16, 64, 256, 1024 (exp2 §D).
* Single cluster (s vertices): rank 2^s (exp1; s ≤ 3). Hence for a faithful single network that holds A (the sub-sum
  over the settings mapping to A equals F(A)): U ≥ Σ_a log rank_a = N, and U ≥ 4·E_S log rank_S for a random half S of the
  nodes (each pair is cut with probability 1/2): at most 0.16N at d = 1024. For T-only networks with profiles: ≤ h_T.
* So cut ranks give U ≥ N at best — exactly the target — and say nothing about N log d. The whole gap between N and
  N log d is the difference between the cut rank of a node (rank of the row-wise Kronecker product of the Ψ_ab) and the
  sum of the pair ranks: one string of a cluster is needed, in different parts, by ~sd other clusters.

## 5. Obstructions for explicit classes of networks
**O1 (zero-pattern rank) [P].** Let M be any real matrix indexed by s-subsets X and s′-subsets Y of [m], s+s′ ≤ m,
with M[X,Y] ≠ 0 iff X ∩ Y = ∅. Then rank M ≥ C(s+s′, s). (Restrict to a ground set of s+s′ elements and to the
pairs (X, complement of X′): non-zero iff X = X′.) If nodes hold exactly (T_a, A_a) and the network is exact on every
(T, A), then Ψ_ab vanishes exactly on the conflicting pairs (put ∅ at the other nodes), so on a matching of m edges
between T-side and A-side with profile (s, s′) = (τm/d, Θ(m)): c_ab ≥ s log(e s′/s) = Θ(s log d). This strengthens Prop. Z
of G9c from the 0/1 conflict test to arbitrary weights/signs (Gottlieb gives C(m, min) for the 0/1 matrix; exp1 table).
When both sides are sparse the bound drops to 2s — but those classes are cheap anyway (§3).
**O2 (gathering) [P under the expansion hypothesis].** T-only networks in which node b learns the raw bits T|_{K_b}
(K_b ⊇ R_b) and every factor u is evaluated at a node with N(u) ⊆ K_b. If |N(U)| ≥ c·min(d|U|, N) for all U ⊆ L
(true whp with an absolute c [L]) and no node knows a constant fraction of R, then Σ_b |K_b ∖ R_b| ≥ (cd − 1)N, and under
the typical profile the cut of b carries ≥ log ∏_a C(k_ab, μ_ab) bits (identity submatrix), so
U ≥ 2(cd−1) N H(τ/d)(1−o(1)) = Omega(τ N log(d/τ)). [H] The same for trees of partial ORs: a tree for u has d leaves in d
distinct clusters, hence ≥ d−1 crossing messages, each a non-constant function of ≥ 1 bit of density ≥ τ/d, and two
different u share at most one neighbour in a girth-6 neighbourhood, so messages are not shared: (d−1) H(τ/d) per u.
exp2 §B: this lower bound equals enumeration at ed = 1 to within 2%, and Lemma BC is within a factor 2.1 of it.
**O3 (hierarchy) [P/H].** Q ≤ d super-clusters: the inter-part messages [N(u)∩T∩part ≠ ∅] total N·Q·H(τ/Q) = O(Nτ log(eQ/τ))
bits, O(N) for Q = O(1), but then the recursion runs at K = Q with e(Q) = Theta(1): cost Theta(N). Fan-in f trees of ORs
cost τ log(ed/(τ f^j)) per level: flat is optimal.
**Not an obstruction:** all of the above are about networks whose settings and messages are restrictions/ORs of T.
Appendix D beats the analogous bound for CNFs by cancellation (parity checks on a signed sum, where the "bad"
assignments contribute 0). In S_t every T contributes a non-zero term, and I found no signed expansion with
O(N) marks of O(1) bits, nor a 2-adic cut-off (the mark expansions 2^{(k−1)^+} = #even subsets have no 2-adic decay).

## 6. Q-hyper restated
Lemma 1 + §6 give [P]: weighted #CSP with k variables of domain 2^λ, each in ≤ r constraints of bounded arity:
exponent e·2(r−1)·kλ (+ o). Q-hyper ⇔ O(e·kλ) for r up to ≈ 1/e; the benchmark hard core is the case λ = 1 (entropy
H(τ/d) per variable), k = N, r = d, tables 1+[all zero]. Necessary: r ≪ eK (table size). Sanity check [H]: a
k′-clique embeds in a degree-r 2-CSP on k variables only for k′ ≲ sqrt(kr/log k), so O(e k) would contradict the k-clique
hypothesis only if r ≫ e² k log k — never in the range r ≤ k^{1/7} that Lemma 1's overhead gP/K allows. No conditional
obstruction from clique-type hypotheses; none found from hyperclique either.

## 7. Experiments (this directory)
* exp1_identities.py / exp1_out.txt: N = 8..12, d = 3, 4, brute force: one-sided form, clause form Σ 2^{|V|}, Fourier Z_{d+1},
  class-wise Holant identity, rank F = rank C for cuts |L1| = 1, 2, 3, N/2; Gottlieb ranks of DISJ and 2^{|X∩Y|} = C(m, min(s,s′)).
* exp2_costs.py / exp2_out.txt: exponents of every encoding at e = 1/d; refined-class max–min; balanced cut-rank estimate.

## 8. Self-review (adversarial)
* Nothing here improves a bound. Tasks 1–2 are answered negatively by counting; these counts are upper-bound formulas of
  specific encodings, not lower bounds on running time.
* "Covered by the paper" (weighted clauses) means: Lemma 1 and the §6 construction apply verbatim; I did not re-verify
  Lemma 1, and the paper's theorems are stated for CNFs only.
* §3 uses leading-order formulas with Theorem B's unknown constant set to 1 and the typical-set feasibility region; the
  location of the hard classes is [H], the three cheap regimes are [P] given G9c's Thm B / Lemma BC (reviewed in G10d).
* Lemma C and O1 are elementary and checked numerically only for N ≤ 12, m ≤ 10. "Single-cluster rank is full" is
  observed for s ≤ 3 (one exception 7 < 8 at N = 8), proved only via the induced-matching argument for s ≪ N/d.
* O1's reduction ("put ∅ at the other nodes") needs, under profiles with fixed |T_c|, a fixed conflict-free completion at the
  other nodes instead of ∅; such completions exist for typical profiles but this is not written out. O1 does not cover
  nodes with auxiliary settings (relays) — O2 covers only raw-bit relays.
* O2 depends on a vertex-expansion constant for random bipartite d-regular graphs that I did not compute, on typical
  profiles, and on a restrictive model (raw bits); the OR-tree extension is a sketch (dependence between messages of
  different u on the same pair is not treated rigorously).
* The balanced-cut bound uses #traces ≤ min over typical sizes; atypical T (clustered) are not counted rigorously.
* §6's "no conditional obstruction" is a failure to find one, not a proof that none exists.
* Open: Q1** (G9c) stands; sharpest form now: compute S_t = Σ_{|T|=τN/d} 2^{−|N(T)|}, τ = Θ(1), in exp2(O(N/d)).
