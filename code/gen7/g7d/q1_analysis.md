# G7d: Q1 (count-preserving sparsification to bounded average density) — analysis

Status: **no proof of 2^{O(k^c n/√log n)}.** What is here: a meta-observation that Q1 without restrictions is
*equivalent* to the goal; a monotonicity lemma that rules out every transformation which only adds variables
and clauses (equality chains / splitting / padding / Tseitin), with exact bounds for route (b); a counting
obstruction for route (a) "sparsify, then re-encode each piece with extension variables" (single piece per
input, any encoding); a note on why an *explicit* obstruction for the extension route would be a superlinear
circuit lower bound; route (c) (inclusion–exclusion / signed sums) bounded; a concrete benchmark family.
Labels: [P] proved here, [P*] proved modulo a cited fact not re-checked, [H] heuristic, [Q] open.

Measure used throughout ("Prop-R cost"): for a k-CNF with n' variables and m' clauses, the Appendix-D bound is
exp2(O(e k n' log D')), D' = max{2, ⌈m'/n'⌉}, e = Θ(1/√log n') (Prop. R of G6c; Thm D(a) of the paper), valid for
k ≤ (log n')^{1/4} and max degree ≤ n'^{δ}. The m-bound is exp2(O(k m'/√log m')) (Thm A(c), α' = 1/2).
All statements about "cost" are about these *upper bounds*, not about the true running time of the algorithm.

## 0. Meta-observation [P]: unrestricted Q1 is the goal itself

If #k-SAT on n variables is computable in time 2^{O(n/√log n)}, the reduction "compute #F, output the O(n)-clause
comparator CNF on n variables whose models are {x : x <_lex N}, N = #F" (or a signed sum of ≤ n such cubes) is a
Q1 reduction with one output of density O(1). Conversely Q1 + Prop. R gives the goal. So Q1 is a *reformulation*,
not a weaker sub-problem: progress requires a reduction of a specific structural kind, and obstructions can only
be proved for restricted classes of reductions. Everything below is of that type.

## 1. Lemma S (adding variables/clauses never helps) [P]

**Lemma S.** Let F have n variables (all occurring) and m clauses, and let F' (any width k' ≤ (log n')^{1/4})
have n' ≥ n variables and m' ≥ m clauses, with n' ≤ n^{O(1)}. Then
n' log max{2, m'/n'} ≥ c · n log max{2, m/n} for an absolute c > 0; hence the Prop-R cost of F' is at least
a constant fraction of that of F (e' = Θ(e) because log n' = Θ(log n)).

*Proof.* If m ≤ 2n the right side is ≤ 2n ≤ 2n'. Otherwise put f(x) = x log(m/x), increasing on (0, m/e].
If n' ≤ m/e: n' log(m'/n') ≥ f(n') ≥ f(n) = n log(m/n). If n' > m/e: left side ≥ n' > m/e, and
m ≥ n log(m/n) (since y ≥ log y), so left side ≥ (1/e) n log(m/n). ∎

Consequences (all count-preserving, all of the form n' ≥ n, m' ≥ m):
* **Route (b), full splitting.** Replace each occurrence of a variable by a fresh copy and link the copies of
  each variable by an equality chain (2 clauses per link). n' = Σ deg ≤ km, m' ≤ m + 2km, max degree 3,
  m'/n' ≤ 3. Every original model extends uniquely. Prop-R cost: e' k' n' log 3 = Θ(e k m) — i.e. the same as
  the m-bound exp2(O(k² m/√log m)) (k factor from n' ≤ km). In terms of n, after Calabro with parameter ε
  (m ≤ 2(8k/ε)^{3k} n): exponent εn + C k² (8k/ε)^{3k} n/√log n; optimal ε ≍ (k^{3k+2}/√log n)^{1/(3k+1)}
  gives exponent n·poly(k)^{O(1)}·(log n)^{-1/(6k+2)} — for k = 3, (log n)^{-1/20}: much **worse** than the
  paper's (log log n)/√log n. Splitting trades n log D for n·D.
* **Partial splitting** (split a variable of degree d into ⌈d/D0⌉ copies): n' ≤ n + km/D0, density O(D0 + k).
  Cost e k (n + km/D0) log D0, minimized at D0 ≍ km/n, giving e k n log(km/n): unchanged (Lemma S).
* **Padding** with fresh variables forced by unit clauses, **Tseitin re-encoding** of existing clauses, any
  encoding that keeps the original clauses: covered by Lemma S.
So extension variables can only help if they let the output have **fewer clauses than F** (m' = O_k(n)), i.e.
genuine compression of the constraint set. This is what the clique example does ("at most one of t": t(t−1)/2
clauses → O(t) clauses via a sequential counter).

## 2. Route (a): sparsify, then compress each piece with extension variables

Required shape: P ≤ 2^{εn} pieces, ε ≈ 1/√log n, each piece ψ_i(x, y) with n'_i = O(n), m'_i = 2^{O(k)} n'_i,
width O(k), with sat F = ⊎_i π_x(sat ψ_i) (or a signed version).

**Prop E (single-piece counting obstruction) [P].** Let 16 ≤ d ≤ √n and let G range over graphs on [n] with
exactly nd/2 edges; φ_G = ∧_{uv∈E}(¬x_u ∨ ¬x_v). For all but a 2^{-Ω(nd log n)} fraction of G, every CNF ψ(x, y)
of width ≤ k' (any number of auxiliary variables, unique extension not even required) with π_x(sat ψ) = sat φ_G
has at least nd/(8k') clauses (for n ≥ n0(k')).
*Proof.* sat φ_G = IS(G) determines G (uv ∈ E iff the indicator of {u,v} is not a model). The number of graphs
is C(C(n,2), nd/2) ≥ (n/(2d))^{nd/2}·... ≥ 2^{(nd/2)(log n − log d − 1)} ≥ 2^{(nd/4) log n − nd/2} for d ≤ √n.
A CNF with s clauses of width ≤ k' uses ≤ k's variables; up to renaming of auxiliaries it is described by
≤ k's·log(2(n + k's)+1) + O(s) bits. If s ≤ nd/(8k') this is ≤ (nd/8)(log n + O(log(dk'))) bits, fewer than
(nd/4) log n − nd/2 − Ω(nd log n)... for n large; so only a 2^{-Ω(nd log n)} fraction of G can be represented. ∎
(I used s ≤ n² w.l.o.g.; constants are not optimized.)

*Consequence for the measure.* Any such single piece has n' ≥ n (projection onto all of x, all x_u occur when
d ≥ 1) and m' ≥ nd/(8k'); by the argument of Lemma S with m replaced by nd/(8k'):
n' log max{2, m'/n'} ≥ c·n log(d/(8k')). So for typical G with d = (log n)^{c0} (an instance of density d, already
counted by Prop. R in exp2(O(e n log d))), **no count-preserving re-encoding into one formula, with or without
extension variables, lowers the Prop-R cost below e·n·Ω(log log n).** In particular "Calabro + per-piece
compression" fails on these instances provided Calabro's algorithm returns φ_G essentially unchanged when all
degrees are below its sunflower threshold (θ ≈ (k/ε)^{Θ(k)} ≫ (log n)^{c0} for small c0) — [P*]: this behaviour
of [Cal09] is plausible from its definition (no sunflower with ≥ θ petals ⇒ no branching) but not re-checked.

**What Prop E does not cover [Q].** With P > 1 pieces the counting gives only P·s ≳ nd/k' (the collection of
pieces is a description of G), which is vacuous at P = 2^{εn}. A positive answer to Q1 on these instances
would need the 2^{εn} pieces themselves to "absorb" the edge information.

**Restricted multi-piece version [P].** If the pieces are the leaves of an adaptive decision tree on the
x-variables of depth ≤ εn (each leaf = subcube, piece = any CNF with π = sat φ_G on that subcube), then the
all-zero leaf alone must encode G[V∖Q] with |Q| ≤ εn (zeros kill nothing). Restrict to graphs of max degree
≤ 2d: the edges touching Q number ≤ 2εnd, so (Q, the leaf's CNF, the ≤ 2εnd edges at Q) determine G; these
cost ≤ n + k's log(3n+3k's) + 2εnd·2 log n bits, so s ≥ (1 − O(ε)) nd/(8k') for typical G [P, outline-level
counting as in Prop E]. Depth-bounded trees are a weak class (degree-branching trees have depth ≈ n along the
0-path), so this is recorded only for completeness.

**Heuristic for variable-branching trees [H].** For a random d-regular G, branching (v ∈ I removes d+1 vertices,
v ∉ I removes 1) until the residual graph has ≤ Dn edges needs residual size ≲ n√(2D/d); the recurrence
L(n) ≥ Σ_i L(n − i − d) over the 0-path gives L ≈ 2^{c n} with 2^{cd} ≈ 1/c, i.e. log L ≈ n log d/d·(1 − √(2D/d)).
With d ≈ (1/ε) log(1/ε) this forces D = Ω(d) = Ω(1/ε): the same D ≳ 1/ε as G6c's Prop. L, now for a random
(not clique-structured) instance where sequential-counter tricks have nothing to compress (Prop E).

**Why an obstruction for the extension route on explicit instances is out of reach [P].** A fan-in-2 circuit of
size s for F gives (Tseitin, auxiliary = gate values, output forced true) a width-3 CNF with s auxiliary variables,
≤ 4s+1 clauses and unique extension. So a statement "explicit F_n (m = n·polylog n clauses) has no count-preserving
single-piece encoding with O(n) clauses" implies an ω(n) circuit lower bound for an explicit function computable by
n·polylog n-size CNFs — beyond current techniques (best explicit bounds ≈ 3.1n). The obstruction of Prop E is
therefore necessarily non-explicit (counting), and the multi-piece case is open.

## 3. Route (c): inclusion–exclusion and signed sums [P, narrow]

* #(A ∧ C) = Σ_{S⊆C} (−1)^{|S|} #(A ∧ "all clauses of S false"): each term fixes the literals of S, so it is a
  restriction of A. With 2^{εn} terms this removes at most εn clauses, while a density-D(ε) piece has
  (D(ε) − 2^{O(k)}) n ≫ εn excess clauses. Mixing with branching gives the Calabro tree again (sign-free).
* #(A ∧ B) = #A − #(A ∧ ¬B): ¬B for a CNF B is a DNF; a count-preserving CNF encoding of ¬B must pick a
  canonical violated clause ("first violated"), which needs a chain over all |B| clauses: m' ≥ |B|, no gain.
* Weighted factors (the algorithm accepts integer tables, W may be negative): a product of local weighted factors
  replacing the dense clauses still has one factor per clause unless the clause set is compressible; Prop E
  applies verbatim to any representation as a product of width-k' factors with auxiliary variables summed out
  *with 0/1 weights*; for general integer weights the counting argument applies to the description length as
  well (weights of ≤ poly(n) bits), [P*] modulo writing out that the weighted sum determines IS(G) — it does,
  since the sum equals the 0/1 indicator. So signed/weighted sums do not escape Prop E for P = 1.

## 4. Benchmark family (where all known routes keep log log n) [P for the listed bounds, Q for the question]

Random d-regular graphs, d = Θ(√log n), as 2-CNF (all-zero is a model, so no unsatisfiability shortcut):
1. Prop. R direct: exp2(O(e n log d)) = exp2(O(n log log n/√log n)).
2. Enumeration: #IS ≥ 2^{α(G)} and α(G) = (2 − o(1)) n ln d/d whp (Frieze–Łuczak), so ≥ 2^{Ω(n log log n/√log n)}.
3. Degree branching to max degree D: 2^{Θ(n log D / D)} leaves; with D = O(1) this is 2^{Θ(n)}; with D ≈ d it is (1).
4. Single-piece re-encoding: Prop E.
**[Q1′]** Count independent sets in n-vertex graphs of average degree d ≈ √log n in time 2^{O(n/d)}.
A yes is necessary for the Appendix-D route to give 2^{O(n/√log n)} without new ideas, and Q1 restricted to
2-CNF IS instances is exactly this. Note random k-CNFs are *not* hard instances for Q1: above density 2^k ln 2 they
are unsatisfiable whp, so the hard instances must be structured (monotone 2-CNF parts like the above).

## 5. Self-review (adversarial)

* Lemma S: elementary; the only subtlety is e' vs e (fine for n' ≤ poly(n)). It concerns the *bound*, not the
  algorithm: a cleverer analysis of the network could, in principle, charge copies of one variable once
  (copies in one node). That just reproduces the original instance, so no gain either way.
* Prop E: counting is rigorous in outline; constants (nd/(8k')) not optimized; requires d ≤ √n. It is a P = 1
  statement; it does **not** refute Q1 (which allows 2^{εn} pieces). The application to "Calabro + compress"
  relies on the unverified [Cal09] behaviour on low-degree 2-CNFs (marked P*).
* The decision-tree heuristic (§2) is a recurrence estimate, not a proof; random-graph structure conditional
  on the tree is not controlled.
* §0 shows any full answer to Q1 is as hard as the goal; I found no construction and no general obstruction.
* Nothing here improves or refutes the paper's theorem; log log n remains.
