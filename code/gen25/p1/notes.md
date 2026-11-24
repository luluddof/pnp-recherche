# G25-P1: satisfiability-preserving (decision-only) sparsification. Idea G2-1

Labels: [P] proved here (complete proof below), [P-cited] taken from a cited note or preprint without re-checking,
[H] heuristic, [Q] open. Throughout, "piece" means a CNF F_i used for DECISION. Soundness means sat(F_i) ⊆ sat(F), or
more weakly F_i SAT ⇒ F SAT. Completeness means F SAT ⇒ some F_i SAT. A family that is sound and complete is called
sat-preserving. A family is "covering" when ∪ sat(F_i) = sat(F) (this is the counting/partition notion).

Summary of outcome.
* Props L and E exclude only covering sparsifiers, and only restricted kinds of them (§1). The brainstorm's check is correct, and a much stronger
  statement holds (Prop 0): for decision, a 1-piece sat-preserving family of density ≤ 1 always EXISTS, in every reasonable
  piece class. So no existence/counting obstruction of Prop L/E type can exist for decision. Unrestricted G2-1 is
  *equivalent* to the goal itself, via search-to-decision. Test (2) of the brainstorm ("Prop L analogue for restriction
  pieces on planted isolated instances") is therefore void as stated: the planted solution is itself a 1-piece family.
* Obstructions for decision can only be statements about restricted PROCEDURES. Proved here:
  - Prop C [P]: a "modular" decision sparsifier is automatically covering. Modular means it stays sat-preserving after any
    bounded-density formula R is conjoined, so it can be run on a dense sub-part and the rest added afterwards.
  - Lemma L' + Prop L' [P]: for restriction+subformula pieces, Prop L holds WITHOUT its component-respecting (product)
    caveat, and for φ_t ∧ M with M any set of monotone clauses.
  - Together: every modular restriction+subformula decision sparsifier needs 2^{n/(19D)} pieces on φ_t, so log log n stays.
* Positive side, Prop P [P]: the "all-false escape" is robust. Adding to Prop L's instance any monotone clauses of max degree
  ≤ t/log t leaves it decidable in 2^{O(n/√log n)}, by one pure-literal step plus the counter. Any hard instance for
  decision needs mixed-sign clauses, or monotone "hub" structure.
* No sat-preserving sparsifier at ε = 1/log n was found (construction: none; prior 6%, now lower). Also [P, elementary]:
  density-adaptive decision algorithms with time 2^{O(n log D / D)} (the best possible for branching on random formulas)
  combined with the counter still leave the log log n factor (§5).
* Density dependence of the preprints (§6). Claude: exp O(k m/√log m) = exp O(k D n/√log(Dn)), linear in D = m/n.
  Appendix D: exp O(e k n log Δ), max degree Δ ≤ (log n)^{4kC0}; with the project's Prop R, exp O(k n log D/√log n) for
  D = ⌈m/n⌉ and Δ ≤ n^{1/(40C0)}. OpenAI: exp O_k(s/√log s), s = n + m (linear in D) for m ≤ n log log n, and an
  n log log n-wire route otherwise.

-----------------------------------------------------------------------------------------------------------------------

## 1. What Props L and E exclude (read from code/gen6/g6c/loglog_analysis.md and code/gen7/g7d/q1_analysis.md)

Prop L (G6c §5) [P-cited, re-read].
Instance: φ_t is a disjoint union of N cliques Q_c of size t, n = Nt, with clauses (¬u ∨ ¬v) for u,v in the same clique.
|sat φ_t| = (t+1)^N.
Class excluded: COVERING sparsifiers (φ_t = union of the outputs' solution sets) whose outputs are CNFs ψ on the SAME
variables with sat ψ ⊆ sat φ_t, ≤ Dn clauses each, and COMPONENT-RESPECTING (ψ = ∧_c ψ_c with ψ_c on the variables of Q_c).
Conclusion: P ≥ ((t+1)/(2+√(2Dt)))^N, so P ≥ 2^{n/(19D)} at t = 32D. Hence P ≤ 2^{εn} forces D ≥ 1/(19ε).
Explicitly NOT covered (the note's own caveats): (i) non-product outputs; (ii) outputs with extension variables (the
sequential counter compresses "at most one of t" to O(t) clauses).

Prop E (G7d §2) [P-cited, re-read].
Instance: φ_G = ∧_{uv ∈ E(G)} (¬x_u ∨ ¬x_v) for a typical graph with nd/2 edges, 16 ≤ d ≤ √n.
Class excluded: SINGLE-piece (P = 1) projection-exact re-encodings: any CNF ψ(x,y), of any width k', with any number of
auxiliary variables, and π_x(sat ψ) = sat φ_G, needs ≥ nd/(8k') clauses (counting argument; G is determined by IS(G)).
Multi-piece version: only for decision trees of depth ≤ εn. The note says explicitly that with P = 2^{εn} pieces the
counting is vacuous.

Both are COUNT- (or projection-) preserving statements. Both arguments need many solutions in every piece: the line lemma
needs two solutions per pair {u,v}, and Prop E needs the full IS(G).

## 2. The check, and Prop 0 (existence is trivial for decision)

Check [P, trivial]. φ_t and φ_G contain only negative clauses, so 0^n satisfies them. The single piece "n unit clauses
¬x_1,...,¬x_n" is sound and complete, has density 1, and is found in polynomial time (every literal is pure). So neither
Prop L nor Prop E says anything about decision. The brainstorm claim is correct.

**Prop 0 [P].** For every CNF F on n variables there is a sat-preserving family with exactly one piece of the form
(partial assignment ρ) ∧ (subformula S ⊆ F|ρ), with |S| ≤ 1:
if F is SAT, take ρ = a full satisfying assignment and S = ∅;
if F is UNSAT, take ρ = 0^n and S = {the empty clause}. Some clause of F is falsified by 0^n, and its restriction is
empty. If the empty clause is not allowed, use the two clauses (x_1),(¬x_1).
*Proof.* Immediate from the definitions. ∎
Consequences.
(a) No information-theoretic obstruction (line lemma, counting, description length) can exist for decision in ANY piece
class that contains full assignments. In particular "restriction + subformula" pieces cannot be obstructed on planted
or isolated instances: the planted solution is a valid 1-piece family.
(b) Equivalence with the goal [P]. If 3-SAT decision runs in time T(n) (non-decreasing), then by self-reducibility
(fix x_1, ..., x_n one at a time; each restriction of a 3-CNF is a ≤3-CNF) a satisfying assignment is found with ≤ 2n calls.
So the 1-piece family of Prop 0 is computed in time O(n·T(n)) poly(n). Conversely, a sat-preserving sparsifier with
2^{O(n/√log n)} pieces of density O(1), running in time 2^{O(n/√log n)}, decides 3-SAT in 2^{O(n/√log n)} via the
counter (§6). So G2-1 without restrictions is a REFORMULATION of the goal (as G7d §0 found for Q1). It is not a weaker
sub-problem. The difficulty is entirely algorithmic: finding which solutions may be dropped.

## 3. What CAN be obstructed: modular sparsifiers, and Prop L without the product caveat

### 3.1 Prop C (modularity forces covering) [P]
Call a decision sparsifier *R-robust* if, for every CNF R on the same variables with at most n clauses, the family
{F_i ∧ R} is sat-preserving for F ∧ R (with soundness sat(F_i) ⊆ sat(F)). This is exactly what one needs in order to
sparsify a dense sub-formula A of F = A ∧ R separately and then add the sparse remainder R back. That is how the
sparsification lemma is normally used, and how Calabro's per-sunflower branching is local.
**Prop C.** An R-robust family is covering: ∪_i sat(F_i) = sat(F).
*Proof.* "⊆" is soundness. Given a ∈ sat(F), take R = the n unit clauses of a. Then F ∧ R is SAT, so some F_i ∧ R is SAT.
Its only possible solution is a, so a ∈ sat(F_i). ∎
So every robust (modular) decision sparsifier falls under the counting obstructions. Any decision-specific gain must
use information about the WHOLE formula (global autarkies, or an actual search for solutions). Local rules such as
"in this sunflower keep only the heart-true branch", decided from the sunflower and its neighbourhood alone, are
robust, and hence give no gain.

### 3.2 Lemma L' (subformula line lemma; no product structure needed) [P]
Let Φ = φ_t ∧ M, where M is ANY set of monotone (all-positive) clauses, ρ is any partial assignment, and S ⊆ Φ|ρ is a
subformula with sat(S) = sat(Φ|ρ) (on the free variables). For a clique Q_c that contains no vertex set to 1 by ρ, let
L_c = {u ∈ Q_c free : some s ∈ sat(Φ|ρ) has s_u = 1}. Then for every u ∈ L_c and every free v ∈ Q_c ∖ {u}, the clause
(¬u ∨ ¬v) belongs to S. Hence |S| ≥ Σ_c C(|L_c|, 2).
*Proof.* Take s ∈ sat(Φ|ρ) with s_u = 1. Then s_v = 0, and u is the only true vertex of Q_c in s. Let s' = s with v set to 1.
The clauses of Φ|ρ are of three kinds:
- restricted monotone clauses: these are still monotone, so s' satisfies them because s does;
- negative clauses not containing v: they are unchanged from s;
- negative clauses (¬v ∨ ¬x) with x ∈ Q_c free: s' violates such a clause only for x = u.
So s' violates exactly one clause of Φ|ρ, namely (¬u ∨ ¬v). Since S ≡ Φ|ρ, s' violates some clause of S. Since S ⊆ Φ|ρ,
that clause must be (¬u ∨ ¬v). ∎
Brute-force check: p1/check_Lprime.py (seeds 1 and 2), 2×300 random trials with N ∈ {2,3}, t ∈ {3,4}, random M of up to 4
monotone clauses of width 1–3, and random ρ of up to 3 variables. All 4908 tested pairs (u,v) are necessary clauses.
0 counterexamples.

### 3.3 Prop L' (covering restriction+subformula pieces on φ_t ∧ M) [P]
Let Φ = φ_t ∧ M as above. Suppose the pieces (ρ_i, S_i) satisfy S_i ⊆ Φ|ρ_i, |S_i| ≤ Dn, sat(ρ_i ∧ S_i) ⊆ sat(Φ), and
∪_i sat(ρ_i ∧ S_i) = sat(Φ). Then P ≥ |sat Φ| / (2 + √(2Dt))^N. For M = ∅: P ≥ ((t+1)/(2+√(2Dt)))^N ≥ 2^{n/(19D)} at t = 32D.
*Proof.* sat(ρ_i ∧ S_i) ⊆ sat(Φ) together with S_i ⊆ Φ|ρ_i forces S_i ≡ Φ|ρ_i, so Lemma L' applies.
Every solution in the piece has at most one true vertex per clique. That vertex lies in L_c, or it is a vertex set to 1
by ρ (then the clique has 1 option). So |sat(piece)| ≤ Π_c (1 + ℓ_c), where ℓ_c = |L_c| and Σ_c ℓ_c(ℓ_c − 1)/2 ≤ DNt.
Let L̄ be the mean of the ℓ_c. By Cauchy–Schwarz, N L̄² ≤ Σ ℓ_c² ≤ 2DNt + N L̄, so L̄ ≤ 1/2 + √(2Dt + 1/4) ≤ 1 + √(2Dt).
By AM–GM, Π(1 + ℓ_c) ≤ (1 + L̄)^N ≤ (2 + √(2Dt))^N. Covering then gives the bound. The numerical part is as in Prop L. ∎
This removes caveat (i) of Prop L for pieces that are restrictions plus subformulas (Calabro's outputs, with or without
the N-clauses of Appendix C, are of this kind up to added unit clauses and added sub-clauses; see the remark below).
Remark (limits of L'). The proof needs every clause of S other than (¬u ∨ ¬v) to be satisfied by s'. This holds for
clauses of Φ|ρ and for any added monotone clauses. It fails for arbitrary added clauses with negative literals, so
Calabro's added hearts (¬u) or (¬u ∨ ¬w) need a separate check: a unit (¬u) kills u from L_c, which only helps the
bound. A general "CNF on the same variables" version remains Prop L's conjecture [Q].

**Corollary [P].** A modular (R-robust) decision sparsifier whose pieces are restriction+subformula pieces needs
2^{n/(19D)} pieces on φ_t (Prop C + Prop L'). At P ≤ 2^{2en}, e = Θ(1/√log n), this means D = Ω(√log n), so the counter
pays log D = Ω(log log n). Within this class, the log log n factor cannot be removed.

## 4. The positive side: the all-false escape is robust (Prop P), and what a hard decision instance must look like

**Prop P [P].** Let Φ = φ_t ∧ M, with M monotone clauses of width ≤ 3 in which every vertex occurs at most Δ_M times.
There is a polynomial-time, decision-preserving reduction of Φ to a single formula Φ' with n' ≤ 3NΔ_M = 3nΔ_M/t free
variables and m' ≤ n't/2 + NΔ_M clauses. With the counter (Prop R), Φ is decided in time exp O(e·n·Δ_M·log(tΔ_M)/t)
(for t, Δ_M ≤ n^{δ}). In particular, it is decided in 2^{O(n/√log n)} whenever Δ_M log(tΔ_M) ≤ t.
*Proof.* Let V+ be the set of vertices occurring in M.
- If Φ has a solution with true set T, then |T| ≤ N, and every clause of M contains a vertex of T. Hence
  |M| ≤ NΔ_M and |V+| ≤ 3|M| ≤ 3NΔ_M. So if |V+| > 3NΔ_M, output UNSAT.
- Otherwise set every vertex outside V+ to 0 (these are pure negative literals).
  * Φ SAT ⇒ Φ' SAT: zero the true vertices outside V+; the negative clauses stay satisfied, and M is untouched.
  * Φ' SAT ⇒ Φ SAT: extend by zeros.
- Φ' has n' = |V+| variables, Σ_c C(|Q_c ∩ V+|, 2) ≤ n'(t−1)/2 negative clauses, and |M| ≤ NΔ_M monotone clauses.
- The function x log(A/x) is increasing for x ≤ A/e, so n' log(m'/n') = O((n' + NΔ_M) log(tΔ_M)) = O(NΔ_M log(tΔ_M)).
  Prop R then gives the time bound. ∎
The same argument works for any negative part whose solutions have a known bound α on the number of true variables (for
example a clique partition into N parts), with N replaced by α.

So "repairing" Prop L's instance against decision by adding monotone clauses needs vertices of M-degree ≳ t/log t. Each
solution's ≤ N true vertices must hit all of M, so M is star-like around a small set of hub vertices.
Example: 3-colouring encodings are exactly φ_G ∧ M (at-most-one per vertex triple, edge clauses, positive "some colour"
3-clauses). There t = 3, the density comes from the edges (Prop E-type), and the decision side has dominating-set
branching plus 2-SAT for graphs of large minimum degree [H/known-type, not re-proved].
Conclusion [H]: hard instances for restricted decision sparsifiers need many MIXED-sign clauses. Monotone/antimonotone
mixtures are absorbed by pure literals. Random dense mixed formulas are easy for decision (unsat: a random sub-formula of
density ≥ 5 is already unsat, a 1-piece subformula family; planted: polynomial-time recovery at large constant density,
Flaxman-type, cited from memory, not checked). Hard instances must therefore be structured mixtures, which no current
technique can obstruct without a model restriction (Prop 0).

Candidate procedural model, left open [Q]. Restriction trees with:
- variable branching (both children kept);
- forced-literal steps;
- autarky steps of width ≤ w (the only general, locally verifiable solution-dropping rule);
- leaves that are equivalent subformulas with ≤ Dn clauses.
Lemma L' gives the leaf cost, and a Prover–Delayer argument (Delayer answers "u = 1" with probability 1/(ℓ+1) on a clique
with ℓ live vertices) reproduces Prop L's bound when w = 0: the per-clique path probability is ≤ (a+1)/(t+1). What is
missing is control of width-w autarkies created by restrictions on an instance with no autarkies at the root. By Prop P,
such an instance needs hub-type monotone parts or mixed clauses. I did not find an instance for which I could prove this.

## 5. Constructions tried for (2): none works
- "Keep one branch per sunflower": by Prop C, a rule decided locally is covering, so Calabro's density bound and Prop L'
  apply. A rule that uses global information is a solution search (Prop 0(b)).
- Heart-true dominance (keep only H true when flipping a heart literal ℓ to true preserves satisfaction). This is valid
  exactly when ℓ is an autarky-type literal relative to the clauses containing ¬ℓ. It is not available in general, and on
  the dense k = 3 case (heart = 1 literal, 2-clause petals) the heart-false branch is the expensive dense 2-CNF part.
- Decision-only extension variables: replacing a biclique of negative 2-clauses A × B with one variable z, via
  (¬a ∨ z), (¬z ∨ ¬b), is projection-exact (π sat = sat), not count-preserving. This is more than the counting route
  allows (it compresses Prop L's cliques AND bicliques). But Prop E already excludes projection-exact single-piece
  re-encodings of random graphs (no large bicliques), so it does not handle Prop E-type density with mixed-sign repair.
- Density-adaptive decision [P, elementary]. Suppose instances of density D could be decided in 2^{O(n log D/D)} (the
  natural branching rate; DPLL achieves 2^{O(n/D)} on random 3-CNF). Combined with the counter's exp O(n log D/√log n),
  the worst case is min(log D/D, log D/√log n) at D ≈ √log n, which equals Θ(log log n/√log n). So log log n stays.
  Removing it needs every density D ≥ C (C a constant) to be handled within 2^{O(n/√log n)}, i.e. a sparsification to
  density O(1), or a counter that does not pay log D.

## 6. Consequence if G2-1 were positive, and the exact density dependence

Claude preprint (claude_3sat.txt):
- Theorem A(c): exp O(k m/(log m)^{1/2}). At m = Dn this is exp O(k D n/√log(Dn)), LINEAR in D (footnote 1:
  "already sparse ⇒ 2^{O(kn/√log n)}").
- Theorem D(a) (Appendix D, Steps 1–4 "proved in full", Step 5 sparsification a sketch): exp O(e k n log Δ), max degree
  Δ ≤ (log n)^{4kC0}, e = Θ(1/√log n). The paper's remark (line ~2440) allows Δ ≤ n^{1/(40C0)}.
- With G6c's Prop R (project, [P-cited]): exp O(e k n log D), D = max{2, ⌈m/n⌉}, Δ ≤ n^{1/(40C0)}.
- After Calabro at ε ∈ [e, 2e]: log Δ ≤ 3k log(8k/e), giving exp O(k² n log log n/√log n) (Theorem D(b)).
- Theorem A(d) (decision): exp O(n/(log n)^{1/20}), from M(ε) = O(ε^{-9}) clauses per piece.
- The paper itself notes (Appendix D remarks): "only log Δ enters" and "we do not know whether the factor log log n can
  be removed".
OpenAI preprint (g14b/dl/paper.txt):
- Theorem 1.1: exp O_k(n log log n/√log n) in general, and exp O_{k,D}(n/√log n) for m ≤ Dn.
- Appendix D.2: exp O_k(s/√log s), s = n + m, i.e. linear in D; used when m ≤ n log log n.
- §2.9 states the same open target: "If sparsification gave B = O_k(1) at σ = Θ(ρ(n)/n), we would already have ...".

Consequence [conditional on the preprints' Theorems A(c)/D(a), resp. Theorem 1.1, and on Prop R for Claude's route].
Suppose a sat-preserving sparsifier runs in time 2^{O(n/√log n)} and outputs ≤ 2^{O(n/√log n)} pieces, each a 3-CNF on
≤ n variables with ≤ C·n clauses (C an absolute constant) and max degree ≤ n^{δ}. Then 3-SAT DECISION (not counting) is in
2^{O(n/√log n)} at every density: time = sparsifier + Σ_i exp O(k n log C/√log n).
The weaker hypothesis "pieces of density (log n)^{o(1)}" gives exp o(n log log n/√log n) through Prop R. With the
OpenAI bounded-density branch, it is enough to have pieces of density ≤ C for an absolute C.
(Max degree ≤ n^{δ} costs nothing extra: a piece with ≤ Cn clauses has at most 3Cn^{1−δ} variables of degree > n^{δ}, and
branching on them costs 2^{O(n^{1−δ})}.)
By Prop 0(b), any such sparsifier is as hard to find as the decision algorithm itself. The realistic forms are
restricted, and the modular (robust) ones are excluded (Prop C + Prop L').

## 7. Self-review
- Prop 0, Prop C: trivial but correct, and they change the framing of G2-1: decision sparsification has no existence
  obstruction; robust sparsification is covering.
- Lemma L': checked by hand and by brute force (small sizes only). Prop L' arithmetic re-done (Cauchy–Schwarz + AM–GM).
  It applies only to subformula pieces; Calabro's added sub-clauses are monotone only if the hearts are positive, so the
  "Calabro is covered" remark is partial (units ¬u only help the bound; added negative 2-clauses ¬u∨¬w inside a clique are
  subsets of φ_t anyway; mixed added clauses are not covered).
- Prop P: rigorous; the cost bound uses Prop R (project) and its max-degree condition.
- §4's [H] statements about 3-COL and planted instances are from memory and not checked.
- Not done: a lower bound in the autarky-tree model (§4, [Q]); a construction (none).
