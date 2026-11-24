# G25-B3 — Brainstorm: ALGORITHMS (routes that avoid the crossing trade-off; lower rungs)

Read: BRIEF.md; PISTES §0, §4, §5; code/spec/requirements.md; g21c, g18c notes; g6c/loglog_analysis.md, g7d/q1_analysis.md
(the loglog question), g9b (decide vs count), g20e (instance adaptivity). Web (Nov 2026): arXiv 2610.41622 (lowered crossover of
Claude's algorithm; says it refutes the Sparse Max-3-SAT Hypothesis), 2611.03880 (Fischer, "The Hyperclique Hypothesis is False",
algebraic), 2611.06813 (Chan-Itakura-Ma-Yu, 4-hyperclique in n^{4-eps}, so Max-3-SAT in 2^{(1-eps/4)n}). arxiv.org is NOT reachable
from this session (DNS blocked); only search snippets were seen. Tags: [P] proved, [H] heuristic, [O] open.

Ideas I looked at and dropped (already closed, no escape found):
 - Fast rectangular MM / Williams polynomial method on the leaves: saves 2^{s/polylog s} on a 2^s leaf, i.e. a lower-order term (G9d).
 - Decide over GF(2) with Valiant-Vazirani isolation: closed by G9b (best GF(2) rho = 0.811).
 - Instance-adaptive pieces for special classes (aim alpha > 1/2): closed by G20e Prop R (the restricted crossing is still a full swap).
 - Noisy wires (approximate counting is enough for decision): the SAT/UNSAT gap forces noise < 2^{-n}, so the noisy identity has full rank. [H, quick count]
 - Grover / quantum: constant in the exponent only.

====================================================================================================
GROUP 1 — other routes toward step 2 or beyond
====================================================================================================

## G1-1  Orientation designs: more directions per gadget (the "D <= log M" resource)
(a) In both algorithms the per-level gain is D^{-1/2} (the LAW of G18c), and D = O(log M) only because each gadget must fill all 2^D
    orientation classes (a packet of L = 2^a slots gives a = log L directions). Proposal: let the gadget's orientation vectors form an
    orthogonal array / code of strength t (L slots, D >> log L directions). If Hermite interpolation only needs the low-order
    marginals (t-wise) to be exact, D can be (log M)^2 and D^{-1/2} gives e = 1/log M = step 2, with NO change to the crossing trade-off.
(b) The jets' price is a graded dimension that depends on how many classes each node has to resolve; a strength-t design might resolve
    the identity using only t-wise class information.
(c) Most likely the identity on the gadget is a full crossing in every direction, so all 2^D classes are needed (a counting bound
    like G9a's rank count). The composition rule of G18c (prod D_j <= #coordinates) would also have to be bypassed: here the
    directions are not coordinates of a product labelling.
(d) Not tried as such. G18h lists "identity with D directions" as a gain object, and G20a tested D = 3 per level (no gain), but
    neither looked at the resource side (directions per gadget size).
(e) Test (< 1 h): reuse g18c/bkernel.py. Build a packet whose slots carry a strength-2 OA of length D = 3 on L = 4 slots, instead of
    the full cube (8 slots). Ask whether the packet identity is a Hermite combination of direction-local kernels (exact rank
    computation mod a 61-bit prime, as in out_B_a1.txt). Any exact identity counts as a "yes".
(f) 2% (value: step 2 if it worked).

## G1-2  Refresh the labels between levels (more than loglog K levels)
(a) Claude's recursion has l = log log K levels because each level turns K = g^2 into g, using up half of the coordinates.
    Proposal: after a level, re-embed the reduced network into a fresh g^2-clique by SPLITTING each vertex into g pieces
    joined by internal wires. Wire count is cheap (K^2 << M). If re-embedding multiplies the wire BITS by c with c*rho < 1,
    the number of levels is no longer capped. The harmonic law of G9d caps this at step 2 (m/log m), so step 2 is the ceiling.
(b) Vertex tensors are arbitrary; the only cost is the cut-width of each vertex's internal tensor once split.
(c) After a level the wires are thick (alphabet = weight window W_E), and splitting a vertex with boundary B probably needs internal
    cuts of about B bits. That gives c >= 2 > 1/rho = sqrt 2: dead. G18d ("depth only buys sqrt loglog; the root dominates") is
    probably this obstruction in another form.
(d) Close to G18d and G4b/G6b (recursion geometries), but the explicit "split and re-embed" count is not written down anywhere.
(e) Test: from claude_3sat.txt sec. 6 (cost of clusters), compute bits per vertex after one level (rho*M/g) and the minimum cut-width
    of a split vertex in the App. B piece tensors (rank of A_s flattenings: out_A_profile.txt already holds the profile). Compare
    c with sqrt 2.
(f) 2%.

## G1-3  Import the Nov-2026 hyperclique algorithms as gadgets
(a) Two new algebraic results (Fischer, 2611.03880; Chan-Itakura-Ma-Yu, 2611.06813) beat a "hyperclique" hypothesis that is
    3-uniform, i.e. three directions at once. If their core is a low-rank or low-border-rank decomposition of a 3-uniform tensor
    (rather than plain fast MM), it may be exactly the "D = 3 direction identity" that G18h and G20a say would beat D^{-1/2}.
(b) It is a new, external algebraic object that appeared after G20a, aimed at the same 3-way structure.
(c) It is most likely a polynomial saving n^{-eps} for fixed k (fast MM type). By G21c's Theorem R, such savings give
    constants only in our setting.
(d) Not tried (the papers are 2-3 weeks old).
(e) Test (< 1 h, needs arxiv.org allowed): read both abstracts and key lemmas. Classify the core identity as "fast MM / rank
    of a fixed tensor" (closed by Theorem R) or "an identity whose terms are cheap in one of 3 directions" (then test it in
    the g20a script).
(f) 4% useful; a clean "no" is about 60% likely and cheap.

## G1-4  Dense Max-3-SAT and weighted counting in the n-scale
(a) Polynomial-valued tensors (one formal z per clause) give the full distribution of satisfied clauses at the same m-scale cost,
    so Max-3-SAT is in 2^{O(m/sqrt log m)} (consistent with 2610.41622). Sparsification does NOT preserve Max-SAT, so dense exact
    Max-3-SAT in n stands at 2^{(1-eps/4)n} (2611.06813). Proposal: sampling/regularity decomposition of dense Max-3-CSP
    (Frieze-Kannan weak regularity) into a low-complexity part plus a sparse residual handled by the 2026 network.
(b) Dense Max-CSPs have low "cut-norm complexity". An exact solver needs the residual handled exactly, and that residual is
    sparse-ish.
(c) The weak-regularity error is eps*n^3 clauses, and exactness needs eps*n^3 < 1, so it fails as stated. A Williams-type hardness
    (Max-3-SAT -> 4-hyperclique) does not forbid 2^{o(n)}, but nothing points to it either.
(d) Not tried in the repo.
(e) Test: write the 1-page reduction "Max-3-SAT with m clauses in 2^{O(m/sqrt log m)}" carefully (degree in z = m, so m+1
    evaluation points), and check the n-scale ceiling: at which density does it beat 2^{(1-eps/4)n}? (Answer: m <= n sqrt(log n).)
(f) 2% for 2^{o(n)} dense; 30% for a clean "Max-3-SAT at density <= sqrt log n" statement (low value, probably folklore after
    2610.41622).

====================================================================================================
GROUP 2 — lower rungs (realistically provable)
====================================================================================================

## G2-1  SATISFIABILITY-preserving sparsification at eps = 1/sqrt(log n)  (decision only)  *** main new idea ***
(a) The loglog factor in n comes from sparsifying with eps ~ 1/sqrt log n, which forces density (1/eps)^{O(1)}. The only proved
    obstruction (G6c Prop L, the line lemma) is for COUNT-preserving sparsifiers that must cover every solution. For DECISION we
    only need pieces F_i with sat(F_i) ⊆ sat(F) and (F SAT => some F_i SAT). The pieces may lose solutions. Each piece is then
    decided by the 2026 counter (count > 0).
(b) [P, trivial] Prop L's hard instance (disjoint "at most one" cliques) and Prop E's (independent sets of random graphs) are both
    satisfied by all-false, so one piece of density 0 suffices. Both proved obstructions disappear. Natural pieces:
    partial assignments ρ_i plus F|ρ_i (fixing variables only strengthens F), Calabro-style branching that keeps only the
    "heart true" side when the petals are implied, and dropping clauses implied in bounded resolution width.
(c) Hard instances must be dense, satisfiable and solution-sparse with no local structure, e.g. isolated solutions (CIKP 2006
    isolation costs only eps*n). For unique-SAT, "keeping a solution" equals "keeping all solutions", but the line-lemma
    argument needs many solutions, so a new obstruction would be needed. An unrestricted construction is "equivalent to the
    goal" in the trivial G7d sense, so only restricted classes can be ruled out.
(d) Not tried: G6c, G7d and the PISTES item "sparsification with extension variables" are all count-preserving.
(e) Test (< 1 h, paper): (1) write the reduction "decision via sat-preserving pieces + Prop R". (2) Try to prove a Prop L analogue
    for restriction-based pieces on "planted dense 3-CNF with an isolated solution" (dense random 3-XOR-like constraints
    encoded in 3-CNF, density polylog n). (3) If that fails, look for a branching rule: high-degree variables whose sunflower
    petals become implied.
(f) Construction 6%. Clean obstruction for restriction-plus-subformula pieces about 25%. Value: decision in 2^{O(n/sqrt log n)}
    at every density, and a provable decision/counting split in the n-scale (with Prop L as the counting side).

## G2-2  Cluster-quotient bound (generalise G21c's Theorem T1)
(a) T1 handles binary #CSP with cost Delta*v*log s / sqrt log v, charging log s per edge. Refinement: partition the variables
    into parts of size <= r, treat each part as a node with an arbitrary tensor, and charge each edge only the bits of the
    shared boundary. Target statement: #3-SAT in 2^{O(beta_r / sqrt log(n/r)) + O(r)} poly(n), where beta_r = number of variable
    occurrences crossing parts in the best r-partition.
(b) Claude's Definition 2 already allows arbitrary node tensors with node-setting bits P_a. It is the same calculation as T1
    with c_{a,b} = shared bits.
(c) Gap: is the gP/K term ~ n*r/K harmless? Also the poly(lambda+B)^{CK^2} factor. This is mostly a corollary, so modest value.
    It removes loglog only for formulas whose dense parts are local (all of Prop L's instances); planar instances already have
    2^{O(sqrt n)}.
(d) Partly done: G21c T1 [C, S, proof sketch only].
(e) Test: redo T1's accounting with shared-bit edges. Evaluate beta_r on the G7d benchmark family and on "expander of
    polylog-size clusters" (the class where it gives 2^{O(n/(log n)^{1.5})}-type bounds).
(f) 50% (statement plus proof at sketch level); value modest (a paper-ready parameterized refinement).

## G2-3  Sub-exponential in treewidth / pathwidth: 3-SAT in 2^{O(pw/sqrt log pw)} poly(n)
(a) A natural rung that ETH used to forbid. Iwata-Yoshida (ESA 2015) show that pathwidth-parameterized SAT, 3-SAT, Max-2-SAT and IS
    are EPNL-complete and equivalent for c^{pw}. A 2^{o(pw)} result for one transfers to all of them.
(b) The 2026 algorithms beat treewidth on expanders (T1). Pathwidth cases need a "streaming" version that compresses the
    2^{pw}-dimensional transfer vector.
(c) The 2026 cost is linear in TOTAL wire bits (about n >> pw), so it does not apply. A 2^{o(pw)} poly(n) algorithm is close to
    simulating k-bit nondeterministic log-space machines in time 2^{o(k)}, which is much stronger than 2^{o(n)}.
(d) Not tried. Search found no statement linking 2^{o(pw)} to ETH-false.
(e) Test: check whether ETH-false implies pw-3-SAT in 2^{o(pw)} poly through any known reduction (for example Iwata-Yoshida's
    EPNL machine model). If not, record it as an open target and give a 1-page reason why exp(e*U) contraction cannot do it.
(f) 3%; value high.

## G2-4  Better k-dependence: k-SAT in 2^{O(k m / sqrt log m)} and n-scale with log k instead of k log(k/eps)
(a) The n-bound carries k^2 log log n (Calabro density (k/eps)^{3k}, so log D = O(k log(k/eps))). Combined with G2-1 (pieces
    that need not cover all solutions) or with the k-uniform sunflower bounds of Alweiss-Lovett-Wu-Zhang, the regime where the
    2026 bound beats PPSZ's 2^{(1-Theta(1/k))n} grows from k ~ (log n)^{1/4} to k ~ (log n)^{1/2-o(1)}.
(b) ALWZ improves the k-dependence of sunflowers, and G21c already notes it improves c(eps).
(c) It cannot touch the eps-dependence (Prop L). The value is a rung "in k", not in n.
(d) Mentioned in G21c §4 (sunflowers), not worked out.
(e) Test: redo G6c's exponent with an ALWZ-based sparsification density and find the largest k(n) with sub-exponential time.
(f) 25%; value low-moderate.

## G2-5  Dense-leaf dichotomy for COUNTING (remove loglog for #3-SAT)
(a) Run Calabro until density D. A leaf that cannot be sparsified further either (i) has a local cluster structure, so G2-2
    absorbs the dense parts into nodes, or (ii) is pseudo-random-dense, so its solution space has entropy about n*polylog(D)/D
    (as for the IS formulas of random graphs, G9c), which a solution-entropy-sensitive counter handles. A regularity-type
    dichotomy would give #3-SAT in 2^{O(n/sqrt log n)} at every density.
(b) Prop L's instances are type (i), Prop E's are type (ii), and G9c already counts type (ii) families in 2^{O(en)}.
(c) Mesoscopic mixtures (clusters of size (log n)^c glued pseudo-randomly) may be neither type. G7d proves single-piece
    re-encodings fail, so pieces must do real work.
(d) Pieces exist (G6c, G7d, G9c); the dichotomy itself has not been tried.
(e) Test: build the mixture family (random dense 2-CNF implication structure on clusters of size t, glued by a sparse random
    3-CNF) and compute both costs, min(beta_r, solution entropy), for n <= 40 by exact counting. If one family keeps both at
    Theta(n log t), the dichotomy is dead.
(f) 5%; value: step 1 exactly in the n-scale for counting.

## G2-6  Explicit constants via PPSZ-style isolation and pre-forcing
(a) Under isolation (CIKP), PPSZ forces a 0.61 fraction of the variables by bounded-width implication. Count only on the free
    0.386n variables, with the forced ones as small decision trees (width k' = O(1)), so the constant becomes c*0.386*poly(k').
(b) Only constants matter for physical sizes (G20b, G20e: at n < 2^{10^5}, PPSZ beats everything).
(c) Forced variables are depth-d functions, so the width blows up as 3^d. The 2026 cost grows like k^2, probably eating the
    0.386. A constant only.
(d) No (G9b mentions PPSZ, G20e compares numbers).
(e) Test: compute the k-dependence of Prop R's cost against the PPSZ forcing depth needed for a 50% forced fraction.
(f) 15%; value low.

====================================================================================================
TOP 3  (probability x value / cost)
====================================================================================================
1. G2-1 sat-preserving sparsification (decision). Both known obstructions vanish for decision. Either outcome (a construction,
   or a Prop L analogue for restriction-based pieces) is a publishable statement about decision vs counting. First test:
   write the reduction, then attack the "planted dense, isolated solution" family with restriction-based pieces.
2. G2-2 cluster-quotient bound. High probability, cheap, extends T1. First test: redo T1's accounting with shared-bit
   edges; evaluate beta_r on G7d's family.
3. G1-3 Nov-2026 hyperclique algorithms as 3-direction gadgets. Very cheap literature check with a large upside; the default
   outcome is a quick, clean closure (Theorem R). Needs arxiv.org allowed.
