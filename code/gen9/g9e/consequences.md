# G9e (track E) — New algorithmic consequences of the claimed 2^{o(n)} #k-SAT algorithm

Date: 2026-11-24. Reliability level of everything below: **P** (written, not re-read by anyone) and
**CONDITIONAL** on the two unreviewed preprints (Anthropic `2610.00577v1`, OpenAI "ETH is false").
If the black box below is wrong, every theorem here is void.

## 0. Black box used (and nothing else)

**(BB)** = Corollary 11 of the Anthropic paper (the Main Theorem with explicit ranges): there is an
absolute constant C and a deterministic algorithm that, given a k-CNF (any k >= 1) with M >= 2 clauses in
which N >= 4 variables occur, outputs its exact number of satisfying assignments in time at most

    min( 2^{C k M / sqrt(log M)} ,  2^{C k^2 N loglog N / sqrt(log N)} ).

Caveats on (BB) itself: the main body of the paper proves only the weak Theorem 2 (2^{delta n} for every
delta); the M-bound is proved in App. A–B, the N-bound in App. D using App. C which is *a proof sketch*.
All my theorems that say "M-bound only" avoid App. C/D. The best theorem (Section 3) uses the M-bound only.

Two bookkeeping rules that decide everything:
* **R1 (M-bound).** A reduction to 3-CNF with M clauses gives 2^{O(M/sqrt(log M))}. A blow-up M = n·g(n)
  keeps a gain over 2^{O(n)} only if g(n) = o(sqrt(log n)). In particular a log n blow-up
  (Cook–Levin, permutation encodings, position labels, weight sums with poly weights) **loses everything**:
  n log n / sqrt(log n) = n sqrt(log n) >> n.
* **R2 (N-bound).** Number of clauses is free, but width costs k^2: with b-bit encodings of a domain of size d
  (b = log d, width 2b, N = nb) the exponent is O(b^3 n loglog n/sqrt(log n)); it beats 2^{nb} only if
  b = o((log n)^{1/4}/sqrt(loglog n)).

## 1. What is already derived (not duplicated here)

Paper §7: Cor. 5 (s_k = 0 for all k; ETH, SETH, NSETH false); E^NP has no linear-size circuits (via
Williams); Cor. 6 (#k-SAT in 2^{delta n}, #ETH false); Cor. 7 (#Independent Set, #Clique, #Vertex Cover,
#Dominating Set of size s, #k-colourings, each in O(2^{delta n}) for every delta — **no explicit rate**).
Explicitly not refuted: Orthogonal Vectors.

OpenAI side (O'Donnell's exposition): permanent of sparse matrices (O(1) non-zeros per row/column,
rational entries) in 2^{O(n/sqrt(log n))}; NTIME[T] simulation via a "schedule-conditioned encoding"
(2^{O(T loglog T/log T)} formulas of size O(T loglog T)).

Follow-ups found by web search (titles/snippets only; arxiv.org and most hosts are blocked for fetching in
this session, so abstracts were NOT read — only the GitHub README below was fetched):
* arXiv 2611.00994 "Subexponential-Time Exact-One Hitting Set".
* arXiv 2611.11832 "Simulating Nondeterministic Turing Machines in Subexponential Deterministic Time".
* arXiv 2611.01774 "Monotone Local Search with Sample-Dependent Running Times".
* arXiv 2611.02548 "Counting Feedback Vertex Sets of All Sizes in Subexponential Time, Conditionally".
* arXiv 2611.02650 "ETH is Dead, Long Live SH"; Szeider's "k x k Clique Hypothesis" note (replacement hypotheses).
* GitHub `xennygrimmato/eth-hamiltonicity` (README fetched; unreviewed, Lean written and reviewed by AI):
  on **sparse** graphs — Hamiltonian cycle parity (deterministic), decision (randomized, via Cut&Count +
  CIKP isolation), claimed exact #HC in 2^{O((m+n)/sqrt(log n))}, sparse TSP/ATSP with weights <= W in
  2^{delta n} poly(W), Eulerian circuits, quantum amplitudes, factoring and number-theoretic items.
  It lists **Graphic TSP, Steiner tree, FVS, longest path as "conjectured, not claimed"**, and says dense
  graphs are not covered.

Consequence for this note: my first idea (Hamiltonian cycle on sparse graphs through counting) is
**already claimed** there. I therefore prove instead items that the README lists as conjectured
(Steiner tree, longest path/cycle), with a self-contained argument and an explicit rate.

## 2. Candidates, with exact bookkeeping

Notation: n vertices / variables, m edges / clauses, L = sqrt(log n).

### C1. #CNF-SAT with unbounded width and m = c·n clauses; Dominating Set and Hitting Set with explicit rate
**Claim.** #CNF-SAT on n variables with m = cn clauses of arbitrary width (c >= 1) is computed in time
2^{O( n · sqrt( c · log(c log n) ) / (log n)^{1/4} )}; this is 2^{o(n)} iff c log c = o(sqrt(log n)).
**Proof sketch.** Count-preserving width reduction of [CIP06] exactly as in the paper's proof of Cor. 7:
while a clause has > w literals pick w of them, branch "delete the others" / "set these w false". A
root–leaf path uses the first branch <= m times and the second <= n/w times, so there are at most
(n/w+1)·(e(1+mw/n))^{n/w} leaves, each a w-CNF with <= m clauses, and each satisfying assignment lies in
exactly one leaf. M-bound per leaf: 2^{C w m/sqrt(log m)}. Total exponent
n·[ (1/w)·log(e(1+cw)) + C c w / L ] + O(log n); take w = sqrt(L·log(e(1+cL))/c).  ∎
**Corollary (Dominating Set).** n domination clauses + an O(n)-gate cardinality circuit (width 3, never
branched on): #dominating sets of size s in 2^{O(n sqrt(loglog n)/(log n)^{1/4})} (w = (log n)^{1/4} sqrt(loglog n)).
Same for Hitting Set / Set Cover with universe n and O(n) sets. (Going through the N-bound instead, as the
paper's proof literally does, gives only the worse (log n)^{1/6}.)
**Previous records:** CNF-SAT with cn clauses 2^{n(1-1/O(log c))} [CIP06]; Dominating Set 1.4864^n (Iwata).
**Novelty:** explicit rate only (the paper has "every delta"); low–medium. Uses M-bound only.
**GAIN LOST:** for m = n^{1+eps} clauses the leaf count forces w >> log n and then w·m/L >> n. So
*"CNF-SAT in (2-eps)^n·poly(m)" does NOT follow*: SETH is refuted only in its k-SAT form; the
CNF-SAT/Set-Cover-type conjectures with polynomially many long clauses survive this method (consistent
with OV not being refuted).

### C2. CSP with growing alphabet, graph homomorphism, k-colouring thresholds
**Claim.** #CSP with n variables over a domain of size d, constraints of arity <= r (any number of
constraints) in time 2^{O(r^2 log^3 d · n · loglog n / sqrt(log n))}. This is d^{o(n)} iff
log d = o( (log n)^{1/4} / (r sqrt(loglog n)) ).
**Proof sketch.** Binary encoding: N = n⌈log d⌉ variables; each forbidden tuple is one clause of width
r⌈log d⌉; unused codes are excluded by <= ⌈log d⌉ clauses of width <= ⌈log d⌉ per variable; bijection with
solutions. Apply the N-bound with k = r⌈log d⌉.  ∎
Special cases: #Hom(G,H) with |V(H)| = h in h^{o(|V(G)|)} for log h = o((log n)^{1/4}/sqrt(loglog n))
(the ETH lower bound "no h^{o(n)}" of Cygan et al. is inverted only in this range);
#k-colourings beat the 2^n of Björklund–Husfeldt–Koivisto iff log k = o((log n)^{1/6}/(loglog n)^{1/3}).
**GAIN LOST:** d or k above these thresholds; hence **Chromatic Number in 2^{o(n)} does not follow**, nor
Subgraph Isomorphism (n log n bits). Uses the N-bound, hence App. C/D (sketch) of the paper.
**Novelty:** threshold computation; medium-low.

### C3. Parameterized: Vertex Cover in 2^{o(k)}; k-Clique NOT in f(k) n^{o(k)}
* Vertex Cover: Nemhauser–Trotter kernel (<= 2k vertices), then #IS with a cardinality circuit: 3-CNF on
  O(k) variables, any number of clauses, N-bound: **time poly(n) + 2^{O(k loglog k/sqrt(log k))}**
  (record: 1.25284^k). Qualitatively known since Cai–Juedes 2003 (VC in 2^{o(k)} iff ETH fails); only the
  rate is new. Works for any problem with a linear-vertex kernel and a local (SNP-type) check.
* k-Clique on N vertices: best we get is 2^{O(N loglog N/sqrt(log N))}, which beats N^k only for
  k >= N loglog N/(log N)^{3/2}. The direct encoding (k slots × log N bits, width 2 log N) gives exponent
  O(k log^3 N ·loglog/sqrt(log(k log N))) >> k log N. **Direction check:** f(k)N^{o(k)} for Clique would
  refute ETH [CHKX06]; the converse is not known, so **no f(k) N^{o(k)} algorithm follows**. GAIN LOST.

### C4. Max-Cut / Max-k-SAT in terms of m; dense case lost
z_e <-> (x_u XOR x_v) (4 clauses of width 3) + an O(m)-gate counter: the number of cuts of each size in
**2^{O(m/sqrt(log m))}** (record 2^{m/5.2..}-type bounds); Max-k-SAT: 2^{O(k m/sqrt(log m))}.
**GAIN LOST** in terms of n for dense graphs (n^2 indicator variables; no sparsification for Max-Cut is
known): Williams' 2^{omega n/3} stays. Novelty: low.

### C5. Other items where the gain is lost (stated so nobody re-derives them)
* Subset Sum / Knapsack with b-bit numbers: adder circuits give M = O(nb), time 2^{O(nb/sqrt(log nb))};
  useful only for b = o(sqrt(log n)), where pseudo-polynomial DP is already polynomial. Lost.
* Hamiltonicity / TSP on dense graphs, Bandwidth, Treewidth/Pathwidth (orderings: n log n bits). Lost.
* TSP / Steiner tree with integer weights <= W through an *unweighted* counter: weight-sum circuit has
  O(m log W) gates, so gain only for log W = o(sqrt(log m)) (the repo's poly(W) claim needs a weighted counter).
* Whole Tutte polynomial of a sparse graph: single integer points (q, v) cost 2^{O(m log(q(v+2))/sqrt(log m))},
  but interpolation needs q, v up to m, i.e. a log m blow-up. Lost.
* Generic NTIME[t] via Cook–Levin/Pippenger–Fischer: size t log t. Lost (2611.11832 must use something else).
* Planar problems: already 2^{O(sqrt n)}; nothing gained.

### C6. (BEST) Connectivity problems on sparse graphs through an exact counter — Section 3
Steiner Tree, ℓ-Cycle / Longest Cycle / Longest Path (and, same lemma: Hamiltonian cycle, connected vertex
cover, max-leaf spanning tree, induced forests/FVS): counting mod p deterministically in
2^{O(m log p/sqrt(log m))}, decision with one-sided error in **2^{O(m loglog m/sqrt(log m))}**.
The log n blow-up of every known connectivity certificate is avoided because *counting* replaces certifying.

## 3. Full proof of the best consequence

Throughout, G = (V,E) is a graph with n vertices and m edges, m >= n-1 >= 2. "Counter circuit": the number
of ones among t bits is computed by a circuit of O(t) two-input gates (full adders; cf. [IPZ01, proof of
Thm 5], as used in the paper's Cor. 7); its Tseitin encoding has one new variable and <= 4 clauses of
width <= 3 per gate, and every input extends to exactly one satisfying assignment. Fixing the output bits
by unit clauses expresses "exactly s ones". Likewise "value(c) <= q" for a b-bit string c and a constant q
is an O(b)-gate circuit, hence O(b) clauses of width <= 3 with unique extension.

**Definition (skeleton formula).** A skeleton is a 3-CNF Phi over variables x_e (e in E), y_v (v in V) and
auxiliary variables such that (i) Phi contains the clauses (¬x_e ∨ y_u), (¬x_e ∨ y_v) for every e = uv, and
(ii) every assignment of (x,y) extends to at most one satisfying assignment. For a satisfying (x,y) write
X = {e : x_e = 1}, Y = {v : y_v = 1}; by (i) (Y,X) is a graph. Let cc(Y,X) be its number of connected
components (isolated vertices of Y count; cc = 0 if Y is empty).

**Lemma 1 (colouring lemma).** Let Phi be a skeleton with M_0 clauses, p >= 2 an integer, b = ⌈log p⌉, and
w : E -> {1,...,p} (multiplicities). There is a 3-CNF Psi = Psi(Phi,p,w) with at most
M_0 + c_1·b·(n+m) clauses (c_1 absolute), computable in polynomial time, such that

    #Psi = sum over (x,y) satisfying Phi of  ( prod_{e in X} w_e ) · p^{cc(Y,X)} .

*Proof.* Add for each vertex v a block c_v of b new variables, and
(a) the O(b) clauses (with O(b) new, functionally determined variables) saying value(c_v) <= p-1;
(b) for each bit i the clause (y_v ∨ ¬c_v^i), i.e. c_v = 0 if v is not in Y;
(c) for each edge e = uv and bit i the two clauses (¬x_e ∨ ¬c_u^i ∨ c_v^i), (¬x_e ∨ c_u^i ∨ ¬c_v^i).
Add for each edge e a block g_e of b new variables, and
(d) for each bit i the clause (x_e ∨ ¬g_e^i), i.e. g_e = 0 if e is not in X;
(e) the O(b) clauses (with functionally determined new variables) saying: x_e -> value(g_e) <= w_e - 1.
Fix (x,y) satisfying Phi. By (a),(b),(c) the admissible (c_v) are exactly the maps Y -> {0,...,p-1} that
are constant on each component of (Y,X) (vertices outside Y are forced to 0): p^{cc(Y,X)} choices. By
(d),(e) the admissible g_e number w_e if e in X and 1 otherwise, independently. Auxiliary variables are
determined. Multiply.  ∎

**Lemma 2 (nonvanishing).** Let p be prime and f in F_p[t_1..t_N] a nonzero polynomial of degree <= 1 in
each variable. Then f is nonzero on at least (p-1)^N points of F_p^N; i.e. a uniform point is a non-root
with probability >= (1-1/p)^N >= e^{-N/(p-1)}.
*Proof.* Induction on N; N = 0 is clear. Write f = g + t_N h with g,h in the first N-1 variables. If h != 0,
for each of the >= (p-1)^{N-1} points a with h(a) != 0 the polynomial f(a,t_N) has degree exactly 1, so
p-1 non-roots. If h = 0 then f = g != 0 has >= (p-1)^{N-1}·p non-roots. Finally ln(1-1/p) >= -1/(p-1).  ∎

**Skeleton for Steiner trees.** Let K ⊆ V be non-empty and 1 <= ℓ <= n-1. Let ST(ℓ) be the set of
subtrees T of G with K ⊆ V(T) and |E(T)| = ℓ. Phi_ST consists of: clauses (i); unit clauses y_v for v in K;
a counter circuit with "sum_v y_v = ℓ+1"; a counter circuit with "sum_e x_e = ℓ". It has
M_0 <= c_2 (n+m) clauses of width <= 3 and is a skeleton.

**Lemma 3.** For every (x,y) satisfying Phi_ST: cc(Y,X) >= 1, with equality iff (Y,X) in ST(ℓ); and every
T in ST(ℓ) arises from exactly one (x,y).
*Proof.* For any graph, cc = |Y| - |X| + beta with beta >= 0 the cyclomatic number; here |Y|-|X| = 1. So
cc = 1 iff beta = 0 and the graph is connected, i.e. a tree on Y ⊇ K with ℓ edges. Conversely a tree with
ℓ edges has ℓ+1 vertices and satisfies all clauses; (x,y) is its indicator.  ∎

**Theorem 1 (counting mod p, deterministic, conditional on (BB)).** For every integer p >= 2, the number
|ST(ℓ)| mod p is computed deterministically in time 2^{O(m log p / sqrt(log m))}.
*Proof.* Psi = Psi(Phi_ST, p, w ≡ 1) has M <= c_3 b m clauses of width 3 (using n <= m+1). By Lemmas 1, 3,
#Psi = p·|ST(ℓ)| + p^2·(integer), so (#Psi / p) mod p = |ST(ℓ)| mod p. (BB), M-bound with k = 3: time
2^{3C M/sqrt(log M)} <= 2^{3C c_3 b m / sqrt(log m)} since M >= m.  ∎

**Theorem 2 (Steiner Tree decision, randomized one-sided, conditional on (BB)).** There is a randomized
algorithm that, given a connected graph with m edges, terminals K and ℓ, decides whether there is a
connected subgraph containing K with at most ℓ edges, in time 2^{O(m loglog m / sqrt(log m))}; it never
errs on NO-instances and errs with probability <= 1/e^2 on YES-instances.
*Proof.* A minimum connected subgraph containing K is a tree, and in a connected graph a tree with ℓ' < n-1
edges extends to one with ℓ'+1 edges; so it suffices to decide ST(ℓ) != ∅ for ℓ := min(ℓ, n-1) (ℓ = 0 is
trivial). Let P = ⌈sqrt(log m)⌉ + 2 and p a prime in [P, 2P]; b = ⌈log p⌉ = O(loglog m).
Trial: draw u_e uniformly in {0,...,p-1} independently. Delete edges with u_e = 0 (add the unit clause
¬x_e); for the others set w_e = u_e. Build Psi = Psi(Phi_ST ∧ deletions, p, w), compute #Psi with (BB),
and accept iff (#Psi / p) mod p != 0. By Lemmas 1 and 3,

    #Psi / p  ≡  f(u) := sum_{T in ST(ℓ)} prod_{e in T} u_e   (mod p)

(trees using a deleted edge are absent, in agreement with u_e = 0). If ST(ℓ) = ∅ then f ≡ 0 and no trial
accepts. Otherwise f is a nonzero multilinear polynomial over F_p in m variables (distinct trees have
distinct edge sets, hence distinct monomials with coefficient 1), so by Lemma 2 a trial accepts with
probability >= e^{-m/(p-1)}. Run R = ⌈2 e^{m/(p-1)}⌉ independent trials: failure probability <= e^{-2}.
Cost: R = 2^{O(m/sqrt(log m))} trials, each one call of (BB) on M <= c_3 b m clauses of width 3, i.e.
2^{O(m loglog m/sqrt(log m))}, plus polynomial overhead.  ∎

**Theorem 3 (ℓ-Cycle, Longest Cycle, Longest Path; same statements).** Replace Phi_ST by Phi_cyc: clauses
(i); for each v a counter on its incident x_e with "y_v -> exactly two" and "¬y_v -> none" (O(deg v) gates,
total O(m)); "sum_v y_v = ℓ". Then every satisfying (Y,X) is 2-regular on ℓ vertices, cc >= 1 with equality
iff it is one cycle of length ℓ. Hence #ℓ-cycles mod p in 2^{O(m log p/sqrt(log m))} deterministically and
existence of an ℓ-cycle (ℓ = n: Hamiltonian cycle; maximise over ℓ: Longest Cycle) with one-sided error in
2^{O(m loglog m/sqrt(log m))}. For paths on ℓ vertices use "sum y = ℓ, sum x = ℓ-1, every degree <= 2":
cc = 1 + beta >= 1, equality iff a tree of maximum degree <= 2, i.e. a path. Proofs are verbatim those of
Theorems 1–2. On graphs of maximum degree Delta the bound reads 2^{O(Delta n loglog n/sqrt(log n))}.

**Comparison with what was known.** Steiner Tree: 2^{|K|} poly (Björklund–Husfeldt–Kaski–Koivisto,
Nederlof), about 1.36^n in n (Fomin et al.); 3^{tw} by Cut&Count. Longest path / Hamiltonicity: 1.657^n
(Björklund), about 1.2^n on cubic graphs. No 2^{o(n)} bound was known on bounded-degree graphs, and under
ETH none exists. Hamiltonian cycle on sparse graphs is already claimed in the `eth-hamiltonicity` repo
(with isolation instead of Lemma 2); **Steiner tree and longest path are listed there as conjectured** —
Theorems 2–3 settle them *on paper*, with an explicit rate, modulo review.

## 4. Adversarial self-review

1. **Everything is conditional** on (BB), i.e. on an unreviewed preprint written by an AI model (the same
   model family as the author of this note — no independence). Section 3 needs only the M-bound (App. A–B),
   not the sketched App. C/D. C2 and C3 need the N-bound and therefore inherit the "proof sketch" status of App. C.
2. (BB) is used with k = 3 and M up to O(m loglog m); it must hold with an absolute constant for all M,
   as Cor. 11 states. If the true statement had an O_k or required bounded density, C1 (growing w) would
   break; Section 3 would not.
3. **Not exact counting, not deterministic decision.** Theorems give counts mod p only for
   log p = o(sqrt(log m)); exact counts would need p > 2^m, i.e. m bits per vertex. Deterministic decision
   fails when the count is divisible by all small moduli. The repo claims exact deterministic #HC at
   2^{O((m+n)/sqrt(log n))}; I could not reconstruct how from an unweighted counter and did not verify it.
4. **Sparse only.** Exponents are in m. For dense graphs nothing is gained (and Steiner tree / Hamiltonicity
   in 2^{o(n)} on dense graphs remains underivable here).
5. Unweighted (edge-count) Steiner tree only; integer weights <= W cost a factor log W in the exponent.
6. Steps not double-checked by a second reader: the leaf count in C1; the constant-free optimisation of w;
   the linear size of counter/comparator encodings (standard); Lemma 2 (elementary, checked by hand).
   Nothing was machine-checked. The counting identity behind Lemmas 1 and 3 (S/p = #trees resp. #cycles mod p)
   was brute-forced on 6 random graphs with n = 7, m = 11, p in {2,3,5} (`check_identity.py`): this detects
   errors in the identity only; it tests neither the CNF encoding nor any running time (level N).
7. **Novelty is not established.** arxiv.org could not be fetched; I saw only titles and search snippets of
   the follow-up papers. Any of C1–C6 may already be in 2611.00994, 2611.02548, 2611.02650 or elsewhere.
   "Previous records" are quoted from memory and were not re-verified online.
8. These are complexity-theoretic statements only: by the paper's own estimate the M-bound beats brute force
   only for log m in the thousands, so no instance that can be written down benefits.

## 5. Sources
* Anthropic paper text: /tmp/claude-0/-home-user/2d12e629-e6d2-542b-bfc0-8b766f57383f/scratchpad/papers/claude_3sat.txt (§1, §7, Thm A, Cor. 11)
* O'Donnell exposition: same directory, odonnell_eth.txt (Thm 6.1 sparse permanent; NTIME remark)
* https://github.com/xennygrimmato/eth-hamiltonicity (README)
* https://arxiv.org/pdf/2611.00994 , https://arxiv.org/pdf/2611.11832 , https://arxiv.org/pdf/2611.01774 ,
  https://arxiv.org/pdf/2611.02548 , https://arxiv.org/pdf/2611.02650 (titles/snippets only)
