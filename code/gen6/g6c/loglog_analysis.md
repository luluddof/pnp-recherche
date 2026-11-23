# G6c: the factor log log n in the bound in n (Appendix D of "A 2^{o(n)} time algorithm for 3-SAT")

Status: **no proof that removes log log n.** What is here: (1) exactly where the factor enters;
(2) a proved sharpening (Prop. R): Appendix D actually pays log of the *average density* m/n,
not log of the maximum degree Δ; (3) proved obstructions for routes (b) and (d), and for route (a)
in the restricted class of "component-respecting" sparsifiers without new variables (Prop. L);
(4) route (c) is *not* refuted; the precise open sub-question is stated (Q2).
Proved statements are marked [P], heuristic ones [H], open ones [Q].

## 1. Where log log n enters

Notation of the paper: e = e_ℓ = Θ((log K)^{-1/2}) with K ≈ n^{1/17}, so e = Θ(1/√log n).
One network with total edge bits U costs 2^{O(eU)} (Step 4, (70)).

* (68): U ≤ 6k(n+log 2Δ) + 28kn·log(96Δ log 2Δ) + 2g².
  - kind 0 (clause marks): at most (k-1)(n̄+1) marks among T0 ≤ (k-1)m tests, cost
    (k-1)(n̄+1)·log(3m/(n̄+1)); the choice n̄ ≈ n/log(2Δ) makes this O(kn). **No log Δ here.**
  - kind 1 (incidence marks in parity checks): kM marks, M ≤ 14n, among T1 ≤ k r s̄ tests,
    cost kM·log(3s̄) with s̄ = 16m/n̄ ≈ 32 Δ log Δ. **This is the only place where log Δ is paid.**
* Step 5: Calabro's sparsification with ε ∈ [e, 2e] gives Δ ≤ 2(8k/ε)^{3k}+2, so
  log Δ = O(k log(k/e)) = O(k log log n). Total exponent 2en + O(e k n · k log(k/e))
  = O(k² n log log n / √log n).

So: log log n = log(1/ε) with ε ≈ e, entering through **log s̄ = log(m/n̄)**, the price of naming
which clause of a parity check carries a mark (≈ 14n such marks).

## 2. Prop. R [P]: only log(m/n) matters, not log Δ

**Prop. R.** Under the hypotheses of Thm D(a), with the Setting changed to: Φ a k-CNF,
3 ≤ k ≤ (log n)^{1/4}, maximum variable degree Δ ≤ n^{1/(40C0)}, and D := max{2, ⌈m/n⌉},
the number of satisfying assignments is computed in time 2^{O(e k n log D)}.

*Proof (check of every inequality of Steps 1–4 with D in place of Δ where Δ is used only through m ≤ Δn).*
Define n̄ as the power of two with n/(2 log 2D) < n̄ ≤ n/log 2D. Then m ≤ Dn gives
m/n̄ < 2D log 2D, s̄ = 16m/n̄ < 32 D log 2D, 3m/(n̄+1) < 6D log 2D ≤ (2D)³ (true for D ≥ 2:
6 log 2D ≤ 8D², e.g. 12 ≤ 32 at D = 2, and the right side grows faster).
Hence by (59) β ≤ n + 3(n̄+1) log 2D + 2 ≤ 4n + 3 log 2D + 2; r̄ ≤ 2.41β+1; (61) r̄ ≤ 10n,
M ≤ 14n for n ≥ n0 (needs only log 2D ≤ log 2Δ = o(n)); pr̄ < 22 log 2D; t̄ ≤ log m + 33 log 2D + 1.
Lemma 13 and Step 2 do not use Δ. (68) becomes U ≤ 6k(n + log 2D) + 28kn log(96 D log 2D) + 2g²
= O(kn log D) + 2g²: kind 0 costs ≤ (k-1)(n̄+1)·3 log 2D ≤ 3(k-1)(n + log 2D) and kind 1 costs
≤ kM log(3s̄) ≤ 14kn log(96 D log 2D). The only place where the true maximum degree is used is
(69): at most Δ⌈n/g⌉ clauses have a variable at a node, so gP = O(nΔ(log Δ + log n)) + O(n s̄),
and Step 4 needs gP/K = o(en), which holds because Δ ≤ n^{1/(40C0)} and s̄ ≤ 64 D log 2D ≤ 64Δ log 2Δ.
Everything else in Step 4 is unchanged. ∎

Consequences.
* The quantity paid is k·n·log(m/n) ≈ k·log C(m,n): the information needed to name ~n clauses
  out of m. Max degree may be as large as n^{δ}.
* The open question is therefore equivalent (for this route) to:
  **[Q1] Is there a count-preserving sparsification into 2^{O(n/√log n)} k-CNFs on the same n
  variables whose average density m/n is 2^{O(k)} (or merely (log n)^{o(1)} for an o(log log n)
  improvement), with max degree ≤ n^{δ}?** Max degree ≤ n^δ is free (Calabro with ε = n^{-δ/(4k)}).
  Calabro gives density (k/ε)^{3k}, i.e. (log n)^{Θ(k)} at ε = 1/√log n.

## 3. Route (b): re-balancing ε  [P, relative to the known density bound]

With density D(ε) = (k/ε)^{Θ(k)} the exponent is f(ε) = εn + c·e·k·n·log D(ε) ≥ n(ε + a log(1/ε))
with a = c' e k². The minimum of ε + a log(1/ε) is at ε = a, value a(1 + ln(1/a))/ln 2·(const)
= Θ(e k² log(1/(e k²))) = Θ(e k² log log n) for fixed k. Using the m-bound on the outputs instead
gives e k D(ε) n, worse. Mixing per output does not help since every output may have density D(ε).
So no choice of ε removes the factor; it is forced by log D(ε) ≥ log(1/ε).

## 4. Route (d): branching on high-degree variables  [P]

Branching on all variables of degree > Δ' costs up to 2^{km/Δ'} and leaves max degree ≤ Δ'.
(i) By Prop. R max degree is not the parameter that matters; average density is.
(ii) It cannot lower average density: take Φ in which every variable has degree exactly
d = km/n (e.g. a disjoint union of complete k-CNF blocks). If Δ' ≥ d nothing is removed and the
density is unchanged; if Δ' < d all n variables are branched on, cost 2^n.
So (d) is useless here; it only matters for the (already free) max-degree condition.

## 5. Route (a): better sparsification. Obstruction for component-respecting sparsifiers [P]

Instance: G_t = disjoint union of N cliques Q_1..Q_N of size t, n = Nt,
φ_t = ∧_{uv ∈ E(G_t)} (¬u ∨ ¬v) (a 2-CNF, so also a 3-CNF). Its solutions: for each clique,
"no vertex" or exactly one vertex true; |sat φ_t| = (t+1)^N.

**Line lemma.** Let ψ be a CNF (any width) on the same variables with sat ψ ⊆ sat φ_t. If two
solutions of ψ differ only in clique c, one with u true and the other with v true (u ≠ v ∈ Q_c,
all else equal = R), then ψ contains a clause whose literals on Q_c are exactly ¬u, ¬v plus
possibly positive literals, and whose other literals are false at R. Distinct pairs {u,v} need
distinct clauses.
*Proof.* ψ must reject w = (R, u=1, v=1). A clause falsified by w and satisfied by (R,u=1,v=0)
must contain ¬v (only v changed); symmetrically it contains ¬u. Any other negative literal ¬x on
Q_c would be true at w (x = 0 there), so there is none; hence the clause determines {u,v}. ∎

**Prop. L.** Call a sparsifier component-respecting if on a formula that is a disjoint union of
variable-disjoint parts each output is a conjunction ψ = ∧_c ψ_c with ψ_c on the variables of part c
(Calabro's algorithm is such: its sunflower branching on φ_t acts inside one clique at a time,
hearts being single literals ¬u). If φ_t is the disjoint union (or just union) of P outputs of
this form on the original variables, each with at most Dn clauses, then
log P ≥ (n/t)·log((t+1)/(2 + √(2Dt))). With t = 32D: log P ≥ n·log(3.3)/(32D) > n/(19D).
*Proof.* For an output, let A_c ⊆ {0} ∪ Q_c be the allowed values of clique c, a_c = |A_c|; the
output has Π a_c solutions. By the line lemma ψ_c has ≥ C(a_c−1, 2) clauses. So
Σ_c (a_c−1)(a_c−2)/2 ≤ Dn = DNt. Π a_c is maximized (log concave vs convex constraint) at equal
a_c = a with (a−1)(a−2) ≤ 2Dt, i.e. a ≤ 2 + √(2Dt) (since (a−2)² ≤ (a-1)(a-2)). The outputs
cover (t+1)^N solutions, so P ≥ ((t+1)/(2+√(2Dt)))^N. For t = 32D, √(2Dt) = 8D and
(32D+1)/(8D+2) ≥ 33/10 for D ≥ 1. ∎

So P ≤ 2^{εn} forces **D ≥ 1/(19ε)**, i.e. log D ≥ log(1/ε) − 4.3, and with ε ≤ 2e (needed so
that 2^{εn} is within the target) log D ≥ ½ log log n − O(1). Within this class of sparsifiers
the Appendix-D pipeline cannot lose the factor.
Caveats (not covered by Prop. L): (i) non-product outputs (2-CNF solution sets are majority-
closed, but correlated outputs might beat the bound; [Q] I conjecture the bound holds for all
outputs on the same variables); (ii) outputs with **new variables**: "at most one of t" has an
O(t)-clause count-preserving encoding with t auxiliary variables (sequential counter with
equivalences), so the clique instance is not hard if extension variables are allowed. A general
count-preserving "sparsification with extension variables" to density 2^{O(k)} with
n' = O(n) variables and 2^{O(n/√log n)} pieces would answer Q1-with-extensions positively; I
found neither a construction nor an obstruction. (Copying a variable along an equality chain does
not help: n' log(m/n') is minimized at n' = n.) Literature search found no lower bound on the
ε-dependence of sparsification density (CIP06, Santhanam–Srinivasan only treat unbounded width).

## 6. Route (c): cheaper encoding of the check marks  [not refuted]

Within Appendix D, kind-1 cost ≈ kM log(I/M), I = Σ_j |Λ_j| (incidences), M ≈ n + r, r = Θ(n)
(X can be ≈ 2^n·C(m, n̄+1), e.g. for an unsatisfiable formula, and one parity check removes at
most a constant fraction). So log log n disappears iff the checks can be chosen with
**I = O(n · 2^{O(k)})** (plus the load condition).
* Hashing the marks does not help directly: Ψ_{a,b} is an exact equality test; two ends must agree
  on *which* incidence is marked, which costs log C(T1, kM) bits.
* [P] Packing bound: if the family {V(ξ) : |V(ξ)| > n̄} contains L pairwise disjoint sets then
  I ≥ L (each needs its own odd incidence). This only gives I ≳ m/n̄, too weak.
* [H] If the violated sets just above the cut-off are numerous (2^{Ω(n)}) and pseudo-random,
  a union bound gives I ≳ n·m/n̄, i.e. exactly the log(m/n̄) per mark that the paper pays.
* But: random dense formulas have no assignment near the cut-off (every ξ violates ≈ m2^{-k}
  clauses, so checks of size ≈ 2^k suffice), and structured instances (e.g. blocks of
  (¬x_i ∨ x_j)) admit checks of size ≈ (block size), not m/n̄. So the heuristic lower bound is not
  realized by obvious instances.
**[Q2]** Is there, for every k-CNF with m clauses on n variables, a list of O(n) checks, found in
time 2^{o(n)}, that kills every ξ with |V(ξ)| > n̄ (cut-off n̄), has total incidence
O(n·2^{O(k)}) and has load O(log m)? A yes (with n̄ = Θ(n))
would remove log log n *and* the sparsification altogether (Prop. R with D replaced by I/n).
Note the kind-0 cost then forces n̄ ≈ n/log(m/n) unless kind-0 marks are also re-encoded; so Q2
must be asked with n̄ = n/log(m/n), where the packing/pseudo-randomness heuristics point to
I ≈ n m/n̄ = n·(m/n)·log(m/n). I regard a positive answer as unlikely but have no proof.

## 7. Self-review
* Prop. R: re-derived each numbered inequality; the one I am least sure of is that nothing in
  Step 3's P/B/τ bookkeeping uses m ≤ Δn beyond (69); I checked (67), (69), B = (k+1)M+1, τ.
* Prop. L: rigorous but only for component-respecting outputs on the original variables; the claim
  that Calabro's algorithm is component-respecting on φ_t is plausible from its definition
  (reductions and sunflower hearts within a component) but I did not re-verify it line by line
  in [Cal09].
* Route (b) is "forced" only relative to the paper's *upper* bound e·k·n·log D for one network;
  no lower bound on network cost is proved.
* Nothing here improves the paper's theorem; the factor log log n remains open.
