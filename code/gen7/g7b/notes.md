# G7b — The DP-type cost model (Limits paper, §5.2 / §9), after G6a

Notation as in G6a/notes.md: point t, groups S_h (a_h = |S_h|), r_h = dim R_{t,h}, ρ_t = Σ r_h,
measured set M_t = {h : r_h > 0}, interaction graph I_t, w_t = tw(I_t), torso G_M (I_t[M] plus a
clique on N(h) for every unmeasured h), Z_t = {o : c_o(t) = 0}, τ = tw(Pr F). For a witness map ω,
X_h = |I^ω_{t,h}|; for uniform random ω: μ_h = E X_h, V_h = Var X_h. c0 = 0.216.

## 0. The DP-type model, quoted precisely

§5.2: "We charge each point t at least κ(t) ≥ max(2^{s_t}, 2^{w_t}, C(n, min(q_t, ⌊n/2⌋)))
(DP-type evaluation)", and product-type additionally 2^{ρ_t}. w_t = tw(I_t), where for each
clause e the groups {α_t(e)} ∪ {h ≠ α_t(e) : P_{e,v}(t)|_{S_h} ≠ 0 for some v ∈ e} form a clique.
Justification given: readings of measured groups are summed by variable elimination along G_M;
"a nonzero measurement row takes at least two values on Boolean inputs, so it costs at least
2^{w_t}". So a group's whole reading vector is ONE elimination variable, charged a factor 2
regardless of r_h; rows in different groups enter a DP table only through bags of G_M.

## 1. Result O1: in the literal DP-type model the gap is closed — trivially, at 2^{Θ(√n)}

**Proposition O1.** For every CNF F on n variables (every R ≥ 0) there is a scheme satisfying
Definition 5.1 with, at every point, s_t = ⌈√n⌉, w_t ≤ ⌈√n⌉ − 1, q_t = 0. Its literal DP-type cost is
≤ 2^{⌈√n⌉}. Together with Theorem 8.1 of the paper, on any family with τ ≥ cn the optimal literal
DP-type cost is 2^{Θ(√n)}: exponent between √((cn+1)/2) − 1 and ⌈√n⌉.

Proof. P_o(T) := e_v (constant; degree ≤ R, P_o(0) = e_v). Then H(T) ≡ #F, deg H ≤ MR trivially.
Recovery: X = any MR+1 distinct nonzero rationals, q_t = 0, Λ = Lagrange interpolation at 0 (exact
for deg ≤ MR). Design (same at every t): s = ⌈√n⌉, consecutive groups of size ≤ s, g = ⌈n/s⌉ ≤ s
groups, α arbitrary, R_{t,h} = Q^{S_h}. (3) holds since every restriction lies in Q^{S_h}. I_t has
g vertices, so w_t ≤ g − 1 ≤ s − 1; C(n, 0) = 1. ∎

Meaning. The literal 2^{w_t} charge ignores how many values a group's reading takes. For this design
a reading is the full assignment of the group; an honest DP table over a bag of w+1 groups has
2^{Σ a_h} ≥ 2^{τ+1−s} entries (computed: W_t = n at n = 12..30). So the literal DP model certifies
nothing beyond 2^{O(√n)} for ANY technique (the paper's remark "our lower-bound technique cannot
exceed 2^{Õ(√n)}" is in fact a property of the model). The DP-type question is only meaningful with a
charge that counts readings. G6a/G6r's "n^{1/4} rows in √n groups" open case is one instance.

## 2. Refined (reading-weighted) DP model

Each measured group h with r_h independent rows takes ≥ 2^{r_h} distinct readings on {0,1}^{S_h}
(choose r_h columns of the row matrix forming an invertible minor; restricted to inputs supported on
those columns the reading map is injective). Following the paper's own elimination argument with
each group reading as one variable, a table whose scope is a bag B costs ≥ 2^{Σ_{h∈B} r_h}. Define
    W_t := min over tree decompositions of G_M of max_bag Σ_{h∈B} r_h,
    κ_ref(t) ≥ max(2^{s_t}, 2^{W_t}, C(n, min(q_t, ⌊n/2⌋))).
W_t ≥ tw(G_M)+1 ≥ w_t (every measured group has r_h ≥ 1), so κ_ref ≥ κ_DP; also W_t ≤ ρ_t, so the
product-type model implies the refined one. [Model proposal; an evaluator that splits a group's
reading into separately eliminated rows is not covered — same caveat class as the paper's.]
Trivial consequence: if every measured group has ≥ r rows, W_t ≥ r·w_t (so √n groups with n^{1/4}
rows each and w_t ≍ √n give W_t ≳ n^{3/4}).

## 3. Theorems in the refined model

**Theorem C (variable width).** At every t: τ ≤ W^a_t + s_t − 1 + |Z_t|, where W^a_t is the min over
tree decompositions of G_M of max_bag Σ_{h∈B} a_h. Moreover, for every tree decomposition of G_M there
is a torso bag with Σ_{h∈B} a_h ≥ τ + 2 − s_t − |Z_t|.
Proof. As in Theorem 7.1(i) of the paper, every edge of Pr(F') (F' = F minus Z_t occurrences) joins
two variables of equal or I_t-adjacent groups (checked exactly). Unmeasured groups are pairwise
non-adjacent (every edge of I_t has a measured end). Take a decomposition of G_M; replace each
measured group by its variables; for each unmeasured h attach the bag S_h ∪ ⋃_{h'∈N(h)} S_{h'} to a bag
containing the clique N(h). All edges of Pr(F') are covered (measured–measured: edge of G_M;
anything touching unmeasured h: its attached bag), running intersection holds, and the attached bags
have size ≤ s + Σ_{N(h)} a ≤ s + (that torso bag). Add the ≤ |Z_t| variables of Z_t to every bag. ∎
**Corollary C1 (full measurement gains nothing).** If every measured group is fully measured
(r_h = a_h), W_t = W^a_t, so at the point of Lemma 6.2 (|Z_t| ≤ q_t): max(s_t, W_t, q_t) ≥ (τ+1)/3,
and κ_ref ≥ 2^{min((τ+1)/3, ⌊n/2⌋)}.

**Theorem D (one point for both budgets).** Let N^ω_t := Σ_h (a_h − r_h − X^ω_h)_+. If every
variable occurs (n ≤ M), there is t ∈ X with |Z_t| + E_ω N^ω_t < 2(q_t + 1).
Proof. G6a's Lemma 1 proof shows nullity Φ^ω(t) ≥ N^ω_t (blockwise rank ≤ min(a_h, r_h + X_h)) and
Σ_t nullity ≤ deg det Φ ≤ nR. Lemma 6.2's G = Π_o c_o is nonzero of degree ≤ MR with ord_t G ≥ |Z_t|,
so Σ_t |Z_t| ≤ MR. Hence Σ_t (|Z_t| + E N_t) ≤ (M+n)R ≤ 2MR < 2Σ_t(q_t+1), using (R1). ∎
(So the treewidth condition of Thm 7.1 — with 2q_t+1 for q_t — and Theorem B hold at the SAME point.)

**Lemma E (per group).** If μ_h ≤ a_h then E(a_h − r_h − X_h)_+ ≥ (c0/2)·min(V_h, √V_h) − r_h.
Proof. (y − r)_+ ≥ y_+ − r; E(a − X)_+ = ½(E|X − a| + a − μ) ≥ ½E|X − a| ≥ (c0/2) min(V, √V) by G6a
Lemma 2's anticoncentration (valid for every real a; reviewed by G6r). ∎

**Theorem F (rows must sit in the bags).** At every t, for every tree decomposition of G_M and every
bag B: Σ_{h∈B} r_h ≥ (c0/2) Σ_{h∈B, μ_h≤a_h} min(V_h, √V_h) − E N_t. Hence, with
Γ_t := min_decomp max_bag Σ_{h∈B, μ_h ≤ a_h} (c0/2) min(V_h, √V_h):  W_t ≥ Γ_t − E N_t, and at the
point of Theorem D, W_t > Γ_t − 2q_t − 2.
Proof. Sum Lemma E over the non-surplus groups of B; the dropped summands E(·)_+ are ≥ 0, so
Σ_{B,ns} E(a_h − r_h − X_h)_+ ≤ E N_t. For W_t take the W-optimal decomposition. ∎

**Corollary G (regular designs: n^{2/3} also in the DP model).** Call the design (γ,β)-regular at t if
in every torso bag, the measured groups h with μ_h ≤ a_h and V_h ≥ β a_h ≥ 1 hold ≥ γ·Σ_{h∈B} a_h
variables. At the point of Theorem D, if the design is (γ,β)-regular there:
    W_t + 2q_t + 2 > (c0/2)·γ·√(β/s_t)·(τ + 1 − s_t − 2q_t).
Proof. Theorem C gives (in the W-optimal decomposition) a bag with Σ_B a_h ≥ τ + 2 − s − |Z_t| ≥
τ + 1 − s − 2q_t; regular groups there have √V_h ≥ √β·√a_h ≥ √β·a_h/√s; apply Theorem F. ∎
On τ ≥ cn (bounded-occurrence expanders, random 3-CNF via Lemma 8.2): every scheme that is
(γ,β)-regular at all points has max(s_t, W_t, q_t) ≥ c'(γ²β)^{1/3} n^{2/3} at some t, so refined
DP-type cost ≥ 2^{Ω(n^{2/3})} — the analogue of G6a Cor. 2 with W_t in place of ρ_t.

## 4. What is NOT proved (self-review)
* The refined model is a proposal (justified like the paper's 2^{w_t}), not the paper's model. In the
  paper's literal DP model nothing beyond 2^{Θ(√n)} can hold (O1) — that part is a theorem.
* Corollary G needs regularity. It FAILS for G6a's natural designs: 32–44% of measured groups are
  surplus (μ_h > a_h) and only 25–100% of measured variables are in regular groups (dp_checks.txt).
  Surplus or low-variance measured groups need no rows by Lemma E's argument; the global identity
  Σ_h(a_h − μ_h) = 0 bounds total surplus only by ρ_t + O(q_t), and ρ_t is not charged here. So the
  unconditional refined-DP statement "max(s, W, q) ≥ n^{1/2+δ}" remains OPEN. Escape route to rule
  out: separators made of surplus "clause-hub" groups whose rows-heavy (deficit) neighbours avoid the
  large bags.
* Unconditionally proved in the refined model: W_t ≥ w_t; W_t ≥ r_min·w_t; Corollary C1; Theorem F
  (W_t ≥ Γ_t − 2q_t − 2) at the same point where τ ≤ (w_t+1)s_t + 2q_t (Theorem D + Thm 7.1(i)).
* Constants are weak: at n ≤ 30, Γ_t ≤ 0.75 while W_t ≥ 9, and Theorem 8.1 is vacuous; the numerics
  are consistency checks, not tightness tests.

## 5. Computation (dp_checks.py → dp_checks.txt; imports G6a's scheme.py/designs.py read-only)
* Part 1 (O1): random 3-regular formulas n = 12,16,20,24,30; trivial design checked exactly
  (condition (3), |Z_t| = 0, q_t = 0, |X| = MR+1); exact tw(I_t) (subset DP) = g−1; literal exponent
  max(s,w) = ⌈√n⌉; exact tw(Pr F) = 5, 6 for n = 12, 16; reading-weighted W_t = n.
* Part 2: G6a's four design families at (n,|X|,R) = (24,8,2), (30,10,2), condition (3) exact:
  strong Lemma-1 form nullity Φ^ω(t) ≥ Σ_h(a_h − r_h − X_h)_+ (3 ω's, exact ranks), Σ_t nullity ≤ nR,
  Σ_t|Z_t| ≤ MR, Theorem D, Lemma E (min slack 1.39–5.82), Theorem F (subset form on all subsets of
  measured groups, and W_t ≥ Γ_t − E N_t with exact weighted treewidth), Theorem C's blow-up step
  (every Pr(F') edge joins equal/adjacent groups; unmeasured groups non-adjacent). All assertions pass.
