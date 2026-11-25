# G25-P3: Appendix D of OpenAI's ETH preprint with explicit constants (tensor-contraction route)

Agent notes (AI; not reviewed by a human). I take App. D as written and do not re-check its correctness.
Tags: [R] read in the paper (paper.txt lines 3007-3655; page 58 rendered to check (D.4)-(D.6));
[P] short proof here (unreviewed); [N] own numerics (scripts in this folder); [C] conditional on an
unproved claim; [O] open.
Scripts: appd_explicit.py / model2.py (cost model), run1.py (legal ε, paper rows), run2.py (rows A-F),
run5.py (rows G-I, constant), lp_kappa.py (taper LP), d5check.py ((D.5) size threshold).
Outputs: run1/run2/run5 *_out.txt (and the transcript numbers below).

Cost model (all rows). A call with M wires enumerates, at the worst node t, at most
exp(π·sqrt(2 n_t p_t / 3)) mode tuples. This is the (D.15) bound 2·sqrt(ζ(2)·n·p). Here
n_t = modes per unit weight and p_t = the precision of the node. A root-to-leaf path costs
Σ_d (π/ln2)·sqrt(2/3·n_t p_t)·[max over t] + M_leaf bits, and the worst path takes the largest M at
every depth (every term increases in M). The factor s^{O(1)} and the packet-entry cost
exp(O(L log s)) with L ≤ s^{1/25} are ignored: they are 2^{O(λ)} or 2^{s^{0.04}λ}, far below 2^{n}.
Instance used for thresholds: random 3-CNF at density α = 4.26 (s = 5.26 n). The target is
E/n < τ, with τ = 1 for brute force and τ = log2 1.307 = 0.386 for PPSZ. PPSZ is a decision
algorithm, so the PPSZ comparison is apples to oranges, as in b4.

## 1. Constraints (D.5)/(D.6) and the legal ε  (task 1)

Paper [R]: λ = ⌈log2 s⌉, r = 64, ε = 2^-17, Δ = ⌈s^{9/10}⌉, H = ⌈s/⌊√λ⌋⌉,
u_d = ⌊ελ/2^d⌋, a_d = 2Σ_{j≤d} u_j, L = 2^a, D_arr = r^{2a} = 2^{12a}.
(D.5) says u = Θ(λ/2^d), λ^{4/5} ≤ u ≤ M, a ≤ 4ελ, L, D_arr ≤ s^{1/1000}.
(D.6) says M ≥ 2 D_arr² Δ, M ≥ L·D_arr, 2L ≤ Δ.

Where each clause is used, and its status at ε = 1/960 [P]:
* M ≥ 2D_arr²Δ is the only binding clause. It is used once, in Lemma D.3 (H_pad ≤ 3M/D_arr, since
  √(2MΔ) ≤ M/D_arr). It must hold for every active call (M ≥ H+1) at the deepest active depth:
  24·a_dmax + 1 + log2⌈s^{0.9}⌉ ≤ log2(H+1), with a_dmax = 2Σ_{j≤dmax}⌊ελ/2^j⌋ < 4ελ.
  The first-order condition is 96ελ < 0.1λ, i.e. ε < 1/960. Because the depth is finite,
  a_dmax = 4ελ(1 - 2^{-dmax-1}). Exact values [N, run1]: ε_max = 1/833 (λ = 10^4), 1/899 (10^6),
  1/930 (10^8), 1/945 (10^10), 1/952.5 (10^12), tending to 1/960 from above.
  **So ε = 1/960 is legal for every λ ≤ 10^12 (and ε < 1/960 is legal for all large λ).
  2^-17 is 136.5x below what (D.6) needs.**
* M ≥ L·D_arr (13a ≤ log2 H) and 2L ≤ Δ: slack at ε = 1/960.
* "L, D_arr ≤ s^{1/1000}" is literally FALSE for ε > 1/48000, because D_arr ≤ 2^{48ελ}. It is used
  only (i) to derive (D.6), already covered; (ii) for the packet-entry cost exp(O(L log s)) ≤
  exp(O(s^{9/10})), which still holds since L ≤ s^{4ε} = s^{1/240}; (iii) for "fewer than s^4
  interpolation points", which still holds. So the clause should read L ≤ s^{1/240},
  D_arr ≤ s^{1/20}. Every lemma survives.
* "λ^{4/5} ≤ u" is used only in (D.11), to absorb inherited precision into O(M·D_t/u). It is a
  running-time device, not a correctness one. Taken literally it is a size condition the paper never
  states [N, d5check.py]:
  - with ε = 2^-17 it needs log2 λ ≥ 85 at d = 0, and log2 λ ≥ 320 at the deepest active depth,
    i.e. s ≥ 2^(2^320);
  - with ε = 1/960 it needs log2 λ ≥ 185.
  The proof's "sufficiently large s" is therefore astronomically larger than the operation-count
  crossover. In the explicit version (row G) this clause is replaced by an explicit slack c_d in a_d,
  and the inherited precision is tracked exactly in the count.
* Lemma D.5 needs 1 ≤ u ≤ a and M ≥ u. With the paper's ε these need λ ≥ 2^{17+d} (u_d ≥ 1), which
  is harmless at the crossover sizes.
Verdict: ε = 1/960 (in fact ε_max(λ) ≈ 1/950) is legal. Every lemma of D.3–D.6 still holds once the
(D.5) wording is fixed as above. The paper's 2^-17 is unnecessary.

## 2. Improvements one at a time  (task 2)
For each step: the lemma that changes, why it stays correct, and the crossover λ* = log2 s*.
λ* scales with the square of the exponent constant.

| row | change (cumulative) | lemma changed | status | λ* vs 2^n | λ* vs PPSZ |
|---|---|---|---|---|---|
| A | paper as written: ε = 2^-17, p_t = p+8hD+8M(a+1), h = ⌈3M/u⌉ (κ = 24), ν = 6, D.2 network (equality trees + subdivision: 2(km+Σ(d_v-3)^+) ≈ 45.1n wires), child ≤ 3M/32 | — | [R] | 7.0e11 | 4.7e12 |
| (A') | b4's model of A (11.78 wires/var, κ = 18) | — | b4 underestimated A by 20x | 3.6e10 | 2.4e11 |
| B | ε = ε_max(λ) ≈ 1/950 | (D.4), (D.5) wording | rigorous [P] | 5.05e9 | 3.4e10 |
| C1 | p_t = p + Λ_i·D_t + (q_t-1)(a+1-i), i.e. the loss the proof of D.5 actually establishes (6hD, not 8hD; evaluation loss (q-1)(a+1-i), not 8M(a+1)) | Lemma D.5 statement | rigorous: the proof shows that this p_t suffices | 3.8e9 | 2.5e10 |
| C2 | h = ⌈M/(u-1)⌉ instead of ⌈3M/u⌉ (κ: 18 → 6) | Lemma D.5 (first line: Σq_t = hL(2u-2+2^{1-u}) ≥ 2ML > M(2L-1), counting all points rather than only the odd ones; q_t ≤ 2M) | rigorous [P] | 1.26e9 | 8.5e9 |
| (C3) | LP-optimal taper | Lemma D.5 (any profile; (D.10) holds for arbitrary q(i)) | rigorous but worthless: LP κ = 5.96 (u = 100), 5.99 (u = 300) against linear 6.06/6.02; gain < 0.3% at u ≥ 300 [N, lp_kappa.py] | ≈ C2 | ≈ C2 |
| D | packing slack: packets ≤ P·M/L and child ≤ 2PM/r, with P = 1 + D_arr√S2/M + L·D_arr/M; modes only in packets whose two symbols differ in coordinate i (ν = 2P(1-1/r) instead of 6); ε re-optimised | Lemma D.3 (1),(3); (D.13) | rigorous: same proof, keeping the exact H_pad ≤ M/D_arr + √S2 + L. (D.6) becomes a soft trade-off | 3.25e8 | 2.2e9 |
| E | network: clause–variable incidence graph (already bipartite; one equality vertex per variable; split variables of degree > Δ0 with binary equality links). M0 = km(1+O(1/Δ0)) = 12.78n instead of 45.1n | D.2 | rigorous [P]: the contraction is #SAT on the occurring variables, as in D.2 | 2.44e7 | 1.6e8 |
| F | degree bound adaptive, not s^{9/10}: S2 = Σ_v deg² ≤ M(Δ0 + 2L_{d-2} + 2L_{d-1}) (old vertices, older packet vertices, new packet vertices) | Lemma D.3 variance bound (the proof counts pairs sharing an endpoint = Σ deg²); (D.6) | rigorous for bounded-degree inputs (the D.2 application); not for the general Thm D.1 with degree s^{9/10} | 2.6e6 | 1.7e7 |
| G | a_d = a_{d-1}+u_d+c_d (c_d ≈ 10-17) instead of a_{d-1}+2u_d; u_d ∝ M_d^{2/3}; depth 2-3 with brute-force leaves; r = 64 (r = 32: +0.6% in c [N, λ=1e6]; r = 128: +1.4% [closed form]) | (D.4); (D.11) re-proved by tracking the inherited precision exactly (D_t ≥ 2^{a_{d-1}+c_d+1}); leaf rule | rigorous [P,N] | **8.97e5** | **6.06e6** |
| H | Smith/Mahler mode basis (n_t → n_t/4: half the exponent) | (D.12)-(D.16): modes Z_k = (V·X_block)_k, weights s_k from the Smith form of the crossing block over the DVR O | algorithm correct [P]; running-time bound **conditional** [C]: needs the Smith profile s_k = (a-i+1)+2k+2v2(k!) (G18d: 23 exact cases, a ≤ 6, unproved) | 2.25e5 | 1.5e6 |
| I | precision κ = 6 → 2 (black-box-optimal recovery instead of Hermite formula (D.9)) | Lemma D.5 replaced by exact lattice recovery | **open** [O]: G18a measures κ_BB = 2.3-6.6 only on toys (a ≤ 3, u ≤ 3); the proved witness floor is ~2 (G17a). κ = 2 is a lower bound on what any recovery achieves, not an achieved value | 7.6e4 | 5.3e5 |

Not done / heuristic: a derandomised Chernoff/Freedman labelling (loads ≤ (1+δ)M/D_arr once
M ≳ D_arr·Δ·log D_arr, instead of D_arr²·Δ). This would change 24a+a' to about 13a and give
about 1.9x in λ* [H].

Corrections to the record (b4, g18d):
1. **b4 row A (4.0e10) underestimates the paper as written by about 20x.** The D.2 network has
   45.1 wires/var (trees + subdivision), not 11.78, and the stated budget is 8hD (κ = 24), not 18.
   The true row A is 7.0e11.
2. **The "LP taper, 1/√3 in the exponent" (g18d/b4) is not a taper-shape effect.** It is the
   constant in h = ⌈3M/u⌉. The paper's own linear taper with h = M/(u-1) already reaches κ → 6,
   which is the asymptotic LP optimum. G18d's 5.84 is a u = 30 value.
3. (D.5) as printed is false at ε = 1/960 (L, D_arr ≤ s^{1/1000}), but only as a convenience bound.
   Its λ^{4/5} ≤ u clause hides a size threshold s ≥ 2^(2^320).

## 3. Final explicit statement  (task 3)

Rigorous part (rows A-G), conditional only on App. D being correct as written [P, N]:
**#k-SAT with n variables and m clauses (s = n+m, λ = log2 s, bounded density) can be counted in time
s^{O(1)} · 2^{E}, with E ≤ c·k·m/√λ, c = 75** (r = 64, ν = 2, κ = 6).
* Closed form: c = (π/ln2)·sqrt(2νκ/3)·sqrt(4·log2 r + 1)·(1-(2/r)^{2/3})^{-3/2}
  = 12.82 · 5 · 1.170 = 75.0. Here 4·log2 r + 1 = 25 comes from 24a_d + a_{d-1} ≤ λ, and
  Σ_d u_d ≈ λ/25.
* Numerically, with every lower-order term included (floors, inherited precision, evaluation loss,
  padding P, leaves), E·√λ/(km) = 73.8 (λ = 10^5), 74.2 (10^6), 74.4 (10^7), 74.5 (10^8), rising
  towards the closed form 75.0.
* Since km ≤ k·s: E ≤ 75·k·s/√(log2 s). For the general Thm D.1 (degree ≤ s^{9/10}, row-D
  hypotheses, closed form only, not run): E ≤ about 232·M/√(log2 s) with M = number of wires.
* Thresholds for random 3-CNF at density 4.26:
  - beats 2^n for log2 s ≥ 9.0·10^5;
  - beats PPSZ 1.307^n for log2 s ≥ 6.1·10^6.
  The paper as written needs 7.0·10^11 and 4.7·10^12: a factor of about 780 in log2 s.
* Conditional on the Smith profile (row H): c = 37.5; λ* = 2.25·10^5 (2^n) and 1.5·10^6 (PPSZ).
* With the open κ = 2 (row I): λ* = 7.6·10^4 and 5.3·10^5. This is a floor for this design, not a result.

## 4. Rigour ledger  (task 4)
* Rigorous given App. D (re-proofs here, unreviewed): B, C1, C2, (C3), D, E, F (bounded-degree
  inputs), G.
* Conditional: H. The algorithm with Smith modes is correct (O is a DVR; the Smith form of each
  crossing block respects the source/target groups). Its running-time bound needs the Smith-profile
  formula, which is numerical only (G18d: a ≤ 6, 23 cases). G18d's C1 (flattening invariants =
  Smith-monomial weights, L ≤ 8) is the matching optimality conjecture. Neither is proved.
* Open: I (κ below 6: needs the black-box recovery precision at large u; only toys) and the
  derandomised-concentration step.
* Not addressed: App. D's correctness itself; the exponent of the s^{O(1)} factor.

## 5. Self-review
* The worst-path argument assumes each per-level term and the child bound increase in M. Checked
  for M·sqrt(P(M)) and 2P(M)M/r.
* The tuple count uses the (D.15) upper bound, not the exact saddle point (exact differs by
  O(n_t log) only).
* The S2 bound for children was first written as M(Δ0+2L_{d-1}). Old packet vertices (2L_{d-2})
  were added afterwards. The effect is below 0.5 bit in a ~λ/2 margin (rows F, G rerun).
* Random-instance numbers: 3m wires and d_v ~ Poisson(12.78) (for the paper's D.2 network). Worst-case
  formulas of the same density give the same c per clause-port.
* PPSZ is a decision algorithm (not counting). The τ values are as in b4.

* Run log: run5.py was stopped by its 590 s limit after λ = 1e8 (the r-scan was not rerun; r = 32 comes from the
  earlier run3.py scan at λ = 1e6: 74.58 vs 74.16 for r = 64). run2_out.txt was rerun after the S2 fix with identical numbers.
* No process left running; nothing written into /home/user/pnp-recherche.
