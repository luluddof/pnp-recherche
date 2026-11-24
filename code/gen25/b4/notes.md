# G25-B4 — "Lower rung": concrete improvements of the two Oct-2026 algorithms (brainstorm + quick numbers)

Tags: [P] short argument here (unreviewed), [N] own numerics (scripts in this folder), [H] heuristic/model, [R] read in paper.
Papers taken as written (untrusted, correctness not assessed). Nothing here changes an exponent of log.
Scripts: appb_rho.py (App. B at the lambda the algorithm uses), claude_cross2.py (ablation of Claude's crossover),
claude_curve.py (best exponent/n vs log2 m), openai_cross.py (App. D crossover model). Outputs are in the text below.

## 0. Key observations (they reorder the caller's list)

(O1) [R+N] **Claude's recursion never meets a small crossing.** In Thm A(b) the pieces are applied at lambda = the
bits per pair of groups, and lambda only GROWS going down (L' ~ rho g^4 lambda / K^2). At the top
lambda = L ~ 2km/K^2 >= 2k m^{7/9} (paper, K <= m^{1/9}); with the overhead terms (c log 2lambda, C0 K^2 log Xi)
lambda must anyway be >> C0 log(Xi)/e ~ 10^4..10^7. App. B's effective rho at crossing size lambda (appb_rho.py,
counting formula of G20b, 4rho = max_s [log2 P(2h+s-1) + log2 P(2h-s)]/lambda, p(h) >= 2^lambda):
   lambda:      8      16     64     256    1024   4096    2^16     2^20
   rho/(1/√2): 1.40   1.25   1.09   1.032  1.010  1.0033  1.0003   1.00002
(Lemma 2's octonion pieces, 4rho = 11/3, beat App. B only for lambda <~ 12.) So the FINITE-SIZE FACT (App. B trivial
for N = 2^lambda <= 256) is about lambda <= 8 and is irrelevant to the algorithm: at every lambda it actually uses,
App. B is within 0.3% of 1/√2. "Best pieces at each size" (Hurwitz sums, JH, exact small solutions) buys ~0.

(O2) [R+N] **The finite-size cost of Claude's algorithm sits in the recursion bookkeeping, not in the pieces.**
Model (Thm A, (40)-(45)): exponent = 2 e_l W (W = #wires), e_{l+1} = (rho+3zeta_l) e_l + 3/K_l,
log K_{l+1} = 2(log K_l + log(g_l/K_l)), top level log m >= 9 log K_l. Reproduces the paper's own remark exactly
(r = 1: l = 4, log2 m* = 8.05e3). Since e ~ e0 sqrt(A/log K) (log K ~ A 2^l), log m* ∝ a_top · A · (e0 W/n / tau)^2:
constants enter SQUARED, and e0 (the base case) and a_top are the cheapest levers.

(O3) [N] **Claude's architecture has a ~20x better finite-size constant than OpenAI's App. D**, after both are
re-optimised: effective exponent ≈ 4.6 W/sqrt(log2 m) (Claude, rows C–E) vs ≈ 108 W/sqrt(log2 s) (OpenAI, row E).
So the lower rung should be built on Claude's recursion.

## 1. Idea list (a what, b gain, c risk, d done?, e first computation < 1 h)

### I1. "Best pieces at each size" in Claude's recursion (the caller's idea 1)
a) Replace App. B by the best known decomposition at the actual crossing size.  b) ≤ 0.3% on rho for lambda >= 4096,
≤ 1% at lambda = 1024; the algorithm uses lambda >= 10^4 (O1). Gain ≈ 0. c) none — it is simply dead.
d) Not done before in this form: G20b/G24a computed the finite-N frontier, not the lambda the algorithm meets.
e) Done here (appb_rho.py). Verdict: drop; record the table above as the reason.

### I2. Cut-width (sweep) base case: e0 = 1/2 -> 1/4 + 1/(2K0)  [P, elementary; NEW as far as grep shows]
a) Thm A(b), l = 0, sums the (K0,lambda)-clique by brute force (2^{K0(K0-1)lambda/2}, e0 = 1/2). Sweep the nodes
1..K0 keeping the partial sum as a table on the cut edges {1..j}x{j+1..K0}: step j costs 2^{lambda(j(K0-j)+K0-j-1)},
max ≈ (K0^2/4 + K0/2) lambda. So e0 = 1/4 + 1/(2K0) (time; space 2^{K0^2 lambda/4}).
b) Worth two levels (rho^2 = 1/2) -> log K /4: density 4.26 vs brute force, log2 m* 9.55e4 -> 2.41e4 (claude_cross2.py
row C -> D); vs PPSZ 5.2e5 -> 1.3e5. Trilinear/fast-matrix-multiplication splits give ≥ 0.25 too (5/18 or ω/9), so 1/4 is
the natural floor of simple schemes.  c) Very low (standard variable elimination; only time matters in Thm A).
d) grep "cutwidth"/"largeur de coupe": nothing in the repo.  e) Write the 5-line lemma + rerun claude_cross2.py with it
as the only change against the paper's other constants (done: rows C/D).

### I3. Exact re-optimisation of the balancing constants (x0, zeta_l, g_l/K_l)  [N]
a) The paper uses x0 = 18, zeta_l = 2^{-l-5}, g_l/K_l = 2^{3l+x0} (crude check of Hoeffding). Use the exact Hoeffding
condition 2 zeta^2 g >= K ln(2gK) with g/K minimal, x0 free (K0 must still beat the 3/K_l term), geometric zeta schedule.
b) density 4.26 vs brute: 1.37e5 -> 9.55e4 (row B -> C); density 1: 3.78e3 -> 2.36e3.  c) low.  d) not done.
e) done here (row C); a clean statement needs the explicit constant c of App. B pieces (I8).

### I4. Top level with m ≈ K^2 instead of K^9 (random clustering + Lemma 8 case (ii))  [H]
a) g = 3K^8 clusters (case (i), arbitrary networks) forces log m >= 9 log K. A random (or enumerated) partition of the
clause chain into g = (g/K)·K clusters makes all pair loads c_ab ≈ 2W/g^2 (concentrated if W/g^2 >> log g), so case (ii)
applies; the remaining overheads (piece assignments (2L)^{cK^2/4}, C0 K^2 log Xi, 2^{2gP/K}) need only
log m >= 2 log K + O(log(C0/e)) (~ +15..25).  b) factor 4.5 in log m*: 2.41e4 -> 5.4e3 (row D -> E, density 4.26).
c) medium: derandomisation (or a randomised counting algorithm), and the unknown constant C0 enters b_top.
d) G4c (a) noted "K <= m^{1/9} instead of m^{1/2}" as a constant sqrt(9/2) in e; here it is a factor 4.5 in log m*.
e) Write Lemma 8' (random clustering, case (ii)) with explicit failure probability; rerun row E with b_top = 25.

### I5. Deterministic "granularity" balancing g/K = O(1/zeta) instead of zeta^{-2} log  [H, speculative]
a) Each node carries ~1/g of the load, so a greedy/LP-rounded partition plausibly balances to (1+zeta) with g/K ~ C/zeta.
b) row F: log2 m* 5.4e3 -> 3.1e3 (density 4.26, brute).  c) high: the max over ALL K^2 pairs must be balanced for every
piece assignment iota; no lemma.  d) not done.  e) exhaustive test: for random iota-loads on the g-node row network,
compute the best partition's max pair load for g/K = 4/zeta vs Hoeffding (small g, K).

### I6. Explicit-constant OpenAI App. D  [N, model; constants from G18a/G18d]
Crossover log2 s* = beta (gamma c w / tau)^2 (u0 = log2 s/beta, c bits per M/sqrt u0, gamma depth sum, w wires/var):
   density 4.26, vs brute:  A paper (eps = 2^-17, c = 38.5 bit, gamma 1.215) 4.0e10  [Dias–Lima via G4c: 5.3e11]
                            B eps at its legal max under (D.6) (u0 = log2 s/960; [H] re-check (D.5))   2.9e8
                            C + LP taper, full packets nu = 2 (c = 12.6 bit)                           3.2e7
                            D + Smith/Mahler modes (c/2, unverified as algorithm)                        7.9e6
                            E + r = 32, Delta = s^{1/2} (beta = 160, gamma = 1.36)                       1.6e6
                            F idealised (K1) beta = 11                                                    1.1e5
   vs PPSZ: x 6.7.  b) up to 2.5e4x in log2 s*, but still ≥ 300x worse than I2–I4 on Claude's side.
c) B: low (pure parameter; ε = 2^-17 is ~136x below what (D.6) needs); D: medium (Smith basis inside the evaluator);
E: medium (Δ smaller needs re-checking the leaf/entry costs). d) G18a/G18d (κ, taper, Smith 4x), G20e (≈1e6), G4c.
e) Re-derive (D.5)/(D.6) with ε = 1/960 symbolically (a ≤ 4ελ ≤ λ/240, r^{4a} s^{0.9} ≤ H) and recompute row B.

### I7. Combine with PPSZ / branching at small n  [P/N, bookkeeping]
a) Output min(PPSZ, algorithm) — there is no interleaving gain (G20e 2c: branch-then-algebra increasing in b).
b) claude_curve.py, re-optimised Claude (rows C–E), exponent/n:
   log2 m:     40     100    400    1000   5000   1e4    3e4    1e5
   r = 4.26:   5.89   5.34   2.99   2.11   1.06   0.75   0.38   0.21
   r = 1:      1.00   0.91   0.51   0.36   0.18   0.13   0.065  0.036
 At every realistic size (log2 m ≤ 40, i.e. m ≤ 10^12) the best the recursion can do is level 0 = the sweep, exponent
 W/2 ≈ 5.9 n at density 4.26: six times worse than brute force, ~15x worse than PPSZ (0.386 n). Paper as written cannot
 run at all below log2 m ≈ 162 (m ≥ K0^9).  Crossovers (density 4.26): vs 2^n 5.4e3, vs #3SAT 1.642^n [Kutzkov, from
 memory] 1.1e4, vs PPSZ 2.9e4 (paper: 2.8e5, 5.5e5, 1.1e6).  c) none.  d) partly G20e (L5).  e) done here.

### I8. Make c (App. B pieces) and C0 explicit -> a fully explicit Theorem A'  [buildable]
a) Def. 7(c) constant c of App. B (piece count 4h+1 ≈ 0.3 lambda^2, coefficient bit sizes) from G20b's appb_exact.py
(h ≤ 7) + the proof; C4, C5 from the operation counts. b) turns rows C–E from a model into a theorem with numbers;
moves b_top by ~10–30 only (constants enter as log). c) low/medium (tedious). d) no.
e) Measure max |coefficient| and |I| in appb_exact.py for h = 2..7, fit c; plug into claude_cross2.py.

### I9. Fewer wires per variable  [H]
a) log m* ∝ (W/n)^2. W = 3m − n (chain) is already a spanning-forest minimum when a variable's occurrences fall in
distinct groups (K >> degree), so for random 3-CNF at 4.26 no gain; for structured inputs (variables concentrated in few
clusters) W can drop a lot. b) W/n 11.8 -> 8 would give log m* / 2.2. c) high for hard instances. d) G20e (S >= phi M).
e) measure distinct-group degree of variables after a min-cut clustering on industrial CNFs.

### I10. d-dimensional grids (D directions per level)  [H]
Same α' = 1/2 by the law ρ_D = D^{-1/2}; overheads per level amortised over D, but no D > 2 pieces are known. Low value.

### I11. Mod-3 evaluator constant 12.6 -> 8.9 (G18r)  — interpolation route only, dominated by the headline. Low value.

## 2. TOP 3 (value x probability / cost)

1. **I2 + I3: explicit, re-optimised Theorem A with the cut-width base** (cheap, safe, factor ~6 in log2 m* at
   density 4.26: 1.37e5 -> 2.4e4 vs brute force; 1.1e6 -> 1.3e5 vs PPSZ). First computation: write the sweep lemma and
   the exact-Hoeffding schedule as a recurrence; rerun claude_cross2.py rows C/D with the explicit c of I8 (1 h).
2. **I4: top level m ≈ K^2 by random clustering** (another factor 4.5: -> 5.4e3 vs brute force, 2.9e4 vs PPSZ).
   First computation: state Lemma 8' with explicit probability, check the pair-load concentration numerically on the
   chain network of random 3-CNF (n = 10^5, g = 100..1000 clusters: max c_ab / mean) (< 1 h).
3. **I6-B/C: OpenAI App. D with legal-max ε and LP taper** (factor ~1000 in log2 s* for free, 4e10 -> 3e7), mainly to
   fix the record and to document that Claude's architecture is ~20x better in the finite-size constant.
   First computation: symbolic re-check of (D.5)/(D.6) at ε = 1/960 (30 min).

Honest bottom line: even with everything above (incl. speculative I5) the lower rung beats brute force only for
log2 m ≳ 3·10^3 (density 4.26) and PPSZ only for log2 m ≳ 1.7·10^4. At m ≤ 10^12 the recursion is ≥ 6x worse than
brute force in the exponent. "Best pieces at each size" is irrelevant because the algorithm's crossings have
lambda ≥ 10^4 bits, where App. B is already at its asymptotic circle.

## 3. Self-review
* The Claude model uses the paper's recurrences (40)/(41)/(45) and reproduces its own r = 1 remark (8.05e3) exactly;
  rows D–F change assumptions (e0, a_top, g/K) that are argued, not proved; C0 and c unknown (enter only b_top).
* The tau targets compare a counting algorithm with decision PPSZ; #3SAT 1.642^n is quoted from memory (as in G9b).
* OpenAI rows use G18d's constants (LP 6.2 sqrt(nu) nats, Smith factor 1/2 unverified); beta = 960 is my reading of
  (D.5)/(D.6), not re-derived line by line. Row F is an idealisation (K1).
* No process left running; nothing written into the repository.
