# G18i — Refuting random 3-CNF at constant density after the October 2026 algorithms

Status: working notes (level P = paper argument, unreviewed; E = experiment). Everything that
uses the two 2026 preprints is CONDITIONAL on them (unrefereed).
Sources read: JOURNAL.md, code/spec/requirements.md, Claude preprint (claude_3sat.txt),
OpenAI preprint (g14b/dl/paper.txt; src_paper.tex is a 404 page), arXiv abstracts in dl/.

## 0. Intermediate facts checked (saved early)

- RRS (Raghavendra–Rao–Schramm, STOC 2017, arXiv:1605.00058, abstract): strong refutation of
  random k-XOR (and k-SAT) at m/n >= Õ(n^{(k/2-1)(1-δ)}) in time exp(Õ(n^δ)), δ in [0,1).
  k=3: m >= Õ(n^{3/2-δ/2}). The Õ hides polylog(n): at m = Δn with Δ constant it gives nothing.
- KMOW (Kothari–Mori–O'Donnell–Witmer, STOC 2017, arXiv:1701.04521, abstract): degree-d SOS with
  d = Θ(n/(Δ^{2/(t-1)} log Δ)) cannot refute; 3-SAT has t = 2, so d = Θ(n/(Δ² log Δ)).
- FKO (Feige–Kim–Ofek, FOCS 2006): poly-size witnesses whp for m > c n^{7/5}; witness found in
  time 2^{O(n^{0.2} log n)} (web summary of the preprint; exponent 0.2 to be re-read in the PDF).
- Müller–Tzameret (arXiv:1101.3970; APAL 2014): poly-size TC^0-Frege refutations at m = Ω(n^{1.4});
  abstract says the only earlier non-trivial upper bound in a non-abstract system was resolution
  at m = Ω(n²/log n) (Beame–Karp–Pitassi–Saks 2002).
- GKM (arXiv:2109.04415), HKM (arXiv:2207.10850): n^{O(ℓ)} time at m >= Õ(n)(n/ℓ)^{k/2-1},
  one log factor in HKM; smoothed instances too; FKO-type witnesses for smoothed instances.
- OpenAI Thm 1.1 / App. D.2: network size s = n + m, exponent O_k(s/√log s); so for m <= Dn the
  constant is LINEAR in (1+D).  Claude Thm 1: 2^{O(min{k² n loglog n/√log n, k m/√log m})}.

Legend: [T] published theorem; [T26] theorem of an unrefereed 2026 preprint; [P] my paper argument,
unreviewed; [E] experiment; [M] cited from memory, statement not re-read today; [?] unproved.

## 1. The map: time versus density for refuting random 3-CNF (n variables, m = Δn clauses)

"Refute" = deterministic algorithm, never says UNSAT on a satisfiable input, says UNSAT whp over F.
Δ0 := any density at which random 3-CNF is UNSAT whp (5.191 by the first moment; 4.4898 rigorous,
Díaz–Kirousis–Mitsche–Pérez-Giménez 2009 [M]; conjectured 4.267).

### 1.1 What the 2026 algorithms give
- Worst case [T26]: OpenAI Thm 1.1 — 2^{O_k(s/√log s)}, s = n+m (App. D.2), so the constant for
  m <= Dn is linear in 1+D; any m: 2^{O_k(n loglog n/√log n)}.  Claude Thm 1 —
  2^{O(min{k² n loglog n/√log n, k m/√log m})}.  Both count; both are purely asymptotic
  (Claude: m = 5n needs n > 10^{10 000}; "Limits" note: log2(n+m) > 10^9 for the OpenAI proof).
- Both are SIZE-OBLIVIOUS: the cost bound depends only on (n, m). Claude §2.2 pads every instance
  to K nodes pairwise joined by λ-bit edges ("only the largest number of bits on an edge matters");
  OpenAI Thm D.1 depends on the wire count only. So on a random formula the bound is exactly the
  worst-case bound: random is neither easier nor harder; no property of randomness is used.
  (Random formulas are expanders, i.e. the generic hard case for every cut-based step: Claude §2.2.)
- **Proposition A [P, conditional on T26] (density-dependent speed-up).** For Δ0 <= Δ <= n^{1-ε},
  random 3-CNF with Δn clauses is refuted whp in time 2^{O(n/√(Δ log n))}.
  Proof. Let n' = ⌈n·√(Δ0(1+ε)/Δ)⌉ and F' = clauses of F inside x_1..x_n'. Each clause lands in F'
  independently with prob. C(n',3)/C(n,3) ~ (n'/n)³ and is then uniform on the n' variables, so F'
  is a random 3-CNF on n' variables with ~Δ0(1+ε)n' clauses, UNSAT whp (n' -> ∞). Run the 2026
  counter on F': m' = Θ(n/√Δ), log m' = Θ(log n). F' ⊆ F, so soundness is trivial. ∎
  Without sub-sampling the bounds are 2^{O(Δ n/√log n)} (m-bound, grows with Δ) or
  2^{O(n loglog n/√log n)} (n-bound). The trick is useless with pre-2026 2^{cn} algorithms
  (2^{cn/√Δ} loses to resolution's 2^{O(n/Δ)}), which is why it is only worth stating now.
- Sign-blind variant [P]: choose S from the hypergraph only (signs unseen) with >= 5.2|S| clauses
  inside; then P(F|S satisfiable) <= 2^{|S|}(7/8)^{5.2|S|} -> 0. Any poly-time densest-subset
  heuristic can be plugged in (tested in §3).

### 1.2 Frontier table (★ = a 2026 algorithm is the best known, conditional on T26)

| density m/n = Δ | best known refutation time | source | lower bounds for restricted systems |
|---|---|---|---|
| Δ0 <= Δ = O(1) | ★ 2^{O(n/√(Δ log n))} (was 2^{O(n/Δ)}, i.e. 2^{Θ_Δ(n)}) | Prop. A on T26; before: Beame–Karp–Pitassi–Saks, SICOMP 31 (2002) [M] | resolution exp(Ω(n/Δ^{2+ε})) (Chvátal–Szemerédi JACM 1988; Ben-Sasson–Wigderson JACM 2001); PC degree Ω(n) every field (Ben-Sasson–Impagliazzo FOCS 1999; Alekhnovich–Razborov FOCS 2001) hence size 2^{Ω(n)} (Impagliazzo–Pudlák–Sgall 1999); SOS degree Ω(n/(Δ² log Δ)) (Grigoriev 2001; Schoenebeck FOCS 2008; KMOW STOC 2017); Res(k), k <= √(log n/loglog n) (Alekhnovich 2005 [M]). Open [M]: cutting planes and bounded-depth Frege on random 3-CNF; no bound for AC0[p]-Frege, TC0-Frege, Frege, EF, IPS. |
| ω(1) <= Δ <= c·log n | ★ 2^{O(n/√(Δ log n))} | Prop. A | same |
| C log n <= Δ <= ~log³ n | 2^{O(n/Δ)} (tree-like resolution, ordered DLL) | BKPS 2002 [M] | resolution exp(Ω(n/Δ^{2+ε})) |
| ~log³ n <= Δ < n^{1/2}, Δ = n^β | exp(Õ(n/Δ²)) = exp(Õ(n^{1-2β})); exp(O(n log³ n/Δ²)) with one log | RRS STOC 2017 (arXiv:1605.00058); GKM (2109.04415); HKM (2207.10850) | SOS degree Ω̃(n/Δ²): tight (KMOW) |
| Δ >= c·n^{0.4} | poly-size NP witness, found in 2^{O(n^{0.2} log n)}; poly-size TC0-Frege proof | Feige–Kim–Ofek FOCS 2006; Müller–Tzameret APAL 2014 (arXiv:1101.3970) | resolution still exponential up to n^{1/2-ε} |
| Δ >= C·n^{1/2} | polynomial (spectral); strong refutation at Õ(n^{1/2}) | Feige–Ofek ICALP 2004 / ToC 2007; Coja-Oghlan–Goerdt–Lanka 2007 [M]; Allen–O'Donnell–Witmer FOCS 2015 (1505.04383) | — |
| Δ >= n/log n | poly-size resolution | BKPS 2002 (quoted in Müller–Tzameret abstract) | — |

Caveats. (i) The crossovers "c log n" and "log³ n" compare O(·) bounds with unknown constants; with
the constants of the 2026 proofs the first crossover is at astronomically large n. (ii) Without Prop. A
the ★ window is only Δ <~ √log n/loglog n. (iii) For m = n^{1+β}, β > 0, the 2026 algorithms are
never the best (2^{n^{1-β/2}/√log n} against RRS exp(Õ(n^{1-2β})) and BKPS 2^{O(n^{1-β})}).
(iv) Feige's hypothesis (poly time, Δ large constant) is untouched; what falls (conditionally) is only
its exponential strengthening "2^{Ω(n)} is needed at constant density".

## 2. Proof-complexity reading

- Trivial fact (Cook–Reckhow 1979): "the run of A on F" is a proof of length T in the ad hoc system
  "A-runs". Not natural, no information.
- **Claude's algorithm -> IPS over Q** (Grochow–Pitassi, arXiv:1404.3820, JACM 2018). [P; same route
  as G4d] Let V_C(x) be the product of the 3 falsified-literal polynomials of clause C (axiom V_C = 0),
  P = Π_C (1 - V_C), and Q_k(x) = Σ_{y_1..y_k ∈ {0,1}} P(y_1..y_k, x_{k+1}..x_n) Π_{i<=k}(x_i y_i + (1-x_i)(1-y_i)).
  Q_0 = P, Q_n = multilinearisation of P, and Q_n ≡ 0 iff F is unsatisfiable. Then
  1 = Σ_j V_{C_j} Π_{l<j}(1 - V_{C_l}) + Σ_k (x_k² - x_k)·D_k + Q_n, with D_k = (Q_{k-1} - Q_k)/(x_k² - x_k)
  (exact division; interpolation in x_k, degree <= m). Each Q_k is a weighted count of the clause
  network with tables over Q[x]; Claude's recursion takes arbitrary tables W_a (§2.2), so it gives a
  circuit of size 2^{O(km/√log m)} for each Q_k. The certificate is linear in the axiom placeholders
  (Hilbert-like IPS).
  [?] needs: the recursion is a table-oblivious arithmetic circuit over Q (no branching on table
  values; division only by fixed integers 2, h!). Believed from the text, not checked line by line.
- **OpenAI's algorithm**: arithmetic in (Z/2^{P0})[V]/(V^ℓ) with 2-adic truncation (Thm D.1). Natural
  home: IPS over the ring Z/2^{n+1} (IPS is defined over any ring), since the output is #F only modulo
  2^{n+1} and a certificate over Q would carry a 2^{n+1}·(garbage) term. [?] Thm D.1 is stated for one
  truncated variable V, not for n free variables x; whether the truncation argument survives over
  (Z/2^P)[x_1..x_n] is not checked. I do not claim a proof system for it.
- **Extended Frege / Frege: nothing follows.** IPS p-simulates EF, and the converse is known only if
  EF proves the PIT axioms (Grochow–Pitassi). A size-poly(T) EF proof needs a PIT-free derivation of
  "Q_n ≡ 0", i.e. a correctness proof of the algorithm that never handles an object of size 2^n
  (2^n = T^{√log n} is super-polynomial in T). The invariant "sub-call = full contraction" is a sum of
  2^{#wires} terms; G8c reduced the question to closed forms of pinned sub-calls (level 1 only) and
  hypothesis (S). Still open. Best EF/Frege bound at density Δ: 2^{O(n/Δ)}, inherited from resolution.
- **Cleanest true statement (conditional on Claude's Thm 1 + obliviousness [?]):**
  every unsatisfiable k-CNF with m clauses has a Hilbert-like IPS_Q refutation of circuit size
  2^{O(km/√log m)}; for random 3-CNF with Δn clauses (Δ0 <= Δ <= n^{1-ε}), whp size
  2^{O(n/√(Δ log n))} (IPS refutations are closed under adding axioms; use F' of Prop. A; more generally
  size 2^{O(μ/√log μ)} where μ = least number of clauses of an unsatisfiable sub-formula).
  For constant Δ this is, to my knowledge (Müller–Tzameret's abstract for the state in 2011; nothing
  found since), the first 2^{o(n)} upper bound for random 3-CNF at constant density in a named
  propositional proof system.
- **Beside, not against, the lower bounds.** Resolution, Res(k), NS, PC/PCR, SA, SOS need 2^{Ω_Δ(n)};
  none of them simulates IPS circuits (PC is p-simulated by IPS; its lower bound is a degree/monomial
  bound, and the certificate above is a deep circuit, not a sparse polynomial). So: a 2^{o(n)} versus
  2^{Ω(n)} average-case separation of IPS from PC, resolution and SOS at constant density.
  It is quantitatively weak (o(n) vs Ω(n) in the exponent) and not new as a separation of systems
  (EF vs PC on PHP is poly vs exp); new is only the formula family. IPS is verified by PIT (coRP), so
  it is a Cook–Reckhow system only if PIT is derandomised.

## 3. The tested idea: "core first" for size-oblivious refuters [E]

Idea. A size-oblivious refuter (cost a function of m' only) should be fed the smallest unsatisfiable
sub-formula, not F. This is the opposite of resolution, which profits from every extra clause.
Question: how small is an unsatisfiable part of a random 3-CNF at density Δ, and what finds it?
Bounds [P]: any unsatisfiable sub-formula has more clauses than variables, and whp every set of
v <= n/(e²Δ) variables contains fewer than v clauses (first moment: (e²Δv/n)^v), so
μ >= n/(e²Δ); a set of xn variables with induced density >= ρ needs x^{2ρ-1} >= (ρ/(eΔ))^ρ/e.
Upper: μ <= Δ0^{3/2}(1+ε)·n/√Δ (Prop. A).  True order of μ between n/Δ and n/√Δ: open [?].

Code: core_exp.py, sound_refuter.py, analyse.py; raw data res_*.jsonl; table analysis.txt.
Generator: 3 distinct variables, uniform signs; cadical 1.5.3 / glucose 4 through pysat.
(a) prefix core: least p with F|{x_1..x_p} UNSAT (bisection; monotone). (b) peel core: same along the
min-degree peeling order (depends on the hypergraph only). (c) MUS: deletion-based clause-minimal
core inside (b), 40 s budget (an upper bound on μ, not the minimum).
Medians (10 formulas per row; 6 for n = 800):

| n | Δ | prefix n'/n | √(4.27/Δ) | peel n'/n | MUS vars/n | MUS clauses/n | n/(e²Δ) bound | CDCL conflicts: F / prefix / peel / MUS |
|---|---|---|---|---|---|---|---|---|
| 200 | 16 | 0.522 | 0.517 | 0.355 | 0.345 | 0.882 | 0.0085 | 86 / 934 / 256 / 534 |
| 400 | 32 | 0.365 | 0.365 | 0.245 | 0.241 | 0.709 | 0.0042 | 76 / 4090 / 832 / 2021 |
| 400 | 64 | 0.260 | 0.258 | 0.165 | 0.160 | 0.404 | 0.0021 | 22 / 848 / 229 / 402 |
| 400 | 128 | 0.184 | 0.183 | 0.114 | 0.109 | 0.249 | 0.0011 | 12 / 227 / 69 / 96 |
| 800 | 128 | 0.185 | 0.183 | 0.106 | 0.106 | 0.292 | 0.0011 | 20 / 3566 / 678 / 1559 |
| 800 | 256 | 0.134 | 0.129 | 0.073 | 0.071 | 0.175 | 0.0005 | 10 / 655 / 140 / 213 |

Log-log slopes in Δ (n = 400, Δ = 24..128): prefix variables -0.50, peel variables -0.55,
MUS variables -0.57, MUS clauses -0.77 (n = 800: -0.49, -0.55, -0.56, -0.81).

Sign-blind refuter (density >= 5.2 inside S, 5 formulas per row, sound_refuter.py):

| n | Δ | prefix |S|/n | √(5.2/Δ) | peel |S|/n | peel/prefix | UNSAT confirmed |
|---|---|---|---|---|---|---|
| 400 | 32 | 0.4075 | 0.4031 | 0.275 | 0.675 | 5/5 both |
| 400 | 128 | 0.205 | 0.2016 | 0.1275 | 0.622 | 5/5 both |
| 1600 | 128 | 0.2037 | 0.2016 | 0.1231 | 0.604 | 5/5 peel (prefix not run: too large for CDCL) |
| 1600 | 256 | 0.1437 | 0.1425 | 0.085 | 0.591 | 5/5 both |
| 1600 | 512 | 0.1019 | 0.1008 | 0.060 | 0.589 | 5/5 both |

(n = 1600, Δ = 64: sizes 0.2863 / 0.1794, UNSAT not checked. First-moment floor for density 5.2 at
Δ = 512: 0.041.)

Reading.
1. Prop. A is confirmed numerically: the prefix core has n√(Δ0/Δ) variables (slope -0.50).
2. Min-degree peeling (poly time, sign-blind) finds a set of 0.54–0.68 times the prefix size, with induced
   density 4.2–4.6 where a random set of that size has density ≈ 1.7. The ratio drifts from 0.675 to
   0.589 over Δ = 32..512: a constant-factor gain, or an exponent slightly below -1/2 (fitted -0.55);
   the data cannot separate the two. The first-moment floor has exponent -0.553 at ρ = 5.2.
3. MUS: 2–3 clauses per variable, clause count 2–3 times below the peel core, apparent slope -0.8,
   BUT MUS clauses/n grows with n at fixed Δ (Δ = 64: 0.325, 0.345, 0.404 at n = 100, 200, 400;
   Δ = 128: 0.240, 0.210, 0.249, 0.292 at n = 100..800), and some small-Δ rows hit the time budget.
   No conclusion on the order of μ. The bound n/(e²Δ) is 100–350 times below what was found.
4. Contrast: CDCL (resolution) needs 10–90 conflicts on F but 3 to 180 times more on any core.
   Resolution uses density as 2^{O(n/Δ)}; a size-oblivious refuter can only use it as 1/√Δ unless
   unsatisfiable cores with o(n/√Δ) clauses exist and can be found.
5. At every size tested, plain CDCL on F is far better than anything the 2026 bounds certify.

Why no structured core can help a separator-based method [P, from T]: a tree decomposition of width w
of an unsatisfiable sub-formula gives a resolution refutation of width O(w); F needs width
Ω(n/Δ^{1+ε}) (Ben-Sasson–Wigderson 2001), and so does every sub-formula. So every unsatisfiable part
of F has treewidth Ω(n/Δ^{1+ε}): choosing the core cannot create small separators.

## 4. Verdict

- Does anything here move the central obstacle? **No.**
  (a) The 2026 algorithms are worst-case and size-oblivious; on random formulas they use nothing
      random, and the only average-case lever (sub-sampling) gains a factor √Δ in the exponent.
  (b) Their saving is a √log factor in the exponent; inside their own framework the floor is about
      2^{m/log m} (ρ -> 1/2; Claude App. B, and G4b/G4c/G5a for the models studied). The obstacle asks
      for polynomial time; even 2^{n^{1-ε}} at constant density is out of reach (RRS gets there only
      for Δ >= n^{ε/2}).
  (c) Their certificates are IPS circuits checked by PIT, a system that never had lower bounds; they
      say nothing about EF or Frege, where the question is exactly the G8c one.
- What does change (conditional): "2^{Ω(n)} at constant density" is no longer a property of every
  known method; random 3-CNF at constant density joins the list of average-case separations of IPS
  from PC/resolution/SOS; the refutation frontier gets a new segment Δ <= c log n.
- Open questions worth a follow-up: (1) order of μ(Δ), the least number of clauses of an unsatisfiable
  sub-formula (between n/Δ and n/√Δ); if μ = O(n/Δ) and findable, core-first would beat resolution's
  2^{O(n/Δ)} at every density by the √log factor. (2) table-obliviousness of both recursions (needed
  for the IPS statements). (3) G8c: a PIT-free proof of Q_n ≡ 0 of size poly(T).

## 5. Adversarial self-review

- Everything marked T26 rests on two unrefereed preprints; I did not re-verify them.
- [M] items (BKPS exponent n/Δ, Res(k) range, openness of CP and bounded-depth Frege on random
  3-CNF, 4.4898) are from memory; arXiv abstracts were read only for RRS, KMOW, GKM, HKM, AOW,
  Müller–Tzameret, Grochow–Pitassi; FKO's 2^{O(n^{0.2} log n)} comes from a web summary.
- "First 2^{o(n)} upper bound in a named system at constant density": negative literature claim, one
  web search plus the 2011 abstract; G4d made the same claim. Could be wrong.
- Prop. A is a two-line corollary; I found no statement of it, but did not search for one specifically.
- Experiment: n <= 800 (1600 for the sign-blind part), 5–20 formulas per point, medians only; the
  generator fixes m exactly; MUS values are local minima and drift with n; fitted slopes are over
  less than a decade and a half of Δ. The "-0.55" is not evidence for an exponent different from -1/2.
- The IPS statement for OpenAI's algorithm is NOT established.
