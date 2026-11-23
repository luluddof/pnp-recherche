# Step 2 — Requirements specification for a polynomial-time 3-SAT algorithm

Purpose: list every algorithm family already **proven** unable to decide 3-SAT in
polynomial time, state the exact scope of each theorem, and derive what a
successful candidate must look like. Written from memory (no web access).
Items marked **[unverified]** are details (exact exponents, years, authors) I am
not certain of; the qualitative statements are standard.

Honesty note: the overwhelming expert consensus is P != NP. Nothing below is a
lower bound against *all* algorithms — no such bound is known (that would be a
proof of P != NP). Every result is a lower bound against a *restricted model*.
The useful content is: "an algorithm of shape X is provably dead; don't retry X".

---

## 0. The master reduction: algorithms -> proof systems (Cook–Reckhow 1979)

* A decider A for 3-SAT that answers UNSAT in time T(n) defines a proof system
  whose refutation of F is "the run of A on F" (length T(n)). If A is poly-time,
  this system is **polynomially bounded**, which implies NP = coNP
  (Cook–Reckhow 1979). So every proof-complexity lower bound below says:
  *the certificates implicitly produced by algorithm family X are proofs in
  system S, and S has no short proofs for family H.*
* Cook 1975 (PV correspondence): if the correctness of A ("A(F)=UNSAT implies F
  unsatisfiable") is provable in the bounded-arithmetic theory PV (or S^1_2),
  then **Extended Frege (EF)** has polynomial-size refutations of all
  unsatisfiable CNFs. Contrapositive: any superpolynomial EF lower bound would
  mean no algorithm with a PV-provable soundness proof works. (No EF lower bound
  is known.)
* Consequence for this project: the earlier finding that **refutation is the
  bottleneck** is not an accident — it is exactly the place where all known
  lower bounds live.

---

## 1. Excluded families (theorem, hard instances, exact scope)

### 1.1 Tree-like resolution = DPLL (incl. lookahead solvers)
* **Theorems.** Every DPLL run on an UNSAT formula is a tree-like resolution
  refutation (folklore; Davis–Logemann–Loveland 1962). Exponential lower bounds:
  Haken 1985 (PHP, general resolution, hence tree-like too); random 3-CNF:
  Chvátal–Szemerédi 1988 (general resolution, see 1.2). Myopic/"drunk" DPLL
  takes exponential time even on **satisfiable** formulas:
  Alekhnovich–Hirsch–Itsykson 2005.
* **Hard instances.** Random 3-CNF at any constant density above threshold
  (incl. 4.26 restricted to UNSAT ones), PHP, Tseitin on expanders.
* **Scope.** Size of tree-like resolution refutations. Applies to every
  branching algorithm whose leaf contradictions are falsified clauses
  (lookahead, failed-literal probing, clause-level propagation are all inside
  resolution).

### 1.2 General resolution (incl. CDCL with restarts)
* **Theorems.** Haken 1985 (PHP_n^{n+1}: 2^{Ω(n)}); Urquhart 1987 (Tseitin on
  expanders); Chvátal–Szemerédi 1988 (random k-CNF, m = cn, c above a constant:
  size 2^{Ω(n)}); Beame–Pitassi 1996 (simplified); Ben-Sasson–Wigderson 2001
  (size >= exp(Ω((w - k)^2 / n)), width w; random 3-CNF width Ω(n / Δ^{1/(k-2)+ε})
  for m = Δn [exponent unverified]); Beame–Karp–Pitassi–Saks 2002 (random 3-CNF
  remains resolution-hard up to m ≈ n^{3/2 - ε} [exact range unverified]).
* **CDCL = resolution.** CDCL with restarts p-simulates resolution and every
  CDCL run is a resolution proof (Pipatsrisawat–Darwiche 2011;
  Atserias–Fichte–Thurley 2011). So no CDCL heuristic, restart policy, learning
  scheme or preprocessing *that stays inside resolution* can be polynomial on
  random UNSAT 3-CNF.
* **Scope.** Size of resolution refutations. Does **not** cover extension
  variables (ER), parity reasoning, cardinality reasoning, symmetry breaking
  (not sound-for-refutation in resolution), or PR/DRAT-with-new-variables.

### 1.3 Bounded width / k-consistency / local lookahead
* **Theorems.** k-consistency (existential k-pebble game) = width-k resolution
  (Atserias–Dalmau 2008). Random 3-CNF at constant density needs width Ω(n)
  (Ben-Sasson–Wigderson 2001). For CSPs: bounded-width algorithms solve exactly
  the CSPs without "ability to count" (Barto–Kozik 2014; Atserias–Bulatov–Dawar
  2009 show k-consistency fails on linear equations over Z_2 — even though they
  are in P via Gaussian elimination).
* **Scope.** Any algorithm that only derives consequences supported on <= k
  variables. Note: Z_2 equations show *bounded width fails on problems in P* —
  a model limitation, not a problem limitation.

### 1.4 Nullstellensatz and Polynomial Calculus (incl. Gaussian elimination)
* **Theorems.** Nullstellensatz degree lower bounds: Beame–Impagliazzo–
  Krajíček–Pitassi–Pudlák 1994. PC degree: Razborov 1998 (PHP);
  Impagliazzo–Pudlák–Sgall 1999 (degree-to-size: size >= exp(Ω((d-k)^2/n)));
  Buss–Grigoriev–Impagliazzo–Pitassi 2001 (Tseitin mod p, char != 2);
  Ben-Sasson–Impagliazzo 1999/2010 and **Alekhnovich–Razborov 2001/2003**:
  random 3-CNF requires PC degree Ω(n) over **every field** (hence size
  2^{Ω(n)}).
* **Hard instances.** Random 3-CNF (any field); Tseitin (fields of char != 2);
  PHP.
* **Scope.** PC/NS over a field, with Boolean axioms x^2 = x. Gaussian
  elimination over GF(2) is degree-1 PC over GF(2) plus linear combination;
  Res(lin)/Res(⊕) is a different system (see §3). PC over **rings that are not
  fields** (Z_m composite, Z) is less covered (partial results exist for NS
  over Z_m in Buss–Impagliazzo–Krajíček–Pudlák–Razborov–Sgall 1996
  [unverified for random 3-CNF]).

### 1.5 Cutting Planes (incl. Pudlák interpolation)
* **Theorems.** Pudlák 1997: feasible interpolation for CP + exponential
  lower bounds for monotone **real** circuits computing clique ⇒ clique–coloring
  formulas need 2^{n^{Ω(1)}} CP size (also Haken–Cook 1999). Rank: CP
  (Chvátal–Gomory) rank of random 3-CNF is Ω(n) (Buresh-Oppenheim–Galesi–Hoory–
  Magen–Pitassi 2003 [authors unverified]). Random **Θ(log n)-CNF** need
  exponential CP size: Fleming–Pankratov–Pitassi–Robere 2017; Hrubeš–Pudlák
  2017 (both via interpolation-style arguments).
* **Open.** Superpolynomial CP size lower bounds for **random 3-CNF** at
  constant density are, to my knowledge, **open** [unverified as of 2026].
  Tseitin is *easy* for CP (Dadush–Tiwari 2020: quasi-poly CP proofs).
* **Scope.** CP proof size (and rank). Covers ILP branch-and-cut with
  CG/Gomory cuts, pseudo-Boolean solvers (RoundingSat etc. are CP-based).

### 1.6 Sherali–Adams, Lovász–Schrijver, SOS / Lasserre
* **Theorems.** Grigoriev 2001 (Positivstellensatz/SOS degree Ω(n) for Tseitin
  and random 3-XOR); **Schoenebeck 2008** (SOS/Lasserre degree Ω(n) for random
  3-XOR and 3-SAT at m = O(n)); Alekhnovich–Arora–Tourlakis 2005 (LS+ rank);
  Buresh-Oppenheim et al. 2003 (LS rank); **Kothari–Mori–O'Donnell–Witmer 2017**
  (degree-d SOS cannot refute random 3-SAT unless m >~ n^{3/2}/d^{1/2}, up to
  logs); matching upper bound Raghavendra–Rao–Schramm 2017. Size: SOS
  size–degree tradeoff, Atserias–Hakoniemi 2019 [unverified] ⇒ exponential SOS
  **size** on random 3-CNF.
* **Optimality of SA/SOS among LPs/SDPs (approximation setting).**
  Chan–Lee–Raghavendra–Steurer 2013/16: poly-size LP relaxations of max-CSP are
  no stronger than O(1)-round Sherali–Adams; Kothari–Meka–Raghavendra 2017:
  sub-exponential size LP needs n^{Ω(...)} rounds [details unverified];
  Lee–Raghavendra–Steurer 2015: same for SDPs vs SOS.
* **Scope.** Degree / rank of the hierarchy (and LP/SDP size for *symmetric/
  CSP-structured* relaxations). Spectral refutation (degree-2 SOS) is included:
  matches the earlier observation m >~ 7.4·n^{1.5}.

### 1.7 Extended formulations (Fiorini et al. 2012, Rothvoss 2014)
* **Theorems.** Yannakakis 1991 (symmetric EFs of TSP/matching are exponential);
  Fiorini–Massar–Pokutta–Tiwary–de Wolf 2012 (cut, TSP, stable-set polytopes:
  extension complexity 2^{Ω(√n)}); **Rothvoss 2014**: perfect matching polytope
  2^{Ω(n)}. Reductions give superpolynomial EFs for "3-SAT polytopes" (Avis–
  Tiwary 2015 [unverified]).
* **Scope — critical.** Bounds the number of facets of *any* polytope projecting
  onto *one specific polytope*. **Perfect matching is in P (Edmonds 1965) yet
  has no small LP.** So these theorems exclude "write 3-SAT as one fixed
  polynomial LP" but say **nothing** about algorithms in general. Not a barrier
  for a combinatorial or adaptive algorithm.

### 1.8 Monotone circuits (Razborov)
* **Theorems.** Razborov 1985 (CLIQUE superpolynomial), Andreev 1985,
  Alon–Boppana 1987 (exponential); Razborov 1985 **perfect matching needs
  n^{Ω(log n)} monotone circuits** (in P!); Tardos 1988 (a monotone function in
  P needs exponential monotone circuits). UNSAT as a function of the clause-
  indicator vector is monotone; monotone CSP-SAT lower bounds:
  Göös–Kamath–Robere–Sokolov 2019 [unverified details].
* **Scope.** Circuits without negation. Again a model restriction: problems in P
  are excluded too.

### 1.9 Bounded-depth Frege / AC^0 / AC^0[p]
* **Theorems.** Circuits: Ajtai 1983, Furst–Saxe–Sipser 1984, Håstad 1986
  (parity ∉ AC^0); Razborov 1987, Smolensky 1987 (MOD_q ∉ AC^0[p]). Since 3-SAT
  is NP-complete under projections, 3-SAT ∉ AC^0, AC^0[p]. Proofs: Ajtai 1988,
  Pitassi–Beame–Impagliazzo 1993, Krajíček–Pudlák–Woods 1995 (PHP exponential
  for depth-d Frege); Ben-Sasson 2002, Pitassi–Rossman–Servedio–Tan 2016
  (Tseitin). Random 3-CNF in bounded-depth Frege: **open** [unverified].
  AC^0[p]-Frege: no superpolynomial lower bounds known.
* **Scope.** Constant-depth circuits / constant-depth line proofs.

### 1.10 OBDDs and knowledge compilation
* **Theorems.** Bryant 1991 (multiplication middle bit: exponential OBDD for
  every order). OBDD proof system (Atserias–Kolaitis–Vardi 2004): PHP and
  others hard (Krajíček 2008; Tveretina–Sinz–Zantema 2010); random CNF need
  exponential OBDD refutations (Friedman–Xu 2013 [unverified]). d-DNNF/
  decision-DNNF exponential lower bounds (Beame–Li–Roy–Suciu 2013;
  Bova–Capelli–Mengel–Slivovsky 2016); decision-DNNF compilers trace
  resolution-like proofs.
* **Scope.** Size of the compiled representation (fixed order / syntactic
  restrictions). Unrestricted Boolean circuits = no restriction (P/poly).

### 1.11 Local search, stable and sequential algorithms (search side)
* **Theorems.** Overlap gap property (OGP): Gamarnik–Sudan 2014/2017 (local
  algorithms on sparse random graphs; balanced sequential local algorithms —
  incl. BP/SP-guided decimation variants — fail for random NAE-k-SAT above
  ~2^{k-1} log^2 k / k); Coja-Oghlan 2011 (BP-guided decimation fails above
  ρ·2^k/k); Achlioptas–Coja-Oghlan 2008 (shattering at ~2^k ln k / k);
  **Bresler–Huang 2021**: low-degree polynomial algorithms fail to find
  solutions of random k-SAT above ~κ*·2^k log k / k; Gamarnik 2021 (PNAS survey).
  Myopic algorithms for 3-SAT fail above 3.26 (Achlioptas–Sorkin 2000).
* **Scope.** (a) Large-k asymptotics — at k = 3 no OGP barrier is proven below
  the threshold 4.267; empirically WalkSAT/SP work near 4.2–4.25. (b) Only the
  **search** problem on **random** instances; says nothing about refutation or
  worst case. (c) Only "stable" algorithms (Lipschitz in the input). Local
  search gives **no UNSAT certificates at all**, so it can never be the full
  algorithm.

### 1.12 Statistical-query and low-degree polynomial algorithms (refutation/
### distinguishing side)
* **Theorems.** Feldman–Perkins–Vempala 2015/2018: SQ algorithms need
  ~n^{r/2} samples to solve planted k-SAT (r = distribution complexity);
  low-degree likelihood-ratio framework (Hopkins–Steurer 2017; Hopkins 2018;
  Kunisky–Wein–Bandeira 2019): degree-D polynomials cannot distinguish random
  3-SAT from quiet-planted 3-SAT for m ≪ n^{3/2}. KMOW 2017 (SOS) as above.
* **Conjecture (not theorem).** Feige's hypothesis 2002: no poly-time
  algorithm refutes (one-sided: never says UNSAT on SAT; says UNSAT on most
  random formulas) random 3-CNF with m = Δn, Δ a large constant.
* **Scope.** Algorithms whose decision is a low-degree polynomial / SQ query of
  the input. Gaussian elimination is the standard counterexample (it beats
  low-degree on random XOR-SAT with planted solution / noiseless parity), so
  these are **heuristic** barriers, but they cover every known random-CSP
  refuter (spectral, SOS, Kikuchi).

### 1.13 Relativizing (and algebrizing) algorithms
* **Theorems.** Baker–Gill–Solovay 1975: oracle B with P^B != NP^B (and A with
  P^A = NP^A). Black-box search needs 2^n queries (Bennett–Bernstein–Brassard–
  Vazirani 1997 give 2^{n/2} even quantumly). Aaronson–Wigderson 2008:
  algebrizing techniques cannot resolve P vs NP in either direction.
* **Meaning for algorithms.** An algorithm (with correctness proof) that only
  accesses the formula through an "evaluate this assignment" oracle, or whose
  proof would go through relative to every oracle, is excluded. The algorithm
  must exploit the **explicit clause structure** (white-box), and its
  correctness proof must be non-relativizing (e.g. use arithmetization-like or
  structural properties of the formula, not just verification).

### 1.14 Levin's universal search (1973) — an explicit algorithm, not a barrier
* **Statement.** Dovetail all programs p_1, p_2, ... (p_i gets 2^{-i} of the
  time); whenever one outputs an assignment, verify it. For satisfiable
  formulas this halts within O(2^{|p|} · T_p(n) + poly) where p is the fastest
  search algorithm. Hence **if P = NP, this explicit algorithm finds
  satisfying assignments in polynomial time** (via search-to-decision self-
  reduction). Adding a brute-force fallback makes it always correct. For
  UNSAT, it needs in addition a proof system with short refutations (exists iff
  NP = coNP), e.g. dovetail over (program, ZFC-proof that the program is a
  correct decider) — polynomial iff P = NP is *provable* in that theory.
* **What it means for this project.**
  1. "Find an algorithm that is polynomial if P = NP" is **already solved**
     (Levin). The open content is entirely a **proof** that some explicit
     algorithm is polynomial — i.e. a mathematical theorem, not code.
  2. Experiments can only give evidence about *growth rates on finite n*; they
     cannot establish polynomiality. A candidate must come with a reason
     (invariant/potential/certificate argument) that could become a proof.
  3. Levin's constants (2^{|p|}) make it useless in practice; P = NP with a
     galactic polynomial would also be useless in practice.

### 1.15 Other relevant facts
* **Cryptography.** P = NP ⇒ no one-way functions (folklore; Impagliazzo 1995
  "five worlds"). The algorithm must invert e.g. Goldreich's local PRG/OWF
  (Goldreich 2000) instances encoded as 3-SAT, SAT-encoded factoring, SHA/AES
  preimage instances — none of which have any known structure exploitable by
  the families above.
* **Short certificates exist below the spectral threshold.** Feige–Kim–Ofek
  2006: random 3-CNF with m ≈ n^{1.4} have **polynomial-size** (NP)
  refutations whp — below n^{1.5}, where no poly-time refuter is known. So
  "short certificates exist but are hard to find" is a real phenomenon.
* **Frege / Extended Frege / IPS.** No superpolynomial lower bounds known.
  IPS lower bounds would imply VP != VNP (Grochow–Pitassi 2014).

---

## 2. Checklist — a polynomial 3-SAT decider proving P = NP must:

1. **Produce poly-size UNSAT certificates in some p-bounded proof system**
   (Cook–Reckhow 1979) — i.e. prove NP = coNP as a side effect.
2. **Not be a DPLL / lookahead / CDCL procedure staying inside resolution**:
   random UNSAT 3-CNF at density 4.26, PHP and Tseitin need 2^{Ω(n)}
   resolution (Chvátal–Szemerédi 1988; Haken 1985; Urquhart 1987;
   Pipatsrisawat–Darwiche 2011 for CDCL ⊆ resolution).
3. **Not be width/consistency-bounded**: needs width Ω(n)
   (Ben-Sasson–Wigderson 2001; Atserias–Dalmau 2008).
4. **Not be low-degree algebra over a field** (incl. Gaussian elimination,
   Gröbner bases with degree cap): PC degree Ω(n) on random 3-CNF over every
   field (Alekhnovich–Razborov 2003; Impagliazzo–Pudlák–Sgall 1999).
5. **Not be a constant/low-degree LP/SDP hierarchy or spectral method**:
   SOS degree Ω(n) at m = O(n) (Schoenebeck 2008; Grigoriev 2001), needs
   m >~ n^{1.5}/√d (KMOW 2017); poly-size LP/SDP for CSPs ≈ SA/SOS
   (Chan–Lee–Raghavendra–Steurer 2016; Lee–Raghavendra–Steurer 2015).
6. **Not be one fixed polynomial-size LP** of a natural polytope (Fiorini et
   al. 2012; Rothvoss 2014) — but note this doesn't exclude adaptive algorithms
   (matching is in P).
7. **Use negation essentially (non-monotone)** and **be non-constant-depth**
   (Razborov 1985; Alon–Boppana 1987; Håstad 1986; Razborov–Smolensky 1987).
8. **Beat cutting planes on clique–coloring** (Pudlák 1997) — i.e. its
   certificates cannot have feasible monotone interpolation; any system
   with feasible interpolation is excluded if one-way functions/clique-hardness
   style separations hold (Krajíček–Pudlák 1998 [unverified]).
9. **Not be an OBDD/d-DNNF compiler** (Bryant 1991; Friedman–Xu 2013;
   Bova et al. 2016).
10. **Refute random 3-CNF with m = Θ(n) clauses**, i.e. refute Feige's
    hypothesis (2002), and therefore beat every SQ/low-degree/SOS method below
    m = n^{1.5} (Feldman–Perkins–Vempala 2018; Hopkins 2018; KMOW 2017).
11. **Find solutions of random k-SAT beyond the OGP/shattering regime** for
    large k (Gamarnik–Sudan 2017; Bresler–Huang 2021) — so it must be
    **unstable** (not Lipschitz/low-degree in the input) for large k.
12. **Be white-box and non-relativizing** (Baker–Gill–Solovay 1975;
    Aaronson–Wigderson 2008): it must use clause structure, not only an
    assignment-evaluation oracle.
13. **Break all cryptography encodable in 3-SAT** (Impagliazzo 1995;
    Goldreich 2000).
14. **Come with a proof**, not just benchmarks — Levin 1973 already gives an
    explicit algorithm that is polynomial iff P = NP. If the soundness proof is
    formalizable in PV, EF becomes p-bounded (Cook 1975), so an EF lower bound
    would kill it.

---

## 3. Least-excluded directions (ranked by plausibility of a simple next experiment)

| Rank | Direction | Why not excluded | Simple experiment |
|---|---|---|---|
| 1 | **Extension variables (ER/EF-style), DRAT/PR clause addition** | No superpolynomial lower bounds for ER/EF/PR-with-new-vars; known to give poly proofs of PHP, Tseitin, clique–coloring-type families where resolution fails (Cook 1976 for PHP in EF; Heule–Kiesl–Biere PR/SDCL short PHP proofs). | On random UNSAT 3-CNF at 4.26 (n = 30…120, 50 instances each), compare CDCL proof size vs CDCL after extension-introducing preprocessing (BVA/SBVA, or definitions y = l1∧l2 / y = x⊕z chosen greedily), and vs SDCL/PR if available. Fit log(size) vs n; look for any change in the exponential base. |
| 2 | **Res(⊕) / branching on linear forms + Gaussian elimination** | Lower bounds known only for tree-like and regular Res(⊕) (Itsykson–Sokolov 2014; Efremenko–Garlík–Itsykson and others 2023–24 [unverified]); general Res(⊕) on random 3-CNF open. | DPLL that branches on parities x_i⊕x_j⊕… (chosen by a heuristic, e.g. sparse vectors from the clause–variable matrix's near-kernel), keeps a GF(2) system, prunes when a clause is falsified under the affine subspace. Measure tree size vs n against plain DPLL. |
| 3 | **Cutting planes / stabbing planes on random 3-CNF** | Size lower bound open for k = 3 (known only for k = Θ(log n)); only *rank* Ω(n) is known. | Run a CP/PB solver (RoundingSat-style: division + saturation) on random UNSAT 3-CNF, measure proof size growth; separately, small exhaustive search for low-size CP refutations at n ≤ 20 to compare with resolution minima. |
| 4 | **Algebra over non-field rings (Z_4, Z_6, Z) / non-Boolean encodings** | Alekhnovich–Razborov is over fields; ring cases only partially covered [unverified]. | Compute minimum Nullstellensatz/PC degree over Z_4 and Z_6 (Hermite/Smith normal form of the degree-d Macaulay matrix) for random 3-CNF n = 10…20; compare with GF(2), GF(3). Any degree drop vs fields would be notable. |
| 5 | **Non-proof-system / global certificates (Feige–Kim–Ofek style)** | FKO 2006 shows poly-size refutations exist at m ≈ n^{1.4}, below the SOS/spectral n^{1.5} barrier; how to *find* them (or extend to m = O(n)) is open. | Hard to test simply: implement FKO's witness (spectral + combinatorial parity-of-cycles certificate) at m = n^{1.4} for n ≤ 200 and check which parts are found by search vs guessed; measure how far below n^{1.4} any similar certificate survives. |

Expectation management: (1)–(3) are the standard "frontier" proof systems;
experts expect exponential growth there too for random 3-CNF, but it is not
proven. Positive evidence at n ≤ 150 would not be a proof; negative evidence
would be expected.
