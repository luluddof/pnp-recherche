# G8a — Limits §9, track 1: (a) regularity in the refined DP model, (b) tightness of n^{2/3}

Notation as in resultats/note_limits_s9.md (G6a) and code/gen7/g7b/notes.md (G7b): point t, groups
S_h (a_h), rows r_h, ρ_t = Σ r_h, witness map ω, X^ω_h = |I^ω_{t,h}|, Y_h = X_h − a_h (Σ_h Y_h = 0),
def^ω_t = Σ_h (−Y_h)_+ = ½Σ_h|Y_h|, interaction graph I_t (clique({α(e)} ∪ {h : some P_{e,v}(t)|_{S_h} ≠ 0})
for every clause e), K := max_t (ρ_t + q_t + 1), N^ω_t = Σ_h (a_h − r_h − X_h)_+ (G7b).
Checks: `checks.py` → `checks.txt` (imports G6a's scheme.py/designs.py read-only).

**Summary.** Neither (a) nor (b) is closed. Proved here:
* (b) Theorem B2 (adversarial witness): K > min(n/4, c_b sep(F)/√|X|). So a scheme reaching n^{2/3}
  needs |X| = Ω(n^{2/3}) points. Proposition S (barrier): Fact R plus Lemma 1 for *every* witness map
  (strong form, at all points at once) can hold with max(s_t, ρ_t, q_t) = O(√(b log n)·n^{2/3}). So
  the √s loss is built into the witness-determinant method, and n^{2/3} is the best that method can give.
  Building an actual scheme that is tight: not achieved.
* (a) Lemma H (boundary rank bound, new and unconditional). Corollary H1: the two sides of every
  bag-separation of the torso are mean-balanced up to W_t + 2q_t + 2. This rules out the *mean-surplus*
  form of the "clause-hub" escape. Proposition P (per-point barrier): a hub–satellite configuration
  satisfies every per-point inequality known so far (G7b's C, D, E, F, Thm 7.1, and Lemma H) with
  s, W = O(√n), q = O(1), E N^ω ≈ 0. Its naive vector realisation violates the *actual*-nullity form
  of Theorem D (a flow/double-counting effect). The next lemma needed is a Rado/flow-type rank bound.
  The unconditional n^{1/2+δ} remains OPEN.

---------------------------------------------------------------------------------------------------
## Part (b)

### Theorem B2 (adversarial witness)
Let F have all d_v ≤ b and M ≥ n. Take any scheme (Def. 5.1) with s_t ≤ n/6 for all t. Then
    K > min( n/4 ,  c_b · sep(F) / √|X| ),     c_b = 0.189/√(2·C(b,2))   (c_3 = 0.0772).

*Proof.* Lemma 1 for a fixed ω gives Σ_t(def^ω_t − ρ_t) ≤ Σ_t(def^ω_t − ρ_t)_+ ≤ nR ≤ MR < Σ_t(q_t+1).
So for every ω we have (★) Σ_t def^ω_t < Σ_t(ρ_t+q_t+1) ≤ |X|·K. Split X = X1 ⊔ X2 by Lemma 3's two cases.
* t ∈ X2 (some |J_{t,h0}| > 2n/3): for *every* ω, X^ω_{h0} ≥ |J_{h0}|, so def^ω_t = Σ_h (Y_h)_+ ≥
  |J_{h0}| − a_{h0} > 2n/3 − s_t ≥ n/2.
* t ∈ X1: Σ_t is a 2/3-separator, so σ'_t ≥ sep(F). Fix signs ε_{t,h} ∈ {±1}. Then
  ½Σ_h|Y_{t,h}| ≥ ½Σ_h ε_{t,h} Y_{t,h}, and
  Σ_{t∈X1,h} ε_{t,h}Y_{t,h} = Σ_v G_v(ω(v)) − C_ε. Here G_v(o) := Σ_{t∈X1} ε_{t,α_t(e_o)}, and C_ε does not
  depend on ω (it collects the a_h and the variables that never occur). Take ε uniform i.i.d. Then
  E C_ε = 0 and E G_v(o) = 0. For each v fix two occurrences o1, o2 that maximise
  N_v := |{t ∈ X1 : α_t(e_{o1}) ≠ α_t(e_{o2})}|. Then
  max_ω Σ_v G_v(ω(v)) ≥ Σ_v max(G_v(o1), G_v(o2)), and E max(G1,G2) = ½E|G1 − G2|.
  G1 − G2 = Σ_{t∈T_v}(ε_{t,h1(t)} − ε_{t,h2(t)}). The terms are independent across t and lie in
  {−2,0,2} with probabilities ¼,½,¼. So G1 − G2 =_d 2(B − B') with B, B' i.i.d. Bin(N_v,½).
  E|B − B'| ≥ 0.378√N_v: for N ≥ 4 this is G6a Lemma 2(ii) with V = N/4 ≥ 1; for N = 1,2,3 the exact
  values are 0.5, 0.75, 0.9375; check 1 confirms the ratio is ≥ 0.5 for N ≤ 1200.
  Since a maximum is at least an average over ε, there is an ω with
  Σ_{t∈X1} def^ω_t ≥ ½·Σ_v 0.378√N_v.
  A split variable has at most C(b,2) pairs of occurrences, so N_v ≥ N'_v/C(b,2), where
  N'_v = |{t ∈ X1 : v ∈ Σ_t}| ≤ |X1|. Then Σ_v √N'_v ≥ Σ_v N'_v/√|X1| = Σ_{t∈X1} σ'_t/√|X1|
  ≥ sep·√|X1|.
Using this ω in (★): |X|K > 0.189·sep·√|X1|/√C(b,2) + |X2|·n/2. If |X2| ≥ |X|/2, then K > n/4.
Otherwise |X1| ≥ |X|/2, and K > 0.189·sep/(√C(b,2)·√(2|X|)). ∎

*Consequences.* Combined with Theorem B:
K > max over the two bounds, i.e. K ≳ max(min(n/√s, n), sep/√|X|). An n^{2/3}-tight scheme on a
bounded-degree expander therefore needs |X| ≥ (c_b sep/K)² = Ω(n^{2/3}), and its design must vary across
points: if a single design were repeated, adversarial ω gets def ≥ σ'/4 at every copy (random-sign
argument at a single t), forcing K = Ω(n). There is no matching upper bound on |X|: Fact R only bounds
|X| from below (|X| ≥ MR/K). This is why B2 does not improve the exponent unconditionally.
*Self-review.* The step from max over ω to the average over ε is valid (for each ε pick the best ω;
averaging over ε shows some ω achieves at least the expectation). Independence across t holds because
each term uses only ε_{t,·}, and h1(t) ≠ h2(t) holds on T_v. Variables that never occur have no choice
and enter C_ε, whose ε-mean is 0. The case-2 bound uses def = Σ_h(Y_h)_+, valid since ΣY = 0.

### Proposition S (the witness-determinant method cannot beat n^{2/3})
Let F be a 3-CNF with d_v ≤ b, M ≥ n, n ≥ n_0(b), and 15 ln n ≤ s ≤ n, with g := ⌈n/s⌉. Then there are
R, X, q, and for each t ∈ X a partition (groups of size ≤ 2s), an α_t and numbers r_{t,h} ≡ r :=
⌈9√(b s ln n)⌉ such that
  (i) Fact R: Σ_t(q+1) ≥ MR+1, with q ≤ 2b+2;
  (ii) for EVERY witness map ω: Σ_t N^ω_t = Σ_t Σ_h (a_h − r − X^ω_{t,h})_+ ≤ nR (strong Lemma 1, hence
       also the def − ρ form and Theorem D's inequality with Z_t = ∅);
  (iii) ρ_t = g·r ≤ 10√(b ln n)·n/√s.
At s = n^{2/3}: max(s_t, ρ_t, q_t) = O(√(b log n)·n^{2/3}). So no argument that uses only {Fact R, Lemma 1
(any ω, all points), Lemma 3} can exceed exponent 2/3. This is what Theorem B extracts, so the √s
loss is a property of the method. Whether a real scheme can do this is not addressed.

*Proof.* At each point, independently: give every variable an i.i.d. uniform label in [g] (the groups),
and set α_t(e) := label of a uniformly random variable u_e ∈ e. Fix ω. For variables that occur, set
k_u := #{v : u_{ω(v)} = u}. A variable that never occurs counts toward its own label (k_v = 1). Then
Σ_u k_u = n, 0 ≤ k_u ≤ 3b, and X_h − a_h = Σ_u 1[lab(u)=h](k_u − 1). Conditioned on the u_e, this is a sum
of independent terms with mean (1/g)Σ(k_u − 1) = 0, range ≤ 3b, and variance ≤ (1/g)Σ(k_u−1)² ≤
(1/g)·3b·2n ≤ 6bs. By Bernstein with r ≤ 6s (true since s ≥ 3b ln n),
P(a_h − X_h ≥ r) ≤ exp(−r²/(24bs)) ≤ n^{−3.37}. Hence E N^ω_t ≤ g·n·n^{−3.37} ≤ 1/n.
Across points, the N^ω_t ∈ [0,n] are independent. Hoeffding gives
P(Σ_t N^ω_t ≥ |X|(1+1/n)) ≤ exp(−2|X|/n²). There are ≤ b^n witness maps, so with |X| := ⌈n³ ln b⌉ + n²,
a union bound shows that with positive probability every ω has Σ_t N^ω_t ≤ |X|(1+1/n). Chernoff plus a
union bound over |X|·g groups gives all group sizes ≤ 2s (needs s ≥ 15 ln n). Set R := ⌈2|X|/n⌉, so
nR ≥ |X|(1+1/n), giving (ii). Set q+1 := ⌈(MR+1)/|X|⌉ ≤ 2M/n + 2 ≤ 2b+2, giving (i). ∎
*Check 3* (n = 3000, 8000; s = n^{2/3}; 200/60 random ω): the deviations of max_h(a_h − X_h) are
≤ 0.73·√(bs ln n), against 9 used in the proof. N^ω = 0 throughout. Mean def is ≈ 0.2 × the bound g·√(6bs)/2. The constant
9 is wasteful; at these n, r exceeds s, so the check only confirms the tail estimate.

### What a tight real scheme must look like (conclusions for (b))
Clause-dependent (Thm A); |X| = Ω(n^{2/3}) with designs that vary across points (B2); and at most points
r_h ≈ fluctuation of X_h for most groups (Thm B). The rows must absorb Θ(s) foreign restrictions per
group in Θ(√s) dimensions. At a single point this is possible (e.g. rank-1 all-ones restrictions:
2×2 minors with nonzero diagonal). The obstruction is polynomial: a dimension count of the constraints
P_o(t) ∈ L_{o,t} (codim ≈ n − s − ρ per point, ≥ 3nR/q points, nR unknowns) suggests Θ(n) *for generic
constraint subspaces*. This is a heuristic, not a proof. Structured families tried mentally (moment /
Vandermonde rows, all-ones rows, cut + measure) all reduce to cuts and give Θ(n). No candidate
n^{2/3} construction was found, and none was verified on small n.

---------------------------------------------------------------------------------------------------
## Part (a)

### Lemma H (boundary rank bound; unconditional)
For every scheme, point t, witness map ω and set A of groups, with ∂A := N_{I_t}(A) \ A:
    nullity Φ^ω(t) ≥ Σ_{h∈A}(X^ω_h − a_h) − Σ_{h∈∂A} r_h.
*Proof.* Let S = {v : α_t(ω(v)) ∈ A}; then |S| = Σ_A X_h (variables that never occur: home = own group).
Take v ∈ S with ω(v) = (e,v) and α(e) = h ∈ A. Its row P_{e,v}(t) restricted to a block h' ≠ h lies in R_{h'}
by (3). If it is nonzero, h' ∈ clique(e) ⊆ N_{I_t}[h] ⊆ A ∪ ∂A. So every row of S lies in
⊕_{h'∈A} Q^{S_{h'}} ⊕ ⊕_{h'∈∂A} R_{h'}, whose dimension is Σ_A a + r(∂A). Then
n − rank Φ ≥ |S| − rank(rows S). ∎
Check 2: exact ranks over Q, G6a's 4 design families, n = 24, 30, all subsets A, 3 ω per point; 5160
instances pass; the bound is attained in 56 of them (random-dense family).

Lemma H complements G6a Lemma 1 (deficit groups need their *own* rows). Surplus regions need rows on
their *boundary*. Averaging: E nullity ≥ Σ_A(μ_h − a_h) − r(∂A). ω-dependent choices of A are allowed.

### Corollary H1 (bag-separations are mean-balanced)
Let t be the point of Theorem D (|Z_t| + E_ω nullity Φ^ω(t) < 2(q_t+1); its proof bounds the actual
nullity). Take any tree decomposition of the torso G_M, any bag B, and any component C of G_M − B. Put
A := C ∪ {unmeasured h with N(h) ∩ C ≠ ∅}. Then ∂A ⊆ B, so
    Σ_{h∈A}(μ_h − a_h) < r(B) + 2q_t + 2 ≤ W_t + 2q_t + 2.
Applying the same bound to the union A'' of the remaining sides and using Σ_all(μ − a) = 0 gives
|Σ_A(μ − a)| ≤ W_t + 2q_t + 2 + |Σ_{h∈B}(μ_h − a_h)|.
*Proof.* Unmeasured groups are pairwise non-adjacent (G7b), and N(h) is a clique of G_M for unmeasured h.
So an unmeasured h meeting C has N(h) ⊆ C ∪ B. A measured neighbour of C outside C lies in B. ∎
Meaning: a clause-hub separator whose surplus is balanced by row-heavy deficit groups on one side is
impossible unless the bag itself carries that surplus in rows. This kills the mean form of the escape
route named in G7b §4. It does NOT handle fluctuations: per ω, Lemma H on a side gains only
√(variance of a whole side) ≈ √n, while the n^{2/3} mechanism needs Σ_h √V_h ≈ n/√s.

### Proposition P (per-point barrier for (a); the remaining obstruction)
On any 3-regular F (n large), take one point with g hubs H_j of size a ≈ 2√n and g pendant groups D_j of
size d = ⌈C√(s log n)⌉ (n = g(a+d)). Choose the D_j so that no clause meets two of them, and give each
D_j's clauses to H_j. Distribute the other clauses so that |C_{H_j}| = a_j + d, which is μ_{H_j} for
3-regular F. Hubs get r = 1; D_j are fully measured with r = d. Interaction: clique(e) = {α(e)} ∪ groups
of the variables of e. Then: s = O(√n), Z = ∅, P(X_{H_j} < a_j − 1) ≤ n^{−C²/4}, so E N^ω ≤ n^{−1} and q = 0
suffices in Theorem D. Lemma H holds in expectation for all A (Σ_A(μ − a) ≤ Σ_{H_j∈A, D_j∉A} d ≤ r(∂A)),
and per ω whenever every hub has ≥ C√(s log n) hub neighbours. W ≤ g + d + 1 (hubs in one bag, pendant
bags {H_j, D_j}). Theorems C and 7.1 hold automatically. So every per-point inequality now known allows
max(s, W, q) = O(√n). [Proof straightforward; Chernoff/Bernstein as in Prop. S; not machine-checked.]
**But** the obvious vector realisation fails. Take P_o(t) = 𝟙_{S_{G(v)}} on foreign occurrences and
e_v + generic on α(e) otherwise. The D_j rows are the only rows that can fill the D_j columns, and they
are also counted in X_{H_j}. Lemma 1 double-counts them. The true nullity is ≈ Σ_j E(a_j − (X_j − d) − 1)_+
≈ g√s = n^{3/4} > 2q + 2, which violates the actual-nullity form of Theorem D. Other realisations route
D-columns through rows of other hubs (generic D-components). These keep W small only if each clause
touches ≤ 1 D-block, and then full rank needs a Hall condition on a random flow. Whether that holds
is UNVERIFIED (generic rank = Rado formula min_S[n − |S| + dim Σ_{v∈S}U_v]).
**Next step suggested.** Replace N^ω by the Rado bound nullity ≥ max_S(|S| − dim Σ_{v∈S}U_{ω(v)}), where
U_o = Q^{S_{α(e)}} ⊕ ⊕_{h∈clique(e)\α(e)} R_h. It contains both Lemma 1 and Lemma H. Then show that on an
expander, the per-group fluctuations of every bag's I_t-neighbourhood force Σ_B r_h ≳ n/√s. This would
remove (γ,β)-regularity. Not proved.

---------------------------------------------------------------------------------------------------
## Status of claims
PROVED: Theorem B2; Proposition S; Lemma H; Corollary H1.
SKETCHED (routine, not checked): Proposition P (the shadow part).
OPEN: unconditional refined-DP n^{1/2+δ}; a real n^{2/3}-tight scheme, or any exponent > 2/3 in the
product-type model (Prop. S shows a new ingredient beyond witness determinants is required).
Model caveat: in the honest DP, a measured group's factor has scope N[h]. G7b's torso adds that clique
only for unmeasured h, so G7b's model may undercharge. Under the N[h]-clique model,
W ≥ max_h r(N[h]); this does not by itself kill Proposition P.
