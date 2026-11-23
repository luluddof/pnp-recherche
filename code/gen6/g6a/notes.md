# G6a — "Can measurements absorb the cut?" (Limits paper, §9)

## 0. The problem, quoted

Definition 5.1 (Limits and Constants of OpenAI's Claimed ETH Refutation). A scheme for a CNF F
(variables [n], occurrences O = {(e,v): v in I_e}, M = |O|) consists of

* (a) Deformation. For each occurrence o=(e,v) a vector P_o(T) in Q[T]^n of degree at most R with
  P_o(0)=e_v; H(T) = sum_x prod_e phi_e(<P_{e,v}(T),x>)_{v in e}. Own coefficient c_o(T)=(P_o(T))_v.
  deg H <= MR.
* (b) Recovery. A finite set X ⊂ Q\{0}, orders q_t >= 0, and a linear functional Λ on Hermite data
  (f^{(j)}(t))_{t in X, j<=q_t} with Λ(f)=f(0) for every f in Q[T] of degree at most MR.
* (c) Designs. For each t in X: a partition of [n] into groups S_{t,1},S_{t,2},...; a group α_t(e) for
  each clause; a subspace R_{t,h} ⊆ Q^{S_{t,h}} (the span of the measurement rows) for each group, such
  that for every clause e, every v in I_e and every group h ≠ α_t(e),
      **P_{e,v}(t)|_{S_{t,h}} ∈ R_{t,h}.**   (3)
  s_t = max_h |S_{t,h}|, ρ_t = Σ_h dim R_{t,h}, μ_t = #measured groups, w_t = tw(interaction graph).

Cost model (§5.2): κ(t) >= max(2^{s_t}, 2^{w_t}, C(n, min(q_t, n/2))) (DP-type), and additionally
κ(t) >= 2^{ρ_t} for the preprint's product-type evaluation.

§9: "From above, a design on an expander with √n groups, O(√n) interaction width, Taylor order O(√n)
and O(1) or polylogarithmically many rows per measured group would bring the scheme to 2^{Õ(√n)}.
From below, one must show that measurements cannot absorb the cut. [...] Turning this [2x2-minor
remark] into a bound for groups with several measurement rows is open."

So the requirement to refute is: for all t, s_t = Õ(√n), q_t = Õ(√n), ρ_t <= (#groups)·polylog = Õ(√n).

Useful exact facts (all elementary):
* (R1) Exactness of Λ up to degree MR forces Σ_t (q_t+1) >= MR+1 (else Π_t (T-t)^{q_t+1} has degree
  <= MR, zero Hermite data, and value ≠ 0 at 0).
* (R2) Given the design, the occurrences are independent: each P_o must be a degree-R curve with
  P_o(0)=e_v and P_o(t) ∈ W_{α_t(e),t} := Q^{S_{t,α_t(e)}} ⊕ ⊕_{h≠α_t(e)} R_{t,h} for all t.
  Cutting (c_o(t)=0) is not a separate choice; it is just one way of satisfying (3).

## 1. Results

Notation for a point t: a_h = |S_{t,h}|. For a "witness map" ω choosing one occurrence ω(v) of each
variable v, put I^ω_{t,h} = {v : α_t(e(ω(v))) = h}, and the deficiency
    def^ω_t = Σ_h (a_h − |I^ω_{t,h}|)_+ .
A variable is **split at t** if its occurrences lie in clauses assigned (by α_t) to >= 2 groups;
σ'_t = number of split variables, d_v = number of occurrences of v.
sep(F) = min |Y| such that every component of Pr(F) − Y has <= 2n/3 vertices.

### Lemma 1 (witness determinant). For every ω, Σ_{t∈X} (def^ω_t − ρ_t)_+ <= nR.
Proof. Let Φ(T) be the n×n matrix with row v equal to P_{ω(v)}(T). Φ(0)=I, entries have degree <= R,
so det Φ is a nonzero polynomial of degree <= nR. Fix t and a group h. For a row v with
α_t(e(ω(v))) ≠ h, (3) gives Φ_v(t)|_{S_h} ∈ R_{t,h}; the other rows are exactly I^ω_{t,h}. Hence the
column block Φ(t)[:,S_h] has rank <= min(a_h, dim R_{t,h} + |I^ω_{t,h}|), and
    nullity Φ(t) >= Σ_h (a_h − dim R_{t,h} − |I^ω_{t,h}|)_+ >= def^ω_t − ρ_t.
For a square polynomial matrix, ord_{T=t} det Φ >= nullity Φ(t) (change basis so the first k columns
span ker Φ(t); those columns of Φ(T)U are divisible by T−t). Summing over the distinct t ∈ X:
Σ_t nullity Φ(t) <= deg det Φ <= nR. ∎

### Lemma 2 (expected deficiency). For ω uniform (independently per variable, uniform over its
occurrences), and s = s_t:
    E def^ω_t = ½ Σ_h E|X_h − a_h| >= max( ½ Σ_h |μ_h − a_h| ,  c0/(2√(2s)) · Σ_{v split} 2(d_v−1)/d_v² ),
with X_h = |I^ω_{t,h}|, μ_h = E X_h, c0 = 0.216.
Proof. (y)_+ = (|y|+y)/2 and Σ_h a_h = Σ_h X_h = n give the identity. X_h is a sum of independent
Bernoulli(π_{v,h}), π_{v,h} = share of v's occurrences assigned to h; V_h = Var X_h <= μ_h.
Anticoncentration: for any real a, E|X−a| >= ½E|X−X'| (X' iid copy). Z = X−X' = Σ ξ_v, ξ_v symmetric
in {−1,0,1}. (i) Given the set of nonzero ξ_v (size k >= 1), P(Z=0) <= ½, so E|Z| >= P(Z≠0) >=
½(1−e^{−2V}); (ii) E Z² = 2V, E Z⁴ <= 2V + 12V², Hölder gives E|Z| >= (2V)^{3/2}/(E Z⁴)^{1/2} >= 0.756√V
for V >= 1. Hence E|X−a| >= 0.216·min(V, √V) (numerically checked: anticonc_check.py, worst ratio 0.74).
If μ_h <= 2s then V_h <= 2s and min(V_h,√V_h) >= V_h/√(2s); if μ_h > 2s then E|X_h−a_h| >= μ_h − s >=
μ_h/2 >= V_h/2. So E|X_h − a_h| >= c0 V_h/√(2s) always, and E|X_h−a_h| >= |μ_h − a_h| (Jensen).
Finally Σ_h V_h = Σ_v (1 − Σ_h π_{v,h}²) and for a split v with d_v occurrences 1−Σπ² >= 2(d_v−1)/d_v². ∎

### Lemma 3 (split variables on separators). For every t:
    either σ'_t >= sep(F), or ½ Σ_h |μ_h − a_h| > (2n/3 − s_t)/2.
Proof. Non-split variables with all clauses in group h form J_h; no clause meets J_h and J_{h'} (h≠h'),
so every component of Pr(F) minus the split set lies in one J_h. If all |J_h| <= 2n/3 the split set
is a 2/3-separator. Otherwise |J_{h0}| > 2n/3, so μ_{h0} >= |J_{h0}| > 2n/3 while a_{h0} <= s_t. ∎

### Theorem B (measurements cannot absorb more than a √s factor of the cut). For every scheme
satisfying Definition 5.1 (any clause-dependent deformation, any partitions/assignments, any
measurement subspaces), on a formula in which every variable occurs and occurs at most b times,
there is a point t ∈ X with
    ρ_t + q_t + 1  >  min( c_b · sep(F) / √s_t ,  (2n/3 − s_t)/2 ),   c_b = c0·(b−1)/(b²√2)  (c_3 ≈ 0.034).
Proof. Average Lemma 1 over ω and use Jensen ((·)_+ convex): Σ_t (E def_t − ρ_t)_+ <= nR <= MR <
Σ_t (q_t+1) by (R1). So some t has E def_t < ρ_t + q_t + 1. Lower-bound E def_t by Lemmas 2–3
(2(d−1)/d² >= 2(b−1)/b² for 2 <= d <= b). ∎
(Unbounded degrees: replace c_b·sep by (c0/(2√2))·Σ_{v split} 2(d_v−1)/d_v²; for random 3-CNF of
density Δ, Markov on degrees still gives Ω_Δ(n) once sep(F) = Ω_Δ(n), which is Lemma 8.2's argument.)

### Theorem A (variable-uniform deformations: no absorption at all). If P_{e,v} = P_v for all e ∋ v
(a global linear change of variables x ↦ Φ(T)x), then there is t with
    ρ_t + q_t + 1 > min( sep(F), 2n/3 − s_t ).
Proof. Same determinant with Φ_v = P_v; now the free rows of block h are J_h (variables all of whose
clauses are in h), so nullity Φ(t) >= Σ_h (a_h − |J_h|)_+ − ρ_t >= n − Σ_h min(a_h,|J_h|) − ρ_t, and
Lemma 3's case split gives n − Σ_h min(a_h,|J_h|) >= min(sep, 2n/3 − s_t). ∎

### Corollaries.
1. **The §9 target design does not exist** on any bounded-occurrence family with sep(F) >= cn (e.g.
   random 3-regular variable–clause incidence; random 3-CNF via Lemma 8.2): with s_t = Õ(√n) it needs
   ρ_t + q_t >= c'·n^{3/4}/polylog at some point, while the target has ρ_t + q_t = Õ(√n).
   Equivalently: with g_t groups, average rows per group must be >= c√(s_t) unless q_t >= c√(n g_t).
2. **Product-type cost >= 2^{Ω(n^{2/3})}** on such formulas (vs 2^{Ω(√n)} in Thm 8.3): at the point of
   Theorem B, κ(t) >= max(2^{s_t}, 2^{ρ_t}, 2^{min(q_t,n/2)}) and max(s, n/√s) >= n^{2/3}.
3. Variable-uniform schemes cost 2^{Ω(n)} (product-type) on expanders; since the preprint's
   construction reaches q, ρ, s = o(n), its deformation must be genuinely clause-dependent.

## 2. What is NOT proved / open
* DP-type model (charges 2^{w_t}, not 2^{ρ_t}): Theorem B forces many rows, but rows spread over
  many groups need not enter one DP bag; the DP-type gap 2^{Ω(√n)} .. 2^{O(ρ(n))} stays open.
* Tightness of the √s loss: the witness map is fixed across points, and a designer who varies α_t
  makes a fixed ω look random — the √s comes only from fluctuations. Whether some design achieves
  ρ ≈ n/√s (product-type cost 2^{Θ̃(n^{2/3})}) is open. Theorem A shows such a design must be
  clause-dependent.
* Theorems use only (a),(b),(c) and the degree bound; they do not cover the preprint's parity form or
  tensor contraction (§5.1 of Limits), which violate (b) or the degree bound.

## 3. Computation (exact, Fractions; files in this directory)
* scheme.py / designs.py / run_checks.py. Random 3-regular formulas (n=24, 30), random partitions per
  point, majority α_t, four valid design families (variable-uniform cut+measure, clause-dependent
  cut+measure, "global-sum" with one all-ones row per group, random dense deformation with minimal
  R_{t,h}). Condition (3) verified exactly; checked for 4 sampled ω per design: nullity Φ^ω(t) >=
  def^ω_t − ρ_t at every t and Σ_t nullity <= nR (equality nR is attained by cut designs, so Lemma 1
  is tight); exact E def_t >= Lemma 2–3 bound. Output: theorem_checks.txt.
* End-to-end (n=12, 3 groups of 4, R=1, |X|=6, q=6): H(T) computed exactly (mod 2^61−1), H(0) = #F =
  977, Hermite recovery from the data at X returns 977, and H(t) computed by the product-type design
  evaluation (enumerating reading tuples, counting inside groups) equals brute-force H(t) exactly in Q
  at every t. Reading-tuple counts reported (64–3072). Output: e2e_checks.txt.
* The naive designs are far from Theorem B's bound (ρ_t ≈ 0.5–1.0·n at these sizes) — consistent with
  the theorem but not a tightness test.
