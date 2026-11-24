# G14c — A sound DP cost model for clause-product schemes (the "interface-rank" model), and a lower bound in it

Follow-up of `resultats/article/main.tex` §6 ("the DP-type model is degenerate"), of G7b (refined model W_t)
and of G10a (W_t is unsound against row-level elimination). Notation of the paper: scheme (Def. 5.1 of
Dias–Lima = DL), point t, groups S_h (a_h, s_t = max a_h), alpha_t, R_{t,h}, rho_t, q_t, w_t = tw(I_t),
witness map omega, Phi^omega, sep(F), sepw(F), theta(d) = 2(d-1)/d^2, c0 = (1-e^{-2})/4 = 0.2162.
Files: `lib.py`, `checks.py` -> `checks.txt` (all assertions pass). Status tags: [PROVED], [SKETCH], [OPEN].

## 0. Summary

* Two models are defined. The **row model** (design-respecting; rows/functionals of single groups are the
  atomic DP variables; charge max(2^{s_t}, 2^{lw_t}, Taylor)) is pointwise >= DL's DP-type charge and is not
  degenerate, but it is sound only for evaluators whose DP variables are functionals of ONE group: an
  evaluator that carries partial sums of readings across groups beats it exponentially (§5, 2^m vs poly).
* The **interface-rank (IR) model** is design-free and sound for a much larger class (§2): charge
  max(2^{bw_t}, C(n, min(q_t, n/2))), where bw_t is the branch-width of the cut-rank function of the form
  matrix A_t = (P_o(t))_o on the leaf set (clauses + variables). It covers in-group DP (which DL exclude),
  row-level elimination and aggregated readings. Non-degenerate: the trivial scheme of the paper's §6 has
  bw_t >= (2/3)(tw(Pr F)+1) at every point, i.e. 2^{Omega(n)} on expanding formulas.
* **Lower bound obtained (unconditional, IR model):** on every formula with n >= 5 (3-CNF), some point has
  bw_t + q_t + 1 > (c0/6)·min(sepw(F), sqrt(sepw(F))). So on the paper's families the IR cost of every
  scheme is 2^{Omega(sqrt n)} — against ALL interface-table DPs, without charging 2^{s_t}.
* **The target n^{1/2+delta} is NOT reached.** What is proved beyond sqrt(n) is a multi-cut inequality
  (Theorem 3): the cut ranks of any y tree edges that fragment the formula into parts of <= s' variables SUM
  to more than min(c·sepw/sqrt(s'), 2n/3 - s') - q_t - 1. It gives n^{2/3} exactly for evaluators that pay
  the sum (product-type, "enumeration nodes"), and only sqrt(n) for true DPs, which pay the max. Prop. M
  shows the averaged single-cut inequality cannot certify more than O(sqrt(n log n)). [OPEN]

## 1. Definitions

Fix a scheme and t in X. **Form matrix** A_t: the M x n matrix over Q with row o = (e,v) equal to P_{e,v}(t).
For a set E of clauses and S of variables, A_t[E,S] is the submatrix (rows = occurrences of clauses in E).
H(t) = sum_{x in {0,1}^n} prod_e phi_e(A_t[e,·] x).

**Leaves** L = clauses ∪ variables. For a bipartition L = (E1 ∪ S1) ⊔ (E2 ∪ S2):
    kappa_t(E1 ∪ S1) := rank A_t[E1, S2] + rank A_t[E2, S1]          (cut rank; symmetric)
(what side-1 clauses read from side-2 variables, plus the converse).
**bw_t** := min over cubic trees T with leaf set L of max over edges f of kappa_t(side of f).
**IR charge:** kappa_IR(t) := max( 2^{bw_t}, C(n, min(q_t, floor(n/2))) );  cost = max_t kappa_IR(t).
(No 2^{s_t} term: groups, alpha, R_{t,h} do not enter; they are one evaluator strategy among others.)

**Row model** (design-respecting). Effective spaces: L_{g->h} := span{P_{e,v}(t)|_{S_h} : alpha(e)=g, v in e}
(h != g), R'_h := sum_{g != h} L_{g->h} ⊆ R_{t,h}, r'_h = dim R'_h, M' = {h : R'_h != 0}. A *row
decomposition* is: finite families 𝔉_h ⊆ R'_h \ {0} spanning R'_h; sets F_{g,h} ⊆ 𝔉_h with
span F_{g,h} ⊇ L_{g->h}; a tree decomposition of the graph on ⊔_h 𝔉_h in which, for every group g,
K_g := 𝔉_g ∪ ⋃_h F_{g,h} is a clique (so K_g lies in one bag, at a node n(g)). Width of a bag B:
sum_h rank(B ∩ 𝔉_h).  **lw_t** := min over row decompositions of the max bag width.
Row charge: max(2^{s_t}, 2^{lw_t}, Taylor). (With 𝔉_h a basis this is tw+1 of G10a's "row graph", minimised
over bases; families may be redundant, the definition is basis-free.) lw_t <= rho'_t <= rho_t (one bag).

## 2. Soundness (c): which evaluators are lower-bounded

**Value-count lemma.** [PROVED] If B is a rational matrix of rank r with s columns, |B·{0,1}^s| >= 2^r.
(Take r independent columns; on inputs supported there x -> Bx is injective.)

**Interface-table DP.** A rooted binary tree with leaf set L. For a node u with leaves (E_u, S_u) the table is
    Tab_u[a, b] = sum_{x in {0,1}^{S_u} : A[E_u^c, S_u]x = a}  prod_{e in E_u} phi_e( A[e,S_u]x + b_e ),
indexed by all achievable a (values exported to outside clauses) and b (values of A[E_u, S_u^c]y imported).
Merge: Tab_u[a,b] = sum over (a1,a2) with (a1+a2)|_{E_u^c} = a of Tab_{u1}[a1, b|_{E1}+a2|_{E1}] ·
Tab_{u2}[a2, b|_{E2}+a1|_{E2}]; H(t) = Tab_root. (Implemented and checked against brute force, Part 1.)

**Proposition 1 (soundness).** [PROVED]
(i) Tab_u has >= 2^{kappa_t(u)} entries (value-count lemma on both index sets; they are independent because
    S_u and S_u^c are disjoint). Hence every interface-table DP stores a table of >= 2^{bw_t} entries.
(ii) Triangle inequality: for a node p with children u, u2: kappa(u) <= kappa(u2) + kappa(p).
(iii) If an evaluator, at each node, either merges two materialised children or computes the node's table by
    brute-force enumeration of the subtree's variables (cost >= 2^{max(|S_p|, rank A[E_p,S_p^c])}), then
    bw_t <= 2·log2(cost). [every node v below a brute-forced p has kappa(v) <= |S_p| + rank A[E_p,S_p^c]:
    rank A[E_v,S_v^c] <= |S_p \ S_v| + rank A[E_p,S_p^c] and rank A[E_v^c,S_v] <= |S_v|]
(iv) bw_t <= lw_t + s_t <= rho_t + s_t. Hence the product-type evaluator, DL's DP-type evaluators
    (group enumeration + elimination of group readings), G7b's group-atomic and G10a's row-level
    eliminators all cost >= 2^{bw_t/2}.
(v) Occurrence-level DPs (clause e split into its k occurrences via phi_e = 1 - prod_v psi_{e,v}, one bit
    c_e carried across every cut that splits e; cut cost rank+rank+#split clauses) have width
    >= bw_t/(k-1).
Proof of (ii): S_u^c = S_{u2} ⊔ S_p^c and E_u^c = E_{u2} ⊔ E_p^c; split both ranks accordingly.
(iv): take an optimal row decomposition, make its tree cubic (copies of bags), hang at n(g) a caterpillar
on {e : alpha(e)=g} ∪ S_g. A tree edge of the row decomposition separates the groups into A | A^c and its
adhesion contains every F_{g,h} with g,h on different sides, so kappa <= sum_{h in A} dim L_{A^c->h} +
sum_{h notin A} dim L_{A->h} <= lw_t. The edge above the caterpillar of g: kappa = r'_g + sum_h dim L_{g->h}
<= width(K_g) <= lw_t. Inside the caterpillar: kappa <= a_g + sum_h dim L_{g->h} <= s_t + lw_t.
(v): keep one occurrence leaf per clause; each split clause moves <= k-1 rows across the cut.

**Not covered** (stated limits of (c)): evaluators that compress tables (state merging, i.e. the algebraic
rank of the 2^{S1} x 2^{S2} cut matrix, which can be smaller than 2^{kappa}; low-degree/moment
representations; Fourier/inclusion–exclusion), fast matrix multiplication at joins, non-tree algorithms, and
any use of Taylor data other than the model's binomial charge. For these the model is a proposal, exactly as
DL's charges are. 2^{rank} is a LOWER bound on table sizes; real tables can be much larger (a rank-1 generic
form takes 2^s values), so small bw_t does not mean fast evaluation.

## 3. Non-degeneracy (b) and comparison with DL (d)

**Proposition 2 (b).** [PROVED] For the trivial scheme of the paper's Prop. 6.1 (P_o ≡ e_v, any groups, any
measurement spaces) and every t:  bw_t >= (2/3)(tw(Pr F)+1). So its IR charge is >= 2^{(2/3)(tw+1)}
= 2^{Omega(n)} on formulas of linear primal treewidth (DL DP-type charge: 2^{ceil(sqrt n)}).
In the row model: max(s_t, lw_t) >= (tw(Pr F)+1)/3 by Prop. 1(iv).
Proof. Here kappa(f) = |C_f|, C_f = {u : f separates the leaf u from the leaf of some clause containing u}.
Bags: B_x = union of C_f over the edges at x (leaf e: I_e; leaf u: {u}). The nodes whose bag contains u are
the nodes of the union of the tree paths from leaf u to the clause leaves containing u: connected. Every
clause is covered by its leaf bag. At an internal node x every u in B_x lies on a path through x, hence in
at least two of the three C_f: |B_x| <= (3/2) bw. So tw(Pr F) <= (3/2) bw_t - 1. ∎

**Proposition 3 (d).** [PROVED]
(d1) Row model: lw_t >= tw(G_{M'}) + 1 >= w_t, for every scheme and point (G_{M'} = torso of I_t on the
     effectively measured groups). So the row charge is pointwise >= DL's DP-type charge, and <= the
     product-type charge.
     Proof. From a row decomposition take bags B'_x = {h : B_x ∩ 𝔉_h != ∅}; |B'_x| <= width(B_x) (nonzero
     functionals). {x : h in B'_x} is a union of subtrees all containing n(h): connected. For a clause e with
     a = alpha(e), the bag at n(a) meets 𝔉_h for every h in clique(e) ∩ M' (F_{a,h} != ∅ since L_{a->h} != 0,
     and 𝔉_a ⊆ K_a). An effectively unmeasured u is never read, so N(u) = {h : L_{u->h} != 0} ⊆ B'_{n(u)}.
     Every edge of I_t has an end in M' and w_t <= tw(G_{M'}) + 1 as in DL §5.2. ∎
(d2) IR model: NOT pointwise >= 2^{w_t}, and it should not be: w_t depends on the design, not on the
     function evaluated. Path formula, trivial scheme, scattered fully measured groups: w_t = 3 (K4), 4 (K5)
     while an explicit DP has width 2 (Part 6). On DL's barrier families (bounded occurrence with linear
     separators; random 3-CNF, Delta > 3) the IR model gives the same order 2^{Omega(sqrt n)} as DL's
     Theorem 8.1 (Theorem 2 below), for a larger evaluator class. DL's Cor. 7.3 (paths cost 2^{Omega(n^{1/3})})
     is not reproduced, correctly: it prices enumeration, and paths are linear-time countable.
     Caveat: DL's Theorem 8.1 needs no occurrence bound (it uses tw); Theorem 2 uses sepw.

## 4. Lower bounds

**Lemma 1 (cut inequality; basis-free form of the paper's Lemma 4.2 and of G8a's Lemma H).** [PROVED]
For every scheme, t, witness map omega and partition of L into parts (E_i ∪ S_i), i = 1..m, with
X_i := #{v : the clause of omega(v) is in E_i} (a never-occurring v counts for the part containing its leaf):
    nullity Phi^omega(t) >= sum_i ( |S_i| - X_i - rank A_t[E_i^c, S_i] )_+ .
For m = 2:  | X_1 - |S_1| | <= nullity Phi^omega(t) + kappa_t(E_1 ∪ S_1); more precisely the deficit
|S_1| - X_1 is paid by rank A[E_2,S_1] and the surplus X_1 - |S_1| by rank A[E_1,S_2].
Proof. Column block S_i of Phi: the rows with home outside E_i are rows of A[E_i^c, S_i] (or zero), so the
block has rank <= min(|S_i|, X_i + rank A[E_i^c,S_i]); sum over the blocks. ∎
It holds verbatim for occurrence-level cuts (O_1 | O_2 any bipartition of occurrences).

**Theorem 2 (sqrt-barrier in the IR model).** [PROVED] Let F be a k-CNF on n >= 3k/2 variables. Every scheme
has a point t with
    bw_t + q_t + 1  >  (c0/6) · g(sepw(F)),     g(y) = min(y, sqrt y),   c0/6 = 0.036.
If d_v <= b: sepw >= theta(b) sep(F). So on the families of the paper's Cor. 5.8 (sep >= cn, d_v <= b; or
random 3-CNF, Delta > 3, w.h.p.) every scheme has IR cost >= 2^{((c0/6) sqrt(eta' n) - 1)/2} = 2^{Omega(sqrt n)}.
Proof. Sum_t nullity Phi^omega(t) <= n_occ R <= MR < sum_t (q_t+1) (paper, Lemmas 4.1, 4.2), so averaging
over uniform omega there is t with E null_t < q_t + 1. Fix an optimal tree at this t. For an edge f with
sides 1,2 let J_i(f) = variables all of whose occurrences are on side i (never-occurring: side of the leaf).
Case 1: some edge has |J_1|, |J_2| <= 2n/3. The split variables Sigma_f (occurrences on both sides) form a
2/3-separator with d_v >= 2, so V := Var X_1 = sum_{Sigma_f} pi(1-pi) >= sum (d_v-1)/d_v^2 >= sepw/2, and by
Lemma 1 and the paper's anticoncentration lemma  bw_t >= kappa(f) >= c0 g(V) - E null >= (c0/2) g(sepw) - E null.
Case 2: otherwise orient every edge towards its side with |J| > 2n/3; a sink exists; it is not a leaf (a leaf
has |J| <= k <= 2n/3); so it is an internal node whose three branches T_1,T_2,T_3 have |J(T_i)| < n/3. The
variables with occurrences in two branches form a 2/3-separator Sigma_3. Lemma 1 with the 3-partition and
sum_i rank A[E_i^c,S_i] <= sum_i kappa(f_i) <= 3 bw_t:  3 bw_t + null >= def = (1/2) sum_i |X_i - |S_i||, and
E def >= (c0/2) sum_i g(V_i) >= (c0/2) g(sum V_i) >= (c0/2) g(sepw)  (g concave, g(0)=0; sum_i V_i >=
sum_{Sigma_3} theta(d_v)). So bw_t >= (c0/6) g(sepw) - E null/3. ∎
With Prop. 1(v) the same bound holds for occurrence-level DPs with constant c0/(6(k-1)).

**Theorem 3 (multi-cut).** [PROVED] For every scheme, t, omega, tree T (any degrees) with leaf set L and set Y of edges, with
parts = leaf sets of the components of T - Y:
    sum_{f in Y} kappa_t(f) + nullity Phi^omega(t) >= def_Y^omega := sum_i (|S_i| - X_i)_+ ,
and if every part has at most s' >= 1 variables,
    E_omega def_Y >= min( c0/(2 sqrt(2 s')) · sepw(F),  2n/3 - s' ).
Hence at the point t of Theorem 2: sum_{f in Y} kappa_t(f) + q_t + 1 > min(c0 sepw/(2 sqrt(2s')), 2n/3 - s').
Proof. Lemma 1; each edge of Y borders exactly one part on each side, and rank A[E_i^c,S_i] is at most the
sum over the edges f bordering part i of rank A[far clauses of f, S_i]. The expectation bound is the paper's
Lemma 4.4 + Lemma 4.5 with parts in place of groups (only |S_i| <= s' and sum_i X_i = n are used). ∎
Consequences.
(a) Star tree whose branches are the groups (the product-type evaluator): sum_f kappa(f) = sum_h (r'_h +
    f_h): this is the paper's Theorem 4.7 again, hence n^{2/3} when the evaluator pays 2^{sum} and 2^{s'}.
(b) "k enumeration nodes": if the evaluator enumerates all interface tuples at k nodes whose removal leaves
    parts of <= s' variables, and enumerates inside parts, its exponent is >= max(s', (c n/sqrt(s') - q)/k):
    >= c'(n/k)^{2/3}. This beats sqrt n iff k <= n^{1/4 - delta}. [PROVED from Thm 3; conditional on the class]
(c) A true DP pays max_f kappa(f) >= (1/|Y|) sum; fragmenting into parts of size s' needs |Y| >= n/s' - 1
    edges, so (c) only gives bw_t >~ sqrt(s'): no gain over Theorem 2.

**Proposition M (method barrier).** [PROVED] For every formula (m clauses) and every cubic tree on the clause
leaves, the variable leaves can be attached so that for every edge f of the resulting tree
E_omega |X_1(f) - |S_1(f)|| <= sqrt(n)/2 + sqrt(n ln(4m)/2). So Lemma 1 with uniform omega, one cut at a time,
cannot give more than O(sqrt(n log n)), whatever the scheme.
Proof. Attach v next to the clause of omega_0(v), omega_0 uniform. Then |S_1(f)| = X_1^{omega_0}(f) on old
edges; Hoeffding + union bound over < 2m edges; E|X - mu| <= sqrt(V) <= sqrt(n)/2; new edges have O(k). ∎

**Why n^{1/2+delta} is open.** A deficit/surplus of a cut is a sum of n independent choices: it fluctuates
by sqrt(n) only; adversarial omega gives |Sigma_f|/2 = Omega(n) for one cut, but the good point t depends on
omega (G8a's Theorem B2 trades this for sep/sqrt|X|). Fluctuations add over a fragmentation (n/sqrt(s')),
but so do the cut ranks; a DP pays their maximum. Any improvement needs either a non-product distribution on
witness maps that is anticoncentrated at scale n^{1/2+delta} on every low-rank cut, or an argument that the
SAME rows cannot serve many nested cuts (a rank inequality across the tree, not cut by cut). Neither is
available. No candidate scheme with IR exponent o(n) on expanders is known either. [OPEN]

## 5. Self-review (adversarial)

1. **Breaking the row model (as G10a did for W_t).** Hub group g = {u}, satellites h_1..h_m of s variables;
   the hub clause reads the single form sum_i 1_{S_{h_i}}. Every row decomposition has the m functionals
   1_{S_{h_i}} in K_g: lw_t >= m, charge 2^m. An evaluator carrying the partial sums z_i = z_{i-1} +
   1_{S_{h_i}}·x (at most ms+1 values) runs in poly(ms). So the row model is sound ONLY for evaluators whose
   DP variables are functionals of single groups. In the IR model the same instance has an explicit tree of
   width <= 3 (Part 4), i.e. the IR model already accounts for this evaluator. This is why the IR model, not
   the row model, is the proposed answer; the row model is kept for (d1).
2. **Breaking the IR model.** (a) Splitting clauses into occurrences: loses at most k-1 (Prop. 1(v)).
   (b) Lazy/symbolic tables (pure-clause subtrees are never tabulated by a sensible evaluator): loses at
   most a factor 2 in the exponent (Prop. 1(ii),(iii)); without this remark (i) alone would overcharge, e.g.
   the cut (all clauses | all variables) has kappa = rank A_t ~ n. (c) State merging: a cut whose side 1 is
   one variable plus many clauses has algebraic rank 2 but kappa large — this is laziness again and is
   absorbed by the minimum over trees up to the factor 2; I have NO proof that the algebraic rank of every
   balanced cut is >= 2^{c·kappa}, so rank-based ("linear-algebraic") DPs are not covered. [OPEN]
   (d) High-degree joins by fast matrix multiplication: not covered (constant in the exponent at most
   omega/3 for table-size bounds of triangles; not analysed).
3. **Theorem 2 checks.** n_occ R <= MR needs no hypothesis; never-occurring variables are isolated in Pr(F)
   and have a deterministic home; Case 2 needs k <= 2n/3. Case 2 never occurred in 800 random trees (Part 5)
   — it is exercised only by the proof, not by the numerics. At n <= 7 the bound is vacuous ((c0/6)sqrt(sepw)
   < 1): the computations test Lemma 1 / Theorem 3 pointwise (exhaustively over all 4096 cuts for n = 6),
   not the asymptotics.
4. **(d) is two statements**, not one: pointwise domination holds only for the row model; for the IR model
   only "same order on the barrier families", and bounded occurrence is needed (sepw), unlike DL's Thm 8.1.
5. The IR model does not charge 2^{s_t}; on formulas of small treewidth it charges little (as it should).
   The paper's n^{2/3} (product-type) is not contradicted: it prices enumeration of reading tuples.
6. Not re-derived here: the paper's Lemmas 4.3 (anticoncentration, Lean-checked) and the expected-deficiency
   computation, which I re-proved for parts only on paper (three lines: if mu_i >= 2a_i use |mu_i - a_i| >=
   V_i/2, else V_i < 2s'; sum_i V_i >= sum_Sigma theta; the 2n/3 - s' branch as in Lemma 4.5).

## 6. Computations (`checks.txt`)
Part 1: interface-table DP implemented; equals brute-force H(t) on G6a's four design families (n = 6, 7),
  1120 node tables, size >= 2^kappa always (ratio 1.00 attained).
Part 2: Lemma 1 on all 4096 cuts x 105 witness maps x 2 points x 8 designs (6.9M cases, tight in 3·10^5);
  Theorem 3 and the partition rank lemma on 27k random (tree, Y, omega); sum_t nullity <= nR for all 729 omega.
Part 3: exact bw_t (3^L subset DP, L <= 14) for 32 (design, point) pairs and the trivial scheme:
  (b) (2/3)(tw+1) <= bw (tw(Pr F) = 4, bw = 4, DL exponent 3); (d1) w <= tw(G_M')+1 <= row-treewidth+1;
  Prop. 1(iv) bw <= rtw + s. bw ranges from 2 (points with no rows) to n (random dense).
Part 4: aggregate example, lw >= 4..12 against IR width 2–3.  Part 5: balanced edge/node; anticoncentration on
cuts; value-count lemma.  Part 6: path formula, w_t = 3, 4 against IR width 2.
