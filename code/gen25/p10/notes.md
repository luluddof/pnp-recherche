# G25-P10 — synthesis: combinations of proved pieces (working notes)

Tags: [P] short proof written here (unreviewed), [P-rel] proof here relative to inputs whose status is given,
[C] conditional on a preprint (or on the "network form" of Claude's Thm A, see S4), [N] own numerics,
[H] heuristic / indicative only, [O] open. Nothing here was written into /home/user/pnp-recherche.
Files here: notes.md; leafscan.py + model2.py/appd_explicit.py/run5.py (copies of P3's model, model2 gets a
`leaffac` switch); window.py (coset-window table); pp_check.py (check on the exact (12,3)@6 point).

Ranking by value (details below):
| # | combination | status | value | goes to |
|---|---|---|---|---|
| 1 | Article 3: DF lower bound A'' (Lean) + explicit App. D (P3) = explicit sandwich for the "sqrt(log) is optimal" claim | P-rel | medium (paper in register A) | Article 3 |
| 2 | explicit constants in the conditional consequences (Max-Cut, Max-3-SAT spectrum, binary #CSP/T1, l-cycle / longest cycle) | P-rel + C | medium-low | B1, B5, B6 |
| 3 | coset window (G24f) = one slice of the (p,q) trade-off; "open window => point strictly inside the circle" | P-rel (dilute bookkeeping of G24f) | medium-low | B2 / coset track |
| 4 | Theorem 1 + flag form + rank perturbation: "pencil + low rank" bound | P (weak form Lean-ready) | low-medium | B2 (endpoint section) |
| 5 | sweep idea in App. D | fails structurally (leaves only) [N] | low | B6 remark |
| 6 | mod-3 evaluator + explicit App. D | no new statement (dominated) | low | Article 2 remark |
| 7 | distance invariant + Theorem 1 (q = 2) | fails; one ingredient (unitary QZ gauge) found | low | B2 remark |

---------------------------------------------------------------------------------------------------------
## 1. Article 3: explicit sandwich for App. D's root  [P-rel]

Inputs. (i) G18e/G18s Theorem A'' (Lean, DFRank.lean) and its corollary (P+): at ONE call with M wires and u
directions, on App. D's labelling (Lemma D.3 loads + (D.6)), any group-factored identity enumerated term by term has
   log2 N >= M (C(2u,u)/4^u - r^{-a}) = (1/sqrt(pi) - o(1)) M/sqrt(u) = (0.564 - o(1)) M/sqrt(u).
(ii) P3: the per-level cost of App. D is (pi/ln2) sqrt(2 nu kappa/3) M/sqrt(u) bits (the (D.15) bound), with
nu = modes factor, kappa = precision factor (P3 code: E += PIB*sqrt(2/3*nufac*best)*M, best ~ kappa/u).
Same M, same u (u = number of label directions of the call, P3 = (D.5)'s u_d; G18s sec. 4 identifies the App. D
root identity as a DF identity with k = (t, zeta, beta)).

**Proposition 1.** For a call of App. D with M wires and u directions (u large), the exponent of the root
enumeration lies between 0.564 M/sqrt(u) (every enumerative DF scheme on this labelling) and c_call M/sqrt(u):
| design | nu | kappa | c_call | ratio upper/lower |
|---|---|---|---|---|
| paper as written (8hD, h = 3M/u) | 6 | 24 | 44.4 | 79 (G18s said "about 70") |
| P3 row C2 (rigorous, (D.6) unchanged) | 6 | 6 | 22.2 | 39 |
| P3 rows D-G (rigorous given App. D; packing slack) | 2 | 6 | 12.8 | 22.7 |
| P3 row H (Smith modes, conditional) | 1/2 | 6 | 6.4 | 11.4 |
| P3 row I (kappa = 2, open) | 1/2 | 2 | 3.7 | 6.6 |
Whole algorithm (r = 64): (D.6) forces u <= a <= log2 M/(4 log2 r), so every enumerative scheme of G18s's family F
pays >= 0.564 sqrt(4 log2 r) M/sqrt(log2 M) = 2.76 M/sqrt(log2 M) at the root, while P3's explicit total is
<= 75 M0/sqrt(log2 s) (M0 = km wires): global gap factor 27 (r = 64). For any even r >= 2 the lower constant is
>= 0.564*2 = 1.13.
Proof: divide the two bounds; c_call = 4.532 sqrt(2 nu kappa/3). QED.
Caveat (exact scope): the lower bound needs (D.6)-type near-uniform loads. P3's rows D-G soften (D.6) (packing slack
P = 1 + D_arr sqrt(S2)/M + L D_arr/M); the lower bound then holds with C(2u,u)/4^u - O(load deviation), not re-derived
here [O, short]. Rows up to C2 satisfy (D.6) as printed, so the ratio 39 is fully covered.
Value: replaces "the gap is a factor ~70 (G18s)" by an explicit 23-39 at the call level and 27 globally; states
precisely what Article 3 can claim: "sqrt(log) optimal within a factor 23 (rigorous designs) / 11 (Smith-conditional)
among enumerative DF schemes with App. D's labelling". Also locates the remaining loss: nu (modes) and kappa
(precision) only; the factor sqrt(4 log2 r + 1) x 1.17 in P3's 75 is a design choice (u vs lambda), not an
enumeration loss.

---------------------------------------------------------------------------------------------------------
## 2. Explicit constants in the conditional consequences  [P-rel + C]

Input: P2's Theorem A': time 2^{(c1+o(1)) W/sqrt(log2 m)}, c1 = 10.6 (deterministic, Lemma T4 top) or 7.47 (Las
Vegas), W = number of 1-bit equality links ("tests"), m = number of nodes (clauses). c1 is a limsup (liminf c1/sqrt2).
**Assumption (NF, "network form")** [C]: Theorem A' holds verbatim for networks whose nodes carry O(1)-size (or
2^{O(Delta b)}-entry) integer tables with B-bit entries, at poly(B) extra cost, W counting the bits on links between
distinct nodes. Reason it is plausible: the paper's objects are networks of tables (Def. 2/3), P2's Lemma S is stated
for arbitrary integer tables with |W_w| < 2^B, and the Lemma T4 variance bound depends only on P/U0, which is
unchanged when links carry b bits (P and U0 both scale by b). G21c T1 used the same assumption. Not re-proved.

(2a) **Max-Cut spectrum.** Graph with n vertices, m edges, max degree Delta = O(1). Vertex nodes; for each edge uv
(u < v) one 1-bit link carrying x_u to v's node; v's table multiplies z^{[x_u != x_v]}. With z = 2^{n+1},
sum_x prod = sum_t N_t z^t gives every N_t = #cuts of size t. W = m:
     time 2^{(c1 + o(1)) m/sqrt(log2 n)} poly(n)        (c1 = 10.6; 7.47 Las Vegas)   [C: NF + Claude's paper]
(G9e C4 had 2^{O(m/sqrt log m)}.) Best known for sparse Max-Cut is about 2^{m/5.2} [from memory]; the asymptotic
crossover c1/sqrt(log2 n) < 1/5.217 is log2 n > 3.1e3 [H: ignores the o(1); P2 shows lower-order terms move
thresholds by large factors].
(2b) **Max-3-SAT spectrum at the price of #3-SAT.** Clause tables weighted z^{[clause satisfied]}, z = 2^{n+1}: no
counter circuit needed, W = 3m - n exactly as for #3-SAT:  time 2^{(c1 + o(1))(3m - n)/sqrt(log2 m)} poly(m) for
all N_t (#assignments satisfying exactly t clauses). [C: weighted clause tables with B = O(mn) bits; the sweep
lemma of P2 is already stated for B-bit tables.] Without weights (pure 3-CNF + adder tree with 2^a clauses merged
per leaf) the extra links are <= (2a+4) m/2^a, i.e. (1+eps) W for every eps > 0.
(2c) **Binary #CSP (T1 explicit).** v variables, domain s, constraint graph H with |E| edges, max degree Delta.
Vertex nodes hold x_u; edge (u,v), u < v: one b-bit link (b = ceil log2 s) carrying x_v to u's node, which applies
f_uv. W = |E| b (better than chaining occurrences, (2|E| - v) b, for Delta >= 3):
     time 2^{(c1 + o(1)) |E| ceil(log2 s)/sqrt(log2 v)} poly(v, s, B)   [C: NF]
Delta-regular: s^{(c1 Delta/2 + o(1)) v/sqrt(log2 v)} (s a power of 2). Beats the treewidth DP s^{tw+1} on a class
with tw >= beta v once log2 v > (c1 Delta/(2 beta))^2; e.g. cubic expanders with beta = 0.1: log2 v > 2.5e4 [H].
(2d) **l-cycle, Hamiltonian cycle, longest cycle / path (B1 explicit).** Network version of G9e Lemma 1 with EDGE
colours: vertex nodes; edge e = uv carries one (b+1)-bit link (x_e, D_e), D_e in [p] (b = ceil log2 p); v's table:
deg_X(v) in {0,2} (cycles) and the chosen edges at v have equal colours, unchosen edges have D_e = 0, weight w_e on
chosen edges (one endpoint), plus a counter of sum_e x_e (adder tree, (eps + 1) m extra link bits). For a disjoint
union of cycles the number of colourings is p^{#cycles}, so #/p mod p = sum over l-cycles of prod w_e mod p
(as in G9e Lemmas 1-3; no isolated-vertex issue because only edges are coloured; for paths use vertex indicators
y_v = [deg_X(v) >= 1] so that no component is an isolated vertex). W = (b + 2 + eps) m. Trials as in G9e Thm 2:
log2 R <= 1.4427 m/(p - 1) + 2 with p prime in (2^{b-1}, 2^b], b = ceil(log2 sqrt(log2 m)):
     decision, one-sided error:  2^{(c1/2 + o(1)) m log2log2 m / sqrt(log2 m)},   c1/2 = 5.3 (3.74 with LV inner)
     #l-cycles mod a fixed prime p (deterministic): 2^{(c1 + o(1)) (ceil log2 p + 2) m/sqrt(log2 m)}
for max degree Delta = O(1) [C: NF + Claude's paper]. Unbounded degree: split vertices into chains of degree-3
nodes (links carry colour-so-far + count): W_b <= b (m + sum_v (deg v - 3)^+) <= 3bm, constant <= 3c1/2 = 15.9.
Pure 3-CNF route (Tseitin + comparators "value <= p-1", weight blocks), no NF needed: W_b <= ~26 b m, constant
~13 c1 = 138: the encoding, not the algorithm, decides the constant (26x).
Indicative crossover vs Bjorklund 1.657^n on cubic graphs (m = 1.5n): log2 m ~ 6.5e4 [H, asymptotic c1 only].
Value: B1 becomes a statement with an explicit constant; (2a)-(2c) are short and natural; all inherit "conditional on
an unverified preprint" and need NF to be checked against the paper's proof (the one real gap).

---------------------------------------------------------------------------------------------------------
## 3. The coset window as one slice of the trade-off  [P-rel; window.py]

Input (G24f sect. 1, 3): dilute regime (fixed E0, H = ha + hb, n -> oo): N = n^{kappa+o(1)}, cost
4 rho = (2E0 + H)/kappa, labels force kappa <= E0, and the Borcherds split gives an EXACT decomposition
T = P_a + P_b of T(x,v,u,y) = delta_xv delta_uy (= the swap: rank N^2 in two flattenings, 1 in the third) with
rank P_a <= D ~ n^{E0+hb}, rank P_b <= D' ~ n^{E0+ha-1} in the two other flattenings. So the window produces feasible
pairs (p,q) = (N^theta, N^T) with theta = (E0+hb)/kappa, T = (E0+ha-1)/kappa.

**Proposition 3.** (i) 4 rho = theta + T + 1/kappa (identity: 2E0 + H = (E0+hb) + (E0+ha-1) + 1).
(ii) [Theorem C] theta + T >= 2, hence 4 rho >= 2 + 2/(2E0 + H - 1) and 4 rho >= 2 + H/E0 (labels).
(iii) Put a = E0+hb, b = E0+ha-1. If (a - b)^2 <= 2(a + b) + 1, then 4 rho < 2 sqrt2 forces (theta, T) STRICTLY
inside Appendix B's circle theta^2 + T^2 = 4. In particular for every H = 0 case, and for every integer window case
(H < 0.83 E0) with E0 <= 7 (E0 = 8..10: all but the most lopsided split hb = H, ha = 0; table window.py).
(iv) Conditional: any bound 2 max(theta,T) + min(theta,T) >= 4 (true in the ST monomial model, Lean; = the
"spreading conjecture" of G25-B2) closes the H = 0 window for E0 <= 5 (at H = 0, 4 rho = 2 theta and
theta >= 4E0/(3E0 - 1), so 4 rho >= 8E0/(3E0 - 1) >= 2 sqrt2 iff E0 <= 5.83) and some H > 0 cases (table), not all.
(The ST bound itself does not apply: the coset decomposition is not monomial.)
Proof. (iii): 4 rho = theta (2E0+H)/a, so 4 rho < 2 sqrt2 iff theta < theta* := 2 sqrt2 a/(a+b+1). On the ray of
slope s = b/a the circle is at theta_c = 2a/sqrt(a^2+b^2). theta* <= theta_c iff 2(a^2+b^2) <= (a+b+1)^2 iff
(a-b)^2 <= 2(a+b)+1. H = 0: a - b = 1. Integer cases: |a - b| <= H + 1 and (H+1)^2 <= 4E0 + 2H - 1 iff H^2 <= 4E0 - 2,
which holds for all H < 0.83 E0 when E0 <= 7. (iv): with b <= a the bound gives theta >= 4a/(2a+b); compare with
theta* (window.py also does b > a). QED
Meaning. A coset construction beating 2 sqrt2 would itself be a feasible family strictly inside the circle at slope
(E0+ha-1)/(E0+hb) - no such point is known at any slope in (0,1) (App. B is ON the circle, Hrubes pieces outside,
monomial tilings at exponent 2). So the window is not an independent escape route: for E0 <= 7 it is contained in
the main construction problem (single-slope version), and "circle optimal at slope s" would close it. On the
lower-bound side the H = 0, E0 = 2 slice needs only theta + T >= 3/sqrt2 = 2.12 at slope 1/2 (6% beyond the line).
Not covered: finite-u interior of the window (only the dilute corner), and a construction outside G24f's split.

---------------------------------------------------------------------------------------------------------
## 4. "Pencil + low rank": Theorem 1 with a rank perturbation  [P]

Setting: K = sum_{l<=r} E_l (x) C_l, rank in (x,y)|(u,v), Swap[x,y,u,v] = [x=v][y=u]; E-span 𝓔 (dim d).
For a 2-dim subspace V of 𝓔 and a complement W (𝓔 = V + W) put rho(V,W) = dim sum_{w in W} im(w)
(= rank of the concatenation of a basis of W). rho_2(𝓔) := min over (V,W).

**Proposition 4.** For every field:
 (weak)  rank(Swap - K) >= N(N-1)/2 - N rho_2(𝓔);
 (flag)  rank(Swap - K) >= sum_{t=1}^{N-1} (N - t - min(rho, (d-2)t))_+ >= (N - rho)(N - rho - 1)/2, rho = rho_2
         (V, W taken over the algebraic closure);
 and the same with the C-span (conjugation by Swap).
Proof (weak). Choose a basis B1, B2 of V and B3..Bd of W; rewrite K = sum_k B_k (x) C'_k. K_V = B1(x)C'_1 + B2(x)C'_2 has
osr <= 2, so rank(Swap - K_V) >= N(N-1)/2 (Theorem 1, Lean: swap_sub_two_kron_rank', any field). K_W = K - K_V has image
inside (sum_W im B_k) (x) K^N, so rank K_W <= rho N. Subadditivity. QED (Lean-ready: three existing lemmas.)
Proof (flag). Lemma P (Lean) makes B1, B2 upper triangular by an equivalence (X, Y), which preserves ranks, rho and
Swap. Then F_t = span(e_1..e_t) satisfies B_i F_t in F_t (i = 1,2), so 𝓔F_t in F_t + sum_W w F_t, of dimension
<= t + min(rho, (d-2)t). Insert in G24g's flag form (sect. 6, [P] there): rank >= sum_t (N - max(t, dim 𝓔F_t)). QED
Corollaries. d = 2: Theorem 1. rho >= N: the G24g bound sum floor((b-1)/(r-1)). d = 3: W = span(E') for ANY nonzero
E' in 𝓔, so rho_2 = min rank of a nonzero element: a decomposition (p, 3) has every nonzero element of the E-span and
of the C-span (over the closure) of rank >= N - (1 + sqrt(1 + 8p))/2; a gain g below N(N-1)/2 at q = 3 needs all
those ranks >= g/N (weak form, over K). Exponent form: p <= N^{2-delta} at any q forces rho_2 >= N - sqrt2 N^{1-delta/2}
on both sides ("the terms beyond a pencil must have almost full joint image").
Check (pp_check.py, exact): on G24g's (12,3)@6 point both 3-spans are singular pencils (det == 0), min rank found 4:
bound (6-4)(5-4)/2 = 1 <= 12, consistent. Honest limit: for generic spans rho_2 = N and the statement is empty; it only
constrains decompositions whose E- or C-span contains low-rank directions (e.g. monomial/ST ones with small classes,
or "identity + sparse correction" ansaetze near the endpoint).

---------------------------------------------------------------------------------------------------------
## 5. Sweep idea in App. D  [N, leafscan.py]  -> fails structurally
Claude's gain came from the BASE level (e0 multiplies every level). App. D has no clique base: its recursion is
root-dominated (mode enumeration, untouched by any contraction order); the only brute-force contraction is at the
leaves (2^{M_leaf}). In P3's row-G model (threshold vs 2^n, density 4.26):
   leaf cost x1: lambda* = 8.97e5 | x0.5: 8.86e5 | x0.25: 7.83e5 | x0 (free leaves): 6.35e5.
Provable universal leaf factor: a random vertex order has every cut <= (n/(2(n-1)))M + O(Delta sqrt(n log n)) (McDiarmid
for permutations), so contraction in 2^{(1/2+o(1))M_leaf}: -1.2% on lambda*. Stronger leaf gains need App. D's leaf
structure (e.g. enumerating equality-vertex values on incidence-graph leaves) [O]. Ceiling -29% (free leaves); the
asymptotic constant 75 is unchanged (leaves vanish as depth grows). Reverse direction (App. D at Claude's base) is
useless: App. D only gains at astronomically large sizes.

## 6. Mod-3 evaluator (G18b) + explicit App. D (P3)  -> nothing new
Bounded-occurrence 3-CNF (m <= n): mod-3 route (12.6+o(1)) n lnln n/sqrt(ln n) nats (decision/mod 3; G18b numerics
10-11 at these sizes; non-rigorous continuous relaxation) vs App. D <= 74 * 3 m/sqrt(lambda) bits. Bits per n:
lambda = 1e5: 0.61 vs 0.70; 1e6: 0.25 vs 0.22: the mod-3 route is (non-rigorously) ahead for log2 n below ~4e5, both
beat 2^n from ~3-5e4. Dominated anyway by Claude's explicit Theorem A' (density 1: log2 m* = 151). No statement worth
a paper line beyond this remark (Article 2, "upper bound side").

## 7. Distance invariant + Theorem 1 (q = 2)  -> fails; one ingredient
Ingredient [P]: Lemma P has a UNITARY version over C: the generalized Schur (QZ) form exists for every square pencil,
singular or not (induction: a common direction x with Ax = lambda Bx always exists). So for q = 2 one may assume
E1, E2 upper triangular by a unitary gauge (X(x)Y)(.)(Y*(x)X*), which preserves Swap and the Frobenius norm; the gauge
group is compact, so in a border sequence the gauges converge.
Obstruction [P]: the Omega-compression (1 - Theta)[Omega,Omega] is unitriangular, but unitriangular matrices can be
arbitrarily close to rank-deficient ([[1,t],[0,1]] is at distance ~1/t from rank 1; [[1,t],[1/t,1]] has rank 1 and
lower part -> I). So neither a distance bound nor the border version of Theorem 1 follows from the pattern alone;
one must show that Theta's strictly-upper Omega-block cannot grow while A compensates (the G24c "B hidden behind A"
regime again). Not done. G24c's product-SVD class does not contain the triangular pencils.

## 8. Small extra remark (B6)
Same instances (random 3-CNF, density 4.26, vs 2^n): explicit Claude log2 m* = 9.4e3 (det) vs explicit OpenAI
log2 s* = 9.0e5: factor ~96; asymptotic exponents 10.6 (3r-1) n/sqrt(log2 m) vs 75*3r n/sqrt(log2 s): ratio 7.7 at r = 4.26
(the threshold ratio is about the square of the constant ratio, plus lower-order terms).

Run log: leafscan.py 63 s; window.py, pp_check.py < 1 min. No process left running.
