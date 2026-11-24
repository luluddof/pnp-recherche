# G18d: the O(log l) gap (Task A) and degree-by-cancellation (Task B)
Tags: [P] proved here (paper proof), [N] numerical (exact mod p unless said), [O] open. Nothing here is in Lean.

## A0. Dual = quadrature on a line arrangement [P]
Commutative class, one point per direction, D' = l-d-1. Lines L_z = {t}, z in {0,1}^l; L_z and L_{z^k} are glued at the
node (z,k) (t = tau_k). Wt := {(W_z)_z : deg W_z <= D', W_z(tau_k) = W_{z^k}(tau_k) for ALL edges}.
  B feasible  <=>  lambda(W) := sum_z W_z(0) lies in span{ev_e : e in B} on Wt   (ev_(z,k)(W) = W_z(tau_k)).
Proof. Fredholm (G8d sec.1): G = E\B infeasible iff exists y on G with sum_k y_{z,k} tau_k^i = 0 (i=1..d, all z), sum y != 0.
A functional on {0,tau_1..tau_l} killing degree <= d is q/omega' with deg q <= D' (omega = t prod(t-tau_k)); so
y_{z,k} = W_z(tau_k)/omega'(tau_k), sum_k y_{z,k} = -W_z(0)/omega'(0); "y lives on edges" = glued; y = 0 on B = W vanishes on B;
sum y != 0 = lambda(W) != 0. QED.   So T(l,d) = min #nodes of a (signed) quadrature rule for lambda exact on Wt.
Theorem 1 of G14d is the Gauss-type lower bound for this rule (moment form with radical: for D'=2 the Gram matrix of
(1,t,u_1..u_l) is diag(1,0,I_l), a space with 1-dim radical needs ambient dim l+3 = bound).

## A1. Pruning lemma [P] and what it gives
On one line L_z the functionals ev_(z,k) span a space of dim <= D'+1 = l-d. Hence every feasible B contains a feasible
row-basis B' with |B'| = rank(ev_B), and  T <= (l-d) * vc(B)  (vc = vertex cover number of B).
Layered sets are NOT minimal: pruned layered (prune.py, rankscan.py; random tau mod 2^31-1; primal re-check l <= 8):
 (l,d)   bound  pruned  layered | max_k: pruned / layered
 (5,2)     8     16      25     |  4 / 5
 (6,2)    16     44      60     |  9 / 10
 (6,3)     9     18      36     |  4 / 6      (previous best known 21)
 (7,3)    18     50      84     | 10 / 12
 (7,4)    10     20      49     |  4 / 7
 (8,4)    20     56     112     | 11 / 14
 (10,5)   80    225     450
 fixed D'=2 (d=l-3), l=7..12: pruned = 2l+6 = 2*bound exactly, max_k = 4 (layered: l(l+1), max_k = l).
 fixed D'=3: pruned = 6l+8 vs bound 2l+4 (ratio -> 3); D'=4: ratio ~2.8-2.85; D'=5,6: ~3.4-3.5 (l <= 12).
So for l-d = O(1) Theorem 1 is sharp up to a constant c(l-d) (2 for l-d=3) and the unpruned layered family is off by ~l/(l-d).
In the window d ~ sqrt(l log l) the lemma gives NOTHING (boundary vertices have ~l/2 < l-d bad edges; rank ~ |B|).

## A2. Small cases: bound | best proved LB | best UB | layered   (total bad edges, generic tau)
 (5,1)  14 | 28 (G8d)                      | 32 (G6d d=1 family; dual MILP 300 s also 32) | 40
 (5,2)   8 | 12 (G8d lazy, generic tau)    | 14 (G6d; dual MILP 300 s also 14)            | 25
 (6,1)  38 | 72 (restriction from T(4,1)=12) | 80 ((l-1)2^{l-2})                          | 96
 (6,2)  16 | 29 (restriction from T(5,2)>=12) | 44 (pruned layered = G6d MILP)            | 60
 (6,3)   9 | 10 (G8d)                      | 18 (pruned layered, NEW)                     | 36
 (7,3)  18 | 24 (restriction from (6,3))   | 50 (pruned layered)                          | 84
EXACT optima at these sizes were NOT obtained (G8d's lazy hitting set needs > 20 min per iteration; the big-M dual MILP
dmilp.py times out with dual bounds 6..23, and its bound is valid only for |w| <= 30). Exact values remain those of
G8d: (3,1) 4, (4,1) 12, (4,2) 4, (5,3) 4. Reading: the optimum is 1.1x-2.3x the bound and 1.25x-3.6x below layered, i.e.
closer to the bound; but at these sizes everything is governed by the 2^{-d} effects, not by the window.

## A3. Two-sided niceness [P]
Conjugating J by coordinate flips x -> x+sigma keeps J(0), the degree and weak niceness and swaps the 01 and 10 blocks
of the directions in sigma; mu is unchanged in distribution. Hence for EVERY orientation sigma
   sum_k rank B_k^{sigma(c(k))} >= bound,   in particular   sum_k min(rank B01_k, rank B10_k) >= bound  (one point per direction),
and sum_k (rank B01_k + rank B10_k) >= 2*bound. This does NOT improve the theorem: min-rank completion gives
rank E_k = max(rank B01_k, rank B10_k) exactly (G11c 1a), and in the commutative class rank B01_k = rank B10_k = #bad
k-edges. Any argument valid for both halves is capped by the commutative optimum, so the gap is not caused by using one half.
A mixed pairing (a0 K00 b0 + a1 K11 b1 must vanish) forces one orientation per point: nothing beyond the above.

## A4. Status of R(l,d) in the window
Proved: (1/l)E[(S_l-d)_+] <= R_avg <= R_max <= layered, ratio 2d^2/l (unchanged). New: sharpness up to a constant when
l-d = O(1) (numerical formulae 2l+6, 6l+8; the lemma itself is proved). The order in the window is still OPEN; no
non-layered set beating layered was found there (IRL1 on the quadrature LP, quad.py: worse than pruned layered for all
(l,d) <= (7,4)). Unproved: the closed forms 2l+6, 6l+8 for all l; c(D') bounded or growing like D'/2.

## B. Degree certified by cancellation
B1 [P] Statement. Prop. 3.2 needs sum_k h_k > deg P, P(t) = per A(t). With a certified bound deg P <= kappa*n*d, Theorem 4
 of G14d becomes Ex >= kappa (1/sqrt2 - o(1)) n (lnln n)^2/sqrt(ln n) (only step (b),(d) change: max h_k >= kappa n d/N).
 A gain needs kappa -> 0. [t^{nd}]P = per(A_d), A_d = L(I (x) J_d)R', J_d = [t^d]J.
B2 [P, sympy] Niceness + low rank do NOT force deg P = nd. l=1, s=t/tau:
   J = [[s, 1-s],[1-s, -s]]: J(0)=X, J(tau)=diag(1,-1) (E=0), per J = 1-2s (degree 1 < nd = 2);
   J = [[s(sqrt2-s), 1-s^2],[1-s^2, s(sqrt2+s)]]: d=2, J(0)=X, J(tau) diagonal, per J == 1 (degree 0 < nd = 4).
 In the actual layered family (l,d)=(3,1), n=4, merged input rows {0,3},{6,5},{1,4},{2,7}, cols {6,0},{2,4},{1,5},{3,7}:
 per A(t) == 2 although two columns have degree 1 (matching degree 2): genuine cancellation (degper2.py, symbolic in tau).
B3 [P] Cancellation is non-generic and is a choice of the designer, not of the hypotheses: adding
 c * t^{d-|x+y|} prod_{k: x_k != y_k}(1-t/tau_k) to J[x,y] for |x+y| <= d-1 changes neither J(0), nor any off-nice block,
 nor any E_k. So the top coefficients on the near-diagonal support {|x+y| <= d-1} are free; for A(t)=J(t) (and whenever the
 near support has a perfect matching with independent parameters) deg P = nd off a proper algebraic subset.
B4 [P] Commutative class, G = J X^{(x)l} a circulant with symbols gamma_w(t). If the certified bound has to hold for the
 selection inputs A(t) = G(t)[U,U] with |U| <= 2 (deg <= kappa |U| d), then gamma_0, gamma_0^2+gamma_w^2 have degree
 <= kappa d, 2 kappa d, so EVERY entry of J has degree <= kappa d: the instance is a (l, kappa d) instance and Theorem 1
 applies with d' = kappa d: no gain. Caveat: this needs the degree certificate to be uniform down to trivial 2x2 inputs;
 it does not cover a certificate that only holds for the large hard inputs, nor the non-commutative class (there only
 deg G[x,y] + deg G[y,x] <= 2 kappa d follows, and 3-cycles can cancel).
 Real field, single-F family, d even: J_d X = lead * (sum_j X_j/tau_j)^d is PSD, so per(G_d[U,U]) >= prod diag > 0 (Marcus):
 deg P = |U| d for every selection input with nonzero diagonal symbol [P modulo the quoted inequality].
B5 [N] degper.py, layered and single-F, (l,d) in {(3,1),(3,2),(4,1),(4,2),(4,3)}, inputs full cube / random selection /
 random 2-to-1 merges: deg P = n d in 34 of 38 runs; exceptions: the B2 example (deg 0) and single-F selections
 (3 of 4, 6 of 8, 6 of 8: structural, no perfect matching in the support of the top coefficient).
Verdict: deg P = Theta(nd) is generic but NOT forced; total cancellation exists in toy cases. Whether a universal J has
deg P = o(nd) on all hard inputs is OPEN; it is excluded only under B4's uniformity assumption (commutative class).

## Self-review
* A0, A1, A3, B1-B4: short proofs, not machine-checked. A0's direction "=>" was also checked by a residue computation.
* All ranks are mod 2^31-1 with random tau (generic tau, w.h.p.); l >= 9 only by the dual test.
* A2 lower bounds for (6,1),(6,2),(7,3) use G6d/G8d restriction monotonicity (commutative class, bad-edge counts) and
  G8d's generic-tau values; they are not all-tau bounds.
* The MILP runs are upper bounds only. "No better non-layered set found" is weak evidence (one heuristic).
* B4 is the only negative result on cancellation and its hypothesis is strong. B3 needs independent parameters:
  with one J shared by all cube copies and merged rows this was not proved in general.
