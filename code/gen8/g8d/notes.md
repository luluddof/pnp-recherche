# G8d (track 4): algebraic constraints beyond the hypercube-Turan relaxation

Setting = G6d/G7a commutative class, one point per direction. P_z(t) = chi(z) g_z(t), deg <= d, P_z(0) = 1;
edge (z,z^k) good iff (P_z + P_{z^k})(tau_k) = 0.  T(l,d) = min total #bad edges; rho_avg = T/(l 2^l).

## 0. Literature (checked by web search, Nov 2026)
* Alon-Krech-Szabo (SIAM J. Discrete Math. 2007): min density c_d of an edge set meeting every Q_d satisfies
  Omega(log d/(d 2^d)) <= c_d <= 4/(d+1)^2 (d odd), 4/(d(d+2)) (d even).
* Ellis-Ivan-Leader (Bull. LMS 2024, arXiv 2401.16289): the vertex-transversal density for Q_d is <= c^d
  (c < 1), "and similar bounds for the edge-Turan densities" (edge version: abstract-level statement only;
  full text not read, arxiv.org is blocked here). So the combinatorial relaxation (Barrier B) admits
  transversals of density 2^{-Theta(d)}: it can NEVER give poly(1/d). G7a's quoted O(1/d^2) is outdated.
* Balogh-Hu-Lidicky-Liu (EJC 35, 2014) and Baber: C4-free density in Q_n <= 0.6068 (later 0.603 is quoted).

## 1. Exact feasibility criterion (PROVED, Fredholm alternative)
G (good set) is feasible iff there is NO edge weighting y: G -> R with
  (i) at every vertex z: sum_{k: (z,z^k) in G} y_{z,k} tau_k^i = 0 for i = 1..d   (y_z in W_tau, dim l-d),
  (ii) sum_e y_e != 0.
The support of such a y is a "certificate": every certificate must contain a bad edge. (d+1)-subcubes with
y_{z,k} = Lagrange weight lambda_k (nodes {0} u tau_T, at 0) are the Barrier-B certificates. At every vertex
the support of y has 0 or >= d+1 edges.

## 2. d = 1: covering-space characterisation (PROVED; generic tau)
With v(z) = (-1)^{|z|} u(z): good k-edge <=> v(z) - v(z^k) = -2 (-1)^{|z|}/tau_k. So G feasible iff every cycle C
of G has sum_k gamma_k(C)/tau_k = 0, where gamma(C) in Z^l is the gain of C (edge z->z^k has gain (-1)^{|z|} e_k).
If the 1/tau_k are Q-linearly independent (generic): feasible <=> all gains vanish <=> every component of G lifts
to the Z^{l-1}-cover Lambda -> Q_l, Lambda = {a in Z^l : sum a in {0,1}}, a ~ a+e_k (sum a = 0) (deck group
{v in 2Z^l : sum v = 0}).  Unbalanced cycles that must be hit besides C4: e.g. the 6-cycle with direction word
1,0,2,0,1,2 in Q_3 and 8-cycles in Q_4 (found automatically, lazy.py). Verified against the rank test on 9000
random good sets, l = 3,4,5 (liftcheck.py: 9000/9000 agree).
For special tau (1/tau_k satisfying integer relations) more cycles are allowed: lower bounds from section 2/4
data are for GENERIC tau only.

## 3. NEW LEMMA (Q_{d+2} lemma) -- PROVED for all d and ALL distinct nonzero tau
Every (d+2)-subcube contains at least 4 bad edges.  (Combinatorial relaxation: 3 suffice.)
Tool: Psi_Q(t) = sum_{z in Q} P_z(t) for subcubes Q: deg <= d, Psi_Q(0) = |Q|, additive under splitting Q into
two halves, and Psi_Q(tau_k) = 0 if Q has no bad k-edge (k a direction of Q) [G7a Barrier B].
Proof. Q has directions S, |S| = d+2. Suppose <= 3 bad edges; n_k = #bad k-edges in Q.
* Some n_a >= 2, or (<=2 bad edges, all directions distinct): take a direction a with n_a >= 1 maximal... in
  detail: if n_a = 3 the facet Q_{a,0} has d+1 clean directions -> d+1 roots, contradiction. If n_a = 2, n_b <= 1
  for one b, all other d directions clean: each a-facet Q_{a,beta} already has d roots, so it must contain a bad
  b-edge; but the single b-edge lies in only one a-facet. Contradiction. Same argument if there are <= 2 bad
  edges in distinct directions a, b (each a-facet needs the unique b-edge).
* n_a = n_b = n_c = 1 (three distinct directions), R(t) := prod_{k in S\{a,b,c}} (1 - t/tau_k) (deg d-1).
  Each a-facet has the d-1 roots of R and cannot also vanish at tau_b and tau_c, so e_b, e_c lie in different
  a-facets, say e_b in Q_{a,0}. Then Psi_{a,0} vanishes at tau_c and the roots of R: Psi_{a,0} = 2^{d+1} R (1-t/tau_c),
  Psi_{a,1} = 2^{d+1} R (1 - t/tau_b), so Psi_Q = 2^{d+1} R (2 - t/tau_b - t/tau_c). Splitting along b instead gives
  Psi_Q = 2^{d+1} R (2 - t/tau_a - t/tau_c). Hence tau_a = tau_b, contradiction.  QED
The step that the combinatorial relaxation lacks: the SAME polynomial Psi_Q is computed through two different
splittings (additivity + rigid root structure).
Tight: T(d+2,d) = 4 (exact, d = 1..3; d = 4 see below).

### Corollaries (all tau, any l >= d+2)
For each S with |S| = d+2 the 2^{l-d-2} disjoint S-subcubes give sum_{k in S} bad_k >= 2^{l-d}. Hence
   rho_max >= rho_avg >= 2^{-d}/(d+2),
versus G6d lemma rho_max >= 2^{-d-1}/(d+1) and Barrier-B averaging rho_avg >= 2^{-d-1}/(d+1):
improvement factor 2(d+1)/(d+2) -> 2. Restriction monotonicity: rho_avg(l,d) >= T(m,d)/(m 2^m) for every m <= l.

## 4. Exact data (lazy certificate generation, lazy.py; generic tau = random tau mod 2^31-1)
Method: hitting-set MILP over certificates found so far; test the optimum algebraically; if infeasible, extract
a minimal certificate (deletion filter) and add all its Aut(Q_l)-images; stop when the MILP optimum is feasible
=> exact optimum (for generic tau, w.h.p. over the random tau). See out_*.txt.
| (l,d) | comb. Barrier B | ALL-tau proved LB (Barrier B + Q_{d+2} lemma, MILP optimal) | generic-tau exact / bounds | G6d lemma (averaged) LB |
|-------|-----|-----|-----|-----|
| (3,1) | 3 | 4 | 4 (exact) | 3 |
| (4,1) | 8 | 11 | 12 (exact) | 8 |
| (4,2) | 3 | 4 | 4 (exact) | 3 |
| (5,3) | 3 | 4 | 4 (exact) | 3 |
| (5,1) | >= 20 | 28 | 28 <= T <= 32 (32 found, optimality run unfinished); max_k bad >= 6 | 20 |
| (5,2) | 8 | 10 | 12 <= T <= 14 (MILP-certified relaxation 12; 14 from G6d) | 7 |
| (6,3) | 8 | 10 | 10 <= T <= 21 (run unfinished) | 6 |
Exhaustive check of the Q_{d+2} lemma (checkQd2.py): d = 1,2,3, every 3-edge bad set infeasible, both for random
tau and for tau_k = k+1; 4-edge bad sets exist (36 for d = 1, 336/352 for d = 2).
Earlier milpmax_*.txt files of G7a are EMPTY (runs never finished).
Logs: out_*.txt (lazy.py), out2_*.txt (lazy2.py = lazy.py + Q_{d+2} inequalities + dual bounds).

## 5. Running-time consequence
* Proved: rho >= 2^{-d}/(d+2) uniformly in l. rho <= l^{-C} then forces only d >= C log2 l - O(loglog l):
  NO new consequence beyond G6d's d-independent bound 2^{Omega(n (loglog n)^3/log n)}.
  The sqrt(log n)-optimality is NOT proved.
* Correction to the target: a bound rho(inf,d) >= c/poly(d) holding uniformly in l is FALSE (OpenAI's
  construction has rho ~ exp(-d^2/2l) sqrt(l)/d < 1/d^2 at d^2 = 2C l ln l, C > 1). A liminf-in-l statement alone
  does not translate into running time. What is needed is a quantitative trade-off, e.g.
  (I') rho(C delta^2, delta) >= delta^{-a}  and  (II) super-multiplicativity rho(l1+l2, d1+d2) >= rho(l1,d1) rho(l2,d2).
  (I')+(II) would give d >= l^{1/2-o(1)}, i.e. time >= 2^{n/(log n)^{1/2+o(1)}} (sqrt(log n) optimal up to
  (log n)^{o(1)}). Both are UNPROVED; (I) with a constant instead of delta^{-a} is false (piecewise construction).
* EIL 2024 shows that every argument that only uses "bad edges hit all (d+1)-subcubes" is capped at 2^{-Theta(d)}.
  Local algebraic constraints such as the Q_{d+2} lemma give only 2^{-d} poly(d) improvements.

## 6. Self-review
* Sections 1, 3 are complete proofs (commutative class, one point per direction). Section 3 holds for all
  distinct nonzero tau; section 2's lift characterisation needs Q-independent 1/tau_k.
* Section 4 numbers: exact optimum of the generic-tau problem with high probability (random tau mod p; a bad
  prime/tau could only make feasibility checks fail on a measure-zero set). They are not lower bounds for
  special tau. MILP solutions verified only mod p (not over Q).
* Ellis-Ivan-Leader edge version: abstract-level claim, not verified in the paper body.
* Nothing here applies to the non-commutative class.
