# G14d: multi-point per-point-error theorem for every placement; factor-2 improvement; rigorous asymptotics

Setting (G6d/G11b, general = non-commutative class). F any field. J(t) in F[t]^{2^l x 2^l}, deg_t <= d,
J(0) = X^{(x)l} (J(0)[x,y] = [y = ~x]). Points tau_1..tau_N in F, distinct, nonzero; point k has a direction c(k) in [l];
m_c = #points of direction c (m_c = 0 allowed), mu(x) = sum_{c: x_c=1} m_c.
B_k := J(tau_k)[ {x: x_c(k)=0}, {y: y_c(k)=1} ]   (the "01 off-nice block"), r_k := rank B_k.
If J(tau_k) = K_k + E_k with K_k weakly c(k)-nice (K_k[x,y] = 0 when x_c != y_c) then rank E_k >= r_k
(B_k is a submatrix of E_k). rho_k = rank(E_k)/2^l.

## 1. THEOREM 1 (improved dual bound; any field; any placement)            [PROVED here; NOT in Lean]
      sum_k rank(E_k) >= sum_k r_k >= sum_{x in {0,1}^l} max(0, N - d - 2 mu(x)).
G11b had sum_x max(0, floor((N-d+1)/2) - mu(x)): Theorem 1 is 2x stronger (exactly 2x for N-d even), needs no
positivity (so no ordered field: G11b's remark "fails over finite fields" concerned its proof only), and uses only
the 01 blocks. One point per direction: sum_k rank E_k >= sum_j C(l,j) max(0, l-d-2j) = 2^l E[(S_l - d)_+],
S_l = sum of l independent +-1 signs.  (Lean's SqrtLog.sum_rank_ge_N is the weaker G11b form N_m.)

Proof. D := N-1-d (if D < 0 nothing to prove). pref_x(t) := prod_{k: x_c(k)=1} (1 - t/tau_k) (degree mu(x), pref_x(0)=1).
A_j := { a = (a_x)_x : a_x = pref_x * r_x, r_x in F[t], deg a_x <= j }  (a_x = 0 if mu(x) > j).
  (i)  dim A_j = sum_x max(0, j+1-mu(x)) = sum_{i=0}^{j} #{x : mu(x) <= i}.
  (ii) a in A_j  =>  a_x(tau_k) = 0 unless x_c(k) = 0.
V_j := { a in A_j : a(tau_k)|_{x_c(k)=0} * B_k = 0 for all k }: at most sum_k r_k linear conditions, so
  (iii) sum_k r_k >= dim A_j - dim V_j.
Pairing. Let a in V_j, b in A_{D-j}, and Phi(t) := sum_{x,y} a_x(t) J(t)[x,y] b_{~y}(t). Then deg Phi <= j+d+(D-j) = N-1.
At t = tau_k (c = c(k)): a_x(tau_k) != 0 forces x_c = 0, b_{~y}(tau_k) != 0 forces (~y)_c = 0 i.e. y_c = 1, so
Phi(tau_k) = a(tau_k)|_0 B_k b~(tau_k)|_1 = 0. N roots, degree <= N-1: Phi = 0. Hence
  0 = Phi(0) = sum_x a_x(0) b_x(0).
For every x with mu(x) <= D-j the vector b = (pref_x at x, 0 elsewhere) lies in A_{D-j}; so
  (iv) a in V_j  =>  a_x(0) = 0 whenever mu(x) <= D-j;   and a_x = 0 whenever mu(x) > j.
So the map V_j -> F^{2^l}, a -> a(0), has image supported on the band {x : D-j < mu(x) <= j}. Its kernel is t*V_{j-1}
(a(0)=0 => all r_x(0)=0 => a/t in A_{j-1}; the conditions survive since tau_k != 0). With V_{-1} = 0, induction gives
  (v) dim V_j <= sum_{i=0}^{j} #{x : D-i < mu(x) <= i}.
(i),(iii),(v) with j = D:  sum_k r_k >= sum_{i=0}^{D} #{x : mu(x) <= min(i, D-i)} = sum_x max(0, D+1-2mu(x)).  QED
(G11b = the case j = floor(D/2), where the band is empty and V_j = 0.)
Checks: l=1 (all N points on one direction): bound = N-d = truth (J[0,1] has <= d roots). d = N-1: bound 1 = truth.

## 2. THEOREM 2 (the adversary's best placement is the balanced one)        [PROVED]
For all placements (m_c) of N points on l directions (unequal, zero, N < l or N >> l allowed):
      (1/N) sum_k rho_k  >=  (1/l) E[(S_l - delta)_+],     delta := d l / N  (real),
hence max_k rho_k >= (1/l) E[(S_l - d l/N)_+]. The exact integer minimum of the Theorem-1 bound over placements
is attained at the balanced placement m_c in {floor(N/l), ceil(N/l)}.
Proof. f(a) := sum_x max(0, N-d-2 sum_c a_c x_c), a in R^l, is convex (max of affine functions, summed) and invariant
under permutations of the coordinates of a (permute x). So f(a) = (1/l!) sum_pi f(pi a) >= f(mean) = f(N/l,..,N/l)
= sum_x max(0, N-d-(2N/l)|x|) = 2^l (N/l) E[(l - delta - 2|x|)_+] and l - 2|x| ~ S_l. Convex + symmetric = Schur-convex,
and the balanced integer vector is majorised by every integer vector of the same sum.  QED
Consequences.
* No loophole from the placement: concentrating points on few directions, or leaving directions empty, makes the
  bound LARGER, not smaller (variance of sum_c m_c eps_c grows). Directions without points give no constraint, but the
  points that carry the constraints then crowd the other directions.
* The only quantity that matters is delta = d l/N = d/(average multiplicity). N <= d gives the bound 0, correctly:
  with N <= d points J can interpolate arbitrary nice values (rho = 0), and then the cost is in h >= n d/N >= n.
* Tight up to the same factor as before: the multi-point layered family g_z = chi(z) + p_B(sum_j s_j pi_j(t)) - p_B(sum_j s_j),
  pi_j = prod_{k on j}(1 - t/tau_k), deg p_B = D, has deg_t = D * max m_c and rho_k = layered(l, D): with equal
  multiplicities delta = D exactly. So multi-point = one-point with d replaced by delta; nothing is gained or lost.
* G6d's corollary (b) (N <= 2 d l/log(1/2rho)) is superseded: it is the regime delta < 1 of Theorem 2.
* Copy-dependent J_c (different J in each copy of the cube): the bound is on sums, so it holds for each copy and for
  the total output error rank charged by the model.
* Complex points are fine (Theorem 1 needs no conjugation).

## 3. Explicit binomial estimate and the degree threshold                    [PROVED, constants not optimised]
(*) For l >= 1, real delta >= 0, real w >= 2 with delta + 2w <= l/sqrt 2:
      E[(S_l - delta)_+] >= (w^2/4) (2l)^{-1/2} exp( -(delta+2w)^2/(2l) - (delta+2w)^4/(6 l^3) ).
Proof. The interval [delta+w, delta+2w] contains >= floor(w/2) >= w/4 integers k = l mod 2; each contributes
(k-delta) Pr[S=k] >= w Pr[S=k]. Pr[S=k] = C(l,(l+k)/2) 2^{-l} >= (2l)^{-1/2} exp(-l D(k/l)) by
C(n,j) >= sqrt(n/(8j(n-j))) e^{n H_e(j/n)} (0<j<n; MacWilliams-Sloane Ch.10 Lemma 7; QUOTED, standard) and j(n-j) <= n^2/4,
with D(s) = ((1+s)ln(1+s)+(1-s)ln(1-s))/2 = sum_{i>=1} s^{2i}/(2i(2i-1)) <= s^2/2 + s^4/(12(1-s^2)) <= s^2/2 + s^4/6 (s^2 <= 1/2),
D increasing in |s|.  QED   (grid check l <= 3000: true/(*) >= 12, adv_asym_out.txt D1.)
COROLLARY 3 (threshold). Put u = delta/sqrt(l), 1 <= u <= sqrt(l)/5, w = l/delta in (*):
      (1/N) sum_k rho_k >= c0 l^{-1/2} u^{-2} exp(-u^2/2 - 13.5 u^4/l),   c0 = 1/(4 sqrt2 e^4)  (> 1/309).
[(u+2/u)^2/2 <= u^2/2 + 4 and (u+2/u)^4 <= 81 u^4 for u >= 1; u + 2/u <= 3u <= 0.6 sqrt(l) < sqrt(l/2).]
Therefore, if max_k rho_k <= l^{-C} with C > 1/2, then either delta > l/5, or delta < sqrt(l) (possible only if
(C-1/2) ln l <= 6.3 + 13.5/l), or
      delta^2 >= (2C-1) l ln l - 2 l ln(2C ln l) - 11.5 l - 108 C^2 ln^2 l,
i.e. delta = d l/N >= (1 - O_C(lnln l/ln l)) sqrt((2C-1) l ln l). (If u^2 > 2C ln l the claim is trivial; otherwise
2 ln u <= ln(2C ln l) and 13.5u^4/l <= 54 C^2 ln^2 l / l.) One point per direction: delta = d.
Upper bound (G7a, one central block): rho_k ~ sqrt(2/(pi l)) e^{-d^2/2l} reaches l^{-C} at d = (1+o(1)) sqrt((2C-1) l ln l).
So d*(l, l^{-C}) = (1+o(1)) sqrt((2C-1) l ln l) with relative error O(lnln l/ln l), for the max AND the average error,
any field, any number of points per direction (delta in place of d).
Exact thresholds (adv_asym_out.txt D2), C=2: l=100: G11b 26, new 28, layered 36; l=1000: 110, 116, 142; l=4000: 252, 262, 312.

## 4. THEOREM 4 (running time, all hypotheses explicit)                       [PROVED as a statement about the cost functional]
Parameters: n; integers q >= 2, l >= 1; M; d; N points with any placement; J as in section 0 (or copy-dependent);
multiplicities h_k >= 1. Cost functional
      Ex := max( n/q ,  max_k rho_k * M * ln n ,  max_k hhat(h_k) ),   hhat(h) = h ln(n/h) for h <= n/e, n/e above.
Hypotheses:
 (H1) shape A(t) = L (I_C (x) J(t)) R', J(0) = X^{(x)l}, deg J <= d; at each used point the structure exploited is
      "weakly c(k)-nice + E_k" for one coordinate direction c(k) of the tensor structure;
 (H2) the error at point k is charged n^{Omega(M rho_k)} (output rank counted as |C| rank E_k), the blocks 2^{Omega(n/q)},
      the derivatives (n/h_k)^{Omega(h_k)};
 (H3) M >= n and l <= ln M/(2 ln q - ln 2)   [true when M >= (q(q-1))^l; implied by G6d's "M >= n, l <= log n/(2 log q)"];
 (H4) sum_k h_k >= n d   (the generic bound deg per A(t) <= n d is the one used).
Then  Ex >= (1/sqrt2 - o(1)) n (lnln n)^2 / sqrt(ln n), for every choice of q, l, M, d, N, placement, J, (h_k).
With unspecified constants in the three terms: exponent = Omega(n (loglog n)^2/sqrt(log n)); the exponent a = 2 is
exact, since O'Donnell/OpenAI's choice (one point per direction, common E, d = O(sqrt(l log l))) gives O(.) in the same
functional.
Proof. L := ln n, L2 := lnln n, T := eps n L2^2/sqrt L, eps < 1/sqrt2; assume Ex <= T, n large.
 (a) q >= n/T, so ln q >= L2/2 - 2 ln L2 = (1/2 - o(1)) L2, and by (H3) l <= (L + Lam)/((1-o(1)) L2), Lam := ln(M/n) >= 0.
 (b) x := d/N. By (H4) max_k h_k >= n x; hhat is nondecreasing, so hhat(n x) <= T: x < 1/5 and x ln(1/x) <= T/n = o(1).
 (c) By Theorem 2, T/(M L) >= max_k rho_k >= (1/l) E[(S_l - x l)_+]. Put u = x sqrt(l).
     l < 25: rho >= 2^{-l}(1-x) >= 2^{-26}, Ex >= 2^{-26} n L > T. u < 1 (l >= 25): monotonicity in delta and Cor. 3 at
     u = 1 give rho >= c/sqrt(l) >= c'/sqrt(ln M), so rho M L >= c' n sqrt(L) > T. u > sqrt(l)/5 is x > 1/5, excluded.
     Otherwise Cor. 3 with 13.5 u^4/l = 13.5 x^2 u^2 = o(u^2):
       (1/2+o(1)) u^2 + 2 ln u >= ln(M L/T) - (1/2) ln l + ln c0 >= Lam/2 + (1-o(1)) L2,
     so u -> infinity and u^2 >= (2-o(1)) L2 + (1-o(1)) Lam;
       x^2 = u^2/l >= (1-o(1)) L2 (2 L2 + Lam)/(L + Lam) >= (2-o(1)) L2^2/L       (as L >= 2 L2).
 (d) x ln(1/x) is increasing on (0,1/e): hhat(n x) >= n (sqrt2-o(1)) (L2/sqrt L)(L2/2)(1-o(1)) = (1/sqrt2-o(1)) n L2^2/sqrt L > T.
 Contradiction.  QED
What the theorem does NOT cover (the exact loopholes, all outside (H1)-(H4)):
 (L1) (H4): a degree bound deg P <= Dg << n d (as in O'Donnell App. B, deg P <= 4Rd) replaces n by Dg/d in the h-term only;
      the theorem then reads Omega(Dg' (lnln)^2/sqrt(ln)) with Dg' = Dg/d. Inside the sparse-permanent reduction nothing
      better than n d is known.
 (L2) (H2): the model CHARGES |C| rank E_k; L, R' could collapse rank, and n^{O(r)} is an upper bound for rank-r permanents.
      The theorem is about this accounting, not about all algorithms.
 (L3) (H1): structures other than coordinate-niceness for one direction per point (blocks from other partitions, a
      different tensor structure per point, points at which two different cubes/labelings are used) are not covered.
 No loophole comes from: the number of points, their placement, unequal h_k, unequal rho_k, the field, copy-dependent J,
 non-commutative J.

## 5. The log gap (Task 2): status
* Improved by exactly 2: R_avg(l,d) >= (1/l) E[(S_l-d)_+] ~ (1/sqrt(2 pi l)) (l/d^2) e^{-d^2/2l}; layered ~ 2/sqrt(2 pi l) e^{-d^2/2l}.
  Ratio layered/bound ~ 2 d^2/l (was 4 d^2/l): still O(log l) in the window. NOT closed.
* Exact commutative optima by exhaustive search (mp_check_out.txt; l <= 2 exact over Q, l = 3 infeasibility mod 2^61-1,
  witnesses exact), 73 cases incl. unequal multiplicities: Theorem 1 holds in all, and is TIGHT in many multi-point cases:
  all placements on one direction (N-d)2^{l-1}; (l=2; m=(2,2); d=2,3); (l=3; m=(2,2,2); N=6; d=4: T=2); (l=3; (3,2,0); d=3: T=4);
  (l=3; (2,2,0); d=2: T=4). Not tight for one point per direction with small l-d: T(3,1)=4 vs 2, T(4,1)=12 vs 7
  (those are governed by the Q_{d+2}/determinant phenomena).
* Diagnostic (hf_rank_out.txt): in the commutative class B is feasible iff lambda: W -> sum_z W(z,0) lies in the span of the
  edge evaluations ev_e (e in B) on the dual space W_{D'} (D' = l-d-1) [verified: true for all layered B]. Any counting
  argument is bounded by rank(ev_B). For layered B: (l,d,|B|,rank,bound) = (8,4,112,56,20), (10,5,450,225,80), (10,4,720,472,136),
  (9,3,504,390,114). So for layered sets the conditions are far from dependent enough to explain the gap: the loss sits in
  step (v) (the symbol space V_j(0) is allowed to fill the whole band D-j < |x| <= j; for layered sets it does not).
  Closing the gap needs either a bad set whose symbol spaces fill the bands (a better construction, necessarily not
  weight-symmetric), or a proof that they cannot. OPEN. Only (iv) is known about V_j(0).
* The bilinear version (a in A_i, c in A_j free) gives nothing beyond Theorem 1; restriction to subcubes gives nothing
  (ranks are not additive in the general class).

## 6. Numerics (all scripts in this directory)
* mp_check.py / mp_check_out.txt: exact commutative optimum vs old/new bound (section 5).
* adv_asym.py / adv_asym_out.txt: B: balanced placement minimises both bounds over ALL compositions, l <= 5, N <= 12, all d;
  C: genuine 2^l x 2^l matrices (multi-point layered, unequal m_c, l <= 4): sum_k rank B_k (exact over Q) >= bound, 15/15
  (not sharp there); D: (*) on a grid, thresholds, ratios, Corollary 3 against exact thresholds.
* hf_rank.py / hf_rank_out.txt: rank of layered bad sets on the dual space.

## 7. Self-review (adversarial)
* Theorem 1: proof complete and elementary; the new step is (iv)-(v). Checked against exhaustive optima (tight cases
  exist, no violation), l=1 by hand, G11c's general-class values (T_gen(3,1)=4 >= 2, T_gen(4,1) >= 10 >= 7). The general
  (non-commutative) class was tested only on commutative-origin matrices: no independent non-commutative low-rank instance.
  NOT formalised in Lean: the Lean theorem SqrtLog.sum_rank_ge_N is the weaker N_m. Formalisation should reuse `core`
  with b free (bilinear) plus the filtration argument (v); positivity (sum of squares) is no longer needed.
* Theorem 2: one-line convexity argument; verified by enumeration. It bounds the AVERAGE of rho_k over points.
* (*) depends on one quoted binomial inequality; constants deliberately crude. Corollary 3's side conditions
  (1 <= u <= sqrt(l)/5) are handled separately in Theorem 4(c).
* Theorem 4 is a theorem about the cost functional under (H1)-(H4). o(1) terms are not made explicit in n. The constant
  1/sqrt2 is meaningful only if the three terms have constant 1 and natural logs.
* For l = 3 the claim "no smaller bad set" in mp_check uses linear algebra mod 2^61-1 (not a proof over Q); irrelevant
  for the theorems, relevant only for the exactness of the listed optima.
* Unproved / open: the O(log l) factor (now 2d^2/l); R_max vs R_avg; small-d regime; (L1)-(L3).
