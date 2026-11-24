# G11b: per-point errors -- the sqrt(log n) loss is optimal (THEOREM, general class of G6d)

## 0. Setting (G6d section 1, general = non-commutative)
J(t): 2^l x 2^l matrix polynomial over R (or C), deg_t J <= d, J(0) = X^{(x)l} (x -> complement).
Points tau_1..tau_N distinct, nonzero; point k has direction c(k) in [l]; m_c = #points of direction c.
J(tau_k) = K_k + E_k with K_k[x,y] = 0 whenever x_{c(k)} != y_{c(k)} (weak niceness), E_k arbitrary.
rho_k = rank(E_k)/2^l. Commutative class, one point per direction: rank E_k = #bad k-edges (G6d sec.3),
R_avg(l,d) = min sum_k rho_k / l, R_max(l,d) = min max_k rho_k >= R_avg.

## 1. THEOREM (dual polynomial method).  Let m = floor((N-d-1)/2) >= 0. Then
      sum_k rank(E_k)  >=  sum_{x in {0,1}^l} max(0, m + 1 - mu(x)),   mu(x) = sum_{c: x_c=1} m_c.
One point per direction (N = l, mu(x) = |x|):
      sum_k rank(E_k) >= N_m := sum_{j=0}^{m} C(l, <=j),  m = floor((l-d-1)/2),
      R_max >= R_avg >= N_m/(l 2^l) >= (1/l) Pr[Bin(l,1/2) <= (l-d-1)/2].
Proof. A_m := { row vectors a(t) in R[t]^{2^l} : a_x(t) = prod_{k: x_{c(k)}=1} (1 - t/tau_k) * r_x(t), deg a_x <= m }.
dim A_m = sum_x max(0, m+1-mu(x)).  For a in A_m: a_x(tau_k) = 0 unless x_{c(k)} = 0.            (1)
Suppose sum_k rank(E_k) < dim A_m. "a(tau_k) E_k = 0" is rank(E_k) linear conditions on a, so some
a != 0 in A_m satisfies a(tau_k) E_k = 0 for all k. If a(0) = 0 (all r_x(0) = 0) replace a by a/t (still in A_m,
same conditions since tau_k != 0); repeat until a(0) != 0.
Put c_y(t) := a_{bar y}(t) and Phi(t) := a(t) J(t) c(t)^T. Then
 * deg Phi <= 2m + d <= N - 1;
 * Phi(tau_k) = a(tau_k) K_k c(tau_k)^T + a(tau_k) E_k c(tau_k)^T = 0: the 2nd term vanishes by choice of a; in the
   1st, a_x(tau_k) != 0 forces x_{c(k)} = 0 and c_y(tau_k) != 0 forces y_{c(k)} = 1 by (1), where K_k[x,y] = 0;
 * Phi(0) = sum_x a_x(0) a_x(0) > 0   (J(0)[x,y] = [y = bar x]).   (Over C: c_y := conj(a_{bar y}).)
A polynomial of degree <= N-1 with N roots is 0, contradicting Phi(0) > 0.  QED
Commutative form (direct): w(z,t) = F(s_1(1-t/tau_1),..,s_l(1-t/tau_l), t), deg F <= m, vanishing at all bad
(z,k); W = w^2 has deg <= l-d-1, W(z,tau_k) = W(z^k,tau_k), W = 0 on bad edges, and
sum_z W(z,t) P_z(t) has degree <= l-1, l roots, value sum_z w(z,0)^2 > 0 at 0. (Smolensky-type squaring trick;
the certificates y_{z,k} = W(z,tau_k)/(tau_k prod_{j!=k}(tau_k - tau_j)) are certificates in the sense of G8d sec.1.)
Remarks. (i) The bound is the per-point analogue of G4c's Hamming-ball lemma (common E: rank >= C(l,<=m)); the
per-point relaxation gains at most a factor ~ 4 d^2/l (see 3). (ii) Tight at d = l-1 (N_0 = 1 = T(l,l-1)).
(iii) Weak for fixed d: R_avg >= ~0.2/sqrt(l) only; the uniform bound 2^{-d}/(d+2) of G8d is better there.
(iv) Equal multiplicities m_c = mu: R_avg >= (1/l) Pr[Bin(l,1/2) <= M-1], M = floor(m/mu) >= (l - (d+1)/mu)/2 - 1,
i.e. d is replaced by d_eff = d/mu = d l/N (the same quantity that controls h >= n d/N).
Unequal multiplicities: the displayed inequality holds, its asymptotic evaluation is NOT done.

## 2. Consequences
(a) Degree threshold. If max_k rho_k <= l^{-C} (C > 1/2 fixed) then d >= (1-o(1)) sqrt((2C-1) l ln l).
    [N_m ~ C(l,m) (l/2d)^2, C(l,m) 2^{-l} ~ sqrt(2/(pi l)) e^{-d^2/2l}; hence R_avg >~ e^{-d^2/2l} sqrt(l)/(c d^2).]
    G7a's layered construction (one central block) has rho_k ~ sqrt(2/(pi l)) e^{-d^2/2l}, i.e. reaches
    rho <= l^{-C} at d = (1+o(1)) sqrt((2C-1) l ln l). So d*(l, l^{-C}) = (1+o(1)) sqrt((2C-1) l ln l), and
    R(l,d) = e^{-d^2/2l} * l^{-1/2} up to a factor O(d^2/l) = O(log l) in the window d = Theta(sqrt(l log l)).
    In particular (I') of G8d is true in the strong form and (II) is not needed.
(b) Running time in the cost model of G6d sec.1 (exponent >= max(n/q, rho_max M log n, h log(n/h)), M >= n,
    l <= log n/(2 log q), h >= n d/N; the h-term read as monotone in h, 2^{Omega(n)} for h >= n/4):
    exponent = Omega(n (loglog n)^2 / sqrt(log n)).
    Proof sketch: if the exponent is <= eps n (loglog n)^2/sqrt(log n), then log q >= (1/2-o(1)) loglog n,
    l <= (1+o(1)) log n/loglog n, rho_max <= (loglog n)^2/(log n)^{3/2}; by the Theorem either d >= l/4
    (h >= n/4) or (d+1)^2/l >= c loglog n, so x = d/l >= c' loglog n/sqrt(log n) and
    h log(n/h) >= n x log(1/x) >= c'' n (loglog n)^2/sqrt(log n). Contradiction for small eps.
    This matches OpenAI's interpolation route 2^{O(n (loglog n)^2/sqrt(log n))} (G4c sec.1) up to the constant.
    Valid for one point per direction and for equal multiplicities (1(iv)), commutative or not.

## 3. Numerics (dual_check.py, bounds_table.py; outputs *_out.txt)
* Test 1: layered feasible configurations, l <= 9, d <= 5: the evaluation map W_m -> R^{bad} is injective (as the
  theorem demands) in all 29 cases (rank mod 2^31-1).
* Test 2: (6,1), (7,2), (6,3): random bad sets of size dim W_m - 1: explicit certificate y built from w^2 and
  verified exactly over Q (edge-consistent, zero on bad edges, vertex moments 1..d vanish, sum y != 0);
  primal system infeasible (mod p). 8/8.
* Test 3: general matrices, (3,0),(4,1),(3,2): a(tau_k) K_k (aX)^T = 0 and Phi(0) = |a(0)|^2, exact, 60/60.
* All known exact values respect the bound: T(4,1)=12>=6, T(5,1)>=28>=7, T(5,2)>=12>=7, T(6,3)>=10>=8.
* UB/LB (layered rho / theorem) ~ 4 d^2/l: e.g. l=1000: d=83: 31.5; d=118: 59.5; d=144: 82.2.
* Threshold degree for rho <= l^{-C}, theorem vs layered: C=2: l=100: 26 vs 36; l=1000: 110 vs 142; l=4000: 252 vs 312;
  C=3: l=4000: 354 vs 404. Ratio -> 1 slowly (relative gap ~ lnln l / ln l from the O(log l) factor).

## 4. Self-review
* Theorem: complete proof, 10 lines, uses only weak niceness, J(0) = X^{(x)l}, deg <= d, distinct nonzero points,
  and positivity of sum a_x(0)^2 (field R or C; FAILS over finite fields). Not yet re-read by a second agent,
  not in Lean (it should be easy to formalise: same ingredients as PerPoint.lean).
* 2(a): asymptotics of N_m by standard binomial estimates, checked numerically, not written out with error terms.
* 2(b) depends on the cost model (output rank = M rho_k charged n^{O(M rho_k)}, deg P <= n d, h-term monotone);
  it is a statement about this accounting, not about all algorithms. Unequal multiplicities: not evaluated.
* What is NOT shown: optimal constant factor of R (gap O(log l)); small-d regime (R(inf,d) between
  2^{-d}/(d+2) and 1/(2(d+1)) remains open but is irrelevant for the running time, because rho <= l^{-C} is needed).
