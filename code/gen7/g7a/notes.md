# G7a: per-point errors, uniform-in-l families, barriers

Setting = G6d, commutative (Hadamard-diagonal) class, one point tau_k per direction k.
Equivalent forms (z in {0,1}^l, s_j = (-1)^{z_j}):
  g_z(t) (deg <= d, g_z(0) = chi(z));  P_z = chi(z) g_z, P_z(0) = 1.
  edge (z, z^k) good  <=>  g_z(tau_k) = g_{z^k}(tau_k)  <=>  (P_z + P_{z^k})(tau_k) = 0.
  rho_k = #bad k-edges / 2^l  (<= 1/2).

## 1. NEW construction: piecewise-OpenAI family, rho -> 1/(2(d+1)) uniformly in l  (PROVED, verified l <= 10)
Partition weights {0..l} into blocks of <= d+1 consecutive weights; p_B in Q[s] of degree <= d with
p_B(l - 2w) = (-1)^w for w in B. For |z| in B:
    g_z(t) = chi(z) + p_B( sum_j (1 - t/tau_j) s_j ) - p_B( sum_j s_j ).
Proof. deg_t <= d, g_z(0) = chi(z). At t = tau_k the first p_B-term does not depend on s_k. For an edge
inside a block, g_z(tau_k) - g_{z^k}(tau_k) = (chi(z) - chi(z^k)) - (p_B(s) - p_B(s')) = 0 because p_B equals
parity on B. Only edges crossing a block boundary can be bad.  QED
Matrix form: J(t) = X^{(x)l} + sum_B Pi_B [p_B(Y(t)) - p_B(Y(0))], Y(t) = sum_j (1 - t/tau_j) X_j,
Pi_B = projector onto the Fourier weight block B (a polynomial in sum_j X_j); all commute.
rho_k = (1/2) Pr[Bin(l-1,1/2) in boundary set], identical for all k:
  * blocks of length d+1 everywhere: rho_k -> 1/(2(d+1)) as l -> infinity for fixed d (exact values
    l <= 10 in piecewise_out.txt, all crossing edges actually bad, no failures inside blocks).
  * one central block: rho_k ~ sqrt(2/(pi l)) exp(-d^2/(2l)) (two boundary levels only), versus OpenAI's
    common error, whose rank is the whole tail Pr[|w - l/2| > d/2] ~ (sqrt(l)/d) exp(-d^2/2l): gain factor ~ l/d.
This refutes the worry in G6d that rho might tend to 1/2 for fixed d >= 2: rho(inf, d) <= 1/(2(d+1)).
Running time: the cost needs rho <= l^{-C} (C ~ 2). Here that means d^2 >= 2 l (C - 1/2) ln l, against
2 l (C - o(1)) ln l for OpenAI. Only the constant in d = Theta(sqrt(l log l)) improves, so the time stays
2^{Theta(n (loglog n)^2 / sqrt(log n))} with a slightly smaller constant. NO exponent improvement.

## 2. Barrier A: single-F (OpenAI-type) per-point families (PROVED)
For any F in R[y_1..y_l] of total degree <= d, g_z(t) = chi(z) + F(s (.) (1 - t/tau)) - F(s), with per-point
errors (bad edges only), every direction has rho_k >= (1/2)(1 - S(l-1, d-1)) = 1/4 - O(d/sqrt l), where
S(m,D) = 2^{-m} sum_{i <= (m+D)/2} C(m,i).
Proof. A k-edge is good iff L_k(s_{-k}) = chi_{-k}(s_{-k}), where L_k = (F(s) - F(s^k))/(2 s_k) has
multilinear degree <= d-1. Smolensky's argument (it works over R): on the agreement set A, chi_T = chi * chi_{T^c}
= L_k * chi_{T^c}, so every function on A is a polynomial of degree <= (m + d - 1)/2. Hence |A| <= S * 2^m.  QED
So per-point errors gain only through several "pieces" (section 1). With one global F, d = Omega(sqrt l) is forced.

## 3. Barrier B: reduction to the hypercube Turan problem (PROVED necessary condition)
Psi(t) = sum_{z in Q} P_z(t) over any (d+1)-subcube Q (directions S) has deg <= d and Psi(0) = 2^{d+1}. It vanishes
at tau_k whenever every k-edge of Q is good. So every (d+1)-subcube contains a bad edge, and the union of
the bad sets is a Q_{d+1}-transversal of Q_l: rho_avg >= (1/2) * tau_l(d+1), where tau_l(m) is the minimum
edge density hitting all m-subcubes. (This is G6d's lemma in combinatorial form.)
Consequences (literature quoted from memory, so treat as unverified):
  * d = 1: tau(2) >= 1 - 0.6068 (Balogh-Hu-Lidicky-Liu, C4-free subgraphs of Q_n). So rho_avg >= 0.196 - o(1),
    better than the lemma's 1/8. Proving rho >= 1/4 combinatorially is Erdos' C4 conjecture; our d=1 and
    piecewise constructions reach 1/4.
  * Uniform rho = 2^{-Theta(d)} would give Q_{d+1}-transversals of density 2^{-Theta(d)}. That would settle the
    Alon-Krech-Szabo gap Omega(log d/(d 2^d)) <= 1 - c_d <= O(1/d^2) in the surprising direction. The best
    known constructions (polychromatic colourings, Offner) have density ~4/d^2.
  * Beating sqrt(log n) needs transversals of density <~ l^{-2} with d = o(sqrt(l log l)). Neither the
    layered nor the polychromatic transversals do this.
Algebra is strictly stronger than the combinatorial condition: the combinatorial optimum is smaller in every
case tested. Comb./alg. totals: (3,1) 3/4, (4,1) 8/12, (4,2) 3/4, (5,2) 8/<=14, (5,3) 3/4, (6,2) 22/<=44,
(6,3) 8/<=21. None of the optimal transversals is algebraically feasible (turan.py). The AKS/Offner
polychromatic colour classes (edge colour (#ones before k mod p, #ones after k mod q), p+q = d+2) are
infeasible at (6,2), (7,2) and (7,3) for every class (poly.py). The layered classes are feasible for every
(l,d) tested.

## 4. Derived-instance recursion (PROVED, structural)
For a good edge write g_z - g_{z^k} = 2 chi(z) (1 - t/tau_k) S_{z,k}(t), with deg S <= d-1 and S(0) = 1.
If the (j,k)-square at z is all good, then (S_{z,j} + S_{z^k,j})(tau_k) = 0. So for each j, y -> S_{(y),j} is an
(l-1, d-1) instance in P-form, whose bad k-edges lie in the (j,k)-squares that contain a bad edge:
rho^{(j)}_k <= 2(rho_j + rho_k). Iterating down to d = 0 (rho = 1/2) gives only rho >= 4^{-d}/2 (weaker than the lemma).

## 5. Data (rho_max = max_k bad_k / 2^l; upper bounds exact over Q)
 (l,d) | piecewise | best found (search)            | lemma LB (rho_avg) 2^{-d-1}/(d+1)
 (6,2) | 0.156     | greedy 0.1875; MILP-max: see milpmax_6_2.txt | 0.0417
 (7,2) | 0.164     | -                                 |
 (6,3) | 0.094     | greedy 0.078 (beats piecewise)    | 0.0156
 (7,3) | 0.094     | MILP-max: see milpmax_7_3.txt     |
 (10,2)| 0.166     | ; (10,3) 0.117; (10,4) 0.070; (10,5) 0.044
The layered good sets are not maximal (extend.py), but the extra edges sit next to the short tail blocks.
For d = 2 their number grows only about linearly (16, 18, 20 for l = 6, 7, 8), so the asymptotics do not change visibly.

## 6. Self-review / unproved items
- Upper bounds (constructions) are proved and verified exactly. Lower bounds: only the lemma or Turan
  condition (exponential in d), plus the Smolensky bound for single-F families. NOT proved: rho(inf,d) >= poly(1/d)
  for the general commutative class. NOT proved: any bound for the non-commutative class beyond G6d's lemma.
- Literature constants in section 3 are quoted from memory.
- G6d's note 2(c) ("disjoint subcube blocks, ranks add") holds for bad-edge counts in the commutative class.
  For general non-commutative J, the rank of a matrix is not superadditive over a block-diagonal pattern of
  submatrices, so that claim needs re-checking there.
- The running-time conclusion uses G6d's cost model (rho <= l^{-C} with C ~ 2 and one point per direction).
