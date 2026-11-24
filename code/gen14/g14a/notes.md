# G14a: the per-point bound over arbitrary fields; remaining loopholes of the running-time consequence

Status tags: [L] Lean-checked, [P] paper proof (here), [N] numerical/SAT, [S] sketch, [O] open.

## 1. Any-field theorem  [L]
F any field, J in F[t]^{2^l x 2^l}, deg <= d, J(0)[x,y] = [y = ~x], tau_1..tau_l distinct nonzero in F,
J(tau_k) = K_k + E_k, K_k[x,y] = 0 if x_k != y_k, 2m+d+1 <= l.  Then  sum_k rank E_k >= N_m = sum_x (m+1-|x|)_+.
Proof = G11b with ONE change: only the LEFT vector a must solve a(tau_k)E_k = 0; the right vector is free.
After dividing by t^v, pick x1 with b_{x1}(0) != 0 (|x1| <= m because b_{x1} != 0 has the |x1| roots tau_k, x1_k=1),
c := pref_{x1} e_{~x1}; Phi = sum_x b_x J[x,~x1] pref_{x1} has deg <= 2m+d < l, l roots, Phi(0) = b_{x1}(0) != 0.
Positivity is not used.  File AnyField.lean (`AnyField.core`, `AnyField.sum_rank_ge_N`, `[Field F]`, axioms:
propext, Classical.choice, Quot.sound; compiled with the project's toolchain via `lake env lean`).
Same as the reviewer's (G13a) argument; found independently here and machine-checked.
Needs |F| > l only to have the points; holds over GF(p^k), C, Q_p; over Z/p^k or Z_p by reduction mod p /
passing to Q_p (rank can only drop), provided the tau_k are distinct nonzero mod p (otherwise Lagrange
interpolation itself is not available without p-adic precision loss = working over Q_p).

### 1b. Refinement (factor ~2), any field  [P; mechanism checked on 2 char-2 solutions, verify_S.py]
L := l-1-d. S := {a in A_L : a(tau_k) E_k = 0 all k}. For a in S of t-valuation j, b = a/t^j: b(0) is supported on
|x| <= L-j, and testing with c = pref_{x0} e_{~x0}, |x0| <= j (deg b + deg c + d <= l-1) gives b_{x0}(0) = 0.
So gr_j S embeds in F^{ {x : j < |x| <= L-j} }, dim S <= sum_{j<L/2} (ball(L-j) - ball(j)), and
   sum_k rank E_k >= N_L - dim S  =  N_m + N_{m-1}  (L = 2m),   2 N_m  (L = 2m+1).
(Stronger than the reviewer's sum_{w<=mc} C(l,w)(ma+1-w).)  Only the blocks B01_k = J(tau_k)[x_k=0, y_k=1] are used.

## 2. Small cases (sum_k r_k, r_k = max(rank B01, rank B10))
| (l,d) | N_m | refined 1b | det lemma (any field) | Q optimum | GF(p) odd, commutative class, ALL tau | char 2, general class |
|---|---|---|---|---|---|---|
| (3,0) | 5 | 6 | 12 | 12 (forced) | 12 | 12 |
| (3,1) | 1 | 2 | 3 | 4 (G11c Thm A) | 4 (p=5,7,11,13) | **3** over GF(4), exact: (1,1,1) SAT; total 2 UNSAT/det lemma |
| (4,2) | 1 | 2 | 3 | 4 | 4 (p=5,7,11) | **3** over GF(8): (1,1,1,0) SAT for 5 tau-sets |
| (4,1) | 6 | 7 | 8 | comm. 12; general 10..11 | 12 for every tau class (p=5: 1, 7: 3, 11: 22, 13: 43 classes) | **8** over GF(8): (2,2,2,2) SAT |
| (5,2) | 7 | 8 | 7 (8 with integrality, G11c) | 9..14 | not run | not run |
* No counterexample to the theorem anywhere (there cannot be one, sec. 1).  C: covered by sec. 1, no search.
* Odd p (commutative class, exact lazy hitting-set MILP + rank mod p, comm_gfp.py): identical to Q's generic value.
* CHARACTERISTIC 2 IS DIFFERENT in small cases [N, SAT + independent re-verification of J, gen_char2.py]:
  the minimum equals the determinant-lemma bound l 2^{l-d-1}/(d+1) in all three cases; G11c's Thm A/B'/Cor 2
  (which divide by 2 resp. use char 0) FAIL in char 2.  rho_max(4,1) = 1/8 in char 2 vs 3/16 over Q.
  This is a constant-size effect: sec. 1 holds in char 2, so in the window rho <= l^{-C} the degree threshold
  d >= (1-o(1)) sqrt((2C-1) l ln l) is the same in every characteristic.  No asymptotic gain, no faster mod-p
  counting, no decision algorithm from this route.  [O]: exact char-2 asymptotics for fixed d (irrelevant for time).

## 3. Output-level theorem (closes loophole (b), any field)  [P, not machine-checked, not stress-tested]
Setting: A(t) in F[t]^{n x n} ARBITRARY (no factorisation assumed), entry degree <= d, A(0) = A; labels
alpha, beta: [n] -> [q]^l; A(tau_k) = B_k + E_k, B_k[u,v] = 0 unless alpha(u)_k = beta(v)_k.
Let pi be a matching (partial is enough) in supp A with alpha(u)_k != beta(pi(u))_k for all k (property (iii) of the
paper: true for EVERY nonzero in the main construction), and assume A[G,pi(G)] nonsingular for all G in dom(pi)
(generic entry values on the support).  Then with m = floor((l-d-1)/2), p = Pr[Bin(l,1/2) <= m]:
        sum_k rank E_k  >=  |dom pi| * p^2.
Proof. Pick S_k subset [q] (k = 1..l); w(u) = #{k: alpha(u)_k notin S_k}, w'(v) = #{k: beta(v)_k in S_k};
G = {u in dom pi: w(u) <= m, w'(pi u) <= m}.  Test space: a_u = prod_{k: alpha(u)_k notin S_k}(1-t/tau_k) r_u, deg <= m,
u in G (dim >= |G|).  If sum rank E_k < |G|: a != 0 with a(tau_k)E_k = 0; divide by t; a(0) != 0 on G; some
v0 in pi(G) has (a(0)A)_{v0} != 0; c = prod_{k: beta(v0)_k in S_k}(1-t/tau_k) e_{v0}.  Phi = a A(t) c^T: deg <= l-1,
Phi(tau_k) = 0 (a_u(tau_k) != 0 => alpha(u)_k in S_k, c(tau_k) != 0 => beta(v0)_k notin S_k => B_k[u,v0] = 0),
Phi(0) != 0.  For uniformly random S_k the two bits per coordinate are independent (the symbols differ), so
E|G| = |dom pi| p^2.  QED.
Consequences.
* (b) How L, R' "help": only the used vertices U_c of each cube matter (sum_c |U_c| = nnz(A) <= Delta n, while
  M = (q(q-1))^l ~ Q^2 >= 16 Delta^2 n; the paper takes M = n log^20 n), and L, R' then merge <= Delta coordinates.
  So G6d's charge M rho_k overcounts by >= M/nnz.  But the theorem above bounds the OUTPUT rank directly:
  max_k rank E_out,k >= n p^2 / l.  Redoing G11b 2(b) with (n p^2/l) log n instead of rho M log n:
  p^2 <= l^{-1/2+o(1)}  =>  d >= (1-o(1)) sqrt((1/2) l ln l)  =>  h log(n/h) = Omega(n (loglog n)^2/sqrt(log n)).
  Same order, worse constant.  Rank collapse through L, R' (and cube-dependent J_c, and any non-factorised A(t))
  does NOT beat sqrt(log n).   Caveat: generic entries (or any A with the stated nonsingular minors).
* Cube-restricted version: sum_k rank E_k[U,~U] >= max_s sum_{x in U}(m+1-|x+s|)_+ >= |U| N_m/2^l  [P].

## 4. Loophole (a): degree bound D < n d
Statement: Prop. 3.2 needs only l h > deg_t per A(t); with a bound D the Taylor order is h = D/l + 1.  The theorem
bounds the ENTRY degree d, not D.  A construction with large d (tiny errors) but D << n d would escape.
* Appendix B: D = 4 R d because every contributing cover uses exactly 4R deformed entries (R = O(n loglog n)
  factors), matrix size N = n polylog.  This replaces n by 4R >= n: it is why App. B LOSES a loglog factor.
  Sec. 3 applies with |dom pi| >= 4R (crossing W-entries of one cover): rank >= 4R p(d)^2 / l, h = 4R d/l.
  Same trade-off in the variable n_eff = D/d.  No gain.
* Non-uniform degrees [P]: if deg A(t)[u,v] <= p_u + q_v with p, q >= 0 (D_pot := sum p + sum q = max weight of a
  matching in the degree pattern, Egervary), the proof of sec. 3 goes through with
  G = {u: w(u) <= M1 - p_u, w'(pi u) <= M2 - q_{pi u}}, M1+M2 = l-1.  All but D_pot/x rows (columns)
  have p_u <= x (q_v <= x); with x = 4 D_pot/n at least n/2 matched pairs remain (pi perfect), so
  sum rank >= (n/2) Pr[Bin(l,1/2) <= (l-1)/2 - 4 D_pot/n]^2.  Hence any degree bound certified by row/column
  degree budgets ("deform only n' rows", "few high-degree rows") gives D_pot/n >= c sqrt(l log l): no gain.
* NOT covered [O]: D certified by cancellation or by global cover structure with free-sign potentials (App. B's
  gadget identity is of this kind), i.e. deg per A(t) << d * (number of crossing entries per cover).  A gain needs
  an encoding of #SAT where contributing covers use o(n sqrt(loglog n / log n)... ) deformed entries, e.g. 1_F as a
  product of R = o(n) sparse low-degree factors.  That is a question about probabilistic polynomials/isolation
  (R >= n seems forced by the 2^{-n} error requirement), not about J(t).  No lower bound proved here.

## 5. Loophole (c): other structured decompositions
* Covered by sec. 3 [P]: any polynomial path A(t), any outer maps, cube-dependent J_c, strong or weak niceness.
* Arbitrary per-point block partitions (R_k, C_k), not coordinates of one product labelling [S]: same proof with
  D(u) = {k : R_k(u) != C_k(pi u)}; coordinates outside D(u) add exactly 1 to w+w'.  Pr[u in G] >= Pr[Bin(2|D(u)|,1/2)
  <= |D(u)|-d-1] / O(sqrt N).  On expanding inputs a constant fraction of matched entries is off-block at a constant
  fraction of the points (blocks of size <= n/q), so again d = Omega(sqrt(N' log N')) per effective point; constants
  and the multi-point bookkeeping (h >= D/N) not worked out.
* Outside the theorem [O]: (i) special points that are easy for another reason than block + low rank (bounded
  treewidth/branchwidth supports, the "parity form" and tensor-contraction routes of G4c sec.1, char-2 per = det);
  (ii) non-polynomial paths or interpolation not through point values; (iii) costs below n^{Theta(r)} for rank-r
  perturbations.  The theorem is a statement about "block-by-coordinate + low rank at l points", nothing more.

## 6. Self-review
* Sec. 1: machine-checked.  Sec. 1b, 3, 4(non-uniform): short paper proofs with the same mechanism; NOT in Lean,
  3 and 4 not numerically stress-tested (random instances are vacuous; an adversarial test needs min-rank search).
* Sec. 3 needs nonsingular minors A[G,pi(G)] (generic weights); for specific 0/1 inputs replace by
  "A[G, C_m] has full row rank"; the bound is for worst-case/generic inputs of the framework.
* p^2 instead of p: for arbitrary labellings I could not align S with orientations; constant in d only.
* Running-time statements remain statements about the cost accounting (rank r charged n^{Theta(r)}, Taylor order
  h >= D/l charged (n/h)^{Theta(h)}), now with D possibly smaller than nd only as in sec. 4.
* Char-2 table entries: SAT (CaDiCaL) solutions re-verified by reconstructing J(t) and recomputing ranks with the
  same GF(2^e) arithmetic class (not an independent implementation); lower bounds there are the det lemma.
  (4,1) over GF(8): one tau-set found SAT; odd-p results are for the commutative class only.
* 4(b) [S]-items and all [O]-items are not results.
