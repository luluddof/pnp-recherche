# G22c — structured algebraic families at large N, exact ranks  (working notes)

Tags: [E] exact computation, [P] short proof (unreviewed), [N] numerical, [O] open.
Setting: R = d*F - sum_l kron(G_l, G'_l) (rows (x,y), cols (u,v)); C(x,y,u,v) = sum_l G_l(x,u) G'_l(y,v); p = rank R, q = #terms.

## 1. Evaluator (ev.py)  [E]
* Builds R sparse (scipy kron), splits into connected components of the bipartite row/column support graph, ranks each block:
  exactly over Z (flint fmpz_mat) for blocks <= 120, else mod two ~62-bit primes (2^62-57, 2^61-1); method flags 'Q' / 'p'
  ('p!' would flag disagreement: never occurred).  Mod-prime rank <= rank over Q, so 'p' values are lower bounds that equal
  rank_Q except with negligible probability; all values at N <= 64 below are 'Q' (exact).
* best_scale: min over scalars c of rank(F - cM) via eigenvalues of F^T M per block, best rational candidates re-ranked exactly.
* dense_check: independent full einsum build + flint rank over Q of the whole N^2 x N^2 matrix (no blocks).
* Validation (validate.py): exact_4_5_3 -> 5, exact_m_4_4_4 -> 4, exact_5_6_6 -> 6, exact_5_7_4 -> 7, exact_6_7_7 -> 7
  (all = claimed, both methods).  Cayley-Dickson structure-tensor ansatz (sect. 2a) gives p = 1,2,4,8 at N = 1,2,4,8 with
  q = N (Hurwitz (8,8)@8 rebuilt from the octonion table), dense check agrees; N = 16: 100 both ways; N = 32: 494 both ways.
* Speed: N = 128 (16384^2) in ~2-14 s, N = 256 (65536^2) in 25 s for CD.

## 2. Families and exact results  (CSV: famA_cd.csv, famA_trunc.csv, scan_*.csv; logs out_*.txt)
### 2a. Cayley-Dickson (composition ansatz)  [E; P for the identity]
Ansatz: T(x; a,c) structure tensor of a bilinear map with m outputs; G_x(a,c) = T(x;a,c), G'_x(b,d) = T(x;d,b)/2.
Identity [P]: with A = M^T M (M(x;(l,u)) = L_l(x,u)), rank A <= m, and D = A + A^swap - 2 Id (the polarised composition
defect L_l^T L_l' + L_l'^T L_l - 2 delta I), one has R = (A - D)/2, so p <= m + rank D, q = m.  Composition formula
[N,N,m] (D = 0) => (m,m) feasible.  Hurwitz = CD at N <= 8.
CD at N = 2^n, q = N, scale 1/2 (the best scale by best_scale):
   N    : 8   16    32    64     128     256
   p    : 8   100   494   2176   9122    37348 (mod p)
   p/N^2: .125 .39  .48   .53    .557    .570
   log p/log N: 1.00 1.66 1.79 1.85 1.88 1.90; local slope log2(p(2N)/p(N)): 3.6, 2.30, 2.14, 2.07, 2.03  -> exponent 2.
Truncations (slices x < q, best of scales 1, 2, 1/2, -1 and best_scale): q -> p is NOT monotone and never below the q = N value
except trivially; N = 16: min max(p,q) = 100 at q = 16; N = 32: 494 at q = 32.  Sedenion defect is Theta(N^2).
Local search (signls.py, annealing over all 225 free signs of a Z_2^4 table, 19000 exact evaluations each):
  from CD start: no improvement on 100 (CD is a local minimum);  from random signs: 245.  [E for the values; search heuristic]
### 2b. Group families
* Representation sums (translations, field multiplications x -> ax): C = sum_g rho(g) (x) rho'(g).  [P] rank(F - C) >=
  N^2 - rank_op(C), and for a group sum rank_op(C) = q * dim(invariants of rho (x) rho') when g -> G'_g is itself a representation (true for G' = G and, for abelian
  groups, G' = G^T); this is O(N) for regular/permutation reps of size N:
  p >= N^2 - O(N).  Measured: field mult (all a, G' = G^T): p = N^2 - N + 1 ... e.g. 271@17, 1639@41, 3421@59, 4969@71;
  squares only: 262@17, 1618@41.  Translations (cyclic, Z_2^n, dihedral; G' = G, G^T, right transl.; fixed scales):
  p = N^2 for all N in 4..64 (non-abelian G^T case is not covered by the bound but measured trivial too).  (field N >= 79: best_scale not run, values trivial N^2.)  All dead: exponent 2, p ~ N^2.
* Structure-tensor slices of groups (the composition ansatz with a group table):
  cyclic f=1: p = N^2 - N exactly (N = 3..64); random signs: worse.  Z_3^n: 6, 72, 702, (6561 fixed scale) for N = 3,9,27,81
  = N^2 - N.  Dihedral: 12, 24, 52, 76, 120, 220, 352, 516, 940, 2172 for N = 4..48 (p/N^2 -> 0.94).  Heisenberg: 52@8, 621@27.
  Clifford cocycle on Z_2^n: 4, 40, 96, 704, 1792 for N = 4..64 (p/N^2 = 0.25..0.69).  Reason: anticommutation of all pairs
  of L_a needs exponent-2 groups AND a composition algebra: only N <= 8.
### 2c. Finite-field / Paley: see 2b (multiplicative = representation sum); affine family not run (q = r(r-1)/2 > N anyway).
### 2d. Tournament / half-of-pairs matching shape (famD.py, G' = +-G^T, coefficient 1):
  xor-matchings M_x = {(a,a^x): a<a^x} with signs 1/CD/Clifford, non-wrapping shifts, Paley shifts: p = N^2 - O(N) for all
  (e.g. N = 16: 216..250, N = 32: 906..1013; N = 4 already 14 although (4,4) exists).  The exact small solutions are not
  "a tournament covered by structured matchings".
### 2e. Recursive: analyse53.py / kron53.py show no N -> 2N rule (only the plain Kronecker square).  Baselines:
  direct sum b(2N) <= 4 b(N) from b(8) = 8: b = 32, 128, 512, 2048 at N = 16..128 (= N^2/8).
  Kronecker square of the 8-square identity (C1 (x) C2 only): (p,q) = (960, 64) at N = 64 [E] (pq/N^2 = 15).
### 2f. Local improvement: only the CD sign search (2a): CD is locally optimal at N = 16.

## 3. Table of balanced values b = min max(p,q) found in each family (exact)
 N   | CD   | Clifford | dihedral | cyclic struct | rep sums | tournament | direct sum | N^1.414 | N^2
 16  | 100  | 96       | 220      | 240           | 256      | 216        | 32         | 50      | 256
 32  | 494  | 704      | 940      | 992           | 1024     | 906        | 128        | 134     | 1024
 64  | 2176 | 1792     | (2172@48)| 4032          | 4096     | -          | 512        | 362     | 4096
 128 | 9122 | -        | -        | -             | -        | -          | 2048       | 955     | 16384
 256 | 37348| -        | -        | -             | -        | -          | 8192       | 2545    | 65536
Every family has local exponent -> 2 and is WORSE than the direct-sum baseline (N^2/8) at every N >= 16.

## 4. Verdict
No structured family tested beats exponent 2; none even beats direct sums of Hurwitz blocks (N^2/8) for N >= 16; nothing
approaches sqrt2 (b(64) = 2176 for CD vs N^1.414 = 362 vs direct sum 512).  Best "family" in closed form remains the
composition ansatz (Hurwitz), exact only for N <= 8; beyond, p = rank of the polarised defect + O(N) ~ 0.57 N^2.
Structural lessons [P]: (i) any C that is a group-representation sum has rank_op O(N) => p >= N^2 - O(N);
(ii) structure-tensor ansatz: p <= m + rank(defect), so a family with b = N^{beta} via this ansatz = integer/real
composition formulas [N,N,N^beta] with small defect — the sum-of-squares composition problem (Hurwitz/Yuzvinsky/HWY).
[O] whether integer composition formulas [N,N,o(N^2)] with explicit structure give a usable family (not tested: no tables
of intercalate matrices implemented).

## 5. Self-review
* All N <= 64 values exact over Q via block decomposition; N = 16, 32 CD also by independent dense full-matrix rank.
  N = 128, 256 CD by two 62-bit primes (agree) — lower bounds on rank_Q; since the claim is "bad", lower bounds suffice.
* Scale: only scalar rescalings optimised (plus 4 fixed scales); per-term weights not optimised.  For N > 77 best_scale
  was skipped; those entries could be beaten by a better scale (not by much: they are rep sums, bound (i)).
* Families are a small sample; "no family" is evidence only about these families.
* Two stray error-log files were accidentally created in the repo root by a mis-cd'd background launch and removed at once.
