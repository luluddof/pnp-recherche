# G4c notes: OpenAI's ETH refutation compared with Claude's

## 1. OpenAI's method (from O'Donnell's exposition of the permanent route)
Goal: per(A) for an n x n matrix with O(1) nonzeros per row and column. Sparse #3SAT reduces to this through the
Valiant / Ben-Dor–Halevi reduction.
- Interpolation: P(t) = per(A(t)) with A(0) = A and deg A(t) <= d. Take l special points t = 1..l. At each one,
  compute P, P', ..., P^(h-1) and recover P(0) by Hermite interpolation. This needs l*h > n*d, so h = n*d/l + 1.
- Easy points: A(i) = B_i + E_out. B_i is block-diagonal with blocks of size s = 3n/q, and rank(E_out) <= r.
  Every square sub-permanent costs 2^{O(s)} n^{O(r)} (Barvinok colour-count DP combined with per-block subset DP).
  Derivatives come from the addition formula, which enumerates pairs |S| = |T| = a < h and costs (n/h)^{O(h)}.
- Construction: rows and columns get word labels in [q]^l, with q = log^10 n and q^l ~ sqrt n. Factor A = L J R'
  with J = I_C ⊗ X^{⊗l}. X^{⊗l} is the swap of the ordered label pair on each Boolean cube {0,1}^l.
  "i-nice" (preserves bit i) makes the middle factor block-diagonal by the i-th symbol.
  J(t) = X^{⊗l} + p(Σ(1 - t/j) X_j) - p(Σ X_j), where p has degree d = O(sqrt(l log l)) and p(Σx) = parity
  except on a fraction l^-100. At t = i the X_i term vanishes, so that point is i-nice. The error is
  E = X^{⊗l} - p(ΣX_j), which a Hadamard change of basis makes diagonal with rank = #disagreements.
- Where the sqrt(log n) loss comes from: exponent ≈ h log(n/h) = n (d/l) log(n/h). Exact niceness forces
  d >= l, because the entry (x, x̄) must vanish at l points. The low-rank error allows d ~ sqrt(l log l), the
  width of the binomial (Hamming-weight) window. With l = Θ(log n / log log n) this gives
  n (log log n)^2 / sqrt(log n). The headline n log log n / sqrt(log n) bound comes from a different route
  (tensor contraction, Appendix D of the preprint).
- Limits paper (Dias–Lima, arXiv 2611.00517):
  - Crossovers with hidden constant C: n > 10^32 (C=2), 10^126 (C=3), 10^3517 (C=10).
  - Constants of the tensor-contraction proof: certified time >= 2^n unless log2(n+m) > 1.03e9
    (K0 ≈ 16077, r = 64 groups, ε = 2^-17); at density 4.26, log2 n > 5.3e11, i.e. C ≈ 6e5.
    The interpolation exponent with constant 1 crosses over at n ≈ 10^2390.
  - Clause-product scheme: the vanishing-product lemma shows that at some point at most q_t occurrences
    lose their own variable. Consequences: tw <= (w+1)s - 1 + q, and the cost is >= 2^{Ω(sqrt tw)},
    2^{Ω(sqrt n)} on random 3-CNF, and 2^{Ω(n^{1/3})} on paths. The parity form and the tensor
    contraction escape this, because they break the degree bound (H2) and condition (H3).
  - Open: the true cost lies somewhere between 2^{sqrt n} and 2^{ρ(n)}.

## 2. Same idea? Yes, the same skeleton
| | Claude | OpenAI (permanent / interpolation) |
|---|---|---|
| object crossed | crossing [x=v][y=u] = swap of the 2 diagonals of a grid rectangle | X^{⊗l} = swap of an ordered word pair on a cube |
| directions | 2 (merge rows / merge columns) | l ≈ log n/loglog n (coordinate i, block by i-th symbol) |
| grid | g × g, K = g^2 nodes | q^l words, Q ≈ sqrt n |
| rectangle | 4 corners/2 diagonals | cube {0,1}^l with l=2 is exactly a rectangle |
| expansion | exact signed sum of poly(λ) pieces, each cheap one way | Hermite combination of l points + low-rank E + Taylor residual |
| per-unit cost | ρ = 11/12 (octonions) or 1/sqrt2 (Heisenberg) | d/l ≈ sqrt(log l / l) |
| sqrt source | Hardy–Ramanujan log p(E) ~ sqrt E | binomial window sqrt l |
| recursion | log K doubles per level, e *= ρ | one shot (one "level" with l directions) |
| effective α in e=(log size)^-α | log(1/ρ)/log 2 = 1/2 | log(l/d)/log l -> 1/2 |
| algebra | Clifford (octonions), Heisenberg/CCR Fock space | Pauli X_j,Z_j (Clifford), Hadamard |
Parameter map: OpenAI q ↔ Claude g; OpenAI l ↔ the number of directions D (Claude D = 2); OpenAI i ↔ the
row/column choice; OpenAI h (Taylor order) ↔ the "residual crossing bits"; OpenAI rank(E_out) ↔ Claude's
extra choice bits (log of the number of pieces); limits paper R = d, |X| = l, q_t = h.

## 3. Combination / transfer; explicit model
Model: D directions per level, gain ρ_D per level, and log K multiplied by D per level, so
e = (log K)^{-α} with α = log(1/ρ_D)/log D. Also K <= sqrt(3m), because K^2 λ ≈ 3m and λ >= 1.
Time is 2^{m e}.
(a) OpenAI interpolation replacing Claude's λ^{O(K^2)} branching: not a bottleneck. The overhead only
    forces K <= m^{1/9} instead of m^{1/2}, which changes the constant (factor sqrt(9/2) ≈ 2.1 at α = 1/2)
    and does not change α.
(b) Claude's recursion on OpenAI's Taylor residual (a clique of q = polylog blocks): at best a factor
    (log q)^{-1/2} = (log log n)^{-1/2}. In addition, the residual edges have degree 2^l, which inflates the
    bits to h(log(n/h) + l). No gain of order sqrt(log n).
(c) Hybrid (OpenAI interpolation across D/2 coordinate pairs, Claude's pieces inside each pair):
    α = [½ log 2 + ½ log D' - ½ log log D'] / [log 2 + log D'] -> 1/2. The result is still
    2^{O(n/sqrt log n)} (up to log log factors).
WALL (exact statement):
- Claude framework: some piece needs c_row >= 2λ - O(log λ) (rank of a permutation matrix). With the uniform
  budget c_row + c_col <= 4ρλ this forces ρ >= 1/2. With D directions the same argument gives ρ_D >= 1/D.
  So α <= 1 for every D, and the time is >= 2^{Ω(m/log m)}.
- OpenAI framework (common E): NEW LEMMA (G4c). If deg_t J <= d, J(0) = X^{⊗l}, and J(i) = K_i + E with K_i
  preserving bit i, then for dist(x,y) > d the entry J(t)[x,y] - E[x,y] has > d roots, so it is identically
  0 and E[x,y] = X^{⊗l}[x,y]. Take A = a Hamming ball of radius floor((l-d-1)/2). Then E[A, Ā] is the
  identity, so rank(E) >= |ball| ≈ 2^l exp(-d^2/2l). OpenAI's p achieves exactly 2 × this bound
  (checked for l <= 100; brute-force matrix ranks for l <= 8 in ball_check.py).
  So d >= sqrt(2 l ln(1/δ)) is forced. Since l <= log n/log log n, OpenAI's sqrt(log n) is intrinsic:
  α = 1/2 is exactly optimal in that framework.
- Common wall: a swap over a product structure with <= log n coordinates; each piece is nice in at most a
  1/(#directions) fraction. Passing 2^{n/log n} requires a non-product source of directions, or algebraic
  structure in the constraints themselves.

## 4. Open questions (ranked by closeness to going below 2^{n/log n})
1. (Limits §9) A design on expanders with sqrt n groups, O(sqrt n) width, Taylor order sqrt n and polylog rows
   would give 2^{Õ(sqrt n)}, which is beyond the wall. By the analysis above, it needs non-product structure.
2. (Limits §5.1/§9) True limits of the parity form and of the 2-adic tensor contraction. Both lie outside
   H1–H3 and have no analysed barrier.
3. (Claude App. H) A graded space with log dim V_E ~ E^α, α < 1/2, would give ρ = 2^{α-1} -> 1/2 and time
   2^{m/(log m)^{1-α}}. This is the most concrete question but is capped at 2^{m/log m}.
4. (G4c, new) A per-point error E_i (allowed by Prop. 3.2) might escape my lemma. Even d = O(1) only gives
   2^{n polyloglog/log n}, which is the same wall.
5. Removing log log n (Claude, n-bound); Naslund's P(N) >= N^{1+ε}; constants and crossovers
   (practically irrelevant: n > 10^{10^8} and log2(n+m) > 10^9).
